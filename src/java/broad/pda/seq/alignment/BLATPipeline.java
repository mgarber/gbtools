package broad.pda.seq.alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeSet;

import broad.core.annotation.PSL;
import broad.core.util.PipelineUtils;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.alignment.sam.SAMRecord;
import broad.pda.seq.fastq.FastqParser;

public class BLATPipeline {
	
	static String path="/seq/mguttman/scripts/GregBLATPipeline/";
	static String blatPath="/seq/mguttman/scripts/BLAT/";
	static int waitTime=60000; //1 minute
	static String queue="week";
	static String saveDir;

	public static void runBLATPipeline(File fastq, String genomeDirectory, String saveDir, int minScore, int minPercentIdentity) throws IOException, InterruptedException{
		new File(saveDir).mkdir();
		BLATPipeline.saveDir=saveDir;
		Runtime run=java.lang.Runtime.getRuntime();
		String faFile=saveDir+"/"+fastq.getName()+".fa";
		
		System.err.println("Converting from fastq to fasta ...");
		
		//Step 1: Convert the fastq to fasta file
		convertFastqToFasta(fastq, faFile);
		
		
		//Step 2: Run BLAT on all reads
		System.err.println("Run BLAT...");
		runBLAT(faFile, genomeDirectory, minScore, minPercentIdentity, run);
		
		//Step 3: Merge all BLAT results
		System.err.println("Merge BLAT File");
		mergeBLATFiles(run); 
		
		//Step 4: Get unique BLAT alignments
		System.err.println("Filter unique and write SAM");
		
		//GetUniquePSLOnly.uniqueToSAM(saveDir+"/R.sort", saveDir+"/"+fastq.getName()+".sam");
		pslToSAM(saveDir+"/R.sort", saveDir+"/"+fastq.getName()+".sam");
	}

	//TODO Classify strata for BLAT alignments
	private static void pslToSAM(String sam, String out) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(sam)));
		FileWriter writer=new FileWriter(out);
		
		int i=0;
		
		String nextLine;
		String name="";
		Collection<PSL> records=new ArrayList<PSL>();
        while ((nextLine = reader.readLine()) != null) {
        	PSL record=new PSL(nextLine);   
        	//SAMRecord record=new SAMRecord(nextLine);
        	if(record.getName().equalsIgnoreCase(name)){records.add(record);}
        	else{ 
        		if(records.size()==1){
        			SAMRecord samRecord=records.iterator().next().toSAM();
        			samRecord.setMappingQuality(255);
        			writer.write(samRecord+"\n");
        		}
        		else{
	        		for(PSL psl: records){
	        			SAMRecord samRecord=psl.toSAM();
	        			samRecord.setWeight((double)1/records.size());
	        			samRecord.setMappingQuality(0);
	        			writer.write(samRecord+"\n");
	        		}
        		}
        		
        		records=new ArrayList<PSL>(); records.add(record);
        	}
        	name=record.getName();
        	//mark the transition in read name
        	i++;
        	if(i%100000 ==0){System.err.println(i);}
        }
		reader.close();		
		writer.close();
	}

	private static void mergeBLATFiles(Runtime run) throws IOException, InterruptedException {
		mergeBLATFiles(new File(saveDir+"/BLAT/").listFiles(), run);
		
	}

	public static void parseBLATFilter(File[] blatFiles, String name, String saveDir) throws IOException, InterruptedException{
		Runtime run=java.lang.Runtime.getRuntime();
		BLATPipeline.saveDir=saveDir;
		
		System.err.println("Merge BLAT File");
		//Merge BLAT files
		mergeBLATFiles(blatFiles, run);
		
		//Step 13: Parse BLAT output
		System.err.println("Parse BLAT output ...");
		parseBLAT(run);

		//Step 15: Convert to SAM file
		System.err.println("Converting to SAM...");
		//convertToSAM(new File(saveDir+"/BlatUnique"), saveDir+"/"+name+".blat.sam");
	}
	

	private static void mergeBLATFiles(File[] blatFiles, Runtime run) throws IOException, InterruptedException {
		concatanate(blatFiles, saveDir+"/R.blat");
		sortBLATFile(saveDir+"/R.blat", run);
	}

	private static void sortBLATFile(String in, Runtime run) throws IOException, InterruptedException {
		String command="sort " +in+" -o "+saveDir+"/R.sort -k 10";
		Process p=run.exec(command);
		p.waitFor();
		writeError(p.getErrorStream());
	}

	private static void concatanate(File[] blatFiles, String save) throws IOException {
		FileWriter writer=new FileWriter(save);

		for(int i=0; i<blatFiles.length; i++){
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(blatFiles[i])));

			String nextLine;
			while ((nextLine = reader.readLine()) != null) {
				PSL psl=PSL.create(nextLine);              	
				if(psl != null) {
					writer.write(psl.toPSL()+"\n");
				}
			}
			reader.close();
		}
		writer.close();

	}

	

	private static void mergeBowtieAndBLAT(Runtime run) throws IOException, InterruptedException {
		String command="perl "+path+"merge_BowtieUnique_and_BlatUnique.pl";
		Process p=run.exec(command);
		p.waitFor();
		writeError(p.getErrorStream());
	}


	private static void parseBLAT(Runtime run) throws IOException, InterruptedException {
		String command="perl "+path+"parse_blat_out.pl  R "+saveDir+"/R.blat "+saveDir+"/R.mdust  35  35 "+saveDir+"/BlatUnique";
		PipelineUtils.bsubProcess(run, command);
	}


	private static void runBLAT(String fasta, String genomeDirectory, int minScore, int minPercentIdentity, Runtime run) throws IOException, InterruptedException {
		//Here we are going to do some fancy work with the farm
		String jobID=BLATAll(genomeDirectory, fasta, minScore, minPercentIdentity, run);
		waitForJobs(jobID);
		//TODO Check whether all jobs completed successfully. If not relaunch
	}
	
	public static void concatanateAndSortBLAT(File[] pslFiles, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		
		
		Map<Integer, Collection<PSL>> map=BEDFileParser.loadPSLBySeqName(pslFiles);
		for(Integer name: map.keySet()){
			Collection<PSL> psls=map.get(name);
			for(PSL psl: psls){writer.write(psl.toPSL()+"\n");}
		}
		
		writer.close();
	}
	
	private static String BLATAll(String genomeDirectory, String fastaFile, int minScore, int minPercentIdentity, Runtime run)throws IOException, InterruptedException{
		new File(saveDir+"/BLAT/").mkdir();
		new File(saveDir+"/bsub/").mkdir();
		
		
		File[] dirs=new File(genomeDirectory).listFiles();
		String UID="U"+System.nanoTime();
		for(int i=0; i<dirs.length; i++){
			if(dirs[i].isDirectory()){
			String chr="chr"+dirs[i].getName();
			String command="bsub -o "+saveDir+"/bsub/"+chr+".bsub"+" -q "+queue+" -J "+UID+" "+blatPath+"blat "+ dirs[i].getAbsolutePath()+"/"+chr+".fa"+" "+fastaFile+" "+saveDir+"/BLAT/"+chr+".blat -minScore="+minScore+" -minIdentity="+minPercentIdentity+" -mask=lower -repeats=lower -out=pslx";
			PipelineUtils.bsubProcess(run, command);

			}
		}
		System.out.println("BLATing Sequences");
		return UID;
	}
	
	private static void waitForJobs(String jobID)throws IOException, InterruptedException{
		int count=5;
		
		while(count>1){
			Thread.sleep(waitTime);
			Process proc = Runtime.getRuntime().exec("bjobs -J "+jobID);
			BufferedReader out = new BufferedReader(new InputStreamReader(proc.getInputStream()));
			count=parseReply(out);
			out.close();
			System.err.println("checking status "+count+" jobs still running");
		}
		System.err.println("done");
		
		checkForFailures(jobID);
	}
	
	private static void checkForFailures(String jobID) throws IOException {
		Process primer3Proc = Runtime.getRuntime().exec("bjobs -a -J "+jobID);
		BufferedReader primer3StdOut = new BufferedReader(new InputStreamReader(primer3Proc.getInputStream()));
		String line = null;
		int i=0;
		while((line = primer3StdOut.readLine()) != null) {
			//split the line
			String[] tokens=line.split(" ");
			String completionStatus=tokens[2];
			if(completionStatus.equalsIgnoreCase("EXIT")){
				System.err.println("WARN: Blat Job "+tokens[0]+" FAILED");
				//To avoid proceeding when the jobs failed and from having to go through every job to determine if its all perfect
				//if it gets here I'm going to throw an exception
				throw new IllegalArgumentException("BLAT jobs failed");
			}
			//System.err.println(i);
			//System.err.println(line);
			i++;
		}
	}

	private static int parseReply(BufferedReader out) throws IOException {
		String line = null;
		int i=0;
		//System.err.println("Parsing LSF output, ready?" + out.ready() + " line: " + out.readLine());
		while((line = out.readLine()) != null) {
			//System.err.println(i);
			//System.err.println(line);
			i++;
		}
		return i;
	}


	private static void runMDust(Runtime run) throws IOException, InterruptedException {
		String file="R";
		System.err.println(new File(file).getAbsolutePath());
		String command="/seq/mguttman/scripts/GregBLATPipeline/mdust/mdust "+file+" >t";
		System.err.println(command);
		Process p=run.exec(command);
		p.waitFor();
		writeError(p.getErrorStream());
		//writeInputStream(p.getInputStream(), "R.mdust");
	}

	private static void writeInputStream(InputStream stream, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		BufferedReader reader=	new BufferedReader(new InputStreamReader(stream));
		String rtrn="";
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			System.err.println(nextLine);
			writer.write(nextLine+"\n");
			rtrn+=nextLine;
		}
		writer.close();
	}

	private static void convertToSAM(File file, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String name=tokens[0];
			String chr=tokens[1];
			Collection<Alignments> blocks=getBlocks(tokens[2], chr);
			String sequence=tokens[3];
			RefSeqGene gene=new RefSeqGene(blocks);
			gene.setSequence(sequence);
			writer.write(gene.toSAM()+"\n");
		}
		reader.close();
		writer.close();
	}


	private static Collection<Alignments> getBlocks(String string, String chr) {
		Collection<Alignments> rtrn=new TreeSet<Alignments>();
		
		String[] tokens=string.split(",");
		
		for(int i=0; i<tokens.length; i++){
			String block=tokens[i];
			int start=new Integer(block.split("-")[0].trim());
			int end=new Integer(block.split("-")[1].trim());
			Alignments exon=new Alignments(chr, start, end);
			rtrn.add(exon);
		}
		
		
		return rtrn;
	}


	private static void getReadsToBLAT(String faFile, Runtime run) throws IOException, InterruptedException {
		String command="perl "+path +"make_R.pl "+faFile+" -single";
		Process p=run.exec(command);
		p.waitFor();
		writeError(p.getErrorStream());
	}


	private static void renameX(Runtime run) throws IOException, InterruptedException {
		//rename
		String command="cp X_sorted X";
		Process p=run.exec(command);
		p.waitFor();
		writeError(p.getErrorStream());
	}

	private static void renameY(Runtime run) throws IOException, InterruptedException {
		//rename
		String command="cp Y_sorted Y";
		Process p=run.exec(command);
		p.waitFor();
		writeError(p.getErrorStream());
	}
	

	private static void mergeGNUandTNUandCNU(Runtime run) throws IOException, InterruptedException {
		String command="perl "+path+" merge_GNU_and_TNU_and_CNU.pl";
		Process p=run.exec(command);
		p.waitFor();
		writeError(p.getErrorStream());
	}


	private static void mergeTUAndGU(Runtime run) throws IOException, InterruptedException {
		String command="perl "+path+" merge_TU_and_GU.pl";
		Process p=run.exec(command);
		p.waitFor();
		writeError(p.getErrorStream());
	}


	private static void mapToGenomicCoordinates(String transcriptomeInfo, Runtime run) throws IOException, InterruptedException {
		String command="perl "+path+" map_transcriptome_bowtie_alignments_to_genomic_coords.pl "+transcriptomeInfo+" -single";
		Process p=run.exec(command);
		p.waitFor();
		writeError(p.getErrorStream());
		
		command="perl "+path+" merge_TU_and_TNU.pl";
		p=run.exec(command);
		p.waitFor();
		writeError(p.getErrorStream());
	}


	private static void sortY(Runtime run) throws IOException, InterruptedException {
		//String command="perl "+path+" sort_Y.pl";
		String command="sort Y -T -n -o Y_sorted";
		Process p=run.exec(command);
		p.waitFor();
		writeError(p.getErrorStream());
		
		//rename
		command="mv Y_sorted Y";
		p=run.exec(command);
		p.waitFor();
		writeError(p.getErrorStream());
	}


	private static void runBowtieToTranscriptome(String faFile, String transcriptomeRef, Runtime run) throws IOException, InterruptedException {
		String command="bowtie -a -f "+transcriptomeRef+" "+faFile+" Y -v 3 --suppress 6,7,8";
		Process p=run.exec(command);
		p.waitFor();
		writeError(p.getErrorStream());
	}


	private static void filterXForConsistentMappers(Runtime run) throws IOException, InterruptedException {
		String command="perl "+path+" filter_X_for_consistent_mappers.pl -single";
		Process p=run.exec(command);
		p.waitFor();
		writeError(p.getErrorStream());
		writeError(p.getInputStream());
		System.err.println(p.exitValue());
	}


	private static void sortX(Runtime run) throws IOException, InterruptedException {
		//String command="perl "+path+" sort_X.pl";
		String command="sort X -T -n -o X_sorted";
		Process p=run.exec(command);
		int completed=p.waitFor();
		System.err.println("Completed "+completed);
		writeError(p.getErrorStream());
		
		//rename
		command="mv X_sorted X";
		p=run.exec(command);
		p.waitFor();
		writeError(p.getErrorStream());
	}


	private static void runBowtieToGenome(String faFile, String genomeRef, Runtime run) throws IOException, InterruptedException {
		String command="bowtie -a -f "+genomeRef+" "+faFile+" X -v 3 --suppress 6,7,8";
		//bowtie-0.12.1/bowtie -a -f /PATH2BOWTIEINDEX/mm9_genes_ucsc_refseq temp3.fa -v 3 --suppress 6,7,8 -p 3 > Y & 
		Process p=run.exec(command);
		p.waitFor();
		writeError(p.getErrorStream());
	}

	
	private static void writeError(InputStream errorStream) throws IOException {
		BufferedReader reader=	new BufferedReader(new InputStreamReader(errorStream));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			System.err.println(nextLine);
		}
		System.err.println();
	}

	public static void convertFastqToGregFasta(File seq, String save) throws IOException{
		FastqParser fastq=new FastqParser(seq);
		fastq.convertToNumberedFasta(save);
	}
	
	public static void convertFastqToFasta(File seq, String save) throws IOException{
		FastqParser fastq=new FastqParser(seq);
		fastq.convertToFasta(save);
	}
	
	/*public static void main(String[] args)throws IOException, InterruptedException{
		if(args.length>2){
			File[] blatFiles=new File(args[0]).listFiles();
			String name=args[1];
			String saveDir=args[2];
			parseBLATFilter(blatFiles, name, saveDir);
		}
		else{System.err.println(usage);}
	}

	private static String usage=" args[0]=BLAT files \n args[1]=name \n args[2]=saveDir";*/
	
	public static void main(String[] args)throws IOException, InterruptedException{
		if(args.length>4){
			File fastq=new File(args[0]);
			String genomeDir=args[1];
			String save=args[2];
			int minScore=new Integer(args[3]);
			int minIdentity=new Integer(args[4]);
			runBLATPipeline(fastq, genomeDir, save, minScore, minIdentity);
		}
		else{System.err.println(usage);}
	}

	private static String usage=" args[0]=fastq file \n args[1]=genome dir \n args[2]=saveDir \n args[3]=minScore \n args[4]=min percent identity";
		
}
