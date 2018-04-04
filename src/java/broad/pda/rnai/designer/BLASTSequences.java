package broad.pda.rnai.designer;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;


public class BLASTSequences {
	static double alpha=.01;
	static String queue="priority";

	public BLASTSequences(String fastaFile, String genomeDirectory, String saveDir)throws IOException{
		Map<String, RNAiScore> map=parseRNAiScore(new File(fastaFile));
		Map<String, Set> kmerBLASTPositions=BLAST(genomeDirectory, fastaFile, saveDir);
		write(saveDir+"/blastResults.txt", kmerBLASTPositions, map);
	}
	
	
	private Map parseRNAiScore(File file)throws IOException{
		 BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		    
		 Map rtrn=new TreeMap();
        String nextLine;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
           String[] tokens=nextLine.split("\t");
			try{
           Alignments position=new Alignments(tokens[0]);
			String kmer=tokens[2];
			double score=new Double(tokens[3]);
			int blast=new Integer(tokens[6]);
			RNAiScore RNAiScore=new RNAiScore(kmer, position, score, blast);
			rtrn.put(kmer, RNAiScore);
			}catch(Exception ex){}
        }
        
        
        reader.close();
        return rtrn;
		
	}
	
	private void write(String save, Map<String, Set> map, Map<String, RNAiScore> scores)throws IOException{
		FileWriter writer=new FileWriter(save);
		for(String kmer: map.keySet()){
			Set<Alignments> set=map.get(kmer);
			RNAiScore rnai=scores.get(kmer);
			rnai.numBlastHits=set.size();
			writer.write(rnai+"\n");
			//for(Alignments align: set){writer.write("\t"+align.toUCSC());}
			//writer.write("\n");
		}
		writer.close();
	}
	
	public static Map BLAST(Map<String, Alignments> kmerPositions, String genomeDirectory, String fastaFile, String saveDir)throws IOException{
		String jobID=BLASTAll(genomeDirectory, fastaFile, saveDir);
		waitForJobs(jobID);
		
		// % identity	 alignment length	 mismatches	 gap openings e-value
		Map<Alignments, double[]> blastOut=parseOutputFiles(saveDir+"/BLAST/"); //make sure it matches the target and only the target
		//Set<String> filtered=filterUnique(blastOut, kmerPositions, genes);
		Map<String, Set> kmerBLASTPositions=findAllThePlacesThatAKMerLands(blastOut);
		return kmerBLASTPositions;
	}
	
	public static Map BLAST(String genomeDirectory, String fastaFile, String saveDir, boolean b)throws IOException{
		String jobID=BLASTAll(genomeDirectory, fastaFile, saveDir);
		waitForJobs(jobID);
		
		// % identity	 alignment length	 mismatches	 gap openings e-value
		Map<Alignments, double[]> blastOut=parseOutputFiles(saveDir+"/BLAST/"); //make sure it matches the target and only the target
		//Set<String> filtered=filterUnique(blastOut, kmerPositions, genes);
		Map<String, Set> kmerBLASTPositions=findAllThePlacesThatAKMerLands(blastOut);
		return kmerBLASTPositions;
	}
	
	public static Map BLAST(String genomeDirectory, String fastaFile, String saveDir, boolean b, int wordSize)throws IOException{
		String jobID=BLASTAll(genomeDirectory, fastaFile, saveDir, wordSize);
		waitForJobs(jobID);
		
		// % identity	 alignment length	 mismatches	 gap openings e-value
		Map<Alignments, double[]> blastOut=parseOutputFiles(saveDir+"/BLAST/"); //make sure it matches the target and only the target
		//Set<String> filtered=filterUnique(blastOut, kmerPositions, genes);
		Map<String, Set> kmerBLASTPositions=findAllThePlacesThatAKMerLands(blastOut);
		return kmerBLASTPositions;
	}
	
	public static Map BLAST(String genomeDirectory, String scoreFile, String saveDir)throws IOException{
		File fastaFile=convert(scoreFile, saveDir);
		
		String jobID=BLASTAll(genomeDirectory, fastaFile.getAbsolutePath(), saveDir);
		waitForJobs(jobID);
		
		// % identity	 alignment length	 mismatches	 gap openings e-value
		Map<Alignments, double[]> blastOut=parseOutputFiles(saveDir+"/BLAST/"); //make sure it matches the target and only the target
		//Set<String> filtered=filterUnique(blastOut, kmerPositions, genes);
		Map<String, Set> kmerBLASTPositions=findAllThePlacesThatAKMerLands(blastOut);
		return kmerBLASTPositions;
	}
	
	
	private static File convert(String scoreFile, String saveDir)throws IOException{
		
		String save=(saveDir+"/kmers.fa"); 
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(scoreFile)));
		FileWriter writer=new FileWriter(save);
		 
         String nextLine;
          while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
            String[] tokens=nextLine.split("\t");
            writer.write(">"+tokens[2]+"\n"+tokens[2]+"\n");					
				
         }
         writer.close();
         
         reader.close();
         return new File(save);
		
	}
	
	private void write(String saveDir, Set<String> filtered, Map<String, Object> scores)throws IOException{
		FileWriter writer=new FileWriter(saveDir+"/scores.txt");
		for(String kmer: filtered){
			writer.write(kmer+"\t"+scores.get(kmer)+"\n");
		}
		writer.close();
	}
	
	private Set filterUnique(Map<Alignments, double[]> blastOut, Alignments align){
		Set rtrn=new TreeSet();
		
		Set good=new TreeSet();
		Set bad=new TreeSet();
		for(Alignments temp: blastOut.keySet()){
			//see if passes criteria
			double[] params=blastOut.get(temp);
			boolean passes=highScoringBLASTHit(params);
			int goodHits=0;
			int badHits=0;
			int geneHits=0;
			if(passes){
				//if so does it overlap alignment
				if(temp.overlapsAtAll(align) || align.overlapsAtAll(temp)){good.add(temp.getName());}
				else{bad.add(temp.getName());}
			}
			
			for(Object name: good){
				if(!bad.contains(name)){rtrn.add(name);}
			}
			
		}
		return rtrn;
	}
	
	private static Map<String, Set> findAllThePlacesThatAKMerLands(Map<Alignments, double[]> blastOut){
		Map rtrn=new TreeMap();
		
		//kmers and all the places it aligns (in addition to the place from which it comes)
		for(Alignments align: blastOut.keySet()){
			double[] params=blastOut.get(align);
			boolean passes=highScoringBLASTHit(params);
			if(passes){
				Set set=new TreeSet();
				if(rtrn.containsKey(align.getName())){set=(Set)rtrn.get(align.getName());}
				set.add(align);
				rtrn.put(align.getName(), set);
			}
		}
		return rtrn;
	}
	
	private static Set filterUnique(Map<Alignments, double[]> blastOut, Map<String, Alignments> kmerPosition, Set<RefSeqGene> genes){
		Set rtrn=new TreeSet();
		//kmers and all the places it aligns (in addition to the place from which it comes)
		
		Set good=new TreeSet();
		Set bad=new TreeSet();
		for(Alignments temp: blastOut.keySet()){
			Alignments align=kmerPosition.get(temp.getName());
			//see if passes criteria
			double[] params=blastOut.get(temp);
			boolean passes=highScoringBLASTHit(params);
			int goodHits=0;
			int genomeHits=0;
			int geneHits=0;
			if(passes){
				//if so does it overlap alignment
				if(temp.overlapsAtAll(align) || align.overlapsAtAll(temp)){good.add(temp.getName()); goodHits++;}
				else if(overlapsGene(temp, genes)){bad.add(temp.getName()); geneHits++;}
				else{genomeHits++;}
			}
			
			for(Object name: good){
				if(!bad.contains(name)){rtrn.add(name);}
			}
			
		}
		return rtrn;
	}
	
	
	private static boolean overlapsGene(Alignments align, Set<RefSeqGene> genes){
		boolean rtrn=false;
		for(RefSeqGene gene: genes){
			if(gene.getAlignment().overlapsAtAll(align) || align.overlapsAtAll(gene.getAlignment())){
				//check if exons overlap
				Alignments[] exons=gene.getExons();
				for(int i=0; i<exons.length; i++){
					if(align.overlapsAtAll(exons[i])|| exons[i].overlapsAtAll(align)){return true;}
				}
			}
		}
		return rtrn;
	}
	
	private static boolean highScoringBLASTHit(double[] params){
		// 0=identity	 1=alignment length	 2=mismatches	 3=gap openings 4=e-value
		if(params[4]<.01 && params[2]<3 && params[1]>19){return true;}
		return false;
	}
	
	
	private static Map parseOutputFiles(String blastDir)throws IOException{
		Map rtrn=new TreeMap();
		File[] files=new File(blastDir).listFiles();
		for(int i=0; i<files.length; i++){
			rtrn.putAll(parse(files[i]));
		}
		return rtrn;
	}
	
	private static Map parse(File file)throws IOException{
		Map rtrn=new TreeMap();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	if(!nextLine.startsWith("#")){
				String[] tokens=nextLine.split("\t");
				if(tokens.length>10){
					int start=new Integer(tokens[8]);
					int end=new Integer(tokens[9]);
					String key=tokens[1]+"\t"+tokens[8]+"\t"+tokens[9]+"\t"+tokens[0];
					if(start>end){key=tokens[1]+"\t"+tokens[9]+"\t"+tokens[8]+"\t"+tokens[0];}
					
					double[] vals={new Double(tokens[2]), new Double(tokens[3]), new Double(tokens[4]), new Double(tokens[5]), new Double(tokens[10])};
					Alignments align=new Alignments(tokens[1]+":"+tokens[9]+"-"+tokens[8]);
					align.setName(tokens[0]);
					rtrn.put(align, vals);
				}
			}
        }		
		return rtrn;
	}
	
	private static void waitForJobs(String jobID)throws IOException{
		int count=5;
		
		while(count>1){
			Process primer3Proc = Runtime.getRuntime().exec("bjobs -J "+jobID);
			BufferedWriter primer3StdIn = new BufferedWriter(new OutputStreamWriter(primer3Proc.getOutputStream()));
			BufferedReader primer3StdOut = new BufferedReader(new InputStreamReader(primer3Proc.getInputStream()));
			BufferedReader primer3StdErr = new BufferedReader(new InputStreamReader(primer3Proc.getErrorStream()));
			count=parseReply(primer3StdOut);
		}
		System.err.println("done");
	}
	
	
	private static int parseReply(BufferedReader primer3StdOut) throws IOException {
		String line = null;
		String [] lineInfo = null;
		int i=0;
		while((line = primer3StdOut.readLine()) != null) {
			//System.err.println(i);
			//System.err.println(line);
			i++;
		}
		return i;
	}
	
	private static String BLASTAll(String genomeDirectory, String fastaFile, String saveDir)throws IOException{
		String junkFile=saveDir+"blast.bsub";
		saveDir+="/BLAST/";
		new File(saveDir).mkdir();
		Runtime run=java.lang.Runtime.getRuntime();
		File[] dirs=new File(genomeDirectory).listFiles();
		String UID="U"+System.nanoTime();
		for(int i=0; i<dirs.length; i++){
			String chr="chr"+dirs[i].getName();
			String command="bsub -o "+junkFile+" -q "+queue+" -J "+UID+" /broad/tools/bin/blastall -p blastn -i "+ fastaFile+" -d "+dirs[i].getAbsolutePath()+"/"+chr+".fa -e "+ alpha+ " -m 9 -o "+saveDir+"/"+chr+".blast";
			run.exec(command);
			//System.err.println(command);
		}
		System.out.println("BLASTing Sequences");
		return UID;
	}
	
	private static String BLASTAll(String genomeDirectory, String fastaFile, String saveDir, int wordSize)throws IOException{
		String junkFile=saveDir+"blast.bsub";
		saveDir+="/BLAST/";
		new File(saveDir).mkdir();
		Runtime run=java.lang.Runtime.getRuntime();
		File[] dirs=new File(genomeDirectory).listFiles();
		String UID="U"+System.nanoTime();
		for(int i=0; i<dirs.length; i++){
			String chr="chr"+dirs[i].getName();
			String command="bsub -o "+junkFile+" -q "+queue+" -J "+UID+" /broad/tools/bin/blastall -p blastn -i "+ fastaFile+" -d "+dirs[i].getAbsolutePath()+"/"+chr+".fa -e "+ alpha+ " -m 9 -o "+saveDir+"/"+chr+".blast "+"-W "+wordSize;
			run.exec(command);
			//System.err.println(command);
		}
		System.out.println("BLASTing Sequences");
		return UID;
	}
	
	public static void writeFasta(Map<String, Object> sequenceScores, String save)throws IOException{
		FileWriter writer=new FileWriter(save);
		int i=1;
		for(String seq: sequenceScores.keySet()){
			writer.write(">"+seq+"\n");
			writer.write(seq+"\n");
			i++;
		}
		
		/*for(String seq: sequenceScores.keySet()){
			writer.write(">"+seq+"\n");
			writer.write(ComputeMIRScore.reverseComplement(seq)+"\n");
			i++;
		}*/
		writer.close();
	}
	
	public static void writeFasta(Set<String> sequences, String save)throws IOException{
		FileWriter writer=new FileWriter(save);
		int i=1;
		for(String seq: sequences){
			writer.write(">"+seq+"\n");
			writer.write(seq+"\n");
			i++;
		}
		
		/*for(String seq: sequenceScores.keySet()){
			writer.write(">"+seq+"\n");
			writer.write(ComputeMIRScore.reverseComplement(seq)+"\n");
			i++;
		}*/
		writer.close();
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>3){
			String fastaFile=args[0];
			String genomeDirectory=args[1];
			String save=args[2];
			new BLASTSequences(fastaFile, genomeDirectory, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=fasta file \n args[1]=genome Directory \n args[2]=save file";
}

