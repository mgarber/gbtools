package broad.pda.ribosome;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;

import broad.core.util.PipelineUtils;

public class CleanFastq {

	public CleanFastq(File[] files, String trimmedDir) throws IOException, InterruptedException{
		Runtime run=Runtime.getRuntime();
		String jobID="U"+System.nanoTime();
		
		//Take raw fastq files and filter and trim according to Ingolia
		Collection<String> fqTrimmed=filterAndTrim(files, trimmedDir, run);
					
	}

	private Collection<String> filterAndTrim(File[] fqFiles, String trimmedDir, Runtime run) throws IOException, InterruptedException {
		//go through each file and trim and filter it
		Collection<String> rtrn=new ArrayList<String>();
		
		for(int i=0; i<fqFiles.length; i++){
			File fq=fqFiles[i];
			System.err.println(fq);
			String save=trimmedDir+"/"+fq.getName();
			rtrn.add(save);
			
			//zcat XXX.fastq.gz |  |  | 
			//step 1
			//String cmd1="/seq/lincRNA/scripts/fastx_toolkit/fastq_illumina_filter --keep N -v -i "+fq+" -o temp.fq";
			//Process p=run.exec(cmd1);
			//p.waitFor();
			
			String cmd2="/seq/lincRNA/scripts/fastx_toolkit/fastx_clipper -Q33 -a CTGTAGGCACCATCAAT -l 25 -c -n -v -i "+fq.getAbsolutePath()+" -o temp2.fq";
			Process p=run.exec(cmd2);
			p.waitFor();
			
			String cmd3="/seq/lincRNA/scripts/fastx_toolkit/fastx_trimmer -Q33 -f 2 -i temp2.fq -o "+save;
			System.err.println(cmd3);
			p=run.exec(cmd3);
			p.waitFor();
			
			String cmd4="rm temp2.fq";
			p=run.exec(cmd4);
			p.waitFor();
		}
		
		return rtrn;
	}

	private void runTophat(String fq, String saveDir, String bamSaveDir, Runtime run, String jobID) throws IOException, InterruptedException {
		File file=new File(fq); 
		String name=file.getName().split("\\.")[0];
		File outDir=new File(saveDir+"/"+name); outDir.mkdir();
		String output=bamSaveDir+"/"+name+".bsub";
		//String command="bsub -q week -o "+bamSaveDir+"/"+name+".bsub ";
		String command1="/seq/mgarber/tools/tophat-1.3.2.Linux_x86_64/tophat -G /seq/lincRNA/RNASeq/scripture/RefSeqUCSCScripture.gtf -o "+ outDir.getAbsolutePath()+"/ /seq/lincRNA/data/mm9.nonrandom.bowtie "+fq;
		String command2="mv "+outDir.getAbsolutePath()+"/accepted_hits.bam "+bamSaveDir+"/"+name+".bam";
		String command3="java -jar /seq/mgarber/tools/picard-tools-1.53/BuildBamIndex.jar I="+ bamSaveDir+"/"+name+".bam O="+bamSaveDir+"/"+name+".bam.bai";
		
		String[] commands={command1, command2, command3};
		PipelineUtils.bsubProcess(run, jobID, commands, output);
		//run.exec(command);
	}

	private void readAndWrite(String input, String output) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		FileWriter writer=new FileWriter(output);
		
		//Keep only last
    	String nextLine;
    	String lastName="";;
        int i=0;
        while ((nextLine = reader.readLine()) != null) {
        	String[] tokens=nextLine.split("\t");
        	String name=tokens[0];
        	String seq=tokens[1];
        	String qual=tokens[2];
        	if(!name.equalsIgnoreCase(lastName)){
        		writer.write("@"+name+"\n"+seq+"\n"+"+"+name+"\n"+qual+"\n");
        	}
        	lastName=name;
        	i++;
        	if(i%1000000 ==0){System.err.println(i);}
        }
        reader.close();
        writer.close();
	}

	public static void main(String[] args)throws IOException, InterruptedException{
		if(args.length>1){
			File[] fastqFiles=new File(args[0]).listFiles();
			String trimmedDir=args[1];
			new CleanFastq(fastqFiles, trimmedDir);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=fastq Files \n args[1]=trimmed fastq directory";


}
