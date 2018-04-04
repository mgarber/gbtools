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

public class ConvertAlignToFastQ {

	public ConvertAlignToFastQ(File[] files, String fastqSaveDir, String tophatSaveDir, String bamSaveDir, String trimmedDir) throws IOException, InterruptedException{
		Runtime run=Runtime.getRuntime();
		String jobID="U"+System.nanoTime();
		
		Collection<String> fqFiles=new ArrayList<String>();
		
		//Convert to Fastq
		for(int i=0; i<files.length; i++){
			String input=files[i].getAbsolutePath();
			String output=fastqSaveDir+"/"+files[i].getName()+".fq";
			readAndWrite(input, output);
			fqFiles.add(output);
		}
		
		
		
		//Take raw fastq files and filter and trim according to Ingolia
		Collection<String> fqTrimmed=filterAndTrim(fqFiles, trimmedDir, run);
				
		
		//Run Tophat
		for(String fq: fqTrimmed){
			runTophat(fq, tophatSaveDir, bamSaveDir, run, jobID);
		}
		
		PipelineUtils.waitForJobs(jobID, run);
		
	}

	private Collection<String> filterAndTrim(Collection<String> fqFiles, String trimmedDir, Runtime run) throws IOException, InterruptedException {
		//go through each file and trim and filter it
		Collection<String> rtrn=new ArrayList<String>();
		
		for(String fq: fqFiles){
			String save=trimmedDir+"/"+new File(fq).getName();
			rtrn.add(save);
			
			//zcat XXX.fastq.gz |  |  | 
			//step 1
			//String cmd1="/seq/lincRNA/scripts/fastx_toolkit/fastq_illumina_filter --keep N -v -i "+fq+" -o temp.fq";
			//Process p=run.exec(cmd1);
			//p.waitFor();
			
			String cmd2="/seq/lincRNA/scripts/fastx_toolkit/fastx_clipper -Q33 -a CTGTAGGCACCATCAAT -l 25 -n -v -i "+fq+" -o temp2.fq";
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
		if(args.length>4){
			File[] input=new File(args[0]).listFiles();
			String fastqOutput=args[1];
			String topHatOutput=args[2];
			String bamSaveDir=args[3];
			String trimmedDir=args[4];
			new ConvertAlignToFastQ(input, fastqOutput, topHatOutput, bamSaveDir, trimmedDir);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=files \n args[1]=fastqOutputDir \n args[2]=tophatOutputDir \n args[3]=directory to save BAM files \n args[4]=trimmed fastq directory";
}
