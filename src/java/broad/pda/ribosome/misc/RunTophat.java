package broad.pda.ribosome.misc;

import java.io.File;
import java.io.IOException;

import broad.core.util.PipelineUtils;

public class RunTophat {

	private static void runTophat(String fq, String saveDir, String bamSaveDir, Runtime run, String jobID) throws IOException, InterruptedException {
		File file=new File(fq); 
		String name=file.getName().split("\\.")[0];
		File outDir=new File(saveDir+"/"+name); outDir.mkdir();
		String output=bamSaveDir+"/"+name+".bsub";
		//String command="bsub -q week -o "+bamSaveDir+"/"+name+".bsub ";
		String command1="~nmcabili/bin/tophat -G /seq/lincRNA/RNASeq/scripture/RefSeqUCSCScripture.gtf -o "+ outDir.getAbsolutePath()+"/ /seq/lincRNA/data/mm9.nonrandom.bowtie "+fq;
		String command2="mv "+outDir.getAbsolutePath()+"/accepted_hits.bam "+bamSaveDir+"/"+name+".bam";
		String command3="java -jar /seq/mgarber/tools/picard-tools-1.66/BuildBamIndex.jar I="+ bamSaveDir+"/"+name+".bam O="+bamSaveDir+"/"+name+".bam.bai";
		
		String[] commands={command1, command2, command3};
		PipelineUtils.bsubProcess(run, jobID, commands, output);
		//run.exec(command);
	}
	
	
	public static void main(String[] args) throws IOException, InterruptedException{
		if(args.length>2){
			File[] fastqFiles=new File(args[0]).listFiles();
			String tophatSaveDir=args[1];
			String bamSaveDir=args[2];
			
			Runtime run=Runtime.getRuntime();
			String jobID="U"+System.nanoTime();
			
			//Run Tophat
			for(int i=0; i<fastqFiles.length; i++){
				File fq=fastqFiles[i];
				runTophat(fq.getAbsolutePath(), tophatSaveDir, bamSaveDir, run, jobID);
			}
			
			PipelineUtils.waitForJobs(jobID, run);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=fastq files \n args[1]=tophat save dir \n args[2]=bam save dir";
	
}
