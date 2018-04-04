package broad.pda.seq.alignment;

import java.io.File;
import java.io.IOException;

public class BowtieJunctionsBatch {
	String queue="priority";
	
	//TODO: Need to add -a if I'm going to use the -m flag
	public BowtieJunctionsBatch(File[] fastqFiles, String junctionsFile, String saveDir) throws IOException{
		Runtime run=java.lang.Runtime.getRuntime();
		for(int i=0; i<fastqFiles.length; i++){
			String command="bsub -q "+queue+" -o "+saveDir+"/"+fastqFiles[i].getName()+".bsub" +" bowtie "+junctionsFile+" "+fastqFiles[i]+" "+saveDir+"/"+fastqFiles[i].getName()+".sam" +" --sam -k 3 -m 1 --un "+saveDir+"/"+fastqFiles[i].getName()+".unaligned.junctions.fq --max "+saveDir+"/"+fastqFiles[i].getName()+".multimappers.junctions.fq";
			run.exec(command);
			System.err.println(command);
		}	
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
		File[] fastqFiles=new File(args[0]).listFiles();
		String junctionsFile=args[1];
		String saveDir=args[2];
		new BowtieJunctionsBatch(fastqFiles, junctionsFile, saveDir);
		}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=fastq files \n args[1]=junctions file \n args[2]=saveDir";
}
