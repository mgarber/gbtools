package broad.pda.miseq;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import broad.core.util.PipelineUtils;
import broad.pda.seq.alignment.Pair;

public class AlignWithBowtie {

	public AlignWithBowtie(File[] fastqFiles, String saveDir, String bowtieIndex, String queue, String options) throws IOException, InterruptedException{
		Map<String, Pair<File>> pairs=getPairs(fastqFiles);
		
		alignWithBowtie(pairs, bowtieIndex, saveDir, queue, options);
	}
	
	
	private static void alignWithBowtie(Map<String, Pair<File>> fastqFiles, String bowtieIndex, String bowtieSaveDir, String queue, String options) throws IOException, InterruptedException {
		Runtime run=Runtime.getRuntime();
		String jobID=PipelineUtils.getJobID();
		for(String name: fastqFiles.keySet()){
			Pair<File> files=fastqFiles.get(name);
			String save=bowtieSaveDir+"/"+name+".sam";
			String output=bowtieSaveDir+"/"+name+".bsub";
			String command="/seq/mgarber/tools/bowtie-0.12.7/bowtie --sam "+bowtieIndex +" "+options;
			if (files.getValue2() != null) {
				command = command +" -1 "+files.getValue1().getAbsolutePath() +" -2 "+files.getValue2().getAbsolutePath();
			} else {
				command = command + " " + files.getValue1().getAbsolutePath();
			}
			command = command + " " + save;
			System.err.println(command);
			PipelineUtils.bsubProcess(run, jobID, command, output, queue);
		}
		PipelineUtils.waitForJobs(jobID, run);
	}
	
	public static Map<String, Pair<File>> getPairs(File[] fastqFiles) {
		Map<String, Pair<File>> rtrn=new TreeMap<String, Pair<File>>();
		
		for(int i=0; i<fastqFiles.length; i++){
			String name=fastqFiles[i].getName();
			String read="";
			String baseName=name.split("\\.")[0];
			if(baseName.endsWith("_1")){
				read="1";
			}
			else if(baseName.endsWith("_2")){
				read="2";
			}
			baseName=baseName.substring(0, baseName.length()-2);
			System.err.println(fastqFiles[i].getName()+" "+baseName);
			
			Pair<File> pair=new Pair<File>();
			if(rtrn.containsKey(baseName)){
				pair=rtrn.get(baseName);
			}
			
			if(read.equalsIgnoreCase("1")){
				pair.setValue1(fastqFiles[i]);
			}
			else if(read.equalsIgnoreCase("2")){
				pair.setValue2(fastqFiles[i]);
			}
			else{
				System.err.println("NO READS");
			}
			rtrn.put(baseName, pair);
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException, InterruptedException{
		if(args.length>3){
			File[] fastq=new File(args[0]).listFiles();
			String tmp=args[1];
			String save=args[2];
			String queue=args[3];
			String options = "";
			if (args.length > 4)
				options =args[4];

			new AlignWithBowtie(fastq, tmp, save, queue, options);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=fastq files \n args[1]=save dir \n args[2]=BAM index \n args[3]=queue";
	
}
