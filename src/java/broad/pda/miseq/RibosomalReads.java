package broad.pda.miseq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.util.PipelineUtils;

public class RibosomalReads {

	String queue="hour";
	
	public RibosomalReads(File[] fastqFiles, String ribosomalBowtieIndex, String bowtieSaveDir, String save) throws IOException, InterruptedException{
		alignWithBowtie(fastqFiles, ribosomalBowtieIndex, bowtieSaveDir);
		parseResults(bowtieSaveDir, save);
	}

	private static void parseResults(String bowtieSaveDir, String save) throws IOException {
		File[] files=new File(bowtieSaveDir).listFiles();
		
		Map<String, Map<String, Integer>> counts=new TreeMap<String, Map<String, Integer>>();
		
		//Go through each SAM file and count
		for(int i=0; i<files.length; i++){
			if(files[i].getName().endsWith("sam")){
				counts.put(files[i].getName(),countSAM(files[i]));
			}
		}
		
		//Make report table
		write(save, counts);	
	}

	private static Map<String, Integer> countSAM(File file) throws IOException {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		 String nextLine;
		 int total=0;
		 while ((nextLine = reader.readLine()) != null) {
			if(!nextLine.startsWith("@")){
			 String[] tokens=nextLine.split("\t");
			String strand=tokens[1];
			String align=tokens[2];
			if(align.equalsIgnoreCase("*")){align="unaligned";}
			String name=align;
			if(strand.equalsIgnoreCase("0")){name=align+"_sense";}
			else if(strand.equalsIgnoreCase("16")){name=align+"_antisense";}
			
			int counter=0;
			if(rtrn.containsKey(name)){
				counter=rtrn.get(name);
			}
			counter++;
			rtrn.put(name, counter);
			
			total++;
			}
		 }
		 reader.close();
		
		//add total
		 rtrn.put("total", total);
		 
		return rtrn;
	}

	private static void write(String save, Map<String, Map<String, Integer>> counts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<String> entries=new TreeSet<String>();
		
		writer.write("Gene");
		for(String file: counts.keySet()){
			writer.write("\t"+file);
			Map<String, Integer> map=counts.get(file);
			entries.addAll(map.keySet());
		}
		writer.write("\n");
		
		for(String gene: entries){
			writer.write(gene);
			for(String file: counts.keySet()){
				Integer count=counts.get(file).get(gene);
				writer.write("\t"+count);
			}
			writer.write("\n");
		}
		
		writer.close();
	}

	private void alignWithBowtie(File[] fastqFiles,	String ribosomalBowtieIndex, String bowtieSaveDir) throws IOException, InterruptedException {
		Runtime run=Runtime.getRuntime();
		String jobID=PipelineUtils.getJobID();
		for(int i=0; i<fastqFiles.length; i++){
			File fastq=fastqFiles[i];
			String save=bowtieSaveDir+"/"+fastq.getName()+".sam";
			String output=bowtieSaveDir+"/"+fastq.getName()+".bsub";
			String command="/seq/mgarber/tools/bowtie-0.12.7/bowtie --sam "+ribosomalBowtieIndex+" "+fastq+" "+save;
			System.err.println(command);
			PipelineUtils.bsubProcess(run, jobID, command, output, queue);
		}
		PipelineUtils.waitForJobs(jobID, run);
	}
	
	public static void main(String[] args)throws IOException, InterruptedException{
		if(args.length>3){
			File[] files=new File(args[0]).listFiles();
			String index=args[1];
			String saveDir=args[2];
			String save=args[3];
			new RibosomalReads(files, index, saveDir, save);
		}
		else{
			System.err.println(usage);
		}
	}
	
	static String usage=" args[0]=fastq files \n args[1]=bowtie index \n args[2]=bowtie save directory \n args[3]=save";
	
}
