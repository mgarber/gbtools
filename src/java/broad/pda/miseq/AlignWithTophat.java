package broad.pda.miseq;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import broad.core.util.PipelineUtils;
import broad.pda.seq.alignment.Pair;

public class AlignWithTophat {

	public AlignWithTophat(File[] fastqFiles, String tempSaveDir, String bamSaveDir, String queue) throws IOException, InterruptedException{
		Runtime run=Runtime.getRuntime();
		String jobID=PipelineUtils.getJobID();
		
		Map<String, Pair<File>> pairs=getPairs(fastqFiles);
		
		for(String pair: pairs.keySet()){
			System.err.println(pair+" "+pairs.get(pair).getValue1()+" "+pairs.get(pair).getValue2());
			Pair<File> files=pairs.get(pair);
			runTophat(files, pair, tempSaveDir, bamSaveDir, run, jobID, queue);
		}
		PipelineUtils.waitForJobs(jobID, run);
		
		
		for(String name: pairs.keySet()){
			File outDir=new File(tempSaveDir+"/"+name); outDir.mkdir();
			String command2="mv "+outDir.getAbsolutePath()+"/accepted_hits.bam "+bamSaveDir+"/"+name+".bam";
			Process p=run.exec(command2);
			p.waitFor();
		}
		
		for(String name: pairs.keySet()){
			String command3="java -jar /seq/mgarber/tools/picard-tools-1.53/BuildBamIndex.jar I="+ bamSaveDir+"/"+name+".bam O="+bamSaveDir+"/"+name+".bam.bai";
			String output=tempSaveDir+"/"+name+".index.bsub";
			PipelineUtils.bsubProcess(run, jobID, command3, output, "hour");
		}
		PipelineUtils.waitForJobs(jobID, run);
		
		
		for(String name: pairs.keySet()){
			String command4="/xchip/igv/tools/igvtools count "+bamSaveDir+"/"+name+".bam "+bamSaveDir+"/"+name+".tdf mm9";
			String output=tempSaveDir+"/"+name+".TDF.bsub";
			PipelineUtils.bsubProcess(run, jobID, command4, output, "hour");
		}
		PipelineUtils.waitForJobs(jobID, run);		
	}
	
	
	public static void AlignWithBowtie(File[] fastqFiles, String saveDir, String bowtieIndex, String queue) throws IOException, InterruptedException{
		Runtime run=Runtime.getRuntime();
		String jobID=PipelineUtils.getJobID();
		
		Map<String, Pair<File>> pairs=getPairs(fastqFiles);
		
		
		alignWithBowtie(pairs, bowtieIndex, saveDir, queue);
		
		
	}
	
	private static Map<String, Pair<File>> getPairs(File[] fastqFiles) {
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


	private void runTophat(Pair<File> pair, String name, String saveDir, String bamSaveDir, Runtime run, String jobID, String queue) throws IOException, InterruptedException {
		File outDir=new File(saveDir+"/"+name); outDir.mkdir();
		String output=bamSaveDir+"/"+name+".bsub";
		String command1="/seq/mgarber/tools/tophat-1.3.2.Linux_x86_64/tophat -G /seq/lincRNA/RNASeq/scripture/RefSeqUCSCScripture.gtf -o "+ outDir.getAbsolutePath()+"/ /seq/lincRNA/data/mm9.nonrandom.bowtie "+pair.getValue1().getAbsolutePath()+" "+pair.getValue2().getAbsolutePath();
		PipelineUtils.bsubProcess(run, jobID, command1, output, queue);
	}
	
	private static void alignWithBowtie(Map<String, Pair<File>> fastqFiles, String bowtieIndex, String bowtieSaveDir, String queue) throws IOException, InterruptedException {
		Runtime run=Runtime.getRuntime();
		String jobID=PipelineUtils.getJobID();
		for(String name: fastqFiles.keySet()){
			Pair<File> files=fastqFiles.get(name);
			String save=bowtieSaveDir+"/"+name+".sam";
			String output=bowtieSaveDir+"/"+name+".bsub";
			String command="/seq/mgarber/tools/bowtie-0.12.7/bowtie --sam "+bowtieIndex+" "+files.getValue1().getAbsolutePath()+" "+files.getValue2().getAbsolutePath()+" "+save;
			System.err.println(command);
			PipelineUtils.bsubProcess(run, jobID, command, output, queue);
		}
		PipelineUtils.waitForJobs(jobID, run);
	}
	
	
	private void runTophat(String fq, String saveDir, String bamSaveDir, Runtime run, String jobID, String queue) throws IOException, InterruptedException {
		File file=new File(fq); 
		String name=file.getName().split("\\.")[0];
		File outDir=new File(saveDir+"/"+name); outDir.mkdir();
		String output=bamSaveDir+"/"+name+".bsub";
		//String command="bsub -q week -o "+bamSaveDir+"/"+name+".bsub ";
		String command1="/seq/mgarber/tools/tophat-1.3.2.Linux_x86_64/tophat -G /seq/lincRNA/RNASeq/scripture/RefSeqUCSCScripture.gtf -o "+ outDir.getAbsolutePath()+"/ /seq/lincRNA/data/mm9.nonrandom.bowtie "+fq;
		String command2="mv "+outDir.getAbsolutePath()+"/accepted_hits.bam "+bamSaveDir+"/"+name+".bam";
		String command3="java -jar /seq/mgarber/tools/picard-tools-1.53/BuildBamIndex.jar I="+ bamSaveDir+"/"+name+".bam O="+bamSaveDir+"/"+name+".bam.bai";
		String command4="/xchip/igv/tools/igvtools count "+bamSaveDir+"/"+name+".bam "+bamSaveDir+"/"+name+".tdf mm9";
		
		String[] commands={command1, command2, command3, command4};
		PipelineUtils.bsubProcess(run, jobID, command1, output, queue);
		//run.exec(command);
	}
	
	public static void main(String[] args) throws IOException, InterruptedException{
		if(args.length>3){
			File[] fastq=new File(args[0]).listFiles();
			String tmp=args[1];
			String save=args[2];
			String queue=args[3];
			new AlignWithTophat(fastq, tmp, save, queue);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=fastq files \n args[1]=temp save dir \n args[2]=BAM save dir \n args[3]=queue";
	
}
