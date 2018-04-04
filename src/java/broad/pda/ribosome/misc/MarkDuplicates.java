package broad.pda.ribosome.misc;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

import broad.core.util.PipelineUtils;

public class MarkDuplicates {

	//Mark duplicates, build new BAM Index, Make TDF
	public MarkDuplicates(Collection<File> files, String saveDir, String genome) throws IOException, InterruptedException{
		
		Runtime run=Runtime.getRuntime();
		String jobID=PipelineUtils.getJobID();
		
		for(File file: files){
			if(file.getName().endsWith(".bam")){
				System.err.println(file.getName());
				
				//Step 1: Mark duplicates
				String saveDuplicates=saveDir+"/"+file.getName().replaceAll(".bam", "")+".noduplicates.bam";
				//String infoDuplicates=saveDir+"/"+file.getName().replaceAll(".bam", "")+".duplicate.info";
				//String cmd1="java -jar -Xmx5000m /seq/mgarber/tools/picard-tools-1.66/MarkDuplicates.jar I="+file.getAbsolutePath()+" O="+saveDuplicates+" M="+infoDuplicates;
				String cmd1="samtools rmdup -s "+file.getAbsolutePath()+" "+saveDuplicates;
				
				
				//Step 2: build BAM index
				String cmd2="java -jar -Xmx5000m /seq/mgarber/tools/picard-tools-1.66/BuildBamIndex.jar I="+saveDuplicates+" O="+saveDuplicates+".bai";
				
				
				//Step 3: Make TDF
				String tdf=saveDir+"/"+file.getName().replaceAll(".bam", "")+".tdf";
				String cmd3="/xchip/igv/tools/igvtools count "+saveDuplicates+" "+tdf+" "+genome;
				
				//Step 4: Run all jobs on a node
				String[] commands={cmd1, cmd2, cmd3};
				String output=saveDir+"/"+file.getName().replaceAll(".bam", "")+".bsub";
				PipelineUtils.bsubProcess(run, jobID, commands, output);
			}
		}
		
	}
	
	
	public static void main(String[] args) throws IOException, InterruptedException{
		if(args.length>2){
			Collection<File> bamFile=getBAM(args[0]);
			String saveDir=args[1];
			String genome=args[2];
			new MarkDuplicates(bamFile, saveDir, genome);
		}
		else{
			System.err.println(usage);
		}
	}
	
	private static Collection<File> getBAM(String string) {
		File[] files=new File(string).listFiles();
		Collection<File> rtrn=new ArrayList<File>();
		
		for(int i=0; i<files.length; i++){
			if(files[i].getName().endsWith(".bam")){rtrn.add(files[i]);}
		}
		
		return rtrn;
	}

	static String usage=" args[0]=BAM files \n args[1]=save directory \n args[2]=genome";
	
}
