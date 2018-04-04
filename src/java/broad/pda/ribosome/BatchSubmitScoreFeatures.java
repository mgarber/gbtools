package broad.pda.ribosome;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;

public class BatchSubmitScoreFeatures {

	
	//batch submit the feature scoring for all of the BAM files
	/*public BatchSubmitScoreFeatures(File[] bamFiles, String genes, String outDir) throws IOException{
		
		Runtime run=Runtime.getRuntime();
		
		for(int i=0; i<bamFiles.length; i++){
			File bamFile=bamFiles[i];
			if(bamFile.getAbsolutePath().endsWith(".bam")){
				String output=outDir+"/"+bamFile.getName()+".bsub";
				String save=outDir+"/"+bamFile.getName()+".scores";
				String bsub="bsub -q week -R rusage[mem=6] -o "+output;
				bsub+=" java -jar -Xmx6000m /seq/lincRNA/scripts/Ribosome/ComputeRibosomeOccupancyByFeature.jar "+bamFile+" /seq/genome/mouse/mouse_Mm9/sizes "+genes+" "+save;
				System.err.println(bsub);
				run.exec(bsub);
			}
		}
	}*/
	
	
	public static void batchSubmitScoreFeatures(File[] bamFiles, File expressionBAM, String genes, String outDir) throws IOException, InterruptedException{
		Runtime run=Runtime.getRuntime();
		
		Map<String, Collection<RefSeqGene>> genesByChr=BEDFileParser.loadDataByChr(new File(genes));
		
		for(int i=0; i<bamFiles.length; i++){
			File bamFile=bamFiles[i];
			if(bamFile.getAbsolutePath().endsWith(".bam")){
				for(String chr: genesByChr.keySet()){
					String output=outDir+"/"+chr+"."+bamFile.getName()+".bsub";
					String save=outDir+"/"+chr+"."+bamFile.getName()+".scores";
					String bsub="bsub -q hour -R rusage[mem=6] -o "+output;
					bsub+=" java -jar -Xmx6000m /seq/lincRNA/scripts/Ribosome/ComputeRibosomeOccupancyByFeature.jar "+bamFile+" "+expressionBAM+" /seq/genome/mouse/mouse_Mm9/ "+genes+" "+save+" "+chr;
					System.err.println(bsub);
					Process p=run.exec(bsub);
					p.waitFor();
				}
			}
		}
		
	}
	
	
	public static void BatchSubmitScoreWindows(File[] bamFiles, File expressionBAM, String genes, String outDir, int windowSize) throws IOException, InterruptedException{
		Runtime run=Runtime.getRuntime();
		
		Map<String, Collection<RefSeqGene>> genesByChr=BEDFileParser.loadDataByChr(new File(genes));
		
		for(int i=0; i<bamFiles.length; i++){
			File bamFile=bamFiles[i];
			if(bamFile.getAbsolutePath().endsWith(".bam")){
				for(String chr: genesByChr.keySet()){
					String output=outDir+"/"+chr+"."+bamFile.getName()+".bsub";
					String save=outDir+"/"+chr+"."+bamFile.getName()+".scores";
					String bsub="bsub -q week -R rusage[mem=6] -o "+output;
					bsub+=" java -jar -Xmx6000m /seq/lincRNA/scripts/Ribosome/ComputeRibosomeOccupancyInWindows.jar "+bamFile+" "+expressionBAM+" /seq/genome/mouse/mouse_Mm9/sizes "+genes+" "+save+" "+windowSize+" "+chr;
					System.err.println(bsub);
					Process p=run.exec(bsub);
					p.waitFor();
				}
			}
		}
		
	}
	
	
	public static void main(String[] args) throws IOException, InterruptedException{
		if(args.length>3){
			File[] bamFiles=new File(args[0]).listFiles();
			File expressionBAM=new File(args[1]);
			String genes=args[2];
			String outDir=args[3];
			if(args.length>4){
				int windowSize=new Integer(args[4]);
				BatchSubmitScoreWindows(bamFiles, expressionBAM, genes, outDir, windowSize);
			}
			else{
				batchSubmitScoreFeatures(bamFiles, expressionBAM, genes, outDir);
			}
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam files \n args[1]=expression BAM \n args[2]=genes \n args[3]=out dir \n args[4]=window size (optional)";
}
