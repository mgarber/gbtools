package broad.pda.ribosome.misc;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import broad.pda.annotation.BEDFileParser;

public class BatchSubmitWindows {
	
	String script="/seq/lincRNA/scripts/Ribosome/FinalScripts/ComputeRibosomeOccupancyInWindows.jar";

	public BatchSubmitWindows(File[] bamFiles, File expressionFile,	String sizeFile, File[] genes, String saveDir, String windowSize) throws IOException, InterruptedException {
		Runtime run=Runtime.getRuntime();
		Map<String, Integer> chrSizes=BEDFileParser.loadChrSizes(sizeFile);
		
		for(int i=0; i<bamFiles.length; i++){
			if(bamFiles[i].getName().endsWith(".bam") && !bamFiles[i].getName().equalsIgnoreCase(expressionFile.getName())){
				String name=bamFiles[i].getName().split("\\.")[0];
				for(String chr: chrSizes.keySet()){
					for(int k=0; k<genes.length; k++){
						String geneName=genes[k].getName().split("\\.")[0];
						String save=saveDir+"/"+chr+"."+geneName+"."+name+".w"+windowSize;
						String cmd="bsub -q week -o "+save+".bsub -R rusage[mem=8]";
						cmd+=" java -jar -Xmx8000m "+script+" "+bamFiles[i].getAbsolutePath()+" "+expressionFile.getAbsolutePath()+" "+sizeFile+" "+genes[k].getAbsolutePath()+" "+save+" "+windowSize+" "+chr;
						System.err.println(cmd);
						Process p=run.exec(cmd);
						p.waitFor();
					}
				}
			}
		}
		
	}

	public static void main(String[] args)throws IOException, InterruptedException{
		if(args.length>5){
			File[] bamFiles=new File(args[0]).listFiles();
			File expressionFile=new File(args[1]);
			String sizeFile=args[2];
			File[] genes=new File(args[3]).listFiles();
			String saveDir=args[4];
			String windowSize=args[5];
			new BatchSubmitWindows(bamFiles, expressionFile, sizeFile, genes, saveDir, windowSize);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=BAM Files \n args[1]=expression file \n args[2]=size file \n args[3]=features \n args[4]=saveDir \n args[5]=window size";
	
	//args[0]=BAM file \n argsp[1]=expression BAM \n args[2]=size file \n args[3]=genes (full BED) \n args[4]=save \n args[5]=window size \n args[6]=chr";
		
		
}
