package broad.pda.ribosome;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;

public class BatchSubmitStartCodon {
	
	static String script="/seq/lincRNA/scripts/Ribosome/FinalScripts/StartCodonAnalysis.jar";

	public BatchSubmitStartCodon(File[] bamFiles, String genes, String outDir) throws IOException, InterruptedException{
		Runtime run=Runtime.getRuntime();
		Map<String, Collection<RefSeqGene>> genesByChr=BEDFileParser.loadDataByChr(new File(genes));
		for(int i=0; i<bamFiles.length; i++){
			File bamFile=bamFiles[i];
			if(bamFile.getAbsolutePath().endsWith(".bam")){
				for(String chr: genesByChr.keySet()){
					String name=bamFile.getName().split("\\.")[0];
					String output=outDir+"/"+chr+"."+name+".bsub";
					String save=outDir+"/"+chr+"."+name+".scores";
					String bsub="bsub -q week -R rusage[mem=5] -o "+output;
					bsub+=" java -jar -Xmx5000m "+script+" "+bamFile+" "+genes+" /seq/genome/mouse/mouse_Mm9/ "+save+" "+chr;
					System.err.println(bsub);
					Process p=run.exec(bsub);
					p.waitFor();
				}
			}
		}	
	}
	
	public static void main(String[] args) throws IOException, InterruptedException{
		if(args.length>2){
			File[] bamFiles=new File(args[0]).listFiles();
			String genes=args[1];
			String outDir=args[2];
			new BatchSubmitStartCodon(bamFiles, genes, outDir);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam files \n args[1]=genes \n args[2]=out dir";
	
}
