package broad.pda.seq.alignment;


import java.io.File;
import java.io.IOException;
import java.util.Map;

import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;

public class TranscriptomeToGenomeLocation {

	public TranscriptomeToGenomeLocation(File sam, Map<String, RefSeqGene> geneMap, String save, Map<String, Integer> chromosomeSizes) throws IOException{
		System.err.println("started");
		RNASeqAlignmentPipeline.transcriptomeToGenome(sam.getAbsolutePath(), geneMap, save, chromosomeSizes);
	}
	

	public static void main(String[] args)throws IOException{
		if(args.length!=3){
			File sam=new File(args[0]);
			Map<String, RefSeqGene> geneMap=BEDFileParser.loadDataByName(new File(args[1]));
			//System.err.println(geneMap.keySet());
			String save=args[2];
			String sizeFile = args[3];
			RNASeqAlignmentPipeline.transcriptomeToGenome(sam.getAbsolutePath(), geneMap, save, BEDFileParser.loadChrSizes(sizeFile));
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=sam file \n args[1]=gene map \n args[2]=save\n args[3]=sizeFile";
	
}
