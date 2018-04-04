package broad.pda.rap;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

import broad.core.annotation.ShortBED;
import broad.core.annotation.ShortBEDReader;
import broad.core.error.ParseException;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.datastructures.Alignments;
import broad.pda.feature.genome.SimpleGenome;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

public class ComputePromoterBias {

	
	private static String USAGE = "java -jar Rapture.jar -in <example.bam> -sizeFile <sizeFile> [-promoters <other.bed>]";
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException, ParseException {

		ArgumentMap argmap = CLUtil.getParameters(args, USAGE);
		String inFile = argmap.getInput();
		String sizeFile = argmap.getMandatory("sizeFile");
		String maskFileDir = argmap.get("maskFileDir", null);
		File[] maskFiles  = maskFileDir != null ?  new File(maskFileDir).listFiles(): null ;
		int extensionFactor = argmap.getInteger("extensionFactor", 0);
		// default file contains 500bp before and 100bp into each gene starting from its TSS
		String regionFile = argmap.get("promoters", "/seq/lincRNA/data/mm9/mm9_RefSeqGenes.TSS_promoter.bed");  
		ShortBEDReader reader = new ShortBEDReader(regionFile);
		Iterator<ShortBED> itr = reader.getAnnotationList().iterator();
		
		long promoterBases = reader.getBaseCoverage();
		SimpleGenome genome = new SimpleGenome(sizeFile, maskFiles);
		long unmaskedGenomeBases = genome.getNumUnmaskedBases();
		double fractionPromoter = (double) promoterBases / (double) unmaskedGenomeBases;  // technically we should mask the promoters the same way ... but oh well

		GenericAlignmentDataModel data = new GenericAlignmentDataModel(inFile, sizeFile);
		double totalCount = data.getTotalNumberOfReads();
		int promoterCount = 0;
		
		while (itr.hasNext()) {
			Alignments curr = new Alignments(itr.next());
			promoterCount += data.getCountsPerAlignment(curr, extensionFactor);
		}
		
		double score = promoterCount * (1000000.0 / totalCount);
		
	
		double enrichment = promoterCount / totalCount / fractionPromoter;
		
		System.out.println(inFile + "\t" + promoterCount + "\t" + totalCount + "\t" + score + "\t" + enrichment);
	}

}
