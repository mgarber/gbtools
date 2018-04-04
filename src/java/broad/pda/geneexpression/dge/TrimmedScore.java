package broad.pda.geneexpression.dge;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Iterator;

import org.broad.igv.Globals;

import broad.core.error.ParseException;
import broad.core.hmm.BadModelException;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.annotation.GTFFileParser;
import broad.pda.gene.RefSeqGene;
import broad.pda.gene.RefSeqGeneWithIsoforms;
import broad.pda.seq.alignment.AlignmentUtils;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;

public class TrimmedScore {
	static String usage="Usage: TrimmedScore -task <task name> "+
	"\n\tscore: Computes expression of a given annotation set by trimming  transcripts for optimal coverage \n\t\t-annotations <Specific regions to segment> \n\t\t-alignment <Alignment (mapped to genome)> \n\t\t-window <Window used to score gene. Default is 75bp> \n\t\t-out <Output file [Defaults to stdout]>\n\t\t-quantile <Quantile to use for trimming, 0.25 by default> " +
		"\n";
	
	public static void main (String [] args) throws IOException, ParseException, BadModelException {
		Globals.setHeadless(true);
		ArgumentMap argMap = CLUtil.getParameters(args,usage , "score");
		 if ("score".equalsIgnoreCase(argMap.getTask())) {
			 String alignmentFile = argMap.getMandatory("alignment");
			 String annotationFile = argMap.getMandatory("annotations");
			 int minimumMappingQuality = argMap.containsKey("minMappingQuality")  ?  argMap.getInteger("minMappingQuality") : -1;
			 int window = argMap.containsKey("window") ? argMap.getInteger("window") : 75;
			 double trimQuantile = argMap.containsKey("quantile") ? argMap.getDouble("quantile") : 0.25;
			 
			 BEDFileParser annotationParser =  annotationFile.endsWith(".gtf") || annotationFile.endsWith(".GTF")? new GTFFileParser(annotationFile) : new BEDFileParser(annotationFile);
			 ContinuousDataAlignmentModel libData = AlignmentUtils.loadAlignmentData(alignmentFile, true, minimumMappingQuality);
			 
			 Iterator<String> chrIt = annotationParser.getChromosomeIterator();
			 BufferedWriter bw = argMap.getOutputWriter();
			 int processed=0;
			 while(chrIt.hasNext()) {
				 String chr = chrIt.next();
				 if(libData.hasDataForChromosome(chr)) {
					 Iterator<RefSeqGeneWithIsoforms> annotationIt = annotationParser.getChrTree(chr).valueIterator();
					 System.err.println("Processing " + chr);
					 while(annotationIt.hasNext()) {	
						 processed++;
						 RefSeqGeneWithIsoforms annotation = annotationIt.next();
						 RefSeqGene trimmedAnnotation = libData.trimGeneEnds(annotation, trimQuantile);
						 if(trimmedAnnotation != null) {	
							 double [] scores = libData.scoreGene(trimmedAnnotation);
							 trimmedAnnotation.setName(annotation.getName());
							 trimmedAnnotation.addExtraField(String.valueOf(scores[0])); //FDR
							 trimmedAnnotation.addExtraField(String.valueOf(scores[1])); //Enrichment
							 trimmedAnnotation.addExtraField(String.valueOf(scores[2])); //total reads
							 trimmedAnnotation.addExtraField(String.valueOf(scores[4])); //RPKM
							 trimmedAnnotation.setBedScore(scores[4]);
							 bw.write(trimmedAnnotation.toBED());
						 } else {
							 annotation.addExtraField("0"); //FDR
							 annotation.addExtraField("0"); //Enrichment
							 annotation.addExtraField("0"); //total reads
							 annotation.addExtraField("0"); //RPKM
							 bw.write(annotation.toBED());
						 }
						 
						 bw.newLine();
						 if(processed % 500 == 0) { System.out.println("\t"+processed + " transcripts");}
					 }
				 }
			 }
			 bw.close();
		 }
	}
	
}
