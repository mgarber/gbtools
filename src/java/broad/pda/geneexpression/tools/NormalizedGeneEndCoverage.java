package broad.pda.geneexpression.tools;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.annotation.GTFFileParser;
import broad.pda.gene.RefSeqGene;
import broad.pda.gene.RefSeqGeneWithIsoforms;
import broad.pda.geneexpression.dge.DGE;
import broad.pda.seq.alignment.AlignmentUtils;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.AlignmentDataModelStats;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;
/**
 * @author skadri
 *
 */
public class NormalizedGeneEndCoverage {

	
	static final String usage = "Usage: NormalizedGeneCoverage -task <task name> "+
			"\n\tcoverage: Computes the normalized coverage along the normalized length of the gene" +
			"\n\t\t-annotations <Specific regions to segment [BED by default]> "+
			"\n\t\t-alignment <Alignment (mapped to genome)> "+
			"\n\t\t-bins <Number of bins [100 by default]>"+
			"\n\t\t-windowSize <Step size/window size [By default: GeneLength/#bins]"+
			"\n\t\t-out <Output Filename>"+
			"\n\t\t-stranded";
	
	static Logger logger = Logger.getLogger(NormalizedGeneEndCoverage.class.getName());
	int numBins;
	double[] coverage;
	
	public NormalizedGeneEndCoverage(String[] args)throws IOException{
	
		Globals.setHeadless(true);
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"coverage");
		
		/*
		 * Read the annotation file
		 */
		String annotationsFile = argMap.getMandatory("annotations");
		/*
		 * Check the format of the annotation files and call the GTF or BED parser accordingly
		 */
		BEDFileParser annotations =  annotationsFile.endsWith(".gtf") || annotationsFile.endsWith(".GTF")? new GTFFileParser(annotationsFile) : new BEDFileParser(annotationsFile);
		
		numBins = argMap.isPresent("bins")? argMap.getInteger("bins") : 100;
		coverage = new double[numBins];
		/*
		 * Read the name of the alignment file
		 */
		String alignmentFile = argMap.getMandatory("alignment");
		
		boolean isStranded = argMap.isPresent("stranded");
		boolean isSecondRead = argMap.isPresent("isSecondRead");
		
		/*
		 * Iterate through list of annotations. Iterator is over chromosome names.
		 */
		Iterator<String> chrom_iter = annotations.getChromosomeIterator();
		
		double[] scores;
		int counter = 0;
		
		if(!isStranded){	
		/*
		 * Initialize the data model using the alignment file
		 * @param: <alignment flieName> <load_chromosome_stats> <minMappingQuality> <remove_PCR_duplicates> <weigh reads by NH flag>
		*/
		ContinuousDataAlignmentModel libDataModel = 	AlignmentUtils.loadAlignmentData(alignmentFile,false,5,false,true);
		/*
		 * For each chromosome
		 */
		while(chrom_iter.hasNext()) {
			String chr = chrom_iter.next();
			/*
			 * If the alignment data has data from that chromosome
			 */
			if(libDataModel.hasDataForChromosome(chr)){
				logger.info("Processing " + chr);
				/*
				 * Obtain an iterator over the interval tree values built from the annotations (genes) on that particular chromosome.
				 * IntervalTree of RefSeqGeneWithIsoforms
				 */
				Iterator<RefSeqGeneWithIsoforms> annotation_iter = annotations.getChrTree(chr).valueIterator();
				/*
				 * Parse annotation tree
				 * While there is an annotated RefSeqGeneWithIsoforms in Interval tree to analyze
				 * Thus, for each gene
				 */
				while(annotation_iter.hasNext()){
					/*
					 * annotation = current node of RefSeqGeneWithIsoforms in interval tree
					 */
					RefSeqGeneWithIsoforms annotation_with_isoforms = annotation_iter.next();
					Iterator<RefSeqGene> isoform_iter = annotation_with_isoforms.getAllIsoforms().iterator();

					while(isoform_iter.hasNext()){

						double max = 0.0;
						RefSeqGene annotation = isoform_iter.next();
						System.out.print(annotation.getName()+"\t");
						double[] c = new double[numBins];
						int annotationLength = annotation.getTranscriptLength();
						//System.out.println(Math.round((float)((annotationLength)/(float)numBins)));
						int step = argMap.isPresent("windowSize")? argMap.getInteger("windowSize") : (int)((float)(annotationLength)/(float)numBins);
						//Step = gene length/ # bins
						int threshold = ((step)*numBins);
						//System.out.println(" threshold: "+(step)*numBins);
						//System.out.print(annotation.getName()+" "+annotation.getSize()+" "+annotation.getTranscriptLength());
						if(annotationLength<threshold || threshold<numBins){
							//DONT PROCESS
							//System.out.println(" does not pass.");
						}
						else{
							//System.out.println(" passes");
							//From gene end, for each bin of size step
							if(annotation.getOrientation().equals("+")){
								for(int i=0;i<numBins;i++){
									int j=i+1;
									RefSeqGene subannotation = null;
									int relativeStart = annotation.getTranscriptLength()-(j*step);
									int relativeEnd = annotation.getTranscriptLength()-(i*step);
									if (relativeStart< 0){
										relativeStart = 0;
									}
									
									subannotation = annotation.trim(relativeStart, relativeEnd);
									//System.out.println("Start: "+relativeStart+" End: "+relativeEnd +" Length: "+subannotation.getSize()+" index = "+(numBins-i-1));
									scores = libDataModel.scoreGene(subannotation);
									//COVERAGE
									c[numBins-i-1] = scores[3];
									if(scores[3]>max)
										max = scores[3];
									//System.out.println("Step = "+step);
									//System.out.println("Start = "+relativeStart+" End = "+relativeEnd);
								}
								
							}
							else{
								for(int i=0;i<numBins;i++){
									int j=i+1;
									RefSeqGene subannotation = null;
									int relativeStart = i*step;
									int relativeEnd = j*step;
									if(relativeEnd>annotation.getTranscriptLength()){
										relativeEnd = annotation.getTranscriptLength();
									}
									subannotation = annotation.trim(relativeStart, relativeEnd);
									//System.out.println(annotation.getName()+" Start: "+relativeStart+" End: "+relativeEnd +" Length: "+subannotation.getSize()+" index = "+(numBins-i-1));
									//System.out.println(annotation.getName()+" "+annotation.getSize()+" Start: "+relativeStart+" End: "+relativeEnd +" index = "+(numBins-i-1));
									scores = libDataModel.scoreGene(subannotation);
									c[numBins-i-1] = scores[3];
									if(scores[3]>max)
										max = scores[3];
								}
								//System.out.println("Step: "+step+" "+annotation.getName()+": Start = "+relativeStart+" End = "+relativeEnd+" AnnotationLength = "+subannotation.getTranscriptLength());
								
								//System.out.println(coverage[i]);	
							}
							if(max>0){
								for(int i=0;i<numBins;i++){
									c[i] = c[i]/max;
									coverage[i] += c[i];
									System.out.print(c[i]+"\t");
								}
								System.out.println();
								counter++;
							}
						}
					}
				}
			}
		}
		}
		else{
			String sizes = argMap.get("sizeFile");
			AlignmentDataModel alignmentsP=new GenericAlignmentDataModel(alignmentFile, sizes, false, 5,false,true);
			AlignmentDataModelStats alignmentDataP = new AlignmentDataModelStats(alignmentsP);
			ContinuousDataAlignmentModel libDataModelP = new ContinuousDataAlignmentModel(alignmentDataP);
			alignmentsP.setPositiveStranded();
			
			AlignmentDataModel alignmentsN=new GenericAlignmentDataModel(alignmentFile, sizes, false, 5,false,true);
			AlignmentDataModelStats alignmentDataN = new AlignmentDataModelStats(alignmentsN);
			ContinuousDataAlignmentModel libDataModelN = new ContinuousDataAlignmentModel(alignmentDataN);
			alignmentsN.setNegativeStranded();

			if(isSecondRead){
				alignmentsP.setSecondRead();
				alignmentsN.setSecondRead();
			}
			/*
			 * For each chromosome
			 */
			while(chrom_iter.hasNext()) {
				String chr = chrom_iter.next();
				/*
				 * If the alignment data has data from that chromosome
				 */
				if(libDataModelP.hasDataForChromosome(chr) || libDataModelN.hasDataForChromosome(chr)){
					logger.info("Processing " + chr);
					/*
					 * Obtain an iterator over the interval tree values built from the annotations (genes) on that particular chromosome.
					 * IntervalTree of RefSeqGeneWithIsoforms
					 */
					Iterator<RefSeqGeneWithIsoforms> annotation_iter = annotations.getChrTree(chr).valueIterator();
					/*
					 * Parse annotation tree
					 * While there is an annotated RefSeqGeneWithIsoforms in Interval tree to analyze
					 * Thus, for each gene
					 */
					while(annotation_iter.hasNext()){
						/*
						 * annotation = current node of RefSeqGeneWithIsoforms in interval tree
						 */
						RefSeqGeneWithIsoforms annotation_with_isoforms = annotation_iter.next();
						Iterator<RefSeqGene> isoform_iter = annotation_with_isoforms.getAllIsoforms().iterator();

						while(isoform_iter.hasNext()){

							double max = 0.0;
							RefSeqGene annotation = isoform_iter.next();
							System.out.print(annotation.getName()+"\t");
							double[] c = new double[numBins];
							int annotationLength = annotation.getTranscriptLength();
							//System.out.println(Math.round((float)((annotationLength)/(float)numBins)));
							int step = argMap.isPresent("windowSize")? argMap.getInteger("windowSize") : (int)((float)(annotationLength)/(float)numBins);
							//Step = gene length/ # bins
							int threshold = ((step)*numBins);
							//System.out.println(" threshold: "+(step)*numBins);
							//System.out.print(annotation.getName()+" "+annotation.getSize()+" "+annotation.getTranscriptLength());
							if(annotationLength<threshold || threshold<numBins){
								//DONT PROCESS
								//System.out.println(" does not pass.");
							}
							else{
								//System.out.println(" passes");
								//From gene end, for each bin of size step
								if(annotation.getOrientation().equals("+")){
									for(int i=0;i<numBins;i++){
										int j=i+1;
										RefSeqGene subannotation = null;
										int relativeStart = annotation.getTranscriptLength()-(j*step);
										int relativeEnd = annotation.getTranscriptLength()-(i*step);
										if (relativeStart< 0){
											relativeStart = 0;
										}
										
										subannotation = annotation.trim(relativeStart, relativeEnd);
										//System.out.println("Start: "+relativeStart+" End: "+relativeEnd +" Length: "+subannotation.getSize()+" index = "+(numBins-i-1));
										scores = DGE.scoreStrandedGene(libDataModelP,libDataModelN,subannotation);
										c[numBins-i-1] = scores[3];
										if(scores[3]>max)
											max = scores[3];
										//System.out.println("Step = "+step);
										//System.out.println("Start = "+relativeStart+" End = "+relativeEnd);
									}
									
								}
								else{
									for(int i=0;i<numBins;i++){
										int j=i+1;
										RefSeqGene subannotation = null;
										int relativeStart = i*step;
										int relativeEnd = j*step;
										if(relativeEnd>annotation.getTranscriptLength()){
											relativeEnd = annotation.getTranscriptLength();
										}
										subannotation = annotation.trim(relativeStart, relativeEnd);
										//System.out.println(annotation.getName()+" Start: "+relativeStart+" End: "+relativeEnd +" Length: "+subannotation.getSize()+" index = "+(numBins-i-1));
										//System.out.println(annotation.getName()+" "+annotation.getSize()+" Start: "+relativeStart+" End: "+relativeEnd +" index = "+(numBins-i-1));
										scores = DGE.scoreStrandedGene(libDataModelP,libDataModelN,subannotation);
										c[numBins-i-1] = scores[3];
										if(scores[3]>max)
											max = scores[3];
									}
									//System.out.println("Step: "+step+" "+annotation.getName()+": Start = "+relativeStart+" End = "+relativeEnd+" AnnotationLength = "+subannotation.getTranscriptLength());
									
									//System.out.println(coverage[i]);	
								}
								if(max>0){
								for(int i=0;i<numBins;i++){
									c[i] = c[i]/max;
									coverage[i] += c[i];
									System.out.print(c[i]+"\t");
								}
								System.out.println();
								counter++;
								}
							}
						}
					}
				}
			}
		}
		
		for(int i=0;i<numBins;i++){
			coverage[i] =coverage[i]/(double)counter;
		}
		
		/*
		 * Output Filename
		 */
		//String outputFileName = argMap.get("out");
		BufferedWriter outBw =argMap.getOutputWriter();
		/*
		 * Write the coverage in a row
		 */
		//BufferedWriter outBw = new BufferedWriter(new FileWriter(outputFileName));
		for(int i=0;i<coverage.length;i++){
			outBw.write(coverage[i]+"\t");
		}
		outBw.newLine();
		outBw.close();
		//writeOutputToFile(outputFileName,coverage);
		
	}
	
	/**
	 * This function writes the output to two files
	 * @param outFileName Name of the output file provided by user
	 * @throws IOException
	 */
	private static void writeOutputToFile (String outFileName,double[] coverage) throws IOException{
		
		/*
		 * Write the coverage in a row
		 */
		BufferedWriter outBw = new BufferedWriter(new FileWriter(outFileName));
		for(int i=0;i<coverage.length;i++){
			outBw.write(coverage[i]+"\t");
		}
		outBw.newLine();
		outBw.close();
	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException{

		NormalizedGeneEndCoverage dummy = new NormalizedGeneEndCoverage(args);
	}

}
