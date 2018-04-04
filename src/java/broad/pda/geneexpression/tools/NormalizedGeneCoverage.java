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
import broad.pda.seq.alignment.AlignmentUtils;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
/**
 * @author skadri
 *
 */
public class NormalizedGeneCoverage {
	
	static final String usage = "Usage: NormalizedGeneCoverage -task <task name> "+
			"\n\tcoverage: Computes the normalized coverage along the normalized length of the gene" +
			"\n\t\t-annotations <Specific regions to segment [BED by default]> "+
			"\n\t\t-alignment <Alignment (mapped to genome)> "+
			"\n\t\t-bins <Number of bins [100 by default]>"+
			"\n\t\t-windowSize <Step size/window size [By default: GeneLength/#bins]"+
			"\n\t\t-out <Output Filename>";
	
	static Logger logger = Logger.getLogger(NormalizedGeneCoverage.class.getName());
	int numBins;
	double[] coverage;
	
	public NormalizedGeneCoverage(String[] args)throws IOException{
	
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
		/*
		 * Initialize the data model using the alignment file
		 * @param: <alignment flieName> <load_chromosome_stats> <minMappingQuality> <remove_PCR_duplicates> <weigh reads by NH flag>
		 */
		ContinuousDataAlignmentModel libDataModel = 	AlignmentUtils.loadAlignmentData(alignmentFile,true,5,false,true);
		
		/*
		 * Iterate through list of annotations. Iterator is over chromosome names.
		 */
		Iterator<String> chrom_iter = annotations.getChromosomeIterator();
		
		double[] scores;
		int counter = 0;
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

						RefSeqGene annotation = isoform_iter.next();
						//System.out.println(annotation.getName());
						int annotationLength = annotation.getTranscriptLength();
						//Step = gene length/ # bins
						if(annotationLength<(2*numBins)){
							//DONT PROCESS
						}
						else{
							int step = argMap.isPresent("windowSize")? argMap.getInteger("windowSize") : (int)((annotationLength)/numBins);
	
							//From gene start, for each bin of size step
							for(int i=0;i<numBins;i++){
								int j=i+1;
								RefSeqGene subannotation = null;
								int relativeStart = 0;
								int relativeEnd   = 0;
								if(annotation.getOrientation().equals("-")){
									relativeStart = annotation.getTranscriptLength()-(j*step)+1;
									relativeEnd = annotation.getTranscriptLength()-(i*step);
									if (relativeStart< 0){
										relativeStart = 0;
									}
									//System.out.println("Step = "+step);
									//System.out.println("Start = "+relativeStart+" End = "+relativeEnd);
								}	
								else{
									relativeStart = i*step;
									relativeEnd = j*step-1;
									if(relativeEnd>annotation.getTranscriptLength()){
										relativeEnd = annotation.getTranscriptLength()-1;
									}
								}
								subannotation = annotation.trim(relativeStart, relativeEnd);
								//System.out.println(annotation.getName()+": Start = "+relativeStart+" End = "+relativeEnd+" GeneLength = "+annotation.getTranscriptLength());
								scores = libDataModel.scoreGene(subannotation);
								coverage[i] = coverage[i]+scores[3];
								//System.out.println(coverage[i]);	
							}
							counter++;
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

		NormalizedGeneCoverage dummy = new NormalizedGeneCoverage(args);
	}

}
