package broad.pda.geneexpression.tools;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

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
public class Gene5endCoverage {

	
	static final String usage = "Usage: Gene5EndCoverage -task <task name> "+
			"\n\tcompute: Computes the normalized coverage along the normalized length of the gene" +
			"\n\t\t-annotations <Specific regions to segment [BED by default]> "+
			"\n\t\t-alignment <Alignment (mapped to genome)> "+
			"\n\t\t-bins <Number of bins [100 by default]>"+
			"\n\t\t-windowSize <Step size/window size [By default: GeneLength/#bins]"+
			"\n\t\t-out <Output Filename>" +
			"\n\t\t-stranded";
	
	static Logger logger = Logger.getLogger(Gene5endCoverage.class.getName());
	int numBins;
	List<Double []> coverage;
	
	public Gene5endCoverage(String[] args)throws IOException{
	
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
		
		boolean isStranded = argMap.isPresent("stranded");
		/*
		 * Check the format of the annotation files and call the GTF or BED parser accordingly
		 */
		BEDFileParser annotations =  annotationsFile.endsWith(".gtf") || annotationsFile.endsWith(".GTF")? new GTFFileParser(annotationsFile) : new BEDFileParser(annotationsFile);
		
		numBins = argMap.isPresent("bins")? argMap.getInteger("bins") : 100;
		coverage = new ArrayList<Double []>();
		
		String outputFileName = argMap.getOutput();
		/*
		 * Read the name of the alignment file
		 */
		String alignmentFile = argMap.getMandatory("alignment");
		
		/*
		 * Iterate through list of annotations. Iterator is over chromosome names.
		 */
		Iterator<String> chrom_iter = annotations.getChromosomeIterator();
		
		/**
		 * To output gene lengths file
		 */
		BufferedWriter bwW = new BufferedWriter(new FileWriter(annotationsFile+".lengths"));
		int genesConsidered = 0;
		
		double[] scores;
		int counter = 0;
		//BufferedWriter bw = new BufferedWriter(new FileWriter(new String(argMap.get("out")+".coverage.log")));
		//BufferedWriter bwAll = new BufferedWriter(new FileWriter(new String(argMap.get("out")+".coverage.all")));
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputFileName+".log"));
		BufferedWriter bwAll = new BufferedWriter(new FileWriter(outputFileName+".all"));
		if(!isStranded){
		/*
		 * For each chromosome
		 */
		while(chrom_iter.hasNext()) {
			/*
			 * Initialize the data model using the alignment file
			 * @param: <alignment flieName> <load_chromosome_stats> <minMappingQuality> <remove_PCR_duplicates> <weigh reads by NH flag>
			 */
			ContinuousDataAlignmentModel libDataModel = AlignmentUtils.loadAlignmentData(alignmentFile,false,5,false,true);
			
			String chr = chrom_iter.next();
			/*
			 * If the alignment data has data from that chromosome
			 */
			if(libDataModel.hasDataForChromosome(chr)){
				logger.info("Processing " + chr+"\n");
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
				logger.info("Entering loop\n");
				while(annotation_iter.hasNext()){
					/*
					 * annotation = current node of RefSeqGeneWithIsoforms in interval tree
					 */
					RefSeqGeneWithIsoforms annotation_with_isoforms = annotation_iter.next();
					Iterator<RefSeqGene> isoform_iter = annotation_with_isoforms.getAllIsoforms().iterator();

					logger.info("Processing " + annotation_with_isoforms.getName()+"\n");
					while(isoform_iter.hasNext()){

						double max = 0.0;
						RefSeqGene annotation = isoform_iter.next();
						bwAll.write(annotation.getName()+"\t"+annotation.getTranscriptLength());
						int annotationLength = annotation.getTranscriptLength();
						//System.out.println(Math.round((float)((annotationLength)/(float)numBins)));
						int step = argMap.isPresent("windowSize")? argMap.getInteger("windowSize") : (int)((float)(annotationLength)/(float)numBins);
						//Step = gene length/ # bins
						int threshold = ((step)*numBins);
						//System.out.println(" threshold: "+(step)*numBins);
						//System.out.print(annotation.getName()+" "+annotation.getSize()+" "+annotation.getTranscriptLength());
						bwW.write(annotation.getName()+"\t");
						bwW.write(new Integer(annotation.getTranscriptLength()).toString());
						bwW.newLine();
						//System.err.println(annotation.getName()+"\t"+annotation.getTranscriptLength());
						
						if(annotationLength<threshold || threshold<numBins){
							//DONT PROCESS
							bw.write(annotation.getName()+" does not pass length threshold. Length: "+annotation.getTranscriptLength()+"\n");
							
						}
						else{
							bw.write(annotation.getName()+" passes. Length: "+annotation.getTranscriptLength()+"\n");
							genesConsidered++;
							//From gene end, for each bin of size step
							Double cov[] = new Double[numBins];
							if(annotation.getOrientation().equals("-")){
								for(int i=0;i<numBins;i++){
									int j=i+1;
									RefSeqGene subannotation = null;
									int relativeStart = annotation.getTranscriptLength()-(j*step);
									int relativeEnd = annotation.getTranscriptLength()-(i*step);
									if (relativeStart< 0){
										relativeStart = 0;
									}
									
									subannotation = annotation.trim(relativeStart, relativeEnd);
									//bw.write("+ strand. Start: "+relativeStart+" End: "+relativeEnd +" Length: "+subannotation.getSize()+" index = "+(numBins-i-1)+"\n");
									scores = libDataModel.scoreGene(subannotation);
									//coverage[numBins-i-1] = coverage[numBins-i-1]+scores[3];
									
									cov[i] = scores[3];
									if(cov[i]>max)
										max = cov[i];
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
									//bw.write("- strand. Start: "+relativeStart+" End: "+relativeEnd +" Length: "+subannotation.getSize()+" index = "+(numBins-i-1)+"\n");
									scores = libDataModel.scoreGene(subannotation);
									//coverage[numBins-i-1] = coverage[numBins-1-i]+scores[3];
									
									cov[i] = scores[3];
									if(cov[i]>max)
										max = cov[i];
								}
								
								//System.out.println("Step: "+step+" "+annotation.getName()+": Start = "+relativeStart+" End = "+relativeEnd+" AnnotationLength = "+subannotation.getTranscriptLength());
								
								//System.out.println(coverage[i]);	
							}
							if(max>0){
							for(int i=0;i<cov.length;i++)
								cov[i] = cov[i]/max;
							
							coverage.add(counter, cov);
							for(int i=0;i<cov.length;i++){
								bwAll.write(cov[i]+"\t");
							}
							counter++;
							}
							//bw.write(counter+"\n");	
						}
						bwAll.newLine();
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
			/*
			 * For each chromosome
			 */
			while(chrom_iter.hasNext()) {
				String chr = chrom_iter.next();
				/*
				 * If the alignment data has data from that chromosome
				 */
				if(libDataModelP.hasDataForChromosome(chr) || libDataModelN.hasDataForChromosome(chr)){
					logger.info("Processing " + chr+"\n");
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
					logger.info("Entering loop\n");
					while(annotation_iter.hasNext()){
						/*
						 * annotation = current node of RefSeqGeneWithIsoforms in interval tree
						 */
						RefSeqGeneWithIsoforms annotation_with_isoforms = annotation_iter.next();
						Iterator<RefSeqGene> isoform_iter = annotation_with_isoforms.getAllIsoforms().iterator();

						logger.info("Processing " + annotation_with_isoforms.getName()+"\n");
						while(isoform_iter.hasNext()){

							double max = 0.0;
							RefSeqGene annotation = isoform_iter.next();
							bwAll.write(annotation.getName()+"\t"+annotation.getTranscriptLength());
							int annotationLength = annotation.getTranscriptLength();
							//System.out.println(Math.round((float)((annotationLength)/(float)numBins)));
							int step = argMap.isPresent("windowSize")? argMap.getInteger("windowSize") : (int)((float)(annotationLength)/(float)numBins);
							//Step = gene length/ # bins
							int threshold = ((step)*numBins);
							//System.out.println(" threshold: "+(step)*numBins);
							//System.out.print(annotation.getName()+" "+annotation.getSize()+" "+annotation.getTranscriptLength());
							bwW.write(annotation.getName()+"\t");
							bwW.write(new Integer(annotation.getTranscriptLength()).toString());
							bwW.newLine();
							//System.err.println(annotation.getName()+"\t"+annotation.getTranscriptLength());
							
							if(annotationLength<threshold || threshold<numBins){
								//DONT PROCESS
								bw.write(annotation.getName()+" does not pass length threshold. Length: "+annotation.getTranscriptLength()+"\n");
								
							}
							else{
								bw.write(annotation.getName()+" passes. Length: "+annotation.getTranscriptLength()+"\n");
								genesConsidered++;
								//From gene end, for each bin of size step
								Double cov[] = new Double[numBins];
								if(annotation.getOrientation().equals("-")){
									for(int i=0;i<numBins;i++){
										int j=i+1;
										RefSeqGene subannotation = null;
										int relativeStart = annotation.getTranscriptLength()-(j*step);
										int relativeEnd = annotation.getTranscriptLength()-(i*step);
										if (relativeStart< 0){
											relativeStart = 0;
										}
										
										subannotation = annotation.trim(relativeStart, relativeEnd);
										//bw.write("+ strand. Start: "+relativeStart+" End: "+relativeEnd +" Length: "+subannotation.getSize()+" index = "+(numBins-i-1)+"\n");
										scores = DGE.scoreStrandedGene(libDataModelP,libDataModelN,subannotation);
										//coverage[numBins-i-1] = coverage[numBins-i-1]+scores[3];
										
										cov[i] = scores[3];
										if(cov[i]>max)
											max = cov[i];
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
										//bw.write("- strand. Start: "+relativeStart+" End: "+relativeEnd +" Length: "+subannotation.getSize()+" index = "+(numBins-i-1)+"\n");
										scores = DGE.scoreStrandedGene(libDataModelP,libDataModelN,subannotation);
										//coverage[numBins-i-1] = coverage[numBins-1-i]+scores[3];
										
										cov[i] = scores[3];
										if(cov[i]>max)
											max = cov[i];
									}
									
									//System.out.println("Step: "+step+" "+annotation.getName()+": Start = "+relativeStart+" End = "+relativeEnd+" AnnotationLength = "+subannotation.getTranscriptLength());
									
									//System.out.println(coverage[i]);	
								}
								if(max>0){
								for(int i=0;i<cov.length;i++)
									cov[i] = cov[i]/max;
								
								coverage.add(counter, cov);
								for(int i=0;i<cov.length;i++){
									bwAll.write(cov[i]+"\t");
								}
								counter++;
								}
								//bw.write(counter+"\n");	
							}
							bwAll.newLine();
						}
					}
				}
			}
		}
		bw.write(counter+"\n");	
		bw.close();
		bwW.close();
		bwAll.close();
		
		logger.info(genesConsidered+" genes were considered for this calculation.");
		
		/*for(int i=0;i<numBins;i++){
			coverage[i] =coverage[i]/(double)counter;
		}*/
		double[] means = new double[numBins];
		for(int j=0;j<numBins;j++){
			means[j] = 0.0;
		}
		for(int i=0;i<coverage.size();i++){
			for(int j=0;j<coverage.get(i).length;j++){
				means[j] = means[j]+coverage.get(i)[j];
			}
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
		/*for(int i=0;i<coverage.length;i++){
			outBw.write(coverage[i]+"\t");
		}*/
		for(int i=0;i<means.length;i++){
			outBw.write((means[i]/(double)coverage.size())+"\t");
		}
		outBw.newLine();
		outBw.close();
		//writeOutputToFile(outputFileName,coverage);
		
	}
		
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException{

		Gene5endCoverage dummy = new Gene5endCoverage(args);
	}

}
