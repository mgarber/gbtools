package broad.pda.geneexpression.tools;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.math.MathException;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.sam.Alignment;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.error.ParseException;
import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.gene.RefSeqGeneWithIsoforms;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.AlignmentDataModelStats;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

public class ReadPileupScanner {
	//[0] positive
	//[1] negative
	private AlignmentDataModel model;
	//private ContinuousDataAlignmentModel modelN;
	int windowSize;
	int extension;
	//Map<String, Integer> chromosomeSizes;
	private int DEFAULT_UPSTREAM_REGION = 2000;
	String outputName;
	private static final int DEFAULT_CHUNK_SIZE = 1000000;
	
	protected static Logger logger = Logger.getLogger(ReadPileupScanner.class.getName());
	

	static final String usage = "Usage: DGE -task <task name> "+
			"\n\tTASK 1: annotate: Identifies the significant 5' gene ends in the end RNA-Seq data for a given annotation set." + 
			"\n\n\tTASK 2: visualize: Outputs wig files to visualize peaks in end RNA-Seq data on a genomic scale." + 
			"\n**************************************************************"+
			"\n\t\tArguments"+
			"\n**************************************************************"+
			"\n\n\t\t-alignment <Alignment (mapped to genome) for which expression must be calculated. Index file MUST be provided.> "+
			"\n\t\t-out <Output file [Defaults to stdout]> "+
			"\n\t\t-annotations <FOR ANNOTATE TASK ONLYAnnotation file for which to calculate expression. [BED by default]> "+

			"\n\n**************************************************************"+
			"\n\t\tOptional Arguments"+
			"\n**************************************************************"+
			"\\nn\t\t-window <Specifies the size of the fixed window use to scan the genome. We recommend a size of 2-5bp for better resolution of gene ends. Default = 5bp> "+
			"\n\t\t-extension <Specifies the size of the fixed window use to scan the genome. We recommend a size of 2-5bp for better resolution of gene ends> "+

			"\n";

	
	public ReadPileupScanner(String[] args) throws IllegalArgumentException, IOException{
		
		
		/*
		 * Gives a log4j error. Check later.
		 */
		Globals.setHeadless(true);
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"annotate");
		
		model=new GenericAlignmentDataModel(argMap.getMandatory("alignment"), null, false, 5,false,true);
		//Default window size = 5bp
		windowSize = argMap.getInteger("window", 5);
		outputName = argMap.getOutput();
		
		if(argMap.getTask().equalsIgnoreCase("annotate")){
			extension = argMap.getInteger("extension", DEFAULT_UPSTREAM_REGION);
			assignGenes(argMap.getMandatory("annotations"));
		}
		else{
			if(argMap.getTask().equalsIgnoreCase("visualize")){
				scan();
			}
			else{
				logger.error(usage);
			}
		}
	}
	public ReadPileupScanner(String sample,int windowS,String sizes,String outFile) throws IOException{
		
		//chromosomeSizes = (BEDFileParser.loadChrSizes(sizes));
		model=new GenericAlignmentDataModel(sample, sizes, false, 5,false,true);
		windowSize = windowS;
		outputName=outFile;
	}
	
	public void assignGenes(String annotationFile) throws IOException{
		
		//Read annotation file
		BEDFileParser annotationParser = new BEDFileParser(annotationFile);		
		//Output files
		BufferedWriter bwP = new BufferedWriter(new FileWriter(outputName+".plus.wig"));
		BufferedWriter bwM = new BufferedWriter(new FileWriter(outputName+".minus.wig"));
		model.setSecondRead();
		
		AlignmentDataModelStats modeldata = new AlignmentDataModelStats(model);
		ContinuousDataAlignmentModel datamodel = new ContinuousDataAlignmentModel(modeldata);
		
		Map<RefSeqGene,Collection<Alignments>> geneToPeakMap = new HashMap<RefSeqGene,Collection<Alignments>>();
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputName));
		BufferedWriter bwPeaks = new BufferedWriter(new FileWriter(outputName+".peaks"));
		//For each chromosome
		for(String chr:model.getChromosomeLengths().keySet()){
			
			double lambda = datamodel.getLambda(chr);
			//If the chromosome is not expressed in the sample
			if(lambda==0.0){
				logger.warn(chr +" is not expressed in sample");
			}
			//If the chromosome is expressed in the sample
			else{
				logger.info("Processing " + chr);
				//Start of new chromosome
				bwP.write("variableStep chrom="+chr+" span="+windowSize+"\n");
				bwM.write("variableStep chrom="+chr+" span="+windowSize+"\n");
				
				/*
				 * Obtain an iterator over the interval tree values built from the annotations (genes) on that particular chromosome.
				 */
				if(annotationParser.containChr(chr)){
				Iterator<RefSeqGeneWithIsoforms> annotation_iter = annotationParser.getChrTree(chr).valueIterator();
				/*
				 * Parse annotation tree
				 * While there is an annotated RefSeqGeneWithIsoforms in Interval tree to analyze
				 * Thus, for each gene
				 */
				while(annotation_iter.hasNext()){
					//For each gene
					RefSeqGene annotation = annotation_iter.next();
					
					Collection<Alignments> peaks = new ArrayList<Alignments>();
					
					List<Double> values = new ArrayList<Double>();
					
					//Get a gene window iterator over the entire gene
					Collection<RefSeqGene> genes = new ArrayList<RefSeqGene>();
					genes.add(annotation);
					
					IntervalTree<Alignment> tree=modeldata.getIntervalTreeCached(chr, annotation.getStart(), annotation.getEnd());
					Iterator<RefSeqGene> giter = new GeneWindowIterator(genes, windowSize, 0);
					//For each window
					while(giter.hasNext()){
						RefSeqGene window = giter.next();
						//logger.info(window.toUCSC());
						double windowCount = 0.0;
						//Get the reads in the window
						//For each exon in the gene
						for(Alignments exon: window.getExonSet()){
						//	long start=System.currentTimeMillis();
							Iterator<Node<Alignment>> overlappers = tree.overlappers(exon.getStart(), exon.getEnd());
						//	long end=System.currentTimeMillis();
							//System.err.println("Interval tree cached Time: "+(end-start));
						/*	start=System.currentTimeMillis();
							CloseableIterator<Alignment> overlappers2 = model.getAlignmentsOverlappingRegion(exon);
							end=System.currentTimeMillis();*/
							//System.err.println("Get alignments "+(end-start));
							//for all reads in the window
							while(overlappers.hasNext()){
								Node<Alignment> read = overlappers.next();
								if(read.getValue().isSecondOfPair()){
									if(readStartFallsInWindow(read.getValue(),exon)){
										//increment count of window
									/*	if((read.getValue().isNegativeStrand() && annotation.isNegativeStrand()) ||
												(!read.getValue().isNegativeStrand() && !annotation.isNegativeStrand()))*/
											windowCount += read.getNumReplicates();
									}
								}
							}
						/*	if(windowCount!=windowCount2)
								System.err.println(window.toUCSC()+" Cached: "+windowCount+" Slow:"+windowCount2);
							overlappers2.close();*/
						}
						//Get the score of each window into an array
						if(windowCount>0.0){
							values.add(windowCount);
						}
					}
					if(values.size()>0){
						//Get the median and standard deviation
						double median = Statistics.median(values);
						double variance = Statistics.variance(values);
						double mean=Statistics.mean(values);
						logger.info("For "+annotation.getName()+" Median = "+median+" Variance = "+variance+ " Mean = "+mean);
						
						//Get a genome window iterator over the genomic region of gene and upstream region
						int start = annotation.getStart();
						int end  = annotation.getEnd();
						if(annotation.isNegativeStrand()){
							end = annotation.getEnd() + getUpstreamNonOverlappingDistance(annotationParser, annotation, extension);
							//logger.info("Start: "+start+" End: "+end);
						}
						else{
							start = annotation.getStart() - getUpstreamNonOverlappingDistance(annotationParser, annotation, extension);
						}
						
						//logger.info("Start = "+start+" End = "+end);
						IntervalTree<Alignment> datatree=modeldata.getIntervalTreeCached(chr, start, end);
						Iterator<Alignments> witer = new GenomeWindowsIterator(chr, windowSize, start, end, 0);
						while(witer.hasNext()){
							Alignments window = witer.next();
							//logger.info(window.toUCSC());
							double windowCount = 0.0;
							//Get the reads in the window
							//CloseableIterator<Alignment> overlappers2 = model.getAlignmentsOverlappingRegion(window);
							Iterator<Node<Alignment>> overlappers = datatree.overlappers(window.getStart(), window.getEnd());
							//for all reads in the window
							while(overlappers.hasNext()){
								Node<Alignment> read = overlappers.next();
//								if(read.getValue().isSecondOfPair()){
									if(readStartFallsInWindow(read.getValue(),window)){
										//logger.info(window.toUCSC()+"\t"+read.getNumReplicates());
										//increment count of window
										windowCount += read.getNumReplicates();
									/*	if((read.getValue().isNegativeStrand() && annotation.isNegativeStrand()) ||
												(!read.getValue().isNegativeStrand() && !annotation.isNegativeStrand()))
											windowCount += read.getNumReplicates();*/
									}
//								}
							}
							//overlappers2.close();
							//Get the z-score of each window
							double zscore = Statistics.zScore(windowCount, mean, variance);
							double pval = Statistics.zscoreToPvalue(zscore);
							
							//Associate the high z-scores with gene
							if(zscore>6){
								peaks.add(window);
							}
							if(annotation.isNegativeStrand()){
								if(windowCount>0.0 && !Double.isInfinite(zscore))
									bwM.write(window.getStart()+"\t"+zscore+"\n");							}
							else{
								if(windowCount>0.0 && !Double.isInfinite(zscore))
									bwP.write(window.getStart()+"\t"+zscore+"\n");
							}
					//		}
							if(windowCount>0.0 && (!Double.isInfinite(zscore)||!Double.isNaN(zscore)))
								bw.write(annotation.getName()+"\t"+window.toUCSC()+"\t"+windowCount+"\t"+zscore+"\n");
						}
						if(peaks.size()>0)
							geneToPeakMap.put(annotation, peaks);
					}
				}
			}
			}
		}
		
		bwP.close();
		bwM.close();
		BufferedWriter bwbed = new BufferedWriter(new FileWriter(outputName+".bed"));
		for(RefSeqGene g:geneToPeakMap.keySet()){
			for(Alignments p:geneToPeakMap.get(g)){
				bwPeaks.write(g.getName()+"\t"+p.toUCSC()+"\n");
				bwbed.write(new RefSeqGene(p).toBED()+"\n");
			}
		}
		bw.close();
		bwbed.close();
		bwPeaks.close();
	}
	
	
	private int getUpstreamNonOverlappingDistance(BEDFileParser annotationCollection,RefSeqGene thisGene,int maxRegion){
		
		if(thisGene.isNegativeStrand()){
			RefSeqGene closestUpstreamGene = annotationCollection.getClosestDownstream(thisGene);
			if(closestUpstreamGene==null){
				return maxRegion;
			}
			else{
				int distanceFromThisGene = closestUpstreamGene.getStart()-thisGene.getEnd();
				if(distanceFromThisGene>maxRegion){
					return maxRegion;
				}
				else{
					return distanceFromThisGene;
				}
			}
		}
		else{
			RefSeqGene closestUpstreamGene = annotationCollection.getClosestUpstream(thisGene);
			if(closestUpstreamGene==null){
				return maxRegion;
			}
			else{
				int distanceFromThisGene = thisGene.getStart() - closestUpstreamGene.getEnd();
				if(distanceFromThisGene>maxRegion){
					return maxRegion;
				}
				else{
					return distanceFromThisGene;
				}
			}
		}
	}
	/**
	 * Scan the genome in bins of window size and outputs wig files where each count represents 
	 * the number of second of pair reads having a start within that bin
	 * @param outFile
	 * @throws IOException
	 */
	public void scan() throws IOException{
		
		BufferedWriter bwP = new BufferedWriter(new FileWriter(outputName+".plus.wig"));
		BufferedWriter bwM = new BufferedWriter(new FileWriter(outputName+".minus.wig"));
		model.setSecondRead();
		AlignmentDataModelStats modeldata = new AlignmentDataModelStats(model);
		ContinuousDataAlignmentModel datamodel = new ContinuousDataAlignmentModel(modeldata);
		//For each chromosome
		for(String chr:model.getChromosomeLengths().keySet()){
			
			double lambda = datamodel.getLambda(chr);
			//If the chromosome is not expressed in the sample
			if(lambda==0.0){
				logger.warn(chr +" is not expressed in sample");
			}
			//If the chromosome is expressed in the sample
			else{
				//load chunks into memory: first chunk
				int chunkNumber=1;
				IntervalTree<Alignment> chunkAlignmentTree=modeldata.getIntervalTreeCached(chr, 0, chunkNumber*DEFAULT_CHUNK_SIZE);
				boolean cached=false;
				
				//Iterator<Alignments> iter = new WindowsIterator(chr, windowSize, 108773064, 108813968, 0);
				Iterator<Alignments> iter = new GenomeWindowsIterator(chr, windowSize, 0, model.getChromosomeLength(chr), 0);
				//Start of new chromosome
				bwP.write("variableStep chrom="+chr+" span="+windowSize+"\n");
				bwM.write("variableStep chrom="+chr+" span="+windowSize+"\n");
				while(iter.hasNext()){
					Alignments window = iter.next();
					//logger.info(window.toUCSC());
					double windowCountPlus = 0.0;
					double windowCountMinus = 0.0;
					//Get the reads in the window
					Iterator<Node<Alignment>> overlappers = chunkAlignmentTree.overlappers(window.getStart(),window.getEnd());
					//for all reads in the window
					while(overlappers.hasNext()){
						Node<Alignment> read = overlappers.next();
						//if(read.isPaired()){
							if(read.getValue().isSecondOfPair()){
								if(readStartFallsInWindow(read.getValue(),window)){
									//increment count of window
									if(read.getValue().isNegativeStrand())
										windowCountMinus += read.getNumReplicates();
									else
										windowCountPlus += read.getNumReplicates();
								}
							}
						//}
					}
					writeToWig(bwP,window,windowCountPlus);
					writeToWig(bwM,window,windowCountMinus);
				}
			}
			
		}	
		bwP.close();
		bwM.close();
	}
	
	
	/**
	 * Scan the genome in bins of window size and outputs wig files where each count represents 
	 * the number of second of pair reads having a start within that bin
	 * @param outFile
	 * @throws IOException
	 */
	public void scan(String chr,int start,int end) throws IOException{
		
		BufferedWriter bwP = new BufferedWriter(new FileWriter(outputName+".plus.wig"));
		BufferedWriter bwM = new BufferedWriter(new FileWriter(outputName+".minus.wig"));
		model.setSecondRead();
		AlignmentDataModelStats modeldata = new AlignmentDataModelStats(model);
		ContinuousDataAlignmentModel datamodel = new ContinuousDataAlignmentModel(modeldata);
		
		double lambda = datamodel.getLambda(chr);
		//If the chromosome is not expressed in the sample
		if(lambda==0.0){
			logger.warn(chr +" is not expressed in sample");
		}
		//If the chromosome is expressed in the sample
		else{
			bwP.write("variableStep chrom="+chr+" span="+windowSize+"\n");
			bwM.write("variableStep chrom="+chr+" span="+windowSize+"\n");
			//load chunks into memory: first chunk
			IntervalTree<Alignment> chunkAlignmentTree=modeldata.getIntervalTreeCached(chr, start,end);
			
			// for every position in chromosome
			Iterator<Alignments> iter = new GenomeWindowsIterator(chr, windowSize, start, end, 0);
			while(iter.hasNext()){
				Alignments window = iter.next();
				//logger.info(window.toUCSC());
				double windowCountPlus = 0.0;
				double windowCountMinus = 0.0;
				//Get the reads in the window
				Iterator<Node<Alignment>> overlappers = chunkAlignmentTree.overlappers(window.getStart(),window.getEnd());
				//for all reads in the window
				while(overlappers.hasNext()){
					Node<Alignment> read = overlappers.next();
					//if(read.isPaired()){
						if(read.getValue().isSecondOfPair()){
							if(readStartFallsInWindow(read.getValue(),window)){
								//increment count of window
								if(read.getValue().isNegativeStrand())
									windowCountMinus += read.getNumReplicates();
								else
									windowCountPlus += read.getNumReplicates();
							}
						}
					//}
				}
				writeToWig(bwP,window,windowCountPlus);
				writeToWig(bwM,window,windowCountMinus);
			}
		}

		bwP.close();
		bwM.close();
	}
	
	
	/**
	 * Returns true if the oriented start of the read falls in the window
	 * @param align
	 * @param window
	 * @return
	 */
	private boolean readStartFallsInWindow(Alignment align,Alignments window){
		
		int start;
		if(align.isNegativeStrand()){
			start = align.getAlignmentEnd();
		}
		else{
			start = align.getAlignmentStart();
		}
		if(start>=window.getStart() && start<=window.getEnd())
			return true;
		else
			return false;
	}
	
	/**
	 * Writes the window score to the wig file
	 * @param bw
	 * @param window
	 * @param windowCount
	 * @throws IOException
	 */
	private void writeToWig(BufferedWriter bw, Alignments window,double windowCount) throws IOException{
		if(windowCount>0.0)
			bw.write(window.getStart()+"\t"+windowCount+"\n");
		else
			;
	}
	
	private static class GenomeWindowsIterator implements Iterator<Alignments> {

		private String chr;
		private int currPosition;
		private int windowSize;
		private int step;	
		private int end;
		
		/**
		 * Constructs an iterator on chromosome chr starting at start, ending at end, of size, windowSize with overlap between consecutive windows 
		 * @param chr
		 * @param windowSize
		 * @param start
		 * @param overlap
		 * @param end
		 */
		public GenomeWindowsIterator(String chr, int windowSize, int start, int end, int overlap){
			if (overlap > windowSize) throw new IllegalArgumentException("Overlap cannot be greater than window size.");
			this.chr = chr;
			this.currPosition = start;
			this.windowSize = windowSize;
			this.step = windowSize - overlap;
			this.end = end;
		}
		
		/**
		 * Returns true if there is another window
		 */
		public boolean hasNext() {
			return (this.end >= (this.currPosition+this.windowSize));
		}

		/**
		 * Returns the next window
		 */
		public Alignments next() {
			Alignments w=new Alignments(this.chr,this.currPosition,this.currPosition+this.windowSize);
			this.currPosition = this.currPosition+this.step;
			return w;
		}

		/**
		 * N/A here
		 */
		public void remove() {
			throw new UnsupportedOperationException(); 
		}
		
	}	
	
	private class GeneWindowIterator implements Iterator<RefSeqGene> {

		int windowSize;
		int step;
		Iterator<? extends RefSeqGene> genes;
		Iterator<? extends RefSeqGene> currentGeneWindows;
		
		GeneWindowIterator(Collection<? extends RefSeqGene> genes, int windowSize, int overlap){
			this.windowSize=windowSize;
			this.step=windowSize-overlap;
			this.genes=genes.iterator();
		}
		
		@Override
		public boolean hasNext() {
			//first, see if there are windows in the current gene
			if(currentGeneWindows!=null && currentGeneWindows.hasNext()){
				return true;
			}
			//else, see if there are still genes
			else if(genes.hasNext()){
				return true;
			}
			return false;
		}

		@Override
		public RefSeqGene next() {
			//if there are windows still in currentGeneWindows send these
			if(currentGeneWindows!=null && currentGeneWindows.hasNext()){
				return currentGeneWindows.next();
			}
			else{
				//else make a new currentGeneWindows
				//make new currentGeneWindows
				RefSeqGene nextGene=genes.next();
				//logger.info(nextGene);
				this.currentGeneWindows=makeWindows(nextGene, this.windowSize, this.step);
				return next();
			}
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException("TODO");
			
		}
		
		private Iterator<? extends RefSeqGene> makeWindows(RefSeqGene gene, int windowSize, int step){
			Collection<? extends RefSeqGene> rtrn=gene.getWindows(windowSize, step, 0);
			return rtrn.iterator();
		}
		
	}

	public static void main (String [] args) throws MathException, ParseException, IOException {
		
		ReadPileupScanner dummy = new ReadPileupScanner(args);
	/*	if(args.length==5){
			//String sample,int window,String sizes,String outFile
			ReadPileupScanner dummy = new ReadPileupScanner(args[0], new Integer(args[1]),args[2],args[3]);
			//annotation file
			dummy.assignGenes(args[4]);
		}
		else if(args.length==7){
			//String sample,int window,String sizes,String outFile
			ReadPileupScanner dummy = new ReadPileupScanner(args[0], new Integer(args[1]),args[2],args[3]);
			//annotation file
			dummy.scan(args[4],new Integer(args[5]),new Integer(args[6]));
		}
		else if(args.length==4){
			//String sample,int window,String sizes,String outFile
			ReadPileupScanner dummy = new ReadPileupScanner(args[0], new Integer(args[1]),args[2],args[3]);
			//annotation file
			dummy.scan();
		}
		else{System.err.println(usage);}*/

	}
	
	//static String usage=" args[0]=input bam file \n\t args[1]=window size \n\t args[2]: Sizes file for genome"
	//		+"\n\t args[3]= outputName \n\t args[4] annotation file";
}
