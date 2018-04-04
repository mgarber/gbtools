package broad.pda.seq.segmentation.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.datastructures.IntervalTree;
import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.annotation.GTFFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.gene.RefSeqGeneWithIsoforms;
import broad.pda.seq.pairedend.PairedEndDistribution;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.AlignmentDataModelStats;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

public class PolyASeqAnalyzer {

	static int DEFAULT_MIN_MAPPING_QUALITY = -1;
	public static int DEFAULT_INSERT_SIZE_FUDGE = 20;
	static int DEFAULT_MAX_3P_EXTENSION = 3000;
	static int STEP = 50;
	static int WINDOW = 200;
	static int USE_SCORE = 1;
	static int DEFAULT_FRAGMENT_SIZE = 500;
	
	static Logger logger = Logger.getLogger(PolyASeqAnalyzer.class.getName());
	
	String outputFileName;
	int minMappingQuality;
	String alignmentFile;
	int[] windows;
	int extensionFactor;
	String sizes;
	boolean printFullScores;
	BEDFileParser annotationParser;
	Map<RefSeqGene, List<RefSeqGene>> peakMap;
	int max3PExtension;
	int fragmentSize;
	
	static String usage="PolyASeqAnalyzer -task score "+
			"\n\t\t -alignment <Alignment file in BAM, SAM or Alignemnt format> "+
			"\n\t\t-annotations <Genes to associate peaks with [BED by default]> "+
			"\n\t\t -maskFileDir <Mask File directory>" +
			"\n\t\t -out <Output file name>"+ 
			"\n\t\t -sizeFile <Chromosome size file> \n"+
			"\n\t\t -chr <Chromsomosome to segment> \n "+
			"-chrSequence <Necessary to filter spliced reads by splice site information. Notice that this is only compatible with region files that contain regions of only one chromosome> "+
			"\n\t\t -minMappingQuality <Minimum quality of the reads. Recommended: 5 [Default:0]>"+
			"\n\t\t -max3PExtension <It will allow the program to search for peaks max3PExtension bases downstream from the 3' end of the gene> "+
			"\n\t\t -windows <Comma separated list of windows to evaluate defaults to contiguous regions of coverage>"+
			"\n\t\t -extensionFactor <Extend reads by this factor (defaults to 0)>"+
			"\n\t\t -fullScores  <Including this flag will write a bed file with all scores for the peaks>"+
			"\n\t\t -findMaxContiguous <Including this flag will find maximum contiguous regions within windows>"+
			"\n\t\t -trim <Include this flag if trimming of the ends of windows based on read coverage  is desired this is expensive>"+
			"\n\t\t -trimQuantile <Coverage quantile below which end bases should be trimmed>"+
			"\n\t\t -alpha <Desired FDR>"+
			"\n\t\t -minLength <>"+
			"\n\t\t -stranded <Setting this flag indicates that the library is stranded>"+
			"\n\t\t -paired <Setting this flag indicates that the library is paired end>"+
			"\n\t\t\t -isSecondRead <>"+
			"\n -stranded <Is the data strand-specific>"+
			"\n -fragmentSize <Fragment Size>";


	public PolyASeqAnalyzer(ArgumentMap argmap) throws IOException {
		
		
		if("score".equalsIgnoreCase(argmap.getTask())) {

		//INITIALIZES IGV
		Globals.setHeadless(true);
		System.out.println("Using Version R4.4");
		logger.debug("DEBUG ON");
		outputFileName = argmap.getOutput();
		minMappingQuality = argmap.containsKey("minMappingQuality") ? argmap.getInteger("minMappingQuality") : DEFAULT_MIN_MAPPING_QUALITY;
		max3PExtension = argmap.containsKey("max3PExtension") ? argmap.getInteger("max3PExtension") : DEFAULT_MAX_3P_EXTENSION;
		alignmentFile = argmap.getMandatory("alignment");
		windows= ContinuousDataAlignmentModel.getWidths(argmap.getMandatory("windows"));
		extensionFactor = argmap.containsKey("extensionFactor") ? argmap.getInteger("extensionFactor") : 0;
		sizes = argmap.get("sizeFile");
		printFullScores = argmap.containsKey("fullScores");
		fragmentSize = argmap.isPresent("fragmentSize") ? argmap.getInteger("fragmentSize") : DEFAULT_FRAGMENT_SIZE;
		String annotationFile = argmap.getMandatory("annotations");
		annotationParser =  annotationFile.endsWith(".gtf") || annotationFile.endsWith(".GTF")? new GTFFileParser(annotationFile) : new BEDFileParser(annotationFile);
		peakMap = new HashMap<RefSeqGene, List<RefSeqGene>>();
		peakFinder(argmap);
		
		}
	}
	
	private void peakFinder(ArgumentMap argmap) throws IOException{
		
		//Optional parameters
		boolean findMaxContiguous = argmap.containsKey("findMaxContiguous");
		boolean trimEnds = argmap.containsKey("trim");
		double trimQuantile = argmap.isPresent("trimQuantile") ? argmap.getDouble("trimQuantile") : 0.25;
		int minRemainingLength = argmap.isPresent("minLength") ? argmap.getInteger("minLength") : Statistics.min(windows);
		double alpha = argmap.containsKey("alpha") ? argmap.getDouble("alpha") : .05;
		File[] maskFiles  = argmap.isPresent("maskFileDir") ?  new File(argmap.get("maskFileDir")).listFiles(): null ;
		Map<String, Integer> maskFileData = ContinuousDataAlignmentModel.parseMaskFiles(maskFiles);
		String chr = argmap.get("chr");
		boolean isStranded = argmap.containsKey("stranded");
		boolean isPaired = argmap.containsKey("paired");
		boolean isSecondRead = false;
		BEDFileParser peaks = new BEDFileParser();
		if(isPaired){
			isSecondRead = argmap.isPresent("isSecondRead")? (Boolean.parseBoolean(argmap.get("isSecondRead"))): true;
		}
		
		if(isStranded){
			
			/*
			 * Run for positive strand first 
			 */
			
			/*
			 * STEP 1: CALL PEAKS 
			 */
			logger.info("Using plus reads:");
			AlignmentDataModel alignments=new GenericAlignmentDataModel(alignmentFile, sizes, false, minMappingQuality,true,false);
					
			if(isPaired){
				if(isSecondRead){
					alignments.setSecondRead();
				}
				else{
					alignments.setFirstRead();
				}
			}
			
			logger.info("AlignmentDataModel loaded, initializing model stats");
			AlignmentDataModelStats alignmentData = new AlignmentDataModelStats(alignments, maskFileData);
			ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(alignmentData);

			/*
			 * SENSE READS FIRST
			 */
			alignments.setPositiveStranded();

			data.setMaskFileData(maskFileData);
			data.setExtensionFactor(extensionFactor);

			List<String> chromosomes = new ArrayList<String>();
			if(chr!= null && chr.length() > 0) {
				chromosomes.add(chr);
			} else {
				Map<String, Integer> chrSizes = data.getChromosomeLengths();
				chromosomes = new ArrayList<String>(chrSizes.keySet());
			}
			data.setTrimEnds(trimEnds);
			data.setMinContguousSegmentSize(minRemainingLength);
			data.setTrimQuantile(trimQuantile);
			data.setFindMaxContiguous(findMaxContiguous);

			Map<Alignments, double[]> scores = new HashMap<Alignments, double[]>();
			Map<RefSeqGene, PairedEndDistribution> insertSizeMap = new HashMap<RefSeqGene,PairedEndDistribution>();
			int totalMappedReadsPlus = 0;
			for(String workChr : chromosomes) {
				logger.info("Processing chromosome " + workChr + ", Scanning windows");
				totalMappedReadsPlus += data.getData().getSum(workChr);
				Collection<Alignments> segments = data.scan(windows, alpha, workChr);
				logger.info("Scoring segments");
				scores.putAll(data.scoreSegments(segments, workChr));
				for(Alignments a:segments){
					if(scores.containsKey(a)){
						RefSeqGene p = new RefSeqGene(a);
						p.setOrientation("+");
					//	System.out.println(p.getName());
					//	System.out.println(" score size:"+scores.get(a).length);
						p.setBedScore(scores.get(a)[USE_SCORE]);
						peaks.addRefSeq(p);
					}
				}
			}
			logger.info("Done. Total mapped reads: " + totalMappedReadsPlus);
			String outName = outputFileName+".plus.bed";
			BEDFileParser.writeSortedBED(outName, scores);
			if(printFullScores) {
				BEDFileParser.writeSortedBEDWithScores(outName+".scores.bed", scores);
			}
			
			/*
			 * ANTISENSE READS NEXT
			 */
			logger.info("Using minus reads:");
//			alignments=new GenericAlignmentDataModel(alignmentFile, sizes, false, minMappingQuality,true,false);
			alignments.setNegativeStranded();

			int totalMappedReadsMinus = 0;
			for(String workChr : chromosomes) {
				logger.info("Processing chromosome " + workChr + ", Scanning windows");
				totalMappedReadsMinus += data.getData().getSum(workChr);
				Collection<Alignments> segments = data.scan(windows, alpha, workChr);
				logger.info("Scoring segments");
				/*
				 * [0] : p-value
				 * [1] : Enrichment
				 * [2] : Sum of reads
				 * [3] : Sum/count
				 */
				scores.putAll(data.scoreSegments(segments, workChr));
				for(Alignments a:segments){
					if(scores.containsKey(a)){
						RefSeqGene p = new RefSeqGene(a);
						p.setOrientation("-");
						p.setBedScore(scores.get(a)[USE_SCORE]);
						peaks.addRefSeq(p);
					}
				}
			}
			
			logger.info("Done. Total mapped reads: " + totalMappedReadsMinus);
			outName = outputFileName+".minus.bed";
			BEDFileParser.writeSortedBED(outName, scores);
			if(printFullScores) {
				BEDFileParser.writeSortedBEDWithScores(outName+".scores.bed", scores);
			}
			
			
			//For each gene in the annotation file
			
			annotationParser.makeGenes();
			/*
			 * MERGE PEAKS TO "MAKE GENES" so that we dont have to check for overlaps later.
			 */
			peaks.merge();
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputFileName+".all.peaks.bed"));
			peaks.writeFullBed(bw);
			bw.close();
			
			/* 
			 *  STEP 2: CALCULATE THE INSERT SIZE DISTRIBUTION FOR EACH PEAK
			 */
			for(RefSeqGene peak: peaks.GetGenes()){
				PairedEndDistribution ped = new PairedEndDistribution(alignments, peak,isSecondRead);
				//ped.ensureNonZeroCounts();
				insertSizeMap.put(peak, ped);
			}
			
			/*
			 * STEP 3: ASSIGN TO GENES
			 */
			Map<RefSeqGene, List<RefSeqGene>> result = assignPeaks(peaks,annotationParser);

			bw = new BufferedWriter(new FileWriter(outputFileName));
			BufferedWriter bws = new BufferedWriter(new FileWriter(outputFileName+".summary.txt"));
			bws.write("GeneName\t#Peaks\n");
			bw.write("GeneName\tPeak\tDistanceFromEndGene\tScore\tInsertSize:Mean\tInsertSize:StdDev\tPolyASize\n");
			
			for(RefSeqGene gene: result.keySet()){
				
				if(result.get(gene)==null){
					bws.write(gene.getName()+"\t");
					bws.write(0+"\n");
				}
				else{
					for(RefSeqGene aPeak: result.get(gene)){
						bw.write(gene.getName()+"\t");
						bw.write(aPeak.toUCSC()+"\t");
						//positive strand
						if(gene.getOrientation().equals("+")){
							bw.write((aPeak.getEnd()-gene.getEnd())+"\t");
						}
						//negative strand
						else{
							bw.write((aPeak.getStart()-gene.getStart())+"\t");
						}
						
						bw.write(aPeak.getBedScore()+"\t");
						double meanInsert = insertSizeMap.get(aPeak).getAccurateSizeDistribution().getMean();
						bw.write(meanInsert+"\t");
						
						bw.write(insertSizeMap.get(aPeak).getAccurateSizeDistribution().getStandardDeviation()+"\t");
						
						bw.write(((double)fragmentSize-meanInsert)+"\n");

					}
					bws.write(gene.getName()+"\t");
					bws.write(result.get(gene).size()+"\n");
				}
				
			}
			bw.close();
			bws.close();
			
		}
		//TODO
		else{
			
		}
	}
	
	
	
	private Map<RefSeqGene, List<RefSeqGene>> assignPeaks(BEDFileParser peaks, BEDFileParser oldAnnotations) throws IOException{
		Map<RefSeqGene, List<RefSeqGene>> result = new HashMap<RefSeqGene, List<RefSeqGene>>();
		Iterator<String> chrIt = oldAnnotations.getChromosomeIterator();
		
		/*
		 * For each chromosome
		 */
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			/*
			 * If the alignment data has data from that chromosome
			 */
			if(peaks.containChr(chr)){
				logger	.info("Processing " + chr);
				/*
				 * Obtain an iterator over the interval tree values built from the annotations (genes) on that particular chromosome.
				 * IntervalTree of RefSeqGeneWithIsoforms
				 */
				Iterator<RefSeqGeneWithIsoforms> annotation_iter = oldAnnotations.getChrTree(chr).valueIterator();
				/*
				 * Parse annotation tree
				 * While there is an annotated RefSeqGeneWithIsoforms in Interval tree to analyze
				 * Thus, for each gene
				 */
				while(annotation_iter.hasNext()){
					/*
					 * annotation = current node of RefSeqGeneWithIsoforms in interval tree
					 */
					RefSeqGeneWithIsoforms gene = annotation_iter.next();
//					Iterator<RefSeqGene> isoform_iter = annotation_with_isoforms.getAllIsoforms().iterator();

//					while(isoform_iter.hasNext()){
						List<RefSeqGene> associatedPeaks = new ArrayList<RefSeqGene>();
								
//						RefSeqGene gene = isoform_iter.next();
						//Get all peaks in the region X downstream of the gene and Y upstream of its gene end
						
						//Get all peaks inside this gene
						IntervalTree<RefSeqGeneWithIsoforms> overlappersTree = peaks.getOverlappers(gene);
						
						if(!overlappersTree.isEmpty()){
							Iterator<RefSeqGeneWithIsoforms> overlappersIter = overlappersTree.valueIterator();
							while(overlappersIter.hasNext()){
								RefSeqGene overlapper = overlappersIter.next();
								//TODO: check if the peak has already been assigned to another gene which is not an isoform
								if(peakPassesTests(gene,overlapper)){
									associatedPeaks.add(overlapper);
								}
							}
						}

						//All peaks downstream of the gene
						int extend = WINDOW;
						while(extend < max3PExtension) {
							Alignments end = null;
							if(gene.getOrientation().equals("-"))
								end = new Alignments(gene.getChr(), gene.getStart() - extend, gene.getStart() - (extend-WINDOW));
							else
								end = new Alignments(gene.getChr(), gene.getEnd() + (extend - WINDOW), gene.getEnd() + extend);
							/*
							 * Get an interval tree for all/any exons that overlap with the extended region
							 */
							IntervalTree<RefSeqGeneWithIsoforms> endOverlappersTree = oldAnnotations.getOverlappers(end);
							/*
							 * If there is an overlap with a gene
							 */
							if(!endOverlappersTree.isEmpty()){
								Iterator<RefSeqGeneWithIsoforms> overlappersIter = endOverlappersTree.valueIterator();
								boolean overlapperIsSameGene = true;
								/*
								 * while the gene is the same gene
								 */
								while(overlappersIter.hasNext() && overlapperIsSameGene){
									RefSeqGene overlapper = overlappersIter.next();
									//compare the end coordiantes of the gene
									if(!(overlapper.getOrientedEnd() == gene.getOrientedEnd())&&(overlapper.getOrientation()==gene.getOrientation())){
										overlapperIsSameGene = false;
									}	
								}
								if(!overlapperIsSameGene)
									break;
								// Because the extended region cannot overlap another annotation
							}
							//No overlap so continue with scoring the region
							//Get all peaks inside the alignment
							overlappersTree = peaks.getOverlappers(end);
							if(!overlappersTree.isEmpty()){
								Iterator<RefSeqGeneWithIsoforms> overlappersIter = overlappersTree.valueIterator();
								while(overlappersIter.hasNext()){
									RefSeqGene overlapper = overlappersIter.next();
									//TODO: check if the peak has already been assigned to another gene which is not an isoform
									if(peakPassesTests(gene,overlapper) && !associatedPeaks.contains(overlapper)){
										associatedPeaks.add(overlapper);
									}
								}
							}
							
							extend += (STEP);
						}
						result.put(gene, associatedPeaks);
//					}
				}
			}
		}
	
		return result;
	}
	
	private boolean peakPassesTests(RefSeqGene gene, RefSeqGene peak){
		
		if(peak.getOrientation().equals(gene.getOrientation())){
			return true;
		}
		else{
			return false;
		}
	}
		
	public static void main(String[] args) throws IOException{
		
	/*	Globals.setHeadless(true);
		System.out.println("Using Version R4.4");
		logger.debug("DEBUG ON");
		
		GenericAlignmentDataModel alignments=new GenericAlignmentDataModel("temp.bam", "dm3.sizes", 5);
		alignments.setPositiveStranded();
		//alignments.setSecondRead();
		AlignmentDataModelStats alignmentData = new AlignmentDataModelStats(alignments);
		ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(alignmentData);
		Alignments al = new Alignments("chr2L",3510283,3517192);
		double XX = data.count(al);
		System.out.println("Positive strand: "+XX);
		
		alignments=new GenericAlignmentDataModel("temp.bam", "dm3.sizes", 5);
		alignments.setNegativeStranded();
		alignmentData = new AlignmentDataModelStats(alignments);
		data = new ContinuousDataAlignmentModel(alignmentData);
		XX = data.count(al);
		System.out.println("Negative strand: "+XX);*/
		ArgumentMap argmap = CLUtil.getParameters(args,usage,"score");
		PolyASeqAnalyzer dummy = new PolyASeqAnalyzer(argmap);
		
	}
}
