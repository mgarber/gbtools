package broad.pda.geneexpression.dge;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.math.MathException;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.annotation.GTFFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.gene.RefSeqGeneWithIsoforms;
import broad.pda.seq.alignment.AlignmentUtils;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.AlignmentDataModelStats;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

public class DGEWindowPeakCaller {
	/*
	 * Output column names
	 */
	static Logger logger = Logger.getLogger(DGEWindowPeakCaller.class.getName());

	static final String usage = "Usage: DGEWindowPeakCaller -task <task name> "+
			"\n\t3P Computes expression of a given annotation set for multiple alignments simultaneously " + 
			"\n\t\t-alignments <Alignments (mapped to genome) listed in a text file one on each line> "+
			"\n\t\t-normalizedOutput <Boolean value. True if user wants an additional matrix file with normalized count values; FALSE by default>"+
			"\n\t\t-conditions <A text file enlisting the biological conditions corresponding to each input alignment file>"+
			
			"\n\t5P: Computes expression of a given annotation set for multiple alignments simultaneously " + 
			"\n\t\t-alignments <Alignments (mapped to genome) listed in  text file one on each line> "+
			"\n\t\t-normalizedOutput <Boolean value. True if user wants an additional matrix file with normalized count values; FALSE by default>"+
			"\n\t\t-conditions <A text file enlisting the biological conditions corresponding to each input alignment file>"+

			"\n\t\t-annotations <Annotation file for which to calculate expression. [BED by default]> "+
			"\n\t\t-window <Window used to score gene. Default is 500bp> "+
			"\n\t\t-maxIntoGene <Maximum distance from annotated 3' end to find the best scoring window within the gene> "+
			"\n\t\t-maxExtension <It will allow the program to search for the best scoring window maxExtension bases from the 3' end of the gene> "+
			"\n\t\t-minMappingQuality <Only use reads mapping quality greater than the specified>" +
			"\n\t\t-weighReadsFlag <If true, the read counts will be penalized for multiple mapping DEFAULT:true>"+
			"\n\t\t-removePCRDuplicatesFlag <If true, PCR duplicates are removed DEFAULT:false>"+
			"\n\t\t-stranded <If true indicates that the library is stranded>"+
			"\n\t\t-reconstruct <If true, will check for reconstructed bed file and reconstruct the ends of annotated genes DEFAULT: FALSE>"+
			"\n\t\t-reconstructions <Bed/GTF file with the reconstructed genes must be supplied IF reconstruct flag is true. 3' & 5' ends will be adjusted based on this data.maxExtension & maxIntoGene will be set to 0.>"+
			"\n\t\t-out <Output file [Defaults to stdout]> "+
			"\n\t\t-scoreFullGene <If true, will return the DGE score as the score over the entire gene (from the best window) DEFAULT: FALSE>"+
			"\n\t\t-collapseIsoforms <If true, the bed file will be collapse with a 40% overlap between isoforms collapsed to single genes. The output will be on genes and on inidividual isoforms>"+
			"\n";

	/*
	 * Indices of important scores
	 */
	static final int COUNT_SCORE = 2;
	static final int RPKM_SCORE = 4;
	static final int ENRICH_SCORE = 1;
	static final int COVERAGE_SCORE = 3;
	static final double MIN_OVERLAP = 0.4;
	/*
	 * We use the counts as the score to compare window scores
	 */
	static final int USE_SCORE = ENRICH_SCORE;
	
	static boolean isPaired;
	static boolean isSecondRead;
	static boolean useBothReads;

	private static String annotationFile;
	private static int maxIntoGene;
	private static int maxExtension;
	private static int minimumMappingQuality;
	private static int window;
	private static boolean weighReadsFlag;
	private static boolean removePCRDuplicatesFlag;
	private static boolean isStranded;
	private static int step;
	static HashMap<RefSeqGene, String> duplicateNameMap;

	public DGEWindowPeakCaller(String[] args) throws IOException, ParseException, MathException {

		/*
		 * Gives a log4j error. Check later.
		 */
		Globals.setHeadless(true);

		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"3P");

		/*
		 * Parameters if not provided are set to defaults:
		 * maxIntoGene : default 0
		 * maxExtension : default 0
		 * window: default 500
		 * minMappingQuality : default -1
		 */
		maxIntoGene = argMap.isPresent("maxIntoGene")? argMap.getInteger("maxIntoGene") : 0;
		maxExtension = argMap.isPresent("maxExtension")? argMap.getInteger("maxExtension") : 0;
		minimumMappingQuality = argMap.isPresent("minMappingQuality")? argMap.getInteger("minMappingQuality") : -1;
		window = argMap.isPresent("window")? argMap.getInteger("window"): 500;
		step = argMap.isPresent("step")? argMap.getInteger("step"): 50;
		//TODO
//		collapseIsoforms = argMap.isPresent("collapseIsoforms")? (Boolean.parseBoolean(argMap.get("collapseIsoforms"))): false;
		/*
		 * FLAG for WEIGHING READS BY NH FLAG
		 * TRUE by default
		 * Convert string to boolean
		 */
		weighReadsFlag = argMap.isPresent("weighReadsFlag")? (Boolean.parseBoolean(argMap.get("weighReadsFlag"))): true;
		
		/*
		 * FLAG for REMOVING PCR DUPLICATES 
		 * FALSE by default
		 * Convert string to boolean
 		 */
		removePCRDuplicatesFlag = argMap.isPresent("removePCRDuplicatesFlag")? (Boolean.parseBoolean(argMap.get("removePCRDuplicatesFlag"))): false;

		isStranded = argMap.isPresent("stranded")? (Boolean.parseBoolean(argMap.get("stranded"))): false;
		
		/*
		 * Read the annotation file
		 */
		annotationFile = argMap.getMandatory("annotations");		
		/*
		 * Check the format of the annotation file and call the GTF or BED parser accordingly
		 */
		BEDFileParser annotationParser =  annotationFile.endsWith(".gtf") || annotationFile.endsWith(".GTF")? new GTFFileParser(annotationFile) : new BEDFileParser(annotationFile);		
		
		/*
		 * HashMap of gene to rowName
		 */
		duplicateNameMap = new HashMap<RefSeqGene, String>();
		
		
		logger.info("Using window: " + window + " max 5' extension: " + maxExtension + " maximum  premature start? "+ maxIntoGene + " minimum mapping quality "+ minimumMappingQuality);

		/*
		 * Read the names of the alignment file into an array
		 */
		BufferedReader br = new BufferedReader(new FileReader(argMap.getMandatory("alignments")));
		Map<String,String> alignmentFiles = new HashMap<String,String>();
		Map<String,Collection<String>> conditionMaps = new HashMap<String,Collection<String>>();
		String s;
		while((s = br.readLine())!= null){
			//Map of sample name to alignment File
			alignmentFiles.put(s.split(whitespaceDelimiter)[0],s.split(whitespaceDelimiter)[1]);
			String cond = s.split(whitespaceDelimiter)[2];
			if(!conditionMaps.containsKey(cond) ) {
				conditionMaps.put(cond, new ArrayList<String>());
			}
			Collection<String> groupSamples = conditionMaps.get(cond);
			groupSamples.add(s.split(whitespaceDelimiter)[0]);
		}
		
		isPaired = argMap.containsKey("paired");
		isSecondRead = true;
		useBothReads = true;
		if(isPaired){
			useBothReads = argMap.isPresent("useBothReads")? (Boolean.parseBoolean(argMap.get("useBothReads"))): true;
			if(!useBothReads)
				isSecondRead = argMap.isPresent("isSecondRead")? (Boolean.parseBoolean(argMap.get("isSecondRead"))): true;
		}
		
		if("3P".equalsIgnoreCase(argMap.getTask())){
			
			//score3PMultiple(annotationParser,alignmentFiles,argMap.getOutput(),conditionMaps);
			
			
		} else if("5P".equalsIgnoreCase(argMap.getTask())){			
			
			score5PMultiple(annotationParser,alignmentFiles,argMap.getOutput(),conditionMaps);
			
		}
	}

	
	public static RefSeqGene getSubAnnotationFromEnd(RefSeqGene annotation, int length, int distanceFromTranscriptEnd) {
		int annotationLength = annotation.getTranscriptLength();
		RefSeqGene subAnnotation = null;
		int relativeStart = 0;
		int relativeEnd   = 0;
		if(annotationLength > length + distanceFromTranscriptEnd) {
			if("+".equals(annotation.getOrientation())) {
				relativeStart = annotationLength - (distanceFromTranscriptEnd + length);
				relativeEnd   = annotationLength - distanceFromTranscriptEnd;
			} else {
				relativeStart = distanceFromTranscriptEnd;
				relativeEnd   = distanceFromTranscriptEnd + length;

			}
//			System.err.println("getSubannotation for " + annotation.getName() + " Annotation_length " + annotationLength + "  Window_size " + length + "("+annotation.getOrientation()+") distFromEnd: " + distanceFromTranscriptEnd + " relstart-relend " + relativeStart+"-"+relativeEnd);
			subAnnotation = annotation.trim(relativeStart, relativeEnd);
			if(subAnnotation!= null) {
				subAnnotation.setChromosome(annotation.getChr());
			}
		}
		return subAnnotation;
	}
	
	public static RefSeqGene getSubAnnotationFromStart(RefSeqGene annotation, int length, int distanceFromTranscriptStart) {
		int annotationLength = annotation.getTranscriptLength();
		RefSeqGene subAnnotation = null;
		int relativeStart = 0;
		int relativeEnd   = 0;
		if(annotationLength > length + distanceFromTranscriptStart) {
			if("-".equals(annotation.getOrientation())) {
				relativeStart = annotationLength - (distanceFromTranscriptStart + length);
				relativeEnd   = annotationLength - distanceFromTranscriptStart;
			} else {
				relativeStart = distanceFromTranscriptStart;
				relativeEnd   = distanceFromTranscriptStart + length;

			}
			//System.err.println("getSubannotation for " + annotation.toUCSC() + " Annotation_length " + annotationLength + "  Window_size " + length + "("+annotation.getOrientation()+") distFromEnd: " + distanceFromTranscriptEnd + " relstart-relend " + relativeStart+"-"+relativeEnd);
			subAnnotation = annotation.trim(relativeStart, relativeEnd);
			if(subAnnotation!= null) {
				subAnnotation.setChromosome(annotation.getChr());
			}
		}
		return subAnnotation;
	}
	
	
	private static void score5PMultiple(BEDFileParser annotationCollection,Map<String,String> alignmentFiles,String outputFile,Map<String,Collection<String>> conditionMaps) throws IOException, MathException{
			
		/*
		 * The columns of the matrix will be each alignment file name
		 */
		List<String> cols = new ArrayList<String>();
		for(String group:conditionMaps.keySet()){
			cols.add(group);
		}
/*		for(String name: conditionMaps.keySet()){
			for(String sample: conditionMaps.get(name))
				cols.add(sample);
		}
	*/	
		// A separate bed file for each condition with the best peak
		Map<String,FileWriter> fw=new HashMap<String,FileWriter>();
		for(String group:conditionMaps.keySet()){
			FileWriter f = new FileWriter(outputFile+"."+group+".dge.peaks.bed");
			fw.put(group, f);
		}
		
		//		For each pair of conditions, a list of genes with alternative promoter usage
		MatrixWithHeaders resultMatrix = new MatrixWithHeaders(cols,cols);
		for(String r:resultMatrix.getRowNames()){
			for(String c:resultMatrix.getColumnNames()){
				resultMatrix.set(r, c, 0.0);
			}	
		}
		//For peaks in one condition not the other
		MatrixWithHeaders resultMatrix2 = new MatrixWithHeaders(cols,cols);
		for(String r:resultMatrix2.getRowNames()){
			for(String c:resultMatrix2.getColumnNames()){
				resultMatrix2.set(r, c, 0.0);
			}	
		}
		FileWriter ggTab = new FileWriter(outputFile+".TSS.switch.txt");
		FileWriter ggTab2 = new FileWriter(outputFile+".summary.txt");
		FileWriter ggMat = new FileWriter(outputFile+".matrix.txt");
		FileWriter ggScores = new FileWriter(outputFile+".scores.txt");
		FileWriter bw = new FileWriter(outputFile);
		
		/*
		 * Iterate through list of annotations. Iterator is over chromosome names.
		 */
		Iterator<String> chrom_iter = annotationCollection.getChromosomeIterator();
		if(!isStranded){
		
		/*
		 * Initialize the data models for all alignment files
		 * @param: <alignment flieName> <load_chromosome_stats> <minMappingQuality> <remove_PCR_duplicates> <weigh reads by NH flag>
		 */
		Map<String,ContinuousDataAlignmentModel> libDataModels = new HashMap<String,ContinuousDataAlignmentModel>();
		for(String ss: alignmentFiles.keySet()){
			ContinuousDataAlignmentModel d = AlignmentUtils.loadAlignmentData(alignmentFiles.get(ss),false,minimumMappingQuality,removePCRDuplicatesFlag,weighReadsFlag);

			//TODO
/*			if(isPaired & !useBothReads){
				if(isSecondRead){
					alignments.setSecondRead();
				}
				else{
					alignments.setFirstRead();
				}
			}		
			libDataModels.put(ss, d);
	*/		
		}	
		/*
		 * For each chromosome
		 */
		while(chrom_iter.hasNext()) {
			String chr = chrom_iter.next();
			System.out.println(chr);
			boolean dataForChr = false;
			for(ContinuousDataAlignmentModel model:libDataModels.values()){
				//If any alignment file has data for the chromosome, 
				//		set flag to true and exit loop
				if(model.hasDataForChromosome(chr)){
					dataForChr = true;
					break;
				}
			}
			/*
			 * If the any alignment file has data from that chromosome
			 */
			if(dataForChr){
				logger.info("Processing " + chr);
				/*
				 * Obtain an iterator over the interval tree values built from the annotations (genes) on that particular chromosome.
				 * IntervalTree of RefSeqGeneWithIsoforms
				 */
				Iterator<RefSeqGeneWithIsoforms> annotation_iter = annotationCollection.getChrTree(chr).valueIterator();
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

						Map<String,RefSeqGene> groupMap = new HashMap<String,RefSeqGene>();
						Map<RefSeqGene, Map<String, Double>> peakMap = new HashMap<RefSeqGene, Map<String, Double>>();
						
						//FOR EACH CONDITION GROUP
						for(String group:conditionMaps.keySet()){
							Map<String,ContinuousDataAlignmentModel> models = new HashMap<String,ContinuousDataAlignmentModel>();
							for(String name:conditionMaps.get(group)){
								models.put(name,libDataModels.get(name));
							}
							//TODO:correct the data structure. Very bad
							Map<RefSeqGene, Map<String, Double>> peakMap1 = getMostEnriched5PWindowForModels(annotation, models, annotationCollection);
							if(peakMap1.keySet().size()>1){
								logger.warn("Peak map has more than one peak. Check");
							}
							peakMap.putAll(peakMap1);
							for(RefSeqGene g:peakMap1.keySet()){
								groupMap.put(group, g);
							}
						}
						
						//OUTPUT:
					//	A separate bed file for each condition with the best peak
						for(String group:conditionMaps.keySet()){
							fw.get(group).write(groupMap.get(group).toBED()+"\n");
						}
						
						//Matrix with gene name then peak name then Y/N in each condition
						for(RefSeqGene g:peakMap.keySet()){
							ggTab.write(annotation+"\t"+g.toUCSC()+"\t");
							for(String group:conditionMaps.keySet()){
								if(groupMap.get(group).overlaps(g)){
									ggTab.write("YES"+"\t");
								}
								else{
									ggTab.write("NO"+"\t");
								}
							}
							ggTab.write("\n");
						}
						
						//	For each pair of conditions, a list of genes with alternative promoter usage
						for(String group1:conditionMaps.keySet()){
							for(String group2:conditionMaps.keySet()){
								if(!group1.equals(group2)){
									if(groupMap.get(group1).overlaps(groupMap.get(group2))){
										ggTab2.write(group1+"_"+group2+"\t"+annotation+"\n");
										double val = (resultMatrix.get(group1, group2))+1.0;
										resultMatrix.set(group1, group2,val);
									}
								}
							}
						}
						
					//	A matrix with conditions across conditions showing # of genes with alternative promoter usage	
						
					//	A list file showing gene to peak (all conditions merged)
						
					}
				}
			}
		}
		}
		else{
			String sizes = null;
			logger.info("Scoring stranded reads: ");
			
			Map<String,ContinuousDataAlignmentModel> libDataModelsP = new HashMap<String,ContinuousDataAlignmentModel>();
			Map<String,ContinuousDataAlignmentModel> libDataModelsN = new HashMap<String,ContinuousDataAlignmentModel>();
			AlignmentDataModel[] alignmentsP = new AlignmentDataModel[alignmentFiles.size()];
			AlignmentDataModelStats[] alignmentDataP = new AlignmentDataModelStats[alignmentFiles.size()];
			AlignmentDataModel[] alignmentsN = new AlignmentDataModel[alignmentFiles.size()];
			AlignmentDataModelStats[] alignmentDataN = new AlignmentDataModelStats[alignmentFiles.size()];
			
			int k=0;
			for(String ss: alignmentFiles.keySet()){
				alignmentsP[k]=new GenericAlignmentDataModel(alignmentFiles.get(ss), sizes, true, minimumMappingQuality,removePCRDuplicatesFlag,weighReadsFlag);
				alignmentsP[k].setPositiveStranded();
				alignmentDataP[k] = new AlignmentDataModelStats(alignmentsP[k]);
				libDataModelsP.put(ss, new ContinuousDataAlignmentModel(alignmentDataP[k]));
			
				alignmentsN[k]=new GenericAlignmentDataModel(alignmentFiles.get(ss), sizes, true, minimumMappingQuality,removePCRDuplicatesFlag,weighReadsFlag);
				alignmentsN[k].setNegativeStranded();
				alignmentDataN[k] = new AlignmentDataModelStats(alignmentsN[k]);
				libDataModelsN.put(ss, new ContinuousDataAlignmentModel(alignmentDataN[k]));
				if(isPaired & !useBothReads){
					if(isSecondRead){
						alignmentsP[k].setSecondRead();
						alignmentsN[k].setSecondRead();
					}
					else{
						alignmentsP[k].setFirstRead();
						alignmentsN[k].setFirstRead();
					}
				}
				k++;
				
			}
			
			//MAP OF PEAK TO MAP OF GROUP NAME TO SCORE
			Map<RefSeqGene, Map<String, Double>> peakToScoreMap = new HashMap<RefSeqGene, Map<String, Double>>();
			Map<RefSeqGene, Collection<RefSeqGene>> geneToPeakMap = new HashMap<RefSeqGene, Collection<RefSeqGene>>();
			//map of peak to set of conditions it is found in
			Map<RefSeqGene, Collection<String>> peakToConditionMap = new HashMap<RefSeqGene,Collection<String>>();
			//map of gene to a map of conditions to peak
			Map<RefSeqGene, Map<String,RefSeqGene>> geneToConditionToPeakMap = new HashMap<RefSeqGene,Map<String,RefSeqGene>>();

			/*
			 * For each chromosome
			 */
			while(chrom_iter.hasNext()) {
				String chr = chrom_iter.next();
				
				boolean dataForChr = false;
				for(String ss:libDataModelsP.keySet()){
					//If any alignment file has data for the chromosome, 
					//		set flag to true and exit loop
					if(libDataModelsP.get(ss).hasDataForChromosome(chr) || libDataModelsN.get(ss).hasDataForChromosome(chr)){
						dataForChr = true;
						break;
					}
				}
				/*
				 * If the any alignment file has data from that chromosome
				 */
				if(dataForChr){
					logger.info("Processing " + chr);
					/*
					 * Obtain an iterator over the interval tree values built from the annotations (genes) on that particular chromosome.
					 * IntervalTree of RefSeqGeneWithIsoforms
					 */
					Iterator<RefSeqGeneWithIsoforms> annotation_iter = annotationCollection.getChrTree(chr).valueIterator();
					/*
					 * Parse annotation tree
					 * While there is an annotated RefSeqGeneWithIsoforms in Interval tree to analyze
					 * Thus, for each gene
					 */
					while(annotation_iter.hasNext()){
						/*
						 * annotation = current node of RefSeqGeneWithIsoforms in interval tree
						 */
						//FOR EACH ANNOTATION
						RefSeqGeneWithIsoforms annotation = annotation_iter.next();
				//		Iterator<RefSeqGene> isoform_iter = annotation_with_isoforms.getAllIsoforms().iterator();

				//		while(isoform_iter.hasNext()){

					
				//			RefSeqGene annotation = isoform_iter.next();

						//FOR THIS GENE
							
							//FOR EACH CONDITION GROUP
							for(String group:conditionMaps.keySet()){
								Map<String,ContinuousDataAlignmentModel> models = new HashMap<String,ContinuousDataAlignmentModel>();
								for(String name:conditionMaps.get(group)){
									if(annotation.isNegativeStrand()){
										models.put(name,libDataModelsN.get(name));
									}
									else{
										models.put(name,libDataModelsP.get(name));
									}
								}
								//TODO:correct the data structure. Very bad
								Map<RefSeqGene, Map<String, Double>> peakMap1 = getMostEnriched5PWindowForModels(annotation, models, annotationCollection);
								if(peakMap1.keySet().size()>1){
									logger.warn("Peak map has more than one peak. Check");
								}
								RefSeqGene peak = null;
								for(RefSeqGene g:peakMap1.keySet()){
									peak = g;
								}
								
								double avg = getAverageScore(peak,peakMap1);
								if(avg>0.0){
									System.out.println(peak.toUCSC()+" passes with "+avg);
									peak.setBedScore(avg);
									
									//ADD TO GENETOPEAK MAP
									if(geneToPeakMap.containsKey(annotation)){
										//The exact peak has already been called
										if(geneToPeakMap.get(annotation).contains(peak)){
											for(RefSeqGene p: geneToPeakMap.get(annotation)){
												if(p.equals(peak)){
													if(p.getBedScore()<peak.getBedScore()){
														p.setBedScore(peak.getBedScore());
														for(RefSeqGene p1: geneToPeakMap.get(annotation)){
															if(p1.equals(peak)){
																if(p1.getBedScore()!=peak.getBedScore()){
																	System.out.println("Bed score is not scored by reference");
																}
															}
														}
													}
												}
											}
										}
										else{
											geneToPeakMap.get(annotation).add(peak);
										}
									}
									else{
										List<RefSeqGene> peaks = new ArrayList<RefSeqGene>();
										peaks.add(peak);
										geneToPeakMap.put(annotation, peaks);
									}
									//ADD TO PEAK TO CONDITION MAP
									//IF PEAK IS IN CONDITION MAP
									//TODO: CHeck that this overlaps with above 
									if(peakToConditionMap.containsKey(peak)){
										peakToConditionMap.get(peak).add(group);
									}
									else{
										Collection<String> conds = new ArrayList<String>();
										conds.add(group);
										peakToConditionMap.put(peak, conds);
									}
									//ADD TO GROUP TO CONDITION MAP
									if(geneToConditionToPeakMap.containsKey(annotation)){
										geneToConditionToPeakMap.get(annotation).put(group, peak);
									}else{
										Map<String,RefSeqGene> temp = new HashMap<String,RefSeqGene>();
										temp.put(group, peak);
										geneToConditionToPeakMap.put(annotation, temp);
									}
									
									//ADD TO PEAK TO CONDITION TO SCORE MAP
									if(peakToScoreMap.containsKey(peak)){
										peakToScoreMap.get(peak).put(group, avg);
									}
									else{
										Map<String,Double> temp2 = new HashMap<String,Double>();
										temp2.put(group,avg);
										peakToScoreMap.put(peak, temp2);
									}
							}
							}
							//OUTPUT:
						//	A matrix with conditions across conditions showing # of genes with alternative promoter usage	
							
						//	A list file showing gene to peak (all conditions merged)
							
						//}
					}
				}

			}
			//	A separate bed file for each condition with the best peak
			for(RefSeqGene peak:peakToConditionMap.keySet()){
				for(String name:peakToConditionMap.get(peak)){
					fw.get(name).write(peak.toBED()+"\n");
				}
			}
			//Matrix with gene name then peak name then Y/N in each condition
			
			double threshold = 40.0;//calculateThreshold(0.25,peakToScoreMap);
			System.out.println("Threshold: "+threshold);
			
			//FOR EVERY CONDITION
			ggTab.write("Gene\tPeak\t");
			ggScores.write("Gene\tPeak\t");
			for(String condition:conditionMaps.keySet()){
				ggTab.write(condition+"\t");
				ggScores.write(condition+"\t");
			}
			ggTab.write("\n");
			ggScores.write("\n");
			
			/**
			 * ALTERNATIVE TSS ONLY
			 */
			//FOR EVERY GENE
			for(RefSeqGene gene:geneToPeakMap.keySet()){
				
				if(geneToConditionToPeakMap.get(gene).size()>0){
					//FOR EVERY CONDITION
					List<String> visited = new ArrayList<String>();
					for(String group1:geneToConditionToPeakMap.get(gene).keySet()){
						visited.add(group1);
						for(String group2:geneToConditionToPeakMap.get(gene).keySet()){
							if(!group1.equals(group2) && !visited.contains(group2)){
								RefSeqGene p1 = geneToConditionToPeakMap.get(gene).get(group1);
								RefSeqGene p2 = geneToConditionToPeakMap.get(gene).get(group2);
								//Both peaks pass threshold
								if(passesPeakThreshold(p1,peakToScoreMap,threshold) && passesPeakThreshold(p2,peakToScoreMap,threshold)){
									if(!p1.overlaps(p2)){
										bw.write(gene.getName()+"\t"+group1+"\t"+group2+"\t"+p1.toUCSC()+"\t"+p2.toUCSC()+"\n");
									}
								}
							}
						}
					}
					
				}
			}
			//FOR EVERY GENE
			for(RefSeqGene gene:geneToPeakMap.keySet()){
				//FOR EVERY PEAK FOR THAT GENE
				for(RefSeqGene peak:geneToPeakMap.get(gene)){
					//IF PEAK PASSES THRESHOLD IN AT LEAST ONE CONDITION
					if(passesPeakThreshold(peak,peakToScoreMap,threshold)){
						ggTab.write(gene.getName()+"\t"+peak.toUCSC()+"\t");
						ggScores.write(gene.getName()+"\t"+peak.toUCSC()+"\t");
						//FOR EVERY CONDITION
						for(String condition:conditionMaps.keySet()){
							
							//CHECK IF THE GENE HAS ANY PEAKS IN THAT CONDITION
							if(geneToConditionToPeakMap.get(gene).containsKey(condition)){
								RefSeqGene other = geneToConditionToPeakMap.get(gene).get(condition);
								if(peak.overlaps(other)){
									ggTab.write("YES"+"\t");
									
									ggScores.write(peakToScoreMap.get(other).get(condition).toString()+"\t");
								}
								else{
									ggTab.write("NO"+"\t");
									ggScores.write("NO-OVERLAP"+"\t");
								}
							}
							else{
								ggTab.write("N/A"+"\t");
								ggScores.write("N/A"+"\t");
							}
						}
						ggTab.write("\n");
						ggScores.write("\n");
						
						for(String group1:conditionMaps.keySet()){
							for(String group2:conditionMaps.keySet()){
								//PEAK in BOTH CONDITIONS
								if(!group1.equals(group2)){
									//If gene has a peak for both conditions
									if(geneToConditionToPeakMap.get(gene).containsKey(group1) 
											&& geneToConditionToPeakMap.get(gene).containsKey(group2)){
										//If the peaks do not overlap: ALTERNATIVE USAGE
										if(geneToConditionToPeakMap.get(gene).get(group1).overlaps(peak) && (!geneToConditionToPeakMap.get(gene).get(group2).overlaps(peak))){
											ggTab2.write(group1+"_"+group2+"\t"+gene.getName()+"\n");
											System.out.println(geneToConditionToPeakMap.get(gene).get(group1).toUCSC()+"\t"+geneToConditionToPeakMap.get(gene).get(group2).toUCSC());
											double val = (resultMatrix.get(group1, group2))+1.0;
											resultMatrix.set(group1, group2,val);
										}
									}
									//PEAK in ONE GROUP NOT THE OTHER
									else{
										if(geneToConditionToPeakMap.get(gene).containsKey(group1) && geneToConditionToPeakMap.get(gene).get(group1).overlaps(peak)){
											double val = (resultMatrix2.get(group1, group2))+1.0;
											resultMatrix2.set(group1, group2,val);
										}
										else if(geneToConditionToPeakMap.get(gene).containsKey(group2) && geneToConditionToPeakMap.get(gene).get(group2).overlaps(peak)){
											double val = (resultMatrix2.get(group2, group1))+1.0;
											resultMatrix2.set(group2, group1,val);
										}
									}
								}
							}
						}
					}
				}
			}
		}

		for(String group1:resultMatrix.getColumnNames()){
				ggMat.write(group1+"\t");
		}
		ggMat.write("\n");
		for(String group1:resultMatrix.getRowNames()){
			for(String group2:resultMatrix.getColumnNames()){
				ggMat.write(resultMatrix.get(group1, group2)+"\t");
			}
			ggMat.write("\n");
		}
		ggMat.write("\n");
		for(String group1:resultMatrix2.getColumnNames()){
			ggMat.write(group1+"\t");
		}
		ggMat.write("\n");
		for(String group1:resultMatrix2.getRowNames()){
			for(String group2:resultMatrix2.getColumnNames()){
				ggMat.write(resultMatrix2.get(group1, group2)+"\t");
			}
			ggMat.write("\n");
		}
		
		for(String group:conditionMaps.keySet()){
			fw.get(group).close();
		}
		ggTab.close();
		ggTab2.close();
		ggMat.close();
		ggScores.close();
		bw.close();
	}
	
	
	private static double calculateThreshold(double d,
			Map<RefSeqGene, Map<String, Double>> peakToScoreMap) {
		
		List<Double> vals = new ArrayList<Double>();
		
		for(RefSeqGene peak:peakToScoreMap.keySet()){
			for(String name:peakToScoreMap.get(peak).keySet()){
				vals.add(peakToScoreMap.get(peak).get(name));
			}
		}
		
		Collections.sort(vals);
		double thresh = Statistics.quantile(vals, d);
		return thresh;
	}


	private static boolean passesPeakThreshold(RefSeqGene peak,
			Map<RefSeqGene, Map<String, Double>> peakToScoreMap,double threshold) {
		
		for(String name:peakToScoreMap.get(peak).keySet()){
				if(peakToScoreMap.get(peak).get(name)>threshold)
					return true;
			}
		return false;
	}


	private static double getAverageScore(RefSeqGene g,
			Map<RefSeqGene, Map<String, Double>> peakMap) {
		if(!peakMap.containsKey(g)){
			System.out.println("peak map does not contain data for "+g.toUCSC());
			return 0.0;
		}
		double avg = 0.0;
		for(String sample:peakMap.get(g).keySet()){
			//System.out.println(g.getChr()+" "+peakMap.get(g).get(sample));
			avg += peakMap.get(g).get(sample);
		}
		avg = avg/(double)peakMap.get(g).size();
		System.out.println("Average for: "+g.toUCSC()+" "+avg);
		return avg;
	}


	/**
	 * Returns a map of the peak to a map of sample name to score
	 * @param annotation
	 * @param libDataModels
	 * @param annotationCollection
	 * @return
	 * @throws IOException
	 */
	private static Map<RefSeqGene, Map<String,Double>> getMostEnriched5PWindowForModels(RefSeqGene annotation,Map<String,ContinuousDataAlignmentModel> libDataModels,BEDFileParser annotationCollection) throws IOException{

		Map<RefSeqGene, Map<String,Double>> peakMap = new HashMap<RefSeqGene,Map<String,Double>>();
		
		RefSeqGene annotationStart = annotation;					
		int annotationLength = annotation.getTranscriptLength();
		
		//If the length of the annotated transcript > window size being analyzed,
		//get sub annotation for window length
		
		if(annotationLength>window){
			annotationStart = getSubAnnotationFromStart(annotation,window,0);
			if(annotationStart == null){
				logger.warn("Annotation end for " + (annotationLength - window) + "-" + annotationLength + " --> " + annotation.toBED() + " was null.");
			}
		}
		//System.out.println(annotationStart.getSequence());
		
		double bestEnrichment = getEnrichmentForWindow(libDataModels,annotationStart);
		RefSeqGene bestWindow = annotationStart;	 
		
		// Use COUNTS as the score
		// The best score is the window score for now
		
		double tmpEnrichment = 0.0;
		int intoGene = step;
		/*
		 * 							 
		 * If downstream extension is allowed, use sliding windows with overlaps of step
		 * while into region length is smaller than (annotation - window), i.e. it lies within the annotation
		 * and it is at distance less than the max region allowed downstream of start of gene
		 */
		while((intoGene<(annotationLength - window)) && (intoGene<maxIntoGene)){

			/*
			 * get annotation for region of length window, "retreat" length from end of transcript
			 */
			annotationStart = getSubAnnotationFromStart(annotation, window, intoGene);
			/*
			 * calculate score for the window and find max score of windows
			 * Uses scanPRate(RefSeqGene, tree) which adds extension factor in AlignmentDataModelStats.
			 * Difference between scoregene and scanPRate from above is that this does not include exons.
			 */
			if(annotationStart !=null){
				tmpEnrichment = getEnrichmentForWindow(libDataModels,annotationStart);
				if(tmpEnrichment>bestEnrichment){
					bestEnrichment = tmpEnrichment;
					bestWindow = annotationStart;
				}
			}
			intoGene += step;
		}

		int extend = step;
		tmpEnrichment = 0.0;
		while(extend < maxExtension) {
			Alignments start = null;
			if(annotation.getOrientation().equals("+"))
				start = new Alignments(annotation.getChr(), annotation.getStart() - extend, annotation.getStart() - (extend-window));
			else
				start = new Alignments(annotation.getChr(), annotation.getEnd() + (extend - window), annotation.getEnd() + extend);
			/*
			 * Get an interval tree for all/any exons that overlap with the extended region
			 */
			IntervalTree<RefSeqGeneWithIsoforms> endOverlappersTree = annotationCollection.getOverlappers(start);
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
					overlapperIsSameGene = (overlapper.getOrientedStart() == annotation.getOrientedStart());
				}
				if(!overlapperIsSameGene)
					break;
				// Because the extended region cannot overlap another annotation
			}
			//No overlap so continue with scoring the region
			annotationStart = new RefSeqGene(start);
			tmpEnrichment = getEnrichmentForWindow(libDataModels,annotationStart);
			if(tmpEnrichment>bestEnrichment){
				bestEnrichment = tmpEnrichment;
				bestWindow = annotationStart;
			}
			extend += (step);
		}
		
		peakMap.put(bestWindow, getScoresForWindow(libDataModels,bestWindow));
		return peakMap;
	}

	
	/**
	 * This function calculates the best enrichment score for a given window
	 * @param libDataModels
	 * @param annotationEnd
	 * @return Vector of size libDataModels with enrichment of each 
	 * @throws IOException
	 */
	private static double getEnrichmentForWindow(Map<String,ContinuousDataAlignmentModel> libDataModels,RefSeqGene annotationEnd) throws IOException{
		
		//double max = 0.0;
		double sum = 0.0;
		for(ContinuousDataAlignmentModel model:libDataModels.values()){
			if(model.hasDataForChromosome(annotationEnd.getChr())){	
				double score = model.scoreGene(annotationEnd)[ENRICH_SCORE];
				/*if(score>max){
					max = score;
				}*/
				sum += score;
			}
		}
		
		//return max;
		return sum;
	}
	

	/**
	 * This function calculates the scores for all alignment models for a given window
	 * @param libDataModels
	 * @param annotationEnd
	 * @return Vector of size libDataModels with enrichment of each 
	 * @throws IOException
	 */
	private static Map<String,Double> getScoresForWindow(Map<String,ContinuousDataAlignmentModel> libDataModels,RefSeqGene annotationEnd) throws IOException{
		
		Map<String,Double> scores = new HashMap<String,Double>();
		for(String ss:libDataModels.keySet()){
			if(libDataModels.get(ss).hasDataForChromosome(annotationEnd.getChr())){
				double score = libDataModels.get(ss).scoreGene(annotationEnd)[COUNT_SCORE];
				//System.out.println(annotationEnd.toUCSC()+" Score: "+score);
				scores.put(ss, score);
			}
			else{
				scores.put(ss, 0.0);
			}
		}
		
		return scores;
	}
	
	/**
	 * This function calculates the score for a gene in a strand-specific manner
	 * @param libDataModelP
	 * @param libDataModelN
	 * @param annotationEnd
	 * @return Vector of size libDataModels with enrichment of each 
	 * @throws IOException
	 */
	public static double[] scoreStrandedGene(ContinuousDataAlignmentModel libDataModelP,ContinuousDataAlignmentModel libDataModelN,RefSeqGene annotationEnd) throws IOException{
		
		double[] scores;
		if(!annotationEnd.isUnoriented()){
			if(annotationEnd.isNegativeStrand())
				scores = libDataModelN.scanPRate(annotationEnd);					 
			else
				scores = libDataModelP.scanPRate(annotationEnd);
		}
		else{
			
			double[] scoresP = libDataModelP.scanPRate(annotationEnd);
			double[] scoresN = libDataModelN.scanPRate(annotationEnd);
			scores = new double[scoresP.length];
			//Everything is not additive
			//[0]c alculatePVal(new Double(sum).intValue(), getLambda(chr), count, getNumberMarkers(chr)) - take the minimum
			//[1] enrich - add
			//[2] sum - add 
			//[3] avgCoverage - add 
			//[4] rpkm - add
			scores[0] = Math.min(scoresP[0], scoresN[0]);
			for(int i=0;i<scores.length;i++){
				scores[i] = scoresP[i]+scoresN[i];
			}
		}

		return scores;
	}

		
	/**
	 * This function converts an array of double to List of type Double
	 * @param double[]
	 * @return
	 */
	public static List<Double> array2List(double[] array){
		List<Double> lt=new ArrayList<Double>();

		for(int i=0;i<array.length;i++){
			lt.add(array[i]);
		}
	
		return lt;
	}
	
	
	public static String whitespaceDelimiter = "\\s++"; //$NON-NLS-1$
	
	public static void main (String [] args) throws MathException, ParseException, IOException {
		
		//AlignmentDataModel alignments=new GenericAlignmentDataModel("DC_6hr_1_dge.th.1.4.1.sorted.bam", "sizes", 5);
		DGEWindowPeakCaller d = new DGEWindowPeakCaller(args);
		//MatrixWithHeaders m = new MatrixWithHeaders("t1");
		//writeNormalizedMatrix("t1.out", m);
		//DGE dummy = new DGE(args);

	}
}