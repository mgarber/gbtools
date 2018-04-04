package broad.pda.geneexpression.dge;

import java.io.BufferedReader;
import java.io.BufferedWriter;
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

import broad.core.annotation.BED;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
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

public class DGE {
	/*
	 * Output column names
	 */
	static final String ANNOTATED_END_EXPR_COL = "AnnotatedEndExpression";
	static final String UPSTREAM_EXPR_COL = "BestUpstreamExpression";
	static final String DOWNSTREAM_EXPR_COL = "BestDownstreamExpression";
	static final String BEST_EXPR_COL = "Expression";
	static final String BEST_PVAL_COL = "p-value";
	//static final String THREE_P_PEXPRESSION_COL = "3PEndExpression";
	static final String BEST_DIST_TO_END_COL = "expressionWindowDistToEnd";
	static Logger logger = Logger.getLogger(DGE.class.getName());

	static final String usage = "Usage: DGE -task <task name> "+
			"\n\tTASK 1: score3P: Computes 3' expression of a given annotation set" + 
			"\n\tTASK 2: score5P: Computes 5' expression of a given annotation set " + 
			"\n**************************************************************"+
			"\n\t\tMANDATORY arguments"+
			"\n**************************************************************"+
			"\n\n\t\t-alignment <Alignment (mapped to genome) for which expression must be calculated. Index file MUST be provided.> "+
			"\n\t\t\t OR"+
			"\n\t\t-alignments <Three column text file listing, one line per sample, for which expression must be calculated.Index files MUST be provided." +
			"\n\t\t\t Sample_Name\t Alignment_file\t Grouping_for_sample (Replicates of a condition have the same grouping)> "+
			"\n\n\t\t-annotations <Annotation file for which to calculate expression. [BED by default]> "+

			"\n\n**************************************************************"+
			"\n\t\tOPTIONAL arguments"+
			"\n**************************************************************"+
			"\\nn\t\t-window <Window used to score gene. Default is 500bp> "+
			"\n\t\t-maxIntoGene <Maximum distance from annotated 3' end to find the best scoring window within the gene> "+
			"\n\t\t-maxExtension <It will allow the program to search for the best scoring window maxExtension bases from the 3' end of the gene> "+
			"\n\t\t-minMappingQuality <Only use reads mapping quality greater than the specified>" +
			"\n\t\t-dontWeighReadsFlag <If provided the read counts will NOT be penalized for multiple mapping >"+
			"\n\t\t-removePCRDuplicatesFlag <If provided, PCR duplicates are removed DEFAULT:false>"+
			"\n\t\t-stranded <If provided indicates that the library is stranded. Default is unstranded.>"+
			"\n\t\t-oppositeStrand <If provided indicates that the read in the opposite direction of transcription should be used for quantification. Default is same strand.>"+
			"\n\t\t-reconstructions <If provided a Bed/GTF file with the reconstructed genes is supplied, 3' & 5' ends will be adjusted based on this data.maxExtension & maxIntoGene will be set to 0.>"+
			"\n\t\t-dontCollapseIsoforms <By default, the bed file will be collapse with a 40% overlap between isoforms collapsed to single genes. The output will be on genes and on inidividual isoforms. " +
			"If this flag is provided, this will NOT be done.>"+
			"\n\t\t-out <Output file [Defaults to stdout]> "+
			
			"\n\n**************************************************************"+
			"\n\t\tArguments specific to multiple file run"+
			"\n**************************************************************"+
			"\n\n\t\t-normalizedOutput <If present, an additional matrix file with normalized count values; FALSE by default>"+
			"\n\t\t-scoreFullGene <If present, will return the DGE score as the score over the entire gene (from the best window) DEFAULT: FALSE>"+
			
			"\n\n\tTASK 3: normalizeMatrix: Takes an input DGE matrix, and outputs a normalized matrix." + 
			"\n\t\t-in <Input matrix> "+
			"\n\t\t-out <Normalized matrix> "+

			"\n";

	/*
	 * Indices of important scores
	 */
	static final int COUNT_SCORE = 2;
	static final int PVAL_SCORE = 0;
	static final int RPKM_SCORE = 4;
	static final int STEP = 50;
	static final double MIN_OVERLAP = 0.4;
	/*
	 * We use the counts as the score to compare window scores
	 */
	static final int USE_SCORE = COUNT_SCORE;
	
	private static String annotationFile;
	private static int maxIntoGene;
	private static int maxExtension;
	private static int minimumMappingQuality;
	private static int window;
	private static boolean weighReadsFlag;
	private static boolean removePCRDuplicatesFlag;
	private static boolean fullGeneScoreFlag;
	private static boolean isStranded;
	private static boolean collapseIsoforms;
	private static boolean oppositeStrand;
	static HashMap<RefSeqGene, String> duplicateNameMap;

	public DGE(String[] args) throws IOException, ParseException, MathException {

		/*
		 * Gives a log4j error. Check later.
		 */
		Globals.setHeadless(true);
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"score3P");

		if ("normalizeMatrix".equalsIgnoreCase(argMap.getTask())) {
			
			writeNormalizedMatrix(argMap.getOutput(),new MatrixWithHeaders(argMap.getInput()));
		}
		else{
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
		collapseIsoforms = !argMap.isPresent("dontCollapseIsoforms");
		oppositeStrand = argMap.isPresent("oppositeStrand");
		/*
		 * FLAG for WEIGHING READS BY NH FLAG
		 * TRUE by default
		 * Convert string to boolean
		 */
		weighReadsFlag = !argMap.isPresent("dontWeighReadsFlag");
		
		/*
		 * FLAG for REMOVING PCR DUPLICATES 
		 * FALSE by default
		 * Convert string to boolean
 		 */
		removePCRDuplicatesFlag = argMap.isPresent("removePCRDuplicatesFlag");
		/*
		 * FLAG for TO RETURN SCORE OVER FULL GENE 
		 * FALSE by default
		 * Convert string to boolean
		 */
		fullGeneScoreFlag = argMap.isPresent("scoreFullGene");

		isStranded = argMap.isPresent("stranded");
		
		//IF user wants the normalized matrix,
		/*
		 * FLAG to return a normalized matrix. FALSE by default. Convert string to boolean
		 */
		boolean normalizedOutput = argMap.isPresent("normalizedOutput");

		/*
		 * Read the annotation file
		 */
		annotationFile = argMap.getMandatory("annotations");		
		/*
		 * Check the format of the annotation file and call the GTF or BED parser accordingly
		 */
		BEDFileParser annotationParser =  annotationFile.endsWith(".gtf") || annotationFile.endsWith(".GTF")? new GTFFileParser(annotationFile) : new BEDFileParser(annotationFile);		

		/**
		 * If user wants to reconstruct gene ends.
		 */
		if(argMap.isPresent("reconstructions")) {
			
			String reconstructionFile = argMap.getMandatory("reconstructions");
			BEDFileParser reconstructions =  reconstructionFile.endsWith(".gtf") || reconstructionFile.endsWith(".GTF")? new GTFFileParser(reconstructionFile) : new BEDFileParser(reconstructionFile);
			logger.info("A reconstructed bed file provided, going to adjust ends of genes and ignore 3' extensions");

			annotationParser = reconstructGeneEnds(reconstructions, annotationParser, annotationFile);
			maxIntoGene = 0;
			maxExtension = 0;
		}
		
		/*
		 * HashMap of gene to rowName
		 */
		duplicateNameMap = new HashMap<RefSeqGene, String>();
		
		
		logger.info("Using window: " + window + " max 5' extension: " + maxExtension + " maximum  premature start? "+ maxIntoGene + " minimum mapping quality "+ minimumMappingQuality);

		/*
		 * If the task assigned was "score"
		 */
		if ("score3P".equalsIgnoreCase(argMap.getTask())) {
			/*
			 * If "alignment" is present, it is single file
			 */
			if(argMap.isPresent("alignment")){
				String alignmentFile = argMap.getMandatory("alignment");
				score3P(annotationParser,alignmentFile,argMap.getOutput());
			}
			/*
			 * if "alignments" is present, it is multiple run
			 */
			 else if(argMap.isPresent("alignments")){
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
					
					score3PMultiple(annotationParser,alignmentFiles,argMap.getOutput(),normalizedOutput,conditionMaps);
					
					
				} 
			 else{
				 throw new IllegalArgumentException("Argument alignment or alignments must be provided\n"+usage);
			 }
			
		}  else if("score5P".equalsIgnoreCase(argMap.getTask())){
			/*
			 * If "alignment" is present, it is single file
			 */
			if(argMap.isPresent("alignment")){
				String alignmentFile = argMap.getMandatory("alignment");
				score5P(annotationParser,alignmentFile,argMap.getOutput());
				
			}
			/*
			 * if "alignments" is present, it is multiple run
			 */
			else if(argMap.isPresent("alignments")){
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
				
				score5PMultiple(annotationParser,alignmentFiles,argMap.getOutput(),normalizedOutput,conditionMaps);	
			}
			 else{
				 throw new IllegalArgumentException("Argument alignment or alignments must be provided\n"+usage);
			 }
		} 
	}
	}

	
	public static RefSeqGene getSubAnnotationFromEnd(RefSeqGene annotation, int length, int distanceFromTranscriptEnd) {
		int annotationLength = annotation.getTranscriptLength();
		RefSeqGene subAnnotation = null;
		int relativeStart = 0;
		int relativeEnd   = 0;
		if(annotationLength >= length + distanceFromTranscriptEnd) {
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
				subAnnotation.setOrientation(annotation.getOrientation());
			}
		}
		else{
			System.err.println("Annotation length "+annotationLength+" is less than length "+length +" distance from end "+ distanceFromTranscriptEnd);
		}
		return subAnnotation;
	}
	
	public static RefSeqGene getSubAnnotationFromStart(RefSeqGene annotation, int length, int distanceFromTranscriptStart) {
		int annotationLength = annotation.getTranscriptLength();
		RefSeqGene subAnnotation = null;
		int relativeStart = 0;
		int relativeEnd   = 0;
		if(annotationLength >= length + distanceFromTranscriptStart) {
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
				subAnnotation.setOrientation(annotation.getOrientation());
			}
		}
		return subAnnotation;
	}
	
	/**
	 * This function writes the output to two files
	 * @param outFileName Name of the output file provided by user
	 * @param resultMatrix matrix to be written to the output file
	 * @param allGenes List of all genes in annotation file
	 * @throws IOException
	 */
	private static void writeOutputFiles(String outFileName,MatrixWithHeaders resultMatrix, List<RefSeqGene> allGenes, HashMap<RefSeqGene, String> duplicateNameMap) throws IOException{
		
		// Name of the bed file is outFileName followed by .bed
		String bedFileName = outFileName+".bed";
		
		/*
		 * Write the expression file
		 */
		BufferedWriter outBw = new BufferedWriter(new FileWriter(outFileName));
		resultMatrix.write(outBw);
		outBw.close();
		/*
		 * Write the bed File
		 */
		BufferedWriter bedBw = new BufferedWriter(new FileWriter(bedFileName));
		for(RefSeqGene gene: allGenes){
			bedBw.write(gene.toBED());
			bedBw.append("\t"+duplicateNameMap.get(gene));
			bedBw.newLine();
		}
		bedBw.close();
		
	}
	
	/**
	 * This function writes the output to three files
	 * @param outFileName Name of the output file provided by user
	 * @param resultMatrix matrix to be written to the output file
	 * @param allGenes List of all genes in annotation file
	 * @throws IOException
	 */
	private static void writeOutputFiles(String outFileName,MatrixWithHeaders resultMatrix, List<RefSeqGene> allGenes, HashMap<RefSeqGene, String> duplicateNameMap,List<BED> windows) throws IOException{
		
		writeOutputFiles(outFileName,resultMatrix,allGenes,duplicateNameMap);
		// Name of the second bed file is outFileName followed by .windows.bed
		String bedFileName = outFileName+".windows.bed";
		
		/*
		 * Write the bed File
		 */
		BufferedWriter bedBw = new BufferedWriter(new FileWriter(bedFileName));
		for(BED gene: windows){
			bedBw.write(gene.toString());
			//bedBw.append("\t"+duplicateNameMap.get(gene));
			bedBw.newLine();
		}
		bedBw.close();
		
	}
	
	public static void writeNormalizedMatrix(String fileName, MatrixWithHeaders resultMatrix) throws IOException{
		
		Map<String,Double> factors = calculateNormalizationFactors(resultMatrix);
		MatrixWithHeaders normalizedMatrix = resultMatrix.multiplyColumnsWithConstants(factors);
		String normalizedFileName = fileName+".normalized";
		String normalizedFactorFile = fileName+".coverage";
		BufferedWriter outBw = new BufferedWriter(new FileWriter(normalizedFileName));
		normalizedMatrix.write(outBw);
		outBw.close();
		
		BufferedWriter fBw = new BufferedWriter(new FileWriter(normalizedFactorFile));
		for(String ss:factors.keySet()){
			fBw.write(ss+"\t"+((double)1.0/factors.get(ss))+"\n");
		}
		fBw.close();
	}
	
	
	private static void score3P(BEDFileParser annotationCollection, String alignmentFile,String outputFile) throws IOException, ParseException, MathException {
		
		/*
		 * Restrict the annotations to one per gene. No duplicate annotations.
		 */
		List<RefSeqGene> allGenes = annotationCollection.GetGenes();
		
		/*
		 * Map of each gene to its peak.
		 */
		Map<String,RefSeqGene> geneToWindowMap = new HashMap<String,RefSeqGene>();
		
		/* each row of the matrix is one RefSeqGene name*/
		List<String> rows = new ArrayList<String>();
		/*
		 * HashMap of gene name to the number of duplicates
		 */
		HashMap<String, Integer> duplicateMap = new HashMap<String, Integer>();
		
		int duplicates=0;
		for(RefSeqGene gene : allGenes) {
			if(!rows.contains(gene.getName())) {
				String name = gene.getName();
				if(duplicateNameMap.containsKey(gene)){
					logger.info("Entry for "+name+" already exists");
				}
				else{
					rows.add(name);
					duplicateMap.put(name, 1);
					duplicateNameMap.put(gene, name);
				}
			} 
			// If the gene name has another annotation
			else {
				
				if(duplicateNameMap.containsKey(gene)){
					logger.info("Entry for "+gene.getName()+" already exists in "+duplicateNameMap.get(gene));
				}
				else{
					//Row name is now the geneName appended with the duplicate number
					duplicateMap.put(gene.getName(), (duplicateMap.get(gene.getName())+1));
					String name = (gene.getName()+"_"+duplicateMap.get(gene.getName()));
					rows.add(name);
					duplicateNameMap.put(gene, name);
					//logger.warn("Duplicated annotation : "+ gene.toBED());
					duplicates++;
				}
			}
		}
		logger.info("Found " + duplicates + " duplicates, ignoring them, going to process " + rows.size() + " annotations");
		
		/*
		 * Initialize the lists for output of expression score
		 * These lists names will also be the columns of the matrix
		 */
		List<String> cols = new ArrayList<String>();
		cols.add(BEST_EXPR_COL);
		cols.add(BEST_PVAL_COL);
		cols.add(ANNOTATED_END_EXPR_COL);
		if(maxIntoGene>0 || maxExtension>0){

			if(maxIntoGene>0)
				cols.add(UPSTREAM_EXPR_COL);
			if(maxExtension>0)
				cols.add(DOWNSTREAM_EXPR_COL);

			cols.add(BEST_DIST_TO_END_COL);
		}

		/* 
		 * Initialize the matrix with row and column names
		 */
		MatrixWithHeaders resultMatrix = new MatrixWithHeaders(rows,cols);

		/*
		 * Iterate through list of annotations. Iterator is over chromosome names.
		 */
		Iterator<String> chrom_iter = annotationCollection.getChromosomeIterator();

		if(!isStranded){
			
		logger.info("Scoring using all reads:");
		/*
		 * Initialize the data model using the alignment file
		 * @param: <alignment flieName> <load_chromosome_stats> <minMappingQuality> <remove_PCR_duplicates> <weigh reads by NH flag>
		 */
		ContinuousDataAlignmentModel libDataModel = AlignmentUtils.loadAlignmentData(alignmentFile,true,minimumMappingQuality,removePCRDuplicatesFlag,weighReadsFlag);
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
						RefSeqGene annotationEnd = annotation;					
						int annotationLength = annotation.getTranscriptLength();
						
						/*
						 * If the length of the annotated transcript > window size being analyzed,
						 * get sub annotation for window length
						 */
						if(annotationLength>window){
							annotationEnd = getSubAnnotationFromEnd(annotation,window,0);
							if(annotationEnd == null){
								logger.warn("Annotation end for " + (annotationLength - window) + "-" + annotationLength + " --> " + annotation.toBED() + " was null.");
							}
						}
						//System.out.println(annotationEnd.getSequence());
						/*
						 * Calculate 
						 * [0] p-value
						 * [1] enrichment score
						 * [2] # aligned reads
						 * [3] average coverage
						 * [4] RPKM value
						 * [5] local lambda
						 * [6] gene length
						 * [7] poisson(#aligned reads)
						 * for the annotated region 	
						 * Finally calls scanPRate(RefSeqGene gene, int extensionFactor) in AlignmentDataModelStats					 
						 */
						double [] bestScores = libDataModel.scanPRate(annotationEnd);					 
						/*
						 * Use COUNTS as the score
						 * The best score is the window score for now
						 */
						double annotatedEndExpr = bestScores[USE_SCORE];
						int bestScoreDistanceFromEnd = 0;
						geneToWindowMap.put(duplicateNameMap.get(annotation), annotationEnd);

						double upstreamExpr = 0;
						double downstreamExpr = 0;

						double [] tmpScores = null;
						int retreat = STEP;
						/*
						 * If upstream extension is allowed, use sliding windows with overlaps of STEP
						 * while retreat region length is smaller than (annotation - window), i.e. it lies within the annotation
						 * and it is at distance less than the max region allowed upstream of end of gene
						 */
						while((retreat<(annotationLength - window)) && (retreat<maxIntoGene)){

							/*
							 * get annotation for region of length window, "retreat" length from end of transcript
							 */
							annotationEnd = getSubAnnotationFromEnd(annotation, window, retreat);
							/*
							 * calculate score for the window and find max score of windows
							 * Uses scanPRate(RefSeqGene, tree) which adds extension factor in AlignmentDataModelStats.
							 * Difference between scoregene and scanPRate from above is that this does not include exons.
							 */
							if(annotationEnd !=null){
								tmpScores = libDataModel.scoreGene(annotationEnd);
								if(tmpScores[USE_SCORE]>bestScores[USE_SCORE]){
									for(int i=0;i<bestScores.length;i++){
										bestScores[i] = tmpScores[i];
									}
									bestScoreDistanceFromEnd = -retreat;
									geneToWindowMap.put(duplicateNameMap.get(annotation), annotationEnd);
								}
								if(tmpScores[USE_SCORE]>upstreamExpr)
									upstreamExpr = tmpScores[USE_SCORE];
							}
							retreat += STEP;
						}

						int extend = STEP;
						tmpScores = null;
						while(extend < maxExtension) {
							Alignments end = null;
							if(annotation.getOrientation().equals("-"))
								end = new Alignments(annotation.getChr(), annotation.getStart() - extend, annotation.getStart() - (extend-window));
							else
								end = new Alignments(annotation.getChr(), annotation.getEnd() + (extend - window), annotation.getEnd() + extend);
							/*
							 * Get an interval tree for all/any exons that overlap with the extended region
							 */
							IntervalTree<RefSeqGeneWithIsoforms> endOverlappersTree = annotationCollection.getOverlappers(end);
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
									if(!(overlapper.getOrientedEnd() == annotation.getOrientedEnd())){
										overlapperIsSameGene = false;
									}	
								}
								if(!overlapperIsSameGene)
									break;
								// Because the extended region cannot overlap another annotation
							}
							//No overlap so continue with scoring the region
							tmpScores = libDataModel.scoreGene(new RefSeqGene(end));
							if(tmpScores[USE_SCORE]>bestScores[USE_SCORE]){
								for(int i=0;i<bestScores.length;i++){
									bestScores[i] = tmpScores[i];
								}
								bestScoreDistanceFromEnd = extend;
								geneToWindowMap.put(duplicateNameMap.get(annotation), new RefSeqGene(end));
							}
							if(tmpScores[USE_SCORE]>downstreamExpr)
								downstreamExpr = tmpScores[USE_SCORE];

							extend += (STEP);
						}
						
						/*
						 * Write the result for this annotation to the matrix
						 */
						resultMatrix.set(duplicateNameMap.get(annotation), BEST_EXPR_COL, bestScores[USE_SCORE]);
						resultMatrix.set(duplicateNameMap.get(annotation), BEST_PVAL_COL, bestScores[PVAL_SCORE]);
						resultMatrix.set(duplicateNameMap.get(annotation), ANNOTATED_END_EXPR_COL, annotatedEndExpr);
						if(maxIntoGene > 0)
							resultMatrix.set(duplicateNameMap.get(annotation), UPSTREAM_EXPR_COL, upstreamExpr);
						if(maxExtension > 0)
							resultMatrix.set(duplicateNameMap.get(annotation), DOWNSTREAM_EXPR_COL, downstreamExpr);
						if(maxIntoGene > 0 || maxExtension > 0)
							resultMatrix.set(duplicateNameMap.get(annotation), BEST_DIST_TO_END_COL, bestScoreDistanceFromEnd);
						
						//System.err.println("Going to score " + annotationEnd.toUCSC());
						annotation.setBedScore(bestScores[USE_SCORE]); //set RPKM as score
						annotation.setExtraFields(bestScores);
						annotation.addExtraField(duplicateNameMap.get(annotation));
						
						//bw.write(annotation.getAllIsoforms().iterator().next().toBED());
						//bw.newLine();
					}
				}
			}
			/**
			 * If there is no data for chromosome in annotation file, 
			 * put zeroes in all fields
			 */
			else{
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
						double [] bestScores = new double[5];
						RefSeqGene annotation = isoform_iter.next();
						resultMatrix.set(duplicateNameMap.get(annotation), BEST_EXPR_COL, 0.0);
						resultMatrix.set(duplicateNameMap.get(annotation), BEST_PVAL_COL, 1.0);
						resultMatrix.set(duplicateNameMap.get(annotation), ANNOTATED_END_EXPR_COL, 0.0);
						if(maxIntoGene > 0)
							resultMatrix.set(duplicateNameMap.get(annotation), UPSTREAM_EXPR_COL, 0.0);
						if(maxExtension > 0)
							resultMatrix.set(duplicateNameMap.get(annotation), DOWNSTREAM_EXPR_COL, 0.0);
						if(maxIntoGene > 0 || maxExtension > 0) 
							resultMatrix.set(duplicateNameMap.get(annotation), BEST_DIST_TO_END_COL, 0.0);
						
						annotation.setBedScore(bestScores[USE_SCORE]); //set RPKM as score
						annotation.setExtraFields(bestScores);
						annotation.addExtraField(duplicateNameMap.get(annotation));
					}
				}
			}
		} 
		}
		/*
		 * STRANDED DATA: SCORES STRANDS SEPARATELY
		 */
		else{
			String sizes = null;
			AlignmentDataModel alignmentsP=new GenericAlignmentDataModel(alignmentFile, sizes, true, minimumMappingQuality,removePCRDuplicatesFlag,weighReadsFlag);
			AlignmentDataModelStats alignmentDataP = new AlignmentDataModelStats(alignmentsP);
			ContinuousDataAlignmentModel libDataModelP = new ContinuousDataAlignmentModel(alignmentDataP);
			alignmentsP.setPositiveStranded();
			AlignmentDataModel alignmentsN=new GenericAlignmentDataModel(alignmentFile, sizes, true, minimumMappingQuality,removePCRDuplicatesFlag,weighReadsFlag);
			AlignmentDataModelStats alignmentDataN = new AlignmentDataModelStats(alignmentsN);
			ContinuousDataAlignmentModel libDataModelN = new ContinuousDataAlignmentModel(alignmentDataN);
			alignmentsN.setNegativeStranded();
			
			logger.info("Scoring stranded reads: ");
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
							RefSeqGene annotationEnd = annotation;					
							int annotationLength = annotation.getTranscriptLength();
							
							/*
							 * If the length of the annotated transcript > window size being analyzed,
							 * get sub annotation for window length
							 */
							if(annotationLength>window){
								annotationEnd = getSubAnnotationFromEnd(annotation,window,0);
								//TODO: Move this check outside this 
								if(annotationEnd == null){
									logger.warn("Annotation end for " + (annotationLength - window) + "-" + annotationLength + " --> " + annotation.toBED() + " was null.");
								}
							}
							//System.out.println(annotationEnd.getSequence());
							/*
							 * Calculate 
							 * [0] p-value
							 * [1] enrichment score
							 * [2] # aligned reads
							 * [3] average coverage
							 * [4] RPKM value
							 * [5] local lambda
							 * [6] gene length
							 * [7] poisson(#aligned reads)
							 * for the annotated region 	
							 * Finally calls scanPRate(RefSeqGene gene, int extensionFactor) in AlignmentDataModelStats					 
							 */
							double[] bestScores = scoreStrandedGene(libDataModelP,libDataModelN,annotationEnd);
							/*
							 * Use COUNTS as the score
							 * The best score is the window score for now
							 */
							double annotatedEndExpr = bestScores[USE_SCORE];
							int bestScoreDistanceFromEnd = 0;
							geneToWindowMap.put(duplicateNameMap.get(annotation), annotationEnd);

							double upstreamExpr = 0;
							double downstreamExpr = 0;

							double [] tmpScores = null;
							int retreat = STEP;
							/*
							 * If upstream extension is allowed, use sliding windows with overlaps of STEP
							 * while retreat region length is smaller than (annotation - window), i.e. it lies within the annotation
							 * and it is at distance less than the max region allowed upstream of end of gene
							 */
							while((retreat<(annotationLength - window)) && (retreat<maxIntoGene)){
								/*
								 * get annotation for region of length window, "retreat" length from end of transcript
								 */
								annotationEnd = getSubAnnotationFromEnd(annotation, window, retreat);
								/*
								 * calculate score for the window and find max score of windows
								 * Uses scanPRate(RefSeqGene, tree) which adds extension factor in AlignmentDataModelStats.
								 * Difference between scoregene and scanPRate from above is that this does not include exons.
								 */
								if(annotationEnd !=null){
									tmpScores = scoreStrandedGene(libDataModelP,libDataModelN,(annotationEnd));
									
									if(tmpScores[USE_SCORE]>bestScores[USE_SCORE]){
										for(int i=0;i<bestScores.length;i++){
											bestScores[i] = tmpScores[i];
										}
										bestScoreDistanceFromEnd = -retreat;
										geneToWindowMap.put(duplicateNameMap.get(annotation), annotationEnd);
									}
									if(tmpScores[USE_SCORE]>upstreamExpr)
										upstreamExpr = tmpScores[USE_SCORE];
								}
								retreat += STEP;
							}

							int extend = STEP;
							tmpScores = null;
							while(extend < maxExtension) {
								Alignments end = null;
								if(annotation.getOrientation().equals("-"))
									end = new Alignments(annotation.getChr(), annotation.getStart() - extend, annotation.getStart() - (extend-window));
								else
									end = new Alignments(annotation.getChr(), annotation.getEnd() + (extend - window), annotation.getEnd() + extend);
								/*
								 * Get an interval tree for all/any exons that overlap with the extended region
								 */
								IntervalTree<RefSeqGeneWithIsoforms> endOverlappersTree = annotationCollection.getOverlappers(end);
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
										if(!(overlapper.getOrientedEnd() == annotation.getOrientedEnd())){
											overlapperIsSameGene = false;
										}	
									}
									if(!overlapperIsSameGene)
										break;
									// Because the extended region cannot overlap another annotation
								}
								RefSeqGene geneEnd = new RefSeqGene(end);
								geneEnd.setOrientation(annotation.getOrientation());
								//No overlap so continue with scoring the region
								tmpScores = scoreStrandedGene(libDataModelP,libDataModelN,geneEnd);

								
								if(tmpScores[USE_SCORE]>bestScores[USE_SCORE]){
									for(int i=0;i<bestScores.length;i++){
										bestScores[i] = tmpScores[i];
									}
									bestScoreDistanceFromEnd = extend;
									geneToWindowMap.put(duplicateNameMap.get(annotation), geneEnd);
								}
								if(tmpScores[USE_SCORE]>downstreamExpr)
									downstreamExpr = tmpScores[USE_SCORE];

								extend += (STEP);
							}
							
							/*
							 * Write the result for this annotation to the matrix
							 */
							resultMatrix.set(duplicateNameMap.get(annotation), BEST_EXPR_COL, bestScores[USE_SCORE]);
							resultMatrix.set(duplicateNameMap.get(annotation), BEST_PVAL_COL, bestScores[PVAL_SCORE]);
							resultMatrix.set(duplicateNameMap.get(annotation), ANNOTATED_END_EXPR_COL, annotatedEndExpr);
							if(maxIntoGene > 0)
								resultMatrix.set(duplicateNameMap.get(annotation), UPSTREAM_EXPR_COL, upstreamExpr);
							if(maxExtension > 0)
								resultMatrix.set(duplicateNameMap.get(annotation), DOWNSTREAM_EXPR_COL, downstreamExpr);
							if(maxIntoGene > 0 || maxExtension > 0)
								resultMatrix.set(duplicateNameMap.get(annotation), BEST_DIST_TO_END_COL, bestScoreDistanceFromEnd);
							
							//System.err.println("Going to score " + annotationEnd.toUCSC());
							annotation.setBedScore(bestScores[USE_SCORE]); //set RPKM as score
							annotation.setExtraFields(bestScores);
							annotation.addExtraField(duplicateNameMap.get(annotation));
							
							//bw.write(annotation.getAllIsoforms().iterator().next().toBED());
							//bw.newLine();
						}
					}
				}
				/**
				 * If there is no data for chromosome in annotation file, 
				 * put zeroes in all fields
				 */
				else{
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
							double [] bestScores = new double[5];
							RefSeqGene annotation = isoform_iter.next();
							resultMatrix.set(duplicateNameMap.get(annotation), BEST_EXPR_COL, 0.0);
							resultMatrix.set(duplicateNameMap.get(annotation), BEST_PVAL_COL, 1.0);
							resultMatrix.set(duplicateNameMap.get(annotation), ANNOTATED_END_EXPR_COL, 0.0);
							if(maxIntoGene > 0)
								resultMatrix.set(duplicateNameMap.get(annotation), UPSTREAM_EXPR_COL, 0.0);
							if(maxExtension > 0)
								resultMatrix.set(duplicateNameMap.get(annotation), DOWNSTREAM_EXPR_COL, 0.0);
							if(maxIntoGene > 0 || maxExtension > 0) 
								resultMatrix.set(duplicateNameMap.get(annotation), BEST_DIST_TO_END_COL, 0.0);
							
							annotation.setBedScore(bestScores[USE_SCORE]); //set RPKM as score
							annotation.setExtraFields(bestScores);
							annotation.addExtraField(duplicateNameMap.get(annotation));
						}
					}
				}
			} 
		
		}
		/*
		 * Write to output File
		 */
		writeOutputFiles(outputFile, resultMatrix, allGenes, duplicateNameMap);
		if(collapseIsoforms){
			writeCollapsedOutputForSingle(annotationCollection,outputFile,resultMatrix,geneToWindowMap);
		}
		//BufferedWriter bw =argMap.getOutputWriter(); 
		//resultMatrix.write(bw);
		//bw.close();	
	}
	
	private static void score3PMultiple(BEDFileParser annotationCollection,Map<String,String> alignmentFiles,String outputFile,boolean normalizedOutput,Map<String,Collection<String>> conditionMaps) throws IOException, MathException{
		
		/*
		 * Restrict the annotations to one per gene. No duplicate annotations.
		 */
		List<RefSeqGene> allGenes = annotationCollection.GetGenes();
		List<BED> windows = new ArrayList<BED>();
		List<String> rows = new ArrayList<String>();
		
		/*
		 * Map of each gene to its peak.
		 */
		Map<String,RefSeqGene> geneToWindowMap = new HashMap<String,RefSeqGene>();

		/*
		 * HashMap of gene name to the number of duplicates
		 */
		HashMap<String, Integer> duplicateMap = new HashMap<String, Integer>();
		
		int duplicates=0;
		for(RefSeqGene gene : allGenes) {
			if(!rows.contains(gene.getName())) {
				String name = gene.getName();
				if(duplicateNameMap.containsKey(gene)){
					logger.info("Entry for "+name+" already exists");
				}
				else{
					rows.add(name);
					duplicateMap.put(name, 1);
					duplicateNameMap.put(gene, name);
				}
			} 
			// If the gene name has another annotation
			else {
				
				if(duplicateNameMap.containsKey(gene)){
					logger.info("Entry for "+gene.getName()+" already exists in "+duplicateNameMap.get(gene));
				}
				else{
					//Row name is now the geneName appended with the duplicate number
					duplicateMap.put(gene.getName(), (duplicateMap.get(gene.getName())+1));
					String name = (gene.getName()+"_"+duplicateMap.get(gene.getName()));
					rows.add(name);
					duplicateNameMap.put(gene, name);
					//logger.warn("Duplicated annotation : "+ gene.toBED());
					duplicates++;
				}
			}
		}
		logger.info("Found " + duplicates + " duplicates, going to process " + rows.size() + " annotations");
		
		/*
		 * The columns of the matrix will be each alignment file name
		 */
		List<String> cols = new ArrayList<String>();
		for(String name: conditionMaps.keySet()){
			for(String sample: conditionMaps.get(name)){
				cols.add(sample);
			}
		}
		
		/* 
		 * Initialize the matrix with row and column names
		 */
		MatrixWithHeaders resultMatrix = new MatrixWithHeaders(rows,cols);
		MatrixWithHeaders pvalueMatrix = new MatrixWithHeaders(rows,cols);
		MatrixWithHeaders fullGeneResultMatrix = null;
		if(fullGeneScoreFlag){
			/* 
			 * Initialize the matrix to store score of full gene with row and column names
			 */
			fullGeneResultMatrix = new MatrixWithHeaders(rows,cols);
		}
		//ELSE DONT INITIALIZE
		
		
		//ContinuousDataAlignmentModel libDataModel = null;
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
			libDataModels.put(ss, AlignmentUtils.loadAlignmentData(alignmentFiles.get(ss),false,minimumMappingQuality,removePCRDuplicatesFlag,weighReadsFlag));
		}		
		/*
		 * For each chromosome
		 */
		while(chrom_iter.hasNext()) {
			String chr = chrom_iter.next();
			
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

						RefSeqGene annotationEnd = annotation;					
						int annotationLength = annotation.getTranscriptLength();
						
						//If the length of the annotated transcript > window size being analyzed,
						//get sub annotation for window length
						
						if(annotationLength>window){
							//logger.info(annotation.getName()+" window: "+window);
							annotationEnd = getSubAnnotationFromEnd(annotation,window,0);
							if(annotationEnd == null){
								logger.warn("Annotation end for " + (annotationLength - window) + "-" + annotationLength + " --> " + annotation.toBED() + " was null.");
							}
						}
						//System.out.println(annotationEnd.getSequence());
						
						double bestEnrichment = getEnrichmentForWindow(libDataModels,annotationEnd);
						int bestWindowEnd = 0;
						boolean insideGene = true;
						RefSeqGene bestWindow = annotationEnd;	 
						
						// Use COUNTS as the score
						// The best score is the window score for now
						
//						int bestScoreDistanceFromEnd = 0;
						double tmpEnrichment = 0.0;
						int retreat = STEP;
						/*
						 * 							 
						 * If upstream extension is allowed, use sliding windows with overlaps of STEP
						 * while retreat region length is smaller than (annotation - window), i.e. it lies within the annotation
						 * and it is at distance less than the max region allowed upstream of end of gene
						 */
						while((retreat<(annotationLength - window)) && (retreat<maxIntoGene)){

							/*
							 * get annotation for region of length window, "retreat" length from end of transcript
							 */
							annotationEnd = getSubAnnotationFromEnd(annotation, window, retreat);
							/*
							 * calculate score for the window and find max score of windows
							 * Uses scanPRate(RefSeqGene, tree) which adds extension factor in AlignmentDataModelStats.
							 * Difference between scoregene and scanPRate from above is that this does not include exons.
							 */
							if(annotationEnd !=null){
								tmpEnrichment = getEnrichmentForWindow(libDataModels,annotationEnd);
								if(tmpEnrichment>bestEnrichment){
									bestEnrichment = tmpEnrichment;
									bestWindow = annotationEnd;
									bestWindowEnd = retreat;
									insideGene = true;
								}
							}
							retreat += STEP;
						}

						int extend = STEP;
						tmpEnrichment = 0.0;
						while(extend < maxExtension) {
							Alignments end = null;
							if(annotation.getOrientation().equals("-"))
								end = new Alignments(annotation.getChr(), annotation.getStart() - extend, annotation.getStart() - (extend-window));
							else
								end = new Alignments(annotation.getChr(), annotation.getEnd() + (extend - window), annotation.getEnd() + extend);
							/*
							 * Get an interval tree for all/any exons that overlap with the extended region
							 */
							IntervalTree<RefSeqGeneWithIsoforms> endOverlappersTree = annotationCollection.getOverlappers(end);
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
									overlapperIsSameGene = (overlapper.getOrientedEnd() == annotation.getOrientedEnd());
								}
								if(!overlapperIsSameGene)
									break;
								// Because the extended region cannot overlap another annotation
							}
							//No overlap so continue with scoring the region
							annotationEnd = new RefSeqGene(end);
							tmpEnrichment = getEnrichmentForWindow(libDataModels,annotationEnd);
							if(tmpEnrichment>bestEnrichment){
								bestEnrichment = tmpEnrichment;
								bestWindow = annotationEnd;
								bestWindowEnd = extend;
								insideGene = false;
							}
							extend += (STEP);
						}
						/*
						 * Write the result for this annotation to the matrix
						 */
						Map<String,Double>[] scores = getScoresForWindow(libDataModels,bestWindow);
						
						String rowName = duplicateNameMap.get(annotation);
						if(fullGeneScoreFlag){
							fullGeneResultMatrix = getScoreFromWindowToGeneStart(libDataModels,bestWindow,annotation,bestWindowEnd,insideGene,fullGeneResultMatrix,rowName);
						}
						
						for(int i=0;i<cols.size();i++){
							//resultMatrix.set(duplicateNameMap.get(annotation), cols.get(i), scores[i]);
							resultMatrix.set(rowName, cols.get(i), scores[0].get(cols.get(i)));
							pvalueMatrix.set(rowName, cols.get(i), scores[1].get(cols.get(i)));
						}
						geneToWindowMap.put(duplicateNameMap.get(annotation), bestWindow);
						windows.add(new BED(duplicateNameMap.get(annotation),annotation.getChr(),bestWindow.getStart(),bestWindow.getEnd()));
						//bw.write(annotation.getAllIsoforms().iterator().next().toBED());
						//bw.newLine();
					}
				}
			}
			else{
				logger.info("No data for " + chr);
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
						for(int i=0;i<cols.size();i++){
							resultMatrix.set(duplicateNameMap.get(annotation), cols.get(i), 0.0);
							pvalueMatrix.set(duplicateNameMap.get(annotation), cols.get(i), 1.0);
							if(fullGeneScoreFlag){
								fullGeneResultMatrix.set(duplicateNameMap.get(annotation), cols.get(i), 0.0);
							}
						}
						
						windows.add(new BED(duplicateNameMap.get(annotation),annotation.getChr(),annotation.getStart(),annotation.getEnd()));
					}
				}
			}
		}
		}
		//STRANDED READS
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
				libDataModelsN.put(ss,new ContinuousDataAlignmentModel(alignmentDataN[k]));
				k++;
			}
			/*
			 * For each chromosome
			 */
			while(chrom_iter.hasNext()) {
				String chr = chrom_iter.next();
				
				boolean dataForChr = false;
				for(String ss: libDataModelsP.keySet()){
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
						RefSeqGeneWithIsoforms annotation_with_isoforms = annotation_iter.next();
						Iterator<RefSeqGene> isoform_iter = annotation_with_isoforms.getAllIsoforms().iterator();

						while(isoform_iter.hasNext()){

							RefSeqGene annotation = isoform_iter.next();

							RefSeqGene annotationEnd = annotation;					
							int annotationLength = annotation.getTranscriptLength();
							
							//If the length of the annotated transcript > window size being analyzed,
							//get sub annotation for window length
							
							if(annotationLength>window){
								annotationEnd = getSubAnnotationFromEnd(annotation,window,0);
								if(annotationEnd == null){
									logger.warn("Annotation end for " + (annotationLength - window) + "-" + annotationLength + " --> " + annotation.toBED() + " was null.");
									break;
								}
							}
							//System.out.println(annotationEnd.getSequence());
							
							double bestEnrichment = getEnrichmentForWindow(libDataModelsP,libDataModelsN,annotationEnd);
							int bestWindowEnd = 0;
							boolean insideGene = true;
							RefSeqGene bestWindow = annotationEnd;	 
							
							// Use COUNTS as the score
							// The best score is the window score for now
							
//							int bestScoreDistanceFromEnd = 0;
							double tmpEnrichment = 0.0;
							int retreat = STEP;
							/*
							 * 							 
							 * If upstream extension is allowed, use sliding windows with overlaps of STEP
							 * while retreat region length is smaller than (annotation - window), i.e. it lies within the annotation
							 * and it is at distance less than the max region allowed upstream of end of gene
							 */
							while((retreat<(annotationLength - window)) && (retreat<maxIntoGene)){

								/*
								 * get annotation for region of length window, "retreat" length from end of transcript
								 */
								annotationEnd = getSubAnnotationFromEnd(annotation, window, retreat);
								/*
								 * calculate score for the window and find max score of windows
								 * Uses scanPRate(RefSeqGene, tree) which adds extension factor in AlignmentDataModelStats.
								 * Difference between scoregene and scanPRate from above is that this does not include exons.
								 */
								if(annotationEnd !=null){
									tmpEnrichment = getEnrichmentForWindow(libDataModelsP,libDataModelsN,annotationEnd);
									if(tmpEnrichment>bestEnrichment){
										bestEnrichment = tmpEnrichment;
										bestWindow = annotationEnd;
										bestWindowEnd = retreat;
										insideGene = true;
									}
								}
								retreat += STEP;
							}

							int extend = STEP;
							tmpEnrichment = 0.0;
							while(extend < maxExtension) {
								Alignments end = null;
								if(annotation.getOrientation().equals("-"))
									end = new Alignments(annotation.getChr(), annotation.getStart() - extend, annotation.getStart() - (extend-window));
								else
									end = new Alignments(annotation.getChr(), annotation.getEnd() + (extend - window), annotation.getEnd() + extend);
								/*
								 * Get an interval tree for all/any exons that overlap with the extended region
								 */
								IntervalTree<RefSeqGeneWithIsoforms> endOverlappersTree = annotationCollection.getOverlappers(end);
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
										overlapperIsSameGene = (overlapper.getOrientedEnd() == annotation.getOrientedEnd());
									}
									if(!overlapperIsSameGene)
										break;
									// Because the extended region cannot overlap another annotation
								}
								//No overlap so continue with scoring the region
								annotationEnd = new RefSeqGene(end);
								tmpEnrichment = getEnrichmentForWindow(libDataModelsP,libDataModelsN,annotationEnd);
								if(tmpEnrichment>bestEnrichment){
									bestEnrichment = tmpEnrichment;
									bestWindow = annotationEnd;
									bestWindowEnd = extend;
									insideGene = false;
								}
								extend += (STEP);
							}
							/*
							 * Write the result for this annotation to the matrix
							 */
							Map<String,Double>[] scores = getScoresForWindow(libDataModelsP,libDataModelsN,bestWindow);
							
							String rowName = duplicateNameMap.get(annotation);
							if(fullGeneScoreFlag){
								//TODO Check for unoriented
								if(annotation.isNegativeStrand())
									fullGeneResultMatrix = getScoreFromWindowToGeneStart(libDataModelsN,bestWindow,annotation,bestWindowEnd,insideGene,fullGeneResultMatrix,rowName);
								else
									fullGeneResultMatrix = getScoreFromWindowToGeneStart(libDataModelsP,bestWindow,annotation,bestWindowEnd,insideGene,fullGeneResultMatrix,rowName);
							}
							
							for(int i=0;i<cols.size();i++){
								//resultMatrix.set(duplicateNameMap.get(annotation), cols.get(i), scores[i]);
								resultMatrix.set(rowName, cols.get(i), scores[0].get(cols.get(i)));
								pvalueMatrix.set(rowName, cols.get(i), scores[1].get(cols.get(i)));
							}
							geneToWindowMap.put(duplicateNameMap.get(annotation), bestWindow);
							windows.add(new BED(duplicateNameMap.get(annotation),annotation.getChr(),bestWindow.getStart(),bestWindow.getEnd()));
							//bw.write(annotation.getAllIsoforms().iterator().next().toBED());
							//bw.newLine();
						}
					}
				}
				else{
					logger.info("No data for " + chr);
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
							for(int i=0;i<cols.size();i++){
								resultMatrix.set(duplicateNameMap.get(annotation), cols.get(i), 0.0);
								pvalueMatrix.set(duplicateNameMap.get(annotation), cols.get(i), 1.0);
								if(fullGeneScoreFlag){
									fullGeneResultMatrix.set(duplicateNameMap.get(annotation), cols.get(i), 0.0);
								}
							}
							
							windows.add(new BED(duplicateNameMap.get(annotation),annotation.getChr(),annotation.getStart(),annotation.getEnd()));
						}
					}
				}
			}
		}
		/*
		 * Write to output File
		 */
		writeOutputFiles(outputFile, resultMatrix, allGenes, duplicateNameMap,windows);
		writeOutputFiles(outputFile+".pvalues.matrix", pvalueMatrix, allGenes, duplicateNameMap,windows);
		if(fullGeneScoreFlag){
			BufferedWriter bw =new BufferedWriter(new FileWriter(outputFile+".fullgene.matrix")); 
			fullGeneResultMatrix.write(bw);
			bw.close();
		}
		
		if(collapseIsoforms){
			writeCollapsedOutputForMultiple(annotationCollection,outputFile,resultMatrix,geneToWindowMap);
		}
		
		if(normalizedOutput){
			
			writeNormalizedMatrix(outputFile, resultMatrix);
			if(fullGeneScoreFlag){
				writeNormalizedMatrix(outputFile+".fullgene.matrix", fullGeneResultMatrix);
			}
		}
		/*
		 * Call R-script to run DESeq
		 */
		//if differential expression has to be performed, conditions file is mandatory
	/*	String command = "RScript runDESeq.R "+argMap.getOutput()+" "+argMap.getMandatory("conditions");
	    Process p = Runtime.getRuntime().exec(command);
	    BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
	    String s;
        while ((s = stdInput.readLine()) != null) {
            logger.info(s);
        }*/
        
	}

	private static void score5P(BEDFileParser annotationCollection, String alignmentFile,String outputFile) throws IOException, ParseException, MathException {
		/*
		 * Restrict the annotations to one per gene. No duplicate annotations.
		 */
		List<RefSeqGene> allGenes = annotationCollection.GetGenes();
		/*
		 * Map of each gene to its peak.
		 */
		Map<String,RefSeqGene> geneToWindowMap = new HashMap<String,RefSeqGene>();
		
		/* each row of the matrix is one RefSeqGene name*/
		List<String> rows = new ArrayList<String>();
		/*
		 * HashMap of gene name to the number of duplicates
		 */
		HashMap<String, Integer> duplicateMap = new HashMap<String, Integer>();
		int duplicates=0;
		for(RefSeqGene gene : allGenes) {
			if(!rows.contains(gene.getName())) {
				String name = gene.getName();
				if(duplicateNameMap.containsKey(gene)){
					logger.info("Entry for "+name+" already exists");
				}
				else{
					rows.add(name);
					duplicateMap.put(name, 1);
					duplicateNameMap.put(gene, name);
				}
			} 
			// If the gene name has another annotation
			else {
				
				if(duplicateNameMap.containsKey(gene)){
					logger.info("Entry for "+gene.getName()+" already exists in "+duplicateNameMap.get(gene));
				}
				else{
					//Row name is now the geneName appended with the duplicate number
					duplicateMap.put(gene.getName(), (duplicateMap.get(gene.getName())+1));
					String name = (gene.getName()+"_"+duplicateMap.get(gene.getName()));
					rows.add(name);
					duplicateNameMap.put(gene, name);
					//logger.warn("Duplicated annotation : "+ gene.toBED());
					duplicates++;
				}
			}
		}
		logger.info("Found " + duplicates + " duplicates, ignoring them, going to process " + rows.size() + " annotations");
		

		/*
		 * Initialize the lists for output of expression score
		 * These lists names will also be the columns of the matrix
		 */
		List<String> cols = new ArrayList<String>();
		cols.add(BEST_EXPR_COL);
		cols.add(BEST_PVAL_COL);
		cols.add(ANNOTATED_END_EXPR_COL);
		if(maxIntoGene>0 || maxExtension>0){

			if(maxIntoGene>0)
				cols.add(DOWNSTREAM_EXPR_COL);
			if(maxExtension>0)
				cols.add(UPSTREAM_EXPR_COL);

			cols.add(BEST_DIST_TO_END_COL);
		}
		/* 
		 * Initialize the matrix with row and column names
		 */
		MatrixWithHeaders resultMatrix = new MatrixWithHeaders(rows,cols);

		/*
		 * Iterate through list of annotations. Iterator is over chromosome names.
		 */
		Iterator<String> chrom_iter = annotationCollection.getChromosomeIterator();

		if(!isStranded){

		/*
		 * Initialize the data model using the alignment file
		 * @param: <alignment flieName> <load_chromosome_stats> <minMappingQuality> <remove_PCR_duplicates> <weigh reads by NH flag>
		 */
		ContinuousDataAlignmentModel libDataModel = 	AlignmentUtils.loadAlignmentData(alignmentFile,true,minimumMappingQuality,removePCRDuplicatesFlag,weighReadsFlag);

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
						RefSeqGene annotationStart = annotation;					
						int annotationLength = annotation.getTranscriptLength();
						
						/*
						 * If the length of the annotated transcript > window size being analyzed,
						 * get sub annotation for window length
						 */
						if(annotationLength>window){
							annotationStart = getSubAnnotationFromStart(annotation,window,0);
							if(annotationStart == null){
								logger.warn("Annotation end for " + (annotationLength - window) + "-" + annotationLength + " --> " + annotation.toBED() + " was null.");
							}
						}
						//System.out.println(annotationStart.getSequence());
						/*
						 * Calculate 
						 * [0] p-value
						 * [1] enrichment score
						 * [2] # aligned reads
						 * [3] average coverage
						 * [4] RPKM value
						 * [5] local lambda
						 * [6] gene length
						 * [7] poisson(#aligned reads)
						 * for the annotated region 	
						 * Finally calls scanPRate(RefSeqGene gene, int extensionFactor) in AlignmentDataModelStats					 
						 */
						double [] bestScores = libDataModel.scanPRate(annotationStart);		
						geneToWindowMap.put(duplicateNameMap.get(annotation), annotationStart);
						/*
						 * Use COUNTS as the score
						 * The best score is the window score for now
						 */
						double annotatedStartExpr = bestScores[USE_SCORE];
						int bestScoreDistanceFromStart = 0;

						double upstreamExpr = 0;
						double downstreamExpr = 0;

						double [] tmpScores = null;
						int intoGene = STEP;
						/*
						 * If downstream extension is allowed, use sliding windows with overlaps of STEP
						 * while intoGene region length is smaller than (annotation - window), i.e. it lies within the annotation
						 * and it is at distance less than the max region allowed downstream from start of gene
						 */
						while((intoGene<(annotationLength - window)) && (intoGene<maxIntoGene)){

							/*
							 * get annotation for region of length window, "intoGene" length from start of transcript
							 */
							annotationStart = getSubAnnotationFromStart(annotation, window, intoGene);
							/*
							 * calculate score for the window and find max score of windows
							 * Uses scanPRate(RefSeqGene, tree) which adds extension factor in AlignmentDataModelStats.
							 * Difference between scoregene and scanPRate from above is that this does not include exons.
							 */
							if(annotationStart !=null){
								tmpScores = libDataModel.scoreGene(annotationStart);
								if(tmpScores[USE_SCORE]>bestScores[USE_SCORE]){
									for(int i=0;i<bestScores.length;i++){
										bestScores[i] = tmpScores[i];
									}
									bestScoreDistanceFromStart = intoGene;
									geneToWindowMap.put(duplicateNameMap.get(annotation), annotationStart);
								}
								if(tmpScores[USE_SCORE]>downstreamExpr)
									downstreamExpr = tmpScores[USE_SCORE];
							}
							intoGene += STEP;
						}

						int extend = STEP;
						tmpScores = null;
						while(extend < maxExtension) {
							Alignments start = null;
							if(annotation.getOrientation().equals("+"))
								start = new Alignments(annotation.getChr(), annotation.getStart() - extend, annotation.getStart() - (extend-window));
							else
								start = new Alignments(annotation.getChr(), annotation.getEnd() + (extend - window), annotation.getEnd() + extend);
							/*
							 * Get an interval tree for all/any exons that overlap with the extended region
							 */
							IntervalTree<RefSeqGeneWithIsoforms> startOverlappersTree = annotationCollection.getOverlappers(start);
							/*
							 * If there is an overlap with a gene
							 */
							if(!startOverlappersTree.isEmpty()){
								Iterator<RefSeqGeneWithIsoforms> overlappersIter = startOverlappersTree.valueIterator();
								boolean overlapperIsSameGene = true;
								/*
								 * while the gene is the same gene
								 */
								while(overlappersIter.hasNext() && overlapperIsSameGene){
									RefSeqGene overlapper = overlappersIter.next();
									//compare the end coordiantes of the gene
									if(!(overlapper.getOrientedStart() == annotation.getOrientedStart())){
										overlapperIsSameGene = false;
									}	
								}
								if(!overlapperIsSameGene)
									break;
								// Because the extended region cannot overlap another annotation
							}
							//No overlap so continue with scoring the region
							tmpScores = libDataModel.scoreGene(new RefSeqGene(start));
							if(tmpScores[USE_SCORE]>bestScores[USE_SCORE]){
								for(int i=0;i<bestScores.length;i++){
									bestScores[i] = tmpScores[i];
								}
								bestScoreDistanceFromStart = -extend;
								geneToWindowMap.put(duplicateNameMap.get(annotation), new RefSeqGene(start));
							}
							if(tmpScores[USE_SCORE]>upstreamExpr)
								upstreamExpr = tmpScores[USE_SCORE];

							extend += (STEP);
						}
						
						/*
						 * Write the result for this annotation to the matrix
						 */
						resultMatrix.set(duplicateNameMap.get(annotation), BEST_EXPR_COL, bestScores[USE_SCORE]);
						resultMatrix.set(duplicateNameMap.get(annotation), BEST_PVAL_COL, bestScores[PVAL_SCORE]);
						resultMatrix.set(duplicateNameMap.get(annotation), ANNOTATED_END_EXPR_COL, annotatedStartExpr);
						if(maxIntoGene > 0)
							resultMatrix.set(duplicateNameMap.get(annotation), DOWNSTREAM_EXPR_COL, downstreamExpr);
						if(maxExtension > 0)
							resultMatrix.set(duplicateNameMap.get(annotation), UPSTREAM_EXPR_COL, upstreamExpr);
						if(maxIntoGene > 0 || maxExtension > 0)
							resultMatrix.set(duplicateNameMap.get(annotation), BEST_DIST_TO_END_COL, bestScoreDistanceFromStart);
						
						//System.err.println("Going to score " + annotationStart.toUCSC());
						annotation.setBedScore(bestScores[USE_SCORE]); //set RPKM as score
						annotation.setExtraFields(bestScores);
						annotation.addExtraField(duplicateNameMap.get(annotation));
												
						//bw.write(annotation.getAllIsoforms().iterator().next().toBED());
						//bw.newLine();
					}
				}
			}
			/**
			 * If there is no data for chromosome in annotation file, 
			 * put zeroes in all fields
			 */
			else{
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
						double [] bestScores = new double[5];
						RefSeqGene annotation = isoform_iter.next();
						resultMatrix.set(duplicateNameMap.get(annotation), BEST_EXPR_COL, 0.0);
						resultMatrix.set(duplicateNameMap.get(annotation), BEST_PVAL_COL, 1.0);
						resultMatrix.set(duplicateNameMap.get(annotation), ANNOTATED_END_EXPR_COL, 0.0);
						if(maxIntoGene > 0)
							resultMatrix.set(duplicateNameMap.get(annotation), DOWNSTREAM_EXPR_COL, 0.0);
						if(maxExtension > 0)
							resultMatrix.set(duplicateNameMap.get(annotation), UPSTREAM_EXPR_COL, 0.0);
						if(maxIntoGene > 0 || maxExtension > 0) 
							resultMatrix.set(duplicateNameMap.get(annotation), BEST_DIST_TO_END_COL, 0.0);
						
						annotation.setBedScore(bestScores[USE_SCORE]); //set RPKM as score
						annotation.setExtraFields(bestScores);
						annotation.addExtraField(duplicateNameMap.get(annotation));
					}
				}
			}
		}
		}
		/*
		 * STRANDED DATA: SCORES STRANDS SEPARATELY
		 */
		else{
			String sizes = null;
			AlignmentDataModel alignmentsP=new GenericAlignmentDataModel(alignmentFile, sizes, true, minimumMappingQuality,removePCRDuplicatesFlag,weighReadsFlag);
			AlignmentDataModelStats alignmentDataP = new AlignmentDataModelStats(alignmentsP);
			ContinuousDataAlignmentModel libDataModelP = new ContinuousDataAlignmentModel(alignmentDataP);
			alignmentsP.setPositiveStranded();
			AlignmentDataModel alignmentsN=new GenericAlignmentDataModel(alignmentFile, sizes, true, minimumMappingQuality,removePCRDuplicatesFlag,weighReadsFlag);
			AlignmentDataModelStats alignmentDataN = new AlignmentDataModelStats(alignmentsN);
			ContinuousDataAlignmentModel libDataModelN = new ContinuousDataAlignmentModel(alignmentDataN);
			alignmentsN.setNegativeStranded();
			
			logger.info("Scoring stranded reads: ");
			
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
							RefSeqGene annotationStart = annotation;					
							int annotationLength = annotation.getTranscriptLength();
							
							/*
							 * If the length of the annotated transcript > window size being analyzed,
							 * get sub annotation for window length
							 */
							if(annotationLength>window){
								annotationStart = getSubAnnotationFromStart(annotation,window,0);
								if(annotationStart == null){
									logger.warn("Annotation end for " + (annotationLength - window) + "-" + annotationLength + " --> " + annotation.toBED() + " was null.");
								}
							}
							//System.out.println(annotationStart.getSequence());
							/*
							 * Calculate 
							 * [0] p-value
							 * [1] enrichment score
							 * [2] # aligned reads
							 * [3] average coverage
							 * [4] RPKM value
							 * [5] local lambda
							 * [6] gene length
							 * [7] poisson(#aligned reads)
							 * for the annotated region 	
							 * Finally calls scanPRate(RefSeqGene gene, int extensionFactor) in AlignmentDataModelStats					 
							 */
							double[] bestScores = scoreStrandedGene(libDataModelP,libDataModelN,annotationStart);	
							/*
							 * Use COUNTS as the score
							 * The best score is the window score for now
							 */
							double annotatedStartExpr = bestScores[USE_SCORE];
							geneToWindowMap.put(duplicateNameMap.get(annotation), annotationStart);
							int bestScoreDistanceFromStart = 0;

							double upstreamExpr = 0;
							double downstreamExpr = 0;

							double [] tmpScores = null;
							int intoGene = STEP;
							/*
							 * If downstream extension is allowed, use sliding windows with overlaps of STEP
							 * while intoGene region length is smaller than (annotation - window), i.e. it lies within the annotation
							 * and it is at distance less than the max region allowed downstream from start of gene
							 */
							while((intoGene<(annotationLength - window)) && (intoGene<maxIntoGene)){
								/*
								 * get annotation for region of length window, "intoGene" length from start of transcript
								 */
								annotationStart = getSubAnnotationFromStart(annotation, window, intoGene);
								/*
								 * calculate score for the window and find max score of windows
								 * Uses scanPRate(RefSeqGene, tree) which adds extension factor in AlignmentDataModelStats.
								 * Difference between scoregene and scanPRate from above is that this does not include exons.
								 */
								if(annotationStart !=null){
									tmpScores = scoreStrandedGene(libDataModelP,libDataModelN,annotationStart);
									if(tmpScores[USE_SCORE]>bestScores[USE_SCORE]){
										for(int i=0;i<bestScores.length;i++){
											bestScores[i] = tmpScores[i];
										}
										bestScoreDistanceFromStart = intoGene;
										geneToWindowMap.put(duplicateNameMap.get(annotation), annotationStart);
									}
									if(tmpScores[USE_SCORE]>downstreamExpr)
										downstreamExpr = tmpScores[USE_SCORE];
								}
								intoGene += STEP;
							}

							int extend = STEP;
							tmpScores = null;
							while(extend < maxExtension) {
								Alignments start = null;
								if(annotation.getOrientation().equals("+"))
									start = new Alignments(annotation.getChr(), annotation.getStart() - extend, annotation.getStart() - (extend-window));
								else
									start = new Alignments(annotation.getChr(), annotation.getEnd() + (extend - window), annotation.getEnd() + extend);
								/*
								 * Get an interval tree for all/any exons that overlap with the extended region
								 */
								IntervalTree<RefSeqGeneWithIsoforms> startOverlappersTree = annotationCollection.getOverlappers(start);
								/*
								 * If there is an overlap with a gene
								 */
								if(!startOverlappersTree.isEmpty()){
									Iterator<RefSeqGeneWithIsoforms> overlappersIter = startOverlappersTree.valueIterator();
									boolean overlapperIsSameGene = true;
									/*
									 * while the gene is the same gene
									 */
									while(overlappersIter.hasNext() && overlapperIsSameGene){
										RefSeqGene overlapper = overlappersIter.next();
										//compare the end coordiantes of the gene
										if(!(overlapper.getOrientedStart() == annotation.getOrientedStart())){
											overlapperIsSameGene = false;
										}	
									}
									if(!overlapperIsSameGene)
										break;
									// Because the extended region cannot overlap another annotation
								}
								RefSeqGene geneStart = new RefSeqGene(start);
								geneStart.setOrientation(annotation.getOrientation());
								//No overlap so continue with scoring the region
								tmpScores = scoreStrandedGene(libDataModelP,libDataModelN,geneStart);
								if(tmpScores[USE_SCORE]>bestScores[USE_SCORE]){
									for(int i=0;i<bestScores.length;i++){
										bestScores[i] = tmpScores[i];
									}
									bestScoreDistanceFromStart = -extend;
									geneToWindowMap.put(duplicateNameMap.get(annotation), geneStart);
								}
								if(tmpScores[USE_SCORE]>upstreamExpr)
									upstreamExpr = tmpScores[USE_SCORE];

								extend += (STEP);
							}
							
							/*
							 * Write the result for this annotation to the matrix
							 */
							resultMatrix.set(duplicateNameMap.get(annotation), BEST_EXPR_COL, bestScores[USE_SCORE]);
							resultMatrix.set(duplicateNameMap.get(annotation), BEST_PVAL_COL, bestScores[PVAL_SCORE]);
							resultMatrix.set(duplicateNameMap.get(annotation), ANNOTATED_END_EXPR_COL, annotatedStartExpr);
							if(maxIntoGene > 0)
								resultMatrix.set(duplicateNameMap.get(annotation), DOWNSTREAM_EXPR_COL, downstreamExpr);
							if(maxExtension > 0)
								resultMatrix.set(duplicateNameMap.get(annotation), UPSTREAM_EXPR_COL, upstreamExpr);
							if(maxIntoGene > 0 || maxExtension > 0)
								resultMatrix.set(duplicateNameMap.get(annotation), BEST_DIST_TO_END_COL, bestScoreDistanceFromStart);
							
							//System.err.println("Going to score " + annotationStart.toUCSC());
							annotation.setBedScore(bestScores[USE_SCORE]); //set RPKM as score
							annotation.setExtraFields(bestScores);
							annotation.addExtraField(duplicateNameMap.get(annotation));

							//bw.write(annotation.getAllIsoforms().iterator().next().toBED());
							//bw.newLine();
						}
					}
				}
				/**
				 * If there is no data for chromosome in annotation file, 
				 * put zeroes in all fields
				 */
				else{
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
							double [] bestScores = new double[5];
							RefSeqGene annotation = isoform_iter.next();
							resultMatrix.set(duplicateNameMap.get(annotation), BEST_EXPR_COL, 0.0);
							resultMatrix.set(duplicateNameMap.get(annotation), BEST_PVAL_COL, 1.0);
							resultMatrix.set(duplicateNameMap.get(annotation), ANNOTATED_END_EXPR_COL, 0.0);
							if(maxIntoGene > 0)
								resultMatrix.set(duplicateNameMap.get(annotation), DOWNSTREAM_EXPR_COL, 0.0);
							if(maxExtension > 0)
								resultMatrix.set(duplicateNameMap.get(annotation), UPSTREAM_EXPR_COL, 0.0);
							if(maxIntoGene > 0 || maxExtension > 0) 
								resultMatrix.set(duplicateNameMap.get(annotation), BEST_DIST_TO_END_COL, 0.0);
							
							annotation.setBedScore(bestScores[USE_SCORE]); //set RPKM as score
							annotation.setExtraFields(bestScores);
							annotation.addExtraField(duplicateNameMap.get(annotation));
						}
					}
				}
			}
		}
		/*
		 * Write to output File
		 */
		writeOutputFiles(outputFile, resultMatrix, allGenes, duplicateNameMap);
		if(collapseIsoforms){
			writeCollapsedOutputForSingle(annotationCollection,outputFile,resultMatrix,geneToWindowMap);
		}
		//BufferedWriter bw =argMap.getOutputWriter(); 
		//resultMatrix.write(bw);
		//bw.close();	
	}
	
	
	private static void score5PMultiple(BEDFileParser annotationCollection,Map<String,String> alignmentFiles,String outputFile,boolean normalizedOutput,Map<String,Collection<String>> conditionMaps) throws IOException, MathException{
		

		/*
		 * Restrict the annotations to one per gene. No duplicate annotations.
		 */
		List<RefSeqGene> allGenes = annotationCollection.GetGenes();
		List<BED> windows = new ArrayList<BED>();
		List<String> rows = new ArrayList<String>();
		/*
		 * Map of each gene to its peak.
		 */
		Map<String,RefSeqGene> geneToWindowMap = new HashMap<String,RefSeqGene>();

		/*
		 * HashMap of gene name to the number of duplicates
		 */
		HashMap<String, Integer> duplicateMap = new HashMap<String, Integer>();
		
		int duplicates=0;
		for(RefSeqGene gene : allGenes) {
			if(!rows.contains(gene.getName())) {
				String name = gene.getName();
				if(duplicateNameMap.containsKey(gene)){
					logger.info("Entry for "+name+" already exists");
				}
				else{
					rows.add(name);
					duplicateMap.put(name, 1);
					duplicateNameMap.put(gene, name);
				}
			} 
			// If the gene name has another annotation
			else {
				
				if(duplicateNameMap.containsKey(gene)){
					logger.info("Entry for "+gene.getName()+" already exists in "+duplicateNameMap.get(gene));
				}
				else{
					//Row name is now the geneName appended with the duplicate number
					duplicateMap.put(gene.getName(), (duplicateMap.get(gene.getName())+1));
					String name = (gene.getName()+"_"+duplicateMap.get(gene.getName()));
					rows.add(name);
					duplicateNameMap.put(gene, name);
					//logger.warn("Duplicated annotation : "+ gene.toBED());
					duplicates++;
				}
			}
		}
		logger.info("Found " + duplicates + " duplicates, going to process " + rows.size() + " annotations");
		
		/*
		 * The columns of the matrix will be each alignment file name
		 */
		List<String> cols = new ArrayList<String>();
		for(String name: conditionMaps.keySet()){
			for(String sample: conditionMaps.get(name))
				cols.add(sample);
		}
		/* 
		 * Initialize the matrix with row and column names
		 */
		MatrixWithHeaders resultMatrix = new MatrixWithHeaders(rows,cols);
		MatrixWithHeaders pvalueMatrix = new MatrixWithHeaders(rows,cols);
		MatrixWithHeaders fullGeneResultMatrix = null;
		if(fullGeneScoreFlag){
			/* 
			 * Initialize the matrix to store score of full gene with row and column names
			 */
			fullGeneResultMatrix = new MatrixWithHeaders(rows,cols);
		}
		//ELSE DONT INITIALIZE
		
		//ContinuousDataAlignmentModel libDataModel = null;
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
			libDataModels.put(ss, AlignmentUtils.loadAlignmentData(alignmentFiles.get(ss),false,minimumMappingQuality,removePCRDuplicatesFlag,weighReadsFlag));
		}	
		/*
		 * For each chromosome
		 */
		while(chrom_iter.hasNext()) {
			String chr = chrom_iter.next();
			
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
						int bestWindowStart = 0;
						boolean insideGene = true;
						RefSeqGene bestWindow = annotationStart;	 
						
						// Use COUNTS as the score
						// The best score is the window score for now
						
//						int bestScoreDistanceFromEnd = 0;
						double tmpEnrichment = 0.0;
						int intoGene = STEP;
						/*
						 * 							 
						 * If downstream extension is allowed, use sliding windows with overlaps of STEP
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
									bestWindowStart = intoGene;
									insideGene = true;
								}
							}
							intoGene += STEP;
						}

						int extend = STEP;
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
								bestWindowStart = extend;
								insideGene = false;
							}
							extend += (STEP);
						}
						/*
						 * Write the result for this annotation to the matrix
						 */
						Map<String,Double>[] scores = getScoresForWindow(libDataModels,bestWindow);
						
						String rowName = duplicateNameMap.get(annotation);
						if(fullGeneScoreFlag){
							fullGeneResultMatrix = getScoreFromWindowToGeneEnd(libDataModels,bestWindow,annotation,bestWindowStart,insideGene,fullGeneResultMatrix,rowName);
						}
						
						for(int i=0;i<cols.size();i++){
							//resultMatrix.set(duplicateNameMap.get(annotation), cols.get(i), scores[i]);
							resultMatrix.set(rowName, cols.get(i), scores[0].get(cols.get(i)));
							pvalueMatrix.set(rowName, cols.get(i), scores[1].get(cols.get(i)));
						}
						geneToWindowMap.put(duplicateNameMap.get(annotation), bestWindow);
						windows.add(new BED(duplicateNameMap.get(annotation),annotation.getChr(),bestWindow.getStart(),bestWindow.getEnd()));
						//bw.write(annotation.getAllIsoforms().iterator().next().toBED());
						//bw.newLine();
					}
				}
			}
			else{
				System.err.println("No data for " + chr);
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
						for(int i=0;i<cols.size();i++){
							resultMatrix.set(duplicateNameMap.get(annotation), cols.get(i), 0.0);
							pvalueMatrix.set(duplicateNameMap.get(annotation), cols.get(i), 1.0);
							if(fullGeneScoreFlag){
								fullGeneResultMatrix.set(duplicateNameMap.get(annotation), cols.get(i), 0.0);
							}
						}
						
						windows.add(new BED(duplicateNameMap.get(annotation),annotation.getChr(),annotation.getStart(),annotation.getEnd()));
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
				k++;
			}
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
						RefSeqGeneWithIsoforms annotation_with_isoforms = annotation_iter.next();
						Iterator<RefSeqGene> isoform_iter = annotation_with_isoforms.getAllIsoforms().iterator();

						while(isoform_iter.hasNext()){

							RefSeqGene annotation = isoform_iter.next();

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
							
							double bestEnrichment = getEnrichmentForWindow(libDataModelsP,libDataModelsN,annotationStart);
							int bestWindowStart = 0;
							boolean insideGene = true;
							RefSeqGene bestWindow = annotationStart;	 
							
							// Use COUNTS as the score
							// The best score is the window score for now
							
//							int bestScoreDistanceFromEnd = 0;
							double tmpEnrichment = 0.0;
							int intoGene = STEP;
							/*
							 * 							 
							 * If downstream extension is allowed, use sliding windows with overlaps of STEP
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
									tmpEnrichment = getEnrichmentForWindow(libDataModelsP,libDataModelsN,annotationStart);
									if(tmpEnrichment>bestEnrichment){
										bestEnrichment = tmpEnrichment;
										bestWindow = annotationStart;
										bestWindowStart = intoGene;
										insideGene = true;
									}
								}
								intoGene += STEP;
							}

							int extend = STEP;
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
								tmpEnrichment = getEnrichmentForWindow(libDataModelsP,libDataModelsN,annotationStart);
								if(tmpEnrichment>bestEnrichment){
									bestEnrichment = tmpEnrichment;
									bestWindow = annotationStart;
									bestWindowStart = extend;
									insideGene = false;
								}
								extend += (STEP);
							}
							/*
							 * Write the result for this annotation to the matrix
							 */
							Map<String,Double>[] scores = getScoresForWindow(libDataModelsP,libDataModelsN,bestWindow);
							
							String rowName = duplicateNameMap.get(annotation);
							if(fullGeneScoreFlag){
								//TODO Check for unoriented
								if(annotation.isNegativeStrand())
									fullGeneResultMatrix = getScoreFromWindowToGeneEnd(libDataModelsN,bestWindow,annotation,bestWindowStart,insideGene,fullGeneResultMatrix,rowName);
								else
									fullGeneResultMatrix = getScoreFromWindowToGeneEnd(libDataModelsP,bestWindow,annotation,bestWindowStart,insideGene,fullGeneResultMatrix,rowName);
							}

							
							for(int i=0;i<cols.size();i++){
								//resultMatrix.set(duplicateNameMap.get(annotation), cols.get(i), scores[i]);
								resultMatrix.set(rowName, cols.get(i), scores[0].get(cols.get(i)));
								pvalueMatrix.set(rowName, cols.get(i), scores[1].get(cols.get(i)));
							}
							geneToWindowMap.put(duplicateNameMap.get(annotation), bestWindow);
							windows.add(new BED(duplicateNameMap.get(annotation),annotation.getChr(),bestWindow.getStart(),bestWindow.getEnd()));
							//bw.write(annotation.getAllIsoforms().iterator().next().toBED());
							//bw.newLine();
						}
					}
				}
				else{
					System.err.println("No data for " + chr);
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
							for(int i=0;i<cols.size();i++){
								resultMatrix.set(duplicateNameMap.get(annotation), cols.get(i), 0.0);
								pvalueMatrix.set(duplicateNameMap.get(annotation), cols.get(i), 1.0);
								if(fullGeneScoreFlag){
									fullGeneResultMatrix.set(duplicateNameMap.get(annotation), cols.get(i), 0.0);
								}
							}
							
							windows.add(new BED(duplicateNameMap.get(annotation),annotation.getChr(),annotation.getStart(),annotation.getEnd()));
						}
					}
				}
			}
		}
		/*
		 * Write to output File
		 */
		writeOutputFiles(outputFile, resultMatrix, allGenes, duplicateNameMap,windows);
		writeOutputFiles(outputFile+".pvalues.matrix", pvalueMatrix, allGenes, duplicateNameMap,windows);
		//BufferedWriter bw =new BufferedWriter(new FileWriter(argMap.getOutput()+".fullgene.matrix")); 
		//fullGeneResultMatrix.write(bw);
		//bw.close();
		if(fullGeneScoreFlag){
			BufferedWriter bw =new BufferedWriter(new FileWriter(outputFile+".fullgene.matrix")); 
			fullGeneResultMatrix.write(bw);
			bw.close();
			if(normalizedOutput){
				writeNormalizedMatrix(outputFile+".fullgene.matrix", fullGeneResultMatrix);
			}
		}
		
		if(normalizedOutput){			
			writeNormalizedMatrix(outputFile, resultMatrix);
		}
		if(collapseIsoforms){
			writeCollapsedOutputForMultiple(annotationCollection,outputFile,resultMatrix,geneToWindowMap);
		}
		/*
		 * Call R-script to run DESeq
		 */
		//if differential expression has to be performed, conditions file is mandatory
	/*	String command = "RScript runDESeq.R "+argMap.getOutput()+" "+argMap.getMandatory("conditions");
	    Process p = Runtime.getRuntime().exec(command);
	    BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
	    String s;
        while ((s = stdInput.readLine()) != null) {
            logger.info(s);
        }*/
        
	}

	
	public static BEDFileParser reconstructGeneEnds(BEDFileParser completeReconstruction, BEDFileParser annotations,String annotationFile) throws IOException, MathException{
		
		/*
		 *  ARGUMENTS:
		 *  double minSpliceFrequency
		 *  int[] fixedWidth
		 *  String chrToSegment
		 *  boolean trimEnds	FALSE FOR NOW
		 *  boolean filterCanonical
		 *  String sequenceFile
		 *  double alpha
		 *  int chrStart
		 *  int chrEnd
		 *  double trimQuantile
		 */
		//BEDFileParser completeReconstruction = completeCDAM.segmentChromosome(0.1, fixedWidth, null, true, true, sequenceFile, 0.05, 0, Integer.MAX_VALUE, 0.1);
		completeReconstruction.makeGenes(); //This avoids having to compare to many different isoforms of different sizes usually present in reconstructed transcriptomes.
		completeReconstruction.writeFullBed("reconstructions.togenes.bed");
		annotations.equalizeTranscriptEnds(completeReconstruction);
		annotations.writeFullBed(annotationFile+".3p5pAdjusted.bed");
		
		return annotations;
	}
	
	/*
	 * This function will reconstruct the gene ends of all genes in the annotation file using a chi-square distance distrbution at the defined gene ends
	 * @return
	 * @throws IOException
	 * @throws MathException
	 
	private static BEDFileParser trimGeneEnds(ContinuousDataAlignmentModel completeCDAM,BEDFileParser annotations,String annotationFile) throws IOException, MathException{
			
			BEDFileParser reconstructedGenes = completeCDAM.reconstructGeneEndsUsingDivergence(annotations,true,0.1,annotationFile);
			annotations.writeFullBed(annotationFile+".3pAdjustedTrimmed.bed");
			return reconstructedGenes;
		}*/

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
				double score = model.scoreGene(annotationEnd)[1];
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
	 * This function calculates the best enrichment score for a given window in a strand specific way
	 * @param libDataModelsP
	 * @param libDataModelsM
	 * @param annotationEnd
	 * @return Vector of size libDataModels with enrichment of each 
	 * @throws IOException
	 */
	private static double getEnrichmentForWindow(Map<String,ContinuousDataAlignmentModel> libDataModelsP,Map<String,ContinuousDataAlignmentModel> libDataModelsN,RefSeqGene annotationEnd) throws IOException{
		
		//double max = 0.0;
		double sum = 0.0;
		if(!annotationEnd.isUnoriented()){
			if(annotationEnd.isNegativeStrand()){
				for(ContinuousDataAlignmentModel model:libDataModelsN.values()){
					if(model.hasDataForChromosome(annotationEnd.getChr())){	
						double score = model.scoreGene(annotationEnd)[1];
						/*if(score>max){
							max = score;
						}*/
						sum += score;
					}
				}
			}
			else{
				for(ContinuousDataAlignmentModel model:libDataModelsP.values()){
					if(model.hasDataForChromosome(annotationEnd.getChr())){	
						double score = model.scoreGene(annotationEnd)[1];
						/*if(score>max){
							max = score;
						}*/
						sum += score;
					}
				}
			}
		}
		//For unoriented genes, sum up across both positive and negative models
		else{
			for(ContinuousDataAlignmentModel model:libDataModelsN.values()){
				if(model.hasDataForChromosome(annotationEnd.getChr())){	
					double score = model.scoreGene(annotationEnd)[1];
					/*if(score>max){
						max = score;
					}*/
					sum += score;
				}
			}
			for(ContinuousDataAlignmentModel model:libDataModelsP.values()){
				if(model.hasDataForChromosome(annotationEnd.getChr())){	
					double score = model.scoreGene(annotationEnd)[1];
					/*if(score>max){
						max = score;
					}*/
					sum += score;
				}
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
	private static Map<String,Double>[] getScoresForWindow(Map<String,ContinuousDataAlignmentModel> libDataModels,RefSeqGene annotationEnd) throws IOException{
		
		//[0] is counts and [1] is pvalues
		Map<String,Double>[] scores = new HashMap[2];
		for(int i=0;i<2;i++){
			scores[i] = new HashMap<String,Double>();
		}
		for(String ss:libDataModels.keySet()){
			if(libDataModels.get(ss).hasDataForChromosome(annotationEnd.getChr())){	
				double[] s = libDataModels.get(ss).scoreGene(annotationEnd);
				scores[0].put(ss, s[USE_SCORE]);
				scores[1].put(ss, s[PVAL_SCORE]);
			}
			else{
				scores[0].put(ss, 0.0);
				scores[1].put(ss, 1.0);
			}
		}
		return scores;
	}
	
	/**
	 * This function calculates the scores for all alignment models for a given window in a strand-specific manner
	 * @param libDataModelsP
	 * @param libDataModelsN
	 * @param annotationEnd
	 * @return Vector of size libDataModels with enrichment of each 
	 * @throws IOException
	 */
	private static Map<String,Double>[] getScoresForWindow(Map<String,ContinuousDataAlignmentModel> libDataModelsP,Map<String,ContinuousDataAlignmentModel> libDataModelsN,RefSeqGene annotationEnd) throws IOException{
		
		Map<String,Double>[] scores = new HashMap[2];
		for(int i=0;i<2;i++){
			scores[i] = new HashMap<String,Double>();
		}
		if(!annotationEnd.isUnoriented()){
			if(annotationEnd.isNegativeStrand()){
				for(String st:libDataModelsN.keySet()){
					if(libDataModelsN.get(st).hasDataForChromosome(annotationEnd.getChr())){
						double[] s=libDataModelsN.get(st).scoreGene(annotationEnd);
						scores[0].put(st,s[USE_SCORE]);
						scores[1].put(st,s[PVAL_SCORE]);
					}
				}
			}
			else{
				for(String st:libDataModelsP.keySet()){
					if(libDataModelsP.get(st).hasDataForChromosome(annotationEnd.getChr())){	
						double[] s=libDataModelsP.get(st).scoreGene(annotationEnd);
						scores[0].put(st, s[USE_SCORE]);
						scores[1].put(st, s[PVAL_SCORE]);
					}
				}
			}
		}
		//For unoriented genes, sum up across both positive and negative models
		//IN THIS CASE, PVALUE WILL BE RECALCULATED
		else{
			String chr = annotationEnd.getChr();
			for(String st:libDataModelsN.keySet()){
				if(libDataModelsN.get(st).hasDataForChromosome(chr) && libDataModelsP.get(st).hasDataForChromosome(chr)){	
					double[] s = libDataModelsN.get(st).scoreGene(annotationEnd);
					double[] t=libDataModelsP.get(st).scoreGene(annotationEnd);
					scores[0].put(st, (t[USE_SCORE]+s[USE_SCORE]));
					double lambda = (libDataModelsP.get(st).getNumberOfReads(chr)+libDataModelsN.get(st).getNumberOfReads(chr))
							/libDataModelsP.get(st).getNumberMarkers(chr);
					double pval=ContinuousDataAlignmentModel.calculatePVal(new Double(t[USE_SCORE]+s[USE_SCORE]).intValue(), lambda, annotationEnd.getSize(), libDataModelsP.get(st).getNumberMarkers(chr));
					scores[1].put(st, pval);
				}
				else{
					if(libDataModelsN.get(st).hasDataForChromosome(chr)){
						double[] s=libDataModelsN.get(st).scoreGene(annotationEnd);
						scores[0].put(st,s[USE_SCORE]);
						scores[1].put(st,s[PVAL_SCORE]);
					}
					if(libDataModelsP.get(st).hasDataForChromosome(chr)){
						double[] s=libDataModelsP.get(st).scoreGene(annotationEnd);
						scores[0].put(st, s[USE_SCORE]);
						scores[1].put(st, s[PVAL_SCORE]);
					}
				}
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
			if(annotationEnd.isNegativeStrand()){
				if(!oppositeStrand){
					scores = libDataModelN.scanPRate(annotationEnd);
				}
				else
					scores = libDataModelP.scanPRate(annotationEnd);
			}
			else
				if(!oppositeStrand)
					scores = libDataModelP.scanPRate(annotationEnd);
				else
					scores = libDataModelN.scanPRate(annotationEnd);
		}
		else{
			String chr = annotationEnd.getChr();
			double[] scoresP = libDataModelP.scanPRate(annotationEnd);
			double[] scoresN = libDataModelN.scanPRate(annotationEnd);
			scores = new double[scoresP.length];
			//Everything is not additive
			//[0]c alculatePVal(new Double(sum).intValue(), getLambda(chr), count, getNumberMarkers(chr)) - take the minimum
			//[1] enrich - add
			//[2] sum - add 
			//[3] avgCoverage - add 
			//[4] rpkm - add 
			double lambda = (libDataModelP.getNumberOfReads(chr)+libDataModelN.getNumberOfReads(chr))
					/libDataModelP.getNumberMarkers(chr);
			double pval=ContinuousDataAlignmentModel.calculatePVal(new Double(scoresP[COUNT_SCORE]+scoresN[COUNT_SCORE]).intValue(), lambda, annotationEnd.getSize(), libDataModelP.getNumberMarkers(chr));
			scores[PVAL_SCORE]= pval;
			//scores[0] = Math.min(scoresP[0], scoresN[0]);
			for(int i=1;i<scores.length;i++){
				scores[i] = scoresP[i]+scoresN[i];
			}
		}

		return scores;
	}

	
	/**
	 * This function calculates the scores for all alignment models from a given window to the start of the gene and adds it to the matrix
	 * @param libDataModels
	 * @param annotationEnd
	 * @param annotation
	 * @param bestWindowEnd : value of this parameter is the end of the best window. If the window is inside the gene of interest, get the annotation upto this distance from the gene end. If the window is outside the gene, get the annotation upto this many bps outside the gene end. 
	 * @return  Vector of size libDataModels with score of each 
	 * @throws IOException
	 */
	private static MatrixWithHeaders getScoreFromWindowToGeneStart(Map<String,ContinuousDataAlignmentModel> libDataModels,RefSeqGene annotationEnd,RefSeqGene annotation,int bestWindowEnd,boolean insideGene,MatrixWithHeaders fullGeneScoreMatrix,String rowName) throws IOException{
		
		RefSeqGene subannotation;
		//1.If the annotationEnd is not within the RefSeq Gene
		if(annotationEnd.percentOverlapping(annotation)<1.0){
			//System.out.println("Annotation End is not inside the gene:False compared to "+new Boolean(insideGene).toString());
			//Then, Get subannotation from gene start to gene End, including the alignment from gene end to annotationEnd
			//subannotation = annotation.getExtended3primeIsoform(annotationEnd.absoluteToRelativePosition(annotation.getOrientedEnd())+annotationEnd.getTranscriptLength());
			subannotation = annotation.extendAnnotation(bestWindowEnd);
		}
		//Else get subannotation from gene start to annotationEnd window
		else{
		//	System.out.println("Annotation End is inside the gene:True compared to "+new Boolean(insideGene).toString());
			//subannotation = getSubAnnotationFromEnd(annotation, annotation.getTranscriptLength()-bestWindowEnd,bestWindowEnd-1);
			subannotation=annotation;
		}
	
		//Score the subAnnotation
		Map<String,Double> scores = new HashMap<String,Double>();
		for(String st:libDataModels.keySet()){
			//double compareScore = libDataModels[i].scoreGene(annotation)[USE_SCORE];
			if(libDataModels.get(st).hasDataForChromosome(subannotation.getChr())){	
				scores.put(st,libDataModels.get(st).scoreGene(subannotation)[USE_SCORE]);
			}
			else{
				scores.put(st, 0.0);
			}
		}
		
		for(int i=0;i<fullGeneScoreMatrix.columnDimension();i++){
			fullGeneScoreMatrix.set(rowName, fullGeneScoreMatrix.getColoumnName(i), scores.get(fullGeneScoreMatrix.getColoumnName(i)));
		}
		
		return fullGeneScoreMatrix;
	}
	
	/**
	 * This function calculates the scores for all alignment models from a given window to the start of the gene and adds it to the matrix
	 * @param libDataModels
	 * @param annotationStart
	 * @param annotation
	 * @param bestWindowEnd : value of this parameter is the end of the best window. If the window is inside the gene of interest, get the annotation upto this distance from the gene end. If the window is outside the gene, get the annotation upto this many bps outside the gene end. 
	 * @return  Vector of size libDataModels with score of each 
	 * @throws IOException
	 */
	private static MatrixWithHeaders getScoreFromWindowToGeneEnd(Map<String,ContinuousDataAlignmentModel> libDataModels,RefSeqGene annotationStart,RefSeqGene annotation,int bestWindowStart,boolean insideGene,MatrixWithHeaders fullGeneScoreMatrix,String rowName) throws IOException{

		RefSeqGene subannotation;
		//1.If the annotationStart is not within the RefSeq Gene
		if(annotationStart.percentOverlapping(annotation)<1.0){
			//System.out.println("Annotation End is not inside the gene:False compared to "+new Boolean(insideGene).toString());
			//Then, Get subannotation from gene start to gene End, including the alignment from gene end to annotationStart
			//subannotation = annotation.getExtended3primeIsoform(annotationStart.absoluteToRelativePosition(annotation.getOrientedEnd())+annotationStart.getTranscriptLength());
			subannotation = annotation.extendAnnotation(bestWindowStart);
		}
		//Else get subannotation from gene start to annotationStart window
		else{
		//	System.out.println("Annotation End is inside the gene:True compared to "+new Boolean(insideGene).toString());
			//subannotation = getSubAnnotationFromStart(annotation, annotation.getTranscriptLength()-bestWindowStart,bestWindowStart-1);
			subannotation=annotation;
		}
	
		//Score the subAnnotation
		Map<String,Double> scores = new HashMap<String,Double>();
		for(String ss:libDataModels.keySet()){
			//double compareScore = libDataModels[i].scoreGene(annotation)[USE_SCORE];
			if(libDataModels.get(ss).hasDataForChromosome(subannotation.getChr())){	
				scores.put(ss, libDataModels.get(ss).scoreGene(subannotation)[USE_SCORE]);
				
			}
			else{
				scores.put(ss, 0.0);
			}
		}
		
		for(int i=0;i<fullGeneScoreMatrix.columnDimension();i++){
			fullGeneScoreMatrix.set(rowName, fullGeneScoreMatrix.getColoumnName(i), scores.get(fullGeneScoreMatrix.getColoumnName(i)));
		}
		
		return fullGeneScoreMatrix;
	}
	
	/**
	 * This function calculates the scores for all alignment models from a given window to the start of the gene and adds it to the matrix
	 * @param libDataModels
	 * @param annotationStart
	 * @param annotation
	 * @param bestWindowEnd : value of this parameter is the end of the best window. If the window is inside the gene of interest, get the annotation upto this distance from the gene end. If the window is outside the gene, get the annotation upto this many bps outside the gene end. 
	 * @return  Vector of size libDataModels with score of each 
	 * @throws IOException
	 */
	private static MatrixWithHeaders getScoreFromWindowToGeneEnd(ContinuousDataAlignmentModel libDataModel,RefSeqGene annotationStart,RefSeqGene annotation,int bestWindowStart,boolean insideGene,MatrixWithHeaders fullGeneScoreMatrix,String rowName) throws IOException{
		
		RefSeqGene subannotation;
		//1.If the annotationStart is not within the RefSeq Gene
		if(annotationStart.percentOverlapping(annotation)<1.0){
			//System.out.println("Annotation End is not inside the gene:False compared to "+new Boolean(insideGene).toString());
			//Then, Get subannotation from gene start to gene End, including the alignment from gene end to annotationStart
			//subannotation = annotation.getExtended3primeIsoform(annotationStart.absoluteToRelativePosition(annotation.getOrientedEnd())+annotationStart.getTranscriptLength());
			subannotation = annotation.extendAnnotation(bestWindowStart);
		}
		//Else get subannotation from gene start to annotationStart window
		else{
		//	System.out.println("Annotation End is inside the gene:True compared to "+new Boolean(insideGene).toString());
			subannotation = getSubAnnotationFromStart(annotation, annotation.getTranscriptLength()-bestWindowStart,bestWindowStart-1);
		}
	
		//Score the subAnnotation
		double scores = 0.0;
		if(libDataModel.hasDataForChromosome(subannotation.getChr())){	
			scores = libDataModel.scoreGene(subannotation)[USE_SCORE];
		}
		else{
			scores = 0.0;
		}
		
		for(int i=0;i<fullGeneScoreMatrix.columnDimension();i++){
			fullGeneScoreMatrix.set(rowName, fullGeneScoreMatrix.getColoumnName(i), scores);
		}
		
		return fullGeneScoreMatrix;
	}

	private static MatrixWithHeaders getScoreFromWindowToGeneStart(ContinuousDataAlignmentModel libDataModel,RefSeqGene annotationEnd,RefSeqGene annotation,int bestWindowEnd,boolean insideGene,MatrixWithHeaders fullGeneScoreMatrix,String rowName) throws IOException{
		
		RefSeqGene subannotation;
		//1.If the annotationEnd is not within the RefSeq Gene
		if(annotationEnd.percentOverlapping(annotation)<1.0){
			//System.out.println("Annotation End is not inside the gene:False compared to "+new Boolean(insideGene).toString());
			//Then, Get subannotation from gene start to gene End, including the alignment from gene end to annotationEnd
			//subannotation = annotation.getExtended3primeIsoform(annotationEnd.absoluteToRelativePosition(annotation.getOrientedEnd())+annotationEnd.getTranscriptLength());
			subannotation = annotation.extendAnnotation(bestWindowEnd);
		}
		//Else get subannotation from gene start to annotationEnd window
		else{
		//	System.out.println("Annotation End is inside the gene:True compared to "+new Boolean(insideGene).toString());
			subannotation = getSubAnnotationFromEnd(annotation, annotation.getTranscriptLength()-bestWindowEnd,bestWindowEnd-1);
		}
	
		//Score the subAnnotation
		double scores = 0.0;
		//double compareScore = libDataModels[i].scoreGene(annotation)[USE_SCORE];
		if(libDataModel.hasDataForChromosome(subannotation.getChr())){	
			scores = libDataModel.scoreGene(subannotation)[USE_SCORE];
		}
		else{
			scores = 0.0;
		}
		for(int i=0;i<fullGeneScoreMatrix.columnDimension();i++){
			fullGeneScoreMatrix.set(rowName, fullGeneScoreMatrix.getColoumnName(i), scores);
		}
		
		return fullGeneScoreMatrix;
			
	}
	
	private static void writeCollapsedOutputForMultiple(BEDFileParser annotationCollection,String outputFile,MatrixWithHeaders resultMatrix,Map<String,RefSeqGene> geneToWindowMap) throws IOException{

		BEDFileParser collapsedAnnotations = annotationCollection.copy();
		String collapsedOutput = outputFile+".collapsed.table";
		BufferedWriter bw1 = new BufferedWriter(new FileWriter(outputFile+".collapsedGenes.bed"));
		BufferedWriter bw = new BufferedWriter(new FileWriter(collapsedOutput));
				
		for(String colName:resultMatrix.getColumnNames())
			bw.write("\t"+colName);
		bw.newLine();
		List<String> genesAdded = new ArrayList<String>();
		
		Map<String, IntervalTree<RefSeqGeneWithIsoforms>> mergedGenes = collapsedAnnotations.getMergedAnnotationMap(MIN_OVERLAP);
		Iterator<String> chrIt=mergedGenes.keySet().iterator();

		while(chrIt.hasNext()){
			String chr=chrIt.next();
			IntervalTree<RefSeqGeneWithIsoforms> tree=mergedGenes.get(chr);
			Iterator<RefSeqGeneWithIsoforms> geneIt=tree.valueIterator();
			while(geneIt.hasNext()) {
				RefSeqGeneWithIsoforms gene = geneIt.next();
				double scores[] = new double[resultMatrix.columnDimension()];
				for(int i=0;i<scores.length;i++){
					scores[i] = 0.0;//resultMatrix.get(duplicateNameMap.get(gene), i);
				}
					
				List<RefSeqGene> prevPeaks = new ArrayList<RefSeqGene>();
				gene.setName("gene");
				gene.cleanIsoforms();
				int isoformsAdded = 0;
				Iterator<RefSeqGeneWithIsoforms> overlappers = collapsedAnnotations.getOverlappers(gene).valueIterator();
				while(overlappers.hasNext()) {
					RefSeqGeneWithIsoforms overlapper = overlappers.next();
					BEDFileParser parser = new BEDFileParser();
					if(parser.isOverlapCompatible(gene, overlapper, MIN_OVERLAP)) {
						//System.err.println("Adding overlapper " + overlapper.toBED());
						Collection<RefSeqGene> isoforms = overlapper.getAllIsoforms();
						if(isoforms.size()>0){
								
							for(RefSeqGene iso: isoforms) {
								if(geneToWindowMap.containsKey(duplicateNameMap.get(iso))){
									gene.setName(gene.getName()+"_"+iso.getName());
									boolean couldAdd = gene.addContainedIsoform(iso);
									if(!couldAdd) {
										//System.err.println("WARN: Could not add isoform " + overlapper.toBED() + " to " + gene.toBED());
									} else {
										if(isoformsAdded==0){
											prevPeaks.add(geneToWindowMap.get(duplicateNameMap.get(iso)));
											for(int i=0;i<scores.length;i++){
												scores[i] = resultMatrix.get(duplicateNameMap.get(iso), i); 
												//System.err.println(scores[i]);
											}
										}
										else{
											//System.err.println(iso.toBED());
											boolean overLapsPeak = false;
											for(RefSeqGene peak:prevPeaks){
												if(peak.overlaps(geneToWindowMap.get(duplicateNameMap.get(iso))))
													overLapsPeak = true;
											}
											if(!overLapsPeak){
												//System.err.println("Adding: "+resultMatrix.get(duplicateNameMap.get(iso), 0));
												prevPeaks.add(geneToWindowMap.get(duplicateNameMap.get(iso)));
												for(int i=0;i<scores.length;i++)
													scores[i] += resultMatrix.get(duplicateNameMap.get(iso), i);
											}
										}
										isoformsAdded++;
										
										
										//System.err.println("Added isoform " + overlapper.getName() + " to " + gene.getName());
									}
								}
							}
						}
					}
				}
				if(isoformsAdded<1 || genesAdded.contains(gene.getName())){
					//Do nothing
				}
				else{
					genesAdded.add(gene.getName());
					bw.write(gene.getName()+"\t");
					for(int i=0;i<scores.length;i++)
						bw.write(scores[i]+"\t");
					bw.newLine();
					bw1.write(gene.toBED());
		  			bw1.newLine();
				}
				/*if(isoformsAdded ==0) {
						//bw.write(gene.getName()+"\t"+resultMatrix.get(duplicateNameMap.get(gene), column))
						//System.err.println("ERROR: Gene " + gene.getName() + " " + gene.toBED() + "  had no overlapping isoforms");
					}
				}if(geneToWindowMap.containsKey(duplicateNameMap.get(gene))){
					prevPeak = geneToWindowMap.get(duplicateNameMap.get(gene));
				}
				else{
					bw.write(gene.getName()+"\t");
					for(int i=0;i<scores.length;i++)
						bw.write("0.0\t");
					bw.newLine();
				}*/
			}
		}	
		bw.close();
		bw1.close();
	}
	
	private static void writeCollapsedOutputForSingle(BEDFileParser annotationCollection,String outputFile,MatrixWithHeaders resultMatrix,Map<String,RefSeqGene> geneToWindowMap) throws IOException{

		BEDFileParser collapsedAnnotations = annotationCollection.copy();
		String collapsedOutput = outputFile+".collapsed.table";
		BufferedWriter bw1 = new BufferedWriter(new FileWriter(outputFile+".collapsedGenes.bed"));
	  	 
		BufferedWriter bw = new BufferedWriter(new FileWriter(collapsedOutput));
		for(String colName:resultMatrix.getColumnNames()){
			if(!colName.equals(BEST_DIST_TO_END_COL))
				bw.write("\t"+colName);
		}
		bw.newLine();
		List<String> genesAdded = new ArrayList<String>();
		
		Map<String, IntervalTree<RefSeqGeneWithIsoforms>> mergedGenes = collapsedAnnotations.getMergedAnnotationMap(MIN_OVERLAP);
		Iterator<String> chrIt=mergedGenes.keySet().iterator();

		while(chrIt.hasNext()){
			String chr=chrIt.next();
			IntervalTree<RefSeqGeneWithIsoforms> tree=mergedGenes.get(chr);
			Iterator<RefSeqGeneWithIsoforms> geneIt=tree.valueIterator();
			while(geneIt.hasNext()) {
				RefSeqGeneWithIsoforms gene = geneIt.next();
				//-1 for the distance thing
				double scores[] = new double[resultMatrix.columnDimension()-1];
				for(int i=0;i<scores.length;i++){
					scores[i] = 0.0;//resultMatrix.get(duplicateNameMap.get(gene), i);
				}
					
				List<RefSeqGene> prevPeaks = new ArrayList<RefSeqGene>();
				gene.setName("gene");
				gene.cleanIsoforms();
				int isoformsAdded = 0;
				Iterator<RefSeqGeneWithIsoforms> overlappers = collapsedAnnotations.getOverlappers(gene).valueIterator();
				while(overlappers.hasNext()) {
					RefSeqGeneWithIsoforms overlapper = overlappers.next();
					BEDFileParser parser = new BEDFileParser();
					if(parser.isOverlapCompatible(gene, overlapper, MIN_OVERLAP)) {
						//System.err.println("Adding overlapper " + overlapper.toBED());
						Collection<RefSeqGene> isoforms = overlapper.getAllIsoforms();
						if(isoforms.size()>0){
								
							for(RefSeqGene iso: isoforms) {
								if(geneToWindowMap.containsKey(duplicateNameMap.get(iso))){
									gene.setName(gene.getName()+"_"+iso.getName());
									boolean couldAdd = gene.addContainedIsoform(iso);
									if(!couldAdd) {
										//System.err.println("WARN: Could not add isoform " + overlapper.toBED() + " to " + gene.toBED());
									} else {
										if(isoformsAdded==0){
											prevPeaks.add(geneToWindowMap.get(duplicateNameMap.get(iso)));
											for(int i=0;i<scores.length;i++){
												scores[i] = resultMatrix.get(duplicateNameMap.get(iso), i); 
												//System.err.println(scores[i]);
											}
										}
										else{
											//System.err.println(iso.toBED());
											boolean overLapsPeak = false;
											for(RefSeqGene peak:prevPeaks){
												if(peak.overlaps(geneToWindowMap.get(duplicateNameMap.get(iso))))
													overLapsPeak = true;
											}
											if(!overLapsPeak){
												//System.err.println("Adding: "+resultMatrix.get(duplicateNameMap.get(iso), 0));
												prevPeaks.add(geneToWindowMap.get(duplicateNameMap.get(iso)));
												for(int i=0;i<scores.length;i++)
													scores[i] += resultMatrix.get(duplicateNameMap.get(iso), i);
											}
										}
										isoformsAdded++;
										
										
										//System.err.println("Added isoform " + overlapper.getName() + " to " + gene.getName());
									}
								}
							}
						}
					}
				}
				if(isoformsAdded<1 || genesAdded.contains(gene.getName())){
					//Do nothing
				}
				else{
					genesAdded.add(gene.getName());
					bw.write(gene.getName()+"\t");
					for(int i=0;i<scores.length;i++)
						bw.write(scores[i]+"\t");
					bw.newLine();
		  			bw1.write(gene.toBED());
		  			bw1.newLine();
				}
				/*if(isoformsAdded ==0) {
						//bw.write(gene.getName()+"\t"+resultMatrix.get(duplicateNameMap.get(gene), column))
						//System.err.println("ERROR: Gene " + gene.getName() + " " + gene.toBED() + "  had no overlapping isoforms");
					}
				}if(geneToWindowMap.containsKey(duplicateNameMap.get(gene))){
					prevPeak = geneToWindowMap.get(duplicateNameMap.get(gene));
				}
				else{
					bw.write(gene.getName()+"\t");
					for(int i=0;i<scores.length;i++)
						bw.write("0.0\t");
					bw.newLine();
				}*/
			}
		}	
		bw.close();
		bw1.close();   
	}

	
	public static Map<String,Double> calculateNormalizationFactors(MatrixWithHeaders mat){
		
		int M = mat.getNumberColumns();
		int N = mat.getNumberRows();
		
		Map<String,Double> factors = new HashMap<String,Double>();
		double[] means = new double[N];
		//CALCULATE GEOMETRIC MEAN FOR ALL N GENES
		for(int i=0;i<N;i++){
			means[i] = 1.0;
			for(int j=0;j<M;j++){
				means[i] = means[i] * mat.get(i, j);
			}
			means[i] = Math.pow(means[i],((double)1.0/(double)M));
		}
		//For each sample
		for(String c:mat.getColumnNames()){
			
			double[] col = new double[N];
			for(int i=0;i<N;i++){
				col[i] = mat.get(i, c)/means[i]; 
			}
			factors.put(c, (double)1.0/calculateMedian(array2List(col)));
		}
		return factors;
		
	}
	
	private static double calculateMedian(List<Double> values)
	{
	    List<Double> nonInfValues = new ArrayList<Double>();
	 
	    for(Double v: values){
	    	
	    	if(v.isInfinite()||v.isNaN()){
	    		//Dont add
	    	}
	    	else{
	    		//System.out.println(v);
	    		nonInfValues.add(v);
	    	}
	    }
	    Collections.sort(nonInfValues);
	    if (nonInfValues.size() % 2 == 1){
	    	double median = nonInfValues.get((nonInfValues.size()+1)/2-1);
	    	if(median==0.0)
	    		return (1.0);
	    	else
	    		return median;
	    }
	    else
	    {
	    	double lower;
	    	double upper;
	    	if(nonInfValues.size()==0){
	    		lower= 1.0;
	    		upper = 1.0;
	    	}
	    	else{
	    		lower = nonInfValues.get(nonInfValues.size()/2-1);
	    		upper = nonInfValues.get(nonInfValues.size()/2);
	    	}
			double median = (lower + upper) / 2.0;
			return median;
	    }
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
		DGE d = new DGE(args);
		//MatrixWithHeaders m = new MatrixWithHeaders("t1");
		//writeNormalizedMatrix("t1.out", m);
		//DGE dummy = new DGE(args);

	}
}
