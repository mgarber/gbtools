package broad.pda.seq.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.reader.AlignmentQueryReader;
import org.broad.igv.sam.reader.BAMQueryReader;

import broad.core.annotation.AnnotationReader;
import broad.core.annotation.AnnotationReaderFactory;
import broad.core.annotation.BED;
import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.error.ParseException;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.gene.RefSeqGeneWithIsoforms;
import broad.pda.seq.fastq.FastqSequence;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;
import net.sf.samtools.util.CloseableIterator;

public class PolyAUtils {
	
	// ** CONSTANTS -- DON'T TOUCH THESE ** //
	private static final String AT_SIGN = "@";
	private static final char A = 'A';
	private static final char T = 'T';
	private static final char C = 'C';
	private static final char G = 'G';
	private static final char N = 'N';
	private static final short READ_STARTS_WITH_A_CODE = 0;
	private static final short READ_STARTS_WITH_T_CODE = 1;
	private static final short READ_ENDS_WITH_A_CODE = 2;
	private static final short READ_ENDS_WITH_T_CODE = 3;
	private static final String NEW_LINE = "\n";
	private static final String READS_WITH_POLYA_START_WITH_A_FILE = "readsWithPolyA.startA.fastq";
	private static final String READS_WITH_POLYA_START_WITH_T_FILE = "readsWithPolyA.startT.fastq";
	private static final String READS_WITH_POLYA_END_WITH_A_FILE = "readsWithPolyA.endA.fastq";
	private static final String READS_WITH_POLYA_END_WITH_T_FILE = "readsWithPolyA.endT.fastq";
	private static final String READS_WITH_POLYA_START_WITH_A_ALIGNMENT_FILE = "readsWithPolyA.startA.sam";
	private static final String READS_WITH_POLYA_START_WITH_T_ALIGNMENT_FILE = "readsWithPolyA.startT.sam";
	private static final String READS_WITH_POLYA_END_WITH_A_ALIGNMENT_FILE = "readsWithPolyA.endA.sam";
	private static final String READS_WITH_POLYA_END_WITH_T_ALIGNMENT_FILE = "readsWithPolyA.endT.sam";
	private static final String READS_WITH_POLYA_START_WITH_A_ALIGNMENT_SORTED_FILE = "readsWithPolyA.startA.sorted.sam";
	private static final String READS_WITH_POLYA_START_WITH_T_ALIGNMENT_SORTED_FILE = "readsWithPolyA.startT.sorted.sam";
	private static final String READS_WITH_POLYA_END_WITH_A_ALIGNMENT_SORTED_FILE = "readsWithPolyA.endA.sorted.sam";
	private static final String READS_WITH_POLYA_END_WITH_T_ALIGNMENT_SORTED_FILE = "readsWithPolyA.endT.sorted.sam";
	private static final String POLYA_READS_FILE = "polyA.fastq";
	private static final String POLYA_ALIGNMENT_FILE = "polyA.sam";
	private static final String POLYA_ALIGNMENT_SORTED_FILE = "polyA.sorted.sam";
	private static final char POLYA_QUALITY_VALUE = 'I';
	private static final String MATES_OF_FULL_POLYA_READS_A_FILE = "matesOfFullPolyAReads.A.fastq";
	private static final String MATES_OF_FULL_POLYA_READS_T_FILE = "matesOfFullPolyAReads.T.fastq";
	private static final String MATES_OF_FULL_POLYA_READS_A_ALIGNMENT_FILE = "matesOfFullPolyAReads.A.sam";
	private static final String MATES_OF_FULL_POLYA_READS_T_ALIGNMENT_FILE = "matesOfFullPolyAReads.T.sam";
	private static final short FULL_POLYA_READ_IS_A_CODE = 0;
	private static final short FULL_POLYA_READ_IS_T_CODE = 1;
	private static final String BED_FORMAT = "BED";
	private static final String COVERAGE_OF_LAST_EXONS_SPLIT_INTO_JOBS_JOB_SEPARATOR = ",";
	private static final String COVERAGE_OF_LAST_EXONS_SPLIT_INTO_JOBS_INFO_SEPARATOR = ":";
	private static final String COVERAGE_OF_LAST_EXONS_JOB_FILE_START = "coverage_";
	private static final String COVERAGE_OF_LAST_EXONS_JOB_BSUB_FILE_END = ".bsub";
	private static final String COVERAGE_OF_LAST_EXONS_JOB_COVERAGE_FILE_END = ".cov";
	private static final String TAB = "\t";
	private static final String COVERAGE_OF_LAST_EXONS_SUM_FROM_ALL_JOBS_FILE = "coverage_sum.cov";
	private static final String COMMA = ",";
	private static final String COLON = ":";
	private static final String SLASH = "/";
	private static final String SPACE = " ";
	private static final String WRITE_ANNOTATIONS_SPLIT_INTO_JOBS_JOB_SEPARATOR = ",";
	private static final String WRITE_ANNOTATIONS_SPLIT_INTO_JOBS_INFO_SEPARATOR = ":";
	private static final String WRITE_ANNOTATIONS_ALL_JOB_FILE_START = "allAlignments_";
	private static final String WRITE_ANNOTATIONS_FILTERED_JOB_FILE_START = "filteredAlignments_";
	private static final String WRITE_ANNOTATIONS_JOB_BSUB_FILE_END = ".bsub";
	private static final String WRITE_ANNOTATIONS_JOB_ANNOTATION_FILE_END = ".bed";
	private static final String WRITE_ANNOTATIONS_BSUB_JOB_FILE_START = "alignments_";
	private static final String WRITE_ANNOTATIONS_LIST_OF_ALIGNMENTS_FILE = "alignments.tmp";
	private static final String ASTERISK = "*";
	private static final String ASTERISK_LETTERS = "ASTRK";
	private static final String PLUS = "+";
	private static final String MINUS = "-";
	private static final String R = "R";
	private static final String L = "L";
	private static final String PARAMETERS_READ_FILE_SPLIT = "=";
	private static final String ANNOTATIONS_JOBS_DIRECTORY = "annotations/";
	private static final int ONE_MINUTE_IN_MILLISECS = 60000;
	private static final int TEN_MILLISECONDS = 10;
	private static final String BSUB_FILE_SUCCESS_INDICATOR = "Successfully completed.";
	private static final String YES = "YES";
	private static final String NO = "NO";
	private static final double PVALUE_SIGNIFICANCE_CUTOFF = 0.05;
	private static final String POLYA_IN_GENOME_FILTERING_METHOD_POLYA = "polya";
	private static final String POLYA_IN_GENOME_FILTERING_METHOD_SIMPLEREPEATS = "simplerepeats";
	private static final String POLYA_IN_GENOME_FILTERING_METHOD_BOTH = "both";
	// ** CONSTANTS ** //
	
	// ** PARAMETERS ** //
	private static String LOCATION_OF_SCRIPTURE_JAR = "/seq/regevlab/polya/runs/scripture.jar";
	private static String PATH_TO_MOST_FILES = "/seq/regevlab/polya/runs/";
	
	private static int MIN_LENGTH_OF_POLYA_WITHIN_READ = 25;
	private static int MIN_LENGTH_OF_READ_REMAINDER_FOR_POLYA_WITHIN_READ = 25;
	private static int MAX_MISMATCHES_ALLOWED_FOR_POLYA_WITHIN_READ = 3;
	
	private static int EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS_WITH_POLYA_OVERLAP = 50;
	private static int EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS_WITH_ORIGINAL_ALIGNMENT_OVERLAP = 1;
	private static int TYPICAL_ORIGINAL_ALIGNMENT_LENGTH = 76;
	private static int MAX_DISTANCE_FROM_ORIGINAL_ALIGNMENT_ENDPOINT_FOR_FILTERING_ANNOTATIONS = 10;
	
	private static int MAX_VALID_ALIGNMENTS_OUTPUT_PER_READ_WITH_POLYA = 3;
	
	private static int LENGTH_OF_POLYA_FOR_ALIGNMENT_OF_POLYA = 20;
	private static boolean ALIGN_POLYA_WITH_V_OPTION = false;
	private static int MAX_MISMATCHES_ALLOWED_FOR_POLYA_ALIGNMENT_V = 3;
	private static int POLYA_ALIGNMENT_SEED_LENGTH = 10;
	private static int MAX_MISMATCHES_FOR_POLYA_ALIGNMENT_IN_SEED = 3; // max number of mismatches in the seed (i.e., the first POLYA_ALIGNMENT_SEED_LENGTH bases) ; values must be 0, 1, 2, or 3
	private static int MAX_MISMATCHES_FOR_POLYA_ALIGNMENT_TOTAL = 6; // max number of mismatches throughout the entire alignment (this value must be greater than or equal to MAX_MISMATCHES_FOR_POLYA_ALIGNMENT_IN_SEED) ; for some reason, we may get a number of mismatches that is equal to one or two more than this value
	private static int MAX_PERMITTED_TOTAL_OF_QUALITY_VALUES = (((int)POLYA_QUALITY_VALUE) - 33) * MAX_MISMATCHES_FOR_POLYA_ALIGNMENT_TOTAL; // note 1: if POLYA_QUALITY_VALUE='I', as it should, then ((int)POLYA_QUALITY_VALUE - 33)=40 ; note 2: ((int)POLYA_QUALITY_VALUE - 33) represents the Phred quality score of POLYA_QUALITY_VALUE ; note 3: ((int)POLYA_QUALITY_VALUE - 33) must be greater than or equal to 30 ; for questions, see bowtie manual for the -e option
	
	private static int MIN_NUM_OF_POLYA_ALIGNMENTS_FOR_FILTERING_ALIGNMENTS = 12;
	private static int MIN_NUM_OF_ORIGINAL_ALIGNMENTS_FOR_FILTERING_ALIGNMENTS = 1;
	
	private static int MAX_MISMATCHES_ALLOWED_FOR_FULL_POLYA_READ = 10;
	
	private static int MAX_VALID_ALIGNMENTS_OUTPUT_PER_MATE_OF_FULL_POLYA_READ = 8;
	
	private static int EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_TRANSCRIPTS = 250;
	
	private static int EXTENSION_FACTOR_AFTER_LAST_EXON_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_LAST_EXON = 1000;
	
	private static int MAX_DISTANCE_BETWEEN_LAST_EXON_AND_NEXT_EXON_FOR_LINKING_TRANSCRIPTS = 500;
		
	private static int COVERAGE_OF_LAST_EXONS_EXTENSION_FACTOR = 0;
	private static int COVERAGE_OF_LAST_EXONS_NUM_OF_BASES_TO_READ = 500;
	private static int COVERAGE_OF_LAST_EXONS_MIN_NUM_OF_BASES = 500;
	private static int COVERAGE_OF_LAST_EXONS_WINDOW_SIZE = 10; // note: this should be chosen such that ((COVERAGE_OF_LAST_EXONS_NUM_OF_BASES_TO_READ + COVERAGE_OF_LAST_EXONS_EXTENSION_FACTOR) % COVERAGE_OF_LAST_EXONS_WINDOW_SIZE == 0)
	
	private static int EXTENSION_FACTOR_FOR_FILTERING_TRANSCRIPT_END_COORDINATES_USING_MATES_OF_FULL_POLYA_READS = 300;
	
	private static int SIMPLE_REPEAT_MIN_PERCENTANGE_OF_A_OR_T_FOR_FILTERING = 40;
	private static int EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS_WITH_SIMPLE_REPEAT_LOCATIONS_OVERLAP = 10;
	private static String POLYA_IN_GENOME_FILTERING_METHOD = "both";
	// ** PARAMETERS ** //
	
	// ** PARAMETERS / CONSTANTS ** //
	private static String READS_FILE = "/seq/regevlab/polya/data/all_reads.fastq";
	private static String BOWTIE_PATH = "/seq/mguttman/scripts/BowTie/bowtie-0.12.1/bowtie";
	private static String BOWTIE_REFERENCE_INDEX = "/seq/genome/mouse/mouse_Mm9/mm9.nonrandom";
	private static String IGVTOOLS_PATH = "/home/radon00/hmetsky/IGVTools/igvtools";
	private static String CHR_SIZES_FILE = "/seq/genome/mouse/mouse_Mm9/sizes";
	private static String ALL_ANNOTATIONS_OUT_FILE = "allAlignments.bed";
	private static String FILTERED_ANNOTATIONS_OUT_FILE = "filteredAlignments.bed";
	private static boolean SORT_AND_INDEX_READS = true;
	private static String SCRIPTURE_TRANSCRIPTS_FILE = "/seq/regevlab/polya/data/all_transcripts.bed";
	private static String MATES_OF_FULL_POLYA_READS_OUT_FILE = "matesOfFullPolyAReads.bed";
	private static String ALIGNMENTS_FILE = "/seq/regevlab/polya/data/all_alignments.sorted.sam";
	private static String END_COORDINATES_OUT_FILE = "endCoordinates.bed";
	private static String END_COORDINATES_ONLY_WITH_POLYA_MATES_OUT_FILE = "endCoordinates_onlyWithPolyAMates.bed";
	private static String TRANSCRIPT_LINKS_OUT_FILE = "transcriptLinks.txt";
	private static String MODIFIED_TRANSCRIPTS_OUT_FILE = "transcripts_new.bed";
	private static String MODIFIED_TRANSCRIPTS_ONLY_WITH_POLYA_MATES_OUT_FILE = "transcripts_new_onlyWithPolyAMates.bed";
	private static String MODIFIED_TRANSCRIPTS_ONLY_WITH_MODIFIED_OUT_FILE = "transcripts_new_onlyWithModified.bed";
	private static int NUM_OF_ALIGNMENTS_PER_JOB_FOR_WRITING_AND_FILTERING_ANNOTATIONS = 2500;
	private static String MEAN_QUALITY_SCORES_FOR_CORRECT_ORIENTATION_OUT_FILE = "meanQualityScoresForCorrectOrientationOfPolyA.txt";
	private static String MEAN_QUALITY_SCORES_FOR_INCORRECT_ORIENTATION_OUT_FILE = "meanQualityScoresForIncorrectOrientationOfPolyA.txt";
	
	private static boolean ANALYZE_AND_WRITE_ANALYSIS_OF_RESULTS = true;
	private static String TRANSCRIPTS_SCORED_FILE = "/seq/regevlab/polya/data/all_transcripts.scored.bed";
	private static String TRANSCRIPT_NAMES_OF_LABELED_TRANSCRIPTS_WITH_NUM_OF_ENDPOINTS_OUT_FILE = "transcriptNamesOfLabeledTranscriptsWithNumOfEndpoints.txt";
	private static String TRANSCRIPT_NAMES_OF_ALL_TRANSCRIPTS_WITH_WHETHER_THEY_ARE_LABELED_WITH_ENDPOINT_OUT_FILE = "transcriptNamesOfAllTranscriptsWithWhetherTheyAreLabeledWithEndpoint.txt";
	private static String PERCENTAGE_OF_ALL_TRANSCRIPTS_LABELED_WITH_ENDPOINT_OUT_FILE = "percentageOfAllTranscriptsLabeledWithEndpoint.txt";
	private static String TRANSCRIPT_NAMES_OF_LABELED_TRANSCRIPTS_WITH_RPKM_OUT_FILE = "transcriptNamesOfLabeledTranscriptsWithRPKM.txt";
	private static String TRANSCRIPT_NAMES_OF_LABELED_TRANSCRIPTS_WITH_PVALUE_OUT_FILE = "transcriptNamesOfLabeledTranscriptsWithPValue.txt";
	private static String PERCENTAGE_OF_LABELED_TRANSCRIPTS_WITH_SIGNIFICANT_PVALUE_OUT_FILE = "percentageOfLabeledTranscriptsWithSignificantPValue.txt";
	private static String ALL_ENDPOINTS_WITH_NUM_OF_READS_OUT_FILE = "allEndpointsWithNumOfReads.txt";
	private static String TRANSCRIPT_NAMES_WITH_NUM_OF_READS_FOR_THE_ENDPOINT_WITH_MOST_READS_OUT_FILE = "transcriptNamesWithNumOfReadsForTheEndpointWithMostReads.txt";
	
	private static String SIMPLE_REPEAT_LOCATIONS_FILE = "/seq/regevlab/polya/data/simpleRepeats_withPercentages.bed";
	private static String SIMPLE_REPEAT_LOCATIONS_FILTERED_FILE = "simpleRepeats.bed";
	
	private static String PATH_TO_CURRENT_DIRECTORY = null;
	// ** PARAMETERS / CONSTANTS ** //
	
	public static void initializePathToDirectory() throws IOException {
		PATH_TO_CURRENT_DIRECTORY = (new File(".")).getCanonicalPath() + SLASH;
	}
	
	public static void initializeParameters(String parametersFile) throws IOException {
		Map<String, String> parameterValues = new HashMap<String, String>();
		
		BufferedReader br = new BufferedReader(new FileReader(parametersFile));
		String line;
		while ((line = br.readLine()) != null) {
			 if (line.trim().length() == 0)
				 continue;
			
			String[] lineSplit = line.split(PARAMETERS_READ_FILE_SPLIT);
			String parameterName = lineSplit[0];
			String parameterValue = lineSplit[1];
			
			parameterValues.put(parameterName, parameterValue);
		}
		
		if (parameterValues.containsKey("LOCATION_OF_SCRIPTURE_JAR"))
			LOCATION_OF_SCRIPTURE_JAR = parameterValues.get("LOCATION_OF_SCRIPTURE_JAR");
		
		if (parameterValues.containsKey("PATH_TO_MOST_FILES")) {
			if (parameterValues.get("PATH_TO_MOST_FILES").equals("PATH_TO_CURRENT_DIRECTORY"))
				PATH_TO_MOST_FILES = PATH_TO_CURRENT_DIRECTORY;
			else
				PATH_TO_MOST_FILES = parameterValues.get("PATH_TO_MOST_FILES");
		}
		
		if (parameterValues.containsKey("MIN_LENGTH_OF_POLYA_WITHIN_READ"))
			MIN_LENGTH_OF_POLYA_WITHIN_READ = Integer.parseInt(parameterValues.get("MIN_LENGTH_OF_POLYA_WITHIN_READ"));
		
		if (parameterValues.containsKey("MIN_LENGTH_OF_READ_REMAINDER_FOR_POLYA_WITHIN_READ"))
			MIN_LENGTH_OF_READ_REMAINDER_FOR_POLYA_WITHIN_READ = Integer.parseInt(parameterValues.get("MIN_LENGTH_OF_READ_REMAINDER_FOR_POLYA_WITHIN_READ"));
		
		if (parameterValues.containsKey("MAX_MISMATCHES_ALLOWED_FOR_POLYA_WITHIN_READ"))
			MAX_MISMATCHES_ALLOWED_FOR_POLYA_WITHIN_READ = Integer.parseInt(parameterValues.get("MAX_MISMATCHES_ALLOWED_FOR_POLYA_WITHIN_READ"));
		
		if (parameterValues.containsKey("EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS_WITH_POLYA_OVERLAP"))
			EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS_WITH_POLYA_OVERLAP = Integer.parseInt(parameterValues.get("EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS_WITH_POLYA_OVERLAP"));
		
		if (parameterValues.containsKey("EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS_WITH_ORIGINAL_ALIGNMENT_OVERLAP"))
			EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS_WITH_ORIGINAL_ALIGNMENT_OVERLAP = Integer.parseInt(parameterValues.get("EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS_WITH_ORIGINAL_ALIGNMENT_OVERLAP"));
	
		if (parameterValues.containsKey("TYPICAL_ORIGINAL_ALIGNMENT_LENGTH"))
			TYPICAL_ORIGINAL_ALIGNMENT_LENGTH = Integer.parseInt(parameterValues.get("TYPICAL_ORIGINAL_ALIGNMENT_LENGTH"));
		
		if (parameterValues.containsKey("MAX_DISTANCE_FROM_ORIGINAL_ALIGNMENT_ENDPOINT_FOR_FILTERING_ANNOTATIONS"))
			MAX_DISTANCE_FROM_ORIGINAL_ALIGNMENT_ENDPOINT_FOR_FILTERING_ANNOTATIONS = Integer.parseInt(parameterValues.get("MAX_DISTANCE_FROM_ORIGINAL_ALIGNMENT_ENDPOINT_FOR_FILTERING_ANNOTATIONS"));
		
		if (parameterValues.containsKey("MAX_VALID_ALIGNMENTS_OUTPUT_PER_READ_WITH_POLYA"))
			MAX_VALID_ALIGNMENTS_OUTPUT_PER_READ_WITH_POLYA = Integer.parseInt(parameterValues.get("MAX_VALID_ALIGNMENTS_OUTPUT_PER_READ_WITH_POLYA"));
		
		if (parameterValues.containsKey("LENGTH_OF_POLYA_FOR_ALIGNMENT_OF_POLYA"))
			LENGTH_OF_POLYA_FOR_ALIGNMENT_OF_POLYA = Integer.parseInt(parameterValues.get("LENGTH_OF_POLYA_FOR_ALIGNMENT_OF_POLYA"));
		
		if (parameterValues.containsKey("ALIGN_POLYA_WITH_V_OPTION"))
			ALIGN_POLYA_WITH_V_OPTION = Boolean.parseBoolean(parameterValues.get("ALIGN_POLYA_WITH_V_OPTION"));
		
		if (parameterValues.containsKey("MAX_MISMATCHES_ALLOWED_FOR_POLYA_ALIGNMENT_V"))
			MAX_MISMATCHES_ALLOWED_FOR_POLYA_ALIGNMENT_V = Integer.parseInt(parameterValues.get("MAX_MISMATCHES_ALLOWED_FOR_POLYA_ALIGNMENT_V"));
		
		if (parameterValues.containsKey("POLYA_ALIGNMENT_SEED_LENGTH"))
			POLYA_ALIGNMENT_SEED_LENGTH = Integer.parseInt(parameterValues.get("POLYA_ALIGNMENT_SEED_LENGTH"));
		
		if (parameterValues.containsKey("MAX_MISMATCHES_FOR_POLYA_ALIGNMENT_IN_SEED"))
			MAX_MISMATCHES_FOR_POLYA_ALIGNMENT_IN_SEED = Integer.parseInt(parameterValues.get("MAX_MISMATCHES_FOR_POLYA_ALIGNMENT_IN_SEED"));
		
		if (parameterValues.containsKey("MAX_MISMATCHES_FOR_POLYA_ALIGNMENT_TOTAL"))
			MAX_MISMATCHES_FOR_POLYA_ALIGNMENT_TOTAL = Integer.parseInt(parameterValues.get("MAX_MISMATCHES_FOR_POLYA_ALIGNMENT_TOTAL"));
		
		if (parameterValues.containsKey("MAX_PERMITTED_TOTAL_OF_QUALITY_VALUES"))
			MAX_PERMITTED_TOTAL_OF_QUALITY_VALUES = Integer.parseInt(parameterValues.get("MAX_PERMITTED_TOTAL_OF_QUALITY_VALUES"));
		
		if (parameterValues.containsKey("MIN_NUM_OF_POLYA_ALIGNMENTS_FOR_FILTERING_ALIGNMENTS"))
			MIN_NUM_OF_POLYA_ALIGNMENTS_FOR_FILTERING_ALIGNMENTS = Integer.parseInt(parameterValues.get("MIN_NUM_OF_POLYA_ALIGNMENTS_FOR_FILTERING_ALIGNMENTS"));
		
		if (parameterValues.containsKey("MIN_NUM_OF_ORIGINAL_ALIGNMENTS_FOR_FILTERING_ALIGNMENTS"))
			MIN_NUM_OF_ORIGINAL_ALIGNMENTS_FOR_FILTERING_ALIGNMENTS = Integer.parseInt(parameterValues.get("MIN_NUM_OF_ORIGINAL_ALIGNMENTS_FOR_FILTERING_ALIGNMENTS"));
		
		if (parameterValues.containsKey("MAX_MISMATCHES_ALLOWED_FOR_FULL_POLYA_READ"))
			MAX_MISMATCHES_ALLOWED_FOR_FULL_POLYA_READ = Integer.parseInt(parameterValues.get("MAX_MISMATCHES_ALLOWED_FOR_FULL_POLYA_READ"));
		
		if (parameterValues.containsKey("MAX_VALID_ALIGNMENTS_OUTPUT_PER_MATE_OF_FULL_POLYA_READ"))
			MAX_VALID_ALIGNMENTS_OUTPUT_PER_MATE_OF_FULL_POLYA_READ = Integer.parseInt(parameterValues.get("MAX_VALID_ALIGNMENTS_OUTPUT_PER_MATE_OF_FULL_POLYA_READ"));
		
		if (parameterValues.containsKey("EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_TRANSCRIPTS"))
			EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_TRANSCRIPTS = Integer.parseInt(parameterValues.get("EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_TRANSCRIPTS"));
		
		if (parameterValues.containsKey("EXTENSION_FACTOR_AFTER_LAST_EXON_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_LAST_EXON"))
			EXTENSION_FACTOR_AFTER_LAST_EXON_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_LAST_EXON = Integer.parseInt(parameterValues.get("EXTENSION_FACTOR_AFTER_LAST_EXON_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_LAST_EXON"));
		
		if (parameterValues.containsKey("MAX_DISTANCE_BETWEEN_LAST_EXON_AND_NEXT_EXON_FOR_LINKING_TRANSCRIPTS"))
			MAX_DISTANCE_BETWEEN_LAST_EXON_AND_NEXT_EXON_FOR_LINKING_TRANSCRIPTS = Integer.parseInt(parameterValues.get("MAX_DISTANCE_BETWEEN_LAST_EXON_AND_NEXT_EXON_FOR_LINKING_TRANSCRIPTS"));
		
		if (parameterValues.containsKey("COVERAGE_OF_LAST_EXONS_EXTENSION_FACTOR"))
			COVERAGE_OF_LAST_EXONS_EXTENSION_FACTOR = Integer.parseInt(parameterValues.get("COVERAGE_OF_LAST_EXONS_EXTENSION_FACTOR"));
		
		if (parameterValues.containsKey("COVERAGE_OF_LAST_EXONS_NUM_OF_BASES_TO_READ"))
			COVERAGE_OF_LAST_EXONS_NUM_OF_BASES_TO_READ = Integer.parseInt(parameterValues.get("COVERAGE_OF_LAST_EXONS_NUM_OF_BASES_TO_READ"));
		
		if (parameterValues.containsKey("COVERAGE_OF_LAST_EXONS_MIN_NUM_OF_BASES"))
			COVERAGE_OF_LAST_EXONS_MIN_NUM_OF_BASES = Integer.parseInt(parameterValues.get("COVERAGE_OF_LAST_EXONS_MIN_NUM_OF_BASES"));
		
		if (parameterValues.containsKey("COVERAGE_OF_LAST_EXONS_WINDOW_SIZE"))
			COVERAGE_OF_LAST_EXONS_WINDOW_SIZE = Integer.parseInt(parameterValues.get("COVERAGE_OF_LAST_EXONS_WINDOW_SIZE"));
		
		if (parameterValues.containsKey("EXTENSION_FACTOR_FOR_FILTERING_TRANSCRIPT_END_COORDINATES_USING_MATES_OF_FULL_POLYA_READS"))
			EXTENSION_FACTOR_FOR_FILTERING_TRANSCRIPT_END_COORDINATES_USING_MATES_OF_FULL_POLYA_READS = Integer.parseInt(parameterValues.get("EXTENSION_FACTOR_FOR_FILTERING_TRANSCRIPT_END_COORDINATES_USING_MATES_OF_FULL_POLYA_READS"));
	
		if (parameterValues.containsKey("READS_FILE"))
			READS_FILE = parameterValues.get("READS_FILE");
		
		if (parameterValues.containsKey("BOWTIE_PATH"))
			BOWTIE_PATH = parameterValues.get("BOWTIE_PATH");
		
		if (parameterValues.containsKey("BOWTIE_REFERENCE_INDEX"))
			BOWTIE_REFERENCE_INDEX = parameterValues.get("BOWTIE_REFERENCE_INDEX");
		
		if (parameterValues.containsKey("IGVTOOLS_PATH"))
			IGVTOOLS_PATH = parameterValues.get("IGVTOOLS_PATH");
		
		if (parameterValues.containsKey("CHR_SIZES_FILE"))
			CHR_SIZES_FILE = parameterValues.get("CHR_SIZES_FILE");
		
		if (parameterValues.containsKey("ALL_ANNOTATIONS_OUT_FILE"))
			ALL_ANNOTATIONS_OUT_FILE = parameterValues.get("ALL_ANNOTATIONS_OUT_FILE");
		
		if (parameterValues.containsKey("FILTERED_ANNOTATIONS_OUT_FILE"))
			FILTERED_ANNOTATIONS_OUT_FILE = parameterValues.get("FILTERED_ANNOTATIONS_OUT_FILE");
		
		if (parameterValues.containsKey("SORT_AND_INDEX_READS"))
			SORT_AND_INDEX_READS = Boolean.parseBoolean(parameterValues.get("SORT_AND_INDEX_READS"));
		
		if (parameterValues.containsKey("SCRIPTURE_TRANSCRIPTS_FILE"))
			SCRIPTURE_TRANSCRIPTS_FILE = parameterValues.get("SCRIPTURE_TRANSCRIPTS_FILE");
		
		if (parameterValues.containsKey("MATES_OF_FULL_POLYA_READS_OUT_FILE"))
			MATES_OF_FULL_POLYA_READS_OUT_FILE = parameterValues.get("MATES_OF_FULL_POLYA_READS_OUT_FILE");
		
		if (parameterValues.containsKey("FILTERED_ANNOTATIONS_OUT_FILE"))
			FILTERED_ANNOTATIONS_OUT_FILE = parameterValues.get("FILTERED_ANNOTATIONS_OUT_FILE");
		
		if (parameterValues.containsKey("ALIGNMENTS_FILE"))
			ALIGNMENTS_FILE = parameterValues.get("ALIGNMENTS_FILE");
		
		if (parameterValues.containsKey("END_COORDINATES_OUT_FILE"))
			END_COORDINATES_OUT_FILE = parameterValues.get("END_COORDINATES_OUT_FILE");
		
		if (parameterValues.containsKey("END_COORDINATES_ONLY_WITH_POLYA_MATES_OUT_FILE"))
			END_COORDINATES_ONLY_WITH_POLYA_MATES_OUT_FILE = parameterValues.get("END_COORDINATES_ONLY_WITH_POLYA_MATES_OUT_FILE");
		
		if (parameterValues.containsKey("TRANSCRIPT_LINKS_OUT_FILE"))
			TRANSCRIPT_LINKS_OUT_FILE = parameterValues.get("TRANSCRIPT_LINKS_OUT_FILE");
		
		if (parameterValues.containsKey("MODIFIED_TRANSCRIPTS_OUT_FILE"))
			MODIFIED_TRANSCRIPTS_OUT_FILE = parameterValues.get("MODIFIED_TRANSCRIPTS_OUT_FILE");
		
		if (parameterValues.containsKey("MODIFIED_TRANSCRIPTS_ONLY_WITH_POLYA_MATES_OUT_FILE"))
			MODIFIED_TRANSCRIPTS_ONLY_WITH_POLYA_MATES_OUT_FILE = parameterValues.get("MODIFIED_TRANSCRIPTS_ONLY_WITH_POLYA_MATES_OUT_FILE");
		
		if (parameterValues.containsKey("MODIFIED_TRANSCRIPTS_ONLY_WITH_MODIFIED_OUT_FILE"))
			MODIFIED_TRANSCRIPTS_ONLY_WITH_MODIFIED_OUT_FILE = parameterValues.get("MODIFIED_TRANSCRIPTS_ONLY_WITH_MODIFIED_OUT_FILE");
		
		if (parameterValues.containsKey("NUM_OF_ALIGNMENTS_PER_JOB_FOR_WRITING_AND_FILTERING_ANNOTATIONS"))
			NUM_OF_ALIGNMENTS_PER_JOB_FOR_WRITING_AND_FILTERING_ANNOTATIONS = Integer.parseInt(parameterValues.get("NUM_OF_ALIGNMENTS_PER_JOB_FOR_WRITING_AND_FILTERING_ANNOTATIONS"));
		
		if (parameterValues.containsKey("MEAN_QUALITY_SCORES_FOR_CORRECT_ORIENTATION_OUT_FILE"))
			MEAN_QUALITY_SCORES_FOR_CORRECT_ORIENTATION_OUT_FILE = parameterValues.get("MEAN_QUALITY_SCORES_FOR_CORRECT_ORIENTATION_OUT_FILE");
		
		if (parameterValues.containsKey("MEAN_QUALITY_SCORES_FOR_INCORRECT_ORIENTATION_OUT_FILE"))
			MEAN_QUALITY_SCORES_FOR_INCORRECT_ORIENTATION_OUT_FILE = parameterValues.get("MEAN_QUALITY_SCORES_FOR_INCORRECT_ORIENTATION_OUT_FILE");
		
		if (parameterValues.containsKey("ANALYZE_AND_WRITE_ANALYSIS_OF_RESULTS"))
			ANALYZE_AND_WRITE_ANALYSIS_OF_RESULTS = Boolean.parseBoolean(parameterValues.get("ANALYZE_AND_WRITE_ANALYSIS_OF_RESULTS"));
		
		if (parameterValues.containsKey("TRANSCRIPTS_SCORED_FILE"))
			TRANSCRIPTS_SCORED_FILE = parameterValues.get("TRANSCRIPTS_SCORED_FILE");
		
		if (parameterValues.containsKey("TRANSCRIPT_NAMES_OF_LABELED_TRANSCRIPTS_WITH_NUM_OF_ENDPOINTS_OUT_FILE"))
			TRANSCRIPT_NAMES_OF_LABELED_TRANSCRIPTS_WITH_NUM_OF_ENDPOINTS_OUT_FILE = parameterValues.get("TRANSCRIPT_NAMES_OF_LABELED_TRANSCRIPTS_WITH_NUM_OF_ENDPOINTS_OUT_FILE");
		
		if (parameterValues.containsKey("TRANSCRIPT_NAMES_OF_ALL_TRANSCRIPTS_WITH_WHETHER_THEY_ARE_LABELED_WITH_ENDPOINT_OUT_FILE"))
			TRANSCRIPT_NAMES_OF_ALL_TRANSCRIPTS_WITH_WHETHER_THEY_ARE_LABELED_WITH_ENDPOINT_OUT_FILE = parameterValues.get("TRANSCRIPT_NAMES_OF_ALL_TRANSCRIPTS_WITH_WHETHER_THEY_ARE_LABELED_WITH_ENDPOINT_OUT_FILE");
		
		if (parameterValues.containsKey("PERCENTAGE_OF_ALL_TRANSCRIPTS_LABELED_WITH_ENDPOINT_OUT_FILE"))
			PERCENTAGE_OF_ALL_TRANSCRIPTS_LABELED_WITH_ENDPOINT_OUT_FILE = parameterValues.get("PERCENTAGE_OF_ALL_TRANSCRIPTS_LABELED_WITH_ENDPOINT_OUT_FILE");
		
		if (parameterValues.containsKey("TRANSCRIPT_NAMES_OF_LABELED_TRANSCRIPTS_WITH_RPKM_OUT_FILE"))
			TRANSCRIPT_NAMES_OF_LABELED_TRANSCRIPTS_WITH_RPKM_OUT_FILE = parameterValues.get("TRANSCRIPT_NAMES_OF_LABELED_TRANSCRIPTS_WITH_RPKM_OUT_FILE");
		
		if (parameterValues.containsKey("TRANSCRIPT_NAMES_OF_LABELED_TRANSCRIPTS_WITH_PVALUE_OUT_FILE"))
			TRANSCRIPT_NAMES_OF_LABELED_TRANSCRIPTS_WITH_PVALUE_OUT_FILE = parameterValues.get("TRANSCRIPT_NAMES_OF_LABELED_TRANSCRIPTS_WITH_PVALUE_OUT_FILE");
		
		if (parameterValues.containsKey("PERCENTAGE_OF_LABELED_TRANSCRIPTS_WITH_SIGNIFICANT_PVALUE_OUT_FILE"))
			PERCENTAGE_OF_LABELED_TRANSCRIPTS_WITH_SIGNIFICANT_PVALUE_OUT_FILE = parameterValues.get("PERCENTAGE_OF_LABELED_TRANSCRIPTS_WITH_SIGNIFICANT_PVALUE_OUT_FILE");
		
		if (parameterValues.containsKey("ALL_ENDPOINTS_WITH_NUM_OF_READS_OUT_FILE"))
			ALL_ENDPOINTS_WITH_NUM_OF_READS_OUT_FILE = parameterValues.get("ALL_ENDPOINTS_WITH_NUM_OF_READS_OUT_FILE");
		
		if (parameterValues.containsKey("TRANSCRIPT_NAMES_WITH_NUM_OF_READS_FOR_THE_ENDPOINT_WITH_MOST_READS_OUT_FILE"))
			TRANSCRIPT_NAMES_WITH_NUM_OF_READS_FOR_THE_ENDPOINT_WITH_MOST_READS_OUT_FILE = parameterValues.get("TRANSCRIPT_NAMES_WITH_NUM_OF_READS_FOR_THE_ENDPOINT_WITH_MOST_READS_OUT_FILE");
		
		if (parameterValues.containsKey("SIMPLE_REPEAT_MIN_PERCENTANGE_OF_A_OR_T_FOR_FILTERING"))
			SIMPLE_REPEAT_MIN_PERCENTANGE_OF_A_OR_T_FOR_FILTERING = Integer.parseInt(parameterValues.get("SIMPLE_REPEAT_MIN_PERCENTANGE_OF_A_OR_T_FOR_FILTERING"));
		
		if (parameterValues.containsKey("POLYA_IN_GENOME_FILTERING_METHOD"))
			POLYA_IN_GENOME_FILTERING_METHOD = parameterValues.get("POLYA_IN_GENOME_FILTERING_METHOD").toLowerCase();
		
		if (parameterValues.containsKey("SIMPLE_REPEAT_LOCATIONS_FILE"))
			SIMPLE_REPEAT_LOCATIONS_FILE = parameterValues.get("SIMPLE_REPEAT_LOCATIONS_FILE");
		
		if (parameterValues.containsKey("SIMPLE_REPEAT_LOCATIONS_FILTERED_FILE"))
			SIMPLE_REPEAT_LOCATIONS_FILTERED_FILE = parameterValues.get("SIMPLE_REPEAT_LOCATIONS_FILTERED_FILE");
		
		if (parameterValues.containsKey("EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS_WITH_SIMPLE_REPEAT_LOCATIONS_OVERLAP"))
			EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS_WITH_SIMPLE_REPEAT_LOCATIONS_OVERLAP = Integer.parseInt(parameterValues.get("EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS_WITH_SIMPLE_REPEAT_LOCATIONS_OVERLAP"));
		
		MAX_PERMITTED_TOTAL_OF_QUALITY_VALUES = (((int)POLYA_QUALITY_VALUE) - 33) * MAX_MISMATCHES_FOR_POLYA_ALIGNMENT_TOTAL;
	}
	
	public static void setPathToCurrentDirectoryAndPathToMostFiles(String currentDirectory) {
		PATH_TO_CURRENT_DIRECTORY = currentDirectory;
		PATH_TO_MOST_FILES = currentDirectory;
	}
	
	public static void findAndProcessReadsWithPolyA() throws IOException {
		System.out.println("*** STARTING TO FIND AND PROCESS READS WITH POLYA ***");
		
		PolyAUtils.findAndWriteReadsWithPolyA(READS_FILE);
		System.out.println("Found and wrote reads with polyA");
		
		PolyAUtils.alignReadsWithPolyA(BOWTIE_PATH, BOWTIE_REFERENCE_INDEX);
		System.out.println("Aligned reads with polyA");
		
		System.out.println("*** DONE FINDING AND PROCESSING READS WITH POLYA ***\n\n++++++++++++++++++++\n\n");
	}
	
	public static void alignPolyAAndSortAndIndexAlignments() throws IOException {
		System.out.println("*** STARTING TO ALIGN POLYA AND SORT AND INDEX ALIGNMENTS ***");
		
		boolean shouldUsePolyAAlignmentFilteringMethod = false;
		boolean shouldUseSimpleRepeatsFilteringMethod = false;
		if (POLYA_IN_GENOME_FILTERING_METHOD.equals(POLYA_IN_GENOME_FILTERING_METHOD_POLYA))
			shouldUsePolyAAlignmentFilteringMethod = true;
		if (POLYA_IN_GENOME_FILTERING_METHOD.equals(POLYA_IN_GENOME_FILTERING_METHOD_SIMPLEREPEATS))
			shouldUseSimpleRepeatsFilteringMethod = true;
		if (POLYA_IN_GENOME_FILTERING_METHOD.equals(POLYA_IN_GENOME_FILTERING_METHOD_BOTH)) {
			shouldUsePolyAAlignmentFilteringMethod = true;
			shouldUseSimpleRepeatsFilteringMethod = true;
		}
		
		if (shouldUsePolyAAlignmentFilteringMethod) {
			PolyAUtils.alignPolyA(BOWTIE_PATH, BOWTIE_REFERENCE_INDEX);
			System.out.println("Aligned polyA");
		} else {
			System.out.println("Not aligning polyA because it is not set in filtering method...");
		}
		
		PolyAUtils.sortAndIndexAlignments(IGVTOOLS_PATH, SORT_AND_INDEX_READS);
		System.out.println("Sorted and indexed alignment files");
		
		System.out.println("*** DONE ALIGNING POLYA AND SORTING AND INDEXING ALIGNMENTS ***\n\n++++++++++++++++++++\n\n");
	}
	
	public static void writeAnnotationsOfReadsWithPolyAAndFilterAndWriteFilteredAnnotationsToBedFileForJob(String alignmentsForJob, int jobNum) throws IOException {
		PolyAUtils.writeAnnotationsOfReadsWithPolyAAndFilterAndWriteFilteredAnnotationsToBedFileForJob(ALIGNMENTS_FILE, CHR_SIZES_FILE, alignmentsForJob, jobNum);
	}
	
	public static void writeAnnotationsOfReadsWithPolyAAndFilterAndWriteFilteredAnnotationsToBedFileBySplittingIntoJobs(String parametersFile) throws IOException {
		System.out.println("*** STARTING TO WRITE AND FILTER ANNOTATIONS OF READS WITH POLYA ***");
		
		PolyAUtils.writeAnnotationsOfReadsWithPolyAAndFilterAndWriteFilteredAnnotationsToBedFileBySplittingIntoJobs(ALIGNMENTS_FILE, CHR_SIZES_FILE, NUM_OF_ALIGNMENTS_PER_JOB_FOR_WRITING_AND_FILTERING_ANNOTATIONS, parametersFile);
		
		System.out.println("*** DONE WRITING AND FILTERING ANNOTATIONS OF READS WITH POLYA ***\n\n++++++++++++++++++++\n\n");
	}
	
	public static void findAlignAndWriteMatesOfFullPolyAReads() throws IOException {
		System.out.println("*** STARTING TO FIND, ALIGN, AND WRITE MATES OF FULL POLYA READS ***");
		
		PolyAUtils.findAndWriteMatesOfFullPolyAReads(READS_FILE);
		System.out.println("Found and wrote mates of full polyA reads");
		PolyAUtils.alignMatesOfFullPolyAReads(BOWTIE_PATH, BOWTIE_REFERENCE_INDEX);
		System.out.println("Aligned mates of full polyA reads");
		PolyAUtils.writeAlignmentsOfMatesOfFullPolyAReadsToBedFile(MATES_OF_FULL_POLYA_READS_OUT_FILE);
		System.out.println("Wrote alignments of mates of full polyA reads");
		
		System.out.println("*** DONE FINDING, ALIGNING, AND WRITING MATES OF FULL POLYA READS ***\n\n++++++++++++++++++++\n\n");

	}
	
	public static void filterSetOfAlingnmentsByOverlapWithTranscripts() throws IOException, ParseException {
		PolyAUtils.filterSetOfAlignmentsByOverlapWithTranscripts(FILTERED_ANNOTATIONS_OUT_FILE, ALL_ANNOTATIONS_OUT_FILE, SCRIPTURE_TRANSCRIPTS_FILE);
	}
	
	public static void findAndProcessCoordinatesOfEndOfTranscripts() throws IOException, ParseException {
		System.out.println("*** STARTING TO FIND AND PROCESS COORDINATES OF END OF TRANSCRIPTS ***");
		
		BEDFileParser transcripts = PolyAUtils.initiateTranscriptsParser(SCRIPTURE_TRANSCRIPTS_FILE);
		Map<String, Set<RefSeqGene>> m = PolyAUtils.findCoordinatesOfEndsOfTranscripts(transcripts, FILTERED_ANNOTATIONS_OUT_FILE, SCRIPTURE_TRANSCRIPTS_FILE, CHR_SIZES_FILE);
		PolyAUtils.processEndOfTranscriptCoordinates(transcripts, SCRIPTURE_TRANSCRIPTS_FILE, m, MATES_OF_FULL_POLYA_READS_OUT_FILE, END_COORDINATES_OUT_FILE, END_COORDINATES_ONLY_WITH_POLYA_MATES_OUT_FILE, TRANSCRIPT_LINKS_OUT_FILE, MODIFIED_TRANSCRIPTS_OUT_FILE, MODIFIED_TRANSCRIPTS_ONLY_WITH_POLYA_MATES_OUT_FILE, MODIFIED_TRANSCRIPTS_ONLY_WITH_MODIFIED_OUT_FILE, ANALYZE_AND_WRITE_ANALYSIS_OF_RESULTS, TRANSCRIPTS_SCORED_FILE, TRANSCRIPT_NAMES_OF_LABELED_TRANSCRIPTS_WITH_NUM_OF_ENDPOINTS_OUT_FILE, TRANSCRIPT_NAMES_OF_ALL_TRANSCRIPTS_WITH_WHETHER_THEY_ARE_LABELED_WITH_ENDPOINT_OUT_FILE, PERCENTAGE_OF_ALL_TRANSCRIPTS_LABELED_WITH_ENDPOINT_OUT_FILE, TRANSCRIPT_NAMES_OF_LABELED_TRANSCRIPTS_WITH_RPKM_OUT_FILE, TRANSCRIPT_NAMES_OF_LABELED_TRANSCRIPTS_WITH_PVALUE_OUT_FILE, PERCENTAGE_OF_LABELED_TRANSCRIPTS_WITH_SIGNIFICANT_PVALUE_OUT_FILE, FILTERED_ANNOTATIONS_OUT_FILE, CHR_SIZES_FILE, ALL_ENDPOINTS_WITH_NUM_OF_READS_OUT_FILE, TRANSCRIPT_NAMES_WITH_NUM_OF_READS_FOR_THE_ENDPOINT_WITH_MOST_READS_OUT_FILE);
	
		System.out.println("*** DONE FINDING AND PROCESSING COORDINATES OF END OF TRANSCRIPTS ***\n\n++++++++++++++++++++\n\n");
	}
	
	public static void writePlotOfQualityScoresForPolyABasedOnMatchingOrientationWithTranscript() throws IOException, ParseException {
		System.out.println("*** STARTING TO WRITE PLOT OF QUALITY SCORES FOR POLYA BASED ON MATCHING ORIENTATION WITH TRANSCRIPT ***");
		
		BEDFileParser transcripts = PolyAUtils.initiateTranscriptsParser(SCRIPTURE_TRANSCRIPTS_FILE);
		PolyAUtils.writePlotOfQualityScoresForPolyABasedOnMatchingOrientationWithTranscript(transcripts, FILTERED_ANNOTATIONS_OUT_FILE, SCRIPTURE_TRANSCRIPTS_FILE, CHR_SIZES_FILE, READS_FILE, MEAN_QUALITY_SCORES_FOR_CORRECT_ORIENTATION_OUT_FILE, MEAN_QUALITY_SCORES_FOR_INCORRECT_ORIENTATION_OUT_FILE);
	
		System.out.println("*** DONE WRITING PLOT OF QUALITY SCORES FOR POLYA BASED ON MATCHING ORIENTATION WITH TRANSCRIPT ***\n\n++++++++++++++++++++\n\n");
	}
	
	/**
	public static void PIPELINE(String readFile, String bowtiePath, String bowtieReferenceIndex, String igvToolsPath, boolean sortAndIndexReads, String chrSizesFile, String allAnnotationsOutFile, String filteredAnnotationsOutFile, String matesOfFullPolyAReadsOutFile, String scriptureTranscriptsFile, String allAlignmentsInFile) throws IOException, ParseException {
		
		// PolyAUtils.findAndWriteReadsWithPolyA(readFile);
		System.out.println("Found and wrote reads with polyA");
		
		// PolyAUtils.alignReadsWithPolyA(bowtiePath, bowtieReferenceIndex);
		System.out.println("Aligned reads with polyA");
		
		PolyAUtils.alignPolyA(bowtiePath, bowtieReferenceIndex);
		System.out.println("Aligned polyA");
		
		PolyAUtils.sortAndIndexAlignments(igvToolsPath, sortAndIndexReads);
		System.out.println("Sorted and indexed alignment files");
		
		// PolyAUtils.writeAnnotationsOfReadsWithPolyAAndFilterAndWriteFilteredAnnotationsToBedFile(allAnnotationsOutFile, filteredAnnotationsOutFile, allAlignmentsInFile, chrSizesFile);
		System.out.println("Wrote (all and filtered) alignments to BED files");
		
		// PolyAUtils.findAndWriteMatesOfFullPolyAReads(readFile);
		System.out.println("Found and wrote mates of full polyA reads");
		// PolyAUtils.alignMatesOfFullPolyAReads(bowtiePath, bowtieReferenceIndex);
		System.out.println("Aligned mates of full polyA reads");
		// PolyAUtils.writeAlignmentsOfMatesOfFullPolyAReadsToBedFile(matesOfFullPolyAReadsOutFile);
		System.out.println("Wrote alignments of mates of full polyA reads");
		
		// PolyAUtils.filterSetOfAlignmentsByOverlapWithTranscripts(filteredAlignmentsOutFile, allAlignmentsOutFile, scriptureTranscriptsFile);
	}
	*/
	
	private static void filterSimpleRepeatLocationsFile(String simpleRepeatsFileIn, String simpleRepeatsFileOut) throws IOException {
		FileWriter fw = new FileWriter(simpleRepeatsFileOut);
		
		BufferedReader br = new BufferedReader(new FileReader(simpleRepeatsFileIn));
		String line;
		while ((line = br.readLine()) != null && line.trim().length() > 0) {
			String[] lineSplit = line.split(TAB);
			String chr = lineSplit[0];
			String start = lineSplit[1];
			String end = lineSplit[2];
			int percentage_A = Integer.parseInt(lineSplit[3]);
			int percentage_T = Integer.parseInt(lineSplit[4]);
			
			if (percentage_A >= SIMPLE_REPEAT_MIN_PERCENTANGE_OF_A_OR_T_FOR_FILTERING || percentage_T >= SIMPLE_REPEAT_MIN_PERCENTANGE_OF_A_OR_T_FOR_FILTERING)
				fw.write(chr + TAB + start + TAB + end + NEW_LINE);
		}
		br.close();
		
		fw.close();
	}
	
	public static int getNumOfLastExonsWithCorrectSize(String transcriptsFile) throws IOException, ParseException {
		AnnotationReader<? extends GenomicAnnotation> transcripts = AnnotationReaderFactory.create(transcriptsFile, BED_FORMAT);
		Map<GenomicAnnotation, Integer> lastExons = AnnotationReader.getLastExons(transcripts.getAnnotationList());
		
		int count = 0;
		for (GenomicAnnotation lastExon : lastExons.keySet()) {
			int startIndex = lastExon.getStart();
			int endIndex = lastExon.getEnd();
			
			if (endIndex - startIndex >= COVERAGE_OF_LAST_EXONS_MIN_NUM_OF_BASES)
				count++;
		}
		
		return count;
	}
	
	public static void sumReadCoverageAndAverageReadCoveragePercentagesNearLastExonsFromAllJobs(int numOfBasesToReadIn, int windowSizeIn) throws IOException {
		int numOfBasesToRead;
		if (numOfBasesToReadIn > 0) {
			numOfBasesToRead = numOfBasesToReadIn;
		} else {
			numOfBasesToRead = COVERAGE_OF_LAST_EXONS_NUM_OF_BASES_TO_READ;
		}
		int windowSize;
		if (windowSizeIn > 0) {
			windowSize = windowSizeIn;
		} else {
			windowSize = COVERAGE_OF_LAST_EXONS_WINDOW_SIZE;
		}
		
		int numOfWindows = (int)((numOfBasesToRead + COVERAGE_OF_LAST_EXONS_EXTENSION_FACTOR) / windowSize);
		
		LinkedList<String> filesToRead = new LinkedList<String>();
		
		File coverageFolder = new File(PATH_TO_CURRENT_DIRECTORY);
		File[] filesInCoverageFolder = coverageFolder.listFiles();
		for (File fileInCoverageFolder : filesInCoverageFolder) {
			if (!fileInCoverageFolder.isFile())
				continue;
			
			String filenameWithPath = fileInCoverageFolder.getCanonicalPath();
			if (!filenameWithPath.endsWith(COVERAGE_OF_LAST_EXONS_JOB_COVERAGE_FILE_END))
				continue;
			
			filesToRead.add(filenameWithPath);
		}
		
		double[] totalCoverageAtWindow = new double[numOfWindows];
		for (int i = 0; i < numOfWindows; i++) {
			totalCoverageAtWindow[i] = 0;
		}
		
		double[] coveragePercentageTotalsAtWindow = new double[numOfWindows];
		for (int i = 0; i < numOfWindows; i++) {
			coveragePercentageTotalsAtWindow[i] = 0;
		}
		
		int numExons = 0;
		
		for (String file : filesToRead) {
			BufferedReader br = new BufferedReader(new FileReader(file));
			String line;
			while ((line = br.readLine()) != null && line.trim().length() > 0) {
				String[] lineSplit = line.split(TAB);
				double totalCoverage = 0;
				for (int i = 0; i < numOfWindows; i++) {
					double coverage = Double.parseDouble(lineSplit[i + 1]);
					totalCoverageAtWindow[i] += coverage;
					totalCoverage += coverage;
				}
				
				if (totalCoverage > 0) {
					for (int i = 0; i < numOfWindows; i++) {
						double coverage = Double.parseDouble(lineSplit[i + 1]);
						double coveragePercentage = coverage / totalCoverage;
						coveragePercentageTotalsAtWindow[i] += coveragePercentage;
					}
				}
				
				numExons++;
			}
			br.close();
		}

		FileWriter fw = new FileWriter(COVERAGE_OF_LAST_EXONS_SUM_FROM_ALL_JOBS_FILE);
		for (int i = 0; i < numOfWindows; i++) {
			double totalCoverage = totalCoverageAtWindow[i];
			double averagePercentage = coveragePercentageTotalsAtWindow[i] / (double)numExons;
			fw.write(String.valueOf(i) + TAB + String.valueOf(totalCoverage) + TAB + String.valueOf(averagePercentage) + NEW_LINE);
		}
		fw.close();
	}
	
	public static void findReadCoverageNearLastExonsAcrossGenomeBySplittingIntoJobs(String transcriptsFile, String chrSizesFile, String alignmentsFile, int numOfLastExonsPerJob, int numOfBasesToReadIn, int windowSizeIn) throws IOException, ParseException {
		int numOfBasesToRead;
		int minNumOfBases;
		if (numOfBasesToReadIn > 0) {
			numOfBasesToRead = numOfBasesToReadIn;
			minNumOfBases = numOfBasesToReadIn;
		} else {
			numOfBasesToRead = COVERAGE_OF_LAST_EXONS_NUM_OF_BASES_TO_READ;
			minNumOfBases = COVERAGE_OF_LAST_EXONS_MIN_NUM_OF_BASES;
		}
		int windowSize;
		if (windowSizeIn > 0) {
			windowSize = windowSizeIn;
		} else {
			windowSize = COVERAGE_OF_LAST_EXONS_WINDOW_SIZE;
		}
		
		AnnotationReader<? extends GenomicAnnotation> transcripts = AnnotationReaderFactory.create(transcriptsFile, BED_FORMAT);
		Map<GenomicAnnotation, Integer> lastExons = AnnotationReader.getLastExons(transcripts.getAnnotationList());
		
		StringBuilder nextJobExons = new StringBuilder();
		int numOfExonsSubmitted = 0;
		int countOnJob = 0;
		int jobNum = 0;
		for (GenomicAnnotation lastExon : lastExons.keySet()) {
			int startIndex = lastExon.getStart();
			int endIndex = lastExon.getEnd();
			
			if (endIndex - startIndex < minNumOfBases)
				continue;

			String chr = lastExon.getChromosome();
			int orientedEnd = lastExons.get(lastExon);
			
			int index1;
			int index2;
			if (startIndex == orientedEnd) {
				index1 = startIndex;
				index2 = endIndex;
			} else {
				index1 = endIndex;
				index2 = startIndex;
			}
			
			String exonInfo = chr + COVERAGE_OF_LAST_EXONS_SPLIT_INTO_JOBS_INFO_SEPARATOR + String.valueOf(index1) + COVERAGE_OF_LAST_EXONS_SPLIT_INTO_JOBS_INFO_SEPARATOR + String.valueOf(index2);
			nextJobExons.append(exonInfo);
			
			countOnJob++;

			if (countOnJob == numOfLastExonsPerJob) {
				jobNum++;
				
				System.out.println("Submitting job #: " + jobNum);
				String[] jobCmd = { "/bin/sh" , "-c" , "bsub -o " + PATH_TO_CURRENT_DIRECTORY + COVERAGE_OF_LAST_EXONS_JOB_FILE_START + String.valueOf(jobNum) + COVERAGE_OF_LAST_EXONS_JOB_BSUB_FILE_END + " -P polya -R \"rusage[mem=8]\" java -jar " + LOCATION_OF_SCRIPTURE_JAR + " -task polya -sizesFile " + chrSizesFile + " -transcripts " + transcriptsFile + " -alignments " + alignmentsFile + " -exons " + nextJobExons.toString() + " -jobnum " + String.valueOf(jobNum) + " -numbases " + String.valueOf(numOfBasesToRead) + " -windowsize " + String.valueOf(windowSize) };
				Process job = Runtime.getRuntime().exec(jobCmd);
				BufferedReader br = new BufferedReader(new InputStreamReader(job.getErrorStream()));
				String line;
				while ((line = br.readLine()) != null) {
					System.out.println("job" + String.valueOf(jobNum) + ": " + line);
				}
				br.close();

				numOfExonsSubmitted += countOnJob;
				
				nextJobExons = new StringBuilder();
				countOnJob = 0;
				
			} else {
				nextJobExons.append(COVERAGE_OF_LAST_EXONS_SPLIT_INTO_JOBS_JOB_SEPARATOR);
			}
		}
		
		if (countOnJob > 0) {
			jobNum++;
			
			System.out.println("Submitting job #: " + jobNum);
			String[] jobCmd = { "/bin/sh" , "-c" , "bsub -o " + PATH_TO_CURRENT_DIRECTORY + COVERAGE_OF_LAST_EXONS_JOB_FILE_START + String.valueOf(jobNum) + COVERAGE_OF_LAST_EXONS_JOB_BSUB_FILE_END + " -P polya -R \"rusage[mem=8]\" java -jar " + LOCATION_OF_SCRIPTURE_JAR + " -task polya -sizesFile " + chrSizesFile + " -transcripts " + transcriptsFile + " -alignments " + alignmentsFile + " -exons " + nextJobExons.toString() + " -jobnum " + String.valueOf(jobNum) + " -numbases " + String.valueOf(numOfBasesToRead) + " -windowsize " + String.valueOf(windowSize) };
			Process job = Runtime.getRuntime().exec(jobCmd);
			BufferedReader br = new BufferedReader(new InputStreamReader(job.getErrorStream()));
			String line;
			while ((line = br.readLine()) != null) {
				System.out.println("job" + String.valueOf(jobNum) + ": " + line);
			}
			br.close();
			
			numOfExonsSubmitted += countOnJob;
		}
		
		System.out.println("TOTAL EXONS SUBMITTED: " + numOfExonsSubmitted);
	}
	
	public static void findReadCoverageNearLastExonsByWindow(String alignmentsFile, String chrSizesFile, String lastExonsIn, int jobNum, int numOfBasesToRead, int windowSize) throws IOException {
		List<Object[]> lastExons = new LinkedList<Object[]>();
		String[] lastExonsInSplit = lastExonsIn.split(COVERAGE_OF_LAST_EXONS_SPLIT_INTO_JOBS_JOB_SEPARATOR);
		for (String lastExonInfo : lastExonsInSplit) {
			String[] lastExonInfoSplit = lastExonInfo.split(COVERAGE_OF_LAST_EXONS_SPLIT_INTO_JOBS_INFO_SEPARATOR);
			String chr = lastExonInfoSplit[0];
			int index1 = Integer.parseInt(lastExonInfoSplit[1]);
			int index2 = Integer.parseInt(lastExonInfoSplit[2]);
			Object[] o = new Object[] { chr, index1, index2 };
			lastExons.add(o);
		}
		
		GenericAlignmentDataModel alignments = new GenericAlignmentDataModel(alignmentsFile, chrSizesFile, false);
		
		int numOfWindows = (int)((numOfBasesToRead + COVERAGE_OF_LAST_EXONS_EXTENSION_FACTOR) / windowSize);
		
		Map<String, Double[]> exonsWithCoverage = new HashMap<String, Double[]>(lastExons.size());
		
		for (Object[] lastExon : lastExons) {
			String chr = (String)lastExon[0];
			int index1 = (Integer)lastExon[1];
			int index2 = (Integer)lastExon[2];
			
			boolean startIndexIsOrientedEnd;
			int startIndex;
			int endIndex;
			if (index1 < index2) {
				startIndex = index1;
				endIndex = index2;
				startIndexIsOrientedEnd = true;
			} else {
				startIndex = index2;
				endIndex = index1;
				startIndexIsOrientedEnd = false;
			}
			
			Double[] coverageAtWindows = new Double[numOfWindows];
			
			int rangeToCoverStart;
			int rangeToCoverEnd;
			if (startIndexIsOrientedEnd) {
				rangeToCoverStart = startIndex - COVERAGE_OF_LAST_EXONS_EXTENSION_FACTOR;
				rangeToCoverEnd = startIndex + numOfBasesToRead;
				int j = 0;
				for (int i = rangeToCoverEnd; i > rangeToCoverStart; i -= windowSize) {
					coverageAtWindows[j++] = alignments.getCountsOnWindow(chr, i - windowSize + 1, i + 1);
				}
			} else {
				rangeToCoverStart = endIndex - numOfBasesToRead;
				rangeToCoverEnd = endIndex + COVERAGE_OF_LAST_EXONS_EXTENSION_FACTOR;
				int j = 0;
				for (int i = rangeToCoverStart; i < rangeToCoverEnd; i += windowSize) {
					coverageAtWindows[j++] = alignments.getCountsOnWindow(chr, i, i + windowSize);
				}
			}
			
			int orientedEnd;
			if (startIndexIsOrientedEnd)
				orientedEnd = startIndex;
			else
				orientedEnd = endIndex;
			String coor = chr + COLON + String.valueOf(orientedEnd);
			
			exonsWithCoverage.put(coor, coverageAtWindows);
		}
		
		FileWriter fw = new FileWriter(PATH_TO_CURRENT_DIRECTORY + COVERAGE_OF_LAST_EXONS_JOB_FILE_START + String.valueOf(jobNum) + COVERAGE_OF_LAST_EXONS_JOB_COVERAGE_FILE_END);
		for (String exonEndCoor : exonsWithCoverage.keySet()) {
			Double[] coverageValues = exonsWithCoverage.get(exonEndCoor);
			
			StringBuilder line = new StringBuilder();
			line.append(exonEndCoor);
			line.append(TAB);
			
			for (int i = 0; i < numOfWindows; i++) {
				line.append(coverageValues[i]);
				if (i != numOfWindows - 1)
					line.append(TAB);
			}
			fw.write(line + NEW_LINE);
		}
		fw.close();
	}
	
	public static void pauseUntilAllAnnotationsAreWrittenAndConcatenate() throws IOException {
		int numOfAlignments = 0;
		BufferedReader brAlignments = new BufferedReader(new FileReader(PATH_TO_CURRENT_DIRECTORY + ANNOTATIONS_JOBS_DIRECTORY + WRITE_ANNOTATIONS_LIST_OF_ALIGNMENTS_FILE));
		String alignmentInfo;
		while ((alignmentInfo = brAlignments.readLine()) != null && alignmentInfo.trim().length() > 0) {
			numOfAlignments++;
		}
		
		int numOfJobs = (int)((double)numOfAlignments / (double)NUM_OF_ALIGNMENTS_PER_JOB_FOR_WRITING_AND_FILTERING_ANNOTATIONS);
		if (numOfAlignments % NUM_OF_ALIGNMENTS_PER_JOB_FOR_WRITING_AND_FILTERING_ANNOTATIONS > 0)
			numOfJobs++;
		
		Set<String> bsubFilesThatAreSuccessfullyWritten = new HashSet<String>();
		
		while (true) {
			File dir = new File(PATH_TO_CURRENT_DIRECTORY + ANNOTATIONS_JOBS_DIRECTORY);
			String[] filesInDirArray = dir.list();
			Set<String> filesInDir = new HashSet<String>();
			for (String fileInDir : filesInDirArray) {
				filesInDir.add(PATH_TO_CURRENT_DIRECTORY + ANNOTATIONS_JOBS_DIRECTORY + fileInDir);
			}
						
			boolean jobContainsError = false;
			boolean allAnnotationsAreSuccessfullyWritten = true;
			for (int i = 1; i <= numOfJobs; i++) {
				String bsubFileToCheck = PATH_TO_CURRENT_DIRECTORY + ANNOTATIONS_JOBS_DIRECTORY + WRITE_ANNOTATIONS_BSUB_JOB_FILE_START + String.valueOf(i) + WRITE_ANNOTATIONS_JOB_BSUB_FILE_END;
								
				if (bsubFilesThatAreSuccessfullyWritten.contains(bsubFileToCheck))
					continue;
				
				if (!filesInDir.contains(bsubFileToCheck)) {
					allAnnotationsAreSuccessfullyWritten = false;
					break;
				}
				
				try {
					Thread.sleep(TEN_MILLISECONDS);
				}
				catch (InterruptedException e) {
					e.printStackTrace();
					break;
				}
				
				BufferedReader br = new BufferedReader(new FileReader(bsubFileToCheck));
				String line;
				boolean lineContainsSuccessIndicator = false;
				while ((line = br.readLine()) != null) {
					if (line.trim().length() == 0)
						continue;
					
					if (line.indexOf(BSUB_FILE_SUCCESS_INDICATOR) >= 0) {
						lineContainsSuccessIndicator = true;
						break;
					}
				}
				
				if (!lineContainsSuccessIndicator) {
					jobContainsError = true;
					allAnnotationsAreSuccessfullyWritten = false;
					break;
				}
				
				bsubFilesThatAreSuccessfullyWritten.add(bsubFileToCheck);
			}
			
			if (jobContainsError) {
				System.out.println("There exists at least one job that encountered an error");
				System.exit(1);
			}
			
			if (allAnnotationsAreSuccessfullyWritten) {
				System.out.println("All annotations are successfully written");
				break;
			}
			
			try {
				Thread.sleep(ONE_MINUTE_IN_MILLISECS);
			}
			catch (InterruptedException e) {
				e.printStackTrace();
				break;
			}
		}
		
		System.out.println("Concatenating all alignments...");
		String[] jobCmdAll = { "/bin/sh" , "-c" , "cat " + PATH_TO_CURRENT_DIRECTORY + ANNOTATIONS_JOBS_DIRECTORY + WRITE_ANNOTATIONS_ALL_JOB_FILE_START + ASTERISK + WRITE_ANNOTATIONS_JOB_ANNOTATION_FILE_END + " > " + PATH_TO_CURRENT_DIRECTORY + ALL_ANNOTATIONS_OUT_FILE };
		Process jobAll = Runtime.getRuntime().exec(jobCmdAll);
		BufferedReader brAll = new BufferedReader(new InputStreamReader(jobAll.getErrorStream()));
		String lineAll;
		while ((lineAll = brAll.readLine()) != null) {
			System.out.println(jobCmdAll[2] + " : " + lineAll);
		}
		brAll.close();
		
		System.out.println("Concatenating filtered alignments...");
		String[] jobCmdFiltered = { "/bin/sh" , "-c" , "cat " + PATH_TO_CURRENT_DIRECTORY + ANNOTATIONS_JOBS_DIRECTORY + WRITE_ANNOTATIONS_FILTERED_JOB_FILE_START + ASTERISK + WRITE_ANNOTATIONS_JOB_ANNOTATION_FILE_END + " > " + PATH_TO_CURRENT_DIRECTORY + FILTERED_ANNOTATIONS_OUT_FILE };
		Process jobFiltered = Runtime.getRuntime().exec(jobCmdFiltered);
		BufferedReader brFiltered = new BufferedReader(new InputStreamReader(jobFiltered.getErrorStream()));
		String lineFiltered;
		while ((lineFiltered = brFiltered.readLine()) != null) {
			System.out.println(jobCmdFiltered[2] + " : " + lineFiltered);
		}
		brFiltered.close();
	}
	
	public static void writeAnnotationsOfReadsWithPolyAAndFilterAndWriteFilteredAnnotationsToBedFileForJob(String allAlignmentsInFile, String chrSizesFile, String alignmentsForJob, int jobNum) throws IOException {
		boolean shouldUsePolyAAlignmentFilteringMethod = false;
		boolean shouldUseSimpleRepeatsFilteringMethod = false;
		if (POLYA_IN_GENOME_FILTERING_METHOD.equals(POLYA_IN_GENOME_FILTERING_METHOD_POLYA))
			shouldUsePolyAAlignmentFilteringMethod = true;
		if (POLYA_IN_GENOME_FILTERING_METHOD.equals(POLYA_IN_GENOME_FILTERING_METHOD_SIMPLEREPEATS))
			shouldUseSimpleRepeatsFilteringMethod = true;
		if (POLYA_IN_GENOME_FILTERING_METHOD.equals(POLYA_IN_GENOME_FILTERING_METHOD_BOTH)) {
			shouldUsePolyAAlignmentFilteringMethod = true;
			shouldUseSimpleRepeatsFilteringMethod = true;
		}
		
		GenericAlignmentDataModel polyAAlignments = null;
		if (shouldUsePolyAAlignmentFilteringMethod)
			polyAAlignments = new GenericAlignmentDataModel(PATH_TO_MOST_FILES + POLYA_ALIGNMENT_SORTED_FILE, chrSizesFile, false);
		
		BEDFileParser simpleRepeatLocations = null;
		if (shouldUseSimpleRepeatsFilteringMethod)
			simpleRepeatLocations = new BEDFileParser(SIMPLE_REPEAT_LOCATIONS_FILTERED_FILE);
		
		GenericAlignmentDataModel allAlignments = new GenericAlignmentDataModel(allAlignmentsInFile, chrSizesFile, false);
		
		FileWriter fwAll = new FileWriter(PATH_TO_CURRENT_DIRECTORY + ANNOTATIONS_JOBS_DIRECTORY + WRITE_ANNOTATIONS_ALL_JOB_FILE_START + String.valueOf(jobNum) + WRITE_ANNOTATIONS_JOB_ANNOTATION_FILE_END);	
		FileWriter fwFiltered = new FileWriter(PATH_TO_CURRENT_DIRECTORY + ANNOTATIONS_JOBS_DIRECTORY + WRITE_ANNOTATIONS_FILTERED_JOB_FILE_START + String.valueOf(jobNum) + WRITE_ANNOTATIONS_JOB_ANNOTATION_FILE_END);
		
		int countAll = 0;
		int countFiltered = 0;
		
		String[] alignmentsForJobSplit = alignmentsForJob.split(WRITE_ANNOTATIONS_SPLIT_INTO_JOBS_JOB_SEPARATOR);
		for (String alignmentForJob : alignmentsForJobSplit) {
			String[] alignmentForJobSplit = alignmentForJob.split(WRITE_ANNOTATIONS_SPLIT_INTO_JOBS_INFO_SEPARATOR);
			String chr = alignmentForJobSplit[0];
			if (chr.equals(ASTERISK_LETTERS))
				chr = ASTERISK;
			int startIndex = Integer.parseInt(alignmentForJobSplit[1]);
			int endIndex = Integer.parseInt(alignmentForJobSplit[2]);
			String polyALoc = alignmentForJobSplit[3];
			boolean polyAIsOnRight = polyALoc.equals(R);
			
			Alignments a1 = new Alignments(chr, startIndex, endIndex);
			BED b = new BED(a1);
			
			String orientation = polyAIsOnRight ? PLUS : MINUS;
			b.setOrientation(orientation);
			
			fwAll.write(b.toString() + NEW_LINE);
			
			countAll++;
			if (countAll % 10000 == 0)
				fwAll.flush();
			
			if (startIndex < 0)
				continue;

			double numOfOriginalAlignmentsTouchingThisAlignmentInGenome = allAlignments.getCountsPerAlignmentWithSameEndpointForPolyA(a1, EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS_WITH_ORIGINAL_ALIGNMENT_OVERLAP, polyAIsOnRight, TYPICAL_ORIGINAL_ALIGNMENT_LENGTH, MAX_DISTANCE_FROM_ORIGINAL_ALIGNMENT_ENDPOINT_FOR_FILTERING_ANNOTATIONS);
			if (numOfOriginalAlignmentsTouchingThisAlignmentInGenome >= MIN_NUM_OF_ORIGINAL_ALIGNMENTS_FOR_FILTERING_ALIGNMENTS)
				continue;
			
			if (shouldUsePolyAAlignmentFilteringMethod) {
				double numOfPolyATouchingThisAlignmentInGenome = polyAAlignments.getCountsPerAlignmentForPolyA(a1, EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS_WITH_POLYA_OVERLAP);
				if (numOfPolyATouchingThisAlignmentInGenome >= MIN_NUM_OF_POLYA_ALIGNMENTS_FOR_FILTERING_ALIGNMENTS)
					continue;
			}
			
			if (shouldUseSimpleRepeatsFilteringMethod) {
				RefSeqGene a2 = new RefSeqGene(chr, startIndex - EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS_WITH_SIMPLE_REPEAT_LOCATIONS_OVERLAP, endIndex + EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS_WITH_SIMPLE_REPEAT_LOCATIONS_OVERLAP);
				int numOfSimpleRepeatsTouchingThisAlignmentInGenome = simpleRepeatLocations.getOverlappers(a2).size();
				if (numOfSimpleRepeatsTouchingThisAlignmentInGenome > 0)
					continue;
			}

			fwFiltered.write(b.toString() + NEW_LINE);
			
			countFiltered++;
			if (countFiltered % 10000 == 0)
				fwFiltered.flush();
		}
		
		fwAll.close();
		fwFiltered.close();
	}
	
	private static boolean polyAForAlignmentIsOnRight(Alignment a, boolean startsAlignment) {
		boolean isNegativeStrand = a.isNegativeStrand(); // is true when 16 (reversed) ; b is false when 0 (correct)
		boolean polyAIsOnRight;
		if (startsAlignment) {
			if (isNegativeStrand)
				polyAIsOnRight = true;
			else
				polyAIsOnRight = false;
		} else {
			if (isNegativeStrand)
				polyAIsOnRight = false;
			else
				polyAIsOnRight = true;
		}
		return polyAIsOnRight;
	}
	
	public static void writeAnnotationsOfReadsWithPolyAAndFilterAndWriteFilteredAnnotationsToBedFileBySplittingIntoJobs(String allAlignmentsInFile, String chrSizesFile, int numOfAlignmentsPerJob, String parametersFile) throws IOException {
		boolean annotationsJobsDirectoryExists = (new File(ANNOTATIONS_JOBS_DIRECTORY)).exists();
		if (!annotationsJobsDirectoryExists) {
			System.out.println("Making annotations jobs directory");
			boolean madeAnnotationsJobsDirectory = (new File(ANNOTATIONS_JOBS_DIRECTORY)).mkdir();
			System.out.println("Made annotations jobs directory: " + madeAnnotationsJobsDirectory);
		}
		
		boolean shouldUsePolyAAlignmentFilteringMethod = false;
		boolean shouldUseSimpleRepeatsFilteringMethod = false;
		if (POLYA_IN_GENOME_FILTERING_METHOD.equals(POLYA_IN_GENOME_FILTERING_METHOD_POLYA))
			shouldUsePolyAAlignmentFilteringMethod = true;
		if (POLYA_IN_GENOME_FILTERING_METHOD.equals(POLYA_IN_GENOME_FILTERING_METHOD_SIMPLEREPEATS))
			shouldUseSimpleRepeatsFilteringMethod = true;
		if (POLYA_IN_GENOME_FILTERING_METHOD.equals(POLYA_IN_GENOME_FILTERING_METHOD_BOTH)) {
			shouldUsePolyAAlignmentFilteringMethod = true;
			shouldUseSimpleRepeatsFilteringMethod = true;
		}
		
		if (shouldUseSimpleRepeatsFilteringMethod) {
			System.out.println("Filtering simple repeats file...");
			filterSimpleRepeatLocationsFile(SIMPLE_REPEAT_LOCATIONS_FILE, SIMPLE_REPEAT_LOCATIONS_FILTERED_FILE);
			System.out.println("Done filtering simple repeats file...");
		}
		
		FileWriter fw = new FileWriter(PATH_TO_CURRENT_DIRECTORY + ANNOTATIONS_JOBS_DIRECTORY + WRITE_ANNOTATIONS_LIST_OF_ALIGNMENTS_FILE);
		int count = 0;
		
		AlignmentQueryReader readerStartA = new BAMQueryReader(new File(PATH_TO_MOST_FILES + READS_WITH_POLYA_START_WITH_A_ALIGNMENT_FILE));  //SamQueryReaderFactory.getReader(new ResourceLocator(PATH_TO_MOST_FILES + READS_WITH_POLYA_START_WITH_A_ALIGNMENT_FILE), false);
		CloseableIterator<Alignment> readerIterStartA = readerStartA.iterator();
		while (readerIterStartA.hasNext()) {
			Alignment a = readerIterStartA.next();
			String chr = a.getChromosome();
			if (chr.equals(ASTERISK))
				chr = ASTERISK_LETTERS;
			int startIndex = a.getAlignmentStart();
			int endIndex = a.getAlignmentEnd();
			String polyALoc = polyAForAlignmentIsOnRight(a, true) ? R : L;
			String alignmentInfo = chr + WRITE_ANNOTATIONS_SPLIT_INTO_JOBS_INFO_SEPARATOR + String.valueOf(startIndex) + WRITE_ANNOTATIONS_SPLIT_INTO_JOBS_INFO_SEPARATOR + String.valueOf(endIndex) + WRITE_ANNOTATIONS_SPLIT_INTO_JOBS_INFO_SEPARATOR + polyALoc;
			fw.write(alignmentInfo + NEW_LINE);
			count++;
			if (count % 10000 == 0)
				fw.flush();
		}
		
		AlignmentQueryReader readerStartT =  new BAMQueryReader(new File(PATH_TO_MOST_FILES + READS_WITH_POLYA_START_WITH_T_ALIGNMENT_FILE));//SamQueryReaderFactory.getReader(new ResourceLocator(PATH_TO_MOST_FILES + READS_WITH_POLYA_START_WITH_T_ALIGNMENT_FILE), false);
		CloseableIterator<Alignment> readerIterStartT = readerStartT.iterator();
		while (readerIterStartT.hasNext()) {
			Alignment a = readerIterStartT.next();
			String chr = a.getChromosome();
			if (chr.equals(ASTERISK))
				chr = ASTERISK_LETTERS;
			int startIndex = a.getAlignmentStart();
			int endIndex = a.getAlignmentEnd();
			String polyALoc = polyAForAlignmentIsOnRight(a, true) ? R : L;
			String alignmentInfo = chr + WRITE_ANNOTATIONS_SPLIT_INTO_JOBS_INFO_SEPARATOR + String.valueOf(startIndex) + WRITE_ANNOTATIONS_SPLIT_INTO_JOBS_INFO_SEPARATOR + String.valueOf(endIndex) + WRITE_ANNOTATIONS_SPLIT_INTO_JOBS_INFO_SEPARATOR + polyALoc;
			fw.write(alignmentInfo + NEW_LINE);
			count++;
			if (count % 10000 == 0)
				fw.flush();
		}
		
		AlignmentQueryReader readerEndA = new BAMQueryReader(new File(PATH_TO_MOST_FILES + READS_WITH_POLYA_END_WITH_A_ALIGNMENT_FILE));//SamQueryReaderFactory.getReader(new ResourceLocator(PATH_TO_MOST_FILES + READS_WITH_POLYA_END_WITH_A_ALIGNMENT_FILE), false);
		CloseableIterator<Alignment> readerIterEndA = readerEndA.iterator();
		while (readerIterEndA.hasNext()) {
			Alignment a = readerIterEndA.next();
			String chr = a.getChromosome();
			if (chr.equals(ASTERISK))
				chr = ASTERISK_LETTERS;
			int startIndex = a.getAlignmentStart();
			int endIndex = a.getAlignmentEnd();
			String polyALoc = polyAForAlignmentIsOnRight(a, false) ? R : L;
			String alignmentInfo = chr + WRITE_ANNOTATIONS_SPLIT_INTO_JOBS_INFO_SEPARATOR + String.valueOf(startIndex) + WRITE_ANNOTATIONS_SPLIT_INTO_JOBS_INFO_SEPARATOR + String.valueOf(endIndex) + WRITE_ANNOTATIONS_SPLIT_INTO_JOBS_INFO_SEPARATOR + polyALoc;
			fw.write(alignmentInfo + NEW_LINE);
			count++;
			if (count % 10000 == 0)
				fw.flush();
		}
		
		AlignmentQueryReader readerEndT = new BAMQueryReader(new File(PATH_TO_MOST_FILES + READS_WITH_POLYA_END_WITH_T_ALIGNMENT_FILE));//SamQueryReaderFactory.getReader(new ResourceLocator(PATH_TO_MOST_FILES + READS_WITH_POLYA_END_WITH_T_ALIGNMENT_FILE), false);
		CloseableIterator<Alignment> readerIterEndT = readerEndT.iterator();
		while (readerIterEndT.hasNext()) {
			Alignment a = readerIterEndT.next();
			String chr = a.getChromosome();
			if (chr.equals(ASTERISK))
				chr = ASTERISK_LETTERS;
			int startIndex = a.getAlignmentStart();
			int endIndex = a.getAlignmentEnd();
			String polyALoc = polyAForAlignmentIsOnRight(a, false) ? R : L;
			String alignmentInfo = chr + WRITE_ANNOTATIONS_SPLIT_INTO_JOBS_INFO_SEPARATOR + String.valueOf(startIndex) + WRITE_ANNOTATIONS_SPLIT_INTO_JOBS_INFO_SEPARATOR + String.valueOf(endIndex) + WRITE_ANNOTATIONS_SPLIT_INTO_JOBS_INFO_SEPARATOR + polyALoc;
			fw.write(alignmentInfo + NEW_LINE);
			count++;
			if (count % 10000 == 0)
				fw.flush();
		}
		
		fw.close();
		
		StringBuilder nextJobAlignments = new StringBuilder();
		int numOfAlignmentsSubmitted = 0;
		int countOnJob = 0;
		int jobNum = 0;
		
		String parameters = "";
		if (!parametersFile.equals(""))
			parameters = "-parameters " + PATH_TO_MOST_FILES + parametersFile + " ";
		
		BufferedReader brAlignments = new BufferedReader(new FileReader(PATH_TO_CURRENT_DIRECTORY + ANNOTATIONS_JOBS_DIRECTORY + WRITE_ANNOTATIONS_LIST_OF_ALIGNMENTS_FILE));
		String alignmentInfo;
		while ((alignmentInfo = brAlignments.readLine()) != null && alignmentInfo.trim().length() > 0) {
			nextJobAlignments.append(alignmentInfo);
			
			countOnJob++;

			if (countOnJob == numOfAlignmentsPerJob) {
				jobNum++;
				
				System.out.println("Submitting job #: " + jobNum);
				String[] jobCmd = { "/bin/sh" , "-c" , "bsub -o " + PATH_TO_CURRENT_DIRECTORY + ANNOTATIONS_JOBS_DIRECTORY + WRITE_ANNOTATIONS_BSUB_JOB_FILE_START + String.valueOf(jobNum) + WRITE_ANNOTATIONS_JOB_BSUB_FILE_END + " -P polya -R \"rusage[mem=5]\" java -Xmx4000m -jar " + LOCATION_OF_SCRIPTURE_JAR + " -task polya " + parameters + "-polyasubtask writeandfilterannotationsofreadswithpolya -alignmentsForJob " + nextJobAlignments.toString() + " -jobnum " + String.valueOf(jobNum) + " -pathtomostfiles " + PATH_TO_MOST_FILES };
				Process job = Runtime.getRuntime().exec(jobCmd);
				BufferedReader br = new BufferedReader(new InputStreamReader(job.getErrorStream()));
				String line;
				while ((line = br.readLine()) != null) {
					System.out.println("job" + String.valueOf(jobNum) + ": " + line);
				}
				br.close();

				numOfAlignmentsSubmitted += countOnJob;
				
				nextJobAlignments = new StringBuilder();
				countOnJob = 0;
			} else {
				nextJobAlignments.append(WRITE_ANNOTATIONS_SPLIT_INTO_JOBS_JOB_SEPARATOR);
			}
		}
		brAlignments.close();
		
		if (countOnJob > 0) {
			jobNum++;
			
			System.out.println("Submitting job #: " + jobNum);
			String[] jobCmd = { "/bin/sh" , "-c" , "bsub -o " + PATH_TO_CURRENT_DIRECTORY + ANNOTATIONS_JOBS_DIRECTORY + WRITE_ANNOTATIONS_BSUB_JOB_FILE_START + String.valueOf(jobNum) + WRITE_ANNOTATIONS_JOB_BSUB_FILE_END + " -P polya -R \"rusage[mem=5]\" java -Xmx4000m -jar " + LOCATION_OF_SCRIPTURE_JAR + " -task polya " + parameters + "-polyasubtask writeandfilterannotationsofreadswithpolya -alignmentsForJob " + nextJobAlignments.toString() + " -jobnum " + String.valueOf(jobNum) + " -pathtomostfiles " + PATH_TO_MOST_FILES };
			Process job = Runtime.getRuntime().exec(jobCmd);
			BufferedReader br = new BufferedReader(new InputStreamReader(job.getErrorStream()));
			String line;
			while ((line = br.readLine()) != null) {
				System.out.println("job" + String.valueOf(jobNum) + ": " + line);
			}
			br.close();
			
			numOfAlignmentsSubmitted += countOnJob;
		}
		
		System.out.println("TOTAL ALIGNMENTS SUBMITTED: " + numOfAlignmentsSubmitted);
	}
	
	/**
	private static void writeAnnotationsOfReadsWithPolyAAndFilterAndWriteFilteredAnnotationsToBedFile(String allAnnotationsOutFile, String filteredAnnotationsOutFile, String allAlignmentsInFile, String chrSizesFile) throws IOException {
		GenericAlignmentDataModel polyAAlignments = new GenericAlignmentDataModel(POLYA_ALIGNMENT_SORTED_FILE, chrSizesFile, false);
		
		GenericAlignmentDataModel allAlignments = new GenericAlignmentDataModel(allAlignmentsInFile, chrSizesFile, false);
		
		FileWriter fwAll = new FileWriter(allAnnotationsOutFile);	
		FileWriter fwFiltered = new FileWriter(filteredAnnotationsOutFile);
		
		int countAll = 0;
		int countFiltered = 0;
		
		AlignmentQueryReader readerStartA = SamQueryReaderFactory.getReader(new ResourceLocator(READS_WITH_POLYA_START_WITH_A_ALIGNMENT_FILE), false);
		CloseableIterator<Alignment> readerIterStartA = readerStartA.iterator();
		while (readerIterStartA.hasNext()) {
			Alignment a = readerIterStartA.next();
			Alignments a1 = new Alignments(a.getChromosome(), a.getAlignmentStart(), a.getAlignmentEnd());
			BED b = new BED(a1);
			
			fwAll.write(b.toString() + NEW_LINE);
			
			countAll++;
			if (countAll % 10000 == 0)
				fwAll.flush();
			
			if (a.getStart() < 0)
				continue;
			
			double numTouchingPolyAInGenome = polyAAlignments.getCountsPerAlignmentForPolyA(a1, EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS);
			if (numTouchingPolyAInGenome > 0)
				continue;
			
			double numTouchingOriginalAlignments = allAlignments.getCountsPerAlignmentWithSameEndpointsForPolyA(a1, EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS);
			if (numTouchingOriginalAlignments > 0)
				continue;
			
			fwFiltered.write(b.toString() + NEW_LINE);
			
			countFiltered++;
			if (countFiltered % 10000 == 0)
				fwFiltered.flush();
		}
		readerIterStartA.close();
		
		AlignmentQueryReader readerStartT = SamQueryReaderFactory.getReader(new ResourceLocator(READS_WITH_POLYA_START_WITH_T_ALIGNMENT_FILE), false);
		CloseableIterator<Alignment> readerIterStartT = readerStartT.iterator();
		while (readerIterStartT.hasNext()) {
			Alignment a = readerIterStartT.next();
			Alignments a1 = new Alignments(a.getChromosome(), a.getAlignmentStart(), a.getAlignmentEnd());
			BED b = new BED(a1);
			
			fwAll.write(b.toString() + NEW_LINE);
			
			countAll++;
			if (countAll % 10000 == 0)
				fwAll.flush();
			
			if (a.getStart() < 0)
				continue;
			
			double numTouchingPolyAInGenome = polyAAlignments.getCountsPerAlignmentForPolyA(a1, EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS);
			if (numTouchingPolyAInGenome > 0)
				continue;
			
			double numTouchingOriginalAlignments = allAlignments.getCountsPerAlignmentWithSameEndpointsForPolyA(a1, EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS);
			if (numTouchingOriginalAlignments > 0)
				continue;
			
			fwFiltered.write(b.toString() + NEW_LINE);
			
			countFiltered++;
			if (countFiltered % 10000 == 0)
				fwFiltered.flush();
		}
		readerStartT.close();
		
		AlignmentQueryReader readerEndA = SamQueryReaderFactory.getReader(new ResourceLocator(READS_WITH_POLYA_END_WITH_A_ALIGNMENT_FILE), false);
		CloseableIterator<Alignment> readerIterEndA = readerEndA.iterator();
		while (readerIterEndA.hasNext()) {
			Alignment a = readerIterEndA.next();
			Alignments a1 = new Alignments(a.getChromosome(), a.getAlignmentStart(), a.getAlignmentEnd());
			BED b = new BED(a1);
			
			fwAll.write(b.toString() + NEW_LINE);
			
			countAll++;
			if (countAll % 10000 == 0)
				fwAll.flush();
			
			if (a.getStart() < 0)
				continue;
			
			double numTouchingPolyAInGenome = polyAAlignments.getCountsPerAlignmentForPolyA(a1, EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS);
			if (numTouchingPolyAInGenome > 0)
				continue;
			
			double numTouchingOriginalAlignments = allAlignments.getCountsPerAlignmentWithSameEndpointsForPolyA(a1, EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS);
			if (numTouchingOriginalAlignments > 0)
				continue;
			
			fwFiltered.write(b.toString() + NEW_LINE);
			
			countFiltered++;
			if (countFiltered % 10000 == 0)
				fwFiltered.flush();
		}
		readerEndA.close();
		
		AlignmentQueryReader readerEndT = SamQueryReaderFactory.getReader(new ResourceLocator(READS_WITH_POLYA_END_WITH_T_ALIGNMENT_FILE), false);
		CloseableIterator<Alignment> readerIterEndT = readerEndT.iterator();
		while (readerIterEndT.hasNext()) {
			Alignment a = readerIterEndT.next();
			Alignments a1 = new Alignments(a.getChromosome(), a.getAlignmentStart(), a.getAlignmentEnd());
			BED b = new BED(a1);
			
			fwAll.write(b.toString() + NEW_LINE);
			
			countAll++;
			if (countAll % 10000 == 0)
				fwAll.flush();

			if (a.getStart() < 0)
				continue;

			double numTouchingPolyAInGenome = polyAAlignments.getCountsPerAlignmentForPolyA(a1, EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS);
			if (numTouchingPolyAInGenome > 0)
				continue;

			double numTouchingOriginalAlignments = allAlignments.getCountsPerAlignmentWithSameEndpointsForPolyA(a1, EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS);
			if (numTouchingOriginalAlignments > 0)
				continue;

			fwFiltered.write(b.toString() + NEW_LINE);
			
			countFiltered++;
			if (countFiltered % 10000 == 0)
				fwFiltered.flush();
		}
		readerEndT.close();
		
		fwAll.close();
		fwFiltered.close();
	}
	*/
	
	private static void writeAlignmentsOfMatesOfFullPolyAReadsToBedFile(String bedOutFile) throws IOException {
		FileWriter fw = new FileWriter(bedOutFile);
		
		int count = 0;
		
		AlignmentQueryReader readerA = new BAMQueryReader(new File(MATES_OF_FULL_POLYA_READS_A_ALIGNMENT_FILE));//  SamQueryReaderFactory.getReader(new ResourceLocator(MATES_OF_FULL_POLYA_READS_A_ALIGNMENT_FILE), false);
		CloseableIterator<Alignment> readerIterA = readerA.iterator();
		while (readerIterA.hasNext()) {
			Alignment a = readerIterA.next();
			Alignments a1 = new Alignments(a.getChromosome(), a.getAlignmentStart(), a.getAlignmentEnd());
			BED b = new BED(a1);
			fw.write(b.toString() + NEW_LINE);
			
			count++;
			if (count % 10000 == 0)
				fw.flush();
		}
		readerIterA.close();
		
		AlignmentQueryReader readerT =  new BAMQueryReader(new File( MATES_OF_FULL_POLYA_READS_T_ALIGNMENT_FILE)); ;//SamQueryReaderFactory.getReader(new ResourceLocator(MATES_OF_FULL_POLYA_READS_T_ALIGNMENT_FILE), false);
		CloseableIterator<Alignment> readerIterT = readerT.iterator();
		while (readerIterT.hasNext()) {
			Alignment a = readerIterT.next();
			Alignments a1 = new Alignments(a.getChromosome(), a.getAlignmentStart(), a.getAlignmentEnd());
			BED b = new BED(a1);
			fw.write(b.toString() + NEW_LINE);
			
			count++;
			if (count % 10000 == 0)
				fw.flush();
		}
		readerIterT.close();
		
		fw.close();
	}
	
	private static void filterSetOfAlignmentsByOverlapWithTranscripts(String filteredAlignmentsOutFile, String allAlignmentsOutFile, String scriptureTranscriptsFile) throws IOException, ParseException {
		PolyAUtils.filterAlignmentsByOverlapWithTranscripts(filteredAlignmentsOutFile, scriptureTranscriptsFile, "filteredAlignmentsByOverlapWithTranscripts_ReadsWithPolyA.bed", false, false);
		System.out.println("Wrote filtered alignments by overlap with transcripts");
		PolyAUtils.filterAlignmentsByOverlapWithTranscripts(filteredAlignmentsOutFile, scriptureTranscriptsFile, "filteredAlignmentsByOverlapWithTranscriptsWithExtension_ReadsWithPolyA.bed", true, false);
		System.out.println("Wrote filtered alignments by overlap with transcripts with extension");
		
		PolyAUtils.filterAlignmentsByOverlapWithTranscripts("matesOfFullPolyAReads.bed", scriptureTranscriptsFile, "filteredAlignmentsByOverlapWithTranscripts_MatesOfFullPolyAReads.bed", false, false);
		System.out.println("Wrote filtered alignments by overlap with transcripts - MATES");
		PolyAUtils.filterAlignmentsByOverlapWithTranscripts("matesOfFullPolyAReads.bed", scriptureTranscriptsFile, "filteredAlignmentsByOverlapWithTranscriptsWithExtension_MatesOfFullPolyAReads.bed", true, false);
		System.out.println("Wrote filtered alignments by overlap with transcripts with extension - MATES");
		
		PolyAUtils.filterAlignmentsByOverlapWithTranscripts(filteredAlignmentsOutFile, scriptureTranscriptsFile, "filteredAlignmentsByOverlapWithTranscriptsLastExon_ReadsWithPolyA.bed", false, true);
		System.out.println("Wrote filtered alignments by overlap with transcripts using only the last exon");
		PolyAUtils.filterAlignmentsByOverlapWithTranscripts(filteredAlignmentsOutFile, scriptureTranscriptsFile, "filteredAlignmentsByOverlapWithTranscriptsWithExtensionLastExon_ReadsWithPolyA.bed", true, true);
		System.out.println("Wrote filtered alignments) by overlap with transcripts using only the last exon with extension");
		
		PolyAUtils.filterAlignmentsByOverlapWithTranscripts("matesOfFullPolyAReads.bed", scriptureTranscriptsFile, "filteredAlignmentsByOverlapWithTranscriptsLastExon_MatesOfFullPolyAReads.bed", false, true);
		System.out.println("Wrote filtered alignments by overlap with transcripts using only the last exon - MATES");
		PolyAUtils.filterAlignmentsByOverlapWithTranscripts("matesOfFullPolyAReads.bed", scriptureTranscriptsFile, "filteredAlignmentsByOverlapWithTranscriptsWithExtensionLastExon_MatesOfFullPolyAReads.bed", true, true);
		System.out.println("Wrote filtered alignments by overlap with transcripts using only the last exon with extension - MATES");
		
		PolyAUtils.filterAlignmentsByOverlapWithTranscripts(allAlignmentsOutFile, scriptureTranscriptsFile, "filteredAlignmentsByOverlapWithTranscripts_ALL_ReadsWithPolyA.bed", false, false);
		System.out.println("Wrote filtered (ALL) alignments by overlap with transcripts");
		PolyAUtils.filterAlignmentsByOverlapWithTranscripts(allAlignmentsOutFile, scriptureTranscriptsFile, "filteredAlignmentsByOverlapWithTranscriptsWithExtension_ALL_ReadsWithPolyA.bed", true, false);
		System.out.println("Wrote filtered (ALL) alignments by overlap with transcripts with extension");
		PolyAUtils.filterAlignmentsByOverlapWithTranscripts(allAlignmentsOutFile, scriptureTranscriptsFile, "filteredAlignmentsByOverlapWithTranscriptsLastExon_ALL_ReadsWithPolyA.bed", false, true);
		System.out.println("Wrote filtered (ALL) alignments by overlap with transcripts using only last exon");
		PolyAUtils.filterAlignmentsByOverlapWithTranscripts(allAlignmentsOutFile, scriptureTranscriptsFile, "filteredAlignmentsByOverlapWithTranscriptsWithExtensionLastExon_ALL_ReadsWithPolyA.bed", true, true);
		System.out.println("Wrote filtered (ALL) alignments by overlap with transcripts using only last exon with extension");
	}
	
	//private static int findCoordinatesOfEndsOfTranscriptsHandleOverlappingAlignment(Alignment readWithPolyAAlignment, String chr, RefSeqGene transcriptAnnot, int start, int end, int overlapStart, int overlapEnd, boolean transcriptAnnotEndIsOnRight, Map<String, Set<RefSeqGene>> coordinatesOfTranscriptEndsWithTranscript, boolean startsAlignment) {
	private static void findCoordinatesOfEndsOfTranscriptsHandleOverlappingAlignment(Alignment readWithPolyAAlignment, String chr, RefSeqGene transcriptAnnot, int start, int end, int overlapStart, int overlapEnd, boolean transcriptAnnotEndIsOnRight, Set<Integer> possibleCoordinatesOfTranscriptEnd, boolean startsAlignment) {
		if (!(readWithPolyAAlignment.getStart() == start && readWithPolyAAlignment.getEnd() == end))
			return;
		
		boolean isNegativeStrand = readWithPolyAAlignment.isNegativeStrand(); // is true when 16 (reversed) ; b is false when 0 (correct)
		boolean polyAIsOnRight;
		if (startsAlignment) {
			if (isNegativeStrand)
				polyAIsOnRight = true;
			else
				polyAIsOnRight = false;
		} else {
			if (isNegativeStrand)
				polyAIsOnRight = false;
			else
				polyAIsOnRight = true;
		}
		
		int coor;
		if (polyAIsOnRight)
			coor = end;
		else
			coor = start;
		
		if (!(coor >= overlapStart && coor <= overlapEnd))
			return;
		
		// ++++++++++
		/**
		if (polyAIsOnRight == transcriptAnnotEndIsOnRight)
			System.out.println("YES ::: " + polyAIsOnRight + " ::: " + transcriptAnnot.getName() + " ::: " + chr + ":" + String.valueOf(coor));
		else
			System.out.println("NO ::: " + polyAIsOnRight + " ::: " + transcriptAnnot.getName() + " ::: " + chr + ":" + String.valueOf(coor));
		*/
		// ++++++++++
				
		if (polyAIsOnRight != transcriptAnnotEndIsOnRight)
			return;
		
		/**
		String completeCoor = chr + COLON + String.valueOf(coor);
		
		if (coordinatesOfTranscriptEndsWithTranscript.containsKey(completeCoor)) {
			coordinatesOfTranscriptEndsWithTranscript.get(completeCoor).add(transcriptAnnot);
		} else {
			Set<RefSeqGene> s = new HashSet<RefSeqGene>();
			s.add(transcriptAnnot);
			coordinatesOfTranscriptEndsWithTranscript.put(completeCoor, s);
		}
		*/
		
		possibleCoordinatesOfTranscriptEnd.add(coor);
	}
	
	private static int[] findOverlapPointsOfTranscriptsForGivenReadWithPolyA(GenomicAnnotation read) {
		int start = read.getStart();
		int end = read.getEnd();
		String orientation = read.getOrientation();
		
		int overlapStart;
		int overlapEnd;
		

		if (orientation.equals(PLUS)) {
			overlapStart = start - EXTENSION_FACTOR_AFTER_LAST_EXON_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_LAST_EXON;
			overlapEnd = end;
		} else if (orientation.equals(MINUS)) {
			overlapStart = start;
			overlapEnd = end + EXTENSION_FACTOR_AFTER_LAST_EXON_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_LAST_EXON;
		} else {
			overlapStart = start - EXTENSION_FACTOR_AFTER_LAST_EXON_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_LAST_EXON;
			overlapEnd = end + EXTENSION_FACTOR_AFTER_LAST_EXON_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_LAST_EXON;
		}
		
		int[] overlapPoints = new int[] { overlapStart, overlapEnd };
		return overlapPoints;
	}
	
	/**
	private static int[] findOverlapPointsOfTranscriptsForGivenReadWithPolyA(GenomicAnnotation read, GenericAlignmentDataModel alignmentsStartA, GenericAlignmentDataModel alignmentsStartT, GenericAlignmentDataModel alignmentsEndA, GenericAlignmentDataModel alignmentsEndT) throws ParseException {
		int start = read.getStart();
		int end = read.getEnd();
		String chr = read.getChromosome();
		Alignments a1 = new Alignments(chr, start, end);

		boolean differingLocOfPolyA = false;
		boolean polyALocIsOnRight = false;
		boolean polyALocHasBeenSet = false;
		
		CloseableIterator<Alignment> overlappingAlignmentsStartA = alignmentsStartA.getAlignmentsOverlappingRegion(a1);
		while(overlappingAlignmentsStartA.hasNext()) {
			Alignment readWithPolyAAlignment = overlappingAlignmentsStartA.next();

			boolean startsAlignment = true;
			
			if (!(readWithPolyAAlignment.getStart() == start && readWithPolyAAlignment.getEnd() == end))
				continue;
			
			boolean isNegativeStrand = readWithPolyAAlignment.isNegativeStrand(); // is true when 16 (reversed) ; b is false when 0 (correct)
			boolean polyAIsOnRight;
			if (startsAlignment) {
				if (isNegativeStrand)
					polyAIsOnRight = true;
				else
					polyAIsOnRight = false;
			} else {
				if (isNegativeStrand)
					polyAIsOnRight = false;
				else
					polyAIsOnRight = true;
			}
			
			if (!polyALocHasBeenSet) {
				polyALocIsOnRight = polyAIsOnRight;
				
				polyALocHasBeenSet = true;
			} else {
				if (polyALocIsOnRight != polyAIsOnRight)
					differingLocOfPolyA = true;
			}
		}
		overlappingAlignmentsStartA.close();
		
		CloseableIterator<Alignment> overlappingAlignmentsStartT = alignmentsStartT.getAlignmentsOverlappingRegion(a1);
		while(overlappingAlignmentsStartT.hasNext()) {
			Alignment readWithPolyAAlignment = overlappingAlignmentsStartT.next();

			boolean startsAlignment = true;
			
			if (!(readWithPolyAAlignment.getStart() == start && readWithPolyAAlignment.getEnd() == end))
				continue;
			
			boolean isNegativeStrand = readWithPolyAAlignment.isNegativeStrand(); // is true when 16 (reversed) ; b is false when 0 (correct)
			boolean polyAIsOnRight;
			if (startsAlignment) {
				if (isNegativeStrand)
					polyAIsOnRight = true;
				else
					polyAIsOnRight = false;
			} else {
				if (isNegativeStrand)
					polyAIsOnRight = false;
				else
					polyAIsOnRight = true;
			}
			
			if (!polyALocHasBeenSet) {
				polyALocIsOnRight = polyAIsOnRight;
				
				polyALocHasBeenSet = true;
			} else {
				if (polyALocIsOnRight != polyAIsOnRight)
					differingLocOfPolyA = true;
			}
		}
		overlappingAlignmentsStartT.close();
		
		CloseableIterator<Alignment> overlappingAlignmentsEndA = alignmentsEndA.getAlignmentsOverlappingRegion(a1);
		while(overlappingAlignmentsEndA.hasNext()) {
			Alignment readWithPolyAAlignment = overlappingAlignmentsEndA.next();

			boolean startsAlignment = false;
			
			if (!(readWithPolyAAlignment.getStart() == start && readWithPolyAAlignment.getEnd() == end))
				continue;
			
			boolean isNegativeStrand = readWithPolyAAlignment.isNegativeStrand(); // is true when 16 (reversed) ; b is false when 0 (correct)
			boolean polyAIsOnRight;
			if (startsAlignment) {
				if (isNegativeStrand)
					polyAIsOnRight = true;
				else
					polyAIsOnRight = false;
			} else {
				if (isNegativeStrand)
					polyAIsOnRight = false;
				else
					polyAIsOnRight = true;
			}
			
			if (!polyALocHasBeenSet) {
				polyALocIsOnRight = polyAIsOnRight;
				
				polyALocHasBeenSet = true;
			} else {
				if (polyALocIsOnRight != polyAIsOnRight)
					differingLocOfPolyA = true;
			}
		}
		overlappingAlignmentsEndA.close();
		
		CloseableIterator<Alignment> overlappingAlignmentsEndT = alignmentsEndT.getAlignmentsOverlappingRegion(a1);
		while(overlappingAlignmentsEndT.hasNext()) {
			Alignment readWithPolyAAlignment = overlappingAlignmentsEndT.next();
			
			boolean startsAlignment = false;
			
			if (!(readWithPolyAAlignment.getStart() == start && readWithPolyAAlignment.getEnd() == end))
				continue;
			
			boolean isNegativeStrand = readWithPolyAAlignment.isNegativeStrand(); // is true when 16 (reversed) ; b is false when 0 (correct)
			boolean polyAIsOnRight;
			if (startsAlignment) {
				if (isNegativeStrand)
					polyAIsOnRight = true;
				else
					polyAIsOnRight = false;
			} else {
				if (isNegativeStrand)
					polyAIsOnRight = false;
				else
					polyAIsOnRight = true;
			}
			
			if (!polyALocHasBeenSet) {
				polyALocIsOnRight = polyAIsOnRight;
				
				polyALocHasBeenSet = true;
			} else {
				if (polyALocIsOnRight != polyAIsOnRight)
					differingLocOfPolyA = true;
			}
		}
		overlappingAlignmentsEndT.close();
		
		int overlapStart;
		int overlapEnd;
		
		if (differingLocOfPolyA) {
			overlapStart = start - EXTENSION_FACTOR_AFTER_LAST_EXON_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_LAST_EXON;
			overlapEnd = end + EXTENSION_FACTOR_AFTER_LAST_EXON_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_LAST_EXON;
		} else {
			if (polyALocIsOnRight) {
				overlapStart = start - EXTENSION_FACTOR_AFTER_LAST_EXON_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_LAST_EXON;
				overlapEnd = end;
			} else {
				overlapStart = start;
				overlapEnd = end + EXTENSION_FACTOR_AFTER_LAST_EXON_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_LAST_EXON;
			}
		}
		
		int[] overlapPoints = new int[] { overlapStart, overlapEnd };
		return overlapPoints;
	}
	*/
	
	public static BEDFileParser initiateTranscriptsParser(String transcriptsFile) throws IOException {
		BEDFileParser transcripts = new BEDFileParser(transcriptsFile);
		return transcripts;
	}
	
	public static Map<String, Set<RefSeqGene>> findCoordinatesOfEndsOfTranscripts(BEDFileParser transcripts, String readsAnnotationFile, String transcriptsFile, String chrSizesFile) throws IOException, ParseException {
		AnnotationReader<? extends GenomicAnnotation> readsAnnotations = AnnotationReaderFactory.create(readsAnnotationFile, BED_FORMAT);
		System.out.println("READS ANNOTATIONS SIZE: " + readsAnnotations.getAnnotationList().size());

		Map<String, Set<RefSeqGene>> coordinatesOfTranscriptEndsWithTranscript = new HashMap<String, Set<RefSeqGene>>();
		
		System.out.println("Starting to find end coordinates...");
		
		System.out.println("Finding transcripts to handle...");
		
		Set<RefSeqGene> transcriptsToHandle = new HashSet<RefSeqGene>();
		List<? extends GenomicAnnotation> readsList = readsAnnotations.getAnnotationList();
		int count = 0;
		for (GenomicAnnotation read : readsList) {
			int[] overlapPoints = findOverlapPointsOfTranscriptsForGivenReadWithPolyA(read);
			int overlapStart = overlapPoints[0];
			int overlapEnd = overlapPoints[1];
			
			IntervalTree<RefSeqGeneWithIsoforms> transcriptsTree = transcripts.getChrTree(read.getChromosome());
			if(transcriptsTree == null) {
				continue;
			}
			
			Iterator<? extends Node<RefSeqGeneWithIsoforms>> transcriptsOverlapperIt = transcriptsTree.overlappers(overlapStart, overlapEnd);
			while (transcriptsOverlapperIt.hasNext()) {
				Node<RefSeqGeneWithIsoforms> transcriptOverlapperNode = transcriptsOverlapperIt.next();
				RefSeqGeneWithIsoforms transcriptsOverlap = transcriptOverlapperNode.getValue();
				Collection<RefSeqGene> transcriptsOverlapCollection = transcriptsOverlap.getAllIsoforms();
				for (RefSeqGene transcriptOverlap : transcriptsOverlapCollection) {
					transcriptsToHandle.add(transcriptOverlap);
				}
			}
			
			count++;
			if (count % 10000 == 0)
				System.out.println(count + " / " + readsList.size());
		}
		
		System.out.println("Done finding transcripts to handle...");

		count = 0;
		for (RefSeqGene transcriptAnnot : transcriptsToHandle){
			count++;
			if (count % 10000 == 0)
				System.out.println(count + " / " + transcriptsToHandle.size());
			
			if (transcriptAnnot.getOrientation().equals(ASTERISK))
				continue;
			
			String chr = transcriptAnnot.getChr();
			
			IntervalTree<? extends GenomicAnnotation> readsAnnotationsTree = readsAnnotations.getChromosomeTree(transcriptAnnot.getChr());
			if(readsAnnotationsTree == null) {
				continue;
			}

			Alignments lastExon = transcriptAnnot.getOrientedLastExon();
			int transcriptAnnotEnd = transcriptAnnot.getOrientedEnd();
			
			boolean transcriptAnnotEndIsOnRight;
			int overlapStart;
			int overlapEnd;
			if (lastExon.getStart() == transcriptAnnotEnd) {
				overlapStart = transcriptAnnotEnd - EXTENSION_FACTOR_AFTER_LAST_EXON_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_LAST_EXON;
				overlapEnd = lastExon.getEnd();
				transcriptAnnotEndIsOnRight = false;
			} else {
				overlapStart = lastExon.getStart();
				overlapEnd = transcriptAnnotEnd + EXTENSION_FACTOR_AFTER_LAST_EXON_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_LAST_EXON;
				transcriptAnnotEndIsOnRight = true;
			}

			Set<Integer> possibleCoordinatesOfTranscriptEnds = new HashSet<Integer>();
			
			Iterator<? extends Node<? extends GenomicAnnotation>> readsAnnotationsOverlaperIt = readsAnnotationsTree.overlappers(overlapStart, overlapEnd);
			while(readsAnnotationsOverlaperIt.hasNext()) {
				Node<? extends GenomicAnnotation> readOverlapperNode = readsAnnotationsOverlaperIt.next();
				GenomicAnnotation readAnnot = readOverlapperNode.getValue();
				
				int start = readAnnot.getStart();
				int end = readAnnot.getEnd();
				String orientation = readAnnot.getOrientation();
				boolean polyAIsOnRight = orientation.equals(PLUS);
				
				int coor;
				if (polyAIsOnRight)
					coor = end;
				else
					coor = start;
				
				if (!(coor >= overlapStart && coor <= overlapEnd))
					continue;
				
				if (polyAIsOnRight != transcriptAnnotEndIsOnRight)
					continue;
				
				possibleCoordinatesOfTranscriptEnds.add(coor);
			}
			
			List<Integer> possibleCoordinatesOfTranscriptEndsList = new ArrayList<Integer>(possibleCoordinatesOfTranscriptEnds);
			Collections.sort(possibleCoordinatesOfTranscriptEndsList);
			
			Set<Integer> possibleCoordinatesToRemove = new HashSet<Integer>();
			
			if (transcriptAnnotEndIsOnRight) {
				for (int i = 0; i < possibleCoordinatesOfTranscriptEndsList.size(); i++) {
					int coor = possibleCoordinatesOfTranscriptEndsList.get(i);
					
					if (i + 1 < possibleCoordinatesOfTranscriptEndsList.size()) {
						int nextCoor = possibleCoordinatesOfTranscriptEndsList.get(i + 1);
						
						if (nextCoor - coor <= TYPICAL_ORIGINAL_ALIGNMENT_LENGTH)
							possibleCoordinatesToRemove.add(coor);
					}
				}
			} else {
				for (int i = 0; i < possibleCoordinatesOfTranscriptEndsList.size(); i++) {
					int coor = possibleCoordinatesOfTranscriptEndsList.get(i);
					
					if (i > 0) {
						int prevCoor = possibleCoordinatesOfTranscriptEndsList.get(i - 1);
						
						if (coor - prevCoor <= TYPICAL_ORIGINAL_ALIGNMENT_LENGTH)
							possibleCoordinatesToRemove.add(coor);
					}
				}
			}

			possibleCoordinatesOfTranscriptEnds.removeAll(possibleCoordinatesToRemove);
			
			for (int coor : possibleCoordinatesOfTranscriptEnds) {
				String completeCoor = chr + COLON + String.valueOf(coor);
				
				if (coordinatesOfTranscriptEndsWithTranscript.containsKey(completeCoor)) {
					coordinatesOfTranscriptEndsWithTranscript.get(completeCoor).add(transcriptAnnot);
				} else {
					Set<RefSeqGene> s = new HashSet<RefSeqGene>();
					s.add(transcriptAnnot);
					coordinatesOfTranscriptEndsWithTranscript.put(completeCoor, s);
				}
			}
		}
		System.out.println("count: " + count);
		System.out.println("TOTAL coordinates: " + coordinatesOfTranscriptEndsWithTranscript.size());

		return coordinatesOfTranscriptEndsWithTranscript;
	}
	
	/**
	public static Map<String, Set<RefSeqGene>> findCoordinatesOfEndsOfTranscripts(BEDFileParser transcripts, String readsAnnotationFile, String transcriptsFile, String chrSizesFile) throws IOException, ParseException {
		GenericAlignmentDataModel alignmentsStartA = new GenericAlignmentDataModel(READS_WITH_POLYA_START_WITH_A_ALIGNMENT_SORTED_FILE, chrSizesFile, false);
		GenericAlignmentDataModel alignmentsStartT = new GenericAlignmentDataModel(READS_WITH_POLYA_START_WITH_T_ALIGNMENT_SORTED_FILE, chrSizesFile, false);
		GenericAlignmentDataModel alignmentsEndA = new GenericAlignmentDataModel(READS_WITH_POLYA_END_WITH_A_ALIGNMENT_SORTED_FILE, chrSizesFile, false);
		GenericAlignmentDataModel alignmentsEndT = new GenericAlignmentDataModel(READS_WITH_POLYA_END_WITH_T_ALIGNMENT_SORTED_FILE, chrSizesFile, false);

		AnnotationReader<? extends GenomicAnnotation> readsAnnotations = AnnotationReaderFactory.create(readsAnnotationFile, BED_FORMAT);
		System.out.println("READS ANNOTATIONS SIZE: " + readsAnnotations.getAnnotationList().size());

		Map<String, Set<RefSeqGene>> coordinatesOfTranscriptEndsWithTranscript = new HashMap<String, Set<RefSeqGene>>();
		
		System.out.println("Starting to find end coordinates...");
		
		System.out.println("Finding transcripts to handle...");
		
		Set<RefSeqGene> transcriptsToHandle = new HashSet<RefSeqGene>();
		List<? extends GenomicAnnotation> readsList = readsAnnotations.getAnnotationList();
		int count = 0;
		for (GenomicAnnotation read : readsList) {
			int[] overlapPoints = findOverlapPointsOfTranscriptsForGivenReadWithPolyA(read, alignmentsStartA, alignmentsStartT, alignmentsEndA, alignmentsEndT);
			int overlapStart = overlapPoints[0];
			int overlapEnd = overlapPoints[1];
			
			IntervalTree<RefSeqGeneWithIsoforms> transcriptsTree = transcripts.getChrTree(read.getChromosome());
			if(transcriptsTree == null) {
				continue;
			}
			
			Iterator<? extends Node<RefSeqGeneWithIsoforms>> transcriptsOverlapperIt = transcriptsTree.overlappers(overlapStart, overlapEnd);
			while (transcriptsOverlapperIt.hasNext()) {
				Node<RefSeqGeneWithIsoforms> transcriptOverlapperNode = transcriptsOverlapperIt.next();
				RefSeqGeneWithIsoforms transcriptsOverlap = transcriptOverlapperNode.getValue();
				Collection<RefSeqGene> transcriptsOverlapCollection = transcriptsOverlap.getAllIsoforms();
				for (RefSeqGene transcriptOverlap : transcriptsOverlapCollection) {
					transcriptsToHandle.add(transcriptOverlap);
				}
			}
			
			count++;
			if (count % 10000 == 0)
				System.out.println(count + " / " + readsList.size());
		}
		
		System.out.println("Done finding transcripts to handle...");

		count = 0;
		for (RefSeqGene transcriptAnnot : transcriptsToHandle){
			if (transcriptAnnot.getOrientation().equals(ASTERISK))
				continue;
			
			String chr = transcriptAnnot.getChr();
			
			IntervalTree<? extends GenomicAnnotation> readsAnnotationsTree = readsAnnotations.getChromosomeTree(transcriptAnnot.getChr());
			if(readsAnnotationsTree == null) {
				continue;
			}

			Alignments lastExon = transcriptAnnot.getOrientedLastExon();
			int transcriptAnnotEnd = transcriptAnnot.getOrientedEnd();
			
			boolean transcriptAnnotEndIsOnRight;
			int overlapStart;
			int overlapEnd;
			if (lastExon.getStart() == transcriptAnnotEnd) {
				overlapStart = transcriptAnnotEnd - EXTENSION_FACTOR_AFTER_LAST_EXON_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_LAST_EXON;
				overlapEnd = lastExon.getEnd();
				transcriptAnnotEndIsOnRight = false;
			} else {
				overlapStart = lastExon.getStart();
				overlapEnd = transcriptAnnotEnd + EXTENSION_FACTOR_AFTER_LAST_EXON_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_LAST_EXON;
				transcriptAnnotEndIsOnRight = true;
			}

			Set<Integer> possibleCoordinatesOfTranscriptEnds = new HashSet<Integer>();
			
			Iterator<? extends Node<? extends GenomicAnnotation>> readsAnnotationsOverlaperIt = readsAnnotationsTree.overlappers(overlapStart, overlapEnd);
			while(readsAnnotationsOverlaperIt.hasNext()) {
				Node<? extends GenomicAnnotation> readOverlapperNode = readsAnnotationsOverlaperIt.next();
				GenomicAnnotation readAnnot = readOverlapperNode.getValue();
				
				int start = readAnnot.getStart();
				int end = readAnnot.getEnd();
				Alignments a1 = new Alignments(chr, start, end);

				CloseableIterator<Alignment> overlappingAlignmentsStartA = alignmentsStartA.getAlignmentsOverlappingRegion(a1);
				while(overlappingAlignmentsStartA.hasNext()) {
					Alignment readWithPolyAAlignment = overlappingAlignmentsStartA.next();

					findCoordinatesOfEndsOfTranscriptsHandleOverlappingAlignment(readWithPolyAAlignment, chr, transcriptAnnot, start, end, overlapStart, overlapEnd, transcriptAnnotEndIsOnRight, possibleCoordinatesOfTranscriptEnds, true);
				}
				overlappingAlignmentsStartA.close();
				
				CloseableIterator<Alignment> overlappingAlignmentsStartT = alignmentsStartT.getAlignmentsOverlappingRegion(a1);
				while(overlappingAlignmentsStartT.hasNext()) {
					Alignment readWithPolyAAlignment = overlappingAlignmentsStartT.next();

					findCoordinatesOfEndsOfTranscriptsHandleOverlappingAlignment(readWithPolyAAlignment, chr, transcriptAnnot, start, end, overlapStart, overlapEnd, transcriptAnnotEndIsOnRight, possibleCoordinatesOfTranscriptEnds, true);
				}
				overlappingAlignmentsStartT.close();
				
				CloseableIterator<Alignment> overlappingAlignmentsEndA = alignmentsEndA.getAlignmentsOverlappingRegion(a1);
				while(overlappingAlignmentsEndA.hasNext()) {
					Alignment readWithPolyAAlignment = overlappingAlignmentsEndA.next();

					findCoordinatesOfEndsOfTranscriptsHandleOverlappingAlignment(readWithPolyAAlignment, chr, transcriptAnnot, start, end, overlapStart, overlapEnd, transcriptAnnotEndIsOnRight, possibleCoordinatesOfTranscriptEnds, false);
				}
				overlappingAlignmentsEndA.close();
				
				CloseableIterator<Alignment> overlappingAlignmentsEndT = alignmentsEndT.getAlignmentsOverlappingRegion(a1);
				while(overlappingAlignmentsEndT.hasNext()) {
					Alignment readWithPolyAAlignment = overlappingAlignmentsEndT.next();
					
					findCoordinatesOfEndsOfTranscriptsHandleOverlappingAlignment(readWithPolyAAlignment, chr, transcriptAnnot, start, end, overlapStart, overlapEnd, transcriptAnnotEndIsOnRight, possibleCoordinatesOfTranscriptEnds, false);
				}
				overlappingAlignmentsEndT.close();
				
			}
			
			List<Integer> possibleCoordinatesOfTranscriptEndsList = new ArrayList<Integer>(possibleCoordinatesOfTranscriptEnds);
			Collections.sort(possibleCoordinatesOfTranscriptEndsList);
			
			Set<Integer> possibleCoordinatesToRemove = new HashSet<Integer>();
			
			if (transcriptAnnotEndIsOnRight) {
				for (int i = 0; i < possibleCoordinatesOfTranscriptEndsList.size(); i++) {
					int coor = possibleCoordinatesOfTranscriptEndsList.get(i);
					
					if (i + 1 < possibleCoordinatesOfTranscriptEndsList.size()) {
						int nextCoor = possibleCoordinatesOfTranscriptEndsList.get(i + 1);
						
						if (nextCoor - coor <= TYPICAL_ORIGINAL_ALIGNMENT_LENGTH)
							possibleCoordinatesToRemove.add(coor);
					}
				}
			} else {
				for (int i = 0; i < possibleCoordinatesOfTranscriptEndsList.size(); i++) {
					int coor = possibleCoordinatesOfTranscriptEndsList.get(i);
					
					if (i > 0) {
						int prevCoor = possibleCoordinatesOfTranscriptEndsList.get(i - 1);
						
						if (coor - prevCoor <= TYPICAL_ORIGINAL_ALIGNMENT_LENGTH)
							possibleCoordinatesToRemove.add(coor);
					}
				}
			}

			possibleCoordinatesOfTranscriptEnds.removeAll(possibleCoordinatesToRemove);
			
			for (int coor : possibleCoordinatesOfTranscriptEnds) {
				String completeCoor = chr + COLON + String.valueOf(coor);
				
				if (coordinatesOfTranscriptEndsWithTranscript.containsKey(completeCoor)) {
					coordinatesOfTranscriptEndsWithTranscript.get(completeCoor).add(transcriptAnnot);
				} else {
					Set<RefSeqGene> s = new HashSet<RefSeqGene>();
					s.add(transcriptAnnot);
					coordinatesOfTranscriptEndsWithTranscript.put(completeCoor, s);
				}
			}
			
			count++;
			if (count % 10000 == 0)
				System.out.println(count + " / " + transcriptsToHandle.size());
		}
		System.out.println("count: " + count);
		System.out.println("TOTAL coordinates: " + coordinatesOfTranscriptEndsWithTranscript.size());

		return coordinatesOfTranscriptEndsWithTranscript;
	}
	*/
	
	/**
	private static LightweightGenomicAnnotation findLastExon(RefSeqGene annot) {
		LightweightGenomicAnnotation lastExon = null;
		int transcriptAnnotEnd = annot.getOrientedEnd();
		List<? extends GenomicAnnotation> transcriptAnnotBlocks = annot.getBlocks();
		Iterator<? extends GenomicAnnotation> transcriptAnnotBlocksIt = transcriptAnnotBlocks.iterator();
		while (transcriptAnnotBlocksIt.hasNext()) {
			GenomicAnnotation block = transcriptAnnotBlocksIt.next();
			if (block.getStart() == transcriptAnnotEnd || block.getEnd() == transcriptAnnotEnd) {
				lastExon = block;
				break;
			}
		}
		return lastExon;
	}
	*/
	
	public static void processEndOfTranscriptCoordinates(BEDFileParser transcripts, String transcriptsFile, Map<String, Set<RefSeqGene>> coordinatesOfTranscriptEndsWithTranscript, String matesOfFullPolyAReadsFile, String endCoordinatesOutFile, String endCoordinatesOnlyWithPolyAMatesOutFile, String transcriptLinksOutFile, String transcriptsModifiedOutFile, String transcriptsModifiedOnlyWithPolyAMatesOutFile, String transcriptsModifiedOnlyWithModifiedOutFile, boolean analyzeAndWriteAnalysisOfResults, String transcriptsScoredFile, String transcriptNamesOfLabeledTranscriptsWithNumOfEndpointsOutFile, String transcriptNamesOfAllTranscriptsWithWhetherTheyAreLabeledWithEndpointOutFile, String percentageOfAllTranscriptsLabeledWithEndpointOutFile, String transcriptNamesOfLabeledTranscriptsWithRPKMOutFile, String transcriptNamesOfLabeledTranscriptsWithPValueOutFile, String percentageOfLabeledTranscriptsWithSignificantPValueOutFile, String readsAnnotationsFile, String chrSizesFile, String allEndpointsWithNumOfReadsOutFile, String transcriptNamesWithNumOfReadsForTheEndpointWithMostReadsOutFile) throws IOException, ParseException {
		Map<RefSeqGene, Set<String>> transcriptEndCoors = new HashMap<RefSeqGene, Set<String>>();
		Map<RefSeqGene, Set<RefSeqGene>> transcriptLinks = new HashMap<RefSeqGene, Set<RefSeqGene>>();
		Set<RefSeqGene> transcriptsThatAreLinkedTo = new HashSet<RefSeqGene>();

		System.out.println("Starting to process end of transcript coordinates...");
				
		int count = 0;
		for (String completeCoor : coordinatesOfTranscriptEndsWithTranscript.keySet()) {			
			Set<RefSeqGene> transcriptsForCoor = coordinatesOfTranscriptEndsWithTranscript.get(completeCoor);
			
			String[] completeCoorSplit = completeCoor.split(COLON);
			String chr = completeCoorSplit[0];
			int coor = Integer.parseInt(completeCoorSplit[1]);
			
			IntervalTree<RefSeqGeneWithIsoforms> transcriptsTree = transcripts.getChrTree(chr);
			if(transcriptsTree == null) {
				continue;
			}
			
			Set<RefSeqGene> transcriptsOverlappingCoorWithLastExon = new HashSet<RefSeqGene>();
			Iterator<? extends Node<RefSeqGeneWithIsoforms>> transcriptsOverlapperIt = transcriptsTree.overlappers(coor, coor + 1);
			while (transcriptsOverlapperIt.hasNext()) {
				Node<RefSeqGeneWithIsoforms> transcriptOverlapperNode = transcriptsOverlapperIt.next();
				RefSeqGeneWithIsoforms transcriptsOverlap = transcriptOverlapperNode.getValue();
				Collection<RefSeqGene> transcriptsOverlapCollection = transcriptsOverlap.getAllIsoforms();
				for (RefSeqGene transcriptOverlap : transcriptsOverlapCollection) {
					Alignments transcriptOverlapLastExon = transcriptOverlap.getOrientedLastExon();
				
					if (coor >= transcriptOverlapLastExon.getStart() && coor <= transcriptOverlapLastExon.getEnd())
						transcriptsOverlappingCoorWithLastExon.add(transcriptOverlap);
				}
			}

			for (RefSeqGene transcript : transcriptsForCoor) {
				Alignments lastExon = transcript.getOrientedLastExon();
				
				LightweightGenomicAnnotation transcriptWithoutLastExon;
				if (transcript.getOrientedEnd() == lastExon.getEnd()) {
					int start = transcript.getStart();
					int end = lastExon.getStart();
					transcriptWithoutLastExon = new BasicLightweightAnnotation(chr, start, end);
				} else {
					int start = lastExon.getEnd();
					int end = transcript.getEnd();
					transcriptWithoutLastExon = new BasicLightweightAnnotation(chr, start, end);
				}
				
				boolean transcriptOverlapsCoor = transcriptsOverlappingCoorWithLastExon.contains(transcript);
								
				if (transcriptOverlapsCoor) {
					if (transcriptEndCoors.containsKey(transcript)) {
						transcriptEndCoors.get(transcript).add(completeCoor);
					} else {
						Set<String> s = new HashSet<String>();
						s.add(completeCoor);
						transcriptEndCoors.put(transcript, s);
					}
				} else {
					if (transcriptsOverlappingCoorWithLastExon.size() > 0) {
						for (RefSeqGene t : transcriptsOverlappingCoorWithLastExon) {		
							if (!t.getOrientation().equals(ASTERISK) && !t.getOrientation().equals(transcript.getOrientation()))
								continue;
							
							if (t.overlapsExon(transcriptWithoutLastExon))
								continue;
														
							if (transcriptLinks.containsKey(transcript)) {
								transcriptLinks.get(transcript).add(t);
							} else {
								Set<RefSeqGene> s = new HashSet<RefSeqGene>();
								s.add(t);
								transcriptLinks.put(transcript, s);
							}
							
							transcriptsThatAreLinkedTo.add(t);
							
							int transcriptsInLinkStart;
							int transcriptsInLinkEnd;
							if (transcript.getOrientation().equals(PLUS)) {
								transcriptsInLinkStart = transcript.getEnd();
								transcriptsInLinkEnd = t.getEnd();
							} else {
								transcriptsInLinkStart = t.getStart();
								transcriptsInLinkEnd = transcript.getStart();
							}
							
							Iterator<? extends Node<RefSeqGeneWithIsoforms>> transcriptsInLinkIt = transcriptsTree.overlappers(transcriptsInLinkStart, transcriptsInLinkEnd);
							while (transcriptsInLinkIt.hasNext()) {
								Node<RefSeqGeneWithIsoforms> transcriptsInLinkNode = transcriptsInLinkIt.next();
								RefSeqGeneWithIsoforms transcriptsInLink = transcriptsInLinkNode.getValue();
								Collection<RefSeqGene> transcriptsInLinkCollection = transcriptsInLink.getAllIsoforms();
								for (RefSeqGene transcriptInLink : transcriptsInLinkCollection) {
									if (transcript.getOrientation().equals(PLUS)) {
										if (transcriptInLink.getStart() >= transcript.getEnd() && transcriptInLink.getEnd() <= t.getEnd())
											transcriptsThatAreLinkedTo.add(transcriptInLink);
									} else {
										if (transcriptInLink.getStart() >= t.getStart() && transcriptInLink.getEnd() <= transcript.getStart())
											transcriptsThatAreLinkedTo.add(transcriptInLink);
									}
								}
							}
							
							if (transcriptEndCoors.containsKey(transcript)) {
								transcriptEndCoors.get(transcript).add(completeCoor);
							} else {
								Set<String> s = new HashSet<String>();
								s.add(completeCoor);
								transcriptEndCoors.put(transcript, s);
							}
						}
					} else {
						if (transcriptEndCoors.containsKey(transcript)) {
							transcriptEndCoors.get(transcript).add(completeCoor);
						} else {
							Set<String> s = new HashSet<String>();
							s.add(completeCoor);
							transcriptEndCoors.put(transcript, s);
						}
					}
				}
			}
			
			count++;
			if (count % 100 == 0)
				System.out.println(count + " / " + coordinatesOfTranscriptEndsWithTranscript.size());
		}

		System.out.println("Starting to process linking of last exons to nearby single exon transcripts...");
		
		Map<RefSeqGene, Set<String>> transcriptEndCoorsToModify = new HashMap<RefSeqGene, Set<String>>();
		
		for (RefSeqGene transcriptWithSingleExon : transcriptEndCoors.keySet()) {
			if (transcriptWithSingleExon.getNumExons() != 1)
				continue;
						
			Set<String> transcriptWithSingleExonEndCoorsToHandle = new HashSet<String>();
			
			for (String completeCoor : transcriptEndCoors.get(transcriptWithSingleExon)) {
				String[] completeCoorSplit = completeCoor.split(COLON);
				String chr = completeCoorSplit[0];
				int coor = Integer.parseInt(completeCoorSplit[1]);
				
				if (coor >= transcriptWithSingleExon.getStart() && coor <= transcriptWithSingleExon.getEnd())
					transcriptWithSingleExonEndCoorsToHandle.add(completeCoor);
			}
			
			if (transcriptWithSingleExonEndCoorsToHandle.size() == 0)
				continue;
			
			String chr = transcriptWithSingleExon.getChr();
			
			IntervalTree<RefSeqGeneWithIsoforms> transcriptsTree = transcripts.getChrTree(chr);
			if(transcriptsTree == null) {
				continue;
			}
			
			int overlapStart;
			int overlapEnd;
			if (transcriptWithSingleExon.getOrientation().equals(PLUS)) {
				overlapStart = transcriptWithSingleExon.getStart() - MAX_DISTANCE_BETWEEN_LAST_EXON_AND_NEXT_EXON_FOR_LINKING_TRANSCRIPTS;
				overlapEnd = transcriptWithSingleExon.getEnd();
			} else if (transcriptWithSingleExon.getOrientation().equals(MINUS)) {
				overlapStart = transcriptWithSingleExon.getStart();
				overlapEnd = transcriptWithSingleExon.getEnd() + MAX_DISTANCE_BETWEEN_LAST_EXON_AND_NEXT_EXON_FOR_LINKING_TRANSCRIPTS;
			} else {
				overlapStart = transcriptWithSingleExon.getStart() - MAX_DISTANCE_BETWEEN_LAST_EXON_AND_NEXT_EXON_FOR_LINKING_TRANSCRIPTS;
				overlapEnd = transcriptWithSingleExon.getEnd() + MAX_DISTANCE_BETWEEN_LAST_EXON_AND_NEXT_EXON_FOR_LINKING_TRANSCRIPTS;
			}
			
			Iterator<? extends Node<RefSeqGeneWithIsoforms>> transcriptsOverlapperIt = transcriptsTree.overlappers(overlapStart, overlapEnd);
			while (transcriptsOverlapperIt.hasNext()) {
				Node<RefSeqGeneWithIsoforms> transcriptOverlapperNode = transcriptsOverlapperIt.next();
				RefSeqGeneWithIsoforms transcriptsOverlap = transcriptOverlapperNode.getValue();
				Collection<RefSeqGene> transcriptsOverlapCollection = transcriptsOverlap.getAllIsoforms();
				for (RefSeqGene transcriptOverlap : transcriptsOverlapCollection) {					
					if (!transcriptWithSingleExon.getOrientation().equals(ASTERISK) && !transcriptWithSingleExon.getOrientation().equals(transcriptOverlap.getOrientation()))
						continue;
					
					if (transcriptOverlap.equals(transcriptWithSingleExon))
						continue;
					
					Alignments transcriptOverlapLastExon = transcriptOverlap.getOrientedLastExon();
					
					if (transcriptOverlap.getOrientation().equals(PLUS)) {
						if (!(transcriptOverlapLastExon.getEnd() >= transcriptWithSingleExon.getStart() - MAX_DISTANCE_BETWEEN_LAST_EXON_AND_NEXT_EXON_FOR_LINKING_TRANSCRIPTS && transcriptOverlapLastExon.getEnd() <= transcriptWithSingleExon.getEnd()))
							continue;
					} else {
						if (!(transcriptOverlapLastExon.getStart() >= transcriptWithSingleExon.getStart() && transcriptOverlapLastExon.getStart() <= transcriptWithSingleExon.getEnd() + MAX_DISTANCE_BETWEEN_LAST_EXON_AND_NEXT_EXON_FOR_LINKING_TRANSCRIPTS))
							continue;
					}
					
					LightweightGenomicAnnotation transcriptOverlapWithoutLastExon;
					if (transcriptOverlap.getOrientedEnd() == transcriptOverlapLastExon.getEnd()) {
						int start = transcriptOverlap.getStart();
						int end = transcriptOverlapLastExon.getStart();
						transcriptOverlapWithoutLastExon = new BasicLightweightAnnotation(chr, start, end);
					} else {
						int start = transcriptOverlapLastExon.getEnd();
						int end = transcriptOverlap.getEnd();
						transcriptOverlapWithoutLastExon = new BasicLightweightAnnotation(chr, start, end);
					}
					
					if (transcriptWithSingleExon.overlapsExon(transcriptOverlapWithoutLastExon))
						continue;
					
					if (transcriptLinks.containsKey(transcriptOverlap)) {
						transcriptLinks.get(transcriptOverlap).add(transcriptWithSingleExon);
					} else {
						Set<RefSeqGene> s = new HashSet<RefSeqGene>();
						s.add(transcriptWithSingleExon);
						transcriptLinks.put(transcriptOverlap, s);
					}
					
					transcriptsThatAreLinkedTo.add(transcriptWithSingleExon);
					
					int transcriptsInLinkStart;
					int transcriptsInLinkEnd;
					if (transcriptOverlap.getOrientation().equals(PLUS)) {
						transcriptsInLinkStart = transcriptOverlap.getEnd();
						transcriptsInLinkEnd = transcriptWithSingleExon.getEnd();
					} else {
						transcriptsInLinkStart = transcriptWithSingleExon.getStart();
						transcriptsInLinkEnd = transcriptOverlap.getStart();
					}
					
					Iterator<? extends Node<RefSeqGeneWithIsoforms>> transcriptsInLinkIt = transcriptsTree.overlappers(transcriptsInLinkStart, transcriptsInLinkEnd);
					while (transcriptsInLinkIt.hasNext()) {
						Node<RefSeqGeneWithIsoforms> transcriptsInLinkNode = transcriptsInLinkIt.next();
						RefSeqGeneWithIsoforms transcriptsInLink = transcriptsInLinkNode.getValue();
						Collection<RefSeqGene> transcriptsInLinkCollection = transcriptsInLink.getAllIsoforms();
						for (RefSeqGene transcriptInLink : transcriptsInLinkCollection) {
							if (transcriptOverlap.getOrientation().equals(PLUS)) {
								if (transcriptInLink.getStart() >= transcriptOverlap.getEnd() && transcriptInLink.getEnd() <= transcriptWithSingleExon.getEnd())
									transcriptsThatAreLinkedTo.add(transcriptInLink);
							} else {
								if (transcriptInLink.getStart() >= transcriptWithSingleExon.getStart() && transcriptInLink.getEnd() <= transcriptOverlap.getStart())
									transcriptsThatAreLinkedTo.add(transcriptInLink);
							}
						}
					}
					
					for (String completeCoor : transcriptWithSingleExonEndCoorsToHandle) {
						if (transcriptEndCoorsToModify.containsKey(transcriptOverlap)) {
							transcriptEndCoorsToModify.get(transcriptOverlap).add(completeCoor);
						} else {
							Set<String> s = new HashSet<String>();
							s.add(completeCoor);
							transcriptEndCoorsToModify.put(transcriptOverlap, s);
						}
					}
				}
			}
		}
				
		for (RefSeqGene transcript : transcriptEndCoorsToModify.keySet()) {
			Set<String> transcriptCoors = transcriptEndCoorsToModify.get(transcript);
			
			if (transcriptEndCoors.containsKey(transcript)) {
				for (String completeCoor : transcriptCoors) {
					transcriptEndCoors.get(transcript).add(completeCoor);					
				}
			} else {
				Set<String> s = new HashSet<String>();
				for (String completeCoor : transcriptCoors) {
					s.add(completeCoor);					
				}
				transcriptEndCoors.put(transcript, s);
			}
		}
		
		System.out.println("Done processing end of transcript coordinates...");

		Set<String> allEndCoordinates = new HashSet<String>();
		for (RefSeqGene g : transcriptEndCoors.keySet()) {
			Set<String> s = transcriptEndCoors.get(g);
			for (String completeCoor : s) {
				allEndCoordinates.add(completeCoor);
			}
		}
		
		Map<String, Integer> coordinatesWithEndOfTranscriptNearMateOfFullPolyARead = findCoordinatesWithEndOfTranscriptNearMateOfFullPolyARead(allEndCoordinates, matesOfFullPolyAReadsFile);
		
		System.out.println("Starting to write information to files...");
		
		writeEndOfTranscriptInformationToFiles(transcripts, transcriptEndCoors, allEndCoordinates, transcriptLinks, transcriptsThatAreLinkedTo, coordinatesWithEndOfTranscriptNearMateOfFullPolyARead, endCoordinatesOutFile, endCoordinatesOnlyWithPolyAMatesOutFile, transcriptLinksOutFile, transcriptsFile, transcriptsModifiedOutFile, transcriptsModifiedOnlyWithPolyAMatesOutFile, transcriptsModifiedOnlyWithModifiedOutFile);
		
		if (analyzeAndWriteAnalysisOfResults)
			analyzeAndWriteAnalysisOfResults(transcripts, transcriptsScoredFile, transcriptEndCoors, transcriptNamesOfLabeledTranscriptsWithNumOfEndpointsOutFile, transcriptNamesOfAllTranscriptsWithWhetherTheyAreLabeledWithEndpointOutFile, percentageOfAllTranscriptsLabeledWithEndpointOutFile, transcriptNamesOfLabeledTranscriptsWithRPKMOutFile, transcriptNamesOfLabeledTranscriptsWithPValueOutFile, percentageOfLabeledTranscriptsWithSignificantPValueOutFile, allEndCoordinates, readsAnnotationsFile, chrSizesFile, allEndpointsWithNumOfReadsOutFile, transcriptNamesWithNumOfReadsForTheEndpointWithMostReadsOutFile);
	}
	
	private static void analyzeAndWriteAnalysisOfResults(BEDFileParser transcripts, String transcriptsScoredFile, Map<RefSeqGene, Set<String>> transcriptEndCoors, String transcriptNamesOfLabeledTranscriptsWithNumOfEndpointsOutFile, String transcriptNamesOfAllTranscriptsWithWhetherTheyAreLabeledWithEndpointOutFile, String percentageOfAllTranscriptsLabeledWithEndpointOutFile, String transcriptNamesOfLabeledTranscriptsWithRPKMOutFile, String transcriptNamesOfLabeledTranscriptsWithPValueOutFile, String percentageOfLabeledTranscriptsWithSignificantPValueOutFile, Set<String> allEndCoordinates, String readsAnnotationsFile, String chrSizesFile, String allEndpointsWithNumOfReadsOutFile, String transcriptNamesWithNumOfReadsForTheEndpointWithMostReadsOutFile) throws IOException, ParseException {
		BEDFileParser transcriptsScored = initiateTranscriptsParser(transcriptsScoredFile);
		
		Set<String> transcriptNamesThatHaveEndpoints = new HashSet<String>();
		
		Map<String, List<Integer>> transcriptNamesWithNumOfEndpointsPerIsoform = new HashMap<String, List<Integer>>();
		Map<String, List<Double>> transcriptNamesWithRPKMPerIsoform = new HashMap<String, List<Double>>();
		Map<String, List<Double>> transcriptNamesWithPValuePerIsoform = new HashMap<String, List<Double>>();
		
		for (RefSeqGene transcript : transcriptEndCoors.keySet()) {
			transcriptNamesThatHaveEndpoints.add(transcript.getName());
			
			int numOfEndCoors = transcriptEndCoors.get(transcript).size();
			
			if (transcriptNamesWithNumOfEndpointsPerIsoform.containsKey(transcript.getName())) {
				transcriptNamesWithNumOfEndpointsPerIsoform.get(transcript.getName()).add(numOfEndCoors);
			} else {
				List<Integer> l = new LinkedList<Integer>();
				l.add(numOfEndCoors);
				transcriptNamesWithNumOfEndpointsPerIsoform.put(transcript.getName(), l);
			}
			
			RefSeqGene transcriptScored = null;
			
			IntervalTree<RefSeqGeneWithIsoforms> transcriptsScoredTree = transcriptsScored.getChrTree(transcript.getChr());
			if(transcriptsScoredTree == null) {
				continue;
			}
			
			boolean shouldBreak = false;
			Iterator<? extends Node<RefSeqGeneWithIsoforms>> transcriptsScoredOverlapperIt = transcriptsScoredTree.overlappers(transcript.getStart(), transcript.getStart() + 1);
			while (transcriptsScoredOverlapperIt.hasNext()) {
				Node<RefSeqGeneWithIsoforms> transcriptScoredOverlapperNode = transcriptsScoredOverlapperIt.next();
				RefSeqGeneWithIsoforms transcriptsScoredOverlap = transcriptScoredOverlapperNode.getValue();
				Collection<RefSeqGene> transcriptsScoredOverlapCollection = transcriptsScoredOverlap.getAllIsoforms();
				for (RefSeqGene transcriptScoredOverlap : transcriptsScoredOverlapCollection) {
					if (transcript.equals(transcriptScoredOverlap)) {
						transcriptScored = transcriptScoredOverlap;
						shouldBreak = true;
						break;
					}
				}
				if (shouldBreak)
					break;
			}
			
			double rpkm = transcriptScored.getRPKM();

			if (transcriptNamesWithRPKMPerIsoform.containsKey(transcript.getName())) {
				transcriptNamesWithRPKMPerIsoform.get(transcript.getName()).add(rpkm);
			} else {
				List<Double> l = new LinkedList<Double>();
				l.add(rpkm);
				transcriptNamesWithRPKMPerIsoform.put(transcript.getName(), l);
			}
			
			double pvalue = transcriptScored.getPValue();

			if (transcriptNamesWithPValuePerIsoform.containsKey(transcript.getName())) {
				transcriptNamesWithPValuePerIsoform.get(transcript.getName()).add(pvalue);
			} else {
				List<Double> l = new LinkedList<Double>();
				l.add(pvalue);
				transcriptNamesWithPValuePerIsoform.put(transcript.getName(), l);
			}
		}
		
		Map<String, Double> transcriptNamesWithMedianNumOfEndpoints = new HashMap<String, Double>();
		for (String transcriptName : transcriptNamesWithNumOfEndpointsPerIsoform.keySet()) {
			List<Integer> numOfEndpointsPerIsoform = transcriptNamesWithNumOfEndpointsPerIsoform.get(transcriptName);
			Collections.sort(numOfEndpointsPerIsoform);
			
			if (numOfEndpointsPerIsoform.size() % 2 == 0) {
				int middle = numOfEndpointsPerIsoform.size() / 2;
				double meanOfMiddleTwo = ( (double)numOfEndpointsPerIsoform.get(middle - 1) + (double)numOfEndpointsPerIsoform.get(middle) ) / 2.0;
				transcriptNamesWithMedianNumOfEndpoints.put(transcriptName, meanOfMiddleTwo);
			} else {
				int middle = (int)((double)numOfEndpointsPerIsoform.size() / 2.0);
				double median = (double)numOfEndpointsPerIsoform.get(middle);
				transcriptNamesWithMedianNumOfEndpoints.put(transcriptName, median);
			}
		}
		
		Map<String, Double> transcriptNamesWithMedianRPKM = new HashMap<String, Double>();
		for (String transcriptName : transcriptNamesWithRPKMPerIsoform.keySet()) {
			List<Double> rpkmPerIsoform = transcriptNamesWithRPKMPerIsoform.get(transcriptName);
			Collections.sort(rpkmPerIsoform);
			
			if (rpkmPerIsoform.size() % 2 == 0) {
				int middle = rpkmPerIsoform.size() / 2;
				double meanOfMiddleTwo = ( (double)rpkmPerIsoform.get(middle - 1) + (double)rpkmPerIsoform.get(middle) ) / 2.0;
				transcriptNamesWithMedianRPKM.put(transcriptName, meanOfMiddleTwo);
			} else {
				int middle = (int)((double)rpkmPerIsoform.size() / 2.0);
				double median = (double)rpkmPerIsoform.get(middle);
				transcriptNamesWithMedianRPKM.put(transcriptName, median);
			}
		}
		
		Map<String, Double> transcriptNamesWithMedianPValue = new HashMap<String, Double>();
		for (String transcriptName : transcriptNamesWithPValuePerIsoform.keySet()) {
			List<Double> pvaluePerIsoform = transcriptNamesWithPValuePerIsoform.get(transcriptName);
			Collections.sort(pvaluePerIsoform);
			
			if (pvaluePerIsoform.size() % 2 == 0) {
				int middle = pvaluePerIsoform.size() / 2;
				double meanOfMiddleTwo = ( (double)pvaluePerIsoform.get(middle - 1) + (double)pvaluePerIsoform.get(middle) ) / 2.0;
				transcriptNamesWithMedianPValue.put(transcriptName, meanOfMiddleTwo);
			} else {
				int middle = (int)((double)pvaluePerIsoform.size() / 2.0);
				double median = (double)pvaluePerIsoform.get(middle);
				transcriptNamesWithMedianPValue.put(transcriptName, median);
			}
		}
		
		FileWriter fwTranscriptNamesOfLabeledTranscriptsWithNumOfEndpoints = new FileWriter(transcriptNamesOfLabeledTranscriptsWithNumOfEndpointsOutFile);
		for (String transcriptName : transcriptNamesWithMedianNumOfEndpoints.keySet()) {
			double numOfEndpoints = transcriptNamesWithMedianNumOfEndpoints.get(transcriptName);
			fwTranscriptNamesOfLabeledTranscriptsWithNumOfEndpoints.write(transcriptName + TAB + String.valueOf(numOfEndpoints) + NEW_LINE);
		}
		fwTranscriptNamesOfLabeledTranscriptsWithNumOfEndpoints.close();
		
		Set<String> allTranscriptNames = new HashSet<String>();
		List<RefSeqGene> allTranscripts = transcripts.GetGenes();
		for (RefSeqGene transcript : allTranscripts) {
			allTranscriptNames.add(transcript.getName());
		}
		
		FileWriter fwTranscriptNamesOfAllTranscriptsWithWhetherTheyAreLabeledWithEndpoint = new FileWriter(transcriptNamesOfAllTranscriptsWithWhetherTheyAreLabeledWithEndpointOutFile);
		for (String transcriptName : allTranscriptNames) {
			String out;
			if (transcriptNamesThatHaveEndpoints.contains(transcriptName))
				out = YES;
			else
				out = NO;
			fwTranscriptNamesOfAllTranscriptsWithWhetherTheyAreLabeledWithEndpoint.write(transcriptName + TAB + out + NEW_LINE);
		}
		fwTranscriptNamesOfAllTranscriptsWithWhetherTheyAreLabeledWithEndpoint.close();
		
		FileWriter fwPercentageOfAllTranscriptsLabeledWithEndpoint = new FileWriter(percentageOfAllTranscriptsLabeledWithEndpointOutFile);
		double percentageOfTranscriptsLabeledWithEndpoint = (double)transcriptNamesThatHaveEndpoints.size() / (double)allTranscriptNames.size();
		fwPercentageOfAllTranscriptsLabeledWithEndpoint.write(String.valueOf(percentageOfTranscriptsLabeledWithEndpoint));
		fwPercentageOfAllTranscriptsLabeledWithEndpoint.close();
		
		FileWriter fwTranscriptNamesOfLabeledTranscriptsWithRPKM = new FileWriter(transcriptNamesOfLabeledTranscriptsWithRPKMOutFile);
		for (String transcriptName : transcriptNamesWithMedianRPKM.keySet()) {
			double rpkm = transcriptNamesWithMedianRPKM.get(transcriptName);
			fwTranscriptNamesOfLabeledTranscriptsWithRPKM.write(transcriptName + TAB + String.valueOf(rpkm) + NEW_LINE);
		}
		fwTranscriptNamesOfLabeledTranscriptsWithRPKM.close();
		
		FileWriter fwTranscriptNamesOfLabeledTranscriptsWithPValue = new FileWriter(transcriptNamesOfLabeledTranscriptsWithPValueOutFile);
		for (String transcriptName : transcriptNamesWithMedianPValue.keySet()) {
			double pvalue = transcriptNamesWithMedianPValue.get(transcriptName);
			fwTranscriptNamesOfLabeledTranscriptsWithPValue.write(transcriptName + TAB + String.valueOf(pvalue) + NEW_LINE);
		}
		fwTranscriptNamesOfLabeledTranscriptsWithPValue.close();
		
		FileWriter fwPercentageOfLabeledTranscriptsWithSignificantPValue = new FileWriter(percentageOfLabeledTranscriptsWithSignificantPValueOutFile);
		int numOfTranscriptsWithSignificantPValue = 0;
		for (String transcriptName : transcriptNamesWithMedianPValue.keySet()) {
			double pvalue = transcriptNamesWithMedianPValue.get(transcriptName);
			if (pvalue <= PVALUE_SIGNIFICANCE_CUTOFF)
				numOfTranscriptsWithSignificantPValue++;
		}
		double percentageOfTranscriptsWithSignificantPValue = (double)numOfTranscriptsWithSignificantPValue / (double)transcriptNamesWithMedianPValue.size();
		fwPercentageOfLabeledTranscriptsWithSignificantPValue.write(String.valueOf(percentageOfTranscriptsWithSignificantPValue));
		fwPercentageOfLabeledTranscriptsWithSignificantPValue.close();
	
		GenericAlignmentDataModel alignmentsStartA = new GenericAlignmentDataModel(READS_WITH_POLYA_START_WITH_A_ALIGNMENT_SORTED_FILE, chrSizesFile, false);
		GenericAlignmentDataModel alignmentsStartT = new GenericAlignmentDataModel(READS_WITH_POLYA_START_WITH_T_ALIGNMENT_SORTED_FILE, chrSizesFile, false);
		GenericAlignmentDataModel alignmentsEndA = new GenericAlignmentDataModel(READS_WITH_POLYA_END_WITH_A_ALIGNMENT_SORTED_FILE, chrSizesFile, false);
		GenericAlignmentDataModel alignmentsEndT = new GenericAlignmentDataModel(READS_WITH_POLYA_END_WITH_T_ALIGNMENT_SORTED_FILE, chrSizesFile, false);
		
		AnnotationReader<? extends GenomicAnnotation> readsAnnotations = AnnotationReaderFactory.create(readsAnnotationsFile, BED_FORMAT);
		
		Map<String, Integer> numOfReadsForAllEndCoordinates = new HashMap<String, Integer>();
		
		for (String completeCoor : allEndCoordinates) {
			int numOfReadsForEndpoint = calculateNumOfReadsForEndpoint(completeCoor, readsAnnotations, alignmentsStartA, alignmentsStartT, alignmentsEndA, alignmentsEndT);
			
			numOfReadsForAllEndCoordinates.put(completeCoor, numOfReadsForEndpoint);
		}
		
		Map<RefSeqGene, Integer> transcriptsWithNumOfReadsForEndpointWithMostReads = new HashMap<RefSeqGene, Integer>();
		for (RefSeqGene transcript : transcriptEndCoors.keySet()) {
			Set<String> completeCoorsWithMaxNumOfReadsEndpoint = new HashSet<String>();
			int maxNumOfReadsForEndpoint = -1;
			
			Set<String> endCoors = transcriptEndCoors.get(transcript);
			for (String completeCoor : endCoors) {
				int numOfReadsForThisEndpoint = numOfReadsForAllEndCoordinates.get(completeCoor);
				
				if (numOfReadsForThisEndpoint > maxNumOfReadsForEndpoint) {
					completeCoorsWithMaxNumOfReadsEndpoint.clear();
					completeCoorsWithMaxNumOfReadsEndpoint.add(completeCoor);
					maxNumOfReadsForEndpoint = numOfReadsForThisEndpoint;
				} else if (numOfReadsForThisEndpoint == maxNumOfReadsForEndpoint) {
					completeCoorsWithMaxNumOfReadsEndpoint.add(completeCoor);
				}
			}
			
			transcriptsWithNumOfReadsForEndpointWithMostReads.put(transcript, maxNumOfReadsForEndpoint);
		}
		
		Map<String, List<Integer>> transcriptNamesWithNumOfReadsForEndpointWithMostReadsPerIsoform = new HashMap<String, List<Integer>>();
		for (RefSeqGene transcript : transcriptsWithNumOfReadsForEndpointWithMostReads.keySet()) {
			int numOfReadsForEndpoint = transcriptsWithNumOfReadsForEndpointWithMostReads.get(transcript);
			
			if (transcriptNamesWithNumOfReadsForEndpointWithMostReadsPerIsoform.containsKey(transcript.getName())) {
				transcriptNamesWithNumOfReadsForEndpointWithMostReadsPerIsoform.get(transcript.getName()).add(numOfReadsForEndpoint);
			} else {
				List<Integer> l = new LinkedList<Integer>();
				l.add(numOfReadsForEndpoint);
				transcriptNamesWithNumOfReadsForEndpointWithMostReadsPerIsoform.put(transcript.getName(), l);
			}
		}
		
		Map<String, Double> transcriptNamesWithMedianNumOfReadsForEndpointWithMostReads = new HashMap<String, Double>();
		for (String transcriptName : transcriptNamesWithNumOfReadsForEndpointWithMostReadsPerIsoform.keySet()) {
			List<Integer> numOfReadsForEndpointPerIsoform =transcriptNamesWithNumOfReadsForEndpointWithMostReadsPerIsoform.get(transcriptName);
			Collections.sort(numOfReadsForEndpointPerIsoform);
			
			if (numOfReadsForEndpointPerIsoform.size() % 2 == 0) {
				int middle = numOfReadsForEndpointPerIsoform.size() / 2;
				double meanOfMiddleTwo = ( (double)numOfReadsForEndpointPerIsoform.get(middle - 1) + (double)numOfReadsForEndpointPerIsoform.get(middle) ) / 2.0;
				transcriptNamesWithMedianNumOfReadsForEndpointWithMostReads.put(transcriptName, meanOfMiddleTwo);
			} else {
				int middle = (int)((double)numOfReadsForEndpointPerIsoform.size() / 2.0);
				double median = (double)numOfReadsForEndpointPerIsoform.get(middle);
				transcriptNamesWithMedianNumOfReadsForEndpointWithMostReads.put(transcriptName, median);
			}
		}
		
		FileWriter fwAllEndpointsWithNumOfReads = new FileWriter(allEndpointsWithNumOfReadsOutFile);
		for (String completeCoor : numOfReadsForAllEndCoordinates.keySet()) {
			int numOfReads = numOfReadsForAllEndCoordinates.get(completeCoor);
			fwAllEndpointsWithNumOfReads.write(completeCoor + TAB + String.valueOf(numOfReads) + NEW_LINE);
		}
		fwAllEndpointsWithNumOfReads.close();
		
		FileWriter fwTranscriptNamesWithNumOfReadsForTheEndpointWithMostReads = new FileWriter(transcriptNamesWithNumOfReadsForTheEndpointWithMostReadsOutFile);
		for (String transcriptName : transcriptNamesWithMedianNumOfReadsForEndpointWithMostReads.keySet()) {
			double numOfReads = transcriptNamesWithMedianNumOfReadsForEndpointWithMostReads.get(transcriptName);
			fwTranscriptNamesWithNumOfReadsForTheEndpointWithMostReads.write(transcriptName + TAB + String.valueOf(numOfReads) + NEW_LINE);
		}
		fwTranscriptNamesWithNumOfReadsForTheEndpointWithMostReads.close();
	}
	
	private static int calculateNumOfReadsForEndpoint(String completeCoor, AnnotationReader<? extends GenomicAnnotation> readsAnnotations, GenericAlignmentDataModel alignmentsStartA, GenericAlignmentDataModel alignmentsStartT, GenericAlignmentDataModel alignmentsEndA, GenericAlignmentDataModel alignmentsEndT)throws IOException {
		String[] completeCoorSplit = completeCoor.split(COLON);
		String chr = completeCoorSplit[0];
		int coor = Integer.parseInt(completeCoorSplit[1]);
		
		int numOfReadsForThisEndpoint = 0;
		
		IntervalTree<? extends GenomicAnnotation> readsAnnotationsTree = readsAnnotations.getChromosomeTree(chr);
		if(readsAnnotationsTree == null) {
			return -1;
		}

		Iterator<? extends Node<? extends GenomicAnnotation>> readsAnnotationsOverlapperIt = readsAnnotationsTree.overlappers(coor - 2, coor + 2);
		while(readsAnnotationsOverlapperIt.hasNext()) {
			Node<? extends GenomicAnnotation> readOverlapperNode = readsAnnotationsOverlapperIt.next();
			GenomicAnnotation readAnnot = readOverlapperNode.getValue();
			
			int start = readAnnot.getStart();
			int end = readAnnot.getEnd();

			if (start != coor && end != coor)
				continue;

			Alignments a1 = new Alignments(chr, start, end);

			CloseableIterator<Alignment> overlappingAlignmentsStartA = alignmentsStartA.getAlignmentsOverlappingRegion(a1);
			while(overlappingAlignmentsStartA.hasNext()) {
				Alignment readWithPolyAAlignment = overlappingAlignmentsStartA.next();

				boolean b = checkIfOrientationOfAlignmentMatchesCoor(readWithPolyAAlignment, start, end, coor, true);
				if (b)
					numOfReadsForThisEndpoint++;
			}
			overlappingAlignmentsStartA.close();
			
			CloseableIterator<Alignment> overlappingAlignmentsStartT = alignmentsStartT.getAlignmentsOverlappingRegion(a1);
			while(overlappingAlignmentsStartT.hasNext()) {
				Alignment readWithPolyAAlignment = overlappingAlignmentsStartT.next();

				boolean b = checkIfOrientationOfAlignmentMatchesCoor(readWithPolyAAlignment, start, end, coor, true);
				if (b)
					numOfReadsForThisEndpoint++;
			}
			overlappingAlignmentsStartT.close();
			
			CloseableIterator<Alignment> overlappingAlignmentsEndA = alignmentsEndA.getAlignmentsOverlappingRegion(a1);
			while(overlappingAlignmentsEndA.hasNext()) {
				Alignment readWithPolyAAlignment = overlappingAlignmentsEndA.next();

				boolean b = checkIfOrientationOfAlignmentMatchesCoor(readWithPolyAAlignment, start, end, coor, false);
				if (b)
					numOfReadsForThisEndpoint++;
			}
			overlappingAlignmentsEndA.close();
			
			CloseableIterator<Alignment> overlappingAlignmentsEndT = alignmentsEndT.getAlignmentsOverlappingRegion(a1);
			while(overlappingAlignmentsEndT.hasNext()) {
				Alignment readWithPolyAAlignment = overlappingAlignmentsEndT.next();
				
				boolean b = checkIfOrientationOfAlignmentMatchesCoor(readWithPolyAAlignment, start, end, coor, false);
				if (b)
					numOfReadsForThisEndpoint++;
			}
			overlappingAlignmentsEndT.close();
		}
		
		return numOfReadsForThisEndpoint;
	}
	
	private static boolean checkIfOrientationOfAlignmentMatchesCoor(Alignment readWithPolyAAlignment, int start, int end, int coor, boolean startsAlignment) {
		if (!(readWithPolyAAlignment.getStart() == start && readWithPolyAAlignment.getEnd() == end))
			return false;

		boolean isNegativeStrand = readWithPolyAAlignment.isNegativeStrand(); // is true when 16 (reversed) ; b is false when 0 (correct)
		boolean polyAIsOnRight;
		if (startsAlignment) {
			if (isNegativeStrand)
				polyAIsOnRight = true;
			else
				polyAIsOnRight = false;
		} else {
			if (isNegativeStrand)
				polyAIsOnRight = false;
			else
				polyAIsOnRight = true;
		}
		
		int coorToCheckAgainst = -1;
		if (polyAIsOnRight)
			coorToCheckAgainst = end;
		else
			coorToCheckAgainst = start;

		if (coorToCheckAgainst == coor)
			return true;
		else
			return false;
	}
	
	private static Map<String, Integer> findCoordinatesWithEndOfTranscriptNearMateOfFullPolyARead(Set<String> coordinatesOfTranscriptEnds, String matesOfFullPolyAReadsFile) throws IOException, ParseException {
		Map<String, Integer> coordinatesWithEndOfTranscriptNearMateOfFullPolyARead = new HashMap<String, Integer>();
		
		AnnotationReader<? extends GenomicAnnotation> matesOfFullPolyAReads = AnnotationReaderFactory.create(matesOfFullPolyAReadsFile, BED_FORMAT);

		for (String completeCoor : coordinatesOfTranscriptEnds) {
			String[] completeCoorSplit = completeCoor.split(COLON);
			String chr = completeCoorSplit[0];
			int coor = Integer.parseInt(completeCoorSplit[1]);
			
			int start = coor - EXTENSION_FACTOR_FOR_FILTERING_TRANSCRIPT_END_COORDINATES_USING_MATES_OF_FULL_POLYA_READS;
			int end = coor + EXTENSION_FACTOR_FOR_FILTERING_TRANSCRIPT_END_COORDINATES_USING_MATES_OF_FULL_POLYA_READS;
			
			IntervalTree<? extends GenomicAnnotation> matesOfFullPolyAReadsTree = matesOfFullPolyAReads.getChromosomeTree(chr);
			if(matesOfFullPolyAReadsTree == null) {
				continue;
			}
			
			int numOfMatesOfFullPolyAReads = 0;
			Iterator<? extends Node<? extends GenomicAnnotation>> matesOfFullPolyAReadsOverlapperIt = matesOfFullPolyAReadsTree.overlappers(start, end);
			while (matesOfFullPolyAReadsOverlapperIt.hasNext()) {
				Node<? extends GenomicAnnotation> mateOfFullPolyAReadOverlapperNode = matesOfFullPolyAReadsOverlapperIt.next();
				numOfMatesOfFullPolyAReads++;
			}
			if (numOfMatesOfFullPolyAReads > 0)
				coordinatesWithEndOfTranscriptNearMateOfFullPolyARead.put(completeCoor, numOfMatesOfFullPolyAReads);
		}
		
		return coordinatesWithEndOfTranscriptNearMateOfFullPolyARead;
	}
	
	private static void writeEndOfTranscriptInformationToFiles(BEDFileParser transcripts, Map<RefSeqGene, Set<String>> transcriptEndCoors, Set<String> allEndCoordinates, Map<RefSeqGene, Set<RefSeqGene>> transcriptLinks, Set<RefSeqGene> transcriptsThatAreLinkedTo, Map<String, Integer> coordinatesWithEndOfTranscriptNearMateOfFullPolyARead, String endCoordinatesOutFile, String endCoordinatesOnlyWithPolyAMatesOutFile, String transcriptLinksOutFile, String transcriptsInFile, String transcriptsModifiedOutFile, String transcriptsModifiedOnlyWithPolyAMatesOutFile, String transcriptsModifiedOnlyWithModifiedOutFile) throws IOException {
		FileWriter fwEndCoordinates = new FileWriter(endCoordinatesOutFile);
		for (String completeCoor : allEndCoordinates) {
			String[] completeCoorSplit = completeCoor.split(COLON);
			String chr = completeCoorSplit[0];
			int coor = Integer.parseInt(completeCoorSplit[1]);
			Alignments a = new Alignments(chr, coor, coor + 1);
			BED b = new BED(a);
			fwEndCoordinates.write(b.toString() + NEW_LINE);
		}
		fwEndCoordinates.close();
		
		FileWriter fwEndCoordinatesWithPolyAMates = new FileWriter(endCoordinatesOnlyWithPolyAMatesOutFile);
		for (String completeCoor : coordinatesWithEndOfTranscriptNearMateOfFullPolyARead.keySet()) {
			String[] completeCoorSplit = completeCoor.split(COLON);
			String chr = completeCoorSplit[0];
			int coor = Integer.parseInt(completeCoorSplit[1]);
			Alignments a = new Alignments(chr, coor, coor + 1);
			BED b = new BED(a);
			fwEndCoordinatesWithPolyAMates.write(b.toString() + NEW_LINE);
		}
		fwEndCoordinatesWithPolyAMates.close();
		
		FileWriter fwTranscriptLinks = new FileWriter(transcriptLinksOutFile);
		for (RefSeqGene rsgFrom : transcriptLinks.keySet()) {
			Set<RefSeqGene> rsgTos = transcriptLinks.get(rsgFrom);
			for (RefSeqGene rsgTo : rsgTos) {
				fwTranscriptLinks.write(rsgFrom.toString() + " ---to--- " + rsgTo.toString() + NEW_LINE);
			}
		}
		fwTranscriptLinks.close();
		
		System.out.println("Done writing basic information to files...");
		
		FileWriter fwTranscriptsModified = new FileWriter(transcriptsModifiedOutFile);
		FileWriter fwTranscriptsModifiedOnlyWithPolyAMates = new FileWriter(transcriptsModifiedOnlyWithPolyAMatesOutFile);
		FileWriter fwTranscriptsModifiedOnlyWithModified = new FileWriter(transcriptsModifiedOnlyWithModifiedOutFile);
				
		System.out.println("Starting to modify transcripts file...");
		System.out.println("Number of transcript end transcripts: " + transcriptEndCoors.size());		
		
		List<RefSeqGene> transcriptsList = transcripts.GetGenes();
		for (RefSeqGene transcript : transcriptsList) {
			if (transcriptsThatAreLinkedTo.contains(transcript))
				continue;
			
			RefSeqGene transcriptClone = new RefSeqGene(transcript.toBED(), false);
						
			if (transcriptEndCoors.containsKey(transcript)) {								
				Set<String> endCoordinates = transcriptEndCoors.get(transcript);
				
				boolean wroteAModifiedTranscriptForOnlyMatesOfPolyAFile = false;
				
				for (String completeCoor : endCoordinates) {
					String[] completeCoorSplit = completeCoor.split(COLON);
					String chr = completeCoorSplit[0];
					int coor = Integer.parseInt(completeCoorSplit[1]);
					
					transcriptClone.updateLastExonWithNewEnd(coor);
					String bed = transcriptClone.toBED();
					
					fwTranscriptsModified.write(bed + NEW_LINE);
										
					if (coordinatesWithEndOfTranscriptNearMateOfFullPolyARead.containsKey(completeCoor)) {
						fwTranscriptsModifiedOnlyWithPolyAMates.write(bed + NEW_LINE);
						wroteAModifiedTranscriptForOnlyMatesOfPolyAFile = true;
					}
					
					fwTranscriptsModifiedOnlyWithModified.write(bed + NEW_LINE);
				}
				
				if (!wroteAModifiedTranscriptForOnlyMatesOfPolyAFile)
					fwTranscriptsModifiedOnlyWithPolyAMates.write(transcript.toBED() + NEW_LINE);				
			} else {
				String bed = transcript.toBED();
				fwTranscriptsModified.write(bed + NEW_LINE);
				fwTranscriptsModifiedOnlyWithPolyAMates.write(bed + NEW_LINE);
			}
		}
		fwTranscriptsModified.close();
		fwTranscriptsModifiedOnlyWithPolyAMates.close();
		fwTranscriptsModifiedOnlyWithModified.close();
	}

	public static void filterAlignmentsByOverlapWithLastExonOrAfterLastExon(String alignmentsFile, String transcriptsFile, String filteredAlignmentsFile) throws IOException, ParseException {
		AnnotationReader<? extends GenomicAnnotation> alignments = AnnotationReaderFactory.create(alignmentsFile, BED_FORMAT);
		AnnotationReader<? extends GenomicAnnotation> transcripts = AnnotationReaderFactory.create(transcriptsFile, BED_FORMAT);
		
		alignments.filterByOverlapOfLastExonOrAfterLastExon(transcripts.getAnnotationList(), EXTENSION_FACTOR_AFTER_LAST_EXON_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_LAST_EXON);

		FileWriter fw = new FileWriter(filteredAlignmentsFile);
		Iterator<String> chrIt = alignments.getChromosomeAnnotationMap().keySet().iterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			IntervalTree<? extends GenomicAnnotation> tree = alignments.getChromosomeTree(chr);
			Iterator<? extends GenomicAnnotation> annotIt = tree.valueIterator();
			while(annotIt.hasNext()) {
				LightweightGenomicAnnotation annot = annotIt.next();
				fw.write(annot.toString() + NEW_LINE);
			}
		}
		
		fw.close();
	}
	
	private static void filterAlignmentsByOverlapWithTranscripts(String alignmentsFile, String transcriptsFile, String filteredAlignmentsFile, boolean extendWhenFindingOverlappers, boolean checkOverlapWithOnlyLastExon) throws IOException, ParseException {
		AnnotationReader<? extends GenomicAnnotation> alignments = AnnotationReaderFactory.create(alignmentsFile, BED_FORMAT);
		AnnotationReader<? extends GenomicAnnotation> transcripts = AnnotationReaderFactory.create(transcriptsFile, BED_FORMAT);
		
		int extensionFactor = 0;
		if (extendWhenFindingOverlappers)
			extensionFactor = EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_TRANSCRIPTS;

		if (checkOverlapWithOnlyLastExon)
			alignments.filterByOverlapOfLastExon(transcripts.getAnnotationList(), extensionFactor);
		else
			alignments.filterByOverlap(transcripts.getAnnotationList(), extensionFactor);

		FileWriter fw = new FileWriter(filteredAlignmentsFile);
		Iterator<String> chrIt = alignments.getChromosomeAnnotationMap().keySet().iterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			IntervalTree<? extends GenomicAnnotation> tree = alignments.getChromosomeTree(chr);
			Iterator<? extends GenomicAnnotation> annotIt = tree.valueIterator();
			while(annotIt.hasNext()) {
				LightweightGenomicAnnotation annot = annotIt.next();
				fw.write(annot.toString() + NEW_LINE);
			}
		}
		
		fw.close();
	}
	
	private static String reverseAndChangeBasePairsInSequence(String sequence) {
		StringBuilder newSequence = new StringBuilder(sequence.length());
		for (int i = sequence.length() - 1; i >= 0; i--) {
			char base = sequence.charAt(i);
			char newBase = 0;
			if (base == A)
				newBase = T;
			if (base == T)
				newBase = A;
			if (base == C)
				newBase = G;
			if (base == G)
				newBase = C;
			if (base == N)
				newBase = N;
			newSequence.append(newBase);
		}
		return newSequence.toString();
	}
	
	private static void checkAndProcessIfOrientationOfPolyAMatchesTranscript(Alignment readWithPolyAAlignment, String chr, RefSeqGene transcriptAnnot, int start, int end, int overlapStart, int overlapEnd, boolean transcriptAnnotEndIsOnRight, boolean startsAlignment, Map<String, String> readsWithCorrectOrientation, Map<String, String> readsWithIncorrectOrientation, Map<String, Set<String>> readsWithCoors) {
		if (!(readWithPolyAAlignment.getStart() == start && readWithPolyAAlignment.getEnd() == end))
			return;
		
		boolean isNegativeStrand = readWithPolyAAlignment.isNegativeStrand(); // is true when 16 (reversed) ; b is false when 0 (correct)
		boolean polyAIsOnRight;
		if (startsAlignment) {
			if (isNegativeStrand)
				polyAIsOnRight = true;
			else
				polyAIsOnRight = false;
		} else {
			if (isNegativeStrand)
				polyAIsOnRight = false;
			else
				polyAIsOnRight = true;
		}
		
		int coor;
		if (polyAIsOnRight)
			coor = end;
		else
			coor = start;
		
		if (!(coor >= overlapStart && coor <= overlapEnd))
			return;
		
		String sequence;
		if (isNegativeStrand)
			sequence = reverseAndChangeBasePairsInSequence(readWithPolyAAlignment.getReadSequence());
		else
			sequence = readWithPolyAAlignment.getReadSequence();
		
		String readName = readWithPolyAAlignment.getReadName();
		if (polyAIsOnRight == transcriptAnnotEndIsOnRight)
			readsWithCorrectOrientation.put(readName, sequence);
		else
			readsWithIncorrectOrientation.put(readName, sequence);
		
		String completeCoor = chr + COLON + String.valueOf(coor);
		if (readsWithCoors.containsKey(readName)) {
			readsWithCoors.get(readName).add(completeCoor);
		} else {
			Set<String> s = new HashSet<String>();
			s.add(completeCoor);
			readsWithCoors.put(readName, s);
		}
	}
	
	public static void writePlotOfQualityScoresForPolyABasedOnMatchingOrientationWithTranscript(BEDFileParser transcripts, String readsAnnotationFile, String transcriptsFile, String chrSizesFile, String readsFile, String meanQualityScoresForCorrectOrientationOutFile, String meanQualityScoresForIncorrectOrientationOutFile) throws IOException, ParseException {
		GenericAlignmentDataModel alignmentsStartA = new GenericAlignmentDataModel(READS_WITH_POLYA_START_WITH_A_ALIGNMENT_SORTED_FILE, chrSizesFile, false);
		GenericAlignmentDataModel alignmentsStartT = new GenericAlignmentDataModel(READS_WITH_POLYA_START_WITH_T_ALIGNMENT_SORTED_FILE, chrSizesFile, false);
		GenericAlignmentDataModel alignmentsEndA = new GenericAlignmentDataModel(READS_WITH_POLYA_END_WITH_A_ALIGNMENT_SORTED_FILE, chrSizesFile, false);
		GenericAlignmentDataModel alignmentsEndT = new GenericAlignmentDataModel(READS_WITH_POLYA_END_WITH_T_ALIGNMENT_SORTED_FILE, chrSizesFile, false);

		AnnotationReader<? extends GenomicAnnotation> readsAnnotations = AnnotationReaderFactory.create(readsAnnotationFile, BED_FORMAT);
		System.out.println("READS ANNOTATIONS SIZE: " + readsAnnotations.getAnnotationList().size());
		
		Map<String, String> readsWithIncorrectOrientation = new HashMap<String, String>();
		Map<String, String> readsWithCorrectOrientation = new HashMap<String, String>();
		
		Map<String, Set<String>> readsWithCoors = new HashMap<String, Set<String>>();
		
		System.out.println("Finding transcripts to handle...");
		
		Set<RefSeqGene> transcriptsToHandle = new HashSet<RefSeqGene>();
		List<? extends GenomicAnnotation> readsList = readsAnnotations.getAnnotationList();
		int count = 0;
		for (GenomicAnnotation read : readsList) {
			int[] overlapPoints = findOverlapPointsOfTranscriptsForGivenReadWithPolyA(read);
			int overlapStart = overlapPoints[0];
			int overlapEnd = overlapPoints[1];
			
			IntervalTree<RefSeqGeneWithIsoforms> transcriptsTree = transcripts.getChrTree(read.getChromosome());
			if(transcriptsTree == null) {
				continue;
			}
			
			Iterator<? extends Node<RefSeqGeneWithIsoforms>> transcriptsOverlapperIt = transcriptsTree.overlappers(overlapStart, overlapEnd);
			while (transcriptsOverlapperIt.hasNext()) {
				Node<RefSeqGeneWithIsoforms> transcriptOverlapperNode = transcriptsOverlapperIt.next();
				RefSeqGeneWithIsoforms transcriptsOverlap = transcriptOverlapperNode.getValue();
				Collection<RefSeqGene> transcriptsOverlapCollection = transcriptsOverlap.getAllIsoforms();
				for (RefSeqGene transcriptOverlap : transcriptsOverlapCollection) {
					transcriptsToHandle.add(transcriptOverlap);
				}
			}

			count++;
			if (count % 10000 == 0)
				System.out.println(count + " / " + readsList.size());
		}
		
		System.out.println("Done finding transcripts to handle...");

		count = 0;
		for (RefSeqGene transcriptAnnot : transcriptsToHandle){
			if (transcriptAnnot.getOrientation().equals(ASTERISK))
				continue;
			
			IntervalTree<? extends GenomicAnnotation> readsAnnotationsTree = readsAnnotations.getChromosomeTree(transcriptAnnot.getChr());
			if(readsAnnotationsTree == null) {
				continue;
			}

			Alignments lastExon = transcriptAnnot.getOrientedLastExon();
			int transcriptAnnotEnd = transcriptAnnot.getOrientedEnd();
			
			boolean transcriptAnnotEndIsOnRight;
			int overlapStart;
			int overlapEnd;
			if (lastExon.getStart() == transcriptAnnotEnd) {
				overlapStart = transcriptAnnotEnd - EXTENSION_FACTOR_AFTER_LAST_EXON_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_LAST_EXON;
				overlapEnd = lastExon.getEnd();
				transcriptAnnotEndIsOnRight = false;
			} else {
				overlapStart = lastExon.getStart();
				overlapEnd = transcriptAnnotEnd + EXTENSION_FACTOR_AFTER_LAST_EXON_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_LAST_EXON;
				transcriptAnnotEndIsOnRight = true;
			}

			Iterator<? extends Node<? extends GenomicAnnotation>> readsAnnotationsOverlaperIt = readsAnnotationsTree.overlappers(overlapStart, overlapEnd);
			while(readsAnnotationsOverlaperIt.hasNext()) {
				Node<? extends GenomicAnnotation> readOverlapperNode = readsAnnotationsOverlaperIt.next();
				GenomicAnnotation readAnnot = readOverlapperNode.getValue();
				
				int start = readAnnot.getStart();
				int end = readAnnot.getEnd();
				String chr = readAnnot.getChromosome();
				Alignments a1 = new Alignments(chr, start, end);

				CloseableIterator<Alignment> overlappingAlignmentsStartA = alignmentsStartA.getAlignmentsOverlappingRegion(a1);
				while(overlappingAlignmentsStartA.hasNext()) {
					Alignment readWithPolyAAlignment = overlappingAlignmentsStartA.next();

					checkAndProcessIfOrientationOfPolyAMatchesTranscript(readWithPolyAAlignment, chr, transcriptAnnot, start, end, overlapStart, overlapEnd, transcriptAnnotEndIsOnRight, true, readsWithCorrectOrientation, readsWithIncorrectOrientation, readsWithCoors);
				}
				overlappingAlignmentsStartA.close();
				
				CloseableIterator<Alignment> overlappingAlignmentsStartT = alignmentsStartT.getAlignmentsOverlappingRegion(a1);
				while(overlappingAlignmentsStartT.hasNext()) {
					Alignment readWithPolyAAlignment = overlappingAlignmentsStartT.next();

					checkAndProcessIfOrientationOfPolyAMatchesTranscript(readWithPolyAAlignment, chr, transcriptAnnot, start, end, overlapStart, overlapEnd, transcriptAnnotEndIsOnRight, true, readsWithCorrectOrientation, readsWithIncorrectOrientation, readsWithCoors);
				}
				overlappingAlignmentsStartT.close();
				
				CloseableIterator<Alignment> overlappingAlignmentsEndA = alignmentsEndA.getAlignmentsOverlappingRegion(a1);
				while(overlappingAlignmentsEndA.hasNext()) {
					Alignment readWithPolyAAlignment = overlappingAlignmentsEndA.next();

					checkAndProcessIfOrientationOfPolyAMatchesTranscript(readWithPolyAAlignment, chr, transcriptAnnot, start, end, overlapStart, overlapEnd, transcriptAnnotEndIsOnRight, false, readsWithCorrectOrientation, readsWithIncorrectOrientation, readsWithCoors);
				}
				overlappingAlignmentsEndA.close();
				
				CloseableIterator<Alignment> overlappingAlignmentsEndT = alignmentsEndT.getAlignmentsOverlappingRegion(a1);
				while(overlappingAlignmentsEndT.hasNext()) {
					Alignment readWithPolyAAlignment = overlappingAlignmentsEndT.next();
					
					checkAndProcessIfOrientationOfPolyAMatchesTranscript(readWithPolyAAlignment, chr, transcriptAnnot, start, end, overlapStart, overlapEnd, transcriptAnnotEndIsOnRight, false, readsWithCorrectOrientation, readsWithIncorrectOrientation, readsWithCoors);
				}
				overlappingAlignmentsEndT.close();
			}
			
			count++;
			if (count % 10000 == 0)
				System.out.println(count + " / " + transcriptsToHandle.size());
		}
		System.out.println("count: " + count);
		
		System.out.println("incorrect count: " + readsWithIncorrectOrientation.size());
		System.out.println("correct count: " + readsWithCorrectOrientation.size());

		Map<String, String> readsWithCoorsAsString = new HashMap<String, String>();
		for (String readName : readsWithCoors.keySet()) {
			Set<String> coors = readsWithCoors.get(readName);
			StringBuilder sb = new StringBuilder();
			boolean first = true;
			for (String completeCoor : coors) {
				if (!first)
					sb.append(COMMA);
				sb.append(completeCoor);
				first = false;
			}
			readsWithCoorsAsString.put(readName, sb.toString());
		}
		
		FileWriter fwCorrect = new FileWriter(meanQualityScoresForCorrectOrientationOutFile);
		FileWriter fwIncorrect = new FileWriter(meanQualityScoresForIncorrectOrientationOutFile);

		BufferedReader br = new BufferedReader(new FileReader(readsFile));
		String line;
		while ((line = br.readLine()) != null && line.trim().length() > 0) {
			if (!line.startsWith(AT_SIGN))
				continue;
			
			String name = line;
			String read = br.readLine();
			String description = br.readLine();
			String quality = br.readLine();
			
			String nameWithoutAtSign = name.replaceFirst(AT_SIGN, "");
			
			if (readsWithCorrectOrientation.containsKey(nameWithoutAtSign)) {
				String readWithPolyAAlignment = readsWithCorrectOrientation.get(nameWithoutAtSign);
				if (read.indexOf(readWithPolyAAlignment) >= 0) {
					double mean = calculateMeanPhredQualityScoreOfPolyAForReadWithPolyA(readWithPolyAAlignment, read, quality);
					if (mean > 0)
						fwCorrect.write(readsWithCoorsAsString.get(nameWithoutAtSign) + TAB + String.valueOf(mean) + NEW_LINE);
				}
			} 
			if (readsWithIncorrectOrientation.containsKey(nameWithoutAtSign)){
				String readWithPolyAAlignment = readsWithIncorrectOrientation.get(nameWithoutAtSign);
				if (read.indexOf(readWithPolyAAlignment) >= 0) {
					double mean = calculateMeanPhredQualityScoreOfPolyAForReadWithPolyA(readWithPolyAAlignment, read, quality);
					if (mean > 0)
						fwIncorrect.write(readsWithCoorsAsString.get(nameWithoutAtSign) + TAB + String.valueOf(mean) + NEW_LINE);
				}
			}
		}

		br.close();
		fwCorrect.close();
		fwIncorrect.close();
	}
	
	/**
	public static void writePlotOfQualityScoresForPolyABasedOnMatchingOrientationWithTranscript(BEDFileParser transcripts, String readsAnnotationFile, String transcriptsFile, String chrSizesFile, String readsFile, String meanQualityScoresForCorrectOrientationOutFile, String meanQualityScoresForIncorrectOrientationOutFile) throws IOException, ParseException {
		GenericAlignmentDataModel alignmentsStartA = new GenericAlignmentDataModel(READS_WITH_POLYA_START_WITH_A_ALIGNMENT_SORTED_FILE, chrSizesFile, false);
		GenericAlignmentDataModel alignmentsStartT = new GenericAlignmentDataModel(READS_WITH_POLYA_START_WITH_T_ALIGNMENT_SORTED_FILE, chrSizesFile, false);
		GenericAlignmentDataModel alignmentsEndA = new GenericAlignmentDataModel(READS_WITH_POLYA_END_WITH_A_ALIGNMENT_SORTED_FILE, chrSizesFile, false);
		GenericAlignmentDataModel alignmentsEndT = new GenericAlignmentDataModel(READS_WITH_POLYA_END_WITH_T_ALIGNMENT_SORTED_FILE, chrSizesFile, false);

		AnnotationReader<? extends GenomicAnnotation> readsAnnotations = AnnotationReaderFactory.create(readsAnnotationFile, BED_FORMAT);
		System.out.println("READS ANNOTATIONS SIZE: " + readsAnnotations.getAnnotationList().size());
		
		Map<String, String> readsWithIncorrectOrientation = new HashMap<String, String>();
		Map<String, String> readsWithCorrectOrientation = new HashMap<String, String>();
		
		Map<String, Set<String>> readsWithCoors = new HashMap<String, Set<String>>();
		
		System.out.println("Finding transcripts to handle...");
		
		Set<RefSeqGene> transcriptsToHandle = new HashSet<RefSeqGene>();
		List<? extends GenomicAnnotation> readsList = readsAnnotations.getAnnotationList();
		int count = 0;
		for (GenomicAnnotation read : readsList) {
			int[] overlapPoints = findOverlapPointsOfTranscriptsForGivenReadWithPolyA(read, alignmentsStartA, alignmentsStartT, alignmentsEndA, alignmentsEndT);
			int overlapStart = overlapPoints[0];
			int overlapEnd = overlapPoints[1];
			
			IntervalTree<RefSeqGeneWithIsoforms> transcriptsTree = transcripts.getChrTree(read.getChromosome());
			if(transcriptsTree == null) {
				continue;
			}
			
			Iterator<? extends Node<RefSeqGeneWithIsoforms>> transcriptsOverlapperIt = transcriptsTree.overlappers(overlapStart, overlapEnd);
			while (transcriptsOverlapperIt.hasNext()) {
				Node<RefSeqGeneWithIsoforms> transcriptOverlapperNode = transcriptsOverlapperIt.next();
				RefSeqGeneWithIsoforms transcriptsOverlap = transcriptOverlapperNode.getValue();
				Collection<RefSeqGene> transcriptsOverlapCollection = transcriptsOverlap.getAllIsoforms();
				for (RefSeqGene transcriptOverlap : transcriptsOverlapCollection) {
					transcriptsToHandle.add(transcriptOverlap);
				}
			}

			count++;
			if (count % 10000 == 0)
				System.out.println(count + " / " + readsList.size());
		}
		
		System.out.println("Done finding transcripts to handle...");

		count = 0;
		for (RefSeqGene transcriptAnnot : transcriptsToHandle){
			if (transcriptAnnot.getOrientation().equals(ASTERISK))
				continue;
			
			IntervalTree<? extends GenomicAnnotation> readsAnnotationsTree = readsAnnotations.getChromosomeTree(transcriptAnnot.getChr());
			if(readsAnnotationsTree == null) {
				continue;
			}

			Alignments lastExon = transcriptAnnot.getOrientedLastExon();
			int transcriptAnnotEnd = transcriptAnnot.getOrientedEnd();
			
			boolean transcriptAnnotEndIsOnRight;
			int overlapStart;
			int overlapEnd;
			if (lastExon.getStart() == transcriptAnnotEnd) {
				overlapStart = transcriptAnnotEnd - EXTENSION_FACTOR_AFTER_LAST_EXON_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_LAST_EXON;
				overlapEnd = lastExon.getEnd();
				transcriptAnnotEndIsOnRight = false;
			} else {
				overlapStart = lastExon.getStart();
				overlapEnd = transcriptAnnotEnd + EXTENSION_FACTOR_AFTER_LAST_EXON_FOR_FILTERING_ALIGNMENTS_BY_OVERLAP_WITH_LAST_EXON;
				transcriptAnnotEndIsOnRight = true;
			}

			Iterator<? extends Node<? extends GenomicAnnotation>> readsAnnotationsOverlaperIt = readsAnnotationsTree.overlappers(overlapStart, overlapEnd);
			while(readsAnnotationsOverlaperIt.hasNext()) {
				Node<? extends GenomicAnnotation> readOverlapperNode = readsAnnotationsOverlaperIt.next();
				GenomicAnnotation readAnnot = readOverlapperNode.getValue();
				
				int start = readAnnot.getStart();
				int end = readAnnot.getEnd();
				String chr = readAnnot.getChromosome();
				Alignments a1 = new Alignments(chr, start, end);

				CloseableIterator<Alignment> overlappingAlignmentsStartA = alignmentsStartA.getAlignmentsOverlappingRegion(a1);
				while(overlappingAlignmentsStartA.hasNext()) {
					Alignment readWithPolyAAlignment = overlappingAlignmentsStartA.next();

					checkAndProcessIfOrientationOfPolyAMatchesTranscript(readWithPolyAAlignment, chr, transcriptAnnot, start, end, overlapStart, overlapEnd, transcriptAnnotEndIsOnRight, true, readsWithCorrectOrientation, readsWithIncorrectOrientation, readsWithCoors);
				}
				overlappingAlignmentsStartA.close();
				
				CloseableIterator<Alignment> overlappingAlignmentsStartT = alignmentsStartT.getAlignmentsOverlappingRegion(a1);
				while(overlappingAlignmentsStartT.hasNext()) {
					Alignment readWithPolyAAlignment = overlappingAlignmentsStartT.next();

					checkAndProcessIfOrientationOfPolyAMatchesTranscript(readWithPolyAAlignment, chr, transcriptAnnot, start, end, overlapStart, overlapEnd, transcriptAnnotEndIsOnRight, true, readsWithCorrectOrientation, readsWithIncorrectOrientation, readsWithCoors);
				}
				overlappingAlignmentsStartT.close();
				
				CloseableIterator<Alignment> overlappingAlignmentsEndA = alignmentsEndA.getAlignmentsOverlappingRegion(a1);
				while(overlappingAlignmentsEndA.hasNext()) {
					Alignment readWithPolyAAlignment = overlappingAlignmentsEndA.next();

					checkAndProcessIfOrientationOfPolyAMatchesTranscript(readWithPolyAAlignment, chr, transcriptAnnot, start, end, overlapStart, overlapEnd, transcriptAnnotEndIsOnRight, false, readsWithCorrectOrientation, readsWithIncorrectOrientation, readsWithCoors);
				}
				overlappingAlignmentsEndA.close();
				
				CloseableIterator<Alignment> overlappingAlignmentsEndT = alignmentsEndT.getAlignmentsOverlappingRegion(a1);
				while(overlappingAlignmentsEndT.hasNext()) {
					Alignment readWithPolyAAlignment = overlappingAlignmentsEndT.next();
					
					checkAndProcessIfOrientationOfPolyAMatchesTranscript(readWithPolyAAlignment, chr, transcriptAnnot, start, end, overlapStart, overlapEnd, transcriptAnnotEndIsOnRight, false, readsWithCorrectOrientation, readsWithIncorrectOrientation, readsWithCoors);
				}
				overlappingAlignmentsEndT.close();
			}
			
			count++;
			if (count % 10000 == 0)
				System.out.println(count + " / " + transcriptsToHandle.size());
		}
		System.out.println("count: " + count);
		
		System.out.println("incorrect count: " + readsWithIncorrectOrientation.size());
		System.out.println("correct count: " + readsWithCorrectOrientation.size());

		Map<String, String> readsWithCoorsAsString = new HashMap<String, String>();
		for (String readName : readsWithCoors.keySet()) {
			Set<String> coors = readsWithCoors.get(readName);
			StringBuilder sb = new StringBuilder();
			boolean first = true;
			for (String completeCoor : coors) {
				if (!first)
					sb.append(COMMA);
				sb.append(completeCoor);
				first = false;
			}
			readsWithCoorsAsString.put(readName, sb.toString());
		}
		
		FileWriter fwCorrect = new FileWriter(meanQualityScoresForCorrectOrientationOutFile);
		FileWriter fwIncorrect = new FileWriter(meanQualityScoresForIncorrectOrientationOutFile);

		BufferedReader br = new BufferedReader(new FileReader(readsFile));
		String line;
		while ((line = br.readLine()) != null && line.trim().length() > 0) {
			if (!line.startsWith(AT_SIGN))
				continue;
			
			String name = line;
			String read = br.readLine();
			String description = br.readLine();
			String quality = br.readLine();
			
			String nameWithoutAtSign = name.replaceFirst(AT_SIGN, "");
			
			if (readsWithCorrectOrientation.containsKey(nameWithoutAtSign)) {
				String readWithPolyAAlignment = readsWithCorrectOrientation.get(nameWithoutAtSign);
				if (read.indexOf(readWithPolyAAlignment) >= 0) {
					double mean = calculateMeanPhredQualityScoreOfPolyAForReadWithPolyA(readWithPolyAAlignment, read, quality);
					if (mean > 0)
						fwCorrect.write(readsWithCoorsAsString.get(nameWithoutAtSign) + TAB + String.valueOf(mean) + NEW_LINE);
				}
			} 
			if (readsWithIncorrectOrientation.containsKey(nameWithoutAtSign)){
				String readWithPolyAAlignment = readsWithIncorrectOrientation.get(nameWithoutAtSign);
				if (read.indexOf(readWithPolyAAlignment) >= 0) {
					double mean = calculateMeanPhredQualityScoreOfPolyAForReadWithPolyA(readWithPolyAAlignment, read, quality);
					if (mean > 0)
						fwIncorrect.write(readsWithCoorsAsString.get(nameWithoutAtSign) + TAB + String.valueOf(mean) + NEW_LINE);
				}
			}
		}

		br.close();
		fwCorrect.close();
		fwIncorrect.close();
	}
	*/
	
	private static double calculateMeanPhredQualityScoreOfPolyAForReadWithPolyA(String readWithPolyAAlignment, String read, String quality) {
		String qualityScoresForPolyA;
		if (read.startsWith(readWithPolyAAlignment))
			qualityScoresForPolyA = quality.substring(readWithPolyAAlignment.length());
		else if (read.endsWith(readWithPolyAAlignment))
			qualityScoresForPolyA = quality.substring(0, read.lastIndexOf(readWithPolyAAlignment));
		else
			return -1;

		int qualityScoresTotal = 0;
		for (int i = 0; i < qualityScoresForPolyA.length(); i++) {
			int qualityScore = (int)qualityScoresForPolyA.charAt(i) - 33;
			qualityScoresTotal += qualityScore;
		}
		
		return ((double)qualityScoresTotal / (double)qualityScoresForPolyA.length());
	}
	
	private static Map<Alignment, Short> getAlignmentsOfMatesOfFullPolyAReads() throws IOException{
		Map<Alignment, Short> alignmentsOfMatesOfFullPolyAReads = new HashMap<Alignment, Short>();
		
		AlignmentQueryReader readerA =  new BAMQueryReader(new File( MATES_OF_FULL_POLYA_READS_A_ALIGNMENT_FILE)  );//SamQueryReaderFactory.getReader(new ResourceLocator(MATES_OF_FULL_POLYA_READS_A_ALIGNMENT_FILE), false);
		CloseableIterator<Alignment> readerIterA = readerA.iterator();
		while (readerIterA.hasNext()) {
			Alignment a = readerIterA.next();
			alignmentsOfMatesOfFullPolyAReads.put(a, FULL_POLYA_READ_IS_A_CODE);
		}
		readerIterA.close();
		
		AlignmentQueryReader readerT =  new BAMQueryReader(new File( MATES_OF_FULL_POLYA_READS_T_ALIGNMENT_FILE));// SamQueryReaderFactory.getReader(new ResourceLocator(MATES_OF_FULL_POLYA_READS_T_ALIGNMENT_FILE), false);
		CloseableIterator<Alignment> readerIterT = readerT.iterator();
		while (readerIterT.hasNext()) {
			Alignment a = readerIterT.next();
			alignmentsOfMatesOfFullPolyAReads.put(a, FULL_POLYA_READ_IS_T_CODE);
		}
		readerIterT.close();
		
		return alignmentsOfMatesOfFullPolyAReads;
	}
	
	private static void alignMatesOfFullPolyAReads(String bowtiePath, String referenceIndex) throws IOException {
		String line;
		
		String[] bowtieCmdA = { "/bin/sh" , "-c" , bowtiePath + " -q -k " + String.valueOf(MAX_VALID_ALIGNMENTS_OUTPUT_PER_MATE_OF_FULL_POLYA_READ) + " -S --sam-nohead " + referenceIndex + " " + MATES_OF_FULL_POLYA_READS_A_FILE + " > " + MATES_OF_FULL_POLYA_READS_A_ALIGNMENT_FILE };
		Process bowtieA = Runtime.getRuntime().exec(bowtieCmdA);
		BufferedReader brA = new BufferedReader(new InputStreamReader(bowtieA.getErrorStream()));
		while ((line = brA.readLine()) != null) {
			System.out.println("bowtieA: " + line);
		}
		brA.close();
		
		String[] bowtieCmdT = { "/bin/sh" , "-c" , bowtiePath + " -q -k " + String.valueOf(MAX_VALID_ALIGNMENTS_OUTPUT_PER_MATE_OF_FULL_POLYA_READ) + " -S --sam-nohead " + referenceIndex + " " + MATES_OF_FULL_POLYA_READS_T_FILE + " > " + MATES_OF_FULL_POLYA_READS_T_ALIGNMENT_FILE };
		Process bowtieT = Runtime.getRuntime().exec(bowtieCmdT);
		BufferedReader brT = new BufferedReader(new InputStreamReader(bowtieT.getErrorStream()));
		while ((line = brT.readLine()) != null) {
			System.out.println("bowtieT: " + line);
		}
		brT.close();
	}
	
	private static void findAndWriteMatesOfFullPolyAReads(String readFile) throws IOException {
		FileWriter fwA = new FileWriter(MATES_OF_FULL_POLYA_READS_A_FILE);
		FileWriter fwT = new FileWriter(MATES_OF_FULL_POLYA_READS_T_FILE);
		
		Set<String> namesOfFullPolyAReads_A = new HashSet<String>();
		Set<String> namesOfFullPolyAReads_T = new HashSet<String>();
		Set<Integer> linesWithFullPolyARead = new HashSet<Integer>();
		
		BufferedReader br1 = new BufferedReader(new FileReader(readFile));
		String line;
		int i = 0;
		while ((line = br1.readLine()) != null && line.trim().length() > 0) {
			if (!line.startsWith(AT_SIGN))
				continue;
			
			String name = line;
			String read = br1.readLine();
			br1.readLine();
			br1.readLine();
			
			String nameWithoutPairNum = null;
			if (name.indexOf(SLASH) > 0) {
				String[] nameSplit = name.split(SLASH);
				nameWithoutPairNum = nameSplit[0];
			} else {
				nameWithoutPairNum = name;
			}
			
			if (readIsWithinMinDistanceOfFullPolyARead_A(read)) {
				namesOfFullPolyAReads_A.add(nameWithoutPairNum);
				linesWithFullPolyARead.add(i);
			}
			if (readIsWithinMinDistanceOfFullPolyARead_T(read)) {
				namesOfFullPolyAReads_T.add(nameWithoutPairNum);
				linesWithFullPolyARead.add(i);
			}
				
			i++;
		}
		br1.close();
		
		BufferedReader br2 = new BufferedReader(new FileReader(readFile));
		i = 0;
		while ((line = br2.readLine()) != null && line.trim().length() > 0) {
			if (!line.startsWith(AT_SIGN))
				continue;
			
			String name = line;
			String read = br2.readLine();
			String description = br2.readLine();
			String quality = br2.readLine();
			
			String nameWithoutPairNum = null;
			if (name.indexOf(SLASH) > 0) {
				String[] nameSplit = name.split(SLASH);
				nameWithoutPairNum = nameSplit[0];
			} else {
				nameWithoutPairNum = name;
			}
			
			boolean lineIsFullPolyARead = linesWithFullPolyARead.contains(i);
			
			if (!lineIsFullPolyARead) {
				if (namesOfFullPolyAReads_A.contains(nameWithoutPairNum)) {
					FastqSequence seq = new FastqSequence(name, read, description, quality);
					fwA.write(seq + NEW_LINE);
				}
				if (namesOfFullPolyAReads_T.contains(nameWithoutPairNum)) {
					FastqSequence seq = new FastqSequence(name, read, description, quality);
					fwT.write(seq + NEW_LINE);
				}
			}
				
			i++;
		}
		br2.close();
		
		fwA.close();
		fwT.close();
	}
	
	private static boolean readIsWithinMinDistanceOfFullPolyARead_A(String read) {
		int AMismatches = 0;
		for (int i = 0; i < read.length(); i++) {
			if (read.charAt(i) != A)
				AMismatches++;
			
			if (AMismatches > MAX_MISMATCHES_ALLOWED_FOR_FULL_POLYA_READ)
				return false;
		}
		
		return true;
	}
	
	private static boolean readIsWithinMinDistanceOfFullPolyARead_T(String read) {
		int TMismatches = 0;
		for (int i = 0; i < read.length(); i++) {
			if (read.charAt(i) != T)
				TMismatches++;
			
			if (TMismatches > MAX_MISMATCHES_ALLOWED_FOR_FULL_POLYA_READ)
				return false;
		}
		
		return true;
	}
	
	public static void printLocOfPolyAs(Map<Alignment, Short> alignments) {
		int startsWithA = 0;
		int startsWithT = 0;
		int endsWithA = 0;
		int endsWithT = 0;
		
		for (Alignment a : alignments.keySet()) {
			short loc = alignments.get(a);
			
			if (loc == READ_STARTS_WITH_A_CODE)
				startsWithA++;
			if (loc == READ_STARTS_WITH_T_CODE)
				startsWithT++;
			if (loc == READ_ENDS_WITH_A_CODE)
				endsWithA++;
			if (loc == READ_ENDS_WITH_T_CODE)
				endsWithT++;
		}
		
		System.out.println("Starts with A: " + startsWithA);
		System.out.println("Starts with T: " + startsWithT);
		System.out.println("Ends with A: " + endsWithA);
		System.out.println("Ends with T: " + endsWithT);
	}
	
	public static void writeAlignmentsToBedFile(Map<Alignment, Short> alignments, String bedOutFile) throws IOException {
		FileWriter fw = new FileWriter(bedOutFile);
		
		for (Alignment a : alignments.keySet()) {
			Alignments a1 = new Alignments(a.getChromosome(), a.getAlignmentStart(), a.getAlignmentEnd());
			BED b = new BED(a1);
			fw.write(b.toString() + NEW_LINE);
		}
		
		fw.close();
	}

	public static Map<Alignment, Short> getFilteredAlignmentsOfReadsWithPolyA(Map<Alignment, Short> alignments, String chrSizesFile) throws IOException{
		Map<Alignment, Short> filteredAlignments = new HashMap<Alignment, Short>();
		
		GenericAlignmentDataModel polyAAlignments = new GenericAlignmentDataModel(POLYA_ALIGNMENT_SORTED_FILE, chrSizesFile, false);
		
		for (Alignment a : alignments.keySet()) {
			if (a.getStart() < 0)
				continue;
			
			Alignments a1 = new Alignments(a.getChromosome(), a.getAlignmentStart(), a.getAlignmentEnd());
			
			double numTouchingPolyAInGenome = polyAAlignments.getCountsPerAlignmentForPolyA(a1, EXTENSION_FACTOR_FOR_FILTERING_ALIGNMENTS_WITH_POLYA_OVERLAP);
			if (numTouchingPolyAInGenome > 0)
				continue;
				
			filteredAlignments.put(a, alignments.get(a));
		}
		
		return filteredAlignments;
	}
	
	public static void sortAndIndexAlignments(String igvToolsPath, boolean sortAndIndexReads) throws IOException {
		String line;
		
		if (sortAndIndexReads) {
			String[] sortCmdStartA = { "/bin/sh" , "-c" , igvToolsPath + " sort " + READS_WITH_POLYA_START_WITH_A_ALIGNMENT_FILE + " " + READS_WITH_POLYA_START_WITH_A_ALIGNMENT_SORTED_FILE };
			Process sortStartA = Runtime.getRuntime().exec(sortCmdStartA);
			BufferedReader brSortStartA = new BufferedReader(new InputStreamReader(sortStartA.getErrorStream()));
			while ((line = brSortStartA.readLine()) != null) {
				System.out.println("sortStartA: " + line);
			}
			brSortStartA.close();
			
			String[] sortCmdStartT = { "/bin/sh" , "-c" , igvToolsPath + " sort " + READS_WITH_POLYA_START_WITH_T_ALIGNMENT_FILE + " " + READS_WITH_POLYA_START_WITH_T_ALIGNMENT_SORTED_FILE };
			Process sortStartT = Runtime.getRuntime().exec(sortCmdStartT);
			BufferedReader brSortStartT = new BufferedReader(new InputStreamReader(sortStartT.getErrorStream()));
			while ((line = brSortStartT.readLine()) != null) {
				System.out.println("sortStartT: " + line);
			}
			brSortStartT.close();
			
			String[] sortCmdEndA = { "/bin/sh" , "-c" , igvToolsPath + " sort " + READS_WITH_POLYA_END_WITH_A_ALIGNMENT_FILE + " " + READS_WITH_POLYA_END_WITH_A_ALIGNMENT_SORTED_FILE };
			Process sortEndA = Runtime.getRuntime().exec(sortCmdEndA);
			BufferedReader brSortEndA = new BufferedReader(new InputStreamReader(sortEndA.getErrorStream()));
			while ((line = brSortEndA.readLine()) != null) {
				System.out.println("sortEndA: " + line);
			}
			brSortEndA.close();
			
			String[] sortCmdEndT = { "/bin/sh" , "-c" , igvToolsPath + " sort " + READS_WITH_POLYA_END_WITH_T_ALIGNMENT_FILE + " " + READS_WITH_POLYA_END_WITH_T_ALIGNMENT_SORTED_FILE };
			Process sortEndT = Runtime.getRuntime().exec(sortCmdEndT);
			BufferedReader brSortEndT = new BufferedReader(new InputStreamReader(sortEndT.getErrorStream()));
			while ((line = brSortEndT.readLine()) != null) {
				System.out.println("sortEndT: " + line);
			}
			brSortEndT.close();
			
			String[] indexCmdStartA = { "/bin/sh" , "-c" , igvToolsPath + " index " + READS_WITH_POLYA_START_WITH_A_ALIGNMENT_SORTED_FILE };
			Process indexStartA = Runtime.getRuntime().exec(indexCmdStartA);
			BufferedReader brIndexStartA = new BufferedReader(new InputStreamReader(indexStartA.getErrorStream()));
			while ((line = brIndexStartA.readLine()) != null) {
				System.out.println("indexStartA: " + line);
			}
			brIndexStartA.close();
			
			String[] indexCmdStartT = { "/bin/sh" , "-c" , igvToolsPath + " index " + READS_WITH_POLYA_START_WITH_T_ALIGNMENT_SORTED_FILE };
			Process indexStartT = Runtime.getRuntime().exec(indexCmdStartT);
			BufferedReader brIndexStartT = new BufferedReader(new InputStreamReader(indexStartT.getErrorStream()));
			while ((line = brIndexStartT.readLine()) != null) {
				System.out.println("indexStartT: " + line);
			}
			brIndexStartT.close();
			
			String[] indexCmdEndA = { "/bin/sh" , "-c" , igvToolsPath + " index " + READS_WITH_POLYA_END_WITH_A_ALIGNMENT_SORTED_FILE };
			Process indexEndA = Runtime.getRuntime().exec(indexCmdEndA);
			BufferedReader brIndexEndA = new BufferedReader(new InputStreamReader(indexEndA.getErrorStream()));
			while ((line = brIndexEndA.readLine()) != null) {
				System.out.println("indexEndA: " + line);
			}
			brIndexEndA.close();
			
			String[] indexCmdEndT = { "/bin/sh" , "-c" , igvToolsPath + " index " + READS_WITH_POLYA_END_WITH_T_ALIGNMENT_SORTED_FILE };
			Process indexEndT = Runtime.getRuntime().exec(indexCmdEndT);
			BufferedReader brIndexEndT = new BufferedReader(new InputStreamReader(indexEndT.getErrorStream()));
			while ((line = brIndexEndT.readLine()) != null) {
				System.out.println("indexEndT: " + line);
			}
			brIndexEndT.close();
		}
		
		boolean shouldUsePolyAAlignmentFilteringMethod = false;
		boolean shouldUseSimpleRepeatsFilteringMethod = false;
		if (POLYA_IN_GENOME_FILTERING_METHOD.equals(POLYA_IN_GENOME_FILTERING_METHOD_POLYA))
			shouldUsePolyAAlignmentFilteringMethod = true;
		if (POLYA_IN_GENOME_FILTERING_METHOD.equals(POLYA_IN_GENOME_FILTERING_METHOD_SIMPLEREPEATS))
			shouldUseSimpleRepeatsFilteringMethod = true;
		if (POLYA_IN_GENOME_FILTERING_METHOD.equals(POLYA_IN_GENOME_FILTERING_METHOD_BOTH)) {
			shouldUsePolyAAlignmentFilteringMethod = true;
			shouldUseSimpleRepeatsFilteringMethod = true;
		}
		
		if (shouldUsePolyAAlignmentFilteringMethod) {
			String[] sortCmdPolyA = { "/bin/sh" , "-c" , igvToolsPath + " sort -t /broad/shptmp/hmetsky/ " + POLYA_ALIGNMENT_FILE + " " + POLYA_ALIGNMENT_SORTED_FILE };
			Process sortPolyA = Runtime.getRuntime().exec(sortCmdPolyA);
			BufferedReader brSortPolyA = new BufferedReader(new InputStreamReader(sortPolyA.getErrorStream()));
			while ((line = brSortPolyA.readLine()) != null) {
				System.out.println("sortPolyA: " + line);
			}
			brSortPolyA.close();
			
			String[] indexCmdPolyA = { "/bin/sh" , "-c" , igvToolsPath + " index " + POLYA_ALIGNMENT_SORTED_FILE };
			Process indexPolyA = Runtime.getRuntime().exec(indexCmdPolyA);
			BufferedReader brIndexPolyA = new BufferedReader(new InputStreamReader(indexPolyA.getErrorStream()));
			while ((line = brIndexPolyA.readLine()) != null) {
				System.out.println("indexPolyA: " + line);
			}
			brIndexPolyA.close();
		}
	}
	
	public static void alignPolyA(String bowtiePath, String referenceIndex) throws IOException {
		String polyA = constructPoly(A, LENGTH_OF_POLYA_FOR_ALIGNMENT_OF_POLYA);
		String polyT = constructPoly(T, LENGTH_OF_POLYA_FOR_ALIGNMENT_OF_POLYA);
		String quality = constructPoly(POLYA_QUALITY_VALUE, LENGTH_OF_POLYA_FOR_ALIGNMENT_OF_POLYA);
		String polyAName = "@PolyA_A length=" + String.valueOf(MIN_LENGTH_OF_POLYA_WITHIN_READ);
		String polyTName = "@PolyA_T length=" + String.valueOf(MIN_LENGTH_OF_POLYA_WITHIN_READ);
		String polyADescription = "+PolyA_A length=" + String.valueOf(MIN_LENGTH_OF_POLYA_WITHIN_READ);
		String polyTDescription = "+PolyA_T length=" + String.valueOf(MIN_LENGTH_OF_POLYA_WITHIN_READ);
		String reads = polyAName + NEW_LINE + polyA + NEW_LINE + polyADescription + NEW_LINE + quality + NEW_LINE + polyTName + NEW_LINE + polyT + NEW_LINE + polyTDescription + NEW_LINE + quality;
		FileWriter fw = new FileWriter(POLYA_READS_FILE);
		fw.write(reads);
		fw.close();
		
		String[] bowtieCmd;
		if (ALIGN_POLYA_WITH_V_OPTION)
			bowtieCmd = new String[] { "/bin/sh" , "-c" , bowtiePath + " -q -v " + String.valueOf(MAX_MISMATCHES_ALLOWED_FOR_POLYA_ALIGNMENT_V) + " -a -S --sam-nohead " + referenceIndex + " " + POLYA_READS_FILE + " > " + POLYA_ALIGNMENT_FILE };
		else
			bowtieCmd = new String[] { "/bin/sh" , "-c" , bowtiePath + " -q -n " + String.valueOf(MAX_MISMATCHES_FOR_POLYA_ALIGNMENT_IN_SEED) + " -e " + String.valueOf(MAX_PERMITTED_TOTAL_OF_QUALITY_VALUES) + " -l " + String.valueOf(POLYA_ALIGNMENT_SEED_LENGTH) + " -a -S --sam-nohead " + referenceIndex + " " + POLYA_READS_FILE + " > " + POLYA_ALIGNMENT_FILE };
		Process bowtie = Runtime.getRuntime().exec(bowtieCmd);
		BufferedReader br = new BufferedReader(new InputStreamReader(bowtie.getErrorStream()));
		String line;
		while ((line = br.readLine()) != null) {
			System.out.println("bowtie: " + line);
		}
		br.close();
	}
	
	public static void alignReadsWithPolyA(String bowtiePath, String referenceIndex) throws IOException {
		String line;
		
		String[] bowtieCmdStartA = { "/bin/sh" , "-c" , bowtiePath + " -q -k " + String.valueOf(MAX_VALID_ALIGNMENTS_OUTPUT_PER_READ_WITH_POLYA) + " -S --sam-nohead " + referenceIndex + " " + READS_WITH_POLYA_START_WITH_A_FILE + " > " + READS_WITH_POLYA_START_WITH_A_ALIGNMENT_FILE };
		Process bowtieStartA = Runtime.getRuntime().exec(bowtieCmdStartA);
		BufferedReader brStartA = new BufferedReader(new InputStreamReader(bowtieStartA.getErrorStream()));
		while ((line = brStartA.readLine()) != null) {
			System.out.println("bowtieStartA: " + line);
		}
		brStartA.close();
		
		String[] bowtieCmdStartT = { "/bin/sh" , "-c" , bowtiePath + " -q -k " + String.valueOf(MAX_VALID_ALIGNMENTS_OUTPUT_PER_READ_WITH_POLYA) + " -S --sam-nohead " + referenceIndex + " " + READS_WITH_POLYA_START_WITH_T_FILE + " > " + READS_WITH_POLYA_START_WITH_T_ALIGNMENT_FILE };
		Process bowtieStartT = Runtime.getRuntime().exec(bowtieCmdStartT);
		BufferedReader brStartT = new BufferedReader(new InputStreamReader(bowtieStartT.getErrorStream()));
		while ((line = brStartT.readLine()) != null) {
			System.out.println("bowtieStartT: " + line);
		}
		brStartT.close();
		
		String[] bowtieCmdEndA = { "/bin/sh" , "-c" , bowtiePath + " -q -k " + String.valueOf(MAX_VALID_ALIGNMENTS_OUTPUT_PER_READ_WITH_POLYA) + " -S --sam-nohead " + referenceIndex + " " + READS_WITH_POLYA_END_WITH_A_FILE + " > " + READS_WITH_POLYA_END_WITH_A_ALIGNMENT_FILE };
		Process bowtieEndA = Runtime.getRuntime().exec(bowtieCmdEndA);
		BufferedReader brEndA = new BufferedReader(new InputStreamReader(bowtieEndA.getErrorStream()));
		while ((line = brEndA.readLine()) != null) {
			System.out.println("bowtieEndA: " + line);
		}
		brEndA.close();
		
		String[] bowtieCmdEndT = { "/bin/sh" , "-c" , bowtiePath + " -q -k " + String.valueOf(MAX_VALID_ALIGNMENTS_OUTPUT_PER_READ_WITH_POLYA) + " -S --sam-nohead " + referenceIndex + " " + READS_WITH_POLYA_END_WITH_T_FILE + " > " + READS_WITH_POLYA_END_WITH_T_ALIGNMENT_FILE };
		Process bowtieEndT = Runtime.getRuntime().exec(bowtieCmdEndT);
		BufferedReader brEndT = new BufferedReader(new InputStreamReader(bowtieEndT.getErrorStream()));
		while ((line = brEndT.readLine()) != null) {
			System.out.println("bowtieEndT: " + line);
		}
		brEndT.close();
	}
	
	public static Map<Alignment, Short> getAlignmentsOfReadsWithPolyA() throws IOException {
		Map<Alignment, Short> alignmentsWithPolyA = new HashMap<Alignment, Short>();
		
		AlignmentQueryReader readerStartA = new BAMQueryReader(new File( READS_WITH_POLYA_START_WITH_A_ALIGNMENT_FILE) );// SamQueryReaderFactory.getReader(new ResourceLocator(READS_WITH_POLYA_START_WITH_A_ALIGNMENT_FILE), false);
		CloseableIterator<Alignment> readerIterStartA = readerStartA.iterator();
		while (readerIterStartA.hasNext()) {
			Alignment a = readerIterStartA.next();
			alignmentsWithPolyA.put(a, READ_STARTS_WITH_A_CODE);
		}
		readerIterStartA.close();
		
		AlignmentQueryReader readerStartT =  new BAMQueryReader(new File( READS_WITH_POLYA_START_WITH_T_ALIGNMENT_FILE));// SamQueryReaderFactory.getReader(new ResourceLocator(READS_WITH_POLYA_START_WITH_T_ALIGNMENT_FILE), false);
		CloseableIterator<Alignment> readerIterStartT = readerStartT.iterator();
		while (readerIterStartT.hasNext()) {
			Alignment a = readerIterStartT.next();
			alignmentsWithPolyA.put(a, READ_STARTS_WITH_T_CODE);
		}
		readerStartT.close();
		
		AlignmentQueryReader readerEndA = new BAMQueryReader(new File( READS_WITH_POLYA_END_WITH_A_ALIGNMENT_FILE));//SamQueryReaderFactory.getReader(new ResourceLocator(READS_WITH_POLYA_END_WITH_A_ALIGNMENT_FILE), false);
		CloseableIterator<Alignment> readerIterEndA = readerEndA.iterator();
		while (readerIterEndA.hasNext()) {
			Alignment a = readerIterEndA.next();
			alignmentsWithPolyA.put(a, READ_ENDS_WITH_A_CODE);
		}
		readerEndA.close();
		
		AlignmentQueryReader readerEndT = new BAMQueryReader(new File( READS_WITH_POLYA_END_WITH_T_ALIGNMENT_FILE));//SamQueryReaderFactory.getReader(new ResourceLocator(READS_WITH_POLYA_END_WITH_T_ALIGNMENT_FILE), false);
		CloseableIterator<Alignment> readerIterEndT = readerEndT.iterator();
		while (readerIterEndT.hasNext()) {
			Alignment a = readerIterEndT.next();
			alignmentsWithPolyA.put(a, READ_ENDS_WITH_T_CODE);
		}
		readerEndT.close();
		
		return alignmentsWithPolyA;
	}
	
	public static void findAndWriteReadsWithPolyA(String readFile) throws IOException {
		FileWriter fwStartA = new FileWriter(READS_WITH_POLYA_START_WITH_A_FILE);
		FileWriter fwStartT = new FileWriter(READS_WITH_POLYA_START_WITH_T_FILE);
		FileWriter fwEndA = new FileWriter(READS_WITH_POLYA_END_WITH_A_FILE);
		FileWriter fwEndT = new FileWriter(READS_WITH_POLYA_END_WITH_T_FILE);
		
		BufferedReader br = new BufferedReader(new FileReader(readFile));
		String line;
		int readsWithPolyAWritten = 0;
		while ((line = br.readLine()) != null && line.trim().length() > 0) {
			if (!line.startsWith(AT_SIGN))
				continue;
			
			String name = line;
			String read = br.readLine();
			String description = br.readLine();
			String quality = br.readLine();
			
			String AStartRemainder = stripPolyAFromReadStart(read, A);
			String TStartRemainder = stripPolyAFromReadStart(read, T);
			String AEndRemainder = stripPolyAFromReadEnd(read, A);
			String TEndRemainder = stripPolyAFromReadEnd(read, T);
			int numTrue = 0;
			if (AStartRemainder != null)
				numTrue++;
			if (TStartRemainder != null)
				numTrue++;
			if (AEndRemainder != null)
				numTrue++;
			if (TEndRemainder != null)
				numTrue++;
			
			if (numTrue != 1)
				continue;
			
			if (AStartRemainder != null) {
				FastqSequence seq = new FastqSequence(name, AStartRemainder, description, quality);
				fwStartA.write(seq + NEW_LINE);
			}
			if (TStartRemainder != null) {
				FastqSequence seq = new FastqSequence(name, TStartRemainder, description, quality);
				fwStartT.write(seq + NEW_LINE);
			}
			if (AEndRemainder != null) {
				FastqSequence seq = new FastqSequence(name, description, AEndRemainder, quality);
				fwEndA.write(seq + NEW_LINE);
			}
			if (TEndRemainder != null) {
				FastqSequence seq = new FastqSequence(name, description, TEndRemainder, quality);
				fwEndT.write(seq + NEW_LINE);
			}
			readsWithPolyAWritten++;
			
			if (readsWithPolyAWritten % 1000 == 0) {
				fwStartA.flush();
				fwStartT.flush();
				fwEndA.flush();
				fwEndT.flush();
			}
		}
		
		br.close();
		fwStartA.close();
		fwStartT.close();
		fwEndA.close();
		fwEndT.close();
	}
	
	private static String stripPolyAFromReadStart(String read, char c) {
		int mismatches = 0;
		int j = read.length();
		for (int i = 0; i < read.length(); i++) {
			if (read.charAt(i) != c)
				mismatches++;
			
			if (mismatches > MAX_MISMATCHES_ALLOWED_FOR_POLYA_WITHIN_READ) {
				if (i + 1 < MIN_LENGTH_OF_POLYA_WITHIN_READ)
					j = -1;
				else
					j = i;
				break;
			}
		}
		
		if (j >= 0 && read.length() - j >= MIN_LENGTH_OF_READ_REMAINDER_FOR_POLYA_WITHIN_READ)
			return read.substring(j);
		else
			return null;
	}
	
	private static String stripPolyAFromReadEnd(String read, char c) {
		int mismatches = 0;
		int j = 0;
		for (int i = read.length() - 1; i >= 0; i--) {
			if (read.charAt(i) != c)
				mismatches++;
			
			if (mismatches > MAX_MISMATCHES_ALLOWED_FOR_POLYA_WITHIN_READ) {
				if (read.length() - i - 1 < MIN_LENGTH_OF_POLYA_WITHIN_READ)
					j = -1;
				else
					j = i;
				break;
			}
		}

		if (j >= 0 && j + 1 >= MIN_LENGTH_OF_READ_REMAINDER_FOR_POLYA_WITHIN_READ)
			return read.substring(0, j + 1);
		else
			return null;
	}
	
	private static String constructPoly(char letter, int n) {
		StringBuilder sb = new StringBuilder(n);
		for (int i = 0; i < n; i++) {
			sb.append(letter);
		}
		return sb.toString();
	}
	
}