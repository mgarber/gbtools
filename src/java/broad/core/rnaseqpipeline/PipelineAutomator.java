package broad.core.rnaseqpipeline;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.math.MathException;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.error.ParseException;
import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;
import broad.core.util.PipelineUtils;
import broad.pda.countreads.DuplicatesAndLibrarySizeAnalysis;
import broad.pda.countreads.LibraryCompositionByRnaClass;
import broad.pda.geneexpression.dge.DGE;

/**
 * This class will use an input list of files and a config file to run an automated pipeline to align and process sequencing data
 * @author skadri 
 * @author prussell
 *
 */
public class PipelineAutomator {

	static Logger logger = Logger.getLogger(PipelineAutomator.class.getName());

	/*
	 * For a file with paired data in the fq list file, there will be 4 columns
	 * Name	left.fq	right.fq	condition
	 * Name is unique
	 * Condition represents a string describing the sample/expt. Replicates have the same string in the condition
	 */
	private static int PAIRED_COLUMNS = 4;
	private static int UNPAIRED_COLUMNS = PAIRED_COLUMNS - 1;
	
	PipelineAutomatorConfigFileParser configP;
	
	/**
	 * The sample names
	 */
	TreeSet<String> sampleNames;
	
	/*
	 * Mapping of name to FastQ file paths
	 */
	Map<String,String> leftFqs;
	Map<String,String> rightFqs;
	/*
	 * Mapping of name to whether there is paired data
	 */
	Map<String,Boolean> pairedData;
	
	/*
	 * Mappings of names to CURRENT fqs being used. 
	 * This is because if the user filters out rRNA, the name of the fqs to align changes.
	 */
	Map<String,String> currentLeftFqs;
	Map<String,String> currentRightFqs;

	/**
	 * The paths of the current bam files being used
	 * This is because we might run extra alignments and merge bam files
	 */
	Map<String,String> currentBamFiles;
	
	/**
	 * The directory containing the current bam files
	 */
	String currentBamDir;
	
	Map<String,String> nameToCondition;

	String queueName;
	/*
	 * LIST OF COMMANDS (in order for ease of coding)
	 */
	static String SPLIT_TRIM_BARCODES = "SPLIT_TRIM_BARCODES";
	static String LIBRARY_STATS = "LIBRARY_STATS";
	static String RNA_CLASSES = "RNA_CLASSES";
	static String FILTER_RRNA = "FILTER_RRNA";
	static String ALIGN = "ALIGN";
	static String RUN_DGE = "RUN_DGE";
	
	/**
	 * Output directories
	 */
	static String TOPHAT_DIRECTORY= "tophat_to_genome";
	static String NOVOALIGN_DIRECTORY = "novoalign_unmapped_to_genome";
	static String MERGED_TOPHAT_NOVOALIGN_DIRECTORY = "merged_alignments_tophat_novoalign";
	static String LIBRARY_STATS_DIRECTORY = "library_stats";
	static String FILTER_RRNA_DIRECTORY = "filter_rRNA";


	
	/**
	 * @param inputListFile The input list of fastq files
	 * @param configFileName The config file
	 * @throws IOException
	 * @throws ParseException
	 * @throws MathException
	 * @throws InterruptedException
	 */
	public PipelineAutomator(String inputListFile,String configFileName) throws IOException, ParseException, MathException, InterruptedException {
		
		Globals.setHeadless(true);
		logger.info("Using Version R4.4");
		logger.debug("DEBUG ON");

		this.configP = new PipelineAutomatorConfigFileParser(configFileName);
		
		//If this flag is false, then the expected input file formats are different
		boolean runBasic = false;
		/*
		 * Format of the inputListFile can depend on the first task in question
		 */
		// IF THE CONFIG FILE CONTAINS A BASIC COMMAND
		if(this.containsBasicCommand()){
			
			runBasic = true;
			
			//Initialize the FQ file lists
			this.sampleNames = new TreeSet<String>();
			this.leftFqs = new TreeMap<String,String>();
			this.rightFqs = new TreeMap<String,String>();
			this.currentLeftFqs = new TreeMap<String,String>();
			this.currentRightFqs = new TreeMap<String,String>();
			this.pairedData = new TreeMap<String,Boolean>();
			this.currentBamFiles = new TreeMap<String,String>();
			this.nameToCondition = new TreeMap<String,String>();
			this.queueName = configP.basicOptions.getQueueName();
			
			
			//Read the Fq list
			readFqList(inputListFile);

			// Split and trim barcodes
			if(this.configP.hasCommandFor(SPLIT_TRIM_BARCODES)){
				splitTrimBarcodes();
			}
			
			// Count reads, quantify duplicates, estimate library size
			if(this.configP.hasCommandFor(LIBRARY_STATS)){
				calculateLibraryStats();
			}
			
			// Characterize library composition by RNA class
			if(this.configP.hasCommandFor(RNA_CLASSES)){
				quantifyRNAClasses();
			}
			
			// Update current fastqs with reads not matching rRNA
			if(this.configP.hasCommandFor(FILTER_RRNA)){
				filterrRNA();
			}
			
			// Align to genome; generate bam and tdf files
			if(this.configP.hasCommandFor(ALIGN)){
				alignToGenome();
			}
			
		}
		
		if(this.configP.hasCommandFor(RUN_DGE)){
			
			String DGEInputFile= null;
			if(!this.configP.hasCommandFor(ALIGN)){
				//Input file format
				//Name	Bam_file	condition
				//TODO: Not needed right now but provision for future flexibility
				DGEInputFile = inputListFile;
				//Make the DGE command using all the options
				
			}
			else{
				/*
				 * Create my own DGE input file using the output of above
				 */
				DGEInputFile = createDGEInputFile(inputListFile);
			}
			
			runDGE(runBasic,DGEInputFile);
			
		}
	}
	
	/**
	 * Returns true if the config file contains one or more basic commands
	 * @return
	 */
	private boolean containsBasicCommand(){
		return (this.configP.hasCommandFor(SPLIT_TRIM_BARCODES) || this.configP.hasCommandFor(LIBRARY_STATS) 
				|| this.configP.hasCommandFor(RNA_CLASSES) || this.configP.hasCommandFor(FILTER_RRNA) || this.configP.hasCommandFor(ALIGN));
	}
	
	/**
	 * 
	 * @param fileName
	 * @throws IOException 
	 */
	private void readFqList(String fileName) throws IOException{
		
		/*	Input file is a list
		*	Name		Left.fq		Right.fq	condition
		*	If unpaired, only Name	Left.fq	condition
		*	The NAME must be unique
		*/
		//TODO: If the name is not unique, assign a unique name
		logger.info("");
		logger.info("Reading list of fastq files from " + fileName + "...");
		
		FileReader r = new FileReader(fileName);
		BufferedReader b = new BufferedReader(r);
		if(!b.ready()) {
			throw new IllegalArgumentException("File " + fileName + " is empty.");
		}
		StringParser p = new StringParser();
		
		boolean first = true;
		boolean unpaired = true;
		while(b.ready()) {
			String line = b.readLine();
			p.parse(line);
			int cols = p.getFieldCount();
			if(cols == 0) continue;
			
			if(first) {
				if(cols < UNPAIRED_COLUMNS) {
					throw new IllegalArgumentException("Illegal number of columns in " + fileName);
				}
				if(cols == PAIRED_COLUMNS) unpaired = false;
				first = false;
			}
			
			String sampleName = p.asString(0);
			this.sampleNames.add(sampleName);
			String leftReads = p.asString(1);
			this.leftFqs.put(sampleName, leftReads);
			
			if(unpaired) {
				if(cols != UNPAIRED_COLUMNS) {
					throw new IllegalArgumentException("Illegal number of columns in " + fileName);
				}
				this.pairedData.put(sampleName, Boolean.valueOf(false));
				this.nameToCondition.put(sampleName, p.asString(2));
			} else {
				if(cols != PAIRED_COLUMNS) {
					throw new IllegalArgumentException("Illegal number of columns in " + fileName);
				}				
				this.rightFqs.put(sampleName, p.asString(2));
				this.pairedData.put(sampleName, Boolean.valueOf(true));
				this.nameToCondition.put(sampleName, p.asString(3));
			}
			
		}
		
		b.close();
		r.close();
		
		if(first) {
			throw new IllegalArgumentException("File " + fileName + " is invalid.");
		}
		
		this.currentLeftFqs = this.leftFqs;
		this.currentRightFqs = this.rightFqs;
		
		int pa = this.currentRightFqs.isEmpty() ? 0 : this.currentRightFqs.keySet().size();
		int u = this.currentLeftFqs.keySet().size() - pa;
		logger.info("There are " + pa + " sets of paired end reads and " + u + " sets of single end reads.");
		
	}
	
	/**
	 * TASK 1: SPLIT_TRIM_BARCODES
	 */
	private void splitTrimBarcodes(){
		// TODO: finish
	}
	
	/**
	 * TASK 2: CALCULATE LIBRARY STATS
	 * For each library, count total reads, unique reads, percent duplicates, and estimated library size
	 * @author prussell
	 * @throws IOException 
	 */
	private void calculateLibraryStats() throws IOException {
		
		logger.info("");
		logger.info("Calculating library stats...");
		
		// Make directory for library stats
		File dir = new File(LIBRARY_STATS_DIRECTORY);
		boolean madeDir = dir.mkdir();
		if(!dir.exists()) {
			throw new IOException("Could not create directory " + LIBRARY_STATS_DIRECTORY);
		}
		String output = LIBRARY_STATS_DIRECTORY + "/" + "library_stats.out";
		File outputFile = new File(output);
		if(outputFile.exists()) {
			logger.info("WARNING: library stats file " + output + " exists. Not calculating new library stats.");
			return;
		}
		
		logger.info("Writing library stats to file " + output + ".");
		FileWriter w = new FileWriter(output);
		
		// Get stats for each sample
		boolean first = true;
		for(String sample : this.leftFqs.keySet()) {
			logger.info("Calculating library stats for sample " + sample + "...");
			if(first && this.pairedData.get(sample).booleanValue()) {
				w.write("Sample\tTotal_read_pairs\tUnique_read_pairs\tPercent_duplicated\tEst_library_size\n");
			}
			if(first && !this.pairedData.get(sample).booleanValue()) {
				w.write("Sample\tTotal_reads\tUnique_reads\tPercent_duplicated\tEst_library_size\n");
			}
			first = false;
			DuplicatesAndLibrarySizeAnalysis d = new DuplicatesAndLibrarySizeAnalysis(this.leftFqs.get(sample), this.pairedData.get(sample).booleanValue() ? this.rightFqs.get(sample) : null, logger);
			w.write(sample + "\t" + d.getTotalReads() + "\t" + d.getNumUniqueReads() + "\t" + d.getPercentDuplicated() + "\t" + d.getEstimatedLibrarySize() + "\n");
		}
		w.close();
		logger.info("");
		logger.info("Done calculating library stats.");
	}
	
	/**
	 * TASK 3: CALCULATE RNA CLASSES
	 * Quantify percentage of reads originating from each RNA class
	 * @author prussell
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	private void quantifyRNAClasses() throws IOException, InterruptedException{
		
		logger.info("");
		logger.info("Quantifying RNA classes...");
		
		// Make directory for RNA classes
		String rnaClassDir = "rna_class_counts";
		File dir = new File(rnaClassDir);
		boolean madeDir = dir.mkdir();
		if(!dir.exists()) {
			throw new IOException("Could not create directory " + rnaClassDir);
		}
		
		// Get options from config file
		Map<String, String> classFiles = this.configP.basicOptions.getRNAClassFileNames();
		String genomeBowtieIndex = this.configP.basicOptions.getGenomeBowtieIndex();
		String bowtie2Executable = this.configP.basicOptions.getBowtie2ExecutablePath();
		String bowtie2BuildExecutable = this.configP.basicOptions.getBowtie2BuildExecutablePath();
		String samtoolsExecutable = this.configP.basicOptions.getSamtoolsPath();
		
		// Align and count reads for each library and RNA class
		LibraryCompositionByRnaClass lcrc = new LibraryCompositionByRnaClass(genomeBowtieIndex, classFiles, this.leftFqs, this.rightFqs, logger);
		Map<String, Integer> totalReadCounts = lcrc.getTotalReadCounts();
		Map<String, Map<String, Integer>> classCounts = lcrc.alignAndGetCounts(samtoolsExecutable, bowtie2Executable, bowtie2BuildExecutable, rnaClassDir);
		
		// Write tables to file
		String countFile = rnaClassDir + "/counts_by_class.out";
		String pctFile = rnaClassDir + "/percentages_by_class.out";
		logger.info("Writing table of counts to file " + countFile);
		logger.info("Writing table of percentages to file " + pctFile);
		FileWriter countWriter = new FileWriter(countFile);
		FileWriter pctWriter = new FileWriter(pctFile);
		
		// Get set of class names
		Iterator<String> iter = classCounts.keySet().iterator();
		Set<String> classNames = classCounts.get(iter.next()).keySet();
		String header = "Sample\t";
		for(String className : classNames) {
			header += className;
			header += "\t";
		}
		header += "\n";
		countWriter.write(header);
		pctWriter.write(header);
		
		// Get counts and make percentages
		for(String sample : classCounts.keySet()) {
			String countLine = sample + "\t";
			String pctLine = sample + "\t";
			for(String className : classNames) {
				countLine += classCounts.get(sample).get(className).toString();
				String pct = Double.valueOf((double)classCounts.get(sample).get(className).intValue() / (double)totalReadCounts.get(sample).intValue()).toString();
				pctLine += pct;
				countLine += "\t";
				pctLine += "\t";
			}
			countLine += "\n";
			pctLine += "\n";
			countWriter.write(countLine);
			pctWriter.write(pctLine);
		}

		countWriter.close();
		pctWriter.close();
		
		logger.info("");
		logger.info("Done quantifying RNA classes.");
				
	}
	
	/**
	 * TASK 4: FILTER OUT RRNA
	 * Align reads to rRNA, retain reads that do not align, and replace current fastq files with files of non aligning reads
	 * @author prussell
	 * @throws InterruptedException 
	 * @throws IOException 
	 */
	private void filterrRNA() throws IOException, InterruptedException{
		
		// Establish paths and software locations
		File outDirFile = new File(FILTER_RRNA_DIRECTORY);
		boolean madeDir = outDirFile.mkdir();
		String outIndex = FILTER_RRNA_DIRECTORY + "/rRNA";
		String rRnaFasta = this.configP.basicOptions.getRrnaSeqsForDepletion();
		String bowtieBuild = this.configP.basicOptions.getBowtie2BuildExecutablePath();
		String bowtie = this.configP.basicOptions.getBowtie2ExecutablePath();
		
		logger.info("");
		logger.info("Filtering ribosomal RNA by removing reads that map to sequences in " + rRnaFasta);
		
		// Make bowtie2 index for ribosomal RNA
		AlignmentUtils.makeBowtie2Index(rRnaFasta, outIndex, bowtieBuild, FILTER_RRNA_DIRECTORY);
		
		// Establish output file names
		Map<String, String> outRibosomalMap = new TreeMap<String, String>();
		Map<String, String> outFilteredUnpairedMap = new TreeMap<String, String>();
		Map<String, String> outFilteredPairedArgMap = new TreeMap<String, String>();
		Map<String, String> outFilteredPaired1Map = new TreeMap<String, String>();
		Map<String, String> outFilteredPaired2Map = new TreeMap<String, String>();
		ArrayList<String> jobIDs = new ArrayList<String>();
		
		// Align each sample to rRNA
		for(String sample : this.sampleNames) {

			String outRibosomal = FILTER_RRNA_DIRECTORY + "/" + sample + "_ribosomal_mappings.sam";
			outRibosomalMap.put(sample, outRibosomal);
			String outFilteredUnpaired = FILTER_RRNA_DIRECTORY + "/" + sample + "_filtered_rRNA.fq";
			outFilteredUnpairedMap.put(sample, outFilteredUnpaired);
			String outFilteredPairedArg = FILTER_RRNA_DIRECTORY + "/" + sample + "_filtered_rRNA_%.fq";
			outFilteredPairedArgMap.put(sample, outFilteredPairedArg);
			String outFilteredPaired1 = FILTER_RRNA_DIRECTORY + "/" + sample + "_filtered_rRNA_1.fq";
			outFilteredPaired1Map.put(sample,outFilteredPaired1);
			String outFilteredPaired2 = FILTER_RRNA_DIRECTORY + "/" + sample + "_filtered_rRNA_2.fq";
			outFilteredPaired2Map.put(sample,outFilteredPaired2);
			
			boolean paired = this.pairedData.get(sample).booleanValue();
			
			// check if filtered fastq files exist
			if(!paired) {
				File filteredFile = new File(outFilteredUnpaired);
				if(filteredFile.exists()) {
					logger.info("WARNING: filtered file for sample " + sample + " already exists: " + outFilteredUnpaired + ". Not creating new file.");
					continue;
				}
			} else {
				logger.info("Filtering sample " + sample + "...");
				File filteredFile1 = new File(outFilteredPaired1);
				File filteredFile2 = new File(outFilteredPaired2);
				if(filteredFile1.exists() && filteredFile2.exists()) {
					logger.info("WARNING: filtered files for sample " + sample + " already exist: " + outFilteredPaired1 + " and " + outFilteredPaired2 + ". Not creating new files.");
					continue;
				}
				
			}
			
			// Align to rRNA and keep unmapped reads
			String jobID = AlignmentUtils.runBowtie(outIndex, this.currentLeftFqs.get(sample), paired ? this.currentRightFqs.get(sample) : null, outRibosomal, 
					outFilteredUnpaired, bowtie, AlignmentUtils.MAX_INSERT_SIZE, FILTER_RRNA_DIRECTORY, paired);
			jobIDs.add(jobID);
			
		}
		
		PipelineUtils.waitForAllJobs(jobIDs, Runtime.getRuntime());
		
		logger.info("");
		logger.info("Done aligning to rRNAs. Updating current fastq files to filtered files.");
		
		// Update current fastq files to the ribosome filtered files
		for(String sample : this.leftFqs.keySet()) {
			if(this.pairedData.get(sample).booleanValue()) {
				this.currentLeftFqs.put(sample, outFilteredPaired1Map.get(sample));
				this.currentRightFqs.put(sample, outFilteredPaired2Map.get(sample));
				logger.info("Current fastq files for sample " + sample + " are " + this.currentLeftFqs.get(sample) + " and " + this.currentRightFqs.get(sample) + ".");
			} else {
				this.currentLeftFqs.put(sample, outFilteredUnpairedMap.get(sample));
				logger.info("Current fastq file for sample " + sample + " is " + this.currentLeftFqs.get(sample) + ".");
			}
		}
		
		logger.info("");
		logger.info("Done filtering rRNA and updating current fastq files.");
	}
	
	/**
	 * TASK 5: ALIGN TO GENOME
	 * 1. Align to genome with tophat
	 * 2. If specified, use Novoalign to align unmapped reads and merge Tophat+Novoalign bam files
	 * 3. Report alignment statistics
	 * 4. Index bam files
	 * 5. Make tdf files
	 * 6. Run Picard metrics
	 * @author prussell
	 * @throws InterruptedException 
	 * @throws IOException 
	 */
	private void alignToGenome() throws IOException, InterruptedException{
		
		logger.info("");
		logger.info("Aligning to genome...\n");
		
		// Establish index files and executables
		logger.info("Aligning current fastq files to genome using Tophat...");
		String genomeIndex = this.configP.basicOptions.getGenomeBowtieIndex();
		String samtools = this.configP.basicOptions.getSamtoolsPath();
		String tophat = this.configP.basicOptions.getTophatPath();
		File tophatDirFile = new File(TOPHAT_DIRECTORY);
		tophatDirFile.mkdir();
		// Get tophat options
		Map<String, String> tophatOptions = this.configP.basicOptions.getTophatOptions();
		// Use a different output directory for each sample
		if(tophatOptions.containsKey("-o") || tophatOptions.containsKey("--output-dir")) {
			logger.info("WARNING: Overriding tophat output directory provided in config file. Creating directories for each sample.");
			tophatOptions.remove("-o");
			tophatOptions.remove("--output-dir");
		}
		String optionsString = "";
		for(String flag : tophatOptions.keySet()) {
			optionsString += flag + " " + tophatOptions.get(flag) + " ";
		}
		
		// Establish output directories
		Map<String, String> tophatDirsPerSample = new TreeMap<String, String>();
		Map<String, String> tophatSubdirsPerSample = new TreeMap<String, String>();
		Map<String, String> tophatBamFinalPath = new TreeMap<String, String>();  // where the bam file will be moved to
		ArrayList<String> tophatJobIDs = new ArrayList<String>();
		// Run tophat
		for(String sample : this.sampleNames) {
			tophatDirsPerSample.put(sample, TOPHAT_DIRECTORY + "/" + TOPHAT_DIRECTORY + "_" + sample);
			tophatSubdirsPerSample.put(sample, TOPHAT_DIRECTORY + "_" + sample);
			File outdir = new File(tophatDirsPerSample.get(sample));
			outdir.mkdir();
			tophatBamFinalPath.put(sample, TOPHAT_DIRECTORY + "/" + sample + ".bam");
			File finalBamFile = new File(tophatBamFinalPath.get(sample));
			// Check if tophat has already been run
			if(finalBamFile.exists()) {
				logger.info("WARNING: alignment file " +tophatBamFinalPath.get(sample) + " already exists. Not rerunning alignment.");
				continue;
			}
			// Make tophat command
			String tophatCmmd = tophat + " ";
			tophatCmmd += optionsString + " ";
			tophatCmmd += "--output-dir " + outdir + " "; 
			tophatCmmd += genomeIndex + " ";
			tophatCmmd += this.currentLeftFqs.get(sample) + " ";
			if(this.pairedData.get(sample).booleanValue()) tophatCmmd += this.currentRightFqs.get(sample);
			logger.info("Writing tophat output for sample " + sample + " to directory " + outdir + ".");
			logger.info("Running tophat command: " + tophatCmmd);
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			tophatJobIDs.add(jobID);
			logger.info("LSF job ID is " + jobID + ".");
			// Submit job
			//PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, tophatCmmd, outdir + "/tophat_" + jobID + ".bsub", queueName, 16);
			PipelineUtils.bsubMultiProcess(Runtime.getRuntime(), jobID, tophatCmmd, outdir + "/tophat_" + jobID + ".bsub", queueName);
		}
		// Wait for tophat jobs to finish
		logger.info("Waiting for tophat jobs to finish...");
		PipelineUtils.waitForAllJobs(tophatJobIDs, Runtime.getRuntime());
		logger.info("All samples done aligning to genome.");
		
		// Move all tophat bam files to one directory
		logger.info("");
		logger.info("Moving all tophat bam files to directory " + TOPHAT_DIRECTORY + "...");
		for(String sample : this.sampleNames) {
			String cmmd = "mv "+TOPHAT_DIRECTORY+ "/" + tophatSubdirsPerSample.get(sample) + "/accepted_hits.bam " + tophatBamFinalPath.get(sample);
			//String cmmd = "mv " + tophatSubdirsPerSample.get(sample) + "/accepted_hits.bam " + tophatBamFinalPath.get(sample);
			Process p = Runtime.getRuntime().exec(cmmd, null);
			p.waitFor();
			// Update current bam files
			this.currentBamFiles.put(sample, tophatBamFinalPath.get(sample));
		}
		// Update current bam directory
		this.currentBamDir = TOPHAT_DIRECTORY;
		
		// Count aligned and unaligned reads
/*		logger.info("");
		logger.info("Counting mapped and unmapped reads...");
		String countFile = TOPHAT_DIRECTORY + "/mapped_unmapped_count.out";
		String pctFile = TOPHAT_DIRECTORY + "/mapped_unmapped_percentage.out";
		FileWriter w = new FileWriter(countFile);
		FileWriter wp = new FileWriter(pctFile);
		String header = "Sample\tMapped\tUnmapped\n";
		w.write(header);
		wp.write(header);
		for(String sample : this.sampleNames) {
			int mapped = AlignmentUtils.countAlignments(samtools, TOPHAT_DIRECTORY + "/" + sample + ".bam", tophatDirsPerSample.get(sample), false);
			int unmapped = AlignmentUtils.countAlignments(samtools, tophatDirsPerSample.get(sample) + "/unmapped.bam", tophatDirsPerSample.get(sample), true, false);
			int total = mapped + unmapped;
			w.write(sample + "\t" + mapped + "\t" + unmapped + "\n");
			wp.write(sample + "\t" + (double)mapped/(double)total + "\t" + (double)unmapped/(double)total + "\n");
		}
		w.close();
		wp.close();
		logger.info("Wrote table of counts to file " + countFile);
		logger.info("Wrote table of percentages to file " + pctFile);*/
		
		String picardDir = this.configP.basicOptions.getPicardDirectory();
		// *** Novoalign steps ***
		// Only run if novoalign path was provided in config file
		if(this.configP.basicOptions.hasNovoalignPath()) {
		
			logger.info("");
			logger.info("Entering steps to align unmapped reads with Novoalign...");
			
			// Make sure genome novoindex was provided in config file
			if(!this.configP.basicOptions.getAllOptionMaps().containsKey(PipelineAutomatorConfigFileParser.OPTION_GENOME_NOVOINDEX)) {
				throw new IllegalArgumentException("Novoalign index for genome is required. Specify in config file with option genome_novoindex.");
			}
			
			String novoalign = this.configP.basicOptions.getNovoalignPath();
			String novoIndex = this.configP.basicOptions.getGenomeNovoindex();
			
			// Convert tophat unmapped.sam files to fastq files
			logger.info("");
			logger.info("Converting tophat unmapped bam files to fastq format...");
			ArrayList<String> convertJobIDs = new ArrayList<String>();
			// Store names of fastq files
			Map<String, String> unmappedFastq1 = new TreeMap<String,String>();
			for(String sample : this.sampleNames) {
				String dir = tophatDirsPerSample.get(sample);
				// Use Picard program SamToFastq
				String cmmd = "java -jar " + picardDir + "/SamToFastq.jar INPUT=" + dir + "/unmapped.bam VALIDATION_STRINGENCY=SILENT ";
				/* Tophat seems to discard pair information when writing unmapped reads to file unmapped.bam
				 * So proceed with unmapped reads as if they were single end
				 */
				String fastq1 = dir + "/unmapped.fq";
				File fastq1file = new File(fastq1);
				// Check if fastq file already exists
				if(fastq1file.exists()) {
					logger.info("WARNING: fastq file " + fastq1 + " already exists. Not rerunning format conversion.");
					unmappedFastq1.put(sample, fastq1);
					continue;
				}
				// Complete command
				cmmd += " FASTQ=" + fastq1;
				logger.info("Running Picard command: " + cmmd);
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				convertJobIDs.add(jobID);
				logger.info("LSF job ID is " + jobID + ".");
				// Submit job
				PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, dir + "/sam_to_fastq_" + jobID + ".bsub", "hour", 1);
			}
			logger.info("Waiting for SamToFastq jobs to finish...");
			PipelineUtils.waitForAllJobs(convertJobIDs, Runtime.getRuntime());
			logger.info("Done converting bam files to fastq format.");
		
			// Run novoalign with unmapped reads
			logger.info("");
			logger.info("Aligning unmapped reads to genome with Novoalign...");
			File novoDirFile = new File(NOVOALIGN_DIRECTORY);
			boolean madeNovDir = novoDirFile.mkdir();
			// Make maps for output locations
			Map<String, String> novoDirsPerSample = new TreeMap<String, String>();
			Map<String, String> novoSubdirsPerSample = new TreeMap<String, String>();
			Map<String, String> novoSamOutput = new TreeMap<String, String>();
			Map<String, String> novoBamOutput = new TreeMap<String, String>();
			Map<String, String> novoSortedBam = new TreeMap<String, String>();
			Map<String, String> novoBamFinalPath = new TreeMap<String, String>(); // the final location for bam file
			// Make directories for each sample
			for(String sample : this.sampleNames) {
				String dir = NOVOALIGN_DIRECTORY + "/" + NOVOALIGN_DIRECTORY + "_" + sample;
				File dirFile = new File(dir);
				boolean madeSampleDir = dirFile.mkdir();
				novoDirsPerSample.put(sample,dir);
				novoSubdirsPerSample.put(sample, NOVOALIGN_DIRECTORY + "_" + sample);
				novoSamOutput.put(sample, dir + "/" + novoalign + "_" + sample + ".sam");
				novoBamOutput.put(sample, dir + "/" + novoalign + "_" + sample + ".bam");
				novoSortedBam.put(sample, dir + "/" + novoalign + "_" + sample + ".sorted.bam");
				novoBamFinalPath.put(sample, NOVOALIGN_DIRECTORY + "/" + sample + ".bam");
			}
			// Get novoalign options
			Map<String, String> novoalignOptions = this.configP.basicOptions.getNovoalignOptions();
			// Override certain flags if provided in config file
			if(novoalignOptions.containsKey("-F")) {
				logger.info("WARNING: Overriding novoalign option -F provided in config file. Using STDFQ.");
				tophatOptions.remove("-F");
			}
			if(novoalignOptions.containsKey("-f")) {
				logger.info("WARNING: Overriding novoalign option -f provided in config file. Using unmapped reads from Tophat.");
				tophatOptions.remove("-f");
			}
			if(novoalignOptions.containsKey("-o")) {
				logger.info("WARNING: Overriding novoalign option -o provided in config file. Using SAM.");
				tophatOptions.remove("-o");
			}
			if(novoalignOptions.containsKey("-d")) {
				logger.info("WARNING: Overriding novoalign option -d provided in config file. Using " + novoIndex);
				tophatOptions.remove("-d");
			}
			// Make string of Novoalign options
			String novoOptionsString = " -d " + novoIndex + " -F STDFQ -o SAM ";
			for(String flag : novoalignOptions.keySet()) {
				novoOptionsString += flag + " " + novoalignOptions.get(flag) + " ";
			}
			
			// Run novoalign
			String cmmdBase = novoalign + novoOptionsString;
			ArrayList<String> novoJobIDs = new ArrayList<String>();
			Map<String, String> novoBsubFiles = new TreeMap<String, String>();
			for(String sample : this.sampleNames) {
				File outdir = new File(novoDirsPerSample.get(sample));
				outdir.mkdir();
				File sam = new File(novoSamOutput.get(sample));
				// Check if Novoalign has already been run
				if(sam.exists()) {
					logger.info("WARNING: alignment file " + sam + " already exists. Not rerunning alignment.");
					continue;
				}
				String cmmd = cmmdBase + " -f " + unmappedFastq1.get(sample);
				/* Tophat seems to discard pair information when writing unmapped reads to file unmapped.bam
				 * So proceed with unmapped reads as if they were single end
				 */
				logger.info("Writing novoalign output for sample " + sample + " to directory " + outdir + ".");
				logger.info("Running novoalign command: " + cmmd);
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				novoJobIDs.add(jobID);
				logger.info("LSF job ID is " + jobID + ".");
				// Submit job
				String bsubFile = outdir + "/novoalign_" + jobID + ".bsub";
				PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, bsubFile, "week", 8);
				novoBsubFiles.put(sample,bsubFile);
			}
			// Wait for novoalign jobs to finish
			logger.info("Waiting for novoalign jobs to finish...");
			PipelineUtils.waitForAllJobs(novoJobIDs, Runtime.getRuntime());
			logger.info("All samples done aligning unmapped reads with novoalign.");
			// Parse bsub output to sam files
			logger.info("");
			logger.info("Parsing LSF output to sam files...");
			StringParser stringparse = new StringParser();
			for(String sample : this.sampleNames) {
				File sam = new File(novoSamOutput.get(sample));
				// Check if Novoalign has already been run
				if(sam.exists()) {
					logger.info("WARNING: alignment file " + sam + " already exists. Not checking for or parsing bsub file.");
					continue;
				}
				String bsub = novoBsubFiles.get(sample);
				String bsub_only = bsub + ".bsub_output_only";
				logger.info("Parsing sam lines from bsub file " + bsub + " to sam file " + sam + "...");
				logger.info("Saving bsub output (excluding sam lines) to file " + bsub_only + "...");
				FileReader fr = new FileReader(bsub);
				BufferedReader br = new BufferedReader(fr);
				FileWriter fw = new FileWriter(sam);
				FileWriter fwb = new FileWriter(bsub_only);
				while(br.ready()) {
					String line = br.readLine();
					if(AlignmentUtils.isSamLine(stringparse, line)) {
						fw.write(line + "\n");
					} else {
						fwb.write(line + "\n");
					}
				}
				fw.close();
				fwb.close();
				fr.close();
				logger.info("Deleting bsub file " + bsub);
				File bsubFile = new File(bsub);
				bsubFile.delete();
			}
			logger.info("Done creating sam files.");
			
			
			
			// Reheader novoalign sam files to match tophat sam header
			// First get tophat header and write to novoalign directory
			Iterator<String> tophatBamIter = tophatBamFinalPath.keySet().iterator();
			String firstTophatBam = tophatBamFinalPath.get(tophatBamIter.next());
			String getHeaderCmmd = samtools + " view -H -o ";
			String tmpHeader = NOVOALIGN_DIRECTORY + "/tophat_sam_header_sorted.sam";
			String tophatHeader = NOVOALIGN_DIRECTORY + "/tophat_sam_header.sam";
			getHeaderCmmd += tmpHeader + " ";
			getHeaderCmmd += firstTophatBam;
			String getHeaderJobID = Long.valueOf(System.currentTimeMillis()).toString();
			logger.info("");
			logger.info("Getting sam header from tophat alignments to reheader novoalign sam files");
			logger.info("Running command: " + getHeaderCmmd);
			logger.info("LSF job ID is " + getHeaderJobID + ".");
			PipelineUtils.bsubProcess(Runtime.getRuntime(), getHeaderJobID, getHeaderCmmd, NOVOALIGN_DIRECTORY + "/get_sam_header_" + getHeaderJobID + ".bsub", "hour", 1);
			logger.info("Waiting for samtools view to finish...");
			PipelineUtils.waitForJobs(getHeaderJobID, Runtime.getRuntime());
			// Change sort order to unsorted
			FileReader r = new FileReader(tmpHeader);
			BufferedReader b = new BufferedReader(r);
			FileWriter wr = new FileWriter(tophatHeader);
			while(b.ready()) {
				String line = b.readLine();
				wr.write(line.replaceAll("coordinate", "unsorted") + "\n");
			}
			r.close();
			b.close();
			wr.close();
			File f = new File(tmpHeader);
			f.delete();
			
			
			logger.info("Done getting header.");
			// Now reheader each novoalign sam
			logger.info("");
			logger.info("Replacing headers in novoalign sam files with header from tophat alignments...");
			ArrayList<String> rhJobIDs = new ArrayList<String>();
			for(String sample : this.sampleNames) {
				String novo = novoSamOutput.get(sample);
				String reheadered = novoBamOutput.get(sample);
				String cmmd = "java -jar " + picardDir + "/ReplaceSamHeader.jar INPUT=" + novo + " HEADER=" + tophatHeader + " OUTPUT=" + reheadered;
				logger.info("Running Picard command: " + cmmd);
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				rhJobIDs.add(jobID);
				logger.info("LSF job ID is " + jobID + ".");
				PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, novoDirsPerSample.get(sample) + "/reheader_" + jobID + ".bsub", "hour", 1);
			}
			logger.info("Waiting for ReplaceSamHeader jobs to finish...");
			PipelineUtils.waitForAllJobs(rhJobIDs, Runtime.getRuntime());
			
			// Convert novoalign sam files to bam format
			/*logger.info("");
			logger.info("Converting novoalign sam files to bam format...");
			// Replace original sam file with reheadered file
			for(String sample : this.sampleNames) {
				File reheadered = new File(novoSamOutput.get(sample) + ".reheadered");
				boolean renamed = reheadered.renameTo(new File(novoBamOutput.get(sample)));
			}
			logger.info("All sam file headers replaced.");
			ArrayList<String> cbJobIDs = new ArrayList<String>();
			for(String sample : this.sampleNames) {
				File bam = new File(novoBamOutput.get(sample));
				// Check if bam files already exist
				if(bam.exists()) {
					logger.info("WARNING: bam file " + bam + " already exists. Not redoing format conversion from sam.");
					continue;
				}
				// Use Picard program SamFormatConverter
				String cmmd = "java -jar " + picardDir + "/SamFormatConverter.jar INPUT=" + novoSamOutput.get(sample) + " OUTPUT=" + novoBamOutput.get(sample);
				logger.info("Running Picard command: " + cmmd);
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				cbJobIDs.add(jobID);
				logger.info("LSF job ID is " + jobID + ".");
				PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, novoDirsPerSample.get(sample) + "/sam_to_bam_" + jobID + ".bsub", "hour", 1);
			}
			// Wait for jobs to finish
			logger.info("Waiting for SamFormatConverter jobs to finish...");
			PipelineUtils.waitForAllJobs(cbJobIDs, Runtime.getRuntime());
			logger.info("All samples done converting to bam format.");*/
			
			
			// Sort the bam files
			logger.info("");
			logger.info("Sorting novoalign bam files...");
			ArrayList<String> sbJobIDs = new ArrayList<String>();
			for(String sample : this.sampleNames) {
				File sorted = new File(novoSortedBam.get(sample));
				// Check if sorted files already exist
				if(sorted.exists()) {
					logger.info("WARNING: sorted bam file " + sorted + " already exists. Not resorting bam file.");
					continue;
				}
				// Use Picard program SortSam
				String cmmd = "java -jar " + picardDir + "/SortSam.jar INPUT=" + novoBamOutput.get(sample) + " OUTPUT=" + novoSortedBam.get(sample) + " SORT_ORDER=coordinate";
				logger.info("Running Picard command: " + cmmd);
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				sbJobIDs.add(jobID);
				logger.info("LSF job ID is " + jobID + ".");
				PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, novoDirsPerSample.get(sample) + "/sort_bam_" + jobID + ".bsub", "hour", 1);
			}
			// Wait for jobs to finish
			logger.info("Waiting for SortSam jobs to finish...");
			PipelineUtils.waitForAllJobs(sbJobIDs, Runtime.getRuntime());
			logger.info("All bam files sorted.");
			
			// Move all novoalign bam files to one directory
			logger.info("");
			logger.info("Moving all sorted novoalign bam files to directory " + NOVOALIGN_DIRECTORY + "...");
			for(String sample : this.sampleNames) {
				String cmmd = "mv " + novoSortedBam.get(sample) + " " + novoBamFinalPath.get(sample);
				Process p = Runtime.getRuntime().exec(cmmd, null);
				p.waitFor();
			}
	
			// Merge tophat and novoalign bam files
			logger.info("");
			logger.info("Merging tophat and novoalign bam files in directory " + MERGED_TOPHAT_NOVOALIGN_DIRECTORY	+ "...");
			File mergedDir = new File(MERGED_TOPHAT_NOVOALIGN_DIRECTORY);
			boolean madeMergedDir = mergedDir.mkdir();
			ArrayList<String> mergeJobIDs = new ArrayList<String>();
			for(String sample : this.sampleNames) {
				String tophatBam = tophatBamFinalPath.get(sample);
				String novoBam = novoBamFinalPath.get(sample); 
				String mergedBam = MERGED_TOPHAT_NOVOALIGN_DIRECTORY + "/" + sample + ".bam";
				File mergedFile = new File(mergedBam);
				// Check if merged files already exist
				if(mergedFile.exists()) {
					logger.info("WARNING: merged bam file " + mergedBam + " already exists. Not re-merging tophat and novoalign files.");
					continue;					
				}
				// Use Picard program MergeSamFiles
				String cmmd = "java -jar " + picardDir + "/MergeSamFiles.jar INPUT=" + tophatBam + " INPUT=" + novoBam + " OUTPUT=" + mergedBam;
				logger.info("Running Picard command: " + cmmd);
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				mergeJobIDs.add(jobID);
				logger.info("LSF job ID is " + jobID + ".");
				// Submit job
				PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, MERGED_TOPHAT_NOVOALIGN_DIRECTORY + "/merge_bams_" + jobID + ".bsub", "hour", 1);
			}
			// Wait for jobs to finish
			logger.info("Waiting for MergeSamFiles jobs to finish...");
			PipelineUtils.waitForAllJobs(mergeJobIDs, Runtime.getRuntime());
			logger.info("All bam files merged.");
			
			// Update current bam files to merged bams
			logger.info("");
			logger.info("Updating current bam files...");
			for(String sample : this.sampleNames) {
				this.currentBamFiles.put(sample, MERGED_TOPHAT_NOVOALIGN_DIRECTORY + "/" + sample + ".bam");
				logger.info("Current bam file for sample " + sample + " is " + this.currentBamFiles.get(sample));
			}
			// Update current bam directory
			this.currentBamDir = MERGED_TOPHAT_NOVOALIGN_DIRECTORY;
	
			// Count mapped and unmapped reads
/*			logger.info("");
			logger.info("Counting mapped and unmapped reads...");
			String mergedCountFile = MERGED_TOPHAT_NOVOALIGN_DIRECTORY + "/mapped_unmapped_count.out";
			String mergedPctFile = MERGED_TOPHAT_NOVOALIGN_DIRECTORY + "/mapped_unmapped_percentage.out";
			FileWriter mw = new FileWriter(mergedCountFile);
			FileWriter mwp = new FileWriter(mergedPctFile);
			String mergedHeader = "Sample\tMapped\tUnmapped\n";
			w.write(mergedHeader);
			wp.write(mergedHeader);
			for(String sample : this.sampleNames) {
				int mapped = AlignmentUtils.countAlignments(samtools, MERGED_TOPHAT_NOVOALIGN_DIRECTORY + "/" + sample + ".bam", MERGED_TOPHAT_NOVOALIGN_DIRECTORY, false, false);
				int unmapped = AlignmentUtils.countAlignments(samtools, MERGED_TOPHAT_NOVOALIGN_DIRECTORY + "/" + sample + ".bam", MERGED_TOPHAT_NOVOALIGN_DIRECTORY, true, false);
				int total = mapped + unmapped;
				mw.write(sample + "\t" + mapped + "\t" + unmapped + "\n");
				mwp.write(sample + "\t" + (double)mapped/(double)total + "\t" + (double)unmapped/(double)total + "\n");
			}
			mw.close();
			mwp.close();
			logger.info("Wrote table of counts to file " + mergedCountFile);
			logger.info("Wrote table of percentages to file " + mergedPctFile);

	*/
	
		} // *** Done with novoalign steps ***
		
		// Index current bam files
		logger.info("");
		logger.info("Indexing bam files...");
		ArrayList<String> indexJobIDs = new ArrayList<String>();
		for(String sample : this.sampleNames) {
			String bam = this.currentBamFiles.get(sample);
			String index = this.currentBamDir + "/" + sample + ".bai";
			File indexfile = new File(index);
			// Check if bai files exist
			if(indexfile.exists()) {
				logger.info("WARNING: index " + index + " already exists. Not re-indexing bam file.");
				continue;
			}
			String cmmd = samtools + " index " + bam + " " + index;
			logger.info("Running samtools command: " + cmmd);
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			indexJobIDs.add(jobID);
			logger.info("LSF job ID is " + jobID + ".");
			// Submit job
			PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, "index_bam_" + jobID + ".bsub", "hour", 1);
		}
		logger.info("Waiting for samtools jobs to finish...");
		PipelineUtils.waitForAllJobs(indexJobIDs, Runtime.getRuntime());
		logger.info("All bam files indexed.");
		
		// Make tdf files from current bam files
		logger.info("");
		logger.info("Making tdf files...");
		String assembly = this.configP.basicOptions.getGenomeAssembly();
		String igvtools = this.configP.basicOptions.getIgvtoolsPath();
		ArrayList<String> tdfJobIDs = new ArrayList<String>();
		for(String sample : this.sampleNames) {
			String bam = this.currentBamFiles.get(sample);
			String tdf = this.currentBamDir + "/" + sample + ".tdf";
			File tdffile = new File(tdf);
			// Check if tdf files exist
			if(tdffile.exists()) {
				logger.info("WARNING: tdf file " + tdf + " already exists. Not remaking tdf file.");
				continue;
			}
			// Use igvtools count
			String cmmd = igvtools + " count -w 3 " + bam + " " + tdf + " " + assembly;
			logger.info("Running igvtools command: " + cmmd);
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			tdfJobIDs.add(jobID);
			logger.info("LSF job ID is " + jobID + ".");
			// Submit job
			PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, "make_tdf_" + jobID + ".bsub", "hour", 1);
			
		}
		logger.info("Waiting for igvtools jobs to finish...");
		// Igvtools count always ends by crashing even though it worked, so catch the exception and check if files were really created
		try {
			PipelineUtils.waitForAllJobs(tdfJobIDs, Runtime.getRuntime());
		} catch(IllegalArgumentException e) {
			boolean ok = true;
			for(String sample : this.sampleNames) {
				String bam = this.currentBamFiles.get(sample);
				String tdf = this.currentBamDir + "/" + sample + ".tdf";
				File tdffile = new File(tdf);
				if(!tdffile.exists()) {
					ok = false;
					break;
				}
			}
			if(!ok) {
				throw new IllegalArgumentException(e.getMessage());
			}
		}
		logger.info("All tdf files created.");
		
		
		// Collect Picard metrics
		logger.info("");
		logger.info("Collecting Picard metrics for bam files in directory " + this.currentBamDir + "...");
		String picardMetricsDir = this.currentBamDir + "/picard_metrics";
		logger.info("Writing all metrics in directory " + picardMetricsDir + "....");
		File picardMetricsDirFile = new File(picardMetricsDir);
		boolean madePicardMetricsDir = picardMetricsDirFile.mkdir();
		ArrayList<String> pmJobIDs = new ArrayList<String>();
		// Inputs to picard metrics
		String refFlat = this.configP.basicOptions.getPicardRefFlat();
		String ribIntervals = this.configP.basicOptions.getPicardRibosomalIntervals();
		String genomeFasta = this.configP.basicOptions.getGenomeFasta();
		String strandSpecificity = this.configP.basicOptions.getPicardStrandSpecificity();
		
		for(String sample : this.sampleNames) {
			
			// Establish output files
			String asMetrics = picardMetricsDir + "/" + sample + "_alignmentSummaryMetrics.out";
			String isMetrics = picardMetricsDir + "/" + sample + "_insertSizeMetrics.out";
			String isHistogram = picardMetricsDir + "/" + sample + "_insertSizeMetrics.histogram";
			String rsMetrics = picardMetricsDir + "/" + sample + "_rnaSeqMetrics.out";
			String rsChart = picardMetricsDir + "/" + sample + "_rnaSeqMetrics_NormalizedPositionCoverage.pdf";
			File asFile = new File(asMetrics);
			File isFile = new File(isMetrics);
			File ishFile = new File(isHistogram);
			File rsFile = new File(rsMetrics);
			File rscFile = new File(rsChart);
			
			// Run alignment summary metrics
			if(asFile.exists()) {
				logger.info("WARNING: alignment summary metrics file " + asMetrics + " already exists. Not rerunning alignment summary metrics.");
			} else {
				// Use Picard program CollectAlignmentSummaryMetrics
				String cmmd = "java -jar " + picardDir + "/CollectAlignmentSummaryMetrics.jar";
				cmmd += " INPUT=" + this.currentBamFiles.get(sample);
				cmmd += " OUTPUT=" + asMetrics;
				cmmd += " REFERENCE_SEQUENCE=" + genomeFasta;
				logger.info("Running Picard command: " + cmmd);
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				pmJobIDs.add(jobID);
				logger.info("LSF job ID is " + jobID + ".");
				// Submit job
				PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, picardMetricsDir + "/picard_alignment_summary_metrics_" + jobID + ".bsub", "hour", 1);
			}
			
			// Run insert size metrics
			if(isFile.exists() && ishFile.exists()) {
				logger.info("WARNING: insert size metrics files " + isMetrics + " and " + isHistogram + " already exist. Not rerunning insert size metrics.");
			} else {
				// Use Picard program CollectInsertSizeMetrics
				String cmmd = "java -jar " + picardDir + "/CollectInsertSizeMetrics.jar";
				cmmd += " INPUT=" + this.currentBamFiles.get(sample);
				cmmd += " OUTPUT=" + isMetrics;
				cmmd += " REFERENCE_SEQUENCE=" + genomeFasta;
				cmmd += " HISTOGRAM_FILE=" + isHistogram;
				logger.info("Running Picard command: " + cmmd);
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				pmJobIDs.add(jobID);
				logger.info("LSF job ID is " + jobID + ".");
				// Submit job
				PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, picardMetricsDir + "/picard_insert_size_metrics_" + jobID + ".bsub", "hour", 1);
			}
			
			// Run RNA-seq metrics
			if(rsFile.exists() && rscFile.exists()) {
				logger.info("WARNING: RNA-seq metrics files " + rsMetrics + " and " + rsChart + " already exist. Not rerunning RNA-seq metrics.");
			} else {
				// Use Picard program CollectRnaSeqMetrics
				String cmmd = "java -jar ";
				// Use user-provided CollectRnaSeqMetrics executable if provided
				if(this.configP.basicOptions.getPicardCollectRnaSeqMetrics() != null) cmmd += this.configP.basicOptions.getPicardCollectRnaSeqMetrics();
				else cmmd += picardDir + "/CollectRnaSeqMetrics.jar";
				cmmd += " INPUT=" + this.currentBamFiles.get(sample);
				cmmd += " OUTPUT=" + rsMetrics;
				cmmd += " REFERENCE_SEQUENCE=" + genomeFasta;
				cmmd += " CHART_OUTPUT=" + rsChart;
				cmmd += " REF_FLAT=" + refFlat;
				if(ribIntervals != null) cmmd += " RIBOSOMAL_INTERVALS=" + ribIntervals;
				cmmd += " STRAND_SPECIFICITY=" + strandSpecificity;
				logger.info("Running Picard command: " + cmmd);
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				pmJobIDs.add(jobID);
				logger.info("LSF job ID is " + jobID + ".");
				// Submit job
				PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, picardMetricsDir + "/picard_rnaseq_metrics_" + jobID + ".bsub", "hour", 1);
			}
			
		}
		
		// Wait for jobs to finish
		logger.info("Waiting for Picard metrics jobs to finish...");
		PipelineUtils.waitForAllJobs(pmJobIDs, Runtime.getRuntime());
		logger.info("All Picard metrics done.");

		logger.info("");
		logger.info("Genome alignments all done.\n");
		
	}
	
	private String createDGEInputFile(String inputListFile) throws IOException{
		
		String newInputFile = inputListFile+".dgeInputFile";
		FileWriter w = new FileWriter(newInputFile);
		for(String name: this.currentBamFiles.keySet()){
			w.write(name+"\t"+this.currentBamFiles.get(name)+"\t"+this.nameToCondition.get(name)+"\n");
		}
		w.close();
		
		return newInputFile;
	}

	/**
	 * TASK 6: RUN DGE
	 * @author skadri
	 * @throws MathException 
	 * @throws ParseException 
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	private void runDGE(boolean runBasic,String DGEInputFile) throws IOException, ParseException, MathException, InterruptedException{
		
		//ALIGNMENTS
		if("multiple".equalsIgnoreCase(this.configP.dgeOptions.getrunType())){
			/*
			 * Construct the DGE argument list
			 */
			ArrayList<String> arguments =  new ArrayList<String>();
			//3'DGE or 5'DGE
			//single/multiple
			//THis determines task
			arguments.add("-task");
			if("3DGE".equalsIgnoreCase(this.configP.dgeOptions.getdgeType())){
				arguments.add("score3PMultiple");
			}
			else if("5DGE".equalsIgnoreCase(this.configP.dgeOptions.getdgeType())){
				arguments.add("score5PMultiple");
			}
			else
				throw new IllegalArgumentException("Illegal DGE type (3DGE/5DGE): "+this.configP.dgeOptions.getdgeType());
			
			//Changing the output to a directory
			this.configP.dgeOptions.setOutputInDirectory();
			
			arguments.add("-alignments");
			arguments.add(DGEInputFile);
			for(String flag:this.configP.getRunSpecificDGEOptions()){
				arguments.add("-"+flag);
				arguments.add(this.configP.dgeOptions.getOption(flag));
			}
			logger.info("Parameters to DGE:" + arguments);
			DGE dge = new DGE(l2a(arguments));
		}
		else if("single".equalsIgnoreCase(this.configP.dgeOptions.getrunType())){
			DGEInput inp = readDGEInputFile(DGEInputFile);
			this.queueName = configP.basicOptions.getQueueName();
			if(queueName==null)
				queueName = "week";
			//DGE arguments
			String optionsString = "java -jar "+this.configP.dgeOptions.getdgePath()+" -task ";
			if("3DGE".equalsIgnoreCase(this.configP.dgeOptions.getdgeType())){
				optionsString += "score3P ";
			}
			else if("5DGE".equalsIgnoreCase(this.configP.dgeOptions.getdgeType())){
				optionsString += "score5P ";
			}
			else{
				throw new IllegalArgumentException("Illegal DGE type (3DGE/5DGE): "+this.configP.dgeOptions.getdgeType());
			}
			for(String flag:this.configP.getRunSpecificDGEOptions()){
				if(flag.equalsIgnoreCase("out"))
					;//dont add
				else{
					optionsString += "-"+flag+" ";
					optionsString += this.configP.dgeOptions.getOption(flag)+" ";
				}
			}
			
			//OUTPUT FILE NAMES ARE DIFF
			for(String sample:inp.getSamples()){
				//MULTIPLE BSUB ROUTINES FOR EACH DGE
				String thisCommand = optionsString;
				//ALIGNMENT
				thisCommand +="-alignment ";
				thisCommand += inp.getBamFileFor(sample)+" ";
				//OTHER FLAGS
				
				//Changing the output to a directory
				String outputName = inp.getBamFileFor(sample)+".dge/"+inp.getBamFileFor(sample)+".dge.exp";
				thisCommand += "-out ";
				thisCommand += outputName;
				
				
				//Submit the job
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				logger.info("LSF job ID is " + jobID + ".");
				// Submit job
				PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, thisCommand, outputName+"_"+ jobID + ".bsub", queueName, 5);
				//PipelineUtils.
			}
		}
		else{
			throw new IllegalArgumentException("Illegal Run type (single/multiple): "+this.configP.dgeOptions.getrunType());
		}
	}

	/**
	 * Reads the "Name		Bam_file	condition" list file
	 * Helper function to runDGE()
	 * @param inputListFile
	 * @throws FileNotFoundException 
	 */
	private static DGEInput readDGEInputFile(String inputListFile) throws FileNotFoundException {
		
		DGEInput dge = new DGEInput();
		Scanner reader = new Scanner(new File(inputListFile));
		while (reader.hasNextLine()) {
			String[] str = StringParser.getTokens(reader.nextLine());
			dge.addInput(str);
		}
		return dge;
	}


	/**
	 * Helper function: Converts list to array and returns the array
	 * @param list
	 * @return
	 */
	private static String[] l2a(List<String> list){
		String[] rtrn=new String[list.size()];
		int i=0;
		for(String val: list){rtrn[i++]=val;}
		return rtrn;
	}

	public static void main (String [] args) throws IOException, ParseException, MathException, InterruptedException{
		
		CommandLineParser p = new CommandLineParser();
		p.setProgramDescription("\n*** Configurable pipeline for RNA-seq read processing and analysis ***");
		p.addStringArg("-r", "File containing list of fastq files. \n\t\tLine format: \n\t\t<sample_name> <left.fq> <right.fq> <condition> \n\t\tOR \n\t\t<sample_name> <unpaired.fq> <condition>", true);
		p.addStringArg("-c", "Config file", true);
		p.parse(args);
		String fastqList = p.getStringArg("-r");
		String configFile = p.getStringArg("-c");
		
		PipelineAutomator PA = new PipelineAutomator(fastqList,configFile);
	}

	public static class DGEInput {
		Map<String,String> nameToBam;
		Map<String,String> nameToCondition;
		
		public DGEInput(){
			this.nameToBam = new HashMap<String,String>();
			this.nameToCondition = new HashMap<String,String>();
		}
		
		/**
		 * 
		 * @param str
		 */
		public void addInput(String[] str){
			//Mandatory 3 columns are required.
			if(str.length<3){
				throw new IllegalArgumentException("Illegal number of columns in input DGE");
			}
			this.nameToBam.put(str[0], str[1]);
			this.nameToCondition.put(str[0], str[2]);
		}
		
		public Set<String> getSamples(){
			return this.nameToBam.keySet();
		}
		
		public String getBamFileFor(String name){
			return this.nameToBam.get(name);
		}
	}
	
	/**
	 * Bowtie command lines needed for this pipeline
	 * @author prussell
	 *
	 */
	public static class AlignmentUtils {
		
		static int MAX_INSERT_SIZE = 1000;
		
		/**
		 * Check whether the string represents a valid sam file line
		 * @param line The string
		 * @return True if the line is a sam header line or a valid alignment line, false otherwise
		 */
		public static boolean isSamLine(String line) {
			StringParser p = new StringParser();
			return isSamLine(p, line);
		}

		/**
		 * Check whether the string represents a valid sam file line using an exisiting StringParser object
		 * @param p The StringParser object
		 * @param line The string
		 * @return True if the line is a sam header line or a valid alignment line, false otherwise
		 */		
		public static boolean isSamLine(StringParser p, String line) {
			p.parse(line);
			if(p.getFieldCount() == 0) return false;
			if(p.getFieldCount() >= 11) {
				// Line might be an alignment line
				// Check if correct fields are ints
				try {
					int i1 = p.asInt(1);
					int i3 = p.asInt(3);
					int i4 = p.asInt(4);
					int i7 = p.asInt(7);
					int i8 = p.asInt(8);
					return true;
				} catch(NumberFormatException e) {
					// The value in a sam int position is not an int
					return false;
				}
			}
			// Check if this is a header line
			String first = p.asString(0);
			if(first.equals("@HD") || first.equals("@SQ") || first.equals("@RG") || first.equals("@PG") || first.equals("@CO")) return true;
			// Line is neither a header line nor an alignment line
			return false;
		}
		

		
		/**
		 * Count total number of mapped reads in a sam or bam file
		 * @param samtoolsExecutable Samtools executable file
		 * @param samFile The sam file
		 * @param logDir Directory to write log information to
		 * @param samFormat True if sam format, false if bam format
		 * @return The alignment count obtained by running the command "samtools view -c -F 4"
		 * @throws IOException
		 * @throws InterruptedException
		 */
		public static int countAlignments(String samtoolsExecutable, String samFile, String logDir, boolean samFormat) throws IOException, InterruptedException {
			return countAlignments(samtoolsExecutable, samFile, logDir, false, samFormat);
		}
		
		/**
		 * Count total number of mapped or unmapped reads in a sam file
		 * @param samtoolsExecutable Samtools executable file
		 * @param samFile The sam file
		 * @param logDir Directory to write log information to
		 * @param countUnmapped True if counting unmapped reads, false if counting mapped reads
		 * @param samFormat True if sam format, false if bam format
		 * @return The alignment count obtained by running the command "samtools view -c -F 4"
		 * @throws IOException
		 * @throws InterruptedException
		 */
		public static int countAlignments(String samtoolsExecutable, String samFile, String logDir, boolean countUnmapped, boolean samFormat) throws IOException, InterruptedException {

			File dir = new File(logDir);
			boolean madeDir = dir.mkdir();
			
			// Use samtools to count alignments
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			String output = logDir + "/count_alignments_" + jobID + ".bsub";
			String cmmd = samtoolsExecutable + " view -c ";
			if(samFormat) cmmd += " -S ";
			if(!countUnmapped) cmmd += " -F 4 ";
			else cmmd += " -f 4 ";
			cmmd += samFile;
			
			// Capture output of samtools process and read into int
			Process p = Runtime.getRuntime().exec(cmmd);
			InputStream is = p.getInputStream();
			InputStreamReader isr = new InputStreamReader(is);
			BufferedReader br = new BufferedReader(isr);
			String line = br.readLine();
			int count = Integer.parseInt(line);
			is.close();
			isr.close();
			br.close();
			
			String type = countUnmapped ? "unmapped" : "mapped";
			logger.info("There are " + count + " " + type + " reads in sam file " + samFile + ".");
			
			return count;
			
		}

		
		
		/**
		 * Make bowtie2 index of fasta file
		 * @param fastaFile The fasta file
		 * @param outBtIndexBase Output index file without .bt2 extension
		 * @param bowtie2BuildExecutable The bowtie2-build executable
		 * @param bsubOutputDir Output directory for bsub file
		 * @throws IOException
		 * @throws InterruptedException
		 */
		public static void makeBowtie2Index(String fastaFile, String outBtIndexBase, String bowtie2BuildExecutable, String bsubOutputDir) throws IOException, InterruptedException {
			logger.info("");
			logger.info("Writing bowtie2 index for file " + fastaFile + " to files " + outBtIndexBase);
			File bsubDir = new File(bsubOutputDir);
			boolean madeDir = bsubDir.mkdir();
			
			// Check if index already exists
			String f1 = outBtIndexBase + ".1.bt2";
			String f2 = outBtIndexBase + ".2.bt2";
			String f3 = outBtIndexBase + ".3.bt2";
			String f4 = outBtIndexBase + ".4.bt2";
			String r1 = outBtIndexBase + ".rev.1.bt2";
			String r2 = outBtIndexBase + ".rev.2.bt2";
			File f1file = new File(f1);
			File f2file = new File(f2);
			File f3file = new File(f3);
			File f4file = new File(f4);
			File r1file = new File(r1);
			File r2file = new File(r2);
			if(f1file.exists() && f2file.exists() && f3file.exists() && f4file.exists() && r1file.exists() && r2file.exists()) {
				logger.info("WARNING: bowtie2 index files already exist. Not writing new files.");
				return;
			}
			
			String cmmd = bowtie2BuildExecutable + " " + fastaFile + " " + outBtIndexBase;
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			String output = bsubOutputDir + "/make_bowtie_index_" + jobID + ".bsub";
			int prc = PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, output, "week", 4);
			PipelineUtils.waitForJobs(jobID, Runtime.getRuntime());
			logger.info("Done creating bowtie2 index for file " + fastaFile);
		}
		
		/**
		 * Run bowtie2 with unpaired reads and default max insert size
		 * @param bowtie2IndexBase Index filename prefix (minus trailing .X.bt2)
		 * @param reads Fastq file containing unpaired reads
		 * @param outSamFile Output sam file
		 * @param outUnalignedFastq Output fastq file for unaligned reads
		 * @param bowtie2Executable Bowtie2 executable
		 * @param bsubOutDir Output directory for bsub file
		 * @throws IOException
		 * @throws InterruptedException
		 * @return The bsub job ID
		 */
		public static String runBowtie2(String bowtie2IndexBase, String reads, String outSamFile, String outUnalignedFastq, String bowtie2Executable, String bsubOutDir) throws IOException, InterruptedException {
			return runBowtie2(bowtie2IndexBase, reads, outSamFile, outUnalignedFastq, bowtie2Executable, MAX_INSERT_SIZE, bsubOutDir);
		}

		
		/**
		 * Run bowtie2 with unpaired reads and specified max insert size
		 * @param bowtie2IndexBase Index filename prefix (minus trailing .X.bt2)
		 * @param reads Fastq file containing unpaired reads
		 * @param outSamFile Output sam file
		 * @param outUnalignedFastq Output fastq file for unaligned reads
		 * @param bowtie2Executable Bowtie2 executable
		 * @param maxInsertSize Max insert size
		 * @param bsubOutDir Output directory for bsub file
		 * @throws IOException
		 * @throws InterruptedException
		 * @return The bsub job ID
		 */
		public static String runBowtie2(String bowtie2IndexBase, String reads, String outSamFile, String outUnalignedFastq, String bowtie2Executable, int maxInsertSize, String bsubOutDir) throws IOException, InterruptedException {
			return runBowtie2(bowtie2IndexBase, reads, null, outSamFile, outUnalignedFastq, bowtie2Executable, maxInsertSize, bsubOutDir, false);
		}
		
		/**
		 * Run bowtie2 with paired reads and default max insert size
		 * @param bowtie2IndexBase Index filename prefix (minus trailing .X.bt2)
		 * @param read1Fastq Fastq file containing read1 if paired or single end reads if unpaired
		 * @param read2Fastq Fastq file containing read2 if paired or null if unpaired
		 * @param outSamFile Output sam file
		 * @param outUnalignedFastq Output fastq file for unaligned reads
		 * @param bowtie2Executable Bowtie2 executable
		 * @param bsubOutDir Output directory for bsub file
		 * @throws IOException
		 * @throws InterruptedException
		 * @return The bsub job ID
		 */
		public static String runBowtie2(String bowtie2IndexBase, String read1Fastq, String read2Fastq, String outSamFile, String outUnalignedFastq, String bowtie2Executable, String bsubOutDir) throws IOException, InterruptedException {
			return runBowtie2(bowtie2IndexBase, read1Fastq, read2Fastq, outSamFile, outUnalignedFastq, bowtie2Executable, MAX_INSERT_SIZE, bsubOutDir, true);
		}
		
		/**
		 * Run bowtie2 with paired reads and specified max insert size
		 * @param bowtie2IndexBase Index filename prefix (minus trailing .X.bt2)
		 * @param read1Fastq Fastq file containing read1 if paired or single end reads if unpaired
		 * @param read2Fastq Fastq file containing read2 if paired or null if unpaired
		 * @param outSamFile Output sam file
		 * @param outUnalignedFastq Output fastq file for unaligned reads
		 * @param bowtie2Executable Bowtie2 executable
		 * @param maxInsertSize Max insert size
		 * @param bsubOutDir Output directory for bsub files
		 * @param readsPaired Whether the reads are paired
		 * @throws IOException
		 * @throws InterruptedException
		 * @return The bsub job ID
		 */
		public static String runBowtie2(String bowtie2IndexBase, String read1Fastq, String read2Fastq, 
				String outSamFile, String outUnalignedFastq, String bowtie2Executable, int maxInsertSize, 
				String bsubOutDir, boolean readsPaired) throws IOException, InterruptedException {

			File bsubOutDirectory = new File(bsubOutDir);
			boolean madeDir = bsubOutDirectory.mkdir();
			
			String executable = bowtie2Executable + " ";
			String index = "-x " + bowtie2IndexBase + " ";
			
			String reads;
			if(read2Fastq != null) reads = "-1 " + read1Fastq + " -2 " + read2Fastq + " ";
			else reads = "-U " + read1Fastq + " ";

			String options = "-q --local -M 0 --met-stderr ";
			options += "--maxins " + maxInsertSize + " ";
			if(!readsPaired) options += "--un " + outUnalignedFastq + " ";
			else options += "--un-conc " + outUnalignedFastq + " ";
			
			String sam = "-S " + outSamFile;
			
			String cmmd = executable + options + index + reads + sam;
			
			logger.info("");
			logger.info("Running bowtie2 command:");
			logger.info(cmmd);
			
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			String output = bsubOutDir + "/run_bowtie_" + jobID + ".bsub";
			logger.info("Writing bsub output to file " + output);
			int prc = PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, output, "week", 4);
			logger.info("Job ID is " + jobID);
			return jobID;
			
		}

		/**
		 * Run bowtie2 with paired reads and specified max insert size
		 * @param bowtie2IndexBase Index filename prefix (minus trailing .X.bt2)
		 * @param read1Fastq Fastq file containing read1 if paired or single end reads if unpaired
		 * @param read2Fastq Fastq file containing read2 if paired or null if unpaired
		 * @param outSamFile Output sam file
		 * @param outUnalignedFastq Output fastq file for unaligned reads
		 * @param bowtie2Executable Bowtie2 executable
		 * @param maxInsertSize Max insert size
		 * @param bsubOutDir Output directory for bsub files
		 * @param readsPaired Whether the reads are paired
		 * @throws IOException
		 * @throws InterruptedException
		 * @return The bsub job ID
		 */
		public static String runBowtie(String bowtieIndexBase, String read1Fastq, String read2Fastq, 
				String outSamFile, String outUnalignedFastq, String bowtie2Executable, int maxInsertSize, 
				String bsubOutDir, boolean readsPaired) throws IOException, InterruptedException {

			File bsubOutDirectory = new File(bsubOutDir);
			boolean madeDir = bsubOutDirectory.mkdir();
			
			String executable = bowtie2Executable + " ";
			String index = bowtieIndexBase + " ";
			
			String reads;
			if(read2Fastq != null) reads = "-1 " + read1Fastq + " -2 " + read2Fastq + " ";
			else reads = "-12 " + read1Fastq + " ";

			String options = "-q -k 1 --best ";
			options += "--maxins " + maxInsertSize + " ";
			options += "--un " + outUnalignedFastq + " ";
			
			String sam = "--sam " + outSamFile;
			
			String cmmd = executable + options + index + reads + sam;
			
			logger.info("");
			logger.info("Running bowtie command:");
			logger.info(cmmd);
			
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			String output = bsubOutDir + "/run_bowtie_" + jobID + ".bsub";
			logger.info("Writing bsub output to file " + output);
			int prc = PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, output, "week", 4);
			logger.info("Job ID is " + jobID);
			return jobID;
			
		}
	}
	
}

