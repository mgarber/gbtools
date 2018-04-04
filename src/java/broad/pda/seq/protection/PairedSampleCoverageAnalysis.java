package broad.pda.seq.protection;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.IntervalTree;
import broad.core.math.Distribution;
import broad.core.parser.CommandLineParser;
import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;

/**
 * @author prussell
 *
 */
public class PairedSampleCoverageAnalysis {

	private SampleCounts backgroundCounts;
	private Map<String, Map<RefSeqGene,double[]>> backgroundExpressionScores;
	private SampleCounts signalCounts;
	private Map<String, Collection<RefSeqGene>> genes;
	private static double EXPRESSION_PVALUE_CUTOFF = 1;
	private static int MIN_WINDOW_COUNT = 10;
	private int windowSize;
	private int stepSize;
	private double pValueCutoffSkellam;
	private double pValueCutoffScan;
	private double trimPeakByQuantile;
		
	/**
	 * Construct object from bam files
	 * @param backgroundAlignmentFile Bam file of background sample alignments
	 * @param signalAlignmentFile Bam file of signal sample alignments
	 * @param bedFile Bed file of annotations
	 * @param chrSizeFile Table of chromosome sizes
	 * @param window Window size
	 * @param step Step size
	 * @param alphaSkellam Max P value for skellam statistic
	 * @param alphaScan Max P value for scan statistic
	 * @param trimQuantile Coverage quantile within untrimmed peak for trimming by max contiguous subsequence
	 * @throws IllegalArgumentException
	 * @throws IOException
	 */
	public PairedSampleCoverageAnalysis(String backgroundAlignmentFile, String signalAlignmentFile, String bedFile, String chrSizeFile, int window, int step, double alphaSkellam, double alphaScan, double trimQuantile) throws IllegalArgumentException, IOException {
		this(new SampleCounts(backgroundAlignmentFile, chrSizeFile), new SampleCounts(signalAlignmentFile, chrSizeFile), bedFile, window, step, alphaSkellam, alphaScan, trimQuantile);
	}
	
	/**
	 * Construct object from ContinuousDataAlignmentModels
	 * @param backgroundData Background sample alignment data
	 * @param signalData Signal sample alignment data
	 * @param bedFile Bed file of annotations
	 * @param window Window size
	 * @param step Step size
	 * @param alphaSkellam Max P value for skellam statistic
	 * @param alphaScan Max P value for scan statistic
	 * @param trimQuantile Coverage quantile within untrimmed peak for trimming by max contiguous subsequence
	 * @throws IOException
	 */
	public PairedSampleCoverageAnalysis(ContinuousDataAlignmentModel backgroundData, ContinuousDataAlignmentModel signalData, String bedFile, int window, int step, double alphaSkellam, double alphaScan, double trimQuantile) throws IOException {
		this(new SampleCounts(backgroundData), new SampleCounts(signalData), bedFile, window, step, alphaSkellam, alphaScan, trimQuantile);
	}
	
	/**
	 * Construct object from SampleCounts objects
	 * @param backgroundSampleCounts Background SampleCounts data
	 * @param signalSampleCounts Signal SampleCounts data
	 * @param bedFile Bed file of annotations
	 * @param window Window size
	 * @param step Step size
	 * @param alphaSkellam Max P value for skellam statistic
	 * @param alphaScan Max P value for scan statistic
	 * @param outFilePrefix Output bed file for significant peaks
	 * @throws IOException
	 */
	private PairedSampleCoverageAnalysis(SampleCounts backgroundSampleCounts, SampleCounts signalSampleCounts, String bedFile, int window, int step, double alphaSkellam, double alphaScan, double trimQuantile) throws IOException {
		this.backgroundCounts = backgroundSampleCounts;
		this.signalCounts = signalSampleCounts;
		this.genes = new TreeMap<String, Collection<RefSeqGene>>();
		this.genes.putAll(BEDFileParser.loadDataByChr(new File(bedFile)));		
		this.backgroundExpressionScores = new TreeMap<String, Map<RefSeqGene,double[]>>();
		this.windowSize = window;
		this.stepSize = step;
		this.pValueCutoffSkellam = alphaSkellam;
		this.pValueCutoffScan = alphaScan;
		this.trimPeakByQuantile = trimQuantile;
		for(String chr : this.genes.keySet()) {
			this.backgroundExpressionScores.put(chr, this.backgroundCounts.getData().scoreGenes(this.genes.get(chr)));
		}
		System.err.println("Instantiated PairedSampleCoverageAnalysis object.");
	}
	
	/**
	 * Compute Skellam P-value of read counts in a region given the two parameters
	 * @param gene The region
	 * @param backgroundLambda Poisson lambda for background sample
	 * @param signalLambda Poisson lambda for signal sample
	 * @return The probability under the null hypothesis of observing a greater difference
	 * @throws IOException
	 */
	private static double getSkellamPvalue(double backgroundLambda, double signalLambda, int backgroundCount, int signalCount) throws IOException {
		return Distribution.skellamRightTail(signalCount - backgroundCount, signalLambda, backgroundLambda);
	}
	
	
	/**
	 * Get the gene set
	 * @return The gene set
	 */
	public Map<String, Collection<RefSeqGene>> getGenes() {
		return this.genes;
	}
	
	/**
	 * Get the window size
	 * @return The window size
	 */
	public int getWindowSize() {return this.windowSize;}
	
	/**
	 * Get the step size
	 * @return The step size
	 */
	public int getStepSize() {return this.stepSize;}
	
	/**
	 * Get the p value cutoff for skellam statistic
	 * @return The p value cutoff
	 */
	public double getPvalueCutoffSkellam() {return this.pValueCutoffSkellam;}
	
	/**
	 * Get the p value cutoff for scan statistic
	 * @return The p value cutoff
	 */
	public double getPvalueCutoffScan() {return this.pValueCutoffScan;}

	
	/**
	 * Get the set of sliding windows in a gene as RefSeqGene objects
	 * @param gene The gene
	 * @return The set of all windows as RefSeqGene objects
	 */
	public Collection<RefSeqGene> getWindows(RefSeqGene gene) {
		return gene.getWindows(this.windowSize, this.stepSize, 0);
	}
	
	
	/**
	 * Get all windows passing skellam score cutoff
	 * @param gene The gene
	 * @return All significant windows, or null if can't make Poisson model due to lack of reads or no significant windows
	 * @throws IOException
	 */
	private Collection<RefSeqGene> getSkellamSignificantWindows(RefSeqGene gene) throws IOException {
			
		try {
			if(this.backgroundExpressionScores.get(gene.getChr()).get(gene)[0] > EXPRESSION_PVALUE_CUTOFF) return null;
		} catch(NullPointerException e) {
			throw new NullPointerException("Can't get expression score for gene " + gene.getName() + " (" + gene.getChr() + ")");
		}
		
		int treeStartPos = Math.min(0,gene.getStart() - 100);
		int treeEndPos = gene.getEnd() + 100;
		IntervalTree<org.broad.igv.sam.Alignment> backgroundIntervalTree = this.backgroundCounts.getData().getData().getIntervalTree(gene.getChr(), treeStartPos, treeEndPos);
		IntervalTree<org.broad.igv.sam.Alignment> signalIntervalTree = this.signalCounts.getData().getData().getIntervalTree(gene.getChr(), treeStartPos, treeEndPos);
		
		Collection<RefSeqGene> windows = this.getWindows(gene);
		Collection<RefSeqGene> sigWindows = new TreeSet<RefSeqGene>();
		
		try {
			double backgroundLambda = this.backgroundCounts.getPoissonLambda(gene, this.windowSize);
			double signalLambda = this.signalCounts.getPoissonLambda(gene, this.windowSize);
			System.err.println("Scoring gene " + gene.getName() + ". Got " + windows.size() + " " + this.windowSize + "bp windows. Background lambda = " + backgroundLambda + ". Signal lambda = " + signalLambda + ".");
			if(signalLambda == 0) return null;

			for(RefSeqGene window : windows) {
				int backgroundCount = (int) Math.round(this.backgroundCounts.getData().getData().getCountsPerAlignment(window, backgroundIntervalTree, 0));
				if(backgroundCount <= 0) continue;
				int signalCount = (int) Math.round(this.signalCounts.getData().getData().getCountsPerAlignment(window, signalIntervalTree, 0));
				if(signalCount < Math.max(MIN_WINDOW_COUNT, signalLambda)) continue; // Skip if signal read count is too low to be interesting
				double pval = PairedSampleCoverageAnalysis.getSkellamPvalue(backgroundLambda, signalLambda, backgroundCount, signalCount);
				if(pval <= this.pValueCutoffSkellam) sigWindows.add(window);
			}
		} catch(IllegalArgumentException e) {
			System.err.println("Skipping gene " + gene.getName() + " because can't make poisson model." );
			return null;
		}

		if(sigWindows.isEmpty()) return null;
		return sigWindows;
		
	}
	
	/**
	 * Get windows passing skellam score cutoff from all genes
	 * @return Map associating genes containing significant windows with their set of significant windows
	 * @throws IOException
	 */
	private Map<RefSeqGene, Collection<RefSeqGene>> getSkellamSigWindowsAllGenes() throws IOException {
		System.err.println("Getting windows that pass skellam threshold...");
		Map<RefSeqGene, Collection<RefSeqGene>> rtrn = new TreeMap<RefSeqGene, Collection<RefSeqGene>>();
		for(String chr : this.genes.keySet()) {
			for(RefSeqGene gene : this.genes.get(chr)) {
				Collection<RefSeqGene> geneSigWindows = this.getSkellamSignificantWindows(gene);
				if(geneSigWindows == null) continue;
				rtrn.put(gene, geneSigWindows);
			}
		}
		return rtrn;
	}
	


	/**
	 * Get all fixed size windows below a skellam P value cutoff and scan P value cutoff
	 * @param alphaSkellam The P value cutoff for skellam statistic
	 * @param alphaScan The P value cutoff for scan statistic
	 * @return The set of significant windows
	 * @throws IOException 
	 */
	public Map<RefSeqGene, Collection<RefSeqGene>> getSignificantWindows(double alphaSkellam, double alphaScan) throws IOException {
		
		System.err.println("Getting significant fixed-width windows...");
		Map<RefSeqGene, Collection<RefSeqGene>> skellamSigWindows = this.getSkellamSigWindowsAllGenes();
		return this.filterScanPvalue(skellamSigWindows, alphaScan);
	}
	
	
	/**
	 * Get merged and trimmed significant peaks
	 * @param alphaSkellam P value cutoff for skellam statistic
	 * @param alphaScan P value cutoff for scan statistic
	 * @return Map associating each gene with its significant windows or empty map if there are none
	 * @throws IOException
	 */
	public Map<RefSeqGene, Collection<RefSeqGene>> getSignificantPeaks(double alphaSkellam, double alphaScan) throws IOException {
		// Get fixed size windows that pass cutoff
		Map<RefSeqGene, Collection<RefSeqGene>> sigWindowsFixedSize = this.getSignificantWindows(alphaSkellam, alphaScan);
		Map<RefSeqGene, Collection<RefSeqGene>> rtrn = new TreeMap<RefSeqGene, Collection<RefSeqGene>>();
		
		System.err.println("Merging significant windows and trimming peaks...");
		for(RefSeqGene gene : sigWindowsFixedSize.keySet()) {
			// Sort the significant windows from the gene
			TreeSet<RefSeqGene> sortedWindows = new TreeSet<RefSeqGene>();
			sortedWindows.addAll(sigWindowsFixedSize.get(gene));
			
			// Merge overlapping windows
			Collection<RefSeqGene> mergedWindows = RefSeqGene.mergeAllExonsAnyOrientation(sortedWindows);
			
			// Trim merged regions
			Collection<RefSeqGene> trimmedWindows = new TreeSet<RefSeqGene>();
			for(RefSeqGene window : mergedWindows) {
				try {
					RefSeqGene trimmed = this.signalCounts.getData().trimMaxContiguous(window, this.trimPeakByQuantile);
					if(trimmed != null) trimmedWindows.add(trimmed);
				} catch(NullPointerException e) {
					continue;
				}
			}
			
			if(!trimmedWindows.isEmpty()) rtrn.put(gene, trimmedWindows);
		}
		return rtrn;
		
	}
	
	/**
	 * Filter a collection of regions for significant scan statistic
	 * @param windows The regions
	 * @param alpha P value cutoff
	 * @return The significant windows in the set
	 * @throws IOException
	 */
	private Map<RefSeqGene, Collection<RefSeqGene>> filterScanPvalue(Map<RefSeqGene, Collection<RefSeqGene>> windows, double alpha) throws IOException {
		Map<RefSeqGene, Collection<RefSeqGene>> rtrn = new TreeMap<RefSeqGene, Collection<RefSeqGene>>();
		for(RefSeqGene gene : windows.keySet()) {
			Collection<RefSeqGene> scanSigWindows = new TreeSet<RefSeqGene>();
			for(RefSeqGene window : windows.get(gene)) {
				double scanPval = this.signalCounts.getData().scanPRate(window)[0];
				if(scanPval < alpha) scanSigWindows.add(window);
			}
			if(!scanSigWindows.isEmpty()) rtrn.put(gene, scanSigWindows);
		}
		return rtrn;
	}
	
	
	/**
	 * Write significant windows to bed file
	 * @param outFile The output file
	 * @param alphaSkellam The P value cutoff for skellam statistic
	 * @param alphaScan The P value cutoff for scan statistic
	 * @throws IOException
	 */
	public void writePeaksAsBed(String outFile, double alphaSkellam, double alphaScan) throws IOException {
		FileWriter w = new FileWriter(outFile);
		Map<RefSeqGene, Collection<RefSeqGene>> significantPeaks = getSignificantPeaks(alphaSkellam, alphaScan);
		System.err.println("\nWriting peaks to bed file " + outFile);
		for(RefSeqGene gene : significantPeaks.keySet()) {
			for(RefSeqGene window : significantPeaks.get(gene))	w.write(window.toBED(255,0,0) + "\n");
		}
		w.close();
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws IllegalArgumentException 
	 */
	public static void main(String[] args) throws IllegalArgumentException, IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Background bam file", true);
		p.addStringArg("-s", "Signal bam file", true);
		p.addStringArg("-g", "Gene bed file", true);
		p.addStringArg("-c", "Chromosome size table", true);
		p.addIntegerArg("-w", "Window size", true);
		p.addIntegerArg("-t", "Step size", true);
		p.addFloatArg("-psk", "P value cutoff for skellam statistic", false, Float.valueOf((float)0.05));
		p.addFloatArg("-psc", "P value cutoff for scan statistic", false, Float.valueOf((float)0.05));
		p.addFloatArg("-q", "Coverage quantile within untrimmed peak for trimming by max contiguous subsequence", false, Float.valueOf((float)0.5));
		p.addStringArg("-o", "Output file prefix", true);
		p.parse(args);
		String backgroundFile = p.getStringArg("-b");
		String signalFile = p.getStringArg("-s");
		String bedFile = p.getStringArg("-g");
		String chrSizeFile = p.getStringArg("-c");
		int windowSize = p.getIntegerArg("-w").intValue();
		int stepSize = p.getIntegerArg("-t").intValue();
		double alphaSkellam = p.getFloatArg("-psk").floatValue();
		double alphaScan = p.getFloatArg("-psc").floatValue();
		double trimQuantile = p.getFloatArg("-q").floatValue();
		String outFilePrefix = p.getStringArg("-o");
		
		
		PairedSampleCoverageAnalysis svbs = new PairedSampleCoverageAnalysis(backgroundFile, signalFile, bedFile, chrSizeFile, windowSize, stepSize, alphaSkellam, alphaScan, trimQuantile);
		
		if(outFilePrefix != null) {
			svbs.writePeaksAsBed(outFilePrefix + ".bed", svbs.pValueCutoffSkellam, svbs.pValueCutoffScan);
		}

		
	}

}
