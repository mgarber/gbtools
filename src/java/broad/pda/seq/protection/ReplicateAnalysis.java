package broad.pda.seq.protection;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.parser.CommandLineParser;
import broad.pda.gene.RefSeqGene;

/**
 * @author prussell
 *
 */
public class ReplicateAnalysis {

	/**
	 * Map associating sample name with paired sample object for the dataset
	 */
	private Map<String,PairedSampleCoverageAnalysis> pairedSamples;
	
	/**
	 * The transcript annotation
	 */
	private Map<String, Collection<RefSeqGene>> transcripts;
	
	/**
	 * The window size
	 */
	private int windowSize;
	
	/**
	 * The step size
	 */
	private int stepSize;
	
	/**
	 * Constructor with two samples being provided as bam files, write individual sample significant windows to file
	 * @param sample1name Sample 1 name
	 * @param backgroundAlignmentFile1 Sample 1 bam alignment file for background data
	 * @param signalAlignmentFile1 Sample 1 bam alignment file for signal data
	 * @param sample2name Sample 2 name
	 * @param backgroundAlignmentFile2 Sample 2 bam alignment file for background data
	 * @param signalAlignmentFile2 Sample 2 bam alignment file for signal data
	 * @param bedFile Bed file of gene annotations
	 * @param chrSizeFile Chromosome size file
	 * @param window Window size
	 * @param step Step size
	 * @param alphaSkellam Max P value to call window significant for skellam statistic
	 * @param alphaScan Max P value for scan statistic
	 * @throws IllegalArgumentException
	 * @throws IOException
	 */
	public ReplicateAnalysis(String sample1name, String backgroundAlignmentFile1, String signalAlignmentFile1, String sample2name, String backgroundAlignmentFile2, String signalAlignmentFile2, String bedFile, String chrSizeFile, int window, int step, double alphaSkellam, double alphaScan, double trimQuantile) throws IllegalArgumentException, IOException {
		this(sample1name, new PairedSampleCoverageAnalysis(backgroundAlignmentFile1, signalAlignmentFile1, bedFile, chrSizeFile, window, step, alphaSkellam, alphaScan, trimQuantile), sample2name, new PairedSampleCoverageAnalysis(backgroundAlignmentFile2, signalAlignmentFile2, bedFile, chrSizeFile, window, step, alphaSkellam, alphaScan, trimQuantile));
	}
	
	/**
	 * Constructor with three samples being provided as bam files, write individual sample significant windows to file
	 * @param sample1name Sample 1 name
	 * @param backgroundAlignmentFile1 Sample 1 bam alignment file for background data
	 * @param signalAlignmentFile1 Sample 1 bam alignment file for signal data
	 * @param sample2name Sample 2 name
	 * @param backgroundAlignmentFile2 Sample 2 bam alignment file for background data
	 * @param signalAlignmentFile2 Sample 2 bam alignment file for signal data
	 * @param sample3name Sample 3 name
	 * @param backgroundAlignmentFile3 Sample 3 bam alignment file for background data
	 * @param signalAlignmentFile3 Sample 3 bam alignment file for signal data 
	 * @param bedFile Bed file of gene annotations
	 * @param chrSizeFile Chromosome size file
	 * @param window Window size
	 * @param step Step size
	 * @param alphaSkellam Max P value for skellam statistic
	 * @param alphaScan Max P value for scan statistic
	 * @throws IllegalArgumentException
	 * @throws IOException
	 */
	public ReplicateAnalysis(String sample1name, String backgroundAlignmentFile1, String signalAlignmentFile1, String sample2name, String backgroundAlignmentFile2, String signalAlignmentFile2, String sample3name, String backgroundAlignmentFile3, String signalAlignmentFile3, String bedFile, String chrSizeFile, int window, int step, double alphaSkellam, double alphaScan, double trimQuantile) throws IllegalArgumentException, IOException {
		this(sample1name, new PairedSampleCoverageAnalysis(backgroundAlignmentFile1, signalAlignmentFile1, bedFile, chrSizeFile, window, step, alphaSkellam, alphaScan, trimQuantile), sample2name, new PairedSampleCoverageAnalysis(backgroundAlignmentFile2, signalAlignmentFile2, bedFile, chrSizeFile, window, step, alphaSkellam, alphaScan, trimQuantile), sample3name, new PairedSampleCoverageAnalysis(backgroundAlignmentFile3, signalAlignmentFile3, bedFile, chrSizeFile, window, step, alphaSkellam, alphaScan, trimQuantile));		
	}

	/**
	 * Constructor with two samples being provided as paired data objects
	 * @param sample1name Sample 1 name
	 * @param sample1pairedCoverageAnalysis Sample 1 paired data object
	 * @param sample2name Sample 2 name
	 * @param sample2pairedCoverageAnalysis Sample 2 paired data object
	 */
	private ReplicateAnalysis(String sample1name, PairedSampleCoverageAnalysis sample1pairedCoverageAnalysis, String sample2name, PairedSampleCoverageAnalysis sample2pairedCoverageAnalysis) {
		Map<String,PairedSampleCoverageAnalysis> replicateData = new TreeMap<String,PairedSampleCoverageAnalysis>();
		replicateData.put(sample1name, sample1pairedCoverageAnalysis);
		replicateData.put(sample2name, sample2pairedCoverageAnalysis);
		this.pairedSamples = replicateData;
		this.stepSize = sample1pairedCoverageAnalysis.getStepSize();
		this.windowSize = sample1pairedCoverageAnalysis.getWindowSize();
		this.transcripts = sample1pairedCoverageAnalysis.getGenes();
	}
	
	/**
	 * Constructor with three samples being provided as paired data objects
	 * @param sample1name Sample 1 name
	 * @param sample1pairedCoverageAnalysis Sample 1 paired data object
	 * @param sample2name Sample 2 name
	 * @param sample2pairedCoverageAnalysis Sample 2 paired data object
	 * @param sample3name Sample 3 name
	 * @param sample3pairedCoverageAnalysis Sample 3 paired data object
	 */
	private ReplicateAnalysis(String sample1name, PairedSampleCoverageAnalysis sample1pairedCoverageAnalysis, String sample2name, PairedSampleCoverageAnalysis sample2pairedCoverageAnalysis, String sample3name, PairedSampleCoverageAnalysis sample3pairedCoverageAnalysis) {
		Map<String,PairedSampleCoverageAnalysis> replicateData = new TreeMap<String,PairedSampleCoverageAnalysis>();
		replicateData.put(sample1name, sample1pairedCoverageAnalysis);
		replicateData.put(sample2name, sample2pairedCoverageAnalysis);
		replicateData.put(sample3name, sample3pairedCoverageAnalysis);
		this.pairedSamples = replicateData;		
		this.stepSize = sample1pairedCoverageAnalysis.getStepSize();
		this.windowSize = sample1pairedCoverageAnalysis.getWindowSize();
		this.transcripts = sample1pairedCoverageAnalysis.getGenes();
	}	
	
	/**
	 * Constructor with arbitrary set of paired data objects
	 * @param replicateData Map associating sample name with paired data object
	 */
	public ReplicateAnalysis(Map<String,PairedSampleCoverageAnalysis> replicateData) {
		this.pairedSamples = replicateData;
		Iterator<String> sampleIter = this.pairedSamples.keySet().iterator();
		PairedSampleCoverageAnalysis sample1pairedCoverageAnalysis = this.pairedSamples.get(sampleIter.next());
		this.stepSize = sample1pairedCoverageAnalysis.getStepSize();
		this.windowSize = sample1pairedCoverageAnalysis.getWindowSize();
		this.transcripts = sample1pairedCoverageAnalysis.getGenes();
	}
	
	/**
	 * Get windows that are significant in all replicates and write individual sample significant sets to bed files
	 * @param alphaSkellam P value cutoff for skellam statistic
	 * @param alphaScan Max P value for scan statistic
	 * @param bedFile Prefix of bed file to write significant windows for each sample
	 * @return The set of windows that are significant in all replicates by this cutoff
	 * @throws IOException 
	 */
	private Collection<RefSeqGene> intersectSignificantWindows(double alphaSkellam, double alphaScan) throws IOException {
		
		System.err.println("\nGetting significant windows for each sample...");
		
		ArrayList<Collection<RefSeqGene>> eachSampleSigWindows = new ArrayList<Collection<RefSeqGene>>();
		Collection<RefSeqGene> rtrn = new TreeSet<RefSeqGene>();
		
		for(String sample : this.pairedSamples.keySet()) {
			Map<RefSeqGene, Collection<RefSeqGene>> sampleSigWindows = this.pairedSamples.get(sample).getSignificantWindows(alphaSkellam, alphaScan);
			for(RefSeqGene gene : sampleSigWindows.keySet()) {
				for(RefSeqGene window : sampleSigWindows.get(gene)) window.setBedScore(0); // Reset bed scores so windows will be recognized as the same
				eachSampleSigWindows.add(sampleSigWindows.get(gene));
			}
		}
		
		System.err.println("\nIntersecting significant window sets.");
		
		Iterator<String> iter = this.pairedSamples.keySet().iterator();
		PairedSampleCoverageAnalysis firstSample = this.pairedSamples.get(iter.next());
		
		for(String chr : this.transcripts.keySet()) {
			for(RefSeqGene gene : this.transcripts.get(chr)) {
				for(RefSeqGene window : firstSample.getWindows(gene)) {
					boolean inAllSamples = true;
					for(Collection<RefSeqGene> sigWindows : eachSampleSigWindows) {
						if(!sigWindows.contains(window)) {
							inAllSamples = false;
							break;
						}
					}
					if(inAllSamples) rtrn.add(window);
				}
			}
		}
		
		System.err.println("Found " + rtrn.size() + " windows that are significant in all " + this.pairedSamples.keySet().size() + " replicates.");
		
		return rtrn;
		
	}
	
	/**
	 * Write windows that are significant in all replicates with FDR correction to a bed file and write individual sample significant window sets to bed files
	 * @param alphaSkellam P value cutoff for skellam statistic
	 * @param alphaScan Max P value for scan statistic
	 * @param outFile Output bed file for intersection
	 * @throws IOException
	 */
	private void writeSignificantWindowsIntersection(double alphaSkellam, double alphaScan, String outFile) throws IOException {
		Collection<RefSeqGene> sigWindows = intersectSignificantWindows(alphaSkellam, alphaScan);
		System.err.println("Writing intersection of significant windows to file " + outFile);
		RefSeqGene.writeBedFile(sigWindows, outFile);
		System.err.println("Wrote file.\n");
	}
	
	/**
	 * Combine position level counts of all background samples
	 * @return SampleCounts object whose underlying data is the combined counts
	 */
	/*public SampleCounts mergeBackgroundSamples() {
		// TODO
	}*/
	
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws IllegalArgumentException 
	 */
	public static void main(String[] args) throws IllegalArgumentException, IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-n1", "Sample 1 name", true);
		p.addStringArg("-b1", "Background 1 bam file", true);
		p.addStringArg("-s1", "Signal 1 bam file", true);
		p.addStringArg("-n2", "Sample 2 name", false);
		p.addStringArg("-b2", "Background 2 bam file", false);
		p.addStringArg("-s2", "Signal 2 bam file", false);
		p.addStringArg("-n3", "Sample 3 name", false);
		p.addStringArg("-b3", "Background 3 bam file", false);
		p.addStringArg("-s3", "Signal 3 bam file", false);
		p.addStringArg("-g", "Gene bed file", true);
		p.addStringArg("-s", "Chromosome size file", true);
		p.addIntegerArg("-w", "Window size", true);
		p.addIntegerArg("-t", "Step size", true);
		p.addStringArg("-o", "Output bed file for intersection of significant windows", true);
		p.addStringArg("-ss", "Output bed file prefix for individual samples", false);
		p.addFloatArg("-psk", "P value cutoff for skellam statistic", false, Float.valueOf((float)0.05));
		p.addFloatArg("-psc", "P value cutoff for scan statistic", false, Float.valueOf((float)0.05));
		p.addFloatArg("-q", "Coverage quantile within untrimmed peak for trimming by max contiguous subsequence", false, Float.valueOf((float)0.5));
		p.parse(args);
		String sample1name = p.getStringArg("-n1");
		String backgroundAlignmentFile1 = p.getStringArg("-b1");
		String signalAlignmentFile1 = p.getStringArg("-s1");
		String sample2name = p.getStringArg("-n2");
		String backgroundAlignmentFile2 = p.getStringArg("-b2");
		String signalAlignmentFile2 = p.getStringArg("-s2");
		String sample3name = p.getStringArg("-n3");
		String backgroundAlignmentFile3 = p.getStringArg("-b3");
		String signalAlignmentFile3 = p.getStringArg("-s3");
		String bedFile = p.getStringArg("-g");
		String chrSizeFile = p.getStringArg("-s");
		double trimQuantile = p.getFloatArg("-q").floatValue();
		int window = p.getIntegerArg("-w").intValue();
		int step = p.getIntegerArg("-t").intValue();
		String outFile = p.getStringArg("-o");
		String indSampleBedPrefix = p.getStringArg("-ss");
		float alphaSkellam = p.getFloatArg("-psk").floatValue();
		float alphaScan = p.getFloatArg("-psc").floatValue();

		Map<String,PairedSampleCoverageAnalysis> replicateData = new TreeMap<String,PairedSampleCoverageAnalysis>();
		
		System.err.println("\nScoring data for sample " + sample1name + ".");
		replicateData.put(sample1name, new PairedSampleCoverageAnalysis(backgroundAlignmentFile1, signalAlignmentFile1, bedFile, chrSizeFile, window, step, alphaSkellam, alphaScan, trimQuantile));
		if(indSampleBedPrefix != null) {
			PairedSampleCoverageAnalysis psca = replicateData.get(sample1name);
			String bed = indSampleBedPrefix + sample1name + ".bed";
			psca.writePeaksAsBed(bed, psca.getPvalueCutoffSkellam(), psca.getPvalueCutoffScan());
		}
		System.err.println("Done scoring sample " + sample1name + ".");
		
		if(sample2name != null && backgroundAlignmentFile2 != null && signalAlignmentFile2 != null) {
			System.err.println("Scoring data for sample " + sample2name + ".");
			replicateData.put(sample2name, new PairedSampleCoverageAnalysis(backgroundAlignmentFile2, signalAlignmentFile2, bedFile, chrSizeFile, window, step, alphaSkellam, alphaScan, trimQuantile));
			if(indSampleBedPrefix != null) {
				PairedSampleCoverageAnalysis psca = replicateData.get(sample2name);
				String bed = indSampleBedPrefix + sample2name + ".bed";
				psca.writePeaksAsBed(bed, psca.getPvalueCutoffSkellam(), psca.getPvalueCutoffScan());
			}
			System.err.println("Done scoring sample " + sample2name + ".");
		}
		
		if(sample3name != null && backgroundAlignmentFile3 != null && signalAlignmentFile3 != null) {
			System.err.println("Reading data for sample " + sample3name + ".");
			replicateData.put(sample3name, new PairedSampleCoverageAnalysis(backgroundAlignmentFile3, signalAlignmentFile3, bedFile, chrSizeFile, window, step, alphaSkellam, alphaScan, trimQuantile));
			if(indSampleBedPrefix != null) {
				PairedSampleCoverageAnalysis psca = replicateData.get(sample3name);
				String bed = indSampleBedPrefix + sample3name + ".bed";
				psca.writePeaksAsBed(bed, psca.getPvalueCutoffSkellam(), psca.getPvalueCutoffScan());
			}
			System.err.println("Done scoring sample " + sample3name + ".");
		}
		
		System.err.println("\nDone scoring individual samples.");
		ReplicateAnalysis par = new ReplicateAnalysis(replicateData);
		par.writeSignificantWindowsIntersection(alphaSkellam, alphaScan, outFile);
		System.err.println("\nAll done.");
		
	}

}
