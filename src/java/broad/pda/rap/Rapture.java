/**
 * 	Jesse Engreitz
 *  August 20, 2012
 */
package broad.pda.rap;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.apache.commons.collections15.CollectionUtils;
import org.apache.commons.collections15.Predicate;
import org.broad.igv.Globals;

import broad.core.annotation.BED;
import broad.core.annotation.BEDReader;
import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.annotation.ShortBEDReader;
import broad.core.error.ParseException;
import broad.core.math.EmpiricalDistribution;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.datastructures.AlignmentsReader;
import broad.pda.seq.alignment.AlignmentUtils;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.slide.Slide;
import broad.pda.seq.slide.SlideAndCountConsumer;
import umontreal.iro.lecuyer.probdist.GammaDist;
import umontreal.iro.lecuyer.probdist.NegativeBinomialDist;

/**
 * @author engreitz
 * Performs RAP-DNA analysis
 */
public class Rapture {

	private ContinuousDataAlignmentModel targetModel, controlModel;
	private Map<String, Integer> chromosomeSizes;
	private File[] maskFiles;
	private static double OVERLAP_FACTOR = 0.8;
	private static double DATA_WEIGHT = 1;
	
	// index of alignments scores
	public static int PVALUE_INDEX = 0;
	public static int FOLD_INDEX = 7;
	
	
	
	public Rapture(String targetFile, String controlFile, int extensionFactor, String sizeFile, String maskFileDir, boolean loadPairsAsFragments) throws IOException, ParseException {
		this.maskFiles  = maskFileDir != null ?  new File(maskFileDir).listFiles(): null ;
		Map<String, Integer> maskFileData = ContinuousDataAlignmentModel.parseMaskFiles(maskFiles);
		
		// loadChromosomeStats = true, minMappingQuality=0, removeDuplicateFlags=true, weighReadCounts=false, strand=null
		this.targetModel = AlignmentUtils.loadAlignmentData(targetFile, true, 0, true, false, null, loadPairsAsFragments);
		this.targetModel.setMaskFileData(maskFileData);
		this.targetModel.setExtensionFactor(extensionFactor);
		
		this.controlModel = AlignmentUtils.loadAlignmentData(controlFile, true, 0, true, false, null, loadPairsAsFragments);
		this.controlModel.setMaskFileData(maskFileData);
		this.controlModel.setExtensionFactor(extensionFactor);
		
		this.chromosomeSizes = BEDFileParser.loadChrSizes(sizeFile);
	}
	
	
	public Rapture(Map<String, Integer> sizes, ContinuousDataAlignmentModel target, ContinuousDataAlignmentModel control, File[] maskFiles) {
		this.chromosomeSizes = sizes;
		this.targetModel = target;
		this.controlModel = control;
		this.maskFiles = maskFiles;
	}
	
	
	/**
	 * 
	 * @param region			either "all" or "chr1".  TODO: make this an array of arbitrary regions
	 * @param windows
	 * @return
	 * @throws IOException
	 * @throws ParseException
	 */
	public Map<Integer, List<Alignments>> scan(String region, int[] windows) throws Exception, IOException, ParseException {
		return scan(region, windows, null);
	}
	
	public Map<Integer, List<Alignments>> scan(String region, int[] windows, Predicate<Alignments> filter) throws Exception, IOException, ParseException {
		List<String> chromosomes = parseRegions(region);
		
		Map<Integer, List<Alignments>> result = new HashMap<Integer, List<Alignments>>();
		for (int window : windows) {
			for (int i = 0; i < chromosomes.size(); i++) {
				String chr = chromosomes.get(i);
				List<Alignments> chrResults = scanChromosome(chr, window);
				
				// Use filter if provided
				if (filter != null) {
					CollectionUtils.filter(chrResults, filter);
				}
				
				Integer wind = new Integer(window);
				if (i == 0)
					result.put(wind, chrResults);
				else
					result.get(wind).addAll(chrResults);
			}
		}
		return result;
	}

	
	/**
	 * Filters Alignments based on their pvalue and fold-difference (target vs control)
	 * Usage:  Predicate<Alignments> filter = new PValueFilter(cutoff, fold);  CollectionUtils.filter(alignments, filter);
	 * @author engreitz
	 *
	 */
	public static class PValueFilter implements Predicate<Alignments> {
		private double cutoff, fold;
		PValueFilter(double cutoff, double fold) { this.cutoff = cutoff; this.fold = fold; }
		public boolean evaluate(Alignments obj) {
			return (obj.getScores().get(PVALUE_INDEX) >= cutoff) && (obj.getScores().get(FOLD_INDEX) >= fold);
		}
	}
	
	
	public List<String> parseRegions(String region) {
		List<String> chromosomes = new ArrayList<String>();
		if (region == null || region.equals("all")) {
			// Skip chrM and chrY
			for (String chr : chromosomeSizes.keySet()) {
				if (!chr.equals("chrM") && !chr.equals("chrY"))
					chromosomes.add(chr);
			}
		} else {
			chromosomes.add(region);
		}
		return chromosomes;
	}
	
	
	/**
	 * Scan a chromosome and return all of the windows and scores
	 * @param chr
	 * @param window
	 * @return
	 * @throws IOException
	 * @throws ParseException
	 */
	private List<Alignments> scanChromosome(String chr, int window) throws Exception, IOException, ParseException {
		System.out.println("Processing " + chr);
		ShortBEDReader maskedRegions = Slide.loadMaskedRegions(maskFiles, chr);
		
		// First scan the control sample
		BasicLightweightAnnotation region = new BasicLightweightAnnotation(chr, 1, chromosomeSizes.get(chr));
		int overlap = (int)Math.round((float)window * OVERLAP_FACTOR);

		// Calculate read counts in each window
		Slide controlSlider = new Slide(window, overlap, maskedRegions);
		SlideAndCountConsumer consumer = new SlideAndCountConsumer(controlModel);
		controlSlider.slide(region, consumer);
		List<Double> controlCounts = consumer.getScores();

		// Estimate prior distribution for Poisson parameter from the control model
		GammaDist gamma = estimateGammaFromControl(controlCounts);
		double controlTotalReads = controlModel.getData().getTotalReads();        // controlModel.getNumberOfReads(chr);
		double targetTotalReads = targetModel.getData().getTotalReads();          // targetModel.getNumberOfReads(chr);
		
		
		// Adjust prior distribution according to number of reads in the target model
		double shape = gamma.getAlpha();
		// change to units of fragments per window per million reads
		double rate = gamma.getLambda() / controlTotalReads * 1000000.0;
		System.out.println("Adjusted rate: " + rate);

		// Now scan the target sample
		Slide targetSlider = new Slide(window, overlap, maskedRegions);
		NBSliderConsumer nbconsumer = new NBSliderConsumer(shape, rate, controlCounts, controlTotalReads, targetTotalReads, targetModel);
		targetSlider.slide(region, nbconsumer);
		List<Alignments> scores = nbconsumer.getWindows();
		
		return scores;
	}
	
	
	/**************************************************************************************************
	 * SlideAndCount consumer that calculates and saves a negative binomial p-value for each window
	 * @author engreitz
	 */
	private class NBSliderConsumer extends SlideAndCountConsumer {
		private List<Alignments> windows;
		private GammaDist prior;
		private List<Double> controlCounts;
		private double controlTotalReads;
		private double targetTotalReads;
		private int counter;
		
		public NBSliderConsumer(double shape, double rate, List<Double> controlCounts, double controlReads, double targetReads, ContinuousDataAlignmentModel d) throws IOException {
			super(targetModel);
			windows = new ArrayList<Alignments>();
			prior = new GammaDist(shape, rate);
			this.controlCounts = controlCounts;
			counter = 0;
			this.controlTotalReads = controlReads;
			this.targetTotalReads = targetReads;
			//System.out.println("ControlTotal\tTargetTotal\tC_Count\tT_Count\tPrior_Shape\tPrior_Rate\tPrior_Mean\tPost_Shape\tPost_Rate\tPost_Mean\tPValue");
		}
		
		@Override
		public int getWindowBatch() {return 100000;}
		
		@Override
		public void consume(List<Alignments> windows) throws Exception, IOException {
			super.consume(windows);
			for(Alignments w : windows) {
				addWindowScores(w);
				processWindow(w);
				counter++;  // tracks the control reads with the incoming windows
			}
		}

		private void addWindowScores(Alignments w) {
			// Get the observed count from the target distribution
			int targetCount = (int)w.getCountScore();
			int controlCount = controlCounts.get(counter).intValue();
			
			// Update the gamma distribution parameters with the control counts
			double posteriorShape = prior.getAlpha() + (controlCount / controlTotalReads * targetTotalReads);
			double posteriorRate = (prior.getLambda() * targetTotalReads / 1000000.0) + DATA_WEIGHT;
			double nbProb = posteriorRate / (posteriorRate + 1);
			double nominalP = -1*Math.log10(1 - NegativeBinomialDist.cdf(posteriorShape, nbProb, targetCount));
			
			double foldDiff = (targetCount/targetTotalReads*1000000 + 1) / (controlCount/controlTotalReads*1000000 + 1);  // Add phantom counts to avoid dividing by zero
			List<Double> scores = new ArrayList<Double>();
			scores.add(nominalP);
			scores.add((double)targetCount);
			scores.add((double)controlCount);
			scores.add(targetTotalReads);
			scores.add(controlTotalReads);
			scores.add((double) targetCount/targetTotalReads*1000000);
			scores.add((double) controlCount/controlTotalReads*1000000);
			scores.add(foldDiff);
			w.setScores(scores); 
			//System.out.println(controlTotalReads + "\t" + targetTotalReads + "\t" + controlCounts.get(counter) + "\t" + observedCount + "\t" + prior.getAlpha() + "\t" + prior.getLambda() * targetTotalReads / 1000000 + "\t" + prior.getAlpha() * prior.getLambda() * targetTotalReads / 1000000 + "\t" + posteriorShape + "\t" + posteriorRate + "\t" + posteriorRate * posteriorShape + "\t" + nominalP);
		}
		
		/**
		 * Override this function to change the saving / outputting behavior of the consumer
		 * @param w
		 */
		private void processWindow(Alignments w) {
			this.windows.add(w);
		}
		
		// Model is:  X = counts observed in a given window
		//			  X ~ Poisson(lambda), lambda ~ Gamma(shape, rate)
		//      Then  f(x) ~ NegativeBinomial(shape, rate/(rate+1))
			
		public List<Alignments> getWindows() { return windows; }
	}
	
	/*
	private class CallPeaksSliderConsumer extends NBSliderConsumer {
		private Predicate<Alignments> filter;
		public CallPeaksSliderConsumer(double shape, double rate, List<Integer> controlCounts, double controlReads, double targetReads, Predicate<Alignments> filter) {
			super(shape, rate, controlCounts, controlReads, targetReads);
		}
	}*/
	
	public void scoreSegments(Collection<Alignments> segments) throws IOException {
		double controlTotalReads = controlModel.getData().getTotalReads();  
		double targetTotalReads = targetModel.getData().getTotalReads();  
		for (Alignments w : segments) {
			scoreSegment(w, controlTotalReads, targetTotalReads);
		}
	}
	
	public void scoreSegments(AlignmentsReader segments) throws IOException {
		double controlTotalReads = controlModel.getData().getTotalReads();  
		double targetTotalReads = targetModel.getData().getTotalReads();  
		Iterator<Alignments> itr = segments.getAnnotationList().iterator();
		while (itr.hasNext())
			scoreSegment(itr.next(), controlTotalReads, targetTotalReads);
	}
	
	public void scoreSegment(Alignments w, double controlTotalReads, double targetTotalReads) throws IOException {
		List<Double> counts = countSegment(w);
		double targetCount = counts.get(0);
		double controlCount = counts.get(1);
		double foldDiff = (targetCount/targetTotalReads*1000000 + 1) / (controlCount/controlTotalReads*1000000 + 1);
		
		List<Double> scores = new ArrayList<Double>();
		scores.add(0.0);
		scores.add(targetCount);
		scores.add(controlCount);
		scores.add(targetTotalReads);
		scores.add(controlTotalReads);
		scores.add(targetCount / targetTotalReads * 1000000.0);
		scores.add(controlCount / controlTotalReads * 1000000.0);
		scores.add(foldDiff);
		w.setScores(scores);
	}
	
	
	public List<Double> countSegment(Alignments w) throws IOException {
		List<Double> counts = new ArrayList<Double>();
		counts.add(targetModel.count(w));
		counts.add(controlModel.count(w));
		return counts;
	}
	
	
	
	/************************************************************************************************
	 * Estimate the parameters of the gamma distribution that serves as the prior distribution
	 * for the Poisson parameter for a window.  Do this by fitting a negative binomial, then 
	 * back-calculating the parameters of the gamma prior.  Note: This gives a fairly similar result 
	 * to just estimating the gamma distribution directly, at least in parameter space.
	 * @param chr
	 * @param windows
	 * @return
	 * @throws IOException
	 */
	private GammaDist estimateGammaFromControl(List<Double> counts) throws IOException {

		int[] countArray = new int[counts.size()];
		for (int i = 0; i < counts.size(); i++) countArray[i] = counts.get(i).intValue();

		double[] params = NegativeBinomialDist.getMLE(countArray, countArray.length);
		System.out.println("Control NB Size = " + params[0] + ". Prob = " + params[1]);	
		
		double gammaShape = params[0];   // gamma shape = negbinomial size
		double gammaRate = params[1] / (1 - params[1]);  // back-calculate gamma rate from negbinomial prob
		System.out.println("Control Gamma Shape = " + gammaShape + ". Rate = " + gammaRate);
		
		return new GammaDist(gammaShape, gammaRate);
	}
	
	
	/**
	 * Old working version that fits a gamma directly instead of via negative binomial
	 * @param chr
	 * @param window
	 * @return
	 * @throws IOException
	private GammaDist estimateGammaFromControl(String chr, int window) throws IOException {
		BasicLightweightAnnotation region = new BasicLightweightAnnotation(chr, 1, chromosomeSizes.get(chr));

		int overlap = (int)Math.round((float)window * OVERLAP_FACTOR);

		// Calculate read counts in each window
		SlideAndCount slider = new SlideAndCount(controlModel, window, overlap, maskedRegions);
		SlideAndCount.SaveCountsConsumer consumer = slider.new SaveCountsConsumer();
		slider.slide(region, consumer);
		List<Integer> counts = consumer.getCounts();

		// Convert to fragments per window per million reads
		double chromosomeReads = controlModel.getNumberOfReads(chr);
		System.out.println("Control sample - number of reads on " + chr + ": " + chromosomeReads);

		double[] fpwm = new double[counts.size()];
		for (int i=0; i<counts.size(); i++) {
			fpwm[i] = counts.get(i) / chromosomeReads * 1000000;
			fpwm[i] = fpwm[i] + 0.5;  // this is necessary for fitting the gamma distribution, which doesn't like zero values 
		}
		System.out.println("counts[100] = " + counts.get(99));
		System.out.println("fpwm[100] = " + fpwm[99]);

		double[] params = GammaDist.getMLE(fpwm, fpwm.length);
		System.out.println("Window = " + window + ".  Shape = " + params[0] + ". Rate = " + params[1]);
		return new GammaDist(params[0], params[1]);
	}
	*/
	
	
	
	/************************************************************************************************
	 * FUNCTIONS FOR BUILDLING EMPIRICAL DISTRIBUTIONS FROM SCANNING COUNTS
	 * @param windows
	 * @return
	 */
	private EmpiricalDistribution edFromWindows(List<Alignments> windows) {
		List<Double> values = new ArrayList<Double>();
		for (Alignments window : windows) {
			values.add(window.getScores().get(PVALUE_INDEX));
		}
		
		EmpiricalDistribution ed = getEmptyEmpiricalDistribution();
		ed.addAll(values);
		return ed;
	}
	
	
	/**
	 * Generate an empty empirical distribution.  Use this function so that all EDs in the class
	 * have a standard number of bins, etc.
	 * @return
	 */
	public static EmpiricalDistribution getEmptyEmpiricalDistribution() {
		return new EmpiricalDistribution(20000, 0, 20);
	}
	
	
	/**
	 * Build a distribution of -log p-values by scanning windows across the entire genome.
	 * @param region
	 * @param windowSize
	 * @return
	 * @throws IOException
	 * @throws ParseException
	 */
	public EmpiricalDistribution buildDistribution(String region, int windowSize) throws Exception, IOException, ParseException {
		List<String> chromosomes = parseRegions(region);
		EmpiricalDistribution ed = null;
		for (String chr : chromosomes) {
			List<Alignments> chrResults = scanChromosome(chr, windowSize);
			EmpiricalDistribution chrDist = edFromWindows(chrResults);
			if (ed == null) {
				ed = chrDist;
			} else {
				ed.addDistribution(chrDist);
			}
		}
		return ed;
	}
	
	public Map<Integer, EmpiricalDistribution> buildDistributions(String region, int[] windows) throws Exception, IOException, ParseException {
		Map<Integer, EmpiricalDistribution> eds = new HashMap<Integer, EmpiricalDistribution>();
		for (int windowSize : windows) {
			eds.put(new Integer(windowSize), buildDistribution(region, windowSize));
		}
		return eds;
	}
	

	
	
	
	/********************************************************************************
	 * COMMAND LINE UTILITIES
	 ********************************************************************************/

	private static String USAGE = "java -jar Rapture.jar [args] TODO"; // TODO
	
	public static Rapture newFromArgs(ArgumentMap argmap) throws IOException, ParseException, Exception {
		Globals.setHeadless(true);
		String targetFile = argmap.getMandatory("target");
		String controlFile = argmap.getMandatory("control");
		//int minMappingQuality = argmap.containsKey("minMappingQuality") ? argmap.getInteger("minMappingQuality") : DEFAULT_MIN_MAPPING_QUALITY;
		int extensionFactor = argmap.getInteger("extensionFactor", 0);
		String sizes   =  argmap.getMandatory("sizeFile");
		String maskFileDir = argmap.get("maskFileDir", null);
		boolean loadPairsAsFragments = argmap.containsKey("pairedEnd") || argmap.containsKey("loadPairsAsFragments");
		Rapture rap = new Rapture(targetFile, controlFile, extensionFactor, sizes, maskFileDir, loadPairsAsFragments);
		return rap;
	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception, IOException, ParseException {
		ArgumentMap argmap = CLUtil.getParameters(args, USAGE);
		Rapture rap = newFromArgs(argmap);
		
		String out = argmap.getOutput();
		String chr = argmap.get("chr", null);
		String task = argmap.getTask();
		

		/**
		 *  TASK = "scan"
		 */
		if (task.equals("scan")) {
			System.out.println("Task = scan");
			int[] windows = ContinuousDataAlignmentModel.getWidths(argmap.getMandatory("windows"));
			rap.scanAndWriteResults(chr, windows, out);
		/**
		 *  TASK = "distribution"
		 */
		} else if (task.equals("distribution")) {
			int[] windows = ContinuousDataAlignmentModel.getWidths(argmap.getMandatory("windows"));
			System.out.println("Task = build null peak distribution");
			Map<Integer, EmpiricalDistribution> eds = rap.buildDistributions(chr, windows);
			for (Integer windowSize : eds.keySet()) {
				writeEdf(out, windowSize.intValue(), eds.get(windowSize));
			}
		/**
		 *  TASK = "peaks"
		 */
		} else if (task.equals("peaks")) {
			System.out.println("Task = call peaks");
			int[] windows = ContinuousDataAlignmentModel.getWidths(argmap.getMandatory("windows"));
			double fold = argmap.getDouble("fold", 0);         // Fold cut-off
			
			if (argmap.containsKey("null")) {
				String nullFile = argmap.getMandatory("null");
				double alpha = argmap.getDouble("alpha", 0.001);   // FDR cut-off
				rap.callPeaksWithNull(chr, windows, out, nullFile, alpha, fold);
			} else {
				double pvalue = -1*Math.log10(argmap.getDouble("pvalue", 0.0001));   // P-value cutoff with no null distribution
				rap.callPeaks(chr, windows, out, pvalue, fold);
			}
		/**
		 *  TASK = "scoresegments"
		 */
		} else if (task.equalsIgnoreCase("scoresegments")) {
			System.out.println("Task = scoreSegments");
			String in = argmap.getInput();
			rap.scoreAndWriteSegments(in, out);
		} else {
			throw new ParseException("Invalid task: " + task);
		} 
		
	}
	
	
	public void scanAndWriteResults(String chr, int[] windows, String out) throws Exception, IOException, ParseException {
		scanAndWriteResults(chr, windows, out, null);
	}
	
	public void scanAndWriteResults(String chr, int[] windows, String out, Predicate<Alignments> filter) throws Exception, IOException, ParseException {
		Map<Integer, List<Alignments>> results = scan(chr, windows, filter);
		for (Integer window : results.keySet()) {
			List<Alignments> currScores = results.get(window);
			Alignments.write(out + "." + window + "_window.bed", currScores);
		
			EmpiricalDistribution ed = edFromWindows(currScores);
			writeEdf(out, window, ed);
		}
	}
	
	public static void writeEdf(String out, int window, EmpiricalDistribution ed) throws IOException {
		ed.write(getEdfFileName(out, window));
	}
	
	public static String getEdfFileName(String base, int window) {
		return base + "." + window + "_window.edf";
	}
	
	public void callPeaksWithNull(String chr, int[] windows, String out, String nullFile, double alpha, double fold) throws Exception, IOException, ParseException {
		EmpiricalDistribution nullDistribution = new EmpiricalDistribution(new File(nullFile));
		double cutoff = nullDistribution.getQuantile(1-alpha);
		callPeaks(chr, windows, out, cutoff, fold);
	}
	
	public void callPeaks(String chr, int[] windows, String out, double cutoff, double fold) throws Exception, IOException, ParseException {
		System.out.println("PValue Cutoff = " + cutoff);
		System.out.println("Fold Cutoff = " + fold);
		Predicate<Alignments> filter = new PValueFilter(cutoff, fold);
		
		Map<Integer, List<Alignments>> results = scan(chr, windows, filter);
		
		// Merge and re-score windows
		Collection<Alignments> mergedResults = new ArrayList<Alignments>();
		for (List<Alignments> curr : results.values()) mergedResults.addAll(curr);
		
		AlignmentsReader reader = new AlignmentsReader(mergedResults);
		reader.merge();
		scoreSegments(reader);
		System.out.println(reader.size());
		reader.write(out + ".bed");
		//scanAndWriteResults(chr, windows, out, filter);
	}
	
	
	public void scoreAndWriteSegments(String in, String out) throws IOException, ParseException {
		List<BED> beds = new BEDReader(in).getAnnotationList();
		List<Alignments> alignments = new LinkedList<Alignments>();
		for (BED b : beds) {
			alignments.add(new Alignments(b));
		}
		Collections.sort(alignments);
		scoreSegments(alignments);
		Alignments.write(out, alignments);
	}
}
