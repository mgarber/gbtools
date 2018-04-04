package broad.pda.rap;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math.MathException;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;

import Jama.Matrix;
import broad.core.annotation.AnnotationReader;
import broad.core.annotation.AnnotationReaderFactory;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.annotation.ShortBEDReader;
import broad.core.datastructures.MatrixWithHeaders;
import broad.pda.datastructures.Alignments;
import broad.pda.enrichment.EnrichmentMetric;
import broad.pda.enrichment.EnrichmentUtils;
import broad.pda.seq.alignment.AlignmentCollection;
import broad.pda.seq.slide.CorrelationConsumer;
import broad.pda.seq.slide.CountBasesCoveredConsumer;
import broad.pda.seq.slide.CountRatioConsumer;
import broad.pda.seq.slide.DistanceToClosestConsumer;
import broad.pda.seq.slide.SavedSlideConsumer;
import broad.pda.seq.slide.ScoreSumConsumer;
import broad.pda.seq.slide.Slide;
import broad.pda.seq.slide.SlideAndCount;
import broad.pda.seq.slide.SlideAndCountConsumer;
import broad.pda.seq.slide.SlideConsumer;
import jsc.independentsamples.SmirnovTest;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;

public class CollectAnnotationCorrelations extends RaptureCommandLineProgram {
    private static final Log log = Log.getInstance(CollectAnnotationCorrelations.class);
	
    @Usage
    public String USAGE = "Runs enrichment tests on a RAP experiment.";
    
	@Option(doc="Control SAM or BAM file for normalization.")
	public File CONTROL;
	
	@Option(doc="File or directory containing genomic annotations (bed or bedGraph format)")
	public File ANNOTATION_DIR;
	
	@Option(doc="Window size")
	public int WINDOW;

	@Option(doc="Overlap between windows")
	public int OVERLAP;
	
	@Option(doc="Ratio file to read/write")
	public File RATIO_FILE = null;
	
	@Option(doc="Percent masked allowable per sliding window", optional=true)
	public double PCT_MASKED_ALLOWED = 20.0;
	
	@Option(doc="Name of the analysis for writing to text file")
	public String ANALYSIS_NAME;
	
	@Option(doc="Number of permutations")
	public int PERMUTATIONS = 1000;
	
	private static final String[] DISCRETE_EXTENSIONS = {"bed"};
	private static final String[] CONTINUOUS_EXTENSIONS = {"bedGraph","bedgraph"};
	
	Slide slider;
	
	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new CollectAnnotationCorrelations().instanceMain(args));
	}
	
	@Override
	protected String[] customCommandLineValidation() {
		// Add custom validation here.
		return super.customCommandLineValidation();
	}
	
	@Override
	protected int doWork() {
		
		try {
			IoUtil.assertFileIsReadable(INPUT);
			IoUtil.assertFileIsReadable(CONTROL);
			IoUtil.assertFileIsWritable(OUTPUT);
			IoUtil.assertDirectoryIsReadable(ANNOTATION_DIR);
			
			Map<String,File> discreteAnnotations = EnrichmentUtils.findAnnotationFiles(ANNOTATION_DIR, DISCRETE_EXTENSIONS);
			Map<String,File> continuousAnnotations = EnrichmentUtils.findAnnotationFiles(ANNOTATION_DIR, CONTINUOUS_EXTENSIONS);
			log.info(continuousAnnotations.size());
			
			slider = new Slide(WINDOW, OVERLAP, getMaskedRegions(), PCT_MASKED_ALLOWED);
			
			// Load continuous data alignment models
			AlignmentCollection target = getContinuousDataModel(INPUT);
			AlignmentCollection control = getContinuousDataModel(CONTROL);
			
			List<GenomicAnnotation> regions = getRegions();
			List<EnrichmentMetric> results = new ArrayList<EnrichmentMetric>();
			
			log.info("Calculating sliding window ratios");
			List<Double> ratios = getAlignmentRatios(target, control, regions);
			
			Map<String, List<Double>> slidingValues = new LinkedHashMap<String, List<Double>>();
			slidingValues.put("ratios", ratios);
			
			
			// DISCRETE ANNOTATIONS
			for (Map.Entry<String, File> entry : discreteAnnotations.entrySet()) {
				log.info("Processing " + entry.getKey());
				double memoryPercent=Runtime.getRuntime().freeMemory()/(double)Runtime.getRuntime().totalMemory();
				log.info("% Memory Used: " + memoryPercent);
				
				// eventually will want to generalize this to accept other types
				AnnotationReader<? extends GenomicAnnotation> reader = AnnotationReaderFactory.create(entry.getValue().getAbsolutePath(), "BED");
				reader.merge();
				
				int numInRegions = reader.getOverlappers(regions).size();
				
				List<Double> densities = getDensities(reader, regions);
				List<Double> baseCoverages = getBaseCoverages(reader, regions);
				
				// now calculate correlations
				double[] densityCorrelation = getCorrelation(ratios, densities);
				EnrichmentMetric m = new EnrichmentMetric(entry.getKey(), "continuous", "discrete", "correlation", "density", densityCorrelation[0], densityCorrelation[1], reader.size(), numInRegions);
				results.add(m);
				slidingValues.put(m.getId(), densities);
				// TODO: Calculate FDR by permutation rather than significance of correlation
				
				double[] baseCorrelation = getCorrelation(ratios, baseCoverages);
				m = new EnrichmentMetric(entry.getKey(), "continuous", "discrete", "correlation", "baseCoverage", baseCorrelation[0], baseCorrelation[1], reader.size(), numInRegions);
				results.add(m);
				slidingValues.put(m.getId(), baseCoverages);
				
				// (useful for sparse BED files): for each window calculate the distance to the closest element and correlate based on this
				List<Double> distances = getDistances(reader, regions);
				double[] distanceCorrelation = getCorrelation(ratios, distances);
				m = new EnrichmentMetric(entry.getKey(), "continuous", "discrete", "correlation", "distanceToClosest", distanceCorrelation[0], distanceCorrelation[1], reader.size(), numInRegions);
				results.add(m);
				slidingValues.put(m.getId(), distances);
				
				// now calculate KS tests of permuted regions
				//double[] ks = getKSTest(reader, target, control, regions);
				//m = new EnrichmentMetric(entry.getKey(), "continuous", "discrete", "KS", "ratios", ks[0], ks[1], reader.size(), numInRegions);
				//results.add(m);
				// TODO: VERY memory intensive at the moment.  Need to try another implementation	
			}
			
			// CONTINUOUS ANNOTATIONS
			for (Map.Entry<String, File> entry : continuousAnnotations.entrySet()) {
				log.info("Processing " + entry.getKey());
				
				// eventually will want to generalize this to accept other types
				// Read BEDGraph file
				AnnotationReader<? extends GenomicAnnotation> reader = AnnotationReaderFactory.create(entry.getValue().getAbsolutePath(), "BEDGraph");
				
				int numInRegions = reader.getOverlappers(regions).size();
				
				// calculate correlations
				List<Double> scores = getSumScores(reader, regions);
				double[] correlation = getCorrelation(ratios, scores);
				results.add(new EnrichmentMetric(entry.getKey(), "continuous","continuous","correlation","sum",correlation[0], correlation[1], reader.size(), numInRegions));

			}
			
			// Output results
			BufferedWriter bw = new BufferedWriter(new FileWriter(OUTPUT));
			for (EnrichmentMetric result : results) 
				bw.write(ANALYSIS_NAME + "\t" + result.toString() + "\n");
			bw.close();
			
			writeResultsMatrix(slidingValues, slider.getWindowNames());
			
		} catch (Exception e) {
			log.error(e);
		}
		
		return 0;
	}
	
	
	private void writeResultsMatrix(Map<String, List<Double>> results, List<String> windowNames) throws IOException {
		//List<List<Double>> values = new ArrayList<List<Double>>();
		ArrayList<String> colNames = new ArrayList<String>();
		
		int ncol = results.size();
		int nrow = windowNames.size();
		
		double [][] values = new double[nrow][ncol];
		
		Iterator<String> itr = results.keySet().iterator();
		for (int i = 0; i < results.size(); i++) {
			String key = itr.next();
			List<Double> curr = results.get(key);
			for (int j = 0; j < curr.size(); j++) {
				values[j][i] = curr.get(j);
			}
			colNames.add(key);
		}

		Matrix matrix = new Matrix(values);
		MatrixWithHeaders m = new MatrixWithHeaders(matrix, windowNames, colNames);
		m.write(OUTPUT.getAbsolutePath() + ".matrix.txt");
	}
	
	
	private double[] getCorrelation(List<Double> l1, List<Double> l2) throws MathException {
		RealMatrix matrix = CorrelationConsumer.formMatrix(l1, l2);
		PearsonsCorrelation p = new PearsonsCorrelation(matrix);
		double[] results = new double[2];
		results[0] = p.getCorrelationMatrix().getEntry(1,0);
		results[1] = p.getCorrelationPValues().getEntry(1,0);
		return results;
	}
	
	
	private double[] getKSTest(AnnotationReader<? extends GenomicAnnotation> reader, AlignmentCollection target, AlignmentCollection control, List<GenomicAnnotation> regions) throws Exception {
		// Note:  this test only makes sense for doing whole chromosomes
		List<Double> regionRatios = getRatiosForRegions(reader, target, control, regions);
		List<Double> nullRatios = new ArrayList<Double>();
		
		ShortBEDReader shuffleReader = reader.cloneRegions();
		
		for (int i = 0; i < PERMUTATIONS; i++) {
			shuffleReader.permuteAnnotationsOnChromosome(sizes);
			nullRatios.addAll(getRatiosForRegions(shuffleReader, target, control, regions));
		}
		
		
		double[] d = new double[regionRatios.size()];
		for (int i=0; i<regionRatios.size();i++) d[i] = regionRatios.get(i).doubleValue();
		double[] e = new double[nullRatios.size()];
		for (int i=0; i<nullRatios.size();i++) e[i] = nullRatios.get(i).doubleValue();
		
		double[] results;
		try {
			SmirnovTest test = new SmirnovTest(d, e);		
			results = new double[] {test.getTestStatistic(), test.getSP()};
		} catch (RuntimeException x) {
			results = new double[] {0,1};
		}
		return results;
	}
	
	private List<Double> getRatiosForRegions(AnnotationReader<? extends GenomicAnnotation> reader, AlignmentCollection target, AlignmentCollection control, List<GenomicAnnotation> regions) throws Exception {
		// WARNING:  this is very slow at the moment.
		// create list of alignment
		List<Alignments> alignments = new ArrayList<Alignments>();
		for (LightweightGenomicAnnotation region : regions) {
			List<? extends GenomicAnnotation> overlappers = reader.getOverlappers(region);
			for (GenomicAnnotation annot : overlappers)
				alignments.add(new Alignments(annot));
		}
		Collections.sort(alignments);
		
		CountRatioConsumer consumer = new CountRatioConsumer(target, control);
		consumer.consume(alignments);
		return consumer.getScores();
	}
	
	
	/**
	 * TODO:  move this to a SlideUtils class
	 * @param target
	 * @param control
	 * @param regions
	 * @return
	 * @throws Exception
	 */
	private List<Double> getAlignmentRatios(AlignmentCollection target, AlignmentCollection control, List<GenomicAnnotation> regions) throws Exception {
		SlideConsumer ratioConsumer;
		
		if (RATIO_FILE != null) {
			if (!RATIO_FILE.exists()) {
				// Calculate and save ratios
				SlideConsumer printer = new SlideAndCount.PrintRatioConsumer(target, control, new BufferedWriter(new FileWriter(RATIO_FILE)));
				slider.slide(regions, printer);
			} 
			// Read ratios back in the slider
			ratioConsumer = new SavedSlideConsumer(RATIO_FILE);				
		} else {
			// otherwise, just calculate on the fly
			ratioConsumer = new CountRatioConsumer(target, control);
		}
			
		slider.slide(regions, ratioConsumer);
		return ratioConsumer.getScores();
	}
	
	
	private List<Double> getDensities(AlignmentCollection reader, List<GenomicAnnotation> regions) throws Exception {
		SlideConsumer counter = new SlideAndCountConsumer(reader);
		slider.slide(regions, counter);
		return counter.getScores();
	}
	
	
	private List<Double> getBaseCoverages(AlignmentCollection reader, List<GenomicAnnotation> regions) throws Exception {
		SlideConsumer counter = new CountBasesCoveredConsumer(reader);
		slider.slide(regions, counter);
		return counter.getScores();
	}
	
	
	private List<Double> getSumScores(AnnotationReader<? extends GenomicAnnotation> reader, List<GenomicAnnotation> regions) throws Exception {
		SlideConsumer scores = new ScoreSumConsumer(reader);
		slider.slide(regions, scores);
		return scores.getScores();
	}
	
	private List<Double> getDistances(AnnotationReader<? extends GenomicAnnotation> reader, List<GenomicAnnotation> regions) throws Exception {
		SlideConsumer distances = new DistanceToClosestConsumer(reader);
		slider.slide(regions, distances);
		return distances.getScores();
	}

}
