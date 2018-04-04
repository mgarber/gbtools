package broad.pda.enrichment;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import broad.core.annotation.AnnotationReader;
import broad.core.annotation.AnnotationReaderFactory;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.ShortBEDReader;
import broad.pda.rap.GenomeCommandLineProgram;
import jsc.independentsamples.SmirnovTest;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;

public class CollectAnnotationEnrichments extends GenomeCommandLineProgram {
	private static final Log log = Log.getInstance(CollectAnnotationEnrichments.class);
	
    @Usage
    public String USAGE = "Runs enrichment tests on a set of regions in a BED-3+ file";
    
	@Option(doc="Input BED file (at least 3 columns: chr	start	end)", shortName="I")
	public File INPUT;
	
	@Option(doc="File to write output results", shortName="O")
	public File OUTPUT;
	
	@Option(doc="File or directory containing genomic annotations (bed or bedGraph format)")
	public File ANNOTATION_DIR;
	
	@Option(doc="Name of the analysis for writing to text file")
	public String ANALYSIS_NAME;
	
	@Option(doc="Number of permutations")
	public int PERMUTATIONS = 1000;
	
	private static final String[] DISCRETE_EXTENSIONS = {"bed"};
	private static final String[] CONTINUOUS_EXTENSIONS = {"bedGraph"};
	
	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new CollectAnnotationEnrichments().instanceMain(args));
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
			IoUtil.assertFileIsWritable(OUTPUT);
			IoUtil.assertDirectoryIsReadable(ANNOTATION_DIR);
			
			Map<String,File> discreteAnnotations = EnrichmentUtils.findAnnotationFiles(ANNOTATION_DIR, DISCRETE_EXTENSIONS);
			Map<String,File> continuousAnnotations = EnrichmentUtils.findAnnotationFiles(ANNOTATION_DIR, CONTINUOUS_EXTENSIONS);
			
			List<GenomicAnnotation> regions = getRegions();
			
			ShortBEDReader peaks = null;
			peaks = new ShortBEDReader(INPUT.getAbsolutePath());
			List<EnrichmentMetric> results = new ArrayList<EnrichmentMetric>();
			
			// TODO:  Should also find nearby genes and write them out for GSEA analysis

			// DISCRETE ANNOTATIONS
			for (Map.Entry<String, File> entry : discreteAnnotations.entrySet()) {
				log.info("Processing " + entry.getKey());
				// eventually will want to generalize this to accept other types
				AnnotationReader<? extends GenomicAnnotation> reader = AnnotationReaderFactory.create(entry.getValue().getAbsolutePath(), "BED");
				reader.merge();
				int numInRegions = reader.getOverlappers(regions).size();
				
				double[] ksDistance = getKSDistance(reader, peaks, regions);
				EnrichmentMetric m = new EnrichmentMetric(entry.getKey(), "discrete", "discrete", "KS", "distanceToClosest", ksDistance[0], ksDistance[1], reader.size(), numInRegions);
				results.add(m);

				double[] overlap = getOverlapFDR(reader, peaks, regions);
				m = new EnrichmentMetric(entry.getKey(), "discrete", "discrete", "q", "overlap", overlap[0], overlap[1], reader.size(), numInRegions);
				results.add(m);
			}
			
			// CONTINUOUS ANNOTATIONS
			for (Map.Entry<String, File> entry : continuousAnnotations.entrySet()) {
				//log.info("Processing " + entry.getKey());				
				//TODO
			}
			
			// Output results
			log.info("Writing " + results.size() + " tests to " + OUTPUT.getPath());
			if (OUTPUT.exists()) OUTPUT.delete();
			BufferedWriter bw = new BufferedWriter(new FileWriter(OUTPUT));
			for (EnrichmentMetric result : results) 
				bw.write(ANALYSIS_NAME + "\t" + result.toString() + "\n");
			bw.close();
			
		} catch (Exception e) {
			log.error(e);
		}
		
		return 0;
	}
	
	

	
	private double[] getKSDistance(AnnotationReader<? extends GenomicAnnotation> reader, AnnotationReader<? extends GenomicAnnotation> peaks, List<GenomicAnnotation> regions) throws IOException {
		// Note: this test only makes sense for doing whole chromosomes, since the AnnotationReader permutes across whole chromosomes
		AnnotationReader<? extends GenomicAnnotation> regionPeaks = peaks.cloneRegions();
		regionPeaks.filterByOverlap(regions);
		
		Collection<Integer> realDistances = regionPeaks.getDistancesToClosest(reader).values();
		List<Integer> nullDistances = new ArrayList<Integer>();
		ShortBEDReader shuffleReader = reader.cloneRegions();

		for (int i = 0; i < PERMUTATIONS; i++) {
			shuffleReader.permuteAnnotationsOnChromosome(sizes);
			Collection<Integer> permutationDistances = regionPeaks.getDistancesToClosest(shuffleReader).values();
			nullDistances.addAll(permutationDistances);
		}
	
		
		double[] d = new double[realDistances.size()];
		for (int i=0; i<realDistances.size();i++) d[i] = new ArrayList<Integer>(realDistances).get(i).doubleValue();
		double[] e = new double[nullDistances.size()];
		for (int i=0; i<nullDistances.size();i++) e[i] = nullDistances.get(i).doubleValue();
		
		double[] results;
		try {
			SmirnovTest test = new SmirnovTest(d, e);		
			results = new double[] {test.getTestStatistic(), test.getSP()};
		} catch (RuntimeException x) {
			results = new double[] {0,1};
		}
		return results;
	}
	
	
	private double[] getOverlapFDR(AnnotationReader<? extends GenomicAnnotation> reader, AnnotationReader<? extends GenomicAnnotation> peaks, List<GenomicAnnotation> regions) {
		AnnotationReader<? extends GenomicAnnotation> regionPeaks = peaks.cloneRegions();
		regionPeaks.filterByOverlap(regions);
		
		int overlappers = peaks.getOverlappers(reader).size();
		List<Integer> nullOverlappers = new ArrayList<Integer>();
		
		ShortBEDReader shuffleReader = reader.cloneRegions();
		int countEqualOrMore = 0;
		for (int i = 0; i < PERMUTATIONS; i++) {
			shuffleReader.permuteAnnotationsOnChromosome(sizes);
			int currOverlap = regionPeaks.getOverlappers(shuffleReader).size();
			nullOverlappers.add(currOverlap);
			if (currOverlap >= overlappers) countEqualOrMore++;
		}
		
		double[] results = new double[2];
		results[0] = overlappers;
		results[1] = (double) countEqualOrMore / (double) nullOverlappers.size();
		return results;
	}
	
}
