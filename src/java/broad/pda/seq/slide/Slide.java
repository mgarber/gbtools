/**
 * 
 */
package broad.pda.seq.slide;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.annotation.ShortBEDReader;
import broad.core.error.ParseException;
import broad.pda.datastructures.Alignments;
import net.sf.picard.util.Log;

/**
 * Jesse Engreitz
 * August 28, 2012
 * Class for sliding across the genome.
 * Optional: skipping regions specified by a mask file
 * TODO:  Incorporate a CoordinateSpace class somehow - maybe that can take care of sliding as well
 * @author engreitz
 */
public class Slide {
    private static final Log log = Log.getInstance(Slide.class);
	
	private int window;
	private int overlap;
	private ShortBEDReader maskedRegions = null;
	private boolean overlapAllowed = false;
	private double pctOverlapAllowed = 0.0;
	ArrayList<String> windowNames = new ArrayList<String>();
	
	public Slide(int window) {
		this(window, (int) Math.floor(window/2));  // default to windows overlapping by 50%
	}

	public Slide(int window, int overlap) {
		if (window < 1 || overlap < 1)
			throw new IllegalArgumentException("Overlap and window should be greater than 0");
		if(overlap > window) {
			throw new IllegalArgumentException("Overlap was larger than window size: overlap " + overlap + ", window " + window);
		}
		this.window = window;
		this.overlap = overlap;
	}
	
	public Slide(int window, int overlap, String maskFileDir, String chr) throws ParseException, IOException, IllegalArgumentException {
		this(window, overlap);
		File[] files = new File(maskFileDir).listFiles();
		this.maskedRegions = loadMaskedRegions(files, chr);
	}

	public Slide(int window, int overlap, String maskFile) throws IOException, ParseException, IllegalArgumentException {
		this(window, overlap);
		this.maskedRegions = loadMaskedRegions(new File(maskFile));
	}
	
	
	public Slide(int window, int overlap, ShortBEDReader maskedRegions) {
		this(window, overlap);
		this.maskedRegions = maskedRegions;
	}
	
	/**
	 * @param window				size of windows (nucleotides)
	 * @param overlap				overlap of windows (nucleotides)
	 * @param maskedRegions			regions to mask
	 * @param pctOverlapAllowed		windows with greater than this % overlap with masked regions will be skipped
	 */
	public Slide(int window, int overlap, ShortBEDReader maskedRegions, double pctOverlapAllowed) {
		this(window, overlap);
		this.maskedRegions = maskedRegions;
		this.pctOverlapAllowed = pctOverlapAllowed;
		this.overlapAllowed = Math.abs(pctOverlapAllowed - 0.0) > 0.0001;
	}
	
	
	/**
	 * Method to scan a region, skipping any windows that overlap masked regions
	 * @param region		region to scan
	 * @param consumer		class to process each window
	 * @throws Exception
	 */
	public void slide(LightweightGenomicAnnotation region, SlideConsumer consumer) throws Exception {
		int pos = region.getStart();
		int cacheWindows = consumer.getWindowBatch();
		List<Alignments> windowList = new ArrayList<Alignments>();  // important that windows are accessible in order
		windowNames = new ArrayList<String>();
		
		consumer.initRegion(region);
		log.info("Starting slide on " + region.toUCSC());
		
		int counter = 0;
		do {
			Alignments windowToScan = new Alignments(region.getChromosome(), pos, pos + window);

			// process if window does not overlap masked region
			if (acceptWindow(windowToScan)) {   
				windowList.add(windowToScan);
				windowNames.add(windowToScan.toUCSC());
				counter++;
				if (windowList.size() == cacheWindows) {
					consumer.consume(windowList);
					windowList.clear();
				}
			}
			
			pos = pos + window - overlap;
		} while(pos + window <= region.getEnd());
		
		consumer.consume(windowList);
		consumer.finishedRegion();
		log.info(counter + " windows accepted.");
	}
	
	public void slide(List<? extends LightweightGenomicAnnotation> regions, SlideConsumer consumer) throws Exception {
		for (LightweightGenomicAnnotation region : regions) 
			slide(region, consumer);
	}
	
	public boolean acceptWindow(Alignments window) throws IOException {
		// Skip windows that overlap with masked regions
		boolean reject = false;
		if (maskedRegions != null) {
			if (!overlapAllowed) {
				reject = (maskedRegions.getOverlappers(window).size() > 0);
			} else {
				double basesCovered = (double) maskedRegions.getBasesCovered(window);
				reject = (basesCovered / (double) window.length())*100.0 > pctOverlapAllowed;
			}
		}
		return !reject;
	}
	
	public final ArrayList<String> getWindowNames() { return windowNames; }
	
	/**
	 * @param file				file to load
	 * @return					AnnotationReader containing masked regions
	 * @throws IOException
	 * @throws ParseException
	 */
	public static ShortBEDReader loadMaskedRegions(File file) throws IOException, ParseException {
		if (file == null) return null;
		return new ShortBEDReader(file.getAbsolutePath());		
	}
	

	/**
	 * @param files				List of files in the mask file directory
	 * @param chr				chromosome to load
	 * @return					AnnotationReader containing masked regions
	 * @throws IOException
	 * @throws ParseException
	 */
	public static ShortBEDReader loadMaskedRegions(File[] files, String chr) throws IOException, ParseException {
		if (files == null) return null;
		
		for (int i=0; i < files.length; i++) {
			String currentChr = files[i].getName().split("\\.")[0];
			if (currentChr.equals(chr)) {
				return new ShortBEDReader(files[i].getAbsolutePath());
			}
		} 
		
		return null; // throw new IOException("Mask file for " + chr + " not found.");
	}
	
	
	private static class EmptyConsumer extends SlideConsumer.AbstractSlideConsumer {
		@Override
		public void consume(List<Alignments> windows) {};
	}
	
	
	public static void main(final String[] args) throws Exception {
		int window = Integer.parseInt(args[0]);
		int overlap = Integer.parseInt(args[1]);
		double pctAllowed = Double.parseDouble(args[2]);
		ShortBEDReader maskedRegions = loadMaskedRegions(new File("/seq/mguttman/ChIPData/MaskFiles/MM9Segments/all.mask.mouse.n36.d2.bin"));
		Slide slider = new Slide(window, overlap, maskedRegions, pctAllowed);
		slider.slide(new Alignments("chr1",1,10000000), new EmptyConsumer());
	}
}
