package broad.pda.seq.slide;

import java.util.List;

import broad.core.annotation.LightweightGenomicAnnotation;
import broad.pda.datastructures.Alignments;

public interface SlideConsumer {
	/**
	 * This function is called by Slide when a new region is about to be processed (i.e. about to slide many windows across the region)
	 * @param region
	 */
	public void initRegion(final LightweightGenomicAnnotation region);		
	
	/**
	 * This function is called for each batch of windows created by Slide.
	 * Use this function to calculate stats, count reads, etc.
	 * @param windows
	 * @throws Exception
	 */
	public void consume(List<Alignments> windows) throws Exception;   
	
	/**
	 * Sets the size of the window batch to process at once
	 * @return
	 */
	public int getWindowBatch();
	
	
	/**
	 * Called after the entire region has been processed
	 */
	public void finishedRegion();
	
	
	/**
	 * Access scores that are generated during sliding (could be optional).
	 * @return
	 */
	public List<Double> getScores();
	
	
	/**
	 * @author engreitz
	 * Extend this class for basic functionality.
	 */
	public abstract class AbstractSlideConsumer implements SlideConsumer {
		public void initRegion(final LightweightGenomicAnnotation region) {}
		public int getWindowBatch() { return 100000; }
		public List<Double> getScores() { throw new UnsupportedOperationException(); }
		public void finishedRegion() {}
	}
}
