/**
 * 
 */
package broad.pda.seq.slide;

import java.util.List;

import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.annotation.ShortBEDReader;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.segmentation.MultipleDataAlignmentModel;
import net.sf.picard.util.Log;

/**
 * @author engreitz
 * Use for sliding over multiple data alignment models at once and calculating an aggregated result from many models
 * Advise against using this for the time being - untested
 */
@Deprecated
public class SlideAndCountMultiple extends Slide {
    private static final Log log = Log.getInstance(SlideAndCountMultiple.class);
    
	public SlideAndCountMultiple(int window, int overlap) {
		super(window, overlap);
	}
	
	public SlideAndCountMultiple(int window, int overlap, ShortBEDReader maskedRegions) {
		super(window, overlap, maskedRegions);
	}
	
	
	public void slide(LightweightGenomicAnnotation region, MultipleConsumer consumer) throws Exception {
		Alignments alignment = new Alignments(region);
		consumer.initRegion(alignment);
		super.slide(region, consumer);
	}
	
	private class MultipleConsumer extends SlideConsumer.AbstractSlideConsumer {
		private List<Integer> regionReads, totalReads;
		MultipleDataAlignmentModel data;
		
		MultipleConsumer(MultipleDataAlignmentModel d) { 
			data = d; 
		}
		
		public void initRegion(LightweightGenomicAnnotation region) {
			try {
				setRegionReads(data.getCount(new Alignments(region)));
				setTotalReads(data.getTotalReads());
			} catch (Exception e) {
				// TODO
			}
		}
		
		public void setRegionReads(List<Integer> regionReads) {this.regionReads = regionReads;}
		public void setTotalReads(List<Integer> totalReads) {this.totalReads = totalReads;}
		public List<Integer> getRegionReads() { return regionReads; }
		public List<Integer> getTotalReads() { return totalReads; }
		
		@Override
		public int getWindowBatch() { return 10000; }
		
		@Override
		public void consume(List<Alignments> windows) throws Exception {
			for (Alignments windowToScan : windows) {
				List<Integer> count = data.getCount(windowToScan);
				windowToScan.setCountScores(count);
			}
		}
	}
}
