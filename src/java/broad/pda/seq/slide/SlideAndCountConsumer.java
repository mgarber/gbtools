package broad.pda.seq.slide;

import java.util.ArrayList;
import java.util.List;

import broad.core.annotation.LightweightGenomicAnnotation;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.alignment.AlignmentCollection;
import net.sf.picard.util.Log;

/**
 * @author engreitz
 * For use with SlideAndCount
 */
public class SlideAndCountConsumer extends SlideConsumer.AbstractSlideConsumer {
    private static final Log log = Log.getInstance(SlideAndCountConsumer.class);
    
	private int regionReads, totalReads;	// useful metrics to have for most SlideAndCount applications
	protected AlignmentCollection d;
	protected List<Double> counts = new ArrayList<Double>();
	
	// flag to actually skip the counting part in the event that the consumer is part of a MultiListenerConsumer
	private boolean skipCount = false;           

	public SlideAndCountConsumer(AlignmentCollection d) {
		this.d = d;
		try {
			setTotalReads(d.size());
		} catch (Exception e) {
			log.error(e);
		}
	}
	
	public SlideAndCountConsumer(AlignmentCollection d, boolean skipCount) {
		this(d);
		this.skipCount = skipCount;
	}
	
	public AlignmentCollection getData() { return d; }
	
	@Override
	public void initRegion(LightweightGenomicAnnotation region) {
		try {
			setRegionReads(d.getCount(new Alignments(region)));
		} catch (Exception e) {
			log.error(e);
		}
	}
	
	public void setRegionReads(int regionReads) {this.regionReads = regionReads;}
	public void setTotalReads(int totalReads) {this.totalReads = totalReads;}
	public int getRegionReads() { return regionReads; }
	public int getTotalReads() { return totalReads; }
	public void setSkipCount(boolean b) { skipCount = b; }
	public boolean getSkipCount() { return skipCount; }
	
	@Override
	public int getWindowBatch() {return 100000;}
	
	@Override
	public void consume(List<Alignments> windows) throws Exception {
		if (!skipCount) {
			for (Alignments windowToScan : windows) {
				double count = d.getCount(windowToScan);
				counts.add(count);
				windowToScan.setCountScore(count);
			}
		}
	}
	
	@Override
	public List<Double> getScores() { return counts; }

}
