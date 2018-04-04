package broad.pda.seq.slide;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import broad.core.annotation.LightweightGenomicAnnotation;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.alignment.AlignmentCollection;


/**
 * Consumer class to have multiple consumers "listen" to the counting results.
 * Use this when you want to have multiple modular calculations on the same count
 * results, but you want to only count once.
 * @author engreitz
 *
 */
public class MultiListenerConsumer extends SlideAndCountConsumer {
	private List<SlideAndCountConsumer> listeners;
	public MultiListenerConsumer(AlignmentCollection d) throws IOException {
		super(d);
		listeners = new ArrayList<SlideAndCountConsumer>();
	}

	public List<SlideAndCountConsumer> getListeners() { return listeners; }
	public void addListener(SlideAndCountConsumer listener) { 
		listeners.add(listener);
		listener.setSkipCount(true);
	}
	
	@Override
	public void initRegion(final LightweightGenomicAnnotation region) {
		super.initRegion(region);
		for (SlideAndCountConsumer c : listeners) {
			c.setRegionReads(getRegionReads());
			c.setTotalReads(getTotalReads());
		}
	}

	@Override
	public void consume(List<Alignments> windows) throws Exception {
		super.consume(windows);         // count the windows in the super class
		for (SlideAndCountConsumer listener : listeners) {
			listener.consume(windows);  // pass the counts to all of the listeners.  
			 							// don't modify the scores in the listeners ... way to encode that?
		}
	}
	
	@Override
	public void finishedRegion() {
		for (SlideAndCountConsumer listener : listeners) 
			listener.finishedRegion();
	}
	
	@Override
	public List<Double> getScores() { throw new UnsupportedOperationException(); }  // override this
}
