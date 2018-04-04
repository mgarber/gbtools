package broad.pda.seq.slide;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.annotation.ShortBED;
import broad.core.annotation.ShortBEDReader;
import broad.pda.datastructures.Alignments;
import net.sf.picard.util.Log;

/**
 * @author engreitz
 * Use this class to load saved sliding data output by other sliders to avoiding recalculating.
 * e.g. Load data saved by SlideAndCount.PrintRatioConsumer
 */
public class SavedSlideConsumer extends SlideConsumer.AbstractSlideConsumer {
	private static final Log log = Log.getInstance(SavedSlideConsumer.class);
	private ShortBEDReader reader = null;
	private Iterator<ShortBED> itr = null;
	protected List<Double> counts = new ArrayList<Double>();
	
	public SavedSlideConsumer(File file) throws IOException {
		this(file.getAbsolutePath());
	}
	
	public SavedSlideConsumer(String filename) throws IOException {
		reader = new ShortBEDReader(filename);
		reset();
	}
	
	private void reset() {
		itr = reader.getAnnotationList().iterator();
	}
	
	@Override
	public void initRegion(LightweightGenomicAnnotation region) {
		 reset(); // send iterator to beginning
	}
	
	@Override 
	public int getWindowBatch() { return 100000; }
	
	@Override
	public void consume(List<Alignments> windows) throws Exception {
		// scan to find first window
		ShortBED curr = null;
		while (itr.hasNext()) {
			ShortBED next = itr.next();
			if (alignmentsEqual(windows.get(0), next)) {
				curr = next;
				break;
			}
		}
		
		for (Alignments w : windows) {
			if (curr == null || !alignmentsEqual(w, curr)) {
				throw new IllegalStateException("SavedSlideConsumer iterator is in illegal state - are you scanning with the same parameters now?");
			}
			
			w.setScore(curr.getScore());
			counts.add(curr.getScore());
			
			if (itr.hasNext()) {
				curr = itr.next();
			} else {
				curr = null;
			}
		}
	}
	
	private boolean alignmentsEqual(Alignments window, ShortBED bed) {
		return (window.getChromosome().equals(bed.getChromosome()) &&
				window.getStart() == bed.getStart() &&
				window.getEnd() == bed.getEnd());
	}
	
	@Override
	public List<Double> getScores() { return counts; }
	
	
	public static void main(String[] args) throws IOException, Exception {
		log.info("For testing purposes only ...");
		SavedSlideConsumer consumer = new SavedSlideConsumer("/seq/lincRNA/RAP/XISTTimecourse/120826_EZH2_ChIP_HiSeq/aligned/scan/6hr_vs_WCE.W100000_O75000.bed");
		Slide slider = new Slide(100000, 75000, Slide.loadMaskedRegions(new File("/seq/mguttman/ChIPData/MaskFiles/MM9Segments/all.mask.mouse.n36.d2.bin")), 20.0);
		slider.slide(new Alignments("chrX", 1, 166650296), consumer);
		System.out.println(consumer.getScores().size());
	}
	
}
