package broad.pda.seq.slide;

import java.util.ArrayList;
import java.util.List;

import broad.pda.datastructures.Alignments;
import broad.pda.seq.alignment.AlignmentCollection;

public class CountBasesCoveredConsumer extends SlideConsumer.AbstractSlideConsumer {
	protected AlignmentCollection d;
	List<Double> basesCovered = new ArrayList<Double>();;
	
	public CountBasesCoveredConsumer(AlignmentCollection d) {
		this.d =d;
	}
	
	@Override
	public void consume(List<Alignments> windows) throws Exception {
		for (Alignments windowToScan: windows) {
			int count = d.getBasesCovered(windowToScan);
			windowToScan.setScore(count);
			basesCovered.add((double)count);
		}
	}
	
	@Override
	public List<Double> getScores() { return basesCovered; }
}
