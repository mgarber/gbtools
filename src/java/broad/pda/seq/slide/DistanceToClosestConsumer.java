package broad.pda.seq.slide;

import java.util.ArrayList;
import java.util.List;

import broad.core.annotation.AnnotationReader;
import broad.core.annotation.GenomicAnnotation;
import broad.pda.datastructures.Alignments;

public class DistanceToClosestConsumer extends SlideConsumer.AbstractSlideConsumer {
	protected AnnotationReader<? extends GenomicAnnotation> d;
	List<Double> scores = new ArrayList<Double>();
	
	public DistanceToClosestConsumer(AnnotationReader<? extends GenomicAnnotation> d) {
		this.d =d;
	}
	
	@Override
	public void consume(List<Alignments> windows) throws Exception {
		for (Alignments windowToScan: windows) {
			double score = (double) d.getDistanceToClosest(windowToScan);
			windowToScan.setScore(score);
			scores.add(score);
		}
	}
	
	@Override
	public List<Double> getScores() { return scores; }
}
