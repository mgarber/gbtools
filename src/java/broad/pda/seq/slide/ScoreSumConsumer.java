package broad.pda.seq.slide;

import java.util.ArrayList;
import java.util.List;

import broad.core.annotation.AnnotationReader;
import broad.core.annotation.GenomicAnnotation;
import broad.pda.datastructures.Alignments;


public class ScoreSumConsumer extends SlideConsumer.AbstractSlideConsumer {
	protected AnnotationReader<? extends GenomicAnnotation> d;
	List<Double> scores = new ArrayList<Double>();
	
	public ScoreSumConsumer(AnnotationReader<? extends GenomicAnnotation> d) {
		this.d =d;
	}
	
	@Override
	public void consume(List<Alignments> windows) throws Exception {
		for (Alignments windowToScan: windows) {
			double score = d.getScoreSum(windowToScan);
			windowToScan.setScore(score);
			scores.add(score);
		}
	}
	
	@Override
	public List<Double> getScores() { return scores; }
}
