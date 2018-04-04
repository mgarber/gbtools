package broad.pda.seq.slide;

import java.util.List;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math.MathException;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;

import broad.pda.datastructures.Alignments;
import broad.pda.seq.alignment.AlignmentCollection;

public class CorrelationConsumer extends SlideConsumer.AbstractSlideConsumer {
	private PearsonsCorrelation p = null;
	private SlideConsumer c1, c2;
	
	/**
	 * @param d
	 * @param annotation
	 */
	public CorrelationConsumer(AlignmentCollection collection1, AlignmentCollection collection2) {
		this(new SlideAndCountConsumer(collection1), new SlideAndCountConsumer(collection2));
	}
	
	public CorrelationConsumer(SlideConsumer c1, SlideConsumer c2) {
		this.c1 = c1;
		this.c2 = c2;
	}
	
	@Override
	public void consume(List<Alignments> windows) throws Exception {
		c1.consume(windows);  // this will save scores in the list
		c2.consume(windows);  // this will overwrite scores in windows but save in list
	}
	
	public double getScore() {
		return getCorrelation();
	}
	
	public double getCorrelation() {
		if (p == null) calculateCorrelation();
		return p.getCorrelationMatrix().getEntry(1,0);
	}
	
	public double getPValue() throws MathException {
		if (p == null) calculateCorrelation();
		return p.getCorrelationPValues().getEntry(1,0);
	}
	
	private void calculateCorrelation() {
		RealMatrix matrix = formMatrix(c1.getScores(), c2.getScores());
		p = new PearsonsCorrelation(matrix);
	}
	
	/*
	public double getCountCorrelation() {
		RealMatrix m = formMatrix(counts, annotCounts);
		PearsonsCorrelation p = new PearsonsCorrelation(m);
		return p.getCorrelationMatrix().getEntry(0,0);
	}
	
	public double getCoverageCorrelation() {
		RealMatrix m = formMatrix(counts, basesCovered);
		PearsonsCorrelation p = new PearsonsCorrelation(m);
		return p.getCorrelationMatrix().getEntry(0,0);
	}
	
	public double getCountCorrelationPValue() throws MathException {
		RealMatrix m = formMatrix(counts, annotCounts);
		PearsonsCorrelation p = new PearsonsCorrelation(m);
		return p.getCorrelationPValues().getEntry(0,0);
	}
	
	public double getCoverageCorrelationPValue() throws MathException {
		RealMatrix m = formMatrix(counts, basesCovered);
		PearsonsCorrelation p = new PearsonsCorrelation(m);
		return p.getCorrelationPValues().getEntry(0,0);
	}
	*/
	
	public static RealMatrix formMatrix(List<Double> l1, List<Double> l2) {
		if (l1.size() != l2.size()) throw new RuntimeException("Two lists must be the same length");
		double[][] d = new double[2][l1.size()];
		d[0] = ArrayUtils.toPrimitive(l1.toArray(new Double[l1.size()]));
		d[1] = ArrayUtils.toPrimitive(l2.toArray(new Double[l2.size()]));
		return new Array2DRowRealMatrix(d).transpose();
	}
	
}
