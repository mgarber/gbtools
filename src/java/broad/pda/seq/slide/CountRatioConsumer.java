package broad.pda.seq.slide;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import broad.core.annotation.LightweightGenomicAnnotation;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.alignment.AlignmentCollection;

public class CountRatioConsumer extends SlideAndCountConsumer {
	protected AlignmentCollection denominator;
	private int dTotalReads;
	public final static double RPKM_OFFSET = 0.1; // this is really Reads per Window per Million
	
	public CountRatioConsumer(AlignmentCollection numerator, AlignmentCollection denominator) throws IOException {
		super(numerator);
		this.denominator = denominator;
		this.dTotalReads = denominator.size();
	}
	
	public AlignmentCollection getNumerator() { return d; }
	public AlignmentCollection getDenominator() { return denominator; }
	
	public int getNumeratorTotalReads() { return getTotalReads(); }
	public int getDenominatorTotalReads() { return dTotalReads; }
	
	@Override
	public int getWindowBatch() { return 10000; }
	
	@Override
	public void initRegion(LightweightGenomicAnnotation region) {}  // override SlideAndCountConsumer
	
	@Override 
	public void consume(List<Alignments> windows) throws Exception {
		for (Alignments windowToScan : windows) {
			double nCount = d.getCount(windowToScan);
			double dCount = denominator.getCount(windowToScan);
			double nFPKM = asFPKM(nCount, getNumeratorTotalReads()) + RPKM_OFFSET;
			double dFPKM = asFPKM(dCount,  getDenominatorTotalReads()) + RPKM_OFFSET;
			double ratio = nFPKM / dFPKM;
			windowToScan.setScore(ratio);
			windowToScan.setCountScore(ratio);
			
			List<Double> scores = new ArrayList<Double>();
			scores.add(nCount);
			scores.add(dCount);
			scores.add((double) getNumeratorTotalReads());
			scores.add((double) getDenominatorTotalReads());
			scores.add(nFPKM);
			scores.add(dFPKM);
			windowToScan.setScores(scores);
			
			getScores().add(ratio);
		}
	}
	
	private double asFPKM(double count, int total) {
		return count / total * 1000000.0;
	}
}
