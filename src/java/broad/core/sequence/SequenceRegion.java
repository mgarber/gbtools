package broad.core.sequence;

import java.util.Collection;
import java.util.List;

import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.annotation.TwoSubjectAnnotation;

public class SequenceRegion extends Sequence implements GenomicAnnotation{
	BasicGenomicAnnotation annotation;
	private String containingSequenceId;
	public static final int INF = 1000000000;
	
	public SequenceRegion(String containingSequenceId) {
		super(null);
		this.containingSequenceId = containingSequenceId;
		annotation = new BasicGenomicAnnotation(containingSequenceId);
		annotation.setSequence(this);
		annotation.setEnd(INF);
	}
	public SequenceRegion(String containingSequenceId, LightweightGenomicAnnotation annotation) {
		super(annotation.getName());
		this.containingSequenceId = containingSequenceId;
		this.annotation = new BasicGenomicAnnotation(annotation);
		
	}
	public int getRegionEnd() {
		return annotation.getEnd() == INF && super.getSequenceBases() != null && getSequenceBases().length() > 0 
			? getSequenceBases().length() + annotation.getStart() - 1
			: annotation.getEnd();
	}
	public void setRegionEnd(int regionEnd) {
		annotation.setEnd(regionEnd);
	}
	public int getRegionStart() {
		return annotation.getStart();
	}
	public void setRegionStart(int regionStart) {
		annotation.setStart(regionStart);
	}
	
	public WindowSlider getSlider(int windowSize, int overlap) {
		return WindowSlider.getSlider(this, windowSize, overlap);
	}
	
	public int length() {
		if(getEnd() == INF && getSequenceBases() == null && getSequenceBases().length() == 0) {
			return INF;
		} else if(getEnd() < INF) {
			return getEnd() - getStart();
		} else {
			return super.getLength();
		}
	}
	
	public String getOrientation() { return inReversedOrientation() ? "-" : "+"; }
	public void setReversedOrientation(boolean reversed) { annotation.setOrientation(!reversed); }
	
	public String getContainingSequenceId() {
		return containingSequenceId;
	}
	
	public String getId() {
		return super.getId() == null ? containingSequenceId + ":" + getStart() + "-" + getEnd() : super.getId();
	}
	
	public int getStart() {
		return annotation.getStart();
	}
	
	public int getEnd() {
		return annotation.getEnd();
	}
	
	public double getScore() {
		return annotation.getScore();
	}
	
	public void setScore(double score) {
		annotation.setScore(score);
	}
	
	public double getExtraScore(int i ) {
		return annotation.getExtraScore(i);
	}
	
	
	public String getName() {
		return getId();
	}
	
	public void setId(String id) {
		super.setId(id);
		annotation.setName(id);
	}
	
	public void setName(String name) {
		setId(name);
	}
	
	public boolean inReversedOrientation() {
		return annotation.inReversedOrientation();
	}
	
	public long getMiddle() {
		return annotation.getMiddle();
	}
	
	public void setStart(int start) {
		setRegionStart(start);
	}
	
	public void setEnd(int end) {
		setRegionEnd(end);
	}
	public Sequence getSequence() {
		return this;
	}
	public void setSequence(Sequence seq) {
		super.setSequenceBases(seq.getSequenceBases());
	}
	
	public String getChromosome() {
		return annotation.getChromosome();
	}
	
	public void setChromosome(String chr) {
		annotation.setChromosome(chr);
	}
	public boolean overlaps(LightweightGenomicAnnotation other, int buffer) {
		return annotation.overlaps(other, buffer);
	}
	public boolean overlaps(LightweightGenomicAnnotation other) {
		return annotation.overlaps(other);
	}
	
	public int getOverlap(LightweightGenomicAnnotation other) { return annotation.getOverlap(other);}
	
	public List<GenomicAnnotation> minus(GenomicAnnotation other) {
		return annotation.minus(other);
	}
	
	public List<GenomicAnnotation> minus(List<? extends GenomicAnnotation> others) {
		return annotation.minus(others);
	}
	
	public void takeIntersection(LightweightGenomicAnnotation other) {
		annotation.takeIntersection(other);
	}
	
	public void takeUnion(LightweightGenomicAnnotation other) {
		annotation.takeUnion(other);
	}
	public void stitchTo(LightweightGenomicAnnotation other) {
		annotation.stitchTo(other);
	}
	
	public boolean isFlankedBy(TwoSubjectAnnotation twoSubjectAnnotation, int buffer) {
		return annotation.isFlankedBy(twoSubjectAnnotation, buffer);
	}
	public int getFivePrimeBases() {
		return annotation.getFivePrimeBases();
	}
	public int getThreePrimeBases() {
		return annotation.getThreePrimeBases();
	}
	
	public String toString() {
		return getContainingSequenceId() + ":" + getRegionStart() + "-" + getRegionEnd();
	}
	
	public int getAbsoluteStart(int absoluteStart){
		return getRegionStart()+absoluteStart;
	}
	
	public int getAbsoluteEnd(int absoluteStart){
		return getRegionEnd()+absoluteStart;
	}
	
	public String getLocationString() { return annotation.getLocationString();}
	public String getChromosomeString() { return annotation.getChromosomeString(); }
	
	public SequenceRegion extractRegionBasedOnGC(float targetGC, int size, int buffer) {
		SequenceRegion theRegion = super.extractRegionBasedOnGC(targetGC, size, buffer);
		if(theRegion != null) {
			theRegion.setRegionStart(getRegionStart() + theRegion.getRegionStart());
			theRegion.setRegionEnd(getRegionStart() + theRegion.getRegionEnd());
			theRegion.setChromosome(annotation.getChromosome());
		}
		return theRegion;
	}
	public boolean contains(LightweightGenomicAnnotation other) {
		return annotation.contains(other);
	}
	public int getDistanceTo(LightweightGenomicAnnotation other) {
		return annotation.getDistanceTo(other);
	}
	public void setOrientation(String orientation) {
		annotation.setOrientation(orientation);
		
	}
	public int compareTo(GenomicAnnotation arg0) {
		return annotation.compareTo(arg0);
	}
	public List<GenomicAnnotation> disect(GenomicAnnotation a) {
		return annotation.disect(a);
	}
	public List<GenomicAnnotation> disect(List<? extends GenomicAnnotation> disectors) {
		return annotation.disect(disectors);
	}
	
	public SequenceRegion getRegion(int start, int end) {
		SequenceRegion region = super.getRegion(start, end);
		if(annotation.getChromosome() != null) {
			region.setChromosome(annotation.getChromosome());
		}
		return region;
	}
	public int getOrientedEnd() {
		return annotation.getOrientedEnd();
	}
	public int getOrientedStart() {
		return annotation.getOrientedStart();
	}
	public void addBlock(String name, int start, int end) {
		//DO nothing as there is nothing to do here.
		
	}
	public List<? extends GenomicAnnotation> getBlocks() {
		return null;
	}
	public boolean mayHaveBlocks() {
		return false;
	}
	
	public String toUCSC() {
		return annotation.toUCSC();
	}
	
	public int getLength() {
		return length();
	}


	public boolean overlaps(Collection<LightweightGenomicAnnotation> others,
			int buffer) {
		return annotation.overlaps(others, buffer);
	}


	public boolean overlaps(Collection<LightweightGenomicAnnotation> others) {
		return annotation.overlaps(others);
	}
	
	public void removeExtraScores() {
		annotation.removeExtraScores();
	}

	public int compareTo(LightweightGenomicAnnotation b) {
		// TODO Auto-generated method stub
		return annotation.compareTo(b);
	}

	public LightweightGenomicAnnotation intersect(LightweightGenomicAnnotation other) {

		return annotation.intersect(other);
	}

	public List<LightweightGenomicAnnotation> intersect(List<? extends LightweightGenomicAnnotation> annotations) {
		return annotation.intersect(annotations);
	}
	public double percentGC() {
		char[] seqChar=this.getSequenceBases().toCharArray();
		double GCCount=0;
		for(int i=0; i<seqChar.length; i++){
			if(seqChar[i]=='G' || seqChar[i]=='C' || seqChar[i]=='g' || seqChar[i]=='c'){GCCount++;}
		}
		double percentGC=GCCount/seqChar.length;
		return percentGC;
	}

	public void addExtraScore(double score) {
		annotation.addExtraScore(score);
		
	}

	public List<Double> getExtraScores() {
		return annotation.getExtraScores();
	}
	
	/**
	 * Checks whether two annotations differ by a small (fudge) factor
	 * @param fudge Maximum difference at either end or start to consider similar
	 * @return
	 */
	public boolean almostEqual(LightweightGenomicAnnotation other, int fudge) {
		return annotation.almostEqual(other, fudge);
	}
	
	
}
