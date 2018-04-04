package broad.core.alignment;

import java.text.DecimalFormat;

import broad.core.annotation.BasicTwoSubjectAnnotation;


public class AlignmentSummary extends BasicTwoSubjectAnnotation {
	float pid;
	float score;
	boolean reversedAlignment;
	private static final DecimalFormat numFormatter = new DecimalFormat("#0.###");
	private int alignmentLength;
	private int missmatches;
	private int gapOpenings;
	private float eValue;

	
	public float getEValue() {
		return eValue;
	}

	public void setEValue(float value) {
		eValue = value;
	}

	public int getGapOpenings() {
		return gapOpenings;
	}

	public void setGapOpenings(int gapOpenings) {
		this.gapOpenings = gapOpenings;
	}

	public int getMissmatches() {
		return missmatches;
	}

	public void setMissmatches(int missmatches) {
		this.missmatches = missmatches;
	}

	protected void setPid(float pid) {
		this.pid = pid;
	}

	public void setQuery(String query) {
		getA().setName(query);
	}

	public void setQueryEnd(int queryEnd) {			
		getA().setEnd(queryEnd);
	}

	protected void setQueryOrientation(String queryOrientation) {
		getA().setOrientation(queryOrientation);
	}

	public void setQueryStart(int queryStart) {
		getA().setStart(queryStart);
	}

	public void setScore(float score) {
		this.score = score;
	}

	public void setSubject(String subject) {
		getB().setName(subject);
	}

	public void setSubjectEnd(int subjectEnd) {
		getB().setEnd(subjectEnd);
	}

	protected void setSubjectOrientation(String subjectOrientation) {
		getB().setOrientation(subjectOrientation);
	}

	public void setSubjectStart(int subjectStart) {
		getB().setStart(subjectStart);
	}
	
	public AlignmentSummary() {}

	public AlignmentSummary(String [] rawData) {
		super();
		pid = rawData[0].length() == 0 ? 0 : Float.parseFloat(rawData[0]);
		setQuery(rawData[1]);
		setSubject(rawData[2]);
		score = rawData[3].length() == 0 ? 0 : Integer.parseInt(rawData[3]);
		int qStart = Integer.parseInt(rawData[6]);
		int qEnd   = Integer.parseInt(rawData[7]);
		if( qStart <= qEnd) {
			setQueryStart(qStart);
			setQueryEnd(qEnd);
		} else {
			setQueryStart(qEnd);
			setQueryEnd(qStart);
		}
		setQueryOrientation(rawData[5]);
		
		setSubjectOrientation(rawData[9]);
		int sStart = Integer.parseInt(rawData[10]);
		int sEnd   = Integer.parseInt(rawData[11]);
		if(sStart <= sEnd) {	
			setSubjectStart(sStart);
			setSubjectEnd(sEnd);
		} else {
			setSubjectStart(sEnd);
			setSubjectEnd(sStart);
		}
	}

	public String toString() {
		StringBuffer buf = new StringBuffer(numFormatter.format(getPid()));
		buf.append("\t").append(getQuery())
			.append("\t").append(getSubject())
			.append("\t").append(getScore())
			.append("\t").append(getA().getOrientation())
			.append("\t").append(getA().getStart())
			.append("\t").append(getA().getEnd())
			.append("\t").append(getB().getOrientation())
			.append("\t").append(getB().getStart())
			.append("\t").append(getB().getEnd());
		
		return buf.toString();
	}

	public float getPid() {
		return pid;
	}

	public String getQuery() {
		return getAName();
	}
	
	public int getQueryEnd() {
		return getA().getEnd();
	}
	
	public int getAEnd(){
		return getA().getEnd();
	}

	public String getQueryOrientation() {
		return getA().getOrientation();
	}

	public int getQueryStart() {
		return getA().getStart();
	}
	
	public int getAStart() {
		return getA().getStart();
	}

	public float getScore() {
		return score;
	}

	public String getSubject() {
		return getB().getName();
	}
	

	public int getSubjectEnd() {
		return getB().getEnd();
	}
	
	public String getSubjectOrientation() {
		return getB().getOrientation();
	}

	public int getSubjectStart() {
		return getB().getStart();
	}
	
	public boolean isReversedOrientation() {
		return reversedAlignment;
	}
	
	public void setReversedOrientation(boolean isReversedOrientation) {
		reversedAlignment = isReversedOrientation;
	}

	public int getAlignmentLength() {
		return alignmentLength;
	}

	public void setAlignmentLength(int alignmentLength) {
		this.alignmentLength = alignmentLength;
	}

	public void merge(AlignmentSummary novel) {
		stitch(novel);
	}
	
	public boolean isStitchable(AlignmentSummary other, int subjectLength) {
		boolean aOverlap = getA().overlaps(other.getA());
		boolean bOverlap = getB().overlaps(other.getB());
		//System.out.println("aOverlap " + aOverlap + ", bOverlap " + bOverlap + " a==b " + (aOverlap==bOverlap) +
		//		" isDirect " + isDirect() + " other.isDirect " + other.isDirect() + " isdirect==other.isDirect " + (isDirect() == other.isDirect()));
		return isDirect() == other.isDirect() 
				&& (getA().getStart() - other.getA().getEnd()) < subjectLength
				&& ( aOverlap == bOverlap);
	}
	
	public void stitch(AlignmentSummary other) {
		int nqStart = Math.min(other.getA().getStart(), getA().getStart());
		int nqEnd   = Math.max(other.getA().getEnd(), getA().getEnd());
		
		int nsStart = Math.min(other.getB().getStart(), getB().getStart());
		int nsEnd   = Math.max(other.getB().getEnd(), getB().getEnd());
		
		int nGapOpen = getGapOpenings() + other.getGapOpenings() + getA().getDistanceTo(other.getA());
		
		int combinedLength = getA().length() + other.getA().length();
		float nPID = (getPid() * getA().length() + other.getPid() * other.getA().length())/(float)combinedLength;
		
		int nMissmatches = getMissmatches() + other.getMissmatches();
		float nEValue     = getEValue() * other.getEValue();
		float nScore      = getScore() + other.getScore();
		int nALignmentLength = nqEnd - nqStart;
		
		//System.out.println ("\tBefore stitching " + toString());
		
		getA().setStart(nqStart);
		getA().setEnd(nqEnd);
		
		getB().setStart(nsStart);
		getB().setEnd(nsEnd);
		
		setGapOpenings(nGapOpen);
		setAlignmentLength(nALignmentLength);
		setMissmatches(nMissmatches);
		setPid(nPID);
		setEValue(nEValue);
		setScore(nScore);
		
		//System.out.println("\t after " + toString());
		
	}
}