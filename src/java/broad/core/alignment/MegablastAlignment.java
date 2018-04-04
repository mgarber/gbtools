package broad.core.alignment;

public class MegablastAlignment extends AlignmentSummary {
	int length;
	int mismatches;
	int gapOpenings;
	float expectation;
	


	MegablastAlignment(String [] rawData) {
		int i = 0;
		setQuery(rawData[i++]);
		setSubject(rawData[i++]);
		setPid(Float.parseFloat(rawData[i++]));
		setLength(Integer.parseInt(rawData[i++]));
		setMismatches(Integer.parseInt(rawData[i++]));
		setGapOpenings(Integer.parseInt(rawData[i++]));
		setQueryStart(Integer.parseInt(rawData[i++]));
		setQueryEnd(Integer.parseInt(rawData[i++]));
		setSubjectStart(Integer.parseInt(rawData[i++]));
		setSubjectEnd(Integer.parseInt(rawData[i++]));
		setReversedOrientation(getSubjectStart() > getSubjectEnd());
		setExpectation(Float.parseFloat(rawData[i++]));
		//setScore(Integer.parseInt(rawData[i++]));
	}
	
	public float getExpectation() {
		return expectation;
	}

	public int getGapOpenings() {
		return gapOpenings;
	}

	public int getLength() {
		return length;
	}

	public int getMismatches() {
		return mismatches;
	}

	protected void setExpectation(float expectation) {
		this.expectation = expectation;
	}

	protected void setLength(int length) {
		this.length = length;
	}

	protected void setMismatches(int mismatches) {
		this.mismatches = mismatches;
	}
	
	
}
