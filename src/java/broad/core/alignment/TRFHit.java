package broad.core.alignment;

import broad.core.sequence.SequenceRegion;

public class TRFHit extends SequenceRegion {
	private double percentA;
	private double percentC;
	private double percentG;
	private double percentT;
	
	private int period;
	private float numberOfCopies;
	private int consensusSize;
	private double percentMatches;
	private double percentIndels;
	
	//alignment score would be used as the hit score;
	private double entropy;
	//consensus would be used as name;
	//Actual hit would be stored as sequence.
	
	public TRFHit (String containingSequence, String [] rawData) {
		super(containingSequence);
		setName(rawData[13]);
		setScore(Double.parseDouble(rawData[7]));
		setStart(Integer.parseInt(rawData[0]));
		setEnd(Integer.parseInt(rawData[1]));
		setSequenceBases(rawData[14]);
		
		int i = 2;
		this.period = Integer.parseInt(rawData[i++]);
		this.numberOfCopies = Float.parseFloat(rawData[i++]);
		this.consensusSize = Integer.parseInt(rawData[i++]);
		this.percentMatches = Double.parseDouble(rawData[i++]);
		this.percentIndels  = Double.parseDouble(rawData[i++]);
		
		i = 8;
		this.percentA = Double.parseDouble(rawData[i++]);
		this.percentC = Double.parseDouble(rawData[i++]);
		this.percentG = Double.parseDouble(rawData[i++]);
		this.percentT = Double.parseDouble(rawData[i++]);
		
		this.entropy = Double.parseDouble(rawData[i]);
	}

	public double getPercentA() {
		return percentA;
	}

	public void setPercentA(double percentA) {
		this.percentA = percentA;
	}

	public double getPercentC() {
		return percentC;
	}

	public void setPercentC(double percentC) {
		this.percentC = percentC;
	}

	public double getPercentG() {
		return percentG;
	}

	public void setPercentG(double percentG) {
		this.percentG = percentG;
	}

	public double getPercentT() {
		return percentT;
	}

	public void setPercentT(double percentT) {
		this.percentT = percentT;
	}

	public int getPeriod() {
		return period;
	}

	public void setPeriod(int period) {
		this.period = period;
	}

	public float getNumberOfCopies() {
		return numberOfCopies;
	}

	public void setNumberOfCopies(float numberOfCopies) {
		this.numberOfCopies = numberOfCopies;
	}

	public int getConsensusSize() {
		return consensusSize;
	}

	public void setConsensusSize(int consensusSize) {
		this.consensusSize = consensusSize;
	}

	public double getPercentMatches() {
		return percentMatches;
	}

	public void setPercentMatches(double percentMatches) {
		this.percentMatches = percentMatches;
	}

	public double getPercentIndels() {
		return percentIndels;
	}

	public void setPercentIndels(double percentIndels) {
		this.percentIndels = percentIndels;
	}

	public double getEntropy() {
		return entropy;
	}

	public void setEntropy(double entropy) {
		this.entropy = entropy;
	}
	

}
