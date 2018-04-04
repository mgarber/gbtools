package broad.pda.nanostring;

import broad.core.annotation.Feature;

public class Observation implements Feature{
	private String transcriptClass;
	private String name;
	private int count;
	private Experiment experiment;
	private double zScore;
	private double confidence;
	private double foldChange;
	
	public Observation(String observationClass, String observationName, int count) {
		this.count = count;
		this.name = observationName.intern();
		this.transcriptClass = observationClass.intern();
	}
	public String getTranscriptClass() {
		return transcriptClass;
	}
	public String getName() {
		return name;
	}
	public double getScore() {
		return count;
	}
	
	public int getCount() {
		return count;
	}
	
	public Experiment getExperiment() {
		return experiment;
	}
	
	public void substractCount(int maxNegativeCount) {
		int newCount = count - maxNegativeCount;
		count = newCount <= maxNegativeCount ? 0 : newCount;
	}
	
	public void addCount(int addCount) {
		count += addCount; 
	}
	
	public void setCount(int count) {
		this.count = count; 
	}
	
	public void normalizeCounts(double normalizingFactor) {
		count = (int) (count * normalizingFactor);
	}
	public void setZScore(double zScore) {
		this.zScore = zScore;
	}
	public double getZScore() {
		return zScore;
	}
	public void setConfidence(double confidence) {
		this.confidence = confidence;
	}
	public double getConfidence() {
		return confidence;
	}
	public void setFoldChange(double d) {
		this.foldChange = d;
	}
	public double getFoldChange() {
		return this.foldChange;
	}
	public void setExperiment(Experiment experiment) {
		this.experiment = experiment;
		
	}
	
	

}
