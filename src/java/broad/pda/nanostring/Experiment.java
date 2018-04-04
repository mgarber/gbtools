package broad.pda.nanostring;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;

import broad.core.math.Statistics;

public class Experiment {

	static final String POSITIVE_CONTROL_SPIKEIN_CLASS = "Positive";
	static final String NEGATIVE_CONTROL_SPIKEIN_CLASS = "Negative";
	private String fileAttributes;
	private String fileName;
	private String id;
	private String owner;
	private int date;
	private float fileVersion;
	private String geneRLF;
	private String comments;
	private int laneId;
	private int fovCount;
	private int fovCounted;
	private String scannerId;
	private int stagePosition;

	private LinkedHashMap<String, Observation> spikeinControls;
	private LinkedHashMap<String, Observation> spikeinNegativeControls;
	private LinkedHashMap<String, Observation> normalizationTranscripts;
	private LinkedHashMap<String, Observation> observations;
	
	public Experiment(String fileAttrs, String fileName) {
		this (fileAttrs, fileName, null, null, null, null, null, null, null, null, null, null, null);
	}
	public Experiment(String fileAttrs, String fileName, String id, String owner,
			String date, String fileVersion, String geneRLF, String comments,
			String laneId, String fovCount, String fovCounted, String scannerId, String stagePosition) {

		this.fileAttributes = fileAttrs;
		this.fileName = fileName;
		this.id = id;
		this.owner = owner != null ? owner.intern() : null;
		this.date = date == null ? 0 : Integer.parseInt(date);
		this.fileVersion = fileVersion != null ? Float.parseFloat(fileVersion) : 0;
		this.geneRLF     = geneRLF != null ? geneRLF.intern() : null;
		this.comments = comments;
		this.laneId = laneId != null ? Integer.parseInt(laneId) : 0;
		this.fovCount = fovCount != null ? Integer.parseInt(fovCount) : 0;
		this.fovCounted = fovCounted != null ? Integer.parseInt(fovCounted) : 0;
		this.scannerId = scannerId;
		this.stagePosition = stagePosition != null ? Integer.parseInt(stagePosition) : 0;
		
		spikeinControls         = new LinkedHashMap<String, Observation>();
		spikeinNegativeControls = new LinkedHashMap<String, Observation>();
		normalizationTranscripts    = new LinkedHashMap<String, Observation>();
		observations            = new LinkedHashMap<String, Observation>();
		
	}

	public String getFileName() {
		return fileName;
	}
	
	public String getName() {
		return fileAttributes == null ? fileName : fileAttributes;
	}
	
	public void setName(String name) {
		fileAttributes = name;
		
	}
	
	public String getFileAttributes() {
		return fileAttributes;
	}

	public String getId() {
		return id;
	}

	public String getOwner() {
		return owner;
	}

	public int getDate() {
		return date;
	}

	public float getFileVersion() {
		return fileVersion;
	}

	public String getGeneRLF() {
		return geneRLF;
	}

	public String getComments() {
		return comments;
	}

	public int getLaneId() {
		return laneId;
	}

	public int getFovCount() {
		return fovCount;
	}

	public int getFovCounted() {
		return fovCounted;
	}

	public String getScannerId() {
		return scannerId;
	}

	public int getStagePosition() {
		return stagePosition;
	}

	public Collection<Observation> getPositiveSpikeControls() {
		return spikeinControls.values();
	}

	public Collection<Observation> getNegativeSpikeControls() {
		return spikeinNegativeControls.values();
	}
	
	public Collection<Observation> observations() {
		return observations != null && !observations.isEmpty() ? observations.values() : new ArrayList<Observation>();
	}
	
	public Collection<Observation> experimentAndControlObservations() {
		ArrayList<Observation> all = new ArrayList<Observation>(observations.size() + normalizationTranscripts.size());
		all.addAll(normalizationTranscripts.values());
		all.addAll(observations.values());
		return all;
	}

	public void addObservation(String observationClass, String observationName, int value) {
		Observation observation = new Observation(observationClass, observationName, value);
		observation.setExperiment(this);
		if(POSITIVE_CONTROL_SPIKEIN_CLASS.equals(observationClass)) {
			spikeinControls.put(observation.getName(), observation);
		}else if(NEGATIVE_CONTROL_SPIKEIN_CLASS.equals(observationClass)) {
			spikeinNegativeControls.put(observation.getName(), observation);
		} else {
			observations.put(observation.getName(), observation);
		}
	}

	public void normalizeByNegativeSpikeIns() {
		int maxNegativeCount = getMaxObservation(getNegativeSpikeControls());
		
		if(observations != null && observations.size() > 0) {
			Iterator<Observation> it = observations.values().iterator();
			while(it.hasNext()) {
				Observation obs = it.next();
				obs.substractCount(maxNegativeCount);
			}
		
		}
		
		if(spikeinControls != null && !spikeinControls.isEmpty()) {
			Iterator<Observation> it = spikeinControls.values().iterator();
			while(it.hasNext()) {
				Observation control = it.next();
				control.substractCount(maxNegativeCount);
			}			
		}
		
		if(normalizationTranscripts != null && !normalizationTranscripts.isEmpty()) {
			Iterator<Observation> it = normalizationTranscripts.values().iterator();
			while(it.hasNext()) {
				Observation control = it.next();
				control.substractCount(maxNegativeCount);
			}			
		}
	}
	
	public void floor() {
		int maxPositiveCount = getMinObservation(getPositiveSpikeControls());
		//int  maxNegativeCount = getMaxObservation(getNegativeSpikeControls());
		if(observations != null && observations.size() > 0) { 
			Iterator<Observation> it = observations.values().iterator();
			while(it.hasNext()) {
				Observation obs = it.next();
				if(obs.getCount() <= maxPositiveCount) {
					obs.setCount(maxPositiveCount); 
				}
			}
		
		}
		
	
	}

	private int getMaxObservation(Collection<Observation> observations) {
		int max = 0;
		if(observations != null && observations.size() > 0) {
			Iterator<Observation> it = observations.iterator();
			while(it.hasNext()) {
				Observation obs = it.next();
				if(obs.getScore() > max) {max = obs.getCount();}
			}
		}
		return max;
		
	}
	
	private int getMinObservation(Collection<Observation> observations) {
		int min = Integer.MAX_VALUE;
		if(observations != null && observations.size() > 0) {
			Iterator<Observation> it = observations.iterator();
			while(it.hasNext()) {
				Observation obs = it.next();
				if(obs.getScore() < min) {min = obs.getCount();}
			}
		}
		
		return Math.max(min, 2);
		
	}

	public int getSumPositiveSpikeIns() {
		int total = 0;
		if(spikeinControls != null && !spikeinControls.isEmpty()) {
			Iterator<Observation> it = spikeinControls.values().iterator();
			while(it.hasNext()) {
				Observation control = it.next();
				total = total + control.getCount();
			}			
		}
		return total;
	}
	
	public int getSumNormalizationTranscripts() {
		int total = 0;
		if(normalizationTranscripts != null && !normalizationTranscripts.isEmpty()) {
			Iterator<Observation> it = normalizationTranscripts.values().iterator();
			while(it.hasNext()) {
				Observation control = it.next();
				total = total + control.getCount();
			}			
		}
		return total;
	}

	public void normalizeToPositiveControlCount(int refCounts) {
		int totalPositiveControlCounts = getSumPositiveSpikeIns();
		double normalizingFactor = refCounts/(double)totalPositiveControlCounts;
		if(observations != null && !observations.isEmpty()) {
			Iterator<Observation> it = observations.values().iterator();
			while(it.hasNext()) {
				Observation obs = it.next();
				obs.normalizeCounts(normalizingFactor);
			}
		}	
		
		if(normalizationTranscripts != null && !normalizationTranscripts.isEmpty()) {
			Iterator<Observation> it = normalizationTranscripts.values().iterator();
			while(it.hasNext()) {
				Observation control = it.next();
				control.normalizeCounts(normalizingFactor);
			}			
		}
		
	}
	
	public void normalizeToTranscripts(Collection<Observation>  normalizalizingObservations) {
		double sumRatios = 0;
		for(Observation refObs : normalizalizingObservations) {
			sumRatios += refObs.getCount()/(double)getObservation(refObs.getName()).getCount();
		}
		double normalizationFactor = sumRatios/normalizalizingObservations.size();
		//System.err.println("Experiment " + getFileName() + " Normalization factor: " + normalizationFactor);
		if(observations != null && !observations.isEmpty()) {
			for(Observation obs : observations.values()) {
				obs.normalizeCounts(normalizationFactor);
			}
		}
		for(Observation obs : normalizationTranscripts.values()) {
			obs.normalizeCounts(normalizationFactor);
		}
		
	}
	
	public void normalizeToNormalizationTranscriptCount(int refCounts) {
		int totalControlCounts = getSumNormalizationTranscripts();
		if(totalControlCounts > 0) {
			double normalizingFactor = refCounts/(double)totalControlCounts;
			if(observations != null && !observations.isEmpty()) {
				Iterator<Observation> it = observations.values().iterator();
				while(it.hasNext()) {
					Observation obs = it.next();
					obs.normalizeCounts(normalizingFactor);
				}
			}	
			Iterator<Observation> it = normalizationTranscripts.values().iterator();
			while(it.hasNext()) {
				Observation obs = it.next();
				obs.normalizeCounts(normalizingFactor);
			}
		}
	}

	public Observation getObservation(String name) {
		return observations.get(name) == null ? normalizationTranscripts.get(name) : observations.get(name);
	}

	public void setNormalizationTranscripts(List<String> controlTranscripts) throws IllegalStateException {
		if(controlTranscripts != null && observations != null) {
			Iterator<String> it = controlTranscripts.iterator();
			while(it.hasNext()) {
				String control = it.next();
				Observation controlObs = observations.get(control);
				if(controlObs == null) {
					throw new IllegalStateException("Control transcript " + control + " is not part of experiment " + getFileName());
				}
				observations.remove(control);
				normalizationTranscripts.put(control, controlObs);
			}
		}
		
	}

	public Collection<Observation> getNormalizationTranscriptsObservations() {
		return normalizationTranscripts.values();
	}

	public Observation getNormalizationTranscriptObservation(String name) {
		return normalizationTranscripts.get(name);
	}


	public void computeConfidenceValues(Collection<Experiment> controls, boolean byRow) {
		List<List<Double>> permuttedPositiveZScores = new ArrayList<List<Double>>(controls.size());
		List<List<Double>> permuttedNegativeZScores = new ArrayList<List<Double>>(controls.size());
		for(int i = 0; i < controls.size(); i++) {
			List<Double> permutationPZScoreList = new ArrayList<Double>(observations.size());
			permuttedPositiveZScores.add(permutationPZScoreList);
			List<Double> permutationNZScoreList = new ArrayList<Double>(observations.size());
			permuttedNegativeZScores.add(permutationNZScoreList);
		}
		Iterator<Observation> observationIt = observations.values().iterator();
		List<Double> observationPositiveZScores = new ArrayList<Double>(observations.size());
		List<Double> observationNegativeZScores = new ArrayList<Double>(observations.size());
		while(observationIt.hasNext()) {
			Observation obs = observationIt.next();
			if(obs.getZScore() >= 0) {
				observationPositiveZScores.add(obs.getZScore());
			} else {
				observationNegativeZScores.add(-obs.getZScore());
			}

			double [] controlCounts = getControlExperimentObservations(obs.getName(), controls);
			double [] permutedZScoresArray = Statistics.computePermutedZScores(obs.getCount(), controlCounts);
			for(int i = 0; i < controls.size(); i++) {
				if(permutedZScoresArray[i] > 0) {
					permuttedPositiveZScores.get(i).add(permutedZScoresArray[i]);
				} else {
					permuttedNegativeZScores.get(i).add(-permutedZScoresArray[i]);
				}
			}
		}
		
		Comparator<Double> comparator = new ReversedDoubleComparator();
		Collections.sort(observationPositiveZScores, comparator);
		Collections.sort(observationNegativeZScores, comparator);
		
		for(int i = 0; i < permuttedPositiveZScores.size(); i++) {
			Collections.sort(permuttedPositiveZScores.get(i), comparator);
			Collections.sort(permuttedNegativeZScores.get(i), comparator);
		}
		
		observationIt = observations.values().iterator();
		while(observationIt.hasNext()) {
			Observation obs = observationIt.next();
			double score = obs.getZScore();
						
			int numOfObservationsAboveZScore = 0;
			if(score < 0) {
				numOfObservationsAboveZScore = countObservationsLargerThan(-score, observationNegativeZScores );
			} else {
				numOfObservationsAboveZScore = countObservationsLargerThan(score, observationPositiveZScores );
			}
			int meanPermuttedObservationsAboveZScore = 0;
			for(int i = 0; i < permuttedPositiveZScores.size(); i++) {
				if(score < 0) {
					meanPermuttedObservationsAboveZScore = meanPermuttedObservationsAboveZScore + countObservationsLargerThan(-score, permuttedNegativeZScores.get(i));
				} else {
					meanPermuttedObservationsAboveZScore = meanPermuttedObservationsAboveZScore + countObservationsLargerThan(score, permuttedPositiveZScores.get(i));
				}
			}
			double fdr = meanPermuttedObservationsAboveZScore/((double) (numOfObservationsAboveZScore * permuttedPositiveZScores.size())); 
			obs.setConfidence(fdr >=1 ? 0 : 1 - fdr);
		}
	}
	
	/**
	 * Assumes sorted list!
	 * @param score
	 * @param list
	 * @return
	 */
	private int countObservationsLargerThan(double score, List<Double> list) {
		int i;
		for(i = 0; i < list.size(); i++) {
			if(score > list.get(i)) {break;}
		}
		return i;
	}

	private double[] getControlExperimentObservations(String name, Collection<Experiment> controls) {
		double [] obsValues = new double [controls.size()];
		Iterator<Experiment> it = controls.iterator();
		int i = 0;
		while(it.hasNext()) {
			Experiment control = it.next();
			Observation obs = control.getObservation(name);
			obsValues[i++] = obs.getCount();
		}
		return obsValues;
	}
	
	private static class ReversedDoubleComparator implements Comparator<Double> {

		public int compare(Double arg0, Double arg1) {
			return (int) ((arg1 - arg0)* 100);
		}
		
	}



}
