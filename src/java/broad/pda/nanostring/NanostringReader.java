package broad.pda.nanostring;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.math.MathUtil;
import broad.core.math.Statistics;



public class NanostringReader {
	private LinkedHashMap<String, Experiment> experiments = new LinkedHashMap<String, Experiment>();
	private LinkedHashMap<String, Experiment> controls    = new LinkedHashMap<String, Experiment>();
	private static final String[] headerRowNames = {"File Attributes", "File name","ID","Owner","Date","File Version","GeneRLF","Comments",
			"Lane ID","FOV Count","FOV Counted","Scanner ID","StagePosition"};
	private static final double FUDGE_FACTOR  = 0.1; // TODO: Implement the ultimate hack to find the fudge factor that maximizes power: http://www.cbil.upenn.edu/PaGE/doc/PaGE_documentation_technical_manual.pdf
	
	public void load (InputStream is) throws IOException, ParseException {
		BufferedReader br = new BufferedReader(new InputStreamReader(is));
		String line = null;
		boolean atExperiments = false;
		
		List<Experiment> experimentsToLoad = new ArrayList<Experiment>();
		List<List<String>> headerRawInfo = new ArrayList<List<String>>();

		while((line = br.readLine()) != null) {
			if(line.trim().length() == 0) {
				continue;
			}
			//System.err.println(line);
			if(!atExperiments ) {
				if(startsWithHeaderLine(line) ) {
					//System.err.println("Line started with header line ");
					List<String> headerItemInfo = new ArrayList<String>();
					String [] rawInfo = line.split("\t+");
					for(int i = 1; i < rawInfo.length; i++) {
						headerItemInfo.add(rawInfo[i]);
					}
					headerRawInfo.add(headerItemInfo);
				} else if(line.startsWith("Code Class")) {
					atExperiments = true;
					// Assume File name must not be empty.
					List<String> fileNameInfo = headerRawInfo.get(1);
					
					for(int i = 0; i < fileNameInfo.size(); i++) {
						//System.err.println("Header info size "+ headerRawInfo.size() + " File Attributes size (headerRawInfo.get(0).size) " + headerRawInfo.get(0).size() + ", fileNameInfo.size() " + fileNameInfo.size() + " --- " + (headerRawInfo.get(0).size() == fileNameInfo.size() ? headerRawInfo.get(0).get(i) : null) );
						Experiment experiment = new Experiment(headerRawInfo.get(0).size() == fileNameInfo.size() ? headerRawInfo.get(0).get(i) : null,
								fileNameInfo.get(i)
								);
						
						if(experiments.get(experiment.getName()) != null) { 
							int replicateNum = 1;
							String name = experiment.getName();
							while(experiments.get(name) != null) {
								name = experiment.getName() + "_" + replicateNum;
								replicateNum++;
							}
							experiment.setName(name);
						}
						experiments.put(experiment.getName(), experiment);
						experimentsToLoad.add( experiment);
					}
					System.err.println("New total " + experiments.size() ); 

				}
			}else {
				String [] rawInfo = line.split("\t+"); 
				//System.err.println(line);
				if(rawInfo.length - 2 != experimentsToLoad.size()) {
					throw new ParseException("Number of Observation for line " + line + " does not match expected number of experiments " + experimentsToLoad.size());
				}
				String observationClass = rawInfo[0];
				String observationName  = rawInfo[1];
				int experimentNumber = 2; //start at 2 since the first two are class and name;
				Iterator<Experiment> experimentIt = experimentsToLoad.iterator();
				while(experimentIt.hasNext()) {
					Experiment experiment = experimentIt.next();
					experiment.addObservation(observationClass, observationName, Integer.parseInt(rawInfo[experimentNumber]));
					experimentNumber++;
				}
				
			}
		}
	}
	
	public void load(String file) throws IOException, ParseException {
		FileInputStream fis = new FileInputStream(file); 
		try {
			load(fis);
		} finally {
			fis.close();
		}
	}
	
	public Collection<String> getControlExperimentNames() {return controls != null ? controls.keySet() : new ArrayList<String>();}
	
	public void normalizeBySpikeIns() {
		if(experiments.isEmpty()) {
			return;
		}
		Collection<Experiment> allExperiments = experiments.values();
		Iterator<Experiment> experimentIt = allExperiments.iterator();
		
		Experiment first = experimentIt.next();
		first.normalizeByNegativeSpikeIns();
		int referenceTotalSpikeIn = first.getSumPositiveSpikeIns();
		first.floor();
		while(experimentIt.hasNext()) {
			Experiment experiment = experimentIt.next();
			experiment.normalizeByNegativeSpikeIns();
			experiment.normalizeToPositiveControlCount(referenceTotalSpikeIn);
			experiment.floor(); 
		}
		
	}
	
	public void normalizeByNormalizationTranscripts() {
		if(experiments == null || experiments.isEmpty()) {
			return;
		}
		Collection<Experiment> allExperiments = experiments.values();
		Iterator<Experiment> experimentIt = allExperiments.iterator();
		
		Experiment first = experimentIt.next();
		//int referenceTotalSpikeIn = first.getSumNormalizationTranscripts();
		
		while(experimentIt.hasNext()) {
			Experiment experiment = experimentIt.next();
			//experiment.normalizeToNormalizationTranscriptCount(referenceTotalSpikeIn);
			experiment.normalizeToTranscripts(first.getNormalizationTranscriptsObservations());
		}
		
		first.normalizeToTranscripts(first.getNormalizationTranscriptsObservations());
	}
	
	public void computeZScores() {
		if( experiments.isEmpty() || controls.isEmpty()) {
			return;
		}
		Collection<Experiment> allExperiments = experiments.values();
		Iterator<Experiment> experimentIt = allExperiments.iterator();
		Experiment ref = experimentIt.next();
		List<Observation> all = new ArrayList<Observation>();
		all.addAll(ref.getNormalizationTranscriptsObservations());
		all.addAll(ref.observations());
		Iterator<Observation> rowIt = all.iterator();
		while(rowIt.hasNext()) {
			Observation obs = rowIt.next();
			double [] controlExperimentObservations = getLogControlExperimentObservations(obs.getName());
			double controlAvg = Statistics.mean(getControlExperimentObservations(obs.getName()));
			Iterator<Experiment> colIt = experiments.values().iterator();
			while(colIt.hasNext()) {
				Experiment experiment = colIt.next();
				Observation experimentObs = experiment.getObservation(obs.getName());
				double zScore = Statistics.zScore(MathUtil.log2(experimentObs.getCount()), controlExperimentObservations, FUDGE_FACTOR);
				experimentObs.setFoldChange(experimentObs.getCount()/controlAvg);
				//double confidence = Statistics.zScoreConfidence(zScore, controlObservations);
				experimentObs.setZScore(zScore);
				//experimentObs.setConfidence(confidence);
			}
		}
		
	}
	
	public MatrixWithHeaders getCountMatrix() {
		MatrixWithHeaders ret = null;
		if(experiments != null && !experiments.isEmpty()) {
			List<String> cols = new ArrayList<String>(experiments.keySet());
			List<String> rows = new ArrayList<String>();
			Experiment first = experiments.values().iterator().next();
			Collection<Observation> observations = first.observations();
			
			for(Observation obs : observations) {
				rows.add(obs.getName());
			}
			ret = new MatrixWithHeaders(rows, cols);
			
			for(Experiment e : experiments.values()) {
				for(String obs : rows) {
					ret.set(obs, e.getName(), e.getObservation(obs).getCount());
				}
			}
		}
		
		return ret;
	}
	
	/**
	 * randomizationType is 0 for row-based, 1 for column (Experiment) based and 2 is for using both.
	 * @param randomizationType
	 * @param minFoldChange 
	 */
	public void computeConfidenceScores(int randomizationType, float minFoldChange) {
		if( experiments.isEmpty() || controls.isEmpty()) {
			return;
		}
		
		Experiment first = experiments.values().iterator().next();		
		Collection<Observation> firstExperimentObservationList = first.experimentAndControlObservations();
		//Initialize datastructures to store observations and permutations
		HashMap<String, List<List<Double>>> negativePermutationsByExperiment = new HashMap<String, List<List<Double>>>(experiments.size());
		HashMap<String, List<List<Double>>> positivePermutationsByExperiment = new HashMap<String, List<List<Double>>>(experiments.size());
		HashMap<String, List<List<Double>>> negativePermutationsByTranscript = new HashMap<String, List<List<Double>>>(firstExperimentObservationList.size());
		HashMap<String, List<List<Double>>> positivePermutationsByTranscript = new HashMap<String, List<List<Double>>>(firstExperimentObservationList.size());
		HashMap<String, List<Double>> negativeObservationsByExperiment = new HashMap<String, List<Double>>();
		HashMap<String, List<Double>> positiveObservationsByExperiment = new HashMap<String, List<Double>>();
		HashMap<String, List<Double>> negativeObservationsByTranscript = new HashMap<String, List<Double>>();
		HashMap<String, List<Double>> positiveObservationsByTranscript = new HashMap<String, List<Double>>();

		Map<String, double[]> controlLogObservations = new HashMap<String, double[]>();
		Iterator<Observation> obsIt = firstExperimentObservationList.iterator();

		//Set up data structures.
		while(obsIt.hasNext()) {
			Observation obs = obsIt.next();
			List<Double> positiveTranscriptObservations = new ArrayList<Double>();
			positiveObservationsByTranscript.put(obs.getName(), positiveTranscriptObservations);
			List<Double> negativeTranscriptObservations = new ArrayList<Double>();
			negativeObservationsByTranscript.put(obs.getName(), negativeTranscriptObservations);
			
			List<List<Double>> positiveObservationPermutedData = new ArrayList<List<Double>>(controls.size());
			List<List<Double>> negativeObservationPermutedData = new ArrayList<List<Double>>(controls.size());
			for(int i = 0; i < controls.size(); i++) {
				List<Double> pList = new ArrayList<Double>();
				positiveObservationPermutedData.add(pList);
				List<Double> nList = new ArrayList<Double>();
				negativeObservationPermutedData.add(nList);
			}
			positivePermutationsByTranscript.put(obs.getName(), positiveObservationPermutedData);
			negativePermutationsByTranscript.put(obs.getName(), negativeObservationPermutedData);
			controlLogObservations.put(obs.getName(),  getLogControlExperimentObservations(obs.getName()));
		}
		
		Iterator<Experiment> experimentIt = experiments.values().iterator();

		//Compute permutations & fill up observed lists
		while(experimentIt.hasNext()) {
			Experiment experiment = experimentIt.next();
			
			List<Double> positiveExperimentObservations = new ArrayList<Double>();
			positiveObservationsByExperiment.put(experiment.getName(), positiveExperimentObservations);
			List<Double> negativeExperimentObservations = new ArrayList<Double>();
			negativeObservationsByExperiment.put(experiment.getName(), negativeExperimentObservations);

			List<List<Double>> negativeExperimentPermutedData = new ArrayList<List<Double>>(controls.size());
			List<List<Double>> positiveExperimentPermutedData = new ArrayList<List<Double>>(controls.size());
			for(int i = 0; i < controls.size(); i++) {
				List<Double> pList = new ArrayList<Double>();
				positiveExperimentPermutedData.add(pList);
				List<Double> nList = new ArrayList<Double>();
				negativeExperimentPermutedData.add(nList);
			}
			positivePermutationsByExperiment.put(experiment.getName(), positiveExperimentPermutedData);
			negativePermutationsByExperiment.put(experiment.getName(), negativeExperimentPermutedData);
			
			obsIt = experiment.experimentAndControlObservations().iterator();
			while(obsIt.hasNext()) {
				Observation obs = obsIt.next();
				//First observations
				double score = obs.getZScore();
				
				if(score > 0) {
					positiveExperimentObservations.add(score);
					positiveObservationsByTranscript.get(obs.getName()).add(score);
				} else {
					negativeExperimentObservations.add(-score);
					negativeObservationsByTranscript.get(obs.getName()).add(-score);
				}
				
				
				//Now permutations
				double [] permutedZScoresArray = Statistics.computePermutedZScores(MathUtil.log2(obs.getCount()), controlLogObservations.get(obs.getName()), FUDGE_FACTOR);
				List<List<Double>> positiveTranscriptPermutedData = positivePermutationsByTranscript.get(obs.getName());
				List<List<Double>> negativeTranscriptPermutedData = negativePermutationsByTranscript.get(obs.getName());
				for(int i = 0; i < controls.size(); i++) {
					if(permutedZScoresArray[i] > 0) {
						positiveExperimentPermutedData.get(i).add(permutedZScoresArray[i]);
						positiveTranscriptPermutedData.get(i).add(permutedZScoresArray[i]);
					} else {
						negativeExperimentPermutedData.get(i).add(-permutedZScoresArray[i]);
						negativeTranscriptPermutedData.get(i).add(-permutedZScoresArray[i]);
					}
				}
			}
		}
		//Sort permutation lists
		ReversedDoubleComparator comparator = new ReversedDoubleComparator();
		Iterator<List<List<Double>>> positivePermutationsByObservationIt = positivePermutationsByTranscript.values().iterator();
		sortLists(positivePermutationsByObservationIt, comparator);
		Iterator<List<List<Double>>> negativePermutationsByObservationIt = negativePermutationsByTranscript.values().iterator();
		sortLists(negativePermutationsByObservationIt, comparator);

		Iterator<List<List<Double>>> positivePermutationsByExperimentIt = positivePermutationsByExperiment.values().iterator();			
		sortLists(positivePermutationsByExperimentIt, comparator);
		Iterator<List<List<Double>>> negativePermutationsByExperimentIt = negativePermutationsByExperiment.values().iterator();			
		sortLists(negativePermutationsByExperimentIt, comparator);
		
		Iterator<List<Double>> positiveTranscriptObservationListIt = positiveObservationsByTranscript.values().iterator();
		while(positiveTranscriptObservationListIt.hasNext()) {
			Collections.sort(positiveTranscriptObservationListIt.next(), comparator);
		}
		Iterator<List<Double>> negativeTranscriptObservationListIt = negativeObservationsByTranscript.values().iterator();
		while(negativeTranscriptObservationListIt.hasNext()) {
			Collections.sort(negativeTranscriptObservationListIt.next(), comparator);
		}
		
		Iterator<List<Double>> positiveExperimentObservationListIt = positiveObservationsByExperiment.values().iterator();
		while(positiveExperimentObservationListIt.hasNext()) {
			Collections.sort(positiveExperimentObservationListIt.next(), comparator);
		}
		Iterator<List<Double>> negativeExperimentObservationListIt = negativeObservationsByExperiment.values().iterator();
		while(negativeExperimentObservationListIt.hasNext()) {
			Collections.sort(negativeExperimentObservationListIt.next(), comparator);
		}
		
		//Now compute confidence scores
		double minLogFoldChange = MathUtil.log2(minFoldChange);
		experimentIt = experiments.values().iterator();
		while(experimentIt.hasNext()) {
			Experiment e = experimentIt.next();
			obsIt = e.experimentAndControlObservations().iterator();
			while(obsIt.hasNext()) {
				Observation obs = obsIt.next();
				double zScore = obs.getZScore();
				if(Math.abs(MathUtil.log2(obs.getFoldChange())) < minLogFoldChange) {
					obs.setConfidence(0);
					continue;
				}
				List<Double> observationListToUse = null;
				switch(randomizationType){
				case 0: 
					if(zScore > 0) {
						observationListToUse = positiveObservationsByTranscript.get(obs.getName());
					} else {
						observationListToUse = negativeObservationsByTranscript.get(obs.getName());
					}
					break;
				case 1:
					if(zScore > 0) {
						observationListToUse = positiveObservationsByExperiment.get(e.getName());
					} else {
						observationListToUse = negativeObservationsByExperiment.get(e.getName());
					}
					break;
				default:
					observationListToUse = new ArrayList<Double>();//used to initialize to the right size, but got lazy when split into positive/negative
					if(zScore > 0) {
						observationListToUse.addAll(positiveObservationsByTranscript.get(obs.getName()));
						observationListToUse.addAll(positiveObservationsByExperiment.get(e.getName()));
					} else {
						observationListToUse.addAll(negativeObservationsByTranscript.get(obs.getName()));
						observationListToUse.addAll(negativeObservationsByExperiment.get(e.getName()));						
					}
					Collections.sort(observationListToUse, comparator);
				}
				
				List<List<Double>> permuttedZScores = null;
				switch(randomizationType){
				case 0: 
					if(zScore > 0) {
						permuttedZScores = positivePermutationsByTranscript.get(obs.getName()) ;
					} else {
						permuttedZScores = negativePermutationsByTranscript.get(obs.getName()) ;
					}
					break;
				case 1:
					if(zScore > 0) {
						permuttedZScores = positivePermutationsByExperiment.get(e.getName());
					} else {
						permuttedZScores = negativePermutationsByExperiment.get(e.getName());
					}
					
					break;
				default:
					permuttedZScores = new ArrayList<List<Double>>();
					if(zScore > 0) {	
						permuttedZScores.addAll(positivePermutationsByTranscript.get(obs.getName()));
						permuttedZScores.addAll(positivePermutationsByExperiment.get(e.getName()));
					} else {
						permuttedZScores.addAll(negativePermutationsByTranscript.get(obs.getName()));
						permuttedZScores.addAll(negativePermutationsByExperiment.get(e.getName()));					
					}
				}
				
				float meanPermuttedObservationsAboveZScore = 0;
				
				double fdr = 1;
				int i = observationListToUse.size() - 1;
				double obsZ = 0;
				//System.out.println("Beging FDR comp for obs " + obs.getName() + " value : " + zScore + " #obs: " + observationListToUse.size());
				//System.out.println("\tFull value list: " + observationListToUse);
				while(Math.abs(obsZ)<= Math.abs(zScore) && i >= 0) {
					obsZ=observationListToUse.get(i--);
					//System.out.print("\tObs: " + obsZ + "- current FDR: " + fdr +" i:"+(i+1));
					double numOfObservationsAboveZScore =	countObservationsLargerThan(Math.abs(obsZ), observationListToUse);
					meanPermuttedObservationsAboveZScore = 0;
					for(int p = 0; p < permuttedZScores.size(); p++) {
						// NOTE: I do not think the following was the right way to do this. The mean is the mean # observations regardless of the size of the negative or positive observations.
						//Since permuted lists may not have the same number of observations than the observed list when we break things into + and -
						// we normalize to the total observations of the same sign made.
						if(permuttedZScores.get(p).size() > 0) {
							meanPermuttedObservationsAboveZScore += 
								 countObservationsLargerThan(Math.abs(obsZ), permuttedZScores.get(p));
								//meanPermuttedObservationsAboveZScore + observationListToUse.size() * countObservationsLargerThan(Math.abs(zScore), permuttedZScores.get(p))/ (float)permuttedZScores.get(p).size();
						}
					}
					double valueFDR = meanPermuttedObservationsAboveZScore/((double) (numOfObservationsAboveZScore * permuttedZScores.size()));
					fdr = Math.min(fdr, valueFDR); 
					//System.out.print(" permuttedZScores size: " + permuttedZScores.size() + " total permutted bos > obs: " + meanPermuttedObservationsAboveZScore +" obs > obs: " + numOfObservationsAboveZScore );
					//System.out.println(" value FDR: " + valueFDR + " end FDR: " + fdr);
				}
				
				obs.setConfidence(fdr >= 1 ? 0 : 1 - fdr);
			}
			//experiment.computeConfidenceValues(controls.values(), byRow);
		}		
	}
	
	public void computeConfidenceScoresUsingTranscripts(String[] transcripts, float minFoldChange) {
		List<List<Double>> controlObservations = new ArrayList<List<Double>>(transcripts.length);
		for(String transcript : transcripts) {
			ArrayList<Double> transcriptObs = new ArrayList<Double>(experiments.size()); 
			controlObservations.add(transcriptObs);
			for (Experiment e : experiments.values()) {
				transcriptObs.add(Math.abs(e.getObservation(transcript).getZScore()));
			}
			Collections.sort(transcriptObs, new Comparator<Double>() {

				public int compare(Double arg0, Double arg1) {
					return (int) (100 * (arg1 - arg0));
				}
				
			});
		}
		
		HashMap<String, List<Double>> transcriptObservations = new HashMap<String, List<Double>>();
		Collection<Observation> obs = experiments.values().iterator().next().observations();
		for(Observation o : obs) {
			ArrayList<Double> transcriptObs = new ArrayList<Double>(experiments.size());
			transcriptObservations.put(o.getName(), transcriptObs);
			for (Experiment e : experiments.values()) {
				transcriptObs.add(Math.abs(e.getObservation(o.getName()).getZScore()));
			}
			Collections.sort(transcriptObs, new Comparator<Double>() {

				public int compare(Double arg0, Double arg1) {
					return (int) (100 * (arg1 - arg0));
				}
				
			});
		}
		double minLogFoldChange = MathUtil.log2(minFoldChange);
		for(Experiment e : experiments.values()) {
			Collection<Observation> experimentObs = e.observations();
			for(Observation o : experimentObs) {
				if(Math.abs(MathUtil.log2(o.getFoldChange())) < minLogFoldChange) {
					o.setConfidence(0);
				} else {
					double oZScore = o.getZScore();
					double absOZScore = Math.abs(oZScore);
					int numEqualOrLargerThanObserved = 0;
					for(List<Double> controlTranscriptObs : controlObservations) {
						numEqualOrLargerThanObserved += countObservationsLargerThan(absOZScore, controlTranscriptObs);
					}
					float expectedControlObsAboveScore = numEqualOrLargerThanObserved /(float)transcripts.length;
					int numObservedAboveScore = countObservationsLargerThan(absOZScore, transcriptObservations.get(o.getName()));
					double confidence =  expectedControlObsAboveScore > numObservedAboveScore ? 0 : 1 - expectedControlObsAboveScore/numObservedAboveScore;
					o.setConfidence(confidence);
				}
			}
		}
	}
	
	private double[] getLogControlExperimentObservations(String name) {
		double [] obsValues = new double [controls.size()];
		Iterator<Experiment> it = controls.values().iterator();
		int i = 0;
		while(it.hasNext()) {
			Experiment control = it.next();
			Observation obs = control.getObservation(name);
			obsValues[i++] = MathUtil.log2(obs.getCount());
		}
		return obsValues;
	}
	
	private double[] getControlExperimentObservations(String name) {
		double [] obsValues = new double [controls.size()];
		Iterator<Experiment> it = controls.values().iterator();
		int i = 0;
		while(it.hasNext()) {
			Experiment control = it.next();
			Observation obs = control.getObservation(name);
			obsValues[i++] = obs.getCount();
		}
		return obsValues;
	}

	private boolean startsWithHeaderLine(String line) {
		boolean startsWithHeader = false;
		for(int i = 0; i < headerRowNames.length; i++) {
			if(line.startsWith(headerRowNames[i])) {return true;}
		}
		return startsWithHeader;
	}
	
	Collection<Experiment> getExperiments() { return experiments.values();}

	public void setControlTranscripts(List<String> controlTranscripts) {
		if(experiments != null && !experiments.isEmpty()) {
			Collection<Experiment> allExperiments = experiments.values();
			Iterator<Experiment> experimentIt = allExperiments.iterator();
			while(experimentIt.hasNext()) {
				Experiment experiment = experimentIt.next();
				experiment.setNormalizationTranscripts(controlTranscripts);
			}
		}
		
	}

	public void setFirstExperiment(int index) {
		if(index < 0 || index >= experiments.size()) {
			throw new IllegalStateException("Tried to set "+ index +" as first experiment but there are only " + experiments.size());
		}
		
		int i = 0;
		Iterator<Experiment> it = experiments.values().iterator();
		while(it.hasNext() && i++ < index){
			it.next();
		}
		Experiment chosen = it.next();	
		experiments.remove(chosen.getName());
		
		LinkedHashMap<String, Experiment> reorderedExperiments = new LinkedHashMap<String, Experiment>(experiments.size());
		
		reorderedExperiments.put(chosen.getName(), chosen);
		
		it = experiments.values().iterator();
		while(it.hasNext()){
			Experiment experiment = it.next();
			reorderedExperiments.put(experiment.getName(), experiment);
		}
		
		experiments = reorderedExperiments;
	}

	/**
	 * 
	 * @param os
	 * @param scoreToWrite 0: Transcript count 1: confidence, other value: z-score
	 * @param minConfidence
	 * @param useGenePatternFormat
	 * @throws IOException
	 */
	public void writeSpikedInNormalized(OutputStream os, int scoreToWrite, boolean useGenePatternFormat, double downConf, double upConf) throws IOException {
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(os));
		if(useGenePatternFormat) {
			bw.write("#1.2\n");
			Experiment first = experiments.values().iterator().next();
			bw.write(String.valueOf(first.observations().size() + first.getNormalizationTranscriptsObservations().size()));
			bw.write("\t");
			bw.write(String.valueOf(experiments.size()));
			bw.newLine();
			bw.write("name\tID");
		}
		if(experiments != null && !experiments.isEmpty()) {
			Collection<Experiment> allExperiments = experiments.values();
			List<Experiment> orderedExperiments = new ArrayList<Experiment>(allExperiments.size());
			Iterator<Experiment> controlExpIt = controls.values().iterator();
			while(controlExpIt.hasNext()) {
				Experiment e = controlExpIt.next();
				bw.write("\t");
				bw.write(e.getName());
				orderedExperiments.add(e);
			}
			Iterator<Experiment> experimentIt = allExperiments.iterator();
			while(experimentIt.hasNext()) {
				Experiment e = experimentIt.next();
				if(!controls.containsKey(e.getName())) {
					bw.write("\t");
					bw.write(e.getName()); 
					orderedExperiments.add(e);
				}
			}
			bw.newLine();
			assert(orderedExperiments.size() == allExperiments.size());
			
			experimentIt = orderedExperiments.iterator();
			Experiment reference = experimentIt.next();
			Iterator<Observation> observationIt = reference.getNormalizationTranscriptsObservations().iterator();
			writeObservations(bw, orderedExperiments, observationIt, scoreToWrite,  useGenePatternFormat, downConf, upConf);
			
			observationIt = reference.observations().iterator();
			writeObservations(bw, orderedExperiments, observationIt, scoreToWrite,  useGenePatternFormat, downConf, upConf);
			
		}
		bw.flush();
	}
	
	public void writeCytoscapeNetwork(String prefix, List<String> ignoreExperiments, List<String> ignoreTranscripts, double downConf, double upConf) throws IOException {
		LinkedHashMap<String, Integer> nodeTypeMap = new LinkedHashMap<String, Integer>();
		Experiment first = experiments.values().iterator().next();
		Iterator<Observation> obsIt = first.observations().iterator();
		while(obsIt.hasNext()) {
			Observation obs = obsIt.next();
			if(ignoreTranscripts != null && !ignoreTranscripts.contains(obs.getName())) {
				nodeTypeMap.put(obs.getName(), 0);
			}
		}

		
		System.out.println("Writing interactions to " + prefix + ".sif file" );
		BufferedWriter netWR = new BufferedWriter(new FileWriter(prefix+".sif"));
		System.out.println("Writing edge z-scores to " + prefix + "_edge_zscores.edo file" );
		BufferedWriter edZSWR = new BufferedWriter(new FileWriter(prefix + "_edge_zscores.eda"));
		System.out.println("Writing edge z-scores to " + prefix + "_edge_type.edo file" );
		BufferedWriter edTyWR = new BufferedWriter(new FileWriter(prefix + "_edge_type.eda"));
		edZSWR.write("EdgeZScore\n");
		edTyWR.write("EdgeType\n");
		Iterator<Experiment>  experimentIt = experiments.values().iterator();
		while(experimentIt.hasNext()) {
			Experiment e = experimentIt.next();
			if(controls.containsKey(e.getName()) ||  (ignoreExperiments != null && ignoreExperiments.contains(e.getName()))) {
				continue;
			}
			nodeTypeMap.put(e.getName(), nodeTypeMap.containsKey(e.getName()) ? 2 : 1);

			obsIt = e.observations().iterator();
			while(obsIt.hasNext()) {
				Observation obs = obsIt.next();
				if(nodeTypeMap.containsKey(obs.getName())  && ( obs.getZScore() > 0 && obs.getConfidence() > upConf || obs.getZScore() < 0 && obs.getConfidence() > upConf)) {
					netWR.append(e.getName()).append(" pp ").append(obs.getName());
					netWR.newLine();
					
					edZSWR.append(e.getName()).append(" (pp) ").append(obs.getName()).append(" = ").append(String.valueOf(obs.getZScore()));
					edZSWR.newLine();
					edTyWR.append(e.getName()).append(" (pp) ").append(obs.getName()).append(" = ").append(obs.getZScore() < 0 ? "Activator" : "Represor");
					edTyWR.newLine();
				}
			}
		}
		edTyWR.close();
		edZSWR.close();
		netWR.close();
		
		System.out.println("Writing node types (Affected, Effector or both) to " + prefix + "_node_type.noa");
		BufferedWriter noTyWR = new BufferedWriter(new FileWriter(prefix + "_node_type.noa"));
		noTyWR.write("NodeType\n");
		Iterator<String> nodeIt = nodeTypeMap.keySet().iterator();
		while(nodeIt.hasNext()) {
			String node = nodeIt.next();
			int type    = nodeTypeMap.get(node);
			noTyWR.write(node);
			noTyWR.write(" = ");
			switch (type) {
			case 0:
				noTyWR.write("Affected");
				break;
			case 1:
				noTyWR.write("Effector");
				break;
			case 2:
				noTyWR.write("Both");
				break;
			}
			noTyWR.newLine();
		}
		noTyWR.close();
	}
	
	public void setControls(int [] cols) {
		if(experiments != null && !experiments.isEmpty()) {
			controls = new LinkedHashMap<String, Experiment>();
			List<Experiment> experimentList = new ArrayList<Experiment>(experiments.values());
			for(int i = 0; i < cols.length; i++) {
				Experiment e = experimentList.get(cols[i]);
				System.err.println("\tExperiment " + e.getName() + " col("+cols[i]+") is set as control ");
				controls.put(e.getName(), e);
			}			
		}	
	}
	
	public void addControls(int [] cols) {
		if(experiments != null && !experiments.isEmpty()) {
			if(controls == null) {
				controls = new LinkedHashMap<String, Experiment>();
			}
			List<Experiment> experimentList = new ArrayList<Experiment>(experiments.values());
			for(int i = 0; i < cols.length; i++) {
				Experiment e = experimentList.get(cols[i]);
				System.err.println("\tExperiment " + e.getName() + " col("+cols[i]+") is set as control ");
				controls.put(e.getName(), e);
			}			
		}		
	}
	

	private void writeObservations(BufferedWriter bw,
			Collection<Experiment> allExperiments,
			Iterator<Observation> observationIt,
			int scoreToWrite,
			boolean useGenePatternFormat,
			double downConf, double upConf) throws IOException {
		Iterator<Experiment> experimentIt;
		while(observationIt.hasNext()) {
			Observation refObs = observationIt.next();
			bw.write(refObs.getName());
			bw.write("\t");			
			if(useGenePatternFormat) {
				bw.write(refObs.getName());
				bw.write("\t");
			}
			
			writeScore(bw, scoreToWrite,  refObs, downConf, upConf);

			experimentIt = allExperiments.iterator();
			experimentIt.next(); //ignore first which we already used
			while(experimentIt.hasNext()) {
				Experiment experiment = experimentIt.next();
				Observation expObs = experiment.getObservation(refObs.getName());
				bw.write("\t");
				writeScore(bw, scoreToWrite, expObs, downConf, upConf);
			}
			bw.newLine();
		}
	}

	private void writeScore(BufferedWriter bw, int scoreToWrite,
			Observation refObs, double downConf, double upConf) throws IOException {
		switch (scoreToWrite) {
		case 0: 
			bw.write(String.valueOf(refObs.getCount()));
			break;
		case 1:
			bw.write(String.valueOf(refObs.getConfidence()));
			break;
		case 3:
			bw.write(String.valueOf((MathUtil.log2(refObs.getFoldChange()))));
			break;
		default:
			double score =  ( refObs.getZScore() > 0 && refObs.getConfidence() >= upConf || refObs.getZScore() < 0 && refObs.getConfidence() >= downConf)  ? refObs.getZScore() : 0;
			bw.write(String.valueOf(score));	
			break;
		}
	}

	private static class ReversedDoubleComparator implements Comparator<Double> {

		public int compare(Double arg0, Double arg1) {
			return (int) ((arg1 - arg0)* 10000);
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
	
	private void sortLists(Iterator<List<List<Double>>> permutationsByObservationIt, Comparator<Double> comparator) {

		while(permutationsByObservationIt.hasNext()) {
			List<List<Double>> permutations = permutationsByObservationIt.next();
			for(int i = 0; i < permutations.size(); i++) {
				Collections.sort(permutations.get(i), comparator);
			}
		}
	}



	
}
