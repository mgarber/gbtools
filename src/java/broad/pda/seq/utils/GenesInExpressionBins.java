package broad.pda.seq.utils;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;

public class GenesInExpressionBins<T extends GenesInExpressionBinsValue> {
	
	private double minExpression;
	private double maxExpression;
	private double minExpression2;
	private double maxExpression2;
	private double[][] expressionIntervals;
	private Collection<RefSeqGene>[] genesInInterval;
	private ArrayList<T> intervalValues;
	private Class<T> valueClass;
	private ContinuousDataAlignmentModel data;
	private String chr;
	
	public GenesInExpressionBins(Collection<RefSeqGene> genes, int numOfBins, Class<T> valueClass, ContinuousDataAlignmentModel data, String chr, double middlePercent) {		
		Collection<RefSeqGene> genesWithMiddleExpression = extractMiddleOfExpressionDataSet(genes, middlePercent);
		
		double binSize = (this.maxExpression2 - this.minExpression2) / (double)numOfBins;
		
		expressionIntervals = new double[numOfBins][2];
		genesInInterval = new Collection[numOfBins];
		intervalValues = new ArrayList<T>(numOfBins);
				
		for (int i = 0; i < numOfBins; i++) {
			expressionIntervals[i][0] = this.minExpression2 + i * binSize;
			expressionIntervals[i][1] = this.minExpression2 + (i + 1) * binSize;

			genesInInterval[i] = new TreeSet<RefSeqGene>();
			intervalValues.add(i, null);
		}
		
		for (RefSeqGene gene : genesWithMiddleExpression) {
			double expression = this.getExpressionOfGene(gene);
			
			for (int i = 0; i < numOfBins; i++) {
				if ((expressionIntervals[i][0] <= expression && expressionIntervals[i][1] > expression) || (expressionIntervals[i][1] == this.maxExpression2 && expression == this.maxExpression2)) {
					genesInInterval[i].add(gene);
					break;
				}
			}
		}
		
		this.valueClass = valueClass;
		this.data = data;
		this.chr = chr;
	}
	
	public T getValueForExpressionGivenBin(int bin) throws IOException{
		if (bin < 0) {
			System.out.println("bin < 0");
			return null;
		}
		
		if (intervalValues.get(bin) == null) {
			Collection<RefSeqGene> genes = genesInInterval[bin];
			
			try {
				T t = valueClass.newInstance();

				t.computeValue(genes, this.data, this.chr);
				
				intervalValues.set(bin, t);
				
				return t;
			}
			catch (IllegalAccessException e) {
				e.printStackTrace();
				return null;
			}
			catch (InstantiationException e) {
				e.printStackTrace();
				return null;
			}
		} else {
			return intervalValues.get(bin);
		}
	}
	
	public T getValueForExpression(double expression) throws IOException{
		int i = -1;
		if (expression <= this.minExpression2) {
			i = 0;
		} else if (expression >= this.maxExpression2) {
			i = expressionIntervals.length - 1;
		} else {
			for (int j = 0; j < expressionIntervals.length; j++) {
				if (expression >= expressionIntervals[j][0] && expression < expressionIntervals[j][1]) {
					i = j;
					break;
				}
			}
		}

		return getValueForExpressionGivenBin(i);
	}
	
	public T getValueForExpressionUsingNearbyExpressionBinsIfNeeded(double expression) throws IOException{
		T initialValue = getValueForExpression(expression);
		if (initialValue.couldComputeValue())
			return initialValue;
		
		Map<Double, Set<Integer>> distancesOfMeanExpressionOfBinFromExpressionWithBinNum = new HashMap<Double, Set<Integer>>();
		List<Double> distancesOfMeanExpressionOfBinFromExpression = new ArrayList<Double>(intervalValues.size());
		for (int j = 0; j < expressionIntervals.length; j++) {
			if (expression >= expressionIntervals[j][0] && expression < expressionIntervals[j][1])
				continue;
			
			double meanExpressionOfBin = (expressionIntervals[j][0] + expressionIntervals[j][1]) / 2.0;
			double distanceOfMeanExpressionOfBinFromExpression = Math.abs(expression - meanExpressionOfBin);
			
			if (!distancesOfMeanExpressionOfBinFromExpressionWithBinNum.containsKey(distanceOfMeanExpressionOfBinFromExpression)) {
				Set<Integer> s = new HashSet<Integer>();
				distancesOfMeanExpressionOfBinFromExpressionWithBinNum.put(distanceOfMeanExpressionOfBinFromExpression, s);
			}
			distancesOfMeanExpressionOfBinFromExpressionWithBinNum.get(distanceOfMeanExpressionOfBinFromExpression).add(j);
			
			distancesOfMeanExpressionOfBinFromExpression.add(distanceOfMeanExpressionOfBinFromExpression);
		}
		
		Collections.sort(distancesOfMeanExpressionOfBinFromExpression);
		
		for (double distanceBetweenExpressions : distancesOfMeanExpressionOfBinFromExpression) {
			Set<Integer> binNums = distancesOfMeanExpressionOfBinFromExpressionWithBinNum.get(distanceBetweenExpressions);
			
			for (int binNum : binNums) {
				T t = getValueForExpressionGivenBin(binNum);
				if (t.couldComputeValue())
					return t;
			}
		}
		
		return null;
	}
	
	private Collection<RefSeqGene> extractMiddleOfExpressionDataSet(Collection<RefSeqGene> genes, double middlePercent) {
		TreeMap<Double, Collection<RefSeqGene>> expressionsWithGenes = new TreeMap<Double, Collection<RefSeqGene>>();
		List<Double> expressions = new ArrayList<Double>(genes.size());
		for (RefSeqGene gene : genes) {
			Double expression = this.getExpressionOfGene(gene);
			expressions.add(expression);
			
			if (!expressionsWithGenes.containsKey(expression)) {
				Collection<RefSeqGene> ts = new TreeSet<RefSeqGene>();
				expressionsWithGenes.put(expression, ts);
			}
			expressionsWithGenes.get(expression).add(gene);
		}
		Collections.sort(expressions);

		int numToRemoveFromEachEnd = (int)(((1 - middlePercent) / 2.0) * genes.size());
		
		int startOfMiddle = numToRemoveFromEachEnd;
		int endOfMiddle = expressions.size() - numToRemoveFromEachEnd - 1;
		
		Collection<RefSeqGene> middleGenes = new TreeSet<RefSeqGene>();
		for (int i = startOfMiddle; i <= endOfMiddle; i++) {
			Double expression = expressions.get(i);
			Collection<RefSeqGene> genesWithExpression = expressionsWithGenes.get(expression);
			for (RefSeqGene gene : genesWithExpression) {
				middleGenes.add(gene);
			}
		}
		
		this.minExpression = expressions.get(0);
		this.maxExpression = expressions.get(expressions.size() - 1);
		
		this.minExpression2 = expressions.get(startOfMiddle);
		this.maxExpression2 = expressions.get(endOfMiddle);
		
		return middleGenes;
	}
	
	/**
	private void findMinAndMaxExpression(Collection<RefSeqGene> genes) {
		this.minExpression = Double.MAX_VALUE;
		this.maxExpression = -Double.MAX_VALUE;
		
		for (RefSeqGene gene : genes) {
			double expression = this.getExpressionOfGene(gene);
			
			this.minExpression = Math.min(this.minExpression, expression);
			this.maxExpression = Math.max(this.maxExpression, expression);
		}
	}
	*/
	
	private double getExpressionOfGene(RefSeqGene gene) {
		return gene.getRPKM();
	}
	
}
