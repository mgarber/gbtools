package broad.pda.arrays.tilling;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.MatrixWithHeaders;
import broad.core.math.EmpiricalDistribution;
import broad.core.math.Statistics;
import broad.pda.datastructures.Alignments;
import broad.pda.datastructures.ContinuousData;
import broad.pda.datastructures.LocationAwareMatrix;

public class WindowPermuteDifferentialExpression {
	private static final int EMPIRICAL_DIST_BINS = 500;
	
	private ContinuousData [] group1Data;
	private ContinuousData [] group2Data;
	
	private Map<String, TreeMap<Alignments, List<Double>>[]> group1Scores;
	private Map<String, TreeMap<Alignments, List<Double>>[]> group2Scores;
	
	private String [] group1;
	private String [] group2;
	
	public WindowPermuteDifferentialExpression(LocationAwareMatrix matrix, String [] group1, String [] group2, String chr, int maxDistance, Collection<? extends LightweightGenomicAnnotation> regions)throws IOException{
		this.group1Scores = new LinkedHashMap<String, TreeMap<Alignments, List<Double>>[]>();
		this.group2Scores = new LinkedHashMap<String, TreeMap<Alignments, List<Double>>[]>();
		group1Data = makeData(matrix, group1);
		group2Data = makeData(matrix, group2);
		this.group1 = group1;
		this.group2 = group2;
		
		//System.err.println(peaks);
		Set<String> chromosomes = new TreeSet<String>();
		if(chr!= null && chr.length()>0) { //Replaced to be compatible with 1.5
			chromosomes.add(chr.replace("chr", ""));
		} else {
			chromosomes = matrix.getChromosomes();
		}
		
		for(String chromosome : chromosomes) {
			Set<Alignments> chrAlignments = regions != null ? 
					group1Data[0].makeContiguousRegions(chromosome, regions) : 
						group1Data[0].makeContiguousRegions(chromosome, maxDistance);
					//System.err.println("Made Contiguous Regions");
					//System.err.println(alignments);
		
			TreeMap<Alignments, List<Double>>[] scores1=getValuesPerTU(group1Data, chrAlignments);
			group1Scores.put(chromosome, scores1);
			
			TreeMap<Alignments, List<Double>>[] scores2=getValuesPerTU(group2Data, chrAlignments);
			group2Scores.put(chromosome, scores2);
		}
	}
	
	public MatrixWithHeaders ttest(double alpha, int numPerm, boolean filterByFDR) {
		MatrixWithHeaders matrix = null;
		
		//Prepare matrix
		List<String> rowNames = new ArrayList<String>();
		List<String> colNames = new ArrayList<String>(4+ group1.length + group2.length);
		colNames.add("PValue");
		colNames.add("PValue(FDR)");
		colNames.add("Group1 MoMs");
		colNames.add("Group2 MoMs");
		for(String col : group1) {
			colNames.add(col);
		}
		for(String col : group2) {
			colNames.add(col);
		}
		
		for (String chr : group1Scores.keySet()) {
			TreeMap<Alignments, List<Double>> chrGroup1Scores = group1Scores.get(chr)[0];
			Set<Alignments> alignments = chrGroup1Scores.keySet();
			for(Alignments align : alignments) {
				rowNames.add(align.toUCSC());
			}
		}
		
		matrix = new MatrixWithHeaders(rowNames, colNames);
		
		List<Double> pvalList = new ArrayList<Double>();
		//Evaluate differential expression
		for (String chr : group1Scores.keySet()) {
			TreeMap<Alignments, List<Double>> chrGroup1Scores = group1Scores.get(chr)[0];
			Set<Alignments> alignments = chrGroup1Scores.keySet();
			//Get null and extreme value distributions
			long startTime = System.currentTimeMillis();
			EmpiricalDistribution[] randomDist = computeRandomScoresMinMax(alignments, group1Data, group2Data, chr, numPerm);
			System.err.println("Finished permutations for chromosome " + chr + ", took " + ((System.currentTimeMillis() - startTime)/1000) + " sec.");
			for(Alignments align : alignments) {
				String alignRowName = align.toUCSC();
				double tstatistic = ttest(align);
				double pMin=randomDist[0].getCummulativeProbability(tstatistic);
				double pMax=1-randomDist[1].getCummulativeProbability(tstatistic);
				double fwerPval = Math.min(pMin, pMax);
				double pVal = randomDist[2].getCummulativeProbability(tstatistic);
				pVal = Math.min(1-pVal, pVal);
				matrix.set(alignRowName, 0, fwerPval);
				matrix.set(alignRowName, 1, pVal);
				double group1MoM=0;
				double group2MoM=0;
				for( int i = 0; i< group1.length; i++) {
					double group1MemberMedian = Statistics.median(group1Scores.get(align.getChr())[i].get(align));
					matrix.set(alignRowName, group1[i], group1MemberMedian);
					group1MoM += group1MemberMedian;
				}
				for(int i = 0; i < group2.length; i++) {
					double group2MemberMedian = Statistics.median(group2Scores.get(align.getChr())[i].get(align));
					matrix.set(alignRowName, group2[i], group2MemberMedian);
					group2MoM += group2MemberMedian;
				}
				if(filterByFDR) {
					pvalList.add(pVal);
				}
				
				group1MoM = group1MoM/(double)group1.length;
				matrix.set(alignRowName,2, group1MoM);
				group2MoM = group2MoM/(double)group2.length;
				matrix.set(alignRowName, 3, group2MoM);
			}
			//Map<Alignments, Double> scores=ttest(alignments);
			//Map<Alignments, double[]> pvalues=computePValues(scores, randomDist);
			//Map<Alignments, Double> FDR=computeFDR(scores, dist, alpha);
			if(filterByFDR) {
				double fdrCutoff = Statistics.FDRCorrect(pvalList, alpha);
				matrix = matrix.filterValuesLargerThanUsingColumn(1, fdrCutoff);
			} else {
				matrix = matrix.filterValuesLargerThanUsingColumn(0, alpha);
			}
		}
		return matrix;
	}
	
	
	private EmpiricalDistribution[] computeRandomScoresMinMax(Set<Alignments> alignments, ContinuousData[] data1, ContinuousData[] data2, String chr, int numRandom){
		ArrayList<Double> minScores=new ArrayList<Double>();
		ArrayList<Double>  maxScores=new ArrayList<Double> ();
		
		ArrayList<Double>  allscores=new ArrayList<Double>();
		
		//double[] temp=this.minMax(scores.values());
		//minScores.add(temp[0]);
		//maxScores.add(temp[1]);
		
		for(int i=0; i<numRandom; i++){
			double min=Double.MAX_VALUE;
			double max=-Double.MAX_VALUE;
			long startTime=System.currentTimeMillis();
			for(Alignments align: alignments){
				double t=getShuffledTStatistc(align, data1, data2);
				if(!new Double(t).equals(Double.NaN)){
					min=Math.min(min, t);
					max=Math.max(max, t);
					allscores.add(t);
				}
				//else{
				//	System.err.println(t);
				//}
			}
			long endTime=System.currentTimeMillis();
			//System.err.println("Finished Permutation "+i+" of "+numRandom+" Took: "+((endTime-startTime)/1000)+"sec");
			minScores.add(min);
			maxScores.add(max);
		}
		
		EmpiricalDistribution maxDist=new EmpiricalDistribution(maxScores, EMPIRICAL_DIST_BINS);
		EmpiricalDistribution minDist=new EmpiricalDistribution(minScores, EMPIRICAL_DIST_BINS);
		EmpiricalDistribution fullDist = new EmpiricalDistribution(allscores, EMPIRICAL_DIST_BINS);
		EmpiricalDistribution[] rtrn={minDist, maxDist, fullDist};
		return rtrn;
	}
	
	/**
	 * This implements our null model: Assume probes have a background intensity distribution,
	 * thus regions that are not differentially expressed would look exactly the same if you swap
	 * the probes with randomly chosen probes across the genome.
	 * @param align
	 * @param data1
	 * @param data2
	 * @return
	 */
	private double getShuffledTStatistc(Alignments align, ContinuousData[] data1, ContinuousData[] data2){
		int[] array=data1[0].getShuffledDataIndex(align.getChr());
				
		int startIndex=align.getStartIndex();
		int endIndex=align.getEndIndex();
		
		
		List<Double> list1=new ArrayList<Double>();
		List<Double> list2=new ArrayList<Double>();
		
		for(int k=0; k<data1.length; k++){
			for(int i=startIndex; i<=endIndex; i++){
				list1.add(data1[k].getDataAtIndex(array[i-startIndex]));
			}
		}
		
		for(int k=0; k<data2.length; k++){
			for(int i=startIndex; i<=endIndex; i++){
				list2.add(data2[k].getDataAtIndex(array[i-startIndex]));
			}
		}
		
		double t = Statistics.tstat(list1, list2);
		//System.err.println("Random: "+list1.size()+" "+list2.size()+" "+t);		
		
		return t;
	}
	
	private TreeMap<Alignments, List<Double>>[] getValuesPerTU(ContinuousData[] data, Set<Alignments> alignments){

		TreeMap<Alignments, List<Double>>[] scores=new TreeMap[data.length];
		for(int i=0; i<scores.length; i++){
			scores[i]=getValuesPerTU(alignments, data[i]);
		}
		return scores;
	}
	
	
	private TreeMap<Alignments, List<Double>> getValuesPerTU(Set<Alignments> alignments, ContinuousData data){
		TreeMap<Alignments, List<Double>> rtrn=new TreeMap<Alignments, List<Double>>();
		for(Alignments align: alignments){
			List<Double> vals = data.getDataForAlignments(align);
			rtrn.put(align, vals);
		}
		return rtrn;
	}	
	
	private ContinuousData[] makeData(LocationAwareMatrix matrix, String [] columns)throws IOException{
		ContinuousData [] data = new ContinuousData[columns.length];
		for(int i=0; i<columns.length; i++){
			data[i]=new ContinuousData(matrix, columns[i]);
		}
		
		return data;
	}
	
	private double ttest(Alignments align) {
		return ttest(align, group1Scores.get(align.getChr()), group2Scores.get(align.getChr()));
	}
	
	private double ttest(Alignments align, Map<Alignments, List<Double>>[] scores1, Map<Alignments, List<Double>>[] scores2 ) {
		List<Double> list1=new ArrayList<Double>();
		List<Double> list2=new ArrayList<Double>();
		
		for(int i=0; i < scores1.length; i++){list1.addAll(scores1[i].get(align));}
		for(int i=0; i < scores2.length; i++){list2.addAll(scores2[i].get(align));}
		
		
		double t = Statistics.tstat(list1, list2);
		return t;
	}
	
	

}
