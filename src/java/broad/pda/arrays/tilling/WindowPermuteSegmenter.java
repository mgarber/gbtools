package broad.pda.arrays.tilling;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.annotation.MaximumContiguousSubsequence;
import broad.core.math.EmpiricalDistribution;
import broad.core.math.MathUtil;
import broad.core.math.Statistics;
import broad.pda.datastructures.Alignments;
import broad.pda.datastructures.ContinuousData;

//need to add observed as 1 permutation

public class WindowPermuteSegmenter {

	int maxDistance;
	
	/**
	 * Keeps contiguous regions as defined by the maxDistance parameter.
	 */
	Map<String, Set<Alignments>> chromosomeContiguousRegions;
	ContinuousData data;
	
	public static final int[] DEFAULT_WINDOWS={5,7,9,10,15,20};
	private static double FUDGE_FACTOR = 0.1;
	public static int MAX_MERGE_DISTANCE = 1000;
		
	public WindowPermuteSegmenter(ContinuousData data,  int maxMergeDistance)throws IOException{
		this.maxDistance = maxMergeDistance;
		chromosomeContiguousRegions = new LinkedHashMap<String, Set<Alignments>>(data.getChromosomes().size());
		this.data = data;
		
		//Map scanPValue=new TreeMap();
		//Map scoreMap=new TreeMap();		
		for(String chr : data.getChromosomes()) {
			System.err.println("Initializing data for chromosome <"+chr+">");
			chromosomeContiguousRegions.put(chr, data.makeContiguousRegions(chr, maxMergeDistance));
			//System.err.println("Made Contiguous Regions");
			//Map<Alignments, double[]> windowsByTU=getValuesPerTU(alignments, data, chr);
			setValuesPerTU(chr);
			//System.err.println("Set Values Per TU");
		}

	}


	public Map<String, Map<Alignments, Double>> segment( int numTails, int [] windowSizes, int numRandom, double alpha) throws IOException {
		Map<String, Map<Alignments, Double>> chrSegments = new LinkedHashMap<String, Map<Alignments,Double>>(data.getChromosomes().size());
		for(String chr : data.getChromosomes()) {
			chrSegments.put(chr, segment(numTails, windowSizes, numRandom, alpha, chr));
		}
		
		return chrSegments;
	}

	public Map<Alignments, Double> segment(int numTails, int[] windowSizes, int numRandom, double alpha, String chr) throws IOException {
		System.err.println("Segmenting " + chr);
		Map<Alignments, Double> segments=new TreeMap<Alignments, Double>();
		long lastTime = System.currentTimeMillis();
		for(int i=0; i<windowSizes.length; i++){	
			System.err.println("\twindow size "+windowSizes[i]+" number " + (i+1) + " out of "+ windowSizes.length);
			//System.setOut(new PrintStream("scan.scores"));
			Map<Alignments, Double > scores = computeScoresPerTU(windowSizes[i], chr);
			long t = System.currentTimeMillis();
			System.err.println("\t\tComputed scores per TU in " + ((t - lastTime)/1000));
			lastTime = t;
			//scoreMap.putAll(scores);
				
			//Map[] randomScoreMap=computeRandomScores(alignments, data, chr, this.numRandom);
			double[][] randomScores = computeRandomScoresMinMax(chr, numRandom, windowSizes[i]);
			t = System.currentTimeMillis();
			System.err.println("\t\tComputed random scores min-max " + ((t - lastTime)/1000));
			lastTime = t;
			System.err.println("Permuted");
			//EmpiricalDistribution[] distributions=getPermutationDistribution(randomScoreMap, data, chr);
			EmpiricalDistribution[] distributions=getPermutationDistribution(randomScores);
			t = System.currentTimeMillis();
			System.err.println("\t\tComputed permutation distributions (" + numRandom + ") in "+ ((t - lastTime)/1000) + " sec.");
			lastTime = t;
			Map<Alignments, double []> pvalMap=computePValue(scores, distributions[0], distributions[1]);
			t = System.currentTimeMillis();
			System.err.println("\t\tComputed pvalues " + ((t - lastTime)/1000));
			lastTime = t;

			//filter by alpha
			Set<Alignments> filtered=filter(pvalMap, numTails, alpha);
			//Merge All overlapping
			t = System.currentTimeMillis();
			System.err.println("\t\tFiltered by alpha " + ((t - lastTime)/1000));
			lastTime = t;			
			Set<Alignments> merged=merge(filtered);
			t = System.currentTimeMillis();
			System.err.println("\t\tMerged segments " + ((t - lastTime)/1000));
			lastTime = t;
			//trim ends
			Map<Alignments, Double> truncatedExons=trimEnds(merged);
			t = System.currentTimeMillis();
			System.err.println("\t\tTruncated exons " + ((t - lastTime)/1000));
			lastTime = t;
			System.err.println("\tGot " + truncatedExons.size() + " segments");
			segments.putAll(truncatedExons);
			t = System.currentTimeMillis();
			System.err.println("\t\tDone with window " + ((t - lastTime)/1000));
			lastTime = t;
		}
		return segments;
	}
	
	
	
	
	private Map<Alignments, Double> trimEnds(Set<Alignments> merged){
		Map<Alignments, Double> rtrn=new TreeMap<Alignments, Double>();
		for(Alignments align: merged){
			Object[] locationsAndData=this.getDataForAlignments(align, data);
			int[] locations=(int[])locationsAndData[0];
			double[] values=(double[])locationsAndData[1];
			double[] coordinates = MaximumContiguousSubsequence.maxSubSum3(values);
			if(coordinates[0]>0 && coordinates[1]!=coordinates[2]){
				Alignments truncAlign=new Alignments(align.getChr(),locations[new Double(coordinates[1]).intValue()], locations[new Double(coordinates[2]).intValue()]);
				double avgExpression=getAvgExpression(values, coordinates);
				rtrn.put(truncAlign, avgExpression);
			}
			else{
				//double avgExpression=getAvgExpression(values);
				//rtrn.put(align, avgExpression);
			}
		}
		return rtrn;
	}
	
	private Set<Alignments> merge(Set<Alignments> filtered){
		Set<Alignments> rtrn=new TreeSet<Alignments>();
		List<Alignments> array=new ArrayList<Alignments>(filtered);
		for(int i=0; i < array.size() - 1; i++){
			Alignments current=array.get(i);
			Alignments next=array.get(i+1);
			Set<Alignments> set=new TreeSet<Alignments>();
			set.add(current);
			while(current.overlaps2(next)){
				set.add(next);
				i++;
				if(i<array.size()-1){
					current=array.get(i);
					next=array.get(i+1);
				}
				else{break;}
			}
			Alignments align=data.collapse(set);
			rtrn.add(align);
		}
		return rtrn;
		
	}
	
	private Set<Alignments> filter(Map<Alignments, double[]> pvalues, int numTails, double alpha){
		Set<Alignments> rtrn=new TreeSet<Alignments>();
		/*try {
			System.setErr(new PrintStream("pval.filter"));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}*/
		for(Alignments align: pvalues.keySet()){
			double[] ps = pvalues.get(align);
			double p = getPValue(ps, numTails);
			//System.err.println("alpha " + alpha + " align pval " + p + " adding it? " + (p>=alpha));
			if(p>=alpha){
				rtrn.add(align);
			}
		}
		return rtrn;
	}
	
	
	
	private double getPValue(double[] ps, int numTails){
		double p=1;
		if(numTails==1){p=ps[0];}
		else{p=Math.max(ps[0],ps[1]);}	
		return p;
	}
	
	private double[][] computeRandomScoresMinMax( String chr, int numRandom, int windowSize){
		double[] minScores=new double[numRandom];
		double[] maxScores=new double[numRandom];
		Set<Alignments> alignments = chromosomeContiguousRegions.get(chr);
		if(alignments != null) {
			for(int i=0; i<numRandom; i++){
				double min=Double.MAX_VALUE;
				double max=-Double.MAX_VALUE;
				//long startTime=System.currentTimeMillis();
				/*try {
					System.setOut(new PrintStream("perm" + i + ".scores"));
				} catch (FileNotFoundException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}*/
				for(Alignments align: alignments){
					List<Double> vals=this.getShuffledDataForAlignments(align);
					Map<Alignments, Double> scores = scan( windowSize, align, vals); //entire permutation distribution
					min=Math.min(min, Statistics.min(scores.values()));
					max=Math.max(max, Statistics.max(scores.values()));
				}
				long endTime=System.currentTimeMillis();
				//System.err.println("Finished Permutation "+i+" of "+numRandom+" Took: "+((endTime-startTime)/1000)+"sec");
				minScores[i]=min;
				maxScores[i]=max;
			}
		}
		
		double[][] rtrn={minScores, maxScores};
		return rtrn;
	}
	
	/**
	 * Shuffles probe values in chromosomes and scans a window of desired size.
	 * @param chr
	 * @param windowSize
	 * @return
	 */
	public List<Double> computeRandomScores(String chr, int windowSize) {
		List<Double> rtrn = new ArrayList<Double>();
		Set<Alignments> alignments = chromosomeContiguousRegions.get(chr);
		for(Alignments align: alignments){
			List<Double> regionData =this.getShuffledDataForAlignments(align);
			Map<Alignments, Double> scores = scan( windowSize, align, regionData); //entire permutation distribution
			rtrn.addAll(scores.values());
		}		
		
		return rtrn;
	}
	
	private Map<Alignments, Double > computeScoresPerTU(int windowSize,  String chr){
		Map<Alignments, Double > rtrn=new TreeMap<Alignments, Double >();
		Set<Alignments> contiguousAlignments = chromosomeContiguousRegions.get(chr);
		if(contiguousAlignments != null) {
			for(Alignments align: contiguousAlignments){
				rtrn.putAll(scan( windowSize, align));
			}
		}
		return rtrn;
	}
	
	private void setValuesPerTU(String chr) {
		//System.err.println("setting values per TU, for chr " + chr + " chrContiguous regions chrs " + chromosomeContiguousRegions.keySet());
		Set<Alignments> alignments = chromosomeContiguousRegions.get(chr);
		for(Alignments align: alignments){
			data.setDataForAlignment(align);
		}
	}
	
	
	private EmpiricalDistribution[] getPermutationDistribution(Map[] randomScoreMap, ContinuousData data, String chr){
		
		double max=-1;
		double min=99;
		
		double minMax=-1;
		double minMin=99;
		
		double[] maxDist=new double[randomScoreMap.length];
		double[] minDist=new double[randomScoreMap.length];
		for(int i=0; i<randomScoreMap.length; i++){
			Map map=randomScoreMap[i];
			maxDist[i]=Statistics.max(map.values());
			minDist[i]=Statistics.min(map.values());
			max=Math.max(maxDist[i], max);
			min=Math.min(maxDist[i], min);
			
			minMax=Math.max(minDist[i], minMax);
			minMin=Math.min(minDist[i], minMin);
		}
		
		
		EmpiricalDistribution distMax=new EmpiricalDistribution(200, min, max);
		EmpiricalDistribution distMin=new EmpiricalDistribution(200, minMin, minMax);
		for(int i=0; i<maxDist.length; i++){
			distMax.add(maxDist[i]);
			distMin.add(minDist[i]);
		}
		
		EmpiricalDistribution[] rtrn={distMax, distMin};
		return rtrn;
		
	}
	
	private EmpiricalDistribution[] getPermutationDistribution(double[][] randomScoresMinMax){
		double[] randomScoresMin=randomScoresMinMax[0];
		double[] randomScoresMax=randomScoresMinMax[1];
		
		double minMax=-Double.MAX_VALUE;
		double minMin=Double.MAX_VALUE;
		
		double maxMax=-Double.MAX_VALUE;
		double maxMin=Double.MAX_VALUE;
		
		
		for(int i=0; i<randomScoresMin.length; i++){
			maxMax=Math.max(randomScoresMax[i], maxMax);
			maxMin=Math.min(randomScoresMax[i], maxMin);
			minMax=Math.max(randomScoresMin[i], minMax);
			minMin=Math.min(randomScoresMin[i], minMin);
		}
		
		
		EmpiricalDistribution distMax=new EmpiricalDistribution(200, maxMin, maxMax);
		EmpiricalDistribution distMin=new EmpiricalDistribution(200, minMin, minMax);
		for(int i=0; i<randomScoresMax.length; i++){
			distMax.add(randomScoresMax[i]);
			distMin.add(randomScoresMin[i]);
		}
		
		EmpiricalDistribution[] rtrn={distMax, distMin};
		return rtrn;
		
	}
	
			
	private Map<Alignments, double []> computePValue(Map<Alignments, Double > scoreMap, EmpiricalDistribution distMax, EmpiricalDistribution distMin){
		Map<Alignments, double []> map=new TreeMap<Alignments, double []>();
		/*try {
			System.setOut(new PrintStream("pvals.txt"));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}*/
		//cern.jet.random.Empirical emp=new cern.jet.random.Empirical(dist, cern.jet.random.Empirical.LINEAR_INTERPOLATION, new cern.jet.random.engine.DRand());
		for(Alignments key: scoreMap.keySet()){
			double val=scoreMap.get(key);
			double pMax=distMax.getCummulativeProbability(val);
			double pMin=1-distMin.getCummulativeProbability(val);
			double[] p={pMax, pMin};
			//System.out.println(key.toUCSC() + "\t" + val + " pMax " + pMax + " pMin " + pMin );
			map.put(key, p);
		}
		return map;
	}
		
	private double getPermutationMax(List<Double> vals, int width, Alignments align){
		vals=MathUtil.shuffle(vals);
		Map<Alignments, Double>  scoreMap=this.scan(width, align,vals);
		double max=0;
		for(Object key: scoreMap.keySet()){
			double val=(Double)scoreMap.get(key);
			//System.err.println(val);
			max=Math.max(val, max);
		}
		return max;
	}
	
	private double[] getPermutationDistribution(List<Double> vals, int width, ContinuousData data, String chr, Alignments align){
		vals=MathUtil.shuffle(vals);
		Map<Alignments, Double>  scoreMap=this.scan(width, align, vals);
		double[] rtrn=new double[scoreMap.size()];
		int i=0;
		for(Object key: scoreMap.keySet()){
			rtrn[i]=(Double)scoreMap.get(key);
			//System.err.println(val);
			//max=Math.max(val, max);
			i++;
		}
		
		
		return rtrn;
		
	}
	
	private Map<Alignments, Double> scan( int width, Alignments region){
		List<Double> vals = region.getScores();
		return scan(width, region, vals);
		
	}
	
	private Map<Alignments, Double> scan( int width, Alignments region, List<Double> vals){
		Map<Alignments, Double> scoreMap = new TreeMap<Alignments, Double>();
		for(int i=0; i<vals.size() - width; i++){
			double score=getWindowScore(i, (i+width), vals);
			//Alignments align=new Alignments(chr, (region.getStartIndex()+i), (region.getStartIndex()+Math.min(i+width, vals.length-1)), data);
			Alignments align=new Alignments(region.getChr(), (region.getStartIndex()+i), (region.getStartIndex()+i+width), data);
			//if(i%100000 == 0){System.out.println(align.toUCSC()+"\t"+score);}
			scoreMap.put(align, score);
		}

		for(int i=vals.size() - 1; i>=width; i--){
			double score=getWindowScore(i-width, i, vals);
			//Alignments align=new Alignments(chr, (region.getStartIndex()+Math.max(0, i-width)), (region.getStartIndex()+i), data);
			Alignments align=new Alignments(region.getChr(), (region.getStartIndex()+ i-width), (region.getStartIndex()+i), data);
			scoreMap.put(align, score);
		}
		
		return scoreMap;
		
	}
	
	private double getWindowScore(int start, int end, List<Double> vals) {
		return getSum( start, end, vals);
		//return getMedianOverVariance(start, end, vals);
		//return getMaxOverVariance(start, end, vals);
	}
	
	private double getMedianOverVariance(int start, int end, List<Double> vals) {
		if(vals.size() <= 3) {return 0;}
		List<Double> valList = new ArrayList<Double>(end -start);
		for(int j=start; j<=end; j++){
			if(j>=0 && j<vals.size()){valList.add(vals.get(j));}
		}
		
		Collections.sort(valList);
		double median = Statistics.median(valList);
		double sd = Statistics.variance(valList);
		
		return median/Math.max(sd,FUDGE_FACTOR);
	}
	
	private double getSum(int start, int end, List<Double> vals){
		double sum=0;
		for(int j=start; j<=end; j++){
			if(j>=0 && j<vals.size()){sum+=vals.get(j);}
		}
		return sum;
	}
	
	private double getMinOverVariance(int start, int end, List<Double> vals) {

		List<Double> valList = new ArrayList<Double>(end -start);
		for(int j=start; j<=end; j++){
			if(j>=0 && j<vals.size()){valList.add(vals.get(j));}
		}
		double sd = Statistics.variance(valList);
		
		return Statistics.min(valList)/Math.max(sd,FUDGE_FACTOR);		
	}
	
	private double getMaxOverVariance(int start, int end, List<Double> vals) {

		List<Double> valList = new ArrayList<Double>(end -start);
		for(int j=start; j<=end; j++){
			if(j>=0 && j<vals.size()){valList.add(vals.get(j));}
		}
		double sd = Statistics.variance(valList);
		
		return Statistics.max(valList)/Math.max(sd,FUDGE_FACTOR);		
	}
	
	
	private double[] getDataForAlignments(Alignments align, ContinuousData data, String chr){
		
		int startIndex=align.getStartIndex();
		int endIndex=align.getEndIndex();
		
		
		double[] rtrn=new double[endIndex-startIndex+1];
		
		for(int i=startIndex; i<=endIndex; i++){
			rtrn[i-startIndex]=data.getData(chr)[i];
		}
		
				
		return rtrn;
	}
	
	private Object[] getDataForAlignments(Alignments align, ContinuousData data){
		
		int startIndex=align.getStartIndex();
		int endIndex=align.getEndIndex();
		
		
		double[] rtrn=new double[endIndex-startIndex+1];
		int[] locations=new int[endIndex-startIndex+1];
		
		for(int i=startIndex; i<=endIndex; i++){
			rtrn[i-startIndex]=data.getData(align.getChr())[i];
			locations[i-startIndex]=data.getLocations(align.getChr())[i];
		}
		
		Object[] array={locations, rtrn};
		return array;
	}
	
	
	private List<Double> getShuffledDataForAlignments(Alignments align){
		List<Double> suffledDataVals =data.getShuffledData(align.getChr());
		int startIndex=align.getStartIndex();
		int endIndex=align.getEndIndex();
		
		
		List<Double> rtrn=new ArrayList<Double>(endIndex-startIndex+1);
		
		for(int i=startIndex; i<=endIndex; i++){
			rtrn.add(suffledDataVals.get(i));
		}
		
				
		return rtrn;
	}
	
	private void write(Map map, Map scoreMap, String save, int numTails)throws IOException{
		FileWriter writer=new FileWriter(save);
		for(Object align: map.keySet()){
			double[] ps=(double[])map.get(align);
			double score=(Double)scoreMap.get(align);
			double p=this.getPValue(ps, numTails);		
			writer.write(align.toString()+"\t"+p+"\t"+score+"\n");
		}
		writer.close();
	}
	
	private void write(Map map, Map scoreMap, Map truncatedExons, String save, int numTails, double alpha)throws IOException{
		FileWriter writer=new FileWriter(save);
		for(Object align: map.keySet()){
			Object truncatedAlign=truncatedExons.get(align);
			double[] ps=(double[])map.get(align);
			double score=(Double)scoreMap.get(align);
			double p=this.getPValue(ps, numTails);		
			if(p>=alpha){writer.write(align.toString()+"\t"+truncatedAlign+"\t"+p+"\t"+score+"\n");}
			
			//for(int i=0; i<vals.length; i++){writer.write("\t"+vals[i]);}
			//writer.write("\n");
		}
		writer.close();
	}
	
	
	public void writeMerged(Map<Alignments, Double> map, String save)throws IOException{
		FileWriter writer=new FileWriter(save);
		for(Alignments align: map.keySet()){
			writer.write(align+"\t"+map.get(align)+"\n"); //add p-values
		}
		writer.close();
	}
	
	private double getAvgExpression(double[] values, double[] coordinates){
		int start=new Double(coordinates[1]).intValue();
		int end=new Double(coordinates[2]).intValue();
		
		double sum=0;
		int count=0;
		for(int i=start; i<end; i++){
			sum+=values[i];
			count++;
		}
		return sum/count;
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>5){
		File file=new File(args[0]); //gff File
		String save=args[1]; //Save File
		int numTails=new Integer(args[2]); //Number of tails just high or high and low
		int numRandom=new Integer(args[3]); //Number of permutations
		double alpha=new Double(args[4]); //P-value cutoff (1-pvalue)
		String chr=args[5]; //chromosome to analyze
		int[] windowSizes=DEFAULT_WINDOWS;
		if(args.length>6){
			windowSizes=new int[args.length-6];
			for(int i=6; i<args.length; i++){
				windowSizes[i-6]=new Integer(args[i]);
			}
		}
		//if(args.length>5){
		
		// Not very efficient, but not worse than before I think, I switch chromosome filtering to data loading.
		WindowPermuteSegmenter wps = new WindowPermuteSegmenter(new ContinuousData(file), MAX_MERGE_DISTANCE);
		Map<Alignments, Double> segments = wps.segment(numTails, windowSizes, numRandom, alpha ,chr);
		wps.writeMerged(segments, save+".segments");
		//}
		//else{new SlidingWindowPermuteByExperimentMemoryEfficient(new ContinuousData(file), save, numTails, numRandom, alpha);}
		}
		else{System.err.println("function: "+function+"\nUSAGE: \n"+usage);}
	}
	
	static String function="Computes enriched expression in sliding windows across the genome for various sized windows";
	static String usage="args[0]=file to analyze \nargs[1]=file to save \nargs[2]=number of tails to test \nargs[3]=number of permutations \nargs[4]=alpha \nargs[5]=chromosome \n args[6]=window Sizes (tab delimited)";
	
}
