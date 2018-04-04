package broad.pda.datastructures;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.annotation.LightweightGenomicAnnotation;

public class ContinuousData {

	private double[] data;
	int[] locations;
	Map<String, int[]> locationMap;
	Map<String, double[]> dataMap;
	ArrayList<Double> allData;
	Map<String, double[]> maxMinByChr;
	
	
	public ContinuousData() {

	}
	
	public ContinuousData(File asciiFile)throws IOException{
		Map<String, Map<Integer, Double>> map=readAndParseAscii(asciiFile);
		this.locationMap=new TreeMap<String, int[]>();
		this.dataMap=new TreeMap<String, double[]>();
		
		for(String chr: map.keySet()){
			System.err.println(chr);
			makeArray(map.get(chr));
			locationMap.put(chr, locations);
			dataMap.put(chr, data);
		}
		
	}
	
	/**
	 * Construct continuous data from a LocationAwareMatrix column
	 * @param matrix
	 * @param column
	 */
	public ContinuousData(LocationAwareMatrix matrix, String column) {
		Map<String, Map<Integer, Double>> map = new HashMap<String, Map<Integer, Double>>();
		allData=new ArrayList<Double>();
		List<String> rows = matrix.getRowNames();
		for(String row : rows) {
			LightweightGenomicAnnotation rowPosition = matrix.getRowPosition(row);
			String chr = rowPosition.getChromosome();
			Map<Integer, Double> chrValMap = map.get(chr);
			if(chrValMap == null) {
				chrValMap = new TreeMap<Integer, Double>();
				map.put(chr, chrValMap);
			}
			chrValMap.put(rowPosition.getStart(), matrix.get(row, column));
			allData.add(matrix.get(row, column));
		}
		
		this.locationMap=new TreeMap<String, int[]>();
		this.dataMap=new TreeMap<String, double[]>();
		
		for(String chr: map.keySet()){
			//System.err.println(chr);
			makeArray(map.get(chr));
			locationMap.put(chr, locations);
			dataMap.put(chr, data);
		}
	}
	
	
	public double[] getMaxMin(String chr, int windowSize){
		double[] data=(double[])this.dataMap.get(chr);
		double[] temp=data.clone();
		Arrays.sort(temp);
		double min=0;
		double max=0;
		for(int i=0; i<windowSize; i++){
			min+=temp[i];
			max+=temp[temp.length-1-i];
		}
		double[] rtrn={min, max};
		return rtrn;
	}
	
	private void makeArray(Map<Integer, Double> map){
		this.locations=new int[map.size()];
		this.data=new double[map.size()];
		
		int i=0;
		for(Integer location: map.keySet()){
			Double value=map.get(location);
			locations[i]=location;
			data[i]=value;
			i++;
		}
	}
	
	private Map<String, Map<Integer, Double>> readAndParseAscii(File file)throws IOException{
		String aLine = new String("");
	    FileInputStream fileInput = new FileInputStream(file);
	    BufferedReader buf = new BufferedReader(new InputStreamReader (fileInput));
	    
	    this.allData=new ArrayList<Double>();
	    
	    Map<String, Map<Integer, Double>> rtrn=new TreeMap<String, Map<Integer, Double>>();

	    while(true) {
	      try {
	    	Map<Integer, Double> temp=new TreeMap<Integer, Double>();
	        aLine = buf.readLine();
	        if(aLine == null) break;
	        String[] tokens=aLine.split("\t");
	        String chr=tokens[0].split(":")[0];
	        if(rtrn.containsKey(chr)){temp=(Map<Integer, Double>)rtrn.get(chr);}
	        int location=Integer.parseInt(tokens[0].split(":")[1].split("-")[0]);
	        double value=new Double(tokens[1]);
	        this.allData.add(value);
	        temp.put(location, value);
	        rtrn.put(chr, temp);
	      } catch (IOException e) {}
	    }
	    
	    return rtrn;
		
	}

	public int[] getBinary(double threshold){
		int[] rtrn=new int[this.data.length];
		for(int i=0; i<rtrn.length; i++){
			rtrn[i]=0;
			if(data[i]>=threshold){rtrn[i]=1;}
		}
		return rtrn;
	}
	
	public Set<String> getChromosomes(){return this.locationMap.keySet();}
	
	public int[] getLocations(String chr){return this.locationMap.get(chr);}
	public double[] getData(String chr){return this.dataMap.get(chr);}
	public long getKey(int index){return (this.locations[index]);}
	
	
	public int getBinaryData(String chr, double threshold, int index){
		double val=((double[])this.dataMap.get(chr))[index];
		if(val<threshold){return 0;}
		else{return 1;}
	}

	
	public double[] resize(double[] array, int newSize) {
        double [] newArray = new double[newSize];
        int elementsToCopy = Math.min(array.length, newSize);
        System.arraycopy(array, 0, newArray, 0, elementsToCopy);
        return newArray;
    }
	
	public double[] getShuffledData(){
		double[] rtrn=this.data;
		int N = rtrn.length;
	    for (int i = 0; i < N; i++) {
	    	int r = i + uniform(N-i);     // between i and N-1
	    	double temp = rtrn[i];
	    	rtrn[i] = rtrn[r];
	    	rtrn[r] = temp;
	       }
	    return rtrn;
	}
	
	
	public List<Double> getShuffledData(String chr){
		int len= dataMap.get(chr).length;
		List<Double> rtrn=new ArrayList<Double>(len);
		Random r = new Random();
		for(int i=0; i<len; i++){
			int idx = r.nextInt(allData.size()-1);
			rtrn.add(allData.get(idx));
		}
		return rtrn;
	}
	
	public int[] getShuffledDataIndex(String chr){
		int len=dataMap.get(chr).length;
		int[] rtrn=new int[len];
		Random r = new Random();
		for(int i=0; i<rtrn.length; i++){
			rtrn[i]=r.nextInt(allData.size() - 1);
		}
		return rtrn;
	}
	
	public Set<Alignments> makeContiguousRegions( String chr, int maxDistanceToMerge){
		Set<Alignments>  rtrn = new TreeSet<Alignments> ();
		int[] locations = getLocations(chr);
		
		for(int i=0; i<locations.length-1; i++){
			long current=locations[i];
			long next=locations[i+1];
			Set<Integer> set=new TreeSet<Integer>();
			set.add(i);
			while((next-current)<maxDistanceToMerge){
				set.add(i+1);
				i++;
				if(i<locations.length-1){
					current=locations[i];
					next=locations[i+1];
				}
				else{break;}
			}
			Alignments align=collapse(set, chr);
			rtrn.add(align);
		}
		return rtrn;
	}
	
	public Set<Alignments> makeContiguousRegions(String chr, Collection<? extends LightweightGenomicAnnotation> regions){
		Set<Alignments> rtrn=new TreeSet<Alignments>();
		
		int[] locations=getLocations(chr);
		
		for(LightweightGenomicAnnotation peak: regions){
			if(peak.getChromosome().equalsIgnoreCase(chr)){
				//System.err.println(peak);
				for(int i=0; i<locations.length; i++){
					int current=locations[i];
					Set<Integer> set=new TreeSet<Integer>();
					Alignments temp=new Alignments(chr, current, current);
					while(peak.overlaps(temp)){
						set.add(i);
						i++;
						if(i<locations.length){
							current=locations[i];
							temp=new Alignments(chr, current, current);
						}
						else{break;}
					}
					if(set.size()>0){
						Alignments align=collapse(set, chr);
						align.setStrand(peak.getOrientation());
						rtrn.add(align);
					}
				}
			}
		}
		return rtrn;
	}
	
	public Alignments collapse(Set<Integer> set, String chr){
		return new Alignments(chr, (Integer)set.toArray()[0], (Integer)set.toArray()[set.size()-1], this);
	}
	
	public Alignments collapse(Set<Alignments> set){
		Set<Integer> startSet=new TreeSet<Integer>();
		Set<Integer> endSet=new TreeSet<Integer>();
		String chr="";
		for(Alignments align: set){
			startSet.add(align.getStartIndex());
			endSet.add(align.getEndIndex());
			chr=align.getChr();
		}
		return new Alignments(chr, (Integer)startSet.toArray()[0], (Integer)endSet.toArray()[endSet.size()-1],this);
	}
	
	public void setDataForAlignment(Alignments align){
		align.setScores(getDataForAlignments(align));
	}
	
	public List<Double> getDataForAlignments(Alignments align){
		
		int startIndex=align.getStartIndex();
		int endIndex=align.getEndIndex();
		
		
		List<Double> values = new ArrayList<Double>(endIndex-startIndex+1);
		//System.err.println("Getting data for align " + align.toUCSC() + " with chr " + align.getChr());
		for(int i=startIndex; i<=endIndex; i++){
			double [] chrData = getData(align.getChr());
			values.add(chrData[i]);
		}			
		return values;
	}
	
	public double getDataAtIndex(int index){return (Double)this.allData.get(index);}
	
	public static int uniform(int N) {
        return (int) (Math.random() * N);
    }
	
	public long [] resize(long[] array, int newSize) {
        long [] newArray = new long[newSize];
        int elementsToCopy = Math.min(array.length, newSize);
        System.arraycopy(array, 0, newArray, 0, elementsToCopy);
        return newArray;
    }
	
}
