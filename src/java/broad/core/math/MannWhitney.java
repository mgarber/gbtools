package broad.core.math;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public class MannWhitney {

	double Z;
	double p;
	
	public MannWhitney(double[] x, double[] y){
		Map<Double, List<String>> map=new TreeMap<Double, List<String>>();
		double[] all=new double[x.length+y.length];
		for(int i=0; i<x.length; i++){
			all[i]=x[i];
			if(map.containsKey(all[i])){
				List<String> list=map.get(all[i]);
				list.add("x");
			}
			else{
				List<String> list=new ArrayList<String>();
				list.add("x");
				map.put(all[i], list);
			}
		}
		for(int i=x.length; i<all.length; i++){
			all[i]=y[i-x.length];
			if(map.containsKey(all[i])){
				List<String> list=map.get(all[i]);
				list.add("y");
			}
			else{
				List<String> list=new ArrayList<String>();
				list.add("y");
				map.put(all[i], list);
			}
		}
		List<Double>[] xyRanked=rankOrder(map);
		this.Z=calculateZ(xyRanked, x, y);
		this.p=calculateP(Z);
	}
	
	
	public MannWhitney(ArrayList<Double> x, ArrayList<Double> y){
		Map<Double, List<String>>  map=new TreeMap<Double, List<String>>();
		double[] all=new double[x.size()+y.size()];
		for(int i=0; i<x.size(); i++){
			all[i]=x.get(i);
			if(map.containsKey(all[i])){
				List<String> list= map.get(all[i]);
				list.add("x");
			}
			else{
				List<String> list=new ArrayList<String>();
				list.add("x");
				map.put(all[i], list);
			}
		}
		for(int i=x.size(); i<all.length; i++){
			all[i]=y.get(i-x.size());
			if(map.containsKey(all[i])){
				List<String> list= map.get(all[i]);
				list.add("y");
			}
			else{
				List<String> list=new ArrayList<String>();
				list.add("y");
				map.put(all[i], list);
			}
		}
		List<Double>[] xyRanked=rankOrder(map);
		this.Z=calculateZ(xyRanked, x, y);
		this.p=calculateP(Z);
	}

	private double calculateP(double Z){
		cern.jet.random.Normal norm=new cern.jet.random.Normal(0,1, new cern.jet.random.engine.DRand());
		double cdf=norm.cdf(Z);
		return Math.min(1, Math.min((1-cdf), cdf)*2);
	}
	
	public double getZ(){return this.Z;}
	public double getP(){return this.p;}
	
	private double calculateZ(List<Double>[] xyRanked, double[] x, double[] y){
		double U=Statistics.sum(xyRanked[0]);
		double mu=(x.length*(x.length+y.length+1))/2.0;
		double var=((x.length*y.length)*(x.length+y.length+1))/12.0;
		double sigma=Math.sqrt(var);
		double Z=(U-mu)/sigma;
		return Z;
	}
	
	
	private double calculateZ(List<Double>[] xyRanked, List<Double> x, List<Double> y){
		double U=Statistics.sum(xyRanked[0]);
		double mu=(x.size()*(x.size()+y.size()+1))/2.0;
		double var=((x.size()*y.size())*(x.size()+y.size()+1))/12.0;
		double sigma=Math.sqrt(var);
		double Z=(U-mu)/sigma;
		return Z;
	}
	
	
	protected static List<Double>[] rankOrder(Map<Double, List<String>> map){
		ArrayList<Double> x=new ArrayList<Double>();
		ArrayList<Double> y=new ArrayList<Double>();
		double count=0;
		for(double key: map.keySet()){
			List<String> list= map.get(key);
			double end = count+list.size();
			double rank = rank(count, list);
			//System.err.println(count+" "+start+" "+end+" "+rank);
			count=end;
			for(String str: list){
				if(str.equalsIgnoreCase("x")){
					x.add(rank);
				}
				if(str.equalsIgnoreCase("y")){
					y.add(rank);
				}
			}
			
		}
		//System.err.println(x);
		List[] rtrn={x,y};
		return rtrn;
	}
	

	protected static double rank(double count, List<String> list){
		//For tie (duplicated) values, returns the mean of the ranks of the duplicated values: 
		// for first rank r and 4 duplicated values it should return  r+ (r+1) + (r+2) + (r+3) = 4*r + 1+2+3 =(in general) N*r + N(N-1)/2
		// But the rank passed (count) is one behind so we add 1.
		double N = list.size();
		return (N*(count+1) + (N*(N - 1))/2)/N;
	}
	
	
}
