package broad.projection.nmf;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.math.Statistics;
import broad.core.util.ParseGCTFile;
import broad.projection.math.EuclidanNMFCost;
import broad.projection.math.Matrix;
import broad.projection.math.NMFCostFunction;
import broad.projection.math.PoissonNMFCost;
import broad.projection.utils.CLSToHinitial;


/**
 * Non-Negative Matrix Factorization
 *
 * @author Daniel Scanfeld
 * @version 1.0
 */

public class NMF
{

	int numIterations=1000;
	int numPermutations=100;
	
	
	// This uses a supervised seeding startegy based on the class labels in the CLS file
	public NMF(File gctFile, File clsFile, String save, int numIterations, double precision, NMFCostFunction cost)throws IOException{
		Array array=new Array(gctFile.getAbsolutePath());
		
		Matrix Hi=CLSToHinitial.parseCLS(clsFile, false);
		
		double[] costSupervised=array.nmf(Hi, numIterations, cost, precision,10, save+".Hi");
		array.writeNMF(save, true);		
	}
	
	//This is the classical NMF decomposition
	public NMF(File gctFile, String save, int k, int numIterations, double precision, NMFCostFunction cost)throws IOException{
		Array array=new Array(gctFile.getAbsolutePath());
		
		array.globalShift();
		array.nmf(k, numIterations, cost, precision);
		
		//cluster genes by group: uses the max metagene and/or pearson distance to metagenes to classify groups
		//Map<String, String> genesByCluster=clusterByMaxMetagene(array);
		//Map<String, String> genesByPearson=clusterByMaxPearson(array);
		
		array.writeNMF(save, true);
	}
	
	private double getCummulativeProbability(ArrayList<Double> vals, double val){
		  double counterTotal=0;
		  double countLessThan=0;
		  
		  for(Double observed: vals){
			  if(observed<val){countLessThan++;}
			  counterTotal++;
		  }
		  
		  return countLessThan/counterTotal;
	  }
	
	
	
	
	private void write(String save, Map<String, Set> genesByCluster, File gctFile)throws IOException{
		FileWriter writer=new FileWriter(save+".gct");
		Map<String, ArrayList> gct=ParseGCTFile.loadData(gctFile);
		Map<String, String> PIDToName=ParseGCTFile.loadPIDName(gctFile);
		
		ArrayList<String> header=gct.get("header");
		
		writer.write("#1.2\n");
		writer.write((gct.size()-1)+"\t"+header.size()+"\n");
		writer.write("PID\tFactor");
		
		for(String desc: header){writer.write("\t"+desc);}
		writer.write("\n");
		
		for(String factor: genesByCluster.keySet()){
			Set<String> list=genesByCluster.get(factor);
			for(String gene: list){
				writer.write(gene+"\t"+factor);
				ArrayList<Double> vals=gct.get(gene);
				for(Double val: vals){writer.write("\t"+val);}
				writer.write("\n");
				
			}
		}	
		
		writer.close();
	}
	
	
	private void write(String save, Map<String, String> genesByCluster, Map<String, String> genesByPearson, File gctFile)throws IOException{
		FileWriter writer=new FileWriter(save+".gct");
		Map<String, ArrayList> gct=ParseGCTFile.loadData(gctFile);
		Map<String, String> PIDToName=ParseGCTFile.loadPIDName(gctFile);
		
		ArrayList<String> header=gct.get("header");
		
		writer.write("#1.2\n");
		writer.write((gct.size()-1)+"\t"+header.size()+"\n");
		writer.write("PID\tFactor\tFactorPearson");
		
		for(String desc: header){writer.write("\t"+desc);}
		writer.write("\n");
		
		for(String gene: gct.keySet()){
			if(!gene.equalsIgnoreCase("header")){
			String factor=genesByCluster.get(gene);
			String factorPearson=genesByPearson.get(gene);
			writer.write(gene+"\t"+factor+"\t"+factorPearson);
			ArrayList<Double> vals=gct.get(gene);
			for(Double val: vals){writer.write("\t"+val);}
			writer.write("\n");
			}
		}
		
		
		writer.close();
	}
	
	
	private Map clusterByMaxMetagene(Array array){
		Map<String, Set> rtrn=new TreeMap();
		
		Map<String, String> map=new TreeMap();
		
		Matrix WMatrix=array.x.W;
		
		//for each gene cluster by max metagene
		for(int i=0; i<WMatrix.nr; i++){
			double[] vals=WMatrix.x[i];
			int index=this.getMaxIndex(vals);
			String name=array.geneName.get(i).toString();
			String cluster="F"+(index+1);
			
			Set set=new TreeSet();
			if(rtrn.containsKey(cluster)){
				set=rtrn.get(cluster);
			}
			//set.add(i);
			set.add(array.geneName.get(i).toString());
			rtrn.put(cluster, set);
			map.put(array.geneName.get(i).toString(), cluster);
		}
		
		return map;
	}
	
	private Map clusterByMaxPearson(Array array){
		Map<String, Set> rtrn=new TreeMap();
		
		Map<String, String> map=new TreeMap();
		
		Matrix WMatrix=array.x.W;
		
		//for each gene cluster by max metagene
		for(int i=0; i<WMatrix.nr; i++){
			//double[] vals=WMatrix.x[i];
			double[] vals=array.x.x[i];
			int index=getMaxPearson(array, vals);
			String name=array.geneName.get(i).toString();
			String cluster="F"+(index+1);
			
			Set set=new TreeSet();
			if(rtrn.containsKey(cluster)){
				set=rtrn.get(cluster);
			}
			//set.add(i);
			set.add(array.geneName.get(i).toString());
			rtrn.put(cluster, set);
			map.put(array.geneName.get(i).toString(), cluster);
		}
		
		return map;
	}
	
	//return list of indexes for top num scores
	private int getMaxIndex(double[]vals){
		
		double[] sorted=(vals.clone());
		Arrays.sort(sorted);
		
		double max=sorted[sorted.length-1];
		
		
		for(int i=0; i<vals.length; i++){
			double val=vals[i];
			if(max==val){return(i);}
		}
		
		return -1;
	}
	
	//return list of indexes for top num scores
	private int getMaxPearson(Array array, double[]vals){
		
		Matrix HMatrix=array.x.H;
		double[] pearson=new double[HMatrix.nr];
		for(int i=0; i<HMatrix.nr; i++){
			double[] HVals=HMatrix.x[i];
			pearson[i]=Statistics.pearsonDistance(vals, HVals);
		}
		
		int index=maxIndex(pearson);
		
		return index;
	}
	
	
	private int maxIndex(double[] vals){
		int maxIndex=0;
		double maxVal=vals[0];
		for(int i=0; i<vals.length; i++){
			if(vals[i]>maxVal){maxIndex=i; maxVal=vals[i];}
		}
		return maxIndex;
	}
	
	private void write(String save, double[] vals)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(int i=0; i<vals.length; i++){
			writer.write(vals[i]+"\n");
		}
		
		writer.close();
	}
	
	
	public static void main(String[] args)throws IOException{
		if(args.length>4){
		File file=new File(args[0]);
		String save=args[1];
		int k=new Integer(args[2]);
		int numIterations=new Integer(args[3]);
		double precision=new Double(args[4]);
		NMFCostFunction cost = args.length > 5 && "euclidian".equals(args[5]) ? new EuclidanNMFCost() : new PoissonNMFCost();
		//if(args.length>6){new NMF(file, new File(args[6]), save, k, numIterations, precision, cost);}
		new NMF(file, save, k, numIterations, precision, cost);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=gct file \n args[1]=save file \n args[2]=k \n args[3]=num iterations \n args[4]=precision \n args[5]=cost function \n args[6]=cls file (optional)"; 
}