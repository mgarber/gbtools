package broad.pda.geneexpression;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.math.EmpiricalDistribution;
import broad.core.math.Statistics;
import broad.core.util.ParseGCTFile;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.rnaseq.misc.GetNeighbors;

public class ComputeNeighborCorrelations {
	
	int numRandom=100;

	public ComputeNeighborCorrelations(File lincExpressionFile, File geneExpressionFile, File geneMappingFile, String save)throws IOException{
		Map lincExpression=ParseGCTFile.loadData(lincExpressionFile);
		Map geneExpression=ParseGCTFile.loadData(geneExpressionFile);
		
		Map<Alignments, String> lincMapping=makeLincMapping(lincExpression);
		Map<Alignments, String> geneMapping=makeGeneMapping(geneMappingFile);
		
		// original code : Map<Alignments, Alignments[]> neighbors=new GetNeighbors(lincMapping.keySet(), geneMapping.keySet()).getNeighbors();
		GetNeighbors NeighborsObj= new GetNeighbors(lincMapping.keySet(), geneMapping.keySet());
		Map<Alignments, Alignments[]> neighbors=NeighborsObj.getNeighbors();
		Map<Alignments, Alignments[]>[] randomNeighbors=generateRandomNeighbors(lincMapping.keySet(), geneMapping.keySet(), numRandom);
		
		Map<Alignments, double[]> correlations=computeCorrelations(neighbors, lincMapping, geneMapping, lincExpression, geneExpression);
		Map<Alignments, double[]>[] randomCorrelations=computeCorrelations(randomNeighbors, lincMapping,  geneMapping, lincExpression, geneExpression);
		
		Map<Alignments, double[]> pvalues=computePValues(correlations, randomCorrelations);
		
		write(save, correlations, randomCorrelations, pvalues);
	}
	
	
	private Map computePValues(Map<Alignments, double[]> correlations, Map<Alignments, double[]>[] randomCorrelations){
		Map rtrn=new TreeMap();
		
		EmpiricalDistribution[] dists=computeDist(randomCorrelations);
		for(Alignments align: correlations.keySet()){
			double[] corrs=correlations.get(align);
			
			double[] pvals={1-dists[0].getCummulativeProbability(corrs[0]), 1-dists[1].getCummulativeProbability(corrs[1])};
			rtrn.put(align, pvals);
		}
				
		return rtrn;
	}
	
	private EmpiricalDistribution[] computeDist(Map<Alignments, double[]>[] randomCorr){
		ArrayList left=new ArrayList();
		ArrayList right=new ArrayList();
		
		for(Alignments align: randomCorr[0].keySet()){
			double leftMax=Double.MAX_VALUE;
			double rightMax=Double.MAX_VALUE;
			for(int i=0; i<randomCorr.length; i++){
				leftMax=Math.min(randomCorr[i].get(align)[0], leftMax);
				rightMax=Math.min(randomCorr[i].get(align)[1], rightMax);
			}
			left.add(leftMax);
			right.add(rightMax);
		}
		
		EmpiricalDistribution[] rtrn={new EmpiricalDistribution(left), new EmpiricalDistribution(right)};
		return rtrn;
	}
	
	private void write(String save, Map<Alignments, double[]> correlations, Map<Alignments, double[]>[] randomCorrelations)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(Alignments linc: correlations.keySet()){
			double[] corrs=correlations.get(linc);
			writer.write(linc.toUCSC()+"\tLeft\t"+corrs[0]);
			for(int i=0; i<randomCorrelations.length; i++){writer.write("\t"+randomCorrelations[i].get(linc)[0]);}
			writer.write("\n");
			
			writer.write(linc.toUCSC()+"\tRight\t"+corrs[1]);
			for(int i=0; i<randomCorrelations.length; i++){writer.write("\t"+randomCorrelations[i].get(linc)[1]);}
			writer.write("\n");
		}
		
		writer.close();
	}
	
	private void write(String save, Map<Alignments, double[]> correlations, Map<Alignments, double[]>[] randomCorrelations, Map<Alignments, double[]> pvalues)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(Alignments linc: correlations.keySet()){
			double[] corrs=correlations.get(linc);
			double[] pVals=pvalues.get(linc);
			writer.write(linc.toUCSC()+"\tLeft\t"+pVals[0]+"\t"+corrs[0]);
			for(int i=0; i<randomCorrelations.length; i++){writer.write("\t"+randomCorrelations[i].get(linc)[0]);}
			writer.write("\n");
			
			writer.write(linc.toUCSC()+"\tRight\t"+pVals[1]+"\t"+corrs[1]);
			for(int i=0; i<randomCorrelations.length; i++){writer.write("\t"+randomCorrelations[i].get(linc)[1]);}
			writer.write("\n");
		}
		
		writer.close();
	}
	
	
	private Map makeLincMapping(Map<String, ArrayList> expression){
		Map rtrn=new TreeMap();
		
		for(String name: expression.keySet()){
			if(!name.equalsIgnoreCase("header")){rtrn.put(new Alignments(name), name);}
		}
		
		return rtrn;
	}
	
	private Map makeGeneMapping(File file)throws IOException{
		Set<Alignments> alignments=BEDFileParser.loadAlignmentData(file);
		Map rtrn=new TreeMap();
		
		for(Alignments align: alignments){
			rtrn.put(align, align.getName());
		}
		
		return rtrn;
	}
	
	private Map[] computeCorrelations(Map<Alignments, Alignments[]>[] randomNeighbors, Map<Alignments, String> lincMapping, Map<Alignments, String> geneMapping, Map lincExpression, Map geneExpression){
		Map[] rtrn=new Map[randomNeighbors.length];
		
		for(int i=0; i<rtrn.length; i++){
			rtrn[i]=computeCorrelations(randomNeighbors[i], lincMapping, geneMapping, lincExpression, geneExpression);
		}
		
		return rtrn;
	}
	
	private Map[] generateRandomNeighbors(Set<Alignments> lincs, Set<Alignments> genes, int numRandom){
		Map[] rtrn=new Map[numRandom];
		
		for(int i=0; i<numRandom; i++){
			rtrn[i]=generateRandomNeighbors(lincs, genes);
		}
		
		return rtrn;
	}
	
	
	private Map generateRandomNeighbors(Set<Alignments> lincs, Set<Alignments> genes){
		Map rtrn=new TreeMap();
		
		Map<String, IntervalTree> trees=GetNeighbors.makeIntervalTree(genes);
		
		for(Alignments linc: lincs){
			IntervalTree tree=trees.get(linc.getChr());
			int index=new Double(Math.random()*(tree.size()-1)).intValue();
			Node<Alignments> node=tree.findByIndex(index);
			Alignments[] array={node.getValue(), node.getNext().getValue()};
			rtrn.put(linc, array);
		}
		
		return rtrn;
	}
	
	private Map computeCorrelations(Map<Alignments, Alignments[]> neighbors, Map<Alignments, String> lincMapping, Map<Alignments, String> geneMapping, Map<String, ArrayList> lincExpression, Map<String, ArrayList> geneExpression){
		Map rtrn=new TreeMap();
		
		for(Alignments linc: neighbors.keySet()){
			String lincName=lincMapping.get(linc);
			Alignments leftName=neighbors.get(linc)[0];
			String geneName1="";
			if(leftName!=null){
				geneName1=geneMapping.get(leftName);
			}
			if(geneName1==null){geneName1="";}
			//System.err.println("Name: "+leftName+" "+geneName1);
			
			Alignments rightName=neighbors.get(linc)[1];
			String geneName2="";
			if(rightName!=null){
				geneName2=geneMapping.get(rightName);
			}
			if(geneName2==null){geneName2="";}
			ArrayList lincExpressionVals=lincExpression.get(lincName);
			ArrayList geneExpressionValsLeft=geneExpression.get(geneName1);
			ArrayList geneExpressionValsRight=geneExpression.get(geneName2);
			double pearsonLeft=pearson(lincExpressionVals, geneExpressionValsLeft);
			double pearsonRight=pearson(lincExpressionVals, geneExpressionValsRight);
			double[] array={pearsonLeft, pearsonRight};
			rtrn.put(linc, array);
		}
		
		return rtrn;
	}
	
	private double pearson(ArrayList list1, ArrayList list2){
		if(list1==null || list2==null){return -999;}
		return Statistics.pearsonDistance(list1, list2);
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>3){
			File lincExpressionFile=new File(args[0]);
			File geneExpressionFile=new File(args[1]);
			File geneMappingFile=new File(args[2]);
			String save=args[3];
			new ComputeNeighborCorrelations(lincExpressionFile, geneExpressionFile, geneMappingFile, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=lincExpression file \n args[1]=gene expression file \n args[2]=gene mapping file \n args[3]=save";
	
}
