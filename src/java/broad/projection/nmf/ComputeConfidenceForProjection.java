package broad.projection.nmf;

import java.io.IOException;
import java.util.ArrayList;

import broad.core.math.EmpiricalDistribution;
import broad.projection.math.NMFCostFunction;


public class ComputeConfidenceForProjection {

	double precision=.000001;
	int numPerms=1000;
	NMFCostFunction cost;
	
	public ComputeConfidenceForProjection(Array geneExpression, Array originalWMatrix, String save, int numPerm, NMFCostFunction cost, double precision){
		this.cost=cost;
		this.precision=precision;
		this.numPerms=numPerm;
		Array projection=new ProjectBack(geneExpression, originalWMatrix, precision, cost).getProjection();
		projection.writeNMF(save+".original", true);
		
		Array pvalues=assignPValuesToProjection(projection, originalWMatrix, geneExpression);
		pvalues.writeNMF(save+".pvaluesByMetagene", true);
		
		//pvalues=this.assignOverallPValues(projection, dist);
		//pvalues.writeNMF(save+".pvaluesOverall", true);
	}
	
	
	
	private Array assignPValuesToProjection(Array projection, Array originalWMatrix, Array geneExpression){
		ArrayList vals=new ArrayList();
		ArrayList[] valsByMetagene=new ArrayList[projection.x.H.nr];
		for(int i=0; i<valsByMetagene.length; i++){valsByMetagene[i]=new ArrayList();}
		
		for(int i=0; i<numPerms; i++){
		//generate random gene expression profiles and project
			Array randomized=generateRandomExpression(geneExpression);
			//randomized.write(save+"."+i+".expression.gct", true);
			Array randomProjection=new ProjectBack(randomized, originalWMatrix, precision, this.cost).getProjection();
			//randomProjection.writeNMF(save+"."+i+".projection", true);
			double val=randomProjection.HMax();
			double[] maxByMetagene=randomProjection.HMaxByMetagene();
			vals.add(val);
			for(int j=0; j<maxByMetagene.length; j++){valsByMetagene[j].add(maxByMetagene[j]);}
		}
		
		//Add the observed to the distributions
		double[] maxByMetagene=projection.HMaxByMetagene();
		for(int j=0; j<maxByMetagene.length; j++){valsByMetagene[j].add(maxByMetagene[j]);}
		
	
		Array pvalues=this.assignPValuesByMetagene(projection, valsByMetagene);
		return pvalues;
	}
	
	private EmpiricalDistribution[] computeDists(ArrayList[] lists){
		EmpiricalDistribution[] rtrn=new EmpiricalDistribution[lists.length];
		
		for(int i=0; i<lists.length; i++){rtrn[i]=new EmpiricalDistribution(lists[i], 200);}
		
		return rtrn;
	}
	
	//assign p-values per metagene so that any devition from random is captured even if there is some dominant effect
	private Array assignPValuesByMetagene(Array projection, EmpiricalDistribution[] distsByMetagene){
		Array pvalues=projection.assignPValues(distsByMetagene);
		return pvalues;
	}
	
	//assign p-values per metagene so that any devition from random is captured even if there is some dominant effect
	private Array assignPValuesByMetagene(Array projection, ArrayList[] distsByMetagene){
		Array pvalues=projection.assignPValues(distsByMetagene);
		return pvalues;
	}
	
	//assign p-values overall values
	private Array assignOverallPValues(Array projection, EmpiricalDistribution dist){
		Array pvalues=projection.assignPValues(dist);
		return pvalues;
	}
	
	private Array generateRandomExpression(Array geneExpression){
		//for each gene randomly swap the expression values by sample
		return geneExpression.randomizeArray();
	}
	
	public static void main(String[] args)throws IOException{
		/*if(args.length>2){
		Array expression=new Array(args[0]);
		Array W=new Array(args[1]);
		String save=args[2];
		int numPerm=new Integer(args[3]);
		new ComputeConfidenceForProjection(expression, W, save, numPerm, cost);
		}
		else{System.err.println(usage);}*/
	}
	
	static String usage=" args[0]=expression file  \n args[1]=original W Matrix \n args[2]=save file \n args[3]=num perm";
}
