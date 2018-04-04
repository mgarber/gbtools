package broad.projection.nmf;

import java.io.File;
import java.io.IOException;

import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.projection.math.EuclidanNMFCost;
import broad.projection.math.NMFCostFunction;
import broad.projection.math.PoissonNMFCost;



public class NMFMain {

	public static void main(String[] args)throws IOException{
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE);
		if (argMap.getTask().equalsIgnoreCase("Decompose")){decomposeMain(argMap);}
		else if(argMap.getTask().equalsIgnoreCase("DecomposeSupervised")){decomposeSupervisedMain(argMap);}
		else if(argMap.getTask().equalsIgnoreCase("Project")){projectMain(argMap);}
		else if(argMap.getTask().equalsIgnoreCase("Significance")){significanceMain(argMap);}
		else{System.err.println(helpString+"\n"+USAGE);}
	}
	
	private static void decomposeMain(ArgumentMap argMap)throws IOException{
		File gctFile=new File(argMap.getMandatory("gctFile"));
		String save=argMap.getMandatory("saveFile");
		int k=new Integer(argMap.getMandatory("numMetagenes"));
		int numIterations=new Integer(argMap.getMandatory("numIterations"));
		double precision=new Double(argMap.getMandatory("precision"));
		String cost=argMap.getMandatory("costFunction");
		
		NMFCostFunction costFunction = "euclidian".equals(cost) ? new EuclidanNMFCost() : new PoissonNMFCost();
		
		new NMF(gctFile, save, k, numIterations, precision, costFunction);
	}
	
	private static void decomposeSupervisedMain(ArgumentMap argMap)throws IOException{
		File gctFile=new File(argMap.getMandatory("gctFile"));
		String save=argMap.getMandatory("saveFile");
		int numIterations=new Integer(argMap.getMandatory("numIterations"));
		double precision=new Double(argMap.getMandatory("precision"));
		String cost=argMap.getMandatory("costFunction");
		File clsFile=new File(argMap.getMandatory("clsFile"));
		
		NMFCostFunction costFunction = "euclidian".equals(cost) ? new EuclidanNMFCost() : new PoissonNMFCost();
		
		new NMF(gctFile, clsFile, save, numIterations, precision, costFunction);
	}
	
	
	private static void projectMain(ArgumentMap argMap)throws IOException{
		Array expression=new Array(argMap.getMandatory("gctFile"));
		Array originalW=new Array(argMap.getMandatory("wFile"));
		String save=argMap.getMandatory("saveFile");
		double precision=new Double(argMap.getMandatory("precision"));
		String cost=argMap.getMandatory("costFunction");
		
		NMFCostFunction costFunction = "euclidian".equals(cost) ? new EuclidanNMFCost() : new PoissonNMFCost();
		ProjectBack projection=new ProjectBack(expression, originalW, precision, costFunction);
		projection.write(save);
	}
	
	private static void significanceMain(ArgumentMap argMap){
		Array expression=new Array(argMap.getMandatory("gctFile"));
		Array originalW=new Array(argMap.getMandatory("wFile"));
		String save=argMap.getMandatory("saveFile");
		double precision=new Double(argMap.getMandatory("precision"));
		String cost=argMap.getMandatory("costFunction");
		int numPerm=new Integer(argMap.getMandatory("numPerm"));
		
		NMFCostFunction costFunction = "euclidian".equals(cost) ? new EuclidanNMFCost() : new PoissonNMFCost();
		
		ComputeConfidenceForProjection confidence=new ComputeConfidenceForProjection(expression, originalW, save, numPerm, costFunction, precision);
	}
	
	private static String usageLine(){
		String USAGE="-TASK \n Task=Decompose \narguments: -gctFile -saveFile -numMetagenes -numIterations -precision -costFunction+\n";
		USAGE+="Task=DecomposeSupervised \narguments: -gctFile -saveFile -clsFile -numIterations -precision -costFunction\n";
		USAGE+="Task=Project \narguments: -gctFile -wFile -saveFile -precision -costFunction\n";
		USAGE+="Task=Significance \nargumnets: -gctFile -wFile -saveFile -precision -costFunction -numPerm\n";
		USAGE+="Task=help";
		return USAGE;
	}
	
	static String helpString="-gctFile=input expression file (gct format) \n -saveFile=output file \n -numMetagenes=number of factors to decompose matrix into \n -numIterations=number of times to ressed the random decomposition \n -precision=the precision of the decomposition at which to break the current iteration \n -costFunction=Euclidian or Poisson cost \n -clsFile=class label file to use for semi-supervised decomposition \n -numPerm=the number of permutations to perform";
	static String USAGE=usageLine();
	
}
