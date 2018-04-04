package broad.projection.nmf;

import java.io.IOException;

import broad.projection.math.Matrix;
import broad.projection.math.NMFCostFunction;
import broad.projection.math.PoissonNMFCost;



public class ProjectBack {


	double precision=0.0001;
	NMFCostFunction DEFAULT_COST_FUNCTION = new PoissonNMFCost();
	Array projection;
	
	/*public ProjectBack(File expressionFile, File originalWMatrix, String save)throws IOException{
		Array a=new Array(expressionFile.getAbsolutePath());
		Array projection = a.project(originalWMatrix.getAbsolutePath());
		projection.write(save, true);
	}*/
	
	public ProjectBack(Array expression, Array originalW, double precision, NMFCostFunction cost){
		this.precision=precision;
		this.projection = projectByFixingNMF(expression, originalW, cost, precision);
	}
		
	public Array getProjection(){return this.projection;}
	
	public void write(String save)throws IOException{projection.writeNMF(save, true);}
	
	/*private Array projectByFixingNMF(File gctFile, File WFile, NMFCostFunction costFunction, double precision){
		Array array=new Array(gctFile.getAbsolutePath());
		array.globalShift();
		
		Array WArray=new Array(WFile.getAbsolutePath());
		Matrix WMatrix=WArray.x;
		
		array.projectFixedW(WArray, precision);
		//System.err.println("Cost: " + cost);
		return array;
	}*/
	
	private Array projectByFixingNMF(Array array, Array WArray, NMFCostFunction costFunction, double precision){
		//array.globalShift();
		
		Matrix WMatrix=WArray.x;
		
		array.projectFixedW(WArray, precision, costFunction);
		//System.err.println("Cost: " + cost);
		return array;
	}

	public static void main(String[] args)throws IOException{
		/*if(args.length>3){
			Array expression=new Array(args[0]);
			Array originalW=new Array(args[1]);
			String save=args[2];
			double precision=new Double(args[3]);
			ProjectBack projection=new ProjectBack(expression, originalW, precision);
			projection.write(save);
			}
			else{System.err.println(usage);}*/
	}
	static String usage=" args[0]=gct file \n args[1]=original W Matrix file \n args[2]=save file \n args[3]=precision"; 
}
