package broad.projection.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Set;
import java.util.TreeSet;
import java.util.Vector;

import broad.projection.math.Matrix;


//Given a cls file construct an initial H matrix
public class CLSToHinitial {

	public CLSToHinitial(File clsFile, String save)throws IOException{
		
		Matrix matrix=parseCLS(clsFile, false);
		//Matrix m=new Matrix(matrix);
		write(save, matrix);
	}
	
	
	
	private void write(String save, Matrix m)throws IOException{
		FileWriter writer=new FileWriter(save);
	
		writer.write(m.toString());
		
		writer.close();
	}
	
	/*public static Matrix parseCLS(File clsFile)throws IOException{
		Vector<String> lines=readFile(clsFile);
		int numPhenotypes=0;
		int numSamples=0;
		try{
			numPhenotypes=new Integer((String)lines.get(0).split(" ")[1]);
			numSamples=new Integer((String)lines.get(0).split(" ")[0]);
		}catch(ArrayIndexOutOfBoundsException ex){
			numPhenotypes=new Integer((String)lines.get(0).split("\t")[1]);
			numSamples=new Integer((String)lines.get(0).split("\t")[0]);
		}
		//double[][] array=new double[numPhenotypes][numSamples];
		
		Object[] phenotypes=getPhenotypes(lines.get(2));
		double[][] matrix=makeMatrix(lines.get(2), phenotypes);
		return new Matrix(matrix);
	}*/
	
	
	public static Matrix parseCLS(File clsFile, boolean possibleCombos)throws IOException{
		Vector<String> lines=readFile(clsFile);
		int numPhenotypes=0;
		int numSamples=0;
		try{
			numPhenotypes=new Integer((String)lines.get(0).split(" ")[1]);
			numSamples=new Integer((String)lines.get(0).split(" ")[0]);
		}catch(ArrayIndexOutOfBoundsException ex){
			numPhenotypes=new Integer((String)lines.get(0).split("\t")[1]);
			numSamples=new Integer((String)lines.get(0).split("\t")[0]);
		}
		
		System.err.println("Number phenotypes "+numPhenotypes+" Number of samples "+numSamples);
		
		//double[][] array=new double[numPhenotypes][numSamples];
		String[] names=split(lines.get(2));
		
		Object[] phenotypes=getPhenotypes(names);
		System.err.println("Phenotypes "+phenotypes.length);
		double[][] matrix=makeMatrix(names, phenotypes, possibleCombos);
		return new Matrix(matrix);
	}
	
	public static Matrix permuteCLS(File clsFile, boolean possibleCombos)throws IOException{
		Vector<String> lines=readFile(clsFile);
		int numPhenotypes=0;
		int numSamples=0;
		try{
			numPhenotypes=new Integer((String)lines.get(0).split(" ")[1]);
			numSamples=new Integer((String)lines.get(0).split(" ")[0]);
		}catch(ArrayIndexOutOfBoundsException ex){
			numPhenotypes=new Integer((String)lines.get(0).split("\t")[1]);
			numSamples=new Integer((String)lines.get(0).split("\t")[0]);
		}
		
		System.err.println("Number phenotypes "+numPhenotypes+" Number of samples "+numSamples);
		
		//double[][] array=new double[numPhenotypes][numSamples];
		String[] names=split(lines.get(2));
		
		Object[] phenotypes=getPhenotypes(names);
		String[] permutedNames=shuffle(names);
		System.err.println("Phenotypes "+phenotypes.length);
		double[][] matrix=makeMatrix(permutedNames, phenotypes, possibleCombos);
		return new Matrix(matrix);
	}
	
	private static String[] shuffle(String[] names){
		String[] rtrn=new String[names.length]; 
		ArrayList<String> list=new ArrayList();
		for(int i=0; i<names.length; i++){list.add(names[i]);}
		
		for(int i=0; i<rtrn.length; i++){
			int index=new Double(Math.random()*list.size()).intValue();
			rtrn[i]=list.remove(index);
		}
		
		return rtrn;
	}
	
	private static Vector readFile(File file)throws IOException{
		
			FileInputStream fileInput;
			BufferedReader buf = null;
			Vector temp = new Vector(2000, 500);
			String aLine = new String("");
		    try{
		      fileInput = new FileInputStream(file);
		      buf = new BufferedReader(new InputStreamReader (fileInput));
		    } catch (IOException ex){
		       return temp;
		    }

		    while(true) {
		      try {
		        aLine = buf.readLine();
		        if(aLine == null) break;
		        temp.add(aLine);
		      } catch (IOException e) {
		        temp.removeAllElements();
		        return temp;
		      }
		    }
		    
		    return temp;
		  
	}
	
	private static double[][] makeMatrix(String[] samples, Object[] phenotypes){
		double[][] rtrn=new double[phenotypes.length][samples.length];
		//for each phenotype
		for(int i=0; i<phenotypes.length; i++){
			double[] vals=new double[samples.length];
			for(int j=0; j<samples.length; j++){
				vals[j] = 0.0;
				if(samples[j].equalsIgnoreCase((String)phenotypes[i])){vals[j]=1;};
			}
			rtrn[i]=vals;
		}
		return rtrn;
	}
	
	private static double[][] makeMatrix(String[] samples, Object[] phenotypes, boolean possibleCombo){
		if(!possibleCombo){return makeMatrix(samples, phenotypes);}
		double[][] rtrn=new double[phenotypes.length][samples.length];
		ArrayList<double[]> list=new ArrayList();
		
		//reference phenotype
		for(int i=0; i<phenotypes.length; i++){
			//additional phenotype
			for(int j=i; j<phenotypes.length; j++){
				double[] vals=new double[samples.length];
				for(int k=0; k<samples.length; k++){
					vals[k] = 0.0;
					if(samples[k].equalsIgnoreCase((String)phenotypes[i])){vals[k]=1;}
					else if(samples[k].equalsIgnoreCase((String)phenotypes[j])){vals[k]=1;}
				}
				list.add(vals);
			}
		}
		
		rtrn=convertToMatrix(list);
		
		return rtrn;
	}
	
	private static double[][] convertToMatrix(ArrayList<double[]> list){
		double[][] rtrn=new double[list.size()][];
		for(int i=0; i<list.size(); i++){
			rtrn[i]=list.get(i);
		}
		return rtrn;
	}
	
	
	private static String[] split(String str){
		String[] temp=str.split(" ");
		if(temp.length==1){temp=str.split("\t");}
		return temp;
	}
	
	private static Object[] getPhenotypes(String[] temp){
		Set set=new TreeSet();
		
		for(int i=0; i<temp.length; i++){set.add(temp[i]);}
		
		System.err.println(set);
		return set.toArray();
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>1){
			File file=new File(args[0]);
			String save=args[1];
			new CLSToHinitial(file, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=file \n args[1]=save";
}
