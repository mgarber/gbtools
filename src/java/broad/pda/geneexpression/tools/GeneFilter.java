package broad.pda.geneexpression.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.geneexpression.dge.DGE;

public class GeneFilter {

	static final String usage = "Usage: GeneFilter -task <task name> "+
			"\n\tselectGenes " + 
			"\n\t\t-in matrix "+
			"\n\t\t-out output"+
			"\n\t\t-quantile "+
			
			"\n\t trimMatrix"+
			"\n\t\t-matrix matrix "+
			"\n\t\t-in genes"+
			"\n\t\t-out ";
	
	public GeneFilter(ArgumentMap argmap) throws IOException{
		
		if ("selectGenes".equalsIgnoreCase(argmap.getTask())){
			
			selectGenes(argmap);
		}
		else if ("trimMatrix".equalsIgnoreCase(argmap.getTask())){
			
			trimMatrix(argmap);
		}
	}
	
	public void selectGenes(ArgumentMap argmap) throws ParseException, IOException{
		/*
		 * Read the matrix
		 */
		MatrixWithHeaders matrix = new MatrixWithHeaders(argmap.getInput());
		
		Set<String> selectedGenes = new HashSet<String>();
		
		//for each column
		for(String c:matrix.getColumnNames()){
			//get the 50% quantile gene names
			List<Double> col = DGE.array2List(matrix.getColumn(c));
			double quant = Statistics.quantile(col, argmap.getDouble("quantile"));
			
			//Now go through entire list of genes, and add to set
			for(String r:matrix.getRowNames()){
				if(matrix.get(r, c)>=quant){
					selectedGenes.add(r);
				}
			}
		}
				
		MatrixWithHeaders newMat = new MatrixWithHeaders(new ArrayList<String>(selectedGenes),matrix.getColumnNames());
		//Make a new matrix with the selected Genes
		for(String r:selectedGenes){
			for(String c:matrix.getColumnNames()){
				newMat.set(r, c, matrix.get(r, c));
			}
		}
		newMat.write(argmap.getOutput());
		
		BufferedWriter bf = new BufferedWriter(new FileWriter(argmap.getOutput()+".geneNames"));
		for(String r:selectedGenes){
			bf.write(r);
			bf.newLine();
		}
		bf.close();
	}
	
	public void trimMatrix(ArgumentMap argmap) throws ParseException, IOException{
		
		/*
		 * Read the matrix
		 */
		MatrixWithHeaders matrix = new MatrixWithHeaders(argmap.get("matrix"));
		
		List<String> genes = new ArrayList<String>();
		//Read the genes
		BufferedReader br = new BufferedReader(new FileReader(argmap.getInput()));
		String str;
		while((str = br.readLine())!=null){
			genes.add(str);
		}
		MatrixWithHeaders newMat = new MatrixWithHeaders(genes,matrix.getColumnNames());
		for(String r:genes){
			if(matrix.getRowNames().contains(r))
				for(String c:matrix.getColumnNames()){
					newMat.set(r, c, matrix.get(r, c));
				}
		}
		newMat.write(argmap.getOutput());
		
	}
	
	public static void main (String [] args) throws IOException{
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"selectGenes");
		GeneFilter dummy = new GeneFilter(argMap);
	}
	
}

