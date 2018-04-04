package broad.pda.seq.joshuaAnalyses;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import broad.core.datastructures.MatrixWithHeaders;

public class BinData {

	public BinData(MatrixWithHeaders fullData, MatrixWithHeaders sampled, int numBins){
		//Determine genes in each bin off of full data
		
		//sort by expression
		MatrixWithHeaders sorted=fullData.sortList(0);
		
		int valsPerBin=sorted.getNumberRows()/numBins;
		
		Map<Integer, List<String>> valsByBin=new TreeMap<Integer, List<String>>();
		
		
		
		int count=0;
		for(String gene: sorted.getRowNames()){
			int binNumber=count%valsPerBin;
			List<String> vals=valsByBin.get(binNumber);
			if(vals==null){vals=new ArrayList();}
			vals.add(gene);
			valsByBin.put(binNumber, vals);
			count++;
		}
		
	}
	
}
