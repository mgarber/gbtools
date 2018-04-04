package broad.pda.rnaseq.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;

import broad.pda.seq.alignment.Pair;

public class FilterRunOnPairsByBLAST {

	public FilterRunOnPairsByBLAST(Collection<Pair<String>> runOnPairs, Collection<Pair<String>> parologPairs, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(Pair<String> runOn: runOnPairs){
			boolean hasParalog=parolog(runOn, parologPairs);
			if(!hasParalog){
				writer.write(runOn.getValue1()+"\t"+runOn.getValue2()+"\n");
			}
		}
		
		writer.close();
	}

	private boolean parolog(Pair<String> runOn,	Collection<Pair<String>> parologPairs) {
		for(Pair<String> paralogs: parologPairs){
			if(runOn.getValue1().equalsIgnoreCase(paralogs.getValue1()) ||runOn.getValue1().equalsIgnoreCase(paralogs.getValue2())){
				if(runOn.getValue2().equalsIgnoreCase(paralogs.getValue1()) ||runOn.getValue2().equalsIgnoreCase(paralogs.getValue2())){
					return true;
				}
			}
		}
		return false;
	}
	
	public static void main(String[] args) throws IOException{
		Collection<Pair<String>> runOns=parsePairs(new File(args[0]));
		Collection<Pair<String>> paralogs=parsePairs(new File(args[1]));
		String save=args[2];
		new FilterRunOnPairsByBLAST(runOns, paralogs, save);
	}

	private static Collection<Pair<String>> parsePairs(File file) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		
    	Collection<Pair<String>> rtrn=new ArrayList();
    	String nextLine;
    	int i=0;
    	while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
    		String[] tokens=nextLine.split("\t");
    		Pair<String> pair=new Pair<String>(tokens[1], tokens[2]);
    		rtrn.add(pair);
    	}
    	
    	reader.close();
    	return rtrn;
	}
}
