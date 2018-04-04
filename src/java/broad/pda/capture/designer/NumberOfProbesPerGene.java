package broad.pda.capture.designer;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import broad.pda.annotation.BEDFileParser;

public class NumberOfProbesPerGene {

	public NumberOfProbesPerGene(String reportFile, String save) throws IOException{
		Collection<String> list=BEDFileParser.loadList(reportFile);
		Map<String, Collection<String>> map=new TreeMap<String, Collection<String>>();
		for(String line: list){
			String name=line.split("\t")[0];
			Collection<String> set=new ArrayList<String>();
			if(map.containsKey(name)){set=map.get(name);}
			set.add(line);
			map.put(name, set);
		}
		
		FileWriter writer=new FileWriter(save);
		for(String name: map.keySet()){
			writer.write(name+"\t"+map.get(name).size()+"\n");
		}
		writer.close();
	}
	
	
	public static void main(String[] args) throws IOException{
		String file=args[0];
		String save=args[1];
		new NumberOfProbesPerGene(file, save);
	}
	
}
