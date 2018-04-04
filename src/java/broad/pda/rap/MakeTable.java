package broad.pda.rap;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import broad.pda.annotation.BEDFileParser;

public class MakeTable {

	
	public static void main(String[] args)throws IOException{
		if(args.length>1){
			File[] files=new File(args[0]).listFiles();
			String save=args[1];
			
			Map<String, String>[] maps=new Map[files.length];
			
			for(int i=0; i<files.length; i++){
				Collection<String> lines=BEDFileParser.loadList(files[i].getAbsolutePath(), true);
				Map<String, String> map=makeMap(lines);
				maps[i]=map;
			}
			
			writeTable(save, maps, files);
			
		}
		else{System.err.println(usage);}
	}
	
	private static Map<String, String> makeMap(Collection<String> lines) {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		for(String line: lines){
			String[] tokens=line.split("\t");
			rtrn.put(tokens[0], line);
		}
		
		return rtrn;
	}

	private static void writeTable(String save, Map<String, String>[] maps, File[] files) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("Chromosome\tLength");
		for(int i=0; i<files.length; i++){
			writer.write("\t"+files[i].getName());
		}
		writer.write("\n");
		
		for(String chr: maps[0].keySet()){
			String[] tokens=maps[0].get(chr).split("\t");
			writer.write(chr+"\t"+tokens[1]);
			for(int i=0; i<maps.length; i++){
				writer.write("\t"+maps[i].get(chr).split("\t")[3]);
			}
			writer.write("\n");
		}
		
		writer.close();
	}

	static String usage="args[0]=dir \n args[1]=save";
	
}
