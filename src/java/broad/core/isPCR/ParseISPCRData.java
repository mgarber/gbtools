package broad.core.isPCR;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

public class ParseISPCRData {

	public ParseISPCRData(File[] files, String save)throws IOException{
		Map[] maps=new Map[files.length];
		for(int i=0; i<maps.length; i++){maps[i]=parse(files[i]);}
		write(save, maps);
	}
	
	private void write(String save, Map<String, Set>[] maps)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(int i=0; i<maps.length; i++){
			for(String gene: maps[i].keySet()){
				Set<String> set=maps[i].get(gene);
				writer.write(gene+"\t"+set.size()+"\t");
				for(String str: set){writer.write(str+",");}
				writer.write("\n");
			}
		}
		
		writer.close();
	}
	
	private void write(String save, Map<String, Set> map)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(String gene: map.keySet()){
			Set<String> set=map.get(gene);
			writer.write(gene+"\t"+set.size()+"\t");
			for(String str: set){writer.write(str+",");}
			writer.write("\n");
		}
		
		writer.close();
	}
	
	private Map parse(File file)throws IOException{
		Map rtrn=new TreeMap();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
          String[] tokens=nextLine.split(" ");
          if(nextLine.startsWith(">")){
        	  String key=tokens[1]+"\t"+tokens[3]+"\t"+tokens[4];
        	  String aligned=tokens[0].replaceAll(">", "");
        	  Set set=new TreeSet();
        	  if(rtrn.containsKey(key)){set=(Set)rtrn.get(key);}
        	  set.add(aligned);
        	  rtrn.put(key, set);
          }
          
          
        }
        reader.close();
        
		
		return rtrn;
	}
	
	public static void main(String[] args)throws IOException{
		File[] files=new File(args[0]).listFiles();
		String save=args[1];
		new ParseISPCRData(files, save);
	}
	
}
