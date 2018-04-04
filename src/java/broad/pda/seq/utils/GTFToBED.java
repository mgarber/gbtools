package broad.pda.seq.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;

public class GTFToBED {

	public GTFToBED(File gtfFile, String save)throws IOException{
		Map<String, String> strandMap=new HashMap<String, String> ();
		Map<String, Collection> gtf=parseGTF(gtfFile, strandMap);
		write(gtf,strandMap, save);
	}
	
	private Map<String, Collection> parseGTF(File file,Map<String, String>  strandMap)throws IOException{
		Map<String, Collection> rtrn=new TreeMap();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		
    	String nextLine;
            while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
               if(nextLine.startsWith("chr")){
					try{
            	   String[] tokens=nextLine.split("\t");
            	// -1 on start point because gff is 1 shifted while BED is 0 shifted. On the other hand, gff is inclusive while bed is not-> do not substract 1 from the end coordinate
					Alignments align=new Alignments(tokens[0], new Integer(tokens[3])-1, new Integer(tokens[4]));
					String name1=tokens[8].split(";")[1];
					String name=name1.split(" ")[2];
					name=name.replaceAll("\\s", "_");
					name=name.replaceAll("\"", "");
					
					Collection c=new TreeSet();
					if(rtrn.containsKey(name)){c=rtrn.get(name);}
					if (tokens[2].equalsIgnoreCase("exon"))
						c.add(align);
					if (!strandMap.containsKey(name))
						strandMap.put(name, tokens[6]);
					if (tokens[2].equalsIgnoreCase("transcript"))
						strandMap.put(name, tokens[6]);
					rtrn.put(name, c);
					}catch(Exception ex){System.err.println(nextLine);}
               }
            }
            
            
            reader.close();
            return rtrn;
	}
	
	private void write(Map<String, Collection> gtf, Map<String , String > strandMap, String save)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(String transcript: gtf.keySet()){
			Collection<Alignments> exons=gtf.get(transcript);
			RefSeqGene gene=new RefSeqGene(exons);
			gene.setName(transcript);
			gene.setOrientation(strandMap.get(transcript));
			writer.write(gene.toBED()+"\n");
		}
		
		writer.close();
	}
	
	public static void main(String[] args)throws IOException{
		File gtfFile=new File(args[0]);
		String save=args[1];
		new GTFToBED(gtfFile, save);
	}
	
}
