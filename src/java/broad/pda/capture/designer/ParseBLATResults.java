package broad.pda.capture.designer;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;

public class ParseBLATResults {

	//Filter cross-hybing probes from the list
	//For each oligo name write all matches
	public ParseBLATResults(File[] pslFiles, String save) throws IOException{
		Map<String, Collection<Alignments>> matches=new TreeMap<String, Collection<Alignments>>();
		
		for(int i=0; i<pslFiles.length; i++){
			System.err.println(pslFiles[i]);
			Map<String, Collection<Alignments>> map=parse(pslFiles[i]);
			merge(matches, map);
		}
		
		write(save, matches);
	}

	private Map<String, Collection<Alignments>> parse(File file) throws IOException {
		Map<String, Collection<Alignments>> rtrn=new TreeMap<String, Collection<Alignments>>();
		Collection<String> lines=BEDFileParser.loadList(file.getAbsolutePath());
		for(String line: lines){
			try{
			RefSeqGene gene=new RefSeqGene(line);
			Collection<Alignments> set=new TreeSet<Alignments>();
			if(rtrn.containsKey(gene.getName())){set=rtrn.get(gene.getName());}
			set.add(gene.getAlignment());
			rtrn.put(gene.getName(), set);
			}catch(Exception ex){System.err.println("Skipping line "+line);}
		}
		return rtrn;
	}

	private void merge(Map<String, Collection<Alignments>> matches, Map<String, Collection<Alignments>> map) {
		for(String name: map.keySet()){
			Collection<Alignments> list=map.get(name);
			if(matches.containsKey(name)){
				Collection<Alignments> all=matches.get(name);
				all.addAll(list);
				matches.put(name, all);
			}
			else{
				matches.put(name, list);
			}
		}
	}

	private void write(String save, Map<String, Collection<Alignments>> matches) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String probe: matches.keySet()){
			writer.write(probe);
			for(Alignments align: matches.get(probe)){
				writer.write("\t"+align.toUCSC());
			}
			writer.write("\n");
		}
		
		writer.close();
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>1){
			File[] files=new File(args[0]).listFiles();
			String save=args[1];
			new ParseBLATResults(files, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=files \n args[1]=save";
	
}
