package broad.core.overlaputils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;

public class CountNumOverlapping {

	public CountNumOverlapping(Map<String, Collection<Alignments>> genesToFilter, Map<String, IntervalTree<Alignments>> genesToFilterBy, String save)throws IOException{
		Map<Alignments, Integer> filtered=new TreeMap();
		
		for(String chr: genesToFilter.keySet()){
			System.err.println(chr);
			Collection<Alignments> lincs=genesToFilter.get(chr);
			IntervalTree<Alignments> proteins=genesToFilterBy.get(chr);
			filtered.putAll(filter(lincs, proteins));
		}
		
		write(save, filtered);
		
	}
	
	private Map<Alignments, Integer> filter(Collection<Alignments> lincs, IntervalTree<Alignments> proteins){
		Map<Alignments, Integer> rtrn=new TreeMap();
		if(proteins==null){return rtrn;}
		
		for(Alignments linc: lincs){
			int num=proteins.numOverlappers(linc.getStart(), linc.getEnd());
			rtrn.put(linc, num);
		}
		
		return rtrn;
	}
	
	private void write(String save, Map<Alignments, Integer> map)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(Alignments gene: map.keySet()){
			writer.write(gene+"\t"+map.get(gene)+"\n");
		}
		
		writer.close();
	}
	
	private static Map<String, IntervalTree<Alignments>> makeTree(Collection<Alignments> genes){
		Map<String, IntervalTree<Alignments>> rtrn=new TreeMap<String, IntervalTree<Alignments>>();
		
		for(Alignments align: genes){
			IntervalTree<Alignments> tree=new IntervalTree<Alignments>();
			String chr=align.getChr();
			if(rtrn.containsKey(chr)){tree=rtrn.get(chr);}
			Node node=tree.find(align.getStart(), align.getEnd());
			if(node==null){tree.put(align.getStart(), align.getEnd(), align);}
			else{
				node.incrementCount();
			}
			rtrn.put(align.getChr(), tree);
		}
		
		return rtrn;
		
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
			System.err.println("started");
			Map<String, Collection<Alignments>> genesToFilter=BEDFileParser.loadAlignmentDataByChr(new File(args[0]));
			System.err.println("loaded lincRNAs");
			Map<String, IntervalTree<Alignments>> genesToFilterBy=makeTree(BEDFileParser.loadAlignmentData(new File(args[1])));
			System.err.println("loaded proteins into interval tree");
			String save=args[2];
			new CountNumOverlapping(genesToFilter, genesToFilterBy, save);
			
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=genes to count by (regions) \n args[1]=genes to count (transcripts- Full BED) \n args[2]=save";
	
	
}
