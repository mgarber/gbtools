package broad.pda.seq.pairedend;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;

public class EstimateInsertSize {

	int cutoff=300;
	
	public EstimateInsertSize(Collection<Alignments> segments, Collection<Alignments> pairs, String save, String chr)throws IOException{
		IntervalTree<Alignments> exons=makeExonTree(segments, chr);
		Map<Alignments, Integer> inferredSize=new TreeMap();
		
		for(Alignments pair: pairs){
			Iterator<Node<Alignments>> iter=exons.overlappers(pair.getStart(), pair.getEnd());
			int size=estimateSize(iter, pair);
			inferredSize.put(pair, size);
		}
		
		write(save, inferredSize);
		write(save+".genomic", inferredSize, cutoff);
	}
	
	
	private IntervalTree<Alignments> makeExonTree(Collection<Alignments> segments, String chr){
		IntervalTree<Alignments> tree=new IntervalTree();
		
		for(Alignments align: segments){
			if(align.getChromosomeString().equalsIgnoreCase(chr)){
			tree.put(align.getStart(), align.getEnd(), align);
			}
		}
		
		return tree;
	}
	
	private int estimateSize(Iterator<Node<Alignments>> iter, Alignments pair){
		int sum=0;
		
		//only count bases in between insert
		while(iter.hasNext()){
			Alignments align=iter.next().getValue();
			Alignments region=new Alignments(align.getChr(), Math.max(pair.getStart(), align.getStart()), Math.min(pair.getEnd(), align.getEnd()));
			sum+=(region.getSize());
		}
		
		return sum;
	}
	
	private void write(String save, Map<Alignments, Integer> map)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(Alignments align: map.keySet()){
			int val=map.get(align);
			writer.write(align+"\t"+val+"\t"+align.getSize()+"\n");
		}
		
		writer.close();
	}
	
	private void write(String save, Map<Alignments, Integer> map, int cutoff)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(Alignments align: map.keySet()){
			int val=map.get(align);
			if(val<cutoff){writer.write(align+"\t"+align.getSize()+"\n");}
		}
		
		writer.close();
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>3){
			Collection<Alignments> segments=BEDFileParser.loadAlignmentData(new File(args[0]));
			Collection<Alignments> inserts=BEDFileParser.loadAlignmentData(new File(args[1]));
			String save=args[2];
			String chr=args[3];
			new EstimateInsertSize(segments, inserts, save, chr);
			
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=segments (BED File) \n args[1]=Inserts (BED File) \n args[2]=save \n args[3]=chr";
}
