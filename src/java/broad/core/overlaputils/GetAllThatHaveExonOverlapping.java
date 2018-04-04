package broad.core.overlaputils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;

public class GetAllThatHaveExonOverlapping {

	public GetAllThatHaveExonOverlapping(Collection<Alignments> exons, Collection<RefSeqGene> genes, String save)throws IOException{
		
		Map<String, IntervalTree<RefSeqGene>> geneTree=makeTree(genes);
		Map<Alignments, Collection<RefSeqGene>> overlapping=getOverlapping(exons, geneTree);
		
		write(save, overlapping);
	}
	
	
		
	private static Map<String, IntervalTree<RefSeqGene>> makeTree(Collection<RefSeqGene> genes){
		Map<String, IntervalTree<RefSeqGene>> rtrn=new TreeMap<String, IntervalTree<RefSeqGene>>();
		
		for(RefSeqGene align: genes){
			IntervalTree<RefSeqGene> tree=new IntervalTree<RefSeqGene>();
			String chr=align.getChr();
			if(rtrn.containsKey(chr)){tree=rtrn.get(chr);}
			Node node=tree.find(align.getAlignment().getStart(), align.getAlignment().getEnd());
			if(node==null){tree.put(align.getAlignment().getStart(), align.getAlignment().getEnd(), align);}
			else{
				node.incrementCount();
			}
			rtrn.put(align.getChr(), tree);
		}
		
		return rtrn;
		
	}
	
	private Map<Alignments, Collection<RefSeqGene>> getOverlapping(Collection<Alignments> exons, Map<String, IntervalTree<RefSeqGene>> geneTree){
		Map<Alignments, Collection<RefSeqGene>> rtrn=new TreeMap();
		
		for(Alignments exon: exons){
			Iterator<Node<RefSeqGene>> genes=geneTree.get(exon.getChr()).overlappers(exon.getStart(), exon.getEnd());
			Collection<RefSeqGene> overlapping=overlappingGenes(exon, genes);
			rtrn.put(exon, overlapping);
		}
		
		return rtrn;
	}
	
	private Collection<RefSeqGene> overlappingGenes(Alignments exon, Iterator<Node<RefSeqGene>> genes){
		Collection<RefSeqGene> rtrn=new TreeSet();
		
		while(genes.hasNext()){
			RefSeqGene gene=genes.next().getValue();
			if(gene.overlapsExon(exon)){rtrn.add(gene);}
		}
		
		return rtrn;
	}
	
	//For now just print all on seperate lines
	private void write(String save, Map<Alignments, Collection<RefSeqGene>> overlapping)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(Alignments exon: overlapping.keySet()){
			Collection<RefSeqGene> genes=overlapping.get(exon);
			for(RefSeqGene gene: genes){writer.write(gene+"\n");}
		}
		
		writer.close();
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
			Collection<Alignments> exons=BEDFileParser.loadAlignmentData(new File(args[0]));
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[1]));
			String save=args[2];
			new GetAllThatHaveExonOverlapping(exons, genes, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=exons (BED File) \n args[1]=genes (ESTs, Full BED) \n args[2]=save file";
	
}
