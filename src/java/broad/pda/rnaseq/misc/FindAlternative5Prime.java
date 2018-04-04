package broad.pda.rnaseq.misc;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.util.CollapseByIntersection;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;

public class FindAlternative5Prime {

	public FindAlternative5Prime(Collection<RefSeqGene> genes, Collection<RefSeqGene> predictions, String save)throws IOException{
		Map<RefSeqGene, Integer> alternativeStart=findAlternativeStart(genes, predictions);
		BEDFileParser.writeFullBED(save, alternativeStart);
	}
		
	
	private Map<RefSeqGene, Integer> findAlternativeStart(Collection<RefSeqGene> genes, Collection<RefSeqGene> predictions) {
		Map<RefSeqGene, Integer> rtrn=new TreeMap();
		
		Map<String, IntervalTree<RefSeqGene>> geneTree=CollapseByIntersection.makeIntervalTreeForGenes(genes);
		
		for(RefSeqGene prediction: predictions){
			//if overlaps but does so on the other strand
			Iterator<Node<RefSeqGene>> overlappers=geneTree.get(prediction.getChr()).overlappers(prediction.getStart(), prediction.getEnd());
			int distance=alternativeStart(prediction, overlappers);
			if(distance>0 && overlapsOnStrand(prediction, geneTree)<=1){rtrn.put(prediction, distance);}
		}
		
		return rtrn;
	}


	private int overlapsOnStrand(RefSeqGene prediction,	Map<String, IntervalTree<RefSeqGene>> geneTree) {
		
		Iterator<Node<RefSeqGene>> iter=geneTree.get(prediction.getChr()).overlappers(prediction.getStart(), prediction.getEnd());
		
		int counter=0;
		
		while(iter.hasNext()){
			RefSeqGene gene=iter.next().getValue();
			//if(gene.getOrientation().equalsIgnoreCase(prediction.getOrientation())){counter++;}
			counter++;
		}
		
		return counter;
	}


	//if >0 then true
	private int alternativeStart(RefSeqGene prediction,	Iterator<Node<RefSeqGene>> overlappers) {
		//If the prediction contains an exon before the annotation starts or ends
		
		IntervalTree<Alignments> exonTree=makeExonTree(overlappers, prediction.getChr());
		
		Alignments region=collapse(exonTree, prediction.getOrientation());
		System.err.println(region+" "+prediction.getAlignment().toUCSC());
		
		if(region==null){return -99;}
		
		for(Alignments exon: prediction.getExonSet()){
			//if exon not in gene (doesnt overlap)
			if(!exonTree.overlappers(exon.getStart(), exon.getEnd()).hasNext()){
				//and the exon is outside the interval
				if(exon.getStart()<region.getStart() && prediction.getOrientation().equalsIgnoreCase("+")){return region.getStart()-exon.getStart();}
				if(exon.getEnd()>region.getEnd() && prediction.getOrientation().equalsIgnoreCase("-")){return exon.getEnd()-region.getEnd();}
			}
			
		}
		
		return -99;
	}


	private Alignments collapse(IntervalTree<Alignments> tree, String orientation) {
		if(tree==null){return null;}
		String chr="";
		int start=Integer.MAX_VALUE;
		int end=-Integer.MAX_VALUE;
		
		int i=0;
		Iterator<Alignments> iter=tree.valueIterator();
		
		while(iter.hasNext()){
			Alignments gene=iter.next();
				if(gene.getOrientation().equalsIgnoreCase(orientation)){
					chr=gene.getChr();
					start=Math.min(start, gene.getStart());
					end=Math.max(end, gene.getEnd());
					i++;
				}
		}
		
		if(i>0){return new Alignments(chr, start, end);}
		return null;
	}


	private IntervalTree<Alignments> makeExonTree(Iterator<Node<RefSeqGene>> iter, String chr) {
		Collection<Alignments> exons=new TreeSet();
		
		while(iter.hasNext()){
			RefSeqGene gene=iter.next().getValue();
			exons.addAll(gene.getExonSet());
		}
		
		Map<String, IntervalTree<Alignments>> exonTree=CollapseByIntersection.makeIntervalTree(exons);
		return exonTree.get(chr);
	}


	public static void main(String[] args)throws IOException{
		if(args.length>2){
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[0]));
			Collection<RefSeqGene> predictions=BEDFileParser.loadData(new File(args[1]));
			String save=args[2];
			new FindAlternative5Prime(genes, predictions, save);
		}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=genes \n args[1]=predictions \n args[2]=save";
}
