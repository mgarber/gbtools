package broad.pda.rnaseq.misc;

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
import broad.core.util.CollapseByIntersection;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;

public class CountIsoformsPerGene {

	public CountIsoformsPerGene(Collection<RefSeqGene> genes, String save)throws IOException{
		Collection<Alignments> regions=CollapseByIntersection.CollapseGenesByIntersection(genes, false);
		Map<String, IntervalTree<RefSeqGene>> geneTree=CollapseByIntersection.makeIntervalTreeForGenes(genes);
		Collection<Alignments> requiredExons=new TreeSet();
		
		Map<Alignments, Integer> counts=new TreeMap<Alignments, Integer>();
		for(Alignments region: regions){
			Iterator<Node<RefSeqGene>> geneOverlappers=geneTree.get(region.getChr()).overlappers(region.getStart(), region.getEnd());
			//This tree has all genes that overlap a region
			//go through and count
			//count only if on same strand
			int count=count(geneOverlappers);
			counts.put(region, count);
		}
		write(save, counts);
	}
	
	private int count(Iterator<Node<RefSeqGene>> geneOverlappers) {
		int counter=0;
		
		while(geneOverlappers.hasNext()){
			RefSeqGene gene=geneOverlappers.next().getValue();
			counter++;
		}
		
		return counter;
	}

	private void write(String save, Map<Alignments, Integer> counts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Alignments region: counts.keySet()){
			writer.write(region+"\t"+counts.get(region)+"\n");
		}
		
		writer.close();
	}

	private boolean isRequired(Alignments exon,	Iterator<Node<RefSeqGene>> geneOverlappers) {
		
		while(geneOverlappers.hasNext()){
			RefSeqGene gene=geneOverlappers.next().getValue();
			if(!gene.hasExon(exon)){return false;}
		}
		return true;
	}

	private void write(String save, Collection<Alignments> requiredExons) throws IOException {
		BEDFileParser.writeBED(save, requiredExons);
		
	}

	public static void main(String[] args)throws IOException{
		if(args.length>1){
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[0]));
			String save=args[1];
			new CountIsoformsPerGene(genes, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=genes \n args[1]=save";
	
}
