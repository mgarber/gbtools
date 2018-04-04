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

public class CountIsoformsPerGeneAcrossTypes {

	public CountIsoformsPerGeneAcrossTypes(Collection<Alignments> regions, Collection<RefSeqGene> genes1, Collection<RefSeqGene> genes2, String save)throws IOException{
		Map<String, IntervalTree<RefSeqGene>> geneTree1=CollapseByIntersection.makeIntervalTreeForGenes(genes1);
		Map<String, IntervalTree<RefSeqGene>> geneTree2=CollapseByIntersection.makeIntervalTreeForGenes(genes2);
	
		Collection<Alignments> requiredExons=new TreeSet();
		
		Map<Alignments, Integer> counts=new TreeMap<Alignments, Integer>();
		for(Alignments region: regions){
			try{
			Iterator<Node<RefSeqGene>> geneOverlappers1=geneTree1.get(region.getChr()).overlappers(region.getStart(), region.getEnd());
			Iterator<Node<RefSeqGene>> geneOverlappers2=geneTree2.get(region.getChr()).overlappers(region.getStart(), region.getEnd());
			
			//This tree has all genes that overlap a region
			//go through and count
			//count only if on same strand
			int count1=count(geneOverlappers1);
			int count2=count(geneOverlappers2);
			counts.put(region, count1);
			System.out.println(region+"\t"+count1+"\t"+count2);
			}catch(NullPointerException ex){}
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
			Collection<Alignments> regions=BEDFileParser.loadAlignmentData(new File(args[0]));
			Collection<RefSeqGene> genes1=BEDFileParser.loadData(new File(args[1]));
			Collection<RefSeqGene> genes2=BEDFileParser.loadData(new File(args[2]));
			String save=args[3];
			new CountIsoformsPerGeneAcrossTypes(regions, genes1, genes2, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=genes \n args[1]=save";
	
}
