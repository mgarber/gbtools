package broad.pda.rnaseq.misc;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeSet;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.util.CollapseByIntersection;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;

public class ExonsAlwaysInIsoform {

	//Goal: If an isoform exists that doesnt have the exon then return false else return true
	public ExonsAlwaysInIsoform(Collection<Alignments> exons, Collection<RefSeqGene> genes, String save)throws IOException{
		Map<String, IntervalTree<RefSeqGene>> geneTree=CollapseByIntersection.makeIntervalTreeForGenes(genes);
		Collection<Alignments> requiredExons=new TreeSet();
		for(Alignments exon: exons){
			Iterator<Node<RefSeqGene>> geneOverlappers=geneTree.get(exon.getChr()).overlappers(exon.getStart(), exon.getEnd());
			boolean required=isRequired(exon, geneOverlappers);
			if(required){requiredExons.add(exon);}
		}
		write(save, requiredExons);
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
		if(args.length>2){
			Collection<Alignments> exons=BEDFileParser.loadAlignmentData(new File(args[0]));
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[1]));
			String save=args[2];
			new ExonsAlwaysInIsoform(exons, genes, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" Goal: Get Required exons (ones in all isoforms) \n args[0]=Additional exons \n args[1]=genes \n args[2]=save";
	
}
