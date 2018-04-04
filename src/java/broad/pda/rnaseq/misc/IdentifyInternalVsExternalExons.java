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

public class IdentifyInternalVsExternalExons {

	
	public IdentifyInternalVsExternalExons(Collection<Alignments> exons, Collection<RefSeqGene> genes, String save) throws IOException {
		
		Map<String, IntervalTree<RefSeqGene>> geneTree=CollapseByIntersection.makeIntervalTreeForGenes(genes);
		Collection<Alignments> internalExons=new TreeSet();
		Collection<Alignments> externalExons=new TreeSet();
		for(Alignments exon: exons){
			Iterator<Node<RefSeqGene>> geneOverlappers=geneTree.get(exon.getChr()).overlappers(exon.getStart(), exon.getEnd());
			boolean internal=isInternal(exon, geneOverlappers);
			if(internal){internalExons.add(exon);}
			else{externalExons.add(exon);}
		}
		write(save+".internal", internalExons);
		write(save+".external", externalExons);
	}

	private boolean isInternal(Alignments exon,	Iterator<Node<RefSeqGene>> geneOverlappers) {
		//if exon overlaps either the first or last exon in a gene then return false
		//else return true
		while(geneOverlappers.hasNext()){
			RefSeqGene gene=geneOverlappers.next().getValue();
			if(gene.getFirstExon().overlaps(exon) || gene.getLastExon().overlaps(exon)){return false;}
		}
		return true;
	}

	private void write(String save, Collection<Alignments> externalExons) throws IOException {
		BEDFileParser.writeBED(save, externalExons);
		
	}

	public static void main(String[] args)throws IOException{
		if(args.length>2){
			Collection<Alignments> exons=BEDFileParser.loadAlignmentData(new File(args[0]));
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[1]));
			String save=args[2];
			new IdentifyInternalVsExternalExons(exons, genes, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" Goal: Identify internal vs external exons \n args[0]=Additional exons \n args[1]=genes \n args[2]=save";
	
}
