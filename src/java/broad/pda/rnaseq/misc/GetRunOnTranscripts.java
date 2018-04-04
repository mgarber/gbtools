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

public class GetRunOnTranscripts {
	int n=1;

	public GetRunOnTranscripts(Collection<RefSeqGene> genes, Collection<RefSeqGene> predictions, String save)throws IOException{
		//find predictions that span multiple genes on same strand
		Collection<RefSeqGene> runon=findRunOns(genes, predictions);
		
		
		write(save, runon);
	}
	
	private Collection<RefSeqGene> findRunOns(Collection<RefSeqGene> genes, Collection<RefSeqGene> predictions) {
		Collection<RefSeqGene> rtrn=new TreeSet();
		
		//split genes by strand
		//collapse into units
		//make interval trees of units
		Map<String, IntervalTree<Alignments>>[] transcriptionalUnits=makeTUs(genes);
		Map<String, IntervalTree<RefSeqGene>>[] geneTrees=makeGeneTrees(genes);
		
		//get overlapping units on strand of prediction
		for(RefSeqGene prediction: predictions){
			if(prediction.getNumExons()>1){
			Map<String, IntervalTree<Alignments>> tree=transcriptionalUnits[0];
			if(prediction.getOrientation().equalsIgnoreCase("-")){tree=transcriptionalUnits[1];}
			int num=tree.get(prediction.getChr()).numOverlappers(prediction.getStart(), prediction.getEnd());
			//if overlaps more than 1 return
				if(num>1){
					//valid predictions should overlap exons on both transcripts
					//check if overlaps at least n exons in both genes
					Iterator<Node<Alignments>> tus=tree.get(prediction.getChr()).overlappers(prediction.getStart(), prediction.getEnd());
					boolean overlaps=overlapsExons(prediction, geneTrees, tus, n);
					if(overlaps){rtrn.add(prediction);}
				}
			}
		}
				
		return rtrn;
	}

	private Map<String, IntervalTree<RefSeqGene>>[] makeGeneTrees(Collection<RefSeqGene> genes) {
		Collection<RefSeqGene> plusGenes=new TreeSet();
		Collection<RefSeqGene> minusGenes=new TreeSet();
		
		for(RefSeqGene gene: genes){
			if(gene.getOrientation().equalsIgnoreCase("+")){plusGenes.add(gene);}
			else if(gene.getOrientation().equalsIgnoreCase("-")){minusGenes.add(gene);}
		}
		
		Map<String, IntervalTree<RefSeqGene>>[] rtrn=new Map[2];
		
		rtrn[0]=CollapseByIntersection.makeIntervalTreeForGenes(plusGenes);
		rtrn[1]=CollapseByIntersection.makeIntervalTreeForGenes(minusGenes);
		
		return rtrn;
	}

	private boolean overlapsExons(RefSeqGene prediction, Map<String, IntervalTree<RefSeqGene>>[] geneTrees, Iterator<Node<Alignments>> tus, int n2) {
		Map<String, IntervalTree<RefSeqGene>> geneTree=geneTrees[0];
		if(prediction.getOrientation().equalsIgnoreCase("-")){geneTree=geneTrees[1];}
		
		int counter=0;
		while(tus.hasNext()){
			Alignments tu=tus.next().getValue();
			IntervalTree<Alignments> exons=getExons(geneTree.get(tu.getChr()).overlappers(tu.getStart(), tu.getEnd()));
			int num=numOverlappingExons(prediction, exons);
			if(num>=n2){counter++;}
		}
		
		if(counter>1){return true;}
		return false;
	}

	private IntervalTree<Alignments> getExons(Iterator<Node<RefSeqGene>> overlappers) {
		IntervalTree<Alignments> rtrn=new IntervalTree();
		
		while(overlappers.hasNext()){
			RefSeqGene gene=overlappers.next().getValue();
			for(Alignments exon: gene.getExonSet()){rtrn.put(exon.getStart(), exon.getEnd(), exon);}
		}
		
		return rtrn;
	}

	private int numOverlappingExons(RefSeqGene prediction, IntervalTree<Alignments> exonTree) {
		int counter=0;
		for(Alignments exon: prediction.getExonSet()){
			if(exonTree.overlappers(exon.getStart(), exon.getEnd()).hasNext()){counter++;}
		}
		return counter;
	}

	private Map<String, IntervalTree<Alignments>>[] makeTUs(Collection<RefSeqGene> genes) {
		//split genes by strand
		Collection<RefSeqGene> plusGenes=split(genes, "+");
		Collection<RefSeqGene> minusGenes=split(genes, "-");
		
		//collapse into units
		Collection<Alignments> plusTUs=CollapseByIntersection.CollapseGenesByIntersection(plusGenes, false);
		Collection<Alignments> minusTUs=CollapseByIntersection.CollapseGenesByIntersection(minusGenes, false);
		
		//make interval trees of units
		Map<String, IntervalTree<Alignments>> plusTree=CollapseByIntersection.makeIntervalTree(plusTUs);
		Map<String, IntervalTree<Alignments>> minusTree=CollapseByIntersection.makeIntervalTree(minusTUs);
	
		Map<String, IntervalTree<Alignments>> [] rtrn=new Map[2];
		
		rtrn[0]=plusTree;
		rtrn[1]=minusTree;
		
		return rtrn;
	}
	

	private Collection<RefSeqGene> split(Collection<RefSeqGene> genes, String string) {
		Collection<RefSeqGene> rtrn=new TreeSet();

		for(RefSeqGene gene: genes){
			if(gene.getOrientation().equalsIgnoreCase(string)){rtrn.add(gene);}
		}
		
		return rtrn;
	}

	private void write(String save, Collection<RefSeqGene> runon) throws IOException {
		BEDFileParser.writeFullBED(save, runon);
	}

	public static void main(String[] args)throws IOException{
		if(args.length>2){
		Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[0]));
		Collection<RefSeqGene> predictions=BEDFileParser.loadData(new File(args[1]));
		String save=args[2];
		new GetRunOnTranscripts(genes, predictions, save);
		}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=genes \n args[1]=predictions \n args[2]=save";
	
}
