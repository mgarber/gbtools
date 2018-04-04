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

public class Alternative3PrimeEnd {

	//Goal: Find transcripts that extend past the annotated 3' end of the given gene
	public Alternative3PrimeEnd(Collection<RefSeqGene> genes, Collection<RefSeqGene> predictions, String save) throws IOException{
		Map<RefSeqGene, Integer> alternative3Prime=getLonger3Prime(genes, predictions);
		write(save, alternative3Prime);
	}
	
	private Map<RefSeqGene, Integer> getLonger3Prime(Collection<RefSeqGene> genes, Collection<RefSeqGene> predictions) {
		Map<RefSeqGene, Integer> rtrn=new TreeMap();
		
		//split genes by strand
		//collapse into units
		//make interval trees of units
		Map<String, IntervalTree<Alignments>>[] transcriptionalUnits=makeTUs(genes);
		Map<String, IntervalTree<RefSeqGene>>[] geneTrees=makeGeneTrees(genes);
		
		for(RefSeqGene prediction: predictions){
			Map<String, IntervalTree<RefSeqGene>> tree=geneTrees[0];
			if(prediction.getOrientation().equalsIgnoreCase("-")){tree=geneTrees[1];}
			Iterator<Node<RefSeqGene>> overlappers=tree.get(prediction.getChr()).overlappers(prediction.getStart(), prediction.getEnd());
			int distance=longer3Prime(prediction, overlappers);
			if(distance>0 && overlaps(prediction, transcriptionalUnits)<=1){rtrn.put(prediction, distance);}
		}
		
		return rtrn;
	}

	private int overlaps(RefSeqGene prediction,	Map<String, IntervalTree<Alignments>>[] transcriptionalUnits) {
		Map<String, IntervalTree<Alignments>> tree=transcriptionalUnits[0];
		if(prediction.getOrientation().equalsIgnoreCase("-")){tree=transcriptionalUnits[1];}
		
		return tree.get(prediction.getChr()).numOverlappers(prediction.getStart(), prediction.getEnd());
		
	}

	private int longer3Prime(RefSeqGene prediction,	Iterator<Node<RefSeqGene>> overlappers) {
		//boolean rtrn=false;
		
		if(prediction.getNumExons()==1 || prediction.getOrientation().equalsIgnoreCase("*")){return -99;}
		Collection<Alignments> exons=new TreeSet();
		
		if(!overlappers.hasNext()){return -99;}
		
		while(overlappers.hasNext()){
			RefSeqGene gene=overlappers.next().getValue();
			exons.addAll(gene.getExonSet());
		}
		
		Alignments end=(Alignments)exons.toArray()[exons.size()-1];
		if(prediction.getOrientation().equalsIgnoreCase("-")){end=(Alignments)exons.toArray()[0];}
		
		
		if(!prediction.get3PrimeExon().overlaps(end)){
		if(prediction.getOrientation().equalsIgnoreCase("+") && prediction.get3PrimeExon().getEnd()>end.getEnd()){return prediction.get3PrimeExon().getEnd()-end.getEnd();}
		if(prediction.getOrientation().equalsIgnoreCase("-") && prediction.get3PrimeExon().getStart()<end.getStart()){return end.getStart()-prediction.getStart();}
		}
		
		
		return -99;
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
	
	private void write(String save, Map<RefSeqGene, Integer> alternative3Prime)throws IOException{
		BEDFileParser.writeFullBED(save, alternative3Prime);
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
		Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[0]));
		Collection<RefSeqGene> predictions=BEDFileParser.loadData(new File(args[1]));
		String save=args[2];
		new Alternative3PrimeEnd(genes, predictions, save);
		}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=genes \n args[1]=predictions \n args[2]=save";
	
}
