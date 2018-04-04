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

public class FindShorterStartOverlappingK4 {

	public FindShorterStartOverlappingK4(Collection<RefSeqGene> genes, Collection<RefSeqGene> predictions, Collection<Alignments> k4, String save)throws IOException{
		Collection<RefSeqGene> fullGenes=getFullGenes(genes, predictions);
		Collection<RefSeqGene> overlappingK4=getTruncatedOverlappingK4(predictions, fullGenes, k4, genes);
		write(save, overlappingK4);
	}
	
	
	private Collection<RefSeqGene> getTruncatedOverlappingK4(Collection<RefSeqGene> predictions,Collection<RefSeqGene> fullGenes, Collection<Alignments> k4, Collection<RefSeqGene> genes) {
		Map<String, IntervalTree<RefSeqGene>> geneTree=CollapseByIntersection.makeIntervalTreeForGenes(genes);
		Collection<RefSeqGene> rtrn=new TreeSet();
		
		Map<String, IntervalTree<Alignments>> k4Tree=CollapseByIntersection.makeIntervalTree(k4);
		
		for(RefSeqGene prediction: predictions){
			if(!fullGenes.contains(prediction) && prediction.getNumExons()>1){
				Alignments fivePrime=prediction.get5PrimeExon();
				if(fivePrime!=null){
					if(k4Tree.get(fivePrime.getChr()).overlappers(fivePrime.getStart(), fivePrime.getEnd()).hasNext() && geneTree.get(fivePrime.getChr()).overlappers(fivePrime.getStart(), fivePrime.getEnd()).hasNext()){rtrn.add(prediction);}
				}
			}
		}
		
		return rtrn;
	}


	private void write(String save, Collection<RefSeqGene> overlappingK4) throws IOException {
		BEDFileParser.writeFullBED(save, overlappingK4);
		
	}


	private Collection<RefSeqGene> getFullGenes(Collection<RefSeqGene> genes, Collection<RefSeqGene> predictions) {
		Collection<RefSeqGene> rtrn=new TreeSet();
		
		Map<String, IntervalTree<RefSeqGene>> geneTree=CollapseByIntersection.makeIntervalTreeForGenes(genes);
		
		for(RefSeqGene prediction: predictions){
			//if overlaps but does so on the other strand
			Iterator<Node<RefSeqGene>> overlappers=geneTree.get(prediction.getChr()).overlappers(prediction.getStart(), prediction.getEnd());
			boolean b=isComplete(prediction, overlappers);
			if(b){rtrn.add(prediction);}
		}
		
		return rtrn;
	}
	

	private boolean isComplete(RefSeqGene prediction, Iterator<Node<RefSeqGene>> overlappers) {
		//test if our prediction matches the annotation well
		while(overlappers.hasNext()){
			RefSeqGene gene=overlappers.next().getValue();
			//if 5' end of prediciton overlaps or extends past the 5' end of gene return true
			boolean prime5=is5PrimeComplete(prediction, gene);
			//if 3' end of prediction overlaps or extends past the 5' end of gene return true
			//boolean prime3=is3PrimeComplete(prediction, gene);
			if(prime5){return true;}
		}
		return false;
	}

	private boolean is3PrimeComplete(RefSeqGene prediction, RefSeqGene gene) {
		//if not on same strand cant be the same
		if(!prediction.getOrientation().equalsIgnoreCase(gene.getOrientation())){return false;}
		
		Alignments prediction3Prime=prediction.get3PrimeExon();
		Alignments known3Prime=gene.get3PrimeExon();
		
		if(prediction3Prime.overlaps(known3Prime)){return true;}
		else if(prediction.getOrientation().equalsIgnoreCase("+") && prediction3Prime.getEnd()>=known3Prime.getEnd()){return true;}
		else if(prediction.getOrientation().equalsIgnoreCase("-")&& prediction3Prime.getStart()<=known3Prime.getStart()){return true;}
	
		return false;
		
	}

	
	//if prediction overlaps or extends past the 5' exon of the gene
	private boolean is5PrimeComplete(RefSeqGene prediction, RefSeqGene gene) {
		//if not on same strand cant be the same
		if(!prediction.getOrientation().equalsIgnoreCase(gene.getOrientation())){return false;}
		
		Alignments prediction5Prime=prediction.get5PrimeExon();
		Alignments known5Prime=gene.get5PrimeExon();
		
		if(prediction5Prime.overlaps(known5Prime)){return true;}
		else if(prediction.getOrientation().equalsIgnoreCase("+") && prediction5Prime.getStart()<=known5Prime.getStart()){return true;}
		else if(prediction.getOrientation().equalsIgnoreCase("-")&& prediction5Prime.getEnd()>=known5Prime.getEnd()){return true;}
	
		return false;
	}
	
	private boolean isComplete(RefSeqGene prediction, RefSeqGene gene){
		Alignments predictionFirstExon=prediction.getFirstExon();
		Alignments predictionLastExon=prediction.getLastExon();
		
		Alignments knownFirstExon=gene.getFirstExon();
		Alignments knownLastExon=gene.getLastExon();
		
		boolean first=false;
		boolean last=false;
		
		if(predictionFirstExon.overlaps(knownFirstExon)){first=true;}
		if(predictionLastExon.overlaps(knownLastExon)){last=true;}
		if(predictionFirstExon.getStart()<=knownFirstExon.getStart()){first=true;}
		if(predictionLastExon.getEnd()>=knownLastExon.getEnd()){last=true;}
				
		return (first && last);
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>3){
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[0]));
			Collection<RefSeqGene> predictions=BEDFileParser.loadData(new File(args[1]));
			Collection<Alignments> k4=BEDFileParser.loadAlignmentData(new File(args[2]));
			String save=args[3];
			new FindShorterStartOverlappingK4(genes, predictions, k4, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=genes \n args[1]=predictions \n args[2]=k4 regions \n args[3]=save";
	
}
