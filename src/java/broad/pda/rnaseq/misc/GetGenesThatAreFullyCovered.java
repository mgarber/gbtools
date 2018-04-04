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

public class GetGenesThatAreFullyCovered {

	Map<RefSeqGene, String> temp;
	Map<String, IntervalTree<RefSeqGene>> geneTree;
	Map<String, IntervalTree<RefSeqGene>> predictionTree;
	
	//Goal: Idnetify genes that are fully covered from 5'-->3' of the known annotation
	public GetGenesThatAreFullyCovered(Collection<RefSeqGene> genes, Collection<RefSeqGene> predictions, String save) throws IOException{		
		temp=new TreeMap();
		Map<RefSeqGene, Collection<RefSeqGene>> complete=getFullGenes(genes, predictions);
		double percent=getPercentFullGenes(genes, predictions);


		Collection<RefSeqGene> uniqueGenes=new TreeSet();
		for(RefSeqGene prediction: complete.keySet()){
			Collection<RefSeqGene> full=complete.get(prediction);
			uniqueGenes.addAll(full);
		}
		
	}
	public GetGenesThatAreFullyCovered(Collection<RefSeqGene> genes, Collection<RefSeqGene> predictions){
		//this.predictions= predictions;
		this.geneTree=CollapseByIntersection.makeIntervalTreeForGenes(genes);
		this.predictionTree = CollapseByIntersection.makeIntervalTreeForGenes(predictions);
	}
	
	public void scoreCompleteness() {		
		for(String chr : geneTree.keySet()) {
			Iterator<RefSeqGene> geneIt = geneTree.get(chr).valueIterator();
			while(geneIt.hasNext()) {
				RefSeqGene gene = geneIt.next();
				scoreCompleteness(gene);
			}
		}
	}
	
	public void scoreCompleteness(RefSeqGene gene) {
		IntervalTree<RefSeqGene> chrtree = predictionTree.get(gene.getChr());
		if(chrtree != null) {
			Iterator<Node<RefSeqGene>> overlappers=chrtree.overlappers(gene.getStart(), gene.getEnd());
			boolean isComplete = false;
			while(overlappers.hasNext() && !isComplete){
				RefSeqGene prediction=overlappers.next().getValue();
				isComplete = isComplete(prediction, gene);
			}
			gene.addExtraField(isComplete ? "1" : "0");
		}
		
	}
	
	public void writeGenes(String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		for(String chr : geneTree.keySet()) {
			Iterator<RefSeqGene> geneIt = geneTree.get(chr).valueIterator();
			while(geneIt.hasNext()){
				writer.write(geneIt.next()+"\n");
			}
		}
		writer.close();
	}
	private double getPercentFullGenes(Collection<RefSeqGene> genes, Collection<RefSeqGene> predictions) {
		
		
		
		Map<String, IntervalTree<RefSeqGene>> geneTree=CollapseByIntersection.makeIntervalTreeForGenes(genes);
		
		for(RefSeqGene prediction: predictions){
			//if overlaps but does so on the other strand
			IntervalTree<RefSeqGene> chrtree = geneTree.get(prediction.getChr());
			if(geneTree.get(prediction.getChr()) != null) {
				Iterator<Node<RefSeqGene>> overlappers=chrtree.overlappers(prediction.getStart(), prediction.getEnd());
				getComplete(prediction, overlappers);
			}
			
		}
		int[] array=this.sumFullPartial(temp);
		double full=array[0];
		double all=array[1];
		double percent=full/(full+all);
		System.err.println("V2 Percent: "+percent+" Full overlaps: "+full+" Any overlap: "+all);
		
		return percent;
	}

	private void getComplete(RefSeqGene prediction,Iterator<Node<RefSeqGene>> overlappers) {
		
		
		
		//test if our prediction matches the annotation well
		while(overlappers.hasNext()){
			RefSeqGene gene=overlappers.next().getValue();
			//if 5' end of prediciton overlaps or extends past the 5' end of gene return true
			//boolean prime5=is5PrimeComplete(prediction, gene);
			//if 3' end of prediction overlaps or extends past the 5' end of gene return true
			//boolean prime3=is3PrimeComplete(prediction, gene);
			if(isComplete(prediction, gene)){temp.put(gene, "full");}
			//if(prime5 && !prime3){temp.put(gene, "full");}
			else{if(!temp.containsKey(gene)){temp.put(gene,"partial");}}
		}
		
		
	}

	private int[] sumFullPartial(Map<RefSeqGene, String> temp) {
		int full=0;
		int partial=0;
		
		for(RefSeqGene gene: temp.keySet()){
			String name=temp.get(gene);
			if(name.equalsIgnoreCase("full")){full++;}
			else{partial++; System.out.println(gene);}
		}
		
		int[] rtrn={full, partial};
		return rtrn;
	}

	private Map<RefSeqGene, Collection<RefSeqGene>> getFullGenes(Collection<RefSeqGene> genes, Collection<RefSeqGene> predictions) {
		Map<RefSeqGene, Collection<RefSeqGene>> rtrn=new TreeMap();
		
		Map<String, IntervalTree<RefSeqGene>> geneTree=CollapseByIntersection.makeIntervalTreeForGenes(genes);
		
		for(RefSeqGene prediction: predictions){
			//if overlaps but does so on the other strand
			IntervalTree<RefSeqGene> chrtree = geneTree.get(prediction.getChr());
			if(chrtree != null) {
				Iterator<Node<RefSeqGene>> overlappers=chrtree.overlappers(prediction.getStart(), prediction.getEnd());
				boolean b=isComplete(prediction, overlappers);
				if(b){rtrn.put(prediction, toCollection(prediction, chrtree.overlappers(prediction.getStart(), prediction.getEnd())));}
			}
		}
		
		return rtrn;
	}
	
	
	private Collection<RefSeqGene> toCollection(RefSeqGene prediction, Iterator<Node<RefSeqGene>> overlappers) {
		Collection<RefSeqGene> rtrn=new TreeSet();
		
		while(overlappers.hasNext()){
			RefSeqGene gene=overlappers.next().getValue();
			if(gene.getOrientation().equalsIgnoreCase(prediction.getOrientation())){rtrn.add(gene);}
		}
		
		return rtrn;
	}

	private boolean isComplete(RefSeqGene prediction, Iterator<Node<RefSeqGene>> overlappers) {
		//test if our prediction matches the annotation well
		while(overlappers.hasNext()){
			RefSeqGene gene=overlappers.next().getValue();
			if(isComplete(prediction, gene)){return true;}
		}
		return false;
	}

	private boolean is3PrimeComplete(RefSeqGene prediction, RefSeqGene gene) {
		//if not on same strand cant be the same
		if(!prediction.getOrientation().equalsIgnoreCase(gene.getOrientation())){return false;}
		
		Alignments prediction3Prime=prediction.get3PrimeExon();
		Alignments known3Prime=gene.get3PrimeExon();
		
		if(prediction3Prime.overlaps(known3Prime)){return true;}
		//else if(prediction.getOrientation().equalsIgnoreCase("+") && prediction3Prime.getEnd()>=known3Prime.getEnd()){return true;}
		//else if(prediction.getOrientation().equalsIgnoreCase("-")&& prediction3Prime.getStart()<=known3Prime.getStart()){return true;}
	
		return false;
		
	}

	
	//if prediction overlaps or extends past the 5' exon of the gene
	private boolean is5PrimeComplete(RefSeqGene prediction, RefSeqGene gene) {
		//if not on same strand cant be the same
		if(!prediction.getOrientation().equalsIgnoreCase(gene.getOrientation())){return false;}
		
		Alignments prediction5Prime=prediction.get5PrimeExon();
		Alignments known5Prime=gene.get5PrimeExon();
		
		if(prediction5Prime.overlaps(known5Prime)){return true;}
		//else if(prediction.getOrientation().equalsIgnoreCase("+") && prediction5Prime.getStart()<=known5Prime.getStart()){return true;}
		//else if(prediction.getOrientation().equalsIgnoreCase("-")&& prediction5Prime.getEnd()>=known5Prime.getEnd()){return true;}
	
		return false;
	}
	
	private boolean isComplete(RefSeqGene prediction, RefSeqGene gene){
		//if 5' end of prediciton overlaps or extends past the 5' end of gene return true
		boolean prime5=is5PrimeComplete(prediction, gene);
		//if 3' end of prediction overlaps or extends past the 5' end of gene return true
		boolean prime3=is3PrimeComplete(prediction, gene);
		return (prime3 && prime5);
	}

	public static void main(String[] args)throws IOException{
		if(args.length>2){
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[0]));
			Collection<RefSeqGene> predictions=BEDFileParser.loadData(new File(args[1]));
			String save=args[2];
			//new GetGenesThatAreFullyCovered(genes, predictions, save);
			GetGenesThatAreFullyCovered worker = new GetGenesThatAreFullyCovered(genes, predictions);
			worker.scoreCompleteness();
			worker.writeGenes(save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=genes \n args[1]=predictions \n args[2]=save";
	
}
