package broad.pda.rnaseq.misc;

import java.io.File;
import java.io.FileWriter;
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

public class IsoformsWithEnds {

	public IsoformsWithEnds(Collection<Alignments> alternative5Prime, Collection<RefSeqGene> predictions, Collection<RefSeqGene> known, String save) throws IOException{
		Map<String, IntervalTree<RefSeqGene>> predictionTree=CollapseByIntersection.makeIntervalTreeForGenes(predictions);
		Map<String, IntervalTree<RefSeqGene>> geneTree=CollapseByIntersection.makeIntervalTreeForGenes(known);
		
		Collection<Alignments> both=new TreeSet();
		Collection<Alignments> onlyNovel=new TreeSet();
		
		//get all predictions overlapping alternative 5 prime
		for(Alignments alternative: alternative5Prime){
			Iterator<Node<RefSeqGene>> predictionOverlappers=predictionTree.get(alternative.getChr()).overlappers(alternative.getStart(), alternative.getEnd());
			Iterator<Node<RefSeqGene>> geneOverlappers=geneTree.get(alternative.getChr()).overlappers(alternative.getStart(), alternative.getEnd());
			//check if the predictions have the annotated start as well
			boolean annotated=hasAnnotated(predictionOverlappers, geneOverlappers);
			if(annotated){both.add(alternative);}
			else{onlyNovel.add(alternative);}
		}
		
		write(save+".both.bed", both);
		write(save+".novel.bed", onlyNovel);
		
	}
	
	private void write(String save, Collection<Alignments> both) throws IOException {
		FileWriter writer=new FileWriter(save);
	
		for(Alignments r: both){
			writer.write(r+"\n");
		}
		
		writer.close();
	}

	private boolean hasAnnotated(Iterator<Node<RefSeqGene>> predictionOverlappers, Iterator<Node<RefSeqGene>> geneOverlappers) {
		Collection<RefSeqGene> predictions=new TreeSet();
		Collection<RefSeqGene> genes=new TreeSet();
		while(predictionOverlappers.hasNext()){
			RefSeqGene prediction=predictionOverlappers.next().getValue();
			predictions.add(prediction);
		}
		while(geneOverlappers.hasNext()){
			RefSeqGene prediction=geneOverlappers.next().getValue();
			genes.add(prediction);
		}
		
		for(RefSeqGene prediction: predictions){
			for(RefSeqGene gene: genes){
				if(prediction.getFirstExon().overlapsAtAll(gene.getFirstExon())){return true;}
			}
		}
		
		return false;
	}

	public static void main(String[] args)throws IOException{
		if(args.length>3){
		Collection<Alignments> alternative5Prime=BEDFileParser.loadAlignmentData(new File(args[0]));
		Collection<RefSeqGene> predictions=BEDFileParser.loadData(new File(args[1]));
		Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[2]));
		String save=args[3];
		new IsoformsWithEnds(alternative5Prime, predictions, genes, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=alternative 5' ends (unique regions) \n args[1]=predictions \n args[2]=genes \n args[3]=save";
}
