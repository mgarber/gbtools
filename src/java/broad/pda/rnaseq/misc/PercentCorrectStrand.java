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
import broad.pda.gene.RefSeqGene;

public class PercentCorrectStrand {

	public PercentCorrectStrand(Collection<RefSeqGene> genes, Collection<RefSeqGene> predictions, String save){
		Map<String, IntervalTree<RefSeqGene>> geneTree=CollapseByIntersection.makeIntervalTreeForGenes(predictions);
		Map<String, IntervalTree<RefSeqGene>> knownGeneTree=CollapseByIntersection.makeIntervalTreeForGenes(predictions);
		Collection<RefSeqGene> correctStrand=new TreeSet();
		Collection<RefSeqGene> incorrect=new TreeSet();
		for(RefSeqGene gene: genes){
			try{
			Iterator<Node<RefSeqGene>> overlappers=geneTree.get(gene.getChr()).overlappers(gene.getStart(), gene.getEnd());
			int known=knownGeneTree.get(gene.getChr()).numOverlappers(gene.getStart(), gene.getEnd());
			if(!overlappers.hasNext() || known>1){}
			else if(sameStrand(overlappers, gene.getOrientation())){correctStrand.add(gene);}
			else{incorrect.add(gene);}
			}catch(NullPointerException ex){System.err.println(gene.getAlignment().toUCSC());}
		}
		System.err.println("correct "+correctStrand.size()+" incorrect "+incorrect.size());
	}
	
	private boolean sameStrand(Iterator<Node<RefSeqGene>> overlappers, String orientation) {
		while(overlappers.hasNext()){
			RefSeqGene gene=overlappers.next().getValue();
			if(gene.getOrientation().equalsIgnoreCase(orientation)){return true;}
		}
		return false;
	}

	public static void main(String[] args)throws IOException{
		if(args.length>2){
		Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[0]));
		Collection<RefSeqGene> predictions=BEDFileParser.loadData(new File(args[1]));
		String save=args[2];
		new PercentCorrectStrand(genes, predictions, save);
		}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=genes \n args[1]=predictions \n args[2]=save";
	
}
