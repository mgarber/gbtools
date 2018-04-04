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
import broad.pda.gene.RefSeqGene;

public class FindOverlappingAntisense {
	
	int minNum=2;
	
	public FindOverlappingAntisense(Collection<RefSeqGene> genes, Collection<RefSeqGene> predictions, String save, int minNum)throws IOException{
		this.minNum=minNum;
		//Goal: Find transcripts that overlap proteins but are on the opposite strand
		Collection<RefSeqGene> antisense=getOverlappingAntisense(genes, predictions);
		//Collection<RefSeqGene> overlappingGenes=getOverlappingGenes(genes, predictions);
		write(save, antisense);
		//write(save+".overlapping.novel", overlappingGenes);
	}
	
	
	private Collection<RefSeqGene> getOverlappingGenes(Collection<RefSeqGene> genes, Collection<RefSeqGene> predictions) {
		Collection<RefSeqGene> rtrn=new TreeSet();
		Map<String, IntervalTree<RefSeqGene>> geneTree=CollapseByIntersection.makeIntervalTreeForGenes(genes);
		for(RefSeqGene prediction: predictions){
			//if overlaps but does so on the other strand
			Iterator<Node<RefSeqGene>> overlappers=geneTree.get(prediction.getChr()).overlappers(prediction.getStart(), prediction.getEnd());
			boolean b=overlappingGenes(prediction, overlappers);
			if(b){rtrn.add(prediction);}
		}
		
		return rtrn;
	}


	private boolean overlappingGenes(RefSeqGene prediction, Iterator<Node<RefSeqGene>> overlappers) {
		if(prediction.getNumExons()==1){return false;}
		
		boolean hasSense=false;
		boolean hasAntisense=false;
		
		while(overlappers.hasNext()){
			RefSeqGene gene=overlappers.next().getValue();
			if(!prediction.getOrientation().equalsIgnoreCase("*") && !gene.getOrientation().equalsIgnoreCase("*") && !gene.getOrientation().equalsIgnoreCase(prediction.getOrientation())){
				if(gene.getOrientation().equalsIgnoreCase("+")){hasSense=true;}
				if(gene.getOrientation().equalsIgnoreCase("-")){hasAntisense=true;}
			}
		}
		
		return (hasSense && hasAntisense);
	}


	private Collection<RefSeqGene> getOverlappingAntisense(Collection<RefSeqGene> genes, Collection<RefSeqGene> predictions) {
		Collection<RefSeqGene> rtrn=new TreeSet();
		
		Map<String, IntervalTree<RefSeqGene>> geneTree=CollapseByIntersection.makeIntervalTreeForGenes(genes);
		
		for(RefSeqGene prediction: predictions){
			//if overlaps but does so on the other strand
			Iterator<Node<RefSeqGene>> overlappers=geneTree.get(prediction.getChr()).overlappers(prediction.getStart(), prediction.getEnd());
			boolean b=novelAntisense(prediction, overlappers);
			if(b){rtrn.add(prediction);}
		}
		
		return rtrn;
	}



	//Gets genes antisense to protein coding genes
	private boolean novelAntisense(RefSeqGene prediction,Iterator<Node<RefSeqGene>> overlappers) {
		if(prediction.getNumExons()<=minNum){return false;}
		
		boolean antisense=false;
		boolean sense=false;
		
		while(overlappers.hasNext()){
			RefSeqGene gene=overlappers.next().getValue();
			if(!prediction.getOrientation().equalsIgnoreCase("*") && !gene.getOrientation().equalsIgnoreCase("*") && !gene.getOrientation().equalsIgnoreCase(prediction.getOrientation())){
				antisense=true;
			}
			else if(!prediction.getOrientation().equalsIgnoreCase("*") && !gene.getOrientation().equalsIgnoreCase("*") && gene.getOrientation().equalsIgnoreCase(prediction.getOrientation())){
				sense=true;
			}
		}
		
		if(!antisense && !sense){return false;}
		else if(antisense && sense){return false;}
		else if(sense && !antisense){return false;}
		
		
		return true;
	}




	private void write(String save, Collection<RefSeqGene> antisense) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(RefSeqGene gene: antisense){
			writer.write(gene+"\n");
		}
		
		writer.close();
	}

	public static void main(String[] args)throws IOException{
		if(args.length>3){
		Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[0]));
		Collection<RefSeqGene> predictions=BEDFileParser.loadData(new File(args[1]));
		String save=args[2];
		int minNum=new Integer(args[3]);
		new FindOverlappingAntisense(genes, predictions, save, minNum);
		}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=genes \n args[1]=predictions \n args[2]=save \n args[3]=min num of exons";
}
