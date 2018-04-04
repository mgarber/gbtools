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

public class NovelExonsInGene {

	public NovelExonsInGene(Collection<RefSeqGene> genes, Collection<RefSeqGene> predictions, String save) throws IOException{
		Map<String, IntervalTree<RefSeqGene>> predictionTree=CollapseByIntersection.makeIntervalTreeForGenes(predictions);
		Map<String, IntervalTree<Alignments>> exonTree=CollapseByIntersection.makeIntervalTreeForGeneExons(genes);
		Collection<Alignments> novelExons=new TreeSet();
		
		for(RefSeqGene gene: genes){
			try{
			Iterator<Node<RefSeqGene>> overlappers=predictionTree.get(gene.getChr()).overlappers(gene.getStart(), gene.getEnd());
			Collection<Alignments> exons=getNovelInternal(gene, overlappers);
			novelExons.addAll(exons);
			}catch(NullPointerException ex){}
		}
		
		write(save, novelExons, exonTree);
	}
	
	
	private Collection<Alignments> getNovelInternal(RefSeqGene gene, Iterator<Node<RefSeqGene>> overlappers) {
		Collection<Alignments> exons=new TreeSet();
		Collection<Alignments> rtrn=new TreeSet();
		IntervalTree<Alignments> geneTree=new IntervalTree();
		IntervalTree<Alignments> endNodes=new IntervalTree();
		
		Collection<Alignments> t=gene.getSortedAndUniqueExons();
		for(Alignments exon: t){geneTree.put(exon.getStart(), exon.getEnd(), exon);}
		
		Alignments first=gene.getFirstExon();
		Alignments last=gene.getLastExon();
		
		while(overlappers.hasNext()){
			RefSeqGene prediction=overlappers.next().getValue();
			Alignments firstNode=prediction.getFirstExon();
			Alignments lastNode=prediction.getLastExon();
			endNodes.put(firstNode.getStart(), firstNode.getEnd(), firstNode);
			endNodes.put(lastNode.getStart(), lastNode.getEnd(), lastNode);
			Collection<Alignments> predictionExons=prediction.getExonSet();
			for(Alignments exon : predictionExons){
				boolean novel=!geneTree.overlappers(exon.getStart(), exon.getEnd()).hasNext();
				if(novel){
					//check if internal
					if(exon.compareTo(first)>0 && exon.compareTo(last)<0){rtrn.add(exon);}
				}
			}
		}
		
		Collection<Alignments> temp=new TreeSet();
		for(Alignments exon: rtrn){
			if(!endNodes.overlappers(exon.getStart(), exon.getEnd()).hasNext()){temp.add(exon);}
		}
		
		return temp;
	}


	private void write(String save, Collection<Alignments> novelExons, Map<String, IntervalTree<Alignments>> exonTree) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Alignments exon: novelExons){
			boolean novel=!exonTree.get(exon.getChr()).overlappers(exon.getStart(), exon.getEnd()).hasNext();
			if(novel){writer.write(exon+"\n");}
		}
		
		writer.close();
	}


	public static void main(String[] args)throws IOException{
		if(args.length>2){
		Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[0]));
		Collection<RefSeqGene> predictions=BEDFileParser.loadData(new File(args[1]));
		String save=args[2];
		new NovelExonsInGene(genes, predictions, save);
		}
		else{System.err.println(usage);}
	}

	static String usage=" args[0]=genes \n args[1]=predictions \n args[2]=save";
}
