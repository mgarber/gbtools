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

public class LongerThanAnnotation {

	public LongerThanAnnotation(Collection<RefSeqGene> annotation, Collection<RefSeqGene> predictions, String save)throws IOException{
		Map<String, IntervalTree<RefSeqGene>> geneTree=CollapseByIntersection.makeIntervalTreeForGenes(annotation);
	
		//Collection<RefSeqGene> bigger=getBigger(geneTree, predictions);
		
		Map<RefSeqGene, Collection<Alignments>> extraExons=getExtraExons(geneTree, predictions);
		
		write(save, extraExons);
	}
	
	private Map<RefSeqGene, Collection<Alignments>> getExtraExons(Map<String, IntervalTree<RefSeqGene>> geneTree, Collection<RefSeqGene> predictions) {
		Map<RefSeqGene, Collection<Alignments>> rtrn=new TreeMap();
		
		int i=0;
		for(RefSeqGene prediction: predictions){
			Iterator<Node<RefSeqGene>> iter=geneTree.get(prediction.getChr()).overlappers(prediction.getStart(), prediction.getEnd());
			Collection<Alignments> extraExons=getExtraExons(prediction, iter);
			if(extraExons!=null && !extraExons.isEmpty()){
				rtrn.put(prediction, extraExons);
			}
			i++;
			if(i% 10000 ==0){System.err.println(i);}
		}
		
		return rtrn;
	}

	private Collection<Alignments> getExtraExons(RefSeqGene prediction, Iterator<Node<RefSeqGene>> iter) {
		if(prediction.getNumExons()==1){return null;}
		Collection<Alignments> exons=new TreeSet();
		Collection<Alignments> rtrn=new TreeSet();
		
		while(iter.hasNext()){
			RefSeqGene gene=iter.next().getValue();
			if(gene.getOrientation().equalsIgnoreCase(prediction.getOrientation())){
				exons.addAll(gene.getExonSet());
			}
		}
		
		if(exons.isEmpty()){return null;}
		
		Map<String, IntervalTree<Alignments>> exonTree=CollapseByIntersection.makeIntervalTree(exons);
		
		if(exonTree==null || exonTree.get(prediction.getChr())==null){return null;}
		
		for(Alignments exon: prediction.getExonSet()){
			if(!exonTree.get(exon.getChr()).overlappers(exon.getStart(), exon.getEnd()).hasNext()){rtrn.add(exon);}
		}
		return rtrn;
	}

	//loop through predictions and get overlapping gene
	//check if bigger than gene on left or right
	//if so report
	private Collection<RefSeqGene> getBigger(Map<String, IntervalTree<RefSeqGene>> geneTree, Collection<RefSeqGene> predictions){
		Collection<RefSeqGene> rtrn=new TreeSet();
		
		int i=0;
		for(RefSeqGene prediction: predictions){
			Iterator<Node<RefSeqGene>> iter=geneTree.get(prediction.getChr()).overlappers(prediction.getStart(), prediction.getEnd());
			if(isBigger(prediction, iter)){rtrn.add(prediction);}
			i++;
			if(i% 10000 ==0){System.err.println(i);}
		}
		
		return rtrn;
	}

	private void write(String save, Collection<RefSeqGene> genes)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(RefSeqGene gene: genes){writer.write(gene+"\n");}
		
		writer.close();
	}
	
	private void write(String save, Map<RefSeqGene, Collection<Alignments>> genes)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		Collection<Alignments> temp=new TreeSet();
		for(RefSeqGene gene: genes.keySet()){
			//temp.addAll(genes.get(gene));
			writer.write(gene+"\n");
		}
		
		
			/*for(Alignments exon: temp){
				writer.write(exon+"\n");
			}*/
		
		
		writer.close();
	}
	
	private boolean isBigger(RefSeqGene prediction, Iterator<Node<RefSeqGene>> iter) {
		
		while(iter.hasNext()){
			RefSeqGene gene=iter.next().getValue();
			if(prediction.getStart()>=gene.getStart() && prediction.getEnd()<=gene.getEnd()){return false;}
		}
		
		return true;
	}

	public static void main(String[] args)throws IOException{
		if(args.length>2){
		Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[0]));
		Collection<RefSeqGene> predictions=BEDFileParser.loadData(new File(args[1]));
		String save=args[2];
		new LongerThanAnnotation(genes, predictions, save);
		}
		else{System.err.println(usage);}
	}

	static String usage=" args[0]=genes \n args[1]=predictions \n args[2]=save";
}
