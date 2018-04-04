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

public class Truncated5Prime {

	public Truncated5Prime(Collection<RefSeqGene> genes, Collection<RefSeqGene> predictions, String save) throws IOException{
		Map<String, IntervalTree<RefSeqGene>> predictionTree=CollapseByIntersection.makeIntervalTreeForGenes(predictions);
		Map<String, IntervalTree<Alignments>> exonTree=CollapseByIntersection.makeIntervalTreeForGeneExons(genes);
		
		Map<RefSeqGene, Collection<RefSeqGene>> truncatedEnds=new TreeMap();
		
		int counter=0;
		for(RefSeqGene gene: genes){
			try{
			Iterator<Node<RefSeqGene>> overlappers=predictionTree.get(gene.getChr()).overlappers(gene.getStart(), gene.getEnd());
			//check each overlapper and see what the predictions are
			Collection<RefSeqGene> truncated=getTruncated5Prime(overlappers, gene, exonTree);
			if(truncated!=null && !truncated.isEmpty()){counter++;}
			truncatedEnds.put(gene, truncated);
			}catch(NullPointerException ex){}
		}
		
		System.err.println("Number of unique truncated genes "+counter);
		
		write(save, truncatedEnds);
	}
	
	private Collection<RefSeqGene> getTruncated5Prime(Iterator<Node<RefSeqGene>> overlappers, RefSeqGene gene, Map<String, IntervalTree<Alignments>> exonTree) {
		Collection<RefSeqGene> rtrn=new TreeSet<RefSeqGene>();
		
		while(overlappers.hasNext()){
			RefSeqGene prediction=overlappers.next().getValue();
			if(prediction.getOrientation().equalsIgnoreCase(gene.getOrientation())){
				//check if prediction is 3' complete
				if(!prediction.get5PrimeExon().overlapsAtAll(gene.get5PrimeExon())){
					int compareTo=prediction.get5PrimeExon().compareTo(gene.get5PrimeExon());
					//System.out.println(gene.get3PrimeExon());
					boolean truncated=false;
					if(gene.getOrientation().equalsIgnoreCase("+") && compareTo>0){truncated=true;}
					else if(gene.getOrientation().equalsIgnoreCase("-") && compareTo<0){truncated=true;}
					
					//if not check if ends early
					if(truncated){
						Alignments lastExon=prediction.get5PrimeExon();
						//make sure doesnt overlap known exon
						boolean overlapsKnown=exonTree.get(prediction.getChr()).overlappers(lastExon.getStart(), lastExon.getEnd()).hasNext();
						if(!overlapsKnown){
							boolean is5PrimeComplete=is3PrimeComplete(prediction, gene);
							if(is5PrimeComplete){rtrn.add(prediction); System.out.println(prediction.get5PrimeExon());}
						}
					}
				}
			}
		}
		
		return rtrn;
	}
	
	//if prediction overlaps or extends past the 5' exon of the gene
	private boolean is3PrimeComplete(RefSeqGene prediction, RefSeqGene gene) {
		//if not on same strand cant be the same
		if(!prediction.getOrientation().equalsIgnoreCase(gene.getOrientation())){return false;}
		
		Alignments prediction5Prime=prediction.get3PrimeExon();
		Alignments known5Prime=gene.get3PrimeExon();
		
		if(prediction5Prime.overlaps(known5Prime)){return true;}
		//else if(prediction.getOrientation().equalsIgnoreCase("+") && prediction5Prime.getStart()<=known5Prime.getStart()){return true;}
		//else if(prediction.getOrientation().equalsIgnoreCase("-")&& prediction5Prime.getEnd()>=known5Prime.getEnd()){return true;}
	
		return false;
	}

	private void write(String save, Map<RefSeqGene, Collection<RefSeqGene>> truncatedEnds) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(RefSeqGene gene: truncatedEnds.keySet()){
			Collection<RefSeqGene> truncated=truncatedEnds.get(gene);
			for(RefSeqGene fragment: truncated){
				writer.write(fragment+"\n");
			}
		}
		
		writer.close();
	}

	public static void main(String[] args) throws IOException{
		if(args.length>2){
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[0]));
			Collection<RefSeqGene> predictions=BEDFileParser.loadData(new File(args[1]));
			String save=args[2];
			new Truncated5Prime(genes, predictions, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=genes \n args[1]=predictions \n args[2]=save";
	
}
