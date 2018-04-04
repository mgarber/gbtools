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

public class ClassifyAllPredictions {

	public ClassifyAllPredictions(Collection<RefSeqGene> genes, Collection<RefSeqGene> predictions, Collection<Alignments> k4, String save) throws IOException{
		//Make interval tree for genes
		Map<String, IntervalTree<RefSeqGene>> predictionTree=CollapseByIntersection.makeIntervalTreeForGenes(predictions);
		Map<String, IntervalTree<RefSeqGene>> geneTree=CollapseByIntersection.makeIntervalTreeForGenes(genes);
		
		FileWriter writer=new FileWriter(save);
		
		//We want to iterate through all genes and figure out what predictions overlap it
		for(RefSeqGene gene: genes){
			//System.err.println(gene);
			
			//get overlapping gene
			if(predictionTree.containsKey(gene.getChr())){
				Iterator<Node<RefSeqGene>> overlappingPredictions=predictionTree.get(gene.getChr()).overlappers(gene.getStart(), gene.getEnd());
				
				//go through the overlappers, and classify the prediction
				String group=classify(gene, overlappingPredictions, geneTree.get(gene.getChr())); //TODO Add gene tree
				
				writer.write(gene+"\t"+group+"\n");
			}
			
		}
		
		
		writer.close();
	}
		
	private String classify(RefSeqGene gene, Iterator<Node<RefSeqGene>> overlappingPredictions, IntervalTree<RefSeqGene> geneTree) {
		String rtrn="";
		Collection<RefSeqGene> predictions=getPredictions(overlappingPredictions, gene.getOrientation());
		
		//case -1: if no overlapping predictions --> not constructed (not expressed)
		boolean isConstructed=isConstructed(predictions, gene);
		
		if(!isConstructed){return "Not constructed";}
		
		//case 0: if all overlapping predictions are single exon or unstranded
		boolean isSingleExon=isSingleExon(predictions, gene);
		if(isSingleExon){return "Single exon";}
		
		//case 1: Overlaps predictions on same strand but they ends before the gene --> truncated 3' end
		boolean truncated3Prime=isTruncated3Prime(predictions, gene);
		
		//case 2: Overlaps predictions on same strand but they ends after the gene --> extended 3' end
		boolean extended3Prime=isExtended3Prime(predictions, gene);
		
		//case 3: Overlaps prediction but they start before the gene --> extended 5' start
		boolean extended5Prime=isExtended5Prime(predictions, gene);
		
		//case 4: Overlaps prediction but they start after the gene --> truncated 5' start
		boolean truncated5Prime=isTruncated5Prime(predictions, gene);
		
		//case 5: Overlaps prediction but they have extra internal exons --> novel internal exons
		//boolean hasNovelInternal=hasNovelInternal(predictions, gene);
		
		//case 6: Overlaps prediction but they have extra external exons --> novel external exons
		//boolean hasNovelExternal=hasNovelExternal(predictions, gene);
		
		//case 7: Overlaps prediction and those predictions overlap another gene on same strand as well --> Run-on
		//boolean isRunOn=isRunOn(predictions, gene, geneTree);
		
		if(truncated3Prime){rtrn+="Truncated3Prime,";}
		if(extended3Prime){rtrn+="Extended3Prime,";}
		if(extended5Prime){rtrn+="Extended5Prime,";}
		if(truncated5Prime){rtrn+="Truncated5Prime,";}
		//if(hasNovelInternal){rtrn+="Novel Internal,";}
		//if(hasNovelExternal){rtrn+="Novel External,";}
		//if(isRunOn){rtrn+="Run on";}
		
		
		//if(truncated3Prime){System.out.println(gene);}
		
		return rtrn;
	}
	
	
	private boolean isTruncated5Prime(Collection<RefSeqGene> predictions, RefSeqGene gene) {
		//is extended 5' if no end goes to the end of the transcript
		Collection<Alignments> exons=new TreeSet();
		for(RefSeqGene prediction: predictions){exons.addAll(prediction.getExonSet());}
		
		
		Alignments gene5PrimeExon=gene.get5PrimeExon();
		/*boolean complete=false;
		for(Alignments exon: exons){
			if(exon.overlapsAtAll(gene3PrimeExon) || gene3PrimeExon.overlapsAtAll(exon)){complete=true;}
		}*/
		
		//if last exon overlaps the last
		
		if(gene.getOrientation().equalsIgnoreCase("+")){
			Alignments firstExon=(Alignments)exons.toArray()[0];
			int order=firstExon.compareTo(gene5PrimeExon);
			if(order>0){return true;}
		}
		else{
			Alignments firstExon=(Alignments)exons.toArray()[exons.size()-1];
			int order=firstExon.compareTo(gene5PrimeExon);
			if(order<0){return true;}
		}
		
		
		return false;
	}

	private boolean isExtended5Prime(Collection<RefSeqGene> predictions, RefSeqGene gene) {
		//is extended 5' if no end goes to the end of the transcript
		Collection<Alignments> exons=new TreeSet();
		for(RefSeqGene prediction: predictions){exons.addAll(prediction.getExonSet());}
		
		
		Alignments gene5PrimeExon=gene.get5PrimeExon();
		/*boolean complete=false;
		for(Alignments exon: exons){
			if(exon.overlapsAtAll(gene3PrimeExon) || gene3PrimeExon.overlapsAtAll(exon)){complete=true;}
		}*/
		
		//if last exon overlaps the last
		
		if(gene.getOrientation().equalsIgnoreCase("+")){
			Alignments firstExon=(Alignments)exons.toArray()[0];
			int order=firstExon.compareTo(gene5PrimeExon);
			if(order<0){return true;}
		}
		else{
			Alignments firstExon=(Alignments)exons.toArray()[exons.size()-1];
			int order=firstExon.compareTo(gene5PrimeExon);
			if(order>0){return true;}
		}
		
		
		return false;
	}

	private boolean isExtended3Prime(Collection<RefSeqGene> predictions,RefSeqGene gene) {
		//is truncated 3' if no end goes to the end of the transcript
		Collection<Alignments> exons=new TreeSet();
		for(RefSeqGene prediction: predictions){exons.addAll(prediction.getExonSet());}
		
		
		Alignments gene3PrimeExon=gene.get3PrimeExon();
		/*boolean complete=false;
		for(Alignments exon: exons){
			if(exon.overlapsAtAll(gene3PrimeExon) || gene3PrimeExon.overlapsAtAll(exon)){complete=true;}
		}*/
		
		//if last exon overlaps the last
		
		if(gene.getOrientation().equalsIgnoreCase("+")){
			Alignments lastExon=(Alignments)exons.toArray()[exons.size()-1];
			int order=lastExon.compareTo(gene3PrimeExon);
			if(order>0){return true;}
		}
		else{
			Alignments firstExon=(Alignments)exons.toArray()[0];
			int order=firstExon.compareTo(gene3PrimeExon);
			if(order<0){return true;}
		}
		
		
		return false;
	}

	private boolean isSingleExon(Collection<RefSeqGene> predictions, RefSeqGene gene) {
		int numExons=1;
		for(RefSeqGene prediction: predictions){
			numExons=Math.max(prediction.getNumExons(), numExons);
		}
		if(numExons>1){return false;}
		return true;
	}

	private boolean isConstructed(Collection<RefSeqGene> predictions, RefSeqGene gene) {
		if(predictions!=null && !predictions.isEmpty()){return true;}
		return false;
	}

	private boolean isTruncated3Prime(Collection<RefSeqGene> predictions, RefSeqGene gene) {
		//is truncated 3' if no end goes to the end of the transcript
		Collection<Alignments> exons=new TreeSet();
		for(RefSeqGene prediction: predictions){exons.addAll(prediction.getExonSet());}
		
		
		Alignments gene3PrimeExon=gene.get3PrimeExon();
		boolean complete=false;
		for(Alignments exon: exons){
			if(exon.overlapsAtAll(gene3PrimeExon) || gene3PrimeExon.overlapsAtAll(exon)){complete=true;}
		}
		
		//if last exon overlaps the last
		if(!complete){
			if(gene.getOrientation().equalsIgnoreCase("+")){
				Alignments lastExon=(Alignments)exons.toArray()[exons.size()-1];
				int order=lastExon.compareTo(gene3PrimeExon);
				if(order<0){return true;}
			}
			else{
				Alignments firstExon=(Alignments)exons.toArray()[0];
				int order=firstExon.compareTo(gene3PrimeExon);
				if(order>0){return true;}
			}
		}
		
		return false;
	}

	private Collection<RefSeqGene> getPredictions(Iterator<Node<RefSeqGene>> overlappingPredictions, String orientation){
		Collection<RefSeqGene> rtrn=new TreeSet();
		
		while(overlappingPredictions.hasNext()){
			RefSeqGene prediction=overlappingPredictions.next().getValue();
			if(prediction.getOrientation().equalsIgnoreCase(orientation)){rtrn.add(prediction);}
		}
		
		return rtrn;
	}

	public static void main(String[] args)throws IOException{
		if(args.length>3){
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[0]));
			Collection<RefSeqGene> predictions=BEDFileParser.loadData(new File(args[1]));
			Collection<Alignments> k4=BEDFileParser.loadAlignmentData(new File(args[2]));
			String save=args[3];
			new ClassifyAllPredictions(genes, predictions, k4, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=genes \n args[1]=predictions \n args[2]=k4 regions \n args[3]=save";
}
