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

public class GetPairedKnownAndNovel5PrimeAndOverlapWithK4 {

	public GetPairedKnownAndNovel5PrimeAndOverlapWithK4(Collection<RefSeqGene> novel5Prime, Collection<RefSeqGene> genes, Collection<Alignments> K4, Collection<Alignments> K4K27, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		Map<String, IntervalTree<RefSeqGene>> trees=CollapseByIntersection.makeIntervalTreeForGenes(genes);
		
		Map<Alignments, Collection<Alignments>> novelKnownMap=new TreeMap();
		Collection<Alignments> all5Prime=new TreeSet();
		
		//get novel 5 prime and known 5 prime
		for(RefSeqGene novel: novel5Prime){
			Iterator<Node<RefSeqGene>> overlappers=trees.get(novel.getChr()).overlappers(novel.getStart(), novel.getEnd());
			Collection<Alignments> known5Primes=getKnownFirstExons(overlappers);
			Alignments novel5PrimeExon=novel.get5PrimeExon();
			//System.out.println(novel5PrimeExon+"\tNovel");
			//for(Alignments t: known5Primes){System.out.println(t+"\tKnown");}
			novelKnownMap.put(novel5PrimeExon, known5Primes);
			all5Prime.addAll(known5Primes);
			//all5Prime.add(novel5PrimeExon);
		}
		
		//get all 5' ends not used
		Map<String, IntervalTree<Alignments>> novelTree=CollapseByIntersection.makeIntervalTreeForGeneExons(novel5Prime);
		Collection<Alignments> notExpressedKnown=new TreeSet();
		for(Alignments known: all5Prime){
			Iterator<Node<Alignments>> novel=novelTree.get(known.getChr()).overlappers(known.getStart(), known.getEnd());
			boolean overlapsNovel=novel.hasNext();
			if(!overlapsNovel){notExpressedKnown.add(known);}
		}
		
		Map<String, IntervalTree<Alignments>> k4Tree=CollapseByIntersection.makeIntervalTree(K4);
		Map<String, IntervalTree<Alignments>> bivalentTree=CollapseByIntersection.makeIntervalTree(K4K27);
		
		//make sure the unused 5' ends dont share a k4 with a used one
		for(Alignments known: notExpressedKnown){
			Iterator<Node<Alignments>> overlappingK4=k4Tree.get(known.getChr()).overlappers(known.getStart(), known.getEnd());
			Iterator<Node<Alignments>> overlappingK27=bivalentTree.get(known.getChr()).overlappers(known.getStart(), known.getEnd());
			int k4Only=count(overlappingK4, novelTree);
			int bivalent=count(overlappingK27, novelTree);
			writer.write(known+"\t"+k4Only+"\t"+bivalent+"\n");
		}
		writer.close();
		
		//get expression
		//Map<Alignments, double[]> allExpressionLevels=data.scoreSegments(all5Prime);
		
		
		//write expression levels
		//writeExpressionLevels(save, novelKnownMap, allExpressionLevels);
	}
	
	private int count(Iterator<Node<Alignments>> overlappingK4,	Map<String, IntervalTree<Alignments>> novelTree) {
		int count=0;
		while(overlappingK4.hasNext()){
			Alignments region=overlappingK4.next().getValue();
			Iterator<Node<Alignments>> over=novelTree.get(region.getChr()).overlappers(region.getStart(), region.getEnd());
			boolean hasNovel=over.hasNext();
			if(!hasNovel){count++;}
		}
		return count;
	}

	private void writeExpressionLevels(String save,	Map<Alignments, Collection<Alignments>> novelKnownMap, Map<Alignments, double[]> allExpressionLevels) throws IOException {
		//write expression of novel only
		FileWriter writer=new FileWriter(save+".novel.expression");
		for(Alignments novel: novelKnownMap.keySet()){
			double[] scores=allExpressionLevels.get(novel);
			writer.write(novel+"\t"+scores[0]+"\t"+scores[1]+"\n");
		}
		writer.close();
		
		
		//write expression of novel and known only
		writer=new FileWriter(save+".novelAndKnown.expression");
	
		for(Alignments novel: novelKnownMap.keySet()){
			Collection<Alignments> knownSet=novelKnownMap.get(novel);
			for(Alignments known: knownSet){
				double[] scoresNovel=allExpressionLevels.get(novel);
				double[] scoresKnown=allExpressionLevels.get(known);
				writer.write(novel+"\t"+known+"\t"+scoresNovel[0]+"\t"+scoresKnown[0]+"\t"+scoresNovel[1]+"\t"+scoresKnown[1]+"\n");
			}
		}
		
		writer.close();
	}

	private Collection<Alignments> getKnownFirstExons(Iterator<Node<RefSeqGene>> overlappers) {
		Collection<Alignments> rtrn=new TreeSet();
		
		while(overlappers.hasNext()){
			RefSeqGene gene=overlappers.next().getValue();
			rtrn.add(gene.get5PrimeExon());
		}
		
		return rtrn;
	}

	public static void main(String[] args)throws IOException{
		if(args.length>4){
			Collection<RefSeqGene> novel5Prime=BEDFileParser.loadData(new File(args[0]));
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[1]));
			String save=args[2];
			Collection<Alignments> k4=BEDFileParser.loadAlignmentData(new File(args[3]));
			Collection<Alignments> k27=BEDFileParser.loadAlignmentData(new File(args[4]));
			new GetPairedKnownAndNovel5PrimeAndOverlapWithK4(novel5Prime, genes, k4, k27, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=Novel 5' ends \n args[1]=genes \n args[2]=save \n args[3]=k4 only \n args[4]=k4-k27";
}
