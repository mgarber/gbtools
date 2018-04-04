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
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

public class GetPairedKnownAndNovel3Prime {

	public GetPairedKnownAndNovel3Prime(Collection<RefSeqGene> novel5Prime, Collection<RefSeqGene> genes, ContinuousDataAlignmentModel data, String save) throws IOException{
		Map<String, IntervalTree<RefSeqGene>> trees=CollapseByIntersection.makeIntervalTreeForGenes(genes);
		
		Map<Alignments, Collection<Alignments>> novelKnownMap=new TreeMap();
		Collection<Alignments> all5Prime=new TreeSet();
		
		//get novel 5 prime and known 5 prime
		for(RefSeqGene novel: novel5Prime){
			Iterator<Node<RefSeqGene>> overlappers=trees.get(novel.getChr()).overlappers(novel.getStart(), novel.getEnd());
			Collection<Alignments> known5Primes=getKnownFirstExons(overlappers);
			Alignments novel5PrimeExon=novel.get3PrimeExon();
			System.out.println(novel5PrimeExon+"\tNovel");
			for(Alignments t: known5Primes){System.out.println(t+"\tKnown");}
			novelKnownMap.put(novel5PrimeExon, known5Primes);
			all5Prime.addAll(known5Primes);
			all5Prime.add(novel5PrimeExon);
		}
		
		//get expression
		Map<Alignments, double[]> allExpressionLevels=data.scoreSegments(all5Prime);
		
		
		//write expression levels
		writeExpressionLevels(save, novelKnownMap, allExpressionLevels);
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
			rtrn.add(gene.get3PrimeExon());
		}
		
		return rtrn;
	}

	public static void main(String[] args)throws IOException{
		if(args.length>4){
			Collection<RefSeqGene> novel5Prime=BEDFileParser.loadData(new File(args[0]));
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[1]));
			String save=args[2];
			ContinuousDataAlignmentModel data=new ContinuousDataAlignmentModel(new GenericAlignmentDataModel(args[3], args[4]));
			new GetPairedKnownAndNovel3Prime(novel5Prime, genes, data, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=Novel 5' ends \n args[1]=genes \n args[2]=save \n args[3]=alignments \n args[4]=sizes";
}
