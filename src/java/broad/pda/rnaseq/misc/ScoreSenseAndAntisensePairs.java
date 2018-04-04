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
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

public class ScoreSenseAndAntisensePairs {

	public ScoreSenseAndAntisensePairs(Collection<RefSeqGene> antisense, Collection<RefSeqGene> allGenes, ContinuousDataAlignmentModel data, String save) throws IOException{
		Map<RefSeqGene, RefSeqGene> antisenseSenseMapping=getSense(antisense, allGenes);
		
		Collection<RefSeqGene> allSenseAntisense=getAll(antisenseSenseMapping);
		Map<RefSeqGene, double[]> scores=data.scoreGenesStranded(allSenseAntisense);
		
		write(save, antisenseSenseMapping, scores);
	}
	
	
	private Map<RefSeqGene, RefSeqGene> getSense(Collection<RefSeqGene> antisense, Collection<RefSeqGene> allGenes) {
		Map<String, IntervalTree<RefSeqGene>> trees=CollapseByIntersection.makeIntervalTreeForGenes(allGenes);
		Map<RefSeqGene, RefSeqGene> rtrn=new TreeMap<RefSeqGene, RefSeqGene>();
		
		for(RefSeqGene gene: antisense){
			Iterator<Node<RefSeqGene>> genes=trees.get(gene.getChr()).overlappers(gene.getStart(), gene.getEnd());
			while(genes.hasNext()){
				RefSeqGene sense=genes.next().getValue();
				if(!sense.getOrientation().equalsIgnoreCase(gene.getOrientation())){rtrn.put(gene, sense);}
			}
		}
		
		return rtrn;
	}


	private Collection<RefSeqGene> getAll(Map<RefSeqGene, RefSeqGene> antisenseSenseMapping) {
		Collection<RefSeqGene> rtrn=new TreeSet<RefSeqGene>();
		
		for(RefSeqGene r1: antisenseSenseMapping.keySet()){
			RefSeqGene r2=antisenseSenseMapping.get(r1);
			rtrn.add(r1);
			rtrn.add(r2);
		}
		
		return rtrn;
	}


	private void write(String save, Map<RefSeqGene, RefSeqGene> antisenseSenseMapping, Map<RefSeqGene, double[]> scores) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(RefSeqGene antisense: antisenseSenseMapping.keySet()){
			RefSeqGene sense=antisenseSenseMapping.get(antisense);
			if(scores.containsKey(antisense) && scores.containsKey(sense)){
				writer.write(antisense.getName()+"\t"+sense.getName()+"\t"+scores.get(antisense)[1]+"\t"+scores.get(sense)[1]+"\t"+scores.get(antisense)[0]+"\t"+scores.get(sense)[0]+"\n");
			}
		}
		
		writer.close();
	}


	public static void main(String[] args)throws IOException{
		if(args.length>4){
		Collection<RefSeqGene> antisense=BEDFileParser.loadData(new File(args[0]));
		Collection<RefSeqGene> allGenes=BEDFileParser.loadData(new File(args[1]));
		ContinuousDataAlignmentModel data=new ContinuousDataAlignmentModel(new GenericAlignmentDataModel(args[2], args[3]));
		String save=args[4];
		new ScoreSenseAndAntisensePairs(antisense, allGenes, data, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=antisense genes \n args[1]=all genes \n args[2]=alignments \n args[3]=sizes \n args[4]=save";
	
}
