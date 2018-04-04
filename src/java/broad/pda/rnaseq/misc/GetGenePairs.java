package broad.pda.rnaseq.misc;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.util.CollapseByIntersection;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;

public class GetGenePairs {

	public GetGenePairs(Collection<Alignments> runOn, Collection<RefSeqGene> genes, String save) throws IOException{
		Map<String, IntervalTree<RefSeqGene>> trees=CollapseByIntersection.makeIntervalTreeForGenes(genes);
		FileWriter writer=new FileWriter(save);
		
		for(Alignments pair: runOn){
			Iterator<Node<RefSeqGene>> over=trees.get(pair.getChr()).overlappers(pair.getStart(), pair.getEnd());
			String pairString=pair.toUCSC();
			int counter=0;
			while(over.hasNext()){
				RefSeqGene gene=over.next().getValue();
				//if(gene.getOrientation().equalsIgnoreCase(pair.getOrientation())){
					pairString+="\t"+gene.getName();
					counter++;
				//}
			}
			if(counter<=2){writer.write(pairString+"\n");}
		}
		writer.close();
	}
	
	public static void main(String[] args) throws IOException{
		Collection<Alignments> runOn=BEDFileParser.loadAlignmentData(new File(args[0]));
		Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[1]));
		String save=args[2];
		new GetGenePairs(runOn, genes, save);
	}
	
}
