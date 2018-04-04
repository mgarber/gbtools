package broad.pda.seq.pairedend;

import java.util.Collection;
import java.util.TreeSet;

import broad.pda.gene.RefSeqGene;
import broad.pda.seq.graph.ChromosomeWithBubbles2;
import broad.pda.seq.graph.Path;

public class FillInWithPairedEnds {

	public static Collection<RefSeqGene> FillInWithPairs(ChromosomeWithBubbles2 graph){
		Collection<RefSeqGene> rtrn=new TreeSet();
		
		Collection<Path> paths=graph.getAllPaths();
		for(Path path: paths){
			RefSeqGene gene=path.toFilledInGene();
			rtrn.add(gene);
		}
		
		return rtrn;
	}
	
}
