package broad.pda.rnaseq.misc;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import org.broad.igv.sam.Alignment;

import broad.core.datastructures.IntervalTree;
import broad.core.util.CollapseByIntersection;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;
import net.sf.samtools.util.CloseableIterator;
public class PairedEndsAcrossGenes {

	public PairedEndsAcrossGenes(AlignmentDataModel data, Collection<Alignments> genes) throws IOException{
		Map<String, IntervalTree<Alignments>> geneTree=CollapseByIntersection.makeIntervalTree(genes);
		double percent=spanningGenes(data, geneTree);
		System.err.println("Percent "+percent);
	}
	
	private double spanningGenes(AlignmentDataModel data, Map<String, IntervalTree<Alignments>> geneTree) throws IOException{
		double sum=0;
		double counter=0;
		
		for(String chr: data.getChromosomeLengths().keySet()){
			System.err.println(chr);
			CloseableIterator<Alignment> iter=data.getAlignmentsOverlappingRegion(new Alignments(chr, 0, data.getChromosomeLengths().get(chr)));
			while(iter.hasNext()){
				Alignment pe=iter.next();
				if(spansGenes(pe, geneTree)){sum+=1;}
				counter+=1;
			}
			iter.close();
		}
		
		return sum/counter;
	}
	
	private boolean spansGenes(Alignment pe, Map<String, IntervalTree<Alignments>> trees){
		return (trees.get(pe.getChromosome()).numOverlappers(pe.getAlignmentStart(), pe.getAlignmentEnd())>1);
	}
	
	
	public static void main(String[] args) throws IOException{
		AlignmentDataModel data=new GenericAlignmentDataModel(args[0], args[1]);
		Collection<Alignments> genes=BEDFileParser.loadAlignmentData(new File(args[2]));
		new PairedEndsAcrossGenes(data, genes);
	}
	
}
