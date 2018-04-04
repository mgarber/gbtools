package broad.pda.rnai;

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
import broad.pda.rnai.designer.RNAiFileFormatUtils;

public class GetOverlappingRefSeq {

	public GetOverlappingRefSeq(Collection<RNAiGeneAnnotation> annotations, Collection<RefSeqGene> genes, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		Map<String, IntervalTree<RefSeqGene>> trees=CollapseByIntersection.makeIntervalTreeForGenes(genes);
		
		for(RNAiGeneAnnotation rnai: annotations){
			Alignments region=rnai.getRegion();
			Iterator<Node<RefSeqGene>> overlappers=trees.get(region.getChr()).overlappers(region.getStart(), region.getEnd());
			writer.write(rnai.toRNAi()+"\t"+commaSep(overlappers)+"\n");
		}
		writer.close();
	}

	private String commaSep(Iterator<Node<RefSeqGene>> overlappers) {
		String rtrn="";
		
		while(overlappers.hasNext()){
			RefSeqGene gene=overlappers.next().getValue();
			rtrn+=gene.getName();
			if(overlappers.hasNext()){rtrn+=",";}
		}
		
		return rtrn;
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
			Collection<RNAiGeneAnnotation> annotations=RNAiFileFormatUtils.parseRNAiReportFile(new File(args[0]));
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[1]));
			String save=args[2];
			new GetOverlappingRefSeq(annotations, genes, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=RNAi Report \n args[1]=RefSeq Full BED \n args[2]=save";
}
