package broad.pda.genetools;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.gene.RefSeqGeneWithIsoforms;

public class GetOverlappingGenes {

	public GetOverlappingGenes(Map<String, IntervalTree<RefSeqGeneWithIsoforms>> genesToFilter, Map<String, Collection<Alignments>> genesToFilterBy, String save)throws IOException{
		
		Collection<RefSeqGeneWithIsoforms> filtered=new ArrayList<RefSeqGeneWithIsoforms>();
		
		for(String chr: genesToFilterBy.keySet()){
			System.err.println(chr);
			IntervalTree<RefSeqGeneWithIsoforms> lincs=genesToFilter.get(chr);
			Collection<Alignments> regions=genesToFilterBy.get(chr);
			for(Alignments region: regions){
				Iterator<Node<RefSeqGeneWithIsoforms>> iter=lincs.overlappers(region.getStart(), region.getEnd());
				addAll(filtered, iter);
			}			
		}
		
		write(save, filtered);
	}
	
	
	private void write(String save, Collection<RefSeqGeneWithIsoforms> filtered) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(RefSeqGeneWithIsoforms gene: filtered){
			Collection<RefSeqGene> isoforms=gene.getAllIsoforms();
			for(RefSeqGene isoform: isoforms){
				writer.write(isoform+"\n");
			}
		}
		
		writer.close();
	}


	private void addAll(Collection<RefSeqGeneWithIsoforms> filtered, Iterator<Node<RefSeqGeneWithIsoforms>> iter) {
		while(iter.hasNext()){
			RefSeqGeneWithIsoforms gene=iter.next().getValue();
			filtered.add(gene);
		}
	}


	private Collection<Alignments> filter(Collection<Alignments> lincs, IntervalTree<Alignments> proteins, boolean exclude){
		Collection<Alignments> rtrn=new ArrayList();
		
		for(Alignments linc: lincs){
			boolean overlaps=proteins.overlappers(linc.getStart(), linc.getEnd()).hasNext();
			if(exclude && !overlaps){rtrn.add(linc);}
			else if(!exclude && overlaps){rtrn.add(linc);}
		}
		
		return rtrn;
	}
	
	
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
			System.err.println("started");
			BEDFileParser bed=new BEDFileParser(args[0]);
			Map<String, IntervalTree<RefSeqGeneWithIsoforms>> genesToFilter=bed.getIntervalTreeWithIsoforoms();
			System.err.println("loaded lincRNAs");
			Map<String, Collection<Alignments>> genesToFilterBy=BEDFileParser.loadAlignmentDataByChr(new File(args[1]));
			System.err.println("loaded proteins into interval tree");
			String save=args[2];
			new GetOverlappingGenes(genesToFilter, genesToFilterBy, save);
			
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=genes to filter (Full BED- lincRNA transcripts) \n args[1]=genes to filter by (K4K36/protein coding genes) \n args[2]=save";
	
}
