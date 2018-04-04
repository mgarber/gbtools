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

public class ExcludeOverlappingGenes {

	public ExcludeOverlappingGenes(Map<String, IntervalTree<RefSeqGeneWithIsoforms>> genesToFilter, Map<String, IntervalTree<Alignments>> genesToFilterBy, String save)throws IOException{
		
		Collection<RefSeqGeneWithIsoforms> filtered=new ArrayList<RefSeqGeneWithIsoforms>();
		
		for(String chr: genesToFilter.keySet()){
			System.err.println(chr);
			Iterator <Node<RefSeqGeneWithIsoforms>> lincs=genesToFilter.get(chr).iterator();
			IntervalTree<Alignments> regions=genesToFilterBy.get(chr);
			while(lincs.hasNext()){
				RefSeqGeneWithIsoforms transcript=lincs.next().getValue();
				Iterator<Node<Alignments>> iter=regions.overlappers(transcript.getStart(), transcript.getEnd());
				if(!iter.hasNext()){filtered.add(transcript);}
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
			Map<String, IntervalTree<Alignments>> genesToFilterBy=BEDFileParser.loadAlignmentDataToTree(new File(args[1]));
			System.err.println("loaded proteins into interval tree");
			String save=args[2];
			new ExcludeOverlappingGenes(genesToFilter, genesToFilterBy, save);
			
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=genes to filter (Full BED- lincRNA transcripts) \n args[1]=genes to filter by (K4K36/protein coding genes) \n args[2]=save";
	
}
