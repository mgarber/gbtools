package broad.pda.rnaseq.misc;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.math.Statistics;
import broad.core.util.CollapseByIntersection;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;


public class AntisenseOverlappingSense {

	public AntisenseOverlappingSense(Collection<RefSeqGene> antisense, Collection<RefSeqGene> genes, String save){
		
		Map<RefSeqGene, RefSeqGene> antisenseSenseMapping=getSense(antisense, genes);
		
		//compute genomic overlap
		double[] genomic=genomicOverlap(antisenseSenseMapping);
		
		//compute transcript overlap
		double[] transcript=transcriptOverlap(antisenseSenseMapping);
		
		System.err.println("Genomic overlap: "+Statistics.average(genomic)+" "+Statistics.stdev(genomic));
		System.err.println("Transcript overlap: "+Statistics.average(transcript)+" "+Statistics.stdev(transcript));
		
	}
	
	
	private double[] genomicOverlap(Map<RefSeqGene, RefSeqGene> antisenseSenseMapping) {
		double[] rtrn=new double[antisenseSenseMapping.size()];
		int i=0;
		for(RefSeqGene antisense: antisenseSenseMapping.keySet()){
			Alignments genomicRegion=antisenseSenseMapping.get(antisense).getAlignment();
			Collection<Alignments> exons=antisense.getExonSet();
			int sum=0;
			for(Alignments exon: exons){
				if(exon.overlapsAtAll(genomicRegion)){sum+=exon.getSize();}
			}
			rtrn[i++]=(double)sum/antisense.getTranscriptLength();
		}
		return rtrn;
	}
	
	private double[] transcriptOverlap(Map<RefSeqGene, RefSeqGene> antisenseSenseMapping) {
		double[] rtrn=new double[antisenseSenseMapping.size()];
		int i=0;
		for(RefSeqGene antisense: antisenseSenseMapping.keySet()){
			IntervalTree<Alignments> tree=antisenseSenseMapping.get(antisense).getExonTree();
			Collection<Alignments> exons=antisense.getExonSet();
			int sum=0;
			for(Alignments exon: exons){
				boolean hasNext=tree.overlappers(exon.getStart(), exon.getEnd()).hasNext();
				if(hasNext){sum+=exon.getSize();}
			}
			rtrn[i++]=(double)sum/antisense.getTranscriptLength();
		}
		return rtrn;
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
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
		Collection<RefSeqGene> antisense=BEDFileParser.loadData(new File(args[0]));
		Collection<RefSeqGene> allGenes=BEDFileParser.loadData(new File(args[1]));
		String save=args[2];
		new AntisenseOverlappingSense(antisense, allGenes, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=antisense genes \n args[1]=all genes \n args[2]=save";
	
}
