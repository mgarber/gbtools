package broad.pda.seq.alignment;

import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.util.CollapseByIntersection;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.graph.ChromosomeWithBubbles2;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

public class MakeJunctionsOffConnectivityGraph {

	int readLength=76;
	
	public MakeJunctionsOffConnectivityGraph(String samFile, String genomeDir, String save) throws IOException{
		AlignmentDataModel data=new GenericAlignmentDataModel(samFile, genomeDir+"/sizes");
		Map<String, Collection<RefSeqGene>> genesByChr=new TreeMap();
		Map<String, Collection<RefSeqGene>> junctionsByChr=new TreeMap();
		
		for(String chr: data.getChromosomeLengths().keySet()){
			System.err.println(chr);
			Collection<RefSeqGene> genes=getGenes(chr, data, genomeDir);
			genesByChr.put(chr, genes);
		}
		
		new MakeSpliceJunctionMaps(genesByChr, junctionsByChr, save, genomeDir);
	}
	
	private Collection<RefSeqGene> getGenes(String chr, AlignmentDataModel data, String genomeDir) throws IOException{
		Collection<RefSeqGene> rtrn=new TreeSet<RefSeqGene>();
		
		String seqFile=genomeDir+"/"+chr.replaceAll("chr", "")+"/"+chr+".fa";
		FastaSequenceIO fsio = new FastaSequenceIO(seqFile);
		List<Sequence> seqs = fsio.loadAll();
		Sequence chrSequence = seqs.get(0);
		
		Alignments chrRegion=new Alignments(chr, 0, data.getChromosomeLengths().get(chr));
		
		Map<Alignments, Integer> introns=data.getSplicedReads(chrRegion, 0, 10000000, 0, chrSequence);
		System.err.println("Got introns");
		Collection<Alignments> exons=data.getOverlappingRegions(chr);
		
		System.err.println("Making 2nd graph");
		
		ChromosomeWithBubbles2 graph=constructGraphs(chr, exons, introns.keySet(),0);
		
		System.err.println("Finished 2nd Graph construction");
		
		System.err.println(graph.getPaths(0).size());
		graph.writeGraph("test.graph");
		Collection<RefSeqGene> allPaths=graph.getGenePaths(0);
		
		return allPaths;
	}
	
	
	private ChromosomeWithBubbles2 constructGraphs(String chr, Collection<Alignments> exons, Collection<Alignments> introns, int minDistance) {
		//Add a collapse step up front
		exons=extendByReadLength(exons);
		exons=CollapseByIntersection.collapseByIntersection(exons, false);
		
		Collection<Alignments> decollapsed=CollapseByIntersection.DecollapseByIntronLocation(exons, introns);
		//for(Alignments exon: introns){System.out.println(exon);}
		
		System.err.println(decollapsed.size()+" "+introns.size());
		ChromosomeWithBubbles2 graph=new ChromosomeWithBubbles2(chr, decollapsed, introns, null, null, 0,0,0,0);
		
		return graph;
	}
	
	private Collection<Alignments> extendByReadLength(Collection<Alignments> exons) {
		Collection<Alignments> rtrn=new TreeSet<Alignments>();
		
		for(Alignments exon: exons){
			Alignments extended=new Alignments(exon.getChr(), exon.getStart()-readLength, exon.getEnd()+readLength);
			rtrn.add(extended);
		}
		
		return rtrn;
	}

	public static void main(String[] args)throws IOException{
		if(args.length>2){
		String samFile=args[0];
		String genomeDir=args[1];
		String save=args[2];
		new MakeJunctionsOffConnectivityGraph(samFile, genomeDir, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=sam file \n args[1]=genomeDir \n args[2]=save";
	
}
