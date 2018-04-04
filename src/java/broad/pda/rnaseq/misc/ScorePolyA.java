package broad.pda.rnaseq.misc;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.motif.SequenceMotif;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.core.util.CollapseByIntersection;
import broad.pda.annotation.BEDFileParser;
import broad.pda.chromosome.Chromosome;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.rnai.ExtractSequence;
import broad.pda.rnai.designer.ComputeMIRScore;

public class ScorePolyA {

	String motif="AATAAA";
	
	public ScorePolyA(Map<String, Collection<RefSeqGene>> genesByChr, String genomeDir, String save) throws Exception{
		Collection<Alignments> polyASites=new TreeSet<Alignments>();
		
		int numGenesWithPolyA=0;
		int randomPolyA=0;
		int totalGenes=0;
		
		Map<Alignments, Integer> map=new TreeMap();
		
		for(String chr: genesByChr.keySet()){
			try{
				totalGenes+=genesByChr.get(chr).size();
				System.err.println(chr);
				String sequenceFile=genomeDir+"/"+chr.replaceAll("chr", "").trim()+"/"+chr+".agp";
				Chromosome chrom=new Chromosome(sequenceFile);
				chrom.loadSequence();
				Collection<RefSeqGene> genes=genesByChr.get(chr);
				Collection<Alignments> sites=matchSites(genes, chrom);
				int numPolyA=countInLastExon(genes, sites, map);
				int random=countInFirstExon(genes, sites);
				randomPolyA+=random;
				numGenesWithPolyA+=numPolyA;
				polyASites.addAll(sites);
			}catch(FileNotFoundException ex){}
		}
		
		System.err.println("\nNumber of genes with polyA in last exon "+numGenesWithPolyA+" random "+randomPolyA+" "+map.size());
		write(save, polyASites);
	}
	
	private int countInFirstExon(Collection<RefSeqGene> genes, Collection<Alignments> sites) {
		Map<String, IntervalTree<Alignments>> sitesTree=CollapseByIntersection.makeIntervalTree(sites);
		
		int counter=0;
		for(RefSeqGene gene: genes){
			try{
			Alignments lastExon=gene.get5PrimeExon();
			Iterator<Node<Alignments>> iter=sitesTree.get(lastExon.getChr()).overlappers(lastExon.getStart(), lastExon.getEnd());
			if(iter.hasNext()){counter++;}
			}catch(Exception ex){}
		}
		
		return counter;
	}

	private int countInLastExon(Collection<RefSeqGene> genes, Collection<Alignments> sites, Map<Alignments, Integer> map) {
		Map<String, IntervalTree<Alignments>> sitesTree=CollapseByIntersection.makeIntervalTree(sites);
		
		int counter=0;
		for(RefSeqGene gene: genes){
			try{
			Alignments lastExon=gene.get3PrimeExon();
			Iterator<Node<Alignments>> iter=sitesTree.get(lastExon.getChr()).overlappers(lastExon.getStart(), lastExon.getEnd());
			if(iter.hasNext()){counter++; map.put(lastExon, 1);}
			}catch(Exception ex){System.err.println(gene);}
		}
		
		return counter;
	}

	private Collection<Alignments> matchSites(Collection<RefSeqGene> genes,	Chromosome chrom) throws Exception {
		Collection<Alignments> rtrn=new TreeSet<Alignments>();
		
		for(RefSeqGene gene: genes){
			rtrn.addAll(this.match(gene, chrom));
		}
		
		return rtrn;
	}

	private void write(String save, Collection<Alignments> polyASites) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(Alignments align: polyASites){
			writer.write(align+"\n");
		}
		
		writer.close();
	}

	private Collection<Alignments> match(RefSeqGene gene, Chromosome chrom) throws Exception{
		Collection<Alignments> rtrn=new TreeSet<Alignments>();
		SequenceMotif m=new SequenceMotif(motif,1);
		Collection<Alignments> exons=gene.getExonSet();
		for(Alignments exon: exons){
			String seq=ExtractSequence.getSequenceUnoriented(exon, chrom, false);
			Sequence sequence=new Sequence(exon.getName());
			if(gene.getOrientation().equalsIgnoreCase("+")){sequence.setSequenceBases(seq);}
			else{sequence.setSequenceBases(ComputeMIRScore.reverseComplement(seq));}
			List<SequenceRegion> regions=m.match(sequence);
			Collection<Alignments> matches=convert(regions, exon);
			rtrn.addAll(matches);
		}
		return rtrn;
	}
	
	private Collection<Alignments> convert(List<SequenceRegion> regions, Alignments exon) {
		Collection<Alignments> rtrn=new TreeSet<Alignments>();
		
		for(SequenceRegion region: regions){
			Alignments align=new Alignments(exon.getChr(), region.getAbsoluteStart(exon.getStart()), region.getAbsoluteEnd(exon.getStart()));
			rtrn.add(align);
		}
		
		return rtrn;
	}

	public static void main(String[] args)throws Exception{
		if(args.length>2){
			Map<String, Collection<RefSeqGene>> genesByChr=BEDFileParser.loadDataByChr(new File(args[0]));
			String genomeDir=args[1];
			String save=args[2];
			new ScorePolyA(genesByChr, genomeDir, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=genes \n args[1]=genome dir \n args[2]=save";
	
}
