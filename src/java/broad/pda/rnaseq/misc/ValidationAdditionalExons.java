package broad.pda.rnaseq.misc;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.motif.SequenceMotif;
import broad.core.primer3.Primer3SequenceInputTags.SequenceRegionCoordinates;
import broad.core.primer3.PrimerPair;
import broad.core.primer3.qPCRPrimerDesigner;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.pda.annotation.BEDFileParser;
import broad.pda.chromosome.Chromosome;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.rnai.ExtractSequence;
import broad.pda.rnai.designer.ComputeMIRScore;

public class ValidationAdditionalExons {
	boolean repeatMask=true;
	int numDesigns=2;
	
	public ValidationAdditionalExons(Collection<RefSeqGene> genes, Collection<Alignments> exons, String genomeDir, String save) throws Exception{
		Map<RefSeqGene, Collection<PrimerPair>> map=new TreeMap<RefSeqGene, Collection<PrimerPair>>();
		for(RefSeqGene gene: genes){
			System.err.println(gene.getAlignment().toUCSC());
			Collection<Alignments> targets=findOverlappingExons(gene, exons);
			//design primers with boundary as target region
			Collection<PrimerPair> primers=designPrimers(gene, targets, genomeDir);
			map.put(gene, primers);
		}
		write(save, map);
	}
	private Collection<Alignments> findOverlappingExons(RefSeqGene gene, Collection<Alignments> exons) {
		Collection<Alignments> rtrn=new TreeSet();
		
		for(Alignments exon: exons){
			if(exon.overlaps(gene.getAlignment())){rtrn.add(exon);}
		}
		
		return rtrn;
	}
	
	private Collection<PrimerPair> designPrimers(RefSeqGene gene, Collection<Alignments> targets, String genomeDir) throws Exception {
		String sequenceFile=genomeDir+"/"+gene.getChr().replaceAll("chr", "").trim()+"/"+gene.getChr()+".agp";
		Chromosome chrom = new Chromosome(sequenceFile);
		chrom.loadSequence();
		
		String seq=ExtractSequence.getSequenceForGene(gene, chrom, repeatMask, new TreeMap());
		Sequence geneSequence=new Sequence(gene.getName());
		geneSequence.setSequenceBases(seq);
		
		
		
		Collection<SequenceRegion> regions=new HashSet();
		for(Alignments exon: targets){
			SequenceMotif motif=new SequenceMotif(ExtractSequence.getSequenceUnoriented(exon, chrom, repeatMask), 1);
			if(gene.getOrientation().equalsIgnoreCase("-")){
				motif=new SequenceMotif(ComputeMIRScore.reverseComplement(ExtractSequence.getSequenceUnoriented(exon, chrom, repeatMask)), 1);
			}
			List<SequenceRegion> list=motif.match(geneSequence);
			System.err.println(list);
			regions.addAll(list);
		}
		
		SequenceRegionCoordinates target=makeUnionRegion(regions);
		Collection<PrimerPair> primers=qPCRPrimerDesigner.designCloningPrimers(geneSequence, numDesigns, target);
		return primers;
	}
	
	private SequenceRegionCoordinates makeUnionRegion(Collection<SequenceRegion> regions) {
		
		int start=Integer.MAX_VALUE;
		int end=-Integer.MAX_VALUE;
		boolean started=false;
		for(SequenceRegion region: regions){
			start=Math.min(region.getStart(), start);
			end=Math.max(region.getEnd(), end);
			started=true;
		}
		
		if(!started){return null;}
		return new SequenceRegionCoordinates(start, end);
	}
	private void write(String save, Map<RefSeqGene, Collection<PrimerPair>> map) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(RefSeqGene gene: map.keySet()){
			Collection<PrimerPair> primers=map.get(gene);
			for(PrimerPair primer: primers){
				writer.write(gene.getName()+"\t"+primer.getLeftPrimer()+"\t"+primer.getRightPrimer()+"\t"+primer.getPrimerPairPenalty()+"\t"+primer.getProductSize()+"\n");
			}
		}
		
		writer.close();
	}
	
	public static void main(String[] args)throws Exception{
		if(args.length>3){
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[0]));
			Collection<Alignments> exons=BEDFileParser.loadAlignmentData(new File(args[1]));
			String genomeDir=args[2];
			String save=args[3];
			new ValidationAdditionalExons(genes, exons, genomeDir, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=genes \n args[1]=additional exons \n args[2]=genome directory \n args[3]=save";

	
}
