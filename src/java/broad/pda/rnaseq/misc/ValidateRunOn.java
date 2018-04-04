package broad.pda.rnaseq.misc;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import broad.core.primer3.Primer3SequenceInputTags.SequenceRegionCoordinates;
import broad.core.primer3.PrimerPair;
import broad.core.primer3.qPCRPrimerDesigner;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.rnai.ExtractSequence;

//Design primers to validate run-on transcripts

public class ValidateRunOn {
	boolean repeatMask=true;
	int numDesigns=2;

	public ValidateRunOn(Collection<RefSeqGene> runonTranscripts, Collection<Alignments> boundaries, String genomeDir, String save)throws Exception{
		Map<RefSeqGene, Collection<PrimerPair>> map=new TreeMap<RefSeqGene, Collection<PrimerPair>>();
		
		for(RefSeqGene gene: runonTranscripts){
			System.err.println(gene.getAlignment().toUCSC());
			Alignments target=findOverlappingBoundary(gene, boundaries);
			//design primers with boundary as target region
			Collection<PrimerPair> primers=designPrimers(gene, target, genomeDir);
			map.put(gene, primers);
		}
		write(save, map);
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

	private Collection<PrimerPair> designPrimers(RefSeqGene gene, Alignments target, String genomeDir) throws Exception {
		String seq=ExtractSequence.getSequenceForGene(gene, genomeDir, repeatMask);
		Sequence geneSequence=new Sequence(gene.getName());
		geneSequence.setSequenceBases(seq);
		
		Collection<Alignments> exons=gene.getSortedAndUniqueExons();
		int relativeStart=0;
		for(Alignments exon: exons){
			if(gene.getOrientation().equalsIgnoreCase("+")){
				if(exon.getStart()<target.getStart() && exon.getEnd()<target.getEnd()){relativeStart+=exon.getSize();}
			}
			else{
				if(gene.getOrientation().equalsIgnoreCase("-")){
					if(exon.getStart()>target.getStart() && exon.getEnd()>target.getEnd()){relativeStart+=exon.getSize();}
				}
			}
		}
		
		SequenceRegionCoordinates targetRegion=new SequenceRegionCoordinates(relativeStart-100, relativeStart+100);
		
		
		Collection<PrimerPair> primers=qPCRPrimerDesigner.designCloningPrimers(geneSequence, numDesigns, targetRegion);
		return primers;
	}

	private Alignments findOverlappingBoundary(RefSeqGene gene,	Collection<Alignments> boundaries) {
		for(Alignments target: boundaries){
			if(target.overlaps(gene.getAlignment())){return target;}
		}
		return null;
	}

	public static void main(String[] args)throws Exception{
		if(args.length>3){
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[0]));
			Collection<Alignments> boundaries=BEDFileParser.loadAlignmentData(new File(args[1]));
			String genomeDir=args[2];
			String save=args[3];
			new ValidateRunOn(genes, boundaries, genomeDir, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=runOn Transcripts \n args[1]=boundary regions \n args[2]=genome directory \n args[3]=save";
}
