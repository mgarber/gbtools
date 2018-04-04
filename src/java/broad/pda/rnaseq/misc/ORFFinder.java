package broad.pda.rnaseq.misc;

import java.io.File;
import java.util.Collection;
import java.util.Map;

import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;
import broad.pda.chromosome.Chromosome;
import broad.pda.chromosome.GenericOrganism;
import broad.pda.gene.RefSeqGene;

public class ORFFinder {
	
	
	public static void main(String[] args)throws Exception{
		if(args.length==2){
			Map<String, Collection<RefSeqGene>> genesByChr=BEDFileParser.loadDataByChr(new File(args[0]));
			String genomeDir=args[1];
			GenericOrganism go = new GenericOrganism(new File(genomeDir));
			for(String chr : genesByChr.keySet()) {
				Chromosome c = go.getChromosome(chr);
				if(c == null) {
					System.err.println(chr + " was not found in build directory " + genomeDir);
					continue;
				}
				System.err.println("Chromosome " + chr + " loading sequence .... ");
				c.loadSequence();
				System.err.println(" done");
				Sequence chrSeq = c.getSequence();
				
				Collection<RefSeqGene> chrGenes = genesByChr.get(chr);
				for(RefSeqGene g : chrGenes) {
					g.setSequenceFromChromosome(chrSeq);
					//System.err.println(">"+g.getName());
					//System.err.println(g.getSequence());
					RefSeqGene gCDS = g.findLongestORF();
					if(gCDS != null) {
						gCDS.addExtraField(String.valueOf(gCDS.getTranscriptLength()));
						gCDS.addExtraField(String.valueOf(gCDS.getTranscriptLength()/(double)g.getTranscriptLength()));
						System.out.println(gCDS.toBED());
					} else {
						RefSeqGene empty = new RefSeqGene(g.getChr(), g.getEnd()-1,g.getEnd());
						empty.setName(g.getName());
						empty.setOrientation(g.getOrientation());
						empty.addExtraField("0");
						empty.addExtraField("0");
						System.out.println(empty.toBED());
					}
				}
				c.unloadSequence();
			}
		}else{System.err.println(usage);}
	}
	
	static String usage="Output goes to standard out.\n args[0]=Annotations in full BED format \n args[1]=genome directory";

}
