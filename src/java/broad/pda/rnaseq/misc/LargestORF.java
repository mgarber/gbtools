package broad.pda.rnaseq.misc;

import java.io.File;
import java.io.IOException;
import java.util.Collection;

import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;

public class LargestORF {

	public LargestORF(Collection<RefSeqGene> genes){
		
			for(RefSeqGene gene: genes){
				/*Collection<Alignments> sorted=gene.getSortedAndUniqueExons();
				int i=0;
				for(Alignments exon: sorted){
					if(i==0 || i==sorted.size()-1){}
					else{System.out.println(exon.length());}
					i++;
				}*/
				System.out.println(gene.getTranscriptLength());
			}
		
	}
	
	public static void main(String[] args) throws IOException{
		Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[0]));
		new LargestORF(genes);
	}
	
}
