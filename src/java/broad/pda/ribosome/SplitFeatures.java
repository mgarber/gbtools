package broad.pda.ribosome;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;

public class SplitFeatures {

	public static void main(String[] args)throws IOException{
		Collection<RefSeqGene> genes=BEDFileParser.loadData(new File("X:/Annotations/RefSeq.bed"));
		FileWriter writer=new FileWriter("X:/Annotations/introns.bed");
		for(RefSeqGene gene: genes){
			RefSeqGene introns=gene.getIntrons();
			if(introns!=null){
				introns.setName(introns.getName()+"_intron");
				writer.write(introns.toBED()+"\n");
			}
		}
		writer.close();
	}
	
}
