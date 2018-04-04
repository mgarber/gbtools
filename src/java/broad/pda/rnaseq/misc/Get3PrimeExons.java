package broad.pda.rnaseq.misc;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;

public class Get3PrimeExons {

	public Get3PrimeExons(Collection<RefSeqGene> genes, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		for(RefSeqGene gene: genes){
			writer.write(gene.get3PrimeExon()+"\n");
		}
		writer.close();
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>1){
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[0]));
			String save=args[1];
			new Get3PrimeExons(genes, save);
		}
	}
	
}
