package broad.pda.rnaseq.misc;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;

import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;

public class GetExonsIntrons {

	public GetExonsIntrons(Collection<RefSeqGene> genes, String save) throws IOException{
		Collection<Alignments> exons=new TreeSet<Alignments>();
		Collection<Alignments> introns=new TreeSet<Alignments>();
		
		for(RefSeqGene gene: genes){
			exons.addAll(gene.getExonSet());
			introns.addAll(gene.getIntronSet());
		}
		
		FileWriter writer=new FileWriter(save+".exons.bed");
		for(Alignments exon: exons){
			writer.write(exon+"\n");
		}
		writer.close();
		
		writer=new FileWriter(save+".introns.bed");
		for(Alignments intron: introns){
			writer.write(intron+"\n");
		}
		writer.close();
		
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>1){
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[0]));
			String save=args[1];
			new GetExonsIntrons(genes, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=full bed \n args[1]=save";
}
