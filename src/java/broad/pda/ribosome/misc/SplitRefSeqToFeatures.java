package broad.pda.ribosome.misc;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;


//Take a RefSeqGene and split it into seperate files for each feature
public class SplitRefSeqToFeatures {

	public SplitRefSeqToFeatures(Collection<RefSeqGene> genes, String save) throws IOException{
		FileWriter cds=new FileWriter(save+".cds.bed");
		FileWriter utr5=new FileWriter(save+".utr5.bed");
		FileWriter utr3=new FileWriter(save+".utr3.bed");
		FileWriter intron=new FileWriter(save+".intron.bed");
		
		for(RefSeqGene gene: genes){
			if(gene.getCDS()!=null){cds.write(gene.getCDS()+"\n");}
			if(gene.get5UTRGene()!=null){utr5.write(gene.get5UTRGene()+"\n");}
			if(gene.get3UTRGene()!=null){utr3.write(gene.get3UTRGene()+"\n");}
			if(gene.getIntrons()!=null){intron.write(gene.getIntrons()+"\n");}
		}
		
		cds.close();
		utr5.close();
		utr3.close();
		intron.close();
	}
	
	
	public static void main(String[] args)throws IOException{
		if(args.length>1){
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[0]));
			String save=args[1];
			new SplitRefSeqToFeatures(genes, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=BED file \n args[1]=save";
	
}
