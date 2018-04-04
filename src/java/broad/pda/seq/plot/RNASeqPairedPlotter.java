package broad.pda.seq.plot;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Set;

import org.jibble.epsgraphics.EpsGraphics2D;

import broad.pda.alignment.PairedEndAlignment;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;


public class RNASeqPairedPlotter {

	public static void plot(Alignments region, String save, Set<RefSeqGene> pairedReads)throws IOException{
		RNASeqDotPlotter plot = new RNASeqDotPlotter(region);
		EpsGraphics2D epsG2D = plot.startPlot(save);
		
		for(RefSeqGene pair: pairedReads){
			//System.err.println(pair.getChr()+" "+pair.getAlignment());
			plot.plot(pair, epsG2D);
		}
		
		plot.endPlot(epsG2D);
	}
	
	public static void plot(Alignments region, String save, Collection<PairedEndAlignment> pairedReads)throws IOException{
		RNASeqDotPlotter plot = new RNASeqDotPlotter(region);
		EpsGraphics2D epsG2D = plot.startPlot(save);
		
		for(PairedEndAlignment pair: pairedReads){
			plot.plot(pair, epsG2D);
		}
		
		plot.endPlot(epsG2D);
	}
	
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
		Alignments region=new Alignments(args[0]);
		String save=args[1];
		Set<RefSeqGene> genes=BEDFileParser.loadData(new File(args[2]), region);
		
		
		//File samFile1=new File(args[2]);
		//File samFile2=new File(args[3]);
		
		//Collection<PairedEndAlignment> pairs=PairedEndMapping.getPairs(samFile1, samFile2, region);
		plot(region, save, genes);
		}
		//else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=region \n args[1]=save \n args[2]=sam file 1 \n args[3]=sam file 2";
	
}
