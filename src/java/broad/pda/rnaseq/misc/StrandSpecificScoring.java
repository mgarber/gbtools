package broad.pda.rnaseq.misc;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

public class StrandSpecificScoring {

	public StrandSpecificScoring(ContinuousDataAlignmentModel data, String save, Collection<RefSeqGene> genes)throws IOException{
		Map<RefSeqGene, double[]> strandScore=data.scoreGenesStranded(genes);

		write(save, strandScore, .05);
		
	}
	
	private void write(String save, Map<RefSeqGene, double[]> strandScore, double alpha) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(RefSeqGene gene: strandScore.keySet()){
			double[] p=strandScore.get(gene);
			if(p[0]<alpha){writer.write(gene+"\n");}
		}
		
		writer.close();
	}

	public static void main(String[] args)throws IOException{
		if(args.length>3){
			ContinuousDataAlignmentModel data=new ContinuousDataAlignmentModel(new GenericAlignmentDataModel(args[0], args[1]));
			String save=args[2];
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[3]));
			new StrandSpecificScoring(data, save, genes);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=alignments \n args[1]=sizes \n args[2]=save \n args[3]=genes";
}
