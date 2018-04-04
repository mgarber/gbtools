package broad.pda.rnaseq.expression;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

public class EstimateExpression {

	public EstimateExpression(AlignmentDataModel data, Collection<RefSeqGene> genes, String save, String chr)throws IOException{
		ContinuousDataAlignmentModel model=new ContinuousDataAlignmentModel(data);
		Map<RefSeqGene, double[]> geneExpression=new TreeMap();
		//for(String chr: data.getChromosomeLengths().keySet()){
		//	System.err.println(chr);
			geneExpression.putAll(model.scoreGenes(genes, chr));
		//}
		write(save, geneExpression);
	}
	
	private void write(String save, Map<RefSeqGene, double[]> geneExpression)throws IOException{
		FileWriter writer=new FileWriter(save);
		for(RefSeqGene gene: geneExpression.keySet()){
			double[] vals=geneExpression.get(gene);
			writer.write(gene.toBED()+"\t"+vals[0]+"\t"+vals[1]+"\n");
		}
		writer.close();
	}
	
	
	public static void main(String[] args)throws IOException{
		if(args.length>4){
			AlignmentDataModel data=new GenericAlignmentDataModel(args[0], args[1]);
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[2]));
			String save=args[3];
			String chr=args[4];
			new EstimateExpression(data, genes, save, chr);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=alignment file \n args[1]=sizes \n args[2]=genes \n args[3]=save \n args[4]=chr";
	
}
