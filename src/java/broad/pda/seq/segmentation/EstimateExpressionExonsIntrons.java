package broad.pda.seq.segmentation;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeSet;

import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;

public class EstimateExpressionExonsIntrons {

	public EstimateExpressionExonsIntrons(ContinuousDataAlignmentModel data, Collection<RefSeqGene> genes, String chr)throws IOException{
		Map<RefSeqGene, double[]> geneScores=data.scoreGenes(genes, chr);
		
		//Collection<RefSeqGene> intronGenes=flipExonIntrons(genes);
		//Map<RefSeqGene, double[]> intronScores=data.scoreGenes(intronGenes);
		//			
		//write(save+".exons.segments", geneScores);
		//write(save+".introns.segments", intronScores);
	}
	
	
	private Collection<RefSeqGene> flipExonIntrons(Collection<RefSeqGene> genes) {
		Collection<RefSeqGene> rtrn=new TreeSet();
		
		for(RefSeqGene gene: genes){
			Collection<Alignments> introns=gene.getIntronSet();
			RefSeqGene intron=new RefSeqGene(gene.getChr(), gene.getStart(), gene.getEnd(), gene.getName(), gene.getOrientation(), introns);
			rtrn.add(intron);
		}
		
		return rtrn;
	}


	private void write(String string, Map<RefSeqGene, double[]> geneScores) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		for(RefSeqGene gene: geneScores.keySet()){
			double[] scores=geneScores.get(gene);
			writer.write(gene.getName()+"\t"+scores[0]+"\t"+scores[1]+"\n");
		}
		
		writer.close();
	}


	public static void main(String[] args)throws IOException{
		if(args.length>3){
			ContinuousDataAlignmentModel model=new ContinuousDataAlignmentModel(new GenericAlignmentDataModel(args[0], args[1]));
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[2]));
			String chr=args[3];
			new EstimateExpressionExonsIntrons(model, genes, chr);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=alignment file \n args[1]=sizes \n args[2]=genes \n args[3]=chr";
	
}
