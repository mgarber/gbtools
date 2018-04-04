package broad.pda.seq.utils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

public class ScoreFullGenes {

	public ScoreFullGenes(Map<String, Collection<RefSeqGene>> regionsByChr, ContinuousDataAlignmentModel data, String save, String chr) throws IOException{
		Map<RefSeqGene, double[]> scores=new TreeMap();
		if(chr!=null){
			scores=data.scoreGenes(regionsByChr.get(chr), chr);
		}
		else{
			for(String chrom: regionsByChr.keySet()){
				scores.putAll(data.scoreGenes(regionsByChr.get(chrom), chrom));
			}
		}
		write(save, scores);
	}
	
	
	private void write(String save, Map<RefSeqGene, double[]> scores) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(RefSeqGene region: scores.keySet()){
			double[] vals=scores.get(region);
			writer.write(region+"\t"+vals[0]+"\t"+vals[1]+"\t"+vals[4]+"\n");
		}
		
		writer.close();
	}


	public static void main(String[] args)throws IOException{
		if(args.length>3){
			Map<String, Collection<RefSeqGene>> regionsByChr=BEDFileParser.loadDataByChr(new File(args[0]));
			ContinuousDataAlignmentModel data=new ContinuousDataAlignmentModel(new GenericAlignmentDataModel(args[1], args[2]));
			String save=args[3];
			String chr=null;
			if(args.length>4){chr=args[4];}
			new ScoreFullGenes(regionsByChr, data, save, chr);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=regions \n args[1]=alignment file \n args[2]=genome sizes \n args[3]=save \n args[4]=chr";
	
}
