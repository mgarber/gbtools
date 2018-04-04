package broad.pda.seq.utils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.Set;

import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

public class ScoreIntervals {

	public ScoreIntervals(Collection<Alignments> regions, ContinuousDataAlignmentModel data, String save) throws IOException{
		Map<Alignments, double[]> scores=data.scoreSegments(regions, null);
		write(save, scores);
	}
	
	
	private void write(String save, Map<Alignments, double[]> scores) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Alignments region: scores.keySet()){
			double[] vals=scores.get(region);
			writer.write(region+"\t"+vals[0]+"\t"+vals[1]+"\t"+vals[4]+"\n");
		}
		
		writer.close();
	}


	public static void main(String[] args)throws IOException{
		if(args.length>3){
			Set<Alignments> regions=BEDFileParser.loadAlignmentData(new File(args[0]));
			ContinuousDataAlignmentModel data=new ContinuousDataAlignmentModel(new GenericAlignmentDataModel(args[1], args[2]));
			String save=args[3];
			new ScoreIntervals(regions, data, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=regions \n args[1]=alignment file \n args[2]=genome sizes \n args[3]=save";
	
}
