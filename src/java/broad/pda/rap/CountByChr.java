package broad.pda.rap;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;

import org.broad.igv.Globals;

import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

public class CountByChr {

	public CountByChr(AlignmentDataModel model, String save, Alignments exclude) throws IOException{
		Globals.setHeadless(true);
		FileWriter writer=new FileWriter(save);
		writer.write("Chromosome\tLength\tcount\tExcludeCount\n");
		
		Map<String, Integer> chrs=model.getChromosomeLengths();
		for(String chr: chrs.keySet()){
			System.err.println(chr);
			double excludeCount=0;
			if(chr.equalsIgnoreCase(exclude.getChr())){
				excludeCount=model.getCountsPerAlignment(exclude, 0);
			}
			double count=model.getCounts(chr);
			
			writer.write(chr+"\t"+chrs.get(chr)+"\t"+count+"\t"+(count-excludeCount)+"\n");
		}
		writer.close();
	}
	
	public static void main(String[] args)throws IOException{
		// NOTE from Jesse - You could use PrecomputeDataAlignmentStats instead
		if(args.length>3){
			File[] files=new File(args[0]).listFiles();
			String sizes=args[1];
			String saveDir=args[2];
			for(int i=0; i<files.length; i++){
				try{
				if(files[i].getName().endsWith("bam")){
					System.err.println(files[i]);
					AlignmentDataModel data=new GenericAlignmentDataModel(files[i].getAbsolutePath(), sizes);
					String save=saveDir+"/"+files[i].getName()+".count";
					Alignments exclude=BEDFileParser.loadAlignmentData(new File(args[3])).iterator().next();
					new CountByChr(data, save, exclude);
				}
			}catch(Exception ex){System.err.println("Skipping "+files[i]);}
			}
		}
		else{System.err.println(usage);}
	}
	
	static String usage="args[0]=alignment files \n args[1]=sizes \n args[2]=save \n args[3]=excludeRegion";
	
}
