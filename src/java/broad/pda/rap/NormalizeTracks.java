package broad.pda.rap;

import java.io.FileWriter;
import java.io.IOException;

import broad.core.parser.CommandLineParser;
import broad.core.sequence.SequenceUtils;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;

public class NormalizeTracks {

	int windowSize=500;
	String genome="mm9";
	
	public NormalizeTracks(ContinuousDataAlignmentModel data1, ContinuousDataAlignmentModel data2, String save, String sizes, String chr, int windowSize) throws IOException, InterruptedException{
		this.windowSize=windowSize;
		//For each window compute median at center point
		windowSlider(data1, data2, windowSize, save+".wig", chr);
				
		//make bigWig
		makeBigWig(save, sizes);
	}
	
	
	private void makeBigWig(String save, String sizes) throws IOException, InterruptedException {
		String cmd="/seq/lincRNA/scripts/UCSC/bedGraphToBigWig "+save+".wig "+sizes+" "+save+".bw";
		Runtime run=Runtime.getRuntime();
		Process p=run.exec(cmd);
		p.waitFor();
	}
	

	private void windowSlider(ContinuousDataAlignmentModel data1, ContinuousDataAlignmentModel data2, int windowSize, String save, String chrToUse) throws IOException {
		if(chrToUse != null) {
			if(!chrToUse.startsWith("chr")){chrToUse=null;}
		}
		
		FileWriter writer=new FileWriter(save);
		
		for(String chr: data1.getChromosomeLengths().keySet()){
			if(chrToUse==null || chr.equalsIgnoreCase(chrToUse)){
				System.err.println(chr);
				//writer.write("variableStep  chrom="+chr+" [span=1]\n");
				data1.scanGenome(data2, windowSize, chr, writer, false);
			}
		}
		
		writer.close();
	}

	 

	public static void main(String[] args)throws IOException, InterruptedException{
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-d1", "data1 bam file (sense)", true);
		p.addStringArg("-d2", "data2 bam file (antisense)", true);
		p.addStringArg("-s", "chr size file", true);
		p.addStringArg("-o", "output bigwig file", true);
		p.addStringArg("-c", "chromosome", false, null);
		p.addIntegerArg("-w", "window size", false, 500);
		p.parse(args);
		String bam1 = p.getStringArg("-d1");
		String bam2 = p.getStringArg("-d2");
		String sizes = p.getStringArg("-s");
		String output = p.getStringArg("-o");
		String chr = p.getStringArg("-c");
		int window = p.getIntegerArg("-w").intValue();
		
		ContinuousDataAlignmentModel data1 = SequenceUtils.getDataModel(bam1, sizes, true);
		ContinuousDataAlignmentModel data2 = SequenceUtils.getDataModel(bam2, sizes, true);
		
		new NormalizeTracks(data1, data2, output, sizes, chr, window);
		
	}
	
	static String usage=" args[0]=bam (sense) \n args[1]=bam (antisense) \n args[2]=sizes \n args[3]=save \n args[4]=chr (optional) \n args[5]=windowSize (optional)";
	
}
