package broad.pda.rap;

import java.io.IOException;
import java.util.Map;

import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;

public class BatchSubmitAlignments {

	int chunkSize=5000000; //1Mb
	
	//Iterate through chromosome and batch out chunks
	public BatchSubmitAlignments(String inputFasta, String genomeDirectory, String save, int windowSize, int overlap) throws IOException, InterruptedException{
		Map<String, Integer> chrSizes=BEDFileParser.loadChrSizes(genomeDirectory+"/sizes");
		Runtime run=Runtime.getRuntime();
		
		for(String chr: chrSizes.keySet()){
			System.err.println(chr);
			int chrSize=chrSizes.get(chr);
			for(int i=0; i<chrSize; i+=chunkSize){
				int start=i;
				int end=Math.min(start+chunkSize, chrSize);
				Alignments align=new Alignments(chr, start, end);
				//submit job
				String local=save+"."+align.toFileName();
				String command="bsub -q hour -W 4:00 -o "+ local+".bsub"+" java -jar /seq/lincRNA/scripts/RAP/ComputeCrossHybConsumer.jar "+windowSize+" "+overlap+" 100 "+inputFasta+" "+genomeDirectory+" "+align.toUCSC()+" "+ local+".bedgraph";
				System.err.println(command);
				Process p=run.exec(command);
				p.waitFor();
			}
		}
	}
	
	public static void main(String[] args)throws IOException, InterruptedException{
		if(args.length>4){
			new BatchSubmitAlignments(args[0], args[1], args[2], new Integer(args[3]), new Integer(args[4]));
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=input fasta \n args[1]=genome directory \n args[2]=save \n args[3]=window size \n args[4]=overlap";
	
}
