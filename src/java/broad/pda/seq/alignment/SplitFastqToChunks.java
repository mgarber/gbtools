package broad.pda.seq.alignment;

import java.io.File;
import java.io.IOException;

import broad.pda.seq.fastq.FastqParser;

public class SplitFastqToChunks {

	int chunkSize=10000000;
	
	public SplitFastqToChunks(File file, String save, int chunkSize)throws IOException{
		this.chunkSize=chunkSize;
		FastqParser fastq=new FastqParser(file);
		fastq.writeChunks(save, chunkSize);
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
			File fastqFile=new File(args[0]);
			String save=args[1];
			int chunkSize=new Integer(args[2]);
			new SplitFastqToChunks(fastqFile, save, chunkSize);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=fastq file \n args[1]=save \n args[2]=chunk size";
}
