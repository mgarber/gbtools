package broad.pda.seq.joshuaAnalyses;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Map;

import broad.core.datastructures.IntervalTree;
import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.alignment.sam.SAMRecord;

public class Get3PrimeReads {
	
	int size=200;

	public Get3PrimeReads(File samFile, File geneFile, String save) throws IOException{
		Map<String, IntervalTree<RefSeqGene>> geneTree=BEDFileParser.makeIntervalTreeFor3Prime(geneFile, size);
		FileWriter writer=new FileWriter(save);
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(samFile)));
		String nextLine;
		int i=0;
		int kept=0;
		while ((nextLine = reader.readLine()) != null) {
			SAMRecord sam=new SAMRecord(nextLine);
			RefSeqGene read=sam.getGene();
			boolean keep=geneTree.get(read.getChr()).overlappers(read.getStart(), read.getEnd()).hasNext();
			if(keep){writer.write(sam.toString()+"\n"); kept++;}
			if(i%1000000 ==0){System.err.println("went through "+i+" kept "+kept);}
			i++;
		}
		
		writer.close();
		
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
			File sam=new File(args[0]);
			File geneFile=new File(args[1]);
			String save=args[2];
			new Get3PrimeReads(sam, geneFile, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=sam file \n args[1]=gene file \n args[2]=save";
}
