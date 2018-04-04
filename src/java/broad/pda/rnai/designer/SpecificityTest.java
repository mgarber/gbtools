package broad.pda.rnai.designer;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.TreeSet;

public class SpecificityTest {

	
	public static void main(String[] args)throws IOException{
		if(args.length>5){
			
		double minScore=new Double(args[5]);
		Collection<String> kmers=parse(args[0], minScore);
		String geneSequence=args[1];
		String saveDir=args[2];
		
		String script=args[3];
		String queue=args[4];
		
		
		
		Runtime run=java.lang.Runtime.getRuntime();
		
		for(String kmer: kmers){
			String save=saveDir+"/"+kmer+".smatch";
			String command="bsub -q "+queue+" -o "+save+".bsub "+" java -jar -Xmx2000m "+script+" "+kmer+" "+geneSequence+" 7 "+save+" 3";
			System.err.println(command);
			run.exec(command);
			System.gc();
		}
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=kmers (all Kmers) \n args[1]=gene Sequence  (fa file) \n args[2]=save directory \n args[3]=script \n args[4]=queue \n args[5]=minScore";

	private static Collection<String> parse(String string, double minScore) throws IOException{
		Collection<String> rtrn=new TreeSet();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(string)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			double score=new Double(nextLine.split("\t")[2]);
			if(score>minScore){
				rtrn.add((nextLine.split("\t")[0]));
			}
		}
		reader.close();
		return rtrn;
	}
	
}
