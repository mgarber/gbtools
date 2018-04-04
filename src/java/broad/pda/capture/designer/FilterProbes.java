package broad.pda.capture.designer;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.TreeSet;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;

public class FilterProbes {

	public FilterProbes(String reportFile, String filterFile, List<Sequence> sequences, String save) throws IOException{
		Collection<String> report=BEDFileParser.loadList(reportFile);
		Collection<String> filter=BEDFileParser.loadList(filterFile);
		
		Collection<String> filterSet=new TreeSet<String>();
		for(String line: filter){
			filterSet.add(line.split("\t")[0]);
		}
		
		Collection<String> filteredSeq=new TreeSet<String>();
		for(Sequence seq: sequences){
			if(filterSet.contains(seq.getId())){
				//System.err.println(seq.getId()+" "+seq.getSequenceBases());
				filteredSeq.add(seq.getSequenceBases().toUpperCase());
			}
		}
		
		FileWriter filteredWriter=new FileWriter(save+".filtered.report");
		FileWriter writer=new FileWriter(save);
		
		for(String line: report){
			try{
			String probeSeq=line.split("\t")[6];
			String subseq=substring(probeSeq, 25, 144);
			//System.err.println(subseq.length()+" "+subseq); //TODO Is this the right size????
			if(filteredSeq.contains(subseq)){
				System.err.println("Found one");
				filteredWriter.write(line+"\t"+subseq+"\n");
			}
			else{writer.write(line+"\t"+subseq+"\n");}
			}catch(Exception ex){System.err.println("Skipping "+line);}
		}
		
		writer.close();
		filteredWriter.close();
		
	}
	
	private String substring(String probeSeq, int start, int end) {
		char[] bases=probeSeq.toCharArray();
		String rtrn="";
		
		for(int i=start; i<=end; i++){rtrn+=bases[i];}
		
		return rtrn.toUpperCase();
	}

	public static void main(String[] args)throws IOException{
		if(args.length>3){
			String reportFile=args[0];
			String filterFile=args[1];
			FastaSequenceIO fsio = new FastaSequenceIO(new File(args[2]));
			List<Sequence> seq= fsio.loadAll();
			String save=args[3];
			new FilterProbes(reportFile, filterFile, seq, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=probe report \n args[1]=blat filter \n args[2]=probe seq \n args[3]=save";
}
