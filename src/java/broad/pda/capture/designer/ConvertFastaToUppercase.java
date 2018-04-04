package broad.pda.capture.designer;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;

public class ConvertFastaToUppercase {

	public static void main (String[] args) throws IOException{
		if(args.length>1){
			FastaSequenceIO fsio = new FastaSequenceIO(new File(args[0]));
			Collection<Sequence> seqs= fsio.loadAll();
			String save=args[1];
			
			FileWriter writer=new FileWriter(save);
			for(Sequence seq: seqs){
				writer.write(">"+seq.getId()+"\n"+seq.getSequenceBases().toUpperCase()+"\n");
			}
			writer.close();
		
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=fasta \n args[1]=save";
	
}
