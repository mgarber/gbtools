package broad.pda.capture.designer;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;

public class FastaToText {

	public FastaToText(File[] files, String saveDir) throws IOException{
		
		for(int i=0; i<files.length; i++){
			FileWriter writer=new FileWriter(saveDir+"/"+files[i].getName()+".txt");
			FastaSequenceIO fsio = new FastaSequenceIO(files[i]);
			Collection<Sequence> seqs= fsio.loadAll();
			for(Sequence seq: seqs){writer.write(seq.getSequenceBases()+"\n");}
			writer.close();
		}
		
	}
	
	public static void main(String[] args)throws IOException{
		File[] files=new File(args[0]).listFiles();
		String save=args[1];
		new FastaToText(files, save);
	}
	
}
