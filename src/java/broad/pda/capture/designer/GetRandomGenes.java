package broad.pda.capture.designer;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;

public class GetRandomGenes {

	
	public GetRandomGenes(List<Sequence> seqs, int num, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		for(int i=0; i<num; i++){
			int random=new Double(Math.random()*seqs.size()).intValue();
			Sequence seq=seqs.remove(random);
			writer.write(">"+seq.getId()+"\n"+seq.getSequenceBases()+"\n");
		}
		writer.close();
	}

	public static void main(String[] args) throws IOException{
		if(args.length>2){
			FastaSequenceIO fsio = new FastaSequenceIO(new File(args[0]));
			List<Sequence> seq= fsio.loadAll();
			int num=new Integer(args[1]);
			String save=args[2];
			new GetRandomGenes(seq, num, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=genes \n args[1]=number random \n args[2]=save";
	
}
