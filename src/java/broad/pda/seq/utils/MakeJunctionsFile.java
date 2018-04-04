package broad.pda.seq.utils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;

import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;

public class MakeJunctionsFile {

	public MakeJunctionsFile(File[] files, String save)throws IOException{
		Collection<Alignments> introns=new TreeSet();
		for(int i=0; i<files.length; i++){
			System.err.println(files[i]);
			Collection<RefSeqGene> genes=BEDFileParser.loadData(files[i]);
			introns.addAll(getIntrons(genes));
		}
		writeJunctions(save, introns);
	}
	
	private Collection<? extends Alignments> getIntrons(Collection<RefSeqGene> genes) {
		Collection<Alignments> rtrn=new TreeSet();
		
		for(RefSeqGene gene: genes){
			rtrn.addAll(gene.getIntronSet());
		}
		
		return rtrn;
	}

	private void writeJunctions(String save, Collection<Alignments> introns) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(Alignments align: introns){writer.write(align+"\t"+align.getStrand()+"\n");}
		
		writer.close();
	}

	public static void main(String[] args)throws IOException{
		if(args.length>1){
			File[] files=new File(args[0]).listFiles();
			String save=args[1];
			new MakeJunctionsFile(files, save);
		}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=files \n args[1]=save";
}
