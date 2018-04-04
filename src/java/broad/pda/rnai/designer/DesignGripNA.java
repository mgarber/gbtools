package broad.pda.rnai.designer;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import broad.core.motif.SearchException;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.core.sequence.WindowSlider;
import broad.pda.annotation.BEDFileParser;

public class DesignGripNA {

	int probeLength=20;
	int seed=7;

	public DesignGripNA(Sequence region, List<Sequence> genes, String save) throws SearchException, IOException{
		FileWriter writer=new FileWriter(save);
		
		WindowSlider slider=WindowSlider.getSlider(region, probeLength, probeLength-1);
		int counter=0;
		while(slider.hasNext()){
			SequenceRegion r=slider.next();
			Sequence probe=r.getSequence();
			System.err.println(counter+" "+probe.getSequenceBases());
			SmatchLike smatch=new SmatchLike(probe.getSequenceBases(), genes, seed);
			Collection<String> matches=smatch.getForwardTargets(2);
			writer.write(Sequence.reverseSequence(probe.getSequenceBases()));
			for(String match: matches){writer.write("\t"+match);}
			writer.write("\n");
			counter++;
			writer.flush();
		}
		
		writer.close();
	}
	
	public DesignGripNA(Sequence region, Collection<String> specificKmers, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		Collection<HairpinKmer> kmers=new ArrayList<HairpinKmer>();
		WindowSlider slider=WindowSlider.getSlider(region, probeLength, probeLength-1);
		while(slider.hasNext()){
			SequenceRegion r=slider.next();
			Sequence probe=r.getSequence();
			String seq=Sequence.reverseSequence(probe.getSequenceBases());
			if(specificKmers.contains(seq)){
				HairpinKmer kmer=new HairpinKmer(seq, r.getStart(), r.getEnd());
				int maxComplement=kmer.maxSelfComplementarity();
				writer.write(kmer.getKmerSequence()+"\t"+kmer.getKmerStartPosition()+"\t"+kmer.getKmerEndPosition()+"\t"+kmer.getPercentGC()+"\t"+maxComplement+"\n");
				kmers.add(kmer);
			}
			
		}
		writer.close();
	}
	
	public static void main(String[] args)throws IOException, SearchException{
		if(args.length>2){
			Sequence seq= new FastaSequenceIO(new File(args[0])).loadAll().iterator().next();
			Collection<String> genes=BEDFileParser.loadList(args[1]);
			String save=args[2];
			new DesignGripNA(seq, genes, save);
		}
		else{
			System.err.println(usage);
		}
	}
	static String usage=" args[0]=seq \n args[1]=good kmers \n args[2]=save";
}
