package broad.pda.seq.pairedend;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.TreeMap;

import org.broad.igv.sam.Alignment;

import broad.pda.datastructures.Alignments;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;
import net.sf.samtools.util.CloseableIterator;

public class PairedEndMappingEfficient {

	private int maxDistance=1000000;

	//We need to load all the reads to get the names and then put them together
	//Since we require them to be on the same chromosome we only need to load 1 chromosomes worth of data
	public PairedEndMappingEfficient(AlignmentDataModel data1, AlignmentDataModel data2, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(String chr: data1.getChromosomeLengths().keySet()){
			System.err.println("Working on "+chr);
			Alignments chrRegion=new Alignments(chr, 0, data1.getChromosomeLength(chr));
			CloseableIterator<org.broad.igv.sam.Alignment> reads1=data1.getAlignmentsOverlappingRegion(chrRegion);
			CloseableIterator<Alignment> reads2=data1.getAlignmentsOverlappingRegion(chrRegion);
			Map<String, Collection<Alignment>> readsByName=getReadsByName(reads1);
			writePairs(readsByName, reads2, writer);
		}
		
		writer.close();
	}

	private Map<String, Collection<Alignment>> getReadsByName(CloseableIterator<Alignment> reads) {
		Map<String, Collection<Alignment>> rtrn=new TreeMap<String, Collection<Alignment>>();
		
		while(reads.hasNext()){
			Alignment read=reads.next();
			String name=read.getReadName();
			Collection<Alignment> list=new ArrayList<Alignment>();
			if(rtrn.containsKey(name)){list=rtrn.get(name);}
			list.add(read);
			rtrn.put(name, list);
		}
		
		return rtrn;
	}

	private void writePairs(Map<String, Collection<Alignment>> readsByName,	CloseableIterator<Alignment> reads2, FileWriter writer) throws IOException {
		while(reads2.hasNext()){
			Alignment right=reads2.next();
			Collection<Alignment> leftReads=readsByName.get(right.getReadName());
			Collection<PairedEndAlignment> pairs=getAllPossiblePairs(leftReads, right);
			for(PairedEndAlignment pair: pairs){writer.write(pair.toString()+"\n");}
		}
		
	}
	
	private Collection<PairedEndAlignment> getAllPossiblePairs(Collection<Alignment> leftReads, Alignment right){
		Collection<PairedEndAlignment> rtrn=new HashSet<PairedEndAlignment>();
		for(Alignment left: leftReads){
			PairedEndAlignment pair=new PairedEndAlignment(left, right);
			if(pair.onSameChromosome() && pair.getDistance()<maxDistance && pair.oppositeStrands()){rtrn.add(pair);}
		}
		return rtrn;
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>3){
			AlignmentDataModel data1=new GenericAlignmentDataModel(args[0], args[2]);
			AlignmentDataModel data2=new GenericAlignmentDataModel(args[1], args[2]);
			String save=args[3];
			new PairedEndMappingEfficient(data1, data2, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]= alignments 1 \n args[1]=alignments 2 \n args[2]=sizes \n args[3]=save";
}
