package broad.pda.seq.joshuaAnalyses;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.broad.igv.sam.Alignment;

import broad.pda.datastructures.Alignments;
import broad.pda.seq.fastq.FastqSequence;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

//Find the number of duplicate mate pairs in a sample
public class CountDuplicateMatePairs {

	public CountDuplicateMatePairs(AlignmentDataModel data, String save)throws IOException{
		
		Map<String, Integer> duplicates=new HashMap();
		Set<String> readNames=new TreeSet();
		for(String chr: data.getChromosomeLengths().keySet()){
			System.err.println(chr);
			Iterator<Alignment> reads=data.getAlignmentsOverlappingRegion(new Alignments(chr, 0, data.getChromosomeLengths().get(chr)));
			while(reads.hasNext()){
				Alignment align=reads.next();
				String name=(align.getReadName());
				String seq=(align.getReadSequence());
				if(!readNames.contains(name)){
					if(duplicates.containsKey(seq)){
						int num=duplicates.get(seq);
						num++;
						duplicates.put(seq, num);
					}
					else{duplicates.put(seq, 1);}
				}
				else{readNames.add(name);}
			}
			
		}
		write(save, duplicates);
	}
	
	
	private void write(String save, Map<String, Integer> duplicates) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		int numDuplicates=0;
		int total=0;
		
		for(String seq: duplicates.keySet()){
			int val=duplicates.get(seq);
			if(val>1){numDuplicates+=val;}
			total+=val;
			writer.write(seq+"\t"+val+"\n");
		}
		
		System.err.println(numDuplicates+" "+total);
		
		writer.close();
	}


	//TODO get duplicate reads by seq
	private int countDuplicates(Collection<FastqSequence> sequences) {
		Collection<String> unique=new TreeSet();
		for(FastqSequence seq: sequences){
			unique.add(seq.getSequence());
		}
		
		return unique.size();
	}


	public static void main(String[] args)throws IOException{
		AlignmentDataModel data=new GenericAlignmentDataModel(args[0], args[1]);
		String save=args[2];
		new CountDuplicateMatePairs(data, save);
	}

	
}
