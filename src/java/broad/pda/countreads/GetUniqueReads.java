package broad.pda.countreads;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;

import org.broad.igv.sam.Alignment;

import broad.pda.datastructures.Alignments;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;
import net.sf.samtools.util.CloseableIterator;


//Take a sam file and flag all regions that are multimappers
public class GetUniqueReads {

	public GetUniqueReads(AlignmentDataModel data, String save, String chr)throws IOException{
		Collection<String> multimappers=filterMultiMappers(data, chr);
		write(save, data, multimappers, chr);
	}
	
	private void write(String save, AlignmentDataModel data, Collection<String> multimappers, String chr) throws IOException {
		FileWriter writer=new FileWriter(save);
		CloseableIterator<Alignment> iter=data.getAlignmentsOverlappingRegion(new Alignments(chr, 0, data.getChromosomeLengths().get(chr)));
		while(iter.hasNext()){
			Alignment read=iter.next();
			if(!multimappers.contains(read.getReadName())){writer.write(read.toString());}
		}
		iter.close();
		writer.close();
	}

	public Collection<String> filterMultiMappers(AlignmentDataModel data, String chr) throws IOException{
		Collection<String> temp=new TreeSet();
		Collection<String> multimapper=new TreeSet();
		CloseableIterator<Alignment> iter=data.getAlignmentsOverlappingRegion(new Alignments(chr, 0, data.getChromosomeLengths().get(chr)));
		while(iter.hasNext()){
			Alignment read=iter.next();
			if(temp.contains(read.getReadName())){multimapper.add(read.getReadName());}
			else{temp.add(read.getReadName());}
		}
		iter.close();
	
		System.err.println("Aligned Reads: "+temp.size());
		temp=null;
		System.gc();
		return multimapper;
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>3){
			AlignmentDataModel data=new GenericAlignmentDataModel(args[0], args[1]);
			String save=args[2];
			String chr=args[3];
			new GetUniqueReads(data, save, chr);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=alignment file \n args[1]=sizes \n args[2]=save \n args[3]=chr";
	
}
