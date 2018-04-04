package broad.pda.seq.alignment;

import java.io.IOException;

import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.reader.AlignmentQueryReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
//import org.broad.igv.sam.reader.SamQueryReaderFactory;
import org.broad.igv.util.ResourceLocator;

import net.sf.samtools.util.CloseableIterator;

public class AlignmentStats {
	private int duplicate;
	private int unique;
	private int paired;
	private int properlyPaired;
	private int [] mapQualDist;
	private int mapped;
	private int totalAlignments;
	private int spliced;
	//TODO: Count number of sequences that are aligned
	//TODO Count the number of sequences that are aligned uniquely
	//TODO count the number of sequences that are aligned to a splice junction
	public AlignmentStats(String file) throws IOException{
		init(file);
	}

	private  void init(String  file) throws IOException {
		AlignmentQueryReader reader = AlignmentReaderFactory.getReader(new ResourceLocator(file), false);
		CloseableIterator<Alignment> alignmentIterator = reader.iterator();
		mapQualDist = new int [256];
		while(alignmentIterator.hasNext()) {
			Alignment alignment = alignmentIterator.next();
			mapQualDist[alignment.getMappingQuality()]++;
			if(alignment.isDuplicate()){  duplicate++; } else { unique++ ;}
			if (alignment.isPaired()){paired++;}
			if (alignment.isProperPair()){ properlyPaired++;}
			if(alignment.isMapped()){mapped++;}
			if(alignment.getAlignmentBlocks().length > 1) { spliced++;}
			totalAlignments++;
		}
		alignmentIterator.close();
		reader.close();

		
	}
	
	public static void main(String[] args)throws IOException{
		AlignmentStats stats = new AlignmentStats(args[0]);
		
		System.out.println("unique: " + stats.unique);
		System.out.println("mapped: " + stats.mapped);
		System.out.println("paired: " + stats.paired);
		System.out.println("properlyPaired: " + stats.properlyPaired);
		System.out.println("duplicated: " + stats.duplicate);
		System.out.println("spliced: " + stats.spliced);
		System.out.println("total: " + stats.totalAlignments);
		for(int i = 0; i < stats.mapQualDist.length; i++) {
			System.out.println(i+"\t"+stats.mapQualDist[i]);
		}
	}
	
}
