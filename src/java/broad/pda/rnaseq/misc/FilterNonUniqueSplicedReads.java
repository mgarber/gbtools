package broad.pda.rnaseq.misc;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.reader.AlignmentQueryReader;
import org.broad.igv.sam.reader.BAMQueryReader;

import broad.pda.gene.RefSeqGene;
import net.sf.samtools.util.CloseableIterator;

//Find all reads that are spliced and placed in multiple locations

public class FilterNonUniqueSplicedReads {

	public FilterNonUniqueSplicedReads(String fileName, String save)throws IOException{
		FileWriter writer=new FileWriter(save);
		AlignmentQueryReader reader = new BAMQueryReader(new File(fileName));//SamQueryReaderFactory.getReader(new ResourceLocator(fileName), false);
		CloseableIterator<Alignment> iter=reader.iterator();
		
		Map<String, Alignment> names=new HashMap();
		
		int i=0;
		while(iter.hasNext()){
			Alignment next=iter.next();
			if(next.getAlignmentBlocks()!=null && next.getAlignmentBlocks().length>1){
				if(names.containsKey(next.getReadName())){names.remove(next.getReadName());}
				else{names.put(next.getReadName(), next);}
			}
			else{writer.write(toSam(next)+"\n");}
			
			if(i%100000 ==0){System.err.println(i);}
			i++;
		}
		
		iter.close();
		
		for(String name: names.keySet()){
			Alignment read=names.get(name);
			writer.write(toSam(read)+"\n");
		}
		
		writer.close();
	}
	
	private String toSam(Alignment read) {
		RefSeqGene gene=new RefSeqGene(read);
		return gene.toSAM();
	}

	public static void main(String[] args)throws IOException{
		if(args.length>1){
			String fileName=args[0];
			String save=args[1];
			new FilterNonUniqueSplicedReads(fileName, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=SAM file \n args[1]=save";
	
}
