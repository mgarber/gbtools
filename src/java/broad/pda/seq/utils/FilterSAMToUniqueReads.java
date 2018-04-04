package broad.pda.seq.utils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;

import org.broad.igv.sam.Alignment;

import broad.core.datastructures.IntervalTree;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;


//This class will take a "Generic SAM" file and remove all reads that overlap the same position in the genome
public class FilterSAMToUniqueReads {

	public FilterSAMToUniqueReads(File SAMFile, String save, String chrSizes)throws IOException{
	 	AlignmentDataModel data=new GenericAlignmentDataModel(SAMFile.getAbsolutePath(), chrSizes);
		filterAndWrite(save, data);
	}
	
	private void filterAndWrite(String save, AlignmentDataModel data)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(String chr: data.getChromosomeLengths().keySet()){
			System.err.println(chr);
			IntervalTree<Alignment> tree=data.getIntervalTree(chr, 0, data.getChromosomeLengths().get(chr));  
			Iterator<Alignment> iter=tree.valueIterator();
			while(iter.hasNext()){
				Alignment next=iter.next();
				writer.write(next+"\n");
			}
		}
		
		writer.close();
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
			File SAMFile=new File(args[0]);
			String save=args[1];
			String chrSizes=args[2];
			//for(int i=0; i<SAMFiles.length; i++){
			//	System.err.println(SAMFiles[i]);
				//String save=saveDir+"/"+SAMFiles[i].getName();
				new FilterSAMToUniqueReads(SAMFile, save, chrSizes);
			//}
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=SAM file \n args[1]=save file \n args[2]=chr sizes";
	
}
