package broad.pda.seq.utils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;

import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;

public class MakeTableOfHits {

	public MakeTableOfHits(File[] files, String save)throws IOException{
		Collection<Alignments>[] regions=new Collection[files.length];
		Collection<Alignments> allRegions=new TreeSet();
		
		for(int i=0; i<files.length; i++){
			regions[i]=BEDFileParser.loadAlignmentData(files[i]);
			allRegions.addAll(regions[i]);
		}
		
		writeTable(save, allRegions, regions, files);
	}
	
	private void writeTable(String save, Collection<Alignments> allRegions, Collection<Alignments>[] regions, File[] files)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		writer.write("chr\tstart\tend");
		for(int i=0; i<files.length; i++){
			writer.write("\t"+files[i].getName());
		}
		writer.write("\n");
		
		for(Alignments align: allRegions){
			writer.write(align.toString());
			for(int i=0; i<regions.length;i++){
				writer.write("\t"+has(regions[i], align));
			}
			writer.write("\n");
		}
		
		writer.close();
	}
	
	private String has(Collection<Alignments> regions, Alignments align){
		if(regions.contains(align)){return "1";}
		else{return "0";}
	}
	
	
	public static void main(String[] args)throws IOException{
		if(args.length>1){
			File[] files=new File(args[0]).listFiles();
			String save=args[1];
			new MakeTableOfHits(files, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=files \n args[1]=save file";
}
