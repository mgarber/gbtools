package broad.pda.seq.utils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;


//Need to fix this to handle full BED files so it doesnt throw out isoforms
public class OverlappingRegions {
	Map info;

	
	public OverlappingRegions(File repeatFile, File intervalFile, String save)throws IOException{
		Map<String, Collection<Alignments>> genes=BEDFileParser.loadAlignmentDataByChr(repeatFile);
		Map<String, Collection<Alignments>> intervals=BEDFileParser.loadAlignmentDataByChr(intervalFile);
		this.info=BEDFileParser.loadDataLine(intervalFile);
		
		Map filtered=new TreeMap();
		for(String chr: genes.keySet()){
			try{
			System.err.println(chr);
			filtered.putAll(filter(genes.get(chr), intervals.get(chr)));
			}catch(Exception ex){}
		}
		
		write(save, filtered);
		
	}
	
	
	private Map filter(Collection<Alignments> knownGenes, Collection<Alignments> intervals){
		Set rtrn=new TreeSet();
		Map map=new TreeMap();
		
		
		for(Alignments interval: intervals){
			boolean overlapsGene=false;
			ArrayList list=new ArrayList();
			map.put(interval, list);
			for(Alignments gene: knownGenes){
				Alignments extendedGene=new Alignments(gene.getChr(), gene.getStart(), gene.getEnd());
				if(interval.overlapsAtAll(extendedGene) ||extendedGene.overlapsAtAll(interval)){
					overlapsGene=true;
					if(map.containsKey(interval)){list=(ArrayList)map.get(interval); list.add(gene);}
					map.put(interval, list);
				}
			}
			if(overlapsGene){rtrn.add(interval);}
		}
		System.err.println(rtrn.size());
		return map;
	}
	
	
	private void write(String save, Map<Alignments, ArrayList> map, Map info)throws IOException{
		FileWriter writer=new FileWriter(save);
		for(Object key: map.keySet()){
			if(map.get(key).size()>0){
				ArrayList<Alignments> list=map.get(key);
				writer.write(key.toString());
				for(Alignments align: list){
					String str=(String)info.get(align);
					writer.write("\t"+str.split("\t")[3]);
				}
				writer.write("\n");
			}
		}
		writer.close();
	}
	
	private void write(String save, Map<Alignments, ArrayList> map)throws IOException{
		FileWriter writer=new FileWriter(save);
		for(Object key: map.keySet()){
			if(map.get(key).size()>0){
				ArrayList<Alignments> list=map.get(key);
				//writer.write(key.toString());
				String str=(String)this.info.get(key);
				writer.write(str+"\n");
				//for(Utils.Alignments align: list){
					//writer.write("\t"+align.toUCSC());
					//String str=(String)info.get(align);
					//writer.write("\t"+str.split("\t")[3]);
				//}
				//writer.write("\n");
			}
		}
		writer.close();
	}
	
	
	public static void main(String[] args)throws IOException{
		//allow for one or the other to be a directory and treat appropriately
		//if both are directories throw and error
		
		if(args.length>0){
			File lincRNA=new File(args[0]);
			File file=new File(args[1]);
			String saveDir=args[2];
			
			if(lincRNA.isDirectory() && file.isDirectory()){System.err.println("Both cant be directories");}
			else if(lincRNA.isDirectory()){
				File[] lincs=lincRNA.listFiles();
				for(int i=0; i<lincs.length; i++){
					System.err.println(lincs[i]);
					String save=saveDir+"/"+lincs[i].getName();
					new OverlappingRegions(lincs[i], file, save);
				}
			}
			else if(file.isDirectory()){
				File[] files=file.listFiles();
				for(int i=0; i<files.length; i++){
					System.err.println(files[i]);
					String save=saveDir+"/"+files[i].getName();
					new OverlappingRegions(lincRNA, files[i], save);
				}
			}
			else{
				new OverlappingRegions(lincRNA, file, saveDir);
			}
			
		}
		
		else{
			System.err.println("Function: "+function+"\n USAGE: \n"+usage);
		}
	}
	
	static String usage="args[0]=File to filter by (lincRNAs) (can be a directory) \nargs[1]=File to filter (ie mRNAs, ESTs)(can be a directory) \nargs[2]=save file (must be a directory if either 0 or 1 are)";
	static String function="Find regions that overlap between the specified files";
	
	
}
