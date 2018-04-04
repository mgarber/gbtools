package broad.core.overlaputils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;

import broad.core.datastructures.IntervalTree;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;

public class FilterGenesOverlappingRegions {

	public FilterGenesOverlappingRegions(Map<String, Collection<Alignments>> genesToFilter, Map<String, IntervalTree<Alignments>> genesToFilterBy, String save, boolean exclude)throws IOException{
		
		Collection<Alignments> filtered=new ArrayList();
		
		for(String chr: genesToFilter.keySet()){
			System.err.println(chr);
			Collection<Alignments> lincs=genesToFilter.get(chr);
			if(genesToFilterBy.containsKey(chr)){
				IntervalTree<Alignments> proteins=genesToFilterBy.get(chr);
				filtered.addAll(filter(lincs, proteins, exclude));
			}
		}
		
		write(save, filtered);
	}
	
	
	private Collection<Alignments> filter(Collection<Alignments> lincs, IntervalTree<Alignments> proteins, boolean exclude){
		Collection<Alignments> rtrn=new ArrayList();
		
		for(Alignments linc: lincs){
			boolean overlaps=proteins.overlappers(linc.getStart(), linc.getEnd()).hasNext();
			if(exclude && !overlaps){rtrn.add(linc);}
			else if(!exclude && overlaps){rtrn.add(linc);}
		}
		
		return rtrn;
	}
	
	private void write(String save, Collection<Alignments> genes)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(Alignments gene: genes){
			writer.write(gene+"\n");
		}
		
		writer.close();
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>3){
			System.err.println("started");
			Map<String, Collection<Alignments>> genesToFilter=BEDFileParser.loadAlignmentDataByChr(new File(args[0]));
			System.err.println("loaded lincRNAs");
			Map<String, IntervalTree<Alignments>> genesToFilterBy=BEDFileParser.loadAlignmentDataToTree(new File(args[1]));
			System.err.println("loaded proteins into interval tree");
			String save=args[2];
			boolean exclude=new Boolean(args[3]);
			new FilterGenesOverlappingRegions(genesToFilter, genesToFilterBy, save, exclude);
			
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=genes to filter (Full BED- lincRNA transcripts) \n args[1]=genes to filter by (protein coding genes) \n args[2]=save \n args[3]=exclude/include based on overlap (true=exclude, false=include)";
	
}
