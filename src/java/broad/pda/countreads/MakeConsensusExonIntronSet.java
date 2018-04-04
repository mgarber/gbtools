package broad.pda.countreads;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.util.CollapseByIntersection;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;

public class MakeConsensusExonIntronSet {

	public MakeConsensusExonIntronSet(Collection<RefSeqGene> genes, String save)throws IOException{
		
		Map<String, Collection<Alignments>> regionsByCategory=categorize(genes);
		write(save, regionsByCategory);
	
	}
	
	private void write(String save, Map<String, Collection<Alignments>> regionsByCategory)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(String category: regionsByCategory.keySet()){
			Collection<Alignments> regions=regionsByCategory.get(category);
			for(Alignments region: regions){writer.write(region+"\t"+category+"\n");}
		}
		
		writer.close();
	}
	
	private static Map<String, Collection<Alignments>> categorize(Collection<RefSeqGene> genes){
		Collection<Alignments> exons=new TreeSet();
		Collection<Alignments> introns=new TreeSet();
		Collection<Alignments> intergenic=new TreeSet();
		Collection<Alignments> tu=new TreeSet();
		
		for(RefSeqGene gene: genes){
			exons.addAll(gene.getExonSet());
			introns.addAll(gene.getIntronSet());
			tu.add(gene.getAlignment());
		}
		
		//In order to make each set unique we need to rank priority each group
		
		//Top rank is exons
		exons=CollapseByIntersection.collapseByIntersection(exons, false); //if any base is exonic call it exon --> Collapse by union
		
		
		//Followed by Introns
		introns=CollapseByIntersection.collapseByIntersection(introns, true); //if all bases are intronic then call it intron --> Collapse by intersection
		
		
		//Followed by intergenic space
		tu=CollapseByIntersection.collapseByIntersection(tu, false); //if any base is in unit then part of TU --> Collapse by union
		intergenic=defineIntergenic(tu);
		
		Map<String, Collection<Alignments>> rtrn=new TreeMap();
		
		rtrn.put("Exons", exons);
		rtrn.put("Introns", introns);
		rtrn.put("Intergenic", intergenic);
		
		return rtrn;
	}
	
	//Assume that the set is sorted!!
	private static Collection<Alignments> defineIntergenic(Collection<Alignments> tus){
		Collection<Alignments> intergenic=new TreeSet();
		
		Object[] array=tus.toArray();
		for(int i=0; i<array.length-1; i++){
			Alignments current=(Alignments)array[i];
			Alignments next=(Alignments)array[i+1];
			if(current.getChr().equalsIgnoreCase(next.getChr())){Alignments intergenicRegion=new Alignments(current.getChr(), current.getEnd(), next.getStart()); intergenic.add(intergenicRegion);}
		}
		
		return intergenic;
	}
	
	
	public static void main(String[] args)throws IOException{
		if(args.length>1){
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[0]));
			String save=args[1];
			new MakeConsensusExonIntronSet(genes, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=genes (Full BED) \n args[1]=save";
}
