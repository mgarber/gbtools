package broad.pda.rnaseq.misc;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.graph.ChromosomeWithBubbles2;
import broad.pda.seq.graph.Path;

public class FilterIsoformsByRemovingSpuriousExons {

	public FilterIsoformsByRemovingSpuriousExons(Collection<RefSeqGene> genes, Collection<Alignments> badExons, String save) throws IOException{
		Collection<Alignments> exons=new TreeSet<Alignments>();
		Collection<Alignments> introns=new TreeSet<Alignments>();
		
		for(RefSeqGene gene: genes){
			exons.addAll(gene.getExonSet());
			introns.addAll(gene.getIntronSet());
		}
		
		for(Alignments badExon: badExons){
			if(!exons.contains(badExon)){System.err.println("Cant find "+badExon.toUCSC());}
			exons.remove(badExon);
		}
		
		Map<String, Collection<Alignments>> exonsByChr=splitByChr(exons);
		Map<String, Collection<Alignments>> intronsByChr=splitByChr(introns);
		
		
		FileWriter writer=new FileWriter(save);
		for(String chr: exonsByChr.keySet()){
			ChromosomeWithBubbles2 graph=new ChromosomeWithBubbles2(chr, exonsByChr.get(chr), intronsByChr.get(chr), null, null,0,0,0,0);
			Collection<Path> paths=graph.getAllPaths();
			for(Path path: paths){
				writer.write(path.toGene()+"\n");
			}
		}
		writer.close();
		
		writer=new FileWriter(save+".exons.bed");
		for(Alignments exon: exons){
			writer.write(exon+"\n");
		}
		writer.close();
		
		writer=new FileWriter(save+".introns.bed");
		for(Alignments intron: introns){
			writer.write(intron+"\n");
		}
		writer.close();
		
	}
	
	private Map<String, Collection<Alignments>> splitByChr(Collection<Alignments> exons) {
		Map<String, Collection<Alignments>> rtrn=new TreeMap<String, Collection<Alignments>>();
		
		for(Alignments exon: exons){
			Collection<Alignments> list=new TreeSet<Alignments>();
			if(rtrn.containsKey(exon.getChr())){
				list=rtrn.get(exon.getChr());
			}
			list.add(exon);
			rtrn.put(exon.getChr(), list);
		}
		
		return rtrn;
	}

	public static void main(String[] args)throws IOException{
		if(args.length>2){
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[0]));
			Collection<Alignments> badExons=BEDFileParser.loadAlignmentData(new File(args[1]));
			String save=args[2];
			new FilterIsoformsByRemovingSpuriousExons(genes, badExons, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=full bed \n args[1]=bad exons \n args[2]=save";
}
