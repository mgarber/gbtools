package broad.pda.seq.utils;

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

//Goal: For genes with too many isoforms filter them down

public class FilterIsoforms {
	
	int maxNumIsoforms=16;

	public FilterIsoforms(File file, String save)throws IOException{
		Map<Alignments, Collection<RefSeqGene>> isoforms=parse(file);
		write(save, isoforms, maxNumIsoforms);
	}
	
	private Map<Alignments, Collection<RefSeqGene>> parse(File file) throws IOException {
		Map<Alignments, Collection<RefSeqGene>> rtrn=new TreeMap();
		
		Collection<RefSeqGene> genes=BEDFileParser.loadData(file);
		for(RefSeqGene gene: genes){
			Alignments align=gene.getAlignment();
			Collection<RefSeqGene> list=new TreeSet();
			if(rtrn.containsKey(align)){list=rtrn.get(align);}
			list.add(gene);
			rtrn.put(align, list);
		}
		
		return rtrn;
	}

	private void write(String save, Map<Alignments, Collection<RefSeqGene>> isoforms, int maxNumIsoforms2)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(Alignments align: isoforms.keySet()){
			Collection<RefSeqGene> list=isoforms.get(align);
			int counter=0;
			if(list.size()>maxNumIsoforms2){System.out.println(align+" "+list.size());}
			for(RefSeqGene gene: list){
				if(counter<maxNumIsoforms2){
					writer.write(gene+"\n");
				}
				counter++;
			}
			
		}
		
		writer.close();
	}

	public static void main(String[] args)throws IOException{
		File file=new File(args[0]);
		String save=args[1];
		new FilterIsoforms(file, save);
	}
}
