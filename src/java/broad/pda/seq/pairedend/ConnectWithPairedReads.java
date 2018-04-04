package broad.pda.seq.pairedend;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;



public class ConnectWithPairedReads {

	//take 2 cufflinks assemblies and link with paired info
	public ConnectWithPairedReads(Set<RefSeqGene> set1, Set<RefSeqGene> set2, Collection<PairedEndAlignment> pairs, String save)throws IOException{
		//get connected assemblies
		Map<PairedEndAlignment, Integer> connections=this.connect(set1, set2, pairs);
		
		write(save, connections);
	}
	
	
	private void write(String save, Map<PairedEndAlignment, Integer> connections)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(PairedEndAlignment pair: connections.keySet()){
			int count=connections.get(pair);
			writer.write(pair.toString()+"\t"+pair.getCDNAInsertSize()+"\t"+count+"\n");
		}
		
		writer.close();
	}
	
	private Map<PairedEndAlignment, Integer> connect(Set<RefSeqGene> set1, Set<RefSeqGene> set2, Collection<PairedEndAlignment> pairs){
		Map<PairedEndAlignment, Integer> rtrn=new TreeMap();
		
		for(PairedEndAlignment pair: pairs){
			Set<RefSeqGene> genes1=this.getOverlapping(set1, pair);
			Set<RefSeqGene> genes2=this.getOverlapping(set2, pair);
			Set<PairedEndAlignment> possiblePairs=enumerateAllPossiblePairs(genes1, genes2, pair.getName());
			rtrn=putAll(rtrn, possiblePairs);
			//enumerate all possibilites
			//keep track of inferred insert size
			//keep track of "genomic" insert size
		}
		
		return rtrn;
	}
	
	private Map putAll(Map<PairedEndAlignment, Integer> map, Set<PairedEndAlignment> possiblePairs){
		Map rtrn=map;
		
		for(PairedEndAlignment pair: possiblePairs){
			int count=0;
			if(rtrn.containsKey(pair)){count=map.get(pair);}
			count++;
			rtrn.put(pair, count);
		}
		
		return rtrn;
	}
	
	private Set<PairedEndAlignment> enumerateAllPossiblePairs(Set<RefSeqGene> genes1, Set<RefSeqGene> genes2, String name){
		Set<PairedEndAlignment> rtrn=new TreeSet();
		
		for(RefSeqGene gene1: genes1){
			for(RefSeqGene gene2: genes2){
				PairedEndAlignment pair=new PairedEndAlignment(gene1, gene2, name);
				rtrn.add(pair);
			}
		}
		
		return rtrn;
	}
	
	private Set<RefSeqGene> getOverlapping(Set<RefSeqGene> set1, PairedEndAlignment pair){
		Set<RefSeqGene> rtrn=new TreeSet();
		for(RefSeqGene gene: set1){
			if(gene.getAlignment().overlapsAtAll(pair.getEncompassingAlignment())){
				rtrn.add(gene);
			}
		}
		return rtrn;
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>5){
			Set<RefSeqGene> genes1=BEDFileParser.loadData(new File(args[0]));
			Set<RefSeqGene> genes2=BEDFileParser.loadData(new File(args[1]));
			File leftFile=new File(args[2]);
			File rightFile=new File(args[3]);
			String chr=args[4];
			
			Collection<PairedEndAlignment> pairs=PairedEndMapping.getAllPairs(leftFile, rightFile, chr);
			String save=args[5];
			new ConnectWithPairedReads(genes1, genes2, pairs, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=cufflinks1 \n args[1]=cufflinks2 \n args[2]=left reads (SAM file) \n args[3]=right reads (SAM file) \n args[4]=chr \n args[5]=save";
	
}
