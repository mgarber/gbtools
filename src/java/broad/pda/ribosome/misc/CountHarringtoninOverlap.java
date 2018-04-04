package broad.pda.ribosome.misc;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeSet;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.GeneTools;
import broad.pda.gene.RefSeqGene;

public class CountHarringtoninOverlap {

	public CountHarringtoninOverlap(Map<String, IntervalTree<RefSeqGene>>[] peaks, String save) throws IOException{
		//Find intersection of all
		Collection<RefSeqGene> collapsed=collapseByIntersection(peaks);
		write(save, collapsed);
	}
	
	private void write(String save, Collection<RefSeqGene> collapsed) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(RefSeqGene gene: collapsed){
			writer.write(gene+"\n");
		}
		
		writer.close();
	}
	
	private Collection<RefSeqGene> collapseByIntersection(Map<String, IntervalTree<RefSeqGene>>[] peaks) {
		Collection<RefSeqGene> rtrn=new TreeSet<RefSeqGene>();
		
		for(String chr: peaks[0].keySet()){
			Iterator<Node<RefSeqGene>> iter=peaks[0].get(chr).iterator();
			while(iter.hasNext()){
				RefSeqGene gene=iter.next().getValue();
				
				//check if its in all the others
				boolean allHave=true;
				Collection<RefSeqGene> overlapperSet=new TreeSet<RefSeqGene>();
				for(int i=0; i<peaks.length; i++){
					Iterator<Node<RefSeqGene>> overlappers=peaks[i].get(chr).overlappers(gene.getStart(), gene.getEnd());
					if(!overlappers.hasNext()){allHave=false;}
					else{
						while(overlappers.hasNext()){
							overlapperSet.add(overlappers.next().getValue());
						}
					}
				}
				if(allHave){
					//collapse the overlapper Set
					//Collection<RefSeqGene> collapsed=collapseSet(overlapperSet);
					Iterator<Node<RefSeqGene>> collapsed=GeneTools.merge(overlapperSet).iterator();
					while(collapsed.hasNext()){rtrn.add(collapsed.next().getValue());}
				}
			}
		}
		
		return rtrn;
	}

	

	public static void main(String[] args)throws IOException{
		if(args.length>1){
			File[] files=new File(args[0]).listFiles();
			String save=args[1];
			
			Map<String, IntervalTree<RefSeqGene>>[] maps=new Map[files.length];
			for(int i=0; i<files.length; i++){maps[i]=BEDFileParser.loadDataByChrToTree(files[i]);}
			
			new CountHarringtoninOverlap(maps, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=files \n args[1]=save";
	
}
