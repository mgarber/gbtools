package broad.pda.rnai;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.gene.RefSeqGeneWithIsoforms;

public class GetGTFForNewLincs {

	public GetGTFForNewLincs(File file, Map<String, IntervalTree<RefSeqGeneWithIsoforms>> scriptureTree, Map<String, IntervalTree<RefSeqGene>> codingTree, String save) throws IOException{
	
		Map<String, Collection<RefSeqGeneWithIsoforms>> genesByName=new TreeMap<String, Collection<RefSeqGeneWithIsoforms>>();
		
		Map<String, Alignments> namesAndRegion=parseNamesAndRegion(file);
		for(String name: namesAndRegion.keySet()){
			Alignments region=namesAndRegion.get(name);
			//get reconstructions
			Iterator<Node<RefSeqGeneWithIsoforms>> iter=scriptureTree.get(region.getChr()).overlappers(region.getStart(), region.getEnd());
			//check if overlaps protein
			Collection<RefSeqGeneWithIsoforms> set=filterNonOverlapping(iter, codingTree);
			if(genesByName.containsKey(name)){
				set.addAll(genesByName.get(genesByName.get(name)));
			}
			genesByName.put(name, set);
		}
		
		write(save, genesByName);
	}
	
	
	private Map<String, Alignments> parseNamesAndRegion(File file) throws IOException {
		Map<String, Alignments> rtrn=new TreeMap<String, Alignments>();
		Collection<String> list=BEDFileParser.loadList(file.getAbsolutePath());
		
		for(String line: list){
			String[] tokens=line.split("\t");
			String name=tokens[0];
			Alignments align=new Alignments(tokens[1], tokens[2], tokens[3]);
			rtrn.put(name, align);
		}
		
		return rtrn;
	}


	private void write(String save,	Map<String, Collection<RefSeqGeneWithIsoforms>> genesByName) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String name: genesByName.keySet()){
			int counter=1;
			Collection<RefSeqGeneWithIsoforms> genes=genesByName.get(name);
			boolean hasMultiExonic=hasMultiExonic(genes);
			for(RefSeqGeneWithIsoforms gene: genes){
				if((hasMultiExonic && gene.getNumExons()>1) || (!hasMultiExonic)){
					Collection<RefSeqGene> isoforms=gene.getAllIsoforms();
					for(RefSeqGene isoforom: isoforms){
						String isoformName=name+".2_"+counter;
						isoforom.setName(isoformName);
						writer.write(isoforom+"\n");
						counter++;
					}
				}
			}
		}
		
		writer.close();
	}


	private boolean hasMultiExonic(Collection<RefSeqGeneWithIsoforms> genes) {
		for(RefSeqGeneWithIsoforms gene: genes){
			if(gene.getNumExons()>1){return true;}
		}
		return false;
	}


	private Collection<RefSeqGeneWithIsoforms> filterNonOverlapping(Iterator<Node<RefSeqGeneWithIsoforms>> iter, Map<String, IntervalTree<RefSeqGene>> codingTree) {
		Collection<RefSeqGeneWithIsoforms> rtrn=new TreeSet<RefSeqGeneWithIsoforms>();
		
		while(iter.hasNext()){
			RefSeqGeneWithIsoforms transcript=iter.next().getValue();
			Iterator<Node<RefSeqGene>> overlappers=codingTree.get(transcript.getChr()).overlappers(transcript.getStart(), transcript.getEnd());
			if(!overlappers.hasNext()){rtrn.add(transcript);}
		}
		
		return rtrn;
	}


	public static void main(String[] args) throws IOException{
		if(args.length>3){
			File RNAiReport=new File(args[0]);
			BEDFileParser parser1=new BEDFileParser(args[1]);
			Map<String, IntervalTree<RefSeqGeneWithIsoforms>> scriptureTree=parser1.getIntervalTreeWithIsoforoms();
			Map<String, IntervalTree<RefSeqGene>> codingTree=BEDFileParser.loadDataByChrToTree(new File(args[2]));
			String save=args[3];
			new GetGTFForNewLincs(RNAiReport, scriptureTree, codingTree, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=RNAi report \n args[1]=scripture genes \n args[2]= RefSeq coding genes \n args[3]=save";
	
}
