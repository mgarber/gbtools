package broad.pda.rnaseq.misc;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.util.CollapseByIntersection;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;

public class FindK4Overlapping5Prime {

	public FindK4Overlapping5Prime(Collection<Alignments> k4Regions, Collection<RefSeqGene> predictions, String save, int extension)throws IOException{
		Map<String, IntervalTree<Alignments>> k4Tree=CollapseByIntersection.makeIntervalTree(k4Regions);
		FileWriter writer=new FileWriter(save);
		
		for(RefSeqGene gene: predictions){
			Alignments firstExon=gene.get5PrimeExon();
			firstExon=new Alignments(firstExon.getChr(), firstExon.getStart()-extension, firstExon.getEnd()+extension);
			Iterator<Node<Alignments>> k4s=k4Tree.get(firstExon.getChr()).overlappers(firstExon.getStart(), firstExon.getEnd());
			write(writer, gene, k4s);
		}
		writer.close();
	}

	private void write(FileWriter writer, RefSeqGene gene, Iterator<Node<Alignments>> k4s) throws IOException {
		
		while(k4s.hasNext()){
			Alignments k4=k4s.next().getValue();
			writer.write(k4+"\t"+gene.getName()+"\n");
		}
		
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>3){
			Collection<Alignments> k4Scores=BEDFileParser.loadAlignmentData(new File(args[0]));
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[1]));
			String save=args[2];
			int extension=new Integer(args[3]);
			new FindK4Overlapping5Prime(k4Scores, genes, save, extension);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]= K4 regions \n args[1]=genes \n args[2]=save \n args[3]=extension";
}
