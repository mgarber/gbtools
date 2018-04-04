package broad.pda.seq.pairedend;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.broad.igv.sam.Alignment;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;
import net.sf.samtools.util.CloseableIterator;

public class ConnectWithPairedEnds {
	
	private int maxDistance=1000000;

	public ConnectWithPairedEnds(Collection<RefSeqGene> genes, AlignmentDataModel data, String save)throws IOException{
		Map<String, IntervalTree<RefSeqGene>> trees=this.makeIntervalTrees(genes);
		Collection<PairedEndAlignment> pairMatches=makePairs(trees, data);
		write(save, pairMatches);
	}
	
	/*private Collection<PairedEndAlignment> makePairs(Map<String, IntervalTree<RefSeqGene>> trees, AlignmentDataModel data){
		Collection<PairedEndAlignment> rtrn=new TreeSet();
		
		for(String chr: trees.keySet()){
			System.err.println(chr);
			IntervalTree<Alignment> chunkTree=data.getIntervalTree(chr, 0, data.getChromosomeLengths().get(chr));
			Iterator<Node<Alignment>> iter=chunkTree.iterator();
			int i=0;
			while(iter.hasNext()){
				Alignment pair=iter.next().getValue();
				Iterator<Node<RefSeqGene>> start=trees.get(pair.getChromosome()).overlappers(pair.getAlignmentStart()-1, pair.getAlignmentStart()+1);
				Iterator<Node<RefSeqGene>> end=trees.get(pair.getChromosome()).overlappers(pair.getAlignmentEnd()-1, pair.getAlignmentEnd()+1);
				//get overlapping with the ends
				while(start.hasNext()){
					while(end.hasNext()){
						RefSeqGene left=start.next().getValue();
						RefSeqGene right=end.next().getValue();
						PairedEndAlignment pairAlign=new PairedEndAlignment(left, right, left.getAlignment().toUCSC()+"_"+right.getAlignment().toUCSC());
						boolean isGood=isGood(pairAlign);
						if(isGood){rtrn.add(pairAlign);}
					}
				}
				if(i%100 ==0){System.err.println(i+" "+pair.getChromosome()+":"+pair.getAlignmentStart()+"-"+pair.getAlignmentEnd());}
				i++;
			}
		}
		
		
		
		return rtrn;
	}*/
	
	private Collection<PairedEndAlignment> makePairs(Map<String, IntervalTree<RefSeqGene>> trees, AlignmentDataModel data) throws IOException{
		Collection<PairedEndAlignment> rtrn=new TreeSet();
		
		for(String chr: trees.keySet()){
			//System.err.println(chr);
			Iterator<RefSeqGene> iter=trees.get(chr).valueIterator();
			while(iter.hasNext()){
				RefSeqGene gene=iter.next();
				System.err.println(gene.getAlignment().toUCSC());
				CloseableIterator<Alignment> pairs=data.getAlignmentsOverlappingRegion(gene.getAlignment());
				//get overlapping paired ends and get the other gene that overlaps it
				while(pairs.hasNext()){
					Alignment pair=pairs.next();
					Collection<PairedEndAlignment> newPairs=makePairs(pair, trees.get(chr));
					rtrn.addAll(newPairs);
				}
			}
		}
		
		
		
		return rtrn;
	}
	
	private Collection<PairedEndAlignment> makePairs(Alignment pairAlign, IntervalTree<RefSeqGene> tree){
		Collection<PairedEndAlignment> rtrn=new TreeSet();
		
		Iterator<Node<RefSeqGene>> leftNodes=tree.overlappers(pairAlign.getAlignmentStart()-1, pairAlign.getAlignmentStart()+1);
		
		
		
		while(leftNodes.hasNext()){
			RefSeqGene left=leftNodes.next().getValue();
			Iterator<Node<RefSeqGene>> rightNodes=tree.overlappers(pairAlign.getAlignmentEnd()-1, pairAlign.getAlignmentEnd()+1);
			while(rightNodes.hasNext()){
				RefSeqGene right=rightNodes.next().getValue();
				PairedEndAlignment pair=new PairedEndAlignment(left, right, "temp");
				if(isGood(pair)){rtrn.add(pair);}
			}
		}
		
		return rtrn;
	}
	
	private Collection<PairedEndAlignment> makePairs(RefSeqGene gene, Iterator<Node<RefSeqGene>> overlappers){
		Collection<PairedEndAlignment> rtrn=new TreeSet();
		
		while(overlappers.hasNext()){
			RefSeqGene left=overlappers.next().getValue();
			PairedEndAlignment pair=new PairedEndAlignment(left, gene, "temp");
			if(isGood(pair)){rtrn.add(pair);}
		}
		
		return rtrn;
	}
	
	private boolean isGood(PairedEndAlignment pair){
		//if left and right dont overlap
		if(pair.onSameChromosome() && pair.getDistance()<maxDistance && !pair.pairOverlaps()){return true;}
		return false;
	}
	
	private void write(String save, Collection<PairedEndAlignment> pairs)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(PairedEndAlignment pair: pairs){
			writer.write(pair.toString()+"\n");
		}
		
		writer.close();
	}
	
	private Map<String, IntervalTree<RefSeqGene>> makeIntervalTrees(Collection<RefSeqGene> alignments){
		Map<String, IntervalTree<RefSeqGene>> rtrn=new TreeMap();
		
		for(RefSeqGene align: alignments){
			IntervalTree<RefSeqGene> tree=new IntervalTree();
			String chr=align.getChr();
			if(rtrn.containsKey(chr)){tree=rtrn.get(chr);}
			tree.put(align.getAlignment().getStart(), align.getAlignment().getEnd(), align);
			rtrn.put(align.getChr(), tree);
		}
		
		return rtrn;
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>4){
		Map<String, Collection<RefSeqGene>> genes=BEDFileParser.loadDataByChr(new File(args[0]));
		AlignmentDataModel pairs=new GenericAlignmentDataModel(args[1], args[2]);
		String save=args[3];
		String chr=args[4];
		new ConnectWithPairedEnds(genes.get(chr), pairs, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=genes (Full BED) \n args[1]=pairs (BED File) \n args[2]=sizes \n args[3]=save file \n args[4]=chr";
	
}
