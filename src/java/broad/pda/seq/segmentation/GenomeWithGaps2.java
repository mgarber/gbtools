package broad.pda.seq.segmentation;

import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.util.CollapseByIntersection;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;

public class GenomeWithGaps2 {

	private IntervalTree<Alignments> cDNASpace; //indexed by relative position
	private int genomeSize;
	private boolean directOrientation;
	static Logger logger = Logger.getLogger(GenomeWithGaps2.class.getName());
	
	public GenomeWithGaps2(Collection<Alignments> gaps, String chr, int genomeSize, String orientation){
		gaps=CollapseByIntersection.collapseByIntersection(gaps, true);
		cDNASpace=convertToCovered(gaps, chr, genomeSize);
		this.genomeSize = genomeSize;
		this.directOrientation = "+".equals(orientation);
	}
	
	
	//Need to convert this relative space indexes
	private IntervalTree<Alignments> convertToCovered(Collection<Alignments> gaps, String chr, int genomeSize){
		IntervalTree<Alignments> rtrn=new IntervalTree();
		
		Map<String, IntervalTree<Alignments>> trees=makeIntervalTree(gaps);
		IntervalTree<Alignments> gapTree=trees.get(chr);
		
		int currentBP=0;
		if(gapTree != null) {
			Iterator<Node<Alignments>> iter=gapTree.iterator();
			
			int counter=0;
			while(iter.hasNext()){
				Alignments gap=iter.next().getValue();
				Alignments block=new Alignments(gap.getChr(), currentBP, gap.getStart());
				currentBP=gap.getEnd();
				rtrn.put(counter, counter+block.getSize(), block);
				counter=counter+block.getSize();
			}
			
			Alignments end=new Alignments(chr, currentBP, genomeSize);
			if(end.getSize()>0){rtrn.put(counter, counter+end.getSize(), end);}
		} else {
			logger.debug("No gaps tree for chr " + chr);
		}
		return rtrn;
	}
	
	private Map<String, IntervalTree<Alignments>> makeIntervalTree(Collection<Alignments> alignments){
		Map<String, IntervalTree<Alignments>> rtrn=new TreeMap<String, IntervalTree<Alignments>>();
		
		for(Alignments align: alignments){
			IntervalTree<Alignments> tree=new IntervalTree<Alignments>();
			String chr=align.getChr();
			if(rtrn.containsKey(chr)){tree=rtrn.get(chr);}
			tree.put(align.getStart(), align.getEnd(), align);
			rtrn.put(align.getChr(), tree);
		}
		
		return rtrn;
	}
	
	public RefSeqGene getRelativeWindow(int relativeStart, int relativeEnd){
		//System.err.println("Calculating relative window for  " + relativeStart  + " - " + relativeEnd + " region length " + genomeSize);
		if(relativeEnd > genomeSize) {
			System.err.println("ERROR " + relativeStart  + " - " + relativeEnd + " have no overlappers in this CDNA of size " + (genomeSize));
			return null;  //TODO: Should throw an illegal argument exception since this means that coordinates are outside of range
		}
	
		Iterator<Node<Alignments>> relativeMapping=cDNASpace.overlappers(relativeStart, relativeEnd);
		Collection<Alignments> exons=new TreeSet<Alignments>();
		//if(!relativeMapping.hasNext()) { System.out.println("start-end: "+ relativeStart+"-"+relativeEnd+" had no exons mapping");}
		while(relativeMapping.hasNext()){
			Node<Alignments> node=relativeMapping.next();
			
			//if rs and re in 1 block then done. Might be able to tell this case from num overlappers
			if(node.getStart()<relativeStart && node.getEnd()>relativeEnd){
				Alignments absAlignment=node.getValue();
				int absStart=relativeStart-node.getStart();
				int absEnd=relativeEnd-node.getStart();
				Alignments firstExon=new Alignments(absAlignment.getChr(), absAlignment.getStart()+absStart, absAlignment.getStart()+absEnd); //end unless fully contained
				exons.add(firstExon);
			}
			
			//if rs is within the block start at the appropriate point
			else if(node.getStart()<relativeStart){
				Alignments absAlignment=node.getValue();
				int absStart=relativeStart-node.getStart();
				Alignments firstExon=new Alignments(absAlignment.getChr(), absAlignment.getStart()+absStart, absAlignment.getEnd()); //end unless fully contained
				exons.add(firstExon);
			}
			
			else if(node.getEnd()>relativeEnd){
				Alignments absAlignment=node.getValue();
				int absEnd=relativeEnd-node.getStart();
				Alignments lastExon=new Alignments(absAlignment.getChr(), absAlignment.getStart(), absAlignment.getStart()+absEnd); //end unless fully contained
				exons.add(lastExon);
			}
			else{
				Alignments absAlignment=node.getValue();
				exons.add(absAlignment);
			}
		}
		
		RefSeqGene gene=null;
		if(!exons.isEmpty()){gene=new RefSeqGene(exons);}
		
		return gene;
	}
	
	public int getRelativeGenomeLength(){
		return this.cDNASpace.max().getEnd();
	}
	
}
