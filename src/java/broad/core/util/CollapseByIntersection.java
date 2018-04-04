package broad.core.util;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;


public class CollapseByIntersection {
	
	static int numOverlap;
	static int numOverlapGenes;
	private static final int MAX_RECURSION = 50;
	
	public static Set<Alignments> collapseByIntersection(Collection<Alignments> alignments, boolean intersection){
		numOverlap=10; //just to start
		
		Set<Alignments> rtrn=new TreeSet<Alignments>(alignments);
		
		int i=0;
		while(numOverlap>1){ //TODO: do not use a static variable for this. 
			numOverlap=0;
			Set<Alignments> set=new TreeSet<Alignments>();
			Map<String, IntervalTree<Alignments>> trees=makeIntervalTree(rtrn);
			//for each region get all overlappign and collapse by intersection
			for(Alignments align: rtrn){
				Iterator<Node<Alignments>> regions=getOverlapping(align, trees);
				Alignments collapsed=collapse(regions, intersection);
				if(collapsed.getSize()>0){set.add(collapsed);}
			}
			i++;
			rtrn=set;
			//System.err.println("Iteration "+i+" Num Overlap "+numOverlap);
		}
		
		return rtrn;
	}
	
	/**
	 * This method returns only exons that are trimmed down by intron location. If its fully within intron it just returns it. If bases overlap intron it trims down
	 * @param exons
	 * @param introns
	 * @param strand
	 * @return
	 */
	public static Collection<Alignments> DecollapseByIntronLocation(Collection<Alignments> exons, Collection<Alignments> introns, String strand){
		if(introns==null || introns.isEmpty()){return exons;}
		
		introns=filter(introns, strand);
		
		TreeSet<Alignments> decollapsedExons=new TreeSet<Alignments>();
		
		Map<String, IntervalTree<Alignments>> intronTree=makeIntervalTree(introns);
		
		for(Alignments exon: exons){			
			//There is a reason for the -1 and +1. If an intron abuts but does not overlap the exon overlaps an intron that truly overlaps the exon,
			//the result of the decollapse will not contain the segment that abuts the original exon. The -1 and +1 ensure that butting introns are
			//included and thus exons that abut them are reported.
			Iterator<Node<Alignments>> iter=intronTree.get(exon.getChr()).overlappers(exon.getStart() - 1, exon.getEnd() + 1); 
			if(!iter.hasNext()) {
				decollapsedExons.add(exon);
			} else {
				LinkedList<Alignments> overlappingIntrons = new LinkedList<Alignments>(); 
				while(iter.hasNext()){
					overlappingIntrons.add(iter.next().getValue());
				}
				//Collections.sort(overlappingIntrons);
				//TreeSet<Alignments> diff = new TreeSet<Alignments>();
				//getRecursiveDifference(exon, overlappingIntrons, rtrn, 1);
				carveExon(exon, overlappingIntrons, decollapsedExons);
				//rtrn.addAll(diff);				
			}
		}
		
		// We are now going to filter out any exon that fully contains any exon-intron-exon combination. These affect highly expressed genes producing large numbers of intronic reads.
		Map<String, IntervalTree<Alignments>> exonTree = makeIntervalTree(decollapsedExons);
		Set<Alignments> filteredExons = new TreeSet<Alignments>();
		
		for(Alignments exon : decollapsedExons) {
			Iterator<Node<Alignments>> intronIter=intronTree.get(exon.getChr()).overlappers(exon.getStart(), exon.getEnd());
			boolean keep = true;
			while(intronIter.hasNext() && keep) {
				Alignments intron = intronIter.next().getValue();
				if(exon.fullyContained(intron)) {
					Alignments intronLeftPoint = new Alignments(intron.getChr(), intron.getStart()-1, intron.getStart());
					Alignments intronRightPoint = new Alignments(intron.getChr(), intron.getEnd(), intron.getEnd()+1);						
					Iterator<Node<Alignments>> leftExons=exonTree.get(exon.getChr()).overlappers(intronLeftPoint.getStart(), intronLeftPoint.getEnd());
					Iterator<Node<Alignments>> rightExons=exonTree.get(exon.getChr()).overlappers(intronRightPoint.getStart(), intronRightPoint.getEnd());
					while(leftExons.hasNext() && keep) {
						Alignments leftExon = leftExons.next().getValue();
						if(leftExon.equals(exon)) {
							continue;
						}
						while(rightExons.hasNext() && keep) {
							Alignments rightExon = rightExons.next().getValue();
							if(rightExon.equals(exon)) {
								continue;
							}
							keep = exon.fullyContained(rightExon) && exon.fullyContained(leftExon);
							
						}
					}
				}
			}
	
			if(keep) {
				filteredExons.add(exon);
			}
		}
		
		return filteredExons;
	}
	
	public static Collection<Alignments> filter(Collection<Alignments> introns, String strand) {
		if(strand.equalsIgnoreCase("+")){
			Collection<Alignments> rtrn=new TreeSet<Alignments>();
			for(Alignments intron: introns){if(intron.getOrientation().equalsIgnoreCase("+")){rtrn.add(intron);}}
			return rtrn;
		}
		else if(strand.equalsIgnoreCase("-")){
			Collection<Alignments> rtrn=new TreeSet<Alignments>();
			for(Alignments intron: introns){if(intron.getOrientation().equalsIgnoreCase("-")){rtrn.add(intron);}}
			return rtrn;
		}
		else{return introns;}
	}

	public static Collection<Alignments> DecollapseByIntronLocation(Collection<Alignments> exons, Collection<Alignments> introns){
		return DecollapseByIntronLocation(exons, introns, "*");
	}
	
	
	private static void carveExon(Alignments exon, LinkedList<Alignments> overlappingIntrons, TreeSet<Alignments> difference) {
		if( difference.contains(exon)) {
			return;
		}
		//System.err.println("Carving " + exon.toUCSC());
		if(overlappingIntrons.size() == 0 ) {
			difference.add(exon);
		} else {
			TreeSet<LightweightGenomicAnnotation> intronLeftEnds = new TreeSet<LightweightGenomicAnnotation>();
			TreeSet<LightweightGenomicAnnotation> intronRightEnds = new TreeSet<LightweightGenomicAnnotation>();
			LightweightGenomicAnnotation enlargedExon = new BasicLightweightAnnotation(exon.getChr(), exon.getStart()-1,exon.getEnd()+1);
		
			for(Alignments intron : overlappingIntrons){
				/*if(exon.fullyContained(intron)){
					LightweightGenomicAnnotation leftEnd = new BasicLightweightAnnotation(intron.getChr(), intron.getStart()-1, intron.getStart());
					LightweightGenomicAnnotation rightEnd = new BasicLightweightAnnotation(intron.getChr(), intron.getEnd(), intron.getEnd()+1);
					intronLeftEnds.add(leftEnd);
					intronRightEnds.add(rightEnd);
				} else {
					LightweightGenomicAnnotation leftEnd = new BasicLightweightAnnotation(intron.getChr(), intron.getStart(), intron.getStart()+1);
					LightweightGenomicAnnotation rightEnd = new BasicLightweightAnnotation(intron.getChr(), intron.getEnd() - 1, intron.getEnd());
					if(enlargedExon.overlaps(leftEnd)) {
						intronLeftEnds.add(leftEnd);
					}
					if(enlargedExon.overlaps(rightEnd)) {
						intronRightEnds.add(rightEnd);
					} 
				}*/
				
				LightweightGenomicAnnotation leftEnd = new BasicLightweightAnnotation(intron.getChr(), intron.getStart(), intron.getStart()+1);
				LightweightGenomicAnnotation rightEnd = new BasicLightweightAnnotation(intron.getChr(), intron.getEnd() - 1, intron.getEnd());
				if(enlargedExon.overlaps(leftEnd)) {
					intronLeftEnds.add(leftEnd);
				}
				if(enlargedExon.overlaps(rightEnd)) {
					intronRightEnds.add(rightEnd);
				}
			}
				
			
			if(intronRightEnds.isEmpty()) {
				intronRightEnds.add(new BasicLightweightAnnotation(exon.getChr(), exon.getStart()-1, exon.getStart()));
			}
			if(intronLeftEnds.isEmpty()) {
				intronLeftEnds.add(new BasicLightweightAnnotation(exon.getChr(), exon.getEnd(), exon.getEnd()+1));
			}

			for(LightweightGenomicAnnotation re : intronRightEnds) {
				for(LightweightGenomicAnnotation le : intronLeftEnds) {
					if(re.getEnd() < le.getStart()) {
						Alignments carving = new Alignments(exon.getChr(), re.getEnd(),le.getStart());
						difference.add(carving);
					} else {
						// carve out the middle of the exon
						if(le.getEnd() > exon.getStart()) {
							Alignments leftPiece = new Alignments (exon.getChr(), exon.getStart(), Math.min(exon.getEnd(),le.getStart()));//new Alignments (exon.getChr(), exon.getStart(), Math.min(exon.getEnd(),le.getEnd()+1));
							difference.add(leftPiece);
						}
						if(re.getStart() < exon.getEnd()) {
							Alignments rightPiece = new Alignments (exon.getChr(), re.getEnd(), exon.getEnd());//new Alignments (exon.getChr(), re.getStart()+1, exon.getEnd());
							difference.add(rightPiece);
						}
					}
				}
			}
			
		}
		

	}
	
	private static void getRecursiveDifference(Alignments exon, LinkedList<Alignments> overlappingIntrons, TreeSet<Alignments> difference, int recursionNum) {
		//TreeSet<Alignments> rtrn = new TreeSet<Alignments>();
		//System.err.println("Recursion: " + recursionNum + ", Exon " + exon.toUCSC() + ", #introns: " + overlappingIntrons.size());
		Alignments lastIntron = null;
		if( difference.contains(exon)) {
			return;
		}
		if(overlappingIntrons.size() == 0 || recursionNum >= MAX_RECURSION ) {
			difference.add(exon);
			if(recursionNum >= MAX_RECURSION) {
				//System.err.println("Maximum number of recursions hit (exon "+exon.toUCSC() +"). Graph will get too complex, this is usually the result of alingment artifacts.");
			}
		} else {
			LinkedList<Alignments> usedIntrons = new LinkedList<Alignments>();
			while(overlappingIntrons.size() > 0) {
				Alignments intron = overlappingIntrons.pop();
				//System.err.println("\tIntron " + intron.toUCSC());
				if(lastIntron != null && !intron.overlaps(lastIntron)) {
					break;
				}
				Collection<Alignments> exonIntronDiff = takeDifference(exon, intron);
				//System.err.println("\texonIntrondiff = " + exonIntronDiff);
				
				for(Alignments exonPart : exonIntronDiff) {
					LinkedList<Alignments> nonOverlappingIntrons = new LinkedList<Alignments>();
					for(Alignments otherIntron : overlappingIntrons) {
						if(!otherIntron.overlaps(intron) && otherIntron.overlaps(exonPart)) {
							nonOverlappingIntrons.add(otherIntron);
						}
					}
					for(Alignments usedIntron : usedIntrons) {
						if(!usedIntron.overlaps(intron) && intron.overlaps(exonPart)) {
							nonOverlappingIntrons.add(usedIntron);
						}
					}
					Collections.sort(nonOverlappingIntrons);
					//System.err.println("\t\tnonoverlapping introns" + nonOverlappingIntrons);
					if(!difference.contains(exonPart)){
						getRecursiveDifference(exonPart, nonOverlappingIntrons, difference , recursionNum++);
					}
				}
				lastIntron = intron;
				usedIntrons.push(intron);
			}
		}
		//System.err.println("Done with " + exon.toUCSC());
		
		//return rtrn;
	}


	private static Collection<Alignments> takeDifference(Alignments exon,  Alignments intron) {
		Collection<Alignments> rtrn = new ArrayList<Alignments>();
		if(!exon.overlaps(intron)){
			rtrn.add(exon);
		}
		else if(exon.fullyContained(intron)){
			//System.err.println("Odd case Exon "+exon.toUCSC()+" Intron: "+intron.toUCSC());
			Collection<Alignments> newExons=exon.excludeRegionFullyContained(intron);
			rtrn.addAll(newExons);
		} else if(!intron.fullyContained(exon)){
			Alignments newExon=exon.excludeRegion(intron);
			//System.err.println(intron.toUCSC()+" "+exon.toUCSC()+" "+newExon.toUCSC());
			rtrn.add(newExon);
		}
		
		return rtrn;
	}


	public static Collection<Alignments> DecollapseByIntronLocationOld(Collection<Alignments> exons, Collection<Alignments> introns){
		Collection<Alignments> rtrn=new TreeSet<Alignments>();
		
		Map<String, IntervalTree<Alignments>> intronTree=makeIntervalTree(introns);
		
		for(Alignments exon: exons){
			rtrn.add(exon);
			Iterator<Node<Alignments>> iter=intronTree.get(exon.getChr()).overlappers(exon.getStart(), exon.getEnd());
			while(iter.hasNext()){
				Alignments intron=iter.next().getValue();
				rtrn.addAll(takeDifference(exon, intron));
				
				
			}
		}
		
		return rtrn;
	}
	
	public static Collection<Alignments> DecollapseByIntronLocationWorking(Collection<Alignments> exons, Collection<Alignments> introns){
		if(introns==null || introns.isEmpty()){return exons;}
		Collection<Alignments> rtrn=new TreeSet<Alignments>();
		
		Map<String, IntervalTree<Alignments>> intronTree=makeIntervalTree(introns);
		
		for(Alignments exon: exons){
			rtrn.add(exon);
			Iterator<Node<Alignments>> iter=intronTree.get(exon.getChr()).overlappers(exon.getStart(), exon.getEnd());
			while(iter.hasNext()){
				Alignments intron=iter.next().getValue();
				//MG: Added this case because crashed when an exon overlapped an intron completely
				if(exon.fullyContained(intron)){
					//System.err.println("Odd case Exon "+exon.toUCSC()+" Intron: "+intron.toUCSC());
					Collection<Alignments> newExons=exon.excludeRegionFullyContained(intron);
					rtrn.addAll(newExons);
				}
				
				else if(!intron.fullyContained(exon)){
					Alignments newExon=exon.excludeRegion(intron);
					//System.err.println(intron.toUCSC()+" "+exon.toUCSC()+" "+newExon.toUCSC());
					rtrn.add(newExon);
				}
				
				
			}
		}
		
		return rtrn;
	}
	
	
	/*public static Set<Alignments> CollapseByUnion(Collection<Alignments> alignments){
		
		Collection<Alignments> accountedFor=new TreeSet();
		Set set=new TreeSet();
		Map<String, IntervalTree> trees=makeIntervalTree(alignments);
		//for each region get all overlappign and collapse by intersection
		for(Alignments align: alignments){
			if(!accountedFor.contains(align)){
			Collection<Alignments> regions=getOverlappingAll(align, trees);
			accountedFor.addAll(regions);
			Alignments collapsed=collapse(regions);
			if(collapsed.getSize()>0){set.add(collapsed);}
			}
		}
		
		
		return set;
	}*/
	
	
	public static Collection<RefSeqGene> updateBoundariesWithoutCrossingIntrons(Collection<RefSeqGene> genes){
		numOverlapGenes=10; //just to start
		
		Collection<RefSeqGene> rtrn=genes;
		
		int i=0;
		while(numOverlapGenes>1){
			numOverlapGenes=0;
			Set set=new TreeSet();
			Map<String, IntervalTree<RefSeqGene>> trees=makeIntervalTreeForGenes(rtrn);
			for(RefSeqGene align: rtrn){
				Iterator<Node<RefSeqGene>> regions=getOverlapping(align, trees);
				RefSeqGene collapsed=collapseGenes(regions);
				if(collapsed!=null){set.add(collapsed);}
			}
			rtrn=set;
			i++;
			//System.err.println("Iteration "+i+" Num Overlap "+numOverlapGenes);
		}
		
		return rtrn;
	}
	
		
	private static RefSeqGene collapseGenes(Iterator<Node<RefSeqGene>> regions){
		//these are all regions overlapping a given interval
		
		//now take all exons and collapse
		int i=0;
		Collection<Alignments> exons=new TreeSet();
		while(regions.hasNext()){
			i++;
			RefSeqGene gene=regions.next().getValue();
			exons.addAll(gene.getExonSet());
		}
		numOverlapGenes=Math.max(i, numOverlapGenes);
		
		exons=collapseByIntersection(exons, false);
		
		//then make into transcript
		RefSeqGene gene=null;
		if(exons!=null && !exons.isEmpty() && exons.size()>0){
			gene=new RefSeqGene(exons);
			//System.err.println(gene.getAlignment().toUCSC());
		}
		
		
		
		return gene;
	}
	
	private static Alignments trimExon(Alignments exon, Alignments intron){
		//if exon starts at or after intron and ends at or before intron -> return null
		if(intron.fullyContained(exon)){return null;}
		
		//if exon starts before intron
		if(exon.getStart()<intron.getStart()){
			return new Alignments(exon.getChr(), exon.getStart(), Math.min(exon.getEnd(), intron.getStart()));
		}
		
		//if exon starts after intron
		if(exon.getEnd()>intron.getEnd()){
			return new Alignments (exon.getChr(), Math.max(exon.getStart(), intron.getEnd()), exon.getEnd());
		}
		
		return null;
		
	}
	
	private static Iterator<Node<Alignments>> getOverlapping(Alignments align, Map<String, IntervalTree<Alignments>> trees){
		IntervalTree<Alignments> tree=trees.get(align.getChr());
		Iterator<Node<Alignments>> iter=tree.overlappers(align.getStart(), align.getEnd());
		return iter;
	}
	
	//can loop through everything in the overlapping list and get all of their overlaps as well
	private static Collection<Alignments> getOverlappingAll(Alignments align, Map<String, IntervalTree> trees){
		IntervalTree tree=trees.get(align.getChr());
		
		Collection<Alignments> rtrn=new TreeSet();
		
		//go through and get overlapping
		Iterator<Node<Alignments>> iter=tree.overlappers(align.getStart(), align.getEnd());
		rtrn=addAll(rtrn, iter);
		
		while(iter.hasNext()){
			Node<Alignments> node=iter.next();
			Iterator overlappers=tree.overlappers(node.getStart(), node.getEnd());
			rtrn=addAll(rtrn, overlappers);
		}
		
		return rtrn;
	}
	
	private static Collection<Alignments> addAll(Collection<Alignments> set, Iterator<Node<Alignments>> overlappers){
		Collection<Alignments> rtrn=set;
		while(overlappers.hasNext()){
			rtrn.add(overlappers.next().getValue());
		}
		return rtrn;
	}
	
	private static Iterator<Node<RefSeqGene>> getOverlapping(RefSeqGene gene, Map<String, IntervalTree<RefSeqGene>> trees){
		IntervalTree<RefSeqGene> tree=trees.get(gene.getChr());
		Iterator<Node<RefSeqGene>> iter=tree.overlappers(gene.getAlignment().getStart(), gene.getAlignment().getEnd());
		return iter;
	}
	
	
	
	
	private static Alignments collapse(Iterator<Node<Alignments>> iter, boolean intersection){
		int start=-1;
		int end=Integer.MAX_VALUE;
		String chr="";
		
		int i=0;
		while(iter.hasNext()){
			Alignments align=iter.next().getValue();
			if(i==0){start=align.getStart(); end=align.getEnd();}
			if(intersection){
				start=Math.max(start, align.getStart());
				end=Math.min(end,  align.getEnd());
			}
			else{
				start=Math.min(start, align.getStart());
				end=Math.max(end,  align.getEnd());
			}
			chr=align.getChr();
			i++;
		}
		
		numOverlap=Math.max(i, numOverlap); //TODO: This is not good. Why modify this class variable here? It makes it horribly non-thread safe for no reason.
		Alignments align=new Alignments(chr, start, end);
		return align;
	}
	
	/*private static Alignments collapse(Collection<Alignments> set){
		int start=-1;
		int end=Integer.MAX_VALUE;
		String chr="";
		
		int i=0;
		for(Alignments align: set){
			if(i==0){start=align.getStart(); end=align.getEnd();}
			start=Math.min(start, align.getStart());
			end=Math.max(end,  align.getEnd());
			chr=align.getChr();
			i++;
		}
		
		numOverlap=Math.max(i, numOverlap);
		Alignments align=new Alignments(chr, start, end);
		return align;
	}*/
	
	private static void write(String save, Collection set)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(Object k: set){writer.write(k.toString()+"\n");}
		
		writer.close();
	}
	
	public static Map<String, IntervalTree<Alignments>> makeIntervalTree(Collection<Alignments> alignments){
		Map<String, IntervalTree<Alignments>> rtrn=new TreeMap<String, IntervalTree<Alignments>> ();
		
		for(Alignments align: alignments){
			IntervalTree<Alignments> tree=new IntervalTree<Alignments>();
			String chr=align.getChr();
			if(rtrn.containsKey(chr)){tree=rtrn.get(chr);}
			//if(align.getStart()>=align.getEnd()){System.err.println("ERROR: " +align.toUCSC());}
			tree.put(align.getStart(), align.getEnd(), align);
			rtrn.put(align.getChr(), tree);
		}
		
		return rtrn;
	}
	
	
	public static Map<String, IntervalTree<RefSeqGene>> makeIntervalTreeForGenes(Collection<RefSeqGene> alignments){
		Map<String, IntervalTree<RefSeqGene>> rtrn=new TreeMap<String, IntervalTree<RefSeqGene>>();
		
		for(RefSeqGene align: alignments){
			IntervalTree<RefSeqGene> tree=new IntervalTree<RefSeqGene>();
			String chr=align.getChr();
			if(rtrn.containsKey(chr)){tree=rtrn.get(chr);}
			Node<RefSeqGene> node = tree.find(align.getAlignment().getStart(), align.getAlignment().getEnd()+1);
			if (node != null)
				node.incrementCount();
			else
				tree.put(align.getAlignment().getStart(), align.getAlignment().getEnd()+1, align);
			rtrn.put(align.getChr(), tree);
		}
		
		return rtrn;
	}
	
	public static void main(String[] args)throws IOException{
		Collection<Alignments> exons=BEDFileParser.loadAlignmentData(new File(args[0]));
		Collection<Alignments> introns = BEDFileParser.loadAlignmentData(new File(args[1]));
		Collection<Alignments> decollapsed = CollapseByIntersection.DecollapseByIntronLocation(exons, introns);
		String save=args[2];
		write(save, decollapsed);
	}


	public static Collection<Alignments> CollapseGenesByIntersection(Collection<RefSeqGene> temp, boolean intersection) {
		Collection<Alignments> alignments=convert(temp);
		
		return collapseByIntersection(alignments, intersection);
	}


	private static Collection<Alignments> convert(Collection<RefSeqGene> temp) {
		Collection<Alignments> rtrn=new TreeSet();
		
		for(RefSeqGene gene: temp){rtrn.add(gene.getAlignment());}
		
		return rtrn;
	}

	public static Map<String, IntervalTree<Alignments>> makeIntervalTreeForGeneExons(Collection<RefSeqGene> genes) {
		Collection<Alignments> exons=new TreeSet();
		
		for(RefSeqGene gene: genes){
			exons.addAll(gene.getSortedAndUniqueExons());
		}
		
		return makeIntervalTree(exons);
	}




	
	
}
