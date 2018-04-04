package broad.core.primer3;

import java.io.File;
import java.io.FileWriter;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.IntervalTree;
import broad.core.motif.SequenceMotif;
import broad.core.sequence.SequenceRegion;
import broad.pda.annotation.BEDFileParser;
import broad.pda.chromosome.Chromosome;
import broad.pda.gene.RefSeqGene;
import broad.pda.rnai.designer.ComputeMIRScore;


public class RNAiPCRPrimerDesigner {

	boolean repeatMask=true;
	int numDesigns=5;
	int min3;
	int min5;
	
	public RNAiPCRPrimerDesigner(Map<String, Collection<RefSeqGene>> alignmentsByChr, String sequenceDirectory, String save, boolean repeatMask, int numDesigns, int min3, int min5)throws Exception{
		this.min3=min3;
		this.min5=min5;
		this.numDesigns=numDesigns;
		this.repeatMask=repeatMask;
		Map<RefSeqGene, Collection<PrimerPair>>[] splicedPrimers=this.designConsecutivePrimers(alignmentsByChr, sequenceDirectory, true);
		Map<RefSeqGene, Collection<PrimerPair>> crossJunctionPrimers=this.designBestPrimers(alignmentsByChr, sequenceDirectory, true);
		Map<RefSeqGene, Collection<PrimerPair>> withinExonPrimers=this.designBestPrimers(alignmentsByChr, sequenceDirectory, false);
		//Map<RefSeqGene, Collection<PrimerPair>> filteredPrimers=filterByCrossingJunction(splicedPrimers);
		write(save, splicedPrimers);
		write(save, crossJunctionPrimers, "Standard");
		write(save, withinExonPrimers, "Within");
	}
	
	private Map<RefSeqGene, Collection<PrimerPair>> designBestPrimers(Map<String, Collection<RefSeqGene>> alignmentsByChr, String sequenceDirectory, boolean crossJunction) throws Exception {
		Map<RefSeqGene, Collection<PrimerPair>> rtrn=new TreeMap();
		
		for(String chr: alignmentsByChr.keySet()){
			//System.err.println("working on " +chr+" ...");
			String sequenceFile=sequenceDirectory+"/"+chr.replaceAll("chr", "").trim()+"/"+chr+".agp";
			Chromosome chrom = new Chromosome(sequenceFile);
			chrom.loadSequence();
			Map<RefSeqGene, Collection<PrimerPair>> map=designBestPrimers(chrom, alignmentsByChr.get(chr), crossJunction);
			rtrn.putAll(map);
		}
		
		return rtrn;
	}

	private Map<RefSeqGene, Collection<PrimerPair>> designBestPrimers(Chromosome chrom,	Collection<RefSeqGene> genes, boolean crossJunction) throws Exception {
		Map<RefSeqGene, Collection<PrimerPair>> rtrn=new TreeMap<RefSeqGene, Collection<PrimerPair>>();
		
		for(RefSeqGene gene: genes){
			Collection<PrimerPair> primers=designBestPrimers(gene, chrom, crossJunction);
			rtrn.put(gene, primers);
		}
		
		return rtrn;
	}

	private Collection<PrimerPair> designBestPrimers(RefSeqGene gene, Chromosome chrom, boolean crossJunction) throws Exception {
		return qPCRPrimerDesigner.designPCRPrimers(chrom, gene, repeatMask, numDesigns, crossJunction);
	}

	private Map<RefSeqGene, Collection<PrimerPair>> filterByCrossingJunction(Map<RefSeqGene, Collection<PrimerPair>> splicedPrimers) throws Exception {
		Map<RefSeqGene, Collection<PrimerPair>> rtrn=new TreeMap();
		
		for(RefSeqGene gene: splicedPrimers.keySet()){
			Collection<PrimerPair> primers=splicedPrimers.get(gene);
			Collection<PrimerPair> filtered=filterNonCrossing(primers, gene);
			rtrn.put(gene, filtered);
		}
		
		return rtrn;
	}

	private Collection<PrimerPair> filterNonCrossing(Collection<PrimerPair> primers, RefSeqGene gene) throws Exception {
		Collection<PrimerPair> rtrn=new HashSet();
		for(PrimerPair primer: primers){
			String leftPrimer=primer.getLeftPrimer();
			String rightPrimer=primer.getRightPrimer();
			if(leftPrimer!=null && rightPrimer!=null){
				boolean isCrossing=isCrossing(leftPrimer, rightPrimer, gene);
				if(isCrossing){rtrn.add(primer);}
				else{System.out.println("REJECTED "+leftPrimer+" "+rightPrimer+" Not crossing boundary");}
			}
		}
		return rtrn;
	}

	private boolean isCrossing(String leftPrimer, String rightPrimer,RefSeqGene gene) throws Exception {
		SequenceMotif leftPrimerMotif=new SequenceMotif(leftPrimer, 1);
		SequenceMotif rightPrimerMotif=new SequenceMotif(ComputeMIRScore.reverseComplement(rightPrimer), 1);
		
		Collection<SequenceRegion> leftMatches=leftPrimerMotif.match(gene.getSequenceObject());
		Collection<SequenceRegion> rightMatches=rightPrimerMotif.match(gene.getSequenceObject());
		
		IntervalTree<Integer> spliceLocations=gene.getSpliceJunctionCoordinatesTree();
		boolean left=false;
		boolean right=false;
		
		//System.err.println(spliceLocations.toCollection());
		
		for(SequenceRegion region: leftMatches){
			//System.err.println(region);
			Iterator iter=spliceLocations.overlappers(region.getStart(), region.getEnd());
			left=iter.hasNext();
		}
		
		for(SequenceRegion region: rightMatches){
			Iterator iter=spliceLocations.overlappers(region.getStart(), region.getEnd());
			right=iter.hasNext();
		}
		
		return left || right;
	}

	private Map<RefSeqGene, Collection<PrimerPair>>[] designConsecutivePrimers(Map<String, Collection<RefSeqGene>> alignmentsByChr, String sequenceDirectory, boolean spliced)throws Exception{
		Map<RefSeqGene, Collection<PrimerPair>> singles=new TreeMap();
		Map<RefSeqGene, Collection<PrimerPair>> doubles=new TreeMap();
		
		for(String chr: alignmentsByChr.keySet()){
			//System.err.println("working on " +chr+" ...");
			String sequenceFile=sequenceDirectory+"/"+chr.replaceAll("chr", "").trim()+"/"+chr+".agp";
			Chromosome chrom = new Chromosome(sequenceFile);
			chrom.loadSequence();
			Map[] map=designPrimers(chrom, alignmentsByChr.get(chr), spliced);
			singles.putAll(map[0]);
			doubles.putAll(map[1]);
		}
		
		Map[] rtrn={singles, doubles};
		return rtrn;
	}
	
	
	
	
	private void write(String save, Map<RefSeqGene, Collection<PrimerPair>>[] splicedPrimers)throws Exception{
		FileWriter writer=new FileWriter(save);
		writer.write("Gene Name\tJunction Spanning Primers\tLeft Primer\tRight Primer\tPrimer Pair Penalty\tPrimer Pair Position (Left,Right)\n");
		for(RefSeqGene align: splicedPrimers[1].keySet()){
			Collection<PrimerPair> primers=splicedPrimers[1].get(align);
			for(PrimerPair primer: primers){
				if(primer.getLeftPrimer()!=null && primer.getRightPrimer()!=null){
					writer.write(align.getName()+"\t2\t"+primer.getLeftPrimer()+"\t"+primer.getRightPrimer()+"\t"+primer.getPrimerPairPenalty()+"\t"+primer.getLeftPrimerPosition()+", "+primer.getRightPrimerPosition()+"\n");
				}
			}
		}
		
		for(RefSeqGene align: splicedPrimers[0].keySet()){
			Collection<PrimerPair> primers=splicedPrimers[0].get(align);
			for(PrimerPair primer: primers){
				if(primer.getLeftPrimer()!=null && primer.getRightPrimer()!=null){
					writer.write(align.getName()+"\t1\t"+primer.getLeftPrimer()+"\t"+primer.getRightPrimer()+"\t"+primer.getPrimerPairPenalty()+"\t"+primer.getLeftPrimerPosition()+", "+primer.getRightPrimerPosition()+"\n");
				}
			}
		}
		
		writer.close();
	}
	
	private void write(String save, Map<RefSeqGene, Collection<PrimerPair>> splicedPrimers, String string)throws Exception{
		FileWriter writer=new FileWriter(save, true);
		//writer.write("Gene Name\tJunction Spanning Primers\tLeft Primer\tRight Primer\tPrimer Pair Penalty\tPrimer Pair Position (Left,Right)\n");
		for(RefSeqGene align: splicedPrimers.keySet()){
			Collection<PrimerPair> primers=splicedPrimers.get(align);
			for(PrimerPair primer: primers){
				if(primer.getLeftPrimer()!=null && primer.getRightPrimer()!=null){
					writer.write(align.getName()+"\t"+string+"\t"+primer.getLeftPrimer()+"\t"+primer.getRightPrimer()+"\t"+primer.getPrimerPairPenalty()+"\t"+primer.getLeftPrimerPosition()+", "+primer.getRightPrimerPosition()+"\n");
				}
			}
		}
		
		
		writer.close();
	}
	

	private RefSeqGene[] junctionsSpanning(PrimerPair primer, RefSeqGene gene, boolean verbose) throws Exception {
		Collection<Integer> junctions=gene.getSpliceJunctionCoordinates();
		RefSeqGene leftPrimer=gene.trim(primer.getLeftPrimerPosition(), primer.getLeftPrimerPosition()+primer.getLeftPrimer().toCharArray().length);
		RefSeqGene rightPrimer=gene.trim(primer.getRightPrimerPosition()-primer.getRightPrimer().toCharArray().length, primer.getRightPrimerPosition());
		
		leftPrimer.setName(primer.getLeftPrimer());
		rightPrimer.setName(primer.getRightPrimer());
		
		if(verbose){
		System.out.println(leftPrimer);
		System.out.println(rightPrimer);
		}
		
		RefSeqGene[] array={leftPrimer, rightPrimer};
		return array;
	}

	private Map[] designPrimers(Chromosome chrom, Collection<RefSeqGene> genes, boolean spliced)throws Exception{
		Map<RefSeqGene, Collection> singles=new TreeMap();
		Map<RefSeqGene, Collection> doubles=new TreeMap();
		for(RefSeqGene gene: genes){
			try{
			//System.err.println("currently on gene "+gene.getName());
			if(gene.getNumExons()>1){
				Collection<PrimerPair>[] pairs=designPrimers(chrom, gene, spliced);
				singles.put(gene, pairs[0]);
				doubles.put(gene, pairs[1]);
			}
			else{System.out.println("Skipping gene "+gene.getName()+" because it has no introns");}
			}catch(Exception ex){System.err.println("Skipped\t"+gene.getName()+"\tbecause of an exception");}
		}
		Map[] rtrn={singles, doubles};
		return rtrn;
	}
	
	private Collection<PrimerPair>[] designPrimers(Chromosome chrom, RefSeqGene gene, boolean spliced)throws Exception{
		TreeSet<PrimerPair> both=new TreeSet();
		TreeSet<PrimerPair> singles=new TreeSet();
		
		Collection<PrimerPair> primers=qPCRPrimerDesigner.designIntronPrimers(chrom, gene, repeatMask, null, null, numDesigns*4, min3, min5);
		
		/**Try to find primers with 2 junctions spanned**/
		for(PrimerPair primer: primers){
			RefSeqGene[] primerPositions=this.junctionsSpanning(primer, gene, false);
				if(primerPositions[0].getNumExons()>1 && primerPositions[1].getNumExons()>1){both.add(primer);}
				else if(primerPositions[0].getNumExons()>1){
					//fix left primer and find right primer
					singles.add(primer);
					Collection<PrimerPair> newPrimers=qPCRPrimerDesigner.designIntronPrimers(chrom, gene, repeatMask, primer.getLeftPrimer(), null, numDesigns, min3, min5);
					newPrimers=getBoth(newPrimers, gene);
					both.addAll(newPrimers);
				}
				else if(primerPositions[1].getNumExons()>1){
					//fix right primer and find left
					singles.add(primer);
					Collection<PrimerPair> newPrimers=qPCRPrimerDesigner.designIntronPrimers(chrom, gene, repeatMask, null, primer.getRightPrimer(), numDesigns, min3, min5);
					newPrimers=getBoth(newPrimers, gene);
					both.addAll(newPrimers);
				}
			
		}
		/*if(primers==null || primers.isEmpty()){
			primers=qPCRPrimerDesigner.designIntronPrimers1End(chrom, gene, spliced);	
		}*/
		//gene.setSequence(gene.getSequence(chrom, repeatMask, false));
		//rtrn.put(gene, primers);
		
		Collection<PrimerPair> topSingles=getTopNum(singles, this.numDesigns);
		Collection<PrimerPair> topDoubles=getTopNum(both, this.numDesigns);
		
		Collection[] rtrn=new Collection[2];
		rtrn[0]=topSingles;
		rtrn[1]=topDoubles;
		return rtrn;
	}

	
	
	
	private Collection<PrimerPair> getTopNum(TreeSet<PrimerPair> primers, int n) {
		Collection<PrimerPair> rtrn=new TreeSet<PrimerPair>();
				
		int counter=0;
		for(PrimerPair pair: primers){
			if(!rtrn.contains(pair)){rtrn.add(pair); counter++;}
			if(counter>=n){return rtrn;}
		}
		
		return rtrn;
	}

	private Collection<PrimerPair> getBoth(Collection<PrimerPair> primers,RefSeqGene gene) throws Exception {
		Collection<PrimerPair> rtrn=new TreeSet<PrimerPair>();
		
		for(PrimerPair primer: primers){
			if(primer.getLeftPrimer()!=null && primer.getRightPrimer()!=null){
			RefSeqGene[] primerPositions=this.junctionsSpanning(primer, gene, false);
			if(primerPositions[0].getNumExons()>1 && primerPositions[1].getNumExons()>1){rtrn.add(primer);}
			}
		}
		
		return rtrn;
	}

	public static void main(String[] args)throws Exception{
		if(args.length>2){
			Map<String, Collection<RefSeqGene>> alignments=BEDFileParser.loadDataByChr(new File(args[0]));
			String seqDir=args[1];
			String save=args[2];
			boolean repeatMask=true;
			int numDesigns=5;
			int min3=4;
			int min5=7;
			
			if(args.length>3){repeatMask=new Boolean(args[3]);}
			if(args.length>4){numDesigns=new Integer(args[4]);}
			if(args.length>5){min3=new Integer(args[5]);}
			if(args.length>6){min5=new Integer(args[6]);}
			new RNAiPCRPrimerDesigner(alignments, seqDir, save, repeatMask, numDesigns, min3, min5);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=alignment file (Full BED) \n args[1]=sequence directory \n args[2]=save file \n args[3]=repeatMask (optional, default=true) \n args[4]=numDesigns (optional, default=5) \n args[5]=min number of bps on the 3' end (optional, default=4) \n args[6]=min num of bps on the 5' end (optional, default=7)";
	
}
