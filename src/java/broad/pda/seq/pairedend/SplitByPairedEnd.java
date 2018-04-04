package broad.pda.seq.pairedend;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;



public class SplitByPairedEnd {

	
	/*public SplitByPairedEnd(AlignmentDataModel data, Collection<RefSeqGene> genes, String save)throws IOException{
		for(String chr: data.getChromosomeLengths().keySet()){
			System.err.println(chr);
			Collection<RefSeqGene> split=splitByPairs(genes, data, chr);
			write(save, split);
		}
	}*/
	
	/*public SplitByPairedEnd(AlignmentDataModel data, File[] files, String saveDir)throws IOException{
		ContinuousDataAlignmentModel model=new ContinuousDataAlignmentModel(data);
		Collection<Alignments> cont=model.getOverlappingRegions();
		Collection<Alignments> collapse=CollapseByIntersection.CollapseByIntersection(cont, false);
		Map<Alignments, double[]> vals=model.scoreSegments(collapse);
		//write(save+".pairs", vals, 1);
		
		for(int i=0; i<files.length; i++){
			String save=saveDir+"/"+files[i].getName()+".split.segments";
			Collection<RefSeqGene> genes=BEDFileParser.loadData(files[i]);
			Collection<RefSeqGene> split=splitByPairs(vals.keySet(), genes);
			write(save, split);
		}
	}*/
	
	/*public static Collection<RefSeqGene> splitByPairs(Collection<RefSeqGene> genes, AlignmentDataModel pairs, String chr) throws IOException{
		if(pairs==null){return genes;}
		
		ContinuousDataAlignmentModel model=new ContinuousDataAlignmentModel(pairs);
		Collection<Alignments> cont=model.getOverlappingRegions(chr);
		Collection<Alignments> collapse=CollapseByIntersection.CollapseByIntersection(cont, false);
		
		Collection<RefSeqGene> split=splitByPairs(collapse, genes, chr);
		return split;
	}*/
	
	/*public SplitByPairedEnd(Collection<Alignments> vals, Collection<RefSeqGene> genes, String save)throws IOException{
		ContinuousDataAlignmentModel model=new ContinuousDataAlignmentModel(data);
		Collection<Alignments> cont=model.getOverlappingRegions();
		Collection<Alignments> collapse=CollapseByIntersection.CollapseByIntersection(cont, false);
		Map<Alignments, double[]> vals=model.scoreSegments(collapse);
		write(save+".pairs", vals, .05);
		
		Collection<RefSeqGene> split=splitByPairs(vals, genes);
		write(save, split);
	}*/
		
	
	private void write(String save, Collection<RefSeqGene> split) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(RefSeqGene gene: split){writer.write(gene+"\n");}
		
		writer.close();
	}

	
	//take an exon, if it spans a paired end segment split it
	/*private Collection<RefSeqGene> splitByPairs(Collection<Alignments> vals,	Collection<RefSeqGene> genes) {
		Collection<RefSeqGene> rtrn=new TreeSet();
		
		Map<String, IntervalTree<Alignments>> trees=makeTree(vals);
		
		for(RefSeqGene gene: genes){
			//only split end nodes
			//TODO Should extend this method to break up any NODE overlapping 2 paired ends (example Unc50)
			Collection<Alignments> newExons=new TreeSet();
			int i=0;
			for(Alignments exon: gene.getSortedAndUniqueExons()){
				//TODO add support if first and last exon then it should split on both ends and add both seperately
				
				//if first exon
				if(i==0){
					Alignments newExon=trim(exon, trees, true);
					newExons.add(newExon);
				}
				//if last exon
				else if(i==gene.getExonSizes().length-1){
					Alignments newExon=trim(exon, trees, false);
					newExons.add(newExon);
				}
				//else
				else{newExons.add(exon);}
				i++;
			}
			RefSeqGene updatedGene=new RefSeqGene(newExons);
			rtrn.add(updatedGene);
		}
		
		return rtrn;
	}*/
	
	//take an exon, if it spans a paired end segment split it
	/*private static Collection<RefSeqGene> splitByPairs(Collection<Alignments> vals, Collection<RefSeqGene> genes, String chr) {
		Collection<RefSeqGene> rtrn=new TreeSet();
		
		Map<String, IntervalTree<Alignments>> trees=makeTree(vals);
		
		for(RefSeqGene gene: genes){
			if(gene.getChr().equalsIgnoreCase(chr)){
			Collection<RefSeqGene> splitGenes=splitGene(gene, trees);
			for(RefSeqGene newGene: splitGenes){
				newGene.setOrientation(gene.getOrientation());
				rtrn.add(newGene);	
			}
			}
		}
		
		return rtrn;
	}*/
	
	/*private static Collection<RefSeqGene> splitGene(RefSeqGene gene, Map<String, IntervalTree<Alignments>> trees) {
		Collection<RefSeqGene> rtrn=new TreeSet();
		
		Collection<Alignments> introns=gene.getIntronSet();
		Collection<Alignments> exons=new TreeSet();
		
		for(Alignments exon: gene.getSortedAndUniqueExons()){
			//if exon overlaps 2 paired end segments
			int numOverlappers=trees.get(exon.getChr()).numOverlappers(exon.getStart(), exon.getEnd());
			if(numOverlappers>1){
				Iterator<Node<Alignments>> overlappers=trees.get(exon.getChr()).overlappers(exon.getStart(), exon.getEnd());
				//then split exon
				Collection<Alignments> splitExons=splitExon(exon, overlappers);
				overlappers=trees.get(exon.getChr()).overlappers(exon.getStart(), exon.getEnd());
				splitExons.addAll(splitExon1(exon,overlappers));
				//add split exons to exon collection
				exons.addAll(splitExons);
			}
			else{exons.add(exon);}
		}
		
		//make graph
		ChromosomeWithBubbles2 graph=makeGraph(exons, introns, gene.getChr());
		
		//return allPaths and orphans
		rtrn.addAll(graph.getGenePaths(0));
		rtrn.addAll(graph.getOrphanNodes());
		
		return rtrn;
	}*/
	
	public static Collection<Alignments> splitExon(Alignments exon, Iterator<Node<Alignments>> overlappers) {
		Collection<Alignments> rtrn=new TreeSet<Alignments>();
		
		while(overlappers.hasNext()){
			Alignments pe=overlappers.next().getValue();
			if(pe.fullyContained(exon)){rtrn.add(exon); }
			else if(exon.fullyContained(pe)){rtrn.add(pe);}
			else if(exon.getEnd()<pe.getEnd()){rtrn.add(new Alignments(exon.getChr(), Math.max(exon.getStart(), pe.getStart()), exon.getEnd()));}
			else if(exon.getStart()> pe.getStart()){rtrn.add(new Alignments(exon.getChr(), exon.getStart(), Math.min(exon.getEnd(), pe.getEnd())));}
			else{rtrn.add(new Alignments(exon.getChr(), Math.max(exon.getStart(), pe.getStart()), Math.min(pe.getEnd(), exon.getEnd())));}
		}
		
		return rtrn;
	}
	
	public static Collection<Alignments> splitExon1(Alignments exon, Iterator<Node<Alignments>> overlappers){
		Collection<Alignments> rtrn=new TreeSet<Alignments>();
		
		Collection<Alignments> pairedEnds=new TreeSet();
		while(overlappers.hasNext()){
			Alignments pe=overlappers.next().getValue();
			pairedEnds.add(pe);
		}
		
		if(pairedEnds.isEmpty()){return rtrn;}
		Object[] array=pairedEnds.toArray();
		Alignments r=new Alignments(exon.getChr(), exon.getStart(), ((Alignments)array[0]).getStart());
		if(r.getStart()<r.getEnd()){rtrn.add(r);}
		
		for(int i=0; i<pairedEnds.size()-1; i++){
			Alignments current=(Alignments)array[i];
			Alignments next=(Alignments)array[i+1];
			r=new Alignments(exon.getChr(), current.getEnd(), next.getStart());
			rtrn.add(r);
		}
		r=new Alignments(exon.getChr(), ((Alignments)array[array.length-1]).getEnd(), exon.getEnd());
		if(r.getStart()<r.getEnd()){rtrn.add(r);}
		
		return rtrn;
	}

	/*private static ChromosomeWithBubbles2 makeGraph(Collection<Alignments> exons, Collection<Alignments> introns, String chr){
		//collapse all exons and  keep all intons
				
		exons=CollapseByIntersection.CollapseByIntersection(exons, false);
		
		if(!introns.isEmpty()){
			exons=CollapseByIntersection.DecollapseByIntronLocation(exons, introns);
		}
		
		//for(Alignments exon: exons){System.out.println(exon);}
		
		//populate graph
		ChromosomeWithBubbles2 bubbles=new ChromosomeWithBubbles2(chr, exons, introns);
		return bubbles;
	}*/

	private static Map<String, IntervalTree<Alignments>> makeTree(Collection<Alignments> genes){
		Map<String, IntervalTree<Alignments>> rtrn=new TreeMap<String, IntervalTree<Alignments>>();
		
		for(Alignments align: genes){
			IntervalTree<Alignments> tree=new IntervalTree<Alignments>();
			String chr=align.getChr();
			if(rtrn.containsKey(chr)){tree=rtrn.get(chr);}
			Node node=tree.find(align.getStart(), align.getEnd());
			if(node==null){tree.put(align.getStart(), align.getEnd(), align);}
			else{
				node.incrementCount();
			}
			rtrn.put(align.getChr(), tree);
		}
		
		return rtrn;
		
	}
 
	private Alignments trim(Alignments exon, Map<String, IntervalTree<Alignments>> trees, boolean first) {
		
		//If has more than one overlapper then we need to trim it down
		int numOverlappers=trees.get(exon.getChr()).numOverlappers(exon.getStart(), exon.getEnd());
		if(numOverlappers>1){
			
			Node<Alignments> pe=null;
			if(!first){
				pe=getMin(trees.get(exon.getChr()).overlappers(exon.getStart(), exon.getEnd()));
				Alignments rtrn=new Alignments(exon.getChr(), exon.getStart(), Math.min(exon.getEnd(), pe.getEnd()));
				
				return rtrn;
			}
			else{
				pe=getMax(trees.get(exon.getChr()).overlappers(exon.getStart(), exon.getEnd()));
				Alignments rtrn=new Alignments(exon.getChr(), Math.max(exon.getStart(), pe.getStart()), exon.getEnd());
				
				return rtrn;
			}	
		}
		
		//else just return itself
		else{return exon;}
	}

	private Node<Alignments> getMin(Iterator<Node<Alignments>> iter){
		Node<Alignments> min=iter.next();
		while(iter.hasNext()){
			Node<Alignments> temp=iter.next();
			if(temp.getStart()< min.getStart()){min=temp;}
		}
		return min;
	}
	
	private Node<Alignments> getMax(Iterator<Node<Alignments>> iter){
		Node<Alignments> max=iter.next();
		while(iter.hasNext()){
			Node<Alignments> temp=iter.next();
			if(temp.getStart()> max.getStart()){max=temp;}
		}
		return max;
	}
	
	
	private void write(String save, Map<Alignments, double[]> cont, double alpha) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Alignments align: cont.keySet()){
			double[] vals=cont.get(align);
			//if(vals[0]<alpha){writer.write(align+"\t"+vals[0]+"\t"+vals[1]+"\n");}
			writer.write(align+"\t"+vals[0]+"\t"+vals[1]+"\n");
		}
		
		writer.close();
	}

	
	/*public static void main(String[] args)throws IOException{
		if(args.length>3){
			AlignmentDataModel data=new GenericAlignmentDataModel(args[0], args[1]);
			String save=args[2];
			File f=new File(args[3]);
			if(f.isDirectory()){
				new SplitByPairedEnd(data, f.listFiles(), save);
			}
			else{
				Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[3]));
				new SplitByPairedEnd(data, genes, save);
			}
			
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=Alignment file (Paired End file) \n args[1]=sizes \n args[2]=save (or savedir) \n args[3]=full bed (or directory)";
	*/
	
	
	/*public static void main(String[] args)throws IOException{
		if(args.length>2){
			Collection<Alignments> vals=BEDFileParser.loadAlignmentData(new File(args[0]));
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[1]));
			String save=args[2];
			new SplitByPairedEnd(vals, genes, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=pairs \n args[1]=genes \n args[2]=save";*/
	
	/*public static void main(String[] args)throws IOException{
		if(args.length>2){
		AlignmentDataModel data=new GenericAlignmentDataModel(args[0], args[1]);
		String save=args[2];
		ContinuousDataAlignmentModel model=new ContinuousDataAlignmentModel(data);
		Collection<Alignments> regions=model.getOverlappingRegions();
		writeAlign(save, regions);
		}
		else{System.err.println(usage);}
	}*/
	
	

	private static void writeAlign(String save, Collection<Alignments> regions) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Alignments align: regions){
			writer.write(align+"\n");
		}
		
		writer.close();
	}



	static String usage=" args[0]=pairs \n args[1]=sizes \n args[2]=save";
}
