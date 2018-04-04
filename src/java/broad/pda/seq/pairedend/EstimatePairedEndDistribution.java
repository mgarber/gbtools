package broad.pda.seq.pairedend;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeSet;

import org.broad.igv.sam.Alignment;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.math.EmpiricalDistribution;
import broad.core.util.CollapseByIntersection;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.graph.Path;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.AlignmentDataModelStats;

public class EstimatePairedEndDistribution {

	Collection<String> possibleChromosomes;
	EmpiricalDistribution distribution;
	
	//Assume that there is a distribution of paired end insert sizes that we can estimate directly from the data
	public EstimatePairedEndDistribution(AlignmentDataModel data, Collection<RefSeqGene> genes) throws IOException{
		//get all simple paths
		Collection<RefSeqGene> simplePaths=findAllSimplePaths(genes);
		//make distribution of paired ends that overlap the simple paths
		EmpiricalDistribution dist=estimatePairedDistribution(simplePaths, data);
		//print estimated distribution
		this.distribution=dist;
	}
	
	public EmpiricalDistribution getSizeDistribution(){return this.distribution;}
	
	private Collection<RefSeqGene> findAllSimplePaths(Collection<RefSeqGene> genes) {
		Collection<RefSeqGene> rtrn=new TreeSet<RefSeqGene>();
		this.possibleChromosomes=new TreeSet<String>();
		
		//Make interval tree of genes
		Map<String, IntervalTree<RefSeqGene>> paths=CollapseByIntersection.makeIntervalTreeForGenes(genes);
				
		for(RefSeqGene gene: genes){
			this.possibleChromosomes.add(gene.getChr());
			//for each gene see if there are overlappers
			int numOverlappers=paths.get(gene.getChr()).numOverlappers(gene.getStart(), gene.getEnd());
			//everything has to overlap with itself so if it has more than 1 count it as not unique
			//if not, add it to the collection
			if(numOverlappers==1){rtrn.add(gene);}
		}
		return rtrn;
	}
	
	public static Collection<RefSeqGene> findSimplePathsRefSeqGene(Collection<RefSeqGene> genes) {
		Collection<RefSeqGene> rtrn=new TreeSet<RefSeqGene>();

		//Make interval tree of genes
		Map<String, IntervalTree<RefSeqGene>> paths=CollapseByIntersection.makeIntervalTreeForGenes(genes);
				
		for(RefSeqGene gene: genes){
			//for each gene see if there are overlappers
			int numOverlappers=paths.get(gene.getChr()).numOverlappers(gene.getStart(), gene.getEnd());
			//everything has to overlap with itself so if it has more than 1 count it as not unique
			//if not, add it to the collection
			if(numOverlappers==1){rtrn.add(gene);}
		}
		return rtrn;
	}
	
	public static Collection<Path> findSimplePaths(Collection<Path> paths){
		Collection<Path> rtrn=new TreeSet<Path>();
		
		Map<String, IntervalTree<Path>> trees=TreeUtils.makeIntervalTreeByPath(paths);
		
		for(Path path: paths){
			//for each gene see if there are overlappers
			int numOverlappers=trees.get(path.getChromosome()).numOverlappers(path.getStart(), path.getEnd());
			//everything has to overlap with itself so if it has more than 1 count it as not unique
			//if not, add it to the collection
			if(numOverlappers==1){rtrn.add(path);}
		}
		
		return rtrn;
	}

	private void write(String save, EmpiricalDistribution dist) throws IOException {
		dist.write(save);
	}

	private EmpiricalDistribution estimatePairedDistribution(Collection<RefSeqGene> simplePaths, AlignmentDataModel pairs) throws IOException {
		//make interval tree of paths
		Map<String, IntervalTree<RefSeqGene>> paths=CollapseByIntersection.makeIntervalTreeForGenes(simplePaths);
		ArrayList<Double> sizes=new ArrayList<Double>();
		
		if(possibleChromosomes==null){this.possibleChromosomes=pairs.getChromosomeLengths().keySet();}
		
		//go through all paired alignments
		for(String chr: possibleChromosomes){
			System.err.println(chr);
			IntervalTree<Alignment> tree=pairs.getIntervalTree(chr, 0, pairs.getChromosomeLengths().get(chr));
			Iterator<Node<Alignment>> iter=tree.iterator();
			while(iter.hasNext()){
				Alignment read=iter.next().getValue();				
				//ask which path it overlaps (has to be unique path!!)
				Iterator<Node<RefSeqGene>>overlappingPaths=paths.get(read.getChromosome()).overlappers(read.getAlignmentStart(), read.getAlignmentEnd());
				//estimate the insert size
				double size=estimateSize(read, overlappingPaths);
				//add it to the distribution
				if(size>0){sizes.add(size);}
			}
		}
		
		//make empirical distribution
		//TODO Consider making a discrete empiricalDistribution
		EmpiricalDistribution dist=new EmpiricalDistribution(sizes, 500);
		//return the distribution
		return dist;
	}
	
	public static EmpiricalDistribution estimatePairedInsertDistributionRefSeqGene(Collection<RefSeqGene> simplePaths, IntervalTree<Alignments> pairedEndTree) {
		//make interval tree of paths
		Map<String, IntervalTree<RefSeqGene>> paths=CollapseByIntersection.makeIntervalTreeForGenes(simplePaths);
		ArrayList<Double> sizes=new ArrayList<Double>();

		//if(possibleChromosomes==null){this.possibleChromosomes=pairs.getChromosomeLengths().keySet();}
		
		//go through all paired alignments
		Iterator<Node<Alignments>> iter = pairedEndTree.iterator();
		while(iter.hasNext()){
			Node<Alignments> readNode = iter.next();
			Alignments read=readNode.getValue();
			int readNodeCount = readNode.getNumReplicates();

			//ask which path it overlaps (has to be unique path!!)
			Iterator<Node<RefSeqGene>>overlappingPaths=paths.get(read.getChromosome()).overlappers(read.getStart(), read.getEnd());
			//estimate the insert size
			double size=estimateInsertSizeRefSeqGene(read, overlappingPaths);
			//add it to the distribution
			if(size>0){
				sizes.add(size);
			}
		}

		//make empirical distribution
		//TODO Consider making a discrete empiricalDistribution
		EmpiricalDistribution dist=new EmpiricalDistribution(sizes, 500);
		//return the distribution
		return dist;
	}
	
	
	/**
	public static EmpiricalDistribution estimatePairedInsertDistribution(Collection<Path> simplePaths, AlignmentDataModelStats pairs) {
		//make interval tree of paths
		Map<String, IntervalTree<Path>> paths=TreeUtils.makeIntervalTreeByPath(simplePaths);
		ArrayList<Double> sizes=new ArrayList<Double>();
		
		//if(possibleChromosomes==null){this.possibleChromosomes=pairs.getChromosomeLengths().keySet();}
		
		//go through all paired alignments
		for(String chr: paths.keySet()){
			try{
			IntervalTree<Alignment> tree=pairs.getIntervalTree(chr, 0, pairs.getChromosomeLengths().get(chr));
			Iterator<Node<Alignment>> iter=tree.iterator();
			while(iter.hasNext()){
				Alignment read=iter.next().getValue();
				
				//ask which path it overlaps (has to be unique path!!)
				Iterator<Node<Path>>overlappingPaths=paths.get(read.getChromosome()).overlappers(read.getAlignmentStart(), read.getAlignmentEnd());
				//estimate the insert size
				double size=estimateInsertSize(read, overlappingPaths);
				//add it to the distribution
				if(size>0){sizes.add(size);}
			}
			}catch(NullPointerException ex){System.err.println("Missing "+chr);}
		}
		
		//make empirical distribution
		//TODO Consider making a discrete empiricalDistribution
		EmpiricalDistribution dist=new EmpiricalDistribution(sizes, 500);
		//return the distribution
		return dist;
	}
	 * @throws IOException 
	*/
	
	public static EmpiricalDistribution estimatePairedInsertDistribution(Collection<Path> simplePaths, AlignmentDataModelStats pairs) throws IOException {
		//make interval tree of paths
		Map<String, IntervalTree<Path>> paths=TreeUtils.makeIntervalTreeByPath(simplePaths);
		ArrayList<Double> sizes=new ArrayList<Double>();
		//if(possibleChromosomes==null){this.possibleChromosomes=pairs.getChromosomeLengths().keySet();}

		//go through all paired alignments
		for(String chr: paths.keySet()){
			try{
			IntervalTree<Alignments> pairedEndTree = pairs.getData().getFullIntervalTreeAsAlignments(chr);
			Iterator<Node<Alignments>> iter = pairedEndTree.iterator();
			while(iter.hasNext()){
				Node<Alignments> readNode = iter.next();
				Alignments read=readNode.getValue();
				int readNodeCount = readNode.getNumReplicates();

				//ask which path it overlaps (has to be unique path!!)
				Iterator<Node<Path>>overlappingPaths=paths.get(read.getChromosome()).overlappers(read.getStart(), read.getEnd());
				//estimate the insert size
				double size=estimateInsertSize(read, overlappingPaths);
				//add it to the distribution
				if(size>0){
					sizes.add(size);
				}
			}
			}catch(NullPointerException ex){System.err.println("Missing "+chr);}
		}
		
		//make empirical distribution
		//TODO Consider making a discrete empiricalDistribution
		EmpiricalDistribution dist=new EmpiricalDistribution(sizes, 500);
		//return the distribution
		return dist;
	}
	
	private double estimateSize(Alignment read, Iterator<Node<RefSeqGene>> overlappingPaths) {
		//first check if the overlapping paths has an overlapper (and ONLY one overlapper)
		if(!overlappingPaths.hasNext()){return -999;}
		
		int i=0;
		double distance=-1;
		while(overlappingPaths.hasNext()){
			RefSeqGene path=overlappingPaths.next().getValue();
			RefSeqGene insert=path.trimAbsolute(read.getAlignmentStart(), read.getAlignmentEnd());
			if(insert==null){distance=-1;}
			else{distance=insert.getSize();}
			i++;
		}
		
		if(i>1){return -99;}
		return distance;
	}
	
	private static double estimateInsertSizeRefSeqGene(Alignments read, Iterator<Node<RefSeqGene>> overlappingPaths) {
		//first check if the overlapping paths has an overlapper (and ONLY one overlapper)
		if(!overlappingPaths.hasNext()){return -999;}

		int i=0;
		double distance=-1;
		while(overlappingPaths.hasNext()){
			Node<RefSeqGene> node = overlappingPaths.next();
			RefSeqGene path=node.getValue();
			RefSeqGene insert=path.trimAbsolute(read.getStart(), read.getEnd());
			if(insert==null){distance=-1;}
			else{distance=insert.getSize();}
			i += node.getNumReplicates();
		}
		
		if(i>1){return -99;}
		return distance;
	}
	
	private static double estimateInsertSize(Alignments read, Iterator<Node<Path>> overlappingPaths) {
		//first check if the overlapping paths has an overlapper (and ONLY one overlapper)
		if(!overlappingPaths.hasNext()){return -999;}

		int i=0;
		double distance=-1;
		while(overlappingPaths.hasNext()){
			Node<Path> node = overlappingPaths.next();
			Path path=node.getValue();
			RefSeqGene insert=path.toGene().trimAbsolute(read.getStart(), read.getEnd());
			if(insert==null){distance=-1;}
			else{distance=insert.getSize();}
			i += node.getNumReplicates();
		}
		
		if(i>1){return -99;}
		return distance;
	}
	
	private static double estimateInsertSize(Alignment read, Iterator<Node<Path>> overlappingPaths) {
		//first check if the overlapping paths has an overlapper (and ONLY one overlapper)
		if(!overlappingPaths.hasNext()){return -999;}

		int i=0;
		double distance=-1;
		while(overlappingPaths.hasNext()){
			Path path=overlappingPaths.next().getValue();
			RefSeqGene insert=path.toGene().trimAbsolute(read.getAlignmentStart(), read.getAlignmentEnd());
			if(insert==null){distance=-1;}
			else{distance=insert.getSize();}
			i++;
		}
		
		if(i>1){return -99;}
		return distance;
	}
	
	/*public static void main(String[] args)throws IOException{
		if(args.length>3){
			AlignmentDataModel data=new GenericAlignmentDataModel(args[0], args[1]);
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[2]));
			String save=args[3];
			new EstimatePairedEndDistribution(data, genes, save);
		}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=paired alignments \n args[1]=sizes \n args[2]=genes \n args[3]=save";
	*/
}
