package broad.pda.seq.segmentation;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.TreeMap;
import java.util.TreeSet;

import org.broad.igv.sam.Alignment;

import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.math.EmpiricalDistribution;
import broad.core.util.CollapseByIntersection;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.graph.ChromosomeWithBubbles2;
import broad.pda.seq.graph.ChromosomeWithBubbles2.BubbleEdge;
import broad.pda.seq.graph.ChromosomeWithBubbles2.EdgeSourceType;
import broad.pda.seq.graph.Path;
import broad.pda.seq.pairedend.EstimatePairedEndDistribution;
import broad.pda.seq.pairedend.PairedEndSignificanceInExpressionBin;
import broad.pda.seq.utils.AlignmentsPair;
import broad.pda.seq.utils.GenesInExpressionBins;
import net.sf.samtools.util.CloseableIterator;

public class AddEdgesUsingPairedEnds {

	private static boolean isStrandSpecific=false;
	private static final double fwer = 0.05; 
	private static final int MAX_PAIRED_END_SIZE = 100000;
	private static final int MIN_PAIRED_END_SUPPORT = 2;
	private static final double MIN_PROBABILITY_IN_PAIRED_INSERT_DISTRIBUTION_TO_USE_WHEN_ADDING_EDGE = 0.1;
	private static final double MAX_PROBABILITY_IN_PAIRED_INSERT_DISTRIBUTION_TO_USE_WHEN_ADDING_EDGE = 0.9;
	private static final int NUM_OF_BINS_FOR_GENES = 200;
	private static final double MIDDLE_PERCENT_OF_EXPRESSION_FOR_EXTRACTING_MIDDLE_GENES = 0.9;
	private static final double MIN_POISSON_FOR_CONNECTING_PATHS = 0.1;

	private Collection<RefSeqGene> filterOrphanNodesNotInGene(Collection<RefSeqGene> genes) {
		Collection<RefSeqGene> rtrn=new TreeSet();

		Map<String, IntervalTree<RefSeqGene>> trees=CollapseByIntersection.makeIntervalTreeForGenes(genes);
		
		for(RefSeqGene gene: genes){
			if(gene.getNumExons()==1){
				Iterator<Node<RefSeqGene>> iter=trees.get(gene.getChr()).overlappers(gene.getStart(), gene.getEnd());
				boolean overlaps=false;
				while(iter.hasNext() && !overlaps){
					overlaps=(iter.next().getValue().getNumExons()>1);
				}
				if(!overlaps){rtrn.add(gene);}
			}
			else{rtrn.add(gene);}
		}
		
		return rtrn;
	}


	
	private static Collection<RefSeqGene> makeAllPaths(Collection<RefSeqGene> posPaths, Collection<RefSeqGene> negPaths, Collection<RefSeqGene> orphans) {
		Collection<RefSeqGene> rtrn=new TreeSet();
	
		for(RefSeqGene gene: posPaths){
			gene.setOrientation("+");
			rtrn.add(gene);
		}
		
		for(RefSeqGene gene: negPaths){
			gene.setOrientation("-");
			rtrn.add(gene);
		}
		
		for(RefSeqGene gene: orphans){
			gene.setOrientation("*");
			rtrn.add(gene);
		}
		
		return rtrn;
	}

	private Collection<RefSeqGene> getNewEdges(Collection<RefSeqGene> genes, Collection<RefSeqGene> posPaths, Collection<RefSeqGene> negPaths) {
		Collection<RefSeqGene> rtrn=new TreeSet();
		Map<String, IntervalTree<RefSeqGene>> geneTree=CollapseByIntersection.makeIntervalTreeForGenes(genes);
		for(RefSeqGene posGene: posPaths){
			Iterator<Node<RefSeqGene>> overlappers=geneTree.get(posGene.getChr()).overlappers(posGene.getStart(), posGene.getEnd());
			boolean b=isSame(posGene, overlappers);
			if(!b){rtrn.add(posGene);}
		}
		return rtrn;
	}

	private boolean isSame(RefSeqGene posGene, Iterator<Node<RefSeqGene>> iter) {
		if(iter==null){return false;}
		while(iter.hasNext()){
			RefSeqGene gene=iter.next().getValue();
			if(posGene.equals(gene)){return true;}
		}
		
		return false;
	}

	private void write(String save, Collection<RefSeqGene> novel) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(RefSeqGene gene: novel){
			writer.write(gene+"\n");
		}
		
		writer.close();
	}

	private void write(String save, Collection<RefSeqGene> posPaths, Collection<RefSeqGene> negPaths, Collection<RefSeqGene> orphans)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(RefSeqGene gene: posPaths){
			gene.setOrientation("+");
			writer.write(gene+"\n");
		}
		
		for(RefSeqGene gene: negPaths){
			gene.setOrientation("-");
			writer.write(gene+"\n");
		}
		
		for(RefSeqGene gene: orphans){
			gene.setOrientation("*");
			writer.write(gene+"\n");
		}
		
		writer.close();
	}
	
	private static Collection<RefSeqGene> getOrphans(Collection<RefSeqGene> orphansPos, Collection<RefSeqGene> orphansNeg){
		Collection<RefSeqGene> rtrn=new TreeSet();
		
		for(RefSeqGene pos: orphansPos){
			if(orphansNeg.contains(pos)){rtrn.add(pos);}
		}
		
		return rtrn;
	}

	
	private void write(String save, ChromosomeWithBubbles2 graph)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		Collection<RefSeqGene> paths=graph.getGenePaths(0);
		
		for(RefSeqGene gene: paths){writer.write(gene+"\n");}
		
		writer.close();
	}
	

	
	public static Collection<Path> usePairedEndsToAddEdgesToPaths(ChromosomeWithBubbles2 graph, ContinuousDataAlignmentModel data, Collection<Path> paths, Map<Path, RefSeqGene> pathsToGenes, AlignmentDataModelStats pairedData, String chr, int start, int end, String saveToFile) throws IOException {
		// Collection<Path> paths=graph.getAllPaths();
		
		if(pairedData==null){return paths;}
		
		//Find simple paths for estimation of size
		Collection<Path> simplePaths=EstimatePairedEndDistribution.findSimplePaths(paths);
		System.err.println("Got Simple paths...");
		
		//make distribution of paired ends that overlap the simple paths
		EmpiricalDistribution pairedDist=EstimatePairedEndDistribution.estimatePairedInsertDistribution(simplePaths, pairedData);
		try {
			pairedDist.write(saveToFile + ".paired.size.dist.txt");
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		System.err.println("Estimated insert size distribution from paired reads mean: " + pairedDist.getMean() );
		
		/*
		try {
			pairedDist.write(saveToFile + ".pairedInsertSizeDistribution1.txt");
		}
		catch (IOException e) { e.printStackTrace(); }
		*/
			
		//set graph paired distribution
		//graph.setEstimatedPairedDistribution( pairedDist);
		
		Collection<RefSeqGene> genes = pathsToGenes.values();
		GenesInExpressionBins<PairedEndSignificanceInExpressionBin> genesInBins = new GenesInExpressionBins<PairedEndSignificanceInExpressionBin>(genes, NUM_OF_BINS_FOR_GENES, PairedEndSignificanceInExpressionBin.class, data, chr, MIDDLE_PERCENT_OF_EXPRESSION_FOR_EXTRACTING_MIDDLE_GENES);
		
		IntervalTree<Collection<Path>> plusPathTree=makePathTreeByExon(paths, "+");
		IntervalTree<Collection<Path>> minusPathTree=makePathTreeByExon(paths, "-");
		System.err.println("Made path trees");
		
		joinPathsEfficient(graph, paths, plusPathTree, pairedData, chr, start, end, "+", pairedDist, pathsToGenes, genesInBins);
		
		joinPathsEfficient(graph, paths, minusPathTree, pairedData, chr, start, end, "-", pairedDist, pathsToGenes, genesInBins);
		try {
			graph.writeGraph(saveToFile + ".graph.post.add.pairs.dot", new Alignments(graph.getName(), start, end));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		//graph.updatePairedEndEdges(pairedData);
		//IntervalTree<Collection<Path>> pathTree=makePathTree(graph.getPaths(0));
		//addPairsToPaths(pathTree, pairedData, chr, start, end);
		
		/*Collection<BubbleEdge> edges=graph.getEdges();
		for(BubbleEdge edge: edges){
			if(edge.getType().equals(EdgeSourceType.SPLICED)){System.out.println(edge.getConnection());}
		}*/
		
		//return toCollection(pathTree);
		Collection<Path> rtrn=graph.getPaths(1);
		rtrn.addAll(graph.getOrphanPaths());
		return rtrn;
	}
	
	//Rewrite of adding paired ends by looping through the exons
	private static void joinPaths(ChromosomeWithBubbles2 graph, Collection<Path> paths, IntervalTree<Collection<Path>> pathTreeByExons,AlignmentDataModelStats pairedData, String chr, int start, int end, String orientation) throws IOException{
		//loop through all exons and get paired end reads overlapping it
		System.err.println("Number of paths: " + paths.size());
		
		Collection<Alignments> exons=new TreeSet<Alignments>();
		for(Path path: paths){
			RefSeqGene gene=path.toGene();
			if(gene.getNumExons()>1){
				exons.addAll(gene.getSortedAndUniqueExons());
			}
		}
		
		
		int i=0;
		//Go through each exon
		for(Alignments exon: exons){
			
			//System.err.println("At exon "+exon.toUCSC()+" "+i+" out of "+exons.size());
			i++;
			
			//get overlapping reads
			IntervalTree<Alignment> pairedReads=pairedData.getIntervalTreeTruncatedCached(exon.getChromosome(), exon.getStart(), exon.getEnd());
			
			//for each read get overlapping exons
			Iterator<Node<Alignment>> reads=pairedReads.overlappers(exon.getStart(), exon.getEnd());
			while(reads.hasNext()){
				Alignment read=reads.next().getValue();
				Collection<Path> paths1=new TreeSet();
				if(read.getAlignmentStart() < exon.getStart() ){
					paths1=getOverlappingPaths(pathTreeByExons.overlappers(read.getStart()-1, read.getStart()+1) );
				}
				else if(read.getAlignmentEnd() >=exon.getEnd()){
					paths1=getOverlappingPaths(pathTreeByExons.overlappers(read.getEnd()-1, read.getEnd()+1));
				}
				
				Iterator<Node<Collection<Path>>> pathOverlappers=pathTreeByExons.overlappers(exon.getStart(), exon.getEnd());
				
				//test with one
				boolean traverse=false;
				try{
					Path testPath1=pathOverlappers.next().getValue().iterator().next();
					Path testPath2=paths1.iterator().next();
					if(!testPath1.toGene().getAlignment().overlaps(testPath2.toGene().getAlignment())){
						traverse=true;
					}
				}catch(Exception ex){}
				
				
				if(traverse){
				while(pathOverlappers.hasNext()){
					Collection<Path> paths2=pathOverlappers.next().getValue();
					int counter=0;
					for(Path path1: paths1){
						for(Path path2: paths2){
							//for each pair see if the djistra distance is such that they are connected
							if(!path1.toGene().getAlignment().overlaps(path2.toGene().getAlignment())){
								//make join
								BubbleEdge edge=getEdge(path1, path2, orientation, 1, graph);
								if(edge != null) {
									boolean b=graph.addEdge(edge, edge.getLeftNode(), edge.getRightNode());
								}
							}
							counter++;
						}
					}
				}
				
			}	
		}	
		}
	}
	
	//Rewrite of adding paired ends to the path
	private static void joinPathsEfficient(ChromosomeWithBubbles2 graph, Collection<Path> paths, IntervalTree<Collection<Path>> pathTree,AlignmentDataModelStats pairedData, String chr, int start, int end, String orientation, EmpiricalDistribution pairedInsertDistribution, Map<Path, RefSeqGene> pathsToGenes, GenesInExpressionBins<PairedEndSignificanceInExpressionBin> genesInBins) throws IOException{
		//loop through all paths and get paired end reads overlapping it
		//System.err.println("Number of paths: " + paths.size());
		Collection<LightweightGenomicAnnotation> visitedExons=new TreeSet<LightweightGenomicAnnotation> ();
		
		//Map<LightweightGenomicAnnotation, Integer> connectionCounts=new TreeMap<LightweightGenomicAnnotation, Integer>();

		int i=0;
		for(Path path: paths){
			if (path.getNumberOfEdges() <=1 || path.isToComplex()) {
				//if(path.isToComplex()){System.err.println("Skipping path "+path.getChromosome()+":"+path.getStart()+"-"+path.getEnd());}
				continue;
			}
			long startTime = System.nanoTime();
			//System.err.println("Finding partners for path edge count("+path.getNumberOfEdges()+") --- " + path.toGene().toBED());
			//Check to see whether path's orientation is compatible with join direction.
			if(path.getOrientation().equalsIgnoreCase(orientation) || path.getOrientation().equalsIgnoreCase("*")){
				long startPartT = System.nanoTime();
				IntervalTree<Alignment> pairedReads=pairedData.getIntervalTreeTruncatedCached(path.getChromosome(), path.getStart(), path.getEnd());
				//System.err.println("\tGot reads -- time("+(System.nanoTime() - startPartT)+ " Cache chunk: " + pairedData.getChunkStart() +"-" + pairedData.getChunkEnd() + " size " + (   pairedData.getChunkEnd()-pairedData.getChunkStart() ) + " #reads: " + pairedReads.size());
				//Map<Path, Integer> overlappingPaths=new TreeMap<Path, Integer>();
				
				//long startTime=System.currentTimeMillis();
				startPartT = System.nanoTime();
				Collection<LightweightGenomicAnnotation> exons=path.getSortedNodes();
				//System.err.println("Current Path "+path.getChromosome()+":"+path.getStart()+"-"+path.getEnd()+" path exons "+exons.size());
								
				//System.err.println("\tUp to getSortedNodes took " + (System.nanoTime() - startTime));
				long exonTime = System.nanoTime();
				//All paths that any exon on this particular path does overlap.
				Collection<Path> allPaths=new TreeSet<Path>();
				Map<Path, Integer> overlappingPaths=new TreeMap<Path, Integer>();
				Map<Path, Map<AlignmentsPair, Set<Alignments>>> pairedReadsForOverlappingPathsByExonPair = new TreeMap<Path, Map<AlignmentsPair, Set<Alignments>>>();
				for(LightweightGenomicAnnotation exon: exons){
					//MG: Added check if we already visited this exon. If so then skip
					if(visitedExons.contains(exon)){
						//System.err.println("Already saw this exon so we wont do anything with it "+exon.toUCSC());
						continue;
					}
					
					visitedExons.add(exon);
					
					//TODO use all paths overlapping exon to add them ALL
					//all paths with this exon

					Iterator<Node<Collection<Path>>> iter=pathTree.overlappers(exon.getStart(), exon.getEnd());
					
	
					//int pathCount=0;
					while(iter.hasNext()){
						allPaths.addAll(iter.next().getValue());
						//pathCount+=iter.next().getValue().size();
					}
					
					//Collection<Alignments> usedReads=new TreeSet<Alignments>();
					startPartT = System.nanoTime();
					Iterator<Node<Alignment>> pairIter=pairedReads.overlappers(exon.getStart(), exon.getEnd());
					//System.err.println("\tRead Iterator time: " + (System.nanoTime() - startPartT));
					//System.err.println("\t\tGoing through path exons, exon " + exon.toUCSC() + " overlapper reads took " + (System.nanoTime() - startPartT)+" has overlapping reads "+pairedReads.numOverlappers(exon.getStart(), exon.getEnd())+" "+allPaths.size()+" paths that contain it");
					int counter=0;
					startPartT = System.nanoTime();
					long totalReadTime = 0;
					while(pairIter.hasNext()){
						
						long internalStarts = 0;
						Alignment read=pairIter.next().getValue();
						//Need to check whether the exon contains one end of the read. Added by MG on 05/25/2010
						if(! (exon.contains(new BasicLightweightAnnotation(read.getChr(), read.getStart(), read.getStart() + 1)) || exon.contains(new BasicLightweightAnnotation(read.getChr(), read.getEnd()-1, read.getEnd())))) {
							continue;
						}
						//Alignments readRegion=new Alignments(read.getChromosome(), read.getAlignmentStart(), read.getAlignmentEnd());
						//if(!usedReads.contains(readRegion)){
							//usedReads.add(readRegion);
						long startRead = System.nanoTime();
						if(Math.abs(read.getAlignmentEnd() - read.getAlignmentStart()) < MAX_PAIRED_END_SIZE ) {
						//TODO if read is in the same path then skip it
							if(read.getAlignmentStart() < path.getStart() ){
								internalStarts = System.nanoTime();
								Stack<Path> paths1=getOverlappingPaths(pathTree.overlappers(read.getStart()-1, read.getStart()+1) );
								//System.err.println("\t\t\tGot overlapping paths, took "+(System.nanoTime() - internalStarts) + " tree length " + pathTree.size());
								for(Path p: paths1){
									if(p.overlaps(path) ) {
										//if(overlappingPaths.containsKey(p)){System.err.println("Skipping path "+p);}
										continue;
									} 
									
									internalStarts = System.nanoTime();
									int count= overlappingPaths.containsKey(p) ? overlappingPaths.get(p) + 1: 1;
									//System.err.print("\t\t\tGoing to add -- time to find in map("+(System.nanoTime()- internalStarts)+")-- " + p.toGene() + " to overlappers " + path.toGene()+ " , current map size "+ overlappingPaths.size() + ", supporting read count " + count);
									internalStarts = System.nanoTime();
									overlappingPaths.put(p, count);
									
									Alignments pExon = p.toGene().getSingleOverlappingExon(new Alignments(read.getChr(), read.getStart()-1, read.getStart()+1));
									AlignmentsPair exonPair = new AlignmentsPair(new Alignments(exon.getChromosome(), exon.getStart(), exon.getEnd()), pExon);
									Alignments pairedReadAlignment = new Alignments(read.getChr(), read.getAlignmentStart(), read.getAlignmentEnd());
									
									if (!pairedReadsForOverlappingPathsByExonPair.containsKey(p)) {
										Map<AlignmentsPair, Set<Alignments>> m = new HashMap<AlignmentsPair, Set<Alignments>>();
										pairedReadsForOverlappingPathsByExonPair.put(p, m);
									}
									
									Map<AlignmentsPair, Set<Alignments>> exonPairsMappedToReadAlignments = pairedReadsForOverlappingPathsByExonPair.get(p);
									if (!exonPairsMappedToReadAlignments.containsKey(exonPair)) {
										Set<Alignments> s = new TreeSet<Alignments>();
										exonPairsMappedToReadAlignments.put(exonPair, s);
									}
									exonPairsMappedToReadAlignments.get(exonPair).add(pairedReadAlignment);
									
									//System.err.println(" size after putting back "+ overlappingPaths.size() + " put back took " + (System.nanoTime()- internalStarts));
								}
							} else if ( read.getAlignmentEnd() >=path.getEnd()) {
								internalStarts = System.nanoTime();
								Collection<Path> paths1=getOverlappingPaths(pathTree.overlappers(read.getEnd()-1, read.getEnd()+1));
								//System.err.println("\t\t\tGot overlapping paths, took "+(System.nanoTime() - internalStarts) + " tree length " + pathTree.size());
								for(Path p: paths1){
									if(p.overlaps(path) /*|| overlappingPaths.containsKey(p)*/) {
										//if(overlappingPaths.containsKey(p)){System.err.println("Skipping path "+p);}
										continue;
									}
									internalStarts = System.nanoTime();
									int count= overlappingPaths.containsKey(p) ? overlappingPaths.get(p) + 1: 1;
									//System.err.print("\t\t\tGoing to add -- time to find in map("+(System.nanoTime()- internalStarts)+")-- " + p.toGene() + " to overlappers path, current map size "+ overlappingPaths.size() + ", supporting read count " + count);
									internalStarts = System.nanoTime();
									overlappingPaths.put(p, count);
									
									Alignments pExon = p.toGene().getSingleOverlappingExon(new Alignments(read.getChr(), read.getEnd()-1, read.getEnd()+1));
									AlignmentsPair exonPair = new AlignmentsPair(new Alignments(exon.getChromosome(), exon.getStart(), exon.getEnd()), pExon);
									Alignments pairedReadAlignment = new Alignments(read.getChr(), read.getAlignmentStart(), read.getAlignmentEnd());
									
									if (!pairedReadsForOverlappingPathsByExonPair.containsKey(p)) {
										Map<AlignmentsPair, Set<Alignments>> m = new HashMap<AlignmentsPair, Set<Alignments>>();
										pairedReadsForOverlappingPathsByExonPair.put(p, m);
									}
									
									Map<AlignmentsPair, Set<Alignments>> exonPairsMappedToReadAlignments = pairedReadsForOverlappingPathsByExonPair.get(p);
									if (!exonPairsMappedToReadAlignments.containsKey(exonPair)) {
										Set<Alignments> s = new TreeSet<Alignments>();
										exonPairsMappedToReadAlignments.put(exonPair, s);
									}
									exonPairsMappedToReadAlignments.get(exonPair).add(pairedReadAlignment);
									
									//System.err.println(" size after putting back "+ overlappingPaths.size()+ " put back took " + (System.nanoTime()- internalStarts));
								}
							}
							counter++;
							//System.err.println("\t\t\tRead handling took " + (System.nanoTime() - startRead));
							totalReadTime += (System.nanoTime() - startRead);
						}
					}
					
					startPartT = System.nanoTime();
					//System.err.println("\tPaths 1 has "+overlappingPaths.size()+" Paths 2 has "+allPaths.size());
					
				}
				//System.err.println("\tFinished exons time: " + (System.nanoTime() - exonTime));
				long connectTime = System.nanoTime();
				for(Path connectedPath: overlappingPaths.keySet()){
					
					Map<AlignmentsPair, Set<Alignments>> exonPairsMappedToReadAlignments = pairedReadsForOverlappingPathsByExonPair.get(connectedPath);
					
					for(Path path1: allPaths){						
						//if(!connectedPath.equals(path1) && !connectedPath.overlaps(path1) && overlappingPaths.get(connectedPath) >= MIN_PAIRED_END_SUPPORT){
						if(!connectedPath.equals(path1) && !connectedPath.overlaps(path1)){
							int pairedCounts = 	overlappingPaths.get(connectedPath);
							//connect the old and new path
							BubbleEdge edge=getEdge(path1, connectedPath, orientation, pairedCounts, graph);
							//System.err.println("\tGoing to add edge (just got it) " + edge.getConnection().toUCSC()+" between " + edge.getLeftNode().toUCSC() + " and " + edge.getRightNode().toUCSC() + " -` "+overlappingPaths.get(connectedPath));														
							//TODO Make a global count map and update here
							/*
							int connectionCount=0;
							if(edge!=null){
								LightweightGenomicAnnotation connection=edge.getConnection();
								if(connectionCounts.containsKey(connection)){
									connectionCount=connectionCounts.get(connection);
								}
								connectionCount++;
								connectionCounts.put(connection, connectionCount);
							}
							else{System.err.println("Skipping a null edge "+path1.toGene().getAlignment().toUCSC()+" "+connectedPath.toGene().getAlignment().toUCSC());}
							*/
							//NOTE: Since we are using paths rather than the graph, we will add connections 
							// between all paths that were disconnected. This introduces alternative isoforms that
							// are only supported by paired end data, most likely artifacts. Should have a clan up step.
							if(edge != null && pairedCounts > MIN_PAIRED_END_SUPPORT) {
								BubbleEdge old = graph.findEdge(edge.getLeftNode(), edge.getRightNode());
								//int counts = connectionCounts.get(edge.getConnection());
								boolean couldAdd = false;
								if(old == null) {
									Path leftPath;
									Path rightPath;
									int compare = path1.compareTo(connectedPath);
									if (compare > 0) {
										rightPath = path1;
										leftPath = connectedPath;
									} else {
										leftPath = path1;
										rightPath = connectedPath;
									}
									// RefSeqGene leftPathAsGene = leftPath.toGene();
									// RefSeqGene rightPathAsGene = rightPath.toGene();
									RefSeqGene leftPathAsGene = pathsToGenes.get(leftPath);
									RefSeqGene rightPathAsGene = pathsToGenes.get(rightPath);
									
									// NOTE: MAYBE REWORK pairedReadsForConnectedPath TO WORK ONE CONNECTION SITE AT A TIME, THEN WE CAN USE GAMMA_P 
									
									for (AlignmentsPair exonPair : exonPairsMappedToReadAlignments.keySet()) {
										Alignments connectedPairExon = exonPair.getAlignment2();
										Alignments path1Exon = exonPair.getAlignment1();
										
										if (!pathsToGenes.get(path1).overlapsExon(path1Exon))
											continue;
										
										Set<Alignments> pairedEndsForExonPair = exonPairsMappedToReadAlignments.get(exonPair);
										
										// HERE
									}
									
									boolean shouldAddEdge = false;
									/*
									for (Alignments pairedRead : pairedReadsForConnectedPath) {
										Alignments pairedReadStart = new Alignments(pairedRead.getChr(), pairedRead.getStart(), pairedRead.getStart() + 1);
										Alignments pairedReadEnd = new Alignments(pairedRead.getChr(), pairedRead.getEnd(), pairedRead.getEnd() + 1);
										
										if (leftPathAsGene.hasExon(pairedReadStart) && rightPathAsGene.hasExon(pairedReadEnd)) {
											int insertSize = AddEdgesUsingPairedEnds.estimateInsertSize(pairedRead, leftPathAsGene, rightPathAsGene);
											double cumulativeProbability = pairedInsertDistribution.getCummulativeProbability((double)insertSize);
											
											System.out.println("insert size: " + insertSize + " ::: " + cumulativeProbability);
											
											// STILL CHECKING THIS..
											
											if (cumulativeProbability >= MIN_PROBABILITY_IN_PAIRED_INSERT_DISTRIBUTION_TO_USE_WHEN_ADDING_EDGE && cumulativeProbability <= MAX_PROBABILITY_IN_PAIRED_INSERT_DISTRIBUTION_TO_USE_WHEN_ADDING_EDGE) {
												shouldAddEdge = true;
												break;
											}
											
										}
									}
									*/
									// shouldAddEdge = true;
									if (shouldAddEdge)
										couldAdd=graph.addEdge(edge, edge.getLeftNode(), edge.getRightNode());
									
									//] System.out.println("Adding edge to graph ::: " + connectedPath.toString() + " ::: " + path1.toString());
									
									//edge.setPairedEndCounts(counts);
								} else {
									//old.setPairedEndCounts(counts);
									couldAdd = true;
								}
							}
						}
						//System.err.println("Are we here??");
					}
				}

				
				long endTime=System.nanoTime();
				//System.err.println("\tTime for path: " + (endTime-startTime) + " for connecting edges " + (endTime - connectTime));
			}
			//System.err.println();
			i++;
		}
		pairedData.resetTreeCache();
	}

	private static boolean thereAreEnoughPairedEndsForConnection(Map<AlignmentsPair, Set<Alignments>> exonPairsMappedToReadAlignments, GenesInExpressionBins<PairedEndSignificanceInExpressionBin> genesInBins, RefSeqGene leftPath, RefSeqGene rightPath) throws IOException{
		for (AlignmentsPair exonPair : exonPairsMappedToReadAlignments.keySet()) {
			Set<Alignments> pairedEndsForExonPair = exonPairsMappedToReadAlignments.get(exonPair);
			
			Set<Alignments> validPairedEndsForExonPair = new TreeSet<Alignments>();
			int connectionSiteStart = Integer.MAX_VALUE;
			int connectionSiteEnd = -Integer.MAX_VALUE;
			for (Alignments pairedEnd : pairedEndsForExonPair) {
				if (leftPath.overlapsExon(new Alignments(pairedEnd.getChr(), pairedEnd.getStart(), pairedEnd.getStart() + 1)) && rightPath.overlapsExon(new Alignments(pairedEnd.getChr(), pairedEnd.getEnd(), pairedEnd.getEnd() + 1))) {
					connectionSiteStart = Math.min(connectionSiteStart, pairedEnd.getStart());
					connectionSiteEnd = Math.max(connectionSiteEnd, pairedEnd.getEnd());
					
					validPairedEndsForExonPair.add(pairedEnd);
				}
			}
			
			if (validPairedEndsForExonPair.size() == 0)
				continue;
			
			int numOfPairedEndsAtConnectionSite = validPairedEndsForExonPair.size();
			
			int lengthOfConnectionSite;
			if (leftPath.hasExon(new Alignments(leftPath.getChr(), connectionSiteStart, connectionSiteStart + 1))) {
				int leftPathLength = leftPath.trimAbsolute(connectionSiteStart, leftPath.getEnd()).getSize();
				int rightPathLength = rightPath.trimAbsolute(rightPath.getStart(), connectionSiteEnd).getSize();
				lengthOfConnectionSite = leftPathLength + rightPathLength;
			} else {
				int leftPathLength = leftPath.trimAbsolute(leftPath.getStart(), connectionSiteEnd).getSize();
				int rightPathLength = rightPath.trimAbsolute(connectionSiteStart, rightPath.getEnd()).getSize();
				lengthOfConnectionSite = leftPathLength + rightPathLength;
			}
	
			double rpkmOfLeftPath = leftPath.getRPKM();
			double rpkmOfRightPath = rightPath.getRPKM();
			
			double rpkmOfBoth = Math.min(rpkmOfLeftPath, rpkmOfRightPath); // THIS IS ONE WAY OF TAKING THE RPKM OF BOTH TRANSCRIPTS IN THE CONNECTION; TO USE THE MEAN OF THE TWO RPKM'S, COMMENT OUT THIS LINE AND USE THE BELOW LINE
			//double rpkmOfBoth = (rpkmOfLeftPath + rpkmOfRightPath) / 2;
			
			PairedEndSignificanceInExpressionBin pesieb = genesInBins.getValueForExpressionUsingNearbyExpressionBinsIfNeeded(rpkmOfBoth);
			double gammaPForExpressionBin = pesieb.getGammaP();
			double expectedNumOfPairedEnds = gammaPForExpressionBin * lengthOfConnectionSite;
			
			if (numOfPairedEndsAtConnectionSite >= expectedNumOfPairedEnds) {
				return true;
			} else {
				double poisson = AlignmentDataModelStats.poisson(numOfPairedEndsAtConnectionSite, expectedNumOfPairedEnds);
				return poisson >= MIN_POISSON_FOR_CONNECTING_PATHS;
			}
		}
		
		return false;
	}
	
	private static int estimateInsertSize(Alignments read, RefSeqGene leftPathAsGene, RefSeqGene rightPathAsGene) {
		RefSeqGene leftPathInsert = leftPathAsGene.trimAbsolute(read.getStart(), leftPathAsGene.getEnd());
		RefSeqGene rightPathInsert = rightPathAsGene.trimAbsolute(rightPathAsGene.getStart(), read.getEnd() + 1);
		
		/**
		if (leftPathInsert == null)
			System.out.println("l: " + leftPathAsGene.toBED() + " ::: " + read.toString());
		if (rightPathInsert == null)
			System.out.println("r: " + rightPathAsGene.toBED() + " ::: " + read.toString());
		*/
		
		if (leftPathInsert == null || rightPathInsert == null)
			return -1;
		
		int leftPathInsertSize = leftPathInsert.getSize();
		int rightPathInsertSize = rightPathInsert.getSize();
		
		return (leftPathInsertSize + rightPathInsertSize);
	}
	

	private static Stack<Path> getOverlappingPaths(Iterator<Node<Collection<Path>>> overlappers) {
		Stack<Path> rtrn=new Stack<Path>();
		
		while(overlappers.hasNext()){
			Collection<Path> over=overlappers.next().getValue();
			rtrn.addAll(over);
		}
		
		return rtrn;
	}

	private static Collection<Path> getOverlappingPaths(Iterator<Node<Collection<Path>>> overlappers, Alignment read) {
		Collection<Path> rtrn=new ArrayList<Path>();
		
		while(overlappers.hasNext()){
			Collection<Path> paths=overlappers.next().getValue();
			for(Path path: paths){
				//check if the path overlaps either the left or right of the read
				boolean overlapsEnds=overlapsEnd(path, read);
				
				//if so add to the collection
				if(overlapsEnds){rtrn.add(path);}
				//else skip it
			}
		}
		
		return rtrn;
	}

	private static boolean overlapsEnd(Path path, Alignment read) {
		if(path.overlapsLeft(read)){return true;}
		if(path.overlapsRight(read)){return true;}
		return false;
	}

	private static BubbleEdge getEdge(Path path, Path connectedPath, String orientation, Integer pairedCount, ChromosomeWithBubbles2 graph) {
		//Need to check which node comes first
		BubbleEdge edge = null;
		double count1=0;
		double count2=0;
		LightweightGenomicAnnotation v1 = null;
		LightweightGenomicAnnotation v2 =  null;
		
		
		
		int compare=path.compareTo(connectedPath);
		if(compare>0){
			//then we need to put connect first then other
			//MG swapped v2 with v1 05/25/2010
			v2=path.getStartVertex();
			v1=connectedPath.getEndVertex();
			count2=path.getNodeCount(v2);
			count1=connectedPath.getNodeCount(v1);
		}
		else{
			v1=path.getEndVertex();
			v2=connectedPath.getStartVertex();
			count1=path.getNodeCount(v1);
			count2=connectedPath.getNodeCount(v2);
		}
		//else{System.err.println("We never get in for "+path.toGene().getAlignment().toUCSC()+" "+connectedPath.toGene().getAlignment().toUCSC());}
		
		//if(v1 != null & v2 != null) {
			//Alignments connection=new Alignments(connectedPath.getChromosome(), Math.min(connectedPath.getEnd(), path.getEnd()), Math.max(connectedPath.getStart(), path.getStart()));
		Alignments connection = new Alignments(connectedPath.getChromosome(), v1.getEnd(), v2.getStart()); 

		v1.setOrientation(orientation);
		v2.setOrientation(orientation);
		connection.setOrientation(orientation);
			
		edge = "-".equals(orientation) 
			? graph.createEdge(connection, v2, v1, 0, EdgeSourceType.PAIRED)
			: graph.createEdge(connection, v1, v2, 0, EdgeSourceType.PAIRED);
			
		edge.setNodeCount(v1, count1);
		edge.setNodeCount(v2, count2);
		edge.setPairedEndCounts(pairedCount.doubleValue());
			//System.out.println("Pair_"+connection.toUCSC()+"_"+count + " connected " + v1.toUCSC() + " and " + v2.toUCSC() + " in " + orientation + " orientation");
		//}
		return edge;
	}

	private static Collection<Path> toCollection( IntervalTree<Collection<Path>> pathTree) {
		Collection<Path> rtrn=new TreeSet<Path>();
		
		Iterator<Node<Collection<Path>>> iter=pathTree.iterator();
		while(iter.hasNext()){
			Collection<Path> paths=iter.next().getValue();
			for(Path path: paths){rtrn.add(path);}
		}
		
		return rtrn;
	}

	private static void joinPaths(ChromosomeWithBubbles2 graph, IntervalTree<Collection<Path>> pathTree,AlignmentDataModelStats pairedData, String chr, int start, int end, String orientation) throws IOException {

		CloseableIterator<Alignment> iter=pairedData.getAlignmentsOverlappingRegion(new Alignments(chr, start, end));
		while(iter.hasNext()){
			Alignment read=iter.next();
			Iterator<Node<Collection<Path>>> overlappers=pathTree.overlappers(read.getAlignmentStart(), read.getAlignmentEnd());
			addToGraph(graph, read, overlappers, orientation, pairedData);
		}
				
	}

	private static void addToGraph(ChromosomeWithBubbles2 graph, Alignment read, Iterator<Node<Collection<Path>>> overlappers, String orientation, AlignmentDataModelStats  pairedData) throws IOException {
		Collection<Path> right=new TreeSet<Path>();
		Collection<Path> left=new TreeSet<Path>();
		while(overlappers.hasNext()){
			Collection<Path> paths=overlappers.next().getValue();
			for(Path path: paths){
				if(path.overlapsLeft(read)){left.add(path);}
				if(path.overlapsRight(read)){right.add(path);}
			}
		}
		
		for(Path rightPath: right){
			for(Path leftPath: left){
				if(!rightPath.equals(leftPath) && !rightPath.overlaps(leftPath)){
					LightweightGenomicAnnotation v1=leftPath.getEndVertex();
					LightweightGenomicAnnotation v2=rightPath.getStartVertex();
					Alignments connection=new Alignments(v1.getChromosome(), v1.getEnd(), v2.getStart());
					connection.setOrientation(orientation);
					//Lets now evaluate how significant it is the number of reads that connect the nodes connected by this one read.
					//1) Get the actual nodes connected by the read
					Iterator<LightweightGenomicAnnotation>leftPathOverlappinNodes = leftPath.getOverlappingNodes(read).iterator();
					Iterator<LightweightGenomicAnnotation>rightPathOverlappinNodes = rightPath.getOverlappingNodes(read).iterator();
					//2) For each pair of nodes in each subgraph, find the number of reads connecting them
					double [] scanData = null;
					while(scanData == null && leftPathOverlappinNodes.hasNext()) {
						LightweightGenomicAnnotation leftPathNode = leftPathOverlappinNodes.next();
						while(scanData == null && rightPathOverlappinNodes.hasNext()) {
							LightweightGenomicAnnotation rightPathNode = rightPathOverlappinNodes.next();
							List<LightweightGenomicAnnotation> leftPlusRight = new ArrayList<LightweightGenomicAnnotation>();
							leftPlusRight.add(leftPathNode);
							leftPlusRight.add(rightPathNode);
							RefSeqGene leftRightNode = new RefSeqGene(leftPlusRight);
							double[] pairScanData = pairedData.scanPRate(leftRightNode, 0);
							if(pairScanData[0] < fwer) {
								scanData = pairScanData;
							}
						}
					}
					
					if(scanData != null ) {
						boolean b=false;
						if(orientation.equalsIgnoreCase("+")){
							BubbleEdge edge=new BubbleEdge(connection, v1, v2, 1, EdgeSourceType.PAIRED);
							edge.setCount((float)scanData[2]);
							b=graph.addEdge(edge, v1, v2);
							//if(b){System.out.println("+connection: " + connection);}
						}
						else if(orientation.equalsIgnoreCase("-")){
							BubbleEdge edge=new BubbleEdge(connection, v2, v1, 1, EdgeSourceType.PAIRED);
							edge.setCount((float)scanData[2]);
							b=graph.addEdge(edge, v2, v1);
							//if(b){System.out.println("-connection: " + connection);}
						}
					}
				}
					
			}
		}
	}

	private static void addPairsToPaths(IntervalTree<Collection<Path>> pathTree,AlignmentDataModelStats pairedData, String chr, int start, int end) throws IOException {
				
		CloseableIterator<Alignment> iter=pairedData.getAlignmentsOverlappingRegion(new Alignments(chr, start, end));
		while(iter.hasNext()){
			Alignment read=iter.next();
			Iterator<Node<Collection<Path>>> overlappers=pathTree.overlappers(read.getAlignmentStart(), read.getAlignmentEnd());
			while(overlappers.hasNext()){
				Collection<Path> overlappingPaths=overlappers.next().getValue();
				//if they are within a path add the paired end to the path
				for(Path overlappingPath: overlappingPaths){
					boolean isAdded=overlappingPath.addPairedEnd(read);
				}
				//if(isAdded){System.out.println(read.getChromosome()+"\t"+read.getAlignmentStart()+"\t"+read.getAlignmentEnd());}
			}
		}
		
	}

	private static IntervalTree<Collection<Path>> makePathTree(Collection<Path> paths, String orientation) {
		IntervalTree<Collection<Path>> tree=new IntervalTree<Collection<Path>>();
		
		for(Path path: paths){
			if(path.getOrientation().equalsIgnoreCase(orientation) || path.getOrientation().equalsIgnoreCase("*")){
				Node<Collection<Path>> currentNode=tree.find(path.getStart(), path.getEnd());
				Collection<Path> c=new TreeSet<Path>();
				if(currentNode!=null){c=currentNode.getValue();}
				c.add(path);
				tree.put(path.getStart(), path.getEnd(), c);
			}
		}
		
		return tree;
	}
	
	private static IntervalTree<Collection<Path>> makePathTreeByExon(Collection<Path> paths, String orientation) {
		IntervalTree<Collection<Path>> tree=new IntervalTree<Collection<Path>>();
		
		for(Path path: paths){
			if(path.getOrientation().equalsIgnoreCase(orientation) || path.getOrientation().equalsIgnoreCase("*")){
				Collection<LightweightGenomicAnnotation> exons=path.getSortedNodes();
				for(LightweightGenomicAnnotation exon: exons){
					Node<Collection<Path>> currentNode=tree.find(exon.getStart(), exon.getEnd());
					Collection<Path> pathSet=new TreeSet();
					if(currentNode!=null){
						pathSet=currentNode.getValue();
					}
					pathSet.add(path);
					tree.put(exon.getStart(), exon.getEnd(), pathSet);
				}
				
				//Node<Collection<Path>> currentNode=tree.find(path.getStart(), path.getEnd());
				//Collection<Path> c=new TreeSet<Path>();
				//if(currentNode!=null){c=currentNode.getValue();}
				//c.add(path);
				//tree.put(path.getStart(), path.getEnd(), c);
			}
		}
		
		return tree;
	}
	
	private static IntervalTree<Collection<Path>> makePathTree(Collection<Path> paths) {
		IntervalTree<Collection<Path>> tree=new IntervalTree<Collection<Path>>();
		
		for(Path path: paths){
			Node<Collection<Path>> currentNode=tree.find(path.getStart(), path.getEnd());
			Collection<Path> c=new TreeSet<Path>();
			if(currentNode!=null){c=currentNode.getValue();}
			c.add(path);
			tree.put(path.getStart(), path.getEnd(), c);
		}
		
		return tree;
	}

	
	//TODO 1) Why are there more than 2? What are these?
	//TODO 2) Make sure isoforms are handled and merged in IntervalTree
	//TODO 3) Make sure that conflicting strands get split up and merged in stranded direction
	private static Path mergePaths(Collection<Path> spanningPaths) {
		if(spanningPaths.size()>2){
			for(Path path: spanningPaths){
				System.out.print(path.toGene()+"\t");
			}
			System.out.println();
		}

		Collection<String> orientations=new ArrayList<String>();
		for(Path path: spanningPaths){orientations.add((path.getOrientation()));}
		if(inConflict(orientations)){System.err.println("Conflict"); return null;}
		
		Path mergedPath=new Path();
		for(Path path: spanningPaths){
			mergedPath.setSpliceWeight(path.getSpliceWeight());
			mergedPath.addEdges(path.getEdges());
			mergedPath.addPairedEndEdges(path.getPairedEndEdges());		
		}
		return mergedPath;
	}

	private static boolean inConflict(Collection<String> orientations) {
		String current="*";
		for(String strand: orientations){
			if(current.equalsIgnoreCase("*")|| strand.equalsIgnoreCase("*")|| strand.equalsIgnoreCase(current)){
				current=updateStrand(current, strand);
			}
			else{return true;}
		}
		return false;
	}

	private static String updateStrand(String current, String strand) {
		if(strand.equalsIgnoreCase("*")){return current;}
		if(current.equalsIgnoreCase("*")){return strand;}
		return strand;
	}

	

	

	
}
