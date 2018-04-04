package broad.pda.seq.graph;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.broad.igv.sam.Alignment;

import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.datastructures.ReversibleIterator;
import broad.core.math.EmpiricalDistribution;
import broad.core.sequence.Sequence;
import broad.core.util.CollapseByIntersection;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.graph.ChromosomeWithBubbles2.BubbleEdge;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.AlignmentDataModelStats;
import broad.pda.seq.segmentation.CannonicalSpliceFilter;
import broad.pda.seq.segmentation.ReadFilter;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraDistance;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.util.Pair;

public class ChromosomeWithBubbles2 extends DirectedSparseGraph<LightweightGenomicAnnotation, BubbleEdge> {

	private static final long serialVersionUID = -6733460455807658588L;
	private static final int MAX_GENE_PATHS = 350;
	private static final String quote="\"";

	private IntervalTree<LightweightGenomicAnnotation> vertices;
	private  Map<LightweightGenomicAnnotation, Double> vertexCounts;
	private int chromosomeLength;
	private IntervalTree<BubbleEdge> edgeTree;
	private DijkstraDistance<LightweightGenomicAnnotation, BubbleEdge> distanceAlgorithm; //TODO: This is not currently taking into account edge types nor edge counts Need a filter
	private IntervalTree<Alignments> gapTree;
	private static final List<EdgeSourceType> defaultValidTypes = new ArrayList<EdgeSourceType>();
	//private static final List<EdgeSourceType> defaultTypesForDegreeCounting = new ArrayList<EdgeSourceType>();
	private static final double scaleFactor = 0.0025;
	private String name;
	private int start;
	private int end;
	private double lambda;
	private double numberOfReads;
	private double numberOfMarkers;
	private double rpkmConstant;
	private Map<Alignments, Double> localRate;
	
	private Map<LightweightGenomicAnnotation, Double> intronScores;
	//private Map<? extends LightweightGenomicAnnotation, Double> exonScores; //TODO Consider putting this back in later
	private EmpiricalDistribution pairedDist;
	private double spliceWeight=1;
	
	static {
		defaultValidTypes.add(EdgeSourceType.PAIRED);
		defaultValidTypes.add(EdgeSourceType.SPLICED);
		//defaultTypesForDegreeCounting.add(EdgeSourceType.SPLICED);
	}
	public ChromosomeWithBubbles2 ( String name) {
		this(name, 0, Integer.MAX_VALUE);
	}
	
	public ChromosomeWithBubbles2(String name, int start, int end) {
		super();
		this.name = name != null ? name.intern() : "noname".intern();
		this.vertexCounts=new TreeMap<LightweightGenomicAnnotation, Double>();
		vertices = new IntervalTree<LightweightGenomicAnnotation>();
		edgeTree = new IntervalTree<BubbleEdge>();
		//exonScores= new TreeMap<LightweightGenomicAnnotation,Double>();
		distanceAlgorithm = new DijkstraShortestPath<LightweightGenomicAnnotation, BubbleEdge>(this, false);
		this.start = start;
		this.end   = end;
		//System.err.println("Initialized graph");
	}
	
	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}

	/**
	 * Constructs a chromosome with bubbles with a given set of nodes and connections between them. The graph will keep track
	 * of how many times a given edge has been added. Nodes in this case are though as exons of putative transcripts.
	 * @param name - Name  of the chromosome.
	 * @param vertices - Collection of nodes, note than duplicate nodes will be ignored
	 * @param edges - Collection of edges connecting nodes. We assume a close-open( [start, end) ) specification of the annotation, 
	 * 			      thus edge connects two nodes if the start of the edge is node1.getEnd() and the end of the edge is node2.getStart() - 1. 
	 * 				  <b>Edges must be sorted</b>. 
	 * @throws IllegalArgumentException if a node or and edge are on a different chromsome than the chromosomeWithBubbles instance.
	 */
	public ChromosomeWithBubbles2 (String name, Collection<? extends LightweightGenomicAnnotation> nodes, Collection<? extends LightweightGenomicAnnotation> edges, Map<? extends LightweightGenomicAnnotation, Double> edgeCounts, Map<? extends LightweightGenomicAnnotation, Double> nodeCounts, double lambda, double numMarkers, double numReads, int minSplices) throws IllegalArgumentException{
		this(name);
		this.lambda=lambda;
		this.numberOfMarkers=numMarkers;
		this.numberOfReads=numReads;
		int counter=0;
		this.vertexCounts=new TreeMap<LightweightGenomicAnnotation, Double>();
		//this.exonScores=nodeCounts;
		//The hack below is to overcome the imposibility of adding a new item to a map<? extends LightweightGenomicAnnotation
		intronScores = new TreeMap<LightweightGenomicAnnotation, Double>(); 
		for(LightweightGenomicAnnotation intron: edgeCounts.keySet()) {
			intronScores.put(intron, edgeCounts.get(intron));
		}
		
		//System.err.println("starting graph construction");
		
		for(LightweightGenomicAnnotation node : nodes) {
			if(! name.equals(node.getChromosome())) {
				throw new IllegalArgumentException("Trying to add node " + node.toUCSC() + " to a chromosome with bubbles " + name + " chromosome must be the same");
			}
			counter++;
			if(nodeCounts!=null){vertexCounts.put(node, nodeCounts.get(node));}
			//if(counter%1 ==0){System.err.println(node);}
			addVertex(node);
		}
		//System.err.println("gone through all nodes");
		
		for(LightweightGenomicAnnotation edge : edges) {

			if(! name.equals(edge.getChromosome())) {
				throw new IllegalArgumentException("Trying to add edge " + edge.toUCSC() + " to a chromosome with bubbles " + name + " chromosome must be the same");
			}
			
			double edgeCount = getSpliceCount(edge);
			
			//System.err.println("Edge " + edge.toUCSC() + " intronScores is null? " +( intronScores==null )+ " if not null, does it contain edge? " + (intronScores != null && intronScores.containsKey(edge) ? intronScores.get(edge) : "edge not in map") );
			if(intronScores != null && !intronScores.isEmpty() && edgeCount <= minSplices) {
				continue;
			}
			
			List<LightweightGenomicAnnotation> lefts  = findAbuttingLeftVertices(edge);
			if(lefts.isEmpty()){
				//throw new IllegalStateException("Can't find left node for " + edge.toUCSC());
				//System.err.println("WARN: Can't find left node for " + edge.toUCSC());
				continue;
			}
			List<LightweightGenomicAnnotation>  rights = findAbuttingRightVertices(edge);
			if(rights.isEmpty()){
				//throw new IllegalStateException("Can't find rigth node for " + edge.toUCSC());
				//System.err.println("WARN: Can't find rigth node for " + edge.toUCSC());
				continue;
			}
			

			for(LightweightGenomicAnnotation left : lefts) {
				for(LightweightGenomicAnnotation right : rights) {
					//System.err.println("Adding edge " + edge.toUCSC() + " between " + left.toUCSC() + " and " + right.toUCSC());
					BubbleEdge e = edge.inReversedOrientation() ? findEdge(right, left) : findEdge(left, right); 
					boolean couldAdd = false;
					if(e == null) {
						if(edge.inReversedOrientation()) {
							e = new BubbleEdge(edge, right, left, new Double(getSpliceCount(edge)).floatValue()); 
							e.setSpliceEdgeCount(edgeCount);
							couldAdd = addEdge(e, right, left);
							e.setNodeCount(right, getCount(right));
							e.setNodeCount(left, getCount(left));
							
						} else {
							e = new BubbleEdge(edge, left, right, new Double(edgeCount).floatValue());
							e.setSpliceEdgeCount(edgeCount);
							couldAdd = addEdge(e, left, right);
							e.setNodeCount(right, getCount(right));
							e.setNodeCount(left, getCount(left));
						}
					} else {
						System.err.println("Something is not right, edge " + e.connection.toUCSC() + " was already in the graph. with count " + intronScores.get(e.connection));
						e.incrementCount();
						couldAdd= true;
					}
					//System.err.println("Added edge " + edge.toUCSC() + " between " + left.toUCSC() + " and " + right.toUCSC());
					if (!couldAdd) {
						System.err.println("Adding of edge " + edge + " between " + left + " and " + right + " failed");
					}
				}
			}
		}
		//System.err.println("gone through all edges");
		
	}

	public ChromosomeWithBubbles2(String name, Collection<? extends LightweightGenomicAnnotation> nodes, Collection<BubbleEdge> edges, boolean fromBE, double lambda, double numMarkers, double numReads){
		this(name);
		this.lambda=lambda;
		this.numberOfMarkers=numMarkers;
		
		this.numberOfReads=numReads;
		for(LightweightGenomicAnnotation node : nodes) {
			if(! name.equals(node.getChromosome())) {
				throw new IllegalArgumentException("Trying to add node " + node.toUCSC() + " to a chromosome with bubbles " + name + " chromosome must be the same");
			}
			addVertex(node);
		}
		
		for(BubbleEdge edge: edges){
			addEdge(edge, edge.getLeftNode(), edge.getRightNode());
		}
	}
	
	public ChromosomeWithBubbles2(String name, Map<? extends LightweightGenomicAnnotation, Double> nodes, Map<BubbleEdge, Double> edges, double lambda, double numMarkers, double numReads, double minEdgeSupport){
		this(name);
		this.lambda=lambda;
		this.numberOfMarkers=numMarkers;
		this.numberOfReads=numReads;
		//this.exonScores=nodes;
		this.intronScores=makeIntronScoresFromBE(edges);
		
		for(LightweightGenomicAnnotation node : nodes.keySet()) {
			if(! name.equals(node.getChromosome())) {
				throw new IllegalArgumentException("Trying to add node " + node.toUCSC() + " to a chromosome with bubbles " + name + " chromosome must be the same");
			}
			//node.setScore(nodes.get(node)); //TODO Orphans might not have scores
			addVertex(node);
			vertexCounts.put(node, nodes.get(node));
		}
		
		for(BubbleEdge edge: edges.keySet()){
			//edge.setCount(edges.get(edge));
			//System.err.println("Edge: " + edge.getConnection().toUCSC() + " type " + edge.getType() + " count " + edge.count + " map count " + edges.get(edge));
			if(edge.getType().ordinal() == EdgeSourceType.SPLICED.ordinal()){
				edge.setSpliceEdgeCount(this.getSpliceCount(edge.getConnection()));
			}
			if(edge.getType().ordinal() != EdgeSourceType.SPLICED.ordinal() || edge.getSplicedCounts() > minEdgeSupport) {
				edge.setParent(this);
				//double nodeCount=nodes.get(edge.getLeftNode())+nodes.get(edge.getRightNode());
				//edge.setNodeCount(nodeCount);
				//edge.setNodeCount(edge.getLeftNode(), this.getCount(edge.getLeftNode()));
				//edge.setNodeCount(edge.getRightNode(), this.getCount(edge.getRightNode()));
				addEdge(edge, edge.getLeftNode(), edge.getRightNode());
			}
		}
	}
	
	public ChromosomeWithBubbles2 (String name, Collection<? extends LightweightGenomicAnnotation> nodes, Collection<? extends LightweightGenomicAnnotation> edges, double lambda, double numMarkers, double numReads, int minSplices) throws IllegalArgumentException{
		this(name, nodes, edges, null, null, lambda, numMarkers, numReads, minSplices);
	}
	
	
	private Map<LightweightGenomicAnnotation, Double> makeIntronScoresFromBE(Map<BubbleEdge, Double> edges2) {
		Map<LightweightGenomicAnnotation, Double> rtrn=new TreeMap<LightweightGenomicAnnotation, Double>();
		
		for(BubbleEdge edge: edges2.keySet()){
			if(edge.getType().equals(EdgeSourceType.SPLICED)){
				rtrn.put(edge.getConnection(), edges2.get(edge));
			}
		}
		
		return rtrn;
	}


	public double getLambda(){return this.lambda;}
	
	public Number getDjistraDistance(LightweightGenomicAnnotation v1, LightweightGenomicAnnotation v2){
		return this.distanceAlgorithm.getDistance(v1, v2);
	}
	
	private List<LightweightGenomicAnnotation> findAbuttingRightVertices(LightweightGenomicAnnotation edge) {
		//Node<LightweightGenomicAnnotation> candidateNode = vertices.min(edge.getEnd() , edge.getEnd()+1);
		List<LightweightGenomicAnnotation> abuttingvertices = new ArrayList<LightweightGenomicAnnotation>();

		Iterator<Node<LightweightGenomicAnnotation>> iter=vertices.overlappers(edge.getEnd(), edge.getEnd()+1);
		while(iter.hasNext()){
			//test if any match exactly, if so set it
			LightweightGenomicAnnotation rightTemp=iter.next().getValue();
			if(rightTemp.getStart() == edge.getEnd() ){
				abuttingvertices.add(rightTemp);
			}
		}
		
		return abuttingvertices;
	}

	private List<LightweightGenomicAnnotation> findAbuttingLeftVertices(LightweightGenomicAnnotation edge) {
		//Node<LightweightGenomicAnnotation> candidateNode = vertices.max(edge.getStart() , edge.getStart()+1);
		List<LightweightGenomicAnnotation> abuttingvertices = new ArrayList<LightweightGenomicAnnotation>();

		Iterator<Node<LightweightGenomicAnnotation>> iter=vertices.overlappers(edge.getStart()-1, edge.getStart());
		
		while(iter.hasNext() ){
			//test if any match exactly, if so set it
			LightweightGenomicAnnotation leftTemp=iter.next().getValue();
			if(leftTemp.getEnd() == edge.getStart() ){
				abuttingvertices.add(leftTemp);
			}
		}
		
		return abuttingvertices;
	}

	public int getRelativeGenomeLength(){return chromosomeLength;}
	
	//return true if position satrts in gap but goes past gap
	public boolean isSpanningGap(int start, int end){
		boolean result = false;
		boolean isInAllIntrons = isInAllGaps(start, end);
		//System.err.print("Testing ("+start+","+end+")" + " is in all introns? " + isInAllIntrons + " ");
		if(isInAllIntrons) {
			Iterator<Node<BubbleEdge>> overlappers=edgeTree.overlappers(start, end);
			while(overlappers.hasNext() && !result){
				Node<BubbleEdge> next=overlappers.next();
				result = end > next.getEnd() || start < next.getStart();
				//System.err.print("["+next.getValue().getConnection().toUCSC()+"-"+result);
			}
		}
		//System.err.println(" ");
		return result;
	}
	
		
	//return true if position satrts in gap but goes past gap
	public boolean isGap(int start, int end){
		/*Iterator<Node<LightweightGenomicAnnotation>> overlappers=edges.overlappers(start, end);
		boolean isGap = overlappers.hasNext();
		while(overlappers.hasNext() && isGap){
			Node<LightweightGenomicAnnotation> next=overlappers.next();
			isGap = start >= next.getStart() && end<=next.getEnd();
		}
		return isGap;*/
		return false;
	}
	
	public void removeSpanningGapNodes() {
		List<LightweightGenomicAnnotation> vertexList = vertices.toList();
		for(LightweightGenomicAnnotation v : vertexList) {;
			if(isSpanningGap(v.getStart(), v.getEnd())) {
				System.err.println("Removing " + v.toUCSC());
				boolean couldRemove = removeVertex(v);
				if(!couldRemove) {System.err.println("Could not remove vertex " + v.toUCSC());} 
			}
		}
		
	}
	
	public boolean removeVertex(LightweightGenomicAnnotation vertex) {
		LightweightGenomicAnnotation removed = vertices.remove(vertex.getStart(), vertex.getEnd());
		boolean allGood = true;
		if(removed != null) {
			Collection<BubbleEdge> incidentEdges = new ArrayList<BubbleEdge>(getInEdges(vertex));
			for(BubbleEdge incidentEdge : incidentEdges) {
				//System.err.println("\tIn edge " + incidentEdge.getConnection().toUCSC());
				boolean couldRemove = removeEdge(incidentEdge);
				if(allGood) {allGood = couldRemove;}
				if(!allGood) {System.err.println("Could not remove in edge " + incidentEdge.getConnection().toUCSC());}
				
			}
			Collection<BubbleEdge> outEdges = new ArrayList<BubbleEdge>(getOutEdges(vertex));
			for(BubbleEdge outEdge : outEdges) {
				//System.err.println("\tOut edge " + outEdge.getConnection().toUCSC());
				boolean couldRemove = removeEdge(outEdge);
				if(allGood) {allGood = couldRemove;}
				if(!allGood) {System.err.println("Could not remove out edge " + outEdge.getConnection().toUCSC());}
			}
			boolean couldRemove = super.removeVertex(vertex);
			//System.err.println("Could remove " + vertex.toUCSC() + " from parent graph? " + couldRemove + " does parent graph contains verte? " + containsVertex(vertex));
			if(allGood) {allGood = couldRemove;}
		}
		
		return allGood;
	}
	
	/**
	 * Creates bubbles that connect the splice junctions specified in the data model. Nodes in the graphs (bubbles) are single base pairs 
	 * indicating the splice junction.
	 * @param model
	 * @param minAligmentGap
	 */
	public void loadBubbles(AlignmentDataModelStats model, Sequence chromosomeSequence, int minAligmentGap, int minNumIntrons, boolean filterSplice) throws IOException {
		chromosomeLength = model.chromosomeLength(name);
		//edges = makeGapTree(model, minAligmentGap);
		edges.clear();
		vertices.clear();
		distanceAlgorithm.reset();
		end = Math.min(end, chromosomeLength);
		
		Alignments chrAnnotation = new Alignments(name, 0, chromosomeLength);
		
		List<ReadFilter> filters = new ArrayList<ReadFilter>();
		if(filterSplice && chromosomeSequence != null && chromosomeSequence.getLength() >0) {
			filters.add(new CannonicalSpliceFilter(chromosomeSequence));
		}
		filters.add(new MinReadFilter(minAligmentGap));
		
		if(start > 0 || end < Integer.MAX_VALUE) {
			filters.add(new RegionReadFilter(start, end));
		}
				
		Map<Alignments, Integer> intronMap = model.getSplicedReads(chrAnnotation, filters, minNumIntrons);
		List<Alignments> introns = new ArrayList<Alignments>(intronMap.keySet());
		Collections.sort(introns, new Comparator<Alignments>() {
			public int compare(Alignments o1, Alignments o2) {
				return o1.getStart() == o2.getStart() ? o1.getEnd() - o2.getEnd() : o1.getStart() - o2.getStart();
			}
		});
		for(LightweightGenomicAnnotation intron : introns) {
			addEdgeAddingMissingVertices(intron, intronMap.get(intron));
		}
		
	}
	
	public void loadBubbles(AlignmentDataModel model) throws IOException {
		chromosomeLength = model.getChromosomeLengths().get(name);
		//edges = makeGapTree(model, minAligmentGap);
		edges.clear();
		vertices.clear();
		distanceAlgorithm.reset();
		end = Math.min(end, chromosomeLength);
		
		Alignments chrAnnotation = new Alignments(name, 0, chromosomeLength);
		
		Map<Alignments, Integer> intronMap = model.getSplicedReads(chrAnnotation);
		List<Alignments> introns = new ArrayList<Alignments>(intronMap.keySet());
		Collections.sort(introns, new Comparator<Alignments>() {
			public int compare(Alignments o1, Alignments o2) {
				return o1.getStart() == o2.getStart() ? o1.getEnd() - o2.getEnd() : o1.getStart() - o2.getStart();
			}
		});
		for(LightweightGenomicAnnotation intron : introns) {
			addEdgeAddingMissingVertices(intron, intronMap.get(intron));
		}
		
	}
	
	/**
	 * Creates bubbles that connect the splice junctions specified in the data model. Nodes in the graphs (bubbles) are single base pairs 
	 * indicating the splice junction.
	 * @param model
	 * @param minAligmentGap
	 */
	public void loadBubbles(AlignmentDataModel model, Sequence chromosomeSequence, int minAligmentGap, boolean spliceFilter) throws IOException{
		loadBubbles(new AlignmentDataModelStats(model), chromosomeSequence, minAligmentGap, 0, spliceFilter);
	}
	
	
	/**
	 * Creates bubbles that connect the splice junctions specified in the data model. Nodes in the graphs (bubbles) are single base pairs 
	 * indicating the splice junction.
	 * @param model
	 * @param minAligmentGap
	 */
	public void loadBubbles(AlignmentDataModelStats model, Sequence chromosomeSequence, int minAligmentGap, boolean spliceFilter) throws IOException
	{
		loadBubbles(model, chromosomeSequence, minAligmentGap, 0, spliceFilter);
	}
	

	
	public void addEdgeAddingMissingVertices(LightweightGenomicAnnotation intron, float count) {
		LightweightGenomicAnnotation start = new BasicLightweightAnnotation(intron.getChromosome(), intron.getStart() - 1, intron.getStart());//TODO: revise, should be [getStart()-1, getStart())
		LightweightGenomicAnnotation end   = new BasicLightweightAnnotation(intron.getChromosome(), intron.getEnd(), intron.getEnd() + 1);
		
		
		BubbleEdge edge = intron.inReversedOrientation() ? findEdge(end, start) : findEdge(start, end);
		if(edge != null) {
			edge.setCount(edge.getAllCounts() + count); //update count if duplicated edge is added.
		} else {
			if(!containsVertex(start)) {
				addVertex(start);
			} 
			
			if(!containsVertex(end) ){
				addVertex(end);
			}
			
			if(intron.inReversedOrientation()) {
				edge = new BubbleEdge(intron, end, start, count);
				edge.setParent(this);
				edge.setNodeCount(end, this.getCount(end));
				edge.setNodeCount(start, this.getCount(start));
				//edge.setSpliceEdgeCount(this.getSpliceCount(intron));//MODIFIED BY MANUEL 1/4/2011 This looked like a bug. The splice count should be the one passed in as the paramter PLEASE REVIEW THIS
				intronScores.put(edge.connection, (double)count);
				addEdge(edge, end, start);
			} else {
				edge = new BubbleEdge(intron, start, end, count);
				edge.setParent(this);
				edge.setNodeCount(start, this.getCount(start));
				edge.setNodeCount(end, this.getCount(end));
				//edge.setSpliceEdgeCount(this.getSpliceCount(intron));//MODIFIED BY MANUEL 1/4/2011 This looked like a bug. The splice count should be the one passed in as the paramter PLEASE REVIEW THIS
				addEdge(edge, start, end);
				intronScores.put(edge.connection, (double)count);
			}		
		}
	}
	
	public boolean addVertex(LightweightGenomicAnnotation v) {
		boolean success = super.addVertex(v);
		if(success) {
			vertices.put(v.getStart(), v.getEnd(), v);
		}
		return success;
	}
	
	public boolean addEdge(BubbleEdge e, LightweightGenomicAnnotation v1, LightweightGenomicAnnotation v2) {
		boolean success = super.addEdge(e, v1, v2);
		if(success) {
			edgeTree.put(e.getConnection().getStart(), e.getConnection().getEnd(),e);
			e.setParent(this);
		}
		
		return success;
	}
	
	public boolean addNewEdge(BubbleEdge e, LightweightGenomicAnnotation v1, LightweightGenomicAnnotation v2){
		//check if edge already exist
		boolean success = false;
		BubbleEdge previous=findEdge(v1, v2);
		if(previous==null){
			previous=findEdge(v1,v2);
		}
		
		if(previous==null){
			success = addEdge(e, v1, v2);
		}
	
		return previous != null || success;
	}

	private Number getOverlapperSetDistance(LightweightGenomicAnnotation v1, LightweightGenomicAnnotation v2) {
		Number d = distanceAlgorithm.getDistance(v1, v2);
		
		if(d == null) {
			Iterator<Node<LightweightGenomicAnnotation>> v1OverlapperNodeIt = vertices.overlappers(v1.getStart(), v1.getEnd());
			Iterator<Node<LightweightGenomicAnnotation>> v2OverlapperNodeIt = vertices.overlappers(v2.getStart(), v2.getEnd());
			
			while(d == null && v1OverlapperNodeIt.hasNext()) {
				LightweightGenomicAnnotation v1Overlapper = v1OverlapperNodeIt.next().getValue();
				while(d == null && v2OverlapperNodeIt.hasNext()) {
					LightweightGenomicAnnotation v2Overlapper = v2OverlapperNodeIt.next().getValue();
					d = distanceAlgorithm.getDistance(v1Overlapper, v2Overlapper);
				}
			}
		}
		return d;
	}

	public LightweightGenomicAnnotation findMaxReachableNode(LightweightGenomicAnnotation node, Collection<EdgeSourceType> validTypes, Collection<LightweightGenomicAnnotation> alreadyTraversed) {
		Collection<BubbleEdge> forwardEdges = getForwardEdges(node, 0, validTypes);
		LightweightGenomicAnnotation maxNode = node;
		alreadyTraversed.add(node);
		if(forwardEdges.size() > 0) {
			for(BubbleEdge e : forwardEdges) {
				LightweightGenomicAnnotation opposite = getOpposite(node, e);
				if(!alreadyTraversed.contains(opposite) ){
					LightweightGenomicAnnotation candidateMax = findMaxReachableNode(opposite, validTypes, alreadyTraversed);
					alreadyTraversed.add(opposite);
					maxNode = candidateMax.getEnd() > maxNode.getEnd() ? candidateMax : maxNode;
				}
			}
		}
		return maxNode;
	}
	
	public LightweightGenomicAnnotation findMinReachableNode(LightweightGenomicAnnotation node, Collection<EdgeSourceType> validTypes, Collection<LightweightGenomicAnnotation> alreadyTraversed) {
		Collection<BubbleEdge> backwardEdges = getBackwardEdges(node, 0, validTypes);
		LightweightGenomicAnnotation minNode = node;
		alreadyTraversed.add(node);
		if(backwardEdges.size() > 0) {
			//System.err.print("\n\t<node " + node.toUCSC() + " bes: " + backwardEdges.size()+">");
			for(BubbleEdge e : backwardEdges) {
				LightweightGenomicAnnotation opposite = getOpposite(node, e);
				if(!alreadyTraversed.contains(opposite) ){
					LightweightGenomicAnnotation candidateMin = findMinReachableNode(opposite, validTypes, alreadyTraversed);
					alreadyTraversed.add(opposite);
					minNode = candidateMin.getStart() < minNode.getStart() ? candidateMin : minNode;
				}
			}
		}
		return minNode;
	}

	public Collection<BubbleEdge> getForwardEdges(LightweightGenomicAnnotation vertex, double alpha, Collection<EdgeSourceType> validTypes) {
		Collection<BubbleEdge> forwardEdges = new ArrayList<BubbleEdge>();
		Collection<BubbleEdge> edges = getIncidentEdges(vertex);
		for(BubbleEdge e : edges) {
			if(getOpposite(vertex, e).getStart() > vertex.getEnd() && isEdgeValid(e,  alpha, validTypes)) {
				forwardEdges.add(e);
			}
		}

		return forwardEdges;			
	}
	
	public Collection<BubbleEdge> getBackwardEdges(LightweightGenomicAnnotation vertex, double alpha, Collection<EdgeSourceType> validTypes) {
		Collection<BubbleEdge> backwardEdges = new ArrayList<BubbleEdge>();
		Collection<BubbleEdge> edges = getIncidentEdges(vertex);
		for(BubbleEdge e : edges) {
			if(getOpposite(vertex, e).getEnd() < vertex.getStart() && isEdgeValid(e, alpha, validTypes)) {
				backwardEdges.add(e);
			}
		}

		return backwardEdges;			
	}

	public String getName() {
		return name;
	}
	
	public Collection<LightweightGenomicAnnotation> getVerticesWithForwardEdges(double alpha, Collection<EdgeSourceType> validEdgeTypes) {
		Collection<LightweightGenomicAnnotation> vertices = new ArrayList<LightweightGenomicAnnotation>();
		
		Collection<LightweightGenomicAnnotation> allVertices = getVertices();
		for(LightweightGenomicAnnotation v : allVertices) {
			if(!getForwardEdges(v, alpha, validEdgeTypes).isEmpty()) {
				vertices.add(v);
			}
		}
		return vertices;
	}
	
	public Collection<LightweightGenomicAnnotation> getSourceVertices(Collection<EdgeSourceType> validTypes) {
		Collection<LightweightGenomicAnnotation> vertices = getVertices();
		Collection<LightweightGenomicAnnotation> sources = new TreeSet<LightweightGenomicAnnotation>();
		for(LightweightGenomicAnnotation v: vertices) {
			//System.err.println(v.toUCSC()+" "+inOrientedDegree(v, validTypes)+" "+outOrientedDegree(v, validTypes) + " " + isSpanningGap(v.getStart(), v.getEnd()));
			//if(inOrientedDegree(v, validTypes) == 0 && outOrientedDegree(v, validTypes) > 0 && !isSpanningGap(v.getStart(), v.getEnd())) {
			if(inOrientedDegree(v, validTypes) == 0 && outOrientedDegree(v, validTypes) >= 0 && !isSpanningGap(v.getStart(), v.getEnd())) {
				//System.out.println("adding source");
				sources.add(v);
			}
		}
		
		return sources;
	}
	
	
	public int inOrientedDegree(LightweightGenomicAnnotation v, Collection<EdgeSourceType> validTypes) {
		int inDegree = 0;
		Collection<BubbleEdge> inEdges = getInEdges(v);
		for(BubbleEdge e : inEdges) {
			if(isEdgeValid(e, 1.1, validTypes)) {
				inDegree++;
			}
		}
		//System.err.println("#incident " + inDegree(v) + " #filtered indegree " + inDegree);
		return inDegree;
	}
	
	public int outOrientedDegree(LightweightGenomicAnnotation v, Collection<EdgeSourceType> validTypes) {
		int outDegree = 0;
		Collection<BubbleEdge> outEdges = getOutEdges(v);
		for(BubbleEdge e : outEdges) {
			if(isEdgeValid(e, 1.1, validTypes)) {
				outDegree++;
			}
		}	
		return outDegree;
	}
	
	public Collection<LightweightGenomicAnnotation> getSinkVertices() {
		return getSinkVertices(null);
	}
	
	public Collection<LightweightGenomicAnnotation> getSinkVertices( Collection<EdgeSourceType> validTypes) {
		Collection<LightweightGenomicAnnotation> vertices = getVertices();
		Collection<LightweightGenomicAnnotation> sinks = new TreeSet<LightweightGenomicAnnotation>();
		for(LightweightGenomicAnnotation v: vertices) {
			if(inOrientedDegree(v, validTypes) > 0 && outOrientedDegree(v, validTypes) == 0) {
				sinks.add(v);
			}
		}
		
		return sinks;
	}
	
	public Collection<BubbleEdge> getOrientedInEdges(LightweightGenomicAnnotation v) {
		Collection<BubbleEdge> inEdges = super.getInEdges(v);
		Collection<BubbleEdge> filteredInEdges = new LinkedList<BubbleEdge>();
		for(BubbleEdge e: inEdges) {
			if(!"*".equals(e.getOrientation())){
				filteredInEdges.add(e);
			}
		}
		
		return filteredInEdges;
	}
	
	public Collection<BubbleEdge> getOrientedOutEdges(LightweightGenomicAnnotation v) {
		Collection<BubbleEdge> outEdges = super.getOutEdges(v);
		Collection<BubbleEdge> filteredOutEdges = new LinkedList<BubbleEdge>();
		for(BubbleEdge e: outEdges) {
			if(!"*".equals(e.getOrientation())){
				filteredOutEdges.add(e);
			}
		}
		
		return filteredOutEdges;
	}
	
	public int inOrientedDegree(LightweightGenomicAnnotation v) {
		return getOrientedInEdges(v).size();
	}
	
	public int outOrientedDegree(LightweightGenomicAnnotation v) {
		return getOrientedOutEdges(v).size();
	}
	
	public boolean isInAllGaps(int start, int end) {
		if( gapTree == null || gapTree.isEmpty() ) {
			loadChromosomeGapTree();
		}
		
		return gapTree.overlappers(start, end).hasNext();
	}
	
	
	private Collection<BubbleEdge> filterValidEdges(Collection<BubbleEdge> edgeList) {
		Collection<BubbleEdge> filtered = new TreeSet<BubbleEdge>();
		for(BubbleEdge e : edgeList) {
			//System.out.println("Considering edge");
			if(isEdgeValid(e, 1.1, getValidTypes())) {
				//System.out.println("Edge is valid");
				filtered.add(e);
			}
		}
		
		return filtered;	
	}
	private Collection<BubbleEdge> getOutOrUnorientedValidEdges(LightweightGenomicAnnotation v) {
		return filterValidEdges(getOutOrUnorientedEdges(v));
		
	}
	
	private Collection<BubbleEdge> getInOrUnorientedValidEdges(LightweightGenomicAnnotation v) {
		return filterValidEdges(getInOrUnorientedEdges(v));
		
	}
	
	public Collection<BubbleEdge> getOutOrUnorientedEdges(LightweightGenomicAnnotation v) {
		Collection<BubbleEdge> rtrn = new TreeSet<BubbleEdge>();
		if(getOutEdges(v) == null) {
			System.err.println("V " + v.toUCSC() + " has NULL out edges! ");
		} else {
			//System.out.println("\tgetOutEdges(v): " + getOutEdges(v).size());
			rtrn.addAll(getOutEdges(v));
		}

		Collection<BubbleEdge> inEdges = getInEdges(v);
		for(BubbleEdge e : inEdges) {
			if("*".equals(e.getOrientation()) ){
					rtrn.add(e);
			}
		}
		return rtrn;
	}
	
	public Collection<BubbleEdge> getInOrUnorientedEdges(LightweightGenomicAnnotation v) {
		Collection<BubbleEdge> rtrn = new TreeSet<BubbleEdge>();
		if(getInEdges(v) == null) {
			System.err.println("V " + v.toUCSC() + " has NULL out edges! ");
		} else {
			rtrn.addAll(getInEdges(v));
		}

		Collection<BubbleEdge> outEdges = getOutEdges(v);
		for(BubbleEdge e : outEdges) {
			if("*".equals(e.getOrientation()) ){
					rtrn.add(e);
			}
		}
		return rtrn;
	}
	
	private void loadChromosomeGapTree() {
		//System.err.println("loading chromosome gap tree");
		Collection<Alignments> edgeAlignments = new ArrayList<Alignments>(edges.size());
		Iterator<BubbleEdge> edgeIt = edgeTree.valueIterator();
		//System.err.println("\tEdges:");
		while(edgeIt.hasNext()) {
			LightweightGenomicAnnotation e = edgeIt.next().getConnection();
			/*BubbleEdge be = findEdge(e);
			if(be.count>5){
				System.err.println("\t\t"+e.toUCSC());
			*/	
				edgeAlignments.add(new Alignments(e.getChromosome(), e.getStart(), e.getEnd()));
			//}
			
		}
		Set<Alignments> colapsedEdges = CollapseByIntersection.collapseByIntersection(edgeAlignments, true);
		//System.err.println("\tCollapsed edges:");
		gapTree = new IntervalTree<Alignments>();
		for(Alignments ce : colapsedEdges) {
			//System.err.println("\t\tAdding gap: " + ce.toUCSC());
			gapTree.put(ce.getStart(), ce.getEnd(), ce);
		}
		//System.err.println("Done loading chromosome gap tree");
	}
	
	private BubbleEdge findEdge(LightweightGenomicAnnotation e) {
		LightweightGenomicAnnotation es = new BasicLightweightAnnotation(name, e.getStart()-1, e.getStart());
		LightweightGenomicAnnotation ee = new BasicLightweightAnnotation(name, e.getEnd(), e.getEnd() + 1);
		
		return e.inReversedOrientation() ?  findEdge(ee, es) : findEdge(es, ee);
	}

	public Collection<RefSeqGene> getPaths(int start, int end) {
		return getPaths(start, end, 1, true);
	}
	
	private List<EdgeSourceType> getValidTypes() {
		return defaultValidTypes;
	
	}
	
	public Collection<RefSeqGene> getPaths(int start, int end, double alpha, boolean forward) {		
		return forward? getPaths(start, end, alpha, getValidTypes()) : getBackwardPaths(start, end, alpha, getValidTypes());
	}
	
	public Collection<RefSeqGene> getPaths(int start, int end, double alpha, Collection<EdgeSourceType> validTypes) {
		List<RefSeqGene> paths = new ArrayList<RefSeqGene>();
		
		Iterator<Node<LightweightGenomicAnnotation>> verticesOverlappingIt = vertices.overlappers(start, end);
		List<LightweightGenomicAnnotation> verticesTraversedList = new ArrayList<LightweightGenomicAnnotation>();
		
		while(verticesOverlappingIt.hasNext()) {
			LightweightGenomicAnnotation v= verticesOverlappingIt.next().getValue();
			//System.err.println("New Vertex " + v.toUCSC() + " traversed list " +  verticesTraversedList + " can be reached? " + canBeReachedFrom(v, verticesTraversedList));
			if(v.getEnd() != end && !canBeReachedFrom(v, verticesTraversedList)) {
				// block one is from start to v 
				Alignments firstBlock = new  Alignments(name, start, v.getEnd()); // Change to support vertices not overlapping edges MG
				// find other blocks by following each edge to the next splice
				Collection<BubbleEdge> edges = getForwardEdges(v, alpha, validTypes);
				for(BubbleEdge edge : edges) {
					LightweightGenomicAnnotation nextV = getOpposite(v, edge);//getOpposite(v, edge);
					Alignments orientedFirstBlock  = new Alignments(name, start, v.getEnd()); // Change to support vertices not overlapping edges MG
					orientedFirstBlock.setOrientation(edge.getOrientation());
					Collection<RefSeqGene> nextPaths =  getPaths(nextV.getStart(), nextV.getStart() + (end - start) - firstBlock.length(), alpha, validTypes); //Recursion
					// append blocks found following each edge to the nascent transcript.
					if(nextPaths.size() == 0) {
						Alignments lastBlock = new Alignments(name, nextV.getStart(), nextV.getStart() + (end - start) - firstBlock.length());
						lastBlock.setOrientation(edge.getOrientation());
						//if(!lastBlock.getOrientation().equals(orientedFirstBlock)) {
							//System.err.println("WARN: Making an inconsistent gene, exon " + orientedFirstBlock.toUCSC() + "("+orientedFirstBlock.getOrientation()+") wheras " + lastBlock.toUCSC() + "("+lastBlock.getOrientation()+")");
						//}
						List<Alignments> exons = new ArrayList<Alignments>();
						exons.add(orientedFirstBlock);
						exons.add(lastBlock);
						RefSeqGene prependedPath = new RefSeqGene(exons);
						prependedPath.setOrientation(orientedFirstBlock.getOrientation());
						paths.add(prependedPath);
					} else {
						for(RefSeqGene newPath : nextPaths ) {
							Collection<Alignments> newPathIntrons = newPath.getIntronSet()	;// Change to support vertices not overlapping edges MG
							if(newPathIntrons.size() == 0) {
								verticesTraversedList.add(new BasicLightweightAnnotation(name, newPath.getStart(), newPath.getStart()+1));
							} else  {
							
								Iterator<Alignments> newPathIntronIt = newPathIntrons.iterator();
								while(newPathIntronIt.hasNext()) {
									Alignments intron = newPathIntronIt.next();
									verticesTraversedList.add(new BasicLightweightAnnotation(name, intron.getStart()-1, intron.getStart()));// Change to support vertices not overlapping edges MG
									verticesTraversedList.add(new BasicLightweightAnnotation(name, intron.getEnd(), intron.getEnd() + 1));// Change to support vertices not overlapping edges MG
								}
							}
							List<Alignments> exons = new ArrayList<Alignments>();
							//if(!newPath.getOrientation().equals(orientedFirstBlock)) {
							//	System.err.println("WARN: Making an inconsistent gene, exon " + orientedFirstBlock.toUCSC() + "("+orientedFirstBlock.getOrientation()+") wheras " + newPath.toString() );
							//}							
							exons.add(orientedFirstBlock); 
							exons.addAll(newPath.getExonSet()); //Changed to support vertices not overlapping edges
							RefSeqGene prependedPath = new RefSeqGene(exons);
							prependedPath.setOrientation(orientedFirstBlock.getOrientation());
							paths.add(prependedPath);
						}
					}
				}

			}
		}
		if(paths.isEmpty()){
			paths.add(new RefSeqGene(name, start, end));
		}
		
		//System.err.println(new RefSeqGene(name, start, end));
		
		return paths;
	}
	
	public Collection<RefSeqGene> getBackwardPaths(int start, int end, double alpha, Collection<EdgeSourceType> validTypes) {
		List<RefSeqGene> paths = new ArrayList<RefSeqGene>();
		
		Iterator<Node<LightweightGenomicAnnotation>> verticesOverlappingIt = vertices.overlappers(start, end);
		List<LightweightGenomicAnnotation> overlappingVertices = new ArrayList<LightweightGenomicAnnotation>();
		while(verticesOverlappingIt.hasNext()) {
			overlappingVertices.add(verticesOverlappingIt.next().getValue());
		}
		Iterator<LightweightGenomicAnnotation> vertexIt = new ReversibleIterator<LightweightGenomicAnnotation>(overlappingVertices.listIterator(), false);
		
		List<LightweightGenomicAnnotation> verticesTraversedList = new ArrayList<LightweightGenomicAnnotation>();
		
		while(vertexIt.hasNext() ) {
			LightweightGenomicAnnotation v= vertexIt.next();
			//System.err.println("New Vertex " + v.toUCSC() + " traversed list " +  verticesTraversedList + " can be reached? " + canBeReachedFrom(v, verticesTraversedList));
			if(v.getStart() != start && !canBeReachedFrom(v, verticesTraversedList)) {
				// block one is from start to v 
				Alignments lastBlock = new  Alignments(name, v.getStart(), end); // Change to support vertices not overlappin edges MG
				// find other blocks by following each edge to the next splice
				Collection<BubbleEdge> edges = getBackwardEdges(v, alpha, validTypes);
				//System.err.println("Vertex " + v.toUCSC() );
				for(BubbleEdge edge : edges) {
					LightweightGenomicAnnotation prevV = getOpposite(v, edge);
					//System.err.println("\tEdge, " + edge.getConnection().toUCSC() + edge.getConnection().getOrientation() + " Oposite v " + prevV.toUCSC());
					Alignments orientedLastBlock  = new Alignments(lastBlock); // Change to support vertices not overlappin edges MG
					if(edge.inReversedOrientation()) {
						orientedLastBlock.setOrientation(edge.inReversedOrientation() ? "-" : "+");
					}
					Collection<RefSeqGene> nextPaths =  getBackwardPaths(prevV.getStart()  - (end - start) - lastBlock.length(), prevV.getEnd(), alpha, validTypes); //Recursion
					// append blocks found following each edge to the nascent transcript.
					if(nextPaths.size() == 0) {
						Alignments firstBlock = new Alignments(name, prevV.getStart()  - (end - start) - lastBlock.length(), prevV.getEnd());
						firstBlock.setOrientation(edge.inReversedOrientation() ? "-" : "+");
						List<Alignments> exons = new ArrayList<Alignments>();
						exons.add(firstBlock);
						exons.add(orientedLastBlock);
						RefSeqGene prependedPath = new RefSeqGene(exons);
						paths.add(prependedPath);
					} else {
						for(RefSeqGene newPath : nextPaths ) {
							Collection<Alignments> newPathIntrons = newPath.getIntronSet()	;// Change to support vertices not overlapping edges MG
							if(newPathIntrons.size() == 0) {
								verticesTraversedList.add(new BasicLightweightAnnotation(name, newPath.getStart(), newPath.getStart()+1));
							} else  {
							
								Iterator<Alignments> newPathIntronIt = newPathIntrons.iterator();
								while(newPathIntronIt.hasNext()) {
									Alignments intron = newPathIntronIt.next();
									verticesTraversedList.add(new BasicLightweightAnnotation(name, intron.getStart()-1, intron.getStart()));// Change to support vertices not overlapping edges MG
									verticesTraversedList.add(new BasicLightweightAnnotation(name, intron.getEnd(), intron.getEnd() + 1));// Change to support vertices not overlapping edges MG
								}
							}
							List<Alignments> exons = new ArrayList<Alignments>();
							exons.addAll(newPath.getExonSet()); //Changed to support vertices not overlapping edges
							exons.add(lastBlock); 
							RefSeqGene prependedPath = new RefSeqGene(exons);
							paths.add(prependedPath);
						}
					}
					//System.err.println("\tpaths: " + paths);
				}

			}
		}
		if(paths.isEmpty()){
			paths.add(new RefSeqGene(name, start, end));
		}
		
		//System.err.println(new RefSeqGene(name, start, end));
		
		return paths;
	}
	

	
	public Collection<RefSeqGene> getPaths(int start, int end, boolean forward) {
		return getPaths(start, end, 0, true);
	}
	
	private boolean canBeReachedFrom(LightweightGenomicAnnotation v, List<LightweightGenomicAnnotation> vertices) {
		boolean canBeReached = vertices.contains(v);
		
		Iterator<LightweightGenomicAnnotation> vIt = vertices.iterator();
		while(!canBeReached && vIt.hasNext()) {
			LightweightGenomicAnnotation traversedVertex = vIt.next();
			if(v.compareTo(traversedVertex) >= 0) {
				Number d = distanceAlgorithm.getDistance(v, traversedVertex);
				canBeReached = d != null;
			}
		}
		return canBeReached;
	}

	public Collection<Path> getPaths(double alpha) {
		return getPaths(alpha, 0);
	}
	public Collection<Path> getPaths(double alpha, double minSpliceFrequency) {
		Collection<Path> allPaths = new TreeSet<Path> ();
		//System.out.println("GETTING SOURCES");
		Collection<LightweightGenomicAnnotation> sources = getSourceVertices(getValidTypes()/*defaultTypesForDegreeCounting*/); // 06/01/2010 paired ends were ignored when looking for sources
		//System.err.println("Source vertices: "+sources.size()+" "+sources);
		for(LightweightGenomicAnnotation s : sources ) {
			//System.err.print("Getting paths starting in " + s.toUCSC() + " ... ");
			//System.out.println("using source: " + s);
			Collection<Path> sourcePaths = getPaths(s, alpha, minSpliceFrequency);
			allPaths.addAll(sourcePaths);
			//System.err.println("done");
		}
		
		return allPaths;
	}
	
		
	private boolean isEdgeValid(BubbleEdge e, double alpha, Collection<EdgeSourceType> validTypes) {
		//System.err.print("\tEdge "+e.getCount()+" "+minReads + " type <" + e.getType()+">");
		
		boolean valid = validTypes == null || validTypes.isEmpty();
		
		if(!valid) {
			Iterator<EdgeSourceType> typeIt = validTypes.iterator();
			while(!valid && typeIt.hasNext()) {
				EdgeSourceType type = typeIt.next();
				valid = type.ordinal() == e.getType().ordinal();
			}
		}
		
		return valid && e.getPvalue() <= alpha;

	}
	
	/*public Collection<RefSeqGene> getPaths() {
		return getPaths(0);
	}*/
	
	public Collection<RefSeqGene> getOrphanNodes(){
		Collection<RefSeqGene> rtrn=new TreeSet<RefSeqGene>();
		
		Iterator<LightweightGenomicAnnotation> iter=vertices.valueIterator();
		while(iter.hasNext()){
			LightweightGenomicAnnotation align=iter.next();
			if(this.getIncidentEdges(align).size()==0){rtrn.add(new RefSeqGene(align));}
		}
		
		return rtrn;
	}
	
	public Collection<Path> getOrphanPaths(){
		Collection<Path> rtrn=new TreeSet<Path>();
		
		Iterator<LightweightGenomicAnnotation> iter=vertices.valueIterator();
		while(iter.hasNext()){
			LightweightGenomicAnnotation align=iter.next();
			
			if(this.getIncidentEdges(align).size()==0){
				BubbleEdge edge=new BubbleEdge(align);
				edge.setParent(this);
				Path path=new Path(this.spliceWeight);
				//edge.setNodeCount(align, getCount(align));
				edge.setSpliceEdgeCount(0);
				//path.addEdge(edge);
				rtrn.add(path);
			}
		}
		
		return rtrn;
	}
	
	public BubbleEdge createEdge(Alignments connection, LightweightGenomicAnnotation v1, LightweightGenomicAnnotation v2, Integer count, EdgeSourceType type) {
		BubbleEdge edge=new BubbleEdge(connection, v1, v2, count, type);
		edge.setParent(this);
		return edge;
	}
	
	private double getCount(LightweightGenomicAnnotation align) {
		//TODO Consider putting back in
		/*if(this.exonScores==null){return 0.0;}
		if(!this.exonScores.containsKey(align)){return 0.0;}
		else{return this.exonScores.get(align);}*/
		
		if(this.vertexCounts==null){return 0.0;}
		else if(!this.vertexCounts.containsKey(align)){return 0.0;}
		else{return this.vertexCounts.get(align);}
	}

	private double getSpliceCount(LightweightGenomicAnnotation intron){
		if(this.intronScores==null || !this.intronScores.containsKey(intron)){return 0.0;}
		else{return this.intronScores.get(intron);}
	}
	
	public Collection<RefSeqGene> getOrphanNodes(Alignments region){
		Collection<RefSeqGene> rtrn=new TreeSet<RefSeqGene>();
		
		Iterator<Node<LightweightGenomicAnnotation>> iter=vertices.overlappers(region.getStart(), region.getEnd());
		while(iter.hasNext()){
			LightweightGenomicAnnotation align=iter.next().getValue();
			if(this.getIncidentEdges(align).size()==0){rtrn.add(new RefSeqGene(align));}
		}
		
		return rtrn;
	}

	public Collection<Path> getPaths(LightweightGenomicAnnotation startNode, double alpha, double minSpliceFreq){
		return getPaths(startNode, alpha, getValidTypes(), minSpliceFreq);
	}
	
	public Collection<Path> getPaths(LightweightGenomicAnnotation startNode, double alpha){
		return getPaths(startNode, alpha, getValidTypes(), 0);
	}
	
	public Collection<Path> getPaths(LightweightGenomicAnnotation startNode, double alpha, Collection<EdgeSourceType> validEdgeTypes) {
		return getPaths(startNode, alpha, validEdgeTypes, new HashMap<LightweightGenomicAnnotation, Collection<Path>>(),1, 0,0);
	}
	public Collection<Path> getPaths(LightweightGenomicAnnotation startNode, double alpha, Collection<EdgeSourceType> validEdgeTypes, double minSpliceFreq) {
		return getPaths(startNode, alpha, validEdgeTypes, new HashMap<LightweightGenomicAnnotation, Collection<Path>>(),1, 0, minSpliceFreq);
		//return getPaths(startNode, alpha, validEdgeTypes, new HashMap<LightweightGenomicAnnotation, Collection<Path>>(),1);

	}
	
	private Collection<Path> getPaths(LightweightGenomicAnnotation startNode, double alpha, Collection<EdgeSourceType> validEdgeTypes, Map<LightweightGenomicAnnotation, Collection<Path>> nodesVisited, int predictedPathsSoFar, int direction, double minFreq) {
		Collection<Path> allPaths = new ArrayList<Path>();

		List<BubbleEdge> allOutEdges = new ArrayList<BubbleEdge>(getOutOrUnorientedValidEdges(startNode));
		/*Collections.sort(allOutEdges, new Comparator<BubbleEdge>()  {
			public int compare(BubbleEdge o1, BubbleEdge o2) {
				return (int) (intronScores.get(o2.connection) - intronScores.get(o1.connection)) ;
			}
		}); */
		Collection<BubbleEdge> filteredOutEdges = filterHiFrequencyEdges(minFreq, allOutEdges);
		//if(minFreq > 0) {System.out.println("Paths from " + startNode.toUCSC() + " following " + filteredOutEdges.size() + " instead of  " + allOutEdges.size());}
		/* Sort edges by their counts */
		double t = System.nanoTime();
		//System.err.println("\tOut edges: "+filteredOutEdges.size()+" "+startNode.toUCSC()+ " took "  + (System.nanoTime() - t));
		//System.out.println("\tAll out edges: " + allOutEdges.size() + " minFreq: " + minFreq);
		int newPredictedPaths = predictedPathsSoFar * filteredOutEdges.size();
		boolean isToComplex = false;
		if(newPredictedPaths < MAX_GENE_PATHS) {
			for(BubbleEdge e : filteredOutEdges) {
				if(e.getPvalue() <= alpha && isEdgeValid(e, alpha, validEdgeTypes))  {
					//System.err.println("\t\tEdge " + ": " + e.getConnection().toUCSC()+" type <" + e.getType()+">");
					//ArrayList<LightweightGenomicAnnotation> visitedSoFar = new ArrayList<LightweightGenomicAnnotation>(nodesVisited);
					//System.err.println(" passed, is allowed");
					LightweightGenomicAnnotation next = getOpposite(startNode, e);
					List<BubbleEdge> allInEdges = new ArrayList<BubbleEdge>(getInOrUnorientedValidEdges(next) );
					double totalCounts = 0;
					for(BubbleEdge inE: allInEdges) {
						if(inE.getType().equals(EdgeSourceType.SPLICED)) { //We find edge frequency for splice edges not others at this point.
							totalCounts += intronScores.get(inE.connection);
						}
					}
					if(e.getType().equals(EdgeSourceType.SPLICED)) {
						double thisInFreq = intronScores.get(e.connection)/totalCounts;
						if(thisInFreq < minFreq) {
							//System.err.println("Edge " + e.getConnection().toUCSC() + " had decent out frequency, but fails in frequency: " + thisInFreq + " total in: " + totalCounts);
							continue;
						}
					}
					int comparison = startNode.compareTo(next) < 0 ? -1 : 1;
					
					if(comparison * direction >= 0) { //Only proceed forward if next node is compatible with the direction in which we are traversing.
						
						Collection<Path> toAppend = null;
						if(nodesVisited.containsKey(next) ) {
							toAppend = nodesVisited.get(next);
						} else {
							//System.err.println("\t\tat " + startNode.toUCSC() +" 2next "+next.toUCSC()  + " going " + comparison+" nodes visited: " + nodesVisited);
						
							next.setOrientation(e.getOrientation());
							
							toAppend = getPaths(next, alpha,validEdgeTypes, nodesVisited, newPredictedPaths, comparison, minFreq);
							nodesVisited.put(next, toAppend);
							//System.err.println("Append: "+toAppend.size());
						}
						for(Path g : toAppend) {
							Path gCopy = g.copy();
							//System.err.println("gene "+g.toGene().toBED());
							//Check whether the paths to append are unnoriented, or the edge is unnoriented or the orientation of the edge and the path to append are the same.
							if("*".equalsIgnoreCase(g.getOrientation()) || "*".equalsIgnoreCase(e.getOrientation()) || g.getOrientation().equalsIgnoreCase(e.getOrientation())){
								//gCopy.addEdge(e);
								//System.err.println("\t\tAppending to all paths which have current size: " + allPaths.size());
								allPaths.add(gCopy);
							}
						}
					} else {
						//System.err.println("Note, inconsistent edge for current path found at "+startNode.toUCSC() +", next node would have been: " + next.toUCSC()+ " not traversing");
						
					}
					/*for(Path g : toAppend) {
						Path gCopy = g.copy();
						//System.err.println("gene "+g.toGene().toBED());
						//Check whether the paths to append are unnoriented, or the edge is unnoriented or the orientation of the edge and the path to append are the same.
						if("*".equalsIgnoreCase(g.getOrientation()) || "*".equalsIgnoreCase(e.getOrientation()) || g.getOrientation().equalsIgnoreCase(e.getOrientation())){
							gCopy.addEdge(e);
							//System.err.println("\t\tAppending to all paths which have current size: " + allPaths.size());
							allPaths.add(gCopy);
						}
					}*/
	
				} else {
					//System.err.println(e.getConnection().toUCSC() + " type " + e.type+ " did not pass, is contained in valid edges? " + validEdgeTypes.contains(e) + " valid edges " + validEdgeTypes);
					//throw new RuntimeException(e.getConnection().toUCSC() + " type " + e.type+ " did not pass, is contained in valid edges? " + validEdgeTypes.contains(e) + " valid edges " + validEdgeTypes);
				}
			}
		} else {
			//System.err.println("WARNING: MAX GENE PATHS "+MAX_GENE_PATHS+" REACHED AT NODE " + startNode.toUCSC() +" will stop recursion, transcript will be truncated.");
			isToComplex = true;
		}
		if(allPaths.isEmpty() && !isSpanningGap(startNode.getStart(), startNode.getEnd())){
			BubbleEdge edge=new BubbleEdge(startNode);
			edge.setParent(this);
			//edge.setNodeCount(startNode, getCount(startNode));
			edge.setSpliceEdgeCount(0);
			Path path=new Path(this.spliceWeight);
			path.setIsToComplex(isToComplex);
			//path.addEdge(edge);
			allPaths.add(path);
		}
		
		//System.err.println("\tAll paths: "+allPaths.size());
		return allPaths;
	}

	private Collection<BubbleEdge> filterHiFrequencyEdges(double minFreq,List<BubbleEdge> allEdges) {
		int total = 0;
		for( BubbleEdge edge : allEdges ) {
			if(edge.getType().equals(EdgeSourceType.SPLICED)) {
				//System.out.println("Adding to total: " + intronScores.get(edge.connection) + " ::: " + total);
				total += intronScores.get(edge.connection); //ONLY FILTER COUNT SPLICE EDGES
			}
		}
		Collection<BubbleEdge> filteredEdges = new TreeSet<BubbleEdge>();
		double freq = 1;
		Iterator<BubbleEdge> outEdgeIt = allEdges.iterator();
		//while(outEdgeIt.hasNext() && freq > minFreq) { // IS THIS CONDITION (freq > minFreq) RIGHT??  -- MAYBE ONLY IF allEdges IS SORTED BY intronScores.get(edge.connection)
		while (outEdgeIt.hasNext()) {
			BubbleEdge edge = outEdgeIt.next();
			if(!edge.getType().equals(EdgeSourceType.SPLICED)) {
				filteredEdges.add(edge); // ONLY FILTER SPLICE EDGES THIS SHOULD BE A TEMPORARY HACK
			}else if(intronScores.get(edge.connection) != null) {
				freq = intronScores.get(edge.connection)/total;
			//if(minFreq > 0) {System.out.println("\tedge " + edge.connection.toUCSC() + edge.connection.getOrientation() + " counts " + intronScores.get(edge.connection) +" freq: " + freq + " min freq " + minFreq);}
				if(freq > minFreq) {
					filteredEdges.add(edge);
				}
			} else {
				System.err.println("ERROR: there is are no counts for edge "+ edge.connection.toUCSC());
			}
		}
		return filteredEdges;
	}
	
	private boolean compatible(LightweightGenomicAnnotation node, List<LightweightGenomicAnnotation> nodesVisited) {
		boolean compatible = !nodesVisited.contains(node) && !node.overlaps(nodesVisited);
		if(compatible && nodesVisited.size() > 2) {
			int direction = nodesVisited.get(0).compareTo(nodesVisited.get(1));
			int thisNodeDirection = nodesVisited.get(nodesVisited.size() - 1).compareTo(node);
			
			compatible = direction * thisNodeDirection > 0 ;
		}
		return compatible;
	}
	


	
	public int getChromosomeLength() {
		return chromosomeLength;
	}

	public void setChromosomeLength(int chromosomeLength) {
		this.chromosomeLength = chromosomeLength;
	}


	List<LightweightGenomicAnnotation> getOverlappingVertices(int start, int end) {
		ArrayList<LightweightGenomicAnnotation> overlappers = new ArrayList<LightweightGenomicAnnotation>();
		Iterator<Node<LightweightGenomicAnnotation>> overlapperIt = vertices.overlappers(start, end);
		
		while(overlapperIt.hasNext()) {
			overlappers.add(overlapperIt.next().getValue());
		}
		return overlappers;
	}
	
	public WindowIterator iterator(int windowSize, int start, int overlap, boolean skipGapSpanningWindows) {
		return new WindowIterator(this, windowSize, start, overlap, skipGapSpanningWindows);
	}
	
	public WindowIterator iterator(int windowSize, int start) {
		return new WindowIterator(this, windowSize, start, windowSize - 1, true);
	}
	
	public WindowIterator reverseIterator(int windowSize, int start) {
		return new ReverseWindowIterator(this, windowSize, start, windowSize - 1, true);
	}
	
	public WindowIterator fowardAndBackIterator(int windowSize) {
		return new ForwardAndBackIterator(this, windowSize, windowSize - 1, true);
	}

	public static class BubbleEdge implements Comparable<BubbleEdge>{
		
		private double count;
		private double pairedEndCounts;
		private double splicedCount;
		//private Map<LightweightGenomicAnnotation, Double> nodeCounts;
		private double pvalue;
		private LightweightGenomicAnnotation connection;
		private LightweightGenomicAnnotation v1;
		private LightweightGenomicAnnotation v2;
		private EdgeSourceType type;
		private Collection<Alignment> pairedEndSupport ;
		private ChromosomeWithBubbles2 parentGraph;
				
		/**
		 * Need to store source and destination vertices to be able to compare edges.
		 * @param connection
		 * @param v1
		 * @param v2
		 * @param connectionCount
		 */
		
		BubbleEdge(LightweightGenomicAnnotation connection, LightweightGenomicAnnotation v1, LightweightGenomicAnnotation v2, double connectionCount) {
			this(connection, v1, v2, connectionCount, EdgeSourceType.SPLICED);
		}


		/**
		 * Need to store source and destination vertices to be able to compare edges.
		 * @param connection
		 * @param v1
		 * @param v2
		 * @param connectionCount
		 * @param type of connection.
		 */
		public BubbleEdge(LightweightGenomicAnnotation connection, LightweightGenomicAnnotation v1, LightweightGenomicAnnotation v2, double connectionCount, EdgeSourceType type) {
			this.connection = connection;
			this.count = connectionCount;
			this.v1=(v1);
			this.v2=(v2);
			this.type = type;
			pairedEndSupport = new ArrayList<Alignment>(); //CANT BE A SET since paired ends may be duplicates
			if(this.connection.getOrientation()==null){this.connection.setOrientation("*");}
		}
		
		
		
		BubbleEdge(BubbleEdge toCopy) {
			this.parentGraph = toCopy.parentGraph;
			this.connection = toCopy.connection;
			this.count = toCopy.count;
			this.v1=(toCopy.v1) ;
			this.v2=(toCopy.v2);
			this.type = toCopy.type;
			pairedEndSupport = toCopy.pairedEndSupport;
		}
		
		public BubbleEdge(LightweightGenomicAnnotation startNode) {
			this(startNode, startNode, startNode, 0, EdgeSourceType.SELF);
		}
		
		public void setParent(ChromosomeWithBubbles2 parent) {
			this.parentGraph = parent;
		}
		
		public void setSpliceEdgeCount(double i) {
			this.splicedCount=i;
		}

		public void setNodeCount(LightweightGenomicAnnotation node, double i) {
			//System.err.println("node: " +node);
			//System.err.println("node " + node + " vertex counts " + (parentGraph == null ? null : parentGraph.vertexCounts) + " parent " + parentGraph);
			parentGraph.setExonScores(node, new Double(i));
			//graph.setNodeCounts(node, i);
		}

		public double getNodeCount(LightweightGenomicAnnotation node){
			if(EdgeSourceType.SELF.equals(getType())) {
				return parentGraph.getCount(connection);
			} else {
				return parentGraph.vertexCounts.get(node);		
			}
		}
			
		public LightweightGenomicAnnotation getOpposite(LightweightGenomicAnnotation v) {
			LightweightGenomicAnnotation opposite = null;
			if(v.equals(v1)) {
				opposite = v2;
			} else if (v.equals(v2)) {
				opposite = v1;
			}
			return opposite;
		}

		public EdgeSourceType getType() {
			return type;
		}

		

		public void incrementCount() {
			count++;
		}
		
		public void incrementCount(float amount) {
			count += amount;
		}

		public String getOrientation(){return connection.getOrientation();}

		
		//TODO: Fix the counts
		public double getAllCounts() {
			return getPairedEndCounts()+getSplicedCounts();
		}
		
		public double getPairedEndCounts(){return this.pairedEndCounts;}
		
		public void setPairedEndCounts(double count){this.pairedEndCounts=count;}
		
		public double getSplicedCounts(){return this.splicedCount;}
		
		public void setCount(double count) {
			this.count = count;
		}
		public LightweightGenomicAnnotation getConnection() {
			return connection;
		}
		
		public void addPairedEndSupport(Alignment a){
			pairedEndSupport.add(a);
		}
		
		public Collection<Alignment> getPairedEndSupport() { return pairedEndSupport;}
		
		public boolean inReversedOrientation() {return connection.inReversedOrientation();}

		public double getPvalue() {
			return pvalue;
		}
		public void setPvalue(double pvalue) {
			this.pvalue = pvalue;
		}

		public int compareTo(BubbleEdge o) {
			int comparison = o.connection.compareTo(connection);
			if(comparison == 0) {
				if(connection.inReversedOrientation() != o.connection.inReversedOrientation()) {
					comparison = connection.inReversedOrientation() ? -1 : 1;
				}
				if(comparison == 0) {
					comparison = getLeftNode().compareTo(o.getLeftNode());
				}
				
				if(comparison == 0) {
					comparison = getRightNode().compareTo(o.getRightNode());
				}
			}
			return comparison;
		}
		
		public int hashCode() {
			return (connection.hashCode() + getLeftNode().hashCode() + getRightNode().hashCode())* (connection.inReversedOrientation() ? -1 : 1);
		}
		
		public boolean equals(Object o) {
			boolean ret = false;
			if(o instanceof BubbleEdge) {
				ret = compareTo((BubbleEdge) o) == 0;
			}
			return ret;
		}

		

		public LightweightGenomicAnnotation getLeftNode() {
			return v1;
		}

		public LightweightGenomicAnnotation getRightNode() {
			return v2;
		}

		/*public double getNodeCounts() {
			return this.v1.getScore()+this.v2.getScore();
		}*/
		

	}
	public static class MinReadFilter implements ReadFilter {
		int minLength;
		public MinReadFilter(int minLength) {
			this.minLength = minLength;
		}

		public boolean passes(Alignments read) {
			return read.length() > minLength;
		}
		
	}
	
	public static class RegionReadFilter implements ReadFilter {
		int start;
		int end;
		
		public RegionReadFilter(int start, int end) {
			this.start = start;
			this.end  = end;
		}
		
		public boolean passes(Alignments read) {
			return read.getStart() < end && read.getEnd() > start; //If read overlaps region.
		}
	}
	public static class WindowIterator implements Iterator<Collection<RefSeqGene>> {
		private ChromosomeWithBubbles2 data;
		private int atPosition;
		private int windowSize;
		private int step;
		private int dataAccess;
		private int numWindowsDone;
		private int winWithNoJumps;
		
		private Collection<RefSeqGene> last;
		private List<LightweightGenomicAnnotation> lastOverlappingVertices;
		
		protected WindowIterator() {
			
		}
		
		//DONT REALLY NEED THE SKIPGAPSPANNINGWINDOWS. NOT USED
		public WindowIterator(ChromosomeWithBubbles2 data, int windowSize, int start, int overlap, boolean skipGapSpanningWindows) {
			this.data = data;
			this.windowSize = windowSize;
			this.atPosition = start;
			this.step = windowSize - overlap;
		}
		public boolean hasNext() {
			return atPosition < data.getEnd();
		}
		
		public int getDataAccess() { return dataAccess;}
		public int getNumWindowsDone() { return numWindowsDone;}
		public int getWinWithNoJumps() { return winWithNoJumps;}

		public Collection<RefSeqGene> next() {
			int end   = atPosition + windowSize;
			//long start = System.nanoTime();
			//long t = System.currentTimeMillis();
			/*
			int farthestLast = 0;
			if(last == null) {
				farthestLast = end;
			} else {
				for(RefSeqGene g : last) {
					farthestLast = Math.max(farthestLast, g.getEnd());
				}
			}
			List<LightweightGenomicAnnotation> overlappingVertices = data.getOverlappingVertices(atPosition, farthestLast+1);
			//System.err.println("Gog vertices " + (System.currentTimeMillis() - t));
			Collection<RefSeqGene> result = null;
			if(overlappingVertices.isEmpty()) {	
				//t = System.currentTimeMillis();
				last = null;
				lastOverlappingVertices = null;
				result = new ArrayList<RefSeqGene>(1);
				result.add(new RefSeqGene(data.getName(), atPosition, end));
				//System.err.print("Trivial window " + (System.currentTimeMillis() - t) + " ");
				winWithNoJumps++;
			} else {
				if(lastOverlappingVertices == null || lastOverlappingVertices.isEmpty()) {
					//t = System.currentTimeMillis();
					result = data.getPaths(atPosition, end);
					//System.err.print("Recursion (no vertices) " + (System.currentTimeMillis() - t)+ " ");
					dataAccess++;
				} else {
					//t = System.currentTimeMillis();
					LightweightGenomicAnnotation firstOverlapper = overlappingVertices.get(0);
					LightweightGenomicAnnotation lastOverlapper = overlappingVertices.get(overlappingVertices.size() - 1);
					
					LightweightGenomicAnnotation firstLastOverlapper = lastOverlappingVertices.get(0);
					LightweightGenomicAnnotation lastLastOverlapper = lastOverlappingVertices.get(lastOverlappingVertices.size() - 1);
					
					if(firstOverlapper.equals(firstLastOverlapper) &&
							lastOverlapper.equals(lastLastOverlapper) &&
							firstOverlapper.getStart() > atPosition + 1 && 
							lastOverlapper.getEnd() < end - step - 1 ) {//TODO: This is not working should use farthestLast but it introduces a bug. FIX!

						result = new ArrayList<RefSeqGene>(last.size());
						for(RefSeqGene g : last) {
							List<Alignments> lastExons = new ArrayList<Alignments>(g.getExonSet());
							Alignments lastFirstExon = lastExons.get(0);
							lastFirstExon.setStart(lastFirstExon.getStart() + step);
							
							Alignments lastLastExon = lastExons.get(lastExons.size() - 1);
							lastLastExon.setEnd(lastLastExon.getEnd() + step);
							RefSeqGene newG = new RefSeqGene(lastExons);
							result.add(newG);
						}
						
						//System.out.print("Used Cache " + (System.currentTimeMillis() - t)+ " ");
					} else {
						//t = System.currentTimeMillis();
						result = data.getPaths(atPosition, end);
						//System.err.println("Recursion (could not use cache) " + (System.currentTimeMillis()- t));
						dataAccess++;
					}
					
				}				
				
			}
			//Move to next position that does not span a gap.
			atPosition = atPosition + step;
			while(skipGapSpanningWindows && data.isSpanningGap(atPosition, windowSize) ) {
				atPosition++;
			}
			
			
			last = result;
			//System.err.println("SIZES: "+result.size()+" "+filteredResults.size());
			lastOverlappingVertices = overlappingVertices;
			numWindowsDone++;
			*/
			//System.err.println("Total time: " + (System.nanoTime()  - start));
			last = data.getPaths(atPosition, end, 0, true);
			for(RefSeqGene window : last) {
				Collection<Alignments> windowIntrons = window.getIntronSet();
					for(Alignments intron: windowIntrons) {
						if(data.findEdge(intron.getStart(), intron.getEnd())==null) {throw new RuntimeException("BAD path intron: \n" + intron + "\npath:\n" + window.toBED());}
					}
			}
			atPosition += step;
			return last;
		}


		public void remove() {
			//
		}
		
		public Collection<RefSeqGene> next(int lowerBound, int higherBound) {
			atPosition = higherBound;
			int end = higherBound + windowSize;
			last = data.getPaths(atPosition, end, 0 , true);
			lastOverlappingVertices = data.getOverlappingVertices(atPosition, end);
			return last;
		}
		
		public void jumpTo(int lowerBound, int higherBound) {
			atPosition = higherBound;
		}
		
	}
	
	public static class ReverseWindowIterator extends WindowIterator {	
		public ReverseWindowIterator(ChromosomeWithBubbles2 data, int windowSize, int start, int overlap, boolean skipGapSpanningWindows) {
			super(data, windowSize, start, overlap, skipGapSpanningWindows);
			//Move to next position that does not span a gap.
		}
		public boolean hasNext() {
			return super.atPosition > super.data.getStart();
		}
		
		public Collection<RefSeqGene> next() {
			int start   = super.atPosition - super.windowSize;
			super.last = super.data.getPaths(start, super.atPosition, 0, false);
			/*System.err.println("At position " + super.atPosition + " paths found: " + super.last.size());
			for(RefSeqGene window : super.last) {
				Collection<Alignments> windowIntrons = window.getIntronSet();
					for(Alignments intron: windowIntrons) {
						if(super.data.findEdge(intron.getStart(), intron.getEnd())==null) {throw new RuntimeException("BAD path intron: \n" + intron + "\npath:\n" + window.toBED());}
					}
			}*/
			super.atPosition -= super.step;
			return super.last;
		}


		public void remove() {
			//
		}
		
		public Collection<RefSeqGene> next(int lowerBound, int higherBound) {
			super.atPosition = lowerBound;
			int start = lowerBound - super.windowSize;
			super.last = super.data.getPaths(start, super.atPosition, 0, false);
			super.lastOverlappingVertices = super.data.getOverlappingVertices(start, super.atPosition);
			return super.last;
		}
		
		public void jumpTo(int lowerBound, int higherBound) {
			super.atPosition = lowerBound;
		}
		
	}

	public enum EdgeSourceType {
	    SPLICED, PAIRED, PAIRED_SUPPORT, SELF; 

	}
	
	public static class ForwardAndBackIterator extends WindowIterator {
		LinkedList<WindowIterator> iteratorQueue;
		
		public ForwardAndBackIterator(ChromosomeWithBubbles2 data, int windowSize,  int overlap, boolean skipGapSpanningWindows) {
			iteratorQueue = new LinkedList<WindowIterator>();
			WindowIterator forwardIt = new WindowIterator(data, windowSize, data.getStart(), overlap, skipGapSpanningWindows);
			iteratorQueue.add(forwardIt);
			WindowIterator backIt = new ReverseWindowIterator(data, windowSize,data.getEnd(), overlap, skipGapSpanningWindows);
			iteratorQueue.add(backIt);
		}

		public boolean hasNext() {
			WindowIterator current = iteratorQueue.peek();
			boolean hasNext = false;
			if(current != null) {
				hasNext = current.hasNext();
				if(!hasNext) {
					iteratorQueue.pop();
					System.err.println("Switching iterator");
					hasNext = hasNext();
				}
			}
			return hasNext;
		}

		public Collection<RefSeqGene> next() {
			return hasNext() ? iteratorQueue.peek().next() : null;
		}

		public void remove() {
			// TODO Auto-generated method stub
			
		}
			
		public Collection<RefSeqGene> next(int lowerBound, int higherBound) {
			return hasNext() ? iteratorQueue.peek().next(lowerBound, higherBound) : null;
		}
	}

	public LightweightGenomicAnnotation findEdge(int start, int end) {
		Node<BubbleEdge> edgeNode = edgeTree.find(start, end);
		return edgeNode != null ? edgeNode.getValue().getConnection() : null;
	}

	public void setExonScores(LightweightGenomicAnnotation node, Double score) {
		this.vertexCounts.put(node, score);
		
	}

	public Collection<RefSeqGene> getGenePaths(double alpha) {
		return getGenePaths(alpha, 0);
	}
	public Collection<RefSeqGene> getGenePaths(double alpha, double minSpliceFrequency) {
		Collection<RefSeqGene> genes=new TreeSet<RefSeqGene>();
		Collection<Path> paths=getPaths(alpha, minSpliceFrequency);
		for(Path path: paths){
			genes.add(path.toGene());
		}
		return genes;
	}
	
	public void updatePairedEndEdges(AlignmentDataModelStats pairedData) throws IOException{
		Collection<BubbleEdge> pairedEdges = getPairedEndEdges();
		System.err.println("Got " + pairedEdges.size() + " paired edges");
		for(BubbleEdge pe : pairedEdges) {
			IntervalTree<Alignment> readTree = pairedData.getIntervalTreeCached(getName(), pe.getConnection().getStart(), pe.getConnection().getEnd());
			Iterator<Node<Alignment>> it =readTree.overlappers(pe.getConnection().getStart(), pe.getConnection().getEnd());
			int start =  pe.getConnection().getStart();
			int end   = pe.getConnection().getEnd();
			int num = 0;
			while(it.hasNext()) {
				Alignment a = it.next().getValue();
				num++;
				//pe.addPairedEndSupport(a);
				start = Math.min(start, a.getAlignmentStart());
				end   = Math.max(end, a.getAlignmentEnd());
			}
			double [] scores = pairedData.scanPRate(new BasicLightweightAnnotation(getName(), start, end), num, 0);
			pe.setCount(num);
			pe.setPvalue(scores[0]);
		}
		
	}
	
	private Collection<BubbleEdge> getPairedEndEdges() {
		Collection<BubbleEdge> pairedEdges = new ArrayList<BubbleEdge>();
		Iterator<BubbleEdge> connectionIt = edgeTree.valueIterator();
		while(connectionIt.hasNext()) {
			BubbleEdge e = connectionIt.next();
			//System.err.println("Edge " + e.getConnection().toUCSC() + " type: " + e.getType());
			if(e.getType().equals(EdgeSourceType.PAIRED)) {
				//System.err.println("Added");
				pairedEdges.add(e);
			}	
		}
		return pairedEdges;
	}
	
	public boolean removeEdge(BubbleEdge e) {
		BubbleEdge removed = edgeTree.remove(e.getConnection().getStart(), e.getConnection().getEnd());
		boolean ret = removed != null;
		if(ret || containsEdge(e)) { //The parent graph may contain duplicate edges according to us. So even if we do not find it, the parent may still have it.
			ret = super.removeEdge(e); 
		}
		
		return ret;
	}

	//Write the graph as a DOT format
	public void writeGraph(String save, boolean setUpShortExonLabels)throws IOException{
		Alignments chrRegion=new Alignments(this.getName(), 0, this.end);
		writeGraph(save, chrRegion, setUpShortExonLabels);
	}
	
	public void writeGraph(String save)throws IOException{
		writeGraph(save, true);
	}
	
	public void writeGraph(String save, Alignments region) throws IOException {
		writeGraph(save, region, true);
	}
	
	public void writeGraph(String save, Alignments region, boolean setUpShortExonLabels, String graphName) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("digraph "+ graphName+" {\n rankdir=LR; \n node [shape = rectangle];\n");
		//go through each and node and write it and its edges
		//get all the graphs edges
	//	Map<BubbleEdge, Pair<LightweightGenomicAnnotation>> edgeMap= edges;
		
		//start by writing an invisible connection between every node and every node in sorted order
	//	IntervalTree<LightweightGenomicAnnotation> nodeTree=this.vertices;
		Iterator<Node<LightweightGenomicAnnotation>> iter=vertices.iterator();
		
		int counter=0;
		while(iter.hasNext()){ //TODO: 1) Consider moving to the bottom 2) Consider not adding invisible edge to edges that already exist
			LightweightGenomicAnnotation node=iter.next().getValue();
			if(node.overlaps(region)){
				if(counter>0){writer.write("->");}
				writer.write(quote+node.toUCSC()+quote);
				counter++;
			}
		}
		
		writer.write("[style=invis];\n"); //TODO: Scale the invisible link by the genomic size
		
		iter=vertices.iterator();
		//Then define the sizes of each node based on genomic sizes
		int exonCounter=0;
		while(iter.hasNext()){
			LightweightGenomicAnnotation exon=iter.next().getValue();
			double exonCount=this.getCount(exon);
			double localRate=this.getLocalRate(exon);
			if(exon.overlaps(region)){
				double normWidth=exon.length()*scaleFactor;
				writer.write(quote+exon.toUCSC()+quote+" [width="+normWidth+", fixedsize=true, label="+exonCount+", comment="+localRate); 
				writer.write("];\n");
				exonCounter++;
			}
		}
		
		//This defines actual edges
		Collection<BubbleEdge> edges = getEdges();
		for(BubbleEdge edge: edges){
			Pair<LightweightGenomicAnnotation> nodes=getEndpoints(edge);
			//System.err.println("WRITEGRAPH - edge" + edge.getConnection().toUCSC());
			LightweightGenomicAnnotation first=nodes.getFirst();
			LightweightGenomicAnnotation second=nodes.getSecond();
			if(first.overlaps(region) || second.overlaps(region)){
				writer.write(quote+first.toUCSC()+quote+"->"+quote+second.toUCSC()+quote);
				writer.write("[label="+edge.getAllCounts()+", comment="+quote+edge.getType()+quote);
				if(edge.getType().equals(EdgeSourceType.PAIRED)){
					writer.write(", style=dotted");
				}
				//writer.write(", minlen="+(edge.getConnection().length()*(scaleFactor/3))); //TODO: Add scaling for intron sizes (default no scaling)
				writer.write("];\n");
			}
		}
					
		writer.write("}");
		
		writer.close();
	}
	
	
	//TODO Write this format
	public void writeGraphFormat(String save, Alignments region) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		//Get all connected components
		
		//For each connected component write a graph containing all exons and connections
		
		//writer.write("digraph "+ graphName+" {\n rankdir=LR; \n node [shape = rectangle];\n");
		//go through each and node and write it and its edges
		//get all the graphs edges
	//	Map<BubbleEdge, Pair<LightweightGenomicAnnotation>> edgeMap= edges;
		
		//start by writing an invisible connection between every node and every node in sorted order
	//	IntervalTree<LightweightGenomicAnnotation> nodeTree=this.vertices;
		Iterator<Node<LightweightGenomicAnnotation>> iter=vertices.iterator();
		
		int counter=0;
		while(iter.hasNext()){ //TODO: 1) Consider moving to the bottom 2) Consider not adding invisible edge to edges that already exist
			LightweightGenomicAnnotation node=iter.next().getValue();
			if(node.overlaps(region)){
				if(counter>0){writer.write("->");}
				writer.write(quote+node.toUCSC()+quote);
				counter++;
			}
		}
		
		writer.write("[style=invis];\n"); //TODO: Scale the invisible link by the genomic size
		
		iter=vertices.iterator();
		//Then define the sizes of each node based on genomic sizes
		int exonCounter=0;
		while(iter.hasNext()){
			LightweightGenomicAnnotation exon=iter.next().getValue();
			double exonCount=this.getCount(exon);
			double localRate=this.getLocalRate(exon);
			if(exon.overlaps(region)){
				double normWidth=exon.length()*scaleFactor;
				writer.write(quote+exon.toUCSC()+quote+" [width="+normWidth+", fixedsize=true, label="+exonCount+", comment="+localRate); 
				writer.write("];\n");
				exonCounter++;
			}
		}
		
		//This defines actual edges
		Collection<BubbleEdge> edges = getEdges();
		for(BubbleEdge edge: edges){
			Pair<LightweightGenomicAnnotation> nodes=getEndpoints(edge);
			//System.err.println("WRITEGRAPH - edge" + edge.getConnection().toUCSC());
			LightweightGenomicAnnotation first=nodes.getFirst();
			LightweightGenomicAnnotation second=nodes.getSecond();
			if(first.overlaps(region) || second.overlaps(region)){
				writer.write(quote+first.toUCSC()+quote+"->"+quote+second.toUCSC()+quote);
				writer.write("[label="+edge.getAllCounts()+", comment="+quote+edge.getType()+quote);
				if(edge.getType().equals(EdgeSourceType.PAIRED)){
					writer.write(", style=dotted");
				}
				//writer.write(", minlen="+(edge.getConnection().length()*(scaleFactor/3))); //TODO: Add scaling for intron sizes (default no scaling)
				writer.write("];\n");
			}
		}
					
		writer.write("}");
		
		writer.close();
	}

	public void writeGraph(String save, Alignments region, boolean setUpShortExonLabels) throws IOException {

		writeGraph(save, region ,setUpShortExonLabels, getName()+"_"+lambda+"_"+this.numberOfMarkers+"_"+this.numberOfReads);
	}
	
	private void computeLocalRates() {
		//Collection<LightweightGenomicAnnotation> sources = getSourceVertices(defaultTypesForDegreeCounting);
		//Collection<LightweightGenomicAnnotation> son
	}

	public void setLocalRate(Collection<Path> paths){
		this.localRate=new TreeMap();
		for(Path path: paths){
			Collection<Alignments> exons=path.toGene().getExonSet();
			for(Alignments exon: exons){
				this.localRate.put(exon, path.getLocalLambda());
			}
		}
	}
	
	public double getLocalRate(LightweightGenomicAnnotation exon){
		if(this.localRate==null || !this.localRate.containsKey(exon)){return this.lambda;}
		return Math.max(this.localRate.get(exon), lambda);
	}
	
	public Map<Alignments, Double> getLocalRate(){return this.localRate;}
	
	public Pair<LightweightGenomicAnnotation> getNodePair(BubbleEdge edge) {
		return super.edges.get(edge);
	}

	public Collection<Path> getAllPaths() {
		return getAllPaths(0);
	}
	public Collection<Path> getAllPaths(double minSpliceFrequency) {
		Collection<Path> rtrn=new TreeSet<Path>();
		
		rtrn.addAll(getOrphanPaths());
		rtrn.addAll(getPaths(1, minSpliceFrequency)); 
		
		return rtrn;
	}

	public double getNumberOfMarkers() {
		return this.numberOfMarkers;
	}

	public double getNumberOfReads() {
		return this.numberOfReads;
	}

	public double getRPKMConstant() {
		return this.rpkmConstant;
	}
	
	public void setRPKMConstant(double num) {
		this.rpkmConstant=num;
	}

	public void setLocalRate(Map<Alignments, Double> localRateMap) {
		this.localRate=localRateMap;		
	}

	public void setEstimatedPairedDistribution(EmpiricalDistribution pairedDist) {
		this.pairedDist=pairedDist;
	}
	
	public EmpiricalDistribution getEstimatedPairedDistribution(){
		return this.pairedDist;
	}

	public void setSpliceWeight(double spliceWeight) {this.spliceWeight=spliceWeight;}
	
	public double getSpliceWeight(){return this.spliceWeight;}
/*
	public mxGraph buildVisualizationGraph(Alignments region) {
		Map<String, mxCell> map=new TreeMap();
		Collection<LightweightGenomicAnnotation> nodes=new TreeSet();
		
		mxGraph graph = new mxGraph();
		
		mxStylesheet stylesheet = graph.getStylesheet();
		Hashtable<String, Object> style = new Hashtable<String, Object>();
		style.put(mxConstants.STYLE_OPACITY, 0);
		stylesheet.putCellStyle("invis", style);
		
		//mxCompactTreeLayout layout=new mxCompactTreeLayout(graph ); 
		mxHierarchicalLayout layout=new mxHierarchicalLayout(graph, SwingConstants.WEST); 
		Object parent = graph.getDefaultParent();
		graph.getModel().beginUpdate();
		try
		{
					
		Iterator<Node<LightweightGenomicAnnotation>> iter=vertices.iterator();
		//Then define the sizes of each node based on genomic sizes
		int exonCounter=0;
		while(iter.hasNext()){
			LightweightGenomicAnnotation exon=iter.next().getValue();
			if(exon.overlaps(region)){
				double normWidth=exon.length()/10;
				mxCell v1 = (mxCell)graph.insertVertex(parent, null, "e"+exonCounter, exon.getStart(), 20, normWidth,30);
				map.put(exon.toUCSC(), v1);
				nodes.add(exon);
				exonCounter++;
			}
		}
		
		//This defines actual edges
		Collection<BubbleEdge> edges = getEdges();
		for(BubbleEdge edge: edges){
			Pair<LightweightGenomicAnnotation> nodePair=getEndpoints(edge);
			//if(edge.getType().equals(EdgeSourceType.SPLICED)){
				//System.err.println("WRITEGRAPH - edge" + edge.getConnection().toUCSC());
				LightweightGenomicAnnotation first=nodePair.getFirst();
				LightweightGenomicAnnotation second=nodePair.getSecond();
				if(first.overlaps(region) && second.overlaps(region)){
					mxCell e1=(mxCell) graph.insertEdge(parent, null, map.get(first.toUCSC()).getValue()+"-"+map.get(second.toUCSC()).getValue(), map.get(first.toUCSC()), map.get(second.toUCSC()));
					System.err.println(map.get(first.toUCSC()).getValue()+" "+map.get(second.toUCSC()).getValue());
				}
			//}
		}
		
		//put invisible edges
		//These might not be in order
		Object[] exons=nodes.toArray();
		for(int i=0; i<exons.length-1; i++){
			String current=((LightweightGenomicAnnotation)exons[i]).toUCSC();
			String next=((LightweightGenomicAnnotation)exons[i+1]).toUCSC();
			mxCell e1=(mxCell) graph.insertEdge(parent, null, "i", map.get(current), map.get(next), "invis");
			//System.err.println(map.get(current).getValue()+ " "+map.get(next).getValue());
		}
		layout.execute(parent);
	
		
		//fix the x-axis so that they are "genome" based
		for(String key: map.keySet()){
			mxCell v1=map.get(key);
			Alignments exon=new Alignments(key);
			//v1.getGeometry().setX(exon.getStart()-region.getStart());
			//System.err.println(v1.getValue().toString()+" "+(exon.getStart()-region.getStart())+" "+exon.toUCSC());
		}
		
		
		}finally
	{
		graph.getModel().endUpdate();
	}
	
	//mxGraphComponent graphComponent = new mxGraphComponent(graph);
	return graph;
	
	
	}
*/


}
