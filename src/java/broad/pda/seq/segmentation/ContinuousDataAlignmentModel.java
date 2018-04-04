package broad.pda.seq.segmentation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.math.MathException;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.AlignmentBlock;
import org.broad.igv.sam.ReadMate;
import org.broad.igv.sam.SamAlignment;
import org.jgrapht.util.VertexPair;

import broad.core.annotation.BED;
import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.annotation.MaximumContiguousSubsequence;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.math.EmpiricalDistribution;
import broad.core.math.Statistics;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.core.util.CollapseByIntersection;
import broad.pda.annotation.BEDFileParser;
import broad.pda.annotation.GTFFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.GeneTools;
import broad.pda.gene.RefSeqGene;
import broad.pda.gene.RefSeqGeneWithIsoforms;
import broad.pda.geneexpression.dge.DGEFullRNASeq;
import broad.pda.graph.dotUtils.ExtractRegionsFromDOT;
import broad.pda.seq.alignment.AlignmentCollection;
import broad.pda.seq.alignment.AlignmentUtils;
import broad.pda.seq.alignment.MapPairedEnds;
import broad.pda.seq.alignment.MapPairedEndsFromSingleFile;
import broad.pda.seq.graph.ChromosomeWithBubblesJGraphT;
import broad.pda.seq.graph.ChromosomeWithBubblesJGraphT.BubbleEdge;
import broad.pda.seq.graph.ChromosomeWithBubblesJGraphT.EdgeSourceType;
import broad.pda.seq.graph.ChromosomeWithBubblesJGraphT.WindowIterator;
import broad.pda.seq.graph.Path;
import broad.pda.seq.pairedend.PairedEndDistribution;
import broad.pda.seq.pairedend.PairedEndDistribution.PairedEndAnalysisResult;
import broad.pda.seq.pairedend.SplitByPairedEnd;
import broad.pda.seq.pairedend.TreeUtils;
import broad.pda.seq.utils.PolyAUtils;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;


//rewrite of my classic method making use of the direct alignment counts to avoid any voodoo counts coming from overlapping regions and unneccessary 25bp indexing
//Now support RNA-Seq by incorporating a graph structure for segmentation
public class ContinuousDataAlignmentModel implements AlignmentCollection {

	static Logger logger = Logger.getLogger(ContinuousDataAlignmentModel.class.getName());
	AlignmentDataModelStats data;
	Map<String, Integer> chromosomeLengths;
	int extensionFactor=0; //ChIP extension factor
	Collection<SAMRecord> introns;
	int chunkSize=DEFAULT_CHUNK_SIZE;

	double minNumberOfReadsAtEnd=1.01;//can also filter by  median of segment and see what happens
	int minNumberOfSplices;

	Map<String, Integer> maskedRegions;
	double alpha=.05;
	private int minIntronSize = 10;
	//private int maxIntronSize=1000000;
	private static final int DEFAULT_CHUNK_SIZE = 1000000;
	AlignmentDataModelStats pairedData;
	AlignmentDataModelStats strandSpecificReads;
	boolean isStrandSpecific=false;
	boolean upWeightSplices;
	private boolean findMaxContiguous;
	private boolean trimEnds;
	private double trimQuantile = 0.25;
	private int minAnnotationSize = 0;
	public static int DEFAULT_MIN_MAPPING_QUALITY = 5;
	public static int DEFAULT_INSERT_SIZE_FUDGE = 20;
	public static double DEFAULT_INS_SIZE_PVAL = 0.05;
	public static double DEFAULT_MAKE_GENE_OVERLAP = 0.3;

	/**
	 * Constructor to build the data alignment model
	 * @param data the alignment data model that stores the read data
	 * @throws IOException
	 */
	public ContinuousDataAlignmentModel(AlignmentDataModel data) throws IOException{
		this(new AlignmentDataModelStats(data), null, 0);
	}

	/**
	 * Constructor to build the data alignment model
	 * @param data the alignment data model that computes the statistics and stores the read data
	 * @throws IOException
	 */
	public ContinuousDataAlignmentModel(AlignmentDataModelStats data) throws IOException{
		this(data, null, 0);
	}

	/**
	 * Constructor to build the data alignment model
	 * @param data the alignment data model that computes the statistics and stores the read data
	 * @param EF an extension factor used to extend reads
	 * @param minSpliceSupport the minimum number of spliced reads needed to add a connection in the graph
	 * @throws IOException
	 */
	public ContinuousDataAlignmentModel(AlignmentDataModelStats data, int EF, int minSpliceSupport ) throws IOException{
		this(data, null, EF, minSpliceSupport);
	}

	/**
	 * Constructor to build the data alignment model
	 * @param data the alignment data model that computes the statistics and stores the read data
	 * @param maskFiles the mask files
	 * @param minSpliceSupport the minimum number of spliced reads needed to add a connection in the graph
	 * @throws IOException
	 */
	public ContinuousDataAlignmentModel(AlignmentDataModelStats data, File[] maskFiles, int minSpliceSupport )throws IOException{
		this(data, maskFiles, 0, minSpliceSupport);
	}


	/**
	 * Constructor to build the data alignment model
	 * @param data the alignment data model that computes the statistics and stores the read data
	 * @param maskFiles the mask files
	 * @param EF an extension factor used to extend reads
	 * @param minSpliceSupport the minimum number of spliced reads needed to add a connection in the graph
	 * @param pairedData the alignment data model that stores the stats for paired end reads and the paired end read data
	 * @param strandSpecificReads an alignment data model for the strand specific reads
	 * @throws IOException
	 */
	public ContinuousDataAlignmentModel(AlignmentDataModelStats data, File [] maskFiles, int EF, int minSpliceSupport, AlignmentDataModelStats pairedData, AlignmentDataModelStats strandSpecificReads) throws IOException{
		this(data, parseMaskFiles(maskFiles), EF, minSpliceSupport, pairedData, strandSpecificReads, DEFAULT_CHUNK_SIZE, false);		
	}

	/**
	 * Constructor to build the data alignment model
	 * @param data the alignment data model that computes the statistics and stores the read data
	 * @param maskFiles the mask files
	 * @param EF an extension factor used to extend reads
	 * @param minSpliceSupport the minimum number of spliced reads needed to add a connection in the graph
	 * @param pairedData the alignment data model that stores the stats for paired end reads and the paired end read data
	 * @param strandSpecificReads an alignment data model for the strand specific reads
	 * @param chunkSize the size of the cache chunk used
	 * @throws IOException
	 */
	public ContinuousDataAlignmentModel(AlignmentDataModelStats data, File [] maskFiles, int EF, int minSpliceSupport, AlignmentDataModelStats pairedData, AlignmentDataModelStats strandSpecificReads, boolean upweight) throws IOException{
		this(data, parseMaskFiles(maskFiles), EF, minSpliceSupport, pairedData, strandSpecificReads, DEFAULT_CHUNK_SIZE, upweight);		
	}


	/**
	 * Constructor to build the data alignment model
	 * @param data the alignment data model that computes the statistics and stores the read data
	 * @param maskFileData a map of all unalignable regions of the genome
	 * @param EF an extension factor used to extend reads
	 * @param minSpliceSupport the minimum number of spliced reads needed to add a connection in the graph
	 * @param pairedData the alignment data model that stores the stats for paired end reads and the paired end read data
	 * @param strandSpecificReads an alignment data model for the strand specific reads
	 * @throws IOException
	 */
	public ContinuousDataAlignmentModel(AlignmentDataModelStats data, Map<String, Integer> maskFileData, int EF, int minSpliceSupport, AlignmentDataModelStats pairedData, AlignmentDataModelStats strandSpecificReads) throws IOException{
		this(data, maskFileData, EF, minSpliceSupport, pairedData, strandSpecificReads, DEFAULT_CHUNK_SIZE, false);
	}

	/**
	 * Constructor to build the data alignment model
	 * @param data the alignment data model that computes the statistics and stores the read data
	 * @param maskFileData a map of all unalignable regions of the genome
	 * @param EF an extension factor used to extend reads
	 * @param minSpliceSupport the minimum number of spliced reads needed to add a connection in the graph
	 * @param pairedData the alignment data model that stores the stats for paired end reads and the paired end read data
	 * @param strandSpecificReads an alignment data model for the strand specific reads
	 * @param upWeightSplices whether to weight spliced reads based on their proportioncoverage scores
	 * @throws IOException
	 */
	public ContinuousDataAlignmentModel(AlignmentDataModelStats data, Map<String, Integer> maskFileData, int EF, int minSpliceSupport, AlignmentDataModelStats pairedData, AlignmentDataModelStats strandSpecificReads, boolean upWeight) throws IOException{
		this(data, maskFileData, EF, minSpliceSupport, pairedData, strandSpecificReads, DEFAULT_CHUNK_SIZE, upWeight);
	}

	/**
	 * Constructor to build the data alignment model
	 * @param data the alignment data model that computes the statistics and stores the read data
	 * @param maskFileData a map of all unalignable regions of the genome
	 * @param EF an extension factor used to extend reads
	 * @param minSpliceSupport the minimum number of spliced reads needed to add a connection in the graph
	 * @param pairedData the alignment data model that stores the stats for paired end reads and the paired end read data
	 * @param strandSpecificReads an alignment data model for the strand specific reads
	 * @param chunkSize the size of the cache chunk used
	 * @param upWeightSplices whether to weight spliced reads based on their proportion when computing coverage scores
	 * @throws IOException
	 */
	public ContinuousDataAlignmentModel(AlignmentDataModelStats data, Map<String, Integer> maskFileData, int EF, int minSpliceSupport, AlignmentDataModelStats pairedData, AlignmentDataModelStats strandSpecificReads, int chunkSize, boolean upWeightSplices) throws IOException{
		this.strandSpecificReads=strandSpecificReads;
		this.pairedData=pairedData;
		this.data=data;
		this.extensionFactor=EF;
		this.chromosomeLengths=data.getChromosomeLengths();
		this.maskedRegions=maskFileData;
		this.minNumberOfSplices = minSpliceSupport;	
		this.chunkSize = chunkSize;
		data.setChunkSize(this.chunkSize);
		this.upWeightSplices=upWeightSplices;
	}

	/**
	 * Constructor to build the data alignment model
	 * @param data the alignment data model that computes the statistics and stores the read data
	 * @param maskFiles the mask files
	 * @param EF an extension factor used to extend reads
	 * @param minSpliceSupport the minimum number of spliced reads needed to add a connection in the graph
	 * @param pairedData the alignment data model that stores the stats for paired end reads and the paired end read data
	 * @throws IOException
	 */
	public ContinuousDataAlignmentModel(AlignmentDataModelStats data, File[] maskFiles, int EF, int minSpliceSupport, AlignmentDataModelStats pairedData) throws IOException{
		this(data, maskFiles, EF, minSpliceSupport, pairedData, null);
	}

	/**
	 * Constructor to build the data alignment model
	 * @param data the alignment data model that computes the statistics and stores the read data
	 * @param maskFiles the mask files
	 * @param EF an extension factor used to extend reads
	 * @param minSpliceSupport the minimum number of spliced reads needed to add a connection in the graph
	 * @throws IOException
	 */
	public ContinuousDataAlignmentModel(AlignmentDataModelStats data, File[] maskFiles, int EF, int minSpliceSupport ) throws IOException{
		this(data, maskFiles, EF, minSpliceSupport, null);
	}



	public void setFindMaxContiguous(boolean findMaxContiguous) {
		this.findMaxContiguous = findMaxContiguous;

	}

	public void setTrimEnds(boolean trimEnds) {
		this.trimEnds = trimEnds;

	}	

	public void setMinContguousSegmentSize(int minSize) {
		this.minAnnotationSize = minSize;

	}
	public void setAlpha(double alpha) {
		this.alpha = alpha;
	
	}
	
	public void setUpWeightSplices(boolean setFlag) {
		this.upWeightSplices = setFlag; 
	}
	
	public void setMinNumberOfSplices(int minSpliceSupport) {
		this.minNumberOfSplices = minSpliceSupport;
	}
	
	public void setExtensionFactor(int extensionFactor) {
		this.extensionFactor = extensionFactor;
		data.setExtensionFactor(extensionFactor);
	}
	
	public void setMaskFileData(Map<String, Integer> maskedRegionMap) {
		this.maskedRegions = maskedRegionMap;
	}

	public void setMinimumMappingQuality(double minMapQual) { this.data.setMinimumMappingQuality(minMapQual);}

	public void setTrimQuantile(double quantile) {this.trimQuantile = quantile;}

	/**
	 * Get the underlying AlignmentDataModelStats object
	 * @return the AlignmentDataModelStats object
	 */
	public AlignmentDataModelStats getData() {return this.data;}
	
	/**
	 * Get an estimate of the expected read coverage per base
	 * @param chr chromosome to retrieve
	 * @return lambda
	 * @throws IOException 
	 */
	public double getLambda(String chr) throws IOException{
		return data.getLambda(chr);
	}

	/**
	 * Gets the total number of reads per chr
	 * @param chr chromosome to retrieve
	 * @return total number of reads
	 */
	public double getSum(String chr)throws IOException{
		return data.getSum(chr);
	}

	/**
	 * Gets the total number of reads per chr
	 * @param chr chromosome to retrieve
	 * @return total number of reads
	 */
	public double getNumberOfReads(String chr)throws IOException{
		return getSum(chr);
	}

	/**
	 * Gets the number of allowable bases per chr
	 * @param chr chromosome to retrieve
	 * @return total number of allowable bases
	 */
	public  double getNumberMarkers(String chr)throws IOException{
		return data.getNumberMarkers(chr);
	}

	/**
	 * Get total number of bases
	 * @return total number of bases
	 */
	public double getNumberOfBPs(){
		double rtrn=0;
		for(String chr: this.chromosomeLengths.keySet()){
			rtrn+=this.chromosomeLengths.get(chr);
		}
		return rtrn;
	}

	public int getCount(Collection<Alignments> alignments) throws IOException {		
		return data.getCounts(alignments, this.extensionFactor);
	
	}

	// For AlignmentCollection interface - should be moved to Generic model eventually. JE
	@Override
	public int getCount(Alignments alignments) throws IOException {
		Collection<Alignments> a = new ArrayList<Alignments>();
		a.add(alignments);
		return this.getCount(a);
	}
	
	public int[] getGeneCounts(Collection<RefSeqGene> genes)throws IOException{
		return data.getGeneCounts(genes, this.extensionFactor);
	}
	
	public int[] getGeneCountsWithoutDatasetTotals(Collection<RefSeqGene> genes)throws IOException{
		return data.getGeneCountsWithoutDatasetTotals(genes, this.extensionFactor);
	}
	
	
	public int[] getGeneCounts(Map<String, IntervalTree<RefSeqGene>> geneTree, Collection<String> multimappers)throws IOException{
	
		int intergenic=0;
		int exonic=0;
		int intronic=0;
		int UTR5=0;
		int UTR3=0;
		int total=0;
		int splicedReads=0;
	
		for(String chr: geneTree.keySet()){
			CloseableIterator<Alignment> iter=data.getData().getAlignmentsOverlappingRegion(new Alignments(chr, 0, this.getChromosomeLength(chr)));
			while(iter.hasNext()){
				Alignment read=iter.next();
				if(multimappers!=null && multimappers.contains(read.getReadName())){}
				else{
					if(read.getAlignmentBlocks()!=null && read.getAlignmentBlocks().length>1){splicedReads++;}
					Iterator<Node<RefSeqGene>> genes=geneTree.get(chr).overlappers(read.getAlignmentStart(), read.getAlignmentEnd());
					if(genes.hasNext()){
						int[] num=group(read, genes);
						exonic+=num[0]; intronic+=num[1]; UTR5+=num[2]; UTR3+=num[3];
					}
					else{intergenic++;}
					total++;
				}
			}
			iter.close();
		}
	
		int[] rtrn={exonic, intronic, UTR5, UTR3, intergenic, total, splicedReads};
		return rtrn;
	}

	public double[] getCountsPerBp(Alignments align, IntervalTree<Alignment> tree){
		double[] vals=new double[align.getSize()];
	
		for(int i=align.getStart(); i<align.getEnd(); i++){
			double counter=data.getCountsPerAlignment(new Alignments(align.getChr(), i, i), tree, this.extensionFactor); //consider removing extension factor
			vals[i-align.getStart()]=counter;
		}
	
		return vals;
	}

	public List<Double> getCountsPerBp(RefSeqGene gene, IntervalTree<Alignment> tree){
		ArrayList<Double> rtrn=new ArrayList<Double>();
		//double[] vals=new double[align.getSize()];
		//String chr=align.getChr();
	
		for(int i=0; i<gene.getExons().length; i++){
			double[] vals=this.getCountsPerBp(gene.getExons()[i], tree);
			for(int j=0; j<vals.length; j++){rtrn.add(vals[j]);}
		}
	
		return rtrn;
	}

	private List<Double> getNormalizedCountsPerBp(RefSeqGene gene, IntervalTree<Alignment> tree) throws IOException{
		ArrayList<Double> rtrn=new ArrayList<Double>();
		//double[] vals=new double[align.getSize()];
		//String chr=align.getChr();
	
		for(int i=0; i<gene.getExons().length; i++){
			double[] vals=getCountsPerBp(gene.getExons()[i], tree);
			for(int j=0; j<vals.length; j++){rtrn.add(vals[j]/getLambda(gene.getChr()));}
		}
	
		return rtrn;
	}

	public Map<RefSeqGene, double[]> getGeneExpression(Collection<RefSeqGene> genes)throws IOException{
		return data.getGeneExpression(genes, extensionFactor);
	}

	public int getCounts(Map<String, IntervalTree<Alignments>> geneTree, Collection<String> multimappers)throws IOException{
	
		int count=0;
	
		for(String chr: geneTree.keySet()){
			CloseableIterator<Alignment> iter=data.getData().getAlignmentsOverlappingRegion(new Alignments(chr, 0, this.getChromosomeLength(chr)));
			while(iter.hasNext()){
				Alignment read=iter.next();
				if(multimappers!=null && multimappers.contains(read.getReadName())){}
				else{
					Iterator<Node<Alignments>> genes=geneTree.get(chr).overlappers(read.getAlignmentStart(), read.getAlignmentEnd());
					if(genes.hasNext()){count++;}
				}
			}
			iter.close();
		}
	
		return count;
	}

	
	
	public Collection<Alignments> getOverlappingRegions(String chr, Collection<Alignments> introns) throws IOException{
		return data.getOverlappingRegions(chr, introns);
	}

	private RefSeqGene get3PrimeBases(RefSeqGene gene, int size) {
	
		if(gene.getOrientation().equalsIgnoreCase("+")){
			int relativeStart=gene.getTranscriptLength()-size;
			int relativeEnd=gene.getTranscriptLength();
			RefSeqGene align=gene;
			if(relativeStart>0){align=gene.trim(relativeStart, relativeEnd);}
			return align;
		}
		else if(gene.getOrientation().equalsIgnoreCase("-")){
			int relativeStart=0;
			int relativeEnd=size;
			RefSeqGene align=gene;
			if(relativeEnd<gene.getTranscriptLength()){align=gene.trim(relativeStart, relativeEnd);}
			return align;
		}
		else{
			//System.err.println(gene.getAlignment().toUCSC()+" "+gene.getOrientation());
			int relativeStart=gene.getTranscriptLength()-size;
			int relativeEnd=gene.getTranscriptLength();
			RefSeqGene align=gene;
			if(relativeStart>0){align=gene.trim(relativeStart, relativeEnd);}
			return align;
		}
	}

	public int getChromosomeLength(String chr){
		if(this.chromosomeLengths.containsKey(chr)){return this.chromosomeLengths.get(chr);}
		return 0;
	}

	public Map<String, Integer> getChromosomeLengths() {
		return data.getChromosomeLengths();
	}

	public static int[] getWidths(String str){
		String[] vals=str.split(",");
		int[] rtrn=new int[vals.length];
	
		for(int i=0; i<vals.length; i++){
			rtrn[i]=new Integer(vals[i]);
		}
	
		return rtrn;
	}

	private Collection<Alignments> getIntrons(Alignment record){
		return getIntronsFromExons(toExons(record));
	}
	
	private Collection<Alignments> getIntronsFromExons(Collection<Alignments> exons){
		RefSeqGene gene=new RefSeqGene(exons);
		Collection<Alignments> rtrn=gene.getIntronSet();
		return rtrn;
	}
	

	private Collection<Alignments> getIntrons(Collection<RefSeqGene> allPaths) {
		Collection<Alignments> introns=new TreeSet<Alignments>();
	
		for(RefSeqGene gene: allPaths){
			Collection<Alignments> s=gene.getIntronSet();
			for(Alignments intron: s){
				intron.setOrientation(gene.getOrientation());
				introns.add(intron);
			}
		}
	
		return introns;
	}

	private Collection<Alignments> getExons(Collection<RefSeqGene> allPaths) {
		Collection<Alignments> exons=new TreeSet<Alignments>();
	
		for(RefSeqGene gene: allPaths){
			exons.addAll(gene.getExonSet());
		}
	
		return exons;
	}

	//TODO Should be local significance
	//TODO only use reads that dont overlap the good exon
	//TODO trim by using all good exons
	private Map<Alignments, double[]> getSignificantPieces(Alignments exon, Map<String, IntervalTree<Alignments>> goodExonTree, int minLength, Map<String, IntervalTree<RefSeqGene>> genes) throws IOException{
		IntervalTree<Alignment> tree=data.getIntervalTreeCached(exon.getChr(), exon.getStart(), exon.getEnd());
		IntervalTree<RefSeqGene> geneTree = genes.get(exon.getChr());
		Map<Alignments, double[]> rtrn=new TreeMap<Alignments, double[]>();
		if(geneTree!= null) {
			Iterator<Node<RefSeqGene>> exonGeneOverlaperIt = geneTree.overlappers(exon.getStart(), exon.getEnd());
			if(exonGeneOverlaperIt.hasNext()) {
				double localLambda=computeLambda(exonGeneOverlaperIt, tree, exon.getChr());
				Collection<Alignments> split=splitAlignments(exon, goodExonTree);
				for(Alignments align: split){
					if(align.getEnd()>align.getStart() && align.getSize()>minLength){
						double[] scanP=this.scanPRate(align, goodExonTree, tree, localLambda);
						if(scanP[0]<alpha){
							scanP=this.scanPRate(align, goodExonTree, tree, getLambda(exon.getChr()));
							if(scanP[0]<alpha){rtrn.put(align, scanP);}
						}
					}
				}
			}
		}
		return rtrn;
	}

	private int getPathCoverageData(Path path, int numIntrons, TreeMap<Alignments, Double> exonMap, TreeMap<Alignments, Double> intronMap, Path overlappingPath) {
		RefSeqGene gene=overlappingPath.toGene();
		Set<Alignments> exons = gene.getExonSet();
		for(Alignments exon : exons) {
			exonMap.put(exon, path.getNodeCount(exon));
		}
		Collection<BubbleEdge> introns = overlappingPath.getEdges();
		for(BubbleEdge intron : introns) {
			intronMap.put(new Alignments(intron.getConnection()), intron.getSplicedCounts());
		}			
		numIntrons=Math.max(numIntrons, gene.getNumExons()-1);
		return numIntrons;
	}

	private Map<Alignments, Double> getNodeCounts(Collection<Alignments> decollapsed, Map<Alignments, Integer> exons) throws IOException{
		Map<Alignments, Double> rtrn=new TreeMap<Alignments, Double>();
	
		Map<String, IntervalTree<Alignments>> trees=CollapseByIntersection.makeIntervalTree(exons.keySet());
	
		Collection<Alignments> remainder=new TreeSet<Alignments>();
	
		for(Alignments exon: decollapsed){
			IntervalTree<Alignments> exonChrTree = trees.get(exon.getChr());
			if(exonChrTree != null) {
				Iterator<Node<Alignments>> iter= exonChrTree.overlappers(exon.getStart(), exon.getEnd());
				int i=0;
				double count=0;
				while(iter.hasNext()){
					Alignments align=iter.next().getValue();
					count=exons.get(align);
					i++;
				}
				if(i>1){
					//compute from scratch
					remainder.add(exon);
					//System.out.println(exon);
				}
				rtrn.put(exon, count);
			}
		}
	
		Map<Alignments, Double> counter=countExons(remainder);
		rtrn.putAll(counter);
	
		return rtrn;
	}

	private Map<Alignments, double[]> getDataForAlignments(Collection<Alignments> alignments) throws IOException{
		Map<Alignments, double[]> rtrn=new TreeMap<Alignments, double[]>();
	
		//cache in memory by chunks
		IntervalTree<Alignment> tree=null;
		int sizeContains=0;
		for(Alignments align: alignments){
	
			if(tree==null){tree=data.getIntervalTree(align.getChr(), align.getStart(), align.getStart()+this.chunkSize); sizeContains=align.getStart()+this.chunkSize;}
			else if(align.getStart()>=sizeContains || align.getEnd()>=sizeContains){tree=data.getIntervalTree(align.getChr(), align.getStart(), align.getEnd()+this.chunkSize); sizeContains=align.getEnd()+this.chunkSize;}
			double[] vals=getCountsPerBp(align, tree);
			rtrn.put(align, vals);
	
		}
	
	
		return rtrn;
	}

	/**
	 * This function returns an array of type double with counts for each base pair across the genomic region of the gene
	 * @author skadri
	 * @param align
	 * @return
	 * @throws IOException
	 */
	public double[] getDataForAlignment(Alignments align) throws IOException{	
		
		IntervalTree<Alignment> tree=data.getIntervalTree(align.getChr(), align.getStart(), align.getEnd()); //sizeContains=align.getStart()+this.chunkSize;
		
		double[] vals=getCountsPerBp(align, tree);
		
		return vals;
	}
	
	public List<Double> getDataForGene(RefSeqGene gene) throws IOException{
		
		IntervalTree<Alignment> tree=data.getIntervalTreeCached(gene.getChr(), gene.getStart(), gene.getEnd());
				
		List<Double> vals=getCountsPerBp(gene, tree);
		
		return vals;
	}	
	
	public Map<RefSeqGene, List<Double>> getDataForGene(Collection<RefSeqGene> alignments) throws IOException{
		Map<RefSeqGene, List<Double>> rtrn=new TreeMap<RefSeqGene, List<Double>>();
	
		//TODO Implement a caching scheme for the tree
	
		long treeTime=0;
		//long sortedTime=0;
		long countTime=0;
	
		//cache in memory by chunks
		for(RefSeqGene align: alignments){
			long start=System.currentTimeMillis();
			IntervalTree<Alignment> tree=data.getIntervalTreeCached(align.getChr(), align.getStart(), align.getEnd());
			long end=System.currentTimeMillis();
			treeTime+=(end-start);
	
			/*start=System.currentTimeMillis();
			Collection<Alignments> exons=align.getSortedAndUniqueExons();
			end=System.currentTimeMillis();
			sortedTime+=(end-start);
	*/
			start=System.currentTimeMillis();
			List<Double> vals=getCountsPerBp(align, tree);
			end=System.currentTimeMillis();
			countTime+=(end-start);
	
			rtrn.put(align, vals);
		}
	
		System.err.println("Tree Time: "+treeTime);
		//System.err.println("Sorted Time: "+sortedTime);
		System.err.println("Count Time: "+ countTime);
	
		//data.resetTreeCache();
	
		return rtrn;
	}

	/**
	 * Return all exons that were originally obtained by read assembly but are extended past intron
	 * @param exons
	 * @param decollapsed
	 * @param graph
	 * @param minDistance
	 * @param minSpliceFreq
	 * @return
	 * @throws IOException
	 */
	private Collection<Alignments> getExtendedExons(Collection<Alignments> exons, Collection<Alignments> decollapsed, ChromosomeWithBubblesJGraphT graph, int minDistance, double minSpliceFreq) throws IOException{
		//int distance=76; //TODO: For now make it the read length, might be able to get rid of altogether
		Map<String, IntervalTree<Alignments>> trees=CollapseByIntersection.makeIntervalTree(decollapsed);
	
		Collection<Alignments> rtrn=new TreeSet<Alignments>();
		
		Map<String, IntervalTree<RefSeqGene>> genes=CollapseByIntersection.makeIntervalTreeForGenes(graph.getGenePaths(1, minSpliceFreq));
	
		for(Alignments exon: exons){
			if(!decollapsed.contains(exon)){
				Iterator<Node<Alignments>> overlappers=trees.get(exon.getChr()).overlappers(exon.getStart(), exon.getEnd());
				if(overlappers==null || !overlappers.hasNext()){rtrn.add(exon);}
				else{
					/*while(overlappers.hasNext()){
						Alignments goodExon=overlappers.next().getValue();
						Map<Alignments, double[]> significantPieces=getSignificantPieces(exon, goodExon, distance);
						rtrn.addAll(significantPieces.keySet());
						//int uniqueLength=goodExon.getUniqueLength(exon);
						//if(uniqueLength<0 || uniqueLength>distance){rtrn.add(exon);}
					}*/
					Map<Alignments, double[]> significantPieces=getSignificantPieces(exon, trees, minDistance, genes);
					//TODO add these pieces to the actual good exon
					Collection<Alignments> additionalNodes=makeAdditionalNodes(significantPieces.keySet(), trees);
					rtrn.addAll(additionalNodes);
	
				}
			}
		}
	
		return rtrn;
	}

	public AlignmentDataModelStats getPairedData() {
		return this.pairedData;
	}

	private Collection<Alignments> scan(int windowSize, double alpha,  String chr)throws IOException{
		double T=getNumberMarkers(chr);
	
		long start=System.currentTimeMillis();
		int criticalValue=calculateCriticalValue(new Double(T).intValue(), windowSize, T, alpha, getLambda(chr));
		long end=System.currentTimeMillis();
	
		logger.info("Computing critical values  for window size "+ windowSize+" took: "+(end-start)/1000.0 + " sec. " +criticalValue);
	
		start=System.currentTimeMillis();
		IntervalTree<Alignments> significantWindows = scanGenome(windowSize, criticalValue, chr);
		end=System.currentTimeMillis();
	
		logger.info("Scanning window NEW took: "+(end-start)/1000.0 +  " sec.");
	
	
		start=System.currentTimeMillis();
		Collection<Alignments> windows = dedup(significantWindows);
	
		logger.info("Going to findMaxContiguous segments " + findMaxContiguous + " trim ends? " + trimEnds);
	
		if(findMaxContiguous) {
			logger.info("Finding max contiguous regions within windows ... ");
			windows =findMaxContiguous(windows);
			logger.info("Done");
		}
	
		if(trimEnds) {
			logger.info("Trimming window ends ... ");
			windows = trimEnds(windows);
			end=System.currentTimeMillis();
			logger.info(" took: "+(end-start)/1000.0);
			logger.info(" Done");
		}
	
	
		return windows;
	}

	public Collection<Alignments> scan(int[] windowSizes, double alpha) throws IOException{
		Collection<Alignments> rtrn=new TreeSet<Alignments>();
	
		for(String chr: this.chromosomeLengths.keySet()){
			Collection<Alignments> set=this.scan(windowSizes, alpha, chr);
			rtrn.addAll(set);
		}
	
		return rtrn;
	}

	public Collection<Alignments> scan(int[] windowSizes, double alpha,  String chr) throws IOException{
		Collection<Alignments> rtrn=new TreeSet<Alignments>();
	
		double lambdaVal=getLambda(chr);
		if(lambdaVal>0){
			for(int i=0; i<windowSizes.length; i++){
				Collection<Alignments> list=scan(windowSizes[i], alpha, chr);
				rtrn.addAll(list);
			}
		}
	
		System.gc();
	
	
		return rtrn;
	}

	
	public void scanGenome(ContinuousDataAlignmentModel data2, int windowSize, String chr, FileWriter writer, boolean filterSignificance) throws IOException{
		int counter=0;
	
		double numMarkers=getNumberMarkers(chr);
	
		
		//load chunks into memory: first chunk
		int chunkNumber=1;
	
		IntervalTree<Alignment> chunkAlignmentTree=data.getIntervalTree(chr, 0, chunkNumber*chunkSize);
		IntervalTree<Alignment> chunkAlignmentTree2=data2.getData().getIntervalTree(chr, 0, chunkNumber*chunkSize);
	
	
	
		double score1=0;
		double score2=0;
		
		Alignments previous=null;
		boolean cached=false;
	
		for(int i=0; i<data.getChromosomeLengths().get(chr); i++){
	
			
			//if score is 0 jump to end of fixed width
			int start=i;
			int end=start+windowSize;
			Alignments current=new Alignments(chr, start, end);
			//System.err.println("Current: " + current.toUCSC());
			double sum=0;
			double sum2=0;
	
			if(score1==0 && score2==0){
				sum=data.getCountsPerAlignment(current, chunkAlignmentTree, extensionFactor); 
				sum2=data2.getData().getCountsPerAlignment(current, chunkAlignmentTree2, extensionFactor);
					
				if(sum ==0 && sum2==0) {
					i=start+windowSize;
					//System.err.println("\tScore was 0 and so was sum, advancing from " + start + " to " +i);
				} else {
					//System.err.println("\tScore was 0, used reads to compute sum: " + sum);
				}
				cached=false;
			}else{
				/*
				if have part of the interval scored and cached then to get next score just need to compute piece not contained
				get score for base not covered and add
				get score for base covered in previous that no longer contained and subtract
				 */		
				cached = true;
				Alignments startPosition = new Alignments(chr, previous.getStart(), current.getStart());
				double subtractVal=data.getCountsOfUniqueOverlappers(startPosition, current, chunkAlignmentTree, extensionFactor);
				double subtractVal2=data2.getData().getCountsOfUniqueOverlappers(startPosition, current, chunkAlignmentTree2, extensionFactor);
				Alignments endPosition = new Alignments(chr, previous.getEnd(), current.getEnd());
				double addVal=data.getCountsOfUniqueOverlappers(endPosition, previous, chunkAlignmentTree, extensionFactor);
				double addVal2=data2.getData().getCountsOfUniqueOverlappers(endPosition, previous, chunkAlignmentTree2, extensionFactor);
				sum=(score1-subtractVal)+addVal;
				sum2=(score2-subtractVal2)+addVal2;
				//System.err.println("\tp: " + previous.toUCSC() + " c: " +current.toUCSC()+ " " +data.getCountsPerAlignment(current, chunkAlignmentTree, extensionFactor)+" "+subtractVal+" "+addVal+" "+data.getCountsPerAlignment(previous, chunkAlignmentTree, extensionFactor)+" score: "+score + " sum: " + sum);
			}
	
			double percent=(i/numMarkers)*100;
			if(counter% 1000000 ==1 || sum < 0){
				double memoryPercent=Runtime.getRuntime().freeMemory()/(double)Runtime.getRuntime().totalMemory();
				//System.err.println(counter+" % of markers done "+percent+", % memory used  "+(memoryPercent*100)+", window: "+window.toUCSC()+", sum: "+sum);
			}
			//assert(sum > 0);
			score1=sum;
			score2=sum2;
			previous=current;
	
			int midPoint=current.getMidPoint();
			
			double ratio=(sum+1)/(sum2+1);
			
			double scaledRatio=ratio;
			if(ratio<1){scaledRatio=(-1.0/ratio);}
			
			//double logratio=Math.log(ratio)/Math.log(2);
			if((sum==0 && sum2==0)){scaledRatio=0;}
			if(sum<0 || sum2<0){}
			//else{writer.write(current.getMidPoint()+"\t"+scaledRatio+"\n");}
			else if(midPoint<data.getChromosomeLengths().get(chr)){
				double p=calculatePVal(new Double(sum2).intValue(), getLambda(chr), windowSize, getNumberMarkers(chr));
				if(!filterSignificance || p<alpha){
					writer.write(current.getChr()+"\t"+midPoint+"\t"+(midPoint+1)+"\t"+scaledRatio+"\n");
				}
			}
			
			
			if(end>=(chunkNumber*chunkSize-1)){
				chunkNumber++; 
				chunkAlignmentTree=data.getIntervalTree(chr, i, Math.max(chunkNumber*this.chunkSize, end));
				chunkAlignmentTree2=data2.getData().getIntervalTree(chr, i, Math.max(chunkNumber*this.chunkSize, end));
			}
	
	
			counter++;
		}
	
	}
	
	// to compute score take current score and subtract base pair before and add base pair after ensuring that none of the added is
	// already overlapping the segment and none of the subtracted is in segment
	// iterate from 0 to numMarkers-w
	// reverse iterate from numMarkers to w
	private IntervalTree<Alignments> scanGenome(int fixedWidth, int critVal, String chr) throws IOException{
		int counter=0;
	
		// chr size minus masked regions
		double numMarkers=getNumberMarkers(chr);
	
		IntervalTree<Alignments> rtrnTree=new IntervalTree<Alignments>();
	
		//load chunks into memory: first chunk
		int chunkNumber=1;
	
		// get interval tree for the first chunk
		// key is mapped coordinates of read
		// value is number of reads mapped to position
		IntervalTree<Alignment> chunkAlignmentTree=data.getIntervalTree(chr, 0, chunkNumber*chunkSize);
		
		double score=0;
		Alignments previous=null;
		boolean cached=false;
	
		// for every position in chromosome
		for(int i=0; i<data.getChromosomeLengths().get(chr); i++) {
				
			long startTime=System.currentTimeMillis();
			
			//if score is 0 jump to end of fixed width
			int start=i;
			int end=start+fixedWidth;
			Alignments current=new Alignments(chr, start, end); // the region of length fixedWidth beginning at start
			//System.err.println("Current: " + current.toUCSC());
			double sum=0;
	
			if(score==0) {
				
				// set sum to the number of mappings in current interval
				//sum=data.getCountsPerAlignment(current, chunkAlignmentTree, extensionFactor); 
				sum=data.getCountsPerAlignment(current, chunkAlignmentTree, 0); //9/9/12 re-did handling of extension factors at the get tree level 
	
				// if there are no reads in the current interval, skip to next interval
				if(sum ==0) {
					i=start+fixedWidth;
					//System.err.println("\tScore was 0 and so was sum, advancing from " + start + " to " +i);
				} else {
					//System.err.println("\tScore was 0, used reads to compute sum: " + sum);
				}
				cached=false; // no part of the interval has been scored
			} else {
				/*
				if have part of the interval scored and cached then to get next score just need to compute piece not contained
				get score for base not covered and add
				get score for base covered in previous that no longer contained and subtract
				 */		
				cached = true; // part of the interval has been scored
				// interval between start position of previous interval and start position of current interval - the part that has been covered
				Alignments startPosition = new Alignments(chr, previous.getStart(), current.getStart());
				// the number of mappings in current that have been counted
				double subtractVal=data.getCountsOfUniqueOverlappers(startPosition, current, chunkAlignmentTree, 0); //9/9/12 re-did handling of extension factors at the get tree level
				// interval between end position of previous interval and end position of current interval - the part that has not been covered
				Alignments endPosition = new Alignments(chr, previous.getEnd(), current.getEnd());
				// the number of mappings in current that have not been counted
				double addVal=data.getCountsOfUniqueOverlappers(endPosition, previous, chunkAlignmentTree, 0);//9/9/12 re-did handling of extension factors at the get tree level
				
				// score is the sum of previous interval
				// sum is the score of current interval
				sum=(score-subtractVal)+addVal;
				//System.err.println("\tp: " + previous.toUCSC() + " c: " +current.toUCSC()+ " " +data.getCountsPerAlignment(current, chunkAlignmentTree, extensionFactor)+" "+subtractVal+" "+addVal+" "+data.getCountsPerAlignment(previous, chunkAlignmentTree, extensionFactor)+" score: "+score + " sum: " + sum);
			}
	
			//double percent=(i/numMarkers)*100;
			//if(counter% 1000000 ==1 || sum < 0) {
				//double memoryPercent=Runtime.getRuntime().freeMemory()/(double)Runtime.getRuntime().totalMemory();
				//System.err.println(counter+" % of markers done "+percent+", % memory used  "+(memoryPercent*100)+", window: "+window.toUCSC()+", sum: "+sum);
			//}
			//assert(sum > 0);
			score=sum;
			previous=current;
	
			/***********BAD For Testing*********************/
			//if(sum>critVal){rtrn.add(align);}
			if(sum>critVal){
				// merge current interval into rtrnTree 
				Iterator<Node<Alignments>> iter=rtrnTree.overlappers(current.getStart(), current.getEnd());
				rtrnTree=mergeAndRemove(iter, current, rtrnTree);
			}
			/**********************************************/
	
			// move on to next chunk
			if(end>=(chunkNumber*chunkSize-1)){chunkNumber++; chunkAlignmentTree=data.getIntervalTree(chr, i, Math.max(chunkNumber*this.chunkSize, end));}
	
	
			//System.out.println("Free memory after iteration: " + i + " " + Runtime.getRuntime().freeMemory());
			counter++;
		}
	
		return rtrnTree;
	}


	private Map<Path, double[]> scanGenome(int fixedWidth, int critVal, String chr, ChromosomeWithBubblesJGraphT graph, Sequence chrSeq, double minimumSpliceFrequencyToFollowEdge)throws IOException{
	
		//TODO Ensure that transcript graph contains all node and edge counts
		ChromosomeWithBubblesJGraphT transcriptGraph=acrossGraph (critVal, fixedWidth, graph, chr, chrSeq, minimumSpliceFrequencyToFollowEdge);
	
		Map<Path, double[]> pathScores=this.augmentAndScoreTranscriptGraph(transcriptGraph, chr,  minimumSpliceFrequencyToFollowEdge);
	
		return pathScores;
	}

	private double[] scanPRate(LightweightGenomicAnnotation first)throws IOException{
		return data.scanPRate(first, 0);////9/9/12 re-did handling of extension factors at the get tree level
	}

	private double[] scanPRate(Alignments align, IntervalTree<Alignment> tree)throws IOException{
		String chr=align.getChr();
		double sum=data.getCountsPerAlignment(align, tree, 0);//consider getting number of unique reads as well //9/9/12 re-did handling of extension factors at the get tree level
		int count=align.getSize();
		double enrich=(sum/count)/getLambda(chr);
		//System.err.println("sum " + sum + "  - regionsize " + count + "  - lamgda Chr " + getLambda(chr) + " -  avg " + (sum/count));
		double[] rtrn={calculatePVal(new Double(sum).intValue(), getLambda(chr), count, getNumberMarkers(chr)), enrich, sum, (sum/count)};
		return rtrn;
	}

	public double[] scanPRate(RefSeqGene gene)throws IOException{
		return data.scanPRate(gene, 0);//9/9/12 re-did handling of extension factors at the get tree level
	}

	public double[] scanPRate(RefSeqGene gene, IntervalTree<Alignment> tree)throws IOException{
		return data.scanPRate(gene, tree, 0); //9/9/12 re-did handling of extension factors at the get tree level
	}

	public double[] scanPRateStranded(RefSeqGene gene, IntervalTree<Alignment> tree)throws IOException{
		String chr=gene.getChr();
		double sum=0;
		int geneLength=0;
	
		sum=data.getData().getCountsPerAlignmentStranded(gene, tree, 0, gene.getOrientation());//9/9/12 re-did handling of extension factors at the get tree level
		geneLength=gene.getTranscriptLength();
	
		double avgCoverage = sum/(double) geneLength;
		double enrich=avgCoverage/getLambda(chr);
		double[] rtrn={calculatePVal(new Double(sum).intValue(), getLambda(chr), geneLength, getNumberMarkers(chr)), enrich, avgCoverage, avgCoverage * data.getRPKMConstant(chr)};
		return rtrn;
	}

	public double[] scanPRate(RefSeqGene gene, IntervalTree<Alignment> tree, double localLambda) throws IOException{
		return data.scanPRate(gene, tree, localLambda, 0);//9/9/12 re-did handling of extension factors at the get tree level
	
	}

	public double[] scanPRate(RefSeqGene gene, IntervalTree<Alignment> tree, double localLambda, int numMarkers) throws IOException{
		return data.scanPRate(gene, tree, localLambda, 0, numMarkers); //9/9/12 re-did handling of extension factors at the get tree level
	}

	private double[] scanPRate(Alignments align,Map<String, IntervalTree<Alignments>> goodExonTree,	IntervalTree<Alignment> tree, double localLambda) throws IOException {
	
	
		String chr=align.getChr();
		//Get scores for align that dont overlap goodExonTree
		double sum=data.getCountsPerAlignment(align, goodExonTree, tree, 0);//9/9/12 re-did handling of extension factors at the get tree level//consider getting number of unique reads as well
		int count=align.getSize();
		double enrich=(sum/count)/localLambda;
	
		double[] rtrn={calculatePVal(new Double(sum).intValue(), localLambda, count, getNumberMarkers(chr)), enrich, sum, (sum/count)};
		return rtrn;
	}

	public Map<Path, double[]> scanFromGraph(int windowSize, double alpha, boolean trimEnds, String chr, int start, int end, Sequence chrSequence, boolean filterSplice,  double minimumSpliceFrequencyToFollowEdge) throws IOException{
		this.alpha=alpha;
	
		double T=getNumberMarkers(chr);
		double lambda = getLambda(chr);
		int criticalValue=calculateCriticalValue(new Double(T).intValue(), windowSize, T, alpha, lambda);
	
		ChromosomeWithBubblesJGraphT graph=new ChromosomeWithBubblesJGraphT(chr, start, end);
		graph.loadBubbles(data, chrSequence, minIntronSize , minNumberOfSplices, filterSplice);
	
		Map<Path, double[]> pathScores = scanGenome(windowSize, criticalValue, chr, graph, chrSequence, minimumSpliceFrequencyToFollowEdge);
	
		return pathScores; 
	}

	public double[] scanRegion(int fixedWidth, LightweightGenomicAnnotation region) throws IOException{
		return data.scanRegion(fixedWidth, region);	
	}

	public double[] scoreSegment(Alignments align) throws IOException{
		//System.err.println("Scoring segment: " + align.toUCSC());
		IntervalTree<Alignment> tree=data.getIntervalTreeCached(align.getChr(), align.getStart(), align.getEnd());
		double[] pRate=scanPRate(align, tree);
	
		return pRate;
	}
	

	public Map<Alignments, double[]> scoreSegments(Collection<Alignments> set) throws IOException{
		Map<Alignments, double[]> rtrn=new TreeMap<Alignments, double[]>();
	
		//Keep interval tree to speed it up
		IntervalTree<Alignment> tree=null;
		String chr="";
	
		int i=0;
		for(Alignments align: set){
			if(!chr.equalsIgnoreCase(align.getChr())){chr=align.getChr(); tree=data.getIntervalTree(chr, 0, this.chromosomeLengths.get(chr));}
			//if(i%10000 ==0){System.err.println(i+" "+align.toUCSC());}
			double[] pRate=scanPRate(align, tree);
			rtrn.put(align, pRate);
			i++;
		}
		return rtrn;
	}

	public Map<Alignments, double[]> scoreSegments(Collection<Alignments> set, String chrToUse) throws IOException { return scoreSegments(set, chrToUse, false); }
	public Map<Alignments, double[]> scoreSegments(Collection<Alignments> set, String chrToUse, boolean ignoreAlpha) throws IOException{
		Map<Alignments, double[]> rtrn=new TreeMap<Alignments, double[]>();
	
		int i=0;
		int j=0;
		for(Alignments align: set){
			j++;
			if(chrToUse==null || chrToUse.equalsIgnoreCase(align.getChr())){
				if(this.chromosomeLengths.containsKey(align.getChr())){
					IntervalTree<Alignment> tree=data.getIntervalTreeCached(align.getChr(), align.getStart()-1, align.getEnd()+1);
					double[] pRate=scanPRate(align, tree);
					//align.setCountScore(pRate[4]);//setting rpkm as gene scores
					if (ignoreAlpha || pRate[0] < alpha) {
						rtrn.put(align, pRate);
						i++;
					}
				}
			}
		}
		//System.err.println("Processed " + j + " but only " + i + " passed alpha " + alpha);
		data.resetTreeCache();
	
		return rtrn;
	}

	/**
	 * 
	 * @param align
	 * @return {calculatePVal(new Double(sum).intValue(), localLambda, geneLength, getNumberMarkers(chr)), enrich, sum, avgCoverage, rpkm, localLambda, geneLength,nominalP, fullyContained}
	 * @throws IOException
	 */
	public double[] scoreGene(RefSeqGene align) throws IOException {
		double [] pRate = null;
		if(this.chromosomeLengths.containsKey(align.getChr())){
			IntervalTree<Alignment> tree=data.getIntervalTreeCached(align.getChr(), align.getStart()-1, align.getEnd()+1);
			pRate=scanPRate(align, tree);
			align.setBedScore(pRate[4]);//setting rpkm as gene scores
		}
		return pRate;
	}
	
	public double[] scoreGene(RefSeqGene align, double localLambda) throws IOException {
		double [] pRate = null;
		if(this.chromosomeLengths.containsKey(align.getChr())){
			IntervalTree<Alignment> tree=data.getIntervalTreeCached(align.getChr(), align.getStart()-1, align.getEnd()+1);
			pRate=scanPRate(align, tree, localLambda);
			align.setBedScore(pRate[4]);//setting rpkm as gene scores
		}
		return pRate;
	}
	
	public GeneScore scoreGene(RefSeqGene align, IntervalTree<Alignment> tree, double localLambda) throws IOException {
		double [] pRate = null;
		if(this.chromosomeLengths.containsKey(align.getChr())){
			pRate=scanPRate(align, tree, localLambda);
			align.setBedScore(pRate[4]);//setting rpkm as gene scores
		}
		return new GeneScore(align, pRate);
	}
	
	public GeneScore scoreGene(RefSeqGene align, IntervalTree<Alignment> tree, double localLambda, int numMarkers) throws IOException {
		double [] pRate = null;
		if(this.chromosomeLengths.containsKey(align.getChr())){
			pRate=scanPRate(align, tree, localLambda, numMarkers);
			align.setBedScore(pRate[4]);//setting rpkm as gene scores
		}
		return new GeneScore(align, pRate);
	}
	
	public Map<RefSeqGene, double[]> scoreGenes(Collection<RefSeqGene> set, String chrToUse, IntervalTree<Alignment> tree) throws IOException{
		Map<RefSeqGene, double[]> rtrn=new TreeMap<RefSeqGene, double[]>();

		int i=0;
		for(RefSeqGene align: set){
			//System.err.println("Scoring "+ align.getName() +" " + align.toUCSC());
			if(chrToUse==null || chrToUse.equalsIgnoreCase(align.getChr())){
				double[] pRate = scoreGene(align, tree);
				if(pRate != null)
				rtrn.put(align, pRate);
				i++;
			}
		}

		data.resetTreeCache();

		return rtrn;
	}

	public GeneScore getPeak(RefSeqGene gene, RefSeqGene startCodon, IntervalTree<Alignment> tree) throws IOException {
		//TODO Lets get the width of the peak over the startCodon
		//Get all reads overlapping this region
		//then define width as Gene going from start of first read to end of last
		//set GeneScore to this count and region
		RefSeqGene peak=data.getPeak(gene, startCodon, tree, 0);
		if(peak!=null){
			double[] scores=this.scoreGene(peak, tree);
			//System.err.println("Peak");
			//System.err.println(peak);
			return new GeneScore(peak, scores);
		}
		else{
			double[] scores=this.scoreGene(startCodon, tree);
			return new GeneScore(startCodon, scores);
		}
	}
	
	
	
	
	public RefSeqGene getPeak(RefSeqGene gene, RefSeqGene startCodon) throws IOException {
		IntervalTree<Alignment> tree=data.getIntervalTree(gene.getChr(), gene.getStart(), gene.getEnd());
		return data.getPeak(gene, startCodon, tree, 0);
	}
	
	
	
	public AlignmentDataModelStats getAlignmentDataModelStats(){return this.data;}
	
	public double[] scoreGene(RefSeqGene align, IntervalTree<Alignment> tree) throws IOException {
		double [] pRate=scanPRate(align, tree);
		align.setBedScore(pRate[4]);//setting rpkm as gene scores
		return pRate;
	}

	public Map<RefSeqGene, double[]> scoreGenes(Collection<RefSeqGene> genes) throws IOException{
		return scoreGenes(genes, null);
	}

	public Map<RefSeqGene, double[]> scoreGenes(Collection<RefSeqGene> set, String chrToUse) throws IOException{
		Map<RefSeqGene, double[]> rtrn=new TreeMap<RefSeqGene, double[]>();
	
		int i=0;
		for(RefSeqGene align: set){
			//logger.info("Scoring "+ align.getName() +" " + align.toUCSC());
			if(chrToUse==null || chrToUse.equalsIgnoreCase(align.getChr())){
				double[] pRate = scoreGene( align);
				if(pRate != null)
					rtrn.put(align, pRate);
				i++;
			}
		}
	
		data.resetTreeCache();
	
		return rtrn;
	}
	
	public Map<RefSeqGene, double[]> scoreGenes(Collection<RefSeqGene> set, String chrToUse, int size) throws IOException{
		Map<RefSeqGene, double[]> rtrn=new TreeMap<RefSeqGene, double[]>();
	
		//Keep interval tree to speed it up
	
		int i=0;
		for(RefSeqGene gene: set){
			RefSeqGene align=get3PrimeBases(gene, size);
			if(chrToUse==null || chrToUse.equalsIgnoreCase(align.getChr())){
				IntervalTree<Alignment> tree=data.getIntervalTreeCached(align.getChr(), align.getStart()-1, align.getEnd()+1);
				double[] pRate=scanPRate(align, tree);
				rtrn.put(align, pRate);
				i++;
			}
		}
	
		return rtrn;
	}

	public void scoreGenes(BEDFileParser genes) throws IOException{
		for (String chr : getChromosomeLengths().keySet()) {
			if(genes.containChr(chr) & genes.getChrTree(chr) != null ) {
				Iterator<RefSeqGeneWithIsoforms> chrGeneIt = genes.getChrTree(chr).valueIterator();
				while(chrGeneIt.hasNext()) {
					RefSeqGeneWithIsoforms gene = chrGeneIt.next();
					Iterator<RefSeqGene> isoIt = gene.getAllIsoforms().iterator();
					while(isoIt.hasNext()) {
						RefSeqGene iso = isoIt.next();
						double [] pRate = scoreGene(iso);
						iso.setBedScore(pRate[4]);
						iso.setExtraFields(pRate);
					}
				}

			} else {
				logger.info("No gene found for chromosome " + chr);
			}
			resetTreeCache();
		}
	}
	
	
	public Map<RefSeqGene, double[]> scoreGenesStranded(Collection<RefSeqGene> genes) throws IOException{
		return scoreGenesStranded(genes, null);
	}

	public Map<RefSeqGene, double[]> scoreGenesStranded(Collection<RefSeqGene> set, String chrToUse) throws IOException{
		Map<RefSeqGene, double[]> rtrn=new TreeMap<RefSeqGene, double[]>();
	
		//Keep interval tree to speed it up
	
		int i=0;
		for(RefSeqGene align: set){
			if(chrToUse==null || chrToUse.equalsIgnoreCase(align.getChr())){
				//System.err.println(align);
				IntervalTree<Alignment> tree=data.getIntervalTreeCached(align.getChr(), align.getStart()-1, align.getEnd()+1);
				double[] pRate=scanPRateStranded(align, tree);
				align.setCountScore(pRate[2]);
				rtrn.put(align, pRate);
				i++;
			}
		}
	
		data.resetTreeCache();
	
		return rtrn;
	}

	public Map<Path, double[]> scorePathsForOnlyRPKM(Collection<Path> paths) {
		Map<Path, double[]> rtrn = new TreeMap<Path, double[]>();
	
		for (Path path : paths) {
			double score = path.getCoverage();
			double size = path.getSize();
			double avgCoverage = score / size;
			double rpkm = avgCoverage * data.getRPKMConstant(path.getChromosome());
	
			double[] array = {0, 0, score, avgCoverage, rpkm, 0, size, 0};
			rtrn.put(path, array);
		}
	
		return rtrn;
	}

	private static Map<Path, double[]> scorePathsForOnlyRPKM(Collection<Path> paths, double rpkmConstant) {
		Map<Path, double[]> rtrn = new TreeMap<Path, double[]>();
	
		for (Path path : paths) {
			double score = path.getCoverage();
			double size = path.getSize();
			double avgCoverage = score / size;
			double rpkm = avgCoverage * rpkmConstant;
	
			double[] array = {0, 0, score, avgCoverage, rpkm, 0, size, 0};
			rtrn.put(path, array);
		}
	
		return rtrn;
	}

	public Map<Path, double[]> scorePaths(Collection<Path> paths, double lambda) throws IOException{
		Map<Path, double[]> rtrn=new TreeMap<Path, double[]>();
	
	
		int i=0;
		for(Path path: paths){
			double score=path.getCoverage();
			int size=path.getSize();
			double scanP=calculateApproximatePVal(new Double(score).intValue(), lambda, size, data.getNumberMarkers(path.getChromosome()), this.alpha);
			double localP=calculateApproximatePVal(new Double(score).intValue(), path.getLocalLambda(), size, data.getNumberMarkers(path.getChromosome()), this.alpha);
	
	
			double avgCoverage = score/size;
	
			double enrich=avgCoverage/lambda;
			double rpkm = avgCoverage * data.getRPKMConstant(path.getChromosome());
			double[] array={scanP, enrich, score, avgCoverage, rpkm, lambda, size, localP};
	
			rtrn.put(path, array);
			i++;
			//if(i%1000 ==0){System.err.println(i+" "+paths.size());}
		}
	
		return rtrn;
	}

	public static Map<Path, double[]> scorePaths(Collection<Path> paths, ChromosomeWithBubblesJGraphT graph, double alpha) {
		return scorePaths(paths, graph.getLambda(), graph.getNumberOfMarkers(), graph.getNumberOfReads(), graph.getRPKMConstant(), graph.getLocalRate(), alpha);
	}

	public static Map<Path, double[]> scorePaths(Collection<Path> paths, double lambda, double numMarkers, double numReads, double RPKMContant, Map<Alignments, Double> localRates, double alpha){
		Map<Path, double[]> rtrn=new TreeMap<Path, double[]>();
	
	
		for(Path path: paths){
	
			double score=path.getCoverage();
			int size=path.getSize();
			double scanP=calculateApproximatePVal(new Double(score).intValue(), lambda, size, numMarkers, alpha);
	
			double localLambda=0;
			if(path.getLocalLambda()==0){
				Collection<Alignments> exons=path.toGene().getExonSet();
				for(Alignments exon: exons){
					if(localRates!=null && localRates.containsKey(exon)){
						localLambda=Math.max(localLambda, localRates.get(exon));
					}
				}
			}
			else{localLambda=path.getLocalLambda();}
	
			double localP=calculateApproximatePVal(new Double(score).intValue(), localLambda, size, numMarkers, alpha);
	
	
			double avgCoverage = score/size;
	
			double enrich=avgCoverage/lambda;
			double rpkm = avgCoverage * RPKMContant;
			double[] array={scanP, enrich, score, avgCoverage, rpkm, lambda, size, localP};
	
			//if(path.toGene().getNumExons()==1){System.err.println(path+"\t"+scanP+"\t"+localP+"\t"+localLambda);}
	
			//if(path.toGene().getNumExons()==1){System.err.println(score);}
	
			rtrn.put(path, array);
		}
	
		return rtrn;
	}

	public double count(Alignments align)throws IOException{
		return data.count(align, 0);//9/9/12 re-did handling of extension factors at the get tree level
	}

	private Map<Alignments, Double> countExons(Collection<Alignments> exons) throws IOException{
		Map<Alignments, Double> rtrn=new TreeMap<Alignments, Double>();
	
		int i=0;
		for(Alignments exon: exons){
			//System.err.println(i+" "+exons.size()+" "+exon);
			double val=data.countWithinExon(exon, 0);
			//int val=new Double(data.count(exon, 0)).intValue();
			rtrn.put(exon, val);
			i++;
		}
	
		return rtrn;
	}

	private Collection<RefSeqGene> collapse(Collection<RefSeqGene> genes){
		Collection<RefSeqGene> rtrn=new TreeSet<RefSeqGene>();
	
		Collection<Alignments> exons=new TreeSet<Alignments>();
		for(RefSeqGene gene: genes){
			exons.addAll(gene.getExonSet());
		}
	
		exons=CollapseByIntersection.collapseByIntersection(exons, false);
	
		for(Alignments exon: exons){rtrn.add(new RefSeqGene(exon));}
	
		return rtrn;
	}

	private Collection<Alignments> collapsePaths(Collection<Path> paths) {
		Collection<Alignments> temp=new TreeSet<Alignments>();
	
		for(Path path: paths){
			temp.add(path.toGene().getAlignment());
		}
	
		return CollapseByIntersection.collapseByIntersection(temp, false);
	}

	private Collection<RefSeqGene> collapseAllSingleExons(Collection<RefSeqGene> c){
		Collection<RefSeqGene> rtrn=new TreeSet<RefSeqGene>();
	
		Collection<Alignments> singleExon=new TreeSet<Alignments>();
	
		for(RefSeqGene gene: c){
			if(gene.getNumExons()>1){rtrn.add(gene);}
			else{singleExon.add(gene.getAlignment());}
		}
	
		singleExon=CollapseByIntersection.collapseByIntersection(singleExon, false);
	
		for(Alignments exon: singleExon){rtrn.add(new RefSeqGene(exon));}
	
		return rtrn;
	}
	
	private Alignments collapseMaxContiguous(Alignments align, double[] array)  throws IOException{
		array=subtract(array, Math.max(getLambda(align.getChr()), this.minNumberOfReadsAtEnd));
	
		//double[] array=this.getDataForAlignments(align);
		double[] maxSum=MaximumContiguousSubsequence.maxSubSum3(array);
		//System.err.println("region " + align.toUCSC() + " dist" + maxSum[0]+" "+maxSum[1]+"-"+maxSum[2]);
		Alignments newAlign=null;
		if(maxSum[0]>0){
			newAlign=new Alignments(align.getChr(),
					new Double(((align.getStart()+(maxSum[1]-1)))).intValue(),
					new Double(align.getStart()+(maxSum[2]-1)).intValue());
		}
		return newAlign;
	}

	/**
	 * Trim a region to a contiguous subregion with max number of positions having coverage above the quantile
	 * @param region The original region
	 * @param quantile The coverage quantile of the original region
	 * @return The contiguous subregion with maximal number of positions having coverage above the quantile of the original region or null if can't trim
	 * @throws IOException 
	 */
	public RefSeqGene trimMaxContiguous(RefSeqGene region, double quantile) throws IOException {
		
		int treeStartPos = Math.min(0,region.getStart());
		int treeEndPos = region.getEnd();
		IntervalTree<org.broad.igv.sam.Alignment> tree = getData().getData().getIntervalTree(region.getChr(), treeStartPos, treeEndPos);

		List<Double> coverageData = getCountsPerBp(region,tree);
		RefSeqGene rtrn = collapseMaxContiguous(region, coverageData, quantile);
		int trimmedStart = rtrn.getStart();
		int trimmedEnd = rtrn.getEnd();
		rtrn.setName(rtrn.getChr() + ":" + trimmedStart + "-" + trimmedEnd);
		rtrn.setOrientation(region.getOrientation());
		rtrn.setCDSRegion(Math.max(trimmedStart, region.getStart()), Math.min(trimmedEnd, region.getEnd()));
		
		return rtrn;
		
	}
	
	private RefSeqGene collapseMaxContiguous(RefSeqGene align, List<Double> data, double quantile) throws IOException{
		double[] array=l2a(data);
		Collections.sort(data);
		double cutoff = Math.max(minNumberOfReadsAtEnd, Statistics.quantile(data, quantile));
		array=subtract(array, Math.max(getLambda(align.getChr()), cutoff));
	
		//double[] array=this.getDataForAlignments(align);
		double[] maxSum=MaximumContiguousSubsequence.maxSubSum3(array);
	
		RefSeqGene newAlign=null;
		if(maxSum[0]>0){
			if(maxSum[1]-1<0 && maxSum[2]-1<0){}
			else{
				newAlign=align.trim(new Double(maxSum[1]-1).intValue(), new Double(maxSum[2]-1).intValue());
			}
		}
		return newAlign;
	}

	public static double calculatePVal(int k, double lambda, double w, double T){
		return AlignmentDataModelStats.calculatePVal(k, lambda, w, T);
	}

	private int calculateCriticalValue(int maxLength, double w, double T, double alpha, double lambda){
		int num=0;
		for(int k=1; k<maxLength; k++){
			num++;
			double p=calculatePVal(k, lambda, w, T); //if this gets too slow can consider precomputing
			//System.err.println(k+" "+p+" "+lambda+" "+w+" "+T);
			if(p<alpha){break;}
		}
		return num;
	}

	private static double calculateApproximatePVal(int k, double lambda, double w, double T, double alpha){
		return AlignmentDataModelStats.calculateApproximatePVal(k, lambda, w, T, alpha);
	}

	/**
	 * Segments locally, computing a local lambda using the overlapper genes
	 * lambda may be different than chromosome lambda only if an overlapping genes is multiexonic.
	 * @param genes
	 * @param tree
	 * @param chr
	 * @return
	 */
	private double computeLambda(Iterator<Node<RefSeqGene>> genes, IntervalTree<Alignment> tree, String chr) throws IOException{
	
		int numIntrons=0;
		Collection<Alignments> geneCollection=new TreeSet<Alignments>();
		while(genes.hasNext()){
			RefSeqGene gene=genes.next().getValue();
			numIntrons=Math.max(numIntrons, gene.getNumExons()-1);
			geneCollection.add(gene.getAlignment());
		}
	
		Collection<Alignments> collapsed=CollapseByIntersection.collapseByIntersection(geneCollection, false);
	
		double reads=0;
		int counts=0;
		//Go through each gene and get reads from start to end
		for(Alignments align: collapsed){
			reads+=data.getCountsPerAlignment(align, tree, 0);
			counts+=align.getSize();
		}
	
		double lambda=getLambda(chr);
		if(numIntrons>0) {
			lambda=Math.max(reads/counts, lambda);
		} 
		return lambda;
	}

	private double computeLambdaFromPath(Iterator<Node<Path>> genes, IntervalTree<Alignment> tree, Path path, IntervalTree<LightweightGenomicAnnotation> cache) throws IOException{
	
		Iterator<Node<LightweightGenomicAnnotation>> pathComputedOverlappers = cache.overlappers(path.getStart(), path.getEnd());
		if( pathComputedOverlappers.hasNext() ) {
			LightweightGenomicAnnotation containigRegion = pathComputedOverlappers.next().getValue();
			if(containigRegion.contains(new BasicLightweightAnnotation(path.getChromosome(), path.getStart(), path.getEnd()))){
				return containigRegion.getScore();
			}
		}
	
		// If the chaced tree did not contain a region that contained the current path, go ahead and compute its local lambda.
	
		//Collection<Alignments> geneCollection=new TreeSet<Alignments>();
		TreeMap<Alignments, Double> exonMap = new TreeMap<Alignments, Double>();
		TreeMap<Alignments, Double> intronMap = new TreeMap<Alignments, Double>();
	
		int numIntrons=getPathCoverageData(path, 0, exonMap, intronMap, path);
		while(genes.hasNext()){
			Path overlappingPath=genes.next().getValue();
			numIntrons = getPathCoverageData(path, numIntrons, exonMap,intronMap, overlappingPath);
			//geneCollection.addAl(gene.getAlignment());
		}
	
		Collection<Alignments> collapsedExons =CollapseByIntersection.collapseByIntersection(exonMap.keySet(), false);
	
		double reads=0;
		int counts=0;
		//Go through each gene and get reads from start to end
		for(Alignments exon: collapsedExons){
			if (exonMap.containsKey(exon)) {
				reads+= exonMap.get(exon); 
			} else {
				reads += data.getCountsPerAlignment(exon, tree, 0);
			}
			counts+=exon.getSize();
		}
	
		//Now that we have all reads within exons, lets count spliced reads.
		Set<Alignments> intronSet = intronMap.keySet();
		for(Alignments intron: intronSet){			
			reads+= intronMap.get(intron); 
		}
	
		double lambda=getLambda(path.getChromosome());
		if(numIntrons>0) {
			lambda=Math.max(reads/counts, lambda);
		} 
	
		LightweightGenomicAnnotation mergedLoci = new BasicLightweightAnnotation(path.getChromosome(), exonMap.firstKey().getStart(), exonMap.lastKey().getEnd());
		mergedLoci.setScore(lambda);
		cache.put(mergedLoci.getStart(), mergedLoci.getEnd(), mergedLoci);
		return lambda;
	}

	private  PairedEndDistribution computeInsertSizeDistribution( Collection<RefSeqGene> genes)
			throws IOException {
		BEDFileParser geneReader = new BEDFileParser(genes);
		geneReader.collapse(DEFAULT_INSERT_SIZE_FUDGE);
		PairedEndDistribution ped = new PairedEndDistribution(data.getData(), geneReader);//AWFULL!!!!
		return ped;
	}

	private  PairedEndDistribution computeInsertSizeDistribution( BEDFileParser genes)
			throws IOException {
		BEDFileParser copy = genes.copy();
		copy.collapse(DEFAULT_INSERT_SIZE_FUDGE);
		PairedEndDistribution ped = new PairedEndDistribution(data.getData(), copy);//AWFULL!!!!
		return ped;
	}

	public BEDFileParser updateReconstructionsWithPairs(BEDFileParser reconstructions, PairedEndDistribution ped) throws IOException, IllegalArgumentException, MathException {
		return updateReconstructionsWithPairs(reconstructions.GetGenes(), ped);
	}

	public BEDFileParser updateReconstructionsWithPairs(Collection<RefSeqGene> reconstructions, PairedEndDistribution ped) throws IOException, IllegalArgumentException, MathException {
		BEDFileParser rslt = new  BEDFileParser();
		BEDFileParser input = new BEDFileParser(reconstructions);
		input.makeGenes(DEFAULT_MAKE_GENE_OVERLAP);
		//Step1: Update ends and remove deviant isoforms.
		Iterator<String> chrIt = input.getChromosomeIterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			Iterator<RefSeqGeneWithIsoforms> geneIt = input.getChrTree(chr).valueIterator();
			while(geneIt.hasNext()) {
				RefSeqGeneWithIsoforms reconstructedGene = geneIt.next();				
				Collection<RefSeqGene> isoforms = reconstructedGene.getAllIsoforms(false);
				Collection<RefSeqGene> passingIsoforms = new ArrayList<RefSeqGene>();
				Collection<RefSeqGene> deviantIsoforms = new ArrayList<RefSeqGene>();
				for(RefSeqGene iso : isoforms) {
					//Iterator<Alignment> readIt = data.getIntervalTreeCached(reconstruction.getChr(), reconstruction.getStart(), reconstruction.getEnd()).valueIterator();
					PairedEndAnalysisResult pear = ped.doPairedEndCompatibilityAnalysis(iso, DEFAULT_INS_SIZE_PVAL);
					RefSeqGene fixedUpIsoform = pear.endFixedTranscript( DEFAULT_INS_SIZE_PVAL);
					fixedUpIsoform.addExtraField(String.valueOf(pear.getTranscriptInsertSizeDistribution().getMean()));
					fixedUpIsoform.addExtraField(String.valueOf(pear.getTranscriptInsertSizeDistribution().getStandardDeviation()));
					fixedUpIsoform.addExtraField( String.valueOf( pear.isDeviant() ? 1 :0));
					if(pear.isDeviant()) {
						passingIsoforms.add(fixedUpIsoform);
					} else {
						deviantIsoforms.add(fixedUpIsoform);
					}
				}
				Collection<RefSeqGene> fixedUpIsoformsToUse = passingIsoforms.size() > 0 ? passingIsoforms : deviantIsoforms;
				rslt.addRefSeqSet(fixedUpIsoformsToUse);
			}
		}
		return rslt;
	}

	private Collection<RefSeqGene> acrossGaps(IntervalTree<Alignment> chunkAlignmentTree, int critVal, int fixedWidth, GenomeWithGaps2 gwg){
		long startTime=System.currentTimeMillis();
		long endTime=0;
	
		Collection<RefSeqGene> rtrn=new TreeSet<RefSeqGene>();
	
		int counter=0;		
		for(int i=0; i<gwg.getRelativeGenomeLength(); i++){
			RefSeqGene window=gwg.getRelativeWindow(i, i+fixedWidth);
			double sum=data.getCountsPerAlignment(window, chunkAlignmentTree, 0); //9/9/12 re-did handling of extension factors at the get tree level
			if(sum>critVal){rtrn.add(window);}
	
			/*******Iteration controller***************/
			if(sum==0){i=i+fixedWidth;}
			/******************************************/
	
			counter++;
			if(counter % 10000 ==0){endTime=System.currentTimeMillis(); logger.debug(window.getAlignment().toUCSC()+" "+window.getExons().length+" "+counter+" "+gwg.getRelativeGenomeLength()+" "+(counter/(double)gwg.getRelativeGenomeLength())+" "+(endTime-startTime)); startTime=System.currentTimeMillis();}
			//TODO Collapse on the fly to avoid huge memory footprint
		}
		return rtrn;
	}

	//TODO Works great just doesnt work efficiently due to Manuel's old implementation
	private Collection<RefSeqGene>[] acrossGraph(int critVal, int fixedWidth, ChromosomeWithBubblesJGraphT graph, Collection<Alignments> regions) throws IOException{
	
		int counterWindows=0;
		long startTime=System.currentTimeMillis();
		long endTime=0;
	
		long pathTime=0;
	
		//TODO This might also be inefficient instead collapse on the fly
		//Map<RefSeqGene, Boolean> rtrn=new TreeMap();
		Collection<RefSeqGene> falseSet=new TreeSet<RefSeqGene>();
		Collection<RefSeqGene> trueSet=new TreeSet<RefSeqGene>();
		int counter=0;
	
		for(Alignments region: regions){
			IntervalTree<Alignment> chunkAlignmentTree=data.getIntervalTree(region.getChr(), region.getStart()-100000, region.getEnd()+100000);
			//System.err.println("Scanning "+region.toUCSC());
			for(int i=region.getStart()-100000; i<region.getEnd()+100000; i++){
				double maxScore=0;
				Collection<RefSeqGene> windows=graph.getPaths(i, i+fixedWidth);
				for(RefSeqGene window: windows){
					Alignments firstExon=window.getFirstExon();
					boolean isStartingInGapButEndingOutside=graph.isSpanningGap(firstExon.getStart(), firstExon.getEnd());
					//if(window.getExons().length > 1)
					//System.out.println(window.toBED());
					//TODO This caching will be wildly inefficient because we're double storing every value
					if(!isStartingInGapButEndingOutside){
						//IntervalTree<Collection<Alignment>> chunkAlignmentTree=data.getIntervalTreeCached(window.getChr(), window.getStart(), window.getEnd());
						long startMem = Runtime.getRuntime().freeMemory();
						long t = System.nanoTime();
						double sum=data.getScorePerAlignmentFromCache(window, chunkAlignmentTree, 0); //9/9/12 re-did handling of extension factors at the get tree level
						//System.err.println("Time " + (System.nanoTime() - t) + " free mem change " + (Runtime.getRuntime().freeMemory()- startMem));
						if(sum>critVal){
							boolean b=graph.isGap(window.getStart(), window.getEnd());
							if(b){trueSet.add(window);}
							else{falseSet.add(window);}
						} // its a boolean telling me whether the significant fragment came from within an intron
						maxScore=Math.max(maxScore, sum);
					}
				}
	
				//if(maxScore==0){i=i+fixedWidth;}
				counter++;
				if(counter % 10000 ==0){endTime=System.currentTimeMillis(); logger.debug(((double)i/graph.getRelativeGenomeLength())+" "+new Alignments(graph.getName(), i, i+fixedWidth).toUCSC()+" "+(endTime-startTime)+" "+counterWindows+" "+pathTime); counterWindows=0; startTime=System.currentTimeMillis(); pathTime=0;}
			}
		}
	
		Collection<RefSeqGene>[] array=new Collection[2];
		array[0]=trueSet;
		array[1]=falseSet;
		return array;
	}

	private ChromosomeWithBubblesJGraphT acrossGraph(int critVal, int fixedWidth, ChromosomeWithBubblesJGraphT graph, String chr, Sequence chrSeq, double minimumSpliceFrequencyToFollowEdge)throws IOException{
		long startTime=System.currentTimeMillis();
		long endTime=0;
	
		long pathTime=0;
		long treeTime=0;
		long sumTime=0;
		long checkTime=0;
		long isGapTime=0;
		long totalLength=0;
	
		IntervalTree<Alignments> exonTreeStar=new IntervalTree<Alignments>();
		IntervalTree<Alignments> exonTreePlus=new IntervalTree<Alignments>();
		IntervalTree<Alignments> exonTreeMinus=new IntervalTree<Alignments>();
		IntervalTree<Alignments> exonTreeSinglets=new IntervalTree<Alignments>();
	
		Collection<Alignments> intronsStar=new TreeSet<Alignments>();
		Collection<Alignments> intronsPlus=new TreeSet<Alignments>();
		Collection<Alignments> intronsMinus=new TreeSet<Alignments>();
	
		int counter=0;
	
		WindowIterator it = graph.iterator(fixedWidth,0);
	
		while(it.hasNext()) {
			int minPosition=Integer.MAX_VALUE;
			int maxPosition=0;
			long pathStart=System.currentTimeMillis();
			Collection<RefSeqGene> windows = it.next();
			pathTime+=(System.currentTimeMillis()-pathStart);
	
			for(RefSeqGene window: windows){
				long treeStart=System.currentTimeMillis();
				IntervalTree<Alignment> chunkAlignmentTree=data.getIntervalTreeCached(window.getChr(), window.getStart(), window.getEnd());
				treeTime+=(System.currentTimeMillis()-treeStart);
				long sumStart=System.currentTimeMillis();
				boolean passes=data.passes(window, chunkAlignmentTree, 0, critVal);
				sumTime+=(System.currentTimeMillis()-sumStart);
				totalLength+=window.getAlignment().getSize();
				long checkStart=System.currentTimeMillis();
	
				if(chunkAlignmentTree.isEmpty()){
					it.jumpTo(data.getChunkStart(), data.getChunkEnd());
					//System.err.println("Jumping "+data.getChunkStart()+" "+data.getChunkEnd()+" "+chunkAlignmentTree.size());
				}
				else if(passes){
					if(chrSeq != null) {
						window.setOrientation(GeneTools.orientationForGene(window, chrSeq));
					}
					if(window.getNumExons()==1){
						exonTreeSinglets=mergeAndRemove(exonTreeSinglets, window.getExonSet());
					}
					else if(window.getOrientation().equalsIgnoreCase("+")){
						intronsPlus.addAll(window.getIntronSet());
						exonTreePlus=mergeAndRemove(exonTreePlus, window.getExonSet());
					}
					else if(window.getOrientation().equalsIgnoreCase("-")){
						intronsMinus.addAll(window.getIntronSet());
						exonTreeMinus=mergeAndRemove(exonTreeMinus, window.getExonSet());
					}
					else{
						intronsStar.addAll(window.getIntronSet());
						exonTreeStar=mergeAndRemove(exonTreeStar, window.getExonSet());
					}	
				} 
				checkTime+=(System.currentTimeMillis()-checkStart);
				minPosition=Math.min(minPosition, window.getEnd());
				maxPosition=Math.max(maxPosition, window.getStart());
			}
			counter++;
			if(counter % 100000 ==0){
				endTime=System.currentTimeMillis(); logger.info(((double)minPosition/graph.getRelativeGenomeLength())+" free memory "+ Runtime.getRuntime().freeMemory() + " "+(endTime-startTime)+" Path Time: "+pathTime+" Tree Time: "+treeTime+" Sum Time: "+sumTime+" Check Time: "+checkTime+" Gap Time: "+isGapTime); startTime=System.currentTimeMillis(); pathTime=0; sumTime=0; treeTime=0; checkTime=0; isGapTime=0;
			}
		}
	
	
	
		Collection<Alignments> exons=new TreeSet<Alignments>();
		Collection<Alignments> introns=new TreeSet<Alignments>();
	
		Collection<Alignments> exonTreeCollectionPlus=exonTreePlus.toCollection();
		exonTreeCollectionPlus.addAll(exonTreeSinglets.toCollection());
		exons.addAll(exonTreeCollectionPlus);
		introns.addAll(intronsPlus);
	
		Collection<Alignments> exonTreeCollectionMinus=exonTreeMinus.toCollection();
		exonTreeCollectionMinus.addAll(exonTreeSinglets.toCollection());
		exons.addAll(exonTreeCollectionMinus);
		introns.addAll(intronsMinus);
	
		Collection<Alignments> exonTreeCollectionStar=exonTreeStar.toCollection();
		exonTreeCollectionStar.addAll(exonTreeSinglets.toCollection());
		exons.addAll(exonTreeCollectionStar);
		introns.addAll(intronsStar);
	
	
		data.resetTreeCache();
	
		ChromosomeWithBubblesJGraphT bubbles=this.makeGraphs(chr, exons, introns, null, minimumSpliceFrequencyToFollowEdge);
	
		return bubbles;
	}

	private Collection<RefSeqGene> mergeIntoTranscripts(Collection<RefSeqGene> genes, String chr) throws IOException{
	
		Collection<Alignments> exons=new TreeSet<Alignments>();
		Collection<Alignments> introns=new TreeSet<Alignments>();
	
		for(RefSeqGene gene: genes){
			exons.addAll(gene.getExonSet());
			introns.addAll(gene.getIntronSet());
		}		
	
		long start=System.currentTimeMillis();
		exons=CollapseByIntersection.collapseByIntersection(exons, false);
	
		long end=System.currentTimeMillis();
		logger.debug("Collapse: "+(end-start));
	
		//for(Alignments exon: exons){System.out.println(exon);}
	
		start=System.currentTimeMillis();
		if(!introns.isEmpty()){
			exons=CollapseByIntersection.DecollapseByIntronLocation(exons, introns);
		}
		end=System.currentTimeMillis();
		logger.info("Decollapse: "+(end-start));
	
	
	
		//writeTest(exons, introns);
	
		logger.debug(exons.size()+" "+introns.size());
		ChromosomeWithBubblesJGraphT bubbles=new ChromosomeWithBubblesJGraphT(chr, exons, introns, data.getLambda(chr), data.getNumberMarkers(chr), data.getNumberOfReads(chr), minNumberOfSplices);
	
		Collection<RefSeqGene> rtrn=new TreeSet<RefSeqGene>();
		rtrn.addAll(bubbles.getGenePaths(0));
		rtrn.addAll(bubbles.getOrphanNodes());
	
		return rtrn;
	}

	// remove intervals in iter
	// replace with one interval consisting of the union of the span of all intervals in iter and align
	private IntervalTree<Alignments> mergeAndRemove(Iterator<Node<Alignments>> iter, Alignments align, IntervalTree<Alignments> tree){
		IntervalTree<Alignments> rtrn=tree;
	
		int start=align.getStart();
		int end=align.getEnd();
	
		while(iter.hasNext()){
			Alignments next=iter.next().getValue();
			start=Math.min(next.getStart(), start);
			end=Math.max(next.getEnd(), end);
			rtrn.remove(next.getStart(), next.getEnd());
		}
	
		Alignments newAlign=new Alignments(align.getChr(), start, end);
	
		rtrn.put(start, end, newAlign);
	
		return rtrn;
	}

	private IntervalTree<Alignments> mergeAndRemove(IntervalTree<Alignments> tree, Collection<Alignments> exons){
		IntervalTree<Alignments> rtrn=tree;
	
		for(Alignments exon: exons){
			//System.out.println(exon);
			int start=exon.getStart();
			int end=exon.getEnd();
			Iterator<Node<Alignments>> overlappers=rtrn.overlappers(exon.getStart(), exon.getEnd());
			while(overlappers.hasNext()){
				Alignments next=overlappers.next().getValue();
				start=Math.min(next.getStart(), start);
				end=Math.max(next.getEnd(), end);
				rtrn.remove(next.getStart(), next.getEnd());
			}
			Alignments newAlign=new Alignments(exon.getChr(), start, end);
			rtrn.put(start, end, newAlign);
		}
		return rtrn;
	}

	//TODO: Might need to consider an iterative addition
	private Collection<Alignments> makeAdditionalNodes(Collection<Alignments> significantPieces, Map<String, IntervalTree<Alignments>> goodExons) {
		Collection<Alignments> rtrn=new TreeSet<Alignments>();
	
		for(Alignments region: significantPieces){
			Iterator<Node<Alignments>> right=goodExons.get(region.getChr()).overlappers(region.getStart(), region.getEnd()+1);
			while(right.hasNext()){
				Alignments addition=right.next().getValue();
				Alignments align=new Alignments(region.getChr(), Math.min(region.getStart(), addition.getStart()), Math.max(region.getEnd(), addition.getEnd()));
				rtrn.add(align);
			}

			Iterator<Node<Alignments>> left=goodExons.get(region.getChr()).overlappers(region.getStart()-1, region.getEnd());
			while(left.hasNext()){
				Alignments addition=left.next().getValue();
				Alignments align=new Alignments(region.getChr(), Math.min(region.getStart(), addition.getStart()), Math.max(region.getEnd(), addition.getEnd()));
				rtrn.add(align);
			} 

			Iterator<Node<Alignments>> full=goodExons.get(region.getChr()).overlappers(region.getStart()-1, region.getEnd());
			//TODO Should add both to the same exon
			while(full.hasNext()){
				Alignments addition=full.next().getValue();
				Alignments align=new Alignments(region.getChr(), Math.min(region.getStart(), addition.getStart()), Math.max(region.getEnd(), addition.getEnd()));
				rtrn.add(align);
			}
		}
	
		return rtrn;
	}

	private ChromosomeWithBubblesJGraphT makeGraphWithCounts(String chr, int minDistance, Sequence chrSeq,  double minimumSpliceFrequencyToFollowEdge) throws IOException {
		logger.debug("Going to get read iterator to make graph with counts");
		CloseableIterator<Alignment> iter=data.getReadIterator(new Alignments(chr, 0, data.getChromosomeLengths().get(chr)));
		logger.debug("Got read iterator");
		IntervalTree<Alignment> intronTree=new IntervalTree<Alignment>();
		double numberOfReads=0;
		double splicedReads=0;
	
		Alignments current=null;
		int counter=0;
		Map<Alignments, Integer> exons=new TreeMap<Alignments, Integer>();
		Map<Alignments, Integer> splicedExons=new TreeMap<Alignments, Integer>();
		Map<Alignments, Double> intronMap=new TreeMap<Alignments, Double>();
	
		while(iter.hasNext()){
			Alignment read=iter.next();
			//System.err.print("read " + read.getChr() + ":" + read.getAlignmentStart()+"-"+read.getAlignmentEnd() + " qual " + read.getMappingQuality() + " cigar " + read.getCigarString());
			/*
			 *  IF READ PASSES MAPPING QUALITY TEST
			 */
			ReadMate mate = read.getMate();
			if(read.getMappingQuality() > data.getMinimumMappingQuality() && (!isPaired() || (mate != null && mate.isMapped()) )){ //This is only to use paired reads .. .not sure if this is right.
				//System.err.print(" PASSES QUAL ");
				
				/*
				 * IF READ IS SPLICED
				 */
				
				Collection<Alignments> alignmentBlocks = getAlignmentBlocks(read);
				if( alignmentBlocks.size() > 1  ){
					//System.err.print(" IS SPLICED ");
					/*
					 * ADD READ TO INTRON TREE
					 */
					Collection<Alignments> counts = toExons(read);
					Collection<Alignments> introns = getIntronsFromExons(counts);
					
					intronTree.put(read.getAlignmentStart(), read.getAlignmentEnd(), read);

					for(Alignments exon: counts){
						int count=0;
						if(splicedExons.containsKey(exon)){count=splicedExons.get(exon);}
						count++;
						splicedExons.put(exon, count);
					}

					for(Alignments intron: introns){
						try{
							if(chrSeq != null) {
								intron.setOrientation(GeneTools.orientationFromSpliceSites(intron, chrSeq));
							}
							double count= intronMap.containsKey(intron) ? intronMap.get(intron) : 0;
							count++;
							intronMap.put(intron, count);
						}catch(Exception ex){ex.printStackTrace(); logger.error(intron);}
					}
					splicedReads++;
				}
				else{
					//this is where we want to aggregate
					//System.err.print(" WAS UNSPLICED  ");
					if(current==null){
						current=new Alignments(read.getChromosome(),read.getAlignmentStart(), read.getAlignmentEnd()); 
						counter=1;
						//System.err.println(" current was null ");
					}
					else if(overlaps(read, current)){
						current=new Alignments(current.getChr(), Math.min(current.getStart(), read.getAlignmentStart()), Math.max(current.getEnd(), read.getAlignmentEnd()));
						counter++;
						//System.err.println(" current incremented " + counter);
					}
					else{
						//add current to collection
						exons.put(new Alignments(current), counter);
						//set current to null
						current=new Alignments(read.getChromosome(),read.getAlignmentStart(), read.getAlignmentEnd());
						counter=1;
						//System.err.println(" not overlapped, added to map ");
					}
					/*if (!iter.hasNext()) {
						//System.err.println(" last read, adding current exon ");
						exons.put(current, counter);
					}*/
				}
				numberOfReads++;
			}
		}
		if(current != null) {
			exons.put(current, counter);
		}
		iter.close();
	
		//Map<Alignments, Double> filteredIntronMap = new TreeMap<Alignments, Double>();
		/*
		 *The filtering that follows is redundant, the chromosome with bubbles constructor will already ignore edges with less than the specified minNumberOfSplices
		 *TODO: Remove this filetering and make sure there is no problem 
		 */
		logger.debug("minimum splices to keep intron: " + minNumberOfSplices);
		Set<Alignments> keys =  new TreeSet<Alignments>(intronMap.keySet());
		for(Alignments a : keys) {
			a.setScore(intronMap.get(a));
			//logger.trace("intron " + a + " score is? " + a.getScore() + " score from map? " + intronMap.get(a));
			if(intronMap.get(a) <= minNumberOfSplices) {
				intronMap.remove(a);
				logger.trace("removed candidate intron, not enough support " + a.toUCSC());
			}
		}
	
	
	
		double numMarkers=data.chromosomeLength(chr);
	
		if(maskedRegions!=null && maskedRegions.containsKey(chr)){
			numMarkers=numMarkers - maskedRegions.get(chr); 
		}
		data.setNumberOfReads(chr, numberOfReads);
		data.setLambda(chr, numberOfReads);
		data.setNumberOfMarkers(chr, numMarkers);
	
		logger.info("Made it through all reads. Total reads " + numberOfReads );
	
		//TODO Merge spliced and unspliced while retaining proper counting
		Collection<Alignments> all=new TreeSet<Alignments>();
		all.addAll(exons.keySet());
		all.addAll(splicedExons.keySet());
		BEDFileParser.writeBED(chr+"_exons.bed", all);
		BEDFileParser.writeBED(chr+"introns.bed", intronMap.keySet());
		//writeVal("intronCounts.bed", intronMap);
		logger.info("All size " + all.size());
		all=CollapseByIntersection.collapseByIntersection(all, false);
		logger.info("Collapsed reads");
		BEDFileParser.writeBED(chr+"_collapsed.bed", all);
		Collection<Alignments> decollapsed=CollapseByIntersection.DecollapseByIntronLocation(all, intronMap.keySet());
		logger.info("Decollapsed by introns");
		BEDFileParser.writeBED(chr+"_decollapsed.bed", decollapsed);
	
		ChromosomeWithBubblesJGraphT graph=new ChromosomeWithBubblesJGraphT(chr, decollapsed, intronMap.keySet(), intronMap, null, data.getLambda(chr), data.getNumberMarkers(chr), data.getNumberOfReads(chr), minNumberOfSplices);
		
		logger.info("Made first graph");
	
		//writeFullBED("firstGraph.bed", graph.getGenePaths(1));
	
		//Second, get exons that were fully within introns
		//Also split up covered exons into pieces and compute significance
		Collection<Alignments> orphans=getExtendedExons(all, decollapsed, graph, minDistance, minimumSpliceFrequencyToFollowEdge);
		//System.err.println("Got extended pieces");
		//BEDFileParser.writeBED("extendedpieces.bed", orphans);
		//Add new nodes to graph
		decollapsed.addAll(orphans);
		//BEDFileParser.writeBED("decollapsed.with.orphas.bed", decollapsed);
		Map<Alignments, Double> nodeCounts=getNodeCounts(decollapsed, exons);
	
		//long start=System.currentTimeMillis();
		//Map<Alignments, Double> nodeCounts=this.countExons(decollapsed);
		//writeVal("exonCounts2.bed", nodeCounts);
		//long end=System.currentTimeMillis();
		//System.err.println("Node scoring took: "+(end-start)/1000.0);
		//start=System.currentTimeMillis();
		graph=new ChromosomeWithBubblesJGraphT(chr, decollapsed, intronMap.keySet(), intronMap, nodeCounts, data.getLambda(chr), data.getNumberMarkers(chr), data.getNumberOfReads(chr), 0); //DO NOT LOOK AT SPLICES NOW
		graph.removeSpanningGapNodes();
		logger.info("Made second graph");
		//writeFullBED("secondGraph.bed", graph.getGenePaths(1, minimumSpliceFrequencyToFollowEdge));
		//end=System.currentTimeMillis();
	
		//System.err.println("Making graph took: "+(end-start)/1000.0);
	
		//writeFullBED("secondGraph.bed", graph.getGenePaths(1));
	
		data.resetTreeCache(chr);
	
		double spliceWeight= upWeightSplices ? numberOfReads/splicedReads : 1;
	
		graph.setSpliceWeight(spliceWeight);
		
		return graph;
	}

	private ChromosomeWithBubblesJGraphT makeGraphWithCounts(String chr, int minDistance, Sequence chrSeq,  double minimumSpliceFrequencyToFollowEdge,String strand) throws IOException {
		logger.debug("Going to get read iterator to make graph with counts");
		CloseableIterator<Alignment> iter=data.getReadIterator(new Alignments(chr, 0, data.getChromosomeLengths().get(chr)));
		logger.debug("Got read iterator");
		IntervalTree<Alignment> intronTree=new IntervalTree<Alignment>();
		double numberOfReads=0;
		double splicedReads=0;
	
		Alignments current=null;
		int counter=0;
		Map<Alignments, Integer> exons=new TreeMap<Alignments, Integer>();
		Map<Alignments, Integer> splicedExons=new TreeMap<Alignments, Integer>();
		Map<Alignments, Double> intronMap=new TreeMap<Alignments, Double>();
	
		while(iter.hasNext()){
			Alignment read=iter.next();
			//System.err.print("read " + read.getChr() + ":" + read.getAlignmentStart()+"-"+read.getAlignmentEnd() + " qual " + read.getMappingQuality() + " cigar " + read.getCigarString());
			/*
			 *  IF READ PASSES MAPPING QUALITY TEST
			 */
			ReadMate mate = read.getMate();
			if(read.getMappingQuality() > data.getMinimumMappingQuality() && (!isPaired() || (mate != null && mate.isMapped()) ) && (read.isNegativeStrand() == "-".equals(strand))){ //This is only to use paired reads .. .not sure if this is right.
				//System.err.print(" PASSES QUAL ");
				/*
				 * IF READ IS SPLICED
				 */			
				Collection<Alignments> alignmentBlocks = getAlignmentBlocks(read);
				if( alignmentBlocks.size() > 1  ){
					//System.err.print(" IS SPLICED ");
					/*
					 * ADD READ TO INTRON TREE
					 */
					Collection<Alignments> counts = toExons(read);
					Collection<Alignments> introns = getIntronsFromExons(counts);
					
					intronTree.put(read.getAlignmentStart(), read.getAlignmentEnd(), read);

					for(Alignments exon: counts){
						int count=0;
						if(splicedExons.containsKey(exon)){count=splicedExons.get(exon);}
						count++;
						splicedExons.put(exon, count);
					}

					for(Alignments intron: introns){
						try{
							if(chrSeq != null) {
								intron.setOrientation(GeneTools.orientationFromSpliceSites(intron, chrSeq));
							}
							double count= intronMap.containsKey(intron) ? intronMap.get(intron) : 0;
							count++;
							intronMap.put(intron, count);
						}catch(Exception ex){ex.printStackTrace(); logger.error(intron);}
					}
					splicedReads++;
				}
				else{
					//this is where we want to aggregate
					//System.err.print(" WAS UNSPLICED  ");
					if(current==null){
						current=new Alignments(read.getChromosome(),read.getAlignmentStart(), read.getAlignmentEnd()); 
						counter=1;
						//System.err.println(" current was null ");
					}
					else if(overlaps(read, current)){
						current=new Alignments(current.getChr(), Math.min(current.getStart(), read.getAlignmentStart()), Math.max(current.getEnd(), read.getAlignmentEnd()));
						counter++;
						//System.err.println(" current incremented " + counter);
					}
					else{
						//add current to collection
						exons.put(new Alignments(current), counter);
						//set current to null
						current=new Alignments(read.getChromosome(),read.getAlignmentStart(), read.getAlignmentEnd());
						counter=1;
						//System.err.println(" not overlapped, added to map ");
					}
					/*if (!iter.hasNext()) {
						//System.err.println(" last read, adding current exon ");
						exons.put(current, counter);
					}*/
				}
				numberOfReads++;
			}
		}
		if(current != null) {
			exons.put(current, counter);
		}
		iter.close();
	
		//Map<Alignments, Double> filteredIntronMap = new TreeMap<Alignments, Double>();
		/*
		 *The filtering that follows is redundant, the chromosome with bubbles constructor will already ignore edges with less than the specified minNumberOfSplices
		 *TODO: Remove this filetering and make sure there is no problem 
		 */
		logger.debug("minimum splices to keep intron: " + minNumberOfSplices);
		Set<Alignments> keys =  new TreeSet<Alignments>(intronMap.keySet());
		for(Alignments a : keys) {
			a.setScore(intronMap.get(a));
			logger.trace("intron " + a + " score is? " + a.getScore() + " score from map? " + intronMap.get(a));
			if(intronMap.get(a) <= minNumberOfSplices) {
				intronMap.remove(a);
				logger.trace("removed candidate intron, not enough support " + a.toUCSC());
			}
		}
	
		double numMarkers=data.chromosomeLength(chr);
	
		if(maskedRegions!=null && maskedRegions.containsKey(chr)){
			numMarkers=numMarkers - maskedRegions.get(chr); 
		}
		
		data.setNumberOfReads(chr, numberOfReads);
		data.setLambda(chr, numberOfReads);
		data.setNumberOfMarkers(chr, numMarkers);
	
		logger.info("Made it through all reads. Total reads " + numberOfReads );
	
		//TODO Merge spliced and unspliced while retaining proper counting
		Collection<Alignments> all=new TreeSet<Alignments>();
		all.addAll(exons.keySet());
		all.addAll(splicedExons.keySet());
		if(strand.equals("+")){
			BEDFileParser.writeBED(chr+"_plus_exons.bed", all);
			BEDFileParser.writeBED(chr+"_plus_introns.bed", intronMap.keySet());
		}
		else{
			BEDFileParser.writeBED(chr+"_minus_exons.bed", all);
			BEDFileParser.writeBED(chr+"_minus_introns.bed", intronMap.keySet());
		}
		//writeVal("intronCounts.bed", intronMap);
		logger.info("All size " + all.size());
		all=CollapseByIntersection.collapseByIntersection(all, false);
		logger.info("Collapsed reads");
		BEDFileParser.writeBED(chr+"_collapsed.bed", all);
		Collection<Alignments> decollapsed=CollapseByIntersection.DecollapseByIntronLocation(all, intronMap.keySet());
		logger.info("Decollapsed by introns");
		BEDFileParser.writeBED(chr+"_decollapsed.bed", decollapsed);
	
		ChromosomeWithBubblesJGraphT graph=new ChromosomeWithBubblesJGraphT(chr, decollapsed, intronMap.keySet(), intronMap, null, data.getLambda(chr), data.getNumberMarkers(chr), data.getNumberOfReads(chr), minNumberOfSplices);
		
		logger.info("Made first graph");
	
		//writeFullBED("firstGraph.bed", graph.getGenePaths(1));
	
		//Second, get exons that were fully within introns
		//Also split up covered exons into pieces and compute significance
		Collection<Alignments> orphans=getExtendedExons(all, decollapsed, graph, minDistance, minimumSpliceFrequencyToFollowEdge);
		//System.err.println("Got extended pieces");
		//BEDFileParser.writeBED("extendedpieces.bed", orphans);
		//Add new nodes to graph
		decollapsed.addAll(orphans);
		//BEDFileParser.writeBED("decollapsed.with.orphas.bed", decollapsed);
		Map<Alignments, Double> nodeCounts=getNodeCounts(decollapsed, exons);
	
		//long start=System.currentTimeMillis();
		//Map<Alignments, Double> nodeCounts=this.countExons(decollapsed);
		//writeVal("exonCounts2.bed", nodeCounts);
		//long end=System.currentTimeMillis();
		//System.err.println("Node scoring took: "+(end-start)/1000.0);
	
	
		//start=System.currentTimeMillis();
		graph=new ChromosomeWithBubblesJGraphT(chr, decollapsed, intronMap.keySet(), intronMap, nodeCounts, data.getLambda(chr), data.getNumberMarkers(chr), data.getNumberOfReads(chr), 0); //DO NOT LOOK AT SPLICES NOW
		graph.removeSpanningGapNodes();
		logger.info("Made second graph");
		//writeFullBED("secondGraph.bed", graph.getGenePaths(1, minimumSpliceFrequencyToFollowEdge));
		//end=System.currentTimeMillis();
	
		//System.err.println("Making graph took: "+(end-start)/1000.0);
	
		//writeFullBED("secondGraph.bed", graph.getGenePaths(1));
	
		data.resetTreeCache(chr);
	
		double spliceWeight= upWeightSplices ? numberOfReads/splicedReads : 1;
	
		graph.setSpliceWeight(spliceWeight);
	
	
		return graph;
	}

	
	
	private List<Alignments> getAlignmentBlocks(Alignment read) {
		Stack<Alignments> blocks = new Stack<Alignments>();
		AlignmentBlock[] originalBlocks = read.getAlignmentBlocks();
		
		for(AlignmentBlock ab : originalBlocks) {
			if(blocks.isEmpty()) {
				blocks.push(new Alignments(read.getChr(), ab.getStart(), ab.getEnd()));
			} else {
				Alignments lastBlock = blocks.pop();
				if(ab.getStart() - lastBlock.getEnd() < minIntronSize) {
					lastBlock.setEnd(ab.getEnd());
					blocks.push(lastBlock);
				} else {
					blocks.push(lastBlock);
					blocks.push(new Alignments(read.getChr(), ab.getStart(), ab.getEnd()));
				}
			}
		}
		
		
		return blocks;
	}

	private ChromosomeWithBubblesJGraphT makeGraphs(String chr, Collection<Alignments> exons, Collection<Alignments> introns, Map<Alignments, Double> intronCounts, double minimumSpliceFrequencyToFollowEdge) throws IOException {
		ChromosomeWithBubblesJGraphT graph=constructGraphs(chr, exons, introns,0, minimumSpliceFrequencyToFollowEdge);
	
		logger.info("Finished 2nd Graph construction");
	
		Collection<RefSeqGene> allPaths=graph.getGenePaths(1);
	
		allPaths.addAll(graph.getOrphanNodes());
	
		//Rebuild the graph
		exons=getExons(allPaths);
		introns=getIntrons(allPaths);
	
		Map<Alignments, Double> exonCounts=countExons(exons);
	
		graph=new ChromosomeWithBubblesJGraphT(graph.getName(), exons, introns, intronCounts, exonCounts, data.getLambda(chr), data.getNumberMarkers(chr), data.getNumberOfReads(chr), 0);
	
		data.resetTreeCache();
	
		return graph;
	}

	private Collection<Alignments> trimEnds(Collection<Alignments> windows) throws IOException{
		Set<Alignments> rtrn = new TreeSet<Alignments>();
	
		//loop through each alignment and compute number of reads for the ends
		Map<Alignments, double[]> counts=this.getDataForAlignments(windows);
		for(Alignments align: windows){
			Alignments trunc=trimEnds(align, counts.get(align));
			if(trunc!=null){rtrn.add(trunc);}
		}
	
		return rtrn;
	}

	private Collection<RefSeqGene> trimEndsForGenes(Collection<RefSeqGene> allSegments, double quantile) throws IOException{
		Set<RefSeqGene> rtrn=new TreeSet<RefSeqGene>();
	
		long start=System.currentTimeMillis();
		//loop through each alignment and compute number of reads for the ends
		logger.debug("Getting reads for annotations");
		//System.err.println("Getting reads for annotations");
		Map<RefSeqGene, List<Double>> counts=getDataForGene(allSegments);
		long end=System.currentTimeMillis();
		logger.debug("Counting data: "+(end-start));
		//System.err.println("Counting data: "+(end-start)+"ms.");
		
		for(RefSeqGene align: allSegments){
			List<Double> alignData =  counts.get(align);
			if(trimEnds){
			//	start=System.currentTimeMillis();
				RefSeqGene trunc = trimEnds(align, alignData, quantile);
	
				if(trunc!=null){rtrn.add(trunc);}
			//	end=System.currentTimeMillis();
			//	System.err.println("Trimming ends for Gene: "+align.getName()+" in "+(end-start)+"ms.");
			}
			else{
			//	start=System.currentTimeMillis();
				RefSeqGene trunc=collapseMaxContiguous(align, alignData, quantile);
				if(trunc != null) {
					rtrn.add(align);
				}
			//	end=System.currentTimeMillis();
			//	System.err.println("Collapsing max contiguous for Gene: "+align.getName()+" in "+(end-start)+"ms.");
			}
		}
	
		return rtrn;
	
	}

	private BEDFileParser trimEndsForGenes(BEDFileParser genes, double quantile) throws IOException{
		BEDFileParser rtrn=new BEDFileParser();
		logger.debug("Enter trimming on quantile");
		long start=System.currentTimeMillis();
		//loop through each alignment and compute number of reads for the ends
		logger.debug("Trimming. TrimEnds ? " + trimEnds +" Getting reads for annotations");
		List<RefSeqGene> genesAsList = genes.GetGenes();
		Map<RefSeqGene, List<Double>> counts=getDataForGene(genesAsList);
		long end=System.currentTimeMillis();
		logger.debug("Counting data: "+(end-start));
	
		for(RefSeqGene align: genesAsList){
			List<Double> alignData =  counts.get(align);
			if(trimEnds){
				RefSeqGene trunc = trimEnds(align, alignData, quantile);
	
				if(trunc!=null){rtrn.addRefSeq(trunc);}
			}
			else{
				RefSeqGene trunc=collapseMaxContiguous(align, alignData, quantile);
				if(trunc != null) {
					rtrn.addRefSeq(trunc);
				}
			}
		}
	
		return rtrn;
	
	}

public BEDFileParser reconstructGeneEndsUsingDivergence(BEDFileParser genes,boolean trimEnds,double quantile,String filename) throws IOException, MathException{
		
		//If trimEnds is set, trim ends using trimQuantile
/*		if(trimEnds) {
			logger.info("Trimming reconstructions ");
			genes = trimEndsForGenes(genes, trimQuantile);
			logger.info("Rescoring trimmed genes");
			scoreGenes(genes);
		}*/
	
		//BufferedWriter bw = new BufferedWriter(new FileWriter(filename+".kl.bin20.empirical.distributions"));
		BEDFileParser rtrn=new BEDFileParser();
	
		long starttime=System.currentTimeMillis();
		// Get all isoforms in the bed file and data for each gene in the annotations file.
		System.out.println("Trimming. TrimEnds ? " + trimEnds +" Getting reads for annotations");
		List<RefSeqGene> genesList = genes.GetGenes();
	
		//Get counts for the last exons of genes
		//Map<RefSeqGene, double[]> counts=new HashMap<RefSeqGene,double[]>();

		/*
		 * EDIT GENE LIST 
		 * STEP 1 : REMOVE ALL GENES WITH LESS THAN 5 reads and less than 25th quantile expression amongst the expressed genes
		 */
		Map<RefSeqGene, List<Double>> countsMap=getDataForGene(genesList);
		List<RefSeqGene> filteredGenes = filterGenesByExpression(countsMap);
		
		/*
		 * OUTPUT THE GENES
		 */
		
		long endtime=System.currentTimeMillis();
		System.out.println("Counting data: "+(endtime-starttime));
		
		//Result Map
		Map<RefSeqGene, List<Double>> resultMap=new HashMap<RefSeqGene,List<Double>>();
		//Window size must be less than the gene length
		int window = 300;
		int step = 50;
		
		List<BED> bedsklDivergence = new ArrayList<BED>();
		List<BED> bedsChiDistance = new ArrayList<BED>();
		List<BED> bedsChiPvalues = new ArrayList<BED>();
		
		BufferedWriter pvalbw = new BufferedWriter(new FileWriter("max.kls.pvalues"));
		BufferedWriter kldsbw = new BufferedWriter(new FileWriter("klds.genes"));
		BufferedWriter pvalsbw = new BufferedWriter(new FileWriter("pvals.genes"));
		BufferedWriter zscoresbw = new BufferedWriter(new FileWriter("zscores.genes"));
		Map<RefSeqGene,Double> pvalues = new HashMap<RefSeqGene, Double>();
		// For each gene or annotation
		for(RefSeqGene align: filteredGenes){

			// DECLARE A LIST FOR ITS KL DIVERGENCES
			if(!((align.getOrientation().equals("+"))||(align.getOrientation().equals("-")))){
				System.out.println("Error getting data for "+align.getName()+" orientation "+align.getOrientation());
			}
			else{
				//System.out.println("Getting data for "+align.getName()+" orientation "+align.getOrientation());
				if(this.hasDataForChromosome(align.getChr())){
					//counts.put(align, getDataForAlignment(align.get3PrimeExon()));
					window = 300;
					// 1. Get the counts for the last exon
					// 2. Calculate the histograms to left and right of the point
					// 3. Calculate the chi-squared distance between the two histograms
					// 4. Store in array
					
					// 1. Get the counts for the last exon
					//double[] alignData = getDataForAlignment(align.get3PrimeExon());	
					/**
					  GET COUNTS FOR ENTIRE GENE AND WINDOW SIZE UPSTREAM AND DOWNSTREAM
					 
					int downExtension = window;
					int upExtension = window;
					
					Alignments end = null;
					Alignments start = null;
					//IRRESPECTIVE OF THE STRAND JUST GET GENOMIC EXTENSIONS
					//GENOMIC SPACE START AND END
					start = new Alignments(align.getChr(), align.getStart() - window, align.getStart());
					end = new Alignments(align.getChr(), align.getEnd(), align.getEnd() + window);
					
					 // Get an interval tree for all/any exons that overlap with the extended region
					 
					IntervalTree<RefSeqGeneWithIsoforms> endOverlappersTree = genes.getOverlappers(end);
					
					 // If there is an overlap with a gene
					 
					if(!endOverlappersTree.isEmpty()){
						IteGeneeqGeneWithIsoforms> overlappersIter = endOverlappersTree.valueIterator();
						boolean overlapperIsSameGene = true;
						
						 //while the gene is the same gene
						 
						while(overlappersIter.hasNext() && overlapperIsSameGene){
							RefSeqGene overlapper = overlappersIter.next();
							//compare the end coordiantes of the gene
							overlapperIsSameGene = (overlapper.getOrientedEnd() == align.getOrientedEnd());
							// If some overlap exists, 
							// Because the extended region cannot overlap another annotation
							if(!overlapperIsSameGene){
								RefSeqGene overlap = align.getOverlap(overlapper);
								if(!(overlap==null) && align.getEnd()<overlap.getStart()){
									downExtension = Math.min((overlap.getStart()-align.getEnd()),downExtension);
								}
								//overlapping gene end is inside align
								else{
									downExtension = 0;
								}
							}
						}
					}
					endOverlappersTree = genes.getOverlappers(start);
					if(!endOverlappersTree.isEmpty()){
						Iterator<RefSeqGeneWithIsoforms> overlappersIter = endOverlappersTree.valueIterator();
						boolean overlapperIsSameGene = true;
						
						 // while the gene is the same gene
						 
						while(overlappersIter.hasNext() && overlapperIsSameGene){
							RefSeqGene overlapper = overlappersIter.next();
							//compare the end coordiantes of the gene
							overlapperIsSameGene = (overlapper.getOrientedEnd() == align.getOrientedEnd());
							// If some overlap exists, 
							// Because the extended region cannot overlap another annotation
							if(!overlapperIsSameGene){
								RefSeqGene overlap = align.getOverlap(overlapper);
								if(!(overlap==null) && align.getStart()>overlapper.getEnd()){
									upExtension = Math.min((align.getStart()-overlap.getEnd()),upExtension);
								}
								//overlapping gene start is inside align
								else{
									upExtension = 0;
								}
							}
						}
					}**/
					int upExtension = 0;
					int downExtension = 0;
					Alignments end = null;
					Alignments start = null;
					//IN GENOMIC SPACE
					double[] genedata = l2a(countsMap.get(align));
					double[] alignData = new double[genedata.length+downExtension+upExtension];
					int pointer=0;
					if(upExtension>0){
						start = new Alignments(align.getChr(), align.getStart() - upExtension, align.getStart());
						double[] updata = getDataForAlignment(start);
						for(int i=0;i<updata.length;i++){
							alignData[pointer] = updata[i];
							pointer++;
						}
					}
					for(int i=0;i<genedata.length;i++){
						alignData[pointer] = genedata[i];
						pointer++;
					}
					if(downExtension>0){
						end = new Alignments(align.getChr(), align.getEnd(), align.getEnd() + downExtension);
						double[] downdata = getDataForAlignment(end);
						for(int i=0;i<downdata.length;i++){
							alignData[pointer] = downdata[i];
							pointer++;
						}
					}
					
					/*
					 * END OF GET EXTENDED GENE COUNTS
					 */
					
					//If the length of the gene is shorter than twice the window, we reset window size
					if(alignData.length<(2*window) && alignData.length>0){
						window = (int)(alignData.length/2.0);
					}
					
					//int bins = window/5;
					int bins = 20;
					List<Double> divergence = new ArrayList<Double>();
					List<Double> distance = new ArrayList<Double>();
					List<Double> Chipvalues = new ArrayList<Double>();
						
					//bw.write("Gene: "+align.getName()+"\n");
					//System.out.println("Gene: "+align.getName()+" Length of total region: "+alignData.length+" Transcript Length: "+align.getTranscriptLength()+" Window: "+window+" Upstream:"+ upExtension+ "Downstream: "+downExtension);
					//System.out.println("Relative positions for KLs: ");
					// For every point in the last "exon" starting at 5'end
					int point = window;
						
					//	System.out.println("Last exon length: "+alignData.length+" Window: "+window);
					while(point<(alignData.length-window)){
						
						//System.out.print((point-upExtension)+" ");
						//bw.write("At point:"+(point-upExtension)+"\n");
						// ENSURING THE SAME RANGE IS USED FOR BOTH DISTRIBUTIONS
						double[] rightArr = Arrays.copyOfRange(alignData, point-window,point);
						double[] leftArr = Arrays.copyOfRange(alignData, point,point+window);
						double[] minMax=minMax(leftArr,rightArr);
						
						EmpiricalDistribution rightDist = new EmpiricalDistribution(rightArr,bins,0.0,minMax[1]);
						EmpiricalDistribution leftDist = new EmpiricalDistribution(leftArr,bins,0.0,minMax[1]);
						
						// ADD PSUEDOCOUNTS
						rightDist.addPsuedocounts(1.0);
						leftDist.addPsuedocounts(1.0);
						
						//KL DIVERGENCE
						//double div = rightDist.KLDivergence(leftDist,bw);
						double KLdiv = leftDist.KLDivergenceSym(rightDist);
						divergence.add(KLdiv);
						// To write an intermediate bed file
						/*
						 * Length of trnascript = T
						 * Length of last exon = L
						 * point = x
						 * T-L+x gives transcript location
						 */	
						//int pos = align.transcriptToGenomicPosition(align.getTranscriptLength()-alignData.length+point);
						int pos;
						if(align.getOrientation().equals("+")){
							pos = align.transcriptToGenomicPosition(point-upExtension);
							if(pos<0){
								System.out.println("Gene: "+align.getName()+" Transcript Length: "+align.getTranscriptLength()+" Point: "+point+" Upextension "+upExtension+" Pos is :("+point+"-"+upExtension+")");
							}
						}
						else{
							pos = align.transcriptToGenomicPosition(align.getTranscriptLength()-(point-upExtension+1));
							if(pos<0){
								System.out.println("Gene: "+align.getName()+" Transcript Length: "+align.getTranscriptLength()+" Point: "+point+" Upextension "+upExtension+" Pos is : "+align.getTranscriptLength()+"-("+point+"-"+upExtension+")");
							}
						}
						
						
						BED entry = new BED((align.getName()),align.getChr(),pos,pos);
						entry.setScore(KLdiv);
						bedsklDivergence.add(entry);
						
						double chiDist = leftDist.chiSquareDistance(rightDist);
						distance.add(chiDist);
						BED entryChi = new BED((align.getName()),align.getChr(),pos,pos);
						entryChi.setScore(chiDist);
						bedsChiDistance.add(entryChi);
						
						double chiPvalue = leftDist.testGoodnessOfFit(rightDist);
						//chiPvalue = -Math.log(chiPvalue);
						Chipvalues.add(chiPvalue);
						BED entryPval = new BED((align.getName()),align.getChr(),pos,pos);
						entryPval.setScore(chiPvalue);
						bedsChiPvalues.add(entryPval);
						
						point += step;
							
					}
					//System.out.println();
					resultMap.put(align, divergence);
			kldsbw.write(align.getName()+"\t");
					for(int i=0;i<divergence.size();i++){
						kldsbw.write(divergence.get(i)+"\t");
					}
					kldsbw.write("\n");
					List<Double> D = new ArrayList<Double>();
					for(int i=0;i<divergence.size();i++){
						D.add(i, divergence.get(i));
					}
					pvalsbw.write(align.getName()+"\t");
					for(int i=0;i<divergence.size();i++){
						double p = Statistics.pvalue(D,divergence.get(i));
						pvalsbw.write(p+"\t");
					}
					pvalsbw.write("\n");
					
					zscoresbw.write(align.getName()+"\t");
					for(int i=0;i<divergence.size();i++){
						double z = Statistics.zScore2(divergence.get(i), l2a(D));
						zscoresbw.write(z+"\t");
					}
					zscoresbw.write("\n");
					
					double max = 0.0;
					//int region = (align.get3PrimeExon().length()-(window))/step;
					int region = divergence.size()/3;
					if(align.getOrientation().equals("-")){
						if(divergence.size()<region)
							max = Statistics.max(divergence);
						else
							max = Statistics.max(divergence.subList(0, region));
					}
					else{
						if((divergence.size()-region)>divergence.size())
							max = Statistics.max(divergence.subList((divergence.size()-region), divergence.size()));
						else
							max = Statistics.max(divergence);
					}
					//buffwr.write("P-values based on raw distribution\n");
					double m = Statistics.pvalue(divergence,max);
					EmpiricalDistribution maxKLsDist = new EmpiricalDistribution(divergence,100);
					double m2 = maxKLsDist.getPValue(max, 0.0);
					double m3 = maxKLsDist.getPValue2(max);
					double zscore = Statistics.zScore2(max, l2a(divergence));
					double pval = 1.0 - Statistics.zscoreToPvalue(zscore);
					double zscore2 = maxKLsDist.getZscore(max);
					double pval2 =1.0 - Statistics.zscoreToPvalue(zscore2);
					pvalbw.write(align.getName()+"\t"+(new Double(m).toString())+"\t"+(new Double(m2).toString())+"\t"+(new Double(m3).toString())+"\t"+(new Double(zscore).toString())+"\t"+(new Double(pval).toString())+"\t"+(new Double(zscore2).toString())+"\t"+(new Double(pval2).toString())+"\n");
					//buffwr.write("P-values based on empirical distribution\n");
					
					pvalues.put(align, pval);
				}
				else
					rtrn.addRefSeq(align);
			}
		}
		pvalbw.close();
		kldsbw.close();
		/*
		 * WRITE THE BED FILES FOR KL SYM DIVERGENCES, CHI-Sq distances AND CHI SQ P-VALUES
		 */
		BufferedWriter bedBw = new BufferedWriter(new FileWriter(filename+".kl.sym.bin20.bed"));
		for(BED entry:bedsklDivergence){
			bedBw.write(entry.toString());
			bedBw.newLine();
		}
		bedBw.close();
		bedBw = new BufferedWriter(new FileWriter(filename+".chi.bin20.bed"));
		for(BED entry:bedsChiDistance){
			bedBw.write(entry.toString());
			bedBw.newLine();
		}
		bedBw.close();
		bedBw = new BufferedWriter(new FileWriter(filename+".chi.pvalues.bin20.bed"));
		for(BED entry:bedsChiPvalues){
			bedBw.write(entry.toString());
			bedBw.newLine();
		}
		bedBw.close();
		//bw.close();
		
		/*
		 * NOW, TO OUTPUT THE DISTIBUTION OF MAX KLD
		 */
		BufferedWriter buffwr = new BufferedWriter(new FileWriter("max.klds"));
		Map<RefSeqGene, Double> maxKLs = new HashMap<RefSeqGene, Double> ();
		for(RefSeqGene gene:filteredGenes){
			double m = Statistics.max(resultMap.get(gene));
			maxKLs.put(gene, m);
			buffwr.write(gene.getName()+"\t"+(new Double(m).toString())+"\n");
		}
		buffwr.close();
		
		Map<RefSeqGene, Double> LnPval = new HashMap<RefSeqGene, Double>();
		for(RefSeqGene g:pvalues.keySet()){
			LnPval.put(g,Math.log(pvalues.get(g))/Math.log(10));
		}
		buffwr = new BufferedWriter(new FileWriter("pvalues.of.pvalues"));
		for(RefSeqGene g:LnPval.keySet()){
			double pp = Statistics.pvalue(new ArrayList(LnPval.values()), LnPval.get(g));
			buffwr.write(g.getName()+"\t"+(new Double(pp).toString())+"\n");
		}
		buffwr.close();
		/*
		 * WRITE THE EMPIRICAL DISTRIBUTION - 100 bins
		 */
		/*buffwr = new BufferedWriter(new FileWriter("max.klds.distribution"));
		EmpiricalDistribution maxKLsDist = new EmpiricalDistribution(maxKLs.values(),100);
		maxKLsDist.write(buffwr);
		buffwr.close();
		*/
		/*
		 * OUTPUT THE P-VALUES FOR ALL GENES (on the maxKLs)
		 */
/*		buffwr = new BufferedWriter(new FileWriter("max.klds.pvalues"));
		//buffwr.write("P-values based on raw distribution\n");
		for(RefSeqGene gene:filteredGenes){
			double m = Statistics.pvalue(new ArrayList(resultMap.get(gene)),maxKLs.get(gene));
			buffwr.write(gene.getName()+"\t"+(new Double(m).toString())+"\n");
		}
		//buffwr.write("P-values based on empirical distribution\n");
		buffwr.close();
		buffwr = new BufferedWriter(new FileWriter("max.klds.distribution.pvalues"));
		for(RefSeqGene gene:filteredGenes){
			double m = maxKLsDist.getPValue(maxKLs.get(gene), 0.0);
			buffwr.write(gene.getName()+"\t"+(new Double(m).toString())+"\n");
		}
		//buffwr.write("P-values based on empirical distribution(P(x>i))\n");
		buffwr.close();
		buffwr = new BufferedWriter(new FileWriter("max.klds.distribution.calculated.pvalues"));
		for(RefSeqGene gene:filteredGenes){
			double m = maxKLsDist.getPValue2(maxKLs.get(gene));
			buffwr.write(gene.getName()+"\t"+(new Double(m).toString())+"\n");
		}
		buffwr.close();*/
		
		return genes;
	}
	
	/**
	 * This function returns the min and max values (combined) in two arrays
	 * @param arr1
	 * @param arr2
	 * @return
	 */
	private double[] minMax(double[] arr1, double[] arr2){
		double min=Double.MAX_VALUE;
  		double max=-Double.MAX_VALUE;
  		
  		for(int i=0; i<arr1.length; i++){
  			min=Math.min(arr1[i], min);
  			max=Math.max(arr1[i], max);
  		}
  		for(int i=0; i<arr2.length; i++){
  			min=Math.min(arr2[i], min);
  			max=Math.max(arr2[i], max);
  		}
  		
  		double[] minMax={min, max};
  		return minMax;
	}
	
	/**
	 * This function will filter the genes in the input Map by a number of constraints on the total number of reads on the genes
	 * 1. All genes with expression < 5.0 read counts will be removed
	 * 2. All genes with expression in the 75th quantile of the remaining expressed genes will be returned 
	 * @param countsMap: A map of RefSeqGene to a List of counts along the length of the gene 
	 * @return
	 */
	public List<RefSeqGene> filterGenesByExpression(Map<RefSeqGene, List<Double>> countsMap){
		
		Map<RefSeqGene, Double> SumCountsMap=new HashMap<RefSeqGene, Double>();
		List<RefSeqGene> genes = new ArrayList<RefSeqGene>();
		/*
		 * STEP 1 : GET A DISTRIBUTION OF TOTAL COUNTS OF EACH GENE
		 */
		int removed = 0;
		for(RefSeqGene gene:countsMap.keySet()){
			/*
			 * STEP 2: DO NOT ADD ANY GENES WITH NUMBER OF COUNTS <5
			 */
			double sum = Statistics.sum(countsMap.get(gene));
			if(sum>=5.0){
				SumCountsMap.put(gene, sum);
			}
			else{
				removed ++;
			}
		}
		
		//TO DO: CHECK the quantile function needs an ordered list
		System.out.println(removed+" genes removed because of expression < 5 reads. Considering remaining "+SumCountsMap.size()+" genes.");
		ArrayList<Double> C = new ArrayList(SumCountsMap.values());
		Collections.sort(C);
		double quant = Statistics.quantile(C, 0.25);
		
		for(RefSeqGene gene:SumCountsMap.keySet()){
			if(SumCountsMap.get(gene)>quant){
				genes.add(gene);
			}
		}
		System.out.println((SumCountsMap.size()-genes.size())+" genes removed because of expression < 25th quantile. Considering remaining "+genes.size()+" genes.");
		return genes;
	}
	

	private Alignments trimEnds(Alignments align, double[] array){
		double min = Math.max(minNumberOfReadsAtEnd, Statistics.quantile(array,trimQuantile));
		//double[] array=this.getDataForAlignments(align);
		int trimStart = align.getStart() + MaximumContiguousSubsequence.contiguousStartSubSequenceOverMin(array, min);
		int trimEnd   = align.getStart() + MaximumContiguousSubsequence.contiguousEndSubSequenceOverMin(array, min);
		Alignments newAlign=null;
		if(trimEnd > trimStart){
			if(trimEnd - trimStart < minAnnotationSize) {
				int toAdd = minAnnotationSize - trimEnd + trimStart + 1;
				trimEnd = trimEnd + toAdd/2;
				trimStart = trimStart - toAdd/2;
			}
			newAlign=new Alignments(align.getChr(),trimStart, trimEnd);
		} 
		//System.err.println("region " + align.toUCSC() + " min coverage " + min + " trimStart " + trimStart + " trimEnd  " + trimEnd);
		return newAlign;
	}

	private RefSeqGene trimEnds(RefSeqGene align, List<Double> data, double quantile){
		double[] array=l2a(data);
	//	long start=System.currentTimeMillis();
		Collections.sort(data);
	//	long end=System.currentTimeMillis();
	//	System.err.println("Sorting in "+(end-start)+"ms.");
		
		double cutoff = Math.max(minNumberOfReadsAtEnd, Statistics.quantile(data, quantile));
		//System.err.print("Trimming " + align.getName() + " (" + align.getOrientation() + ") cutoffs " + cutoff);
		//double[] array=this.getDataForAlignments(align);
		int trimStart = MaximumContiguousSubsequence.contiguousStartSubSequenceOverMin(array, cutoff);
		int trimEnd   =  MaximumContiguousSubsequence.contiguousEndSubSequenceOverMin(array, cutoff);
		logger.debug(" trimStart " + trimStart + " trimEdn " + trimEnd + ", transcript length " + align.getTranscriptLength());


		// We only want to trim the last exons not the cut into spliced ones
		if(align.getNumExons() > 1 ) {
			//int genomicTrimStart = align.transcriptToGenomicPosition(trimStart);
			//int genomicTrimEnd   = align.transcriptToGenomicPosition(trimEnd);
			
			int transcriptFirstExonEnd = "-".equals(align.getOrientation())  
					? align.genomicToTranscriptPosition(align.getLastExon().getStart())
					: align.genomicToTranscriptPosition(align.getFirstExon().getEnd() - 1);

			int transcriptLastExonStart = "-".equals(align.getOrientation()) 
					? align.genomicToTranscriptPosition(align.getFirstExon().getEnd() - 1)
					: align.genomicToTranscriptPosition(align.getLastExon().getStart() );
			
			logger.trace("first exon end in transcript " + transcriptFirstExonEnd + " last exon start " + transcriptLastExonStart);

			if(trimStart > transcriptFirstExonEnd) {
				trimStart = Math.max(0, transcriptFirstExonEnd - 50);
			}
			
			if(trimEnd < transcriptLastExonStart) {
				trimEnd = Math.min(align.getTranscriptLength(), transcriptLastExonStart + 50);
			}

			if("-".equals(align.getOrientation()) ){   
				int tmpTrimStart = trimStart;
				trimStart = align.getTranscriptLength() - trimEnd;
				trimEnd = align.getTranscriptLength() - tmpTrimStart;
			}
			
			logger.trace("Reset trimStart and TrimEnd to  " + trimStart + " - " + trimEnd);
			
			//trimStart = genomicTrimStart > align.getFirstExon().getEnd() ? 1 : trimStart;
			//trimEnd = genomicTrimEnd < align.getLastExon().getStart() ? align.getTranscriptLength() : trimEnd;
			//logger.debug("genomic trim end: " + genomicTrimEnd + " genomic trim start " + genomicTrimStart +/* " lastExon start: " + lastExonStart + " firstExonEnd: " + firstExonEnd +*/ " Reset trimStart and TrimEnd to  " + trimStart + " - " + trimEnd);
		}
	
		
		RefSeqGene newAlign=align;
		
		if(trimEnd > trimStart && (1 > trimStart || trimEnd < align.getTranscriptLength()-1)){
			newAlign=align.trim(trimStart, trimEnd);
			logger.debug("trimming ("+trimStart +" - "+ trimEnd+") gene was: " + align.toBED() + " and now is: " +newAlign.toBED());
		} 
		return newAlign;
	}

	private Collection<RefSeqGene> trimGenes(Collection<RefSeqGene> set, double quantile) throws IOException{
		return trimEndsForGenes(set, quantile);
	}

	public RefSeqGene trimGeneEnds(RefSeqGene gene, double quantile) throws IOException{
		long start=System.currentTimeMillis();
		//loop through each alignment and compute number of reads for the ends
		List <RefSeqGene> oneGeneList = new ArrayList<RefSeqGene>(1);
		oneGeneList.add(gene);
		List<Double> counts=getDataForGene(oneGeneList).get(gene);
		long end=System.currentTimeMillis();
	
		RefSeqGene trunc = trimEnds(gene, counts, quantile);
		return trunc;
	
	}

	private static void writeFullBED(String save, Map<Path, double[]> segments)throws IOException{
		FileWriter writer=new FileWriter(save);
	
		for(Path align: segments.keySet()){
			RefSeqGene gene=align.toGene();
			double[] vals=segments.get(align);
			for(double val : vals) {
				gene.addExtraField(String.valueOf(val));
			}
			//gene.setName(align.getLocalLambda()+"_"+vals[0]+"_"+vals[7]);
			//gene.setName(""+vals[2]);
			writer.write(gene+"\n");
		}
	
		writer.close();
	}

	public static void writeFullBED(String save, Map<RefSeqGene, double[]> map, double alpha)throws IOException{
		FileWriter writer=new FileWriter(save);
	
		for(RefSeqGene align: map.keySet()){
			double[] ps=map.get(align);
			if(ps[0]<alpha){
				writer.write(align.toBED());
				for(double p : ps) {
					writer.write("\t"+p);
				}
				writer.write("\n");
			}else {
				logger.debug("Gene " + align.getName() + " in " + align.toUCSC() + " alpha " + ps[0] +" not printed");
			}
	
		}
	
		writer.close();
	}

	public static void writeFullBED(String save, Collection<RefSeqGene> segments)throws IOException{
		BufferedWriter writer=new BufferedWriter(new FileWriter(save));
	
		writeFullBED(writer, segments);
	
		writer.close();
	}

	private static void writeFullBED(BufferedWriter writer, Collection<RefSeqGene> segments)throws IOException{
		for(RefSeqGene align: segments){
			writer.write(align.toBED());
			writer.newLine();
		}
	}

	public void writeConnections(String save, ChromosomeWithBubblesJGraphT graph)throws IOException{
		FileWriter writer=new FileWriter(save);
	
		IntervalTree<BubbleEdge> edges=new IntervalTree<BubbleEdge>();
	
		for(BubbleEdge edge: graph.edgeSet()){
			edges.put(edge.getConnection().getStart(), edge.getConnection().getEnd(), edge);
		}
	
		Collection<Path> paths=graph.getPaths(1);
		Collection<Alignments> transcribedRegions=collapsePaths(paths);
	
		for(Alignments transcript: transcribedRegions){
			writer.write("GENE: "+transcript.toUCSC()+"\n");
			Iterator<Node<BubbleEdge>> edgeIter=edges.overlappers(transcript.getStart(), transcript.getEnd());
			Collection<LightweightGenomicAnnotation> exons=new TreeSet<LightweightGenomicAnnotation>();
			while(edgeIter.hasNext()){
				BubbleEdge edge=edgeIter.next().getValue();
				if(edge.getType().equals(EdgeSourceType.SPLICED)){
					VertexPair<LightweightGenomicAnnotation> nodes=graph.getNodePair(edge);
					LightweightGenomicAnnotation first=nodes.getFirst();
					LightweightGenomicAnnotation second=nodes.getSecond();
					writer.write(first.toUCSC()+"\t"+second.toUCSC()+"\t"+edge.getSplicedCounts()+"\t"+data.getSpliceWeightFactor()+"\n");
					exons.add(first);
					exons.add(second);
				}
			}
			for(LightweightGenomicAnnotation exon: exons){
				writer.write(exon.toUCSC()+"\t"+scanPRate(exon)[2]+"\n");
			}
		}
	
		writer.close();
	}

	private static Map<RefSeqGene, double[]> PathsToGenes(Map<Path, double[]> pathsScores) {
	
		Map<RefSeqGene, double[]> rtrn=new HashMap<RefSeqGene, double[]>();
		for (Path p: pathsScores.keySet()){
			RefSeqGene g = p.toGene();
			g.setCountScore(pathsScores.get(p)[4]);
			rtrn.put(g, pathsScores.get(p));
		}
		return rtrn;
	}

	private static Map<Path, RefSeqGene> pathsToGenesMap(Map<Path, double[]> pathsScores) {
		Map<Path, RefSeqGene> rtrn=new TreeMap<Path, RefSeqGene>();
	
		for(Path gene: pathsScores.keySet()){
			double[] p=pathsScores.get(gene);
			//Gene g=gene.toGene();
			RefSeqGene g=gene.toGene(p);
			g.setBedScore(p[4]);
			rtrn.put(gene, g);
	
		}
	
		return rtrn;
	}

	private static Collection<RefSeqGene> pathToGenesCollection(Map<Path, double[]> pathsScores) {
		Collection<RefSeqGene> rtrn = new TreeSet<RefSeqGene>();
	
		for (Path gene : pathsScores.keySet()) {
			double[] p=pathsScores.get(gene);
			//Gene g=gene.toGene();
			RefSeqGene g=gene.toGene(p);
			g.setBedScore(p[4]);
			rtrn.add(g);
		}
	
		return rtrn;
	}

	public boolean hasDataForChromosome(String chr) { return this.chromosomeLengths.containsKey(chr);}

	private boolean overlaps(Alignment read, Alignments current) {
		Alignments t=new Alignments(read.getChromosome(), read.getAlignmentStart(), read.getAlignmentEnd());
		return current.overlaps(t);
	}

	private boolean isPaired() {
		return this.data.isPaired();
	}

	private void localSegment(Collection<Path> allPaths, double alpha, boolean onlyNonSpliced) throws IOException{
		Map<String, IntervalTree<Path>> trees=TreeUtils.makeIntervalTreeByPath(allPaths);
		IntervalTree<LightweightGenomicAnnotation> cache = new IntervalTree<LightweightGenomicAnnotation>();
		//for each gene
		for(Path gene: allPaths){
			//get all overlapping pieces
			Iterator<Node<Path>> overlappers=trees.get(gene.getChromosome()).overlappers(gene.getStart(), gene.getEnd());
			IntervalTree<Alignment> tree=data.getIntervalTreeCached(gene.getChromosome(), gene.getStart(), gene.getEnd());
			double localLambda=Math.max(data.getLambda(gene.getChromosome()), computeLambdaFromPath(overlappers, tree, gene, cache));
			gene.setLocalLambda(localLambda);
	
		}
	
	}

	/*
	private static Collection<EdgeSourceType> typesToFollow(boolean followPairedEnds) {
		ArrayList<EdgeSourceType> types = new ArrayList<EdgeSourceType>();
		types.add(EdgeSourceType.SPLICED);
		if(followPairedEnds) {
			types.add(EdgeSourceType.PAIRED);
		}
		return types;
	}
	 */
	
	public static void main(String[] args) throws IOException, ParseException, InterruptedException, broad.core.error.ParseException, IllegalArgumentException, MathException {
		Globals.setHeadless(true);
		System.out.println("Using Version R4.4");
		logger.debug("DEBUG ON");
		ArgumentMap argmap = CLUtil.getParameters(args, usage, "full");
		double lambda = argmap.isPresent("lambda") ? argmap.getDouble("lambda") : 0;
		if("full".equalsIgnoreCase(argmap.getTask())) {
			//Mandatory parameters
			String alignmentFile = argmap.getMandatory("alignment");
			int [] defaultWindows = {0};
			double minSpliceFrequency = argmap.containsKey("minSpliceFrequency") ? argmap.getDouble("minSpliceFrequency") : 0.1;
			int[] fixedWidth=argmap.containsKey("windows")? getWidths(argmap.getMandatory("windows")): defaultWindows;
			String sizes =  argmap.get("sizeFile");			
			int minMappingQuality = argmap.containsKey("minMappingQuality") ? argmap.getInteger("minMappingQuality") : DEFAULT_MIN_MAPPING_QUALITY;
			//Optional parameters
			File[] maskFiles  = argmap.isPresent("maskFileDir") ?  new File(argmap.get("maskFileDir")).listFiles(): null ;
			Map<String, Integer> maskFileData = parseMaskFiles(maskFiles);
			String chr =  argmap.get("chr");
	
			boolean trimEnds = argmap.containsKey("trim");
			boolean filterCanonical=!argmap.containsKey("dontFilterCanonicalSplice");
			double trimQuantile = argmap.isPresent("trimQuantile") ? argmap.getDouble("trimQuantile") : 0.25;
			String sequenceFile =  argmap.get("chrSequence");
			double alpha = argmap.containsKey("alpha") ? argmap.getDouble("alpha") : .05;
			int chrStart = argmap.containsKey("start") ? argmap.getInteger("start") :0;
			int chrEnd   = argmap.containsKey("end") ? argmap.getInteger("end") : Integer.MAX_VALUE;
			int minimumSpliceReadSupport = argmap.containsKey("minSpliceSupport") ? argmap.getInteger("minSpliceSupport") : 2;

			ContinuousDataAlignmentModel data = AlignmentUtils.loadAlignmentData(alignmentFile, true, minMappingQuality, true, false ); //Right now making default to not weight reads by mapping hits and to ignore duplicated reads.
			data.setMinNumberOfSplices(minimumSpliceReadSupport);
			data.setMaskFileData(maskFileData);
			data.setExtensionFactor(0);
			/*AlignmentDataModel alignments=new GenericAlignmentDataModel(alignmentFile, sizes, false, minMappingQuality); //Does this solve the scoring issue?
			logger.info("AlignmentDataModel loaded, initializing model stats");
			AlignmentDataModelStats alignmentData = new AlignmentDataModelStats(alignments, maskFileData, chr, lambda, true);
			logger.info("model stats loaded, initializing model");
			long start=System.currentTimeMillis();
			ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(alignmentData, maskFileData, 0, minimumSpliceReadSupport, pairedEndData, strandSpecificData, upweightSplices);
			logger.info("Built the model: "+(System.currentTimeMillis()-start)/1000.0 + " free memory: " + Runtime.getRuntime().freeMemory());
			*/
			data.trimEnds = trimEnds;
			BEDFileParser segmentation = data.segmentChromosome(minSpliceFrequency, fixedWidth, chr, trimEnds, filterCanonical, sequenceFile,
					alpha, chrStart, chrEnd, trimQuantile);
			segmentation.merge();
			BufferedWriter bw = argmap.getOutputWriter();
			segmentation.writeFullBed(bw);
			bw.close();
			
		}  else if ("addpairs".equalsIgnoreCase(argmap.getTask())) {
			String chr =  argmap.get("chr");
			double minMappingQuality = argmap.containsKey("minMappingQuality") ? argmap.getDouble("minMappingQuality") : DEFAULT_MIN_MAPPING_QUALITY;
			double alpha = argmap.isPresent("alpha") ? argmap.getDouble("alpha") : 1.5; //Default is NOT to filter by significance
			String annotationFile = argmap.getMandatory("annotations");
			BEDFileParser annotations = new BEDFileParser(annotationFile);
			String alignmentFile = argmap.getMandatory("alignment");
			logger.info("Constructing paired end data model");
			String out = argmap.getOutput();
			String sizes =  argmap.get("sizeFile");
			File[] maskFiles  = argmap.isPresent("maskFileDir") ?  new File(argmap.get("maskFileDir")).listFiles(): null ;
			Map<String, Integer> maskFileData = parseMaskFiles(maskFiles);
			AlignmentDataModel alignmentDataModel=new GenericAlignmentDataModel(alignmentFile, sizes, false, minMappingQuality); //Does this solve the scoring issue?
			logger.info("AlignmentDataModel loaded, initializing model stats");
			AlignmentDataModelStats alignmentDataStats = new AlignmentDataModelStats(alignmentDataModel, maskFileData, chr, lambda, true/*false*/);
	
			ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(alignmentDataStats, maskFileData, 0, 1, null, null, false);
	
			logger.info("Updating graphs with Paired end information ...");
			PairedEndDistribution ped = data.computeInsertSizeDistribution(annotations);
			//BEDFileParser updatedAnnotations = data.updateReconstructionsWithPairs(annotations, ped);
			//BEDFileParser mergedAnnotations = data.mergeReconstructionWithPairs(updatedAnnotations, ped);
			BEDFileParser mergedAnnotations = ped.mergeReconstructionWithPairs(annotations);
			logger.info("Done adding paired end information");
			Map<RefSeqGene, double[]> rescoredAnnotations = data.scoreGenes(mergedAnnotations.GetGenes());
			writeFullBED(out, rescoredAnnotations, alpha); 
		}   
		else if("fastScore".equalsIgnoreCase(argmap.getTask())) {
			BEDFileParser annotations = new BEDFileParser(argmap.getMandatory("annotations"));
			String chr =  argmap.getMandatory("chr");
			double alpha = argmap.isPresent("alpha") ? argmap.getDouble("alpha") : 0.05;
			String outFile = argmap.getOutput();
			BufferedReader br = argmap.getInputReader();
			ChromosomeWithBubblesJGraphT graph = ExtractRegionsFromDOT.buildFromDOT(br);
			Collection<Path> paths= refSeqGeneToPath (annotations.getChrTree(chr), graph);
			Map<Path, double[]> pathsScores = ContinuousDataAlignmentModel.scorePaths(paths, graph.getLambda(), graph.getNumberOfMarkers(), graph.getNumberOfReads(), graph.getRPKMConstant(), graph.getLocalRate(), alpha);
			Map<RefSeqGene, double[]> geneScores = PathsToGenes(pathsScores);
			writeFullBED(outFile, geneScores,1.1);
	
		}	else if ("score".equalsIgnoreCase(argmap.getTask())) {
			String alignmentFile = argmap.getMandatory("alignment");
			boolean useConstituentExons = argmap.containsKey("useConstituentExons");
			boolean useConstituentIntrons = argmap.containsKey("useConstituentIntrons");
			double minMappingQuality = argmap.containsKey("minMappingQuality") ? argmap.getDouble("minMappingQuality") : DEFAULT_MIN_MAPPING_QUALITY;
	
			//Map<String, Collection<RefSeqGene>> annotations = BEDFileParser.loadDataByChr(new File(argmap.getInput()));
			String annotationFile = argmap.getInput();
			BEDFileParser annotationParser =  annotationFile.endsWith(".gtf") || annotationFile.endsWith(".GTF")? new GTFFileParser(annotationFile) : new BEDFileParser(annotationFile);
			Map<String, Collection<RefSeqGene>> annotations = null;
			if(useConstituentExons) {
				annotationParser.makeGenes(0.1);
				annotations = annotationParser.toConstituentIsoformMap();
				//annotationParser.writeFullBed("Madegenes.bed");
				BufferedWriter ciw = new BufferedWriter(new FileWriter("constituentExons.bed"));
				for(String chr : annotations.keySet()) {
					Collection<RefSeqGene> constituentIsoforms = annotations.get(chr);
					for (RefSeqGene g : constituentIsoforms) {
						ciw.write(g.toBED());
						ciw.newLine();
					}
				}
				ciw.close();
			} else if (useConstituentIntrons){
				annotationParser.makeGenes(0.1);
				//annotationParser.writeFullBed("Madegenes.bed");
				annotations = annotationParser.toConstituentIntroformMap();
			}else {
				annotations = annotationParser.toMap();
			}
			Collection<RefSeqGene> annotationCollection = new ArrayList<RefSeqGene>();
			for(Collection<RefSeqGene> chrAnnotations : annotations.values()) {
				annotationCollection.addAll(chrAnnotations);
			}
			//BEDFileParser.writeFullBED("constituentIsoforms.bed", annotationCollection);
			String save = argmap.getOutput();
			String sizes = argmap.get("sizeFile");
			File[] maskFiles  = argmap.isPresent("maskFileDir") ?  new File(argmap.get("maskFileDir")).listFiles(): null ;
			Map<String, Integer> maskFileData = parseMaskFiles(maskFiles);
			logger.info("Minimum mapping quality to count reads: " + minMappingQuality);
			boolean isStranded = argmap.containsKey("stranded");
	
			if(!isStranded) {
				logger.info("Scoring using all reads ");
				AlignmentDataModel alignments=new GenericAlignmentDataModel(alignmentFile, sizes, false, minMappingQuality);
				Map<RefSeqGene, double[]> scores=new TreeMap<RefSeqGene, double[]>();				
				runScore(annotations, save, maskFileData, alignments, scores);
				writeFullBED(save, scores, 1.1);
			} else {
				logger.info("Scoring minus reads");
				AlignmentDataModel alignments=new GenericAlignmentDataModel(alignmentFile, sizes, false, minMappingQuality);
				alignments.setNegativeStranded();
				Map<RefSeqGene, double[]> scores=new TreeMap<RefSeqGene, double[]>();				
				runScore(annotations, save, maskFileData, alignments, scores);
				writeFullBED(save+".minus", scores, 1.1);
	
				logger.info("Scoring plus reads");
				alignments.setPositiveStranded();
				scores=new TreeMap<RefSeqGene, double[]>();				
				runScore(annotations, save, maskFileData, alignments, scores);
				writeFullBED(save+".plus", scores, 1.1);
			}
		}	else if ("trim".equalsIgnoreCase(argmap.getTask())) {
			String alignmentFile = argmap.getMandatory("alignment");
			double quantile = argmap.containsKey("quantile") ? argmap.getDouble("quantile") : 0.25;
			Map<String, Collection<RefSeqGene>> annotations = BEDFileParser.loadDataByChr(new File(argmap.getInput()));
			String save = argmap.getOutput();
			String sizes = argmap.getMandatory("sizeFile");
			File[] maskFiles  = argmap.isPresent("maskFileDir") ?  new File(argmap.get("maskFileDir")).listFiles(): null ;
			Map<String, Integer> maskFileData = parseMaskFiles(maskFiles);
			AlignmentDataModel alignments=new GenericAlignmentDataModel(alignmentFile, sizes, false);
			Collection<RefSeqGene> trimmed=new TreeSet<RefSeqGene>();
			for(String chr : annotations.keySet()) {
				try{
					logger.info("processing " + chr);
					AlignmentDataModelStats alignmentData = new AlignmentDataModelStats(alignments, maskFileData, chr, false); //Setting lambda greater than one to avoid computing chromosome stats.
					ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(alignmentData, maskFileData, 0, 1, null, null);
					data.setTrimEnds(true);
					Collection<RefSeqGene> chrAnnotations = annotations.get(chr);
					logger.info("Scanning annotations");
					trimmed.addAll(data.trimEndsForGenes(chrAnnotations, quantile));
				}catch(Exception ex){logger.warn("Skipping "+chr);}
			}
			writeFullBED(save, trimmed);
		} else if ("filterbysplice".equalsIgnoreCase(argmap.getTask())) {
			int minSplices = argmap.getInteger("minSpliceSupport");
			double alpha = argmap.getDouble("alpha");
			String outFile = argmap.getOutput();
			BufferedReader br = argmap.getInputReader();
	
			ChromosomeWithBubblesJGraphT graph = ExtractRegionsFromDOT.buildFromDOT(br, minSplices);
			br.close();
	
			Collection<Path> paths = graph.getAllPaths();
	
			Map<Path, double[]> pathsScores =  ContinuousDataAlignmentModel.scorePaths(paths, graph, alpha );
			Collection<RefSeqGene> genes = filterPaths(pathsScores, alpha);
			graph.writeGraph(outFile+".dot");
			writeFullBED(outFile, pathsScores);
	
		} else if (argmap.getTask().toLowerCase().contains("extract") ) {
			String chr  = argmap.getMandatory("chr");
			int start   = argmap.getInteger("start");
			int end     = argmap.getInteger("end");
	
			File dotFile = new File(argmap.getInput());
			String save = argmap.getOutput();
			Alignments region=new Alignments(chr, start, end);
			new ExtractRegionsFromDOT(dotFile, region, save);
	
		}else if ("getIdenticalGappedReadsTranscripts".equalsIgnoreCase(argmap.getTask())) {
			String alignmentFile = argmap.getMandatory("alignment");
			Map<String, Collection<RefSeqGene>> annotations = BEDFileParser.loadDataByChr(new File(argmap.getInput()));
			String save = argmap.getOutput();
			String sizes = argmap.getMandatory("sizeFile");
			String name = argmap.getMandatory("name");
			File[] maskFiles  = argmap.isPresent("maskFileDir") ?  new File(argmap.get("maskFileDir")).listFiles(): null ;
			Map<String, Integer> maskFileData = parseMaskFiles(maskFiles);
			GenericAlignmentDataModel alignments=new GenericAlignmentDataModel(alignmentFile, sizes, false);
			ArrayList<RefSeqGene> transcripts =new ArrayList<RefSeqGene>();
	
	
			//TODO: support filters and sequence files
			/*
			Sequence chrSequence = null;
			if(chrSequenceFile != null) {
				FastaSequenceIO fsio = new FastaSequenceIO(chrSequenceFile);
	
				List<Sequence> seqs = fsio.loadAll();
	
				if(!seqs.isEmpty()) {
					chrSequence = seqs.get(0);
				} else {
					System.err.println("Sequence for " + chr + " was not found in file " + chrSequenceFile + " continuing without splice site analisys");
				}
				System.err.println("Loaded chromosome Sequence");
			}
			 */
			double total=0; 
			double passed=0;
			double biexon=0;
			for(String chr : annotations.keySet()) {
				//System.err.println("processig " + chr);
				for (RefSeqGene gene:annotations.get(chr) )
				{
					Alignments region =new Alignments(gene.getChr(), gene.getStart(), gene.getEnd());
					//boolean pass=alignments.AreSplicedReadsIdentical(region, Collection<ReadFilter> filters,  Sequence chrSeq);
					boolean pass=alignments.AreSplicedReadsIdentical(region, null,  null);
					if (pass) {transcripts.add(gene); passed++;}
					total++;
					if(gene.getNumExons()==2) {biexon++;}
				}
				//System.err.println("done. saving data  to " + save+" for "+chr);
			}	
			writeFullBED(save, transcripts);
			System.out.println(name+"\t"+passed/total+"\t"+passed/biexon);
	
		} else if (argmap.getTask().toUpperCase().contains("PAIREDFILE")){
			String pair1 = argmap.getMandatory("pair1");
			int minMappingQuality = argmap.containsKey("minMappingQuality") ? argmap.getInteger("minMappingQuality") : DEFAULT_MIN_MAPPING_QUALITY;
			boolean forChIP = argmap.containsKey("forChIP");
			String out   = argmap.getOutput();
			boolean isSorted = argmap.containsKey("sorted");
			boolean usePair2Orientation = argmap.containsKey("usePair2Orientation");
			if (argmap.containsKey("pair2")){
				String pair2 = argmap.get("pair2");
				new MapPairedEnds(new File(pair1), new File(pair2), out, isSorted);
			}
			else{
				String sizes = argmap.get("sizeFile");
				new MapPairedEndsFromSingleFile (pair1,sizes,out,"SCRIPTURE", minMappingQuality, usePair2Orientation, forChIP);
			}
	
	
		}else if ("togff".equalsIgnoreCase(argmap.getTask())) {
			String source = argmap.getMandatory("source");
			boolean toCufflinks = argmap.containsKey("cufflinks");
			String prefix = argmap.get("prefix") ;
			boolean keepids = argmap.containsKey("keepIds");
			String file = argmap.getInput();
			Map<String, Collection<RefSeqGene>> data = BEDFileParser.loadDataByChr(new File(file));
	
	
			BufferedWriter bw = argmap.getOutputWriter();
	
			for(String  chr : data.keySet() ) {
				int id = 0;
				for(RefSeqGene g : data.get(chr)) {
					if(g.getExtraFields() != null && g.getExtraFields().length > 4) {
						g.addAttribute("RPKM", g.getExtraFields()[4]);
					}
					String name = keepids ? g.getName() : (prefix != null ? prefix  : "") +  "SCRPTR."+g.getChr() + "." + id++;
					g.setName(name);
					bw.write(toCufflinks ? g.toCufflinksGTF(source, name, name, "") :  g.toGTF(source));
				}
			}
			bw.close();
		}else if ("toMISO".equalsIgnoreCase(argmap.getTask())) {
			String source = argmap.getMandatory("source");
			String annotationFile = argmap.getInput();
			BEDFileParser data = annotationFile.endsWith(".gtf") || annotationFile.endsWith(".GTF")? new GTFFileParser(annotationFile) : new BEDFileParser(annotationFile);
			data.makeGenes(0);
	
			BufferedWriter bw = argmap.getOutputWriter();
			Iterator<String> chrIt = data.getChromosomeIterator();
			while(chrIt.hasNext()) {
				String chr = chrIt.next();
				Iterator<RefSeqGeneWithIsoforms> geneIt = data.getChrTree(chr).valueIterator();
				int id = 0;
				while(geneIt.hasNext()) {
					RefSeqGeneWithIsoforms g = geneIt.next();
					bw.write(g.toMISO(source));
				}
			}
			bw.close();
		}else if (argmap.getTask().toUpperCase().contains("CHIP")){
			String out = argmap.getOutput();
			int minMappingQuality = argmap.containsKey("minMappingQuality") ? argmap.getInteger("minMappingQuality") : DEFAULT_MIN_MAPPING_QUALITY;
			String alignmentFile = argmap.getMandatory("alignment");
			int[] windows= getWidths(argmap.getMandatory("windows"));
			int extensionFactor = argmap.containsKey("extensionFactor") ? argmap.getInteger("extensionFactor") : 0;
			String sizes = argmap.get("sizeFile");
			String chr = argmap.get("chr");
			boolean printFullScores = argmap.containsKey("fullScores");
			boolean loadPairsAsFragments = argmap.containsKey("loadPairsAsFragments") || argmap.containsKey("pairedEnd");
	
			//Optional parameters
			boolean findMaxContiguous = argmap.containsKey("findMaxContiguous");
			boolean trimEnds = argmap.containsKey("trim");
			double trimQuantile = argmap.isPresent("trimQuantile") ? argmap.getDouble("trimQuantile") : 0.25;
			int minRemainingLength = argmap.isPresent("minLength") ? argmap.getInteger("minLength") : Statistics.min(windows);
			double alpha = argmap.containsKey("alpha") ? argmap.getDouble("alpha") : .05;
			File[] maskFiles  = argmap.isPresent("maskFileDir") ?  new File(argmap.get("maskFileDir")).listFiles(): null ;
			Map<String, Integer> maskFileData = parseMaskFiles(maskFiles);
	
			ContinuousDataAlignmentModel data = AlignmentUtils.loadAlignmentData(alignmentFile, true, minMappingQuality, true, false, null, loadPairsAsFragments);	
			data.setMaskFileData(maskFileData);
			data.setExtensionFactor(extensionFactor);
			//AlignmentDataModel alignments=new GenericAlignmentDataModel(alignmentFile, sizes, false, minMappingQuality); //Does this solve the scoring issue?
	
			//logger.info("AlignmentDataModel loaded, initializing model stats");
			//AlignmentDataModelStats alignmentData = new AlignmentDataModelStats(alignments, maskFileData, chr, lambda, false);
			//ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(alignmentData, maskFiles, extensionFactor, 0);
			List<String> chromosomes = new ArrayList<String>();
			if(chr!= null && chr.length() > 0) {
				chromosomes.add(chr);
			} else {
				Map<String, Integer> chrSizes = data.getChromosomeLengths();
				chromosomes = new ArrayList<String>(chrSizes.keySet());
			}
			data.setTrimEnds(trimEnds);
			data.setMinContguousSegmentSize(minRemainingLength);
			data.setTrimQuantile(trimQuantile);
			data.setFindMaxContiguous(findMaxContiguous);
	
			Map<Alignments, double[]> scores = new HashMap<Alignments, double[]>();
			int totalMappedReads = 0;
			for(String workChr : chromosomes) {
				logger.info("Processing chromosome " + workChr + ", Scanning windows");
				totalMappedReads += data.data.getSum(workChr);
				Collection<Alignments> segments = data.scan(windows, alpha, workChr);
				logger.info("Scoring segments");
				scores.putAll(data.scoreSegments(segments, workChr));
				logger.info("Printing results");
			}
			logger.info("Done. Total mapped reads: " + totalMappedReads);
			BEDFileParser.writeSortedBED(out, scores);
			if(printFullScores) {
				BEDFileParser.writeSortedBEDWithScores(out+".scores", scores);
			}
	
		} else if ("scoresegments".equalsIgnoreCase(argmap.getTask())) {
			String alignmentFile = argmap.getMandatory("alignment");
			Map<String, Collection<Alignments>> annotations = BEDFileParser.loadAlignmentDataByChr(new File(argmap.getInput()));
			String save = argmap.getOutput();
			String sizes = argmap.getMandatory("sizeFile");
			File[] maskFiles  = argmap.isPresent("maskFileDir") ?  new File(argmap.get("maskFileDir")).listFiles(): null ;
			double alpha = argmap.isPresent("alpha") ? argmap.getDouble("alpha") : 0.05;
			boolean ignoreAlpha = argmap.isPresent("ignoreAlpha");
			
			if (ignoreAlpha) {
				logger.info("Ignoring alpha - outputting all scores"); 
			} else {
				logger.info("Using alpha = " + alpha);
			}
			Map<String, Integer> maskFileData = parseMaskFiles(maskFiles);
			AlignmentDataModel alignments=new GenericAlignmentDataModel(alignmentFile, sizes, false);
			Map<Alignments, double[]> scores=new TreeMap<Alignments, double[]>();
			AlignmentDataModelStats alignmentData = new AlignmentDataModelStats(alignments, maskFileData);
			ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(alignmentData, maskFileData, 0, 1, null, null);
			data.alpha = alpha;
			for(String chr : annotations.keySet()) {
				logger.info("processing " + chr);
				Collection<Alignments> chrAnnotations = annotations.get(chr);
				scores.putAll(data.scoreSegments(chrAnnotations, chr, ignoreAlpha));
			}
			BEDFileParser.writeBEDWithScores(save, scores);
		}else if("trimSegments".equalsIgnoreCase(argmap.getTask())) { 
			String alignmentFile = argmap.getMandatory("alignment");
			double trimQuantile = argmap.isPresent("trimQuantile") ? argmap.getDouble("trimQuantile") : 0.25;
			int minRemainingLength = argmap.isPresent("minLength") ? argmap.getInteger("minLength") : 200;
			Map<String, Integer> maskFileData = null;
			if(argmap.isPresent("maskFileDir")) {
				maskFileData=parseMaskFiles(new File(argmap.get("maskFileDir")).listFiles());
			}
			Map<String, Collection<Alignments>> annotations = BEDFileParser.loadAlignmentDataByChr(new File(argmap.getInput()));
			String save = argmap.getOutput();
			String sizes = argmap.get("sizeFile");
			boolean findMaxContiguous = argmap.containsKey("findMaxContiguous");
			AlignmentDataModel alignments=new GenericAlignmentDataModel(alignmentFile, sizes, false);
			Map<Alignments, double[]> scores=new TreeMap<Alignments, double[]>();
			AlignmentDataModelStats alignmentData = new AlignmentDataModelStats(alignments, maskFileData);
			ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(alignmentData, maskFileData, 0, 1, null, null);
			data.setMinContguousSegmentSize(minRemainingLength);
			data.setTrimQuantile(trimQuantile);
			data.alpha = 1.1;
			for(String chr : annotations.keySet()) {
				logger.info("processig " + chr);
				Collection<Alignments> chrAnnotations = annotations.get(chr);
				if(findMaxContiguous) {
					chrAnnotations = data.findMaxContiguous(chrAnnotations);
				} else {
					chrAnnotations = data.trimEnds(chrAnnotations);
				}
	
				scores.putAll(data.scoreSegments(chrAnnotations, chr));
				BEDFileParser.writeBEDWithScores(save, scores);
			}
		}else if("adjustEnds".equalsIgnoreCase(argmap.getTask())) { 
			String alignmentFile = argmap.getMandatory("alignment");
			double trimQuantile = argmap.isPresent("trimQuantile") ? argmap.getDouble("trimQuantile") : 0.25;
			int minRemainingLength = argmap.isPresent("minLength") ? argmap.getInteger("minLength") : 200;
			Map<String, Integer> maskFileData = null;
			if(argmap.isPresent("maskFileDir")) {
				maskFileData=parseMaskFiles(new File(argmap.get("maskFileDir")).listFiles());
			}
			Map<String, Collection<Alignments>> annotations = BEDFileParser.loadAlignmentDataByChr(new File(argmap.getInput()));
			String save = argmap.getOutput();
			String sizes = argmap.get("sizeFile");
			boolean findMaxContiguous = argmap.containsKey("findMaxContiguous");
			AlignmentDataModel alignments=new GenericAlignmentDataModel(alignmentFile, sizes, false);
			Map<Alignments, double[]> scores=new TreeMap<Alignments, double[]>();
			AlignmentDataModelStats alignmentData = new AlignmentDataModelStats(alignments, maskFileData);
			ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(alignmentData, maskFileData, 0, 1, null, null);
			data.setMinContguousSegmentSize(minRemainingLength);
			data.setTrimQuantile(trimQuantile);
			data.alpha = 1.1;
			for(String chr : annotations.keySet()) {
				logger.info("processig " + chr);
				Collection<Alignments> chrAnnotations = annotations.get(chr);
				if(findMaxContiguous) {
					chrAnnotations = data.findMaxContiguous(chrAnnotations);
				} else {
					chrAnnotations = data.trimEnds(chrAnnotations);
				}
	
				scores.putAll(data.scoreSegments(chrAnnotations, chr));
				BEDFileParser.writeBEDWithScores(save, scores);
			}
		}else if ("polya".equalsIgnoreCase(argmap.getTask())) {
			PolyAUtils.initializePathToDirectory();
	
			if (argmap.containsKey("parameters")) {
				String parametersFile = argmap.get("parameters");
	
				PolyAUtils.initializeParameters(parametersFile);
			}
	
			Set<String> polyasubtasks = new HashSet<String>();
			String[] polyasubtasksSplit = argmap.getMandatory("polyasubtask").split(",");
			for (String polyasubtask : polyasubtasksSplit) {
				polyasubtasks.add(polyasubtask.toLowerCase());
			}
	
			if (polyasubtasks.contains("primarypipeline")) {
				polyasubtasks.add("findandprocessreadswithpolya");
				polyasubtasks.add("alignpolyaandsortandindexalignments");
				polyasubtasks.add("writeandfilterannotationsofreadswithpolya");
				polyasubtasks.add("pauseuntilallannotationsarewritten");
				polyasubtasks.add("findalignandwritematesoffullpolyareads");
				polyasubtasks.add("findandprocesscoordinatesofendoftranscripts");
				polyasubtasks.add("writeplotofqualityscoresforpolyabasedonmatchingorientationwithtranscript");
			}
	
			if (polyasubtasks.contains("findandprocessreadswithpolya"))
				PolyAUtils.findAndProcessReadsWithPolyA();
	
			if (polyasubtasks.contains("alignpolyaandsortandindexalignments"))
				PolyAUtils.alignPolyAAndSortAndIndexAlignments();
	
			if (polyasubtasks.contains("writeandfilterannotationsofreadswithpolya")) {
				if (argmap.containsKey("alignmentsForJob")) {
					String alignmentsForJob = argmap.get("alignmentsForJob");
					int jobNum = argmap.getInteger("jobnum");
					String currentDirectory = argmap.get("pathtomostfiles");
	
					PolyAUtils.setPathToCurrentDirectoryAndPathToMostFiles(currentDirectory);
	
					PolyAUtils.writeAnnotationsOfReadsWithPolyAAndFilterAndWriteFilteredAnnotationsToBedFileForJob(alignmentsForJob, jobNum);
				} else {
					String parametersFile = "";
					if (argmap.containsKey("parameters"))
						parametersFile = argmap.get("parameters");
					PolyAUtils.writeAnnotationsOfReadsWithPolyAAndFilterAndWriteFilteredAnnotationsToBedFileBySplittingIntoJobs(parametersFile);
				}
			}
	
			if (polyasubtasks.contains("pauseuntilallannotationsarewritten"))
				PolyAUtils.pauseUntilAllAnnotationsAreWrittenAndConcatenate();
	
			if (polyasubtasks.contains("filteralignmentsbyoverlapwithtranscripts"))
				PolyAUtils.filterSetOfAlingnmentsByOverlapWithTranscripts();
	
			if (polyasubtasks.contains("findalignandwritematesoffullpolyareads"))
				PolyAUtils.findAlignAndWriteMatesOfFullPolyAReads();
	
			if (polyasubtasks.contains("findandprocesscoordinatesofendoftranscripts"))
				PolyAUtils.findAndProcessCoordinatesOfEndOfTranscripts();
	
			if (polyasubtasks.contains("writeplotofqualityscoresforpolyabasedonmatchingorientationwithtranscript"))
				PolyAUtils.writePlotOfQualityScoresForPolyABasedOnMatchingOrientationWithTranscript();
	
			// bsub -o out1.bsub -R "rusage[mem=18]" -q priority -P polya java -jar /seq/regevlab/polya/runs/scripture.jar -task polya -parameters parameters.txt -polyasubtask primarypipeline
	
			/**
			//String readFile = argmap.getMandatory("reads");
			//String bowtiePath = argmap.getMandatory("bowtie");
			//String bowtieReferenceIndex = argmap.getMandatory("bowtieReferenceIndex");
			//String igvToolsPath = argmap.getMandatory("igvtools");
			//String chrSizesFile = argmap.getMandatory("sizesFile");
			//String allAnnotationsOutFile = argmap.getMandatory("allOut");
			//String filteredAnnotationsOutFile = argmap.getMandatory("filteredOut");
			//boolean sortAndIndexReads = argmap.containsKey("sortAndIndexReads");
			//String scriptureTranscriptsFile = argmap.getMandatory("transcripts");
			//String matesOfFullPolyAReadsOutFile = argmap.getMandatory("matesOut");
			//String alignmentsFile = argmap.getMandatory("alignments");
			//String endCoordinatesOutFile = argmap.getMandatory("endCoorOut");
			//String endCoordinatesOnlyWithPolyAMatesOutFile = argmap.getMandatory("endCoorOnlyWithMatesOut");
			//String transcriptLinksOutFile = argmap.getMandatory("transcriptLinksOut");
			//String modifiedTranscriptsOutFile = argmap.getMandatory("modifiedTranscriptsOut");
			//String modifiedTranscriptsOnlyWithPolyAMatesOutFile = argmap.getMandatory("modifiedTranscriptsOnlyWithMatesOut");
			//String modifiedTranscriptsOnlyWithModifiedOutFile = argmap.getMandatory("modifiedTranscriptsOnlyWithModified");
	
			//PolyAUtils.PIPELINE(readFile, bowtiePath, bowtieReferenceIndex, igvToolsPath, sortAndIndexReads, chrSizesFile, allAnnotationsOutFile, filteredAnnotationsOutFile, matesOfFullPolyAReadsOutFile, scriptureTranscriptsFile, alignmentsFile);
			// COMMAND: java -jar scripture.jar -task polya -reads ../tmp/all_reads.fastq -bowtie /seq/mguttman/scripts/BowTie/bowtie-0.12.1/bowtie -bowtieReferenceIndex /seq/genome/mouse/mouse_Mm9/mm9.nonrandom -igvtools /home/radon00/hmetsky/IGVTools/igvtools -sizesFile /seq/genome/mouse/mouse_Mm9/sizes -allOut allAlignments.bed -filteredOut filteredAlignments.bed -transcripts /home/radon00/hmetsky/work/tmp/all_transcripts.bed -matesOut matesOfFullPolyAReads.bed -alignments /home/radon00/hmetsky/work/tmp/all_alignments.sorted.sam
	
			//BEDFileParser transcripts = PolyAUtils.initiateTranscriptsParser(scriptureTranscriptsFile);
			//Map<String, Set<RefSeqGene>> m = PolyAUtils.findCoordinatesOfEndsOfTranscripts(transcripts, filteredAnnotationsOutFile, scriptureTranscriptsFile, chrSizesFile);
			//PolyAUtils.processEndOfTranscriptCoordinates(transcripts, scriptureTranscriptsFile, m, matesOfFullPolyAReadsOutFile, endCoordinatesOutFile, endCoordinatesOnlyWithPolyAMatesOutFile, transcriptLinksOutFile, scriptureTranscriptsFile, modifiedTranscriptsOutFile, modifiedTranscriptsOnlyWithPolyAMatesOutFile, modifiedTranscriptsOnlyWithModifiedOutFile);
			// COMMAND: java -jar scripture.jar -task polya -sizesFile /seq/genome/mouse/mouse_Mm9/sizes -transcripts /home/radon00/hmetsky/work/tmp/all_transcripts.bed -matesOut /home/radon00/hmetsky/work/polya/matesOfFullPolyAReads.bed -filteredOut /home/radon00/hmetsky/work/polya/filteredAlignments.bed -endCoorOut endCoordinates.bed -endCoorOnlyWithMatesOut endCoordinates_onlyWithPolyAMates.bed -transcriptLinksOut transcriptLinks.txt -modifiedTranscriptsOut transcripts_new.bed -modifiedTranscriptsOnlyWithMatesOut transcripts_new_onlyWithPolyAMates.bed -modifiedTranscriptsOnlyWithModified transcripts_new_onlyWithModified.bed
	
			//PolyAUtils.filterAlignmentsByOverlapWithLastExonOrAfterLastExon("filteredAlignments.bed", scriptureTranscriptsFile, "filteredAlignments_lastexon.bed");
			 */
	
			/**
			if (argmap.containsKey("alignmentsForJob")) {
				String alignmentsForJob = argmap.get("alignmentsForJob");
				int jobNum = argmap.getInteger("jobnum");
	
				PolyAUtils.writeAnnotationsOfReadsWithPolyAAndFilterAndWriteFilteredAnnotationsToBedFileForJob(alignmentsFile, chrSizesFile, alignmentsForJob, jobNum);
			} else {
				int numOfAlignmentsPerJob = argmap.getInteger("numalignmentsperjob");
				PolyAUtils.writeAnnotationsOfReadsWithPolyAAndFilterAndWriteFilteredAnnotationsToBedFileBySplittingIntoJobs(alignmentsFile, chrSizesFile, numOfAlignmentsPerJob);
	
				// COMMAND: java -jar ../scripture.jar -task polya -reads /home/radon00/hmetsky/work/tmp/all_reads.fastq -sizesFile /seq/genome/mouse/mouse_Mm9/sizes -alignments /home/radon00/hmetsky/work/tmp/all_alignments.sorted.sam -numalignmentsperjob 5000
			}
			 */
	
	
			/**
			if (argmap.containsKey("sumcoverage")) {
				int numOfBases = -1;
				if (argmap.containsKey("numbases"))
					numOfBases = argmap.getInteger("numbases");
				int windowSize = -1;
				if (argmap.containsKey("windowsize"))
					windowSize = argmap.getInteger("windowsize");
	
				PolyAUtils.sumReadCoverageAndAverageReadCoveragePercentagesNearLastExonsFromAllJobs(numOfBases, windowSize);
	
				// COMMAND: java -jar ../../scripture1.jar -task polya -sizesFile /seq/genome/mouse/mouse_Mm9/sizes -transcripts /seq/rinnscratch/mgarber/mouseES/top_expressed_coding.bed -alignments /home/radon00/hmetsky/work/tmp/blat_alignment/sam_alignments/all_alignments.sorted.sam -sumcoverage -numbases 500 -windowsize 10
			} else if (argmap.containsKey("countlastexons")) {
				System.out.println(PolyAUtils.getNumOfLastExonsWithCorrectSize(scriptureTranscriptsFile));
	
				// COMMAND: java -jar ../scripture.jar -task polya -sizesFile /seq/genome/mouse/mouse_Mm9/sizes -transcripts /home/radon00/hmetsky/work/tmp/all_transcripts.bed -alignments /home/radon00/hmetsky/work/tmp/blat_alignment/sam_alignments/all_alignments.sorted.sam -countlastexons
			} else if (argmap.containsKey("exons")) {
				String exons = argmap.get("exons");
				int jobNum = argmap.getInteger("jobnum");
				int numOfBases = argmap.getInteger("numbases");
				int windowSize = argmap.getInteger("windowsize");
	
				PolyAUtils.findReadCoverageNearLastExonsByWindow(alignmentsFile, chrSizesFile, exons, jobNum, numOfBases, windowSize);
	
				// COMMAND: java -jar ../scripture.jar -task polya -sizesFile /seq/genome/mouse/mouse_Mm9/sizes -transcripts /home/radon00/hmetsky/work/tmp/all_transcripts.bed -alignments /home/radon00/hmetsky/work/tmp/blat_alignment/sam_alignments/all_alignments.sorted.sam -exons ... -jobnum ...
			} else if (argmap.containsKey("exonsperjob")) {
				int numOfLastExonsPerJob = argmap.getInteger("exonsperjob");
				int numOfBases = -1;
				if (argmap.containsKey("numbases"))
					numOfBases = argmap.getInteger("numbases");
				int windowSize = -1;
				if (argmap.containsKey("windowsize"))
					windowSize = argmap.getInteger("windowsize");
	
				PolyAUtils.findReadCoverageNearLastExonsAcrossGenomeBySplittingIntoJobs(scriptureTranscriptsFile, chrSizesFile, alignmentsFile, numOfLastExonsPerJob, numOfBases, windowSize);
	
				// COMMAND: java -jar ../../scripture1.jar -task polya -sizesFile /seq/genome/mouse/mouse_Mm9/sizes -transcripts /seq/rinnscratch/mgarber/mouseES/top_expressed_coding.bed -alignments /home/radon00/hmetsky/work/tmp/blat_alignment/sam_alignments/all_alignments.sorted.sam -exonsperjob 100 -numbases 500 -windowsize 10
			}
			 */
	
		}else if("scoreMultiple".equalsIgnoreCase(argmap.getTask())){
			DGEFullRNASeq dge = new DGEFullRNASeq(argmap);
		}else if("polyAseq".equalsIgnoreCase(argmap.getTask())){
			GenericAlignmentDataModel alignments=new GenericAlignmentDataModel("temp.bam", "sizes", 5);
			alignments.setNegativeStranded();
			ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(alignments);
			Alignments al = new Alignments("chrX","1","2600737");
			double XX = data.count(al);
			System.out.println("Negative strand: "+XX);
			alignments.setPositiveStranded();
			data = new ContinuousDataAlignmentModel(alignments);
			XX = data.count(al);
			System.out.println("Positive strand: "+XX);
			
		}
		else{
			System.err.println("Invalid task " + argmap.getTask() + "\n" + usage);
		}
	}

	public void resetTreeCache() {
		data.resetTreeCache();
	}

	private static void runScore(Map<String, Collection<RefSeqGene>> annotations, String save,
			Map<String, Integer> maskFileData, AlignmentDataModel alignments,	Map<RefSeqGene, double[]> scores) throws IOException {
		AlignmentDataModelStats alignmentData = new AlignmentDataModelStats(alignments, maskFileData);
		ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(alignmentData, maskFileData, 0, 1, null, null);
		for(String chr : annotations.keySet()) {
			logger.info("processig " + chr);
			Collection<RefSeqGene> chrAnnotations = annotations.get(chr);
			scores.putAll(data.scoreGenes(chrAnnotations, chr));
			//System.err.print("done. saving data  to " + save+"."+chr);
		}
	}

	private static int sum(Collection<Alignments> maskedRegions){
		int sum=0;
		for(Alignments align: maskedRegions){sum+=align.getSize();}
		return sum;
	}
	
	@Override  // For AlignmentCollection interface.  Should be moved to GenericDataModel eventually.  -JE
	public int size() {
		int result = 0;
		try {
			result = data.getTotalReads();
		} catch (Exception e) {
			System.err.println(e.getMessage());
			System.exit(1);
		}
		return result;
	}

	public double CV(RefSeqGene gene, IntervalTree<Alignment> tree) throws IOException{
		List<Double> vals=getNormalizedCountsPerBp(gene, tree);
		double avg=Statistics.average(vals);
		double sd=Statistics.stdev(vals);
		return sd/avg;
	}

	private double[] l2a(List<Double> list){
		double[] rtrn=new double[list.size()];
	
		int i=0;
		for(Double val: list){rtrn[i++]=val;}
	
		return rtrn;
	}

	private double[] subtract(double[] array, double factor){
		double[] rtrn=new double[array.length];
	
		for(int i=0; i<array.length; i++){
			rtrn[i]=array[i]-factor;
		}
	
		return rtrn;
	}

	private int[] group(Alignment read, Iterator<Node<RefSeqGene>> genes){
		int exonic=0;
		int intronic=0;
		int UTR5=0;
		int UTR3=0;
	
		Alignments align=new Alignments (read.getChromosome(), read.getAlignmentStart(), read.getAlignmentEnd());
	
		while(genes.hasNext()){
			RefSeqGene gene=genes.next().getValue();
			if(align.overlaps(gene.get3UtrIncludingIntrons())){UTR3++;}
			if(align.overlaps(gene.get5UtrIncludingIntrons())){UTR5++;}
			if(align.overlapsCollection(gene.getExonSet())){exonic++;}
			if(align.overlapsCollection(gene.getIntronSet())){intronic++;}
		}
	
	
		if(UTR5>0){int[] rtrn={0,0,1,0}; return rtrn;}
		if(UTR3>0){int[] rtrn={0,0,0,1}; return rtrn;}
		if(exonic>0){int[] rtrn={1,0,0,0}; return rtrn;}
		if(intronic>0){int[] rtrn={0,1,0,0}; return rtrn;}
	
		int[] rtrn={0,0,0,0};
		return rtrn;
	}

	public Collection<String> filterMultiMappers()throws IOException{
		Collection<String> rtrn=new TreeSet<String>();
		Collection<String> temp=new TreeSet<String>();
		for(String chr: data.getChromosomeLengths().keySet()){
			//System.err.println(chr);
			CloseableIterator<Alignment> iter=data.getData().getAlignmentsOverlappingRegion(new Alignments(chr, 0, this.getChromosomeLength(chr)));
			while(iter.hasNext()){
				Alignment read=iter.next();
				if(temp.contains(read.getReadName())){rtrn.add(read.getReadName());}
				else{temp.add(read.getReadName());}
			}
			iter.close();
		}
		return rtrn;
	}

	private Collection<Alignments> toExons( Alignment record){
		Collection<Alignments> rtrn=new ArrayList<Alignments>();
		//AlignmentBlock [] blocks = record.getAlignmentBlocks();
		List<Alignments> blocks = getAlignmentBlocks(record);
		char [] gapTypes =  record.getGapTypes();
		
		if(blocks==null || blocks.size() <= 1){
			rtrn.add(new Alignments(record.getChromosome(), record.getAlignmentStart(), record.getAlignmentEnd())); 
			return rtrn;
		}
		
		
		Alignments lastExon = new Alignments(record.getChromosome(), blocks.get(0).getStart(), blocks.get(0).getEnd());;
		for(int i=0; i<blocks.size(); i++){
			if(i == blocks.size() - 1 || gapTypes[i] == SamAlignment.SKIPPED_REGION) {
				Alignments finalExon = new Alignments(record.getChromosome(), lastExon.getStart(), blocks.get(i).getEnd());
				rtrn.add(finalExon);
				if(i < blocks.size()-1) {
					lastExon = new Alignments(record.getChromosome(), blocks.get(i+1).getStart(), blocks.get(i+1).getEnd());
				}
			}
		}
	
		return rtrn;
	}

	private Collection<Alignments> dedup(IntervalTree<Alignments> allSegments){
		Set<Alignments> set=new TreeSet<Alignments>();
	
		Iterator<Node<Alignments>> iter=allSegments.iterator();
	
		while(iter.hasNext()){
			Alignments align=iter.next().getValue();
			set.add(align);
		}	
	
		return set;
	}

	private Collection<Alignments> findMaxContiguous(Collection<Alignments> windows) throws IOException{
		Set<Alignments> rtrn = new TreeSet<Alignments>();
	
		//loop through each alignment and compute number of reads for the ends
		Map<Alignments, double[]> counts=this.getDataForAlignments(windows);
		for(Alignments align: windows){
			Alignments trunc=collapseMaxContiguous(align, counts.get(align));
			if(trunc!=null){rtrn.add(trunc);}
		}
	
		return rtrn;
	}

	private Collection<Alignments> splitAlignments(Alignments exon, Map<String, IntervalTree<Alignments>> trees){
		Collection<Alignments> rtrn=new TreeSet<Alignments>();
		Iterator<Node<Alignments>> overlappers=trees.get(exon.getChr()).overlappers(exon.getStart(), exon.getEnd());
		//then split exon
		Collection<Alignments> splitExons=SplitByPairedEnd.splitExon(exon, overlappers);
		overlappers=trees.get(exon.getChr()).overlappers(exon.getStart(), exon.getEnd());
		splitExons.addAll(SplitByPairedEnd.splitExon1(exon,overlappers));
		//add split exons to exon collection
		rtrn.addAll(splitExons);
		return rtrn;
	}

	public Collection<Alignments> permuteRegions(Alignments align) throws IOException{
		Collection<Alignments> rtrn=new TreeSet<Alignments>();
		//get records overlapping this region
		Iterator<Alignments> iter=data.getExonAlignmentsOverlappingRegion(align);
	
		int regionSize=align.getSize();
	
		while(iter.hasNext()){
			Alignments record=iter.next();
			//randomly place in the interval
			int start=align.getStart()+new Double(Math.random()*(regionSize)).intValue();
			Alignments r=new Alignments(align.getChr(), start, start+(record.getSize()));
			r.setStrand(record.getStrand());
			rtrn.add(r);
		}
	
		return rtrn;
	}

	public Collection<Alignments> rawRegions(Alignments align) throws IOException{
		Collection<Alignments>  rtrn=new TreeSet<Alignments> ();
		//get records overlapping this region
		Iterator<Alignments> iter=data.getExonAlignmentsOverlappingRegion(align);
	
		while(iter.hasNext()){
			Alignments record=iter.next();
			rtrn.add(record);
		}
	
		return rtrn;
	}

	private static Collection<RefSeqGene> filter(Map<RefSeqGene, double[]> scores, double alpha2) {
		Collection<RefSeqGene> rtrn=new TreeSet<RefSeqGene>();
	
		for(RefSeqGene gene: scores.keySet()){
			double[] p=scores.get(gene);
			if(p[0]<alpha2){
				rtrn.add(gene);
			}
		}
	
		return rtrn;
	}

	public static Collection<RefSeqGene> filterPaths(Map<Path, double[]> scores, double alpha2) {
		Collection<RefSeqGene> rtrn=new TreeSet<RefSeqGene>();
	
		/*
		 * [0] = P-value
		 * [1] = enrichment 
		 * [2] = score (counts?)
		 * [3] = average coverage
		 * [4] = rpkm
		 * [5] = lambda
		 * [6] = size
		 * [7] = localPvalue
		 */
		for(Path gene: scores.keySet()){
			double[] p=scores.get(gene);
			
			//Gene g=gene.toGene();
			RefSeqGene g=gene.toGene(p);
			//System.out.println("For "+g.toBED()+": ");
			//System.out.println("P-value: "+p[0] +", RPKM: "+p[4]+", Score: "+p[2]);
			g.setBedScore(p[4]);
			if(g.getNumExons()>1){
				if(p[0]<alpha2){rtrn.add(g);}
			}
			else if(p.length>7){
	
				if(p[0]<alpha2 && p[7]<alpha2){rtrn.add(g);}
			}
			else{
	
				if(p[0]<alpha2){rtrn.add(g);}
			}
	
		}
	
		return rtrn;
	}
	
	// For AlignmentCollection interface - should be moved to Generic model eventually. JE
	@Override
	public int getBasesCovered(Alignments region) throws IOException {
		return data.getData().getBasesCovered(region, extensionFactor); 
	}
	
	// For AlignmentCollection interface - should be moved to Generic model eventually. JE
	@Override
	public int getBasesCovered(Alignments region, int EF) throws IOException {
		return data.getData().getBasesCovered(region, EF); 
	}

	private static Collection<Path> refSeqGeneToPath( IntervalTree<RefSeqGeneWithIsoforms> chrTree, ChromosomeWithBubblesJGraphT graph) {
		Collection<Path> rtrn = new LinkedList<Path>();
		Iterator<RefSeqGeneWithIsoforms> it =chrTree.valueIterator();
		while(it.hasNext()){
			Collection <RefSeqGene> allIso=it.next().getAllIsoforms();
			for (RefSeqGene iso: allIso){
				Path p=new Path(iso, graph);
				rtrn.add(p);
			}
		}
		return rtrn;
	}

	public static Map<String, Integer> parseMaskFiles(File[] files) throws IOException{
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		if(files==null){return rtrn;}
	
		for(int i=0; i<files.length; i++){
			String chr=files[i].getName().split("\\.")[0];
			Collection<Alignments> data = BEDFileParser.loadAlignmentData(files[i]);
			rtrn.put(chr,sum(data));
		}	
		return rtrn;
	}

	private Map<Path, double[]> contiguousAcrossGraph(String chr,	Sequence chrSequence, int start, int end, double minimumSpliceFrequencyToFollowEdge) throws IOException {
		//Segment and make transcript graph
		ChromosomeWithBubblesJGraphT graph=makeGraphWithCounts(chr,this.minNumberOfSplices, chrSequence,  minimumSpliceFrequencyToFollowEdge);
		graph.writeGraph(chr+"graph.dot");
		logger.info("Done making spliced graph");
		Map<Path, double[]> pathScores=augmentAndScoreTranscriptGraph(graph, chr,  minimumSpliceFrequencyToFollowEdge);
		return pathScores;
	}

	private Map<Path, double[]> contiguousAcrossGraphStrandSpecific(String chr,	Sequence chrSequence, int start, int end, double minimumSpliceFrequencyToFollowEdge) throws IOException {
		//Segment and make transcript graph
		ChromosomeWithBubblesJGraphT graphPlus=makeGraphWithCounts(chr,0, chrSequence,  minimumSpliceFrequencyToFollowEdge,"+");
		graphPlus.writeGraph(chr+"graph.plus.dot");
		ChromosomeWithBubblesJGraphT graphMinus=makeGraphWithCounts(chr,0, chrSequence,  minimumSpliceFrequencyToFollowEdge,"-");
		graphPlus.writeGraph(chr+"graph.minus.dot");
		logger.info("Done making spliced graph");
		Map<Path, double[]> pathScores=augmentAndScoreTranscriptGraph(graphPlus, chr,  minimumSpliceFrequencyToFollowEdge);
		pathScores.putAll(augmentAndScoreTranscriptGraph(graphMinus, chr,  minimumSpliceFrequencyToFollowEdge));
		return pathScores;
	}
	
	private Map<Path, double[]> augmentAndScoreTranscriptGraph(ChromosomeWithBubblesJGraphT graph, String chr, double minimumSpliceFrequencyToFollowEdge) throws IOException{

		logger.info("Getting all paths minumum frequency " + minimumSpliceFrequencyToFollowEdge);
		/*
		 * Get all paths 
		 */
		Collection<Path> pairedEndPaths= graph.getAllPaths(minimumSpliceFrequencyToFollowEdge);
		logger.info("Done getting paths. Total: " + pairedEndPaths.size());

		//Get local rate and build insert distribution
		//TODO Local segment should use path scores as well
		localSegment(pairedEndPaths, alpha, true);
		logger.info("Done with local segmentation");
		graph.setLocalRate(pairedEndPaths);
		logger.info("Done setting local rate in graph");

		//write the graph and all its info to a file
		//if(saveToFile!=null){
		//	graph.writeGraph(saveToFile+"."+chr+".dot");
		//}

		//score each path and get significance
		Map<Path, double[]> scores=scorePaths(pairedEndPaths, data.getLambda(chr));
		logger.info("Done scoring paths.");
		//return paths and scores
		return scores;
	}




	private ChromosomeWithBubblesJGraphT constructGraphs(String chr, Collection<Alignments> exons, Collection<Alignments> introns, int minDistance, double minSpliceFreq) throws IOException{
		//Add a collapse step up front
		exons=CollapseByIntersection.collapseByIntersection(exons, false);
		logger.info("Finished collapsing");

		//First, decollapse exons
		/*Collection<Alignments> decollapsedPlus=CollapseByIntersection.DecollapseByIntronLocation(exons, introns, "+");
		Collection<Alignments> decollapsedMinus=CollapseByIntersection.DecollapseByIntronLocation(exons, introns, "-");

		Collection<Alignments> decollapsed=new TreeSet();
		decollapsed.addAll(decollapsedPlus);
		decollapsed.addAll(decollapsedMinus);
		 */

		Collection<Alignments> decollapsed=CollapseByIntersection.DecollapseByIntronLocation(exons, introns);
		logger.info("Finished decollapse");

		ChromosomeWithBubblesJGraphT graph=new ChromosomeWithBubblesJGraphT(chr, decollapsed, introns, data.getLambda(chr), data.getNumberMarkers(chr), data.getNumberOfReads(chr), minNumberOfSplices);
		logger.info("Made graph");

		//Second, get exons that were fully within introns
		//Also split up covered exons into pieces and compute significance
		Collection<Alignments> orphans=getExtendedExons(exons, decollapsed, graph, minDistance, minSpliceFreq);
		logger.info("Got extended pieces");

		//Add new nodes to graph
		decollapsed.addAll(orphans);
		graph=new ChromosomeWithBubblesJGraphT(chr, decollapsed, introns, data.getLambda(chr), data.getNumberMarkers(chr), data.getNumberOfReads(chr), 0);
		logger.info("Made graph");

		data.resetTreeCache(chr);

		return graph;
	}

	

	

	/*
	private static Collection<EdgeSourceType> typesToFollow(boolean followPairedEnds) {
		ArrayList<EdgeSourceType> types = new ArrayList<EdgeSourceType>();
		types.add(EdgeSourceType.SPLICED);
		if(followPairedEnds) {
			types.add(EdgeSourceType.PAIRED);
		}
		return types;
	}
	 */

	public Map<Path, double[]> scanFromGraph(int windowSize, double alpha, boolean trimEnds, String chr, int start, int end, Sequence chrSequence, boolean filterSplice, String save, double minimumSpliceFrequencyToFollowEdge) throws IOException{
		ChromosomeWithBubblesJGraphT graph=new ChromosomeWithBubblesJGraphT(chr, start, end);
		graph.loadBubbles(data, chrSequence, minIntronSize , minNumberOfSplices, filterSplice);
		return scanFromGraph(graph, windowSize, alpha, trimEnds, chr, start, end, chrSequence, filterSplice, save, minimumSpliceFrequencyToFollowEdge);
	}

	public Map<Path, double[]> scanFromGraph(ChromosomeWithBubblesJGraphT graph, int windowSize, double alpha, boolean trimEnds, String chr, int start, int end, Sequence chrSequence, boolean filterSplice, String save, double minimumSpliceFrequencyToFollowEdge) throws IOException{
		this.alpha=alpha;

		double T=getNumberMarkers(chr);
		double lambda = getLambda(chr);
		int criticalValue=calculateCriticalValue(new Double(T).intValue(), windowSize, T, alpha, lambda);

		Map<Path, double[]> pathScores = scanGenome(windowSize, criticalValue, chr, start, end, graph, chrSequence, save, minimumSpliceFrequencyToFollowEdge);

		return pathScores; 
	}
	
	private Map<Path, double[]> scanGenome(int fixedWidth, int critVal, String chr, int start, int end, ChromosomeWithBubblesJGraphT graph, Sequence chrSeq, String saveToFile, double minimumSpliceFrequencyToFollowEdge)throws IOException{

		//TODO Ensure that transcript graph contains all node and edge counts
		ChromosomeWithBubblesJGraphT transcriptGraph=acrossGraph (critVal, fixedWidth, graph, chr, chrSeq, minimumSpliceFrequencyToFollowEdge);

		Map<Path, double[]> pathScores=this.augmentAndScoreTranscriptGraph(transcriptGraph, chr, start, end, saveToFile, minimumSpliceFrequencyToFollowEdge);

		return pathScores;
	}
	
	private Map<Path, double[]> augmentAndScoreTranscriptGraph(ChromosomeWithBubblesJGraphT graph, String chr, int start, int end, String saveToFile, double minimumSpliceFrequencyToFollowEdge) throws IOException{

		//Add paired end edges to the graph
		if(pairedData != null) {
			//Collection<Path> currentPaths = graph.getAllPaths();
			//Map<Path, double[]> scoresWithOnlyRPKM = scorePathsForOnlyRPKM(currentPaths);
			//Map<Path, RefSeqGene> currentGenes = pathsToGenesMap(scoresWithOnlyRPKM);
			
			AddPairedEndEdges.addPairedEndEdgesToPaths(graph, pairedData, chr, start, end, saveToFile);
			//AddEdgesUsingPairedEnds.usePairedEndsToAddEdgesToPaths(graph, this, currentPaths, currentGenes, pairedData, chr, start, end, saveToFile);
			logger.info("Done adding paired ends (if available)");				
		}
		logger.info("Getting all paths minumum frequency " + minimumSpliceFrequencyToFollowEdge);
		//Get all paths
		Collection<Path> pairedEndPaths= graph.getAllPaths(minimumSpliceFrequencyToFollowEdge);
		logger.info("Done getting paths. Total: " + pairedEndPaths.size());

		//Get local rate and build insert distribution
		//TODO Local segment should use path scores as well
		localSegment(pairedEndPaths, alpha, true);
		logger.info("Done with local segmentation");
		graph.setLocalRate(pairedEndPaths);
		logger.info("Done setting local rate in graph");

		//write the graph and all its info to a file
		if(saveToFile!=null){
			graph.writeGraph(saveToFile+".dot");
		}

		//score each path and get significance
		Map<Path, double[]> scores=scorePaths(pairedEndPaths, data.getLambda(chr));
		logger.info("Done scoring paths.");
		//return paths and scores
		return scores;
	}
	
	
	public  BEDFileParser segmentChromosome(double minSpliceFrequency,
			int[] fixedWidth, String chrToSegment,  boolean trimEnds,
			boolean filterCanonical, String sequenceFile, double alpha,
			int chrStart, int chrEnd, double trimQuantile) throws IOException, IllegalArgumentException, MathException {

		double memoryPercent=Runtime.getRuntime().freeMemory()/(double)Runtime.getRuntime().totalMemory();
		System.err.println("At the start of segmenting, Total memory: "+ Runtime.getRuntime().totalMemory()+ "Free memory: " +Runtime.getRuntime().freeMemory()+ "% memory available  "+(memoryPercent*100));
		BEDFileParser result = new BEDFileParser();
		Sequence chrSequence = null;
		List<String> chromosomes = new ArrayList<String>();
		//IF NO CHROMOSOME NAME IS PROVIDED, ADD ALL CHROMOSOMES TO LIST
		if(chrToSegment == null || chrToSegment.trim().length() == 0) {
			chromosomes.addAll(data.getChromosomeLengths().keySet());
		} 
		//ELSE, ADD THE SPECIFIED CHROMOSOME TO LIST
		else {
			chromosomes.add(chrToSegment);
		}

		//FOR EACH CHROMOSOME UNDER CONSIDERATION
		for(String chr : chromosomes) {
			logger.info("Processing chromosome " + chr + " sequence file " + sequenceFile);
			// CHECK THAT GENOME SEQUENCE FILE IS SUPPLIED
			   //EXTRACT RECORDS FOR THIS CHROMOSOME FROM SEQUENCE FILE
			if(sequenceFile != null) {
				FastaSequenceIO fsio = new FastaSequenceIO(sequenceFile);
				List<String> chrIds = new ArrayList<String>(1);
				chrIds.add(chr);
				List<Sequence> seqs = fsio.extractRecords(chrIds);
				if(seqs.isEmpty() && !chr.contains("chr")) {
					logger.info("Reference sequence for chromosome " + chr + " was not found will try and find sequence for chr"+chr +" in genome reference fasta file");
					chrIds.clear();
					chrIds.add("chr"+chr); //Since many fasta refrence files name chromosomes with chr, add this id to maximize prob to find the ref sequence
					seqs = fsio.extractRecords(chrIds);
				}

				if(!seqs.isEmpty()) {
					chrSequence = seqs.get(0);
					logger.info("Loaded chromosome Sequence with ID: " + chrSequence.getId()  );
				} else {
					logger.warn("Sequence for " + chr + " was not found in file " + sequenceFile + " skipping chromosome");
					continue;
					//throw new RuntimeException("Sequence for " + chr + " was not found in file " + sequenceFile + " continuing without splice site analisys");  //TODO: handle more gracefully.
				}
			}

			memoryPercent=Runtime.getRuntime().freeMemory()/(double)Runtime.getRuntime().totalMemory();
			System.err.println("After loading sequence, Total memory: "+ Runtime.getRuntime().totalMemory()+ "Free memory: " +Runtime.getRuntime().freeMemory()+ "% memory available  "+(memoryPercent*100));
			
			Collection<RefSeqGene> genes;

			if(fixedWidth[0]==0) {
				logger.info("Segmenting across graph");
				Map<Path, double[]> scores;
				scores = contiguousAcrossGraph(chr, chrSequence, chrStart, chrEnd, minSpliceFrequency);
				//scores = contiguousAcrossGraphStrandSpecific(chr, chrSequence, chrStart, chrEnd, minSpliceFrequency);
				//writeFullBED("allScores.bed", scores);
				genes=filterPaths(scores, alpha);
			} else {
				logger.info("Segmenting using windows");
				Map<Path, double[]> scores;
				scores=scanFromGraph(fixedWidth[0], alpha, trimEnds, chr, chrStart, chrEnd, chrSequence, filterCanonical, minSpliceFrequency);
				//writeFullBED("allScores.bed", scores);
				genes=filterPaths(scores, alpha);
			}

			memoryPercent=Runtime.getRuntime().freeMemory()/(double)Runtime.getRuntime().totalMemory();
			System.err.println("After segmenting across graph, Total memory: "+ Runtime.getRuntime().totalMemory()+ "Free memory: " +Runtime.getRuntime().freeMemory()+ "% memory available  "+(memoryPercent*100));
			
			if(data.isPaired()) {
				String unpairedData = "run_" + chr+"_"+ System.currentTimeMillis() + ".unpaired.reconstruction.bed"; 
				writeFullBED(unpairedData, genes);

				//String unpairedData = "run_" + chr+"_"+ System.currentTimeMillis() + ".unpaired.reconstruction.bed"; 
				//writeFullBED(unpairedData, genes);
				logger.info("PAIRED-END Stem\nData is paired ended. Computing insert size distribution");
				PairedEndDistribution ped = computeInsertSizeDistribution(genes);
				ped.ensureNonZeroCounts();
				//ped.getSizeDistribution().write(save + ".insertsize.dist");

				logger.info("Finished computing insert size distribution. Mean insert size: " +ped.getAccurateSizeDistribution().getMean() + " sd: "+ ped.getAccurateSizeDistribution().getStandardDeviation() + " Revisitting annotations for their compatibility with inserts ");

				//BEDFileParser updatedReconstructions = data.updateReconstructionsWithPairs(genes, ped);
				//BEDFileParser mergedAnnotations = data.mergeReconstructionWithPairs(updatedAnnotations, ped);
				BEDFileParser mergedAnnotations = ped.mergeReconstructionWithPairs(genes);
				
				logger.info("Generating data for connecting last exons to single exons that follow by using paired ends");
				//ConnectLastExonToSingleExonsWithPairedEnds cletsewpe = new ConnectLastExonToSingleExonsWithPairedEnds(chr, genes, data, 1000, 0.1, 0.9, 200, 0.02, 0.9);
				//cletsewpe.writePairedInsertDistributionToFile(save+".insert.size.dist");
				logger.info("Connecting last exons to single exons that follow by using paired ends");
				//genesWithConnectionsToSingleExons = cletsewpe.connectLastExonToSingleExonsWithPairedEnds();
				
				if(trimEnds) {
					logger.info("Trimming reconstructions");
					mergedAnnotations = trimEndsForGenes(mergedAnnotations, trimQuantile);
				}
				logger.info("rescoring merged (and trimmed) reconstructions");
				scoreGenes(mergedAnnotations);
				
				result.addRefSeqSet(mergedAnnotations.GetGenes());
				
			} /*else {
				//writeFullBED(save, genes);
				if(trimEnds) {
					logger.info("Trimming reconstructions ");
					genes = trimEndsForGenes(genes, trimQuantile);
					logger.info("Rescoring trimmed genes");
					scoreGenes(genes); 
				}
				result.addRefSeqSet(genes);
			}*/
			memoryPercent=Runtime.getRuntime().freeMemory()/(double)Runtime.getRuntime().totalMemory();
			System.err.println("Paired End done, Total memory: "+ Runtime.getRuntime().totalMemory()+ "Free memory: " +Runtime.getRuntime().freeMemory()+ "% memory available  "+(memoryPercent*100));
			
			if(trimEnds) {
				logger.info("Trimming reconstructions ");
				genes = trimEndsForGenes(genes, trimQuantile);
				logger.info("Rescoring trimmed genes");
				scoreGenes(genes); 
				memoryPercent=Runtime.getRuntime().freeMemory()/(double)Runtime.getRuntime().totalMemory();
				System.err.println("After trimming, Total memory: "+ Runtime.getRuntime().totalMemory()+ "Free memory: " +Runtime.getRuntime().freeMemory()+ "% memory available  "+(memoryPercent*100));
			}
			
			result.addRefSeqSet(genes);
		}
		
		return result;
	}



	//static String usage="version=v6 \n args[0]=BAM File \n args[1]=mask file \n args[2]=save file \n args[3]=sizes \n args[4]=window size \n args[5]=trim ends?";
	static String usage="\nParameters \n -alignment <Alignment file in BAM, SAM or Alignemnt format> \n -maskFileDir <Mask File directory> \n -out <Output file name>"+ 
			"\n -sizeFile <Chromosome size file> \n -chr <Chromsomosome to segment> \n -chrSequence <Necessary to filter spliced reads by splice site information. Notice that this is only compatible with region files that contain regions of only one chromosome> "+
			"\n Optional arguments: \n -windows <Comma separated list of windows to evaluate defaults to contiguous regions of coverage>\n -trim <Include this flag if trimming of the ends of windows based on read coverage  is desired this is expensive> \n -alpha <Desired FDR>" +
			"\n  -dontFilterCanonicalSplice" +
			"\n -start <To segment only a subregion of the chromosome include its start> -end <To segment only a subregion of the chromosome include its end> " + 
			"\n -minSpliceSupport <Minimum count to support splice reads, default is 1> \n-minSpliceFrequency <When there are more than one splice junction, junctions that account for less than the specified portion of junctions are ignored>\n -pairedEnd <Paired end alignment files> -strandSpecificReads <Strand specific alignment file>\n\t-scoreRegions <Full BED to score> -upWeightSplices -lambda <If a prior background expectation for number of reads per base exists> -exons <BED file of exons> -introns <Introns and counts>" +
			"\n\nTask: AddPairs -  Uses a paired end alignment to tune graph \n\t-in <Graph in .dot format. Standard input is assumed> \n\t-pairedEnd <Paired end information (as in previous task), in single line BED format>\n\t -maskFileDir <Directory containing mask files for the genome> \n\t-chr <Chromosome (only a chromosome at a time is supported at this point)> \n\t-sizeFile <Chromosome size file> \n\t-out <Output file name>"+
			"\n\nTask: fastScore -  Computes several expression related scores for a set of annotations using a the graph .dot file \n\t-in <chr.dot file> \n\t-annotations <BED file with annotation to score>\n\t -chr <chr e.g: chrZ>\n\t -alpha <optional>\n\t -out\n" +
			"\n\nTask: score -  Computes several expression related scores for a set of annotations -in <Full BED file with annotations to score> \n\t-alignment <Alignment file in BAM, SAM or Alignemnt format> \n\t-sizeFile <Chromosome size file> \n\t-out <Output file name> \n\t -maskFileDir <Mask File directory>\n\t -useConstituentExons <For each gene a set of contituent exons will be chosen and scored>"+
			"\n\nTask: extractDot - Extracts a graph for the specified region \n\t-in <Dot file from a previous Scripture run \n\t-chr <Chromosome> \n\t-start <Start of region> \n\t-end<End of region> \n\t-out <output file>"+
			"\n\nTask: getIdenticalGappedReadsTranscripts -  Report all transcripts that ALL their introns are spanned by identical gapped reads -in <Full BED file with annotations to score> \n\t-alignment <Alignment file in BAM, SAM or Alignemnt format> \n\t-sizeFile <Chromosome size file> \n\t-out <Output file name> \n\t -maskFileDir <Mask File directory>"+
			"\n\nTask: makePairedFile Makes a paired end alignment file from two sets of independtly aligned left and right ends, ideally the files should be name-sorted \n\t-pair1 <First pair alignments> \n\t-pair2 <Second pair alignments> \n\t-out <output consolidated paired end alignment> \n\t-sorted  <Include this flag if the data is already read name sorted, ideally both input files should be sorted by read name using unix sort for example> \n\t-usePair2Orientation <If the second paired rather than the first should be used to orient insert like for dUTP libraries>\n\t-forChIP <If the alignment if for ChIP rather than RNAseq then Ms will be used instead of Ns>" + 
			"\n\nTask: chipScan - Segment the genome assuming contiguous data. Similar to the default task but optimized for contiguous data. \n -alignment <Alignment file in BAM, SAM or Alignemnt format> \n -extensionFactor <Extend reads by this factor (defaults to 0)> \n -maskFileDir <Mask File directory> \n -out <Output file name>"+ 
			"\n -chr <Chromosome to segment>\n -sizeFile <Chromosome size file> \n -windows <Comma separated list of windows to evaluate defaults to contiguous regions of coverage> \n Optional arguments:\n -findMaxContiguous <Each significant window is trimmed by finding the contiguous sub region with coverage over a predefined threshold> -trim <Include this flag if trimming of the ends of windows based on read coverage  is desired this is expensive> \n -alpha <Desired FDR>" +
			"\n\nTask: trim -  Trims end of transcripts by removing all bases whose coverage is below the specified quantile of transcript expression -in <Full BED file with annotations to trim> \n\t-alignment <Alignment file in BAM, SAM or Alignemnt format> \n\t-sizeFile <Chromosome size file> \n\t-out <Output file name> \n\t -maskFileDir <Mask File directory>\n\t-quantile <Coverage quantile below which end bases should be trimmed>"+
			"\n\nTask: trimSegments -  Trims end of continuous segments by removing all bases whose coverage is below the specified quantile of transcript expression -in <Full BED file with annotations to trim> \n\t-alignment <Alignment file in BAM, SAM or Alignemnt format> \n\t[-sizeFile <Chromosome size file>] \n\t-out <Output file name> \n\t -maskFileDir <Mask File directory>\n\t-quantile <Coverage quantile below which end bases should be trimmed> \n\t\t-findMaxContiguous <To break up segments that have peak/valley shapes>"+
			"\n\nTask: adjustEnds -  Takes an annotation set and adjust transcript ends based on the given alignment \n\t-in <Full BED file with annotations to trim> \n\t-alignment <Alignment file in BAM, SAM or Alignemnt format> \n\t[-sizeFile <Chromosome size file>] \n\t-out <Output file name> \n\t -maskFileDir <Mask File directory>\n\t=trimQuantile <Coverage quantile below which end bases should be trimmed> \n\t\t-findMaxContiguous <To break up segments that have peak/valley shapes>"+
			"\n";



}
