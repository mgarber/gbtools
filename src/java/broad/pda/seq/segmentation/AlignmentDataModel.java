package broad.pda.seq.segmentation;

import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;

import org.broad.igv.sam.Alignment;

import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.sequence.Sequence;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import net.sf.samtools.util.CloseableIterator;

//TODO Add updateInterval tree to interface

public interface AlignmentDataModel{

	//Get all alignments overlapping a specified region
	public CloseableIterator<Alignment> getAlignmentsOverlappingRegion(Alignments align) throws IOException;
	
	//Get all represented chromosomes and their lengths in the genome
	public Map<String, Integer> getChromosomeLengths();

	/**
	 * Get length of a chromosome
	 * @param chr
	 * @return
	 */
	public int getChromosomeLength(String chr);

	/**
	 * Counts the number of reads for the specified chromosome
	 * @return
	 * @throws IOException 
	 */
	public double getCounts(String chr) throws IOException;
	
	//get the number of alignments overlapping a given region
	public double getCountsPerAlignment(LightweightGenomicAnnotation first, int EF) throws IOException;

	//get the number of alignments overlapping a given region with a cached interval tree
	public double getCountsPerAlignment(Alignments align, IntervalTree<Alignment> tree, int EF);
	
	//get the number of alignments overlapping a given region with a cached interval tree
	public double getCountsPerAlignment(RefSeqGene gene, IntervalTree<Alignment> tree, int EF);

	//get the number of alignments overlapping a given region with a cached interval tree
	public double getCountsPerAlignmentStranded(RefSeqGene gene, IntervalTree<Alignment> tree, int EF, String orientation);

	//given an array of alignments generate the counts over the whole region
	public double getCountsPerAlignment(Alignments[] alignments, IntervalTree<Alignment> tree, int EF);
	
	//get the number of reads that are fully contained within the region (doesn't start before or extend past)
	public double getCountsPerAlignmentFullyContained(RefSeqGene gene, IntervalTree<Alignment> tree, int EF);
	
	/**
	 * Similar to @see getCountsPerAlignment(Alignments[] alignments, IntervalTree<Alignment> tree, int EF)
	 * but delegates the management of the alignment interval tree to the implementing class rather than the caller
	 * @param alignments
	 * @param EF
	 * @return
	 * @throws IOException 
	 */
	public double getCountsPerAlignment(Alignments[] alignments,  int EF) throws IOException;

	public double getCountsPerAlignment(Alignments align,
	Map<String, IntervalTree<Alignments>> goodExonTree,
	IntervalTree<Alignment> tree, int extensionFactor);

	public int getCountsWithinExons(Alignments align, Iterator<Alignment> iter, int EF);

	public int getCountsWithinExons(Alignments align, IntervalTree<Alignment> tree, int EF);
	
	public double getCountsOfUniqueOverlappers(Alignments target, Alignments exclude, IntervalTree<Alignment> tree, int EF);
	
	public double getCountsOfUniqueOverlappers(Alignments target, RefSeqGene exclude, IntervalTree<Alignment> tree, int EF);
	
	public double getCountsOnBase(String chr, int index) throws IOException;
	
	public int getBasesCovered(Alignments record, int EF) throws IOException;
	
	//get the number of spliced reads overlapping a given region
	public double getNumberOfSplicedReads(Alignments align) throws IOException;
	
	public GeneCounts getGeneCounts(RefSeqGene gene, int extensionFactor) throws IOException;

	//get the actual number of bases covered by alignments overlapping a given region with a cached interval tree (if a single 76 base read spans the gene, return 76)
	public double getBasesCoveredPerAlignment(RefSeqGene gene, IntervalTree<Alignment> tree, int EF);
			
	//get the exonic regions overlapping a given region
	public Iterator<Alignments> getExonAlignmentsOverlappingRegion(Alignments align) throws IOException;
	
	public RefSeqGene getPeak(RefSeqGene gene, RefSeqGene startCodon, IntervalTree<Alignment> tree, int EF);
	
	public Collection<Alignments> getOverlappingRegions(String chr) throws IOException;
	
	//compute a cached interval tree for a given region
	//Replace with counter per node
	public IntervalTree<Alignment> getIntervalTree(String chr, int start, int end) throws IOException;
	//public IntervalTree<Alignments> getIntervalTreeCachedAlignments(String chr, int start, int end);
	
	
	public IntervalTree<Alignment> getIntervalTreeCached(String chr, int start, int end) throws IOException;

	public RefSeqGene updateGeneByFirstCounts(RefSeqGene gene, IntervalTree<Alignment> tree, int EF);
	
	//get all spliced reads overlapping a region
	//break into intronic regions and store them as introns
	public Map<Alignments, Integer> getSplicedReads(Alignments align, int minIntronSize, int maxIntronSize) throws IOException;
	
	public Map<Alignments, Integer> getSplicedReads(Alignments region, Collection<ReadFilter> filters) throws IOException;
	public Map<Alignments, Integer> getSplicedReads(Alignments region) throws IOException;
	public Map<Alignments, Integer> getSplicedReads(Alignments align, final int minIntronSize, final int maxIntronSize , int minNumIntrons, Sequence chrSeq) throws IOException;
	
	public Map<Alignments, Integer> getSplicedReads(Alignments region, Collection<ReadFilter> filters, int minNumIntrons) throws IOException;

	public Map<Alignments, Integer> getSplicedReads(Alignments region, Collection<ReadFilter> filters, int minNumIntrons, Sequence chrSeq) throws IOException;

	public Map<Alignments, Integer> getSplicedReads(Alignments region, ReadFilter filter) throws IOException;

	public Map<Alignments, Integer> getSplicedReads(int minIntronSize, int maxIntronSize) throws IOException;

	public double getSpliceWeightFactor();

	public Collection<? extends Alignments> getSplicedReadExons(String chr) throws IOException;

	public CloseableIterator<Alignment> getReadIterator();
	public CloseableIterator<Alignment> getReadIterator(Alignments region) throws IOException;
	
	
	
	//public int getFirstReads(String chr);
	
	/**
	 * @param the chromosome over which to iterate
	 * @return CloseableIterator over all reads mapped to the given Chromosome
	 */
	public CloseableIterator<Alignment> getChromosomeReadIterator(String chromosome) throws IOException;

	public double getTotalNumberOfStrandedReads();

	public double getScorePerAlignmentFromCache(RefSeqGene window, IntervalTree<Alignment> chunkAlignmentTree, int extensionFactor);
	public void resetTreeCache();

	public int getChunkStart();

	public int getChunkEnd();
	
	public IntervalTree<Alignment> getIntervalTreeTruncatedCached(String chromosome, int start, int end) throws IOException;

	public IntervalTree<Alignments> getFullIntervalTreeAsAlignments(String chr) throws IOException;

	public Alignments getNextExon(String chr) throws IOException;

	public double getMinimumMappingQuality();

	public void setChunkSize(int chunkSize);

	public void setNegativeStranded() ;

	public void setPositiveStranded() ;
	
	public void setSecondRead() ;
	
	public void setFirstRead() ;

	void setNormalizationStrategy(ReadCountNormalizer normalizer);

	void setMinimumMappingQuality(double minimumMapQuality);

	//get the number of alignments overlapping a given region with a cached interval tree
	public boolean passes(RefSeqGene gene, IntervalTree<Alignment> tree, int EF, double cutoff);

	//get the number of alignments overlapping a given region with a cached interval tree
	public long passesCutoff(Alignments align, IntervalTree<Alignment> tree, int EF, double threshold);

	public boolean hasNextExon(String chr) throws IOException;
	public void restartIterator();

	public void resetTreeCache(String chr);
	
	public void clearFullIntervalTreeAsAlignmentsCached(String chr);
	
	public boolean isStranded();
	
	public void unsetStranded();
	public boolean isPositiveStranded();
	public boolean isNegativeStranded() ;
	public String getModelFilePath();
	void setExtensionFactor(int factor);
}
