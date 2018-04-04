/**
 * Data processing and analysis tools for sequencing of crosslinked samples
 */
package broad.pda.seq.protection;

import java.io.IOException;
import java.util.Collection;

import org.broad.igv.sam.Alignment;

import broad.core.sequence.SequenceUtils;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.AlignmentDataModelStats;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;
import net.sf.samtools.util.CloseableIterator;

/**
 * @author prussell
 * Analyze sequencing read depth by region
 */
public final class SampleCounts {

	private GenericAlignmentDataModel genericAlignmentDataModel;
	private AlignmentDataModelStats alignmentDataModelStats;
	private ContinuousDataAlignmentModel alignmentData;
	
	/**
	 * Construct object and store alignment data
	 * @param data Alignment data
	 */
	public SampleCounts(ContinuousDataAlignmentModel data) {
		
		this.alignmentData = data;
		this.alignmentDataModelStats = this.alignmentData.getData();
		this.genericAlignmentDataModel = (GenericAlignmentDataModel) this.alignmentDataModelStats.getData();

	}
	
	/**
	 * Construct object and store alignment data from bam file
	 * @param alignmentFile Bam alignment file
	 * @param sizeFile Table of chromosome sizes
	 * @throws IllegalArgumentException
	 * @throws IOException
	 */
	public SampleCounts(String alignmentFile, String sizeFile) throws IllegalArgumentException, IOException {
		
		this.alignmentData = SequenceUtils.getDataModel(alignmentFile, sizeFile, false);
		this.alignmentDataModelStats = this.alignmentData.getData();
		this.genericAlignmentDataModel = (GenericAlignmentDataModel) this.alignmentDataModelStats.getData();
		
	}
	
	/**
	 * Get total read count in region
	 * @param alignments The region
	 * @return The number of reads mapping to the region
	 * @throws IOException
	 */
	private int getCount(Alignments alignments) throws IOException{		
		return this.alignmentData.getCount(alignments);
	
	}

	/**
	 * Get total read count in set of regions
	 * @param alignments The regions
	 * @return The number of reads mapping to the regions
	 * @throws IOException
	 */
	private int getCount(Collection<Alignments> alignments) throws IOException{		
		return this.alignmentData.getCount(alignments);
	
	}

	/**
	 * Get total read count in gene
	 * @param gene The gene of interest
	 * @return The number of reads mapping to the gene
	 * @throws IOException
	 */
	public int getCount(RefSeqGene gene) throws IOException {
		return getCount(gene.getExonSet());
	}
	
	/**
	 * Get the alignment data
	 * @return The alignment data
	 */
	public ContinuousDataAlignmentModel getData() {
		return this.alignmentData;
	}
	
	/**
	 * For a RefSeqGene object get the Poisson parameter describing the null model for the number of reads falling in a window
	 * @param region The region of interest
	 * @param windowSize The window size
	 * @return Lambda, the expected number of reads in the window under the null Poisson model
	 * @throws IOException
	 */
	public double getPoissonLambda(RefSeqGene region, int windowSize) throws IOException {
		return getPoissonLambda(region.getExonSet(), windowSize);
	}
	
	/**
	 * For a collection of Alignments objects get the Poisson parameter describing the null model for the number of reads falling in a window
	 * @param regions The regions of interest
	 * @param windowSize The window size
	 * @return Lambda, the expected number of reads in the window under the null Poisson model
	 * @throws IOException
	 */
	private double getPoissonLambda(Collection<Alignments> regions, int windowSize) throws IOException {
		int totalSize = 0;
		for(Alignments region : regions) totalSize += region.getSize();
		return getPoissonLambda(totalSize, getCount(regions), windowSize);
	}
	
	/**
	 * Given a total number of reads, a region length, and a window size, get the parameter of the Poisson distribution modeling the number of reads in a window of this size
	 * @param regionLength The total region length
	 * @param numReads The number of reads
	 * @param windowSize The window size
	 * @return Lambda, the expected number of reads in the window under the null Poisson model
	 */
	private static double getPoissonLambda(int regionLength, int numReads, int windowSize) {
		if(numReads <= 0) {
			throw new IllegalArgumentException("Can't get poisson expected number of reads for " + numReads + " in window.");
		}
		return (double)windowSize * (double)numReads / regionLength;
	}

	/**
	 * Get total number of bases in reads mapped to region
	 * @param region The region
	 * @return The total number of bases in reads mapped to the region
	 * @throws IOException
	 */
	public int getTotalBasesInAlignments(Alignments region) throws IOException {
		int rtrn = 0;
		CloseableIterator<Alignment> iter = this.genericAlignmentDataModel.getReadIterator(region);
		while(iter.hasNext()) {
			Alignment a = iter.next();
			rtrn += Math.abs(a.getAlignmentEnd() - a.getAlignmentStart()) + 1;
		}
		return rtrn;
	}

	/**
	 * Get the average length of reads mapping to the region
	 * @param region The region
	 * @return The average read length
	 * @throws IOException 
	 */
	public double getAverageFragmentLength(Alignments region) throws IOException {
		return (double)this.getTotalBasesInAlignments(region) / (double)this.getCount(region);
	}

	
}
