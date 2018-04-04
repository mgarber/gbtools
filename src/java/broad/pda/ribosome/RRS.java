/**
 * 
 */
package broad.pda.ribosome;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;
import java.util.TreeSet;

import org.broad.igv.Globals;

import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;
import broad.core.sequence.SequenceUtils;
import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;

/**
 * @author prussell
 *
 */
public class RRS {
	
	/**
	 * Read alignments of ribosome-protected fragments to genome
	 */
	private ContinuousDataAlignmentModel ribosomeData;
	
	/**
	 * Background expression
	 */
	private ContinuousDataAlignmentModel expressionData;
	
	/**
	 * Map associating chromosome name with length
	 */
	private Map< String, Integer > chromosomeLengths;
	
	/**
	 * The chromosome names
	 */
	private Collection< String > chromosomeNames;
	
	/**
	 * Gene set
	 */
	private Map< String, Collection<RefSeqGene> > genes;
		
	/**
	 * ComputeRibosomeOccupancyByFeature object to generate gene scores
	 */
	private ComputeRibosomeOccupancyByFeature crof;
	
	/**
	 * The scores on genes
	 */
	private Map< String, Collection<GeneComponent> > geneScores;
	
	
	public static RRS NewInstance(ContinuousDataAlignmentModel ribosomeData, ContinuousDataAlignmentModel expressionData, Map< String, Collection<RefSeqGene> > genes) throws IOException {
		
		Globals.setHeadless(true);
		
		RRS ovu = new RRS();
		
		ovu.ribosomeData = ribosomeData;
		ovu.expressionData = expressionData;
		ovu.genes = genes;
		ovu.chromosomeLengths = ovu.ribosomeData.getChromosomeLengths();
		ovu.chromosomeNames = ovu.chromosomeLengths.keySet();
		ovu.crof = new ComputeRibosomeOccupancyByFeature(ovu.ribosomeData, ovu.expressionData);
		
		ovu.geneScores = new TreeMap< String, Collection<GeneComponent> >();
		
		for(String chr : ovu.chromosomeNames) {
			if(ovu.genes.containsKey(chr)) {
				if(!ovu.genes.get(chr).isEmpty()) ovu.geneScores.put(chr , ovu.crof.scoreFeatures(ovu.genes.get(chr)));
			}
		}

		return ovu;
		
	}
	
	/**
	 * Get a RibosomeProfilingAnalysis object with only one chromosome
	 * @param ribosomeBamFile Bam file of read alignments from ribosome-protected fragments
	 * @param expressionFile RNA-seq bam file
	 * @param chrSizeFile Chromosome size file
	 * @param geneFile Bed file of coding genes or null
	 * @param chr the chromosome to use
	 * @return the RibosomeProfilingAnalysis object with underlying alignments on the given chromosome only
	 * @throws IOException
	 */
	public static RRS NewInstance(String ribosomeBamFile, String expressionFile, String chrSizeFile, String geneFile, String chr) throws IOException {
		
		// Make a new chromosome size file with only one chromosome, to avoid computing global stats on all chromosomes when instantiating AlignmentDataModelStats object
		File dir = new File("tmp");
		boolean madedir = dir.mkdir();
		Random r = new Random();
		String tmpChrSizeFile = dir + "/tmp_" + Integer.valueOf(r.nextInt()).toString();
		FileWriter w = new FileWriter(tmpChrSizeFile);
		StringParser sp = new StringParser();
		FileReader r1 = new FileReader(new File(chrSizeFile));
		BufferedReader b = new BufferedReader(r1);
		boolean found = false;
		while(b.ready()) {
			String line = b.readLine();
			sp.parse(line);
			if(sp.asString(0).equals(chr)) {
				w.write(line + "\n");
				w.close();
				found = true;
				break;
			}
		}
		if(!found) {
			throw new IllegalArgumentException("Chromosome " + chr + " is not present in file " + chrSizeFile);
		}
		
		// Get the ContinuousDataAlignmentModels
		ContinuousDataAlignmentModel radm = SequenceUtils.getDataModel(ribosomeBamFile, tmpChrSizeFile, chr, false);
		ContinuousDataAlignmentModel eadm = SequenceUtils.getDataModel(expressionFile, tmpChrSizeFile, chr, false);
		
		// Read the genes
		Map< String, Collection<RefSeqGene> > genes = new TreeMap< String, Collection<RefSeqGene> >();
		if(geneFile != null) genes.putAll(BEDFileParser.loadDataByChr(new File(geneFile)));
		
		return NewInstance(radm, eadm, genes);
	}
	
	/**
	 * For each gene on the given chromosome, write RPKM scores of ORFs and UTRs, and the ratio. Write separate files for coding, noncoding, lincRNAs. Skip gene class if has not been provided.
	 * @param outfilePrefix Prefix for output files, to be used for separate files for coding, noncoding, lincRNAs.
	 * @param append Whether to append to existing output files.
	 * @param chr Chromsome name
	 * @throws IOException
	 */
	private void writeRRS(String outfilePrefix, boolean append, String chr, boolean normalizeByExpression, boolean fullyContainedReadsInCds) throws IOException {
		if(!this.genes.isEmpty()) ComputeRibosomeOccupancyByFeature.writeRRS(outfilePrefix , this.geneScores.get(chr), append, normalizeByExpression, fullyContainedReadsInCds);
	}
	
	
	/**
	 * For each gene, write RPKM scores of ORFs and UTRs, and the ratio. Write separate files for coding, noncoding, lincRNAs. Skip gene class if has not been provided.
	 * @param outfilePrefix Prefix for output files, to be used for separate files for coding, noncoding, lincRNAs.
	 * @throws IOException
	 */
	private void writeRRS(String outfilePrefix, boolean normalizeByExpression, boolean fullyContainedReadsInCds) throws IOException {
		boolean started = false;
		for(String chr : this.chromosomeNames) {
			this.writeRRS(outfilePrefix, started, chr, normalizeByExpression, fullyContainedReadsInCds);
			started = true;
		}
	}
	
	/**
	 * For each gene on the given chromosome, write max ratio of RPKM scores of ORFs and UTRs. Write separate files for coding, noncoding, lincRNAs. Skip gene class if has not been provided.
	 * @param outfilePrefix Prefix for output files, to be used for separate files for coding, noncoding, lincRNAs.
	 * @param append Whether to append to existing output files.
	 * @param chr Chromsome name
	 * @param translationalEfficiencyRatio whether to calculate the ratio of translational efficiency as opposed to RPKM
	 * @throws IOException
	 */
	private void writeMaxRRS(String outfilePrefix, boolean append, String chr, boolean normalizeByExpression, boolean fullyContainedReadsInCds) throws IOException {
		if(!this.genes.isEmpty()) RibosomeScoring.writeMaxRRSPerGene(this.geneScores.get(chr), ComputeRibosomeOccupancyByFeature.ALPHA, outfilePrefix, append, normalizeByExpression, fullyContainedReadsInCds);
	}
	
	
	/**
	 * For each gene, write max ratio of RPKM scores of ORFs and UTRs. Write separate files for coding, noncoding, lincRNAs. Skip gene class if has not been provided.
	 * @param outfilePrefix Prefix for output files, to be used for separate files for coding, noncoding, lincRNAs.
	 * @param translationalEfficiencyRatio whether to calculate the ratio of translational efficiency as opposed to RPKM
	 * @throws IOException
	 */
	private void writeMaxRRS(String outfilePrefix, boolean normalizeByExpression, boolean fullyContainedReadsInCds) throws IOException {
		boolean append = false;
		for(String chr : this.chromosomeNames) {
			this.writeMaxRRS(outfilePrefix, append, chr, normalizeByExpression, fullyContainedReadsInCds);
			append = true;
		}
	}
	
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		// Set up the command line parser
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-br", "bam file of read alignments from ribosome-protected fragments", true);
		p.addStringArg("-be", "expression bam file", true);
		p.addStringArg("-s", "chromosome size file", true);
		p.addStringArg("-g", "bed file of genes", true);
		p.addStringArg("-o", "output directory for analysis", true);
		p.addStringArg("-chr", "chromosome to restrict analysis to", false);
		p.addBooleanArg("-n", "normalize RRS by expression", true);
		p.addBooleanArg("-f", "only count fully contained reads in CDS", false, Boolean.valueOf(true));

		// Parse the command line and get argument values
		p.parse(args);
		String ribosomeBamFile = p.getStringArg("-br");
		String expressionBamFile = p.getStringArg("-be");
		String geneFile = p.getStringArg("-g");
		String outdir = p.getStringArg("-o");
		String sizeFile = p.getStringArg("-s");
		String chr = p.getStringArg("-chr");
		boolean normalizeByExpression = p.getBooleanArg("-n").booleanValue();
		boolean fullyContainedReadsInCds = p.getBooleanArg("-f").booleanValue();
		
		// Make output directory
		File dir = new File(outdir);
		boolean madedir = dir.mkdir();
		
		// Read chromosome names
		Collection<String> chrsToUse = new TreeSet<String>();
		if(chr != null) chrsToUse.add(chr);
		else {
			StringParser sp = new StringParser();
			FileReader r = new FileReader(new File(sizeFile));
			BufferedReader b = new BufferedReader(r);
			while(b.ready()) {
				String line = b.readLine();
				sp.parse(line);
				chrsToUse.add(sp.asString(0));
			}
			
		}
		
		// Write data
		boolean append = false;
		for(String chromosome : chrsToUse) {
			// Get a RibosomeProfilingAnalysis object
			RRS rpa = RRS.NewInstance(ribosomeBamFile, expressionBamFile, sizeFile, geneFile, chromosome);
		
			// Write RRS scores
			String rpkmRatioFilePrefix = outdir + "/RRS";
			rpa.writeRRS(rpkmRatioFilePrefix, append, chromosome, normalizeByExpression, fullyContainedReadsInCds);
			
			// Write max RRS per gene
			String maxRpkmRatioFilePrefix = outdir + "/MaxRRS";
			rpa.writeMaxRRS(maxRpkmRatioFilePrefix, append, chromosome, normalizeByExpression, fullyContainedReadsInCds);
					
			append = true;
		}
		
		
		
		
		
		
	}

}
