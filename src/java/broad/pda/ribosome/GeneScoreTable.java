/**
 * 
 */
package broad.pda.ribosome;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeSet;

import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;
import broad.core.sequence.SequenceUtils;
import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GeneScore;

/**
 * @author prussell
 *
 */
public class GeneScoreTable {
	
	
	/**
	 * Background expression
	 */
	private ContinuousDataAlignmentModel expressionData;

	/**
	 * Read alignments of ribosome-protected fragments to genome
	 */
	private ContinuousDataAlignmentModel ribosomeData;
	
	/**
	 * Read alignments of cycloheximide experiment to genome
	 */
	private ContinuousDataAlignmentModel cycloheximideData;
	
	/**
	 * Read alignments of harringtonin90s experiment to genome
	 */
	private ContinuousDataAlignmentModel harringtonin90Data;
	
	/**
	 * Read alignments of harringtonin120s experiment to genome
	 */
	private ContinuousDataAlignmentModel harringtonin120Data;
	
	/**
	 * Read alignments of harringtonin150s experiment to genome
	 */
	private ContinuousDataAlignmentModel harringtonin150Data;
	
	/**
	 * Read alignments of harringtonin180s experiment to genome
	 */
	private ContinuousDataAlignmentModel harringtonin180Data;
	
	/**
	 * Compute gene scores for ribosome dataset
	 */
	private ComputeRibosomeOccupancyByFeature ribosomeCrof;
	
	/**
	 * Compute gene scores for cycloheximide dataset
	 */
	private ComputeRibosomeOccupancyByFeature cycloheximideCrof;
	
	/**
	 * Output table file
	 */
	private String outfile;

	/**
	 * Map associating chromosome name with length
	 */
	private Map< String, Integer > chromosomeLengths;
	
	/**
	 * Name of gene set
	 */
	private String geneSet;
	
	/**
	 * @param geneSetFile file with two columns: gene set name and file path for bed file
	 * @param chrSizeFile file with two columns: chromsome name and size
	 * @param expressionFile expression bam file
	 * @param ribosomeFile ribosome bam file
	 * @param cycloheximideFile cycloheximide bam file
	 * @param harringtonin90file harringtonin90s bam file
	 * @param harringtonin120file harringtonin120s bam file
	 * @param harringtonin150file harringtonin150s bam file
	 * @param harringtonin180file harringtonin180s bam file
	 * @throws IOException
	 */
	public GeneScoreTable(String geneSetFile, String geneSetName, String chrSizeFile, String expressionFile, String ribosomeFile, String cycloheximideFile, String harringtonin90file, String harringtonin120file, String harringtonin150file, String harringtonin180file, String outFile, boolean normalizeByExpression, boolean fullyContainedReadsInCds) throws IOException {
		
		this.outfile = outFile;
		this.geneSet = geneSetName;
		
		System.err.println("Reading gene sets from file " + geneSetFile);
		loadGenes(geneSetFile);
		
		System.err.println("Reading chromosome information from file " + chrSizeFile);
		readChrsFromFile(chrSizeFile);
		
		System.err.println("Reading expression data from file " + expressionFile);
		this.expressionData = SequenceUtils.getDataModel(expressionFile, chrSizeFile, false);
		
		System.err.println("Reading ribosome data from file " + ribosomeFile);
		this.ribosomeData = SequenceUtils.getDataModel(ribosomeFile, chrSizeFile, false);
		this.ribosomeCrof = new ComputeRibosomeOccupancyByFeature(this.ribosomeData, this.expressionData);
		
		System.err.println("Reading cycloheximide data from file " + cycloheximideFile);
		this.cycloheximideData = SequenceUtils.getDataModel(cycloheximideFile, chrSizeFile, false);
		this.cycloheximideCrof = new ComputeRibosomeOccupancyByFeature(this.cycloheximideData, this.expressionData);
		
		System.err.println("Reading harringtonin90s data from file " + harringtonin90file);
		this.harringtonin90Data = SequenceUtils.getDataModel(harringtonin90file, chrSizeFile, false);
		
		System.err.println("Reading harringtonin120s data from file " + harringtonin120file);
		this.harringtonin120Data = SequenceUtils.getDataModel(harringtonin120file, chrSizeFile, false);
		
		System.err.println("Reading harringtonin150s data from file " + harringtonin150file);
		this.harringtonin150Data = SequenceUtils.getDataModel(harringtonin150file, chrSizeFile, false);
		
		System.err.println("Reading harringtonin180s data from file " + harringtonin180file);
		this.harringtonin180Data = SequenceUtils.getDataModel(harringtonin180file, chrSizeFile, false);
		
		System.err.println("\nDone reading data.\n");
		
		System.err.println("Scoring genes...\n");
		this.scoreGenes();		
		
		System.err.println("Writing table...\n");
		this.writeTable(this.geneSet, this.outfile, false, normalizeByExpression, fullyContainedReadsInCds);
		
	}
	
	
	
	/**
	 * Gene sets
	 * Key is gene set name
	 * Inner key is chromosome name
	 */
	private Map< String, Collection<RefSeqGene> > genes;
		

	/**
	 * The scores on genes for ribosome data
	 * Key is gene set name
	 * Inner key is chromosome name
	 */
	private Map< String, Map< String, Collection<GeneComponent> > > geneScoresRibosome;
	
	/**
	 * The scores on genes for cycloheximide data
	 * Key is gene set name
	 * Inner key is chromosome name
	 */
	private Map< String, Map< String, Collection<GeneComponent> > > geneScoresCycloheximide;

	/**
	 * The scores on genes for harringtonin90s data
	 * Key is gene set name
	 * Inner key is chromosome name
	 */
	private Map< String, Map< String, Collection<GeneScore> > > geneScoresHarringtonin90s;

	/**
	 * The scores on genes for harringtonin120s data
	 * Key is gene set name
	 * Inner key is chromosome name
	 */
	private Map< String, Map< String, Collection<GeneScore> > > geneScoresHarringtonin120s;

	/**
	 * The scores on genes for harringtonin150s data
	 * Key is gene set name
	 * Inner key is chromosome name
	 */
	private Map< String, Map< String, Collection<GeneScore> > > geneScoresHarringtonin150s;

	/**
	 * The scores on genes for harringtonin180s data
	 * Key is gene set name
	 * Inner key is chromosome name
	 */
	private Map< String, Map< String, Collection<GeneScore> > > geneScoresHarringtonin180s;

	
	/**
	 * Read file specifying chromosome names and sizes
	 * @param fileName file with two columns: chromosome name and size
	 * @throws IOException
	 */
	private void readChrsFromFile(String fileName) throws IOException {

		this.chromosomeLengths = new HashMap<String, Integer >();
		FileReader r = new FileReader(fileName);
		BufferedReader b = new BufferedReader(r);
		StringParser s = new StringParser();
		while(b.ready()) {
			String line = b.readLine();
			s.parse(line);
			if(s.getFieldCount() != 2) {
				throw new IllegalArgumentException("Chromosome file must have two columns: chromosome name, length");
			}
			this.chromosomeLengths.put(s.asString(0), Integer.valueOf(s.asInt(1)));
		}
		
	}
	
	/**
	 * Establish set of RefSeqGene objects
	 * @param fileName bed file
	 * @throws IOException 
	 */
	private void loadGenes(String fileName) throws IOException {
		
		this.genes = new HashMap< String, Collection<RefSeqGene >>();
		this.genes.putAll(BEDFileParser.loadDataByChr(new File(fileName)));
		
	}
	
	/**
	 * Score genes and write data to table
	 * @throws IOException
	 */
	private void scoreGenes() throws IOException {
		
		this.geneScoresCycloheximide = new HashMap<String, Map<String, Collection<GeneComponent>>>();
		this.geneScoresRibosome = new HashMap<String, Map<String, Collection<GeneComponent>>>();
		this.geneScoresHarringtonin90s = new HashMap<String, Map<String, Collection<GeneScore>>>();
		this.geneScoresHarringtonin120s = new HashMap<String, Map<String, Collection<GeneScore>>>();
		this.geneScoresHarringtonin150s = new HashMap<String, Map<String, Collection<GeneScore>>>();
		this.geneScoresHarringtonin180s = new HashMap<String, Map<String, Collection<GeneScore>>>();
			
		System.err.println("Scoring ribosome data for gene set " + this.geneSet);
		Map< String, Collection<GeneComponent> > tmpRibosomeScores = new HashMap< String, Collection<GeneComponent>>();
		for(String chr : this.genes.keySet()) {
			System.err.println(chr);
			tmpRibosomeScores.put(chr, this.ribosomeCrof.scoreFeatures(this.genes.get(chr)));
		}
		this.geneScoresRibosome.put(this.geneSet, tmpRibosomeScores);
			
		System.err.println("Scoring cycloheximide data for gene set " + this.geneSet);
		Map< String, Collection<GeneComponent> > tmpCycloheximideScores = new HashMap< String, Collection<GeneComponent>>();
		for(String chr : this.genes.keySet()) {
			System.err.println(chr);
			tmpCycloheximideScores.put(chr, this.cycloheximideCrof.scoreFeatures(this.genes.get(chr)));
		}
		this.geneScoresCycloheximide.put(this.geneSet, tmpCycloheximideScores);
			
		System.err.println("Scoring harringtonin90s data for gene set " + this.geneSet);
		Map< String, Collection<GeneScore> > tmp90Scores = new HashMap< String, Collection<GeneScore>>();
		for(String chr : this.genes.keySet()) {
			System.err.println(chr);
			Collection<GeneScore> tmpScores = new ArrayList<GeneScore>();
			for(RefSeqGene gene: this.genes.get(chr)){
					
				try {
					GeneComponent fullScore = RibosomeScoring.scoreFeature(gene, this.harringtonin90Data, this.harringtonin90Data.getAlignmentDataModelStats());
					double lambda = fullScore.getGeneScore().getNumberOfReads()/gene.getSize();
					RefSeqGene start = gene.getStartCodon();
					if(start == null) continue;
					RefSeqGene peak = this.harringtonin90Data.getPeak(gene, start);
					
					if(peak!=null){
						GeneScore startScore = new GeneScore(peak, this.harringtonin90Data.scoreGene(peak, lambda));
						tmpScores.add(startScore);
					}
				} catch(NullPointerException e) {
					continue;
				}
					
			}
			tmp90Scores.put(chr, tmpScores);
		}
		this.geneScoresHarringtonin90s.put(this.geneSet, tmp90Scores);
			
			
		System.err.println("Scoring harringtonin120s data for gene set " + this.geneSet);
		Map< String, Collection<GeneScore> > tmp120Scores = new HashMap< String, Collection<GeneScore>>();
		for(String chr : this.genes.keySet()) {
			System.err.println(chr);
			Collection<GeneScore> tmpScores = new ArrayList<GeneScore>();
			for(RefSeqGene gene: this.genes.get(chr)){
				
				try {
					GeneComponent fullScore = RibosomeScoring.scoreFeature(gene, this.harringtonin120Data, this.harringtonin120Data.getAlignmentDataModelStats());
					double lambda = fullScore.getGeneScore().getNumberOfReads()/gene.getSize();
					RefSeqGene start = gene.getStartCodon();
					if(start == null) continue;
					RefSeqGene peak = this.harringtonin120Data.getPeak(gene, start);
					
					if(peak!=null){
						GeneScore startScore = new GeneScore(peak, this.harringtonin120Data.scoreGene(peak, lambda));
						tmpScores.add(startScore);
					}
				} catch(NullPointerException e) {
					continue;
				}
			}
			tmp120Scores.put(chr, tmpScores);
		}
		this.geneScoresHarringtonin120s.put(this.geneSet, tmp120Scores);
			
			
		System.err.println("Scoring harringtonin150s data for gene set " + this.geneSet);
		Map< String, Collection<GeneScore> > tmp150Scores = new HashMap< String, Collection<GeneScore>>();
		for(String chr : this.genes.keySet()) {
			System.err.println(chr);
			Collection<GeneScore> tmpScores = new ArrayList<GeneScore>();
			for(RefSeqGene gene: this.genes.get(chr)){
					
				try {
					GeneComponent fullScore = RibosomeScoring.scoreFeature(gene, this.harringtonin150Data, this.harringtonin150Data.getAlignmentDataModelStats());
					double lambda = fullScore.getGeneScore().getNumberOfReads()/gene.getSize();
					RefSeqGene start = gene.getStartCodon();
					if(start == null) continue;
					RefSeqGene peak = this.harringtonin150Data.getPeak(gene, start);
					
					if(peak!=null){
						GeneScore startScore = new GeneScore(peak, this.harringtonin150Data.scoreGene(peak, lambda));
						tmpScores.add(startScore);
					}
				} catch(NullPointerException e) {
					continue;
				}
			}
			tmp150Scores.put(chr, tmpScores);
		}
		this.geneScoresHarringtonin150s.put(this.geneSet, tmp150Scores);
			
			
		System.err.println("Scoring harringtonin180s data for gene set " + this.geneSet);
		Map< String, Collection<GeneScore> > tmp180Scores = new HashMap< String, Collection<GeneScore>>();
		for(String chr : this.genes.keySet()) {
			System.err.println(chr);
			Collection<GeneScore> tmpScores = new ArrayList<GeneScore>();
			for(RefSeqGene gene: this.genes.get(chr)){
					
				try {
					GeneComponent fullScore = RibosomeScoring.scoreFeature(gene, this.harringtonin180Data, this.harringtonin180Data.getAlignmentDataModelStats());
					double lambda = fullScore.getGeneScore().getNumberOfReads()/gene.getSize();
					RefSeqGene start = gene.getStartCodon();
					if(start == null) continue;
					RefSeqGene peak = this.harringtonin180Data.getPeak(gene, start);
					
					if(peak!=null){
						GeneScore startScore = new GeneScore(peak, this.harringtonin180Data.scoreGene(peak, lambda));
						tmpScores.add(startScore);
					}
				} catch(NullPointerException e) {
					continue;
				}
			}
			tmp180Scores.put(chr, tmpScores);
		}
		this.geneScoresHarringtonin180s.put(this.geneSet, tmp180Scores);
						
		System.err.println("Done scoring gene set " + this.geneSet + ".\n");
			
						
		

	}
	
	private void writeTable(String geneSet, String outfile, boolean append, boolean normalizeByExpression, boolean fullyContainedReadsInCds) throws IOException {
		
		System.err.println("Writing table to file " + outfile);
		FileWriter w = new FileWriter(outfile, append);
		
		String fields = "GeneSet\t";
		fields += "GeneName\t";
		fields += "CDS\t";
		//fields += "Expression_GlobalEnrichment\t";
		//fields += "Expression_GlobalPvalue\t";
		fields += "RRS\t";
		fields += "Ribosome_GlobalEnrichment\t";
		fields += "Ribosome_GlobalPvalue\t";
		fields += "Cycloheximide_GlobalEnrichment\t";
		fields += "Cycloheximide_GlobalPvalue\t";
		fields += "Harringtonin90s_StartCodonLocalEnrichment\t";
		fields += "Harringtonin90s_StartCodonLocalPvalue\t";
		fields += "Harringtonin120s_StartCodonLocalEnrichment\t";
		fields += "Harringtonin120s_StartCodonLocalPvalue\t";
		fields += "Harringtonin150s_StartCodonLocalEnrichment\t";
		fields += "Harringtonin150s_StartCodonLocalPvalue\t";
		fields += "Harringtonin180s_StartCodonLocalEnrichment\t";
		fields += "Harringtonin180s_StartCodonLocalPvalue\t";
		fields += "\n";
		
		TreeSet<String> allIDs = new TreeSet<String>();
		HashMap<String,String> ribLines = new HashMap<String,String>();
		HashMap<String,String> cycLines = new HashMap<String,String>();
		HashMap<String,String> h90Lines = new HashMap<String,String>();
		HashMap<String,String> h120Lines = new HashMap<String,String>();
		HashMap<String,String> h150Lines = new HashMap<String,String>();
		HashMap<String,String> h180Lines = new HashMap<String,String>();
		
		System.err.println("Gathering data.");
		
		if(!append) w.write(fields);
		
		for(String chr : this.geneScoresRibosome.get(geneSet).keySet()) {
			for(GeneComponent score : this.geneScoresRibosome.get(geneSet).get(chr)) {
				String id = geneSet + "\t";
				id += score.getGeneName() + "\t";
				id += score.getGene().getCDSRegion().getChr() + ":" + score.getGene().getCDSRegion().getStart() + "-" + score.getGene().getCDSRegion().getEnd() + "\t";
				//String line = score.getGeneRNASeqScore().getEnrichment() + "\t";
				//line += score.getGeneRNASeqScore().getScanPvalue() + "\t";
				String line = RibosomeScoring.getRRS(score, 0.01, normalizeByExpression, fullyContainedReadsInCds) + "\t";
				line += score.getGeneScore().getEnrichment() + "\t";
				line += score.getGeneScore().getScanPvalue() + "\t";
				allIDs.add(id);
				ribLines.put(id, line);
			}
		}
			
		for(String chr : this.geneScoresCycloheximide.get(geneSet).keySet()) {
			for(GeneComponent score : this.geneScoresCycloheximide.get(geneSet).get(chr)) {
				String id = geneSet + "\t";
				id += score.getGeneName() + "\t";
				id += score.getGene().getCDSRegion().getChr() + ":" + score.getGene().getCDSRegion().getStart() + "-" + score.getGene().getCDSRegion().getEnd() + "\t";
				String line = score.getGeneScore().getEnrichment() + "\t";
				line += score.getGeneScore().getScanPvalue() + "\t";
				allIDs.add(id);
				cycLines.put(id, line);
			}
		}

			
		for(String chr : this.geneScoresHarringtonin90s.get(geneSet).keySet()) {
			for(GeneScore score : this.geneScoresHarringtonin90s.get(geneSet).get(chr)) {
				String id = geneSet + "\t";
				id += score.getGene().getName() + "\t";
				id += score.getGene().getCDSRegion().getChr() + ":" + score.getGene().getCDSRegion().getStart() + "-" + score.getGene().getCDSRegion().getEnd() + "\t";
				String line = score.getEnrichment() + "\t";
				line += score.getScanPvalue() + "\t";
				allIDs.add(id);
				h90Lines.put(id, line);
			}
		}

			
		for(String chr : this.geneScoresHarringtonin120s.get(geneSet).keySet()) {
			for(GeneScore score : this.geneScoresHarringtonin120s.get(geneSet).get(chr)) {
				String id = geneSet + "\t";
				id += score.getGene().getName() + "\t";
				id += score.getGene().getCDSRegion().getChr() + ":" + score.getGene().getCDSRegion().getStart() + "-" + score.getGene().getCDSRegion().getEnd() + "\t";
				String line = score.getEnrichment() + "\t";
				line += score.getScanPvalue() + "\t";
				allIDs.add(id);
				h120Lines.put(id, line);
			}
		}

		for(String chr : this.geneScoresHarringtonin150s.get(geneSet).keySet()) {
			for(GeneScore score : this.geneScoresHarringtonin150s.get(geneSet).get(chr)) {
				String id = geneSet + "\t";
				id += score.getGene().getName() + "\t";
				id += score.getGene().getCDSRegion().getChr() + ":" + score.getGene().getCDSRegion().getStart() + "-" + score.getGene().getCDSRegion().getEnd() + "\t";
				String line = score.getEnrichment() + "\t";
				line += score.getScanPvalue() + "\t";
				allIDs.add(id);
				h150Lines.put(id, line);
			}
				
		}

		for(String chr : this.geneScoresHarringtonin180s.get(geneSet).keySet()) {
			for(GeneScore score : this.geneScoresHarringtonin180s.get(geneSet).get(chr)) {
				String id = geneSet + "\t";
				id += score.getGene().getName() + "\t";
				id += score.getGene().getCDSRegion().getChr() + ":" + score.getGene().getCDSRegion().getStart() + "-" + score.getGene().getCDSRegion().getEnd() + "\t";
				String line = score.getEnrichment() + "\t";
				line += score.getScanPvalue() + "\n";
				allIDs.add(id);
				h180Lines.put(id, line);
			}
		}

		System.err.println("Writing table.");
				
		for(String id : allIDs) {
			
			String line = id;
			
			if(ribLines.get(id) == null ) {
				//line += "-\t-\t-\t-\t-\t";
				line += "-\t-\t-\t";
			} else {
				line += ribLines.get(id);
			}
			
			if(cycLines.get(id) == null) {
				line += "-\t-\t";
			} else {
				line += cycLines.get(id);
			}

			if(h90Lines.get(id) == null) {
				line += "-\t-\t";
			} else {
				line += h90Lines.get(id);
			}

			if(h120Lines.get(id) == null) {
				line += "-\t-\t";
			} else {
				line += h120Lines.get(id);
			}

			if(h150Lines.get(id) == null) {
				line += "-\t-\t";
			} else {
				line += h150Lines.get(id);
			}

			if(h180Lines.get(id) == null) {
				line += "-\t-\n";
			} else {
				line += h180Lines.get(id);
			}

			w.write(line);
			
		}
		
		w.close();
		
		System.err.println("All done.");
		
	}
	
	
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		// Set up the command line parser
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-g", "Gene set bed file", true);
		p.addStringArg("-gn", "Gene set name", true);
		p.addStringArg("-c","Chromosome file",true);
		p.addStringArg("-e","Expression bam file",true);
		p.addStringArg("-r","Ribosome bam file", true);
		p.addStringArg("-cy","Cycloheximide bam file",true);
		p.addStringArg("-h90","Harringtonin90s bam file",true);
		p.addStringArg("-h120","Harringtonin120s bam file",true);
		p.addStringArg("-h150","Harringtonin150s bam file", true);
		p.addStringArg("-h180", "Harringtonin180s bam file",true);
		p.addBooleanArg("-n","Normalize RRS by expression",true);
		p.addStringArg("-o", "Output table file", true);
		p.addBooleanArg("-f", "Only count fully contained reads in CDS", false, Boolean.valueOf(true));
		
		// Parse the command line and get argument values
		p.parse(args);

		String geneFile = p.getStringArg("-g");
		String geneSetName = p.getStringArg("-gn");
		String chrFile = p.getStringArg("-c");
		String expressionFile = p.getStringArg("-e");
		String ribosomeFile = p.getStringArg("-r");
		String cycFile = p.getStringArg("-cy");
		String h90File = p.getStringArg("-h90");
		String h120File = p.getStringArg("-h120");
		String h150File = p.getStringArg("-h150");
		String h180File = p.getStringArg("-h180");
		boolean normalizeByExpression = p.getBooleanArg("-n").booleanValue();
		String outfile = p.getStringArg("-o");
		boolean fullyContainedReadsInCds = p.getBooleanArg("-f").booleanValue();
		
		GeneScoreTable gst = new GeneScoreTable(geneFile, geneSetName, chrFile, expressionFile, ribosomeFile, cycFile, h90File, h120File, h150File, h180File, outfile, normalizeByExpression, fullyContainedReadsInCds);
		

	}

}
