/**
 * Collapse redundant probes and filter for various things
 * This class is designed to be used by ArrayDesigner by submitting LSF jobs calling the main method
 */
package broad.pda.capture.designer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import broad.core.motif.SearchException;
import broad.core.parser.StringParser;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.capture.designer.ArrayDesigner.Probe;
import broad.pda.capture.designer.ArrayDesigner.ProbePurpose;
import broad.pda.capture.designer.ArrayDesigner.SpeciesGene;
import broad.pda.rnai.designer.SmatchLike;
import jaligner.Alignment;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixGenerator;

/**
 * @author prussell
 *
 */
public class FilteredGeneProbes {
	
	// An ArrayDesigner object needed to instantiate some classes
	private ArrayDesigner arrayDesigner;
	// The set of isoform sequences and all their probes
	private Map< Sequence, TreeSet<Probe> > isoformProbeSet;
	// The collapsed filtered set of probes
	private TreeSet<Probe> geneProbes;
	// The transcriptome to avoid
	private Collection<Sequence> genes;
	// The set of qPCR primer sequences
	private Collection<Sequence> qpcrPrimers;
	// The SpeciesGene identifier of this set of isoforms
	private SpeciesGene speciesgene;
	/**
	 * Constructor requires information about the parent gene, probes that have already been generated, transcriptome, and existing qPCR primers
	 * @param speciesGene the SpeciesGene object whose set of probes this is
	 * @param isoformProbesFile the file of probes against all isoforms as generated and writtend from within ArrayDesigner class
	 * @param geneFile fasta file of the transcriptome
	 * @param qpcrFile fasta file of qPCR primer sequences
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public FilteredGeneProbes( String speciesGene, String isoformProbesFile, String geneFile, String qpcrFile) throws FileNotFoundException, IOException {
		this.arrayDesigner = new ArrayDesigner();
		this.isoformProbeSet = new HashMap< Sequence, TreeSet<Probe>>();
		this.loadIsoformProbes(isoformProbesFile);
		this.genes = new HashSet<Sequence>();
		this.loadGenes(geneFile);
		this.qpcrPrimers = new HashSet<Sequence>();
		this.loadQpcrPrimers(qpcrFile);
		this.speciesgene = this.arrayDesigner.new SpeciesGene();
		this.speciesgene.getFromStringRepresentation(speciesGene);
		this.geneProbes = new TreeSet<Probe>();
	}

	/**
	 * Remove gene probes that hybridize to a different gene
	 * @throws SearchException
	 */
	private void filterCrossHybs() throws SearchException {
		
		if(this.genes.isEmpty()) {
			throw new IllegalStateException("Can't filter gene probes for cross hybs because gene set is empty.");
		}
		
		if(this.geneProbes.isEmpty()) {
			throw new IllegalStateException("Gene probe set is empty.  Need to call makeInitialGeneProbes.");
		}
		
		// Copy probes into new set to filter
		TreeSet<Probe> filteredGeneProbes = new TreeSet<Probe>(this.geneProbes);
		// Check each probe
		for(Probe probe : this.geneProbes) {
			// Filter if cross hybridizes to another gene
			if(this.probeCrossHybsToTranscriptome(probe,this.speciesgene)) {
				filteredGeneProbes.remove(probe);
			}
		}
		// Replace probe set with filtered set
		this.geneProbes.clear();
		this.geneProbes.addAll(filteredGeneProbes);
	}

	/**
	 * Filter probes with more than the allowable number of soft masked bases
	 */
	private void filterRepetetiveProbes() {
		// Only repeat mask pulldown probes
		if(this.speciesgene.getPurpose() == ProbePurpose.PULLDOWN) {
			if(this.isoformProbeSet.isEmpty()) {
				throw new IllegalStateException("Can't filter repetetive probes because initial isoform probe set is empty.");
			}
			for(Sequence seq : this.isoformProbeSet.keySet()) {
				// Copy allowable probes into new set
				TreeSet<Probe> filtered = new TreeSet<Probe>();
				for(Probe probe : this.isoformProbeSet.get(seq)) {
					if(probe.getSequence().countSoftMaskedBases() <= ArrayDesigner.REPEAT_CUTOFF) {
						// If number of repeat bases is tolerable, make entire probe uppercase and keep
						Sequence upper = probe.getSequence();
						upper.uppercase();
						probe.setSequence(upper);
						filtered.add(probe);
					}
				}
				// Replace the isoform probes with the filtered set
				this.isoformProbeSet.put(seq, filtered);
			}
		}
	}

	/**
	 * Filter redundant probes
	 * @return the filtered set
	 */
	private TreeSet<Probe> getUniqueOptimalTilingPath() {
		
		System.out.println("Collapsing probes...");
		
		TreeSet<Probe> returnTilingPath = new TreeSet<Probe>();
		
		// Throw exception if isoform set is empty
		if(this.isoformProbeSet.keySet().isEmpty()) {
			throw new IllegalArgumentException("Can't collapse probes for empty set of isoforms.");
		}
		
		// If gene has only one isoform, return the original probe set for the single isoform
		if(this.isoformProbeSet.keySet().size() == 1) {
			Iterator<Sequence> iter = this.isoformProbeSet.keySet().iterator();
			Sequence onlyIsoform = iter.next();
			returnTilingPath = this.isoformProbeSet.get(onlyIsoform);
			//System.out.println("There is only one isoform " + onlyIsoform.getId());
			return returnTilingPath;
		}
		
		// Pick the longest isoform as the reference
		Sequence reference = new Sequence("");
		reference.setSequenceBases("");
		for(Sequence seq : this.isoformProbeSet.keySet()) {
			if(seq.getLength() > reference.getLength()) reference = seq;
		}
		
		//System.out.println("For reference picked sequence " + reference.getId() + " with bases " + reference.getSequenceBases());
		
		// Initialize return set to the probes of the reference isoform
		returnTilingPath = this.isoformProbeSet.get(reference);
		
		for(Sequence isoform : this.isoformProbeSet.keySet()) {
			// Ignore the reference isoform
			if(isoform == reference) continue;
			// Recall the isoform probes for the current isoform
			int perfectMatches = 0;
			int partialMatches = 0;
			int noMatch = 0;
			TreeSet<Probe> currIsoformProbes = this.isoformProbeSet.get(isoform);
			for(Probe probe : currIsoformProbes) {
				// Check for matches among probes currently stored in returnTilingPath
				double bestPercentMatch = -1;
				// Identify best matching probe
				for(Probe refprobe : returnTilingPath) {
					Alignment align = localAlign(probe.getSequenceBases(),refprobe.getSequenceBases());
					// align.getSequence1() is the aligned part of probe
					double percentMatch = (double)align.getSequence1().length/(double)ArrayDesigner.OLIGO_SIZE;
					if(percentMatch > bestPercentMatch) {
						bestPercentMatch = percentMatch;
					}
				}
				// Pass over probes with a perfect match already in set
				if(bestPercentMatch > 0.9) {
					//System.out.println("Skipping probe due to perfect match " + probe.getSequenceBases());
					perfectMatches++;
				}
				else if(bestPercentMatch > 0.3 && bestPercentMatch <= 0.9) {
					// Keep the semi-redundant probe but erase coordinates and add to the special extra tiling path
					probe.setStartPos(-1);
					probe.setEndPos(-1);
					probe.setTilingPath(-1);
					probe.setEvenOdd(-1);
					probe.setDomain(-1);
					returnTilingPath.add(probe);
					//System.out.println("Keeping semi-redundant probe in tiling path -1 due to partial match " + probe.getSequenceBases());
					partialMatches++;
				}
				else {
					// Keep the new probe but erase coordinates and add to the special extra tiling path
					probe.setStartPos(-1);
					probe.setEndPos(-1);
					probe.setTilingPath(-1);	
					probe.setEvenOdd(-1);
					probe.setDomain(-1);
					returnTilingPath.add(probe);
					//System.out.println("Keeping probe in tiling path -1 " + probe.getSequenceBases());
					noMatch++;
				}
			}
			currIsoformProbes = null;
			// Throw away isoform probes to save memory
			// this.isoformProbeSet.remove(isoform);
			
		}
		
		return returnTilingPath;
	}

	/**
	 * Load transcriptome from file
	 * @throws IOException thrown by FastaSequenceIO instance
	 */
	private void loadGenes(String fileName) throws IOException {
		System.out.println("Loading genes from file " + fileName + "...");
		this.genes.clear();
		FastaSequenceIO fsio = new FastaSequenceIO(new File(fileName));
		this.genes= fsio.loadAll();
		System.out.println("Loaded " + this.genes.size() + " genes.\n");
	}
	
	/**
	 * Read isoform probes from a file written by ArrayDesigner
	 * @param fileName
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	private void loadIsoformProbes(String fileName) throws FileNotFoundException, IOException {
		
		System.out.println("Loading isoform probes from file " + fileName);
		
		this.isoformProbeSet.clear();
		
		FileReader reader = new FileReader(fileName);
		BufferedReader b = new BufferedReader(reader);
				
		// Each line of the isoform file has the full isoform sequence in the first column and information for one probe in remaining columns
		// In order to identify the set of probes for each isoform, we will read in the file and combine probes that share an isoform sequence
		Map<String, TreeSet<Probe>> basesToProbes = new HashMap<String, TreeSet<Probe>>();
		
		while(b.ready()) {
			
			String line = b.readLine();
			StringParser p = new StringParser();
			p.parse(line);
			
			if(p.getFieldCount() != 8) {
				throw new IllegalArgumentException("File " + fileName + " is in wrong format, must have 8 fields per line.");
			}
			
			TreeSet<Probe> tempprobes = new TreeSet<Probe>();
			String bases = p.asString(0);
			
			// Copy the set of exising probes for the isoform if it exists
			if(basesToProbes.containsKey(bases)) tempprobes.addAll(basesToProbes.get(bases));
			Probe probe = this.arrayDesigner.new Probe();
			// Read the rest of the fields to get the new probe
			probe.readFromString(p.asString(1) + " " + p.asString(2) + " " + p.asString(3) + " " + p.asString(4) + " " + p.asString(5) + " " + p.asString(6) + " " + p.asString(7));
			// Add the new probe to the isoform set and put back into the map
			tempprobes.add(probe);
			basesToProbes.put(bases, tempprobes);
		}
		
		// Count the number of probes loaded
		int numProbes = 0;
		for(String bases : basesToProbes.keySet()) {
			Sequence seq = new Sequence(bases);
			seq.setSequenceBases(bases);
			this.isoformProbeSet.put(seq,basesToProbes.get(bases));
			numProbes += this.isoformProbeSet.get(seq).size();
		}
		
		System.out.println("Loaded " + numProbes + " probes.\n");
		
	}

	/**
	 * Load qPCR primers from file
	 * @param fileName
	 * @throws IOException
	 */
	private void loadQpcrPrimers(String fileName) throws IOException {
		System.out.println("Loading qPCR primers from file " + fileName + "...");
		this.qpcrPrimers.clear();
		FastaSequenceIO fsio = new FastaSequenceIO(new File(fileName));
		this.qpcrPrimers= fsio.loadAll();
		System.out.println("Loaded " + this.qpcrPrimers.size() + " qPCR primers.\n");
	}

	/**
	 * Align two sequences
	 * @param seq1 Sequence 1
	 * @param seq2 Sequence 2
	 * @return the Alignment object
	 */
	private static Alignment localAlign(String seq1, String seq2) {
		jaligner.Sequence s1 = new jaligner.Sequence(seq1);
		jaligner.Sequence s2 = new jaligner.Sequence(seq2);
		Matrix matrix = MatrixGenerator.generate(1.0f, -100000.0f); // To get the perfect match which we expect here
		jaligner.Alignment alignment = SmithWatermanGotoh.align(s1, s2, matrix, 1000000, 1000000);
		return alignment;
	}

	/**
	 * The main method is only designed to be called from a FilteredGeneProbes object which in turn is instantiated within an ArrayDesigner object
	 * @param args
	 * @throws IOException
	 * @throws FileNotFoundException
	 * @throws SearchException
	 */
		public static void main(String[] args) throws IOException, FileNotFoundException, SearchException {
			
			if(args.length != 5) {
				System.err.println("SpeciesGene Probes GeneFile OutFile QpcrFile");
				System.exit(-1);
			}
			
			String speciesGeneString = args[0];
			String probeFile = args[1];
			String geneFile = args[2];
			String outFile = args[3];
			String qpcrFile = args[4];
			
			FilteredGeneProbes f = new FilteredGeneProbes(speciesGeneString, probeFile, geneFile, qpcrFile);
			
			// Filter repetetive probes
			f.filterRepetetiveProbes();
			
			// Collapse to genes and eliminate redundancy
			f.makeInitialGeneProbes();
			
			// Remove probes that cross-hybridize in transcriptome
			// Takes forever, need to fix or just blast probes to transcriptome later
			//f.filterCrossHybs();
			
			// Determine if probes contain a qPCR primer
			f.markQpcrPrimers(qpcrFile);
			
			// Write gene probes to file
			f.writeGeneProbes(outFile);
	
		}

	/**
	 * Clear gene probe set and replace with initial collapsed set of probes
	 */
	private void makeInitialGeneProbes() {

		// Count the number of initial probes to print comparison later
		int initSize = 0;
		for(Sequence seq : this.isoformProbeSet.keySet()) {
			initSize += this.isoformProbeSet.get(seq).size();
		}
		
		if(this.isoformProbeSet.isEmpty()) {
			throw new IllegalStateException("Cannot collapse probes to gene level because isoform probe set is empty.");
		}

		this.geneProbes.clear();

		// Collapse probe sets of all isoforms and get one optimal tiling path
		TreeSet<Probe> filtered = this.getUniqueOptimalTilingPath();
		this.geneProbes.addAll(filtered);
		
		System.out.println("Collapsed from " + initSize + " to " + this.geneProbes.size() + " probes.\n");
		
	}

	/**
	 * Identify probes containing a qPCR primer site and mark them
	 * @param file the fasta file of qPCR primers
	 * @throws SearchException
	 */
	private void markQpcrPrimers(String file) throws SearchException {
		if(this.geneProbes.isEmpty()) {
			throw new IllegalStateException("Can't mark qPCR primer sites because gene probe set is empty.");
		}
		if(this.qpcrPrimers.isEmpty()) {
			throw new IllegalStateException("qPCR primer set is empty.");
		}
		System.out.println("Marking probes that contain a qPCR primer...");
		int numFound = 0;
		// New probe set to modify and replace the existing one
		TreeSet<Probe> filteredProbes = new TreeSet<Probe>();
		for(Probe probe : this.geneProbes) {
			//System.out.println("Checking probe for qPCR primer sites: " + probe);
			if(probeContainsQpcrPrimer(probe)) {
				numFound++;
				//System.out.println("Probe contains qPCR primer site: " + probe);
				// Replace with a probe with same fields except qPCR primer set to 1
				Probe newprobe = this.arrayDesigner.new Probe();
				newprobe.setSequence(probe.getSequence());
				newprobe.setStartPos(probe.getStartPos());
				newprobe.setEndPos(probe.getEndPos());
				newprobe.setDomain(probe.getDomain());
				newprobe.setTilingPath(probe.getTilingPath());
				newprobe.setEvenOdd(probe.getEvenOdd());
				newprobe.setQpcrPrimer(1);
				filteredProbes.add(newprobe);
			} else filteredProbes.add(probe);
			probe = null;
		}
		
		//System.out.println("Found " + numFound + " probes containing qPCR primer sites.\n");
		// Replace geneProbes with the modified set
		this.geneProbes.clear();
		this.geneProbes.addAll(filteredProbes);
	}

	/**
	 * Whether the probe contains sequence matching a qPCR primer
	 * @param probe the probe
	 * @return whether a match is found
	 * @throws SearchException
	 */
	private boolean probeContainsQpcrPrimer(Probe probe) throws SearchException {
		if(this.qpcrPrimers.isEmpty()) {
			throw new IllegalStateException("qPCR primer set is empty.");
		}
		// Make a list of just this probe because SmatchLike takes a list
		List<Sequence> shortlist = new ArrayList<Sequence>();
		shortlist.add(probe.getSequence());
		// Check each primer as query with the probe as target
		for(Sequence q : this.qpcrPrimers) {
			SmatchLike smatch = new SmatchLike(q.getSequenceBases(), shortlist, q.getLength());
			// The collection of matches which is either the single probe or empty
			Collection<String> crossHybs = smatch.getAllPossibleTargets(0);
			// Probe matches a primer
			if(!crossHybs.isEmpty()) {
				shortlist = null;
				smatch = null;
				crossHybs = null;
				return true;
			}
			crossHybs = null;
			smatch = null;
		} 
		shortlist = null;
		return false;
	}

	/**
	 * Test if probe sequence cross hybridizes to any other gene in transcriptome
	 * @param probe the probe
	 * @param probeSpeciesGene the name of the gene the probe came from
	 * @return whether the probe sequence hybridizes to a different gene
	 * @throws SearchException if SmatchLike instance throws SearchException
	 */
	private boolean probeCrossHybsToTranscriptome(Probe probe, SpeciesGene probeSpeciesGene) throws SearchException {
		String probeSequence = probe.getSequenceBases();
		//System.out.println("Looking for cross-hybs for probe " + probeSpeciesGene + " " + probeSequence);
		SmatchLike smatch = new SmatchLike(probeSequence,(List<Sequence>) this.genes,ArrayDesigner.MATCH_CUTOFF);
		// The sequence IDs for cross hybs
		Collection<String> crossHybs = smatch.getAllPossibleTargets(0);
		// Check each cross hybridizing sequence to see if it's the same gene
		for(String sequenceName : crossHybs) {
			Sequence sequence = new Sequence(sequenceName);
			String sequenceSpecies = sequence.getRefseqSpeciesName();
			String sequenceGene = sequence.getRefseqGeneName();
			SpeciesGene sequenceSpeciesGene = this.arrayDesigner.new SpeciesGene(sequenceSpecies,sequenceGene);
			//System.out.println("Found match to gene " + sequenceSpeciesGene.toString());
			// If cross hybridizes to another gene, return true
			if(!sequenceSpeciesGene.sameGeneAndSpecies(probeSpeciesGene)) {
				System.out.println("Cross-hyb to different gene found: " + sequenceSpeciesGene.toString());
				return true;
			} //else System.out.println("But it's the same gene.");
			
		}
		//System.out.println("No cross hybs found.");
		return false;
	}

	
	/**
	 * Write gene probe set to file
	 * @param file
	 * @throws IOException
	 */
	private void writeGeneProbes(String file) throws IOException {
		if(this.geneProbes.isEmpty()) throw new IllegalStateException("Gene probe set is empty. Need to make gene probes first.");
		FileWriter writer = new FileWriter(file);
		// Write the string representation of probe which can be read back from within the Probe class
		for(Probe probe : this.geneProbes) {
			writer.write(probe + "\n");
		}
		writer.close();
		System.out.println("Wrote filtered and collapsed probes to file " + file + "\n");
	}

}


