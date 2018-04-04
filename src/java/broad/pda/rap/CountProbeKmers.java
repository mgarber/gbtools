/**
 *  Jesse Engreitz
 *  June 14, 2012
 *  
 *  Count the occurrence of kmers in the pulldown probes in the regions contained in the input annotation files.
 */
package broad.pda.rap;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import Jama.Matrix;
import broad.core.annotation.BED;
import broad.core.annotation.BEDReader;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.util.SequenceUtil;

/**
 * @author engreitz
 *
 */
public class CountProbeKmers {

	public CountProbeKmers(String probeFasta, String[] beds, int[] k, boolean allowMismatch, String genomeFasta, String out) throws Exception, IOException {
		for (int i = 0; i < k.length; i++) {
			countProbeKmers(probeFasta, beds, k[i], allowMismatch, genomeFasta, out);
		}
	}
	
	public void countProbeKmers(String probeFasta, String[] beds, int k, boolean allowMismatch, String genomeFasta, String out) throws Exception, IOException {
		Set<String> kmers = parseFastaToKmer(probeFasta, k, allowMismatch);
		Map<String, Integer> kmerCounts = initializeKmerCounts(kmers);
		
		out = out + ".k_" + k;
		BufferedWriter log = new BufferedWriter(new FileWriter(new File(out + ".log")));
		log.write("Found " + kmers.size() + " unique " + k + "-mers.\n");
		
		// Load the genome FASTA file for extracting sequences
		final ReferenceSequenceFile ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(new File(genomeFasta));
		
		// Read the first BED file and set up a matrix
		BEDReader b = new BEDReader(beds[0]);
		List<String> chrs = getChromosomesFromBED(b);
		List<String> bedList = new ArrayList<String>();
		for (String bed : beds) bedList.add(bed);
		
		MatrixWithHeaders m = new MatrixWithHeaders(new Matrix(new double[bedList.size()][chrs.size()]), bedList, chrs);
		
		for (int i=0; i<bedList.size(); i++) {
			System.out.println("Scanning " + beds[i] + " ...");
			MatrixWithHeaders toAdd = new MatrixWithHeaders(new Matrix(new double[1][chrs.size()]), bedList.subList(0,0), chrs);
			if (i>0) b = new BEDReader(beds[i]);
			countKmerMatchesInBED(toAdd, kmerCounts, b, k, allowMismatch, ref, log);
			m.setRow(bedList.get(i), toAdd.getRow(0));
		}
		
		m.write(out);
		writeKmerCounts(kmerCounts, out);
		log.close();
	}
	
	
	private void writeKmerCounts(Map<String, Integer> kmerCounts, String out) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(out + ".count")));
		
		List<Map.Entry<String, Integer>> list = new LinkedList(kmerCounts.entrySet());
		Collections.sort(list, new Comparator() {
			public int compare(Object o1, Object o2) {
	               return ((Comparable) ((Map.Entry) (o1)).getValue()).compareTo(((Map.Entry) (o2)).getValue());
	        }
		});
		Collections.reverse(list);
		
		for (Map.Entry<String, Integer> entry : list) {
			if (entry.getValue() > 0)
				bw.write(entry.getKey() + "\t" + entry.getValue() + "\n");
		}
		bw.close();
	}
	
	
	private Map<String, Integer> initializeKmerCounts(Set<String> kmers) {
		Map<String, Integer> kmerCounts = new HashMap<String, Integer>(); 
		for (String kmer : kmers) {
			Integer newInt = new Integer(0);
			kmerCounts.put(kmer, newInt);
		}
		return kmerCounts;
	}
	
		
	/**
	 * @param   BEDReader b
	 * @return	list of chromosomes present in this BED file
	 */
	private List<String> getChromosomesFromBED(BEDReader b) {
		Iterator<String> chrIt = b.getChromosomeIterator();
		List<String> chrs = new ArrayList<String>();
		while (chrIt.hasNext()) {
			chrs.add(chrIt.next());
		}
		return chrs;
	}
	
	
	/**
	 * Loop through a BED file and record how many intervals match the library of kmers
	 * @param counts
	 * @param kmers
	 * @param b
	 * @param k
	 * @param allowMismatch
	 * @param ref
	 * @throws Exception
	 * @throws IOException
	 * @throws ParseException
	 */
	private void countKmerMatchesInBED(MatrixWithHeaders counts, Map<String,Integer> kmerCounts, BEDReader b, int k, boolean allowMismatch, final ReferenceSequenceFile ref, BufferedWriter log) throws Exception, IOException, ParseException {
		Iterator<String> chrIt = b.getChromosomeIterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			log.write("Starting " + chr + " ...\n");
			counts.set(0, chr, 0);  // Initialize the count to zero
			IntervalTree<BED> tree = b.getChromosomeTree(chr);
			Iterator<BED> annotIt = tree.valueIterator();
			while(annotIt.hasNext()) {
				BED curr = annotIt.next();
				String seq = new String(ref.getSubsequenceAt(curr.getChromosome(), curr.getStart(), curr.getEnd()).getBases());
				String kmerMatch = stringMatchesKmers(kmerCounts, seq, k, allowMismatch);
				if (kmerMatch != null) {
					counts.set(0, chr, counts.get(0, chr) + 1);
					log.write(curr.toString() + "\t" + seq + "\t" + kmerMatch + "\n");
				}
			}
			log.flush();
		}
	}
	
	
	/**
	 * Searches a string for the kmers contained in kmerCounts
	 * @param kmers				Map of kmers and observation counts.  Map values are incremented in this function
	 * @param str				query string
	 * @param k					length of strings in 'kmers'
	 * @param allowMismatch		whether to allow one mismatch in the kmer
	 * @return					returns the first kmer match identified
	 */
	private String stringMatchesKmers(Map<String,Integer> kmerCounts, String str, int k, boolean allowMismatch) {
		str = str.toUpperCase();
		String found = null;
		
		if (!allowMismatch) {
			for (int i=0; i<=str.length()-k; i++) {
				String sub = str.substring(i, i+k);
				if (kmerCounts.keySet().contains(sub)) {
					found = sub;
					kmerCounts.put(sub, kmerCounts.get(sub) + 1);
					break;
				}
			}
		} else {
			for (int i=0; i<=str.length()-k; i++) {
				Set<String> mismatchStrs = getMismatchStrings(str.substring(i,i+k));
				for (String curr : mismatchStrs) {
					if (kmerCounts.keySet().contains(curr)) {
						found = curr;
						kmerCounts.put(curr, kmerCounts.get(curr) + 1);
						break;
					}
				}
			}
		}
		return found;
	}
	
	
	/**
	 * Read a FASTA file and return a set of all kmers (and their reverse complements).
	 * @param fasta				input FASTA file
	 * @param k					length of kmers
	 * @param allowMismatch		if true, generate reference set including degenerate kmers
	 * @return					set of kmers and their reverse complements
	 * @throws IOException
	 */
	private Set<String> parseFastaToKmer(final String fasta, int k, boolean allowMismatch) throws IOException {
		FastaSequenceIO fio = new FastaSequenceIO(fasta);
		List<Sequence> sequences = fio.loadAll();
		Set<String> rtrn = new HashSet<String>();
		for (int i=0; i<sequences.size(); i++) {
			rtrn.addAll(generateKmers(sequences.get(i).getSequenceBases(), k, allowMismatch));
			rtrn.addAll(generateKmers(SequenceUtil.reverseComplement(sequences.get(i).getSequenceBases()), k, allowMismatch));
		}
		return rtrn;
	}
	
	
	public Set<String> generateKmers(final String str, int k) {
		return generateKmers(str, k, false);
	}
	

	/**
	 * Generate all Kmers from a given string. 
	 * @param str			input string
	 * @param k				length of kmers to generate
	 * @param allowMismatch	If true, generate a set with degenerate kmers (add X's at each position)
	 * @return				set of kmers
	 */
	public Set<String> generateKmers(String str, int k, boolean allowMismatch) {
		str = str.toUpperCase();
		HashSet<String> set = new HashSet<String>();
		for (int i = 0; i <= str.length() - k; i++) {
			String sub = str.substring(i, i+k);
			if (!allowMismatch) {
				set.add(sub);
			} else {
				set.addAll(getMismatchStrings(sub));
			}
		}
		return set;
	}

	
	/**
	 * @param sub	input string
	 * @return		Set of strings where each character in the input is replaced with 'X'
	 */
	public static Set<String> getMismatchStrings(final String sub) {
		Set<String> set = new HashSet<String>();
		for (int j = 0; j < sub.length(); j++) {
			char[] base = sub.toCharArray();
			base[j] = 'X';
			set.add(new String(base));
		}
		return set;
	}
	
	
	/**
	 * @param str	input string containing ints delimited by commas
	 * @return		array of ints
	 */
	private static int[] splitStringToInts(final String str){
		String[] vals=str.split(",");
		int[] rtrn=new int[vals.length];
		for(int i=0; i<vals.length; i++){
			rtrn[i]=new Integer(vals[i]);
		}
		return rtrn;
	}
	
	
	static private String USAGE = "java -XX:MaxPermSize=1g -jar CountProbeKmers.jar -probe [ProbeFastaFile] -in [bed1,...] -k [kmer length] -genome [mm9.fa, must be indexed] -out [datamatrix] > [log.file]";
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception, IOException, ParseException {
		ArgumentMap argmap = CLUtil.getParameters(args, USAGE, "full");
		// TODO Auto-generated method stub

		String probeFasta = argmap.getMandatory("probe");
		String[] beds = argmap.getInput().split(",");
		int[] k = splitStringToInts(argmap.getMandatory("k"));
		boolean allowMismatch = argmap.containsKey("allowMismatch");
		String genomeFasta = argmap.getMandatory("genome");
		String out = argmap.getOutput();
		new CountProbeKmers(probeFasta, beds, k, allowMismatch, genomeFasta, out);
	}

}
