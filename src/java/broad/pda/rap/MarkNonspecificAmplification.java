/*
 * Jesse Engreitz
 * February 23, 2012
 * 
 * Goal: Identify reads that resulted from non-specific amplification during either the RT or PCR enrichment steps.
 *       In this experiment, RNA was ligated with a random barcode-containing adaptor, which should represent the
 *       first 10 base pairs of read 1 (ANNNNNNNNT) or (TNNNNNNNNA).  In order to align these reads (R1), I trimmed
 *       10 bp from the front of the read and aligned the rest.  Since it's a random barcode, what we need 
 *       to do is locate where the read mapped to, then compare the 10 adjacent base pairs to the discarded read
 *       sequence.  If it's the same, then this read must be the result of non-specific RT or amplification.  If 
 *       it is not the expected RNA sequence, and it contains an A or T in position 10, then it's likely specific
 *       amplification.
 *       
 *       Sort reads into 3 bins:  specific (contains A/T), non-specific (matches RNA) or other.
 *       
 *       Writes a new SAM file with nonspecific RT products marked with the 0x800 flag (not usually used by SAM).
 *       Outputs the number of reads at each strand as well.
 *     
 */

package broad.pda.rap;


import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.io.FilenameUtils;

import broad.pda.seq.alignment.Pair;
import broad.util.ExtensionFilter;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.reference.FastaSequenceIndex;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.SequenceUtil;


/**
 * For marking nonspecific amplification (by identifying the RT primer?) in RAP-RNA sequencing
 * @author engreitz
 *
 */
public class MarkNonspecificAmplification {

	public MarkNonspecificAmplification(File[] fastqFiles, File[] samFiles, String reference, String referenceIndex, String saveDir) throws IOException, InterruptedException {
		System.out.println("Found "+fastqFiles.length+" FASTQ files.");
		System.out.println("Found "+samFiles.length+" SAM files.");
		Map<String, Pair<File>> pairs = getPairs(fastqFiles, samFiles);
		System.out.println("Found "+pairs.size()+" pairs.\n");
		markNonspecificAmplification(pairs, reference, referenceIndex, saveDir);
	}
	
	private static Map<String, Pair<File>> getPairs(File[] fastqFiles, File[] samFiles) throws InterruptedException {
		Map<String, File> hold = new TreeMap<String, File>();
		Map<String, Pair<File>> rtrn=new TreeMap<String, Pair<File>>();
		
		for (int i=0; i < samFiles.length; i++) {
			String name = samFiles[i].getName();
			String baseName = FilenameUtils.removeExtension(name);
			hold.put(baseName, samFiles[i]);
		}
		
		for (int i=0; i < fastqFiles.length; i++) {
			String name = fastqFiles[i].getName();
			String baseName = FilenameUtils.removeExtension(name);
			baseName = baseName.substring(0, baseName.length()-2);
			
			Pair<File> pair = new Pair<File>();
			if (hold.containsKey(baseName)) {
				pair.setValue1(hold.get(baseName));
				pair.setValue2(fastqFiles[i]);
				rtrn.put(baseName, pair);
				System.out.println("Found pair: "+pair.getValue1().getName()+" and "+pair.getValue2().getName());
			}
		}
		return rtrn;
	}
	
	private static void markNonspecificAmplification(Map<String, Pair<File>> files, String reference, String referenceIndex, String saveDir) throws InterruptedException, IOException {
		FastaSequenceIndex index = new FastaSequenceIndex(new File(referenceIndex));
		IndexedFastaSequenceFile ref = new IndexedFastaSequenceFile(new File(reference), index);
		
		// Iterate through all pairs of SAM-FASTQ files
		for (Pair<File> pair : files.values()) {
			processOneFilePair(pair, ref, saveDir);
		}
	}
	
	private static void processOneFilePair(Pair<File> pair, IndexedFastaSequenceFile ref, String saveDir) throws InterruptedException {
		FastqReader fq = new FastqReader(pair.getValue2());
		SAMFileReader sam = new SAMFileReader(pair.getValue1());
		SAMRecordIterator samItr = sam.iterator();

		File outFile = new File(FilenameUtils.concat(saveDir, pair.getValue1().getName()));
		System.err.println(outFile.getPath());

		final SAMFileWriter samWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(sam.getFileHeader(), true, outFile);
		//SAMTextWriter samWriter = new SAMTextWriter(outFile);
		//System.err.println(sam.getFileHeader().getTextHeader());
		//samWriter.writeHeader(sam.getFileHeader().getTextHeader());

		// first dimension represents specific (1) / nonspecific (0)
		// second dimension represents + strand (1) / - strand (0) counts
		int[][] counts = { {0,0}, {0,0} };
		int other = 0;

		while (fq.hasNext()) {
			//for (int i=0; i<100; i++) {
			if (!samItr.hasNext()) {
				System.err.println("Unexpectedly reached the end of the SAM file.");
				throw new InterruptedException();
			}

			FastqRecord seq = fq.next();
			
			// Read in records for each of the paired reads
			SAMRecord rec;
			SAMRecord rec1 = samItr.next();		
			SAMRecord rec2 = samItr.next();
			if (!rec1.getFirstOfPairFlag())
				rec = rec2;
			else
				rec = rec1;  			// We want the first of pair

			// This means that the SAM read was mapped successfully; otherwise
			// the name will have a "NM:i:0" tag at the end
			String fastaName = seq.getReadHeader().split(" ")[0];
			if (fastaName.equals(rec.getReadName()) && !rec.getReadUnmappedFlag()) {	
				Pair<Boolean> result = checkOneRead(seq, rec, ref);
				int specific = result.getValue1() ? 1 : 0;
				int strand = result.getValue2() ? 1 : 0;
				counts[specific][strand]++;

				// Update the flag - use the 0x800 bit to record if the read is aligned but is
				// a non-specific RT artifact
				if (specific == 0) {
					rec1.setFlags(rec1.getFlags() + 0x800);
					rec2.setFlags(rec2.getFlags() + 0x800);
				}
			} else {
				other++;
			}

			samWriter.addAlignment(rec1);
			samWriter.addAlignment(rec2);
		}

		samWriter.close();
		System.out.printf("%s:\n S+ %d\n S- %d\n NS+ %d\n NS- %d\n Oth %d\n", pair.getValue1().getName(), counts[1][1], counts[1][0], counts[0][1], counts[0][0], other);
	}
	
	private static Pair<Boolean> checkOneRead(final FastqRecord seq, final SAMRecord rec, final IndexedFastaSequenceFile ref) {		
		// Get the sequence that should precede the read.
		// First boolean represents specific / nonspecific, and second boolean represents strand
		// Returning true means it's a specific read (the adjoining sequence in the aligned
		//   gene does not match the read), and returning false means it's a nonspecific read
		// Second boolean is true if the strand is +
		
		ReferenceSequence refgene = ref.getSequence(rec.getReferenceName());
		Pair<Boolean> result = new Pair<Boolean>();
		
		// Note: In our ligation scheme, specific sequences should also be _negative_ strand only
		long start, end;
		if (rec.getReadNegativeStrandFlag()) {
			//System.out.println("-");
			result.setValue2(false);
			start = rec.getAlignmentEnd() + 2;
			end = start + 7;
			if (end >= refgene.length()) {
				result.setValue1(true);
				return result;		
			}
		} else {   // Should never have a specific read that is positive strand
			//System.out.println("+");
			result.setValue2(true);
			//end = rec.getAlignmentStart() - 2;
			//start = end - 7;
			//if (start < 0) {
			result.setValue1(false);
			return result;
		}
		
		//System.out.printf("%d %d %d %d %d\n",refgene.length(), start, end, rec.getAlignmentStart(), rec.getAlignmentEnd());
		//String extbases = new String(ref.getSubsequenceAt(rec.getReferenceName(), Math.max(0, start-10), Math.min(refgene.length(), end+10)).getBases());
		String bases = new String(ref.getSubsequenceAt(rec.getReferenceName(), start, end).getBases());
		
		/*
		System.out.println("         "+seq.getReadString());
		if (rec.getReadNegativeStrandFlag()) {
			System.out.println(SequenceUtil.reverseComplement(extbases));	
			System.out.println("          "+SequenceUtil.reverseComplement(bases));
		} else {
			System.out.println(extbases);
			System.out.println("          "+bases);
		}*/
		
		result.setValue1(!equalsBarcodeSequence(seq.getReadString(), bases, rec.getReadNegativeStrandFlag()));
		return result;
	}
	
	static private boolean equalsBarcodeSequence(String barcode, String bases, boolean isNegativeStrand) {
		if (isNegativeStrand) {
			bases = SequenceUtil.reverseComplement(bases);
		}
		return barcode.substring(1,8).equals(SequenceUtil.reverseComplement(bases));
	}
	
	static private String usage="\nargs[0]=directory containing untrimmed FastQ files \n args[1]=directory containing Bowtie-aligned SAM files \n args[2]=Alignment reference FASTA \n args[3]=saveDir \n\n ";

	public static void main(String[] args) throws IOException, InterruptedException {
		if (args.length == 5) {
			File[] fastq=new File(args[0]).listFiles(new ExtensionFilter(".fq"));
			File[] sams=new File(args[1]).listFiles(new ExtensionFilter(".sam"));
			String reference = args[2];
			String referenceIndex = args[3];
			String outdir = args[4];
			
			new MarkNonspecificAmplification(fastq, sams, reference, referenceIndex, outdir);
		} else {
			System.err.println(usage);
		}
	}

}
