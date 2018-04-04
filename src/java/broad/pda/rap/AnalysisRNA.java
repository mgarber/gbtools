/*
 * Jesse Engreitz
 * February 24, 2012
 * 
 * Code for analysis of RNA RAP data - pull down RNA and sequence RNA.
 * Originally written to analyze data in /seq/lincRNA/Jesse/120221_rRNAPulldownMiSeq/
 */

package broad.pda.rap;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import Jama.Matrix;
import broad.core.datastructures.MatrixWithHeaders;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceRecord;

public class AnalysisRNA {

	private static ArrayList<String> statNames = new ArrayList<String>() {{
		this.add("total");			// total read pairs
		this.add("mapped");			// mapped read pairs
		this.add("pcr");			// PCR artifacts marked by Picard
		this.add("nonspecific");	// nonspecific RT reads marked by MarkNonspecificAmplification
		this.add("probe");			// number of + strand reads that map to the target
		this.add("target");			// number of - strand reads that map to the target
	}};
	
	public AnalysisRNA(String in, File out, File statFile) throws IOException {
		Map<String, RAPSample> samples = readSampleSheet(in);
		
		String onekey = samples.keySet().iterator().next();
		Collection<String> refseqIds = getRefseqIds(samples.get(onekey).file);

		MatrixWithHeaders counts = initDataMatrix(refseqIds, samples.keySet());		
		MatrixWithHeaders stats = initDataMatrix(statNames, samples.keySet());
		
		for (String sampleName : samples.keySet()) {
			countFile(counts, stats, sampleName, samples.get(sampleName));
		}
		writeOutput(counts, stats, out, statFile);
	}

	
	/* Assumes that all SAM input files were generated with the same reference */
	private Collection<String> getRefseqIds(File file) {
		Collection<String> refseqIds = new HashSet<String>();
		
		SAMFileHeader header = new SAMFileReader(file).getFileHeader();
		for (SAMSequenceRecord record : header.getSequenceDictionary().getSequences()) {
			refseqIds.add(record.getSequenceName());
		}
		System.err.println("Read "+refseqIds.size()+" RefSeq ids.");
		return refseqIds;
	}
	
	private MatrixWithHeaders initDataMatrix(Collection<String> rowNames, Collection<String> colNames) {
		List<String> rowList = new ArrayList<String>(rowNames);
		List<String> colList = new ArrayList<String>(colNames);
		Matrix data = new Matrix(new double[rowList.size()][colList.size()]);
		MatrixWithHeaders mat = new MatrixWithHeaders(data, rowList, colList);
		return mat;
	}
	
	
	private static void countFile(MatrixWithHeaders counts, MatrixWithHeaders stats, final String sampleName, RAPSample sample) {
		System.err.println("Reading "+sampleName);
		SAMRecordIterator samFile = new SAMFileReader(sample.file).iterator();
		
		
		while (samFile.hasNext()) {
			SAMRecord read = samFile.next();
			if (!read.getFirstOfPairFlag())   // This SAM file is unsorted, so to count each read pair once, let's look at Read 1
				continue;
			
			if (!read.getReadUnmappedFlag()) {
				stats.set("mapped",sampleName, stats.get("mapped", sampleName) + 1);
				
				if ((read.getFlags() & 0x800) == 0) { // only add if the read isn't nonspecific
					String refseq = read.getReferenceName();
					counts.set(refseq, sampleName, counts.get(refseq, sampleName) + 1);
					
					String refseqBase = null;
					try {
						refseqBase = refseq.split("\\|")[3].split("\\.")[0];
					} catch (java.lang.ArrayIndexOutOfBoundsException e) {
						System.err.println(refseq.split("\\|")[0]);
						System.err.println(refseq);
					}
					
					if (refseqBase.equals(sample.probe.refseq)) {
						if (read.getReadNegativeStrandFlag() && !sample.probe.antisense)
							stats.set("target",sampleName, stats.get("target", sampleName) + 1);
						else
							stats.set("probe",sampleName, stats.get("probe", sampleName) + 1);
					}
				}
			}
			if (read.getDuplicateReadFlag()) 
				stats.set("pcr", sampleName, stats.get("pcr", sampleName) + 1);
			if ((read.getFlags() & 0x800) > 0)
				stats.set("nonspecific", sampleName, stats.get("nonspecific", sampleName) + 1);
			
			stats.set("total", sampleName, stats.get("total", sampleName) + 1);
		}
	}
	
	private static Map<String, RAPSample> readSampleSheet(String in) throws IOException {
		Map<String, RAPSample> samFiles = new TreeMap<String, RAPSample>();
		
		Collection<String> lines=BEDFileParser.loadList(in, true);
		for (String line: lines) {
			String[] tokens=line.split("\t");
			RAPProbe probe = new RAPProbe(tokens[2], Integer.parseInt(tokens[3]) == 1);
			String name = tokens[1];
			samFiles.put(name, new RAPSample(probe, new File(tokens[0])));
		}
		
		/*
		for (int i=0; i<in.length; i++) {
			String filename = in[i].getName();
			SAMRecordIterator sam = new SAMFileReader(in[i]).iterator();

			samFiles.put(filename, sam);
		}*/
		
		return samFiles;
	}
	
	
	private void writeOutput(MatrixWithHeaders counts, MatrixWithHeaders stats, File out, File statFile) throws IOException {
		BufferedWriter bwOut = new BufferedWriter(new FileWriter(out));
		BufferedWriter bwStats = new BufferedWriter(new FileWriter(statFile));
		
		writeCounts(counts, bwOut);
		writeCounts(stats, bwStats);
		//counts.write(bwOut);
		//write("row");
		
		bwOut.close();
		bwStats.close();
	}
	
	
	// Write the data matrix as integers
	private void writeCounts(MatrixWithHeaders stats, BufferedWriter bw) throws IOException {
		for(String columnName : stats.getColumnNames()) {
			bw.write("\t");
			bw.write(columnName);
		}
		bw.newLine();
		bw.flush();
		for(String rowName : stats.getRowNames()) {
			bw.write(rowName);
			for(String colName : stats.getColumnNames()) {
				bw.write("\t");
				bw.write(String.valueOf((int)stats.get(rowName, colName)));
			}
			bw.newLine();
		}
	}
	
	private static class RAPSample {
		public RAPProbe probe;
		File file;
		
		public RAPSample(RAPProbe probe, File file) {
			this.probe = probe;
			this.file = file;
		}
	}
	
	public static class RAPProbe {
		public String refseq;
		public boolean antisense;
		
		public RAPProbe(String refseq, boolean antisense) {
			this.refseq = refseq;
			this.antisense = antisense;
		}
	}
	
	
	private static String USAGE = "AnalysisRNA parameters:\n\n\t-in\tInput directory\n\t-out\tBasename for output files\n"+
								   "\t-filter\tOptional suffix filter for input files\n\t-targets\tTabbed file containing information about the RNA targets\n\n";
	
	public static void main (String [] args) throws IllegalArgumentException, IOException {
		ArgumentMap argmap = CLUtil.getParameters(args, USAGE, "analyze");
		
		final String in = argmap.getInput();
		final String out = argmap.getOutput();
		
		/*
		String filter = argmap.get("filter");
		String targets = argmap.get("targets");
		
		File[] inFiles;
		if (filter == null) {
			inFiles = new File(in).listFiles();
		} else {
			inFiles = new File(in).listFiles(new ExtensionFilter(filter));
		}
		
		if (inFiles.length == 0) {
			System.err.println("No input files were found.");
			throw new IOException();
		}
		for (int i=0;i<inFiles.length;i++)
			System.err.println(inFiles[i].getName());
			
		File targetFile;
		if (targets != null)
			targetFile = new File(targets);
		else
			targetFile = null;
		*/
		
		File outFile = new File(out + ".txt");
		File statFile = new File(out + ".stats.txt");
				
		new AnalysisRNA(in, outFile, statFile);
	}
	
	
}
