/**
 * 
 */
package broad.pda.countreads;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import net.sf.picard.sam.DuplicationMetrics;

/**
 * @author prussell
 *
 */
public class DuplicatesAndLibrarySizeAnalysis {

	/**
	 * Constructor for single end reads
	 * @param singleReadFastq
	 * @param log Logger object
	 * @throws IOException 
	 */
	public DuplicatesAndLibrarySizeAnalysis(String singleReadFastq, Logger log) throws IOException {
		this(singleReadFastq, null, log);
	}
	
	/**
	 * Constructor for paired end reads
	 * @param read1fastq
	 * @param read2fastq
	 * @param log Logger object
	 * @throws IOException 
	 */
	public DuplicatesAndLibrarySizeAnalysis(String read1fastq, String read2fastq, Logger log) throws IOException {
		this.read1file = read1fastq;
		this.read2file = read2fastq;
		this.uniqueReads = new TreeSet<String>();
		this.duplicatedReads = new TreeSet<String>();
		this.estLibrarySize = -1;
		this.numUniqueReads = -1;
		this.totalReads = -1;
		this.pctDup = -1;
		this.logger = log;
		this.readsPaired = false;
		if(read2fastq != null) this.readsPaired = true;
		collapseReadsAndCount();
	}
	
	private String read1file;
	private String read2file;
	private TreeSet<String> uniqueReads;
	private TreeSet<String> duplicatedReads;
	private int numUniqueReads;
	private double pctDup;
	private int totalReads;
	private long estLibrarySize;
	private boolean readsPaired;
	private Logger logger;
	
	/**
	 * Get the number of unique reads
	 * @return The number of unique reads
	 */
	public int getNumUniqueReads() {
		return this.numUniqueReads;
	}
	
	/**
	 * Get the proportion of reads that are duplicates
	 * @return The proportion of reads that are duplicated
	 */
	public double getPercentDuplicated() {
		return this.pctDup;
	}
	
	/**
	 * Get the total number of reads
	 * @return The total number of reads
	 */
	public int getTotalReads() {
		return this.totalReads;
	}
	
	/**
	 * Get the estimated library size based on number of reads and percent duplicates
	 * @return The estimated library size
	 */
	public long getEstimatedLibrarySize() {
		return this.estLibrarySize;
	}
	
	/**
	 * Count unique reads by collapsing identical reads
	 * @param outFile Output file for statistics
	 * @param readsPaired Whether the reads are paired
	 * @throws IOException
	 */
	private void collapseReadsAndCount() throws IOException {
		
		// If paired reads, read read2 file. Otherwise just read read1 file.
		String read2 = this.readsPaired ? this.read2file : this.read1file;
		
		FileReader reader1 = new FileReader(this.read1file);
		BufferedReader buffered1 = new BufferedReader(reader1);
				
		FileReader reader2 = new FileReader(read2);
		BufferedReader buffered2 = new BufferedReader(reader2);
		
		int linesRead = 0;
		this.totalReads = 0;
		
		while(buffered1.ready()) {
			String read1line = buffered1.readLine();
			String read2line = buffered2.readLine();
			linesRead++;
			if(linesRead % 4 == 2) {
				this.totalReads++;
				String combined = read1line + "_" + read2line;
				if(this.uniqueReads.contains(combined)) this.duplicatedReads.add(combined);
				this.uniqueReads.add(combined);
			}
		}
		
		buffered1.close();
		buffered2.close();
		
		this.numUniqueReads = this.uniqueReads.size();
		try {
			this.estLibrarySize = DuplicationMetrics.estimateLibrarySize(this.totalReads, this.numUniqueReads).longValue();
		} catch (NullPointerException e) {
			this.logger.info("Warning: caught NullPointerException. Total reads = " + this.totalReads + ". Unique reads = " + this.numUniqueReads + ".");
		}
		this.pctDup = ((double)this.totalReads - (double)this.numUniqueReads)/this.totalReads;
	}
	
	/**
	 * Write fastq files of unique reads and duplicated reads
	 * @param outUniquePrefix Output fastq file of unique reads
	 * @param outDupPrefix Output fastq file of duplicated reads
	 * @param pairedReads Whether the reads are paired
	 * @throws IOException 
	 */
	public void writeSeparateFiles(String outUniquePrefix, String outDupPrefix, boolean pairedReads) throws IOException {
		
		if(this.duplicatedReads.isEmpty()) {
			throw new IllegalStateException("Duplicated read pair set is empty. Try calling collapseReads() first.");
		}
		
		// If paired reads, read read2 file. Otherwise just read read1 file.
		String read2;
		if(pairedReads) read2 = this.read2file;
		else read2 = this.read1file;
		
		FileReader reader1 = new FileReader(this.read1file);
		FileReader reader2 = new FileReader(read2);
		BufferedReader buffered1 = new BufferedReader(reader1);
		BufferedReader buffered2 = new BufferedReader(reader2);
		
		if(pairedReads) {
		
			FileWriter ou1 = new FileWriter(outUniquePrefix + "_1.fq");
			FileWriter ou2 = new FileWriter(outUniquePrefix + "_2.fq");
			FileWriter od1 = new FileWriter(outDupPrefix + "_1.fq");
			FileWriter od2 = new FileWriter(outDupPrefix + "_2.fq");
		
			this.logger.info("Writing unique reads to " + outUniquePrefix + "_1.fq" + " and " + outUniquePrefix + "_2.fq" + ".");
			this.logger.info("Writing duplicated reads to " + outDupPrefix + "_1.fq" + " and " + outDupPrefix + "_2.fq" + ".");
			
			while(buffered1.ready()) {
				String read1Line1 = buffered1.readLine();
				String read2Line1 = buffered2.readLine();
				String read1Line2 = buffered1.readLine();
				String read2Line2 = buffered2.readLine();
				String read1Line3 = buffered1.readLine();
				String read2Line3 = buffered2.readLine();
				String read1Line4 = buffered1.readLine();
				String read2Line4 = buffered2.readLine();
			
				String combined = read1Line2 + "_" + read2Line2;
				if(this.duplicatedReads.contains(combined)) {
					od1.write(read1Line1 + "\n");
					od1.write(read1Line2 + "\n");
					od1.write(read1Line3 + "\n");
					od1.write(read1Line4 + "\n");
					od2.write(read2Line1 + "\n");
					od2.write(read2Line2 + "\n");
					od2.write(read2Line3 + "\n");
					od2.write(read2Line4 + "\n");
				} else {
					ou1.write(read1Line1 + "\n");
					ou1.write(read1Line2 + "\n");
					ou1.write(read1Line3 + "\n");
					ou1.write(read1Line4 + "\n");
					ou2.write(read2Line1 + "\n");
					ou2.write(read2Line2 + "\n");
					ou2.write(read2Line3 + "\n");
					ou2.write(read2Line4 + "\n");
				}
			
			}
		
			buffered1.close();
			buffered2.close();
			ou1.close();
			ou2.close();
			od1.close();
			od2.close();
		
		}
		
		else {
			
			FileWriter ou1 = new FileWriter(outUniquePrefix + ".fq");
			FileWriter od1 = new FileWriter(outDupPrefix + ".fq");
		
			this.logger.info("Writing unique reads to " + outUniquePrefix + ".fq.");
			this.logger.info("Writing duplicated reads to " + outDupPrefix + ".fq.");
			
			while(buffered1.ready()) {
				String read1Line1 = buffered1.readLine();
				String read2Line1 = buffered2.readLine();
				String read1Line2 = buffered1.readLine();
				String read2Line2 = buffered2.readLine();
				String read1Line3 = buffered1.readLine();
				String read2Line3 = buffered2.readLine();
				String read1Line4 = buffered1.readLine();
				String read2Line4 = buffered2.readLine();
			
				String combined = read1Line2 + "_" + read2Line2;
				if(this.duplicatedReads.contains(combined)) {
					od1.write(read1Line1 + "\n");
					od1.write(read1Line2 + "\n");
					od1.write(read1Line3 + "\n");
					od1.write(read1Line4 + "\n");
				} else {
					ou1.write(read1Line1 + "\n");
					ou1.write(read1Line2 + "\n");
					ou1.write(read1Line3 + "\n");
					ou1.write(read1Line4 + "\n");
				}
			
			}
		
			buffered1.close();
			buffered2.close();
			ou1.close();
			od1.close();

		}
		
	}

}
