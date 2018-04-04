package broad.pda.geneexpression.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.core.parser.StringParser;

/**
 * This class makes a duplicate marked bam file using paired information from the fastq files
 * Samtools are not very efficient for 3'DGE mapped data since it can have a lot of polyTs
 * @author skadri
 *
 */
public class MarkDuplicates {

	static Logger logger = Logger.getLogger(MarkDuplicates.class.getName());
	
	private static int PAIRED_COLUMNS = 4;
	private static int UNPAIRED_COLUMNS = PAIRED_COLUMNS - 1;
	/*
	 * Mapping of name to FastQ file paths
	 */
	Map<String,String> leftFqs;
	Map<String,String> rightFqs;
	/**
	 * The sample names
	 */
	TreeSet<String> sampleNames;
	/*
	 * Mapping of name to whether there is paired data
	 */
	Map<String,Boolean> pairedData;
	
	/*
	 * Mappings of names to CURRENT fqs being used. 
	 * This is because if the user filters out rRNA, the name of the fqs to align changes.
	 */
	Map<String,String> currentLeftFqs;
	Map<String,String> currentRightFqs;

	/**
	 * The paths of the current bam files being used
	 * This is because we might run extra alignments and merge bam files
	 */
	Map<String,String> currentBamFiles;
	
	/**
	 * The directory containing the current bam files
	 */
	String currentBamDir;
	
	Map<String,String> nameToCondition;

	String queueName;
	
	public MarkDuplicates(String inputListFile) throws IOException{
		
		this.sampleNames = new TreeSet<String>();
		this.leftFqs = new TreeMap<String,String>();
		this.rightFqs = new TreeMap<String,String>();
		this.currentLeftFqs = new TreeMap<String,String>();
		this.currentRightFqs = new TreeMap<String,String>();
		this.pairedData = new TreeMap<String,Boolean>();
		this.currentBamFiles = new TreeMap<String,String>();
		this.nameToCondition = new TreeMap<String,String>();
		//Read the Fq list
		readFqList(inputListFile);
		
		logger.info("Marking duplicates...");
		// Make directory for duplicate marked bam files
		File dir = new File("duplicate_marked_files");
		boolean madeDir = dir.mkdir();
		for(String sample : this.sampleNames) {
			String read1file = this.leftFqs.get(sample);
			String read2file = this.rightFqs.get(sample);
			Map<String,String> uniqueReads = new HashMap<String,String>();
			Set<String> duplicatedReadNames = new TreeSet<String>();
			
			//READ FASTQ FILES
			FileReader reader1 = new FileReader(read1file);
			BufferedReader buffered1 = new BufferedReader(reader1);
					
			FileReader reader2 = new FileReader(read2file);
			BufferedReader buffered2 = new BufferedReader(reader2);
			
			int linesRead = 0;
			String readName=null;
			//STORE DUPLICATE NAMES IN A LIST
			while(buffered1.ready()) {
				String read1line = buffered1.readLine();
				String read2line = buffered2.readLine();
				linesRead++;
				//Read names
				if(linesRead % 4 == 1) {
					readName = read1line.substring(1,read1line.length()-2);
				}
				if(linesRead % 4 == 2) {
					String combined = read1line + "_" + read2line;
					//System.out.println("Read name is: "+readName);
					if(uniqueReads.containsKey(combined)){ 
						duplicatedReadNames.add(readName);
						duplicatedReadNames.add(uniqueReads.get(combined));
					}
					else{
						uniqueReads.put(combined,readName);
					}
				}
			}
			/*for (String s:duplicatedReadNames){
				System.out.println(s);
			}*/
			//GO THROUGH THE BAM FILE
			
		}
	}
	
	/**
	 * 
	 * @param fileName
	 * @throws IOException 
	 */
	private void readFqList(String fileName) throws IOException{
		
		/*	Input file is a list
		*	Name		Left.fq		Right.fq	condition
		*	If unpaired, only Name	Left.fq	condition
		*	The NAME must be unique
		*/
		//TODO: If the name is not unique, assign a unique name
		logger.info("");
		logger.info("Reading list of fastq files from " + fileName + "...");
		
		FileReader r = new FileReader(fileName);
		BufferedReader b = new BufferedReader(r);
		if(!b.ready()) {
			throw new IllegalArgumentException("File " + fileName + " is empty.");
		}
		StringParser p = new StringParser();
		
		boolean first = true;
		boolean unpaired = true;
		while(b.ready()) {
			String line = b.readLine();
			p.parse(line);
			int cols = p.getFieldCount();
			if(cols == 0) continue;
			
			if(first) {
				if(cols < UNPAIRED_COLUMNS) {
					throw new IllegalArgumentException("Illegal number of columns in " + fileName);
				}
				if(cols == PAIRED_COLUMNS) unpaired = false;
				first = false;
			}
			
			String sampleName = p.asString(0);
			this.sampleNames.add(sampleName);
			String leftReads = p.asString(1);
			this.leftFqs.put(sampleName, leftReads);
			
			if(unpaired) {
				if(cols != UNPAIRED_COLUMNS) {
					throw new IllegalArgumentException("Illegal number of columns in " + fileName);
				}
				this.pairedData.put(sampleName, Boolean.valueOf(false));
				this.nameToCondition.put(sampleName, p.asString(2));
			} else {
				if(cols != PAIRED_COLUMNS) {
					throw new IllegalArgumentException("Illegal number of columns in " + fileName);
				}				
				this.rightFqs.put(sampleName, p.asString(2));
				this.pairedData.put(sampleName, Boolean.valueOf(true));
				this.nameToCondition.put(sampleName, p.asString(3));
			}
			
		}
		
		b.close();
		r.close();
		
		if(first) {
			throw new IllegalArgumentException("File " + fileName + " is invalid.");
		}
		
		this.currentLeftFqs = this.leftFqs;
		this.currentRightFqs = this.rightFqs;
		
		int pa = this.currentRightFqs.isEmpty() ? 0 : this.currentRightFqs.keySet().size();
		int u = this.currentLeftFqs.keySet().size() - pa;
		logger.info("There are " + pa + " sets of paired end reads and " + u + " sets of single end reads.");
		
	}
	
	public static void main(String[] args) throws IOException{
		
		MarkDuplicates dummy = new MarkDuplicates("xyz.list");
	}
}
