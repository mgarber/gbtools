package broad.pda.seq.alignment;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DateFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.core.util.CollapseByIntersection;
import broad.core.util.PipelineUtils;
import broad.pda.alignment.PairedEndAlignment;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.GeneTools;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.alignment.sam.SAMRecord;
import broad.pda.seq.alignment.sam.SAMUtils;
import broad.pda.seq.fastq.FastqParser;
import broad.pda.seq.graph.ChromosomeWithBubbles2;
import broad.pda.seq.graph.Path;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

public class RNASeqAlignmentPipeline {

	//TODO Have optional GSNAP step
	//TODO Fix the paired end counting for reporting
	//TODO Cleanup all the extra files
	//TODO Report proper paired SAM file
	
	//String blatScript="/seq/mguttman/scripts/BLATPipeline/BLATPipeline.jar";
	private static final Logger logger =  Logger.getLogger(RNASeqAlignmentPipeline.class.getName());
	private static String blatScriptClass="broad.pda.seq.alignment.BLATPipeline";
	private static String bowtie="/seq/mguttman/scripts/TopHat/bin/bowtie";
	private static String tophat="/seq/mguttman/scripts/TopHat/bin/tophat";
	private static String TMP_DIR = "/broad/shptmp/scripturetmp/";
	String queue="week";
	int waitTime=60000; //1 minute
	int chunkSize=1000000;
	private static final int MAX_GAP = 100000;
	
	File novelJunctionSequenceFile;
	File novelJunctionSequenceMap;
	String novelJunctionBowtieIndex;
	
	File knownJunctionSequenceMap;
	String knownJunctionBowtieIndex;
	
	String bowtieGenomeIndex;
	String genomeDir;
	String saveDir;
	
	Pair<File> fastqFiles;
	Pair<Collection<File>> blatOutputDirectories;
	Pair<Collection<File>> tophatOutputDirectories;
	
	int blatMinScore=50;
	int blatMinPercentIdentity=93;
	private int readLength=76;
	boolean hasPairs=false;
	boolean filterBySplice=true;
	int minIntronSize=30;
	double minPercentMapped=.9;
	int minCoverageOnEachSide=10;
	
	Pair<String> alignedGenome;
	Pair<String> alignedKnown;
	Pair<String> alignedNovel;
	Pair<File> tophatAligned;
	Pair<File> blatAligned;
	Pair<File> bowtieUniqueMappers;
	Pair<File> bowtieNonUniqueMappers;
	Pair<File> salvagedReads;
	String save;
	
	int maxNumberOfPlacements=40; //Tophat default
	private int maxPairDistance=1000000;
	private int minPairDistance=50;
	
	private boolean runTophat;
	private boolean useKnownJunctions;
	
	//For reporting
	private Map<File, Integer> readCounts;
	private Map<String, Integer> novelJunctionCounts;
	private Pair<Integer> salvagedCounts;
	private Pair<Integer> bowtieCounts;
	private int uniqueCounter;
	private int multimapperCounter;
	private int multimappingPairs;
	private int uniquePairs;
	private int both;
	private int spliced;
	private int unspliced;
	private int multimapperCounterPair2;
	private int uniqueCounterPair2;
	private int multimapperCounterPair1;
	private int uniqueCounterPair1;
	private int bothPair2;
	private int splicedPair2;
	private int unsplicedPair2;
	private int unsplicedPair1;
	private int splicedPair1;
	private int bothPair1;
	private int alignedCount1;
	private int alignedCount2;
	private Pair<String> bowtieAligned;
	private int maxAlignmentGap;
	private String pairedEndInsert;
	private boolean onlySplicingGaps;
	private Pair<File[]> fastqChunks;
	private boolean rechunk;
	
	Map<String, Integer> chromosomeSizes;
	
	private static final String quote="\"";
	
	public RNASeqAlignmentPipeline(Pair<File> fastqFiles, String genomeDir, String saveDir, String bowtieGenomeIndex, String knownJunctionIndex, File knownJunctionMap, String save, Pair<File> tophatAligned, Pair<File> blatAligned, boolean runTophat, boolean useKnownJunctions, Pair<String> bowtieAligned, int maxAlignmentGapSize, boolean onlySplicinGaps, boolean rechunk) throws IOException, InterruptedException{
		//Initialization
		this.fastqFiles=fastqFiles;
		if(fastqFiles.hasValue2()){this.hasPairs=true;}
		this.bowtieAligned=bowtieAligned;
		this.blatAligned=blatAligned;
		this.tophatAligned=tophatAligned;
		this.knownJunctionBowtieIndex=knownJunctionIndex;
		this.knownJunctionSequenceMap=knownJunctionMap;
		this.bowtieGenomeIndex=bowtieGenomeIndex;
		this.genomeDir=genomeDir;
		File sizeFileFile = new File(genomeDir + "/sizes");
		if(sizeFileFile.exists()) {
			chromosomeSizes = BEDFileParser.loadChrSizes(sizeFileFile.getAbsolutePath());
		} else {
			chromosomeSizes = new HashMap<String, Integer>();
		}
		
		this.saveDir=saveDir;
		this.save=save;
		this.runTophat=runTophat;
		this.useKnownJunctions=useKnownJunctions;
		this.maxAlignmentGap = maxAlignmentGapSize;
		this.onlySplicingGaps = onlySplicinGaps;
		this.rechunk = rechunk;
		
		//For reporting
		this.readCounts=new TreeMap<File, Integer>();
		this.novelJunctionCounts=new TreeMap<String, Integer>();
		
		File tmpdir = new File(TMP_DIR);
		if(!tmpdir.exists()) {
			tmpdir.mkdir();
		}
	}
	
	public RNASeqAlignmentPipeline(String saveDir) {
		this.saveDir = saveDir;
		
		File tmpdir = new File(TMP_DIR);
		if(!tmpdir.exists()) {
			tmpdir.mkdir();
		}
	}
	
	public void makePairs(Pair<File> alignmentPair) throws IOException, InterruptedException {
		hasPairs = true;
		String mergedAlignment = saveDir +"/merged.sam";
		String pairedInsertAlignment = saveDir+"/pairedInsert.sam";
		resolveMapping(alignmentPair, mergedAlignment, pairedInsertAlignment);
		writeToLog("DONE mapping paires");
		//Step 7: Sort and index
		System.err.println("Sorting and indexing...");
		
		Runtime run=java.lang.Runtime.getRuntime();
		sortAndIndex( mergedAlignment, run);
		sortAndIndex( pairedInsertAlignment, run);
	}
	
	public void runAlignment(boolean spliceBySplice) throws IOException, InterruptedException{
		Runtime run=java.lang.Runtime.getRuntime();
		writeToLog("START OF RUN " + DateFormat.getDateTimeInstance().format(new Date()));
		//Step 0: Split files into chunks
		fastqChunks=new Pair<File[]>();
		if(this.blatAligned==null || (this.tophatAligned==null && runTophat)){
			writeToLog("chunking ");			
			System.err.println("Blat aligning? "+blatAligned);
			fastqChunks.setValue1(splitIntoChunks(fastqFiles.getValue1()));
			if(fastqFiles.hasValue2()){
				fastqChunks.setValue2(splitIntoChunks(fastqFiles.getValue2()));
			}
			
			writeToLog("DONE chunking");
			writeToLog("Pair 1 chunks",fastqChunks.getValue1());
			writeToLog("Pair 2 chunks", fastqChunks.getValue2());
		}

		if(this.bowtieAligned==null){
			//Step 1: Get putative junctions
			System.err.println("Find novel junctions...");
			findNovelJunctions(fastqChunks, run);
			writeToLog("DONE findNovelJunctions");
			//Step 2: Align all reads to novel junctions, known junctions, and genome
			System.err.println("Align all reads to novel junctions using Bowtie...");
			this.bowtieAligned=alignReadsWithBowtie(fastqFiles, run);
			writeToLog("DONE alignReadsWithBowtie " + DateFormat.getDateTimeInstance().format(new Date()));
		}
		
		//Step 3: Salvage reads aligned in step 1 but not the remainder
		System.err.println("Salvage reads not aligned by Bowtie...");
		Pair<File> salvagedReads=salvageReads(bowtieAligned);
		writeToLog("DONE salvageReads " + DateFormat.getDateTimeInstance().format(new Date()));
		
		//Step 4: Merge all bowtie and salvaged reads
		System.err.println("Merge Bowtie aligned and salvaged reads...");
		Pair<File> allAligned=mergeBowtieAndSalvaged(bowtieAligned, salvagedReads, run);
		writeToLog("DONE mergeBowtieAndSalvaged " + DateFormat.getDateTimeInstance().format(new Date()));
		
		//Step 5: Resolve unique mappers, unique pairs, and lowest mismatch pairs
		System.err.println("Resolving mapping based on pairs...");
		//TODO Add a tag to flag multimapping junctions vs multimapping non-junctions
		resolveMapping(allAligned, save, save+".PairedInsert.sam");
		writeToLog("DONE mapping paires");
		//Step 7: Sort and index
		System.err.println("Sorting and indexing...");
		sortAndIndex( save, run);
		sortAndIndex( save+".PairedInsert.sam", run);
		writeToLog("DONE Sort & Index " + DateFormat.getDateTimeInstance().format(new Date()));
		//Step 8: Compute basic statistics about the run
		alignmentReporting();
		writeToLog("RUN FINISHED WITHOUT KNOWN ERRORS " + DateFormat.getDateTimeInstance().format(new Date()));
		//Step 8: Try to resolve multimappers using paired ends
		
		//Step 9: Flag poly-A reads and their mates and partial polyA reads
	}
	
	

	private void writeToLog(String string, Object[] values)throws IOException  {
		if(values == null) { return;}
		BufferedWriter bw = new BufferedWriter(new  FileWriter(saveDir+"/log",true));
		bw.write(string);
		bw.newLine();
		for(Object val : values) {
			if(val != null) {
				bw.write(val.toString());
				bw.newLine();
			}
		}
		bw.close();
		
	}

	private void writeToLog(String string) throws IOException {
		BufferedWriter bw = new BufferedWriter(new  FileWriter(saveDir+"/log",true));
		bw.write(string);
		bw.newLine();
		bw.close();
	}

	private void alignmentReporting() {
		//Number of reads
		for(File file: readCounts.keySet()){
			System.out.println("Read Counts\t"+file+"\t"+readCounts.get(file));
		}
		
		//Number of novel junctions identified by BLAT and Tophat
		//Number of novel junctions identified by BLAT
		//Number of novel junctions identified by Tophat
		for(String name: this.novelJunctionCounts.keySet()){
			System.out.println("Novel Junction Counts: "+name+"\t"+this.novelJunctionCounts.get(name));
		}
		
		//Number of aligned reads by Bowtie
		System.out.println("Bowtie aligned\t"+(this.alignedCount1-this.salvagedCounts.getValue1())+"\t"+(this.alignedCount2-this.salvagedCounts.getValue2()));
		
		//Number of salvaged
		System.out.println("Salvaged aligned\t"+this.salvagedCounts.getValue1()+"\t"+this.salvagedCounts.getValue2());
		
		//Total number of aligned reads
		System.out.println("Total aligned reads\t"+this.alignedCount1+"\t"+this.alignedCount2);
		
		//Percent of reads that are unique alignments
		if(!this.hasPairs){
			System.out.println("Uniquely aligned reads-Single\t"+this.uniqueCounter);
			System.out.println("Multimapping aligned reads-Single\t"+this.multimapperCounter);
		}
		
		if(this.hasPairs){
			System.out.println("Uniquely aligned reads-Paired\t"+this.uniqueCounterPair1+"\t"+this.uniqueCounterPair2);
			System.out.println("Multimapping aligned reads-Paired\t"+this.multimapperCounterPair1+"\t"+this.multimapperCounterPair2);
			System.out.println("Unique pairs\t"+this.uniquePairs);
			System.out.println("Multimapping pairs\t"+this.multimappingPairs);
		}
		
		//Percent of reads that are spliced
		if(!this.hasPairs){
			System.out.println("Spliced reads-Single\t"+this.spliced);
			System.out.println("Unspliced reads-Single\t"+this.unspliced);
			System.out.println("Ambiguous spliced/unspliced-Single\t"+this.both);
		}
		
		if(this.hasPairs){
			System.out.println("Spliced reads-Paired\t"+this.splicedPair1+"\t"+this.splicedPair2);
			System.out.println("Unspliced reads-Paired\t"+this.unsplicedPair1+"\t"+this.unsplicedPair2);
			System.out.println("Ambiguous spliced/unspliced-Paired\t"+this.bothPair1+"\t"+this.bothPair2);
		}
	}

	//These alignment files are assumed to be sorted by name
	private void resolveMapping(Pair<File> allAligned, String save, String pairedAlignmentSave) throws IOException {		
		if(!this.hasPairs){
			//just make unique
			getSingle(allAligned.getValue1(), save);
		}
		else{
			//Step 1: Get all pairs
			getPairs(allAligned, save, pairedAlignmentSave);		
		}
	}

	private Pair<File> getSingle(File source, String save) throws IOException {
		FileWriter writer1=new FileWriter(save);
		
		Iterator<Collection<SAMRecord>> pair1Iter=new SortedSAMIterator(source);
		
		this.uniqueCounter=0;
		this.multimapperCounter=0;
		while(pair1Iter.hasNext()){
			Collection<SAMRecord> pair1=pair1Iter.next();
			printSingleMapping(pair1, writer1);
			if(pair1.size()==1){uniqueCounter++;}
			else{multimapperCounter++;}
			alignedCount1++;
		}
		
		writer1.close();
		return new Pair<File> (new File(save), null);
	}

	private Pair<File> getPairsOld(Pair<File> allAligned, String save) throws IOException {
		String save1=save+".Pair1.sam";
		String save2=save+".Pair2.sam";
		String saveInsert=save+".PairedInsert.sam";
		FileWriter writer1=new FileWriter(save1);
		FileWriter writer2=new FileWriter(save2);
		FileWriter writerPairs=new FileWriter(saveInsert);
		this.pairedEndInsert=saveInsert;
		
		//Given 2 files that are sorted get the pairing relationship
		Iterator<Collection<SAMRecord>> pair1Iter=new SortedSAMIterator(allAligned.getValue1());
		Iterator<Collection<SAMRecord>> pair2Iter=new SortedSAMIterator(allAligned.getValue2());
		
		System.err.println("File 1 "+allAligned.getValue1());
		System.err.println("File 2 "+allAligned.getValue2());
		
		boolean done=false;
		
		Collection<SAMRecord> pair1=new TreeSet();
		Collection<SAMRecord> pair2=new TreeSet();
		boolean started=false;
		while((pair1Iter.hasNext() || pair2Iter.hasNext()) && !done){
			if(!started){
				pair1=pair1Iter.next();
				System.err.println("Pair1 "+pair1);
				pair2=pair2Iter.next();
				System.err.println("Pair2 "+pair2);
				this.alignedCount1++;
				this.alignedCount2++;
				//System.err.println(getName(pair1)+" "+getName(pair2));
				started=true;
			}
			String name1=getName(pair1);
			String name2=getName(pair2);
						
			if(pair1!=null && pair1.size()>1){System.err.println("1 "+name1+" "+pair1.size());}
			else{System.err.println("1 is Null "+name1);}
			if(pair2!=null && pair2.size()>1){System.err.println("2 "+name2+" "+pair2.size());}
			else{System.err.println("2 is Null "+name2);}
			if(name1.equalsIgnoreCase(name2)){
				//valid pair
				//add both
				printPairMappingOld(pair1, pair2, writer1, writer2, writerPairs);
				//proceed both iters
				if(pair1Iter.hasNext()){pair1=pair1Iter.next(); this.alignedCount1++;}
				if(pair2Iter.hasNext()){pair2=pair2Iter.next(); this.alignedCount2++;}
			}
			else{
				//2 possibilities
				if(name1.compareTo(name2)<0){
					//then read1 is before read2
					//make pair with read1 and null
					printSingleMapping(pair1, writer1);
					//advance read1 
					if(!pair1Iter.hasNext()){done=true;}
					pair1=pair1Iter.next();
					this.alignedCount1++;
				}
				else{
					//then read1 is after read2
					//make pair with null and read2
					printSingleMapping(pair2, writer2);
					//advance read2 
					if(!pair2Iter.hasNext()){done=true;}
					pair2=pair2Iter.next();
					this.alignedCount2++;
				}
			}
		}
		writer1.close();
		writer2.close();
		writerPairs.close();
		return new Pair<File> (new File(save1), new File(save2));
	}
	
	private void getPairs(Pair<File> allAligned, String save, String pairedAlignmentSave) throws IOException {
		String saveAln=save;
		String saveInsert= pairedAlignmentSave;
		FileWriter alnWriter=new FileWriter(saveAln);
		FileWriter writerPairs=new FileWriter(saveInsert);
		this.pairedEndInsert=saveInsert;
		
		//Given 2 files that are sorted get the pairing relationship
		Iterator<Collection<SAMRecord>> pair1Iter=new SortedSAMIterator(allAligned.getValue1());
		Iterator<Collection<SAMRecord>> pair2Iter=new SortedSAMIterator(allAligned.getValue2());
		
		System.err.println("File 1 "+allAligned.getValue1());
		System.err.println("File 2 "+allAligned.getValue2());
		
		boolean done=false;
		
		Collection<SAMRecord> pair1=new TreeSet<SAMRecord>();
		Collection<SAMRecord> pair2=new TreeSet<SAMRecord>();
		boolean started=false;
		while((pair1Iter.hasNext() || pair2Iter.hasNext()) && !done){
			if(!started){
				pair1=pair1Iter.next();
				//System.err.println("Pair1 "+pair1);
				pair2=pair2Iter.next();
				//System.err.println("Pair2 "+pair2);
				this.alignedCount1++;
				this.alignedCount2++;
				//System.err.println(getName(pair1)+" "+getName(pair2));
				started=true;
			}
			String name1=getName(pair1);
			String name2=getName(pair2);
						
			//if(pair1!=null && pair1.size()>1){System.err.println("1 "+name1+" "+pair1.size());}
			//else{System.err.println("1 is Null "+name1);}
			//if(pair2!=null && pair2.size()>1){System.err.println("2 "+name2+" "+pair2.size());}
			//else{System.err.println("2 is Null "+name2);}
			if(name1.equalsIgnoreCase(name2)){
				//valid pair
				//add both
				//System.err.println("BINGO. got match reads: " + name1);
				printPairMapping(pair1, pair2, alnWriter, writerPairs);
				//proceed both iters
				if(pair1Iter.hasNext()){pair1=pair1Iter.next(); this.alignedCount1++;}
				if(pair2Iter.hasNext()){pair2=pair2Iter.next(); this.alignedCount2++;}
			}
			else{
				//2 possibilities
				if(name1.compareTo(name2)<0){
					//then read1 is before read2
					//make pair with read1 and null
					printSingleMapping(pair1, alnWriter);
					//advance read1 
					if(!pair1Iter.hasNext()){done=true;}
					pair1=pair1Iter.next();
					this.alignedCount1++;
				}
				else{
					//then read1 is after read2
					//make pair with null and read2
					printSingleMapping(pair2, alnWriter);
					//advance read2 
					if(!pair2Iter.hasNext()){done=true;}
					pair2=pair2Iter.next();
					this.alignedCount2++;
				}
			}
		}
		alnWriter.close();
		writerPairs.close();
	}

	private String getName(Collection<SAMRecord> pair) {
		if(pair==null || pair.isEmpty()){return "";}
		else{return pair.iterator().next().getName();}
	}

	private void printSingleMapping(Collection<SAMRecord> pair1, FileWriter writer1) throws IOException{
		int spliced=0;
		int unspliced=0;
		
		Collection<SAMRecord> bestRecords=getBest(pair1);
		
		for(SAMRecord sam: bestRecords){
			//if(sam!=null){
				if(hasPairs) {
					sam.setPairedSequenceFlag();
					sam.setMateUnmappedFlag();
				}
				sam.setWeight(1.0/bestRecords.size());
				if(bestRecords.size()>1){sam.setMappingQuality(0);}
				else{sam.setMappingQuality(255);}
				writer1.write(sam+"\n");
			//}
				if(sam.getGene().getNumExons()>1){spliced++;}
				else{unspliced++;}
		}
		if(spliced >0 && unspliced>0){this.both++;}
		else if(spliced>0){this.spliced++;}
		else{this.unspliced++;}
	}
	
	private Collection<SAMRecord> getBest(Collection<SAMRecord> pair1) {
		Collection<Pair<SAMRecord>> rtrn=new HashSet<Pair<SAMRecord>>();
		
		for(SAMRecord record: pair1){
			Pair<SAMRecord> pair=new Pair<SAMRecord>(record, null);
			rtrn.add(pair);
		}
		
		Collection<Pair<SAMRecord>> best=this.getBestPairs(rtrn);
		Collection<SAMRecord> bestRecords=new TreeSet<SAMRecord>();
		for(Pair<SAMRecord> pair: best){bestRecords.add(pair.getValue1());}
		return bestRecords;
	}

	private void printPairMapping(Collection<SAMRecord> pair1,Collection<SAMRecord> pair2, FileWriter alnWriter, FileWriter writerPairs) throws IOException {
		if(pair1==null && pair2==null){return;}
		else if(pair1==null){
			printSingleMapping(pair2, alnWriter);
			return;
		}
		else if(pair2==null){
			reverseOrientation(pair2); //Reverse pair orientation
			printSingleMapping(pair1, alnWriter);
			return;
		}
		
		
		//go through pairs
		Collection<Pair<SAMRecord>> pairs=getBestPairMappings(pair1, pair2); //make sure this gets the "best" pairings
		
		//No valid pairs then treat the two sets as singlets
		if(pairs.isEmpty()){
			printSingleMapping(pair1, alnWriter);
			reverseOrientation(pair2); //Reverse pair orientation
			printSingleMapping(pair2, alnWriter);
		}
		else{
			for(Pair<SAMRecord> pair: pairs){
				SAMRecord record1=pair.getValue1();
				record1.setFirstOfPairFlag();
				record1.getProperlyPairedMappedFlag();
				SAMRecord record2=pair.getValue2();
				if(record2.getReadNegativeStrandFlag()) {
					record1.setMateNegativeStrandFlag();
				}
				record2.setSecondOfPairFlag();
				record2.getProperlyPairedMappedFlag();
				if(record1.getReadNegativeStrandFlag()) {
					record2.setMateNegativeStrandFlag();
				}
				record2.reverseOrientation();
				SAMRecord insert=new SAMRecord(pair.getValue1().getName(), new RefSeqGene(new PairedEndAlignment(pair).getInsertAlignment()), "*");
				record1.setWeight(1.0/pairs.size());
				record2.setWeight(1.0/pairs.size());
				insert.setWeight(1.0/pairs.size());
				if(pairs.size()>1){
					record1.setMappingQuality(0);
					record2.setMappingQuality(0);
					insert.setMappingQuality(0);
				}
				else{
					record1.setMappingQuality(255);
					record2.setMappingQuality(255);
					insert.setMappingQuality(255);
				}
				alnWriter.write(record1+"\n");
				alnWriter.write(record2+"\n");
				writerPairs.write(insert+"\n");
			}
		}
	}
	private void reverseOrientation(Collection<SAMRecord> records) {
		for(SAMRecord sam : records) {
			sam.reverseOrientation();
		}
		
	}

	private void printPairMappingOld(Collection<SAMRecord> pair1,Collection<SAMRecord> pair2, FileWriter writer1, FileWriter writer2, FileWriter writerPairs) throws IOException {
		if(pair1==null && pair2==null){return;}
		else if(pair1==null){
			printSingleMapping(pair2, writer2);
			return;
		}
		else if(pair2==null){
			printSingleMapping(pair1, writer1);
			return;
		}
		
		
		//go through pairs
		Collection<Pair<SAMRecord>> pairs=getBestPairMappings(pair1, pair2); //make sure this gets the "best" pairings
		
		//No valid pairs then treat the two sets as singlets
		if(pairs.isEmpty()){
			printSingleMapping(pair1, writer1);
			printSingleMapping(pair2, writer2);
		}
		else{
			for(Pair<SAMRecord> pair: pairs){
				SAMRecord record1=pair.getValue1();
				SAMRecord record2=pair.getValue2();
				SAMRecord insert=new SAMRecord(pair.getValue1().getName(), new RefSeqGene(new PairedEndAlignment(pair).getInsertAlignment()), "*");
				record1.setWeight(1.0/pairs.size());
				record2.setWeight(1.0/pairs.size());
				insert.setWeight(1.0/pairs.size());
				if(pairs.size()>1){
					record1.setMappingQuality(0);
					record2.setMappingQuality(0);
					insert.setMappingQuality(0);
				}
				else{
					record1.setMappingQuality(255);
					record2.setMappingQuality(255);
					insert.setMappingQuality(255);
				}
				writer1.write(record1+"\n");
				writer2.write(record2+"\n");
				writerPairs.write(insert+"\n");
			}
		}
		
		
		//go through pairs
		/*Pair<Collection<SAMRecord>> alignments=getPairMappings(pair1, pair2);
		
		if(alignments.getValue1().size()==1){this.uniqueCounterPair1++;}
		else{this.multimapperCounterPair1++;}
		
		int spliced=0;
		int unspliced=0;
		for(SAMRecord sam: alignments.getValue1()){
			if(sam!=null){
				sam.setWeight(1.0/alignments.getValue1().size());
				if(alignments.getValue1().size()>1){sam.setMappingQuality(0);}
				writer1.write(sam+"\n");
				if(sam.getGene().getNumExons()>1){spliced++;}
				else{unspliced++;}
			}
		}
		if(spliced>0 && unspliced>0){bothPair1++;}
		else if(spliced>0){splicedPair1++;}
		else{unsplicedPair1++;}
	
		if(alignments.getValue2().size()==1){this.uniqueCounterPair2++;}
		else{this.multimapperCounterPair2++;}
		
		spliced=0;
		unspliced=0;
		for(SAMRecord sam: alignments.getValue2()){
			if(sam!=null){
				sam.setWeight(1.0/alignments.getValue2().size());
				if(alignments.getValue2().size()>1){sam.setMappingQuality(0);}
				writer2.write(sam+"\n");
				if(sam.getGene().getNumExons()>1){spliced++;}
				else{unspliced++;}
			}
		}
		if(spliced>0 && unspliced>0){bothPair2++;}
		else if(spliced>0){splicedPair2++;}
		else{unsplicedPair2++;}
		
		if(alignments.getValue1()!=null && alignments.getValue2()!=null && !alignments.getValue1().isEmpty() && !alignments.getValue2().isEmpty()){
			Collection<Pair<SAMRecord>> pairs=this.getValidPairs(pair1, pair2);
			if(pairs.size()==1){this.uniquePairs++;}
			else{this.multimappingPairs++;}
			
			
			pairs=this.getBestPairs(pairs);
			for(Pair<SAMRecord> pair: pairs){
				PairedEndAlignment paired=new PairedEndAlignment(pair);
				SAMRecord r=new SAMRecord(pair.getValue1().getName(), new RefSeqGene(paired.getInsertAlignment()), "*");
				if(pairs.size()>1){r.setMappingQuality(0); r.setWeight(1.0/pairs.size());}
				else{r.setMappingQuality(255);}
				//Alignments r=paired.getInsertAlignment();
				writerPairs.write(r+"\n");
			}*/
			/*for(SAMRecord r1: alignments.getValue1()){
				for(SAMRecord r2: alignments.getValue2()){
					if(r1!=null && r2!=null){
						PairedEndAlignment paired=new PairedEndAlignment(r1.getGene(), r2.getGene(), r1.getName());
						Pair<SAMRecord> pair=new Pair<SAMRecord>(r1, r2);
						if(paired!=null && paired.getInsertAlignment()!=null && isValid(pair)){
							SAMRecord r=new SAMRecord(r1.getName(), new RefSeqGene(paired.getInsertAlignment()), "*");
							if(alignments.getValue1().size()>1 || alignments.getValue2().size()>1){r.setMappingQuality(0);}
							else{r.setMappingQuality(255);}
							//Alignments r=paired.getInsertAlignment();
							writerPairs.write(r+"\n");
						}
					}
				}
			}
		}*/
	}
	
	private Collection<Pair<SAMRecord>> getBestPairMappings(Collection<SAMRecord> pair1, Collection<SAMRecord> pair2){
		
		//Get all pairs
		Collection<Pair<SAMRecord>> validPairs=getValidPairs(pair1, pair2);
		
		//Get best pairs
		Collection<Pair<SAMRecord>> bestPairs=getBestPairs(validPairs);
		
		return bestPairs;
	}

	/*private Pair<Collection<SAMRecord>> getPairMappings(Collection<SAMRecord> pair1, Collection<SAMRecord> pair2) {
		Collection<SAMRecord> rtrn1=new TreeSet<SAMRecord>();
		Collection<SAMRecord> rtrn2=new TreeSet<SAMRecord>();
		
		Collection<Pair<SAMRecord>> validPairs=getValidPairs(pair1, pair2);
		

		//Step 4: If both pairs match and only one unique pairing keep it and mark as unique
		if(validPairs.size()==1){
			rtrn1.add(validPairs.iterator().next().getValue1());
			rtrn2.add(validPairs.iterator().next().getValue2());
		}
		else if(pair1==null || pair2==null){
			//Step 2: If only one part of a pair but it is unique, keep it and mark as unique
			if(((pair1!=null &&pair1.size()==1) && (pair2==null || pair2.isEmpty()))){
				rtrn1.add(pair1.iterator().next());
			}
			
			else if((((pair2!=null && pair2.size()==1) && (pair1==null || pair1.isEmpty())))){
				rtrn2.add(pair2.iterator().next());
			}
		}
		//Step 1: If only one possible pair keep it and mark as unique
		else if(pair1.size()==1 && pair2.size()==1){
			rtrn1.add(pair1.iterator().next());
			rtrn2.add(pair2.iterator().next());
		}
		
		else if(validPairs.isEmpty()){
			//resolve each pair individually
			Collection<Pair<SAMRecord>> pairs=new ArrayList();
			for(SAMRecord r: pair1){
				Pair<SAMRecord> pair=new Pair<SAMRecord>(r, null);
				pairs.add(pair);
			}
			for(SAMRecord r: pair2){
				Pair<SAMRecord> pair=new Pair<SAMRecord>(null, r);
				pairs.add(pair);
			}
			Collection<Pair<SAMRecord>> bestPairs=getBestPairs(pairs);
			if(bestPairs.size()==1){
				rtrn1.add(bestPairs.iterator().next().getValue1());
				rtrn2.add(bestPairs.iterator().next().getValue2());
			}
			else{
				for(Pair<SAMRecord> best: bestPairs){
					if(best.getValue1()!=null){rtrn1.add(best.getValue1());}
					if(best.getValue2()!=null){rtrn2.add(best.getValue2());}
				}
			}
		}
		else{
			Collection<Pair<SAMRecord>> bestPairs=getBestPairs(validPairs);
			//Step 5: If both pairs match and non-unique try to resolve best pair by lowest mismatches
			if(bestPairs.size()==1){
				rtrn1.add(bestPairs.iterator().next().getValue1());
				rtrn2.add(bestPairs.iterator().next().getValue2());
			}
			else{
				for(Pair<SAMRecord> best: bestPairs){
					rtrn1.add(best.getValue1());
					rtrn2.add(best.getValue2());
				}
			}
		}
		
		return new Pair<Collection<SAMRecord>>(rtrn1, rtrn2);
	}*/

	private Collection<Pair<SAMRecord>> getBestPairs(Collection<Pair<SAMRecord>> validPairs) {
		
		Map<Pair<SAMRecord>, Double> scores=new HashMap<Pair<SAMRecord>, Double>();
		
		for(Pair<SAMRecord> pair: validPairs){
			double score=scorePairs(pair);
			scores.put(pair, score);
		}
		
		Collection<Pair<SAMRecord>> bestScores=getMaxScores(scores);
		return bestScores;
	}
	
	private Collection<Pair<SAMRecord>> getMaxScores( Map<Pair<SAMRecord>, Double> scores) {
		Collection<Pair<SAMRecord>> rtrn=new HashSet<Pair<SAMRecord>>();
		
		double max=-Double.MAX_VALUE;
		
		for(Pair<SAMRecord> record: scores.keySet()){max=Math.max(scores.get(record), max);}
		
		for(Pair<SAMRecord> record: scores.keySet()){
			double score=scores.get(record);
			if(score==max){rtrn.add(record);}
		}
		
		return rtrn;
	}

	private double scorePairs(Pair<SAMRecord> pair) {
		double score1=0;
		if(pair.getValue1()!=null){score1=pair.getValue1().getScore();}
		double score2=0;
		if(pair.hasValue2()){score2=pair.getValue2().getScore();}
		return score1+score2; 
	}

	private Collection<Pair<SAMRecord>> getValidPairs(Collection<SAMRecord> pair1, Collection<SAMRecord> pair2) {
		//Collection<Pair<SAMRecord>> rtrn=new TreeSet<Pair<SAMRecord>>();
		Map<SAMRecord, Pair<SAMRecord>> map=new TreeMap<SAMRecord, Pair<SAMRecord>>();
		for(SAMRecord r1: pair1){
			for(SAMRecord r2: pair2){
				Pair<SAMRecord> pair=new Pair<SAMRecord>(r1, r2);
				if(isValid(pair)){
					PairedEndAlignment paired=new PairedEndAlignment(pair);
					SAMRecord r=new SAMRecord(pair.getValue1().getName(), new RefSeqGene(paired.getInsertAlignment()), "*");
					map.put(r, pair);
					//rtrn.add(pair);
				}
			}
		}
		return map.values();
	}

	private boolean isValid(Pair<SAMRecord> pair) {
		SAMRecord r1=pair.getValue1();
		SAMRecord r2=pair.getValue2();
		
		//System.err.println(r1.getGene().getAlignment().toUCSC()+"\t"+r2.getGene().getAlignment().toUCSC()+" "+dist);
		//is on same chr
		if(!r1.getGene().getChr().equalsIgnoreCase(r2.getGene().getChr())){return false;}
		
		int dist=new PairedEndAlignment(r1.getGene(), r2.getGene(), r1.getName()).getInsertAlignment().length();
		
		
		//is valid if less than xbp from each other
		if(dist>maxPairDistance){return false;}
		
		//aligned on opposite strands
		if(r1.getGene().getOrientation().equalsIgnoreCase(r2.getGene().getOrientation())){return false;}
		
		//is valid if the pairs dont overlap
		if(dist<=minPairDistance){return false;}
		
		return true;
	}

	private Pair<File> mergeBowtieAndSalvaged(Pair<String> sortedBowtie, Pair<File> salvagedReads, Runtime run) throws IOException, InterruptedException {
		//Step 1: Merge files keeping track of the program names that generated it
		int bowtieCount1=0;
		int bowtieCount2=0;
		
		int salvagedCount1=0;
		int salvagedCount2=0;
		
		File save1=new File(this.saveDir+"/Nonunique.Pair1.sam");
		FileWriter writer1=new FileWriter(save1);
		bowtieCount1=mergeWrite(new File(sortedBowtie.getValue1()), writer1, "bowtie");
		salvagedCount1=mergeWrite(salvagedReads.getValue1(), writer1, "salvaged");
				
		File save2=new File(this.saveDir+"/Nonunique.Pair2.sam");
		if(this.hasPairs && sortedBowtie.hasValue2()){
			FileWriter writer2=new FileWriter(save2);
			bowtieCount2=mergeWrite(new File(sortedBowtie.getValue2()), writer2, "bowtie");
			salvagedCount2=mergeWrite(salvagedReads.getValue2(), writer2, "salvaged");
		}
		
		this.salvagedCounts=new Pair<Integer>(salvagedCount1, salvagedCount2);
		this.bowtieCounts=new Pair<Integer>(bowtieCount1, bowtieCount2);
		
		String jobID="U"+System.nanoTime();
		this.sortByReadName(save1.getAbsolutePath(), jobID, save1+".temp", run);
		if(this.hasPairs){this.sortByReadName(save2.getAbsolutePath(), jobID, save2+".temp", run);}
		waitForAlignmentJobs(jobID, run);
		
		this.renameFile(save1+".temp", save1.getAbsolutePath(), run);
		if(this.hasPairs){this.renameFile(save2+".temp", save2.getAbsolutePath(), run);}
		
		return new Pair<File>(save1, save2);
	}
	
	private int mergeWrite(File sam, FileWriter writer, String programName) throws IOException {
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(sam)));
			
			int counter=0;
			int count=0;
			String nextLine;
			while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
				try{
					SAMRecord record=new SAMRecord(nextLine);
					record.setName(record.getName().split("/")[0]);
					record.setProgramName(programName);
					//System.err.println(programName);
					//if(record.getProgramName()==null){record.setProgramName(programName);}
					writer.write(record+"\n");
					counter++;
				}catch(Exception ex){}
			}
			return counter;
	}

	private Pair<File> salvageReads(Pair<String> sortedBowtie) throws IOException, InterruptedException {
		String save1=saveDir+"/salvaged.Pair1.sam";
		int count1=0;
		int count2=0;
		if(runTophat){count1=new SalvageUnalignedReads(this.blatAligned.getValue1(), this.tophatAligned.getValue1(), new File(sortedBowtie.getValue1()), save1).getSalvagedCount();}
		else{count1=new SalvageUnalignedReads(this.blatAligned.getValue1(), null, new File(sortedBowtie.getValue1()), save1).getSalvagedCount();}
		String save2=null;
		if(this.hasPairs && sortedBowtie.hasValue2()){
			save2=saveDir+"/salvaged.Pair2.sam";
			if(runTophat){count2=new SalvageUnalignedReads(this.blatAligned.getValue2(), this.tophatAligned.getValue2(), new File(sortedBowtie.getValue2()), save2).getSalvagedCount();}
			else{count2=new SalvageUnalignedReads(this.blatAligned.getValue2(), null, new File(sortedBowtie.getValue2()), save2).getSalvagedCount();}
			return new Pair<File>(new File(save1), new File(save2)); 
		}
		
		return new Pair<File>(new File(save1), null); 
	}

	private void sortAndIndexOld(Pair<File> pairedMapping, String save, Runtime run) throws IOException, InterruptedException {
		if(pairedMapping.hasValue2()){
			merge(pairedMapping, save);
			
			String jobID="U"+System.nanoTime();
			
			sortAlignmentFile(save, save+".temp", jobID, run);
			sortAlignmentFile(pairedMapping.getValue1().getAbsolutePath(), pairedMapping.getValue1().getAbsolutePath()+".temp", jobID, run);
			sortAlignmentFile(pairedMapping.getValue2().getAbsolutePath(), pairedMapping.getValue2().getAbsolutePath()+".temp", jobID, run);
			
			this.waitForAlignmentJobs(jobID, run);
			
			renameFile(pairedMapping.getValue1().getAbsolutePath()+".temp", pairedMapping.getValue1().getAbsolutePath(), run);
			renameFile(pairedMapping.getValue2().getAbsolutePath()+".temp", pairedMapping.getValue2().getAbsolutePath(), run);
			renameFile(save+".temp", save, run);
			
			jobID="U"+System.nanoTime();
			
			indexAlignmentFile(save, jobID, run);
			indexAlignmentFile(pairedMapping.getValue1().getAbsolutePath(), jobID, run);
			indexAlignmentFile(pairedMapping.getValue2().getAbsolutePath(), jobID, run);
			
			this.waitForAlignmentJobs(jobID, run);
		}
		else{
			save=pairedMapping.getValue1().getAbsolutePath();
			String jobID="U"+System.nanoTime();
			sortAlignmentFile(save, save+".temp", jobID, run);
			this.waitForAlignmentJobs(jobID, run);
			renameFile(save+".temp", save, run);
			jobID="U"+System.nanoTime();
			indexAlignmentFile(save, jobID, run);
		}
	}
	
	private void sortAndIndex(String alignFile, Runtime run) throws IOException, InterruptedException {
		String jobID="U"+System.nanoTime();
		sortAlignmentFile(alignFile, alignFile+".temp", jobID, run);
		this.waitForAlignmentJobs(jobID, run);
		renameFile(alignFile+".temp", alignFile, run);
		jobID="U"+System.nanoTime();
		indexAlignmentFile(alignFile, jobID, run);
	
	}

	private void merge(Pair<File> pairedMapping, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		writeToFile(pairedMapping.getValue1(), writer);
		writeToFile(pairedMapping.getValue2(), writer);
		writer.close();
		
	}

	private int renameFile(String oldName, String newName, Runtime run) throws IOException, InterruptedException {
		String command="mv "+oldName+" "+newName;
		Process p=run.exec(command);
		int completed=p.waitFor();
		int exitVal=p.exitValue();
		PipelineUtils.writeError(p.getErrorStream());
		return exitVal;
	}

	private int sortAlignmentFile(String alignmentFile, String save, String jobID, Runtime run) throws IOException, InterruptedException{
		String command="bsub -R rusage[mem=5] -q "+queue+" -J "+jobID+ " -o " + saveDir + "/lsf.bsub /xchip/igv/tools/igvtools sort "+alignmentFile+" "+save;
		return PipelineUtils.bsubProcess(run, command);
	}
	
	private int indexAlignmentFile(String alignmentFile, String jobID, Runtime run) throws IOException, InterruptedException{
		String command="bsub -R rusage[mem=5] -q "+queue+" -J "+jobID+ " -o "+ saveDir + "/lsf.bsub /xchip/igv/tools/igvtools index "+alignmentFile;
		return  PipelineUtils.bsubProcess(run, command);
	}
		
	/*private void mergeAlignments(String save) throws IOException {
		//merge the bowtie results and salvaged stuff
		FileWriter writer=new FileWriter(save);
		writeToFile(this.bowtieNonUniqueMappers, writer);
		writeToFile(this.salvagedReads, writer);
		writer.close();
	}*/

	private void writeToFile(File file, FileWriter writer) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		
    	String nextLine;
    	while ((nextLine = reader.readLine()) != null) {
    		writer.write(nextLine+"\n");
    	}
    	reader.close();
	}

	/*private File salvageReads() throws IOException, InterruptedException {
		String save=saveDir+"/salvaged.All.sam";
		String saveFiltered=saveDir+"/salvaged.Filtered.sam";
		new SalvageUnalignedReads(this.blatAligned, this.tophatAligned, this.bowtieNonUniqueMappers, save);
		//Filter the salvaged reads
		filterSalvagedReads(save, saveFiltered);
		return new File(saveFiltered);
	}*/

	private void filterSalvagedReads(String in, String out) throws IOException {
		new FilterSalvagedReads(new File(in), out);
	}

	/*private void getUniqueMappers(String sortedBowtie) throws IOException {
		//Get unique mappers and resolvable non-unique mappers
		this.bowtieUniqueMappers=new File(saveDir+"/BowtieUniqueMappers.sam");
		this.bowtieNonUniqueMappers=new File(saveDir+"/GenomeAndTranscriptome.NU.sam");
		new CombineBowtieFromVariousReferences(sortedBowtie, bowtieUniqueMappers.getAbsolutePath(), bowtieNonUniqueMappers.getAbsolutePath());
	}*/

	private void mergeValidSAMRecords(String save, String[] samFiles) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(int i=0; i<samFiles.length; i++){
			if(samFiles[i]!=null){
				BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(samFiles[i])));
				
				int counter=0;
				int count=0;
				String nextLine;
				while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
					if(!nextLine.startsWith("@")){
						//try{
						String[] tokens=nextLine.split("\t");
						String match=tokens[2];
						if(!match.equalsIgnoreCase("*")){
							writer.write(nextLine+"\n");
							count++;
						}
						//}catch(Exception ex){}
					}
					counter++;
					//if(counter%100000 ==0){System.err.println(counter+" Good: "+count);}
				}
			}
		}
		
		writer.close();
	}

	private String alignReadsWithBowtie(File fastqFile, Runtime run) throws IOException, InterruptedException {
		Pair<File> pair=new Pair<File>(fastqFile, null);
		return alignReadsWithBowtie(fastqFile, run);
	}
		
	
	private Pair<String> alignReadsWithBowtie(Pair<File> fastqFiles, Runtime run) throws IOException, InterruptedException {
		String jobID="U"+System.nanoTime();
		String outDir=saveDir+"/Bowtie/";
		new File(outDir).mkdir();
		this.alignedGenome=bowtie(fastqFiles, bowtieGenomeIndex, outDir, jobID, run);
		Pair<String> novel=bowtieToJunctions(fastqFiles, novelJunctionBowtieIndex, outDir, jobID, run);
		Pair<String> known=null;
		if(this.useKnownJunctions){
			known=bowtieToJunctions(fastqFiles, knownJunctionBowtieIndex, outDir, jobID, run);
		}
		waitForAlignmentJobs(jobID, run);
		writeToLog("BOWTIE Alignments done " + DateFormat.getDateTimeInstance().format(new Date()));
		//convert known and novel junctions to genome coordinate space
		this.alignedNovel=transcriptomeToGenome(novel, this.novelJunctionSequenceMap);
		writeToLog("Mapped novel junction alignments to genome ");
		//Reads are in order of read names
		//Gets all best strata mappers
		filterUnique(alignedNovel, run);
		writeToLog("FILTER unique done ");
		if(this.useKnownJunctions){
			this.alignedKnown=transcriptomeToGenome(known, this.knownJunctionSequenceMap);
			filterUnique(alignedKnown, run);
		}
		
		
		//merge all valid lines into a single file
		//String bowtieAll=outDir+"/bowtieAll.sam";
		//mergeValidSAMRecords(bowtieAll);
		Pair<String> bowtieAll=mergeValidSAMRecords(outDir);
		writeToLog("MERGE ALIGNMENTS DONE " + DateFormat.getDateTimeInstance().format(new Date()));
		//sort by read name
		Pair<String> sortedBowtie=sortByReadName(bowtieAll, run);
		
		return sortedBowtie;
	}

	private Pair<String> sortByReadName(Pair<String> bowtieAll, Runtime run) throws IOException, InterruptedException {
		//send this to the farm
		String jobID="U"+System.nanoTime();
		String save1=bowtieAll.getValue1()+".nameSort.sam";
		sortByReadName(bowtieAll.getValue1(), jobID, save1, run);
		String save2=null;
		if(this.hasPairs && bowtieAll.hasValue2()){
			save2=bowtieAll.getValue2()+".nameSort.sam";
			sortByReadName(bowtieAll.getValue2(), jobID, save2, run);
		}
		waitForAlignmentJobs(jobID, run);
		
		return new Pair<String>(save1, save2);
	}

	private Pair<String> mergeValidSAMRecords(String outDir) throws IOException {
		String save1=outDir+"/bowtieAllPair1.sam";
		String save2=outDir+"/bowtieAllPair2.sam";
		
		String[] samFiles1={this.alignedGenome.getValue1(), null, this.alignedNovel.getValue1()};
		if(this.useKnownJunctions){
			samFiles1[1]=this.alignedKnown.getValue1();
		}
		mergeValidSAMRecords(save1, samFiles1);
		
		if(this.hasPairs){
			String[] samFiles2={this.alignedGenome.getValue2(), null, this.alignedNovel.getValue2()};
			if(this.useKnownJunctions){
				samFiles2[1]=this.alignedKnown.getValue2();
			}	
			mergeValidSAMRecords(save2, samFiles2);
		}
		
		return new Pair<String>(save1, save2);
	}

	private void filterUnique(Pair<String> alignedNovel2, Runtime run) throws IOException, InterruptedException {
		filterUnique(alignedNovel2.getValue1(), run);
		if(alignedNovel2.hasValue2()){filterUnique(alignedNovel2.getValue2(), run);}
	}

	private Pair<String> transcriptomeToGenome(Pair<String> novel, File novelJunctionSequenceMap2) throws IOException {
		String save1=novel.getValue1()+".genome.sam";
		String save2=null;
		
		Map<String, RefSeqGene> geneMap=BEDFileParser.loadDataByName(novelJunctionSequenceMap2);
		
		transcriptomeToGenome(novel.getValue1(), geneMap, save1, chromosomeSizes);
		if(novel.hasValue2()){
			save2=novel.getValue2()+".genome.sam";
			transcriptomeToGenome(novel.getValue2(), geneMap, save2,chromosomeSizes);
		}
		
		return new Pair<String>(save1, save2);
	}

	private Pair<String> bowtieToJunctions(Pair<File> fastqFiles, String novelJunctionBowtieIndex2, String outDir, String jobID, Runtime run) throws IOException, InterruptedException {
		String sam1=bowtieToJunctions(fastqFiles.getValue1(), novelJunctionBowtieIndex2, outDir, jobID, run);
		String sam2=null;
		if(fastqFiles.hasValue2()){
			sam2=bowtieToJunctions(fastqFiles.getValue2(), novelJunctionBowtieIndex2, outDir, jobID, run);
		}
		return new Pair<String>(sam1, sam2);
	}

	private Pair<String> bowtie(Pair<File> fastqFiles, String bowtieGenomeIndex2, String outDir, String jobID, Runtime run) throws IOException, InterruptedException {
		String sam1=bowtie(fastqFiles.getValue1(), bowtieGenomeIndex2, outDir, jobID, run);
		String sam2=null;
		if(fastqFiles.hasValue2()){
			sam2=bowtie(fastqFiles.getValue2(), bowtieGenomeIndex2, outDir, jobID, run);
		}
		return new Pair<String>(sam1, sam2);
	}

	private void filterUnique(String fileName, Runtime run) throws IOException, InterruptedException {
		CombineBowtieFromVariousReferences.writeUniqueResolvableMappers(fileName, fileName+".tmp");
		//rename the file
		int exitCode=renameFile(fileName+".tmp", fileName, run);
		if(exitCode!=0){throw new IllegalArgumentException("Renamingxg of the file failed");}
	}

	private String sortByReadName(String bowtieAll, String jobID, String save, Runtime run) throws IOException, InterruptedException {
		String command="bsub -q "+queue+" -R rusage[mem=8] -o "+saveDir+"/lsf.bsub -J "+jobID+" sort "+bowtieAll+" -o "+save+" -T "+ TMP_DIR;
		PipelineUtils.bsubProcess(run, command);
		
		return save;
	}



	//TODO This still may not work right for truncated alignments
	static String transcriptomeToGenome(String samFile, Map<String, RefSeqGene> geneMap, String save, Map<String, Integer> chromosomeSizes) throws IOException {
		//FileWriter writer=new FileWriter(save);
		final SAMFileReader inputSam = new SAMFileReader(new File(samFile));
		inputSam.setValidationStringency(ValidationStringency.SILENT);
		SAMFileHeader inputHeader = inputSam.getFileHeader();
		//final SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(inputSam.getFileHeader(), false, new File(save));
		SAMFileHeader header = inputSam.getFileHeader().clone();
		SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
		for(String chr : chromosomeSizes.keySet()) {
			dictionary.addSequence(new SAMSequenceRecord(chr, chromosomeSizes.get(chr)));
		}
		header.setSequenceDictionary(dictionary);
		header.addProgramRecord(new SAMProgramRecord("SCRIPTURE"));
		header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
		//MORAN BUG FIX ; change iterator
		//SAMRecordIterator it = inputSam.iterator();
		net.sf.samtools.util.CloseableIterator<net.sf.samtools.SAMRecord> it = inputSam.iterator();
		
		//BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(samFile)));
		SAMFileWriterFactory alnWriterFactory = new SAMFileWriterFactory();
		SAMFileWriter alignmentWriter = alnWriterFactory.makeSAMOrBAMWriter(header, false, new File(save));
		//BAMFileWriter bamWriter = new BAMFileWriter(new File(save));
		//bamWriter.setHeader(header);
		int i=0;
		int count=0;
		//String nextLine;
		//while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
		int unMappableToGenome = 0;
		int cantFindAnnotation = 0;
		while(it.hasNext()) {
			net.sf.samtools.SAMRecord sam = it.next();
			//RefSeqGene relativeSAM=SAMUtils.SAMFormatFullBED(nextLine);
			String chr = sam.getReferenceName();
			RefSeqGene gene=geneMap.get(chr.toUpperCase());
			if(gene==null){
				gene=geneMap.get(chr);
			} 
			if(gene == null) {
				//System.err.println("Looing for "+ relativeSAM.getChr().split("\\.")[0] + " instead of "+ relativeSAM.getChr() + " referenceName: " + sam.getReferenceName());
				gene=geneMap.get(chr.split("\\.")[0].toUpperCase());
			}
			if(gene == null) {
				gene=geneMap.get(chr.split("\\.")[0]);
			}
			
			if(gene==null && !"*".equals(sam.getReferenceName())){
				cantFindAnnotation++;
				System.err.println("Gene for " + sam.getReferenceName() + " was null");
				continue;
			}
			
			net.sf.samtools.SAMRecord  mappedSam = mapSamRecord(header, sam, gene);;

			if(mappedSam != null) {
				alignmentWriter.addAlignment(sam);
				//writer.write(samRecord.toSAM(nextLine, flip)+"\n");
				count++;
			} else {
				unMappableToGenome++;
			}
			
		}
		i++;
		if(i%100000 ==0){System.err.println(i+" Good: "+count);}
		alignmentWriter.close();
		it.close();
		inputSam.close();
		System.err.println("Unmapped: " + unMappableToGenome + " mapped: " + count + " can't find annotation " + cantFindAnnotation);
		return save;
	}

	/**
	 * Maps a single SAM record from the transcriptome to the genome TODO: TEST THIS REFACTORING.
	 * @param header
	 * @param sam
	 * @param gene
	 * @return
	 */
	public static net.sf.samtools.SAMRecord mapSamRecord(SAMFileHeader header,	net.sf.samtools.SAMRecord sam, RefSeqGene gene) {
		net.sf.samtools.SAMRecord  mappedSam = null;
		RefSeqGene samToBed = SAMUtils.SAMFormatFullBED(sam);
		if(gene!=null && samToBed.getNumExons()==1 && header.getSequence(gene.getChr())!= null){
			//System.err.print("Reference " + sam.getReferenceName() ); 
			//System.err.println(" length " + reference.getSequenceLength() + ", difference gene, ref "+(gene.getSize()-reference.getSequenceLength()));
			boolean flip=false;
			int relativeStart=samToBed.getStart();
			int relativeEnd=samToBed.getEnd();
			if(gene.getOrientation().equalsIgnoreCase("-")){
				relativeStart=(gene.getSize())-(samToBed.getEnd());
				relativeEnd=(gene.getSize())-(samToBed.getStart()); //+1
				flip=true;
			}
			RefSeqGene samRecord=gene.trim(relativeStart, relativeEnd); //put back end-1
			if(samRecord != null) {
				samRecord.setName(samToBed.getName());
				samRecord.setSequence(samToBed.getSequence());
				int seqIdx = header.getSequenceIndex(samRecord.getChr());
				
				logger.debug("gene chr: " + samRecord.getChr() + " header idx: " + seqIdx);
				try {
					mappedSam = (net.sf.samtools.SAMRecord) sam.clone();
					mappedSam.setHeader(header);
					samRecord.updatePicardSAM(mappedSam, flip, seqIdx);
				} catch (CloneNotSupportedException e) {
					// TODO Handle this error!
					e.printStackTrace(System.err);
				}
			} else {
				SAMSequenceRecord reference = header.getSequence(sam.getReferenceName());
				logger.error("Could not map record gene was: " + gene.toBED() + " original record: " + sam.toString());
				//logger.info("Could not map back to genome read relative start-end: "+ relativeStart + "-" +relativeEnd + " ---- " + sam.getReadName() + " mapped to " + reference.getSequenceName() + " aligned reference length was " + reference.getSequenceLength() + " while annotation length was " + gene.getSize());
			}
		}
		else if(samToBed.getNumExons()>1){System.err.println("ERROR Too many exons: "+sam.toString());}
		else { logger.error("Could not map sam record to genome: gene " + gene + " header has chr? " + header.getSequence(gene.getChr()) + " record exons: " + samToBed.getNumExons());}
		return mappedSam;
	}

	

	private String bowtie(File fastqFile, String bowtieIndex, String outDir, String jobID, Runtime run) throws IOException, InterruptedException {
		String samFile=outDir+"/"+fastqFile.getName()+"."+new File(bowtieIndex).getName()+".sam";
		String unaligned=outDir+"/"+fastqFile.getName()+"."+new File(bowtieIndex).getName()+".unaligned.fq";
		String multimappers=outDir+"/"+fastqFile.getName()+"."+new File(bowtieIndex).getName()+".multimappers.fq";
		
		String command="bsub -q "+queue+" -R rusage[mem=10] -J "+jobID+" -o "+outDir+"/"+new File(bowtieIndex).getName()+".bsub " +bowtie+" "+bowtieIndex+" "+fastqFile.getAbsolutePath()+" "+samFile +" --sam -a --best --strata -m "+this.maxNumberOfPlacements+" --un "+unaligned+" --max "+multimappers;
		PipelineUtils.bsubProcess(run , command);

		return samFile;
	}

	private String bowtieToJunctions(File fastqFile, String bowtieIndex, String outDir, String jobID, Runtime run) throws IOException, InterruptedException {
		String samFile=outDir+"/"+fastqFile.getName()+"."+new File(bowtieIndex).getName()+".sam";
		String unaligned=outDir+"/"+fastqFile.getName()+"."+new File(bowtieIndex).getName()+".unaligned.fq";
		String multimappers=outDir+"/"+fastqFile.getName()+"."+new File(bowtieIndex).getName()+".multimappers.fq";
		
		String command="bsub -q "+queue+" -R rusage[mem=10] -J "+jobID+" -o "+outDir+"/"+new File(bowtieIndex).getName()+".bsub "+bowtie +" "+bowtieIndex+" "+fastqFile.getAbsolutePath()+" "+samFile +" --sam -a --best --strata --un "+unaligned+" --max "+multimappers;
		
		PipelineUtils.bsubProcess(run, command);
		return samFile;
	}
	
	private void findNovelJunctions(Pair<File[]> fastqChunks, Runtime run) throws IOException, InterruptedException {

		
		String jobID="U"+System.nanoTime();
		
		//Step 1: BLAT to genome and get all junctions
		if(this.blatAligned==null){
			System.err.println("Running BLAT jobs");
			this.blatOutputDirectories=blatToGenome(fastqChunks, jobID, run);
		}
		
		//Step 2: Tophat to genome and get all junctions
		if(this.tophatAligned==null && this.runTophat){
			System.err.println("Running Tophat");
			this.tophatOutputDirectories=tophatToGenome(fastqChunks, jobID, run);
		}
		
		//Step 3: GSNAP to genome and get all junctions (optional for later)
		//String gsnapIDs=gsnapToGenome(fastqFiles, saveDir);
		
		//Step 4: Wait for the jobs to finish
		if(this.blatAligned==null || this.tophatAligned==null){
			waitForAlignmentJobs(jobID, run);
		}
		
		//Step 5: Merge the BLAT results and filter output to novel junctions
		if(this.blatAligned==null){
			System.err.println("Merging BLAT results");
			this.blatAligned=mergeBLATResults();
		}
		
		//Step 6: Merge the Tophat results and get junctions
		if(this.tophatAligned==null && this.runTophat){
			this.tophatAligned=mergeTophatResults();
		}
		
		//Step 8: Merge the junction files and make a junction sequence file
		System.err.println("Merging junctions");
		mergeJunctions(this.blatAligned, this.tophatAligned);
		
		System.err.println("Building bowtie index for novel junctions");
		//Step 9: Build a bowtie index of the novel junctions
		buildBowtieIndex(run);
	}
	
	private void buildBowtieIndex(Runtime run) throws IOException, InterruptedException {
		String command="bowtie-build "+this.novelJunctionSequenceFile.getAbsolutePath()+" "+saveDir+"/novelJunctions";
		this.novelJunctionBowtieIndex=saveDir+"/novelJunctions";
		Process p=run.exec(command);
		p.waitFor();
		PipelineUtils.writeError(p.getErrorStream());
	}
	

	//TODO Merge all these results into a single unit based on the graph and decomposition into transcripts
	//TODO use the full segments and the junction pairs
	private void mergeJunctions(Pair<File> blatResults, Pair<File> tophatResults) throws IOException {
		String bedFile=saveDir+"/novelJunctions.bed";
		String seqFile=saveDir+"/novelJunctions.fa";
		
		Map<String, Collection<RefSeqGene>> junctionsByChr=new TreeMap<String, Collection<RefSeqGene>>();
		
		//Take all tophat junctions because we know Cole is good at filtering crap
		if(this.runTophat){
			junctionsByChr=getAllJunctions(tophatResults.getValue1(), junctionsByChr);
			if(tophatResults.hasValue2()){
				junctionsByChr=getAllJunctions(tophatResults.getValue2(), junctionsByChr);
			}
		}
		
		//Take all BLAT junctions that pass a series of constraints
		junctionsByChr=getFilteredJunctions(blatResults.getValue1(), junctionsByChr);
		if(blatResults.hasValue2()){
			junctionsByChr=getFilteredJunctions(blatResults.getValue2(), junctionsByChr);
			//System.err.println("")
		}
		for(String chr : junctionsByChr.keySet()) {
			System.err.println("chr " + chr + " has " + junctionsByChr.get(chr).size());
		}
		//write all junctions
		writeJunctions(junctionsByChr, bedFile, seqFile);
		
		this.novelJunctionSequenceFile=new File(seqFile);
		this.novelJunctionSequenceMap=new File(bedFile);
	}

	private void writeJunctions(Map<String, Collection<RefSeqGene>> junctionsByChr, String bedFile, String seqFile) throws IOException {
		FileWriter bedWriter=new FileWriter(bedFile);
		FileWriter seqWriter=new FileWriter(seqFile);
		
		int total=0;
		int filtered=0;
		System.err.println("Writing junctions");
		for(String chr: junctionsByChr.keySet()){
			
			System.err.println("\tJunctions for "+ chr);
			Collection<RefSeqGene> junctions=junctionsByChr.get(chr);
			String chrSequenceFile=genomeDir+"/"+chr.replaceAll("chr", "")+"/"+chr+".fa";
			FastaSequenceIO fsio = new FastaSequenceIO(chrSequenceFile);
			List<Sequence> seqs = fsio.loadAll();
			System.err.println("\t\tLoaded chromosome sequence "+ chrSequenceFile + ", making graph" );
			Sequence chrSequence = seqs.get(0);
			
			ChromosomeWithBubbles2 graph=makeGraph(junctions, chr, chrSequence);
			System.err.println("\t\tMade graph, extracting paths");
			Collection<Path> paths=graph.getAllPaths();
			System.err.println("\t\tPaths extracted, writing junction file");
			int pathNum=0;
			for(Path path: paths){
				RefSeqGene junction=path.toGene();
				junction.setName(chr+"_"+pathNum);
				//System.err.println(junction);
				//System.err.println(junction);
				//Alignments intron=junction.getIntronSet().iterator().next();
				//String orientation="+";
				//if(this.filterBySplice){orientation=GeneTools.orientationFromSpliceSites(intron, chrSequence);}
				//if(!orientation.equalsIgnoreCase("*")){
					//filtered++;
				//TODO Extend the read by read length at ends	
				bedWriter.write(junction.toBED()+"\n");
				seqWriter.write(">"+junction.getName()+"\n"+junction.getSequence(chrSequence).getSequenceBases()+"\n");
				//}
				total++;
				pathNum++;
			}
			//System.err.println("Total number of reads "+total+" filtered number "+filtered);
		}
			
		bedWriter.close();
		seqWriter.close();
	}
	
	private ChromosomeWithBubbles2 makeGraph(Collection<RefSeqGene> junctions, String chr, Sequence chrSequence) {
		Collection<Alignments> exons=new TreeSet<Alignments>();
		Collection<Alignments> introns=new TreeSet<Alignments>();
		
		for(RefSeqGene gene: junctions){
			String orientation=GeneTools.orientationForGene(gene, chrSequence);
			
			gene.setOrientation(orientation);
			if (!"*".equals(gene.getOrientation()) || !filterBySplice){
				exons.addAll(gene.getExonSet());
				introns.addAll(gene.getIntronSet());
			}
		}
		System.err.println("\t\t\tmakeGraph, #exons: " + exons.size() + " #introns: " + introns.size());
		exons=CollapseByIntersection.collapseByIntersection(exons, false);
		System.err.println("\t\t\tmakeGraph CollapseByIntersection is done.");
		/*try {
			BEDFileParser.writeBED("exons.bed", exons);
			BEDFileParser.writeBED("introns.bed", introns);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}*/
		exons=CollapseByIntersection.DecollapseByIntronLocation(exons, introns);
		System.err.println("\t\t\tmakeGraph DecollapseByIntronLocation is done.");
		ChromosomeWithBubbles2 graph=new ChromosomeWithBubbles2(chr, exons, introns, null, null, 0, 0, 0,0);
		System.err.println("\t\t\tmakeGraph built chromosome with bubbles");
		return graph;
	}

	private void write(String save, Map<String, Collection<RefSeqGene>> junctionsByChr, String genomeDir) throws IOException {
		FileWriter bedWriter=new FileWriter(save+".bed");
		FileWriter seqWriter=new FileWriter(save+".fa");
		
		int total=0;
		int filtered=0;
		
		for(String chr: junctionsByChr.keySet()){
			//System.err.println(chr);
			Collection<RefSeqGene> junctions=junctionsByChr.get(chr);
			String chrSequenceFile=genomeDir+"/"+chr.replaceAll("chr", "")+"/"+chr+".fa";
			FastaSequenceIO fsio = new FastaSequenceIO(chrSequenceFile);
			List<Sequence> seqs = fsio.loadAll();
			Sequence chrSequence = seqs.get(0);
			
			for(RefSeqGene junction: junctions){
				Alignments intron=junction.getIntronSet().iterator().next();
				String orientation=GeneTools.orientationFromSpliceSites(intron, chrSequence);
				
				if(!orientation.equalsIgnoreCase("*") && intron.getSize()>this.minIntronSize && intron.getSize()<this.maxAlignmentGap){
					filtered++;
					bedWriter.write(junction.toBED()+"\n");
					seqWriter.write(">"+junction.getName()+"\n"+junction.getSequence(chrSequence).getSequenceBases()+"\n");
				}
				total++;
			}
			//System.err.println("Total number of reads "+total+" filtered number "+filtered);
		}
			
		bedWriter.close();
		seqWriter.close();
	}

	private Map<String, Collection<RefSeqGene>> getFilteredJunctions(File blatResults, Map<String, Collection<RefSeqGene>> junctionsByChr) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(blatResults)));
		String nextLine;
		int counter=0;
		int total=0;
		//Collection<RefSeqGene> temp=new TreeSet<RefSeqGene>();
		
		while ((nextLine = reader.readLine()) != null) {
			try{
				RefSeqGene gene=SAMUtils.SAMFormatFullBED(nextLine);
				if(getPercentMapped(gene, readLength)>this.minPercentMapped && getMinCoverage(gene)>this.minCoverageOnEachSide){
					Collection<RefSeqGene> junc=getFilteredJunctions(gene, readLength);
					Collection<RefSeqGene> junctions=new TreeSet<RefSeqGene>();
					if(junctionsByChr.containsKey(gene.getChr())){
						junctions=junctionsByChr.get(gene.getChr());
					}
					junctions.addAll(junc);
					junctionsByChr.put(gene.getChr(), junctions);
					//temp.addAll(junc);
					counter++;
					//writer.write(nextLine+"\n");
				}
				total++;
				if(total %100000 ==0){System.err.println("filtering BLAT file "+counter+" "+total);}
			}catch(Exception ex){
				System.err.println("Skipping: "+nextLine);
			}
		}
		
		int junctionNum = 0;
		for(Collection<RefSeqGene> junctions: junctionsByChr.values()) {
			junctionNum += junctions.size();
		}
		this.novelJunctionCounts.put(blatResults.getAbsolutePath(), junctionNum);
		reader.close();
		return junctionsByChr;
	}
	
	private int getMinCoverage(RefSeqGene gene) {
		int min=Integer.MAX_VALUE;
		Collection<Alignments> exons=gene.getSortedAndUniqueExons();
		for(Alignments exon: exons){min=Math.min(min, exon.length());}
		if(min==Integer.MAX_VALUE){min=0;}
		return min;
	}

	private double getPercentMapped(RefSeqGene gene, int readLength) {
		return (double)gene.getTranscriptLength()/readLength;
	}
	
	private Collection<RefSeqGene> getJunctions(RefSeqGene gene, int readLength2) {
		if(gene.getNumExons()<2){return new TreeSet();}
		
		Collection<RefSeqGene> junctions=new TreeSet();
		junctions.add(gene);
		
		Collection<Alignments> introns=gene.getIntronSet();
		Collection<Alignments> exons=gene.getSortedAndUniqueExons();
		
		Collection<Alignments> newExons=new TreeSet();
		int i=0;
		for(Alignments exon: exons){
			if(i==0){
				int start=Math.min(exon.getStart(), exon.getEnd()-readLength2);
				Alignments newExon=new Alignments(exon.getChr(), start, exon.getEnd());
				newExons.add(newExon);
			}//first
			else if(i== exons.size()-1){
				int end=Math.max(exon.getEnd(), exon.getStart()+readLength2);
				Alignments newExon=new Alignments(exon.getChr(), exon.getStart(), end);
				newExons.add(newExon);
			}//last		
			else{newExons.add(exon);}
			i++;
		}
		
		RefSeqGene junction=new RefSeqGene(newExons);
		junctions.add(junction);
		
		/*if(introns.size()>2){junctions.add(gene);}
		
		for(Alignments intron: introns){
			Alignments exon1=new Alignments(intron.getChr(), intron.getStart()-readLength2, intron.getStart());
			Alignments exon2=new Alignments(intron.getChr(), intron.getEnd(), intron.getEnd()+readLength2);
			Collection<Alignments> exons=new TreeSet<Alignments>();
			exons.add(exon1);
			exons.add(exon2);
			RefSeqGene jun=new RefSeqGene(exons);
			junctions.add(jun);
		}*/
		
		return junctions;
	}

	private Collection<RefSeqGene> getFilteredJunctions(RefSeqGene gene, int readLength2) {
		if(gene.getNumExons()<2 ){
			return new TreeSet<RefSeqGene>();
		}
		
		Collection<RefSeqGene> junctions=new TreeSet<RefSeqGene>();
		Collection<Alignments> introns=gene.getIntronSet();
		
		/*boolean good=true;
		for(Alignments intron: introns){
			if(intron.getSize()>this.minIntronSize && intron.getSize()<this.maxIntronSize){}
			else{good=false;}
		}
		
		if(good){junctions.add(gene);}*/
		
		Collection<Alignments> exons=gene.getSortedAndUniqueExons();
		
		Collection<Alignments> newExons=new TreeSet<Alignments>();
		int i=0;
		for(Alignments exon: exons){
			if(i==0){
				int start=Math.min(exon.getStart(), exon.getEnd()-readLength2);
				Alignments newExon=new Alignments(exon.getChr(), start, exon.getEnd());
				newExons.add(newExon);
			}//first
			else if(i== exons.size()-1){
				int end=Math.max(exon.getEnd(), exon.getStart()+readLength2);
				Alignments newExon=new Alignments(exon.getChr(), exon.getStart(), end);
				newExons.add(newExon);
			}//last		
			else{newExons.add(exon);}
			i++;
		}
		
		boolean good=true; 
		Iterator<Alignments> intronIt = introns.iterator();
		while(intronIt.hasNext() && good){
			Alignments intron = intronIt.next();
			good = intron.getSize()>this.minIntronSize && intron.getSize()<this.maxAlignmentGap;
		}
		
		
		if(good){
			RefSeqGene junction=new RefSeqGene(newExons);
			junctions.add(junction);
		}
		
		/*if(introns.size()>2){junctions.add(gene);}
		
		for(Alignments intron: introns){
			if(intron.getSize()>this.minIntronSize && intron.getSize()<this.maxIntronSize){
				Alignments exon1=new Alignments(intron.getChr(), intron.getStart()-readLength2, intron.getStart());
				Alignments exon2=new Alignments(intron.getChr(), intron.getEnd(), intron.getEnd()+readLength2);
				Collection<Alignments> exons=new TreeSet<Alignments>();
				exons.add(exon1);
				exons.add(exon2);
				RefSeqGene jun=new RefSeqGene(exons);
				junctions.add(jun);
			}
		}*/
		
		return junctions;
	}
	
	private Map<String, Collection<RefSeqGene>> getAllJunctions(File samFile, Map<String, Collection<RefSeqGene>> junctionsByChr) throws IOException {
		//Map<String, Collection<RefSeqGene>> junctionsByChr=new TreeMap<String, Collection<RefSeqGene>>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(samFile)));
		String nextLine;
		int counter=0;
		int total=0;
		while ((nextLine = reader.readLine()) != null) {
			try{
				RefSeqGene gene=SAMUtils.SAMFormatFullBED(nextLine);
				if(gene.getNumExons()>1){
					Collection<RefSeqGene> junc=getJunctions(gene, readLength);
					Collection<RefSeqGene> junctions=new TreeSet();
					if(junctionsByChr.containsKey(gene.getChr())){
						junctions=junctionsByChr.get(gene.getChr());
					}
					junctions.addAll(junc);
					junctionsByChr.put(gene.getChr(), junctions);
					total++;
					//if(total %100000 ==0){System.err.println(counter+" "+total);}	
				}
			}catch(Exception ex){System.err.println("Skipping: "+nextLine);}
		}
			
		this.novelJunctionCounts.put(samFile.getAbsolutePath(), total);
		reader.close();
		return junctionsByChr;
	}

	private Pair<File> mergeBLATResults() throws IOException {
		String save1=saveDir+"/blatMergedPair1.sam";
		String save2=saveDir+"/blatMergedPair2.sam";
		mergeBLATResults(this.blatOutputDirectories.getValue1(), save1);
		if(this.blatOutputDirectories.hasValue2()){
			mergeBLATResults(this.blatOutputDirectories.getValue2(), save2);
			return new Pair<File>(new File(save1), new File(save2));
		}
		return new Pair<File>(new File(save1), null);
	}
	
	private String mergeBLATResults(Collection<File> blatFiles, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		for(File blatDir: blatFiles){
			File[] files=blatDir.listFiles();
			for(int i=0; i<files.length; i++){
				if(files[i].getName().endsWith(".sam")){
					System.err.println(files[i]);
					BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(files[i])));
					String nextLine;
			        while ((nextLine = reader.readLine()) != null) {
			        	writer.write(nextLine+"\n");
			        }
			        reader.close();
				}
			}
		}
		writer.close();
		return save;
	}

	private Pair<File> mergeTophatResults() throws IOException{
		String save1=saveDir+"/tophatMergedPair1.sam";
		String save2=saveDir+"/tophatMergedPair2.sam";
		mergeTophatResults(this.tophatOutputDirectories.getValue1(), save1);
		if(this.tophatOutputDirectories.hasValue2()){
			mergeTophatResults(this.tophatOutputDirectories.getValue2(), save2);
		}
		return new Pair<File>(new File(save1), new File(save2));
	}
	
	private String mergeTophatResults(Collection<File> topHatFiles, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		for(File blatDir: topHatFiles){
			File[] files=blatDir.listFiles();
			for(int i=0; i<files.length; i++){
				if(files[i].getName().equalsIgnoreCase("accepted_hits.sam")){
					System.err.println(files[i]);
					BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(files[i])));
					String nextLine;
			        while ((nextLine = reader.readLine()) != null) {
			        	writer.write(nextLine+"\n");
			        }
			        reader.close();
				}
			}
		}
		writer.close();
		return save;
	}
	

	private Collection<File> tophatToGenome(File[] fastqFiles, String jobID, Runtime run) throws IOException, InterruptedException {
		Collection<File> rtrn=new ArrayList<File>();
		for(int i=0; i<fastqFiles.length; i++){
			File save=new File(saveDir+"/"+fastqFiles[i].getName()+"/Tophat/");
			save.mkdir();
			String command="bsub -q "+queue+" -J "+jobID+" -o "+saveDir+"/"+fastqFiles[i].getName()+"/tophatJunk.bsub "+tophat+" -o "+save.getAbsolutePath() +" "+bowtieGenomeIndex+" "+fastqFiles[i].getAbsolutePath();
			PipelineUtils.bsubProcess(run, command);
			rtrn.add(save);
		}
		return rtrn;
	}
	
	private Pair<Collection<File>> tophatToGenome(Pair<File[]> fastqFiles, String jobID, Runtime run) throws IOException, InterruptedException {
		Collection<File> files1=tophatToGenome(fastqFiles.getValue1(), jobID, run);
		Collection<File> files2=null;
		if(fastqFiles.hasValue2()){files2=tophatToGenome(fastqFiles.getValue2(), jobID, run);}
		return new Pair<Collection<File>>(files1, files2);
	}

	private ArrayList<File> blatToGenome(File[] fastqFiles, String jobID, Runtime run) throws IOException, InterruptedException {
		return runAllBLAT(fastqFiles, genomeDir, saveDir, jobID, run);
	}
	
	private Pair<Collection<File>> blatToGenome(Pair<File[]> fastqFiles, String jobID, Runtime run) throws IOException, InterruptedException {
		Pair<Collection<File>> rtrn=new Pair<Collection<File>>();
		ArrayList<File> blatFiles1=blatToGenome(fastqFiles.getValue1(), jobID, run);
		ArrayList<File> blatFiles2=null;
		if(fastqFiles.hasValue2()){blatFiles2=blatToGenome(fastqFiles.getValue2(), jobID, run);}

		rtrn.setValue1(blatFiles1);
		rtrn.setValue2(blatFiles2);
		
		return rtrn;
	}
	
	private ArrayList<File> runAllBLAT(File[] fastqFiles, String genomeDir, String saveDir, String jobID, Runtime run) throws IOException, InterruptedException {
		ArrayList<File> rtrn=new ArrayList<File>();
		String classPath = System.getProperty("java.class.path");
		System.err.println("Class path: " + classPath);
		for(int i=0; i<fastqFiles.length; i++){
			new File(saveDir+"/"+fastqFiles[i].getName()+"/").mkdir();
			String blatSaveDir=saveDir+"/"+fastqFiles[i].getName()+"/BLAT/";
			new File(blatSaveDir).mkdir();
			File blatPipelineOut = new File(blatSaveDir+"/"+fastqFiles[i].getName()+".sam");
			System.err.println("Looking for " + blatPipelineOut);
			if(blatPipelineOut.exists() && blatPipelineOut.length() > 0 && !rechunk) {
				writeToLog("found non empty BLAT generated SAM file for CHUNK"+ fastqFiles[i].getName() +"  NOT RUNNING BLAT for ");
			} else { 
				String command="bsub -q "+queue+" -J "+jobID+" -o "+saveDir+"/"+fastqFiles[i].getName()+"/"+"blatJunk.bsub"+" java -cp "+ classPath + " " + blatScriptClass;
				command+=" "+fastqFiles[i].getAbsolutePath()+" "+genomeDir+" "+blatSaveDir+" "+this.blatMinScore+" "+this.blatMinPercentIdentity;
				PipelineUtils.bsubProcess(run, command);
			}
			rtrn.add(new File(blatSaveDir));
		}
		return rtrn;
	}

	private void waitForAlignmentJobs(String jobID, Runtime run) throws IOException, InterruptedException {	
		int running=5;
		
		while(running>1){
			Thread.sleep(waitTime);
			Process blatProc = run.exec("bjobs -J "+jobID);
			BufferedReader out = new BufferedReader(new InputStreamReader(blatProc.getInputStream()));
			running=parseReply(out);
			System.err.println("checking status "+running+" jobs still running");
		}
		System.err.println("done");
		PipelineUtils.checkForLSFFailures(jobID, run);
	}
	

	
	private static int parseReply(BufferedReader lsfOut) throws IOException {
		String line = null;
		String [] lineInfo = null;
		int i=0;
		while((line = lsfOut.readLine()) != null) {i++;}
		return i;
	}

	private File[] splitIntoChunks(final File fastqFile) throws IOException {
		File dir=new File(saveDir+"/chunks/");
		if(!dir.exists()) {
			dir.mkdir();
		}
		
		String save=saveDir+"/chunks/"+fastqFile.getName();
		
		String [] existingChunks = dir.list(new FilenameFilter() {

			public boolean accept(File arg0, String arg1) {
				return arg1.contains(fastqFile.getName());
			}
			
		});
		File [] files = null;
		if(existingChunks.length > 0 && !rechunk) {
			writeToLog("found non empty chunk directory and the rechunk flag is unset USING EXISTING CHUNKS");
			files = new File [existingChunks.length];
			for(int i = 0; i < existingChunks.length; i++) {
				files[i] = new File(dir.getAbsolutePath() + "/" + existingChunks[i]);
			}
		} else {
			FastqParser fastq=new FastqParser(fastqFile);
			files=fastq.writeChunks(save, chunkSize);
			readCounts.put(fastqFile, fastq.getNumberOfSequences());
		}
		return files;
	}
	
	/*public static void main(String[] args)throws IOException, InterruptedException{
		Pair<File> salvagedReads=new Pair<File>(new File(args[0]), new File(args[1]));
		String save=args[2];
		
		RNASeqAlignmentPipeline aligner=new RNASeqAlignmentPipeline();
		aligner.runAlignment(salvagedReads, save);
	}*/
	
	public static void main(String[] args)throws IOException, InterruptedException{		
		ArgumentMap argmap = CLUtil.getParameters(args, usage, "full");
		if("full".equals(argmap.getTask())) {
			int maxGapSize = argmap.containsKey("maxAlignmentGapSize") ? argmap.getInteger("maxAlignmentGapSize") : MAX_GAP;
			boolean filterBySplice = argmap.containsKey("filterBySplice");
			File fastqFile1=new File(argmap.getMandatory("fastq1"));
			File fastqFile2=argmap.isPresent("fastq2") ? new File(argmap.get("fastq2")) : null;
			Pair<File> fastqFiles=new Pair<File>(fastqFile1, fastqFile2);
			
			String genomeDir=argmap.getMandatory("genomeDir");
			String saveDir=argmap.getMandatory("saveDir");
			
			String bowtieGenomeIndex=argmap.getMandatory("bowtieGenomeIndex");
			boolean rechunk = argmap.containsKey("rechunk");
			
			String bowtieJunctionIndex=argmap.isPresent("bowtieJunctionIndex")? argmap.get("bowtieJunctionIndex"): null;
			File junctionMap=argmap.isPresent("junctionMap")? new File(argmap.getMandatory("junctionMap")): null;
			
			boolean useKnown=bowtieJunctionIndex!=null && junctionMap!=null;
			
			String save=argmap.getMandatory("save");
			
			Pair<File> tophatFiles=null;
			Pair<File> blatFiles=null;
			Pair<String> bowtieFiles=null;
			
			boolean runTophat=argmap.containsKey("runTophat");
			if(argmap.isPresent("tophat1")){
				File tophat1=new File(argmap.get("tophat1"));
				File tophat2=argmap.isPresent("tophat2") ? new File(argmap.get("tophat2")) : null;
				tophatFiles=new Pair<File>(tophat1, tophat2);
				runTophat=true;
			}
			
			if(argmap.isPresent("blat1")){
				File blat1=new File(argmap.get("blat1"));
				File blat2=argmap.isPresent("blat2") ? new File(argmap.get("blat2")) : null;
				blatFiles=new Pair<File>(blat1, blat2);
			}
			
			if(argmap.isPresent("bowtie1")){
				String bowtie1=argmap.get("bowtie1");
				String bowtie2=argmap.isPresent("bowtie2") ? (argmap.get("bowtie2")) : null;
				bowtieFiles=new Pair<String>(bowtie1, bowtie2);
			}
			
			RNASeqAlignmentPipeline aligner=new RNASeqAlignmentPipeline(fastqFiles, genomeDir, saveDir, bowtieGenomeIndex, bowtieJunctionIndex, junctionMap, save, tophatFiles, blatFiles, runTophat, useKnown, bowtieFiles, maxGapSize, filterBySplice, rechunk);
			aligner.runAlignment(filterBySplice);
		} else if ("makepairs".equalsIgnoreCase(argmap.getTask())) {
			String pair1 = argmap.getMandatory("pair1");
			String pair2 = argmap.getMandatory("pair2");
			Pair<File> alnPair = new Pair<File>(new File(pair1), new File(pair2));
			String outdir = argmap.getOutputDir();
			RNASeqAlignmentPipeline aligner=new RNASeqAlignmentPipeline(outdir);
			aligner.makePairs(alnPair);
		}
	}
	
	static String usage="\nParameters \n -fastq1 <Sequence reads corresponding to the first mate pair> \n -fastq2 <Sequence reads corresponding to the second mate pair> \n -save <Output file name>\n -genomeDir <Genome sequence directory>"+ 
	"\n -saveDir <save directory>  \n -bowtieGenomeIndex <bowtie genome index> \n -bowtieJunctionIndex <bowtie known junction index> \n -junctionMap <BED file mapping junction locations to genome> \n -tophat1 <SAM file of tophat alignment for pair1> \n -tophat2 <SAM file of tophat alignment for pair2> \n -blat1 <SAM file of BLAT alignments for pair1> \n -blat2 <SAM file for BLAT alignments for pair2> \n -runTophat \n -bowtie1 <SAM file of bowtie aligned for pair1> \n -bowtie2 <SAM file of bowtie aligned for pair2>" +
	"\n -maxAlignmentGapSize <Maximum gap size to consider. The higher this number the more complex the internal graph building becomes>\n -filterBySplice <Include this flag to only consider alignments flanked by splice sites>\n";
}
