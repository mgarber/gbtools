package broad.pda.seq.alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.pda.alignment.PairedEndAlignment;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.alignment.sam.SAMRecord;

public class MapPairedEnds {
	private int maxPairDistance=1000000;
	private int minPairDistance=50;
	
	int waitTime=60000; //1 minute
	static Logger logger = Logger.getLogger(MapPairedEnds.class.getName());
	String queue="long";
	
	public MapPairedEnds(File file1, File file2, String save, boolean isSorted) throws IOException, InterruptedException{
		Pair<File> pair=new Pair<File>(file1, file2);
		
		if(!isSorted){
			logger.info("File is not sorted by name attempting to sort now");
			pair=sortByName(pair);
		}
		
		this.getPairs(pair, save);
	}
	
	
	public MapPairedEnds(String fileName, String save, boolean isSorted,boolean sortBam) throws IOException, InterruptedException{
		
		String samFileName=fileName;
		if (sortBam){
			System.err.println("File is a BAM,  attempting to convert to SAM now");
			String tmpName= fileName+".sam";
			convertToSam(fileName,tmpName);
			samFileName=tmpName;
			isSorted=false;
		}
		if(!isSorted){
			System.err.println("File is not sorted by name attempting to sort now");
			File newPairFile=sortByName(new File (samFileName));
			this.getPairs(newPairFile, save);
		}
		else
			this.getPairs(new File (samFileName), save);
		
	}
	
	


	private void convertToSam(String bamName, String samName) throws IOException, InterruptedException {

		Runtime run=java.lang.Runtime.getRuntime();
		String jobID="U"+System.nanoTime();
		
		String command="bsub -q "+queue+" -o junk.bsub -J "+jobID+" samtools  view  "+ bamName+" -o "+samName;
		Process p=run.exec(command);
		p.waitFor();
		writeError(p.getErrorStream());
		waitForAlignmentJobs(jobID, run);
		
	}


	private Pair<File> sortByName(Pair<File> pair) throws IOException, InterruptedException {
		Runtime run=java.lang.Runtime.getRuntime();
		
		String jobID="U"+System.nanoTime();
		String val1=sortByReadName(pair.getValue1().getAbsolutePath(), jobID, pair.getValue1()+".nameSorted.sam", run);
		String val2=sortByReadName(pair.getValue2().getAbsolutePath(), jobID, pair.getValue1()+".nameSorted.sam", run);
		
		Pair<File> rtrn=new Pair<File>(new File(val1), new File(val2));
		
		waitForAlignmentJobs(jobID, run);
		return rtrn;
	}
	
	private File sortByName(File file1) throws IOException, InterruptedException {
		Runtime run=java.lang.Runtime.getRuntime();
		String jobID="U"+System.nanoTime();
		String val1=sortByReadName(file1.getAbsolutePath(), jobID, file1+".nameSorted.sam", run);
		File rtrn=new File(val1);
		waitForAlignmentJobs(jobID, run);
		return rtrn;
	}
	
	
	private void waitForAlignmentJobs(String jobID, Runtime run) throws IOException, InterruptedException {	
		int running=5;
		
		while(running>1){
			Thread.sleep(waitTime);
			Process blatProc = run.exec("bjobs -J "+jobID);
			BufferedReader primer3StdOut = new BufferedReader(new InputStreamReader(blatProc.getInputStream()));
			running=parseReply(primer3StdOut);
			System.err.println("checking status "+running+" jobs still running");
		}
		System.err.println("done");
		checkForFailures(jobID, run);
	}
	
	private static void checkForFailures(String jobID, Runtime run) throws IOException {
		Process primer3Proc = run.exec("bjobs -a -J "+jobID);
		BufferedReader primer3StdOut = new BufferedReader(new InputStreamReader(primer3Proc.getInputStream()));
		String line = null;
		int i=0;
		while((line = primer3StdOut.readLine()) != null) {
			//split the line
			String[] tokens=line.split(" ");
			String completionStatus=tokens[2];
			if(completionStatus.equalsIgnoreCase("EXIT")){System.err.println("WARN: Blat Job "+tokens[0]+" FAILED"); throw new IllegalArgumentException("Job failed");}
			//System.err.println(i);
			//System.err.println(line);
			i++;
		}
	}
	
	private static int parseReply(BufferedReader primer3StdOut) throws IOException {
		String line = null;
		String [] lineInfo = null;
		int i=0;
		while((line = primer3StdOut.readLine()) != null) {i++;}
		return i;
	}
	
	private String sortByReadName(String bowtieAll, String jobID, String save, Runtime run) throws IOException, InterruptedException {
		String command="bsub -q "+queue+" -o junk.bsub -J "+jobID+" sort "+bowtieAll+" -o "+save+" -T /broad/shptmp/";
		Process p=run.exec(command);
		p.waitFor();
		writeError(p.getErrorStream());
		
		return save;
	}
	
	private void writeError(InputStream errorStream) throws IOException {
		BufferedReader reader=	new BufferedReader(new InputStreamReader(errorStream));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			System.err.println(nextLine);
		}
		System.err.println();
	}

	private String getName(Collection<SAMRecord> pair) {
		if(pair==null || pair.isEmpty()){return "";}
		else{return pair.iterator().next().getName();}
	}
	
	private void getPairs(Pair<File> allAligned, String save) throws IOException {
		FileWriter writerPairs=new FileWriter(save);
		
		int properPairs=0;
		
		//Given 2 files that are sorted get the pairing relationship
		Iterator<Collection<SAMRecord>> pair1Iter=new SortedSAMIterator(allAligned.getValue1());
		Iterator<Collection<SAMRecord>> pair2Iter=new SortedSAMIterator(allAligned.getValue2());
		
		boolean done=false;
		
		Collection<SAMRecord> pair1=new TreeSet();
		Collection<SAMRecord> pair2=new TreeSet();
		boolean started=false;
		int counter=0;
		while((pair1Iter.hasNext() || pair2Iter.hasNext()) && !done){
			if(!started){
				pair1=pair1Iter.next();
				pair2=pair2Iter.next();
				started=true;
			}
			String name1=getName(pair1);
			String name2=getName(pair2);
						
			//if(pair1.size()>1){System.err.println("1 "+name1+" "+pair1.size());}
			//if(pair2.size()>1){System.err.println("2 "+name2+" "+pair2.size());}
			if(name1.equalsIgnoreCase(name2)){
				//valid pair
				//add both
				printPairMapping(pair1, pair2, writerPairs);
				//proceed both iters DON'T WE MISS THE LAST PAIR IN THE FILE? (WILL NOT HAVE HASnEXT WHEN WE GO TO THE WHILE LOOP)
				if(pair1Iter.hasNext()){pair1=pair1Iter.next();}
				if(pair2Iter.hasNext()){pair2=pair2Iter.next();}
				properPairs++;
			}
			else{
				//2 possibilities
				if(name1.compareTo(name2)<0){
					//then read1 is before read2
					//make pair with read1 and null
					//printPairMapping(pair1, null, writer1, writer2, writerPairs);
					//advance read1 
					if(!pair1Iter.hasNext()){
						done=true;
					} else { 
						pair1=pair1Iter.next();
					}
				}
				else{
					//then read1 is after read2
					//make pair with null and read2
					//printPairMapping(null, pair2, writer1, writer2,writerPairs);
					//advance read2 
					if(!pair2Iter.hasNext()){
						done=true;
					} else {
						pair2=pair2Iter.next();
					}
				}
			}
			counter++;
			if(counter % 100000 ==0){System.err.println(counter+" proper pairs: "+properPairs);}
		}
		
		System.err.println("Total number of proper pairs "+properPairs);
		writerPairs.close();
	}
	
	
	private void getPairs(File pairFile, String save) throws IOException {
		
		FileWriter writerPairs=new FileWriter(save);
		int properPairs=0;
		//Given a sorted by name file get the pairing relationship
		Iterator<Collection<SAMRecord>> pairIter=new SortedSAMIterator(pairFile);
		Collection<SAMRecord> samSet=new TreeSet();
		int counter=0;
		while(pairIter.hasNext() ){
			//get all Sam records with the same name (this set includes records from both ends 
			//of a fragment)
			samSet=pairIter.next();
			//Sort them into the matched mates
			 Map <SAMRecord,SAMRecord>  samMap=splitSamRecordsByEnds (samSet);
			 
			 for (SAMRecord r:samMap.keySet()){
				 Collection<SAMRecord> pair1=new TreeSet();
				 Collection<SAMRecord> pair2=new TreeSet();
				 pair1.add(r);
				 pair2.add(samMap.get(r));
				 printPairMapping(pair1,pair2, writerPairs);
				 properPairs++;
			 }
			 
			 counter++;
			if(counter % 100000 ==0){System.err.println(counter+" proper pairs: "+properPairs);}
			
		}
		
		System.err.println("Total number of proper pairs "+properPairs);
		writerPairs.close();
	}
	
	
	
	private Map<SAMRecord, SAMRecord> splitSamRecordsByEnds(Collection<SAMRecord> samSet) {

		Map <Integer,SAMRecord> end1= new HashMap <Integer,SAMRecord>();
		Map <Integer,SAMRecord> end2=new HashMap <Integer,SAMRecord>();
		Map<SAMRecord, SAMRecord> rtrn= new HashMap<SAMRecord, SAMRecord>();
		
		for (SAMRecord sam: samSet){
			
			if (sam.getFirstOfPairFlag()) {
				end1.put(sam.getPosition()+1,sam);}
			else if(sam.getSecondOfPairFlag())
				{end2.put(sam.getPosition()+1,sam);}
		}
		for (SAMRecord sam: end1.values()){
			if (end2.containsKey(sam.getMatePosition()))
				rtrn.put(sam, end2.get(sam.getMatePosition()));
		}
		
		return rtrn;
	}


	private void printPairMapping(Collection<SAMRecord> pair1,Collection<SAMRecord> pair2,FileWriter writerPairs) throws IOException {
		//go through pairs
		Pair<Collection<SAMRecord>> alignments=getPairMappings(pair1, pair2);
		
		
		if(alignments.getValue1()!=null && alignments.getValue2()!=null && !alignments.getValue1().isEmpty() && !alignments.getValue2().isEmpty()){
			Collection<Pair<SAMRecord>> pairs=this.getValidPairs(pair1, pair2);
			pairs=this.getBestPairs(pairs);
			for(Pair<SAMRecord> pair: pairs){
				PairedEndAlignment paired=new PairedEndAlignment(pair);
				SAMRecord r=new SAMRecord(pair.getValue1().getName(), new RefSeqGene(paired.getInsertAlignment()), "*");
				if(pairs.size()>1){r.setMappingQuality(0); r.setWeight(1.0/pairs.size());}
				else{r.setMappingQuality(255);}
				//Alignments r=paired.getInsertAlignment();
				writerPairs.write(r+"\n");
			}
			
		}
	}
	
	private Pair<Collection<SAMRecord>> getPairMappings(Collection<SAMRecord> pair1, Collection<SAMRecord> pair2) {
		Collection<SAMRecord> rtrn1=new TreeSet<SAMRecord>();
		Collection<SAMRecord> rtrn2=new TreeSet<SAMRecord>();
		
		Collection<Pair<SAMRecord>> validPairs=new TreeSet();
		if(pair1!=null && pair2!=null){
			validPairs=getValidPairs(pair1, pair2);
		}
		

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
	}

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
	
	
	public static void main(String[] args)throws IOException, InterruptedException{
		if(args.length>3){
			File file1=new File(args[0]);
			File file2=new File(args[1]);
			String save=args[2];
			boolean isSorted=new Boolean(args[3]);
			new MapPairedEnds(file1, file2, save, isSorted);
		}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=file1 \n args[1]=file2 \n args[2]=save \n args[3]=is sorted by name?";
}
