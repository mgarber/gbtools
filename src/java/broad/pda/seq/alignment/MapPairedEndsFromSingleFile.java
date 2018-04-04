package broad.pda.seq.alignment;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;


public class MapPairedEndsFromSingleFile {

	private int maxPairDistance=1000000;
	private int minPairDistance=-1;


	HashMap<String,SAMRecord> cacheMap;
	static Logger logger = Logger.getLogger(MapPairedEndsFromSingleFile.class.getName());

	public  MapPairedEndsFromSingleFile(String alignmentFile, String sizes, String out, String applicationName, int minimumMappingScore, boolean usePair2Orientation,  boolean forChIP) throws IOException {

		cacheMap=new HashMap<String,SAMRecord>();
		//@Option(shortName=StandardOptionDefinitions.SORT_ORDER_SHORT_NAME, doc="Sort order of output file")         public SAMFileHeader.SortOrder SORT_ORDER;


		final SAMFileReader inputSam = new SAMFileReader (new File(alignmentFile));
		inputSam.setValidationStringency(ValidationStringency.STRICT);
		if(inputSam.getFileHeader().getSequenceDictionary().isEmpty() && sizes != null) {
			SAMBAMUtils.createSequenceDictionary(inputSam, sizes);
		}
		CloseableIterator <SAMRecord> readIter =inputSam.iterator();
		SAMFileHeader samHeader=inputSam.getFileHeader(); 
		//GenericAlignmentDataModel alignments=new GenericAlignmentDataModel(alignmentFile, sizes, false);
		//Iterator<Alignment> readIter=alignments.getReadIterator();
		RecordWriter writer = new RecordWriter(out, samHeader.clone(), applicationName);
		writer.setMinimumMappingScore(minimumMappingScore);
		Integer currIndex =-1;
		String currRefName = "";
		int recordsAdded = 0; 
		int unmapped=0;
		int mateUnmapped=0;
		int diffChrs=0;
		int notProcUnknown=0;
		int totalBAMRecords = 0;
		int notPassingQual = 0;
		int timesCacheToBig = 1;
		int cachedUnused = 0;
		while (readIter.hasNext()) {
			totalBAMRecords++;
			//System.err.println("RECORD " + currIndex);
			SAMRecord record = readIter.next();

			//clean the cache map when a record from a new chr appears (file is assumed to be sorted by position)
			if (currIndex==-1) {
				currIndex=record.getReferenceIndex(); 
				currRefName = record.getReferenceName();
			}
			if (!record.getReadUnmappedFlag() &&  record.getReferenceIndex() != currIndex){
				cachedUnused += cacheMap.size();
				cleanCashMap();
				logger.debug("Cleaned cache! Reference was: " + currRefName + " new ref is: " + record.getReferenceName()+"  \n");
				currIndex=record.getReferenceIndex();
				currRefName = record.getReferenceName();
			}

			//System.err.println("Mate flag= "+ record.getMateUnmappedFlag() + " read flag: "+ record.getReadUnmappedFlag() + " mate chr: "  +record.getMateReferenceIndex() + " readChr: "+ record.getReferenceIndex());
			if (record.getMappingQuality() > minimumMappingScore && record.getReadPairedFlag() && !record.getMateUnmappedFlag() && !record.getReadUnmappedFlag() && record.getMateReferenceIndex() == record.getReferenceIndex() ){
				String mateName = record.getAlignmentStart()+"_"+record.getReadName();
				//System.err.println(mateName);
				if (cacheMap.containsKey(mateName)){
					//both mates are found- report to output and remove mate from cache
					SAMRecord mateAlignment=cacheMap.get(mateName);
					//if (mateAlignment != null) //SHould not need to check if we already checked that the key is in. 
					//{
						//System.err.println(	"Read neg strand :" + record.getReadNegativeStrandFlag() +" mate negative strand "+ mateAlignment.getReadNegativeStrandFlag());

					Map<SAMRecord,Integer> insertInfo=createInsertRecordFromFragmentPair(samHeader,record,mateAlignment, usePair2Orientation, forChIP);
					int insertSize =  insertInfo.values().iterator().next();
					if (insertSize <maxPairDistance && insertSize> minPairDistance ){
						SAMRecord pairedRecord = insertInfo.keySet().iterator().next();
						//System.err.println("Original record "+ record.format());
						//System.err.println("Adding insert "+ pairedRecord.format());	
						writer.addAlignment(pairedRecord);
						// System.err.println("Added insert ");
						recordsAdded++;
						if(recordsAdded % 500000 == 0) {
							logger.debug(recordsAdded + " paired end processed, chache size " + cacheMap.size());
						}
					}
				//	}
					cacheMap.remove(mateName);
				}
				else if (record.getAlignmentStart() <= record.getMateAlignmentStart()) //first mate in the fragment; store in cache until second mate is found
				{
					String recordName= record.getMateAlignmentStart()+"_"+record.getReadName();
					cacheMap.put(recordName, record);
				}

				//If hash is too big clean all cashed records that their mate is has a start value smaller than the current record
				if (cacheMap.size()>= 500000){
					cleanCashMap(record);
				}
			} else {
				//logger.trace("Read " + record.getReadName() + " is ignored, is unmapped? " + record.getReadUnmappedFlag() + " is mate unmapped? " + record.getMateUnmappedFlag() + " is mate on a different chr? " + (record.getMateReferenceIndex() != record.getReferenceIndex()));
				if(record.getMappingQuality() < minimumMappingScore) {notPassingQual++;}
				else if(record.getReadUnmappedFlag()) {unmapped++;}
				else if(record.getMateUnmappedFlag()) {mateUnmapped++;}
				else if(record.getMateReferenceIndex() != record.getReferenceIndex()) {diffChrs++;}
				else {notProcUnknown++;}
			}

			//debug
			if (cacheMap.size()>= 500000){
				FileWriter debugOut=new FileWriter("DebugMap_"+timesCacheToBig+".sam");
				logger.debug("cache size is "+ cacheMap.size() +"\n");
				int i=0;
				for (SAMRecord s: cacheMap.values()) 
				{
					debugOut.write(s.format() + "\n");
					i++;
					if (i>20000)
						break;
				}
				debugOut.close();
				cachedUnused += cacheMap.size();
				cleanCashMap();
				timesCacheToBig++;
				logger.error("The cache is too big, that means that pair mates are not being found at the expected distance. Something is wrong with the data. Ignoring reads, cleaning up cache and continuing" );
				//break;

			}

		}

		readIter.close();
		logger.info("End of mapping pairs , \n\ttotal BAM recdords " + totalBAMRecords + " \n\trecords processed " + recordsAdded + " \n\tcache size is "+ cacheMap.size() +" \n\tlow quality alignments " + notPassingQual + " \n\tunmapped reads: " + unmapped + " \n\tmateUnmapped: " + mateUnmapped + " \n\tmate in diff char " + diffChrs  + " \n\talignments left in cache after chromosome scan finished: " + cachedUnused + " \n\tunprocessed unknown " + notProcUnknown + "\n");

		writer.close();
		inputSam.close();




	}





	private void cleanCashMap(SAMRecord record) {

		int i=0;
		int currPos=record.getAlignmentStart();
		List<String> rmlst=new LinkedList<String>();
		for(String key : cacheMap.keySet()){
			if (cacheMap.get(key).getMateAlignmentStart() < currPos)
				rmlst.add(key);
		}
		for (String key: rmlst){
			cacheMap.remove(key);
			i++;
		}
		logger.debug("Cleared "+ i +" records with mate pair that maps before " + currPos + "  , cache size is "+ cacheMap.size() +"\n");

	}





	private void cleanCashMap() {
		cacheMap.clear();
	}





	private Map<SAMRecord, Integer> createInsertRecordFromFragmentPair(SAMFileHeader samHeader,SAMRecord r1, SAMRecord r2, boolean usePair2Orientation, boolean forChIP){

		SAMRecord pair1 = r1.getFirstOfPairFlag() ? r1 : r2;
		SAMRecord pair2 = r1.getFirstOfPairFlag() ? r2 : r1;
		LightweightGenomicAnnotation mappedPair1 = new BasicLightweightAnnotation(r1.getReferenceName(), r1.getAlignmentStart(), r1.getAlignmentEnd());
		LightweightGenomicAnnotation mappedPair2 = new BasicLightweightAnnotation(r2.getReferenceName(), r2.getAlignmentStart(), r2.getAlignmentEnd());
		if(mappedPair1.contains(mappedPair2)) {
			//This is a bad case in which one mapped read engulfs the other. This can't happen. However there are cases that mapped like this and thus are mapping artifacts.
			//When this happens we ignore the reads.
			HashMap<SAMRecord, Integer> recordInfo = new HashMap<SAMRecord, Integer>();
			recordInfo.put(pair1, minPairDistance);
			return recordInfo ;
			
		} else if ( mappedPair2.contains(mappedPair1))  {
			HashMap<SAMRecord, Integer> recordInfo = new HashMap<SAMRecord, Integer>();
			recordInfo.put(pair2, minPairDistance);
			return recordInfo ;
		}
		
		
		//OK all seems good lets handle the feasable cases.
		SAMRecord first = r1.getAlignmentStart() < r2.getAlignmentStart() ? r1 : r2;
		SAMRecord second = r1.getAlignmentStart() < r2.getAlignmentStart() ? r2 : r1;
		int readDist =  second.getAlignmentStart()- first.getAlignmentEnd() ;
		
		hardClipEndSoftClipAtEnd(first);
		hardClipEndSoftClipAtStart(second);
		
		final SAMRecord srec = new SAMRecord(samHeader);
		int start= Math.min(r1.getAlignmentStart(), r2.getAlignmentStart());
		int end=Math.max(r1.getAlignmentEnd(), r2.getAlignmentEnd());
		//int size=Math.abs(end-start);
		srec.setReadName(r1.getReadName());
		srec.setNotPrimaryAlignmentFlag(true);
		srec.setReadNegativeStrandFlag(usePair2Orientation ? pair2.getReadNegativeStrandFlag() : pair1.getReadNegativeStrandFlag());
		srec.setReadPairedFlag(false);
		srec.setAlignmentStart(start); 
		List<CigarElement> firstElements = null;
		List<CigarElement> secondElements = null;
		

		
		Cigar pairedCigar = null;
		if(readDist>0) { //No overlap
			srec.setReadString(first.getReadString() + second.getReadString());
			srec.setBaseQualityString(first.getBaseQualityString()+ second.getBaseQualityString());
			firstElements = first.getCigar().getCigarElements();
			secondElements = second.getCigar().getCigarElements();
			pairedCigar = new Cigar();
			//int elementNum = 0;
			for(CigarElement e : firstElements) {
				//if(e.getOperator().equals(CigarOperator.S) && elementNum > 0) {
				//	e = new CigarElement(e.getLength(), CigarOperator.I);
				//}
				pairedCigar.add(e);
				//elementNum++;
			}
			if(readDist > 1) {
				pairedCigar.add(new CigarElement(readDist-1, forChIP ? CigarOperator.D : CigarOperator.N));
			}
			
			//elementNum = secondElements.size() - 1;;
			for(CigarElement e : secondElements) {
				//if(e.getOperator().equals(CigarOperator.S) && elementNum > 0) {
				//	e = new CigarElement(e.getLength(), CigarOperator.I);
				//}
				pairedCigar.add(e);
				//elementNum--;
			} 
		} else { //Overlap
			pairedCigar = new Cigar();
			StringBuilder readBases = new StringBuilder();
			StringBuilder qualityBases = new StringBuilder();
			
			//Cigar firstReadCigar = getTrimmedCigar(genomicPosition, sam)
			//System.err.println("\tUnique first read bases:" + firstUniquePortion);
			for(int i = 0; i < first.getReadLength(); i++) {
				readBases.append(first.getReadString().charAt(i));
				qualityBases.append(first.getBaseQualities()[i]);
			}
			
			//TODO: go back and examine the routine commented out to pick the best base among overlapping bases and avoid using soft clipped bases 
			List<CigarElement> pairedCigarElements =  new ArrayList<CigarElement>(first.getCigar().getCigarElements());
			//First change soft clipping at end of read to missmatch
			//CigarElement lastReadLastElement = pairedCigarElements.get(pairedCigarElements.size() - 1);
			//if(lastReadLastElement.getOperator().equals(CigarOperator.S)) {
			//	pairedCigarElements.set(pairedCigarElements.size() - 1,  new CigarElement(lastReadLastElement.getLength(), CigarOperator.I));
				
			//}

			
			/*
			System.err.println("Doing non-unique bases: from " + firstUniquePortion + " to read length: " + first.getReadLength());
			for(int  i = firstUniquePortion; i < first.getReadLength(); i++) {
				if(first.getBaseQualities()[i] >  second.getBaseQualities()[i]) {
					readBases.append(first.getReadString().charAt(i));
					qualityBases.append(first.getBaseQualities()[i]);
				} else {
					readBases.append(second.getReadString().charAt(i));
					qualityBases.append(second.getBaseQualities()[i]);
				}
			}
			*/
			
			Cigar secondReadPortion = getEndCigar(first.getAlignmentEnd(), second);
			List<CigarElement> secondReadCigarElements = secondReadPortion.getCigarElements();
			//Need to convert soft clips of first element to missmatches as they will no longer be at the beggining 
			/*if(!secondReadCigarElements.isEmpty()) {
				CigarElement firstSecondReadElment = secondReadCigarElements.get(0);
				if(firstSecondReadElment.getOperator().equals(CigarOperator.S)) {
					CigarElement newfirstSecondReadElment = new CigarElement(firstSecondReadElment.getLength(), CigarOperator.I);
					secondReadCigarElements.set(0, newfirstSecondReadElment);
				}				
			}*/


			Iterator<CigarElement> secondReadCigarElementIt = secondReadCigarElements.iterator();
			if(secondReadCigarElementIt.hasNext()) {
				CigarElement e = secondReadCigarElementIt.next();
				CigarElement currLastE = pairedCigarElements.get(pairedCigarElements.size() - 1);
				//System.err.println("Op1: " + e.getOperator() + " op2: " + currLastE.getOperator() + " are they equal? " + (e.getOperator().equals(currLastE.getOperator())));
				if(e.getOperator().equals(currLastE.getOperator())){
					pairedCigarElements.remove(pairedCigarElements.size() - 1);
					CigarElement newLast = new CigarElement(currLastE.getLength()+ e.getLength(), e.getOperator());
					pairedCigarElements.add(newLast);
					//System.err.println("Removed " + e.getOperator()+e.getLength() + " and added " + newLast.getOperator() + newLast.getLength());
				} else {
					pairedCigarElements.add(e);
					//System.err.println(" added " + e.getOperator() + e.getLength());
				}
			}
			while(secondReadCigarElementIt.hasNext()) {
				pairedCigarElements.add(secondReadCigarElementIt.next());
			}
			
			for(CigarElement ce : pairedCigarElements) {
				pairedCigar.add(ce);
			}
			int numOfBases = getNumberOfReadBases(secondReadPortion);
			//logger.debug("Doing second read unique portion: from " + (second.getReadLength() - numOfBases) + " to read length: " + second.getReadLength());
			for(int i = second.getReadLength() - numOfBases ; i < second.getReadLength(); i++) {
				readBases.append(second.getReadString().charAt(i));
				qualityBases.append(second.getBaseQualities()[i]);
			}
			
		}
		srec.setCigar(pairedCigar);
		srec.setReferenceName(r1.getReferenceName());
		srec.setMappingQuality(Math.round((r1.getMappingQuality() +r2.getMappingQuality())/2));
		srec.setDuplicateReadFlag(r1.getDuplicateReadFlag() && r2.getDuplicateReadFlag());
		int r1NH = r1.getAttribute("NH") ==null ? 1 : Integer.valueOf(r1.getIntegerAttribute("NH"));
		int r2NH = r2.getAttribute("NH") ==null ? 1 : Integer.valueOf(r2.getIntegerAttribute("NH"));
		int r1NM = r1.getAttribute("NM") ==null ? 0 : Integer.valueOf(r1.getIntegerAttribute("NM"));
		int r2NM = r2.getAttribute("NM") ==null ? 0 : Integer.valueOf(r2.getIntegerAttribute("NM"));
		srec.setAttribute("NH",(int) Math.min(r1NH,  r2NH)); //I am assuming that the aligner take into account 
		srec.setAttribute("NM",r1NM + r2NM);
		int dist=Math.abs(end-start);
		HashMap<SAMRecord, Integer> recordInfo = new HashMap<SAMRecord, Integer>();
		recordInfo.put(srec, dist);
		return recordInfo ;

	}

	private void hardClipEndSoftClipAtEnd(SAMRecord e) {
		List<CigarElement> cigarElements = e.getCigar().getCigarElements();
		
		CigarElement lastElement = cigarElements.get(cigarElements.size() - 1);
		if(lastElement.getOperator().equals(CigarOperator.S)) {
			//logger.debug("aln was: " + e.getCigarString() + " - " + e.getReadString() + " - " + e.getBaseQualityString());
			int toTrim = lastElement.getLength();
			List<CigarElement> hardClippedElements = new ArrayList<CigarElement>(cigarElements.size() - 1);
			for(int i = 0 ; i< cigarElements.size() - 1; i++) {
				hardClippedElements.add(cigarElements.get(i));
			}
			Cigar trimmedCigar = new Cigar(hardClippedElements);
			String sequence = e.getReadString();
			String quals    = e.getBaseQualityString();
			
			sequence = sequence.substring(0, sequence.length() - toTrim);
			quals    = quals.substring(0, quals.length() - toTrim);
			
			e.setBaseQualityString(quals);
			e.setCigar(trimmedCigar);
			e.setReadString(sequence);
			//logger.debug("end trimmed aln new is: " + e.getCigarString() + " - " + e.getReadString() + " - " + e.getBaseQualityString());
			
		}
		
	}
	
	private void hardClipEndSoftClipAtStart(SAMRecord e) {
		List<CigarElement> cigarElements = e.getCigar().getCigarElements();
		
		CigarElement firstElement = cigarElements.get(0);
		if(firstElement.getOperator().equals(CigarOperator.S)) {
			//logger.debug("aln was: " + e.getCigarString() + " - " + e.getReadString() + " - " + e.getBaseQualityString());
			int toTrim = firstElement.getLength();

			List<CigarElement> hardClippedElements = new ArrayList<CigarElement>(cigarElements.size() - 1);
			for(int i = 1 ; i< cigarElements.size(); i++) {
				hardClippedElements.add(cigarElements.get(i));
			}
			Cigar trimmedCigar = new Cigar(hardClippedElements);
			String sequence = e.getReadString();
			String quals    = e.getBaseQualityString();
			
			sequence = sequence.substring( toTrim  );
			quals    = quals.substring( toTrim );
			
			e.setBaseQualityString(quals);
			e.setCigar(trimmedCigar);
			e.setReadString(sequence);
			//logger.debug("start trimmed aln is: " + e.getCigarString() + " - " + e.getReadString() + " - " + e.getBaseQualityString());
		}
		
	}





	public void setMaxPairDistance (int maxPairedDist) {
		this.maxPairDistance = maxPairedDist;
	}



	private int getNumberOfReadBases(Cigar secondReadPortion) {
		int numOfReadBases = 0;
		Iterator<CigarElement> cigarElements = secondReadPortion.getCigarElements().iterator();

		while ( cigarElements.hasNext()) {
			CigarElement e = cigarElements.next();
			CigarOperator eo = e.getOperator();
			int length = e.getLength();
			switch (eo) {
			case D: 
			case H:
			case S:
			case EQ:
			case X:
			case M:
				numOfReadBases += length;
				break;
			case I:
				break;
			case N:
				break;
			case P:
				break;
			}

		}
		return numOfReadBases;
	}


	/**
	 * Trims the cigar element to the ones that ends in the given genomic position
	 * @param genomicPosition
	 * @param sam
	 * @return
	 */
	private Cigar getEndCigar(int fromGenomicPosition,  SAMRecord sam) {
		Cigar cigar = sam.getCigar();
		Cigar trimmedCigar = new Cigar();
		List<CigarElement> cigarElements = cigar.getCigarElements();
		int currGenomicPosition = sam.getAlignmentEnd();
		List<CigarElement> reversedTrimmedElements = new ArrayList<CigarElement>();
		CigarElement e = null;
		int i = cigarElements.size() - 1;
		//System.err.println("fromGenomicPosition: " + fromGenomicPosition);
		while ( i >=0 && currGenomicPosition > fromGenomicPosition) {
			e = cigarElements.get(i);
			CigarOperator eo = e.getOperator();
			int length = e.getLength();
			switch (eo) {
			case EQ:
			case X:
			case M:
				int operatorLength = length;
				currGenomicPosition -= length;
				if(currGenomicPosition <= fromGenomicPosition) {
					operatorLength = length - (fromGenomicPosition - currGenomicPosition) ;
				}
				//System.err.println("Added "+operatorLength+eo);
				reversedTrimmedElements.add(new CigarElement(operatorLength, eo));
				break;
			case N:
			case I:
				currGenomicPosition -= length;
				break;
			case D:
				break;
			case H:
				break;
			case P:
				break;
			case S:
				break;
			}
			//System.err.println("i " + i + " -- genomicPosition: " + currGenomicPosition);
			i--;
		}
		for(int j = reversedTrimmedElements.size() - 1 ; j >=0 ; j--) {
			trimmedCigar.add(reversedTrimmedElements.get(j));
		}
		//System.err.println("returning " + trimmedCigar.toString() + " was " + cigar.toString());
		return trimmedCigar;
	}
	

}


