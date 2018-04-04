package broad.pda.seq.alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.pda.seq.alignment.sam.SAMRecord;

public class CombineBowtieFromVariousReferences {

	public CombineBowtieFromVariousReferences(String samFile, String save, String saveNU) throws IOException{
		//write all resolvable mappers
		writeUniqueAndResolvableMappers(samFile, save, saveNU);
	}
	
	private void writeUniqueAndResolvableMappers(String samFile, String save, String saveNU) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(samFile)));
		FileWriter writer=new FileWriter(save);
		FileWriter writerNU=new FileWriter(saveNU);
		
		int i=0;
		
		String nextLine;
		String name="";
		Collection<SAMRecord> records=new TreeSet<SAMRecord>();
        while ((nextLine = reader.readLine()) != null) {
        	SAMRecord record=new SAMRecord(nextLine);
        	if(record.getName().equalsIgnoreCase(name)){records.add(record);}
        	else{
        		SAMRecord best=determineBestMapper(records);
        		if(best!=null){writer.write(best+"\n"); writerNU.write(best+"\n");}
        		else{
        			Collection<SAMRecord> bestRecords=getMinValueRecords(records);
        			for(SAMRecord r: bestRecords){
        				r.setWeight((double)1/bestRecords.size());
        				r.setMappingQuality(0);
        				writerNU.write(r+"\n");
        			}
        			//TODO write all records with a weight flag set to 1/the number of alignments
        			//TODO For now give it a mapping quality of 0 for visualization purposes
        		}
        		records=new TreeSet<SAMRecord>(); records.add(record);
        	}
        	name=record.getName();
        	//mark the transition in read name
        	i++;
        	if(i%100000 ==0){System.err.println(i);}
        }
		reader.close();
		writer.close();
		writerNU.close();
	}

	private static Collection<SAMRecord> getMinValueRecords(Collection<SAMRecord> records) {
		Map<SAMRecord, Integer> scores=new TreeMap<SAMRecord, Integer>();
		
		for(SAMRecord record: records){
			int mm=record.getNumMismatches();
			scores.put(record, mm);
		}
		
		//if there is one with fewest mismatches return this
		Collection<SAMRecord> minScores=getMinimumAlignments(scores);
		return minScores;
	}

	public static void writeUniqueMappers(File samFile, String save) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(samFile)));
		FileWriter writer=new FileWriter(save);
		
		//FileWriter writerNU=new FileWriter(save+".NU.sam");
		
		int good=0;
		int all=0;
		
		String nextLine;
		String name="";
		Collection<SAMRecord> records=new TreeSet<SAMRecord>();
        while ((nextLine = reader.readLine()) != null) {
        	SAMRecord record=new SAMRecord(nextLine);
        	if(record.getName().equalsIgnoreCase(name)){records.add(record);}
        	else{
        		if(records.size()==1){writer.write(records.iterator().next().toString()+"\n"); good++;}
        		/*else{
        			for(SAMRecord rec: records){
        				writerNU.write(rec.toString()+"\n");
        			}
        		}*/
        		all++;
        		//if(all % 10000 ==0){System.err.println(good+" "+all +" "+((double)good/all));}
        		records=new TreeSet<SAMRecord>(); 
        		records.add(record);
        	}
        	name=record.getName();
        	//mark the transition in read name
        }
		reader.close();
		writer.close();
		//writerNU.close();
	}

	public static void writeUniqueResolvableMappers(String samFile, String save) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(samFile)));
		FileWriter writer=new FileWriter(save);
		//FileWriter writerNU=new FileWriter(saveNU);
		
		int i=0;
		
		String nextLine;
		String name="";
		Collection<SAMRecord> records=new TreeSet<SAMRecord>();
        while ((nextLine = reader.readLine()) != null) {
        	SAMRecord record=new SAMRecord(nextLine);
        	if(record.getName().equalsIgnoreCase(name)){records.add(record);}
        	else{
        		SAMRecord best=determineBestMapper(records);
        		if(best!=null){writer.write(best+"\n");}
        		else{
        			Collection<SAMRecord> bestRecords=getMinValueRecords(records);
        			for(SAMRecord r: bestRecords){
        				r.setWeight((double)1/bestRecords.size());
        				r.setMappingQuality(0);
        				writer.write(r+"\n");
        			}
        			//TODO write all records with a weight flag set to 1/the number of alignments
        			//TODO For now give it a mapping quality of 0 for visualization purposes
        		}
        		records=new TreeSet<SAMRecord>(); records.add(record);
        	}
        	name=record.getName();
        	//mark the transition in read name
        	i++;
        	if(i%100000 ==0){System.err.println(i);}
        }
		reader.close();
		writer.close();
	}
	
	private Map<String, SAMRecord> resolveMultimappers(Map<String, Collection<SAMRecord>> samRecords) {
		Map<String, SAMRecord> rtrn=new TreeMap<String, SAMRecord>();
		
		int i=0;
		for(String readName: samRecords.keySet()){
			Collection<SAMRecord> records=samRecords.get(readName);
			SAMRecord bestRecord=determineBestMapper(records);
			if(bestRecord!=null){rtrn.put(readName, bestRecord);}
			i++;
			if(i%10000 ==0){System.err.println(i+" "+readName);}
		}
		
		return rtrn;
	}

	private static SAMRecord determineBestMapper(Collection<SAMRecord> records) {
		if(records.size()==1){return records.iterator().next();}
		
		Collection<SAMRecord> minScores=getMinValueRecords(records);
		//if there is one with fewest mismatches return this
		if(minScores.size()==1){return minScores.iterator().next();}
		
		//if they all have the same mismatches and overlap then return genomic path
		//if there are 2 min scores, 1 junction, 1 genomic
		/*if(minScores.size()==2){
			SAMRecord record1=(SAMRecord)minScores.toArray()[0];
			SAMRecord record2=(SAMRecord)minScores.toArray()[1];
			if((record1.getGene().getNumExons()>1 && record2.getGene().getNumExons()==1) || (record1.getGene().getNumExons()==1 && record2.getGene().getNumExons()>1)){
				if(record1.overlaps(record2)){
					if(record1.getGene().getNumExons()==1){return record1;}
					return record2;
				}
			}
		}*/
		
		
		//TODO Get all and return weighted values
		//if they have sam mismatches and dont overlap then return both with weight divided by the number of placements
		return null;
	}

	private static Collection<SAMRecord> getMinimumAlignments( Map<SAMRecord, Integer> scores) {
		Collection<SAMRecord> rtrn=new TreeSet<SAMRecord>();
		
		int min=Integer.MAX_VALUE;
		
		for(SAMRecord record: scores.keySet()){min=Math.min(scores.get(record), min);}
		
		for(SAMRecord record: scores.keySet()){
			double score=scores.get(record);
			if(score==min){rtrn.add(record);}
		}
		
		return rtrn;
	}

	private void write(String save,	Map<String, SAMRecord> resolved) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String name: resolved.keySet()){
			SAMRecord record=resolved.get(name);
			writer.write(record+"\n");
		}
		
		writer.close();
	}

	//For now will return all
	//If memory becomes a problem then just get name chunks and do all the work on that
	//Then unload and proceed
	private Map<String, Collection<SAMRecord>> getMultimappers(File file) throws IOException{
		Map<String, Collection<SAMRecord>> rtrn=new HashMap<String, Collection<SAMRecord>>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		
		int i=0;
		
		String nextLine;
		String name="";
		Collection<SAMRecord> records=new TreeSet<SAMRecord>();
        while ((nextLine = reader.readLine()) != null) {
        	SAMRecord record=new SAMRecord(nextLine);
        	if(record.getName().equalsIgnoreCase(name)){records.add(record);}
        	else{rtrn.put(name, records); records=new TreeSet<SAMRecord>(); records.add(record);}
        	name=record.getName();
        	//mark the transition in read name
        	i++;
        	if(i%100000 ==0){System.err.println(i);}
        }
		reader.close();
		return rtrn;
	}
	
	/*private String countUnique(Collection<SAMRecord> records) {
		if(records.size()>0){
		String rtrn="";
		SAMRecord first=records.iterator().next();
		for(SAMRecord record: records){
			rtrn+=" "+(first.getName()+" "+record.getName()+" "+first.equals(record));
		}
		return rtrn;
		}
		return "";
	}*/

	/*
	//Requires the input to be sorted by read name
	public CombineBowtieFromVariousReferences(File file, String save) throws IOException{
		String save2=save+".filteredSame.sam";
		//For now asusme we have a SAM file
		//write non-unique aligners
		writeNonUniqueAlignments(file, save);
		writeUniqueAlignments(save, save2);
	}

	private void writeUniqueAlignments(String file, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		
		RefSeqGene previousPSL=null;
		boolean flagged=false;
		
		int alignedCounter=0;
		int uniqueAlignedCounter=0;
		int splicedReads=0;
		int splicedReadsUnique=0;
		int counter=0;
		
		String nextLine;
        while ((nextLine = reader.readLine()) != null) {
        	try{
        		RefSeqGene psl=SAMUtils.SAMFormatFullBED(nextLine);              	
				if(previousPSL!=null && !psl.equals(previousPSL)){
					//write previous PSL
					if(!flagged){writer.write(previousPSL.toSAM()+"\n");}
					alignedCounter++;
					if(psl.getNumExons()>1){splicedReads++;}
					flagged=false;
				}
				else{flagged=true;}
				previousPSL=psl;
				counter++;
				if(counter%100000 ==0){System.err.println("Part 2 "+counter);}
        	}catch(Exception ex){System.err.println("Skipping: "+nextLine); ex.printStackTrace();}
				
          }
		writer.close();
	}

	private void writeNonUniqueAlignments(File file, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		
		RefSeqGene previousPSL=null;
		boolean flagged=false;
		
		int alignedCounter=0;
		int uniqueAlignedCounter=0;
		int splicedReads=0;
		int splicedReadsUnique=0;
		int counter=0;
				
		String nextLine;
        while ((nextLine = reader.readLine()) != null) {
        	try{
        		RefSeqGene psl=SAMUtils.SAMFormatFullBED(nextLine);              	
				//if unique
        		if(previousPSL!=null && (!psl.getName().equalsIgnoreCase(previousPSL.getName()))  ){
					//write previous PSL
					if(!flagged){uniqueAlignedCounter++; if(psl.getNumExons()>1){splicedReadsUnique++;}}
					else{writer.write(previousPSL.toSAM()+"\n");}
					alignedCounter++;
					if(psl.getNumExons()>1){splicedReads++;}
					flagged=false;
				}
        		else if(previousPSL!=null){flagged=true; writer.write(previousPSL.toSAM()+"\n");}
        		else if(previousPSL==null){flagged=false;}
        		else{flagged=true;}
				previousPSL=psl;
				counter++;
				if(counter % 100000 ==0){System.err.println(counter);}
        	}catch(Exception ex){System.err.println("Skipping: "+nextLine); ex.printStackTrace();}
          }
		
		System.err.println("Total aligned "+alignedCounter+" Unique aligned " +uniqueAlignedCounter+" percent "+((double)uniqueAlignedCounter/alignedCounter));
		writer.close();
	}*/

	/*public static void main(String[] args)throws IOException{
		if(args.length>1){
			String file=(args[0]);
			String save=args[1];
			new CombineBowtieFromVariousReferences(file, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=sam file (sorted by name) \n args[1]=save file";
	*/
}
