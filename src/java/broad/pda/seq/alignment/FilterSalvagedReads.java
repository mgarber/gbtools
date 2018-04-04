package broad.pda.seq.alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.pda.datastructures.Alignments;
import broad.pda.seq.alignment.sam.SAMRecord;

public class FilterSalvagedReads {

	double minScore=65;
	private int minIntronSize=77; //1% of all ref seq introns
	private int maxIntronSize=73204; //1% of all ref seq introns
	private int minExonSize=10; //empirical
	
	
	public FilterSalvagedReads(File sam, String save) throws IOException{
		//FileWriter writer=new FileWriter(save);
		//writeAndFilter(sam, writer, "salvaged");
		//writer.close();
		getReadsByName(sam.getAbsolutePath(), save);
	}
	
	//TODO Implement so we can get unique best mapper if exists
	private void getReadsByName(String samFile, String save) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(samFile)));
		FileWriter writer=new FileWriter(save);
		
		int i=0;
		
		String nextLine;
		String name="";
		Collection<SAMRecord> records=new ArrayList<SAMRecord>();
        while ((nextLine = reader.readLine()) != null) {
        	SAMRecord record=new SAMRecord(nextLine);
        	if(record.getName().equalsIgnoreCase(name)){records.add(record);}
        	else{ 
        		//we have a collection by name do something!!
        		Collection<SAMRecord> salvaged=getMinValueRecords(records);
        		for(SAMRecord s: salvaged){
        			s.setWeight((double)1/salvaged.size());
        			if(salvaged.size()>1){s.setMappingQuality(0);}
        			else{s.setMappingQuality(255);} //TODO This will only write unique
        			writer.write(s+"\n");
        		}
        		
        		records=new ArrayList<SAMRecord>(); records.add(record);
        	}
        	name=record.getName();
        	//mark the transition in read name
        	i++;
        	if(i%100000 ==0){System.err.println(i);}
        }
		reader.close();		
		writer.close();
	}
	
	private Collection<SAMRecord> getMinValueRecords(Collection<SAMRecord> records) {
		//first lets filter reads
		Collection<SAMRecord> temp=new ArrayList();
		for(SAMRecord record: records){
			boolean filtered=filterRecord(record);
			if(!filtered){temp.add(record);}
		}
		records=temp;
		
		
		Map<SAMRecord, Double> scores=new TreeMap<SAMRecord, Double>();
		
		for(SAMRecord record: records){
			double score=(record.getLength()-record.getNumMismatches())-record.getNumMismatches();
			scores.put(record, score);
		}
		
		//if there is one with fewest mismatches return this
		Collection<SAMRecord> minScores=getMaxAlignments(scores);
		return minScores;
	}
	
	private Collection<SAMRecord> getMaxAlignments( Map<SAMRecord, Double> scores) {
		Collection<SAMRecord> rtrn=new TreeSet<SAMRecord>();
		
		double max=-Double.MAX_VALUE;
		
		for(SAMRecord record: scores.keySet()){max=Math.max(scores.get(record), max);}
		
		for(SAMRecord record: scores.keySet()){
			double score=scores.get(record);
			if(score==max){rtrn.add(record);}
		}
		
		return rtrn;
	}
	
	//If has gap then has to be minimum size and less than maximum size
	//score has to be larger than cutoff
	private boolean filterRecord(SAMRecord record) {
		double score=(record.getLength()-record.getNumMismatches())-record.getNumMismatches();
		if(score<this.minScore){return true;}
		
		Collection<Alignments> introns=record.getGene().getIntronSet();
		//all introns have to be larger than min
		for(Alignments intron: introns){
			if(intron.getSize()<this.minIntronSize || intron.getSize()>this.maxIntronSize){return true;}
		}
		
		//all exons need to be larger than min exon
		Collection<Alignments> exons=record.getGene().getExonSet();
		for(Alignments exon: exons){
			if(exon.getSize()<this.minExonSize){return true;}
		}
		
		return false;
	}
	
	private void writeAndFilter(File sam, FileWriter writer,String string) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(sam)));
		
		int passed=0;
		int counter=0;
		int count=0;
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			try{
				SAMRecord record=new SAMRecord(nextLine);
				record.setName(record.getName().split("/")[0]);
				record.setProgramName(string);
				boolean filtered=filterRecord(record);
				if(!filtered){writer.write(record+"\n"); passed++;}
				//else{System.err.println(record);}
				counter++;
				if(counter % 100000 ==0){System.err.println(counter+" "+(double)passed/counter);}
			}catch(Exception ex){}
		}
					
		//parse SAM file
		//if valid SAM line write it
		//else skip it
		
	}
	
	public static void main(String[] args) throws IOException{
		File in=new File(args[0]);
		String save=args[1];
		new FilterSalvagedReads(in, save);
	}
	
}
