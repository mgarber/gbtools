package broad.pda.seq.alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Iterator;
import java.util.TreeSet;

import broad.pda.seq.alignment.sam.SAMRecord;

public class SortedSAMIterator implements Iterator<Collection<SAMRecord>> {

	BufferedReader reader;
	boolean done;
	String currentName;
	Collection<SAMRecord> records;
	boolean started;
	
	//TODO The records will eventually need to be a set
	public SortedSAMIterator(File samFile) throws FileNotFoundException{
		reader=new BufferedReader(new InputStreamReader(new FileInputStream(samFile)));
		done=false;
		this.currentName="";
		Collection<SAMRecord> records=new TreeSet<SAMRecord>();
		started=false;
	}
	
	public boolean hasNext() {
		return !done;
	}

	/*public Collection<SAMRecord> next() {
		if(done){return null;}
		Collection<SAMRecord> temp=new TreeSet();
		try {
			String nextLine;
			while ((nextLine = reader.readLine()) != null) {
				SAMRecord record=new SAMRecord(nextLine);
				System.err.println(record.getName()+" "+currentName+" "+record.getName().equalsIgnoreCase(currentName));
	        	if(record.getName().equalsIgnoreCase(currentName)){records.add(record); }
	        	else{temp=records; records=new ArrayList<SAMRecord>(); records.add(record); currentName=record.getName(); System.err.println("AM I HERE?"); return temp;}
	        }
			if(nextLine==null){done=true; System.err.println("HIT END"); close(); return records;}
		} catch (IOException e) {
			System.err.println("HIT END");
			done=true;
			close();
		}
		if(records.size()>1){System.err.println(records.size());}
		System.err.println(temp);
		return temp;
	}*/
	
	public Collection<SAMRecord> next(){
		String nextLine;
		
		Collection<SAMRecord> temp=new TreeSet<SAMRecord>();
        try {
			while ((nextLine = reader.readLine()) != null) {
				if(nextLine.startsWith("@"))  {
					// read until you pass header.
					continue;
				}
				//System.out.println(nextLine);
				SAMRecord record=new SAMRecord(nextLine);
				if(this.currentName==null || this.currentName.isEmpty()){
					this.currentName=record.getName();
					records=new TreeSet<SAMRecord>(); records.add(record);
				}
				else if(record.getName().equalsIgnoreCase(this.currentName)){records.add(record);}
				else{ 
					//we have a collection by name do something!!
					temp=copy(records);
					
					//reset for next run
					this.currentName=record.getName();
					records=new TreeSet<SAMRecord>(); records.add(record);
					
					//return the previous collection
					return temp;
				}
			}
			if(nextLine==null){done=true; close(); System.err.println("HIT END"); return records;}
		} catch (IOException e) {
			//Do nothing
			e.printStackTrace();
		}
		return temp;
	}

	private Collection<SAMRecord> copy(Collection<SAMRecord> records2) {
		if(records2==null){return null;}
		Collection<SAMRecord> rtrn=new TreeSet();
		
		for(SAMRecord record: records2){rtrn.add(record);}
		
		return rtrn;
	}

	private void close() {
		try {
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public void remove() {}

}
