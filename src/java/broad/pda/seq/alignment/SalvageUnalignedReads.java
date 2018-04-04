package broad.pda.seq.alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.TreeSet;

import broad.core.util.PipelineUtils;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.alignment.sam.SAMRecord;

public class SalvageUnalignedReads {

	double minScore=65;
	private int minIntronSize=77; //1% of all ref seq introns
	private int maxIntronSize=73204; //1% of all ref seq introns
	private int minExonSize=10;
	private int salvagedCount;
	
	public SalvageUnalignedReads(File blatAlignedSAM, File tophatAlignedSAM, File bowtieAlignedSAM, String save) throws IOException, InterruptedException{
		Runtime run=java.lang.Runtime.getRuntime();
		
		System.err.println("Concatanating....");
		//make a merged file with all records concatanated
		//append method description to row
		FileWriter writer=new FileWriter(save+".temp");
		//TODO Filter BLAT Reads
		//TODO Filter unique strata if possible
		writeAndFilter(blatAlignedSAM, writer, "BLAT");
		if(tophatAlignedSAM!=null){write(tophatAlignedSAM, writer, "Tophat");}
		write(bowtieAlignedSAM, writer, "Bowtie");
		writer.close();
		
		System.err.println("sorting ....");
		//sort by name
		String command="sort "+save+".temp"+" -o "+save+".temp.sorted.sam -T /broad/shptmp/";
		PipelineUtils.bsubProcess(run, command);
		System.err.println("writing....");
		//get collection by name
		this.salvagedCount=getReadsByName(save+".temp.sorted.sam", save);
		
		//if no bowtie then work on it
		
	}

	public int getSalvagedCount(){return this.salvagedCount;}
	
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
				if(counter % 100000 ==0){System.err.println((double)passed/counter);}
			}catch(Exception ex){}
		}
					
		//parse SAM file
		//if valid SAM line write it
		//else skip it
		
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

	private int getReadsByName(String samFile, String save) throws IOException {
		int counter=0;
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
        		Collection<SAMRecord> salvaged=salvage(records);
        		for(SAMRecord s: salvaged){
        			s.setWeight((double)1/salvaged.size());
        			if(salvaged.size()>1){s.setMappingQuality(0);}
        			else{s.setMappingQuality(255); counter++;}
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
		return counter;
	}

	//This is the method to salvage reads
	private Collection<SAMRecord> salvage(Collection<SAMRecord> records) {
		//if any of these reads came from bowtie then return an empty list;
		for(SAMRecord r: records){
			//System.err.println(r.getProgramName());
			if(r.getProgramName().equalsIgnoreCase("Bowtie")){return new TreeSet();}
		}
		//for now just report them all
		return records;
	}

	private void write(File sam, FileWriter writer, String string) throws IOException {
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(sam)));
			
			int counter=0;
			int count=0;
			String nextLine;
			while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
				try{
				//if(!nextLine.startsWith("@")){
					//try{
					SAMRecord record=new SAMRecord(nextLine);
					record.setName(record.getName().split("/")[0]);
					record.setProgramName(string);
					writer.write(record+"\n");
					/*String[] tokens=nextLine.split("\t");
					String match=tokens[2];
					if(!match.equalsIgnoreCase("*")){
						writer.write(nextLine+"\tPG:Z:"+string+"\n");
						count++;
					}*/
					//}catch(Exception ex){}
				//}
				counter++;
				//if(counter%100000 ==0){System.err.println(counter+" Good: "+count);}
				}catch(Exception ex){}
			}
						
			//parse SAM file
			//if valid SAM line write it
			//else skip it
		
		
	}
	
	public static void main(String[] args)throws IOException, InterruptedException{
		if(args.length>3){
		File blat=new File(args[0]);
		File tophat=new File(args[1]);
		File bowtie=new File(args[2]);
		String save=args[3];
		new SalvageUnalignedReads(blat, tophat, bowtie, save);
		}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=blat \n args[1]=tophat \n args[2]=bowtie \n args[3]=save";
}
