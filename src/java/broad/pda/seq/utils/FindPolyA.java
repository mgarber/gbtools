package broad.pda.seq.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;

import broad.pda.seq.fastq.FastqSequence;



public class FindPolyA {
	
	int polyN=10;
	int minLength=30;
	
	//find reads that end with PolyA or PolyT
	public FindPolyA(File fastq, String saveDir)throws IOException{
		parseAndWrite(fastq, saveDir);
		//Collection<FastqSequence> polyA=getAllSequencesEnding(sequences, polyN, "A");
		//Collection<FastqSequence> polyT=getAllSequencesEnding(sequences, polyN, "T");
		
	}
	
	private void parseAndWrite(File file, String saveDir)throws IOException{
		FileWriter lastT=new FileWriter(saveDir+"/"+file.getName()+".LastTs.fq");
		FileWriter lastA=new FileWriter(saveDir+"/"+file.getName()+".LastAs.fq");
		FileWriter firstT=new FileWriter(saveDir+"/"+file.getName()+".FirstTs.fq");
		FileWriter firstA=new FileWriter(saveDir+"/"+file.getName()+".FirstAs.fq");
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
    	int i=0;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	if(nextLine.startsWith("@")){
        		String firstLine=nextLine;
        		String secondLine=reader.readLine();
        		String thirdLine=reader.readLine();
        		String fourthLine=reader.readLine();
        		
        		FastqSequence seq=new FastqSequence(firstLine, secondLine, thirdLine, fourthLine);
        		String sequence=seq.getSequence();
        		String lastNBps=getLastBps(sequence, polyN);
        		String firstNBps=getFirstBps(sequence, polyN);
    			String polyT=polyN("T", polyN);
    			String polyA=polyN("A", polyN);
    			if(lastNBps.equalsIgnoreCase(polyT)){
    				FastqSequence trimmed=seq.trimEnds('T');
    				if(trimmed.getLength()>minLength){lastT.write(trimmed+"\n");}
    			}
    			if(lastNBps.equalsIgnoreCase(polyA)){
    				FastqSequence trimmed=seq.trimEnds('A');
    				if(trimmed.getLength()>minLength){lastA.write(trimmed+"\n");}
    			}
    			if(firstNBps.equalsIgnoreCase(polyT)){
    				FastqSequence trimmed=seq.trimBeginning('T');
    				if(trimmed.getLength()>minLength){firstT.write(trimmed+"\n");}
    			}
    			if(firstNBps.equalsIgnoreCase(polyA)){
    				FastqSequence trimmed=seq.trimBeginning('A');
    				if(trimmed.getLength()>minLength){firstA.write(trimmed+"\n");}
    			}
        	}
        	
        	
        }
        lastT.close();
        lastA.close();
        firstA.close();
        firstT.close();
	}
	
	
	private Collection<FastqSequence> getAllSequencesEnding(Collection<FastqSequence> sequences, int num, String letter){
		Collection rtrn=new ArrayList();
		
		for(FastqSequence seq: sequences){
			String sequence=seq.getSequence();
			String lastNBps=getLastBps(sequence, num);
			String polyN=polyN(letter, num);
			System.err.println(lastNBps.toCharArray().length+" "+polyN.toCharArray().length);
			if(lastNBps.equalsIgnoreCase(polyN)){rtrn.add(seq);}
		}
		
		return rtrn;
	}
	
	private String getLastBps(String sequence, int num){
		return sequence.substring(sequence.toCharArray().length-num);
	}
	
	private String getFirstBps(String sequence, int num){
		return sequence.substring(0, num);
	}
	
	private String polyN(String letter, int num){
		String rtrn="";
		for(int i=0; i<num; i++){rtrn=rtrn+letter;}
		return rtrn;
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>1){
			File fastq=new File(args[0]);
			String save=args[1];
			new FindPolyA(fastq, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=fastq file \n args[1]=save directory";
}
