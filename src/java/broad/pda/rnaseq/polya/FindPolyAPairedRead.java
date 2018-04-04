package broad.pda.rnaseq.polya;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.TreeSet;

import broad.pda.seq.fastq.FastqParser;
import broad.pda.seq.fastq.FastqSequence;

public class FindPolyAPairedRead {

	private int polyN=20;

	//TODO Flag polyAs and put them in a seperate file
	//Make 2 files Pairs of PolyA and trimmed polyA
	public FindPolyAPairedRead(File fastqFileP1, File fastqFileP2, String save, int polyN) throws IOException{
		this.polyN=polyN;
		FastqParser fastq=new FastqParser(fastqFileP1);
		
		Collection<String> polyAReads=getPolyASeq(fastq);
		fastq=new FastqParser(fastqFileP2);
		polyAReads.addAll(getPolyASeq(fastq));
		
		writeReads(save, fastqFileP1, fastqFileP2, polyAReads);
		//Collection<String> alignmentsOfPairs=getAlignedPairs(samFile, polyAReads);
		//write(save, alignmentsOfPairs);
	}
		
	private void writeReads(String save, File fastqFileP1, File fastqFileP2, Collection<String> polyAReads) throws IOException {
		FileWriter writer1=new FileWriter(save+".Pair1.fq");
		FileWriter writer2=new FileWriter(save+".Pair2.fq");
		FastqParser fastq=new FastqParser(fastqFileP1);
		while(fastq.hasNext()){
			FastqSequence seq=fastq.next();
			if(seq!=null){
				if(polyAReads.contains(seq.getName().split("/")[0].replaceAll("@", ""))){
					if(seq.isPartialPolyA(polyN)){
						seq=seq.trimPolyA();
					}
					writer1.write(seq.toFastq()+"\n");
				}
			}
		}
		fastq=new FastqParser(fastqFileP2);
		while(fastq.hasNext()){
			FastqSequence seq=fastq.next();
			if(seq!=null){
			if(polyAReads.contains(seq.getName().split("/")[0].replaceAll("@", ""))){writer2.write(seq.toFastq()+"\n");}
			}
		}
		writer1.close();
		writer2.close();
	}

	private Collection<String> getAlignedPairs(String samFile, Collection<String> polyAReads) throws IOException {
		Collection<String> rtrn=new TreeSet<String>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(samFile)));
		String nextLine;
		int total=0;
		int polyA=0;
		while ((nextLine = reader.readLine()) != null) {
			String name=nextLine.split("\t")[0];
			if(polyAReads.contains(name)){rtrn.add(nextLine); polyA++; System.out.println(nextLine);}
			total++;
			if(total%100000 ==0){System.err.println(total+" "+polyA+" "+((double)polyA/total));}
		}
		reader.close();
		return rtrn;
	}

	private void write(String save, Collection<String> alignmentsOfPairs) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String seq: alignmentsOfPairs){
			writer.write(seq+"\n");
		}
		
		writer.close();
	}

	private Collection<String> getPolyASeq(FastqParser fastq){
		Collection<String> polyA=new TreeSet<String>();
		
		int total=0;
		int polyANum=0;
		while(fastq.hasNext()){
			FastqSequence seq=fastq.next();
			if(seq!=null){
				boolean polyASeq=seq.isPolyA();
				boolean partialPolyASeq=seq.isPartialPolyA(polyN);
				if(polyASeq || partialPolyASeq){polyA.add(seq.getName().split("/")[0].replaceAll("@", "")); polyANum++;}
			}
			else{System.err.println("NULL");}
			total++;
			if(total%100000 == 0){System.err.println(total+" "+polyANum+" "+((double)polyANum/total));}
		}
		return polyA;
	}
	
	
	public static void main(String[] args)throws IOException{
		if(args.length>3){
			File fastq1=new File(args[0]);
			File fastq2=new File(args[1]);
			String save=args[2];
			int polyN=new Integer(args[3]);
			new FindPolyAPairedRead(fastq1, fastq2, save, polyN);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=pair1 \n args[1]=pair2 \n args[2]=save \n args[3]=poly N";
	
}
