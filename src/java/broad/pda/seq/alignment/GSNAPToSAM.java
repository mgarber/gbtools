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

import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.alignment.sam.SAMRecord;

public class GSNAPToSAM {

	private int maxNumHits=10;

	public GSNAPToSAM(File file, String save) throws IOException{
		int i=0;
		FileWriter writer=new FileWriter(save);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
        while ((nextLine = reader.readLine()) != null) {
        	if(nextLine.startsWith(">")){
        		String[] tokens=nextLine.split("\t");
        		String seq=tokens[0].replaceAll(">", "");
        		int numHits=new Integer(tokens[1]);
        		String readName=tokens[2];
        		ArrayList<String> lines=new ArrayList<String>();
        		
        		//keep going until hit a space
        		while((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)){
        			lines.add(nextLine);
        		}
        		
        		Collection<SAMRecord> sam=makeSAMRecord(seq, numHits, readName, lines);
        		if(numHits<maxNumHits){
        			for(SAMRecord record: sam){
        				writer.write(record+"\n");
        			}
        		}
        		i++;
        		if(i%100000 ==0){System.err.println(i);}
        	}
        }
        reader.close();
        writer.close();
	}

	private Collection<SAMRecord> makeSAMRecord(String seq, int numHits, String readName, ArrayList<String> lines) {
		Collection<SAMRecord> rtrn=new ArrayList<SAMRecord>();
		
		double weight=1.0/numHits;
		
		for(int i=0; i<lines.size(); i++){
			String line=lines.get(i);
			Collection<Alignments> exons=new TreeSet<Alignments>();
			Alignments align=makeAlignment(line);
			exons.add(align);
			boolean add=true;
			int sub=getSub(line);
			if(i+1<lines.size() && lines.get(i+1).startsWith(",")){
				//then we have a splice
				Alignments exon1=makeAlignment(lines.get(i+1));
				if(exon1.getChr().equalsIgnoreCase(align.getChr())){
					exons.add(exon1);
					sub+=getSub(lines.get(i+1));
				}
				else{add=false;}
				i=i+1;
			}
			RefSeqGene gene=new RefSeqGene(exons);
			gene.setName(readName);
			SAMRecord record=new SAMRecord(readName, gene, seq);
			record.setNumMismatches(sub);
			record.setWeight(weight);
			if(numHits==1){record.setMappingQuality(255);}
			else{record.setMappingQuality(0);}
			if(add){rtrn.add(record);}
		}
		
		/*tokens=nextLine.split("\t");
		String alignmentString=tokens[0];
		String match=tokens[1];
		String score=tokens[2];
		String position=tokens[3];*/
		
		
		return rtrn;
	}

	private int getSub(String string) {
		String[] tokens=string.split(",");
		for(int i=0; i<tokens.length; i++){
			if(tokens[i].startsWith("sub:")){
				return new Integer(tokens[i].split(":")[1]);
			}
		}
		return 0;
	}

	private Alignments makeAlignment(String line) {
		
		String info=line.split("\t")[3];
		String chr="chr"+info.substring(1).split(":")[0];
		
		
		//try{
		//String chr="chr"+info.split(":")[0].replaceAll("+", "").replaceAll("-", "");
		int start=new Integer(info.split(":")[1].split("\\.\\.")[0]);
		int end=new Integer(info.split(":")[1].split("\\.\\.")[1]);
		//System.err.println(info+" "+chr+" "+start+" "+end);
		return new Alignments(chr, Math.min(start, end)-1, Math.max(start, end));
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>1){
			File in=new File(args[0]);
			String save=args[1];
			new GSNAPToSAM(in, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=in \n args[1]=out";
	
}
