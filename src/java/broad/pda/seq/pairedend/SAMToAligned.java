package broad.pda.seq.pairedend;

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

import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;

public class SAMToAligned {
	static int strandMask=16;

	public SAMToAligned(File file, String save)throws IOException{
		readAndWrite(file, save);
	}
	
	private void readAndWrite(File file, String save)throws IOException{
		FileWriter writer=new FileWriter(save);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		
    	String nextLine;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	String[] tokens=nextLine.split("\t");
        	String chr=tokens[2];
        	int start=new Integer(tokens[3]);
        	int end=start+tokens[9].toCharArray().length;
        	String strand=getStrand(new Integer(tokens[1]));
        	String cigar=tokens[5]; //parse cigar and put the breaks into their component parts, in parallel keep track of the splices
        	Collection<Alignments> aligns=this.parseCigar(cigar, chr, start);
        	for(Alignments align: aligns){
        		writer.write(align+"\t"+strand+"\n");
        	}
        }
        writer.close();
	}
	
	private Collection<Alignments> parseCigar(String cigar, String chr, int start){
		ArrayList[] array=parseCigar(cigar); //type , numbers
		
		ArrayList<Character> type=array[0];
		ArrayList<Integer> numbers=array[1];
		
		Collection<Alignments> rtrn=new ArrayList();
		
		int currentStart=start;
		boolean shouldReturn=true;
		
		for(int i=0; i<type.size(); i++){
			char flag=type.get(i);
			Integer num=numbers.get(i);
			if(flag==('M')){Alignments align=new Alignments(chr, currentStart, currentStart+num); currentStart=currentStart+num; rtrn.add(align);}
			else if(flag==('N')){currentStart=currentStart+num;}
			else{shouldReturn=false;}
		}
		if(!shouldReturn){return new ArrayList();}
		return rtrn;
	}
	
	private static ArrayList[] parseCigar(String cigar){
		char[] chars=cigar.toCharArray();
		
		ArrayList type=new ArrayList();
		ArrayList numbers=new ArrayList();
		
		ArrayList ordered=new ArrayList();
		String str="";
		for(int i=0; i<chars.length; i++){
			if(chars[i]=='M' || chars[i]=='N' || chars[i]=='I' || chars[i]=='D' || chars[i]=='S' || chars[i]=='H' || chars[i]=='P'){numbers.add(new Integer(str)); type.add(chars[i]); str=str+chars[i]; ordered.add(str);  str="";}
			else{str=str+chars[i];}
		}
		
		//System.err.println(ordered);
		
		ArrayList[] rtrn={type, numbers};
		return rtrn;
	}
	
	private String getStrand(int flag){
		int val=flag&strandMask;
		if(val==flag){return "-";}
		return "+";
	}
	
	public static Map<String, Collection> parseSAMByChr(File file)throws IOException{
		Map<String, Collection> rtrn=new TreeMap();
			
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
    	int i=0;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	Collection<Alignments> aligns=GenericRNASeqLineParser.SAMFormat(nextLine);
           	String chr=((Alignments)aligns.toArray()[0]).getChr();
           	Collection all=new ArrayList();
        	if(rtrn.containsKey(chr)){all=rtrn.get(chr);}
        	all.addAll(aligns);
        	rtrn.put(chr, all);
        	i++;
        	if(i%10000 ==0){System.err.println(i);}
        }
                
       return rtrn;
	}
	
	public static Map<String, Collection> parseSAMByName(File file)throws IOException{
		Map<String, Collection> rtrn=new TreeMap();
			
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
    	int i=0;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	String[] tokens=nextLine.split("\t");
        	String name=tokens[0];
        	Collection<Alignments> aligns=GenericRNASeqLineParser.SAMFormat(nextLine);
           	Collection all=new ArrayList();
        	if(rtrn.containsKey(name)){all=rtrn.get(name);}
        	all.addAll(aligns);
        	rtrn.put(name, all);
        	i++;
        	if(i%10000 ==0){System.err.println(i);}
        }
                
       return rtrn;
	}
	
	
	public static Map<String, Collection> parseSAMByName(File file, String chrFilter)throws IOException{
		Map<String, Collection> rtrn=new TreeMap();
			
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
    	int i=0;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	String[] tokens=nextLine.split("\t");
        	String chr=tokens[2];
        	if(chr.equalsIgnoreCase(chrFilter)){
        	String name=tokens[0];
        	RefSeqGene gene=GenericRNASeqLineParser.SAMFormatFullBED(nextLine);
        	Collection all=new ArrayList();
        	if(rtrn.containsKey(name)){all=rtrn.get(name);}
        	all.add(gene);
        	rtrn.put(name, all);
        	}
        	i++;
        	if(i%1000000 ==0){System.err.println(i);}
        }
                
       return rtrn;
	}
	
	public static Collection parseSAMSplices(File file)throws IOException{
		Collection<Alignments> splices=new ArrayList();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		
    	String nextLine;
    	int i=0;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	Collection<Alignments> splicedAlignments=GenericRNASeqLineParser.SAMFormatSplices(nextLine);
        	splices.addAll(splicedAlignments);
        	i++;
        	if(i%10000 ==0){System.err.println(i);}
        }
		
		
		return splices;
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>1){
		File file=new File(args[0]);
		String save=args[1];
		new SAMToAligned(file, save);
		}
		else{System.err.println(usage);}
		
		
	}
	
	static String usage=" args[0]=file \n args[1]=save";
}
