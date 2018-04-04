package broad.pda.seq.pairedend;

import java.util.ArrayList;
import java.util.Collection;

import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;


public class GenericRNASeqLineParser {

	static int strandMask=16;
	
	public static Collection<Alignments> alignedFormat(String str){
		Collection c=new ArrayList();
		String[] tokens=str.split("\t");
		String chr=tokens[0];
		int start=new Integer(tokens[1]);
		int end=new Integer(tokens[2]);
		String strand=tokens[3];
		Alignments align=new Alignments(chr, start, end, strand);
		c.add(align);
		return c;
	}
	
	public static Collection<Alignments> SAMFormat(String str){
		String[] tokens=str.split("\t");
    	String chr=tokens[2];
    	int start=new Integer(tokens[3]);
    	int end=start+tokens[9].toCharArray().length;
    	String strand=getStrand(new Integer(tokens[1]));
    	String cigar=tokens[5]; //parse cigar and put the breaks into their component parts, in parallel keep track of the splices
    	Collection<Alignments> aligns=parseCigar(cigar, chr, start);
		return aligns;
	}
	
	public static RefSeqGene SAMFormatFullBED(String str){
		String[] tokens=str.split("\t");
    	String chr=tokens[2];
    	int start=new Integer(tokens[3]);
    	int end=start+tokens[9].toCharArray().length;
    	String strand=getStrand(new Integer(tokens[1]));
    	String cigar=tokens[5]; //parse cigar and put the breaks into their component parts, in parallel keep track of the splices
    	Collection<Alignments> aligns=parseCigar(cigar, chr, start);
		String sequence=tokens[9];
    	RefSeqGene gene=new RefSeqGene(chr, start, end, tokens[0], strand, aligns);
		gene.setSequence(sequence);
    	gene.setSAMString(str);
    	return gene;
	}
	
	public static Collection<Alignments> SAMFormatSplices(String str){
		String[] tokens=str.split("\t");
    	String chr=tokens[2];
    	int start=new Integer(tokens[3]);
    	int end=start+tokens[9].toCharArray().length;
    	String strand=getStrand(new Integer(tokens[1]));
    	String cigar=tokens[5]; //parse cigar and put the breaks into their component parts, in parallel keep track of the splices
    	Collection<Alignments> aligns=parseCigarSplices(cigar, chr, start);
		return aligns;
	}
	
	private static Collection<Alignments> parseCigarSplices(String cigar, String chr, int start){
		ArrayList[] array=parseCigar(cigar); //type , numbers
		
		ArrayList<Character> type=array[0];
		ArrayList<Integer> numbers=array[1];
		
		ArrayList<Alignments> rtrn=new ArrayList();
		
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
		
		Collection splices=new ArrayList();
		for(int i=0; i<rtrn.size()-1; i++){
			Alignments first=rtrn.get(i);
			Alignments next=rtrn.get(i+1);
			Alignments junction=new Alignments(first.getChr(), first.getEnd(), next.getStart());
			splices.add(junction);
		}
		
		return splices;
	}
	
	private static Collection<Alignments> parseCigar(String cigar, String chr, int start){
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
	
	private static String getStrand(int flag){
		if(flag==16){return "-";}
		else if(flag==0){return "+";}
		else{return "*";}
	}
	
}
