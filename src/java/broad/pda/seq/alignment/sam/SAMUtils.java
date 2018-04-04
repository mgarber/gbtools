package broad.pda.seq.alignment.sam;

import java.util.ArrayList;
import java.util.Collection;

import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;

public class SAMUtils {

	public static RefSeqGene SAMFormatFullBED(String str){
		String[] tokens=str.split("\t");
    	String chr=tokens[2];
    	int start=new Integer(tokens[3])-1; //no +1
    	String strand=getStrand(new Integer(tokens[1]));
    	String cigar=tokens[5]; //parse cigar and put the breaks into their component parts, in parallel keep track of the splices
    	Collection<Alignments> aligns=parseCigar(cigar, chr, start);
		RefSeqGene gene=new RefSeqGene(aligns);
		gene.setName(tokens[0]);
		gene.setOrientation(strand);
		if(tokens.length>9){gene.setSequence(tokens[9]);}
    	return gene;
	}
	
	public static RefSeqGene SAMFormatFullBED(net.sf.samtools.SAMRecord picardSAMRecord){
		String chr=picardSAMRecord.getReferenceName();
    	int start=picardSAMRecord.getAlignmentStart()-1; //TODO: Do we need the -1, it is a left over from the method above.
    	String strand=picardSAMRecord.getReadNegativeStrandFlag() ? "-" : "+";
    	String cigar= picardSAMRecord.getCigarString(); //parse cigar and put the breaks into their component parts, in parallel keep track of the splices
    	Collection<Alignments> aligns=parseCigar(cigar, chr, start);
		RefSeqGene gene=new RefSeqGene(aligns);
		gene.setName(picardSAMRecord.getReadName());
		gene.setOrientation(strand);
		gene.setChromosome(chr);
		if(picardSAMRecord.getReadString() != null) {gene.setSequence(picardSAMRecord.getReadString()); }
    	return gene;
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
			//end is currentStart+num-1;
			if(flag==('M')){Alignments align=new Alignments(chr, currentStart, currentStart+num); currentStart=currentStart+num; rtrn.add(align);}
			else if(flag==('N')){currentStart=currentStart+num;} //currentStart=currentStart+num-1;
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
		//Old code that doesn't decode a more complex bit flag:
		
		//if(flag==16){return "-";}
		//else if(flag==0){return "+";}
		//else{return "*";}
		
		
		//proper way is:
		if ((flag & SAMRecord.getReadStrandFlag ())!= 0) {return "-";}
		else{return "+";}
	}

	public static boolean isValid(String nextLine) {
		if(nextLine.startsWith("@") || nextLine.split("\t")[2].equalsIgnoreCase("*")){return false;}
		return true;
	}

	
}
