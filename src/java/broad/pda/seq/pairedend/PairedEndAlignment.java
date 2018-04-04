package broad.pda.seq.pairedend;

import java.util.Set;
import java.util.TreeSet;

import org.broad.igv.sam.Alignment;

import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;

public class PairedEndAlignment implements Comparable{

	private RefSeqGene leftAlign;
	private RefSeqGene rightAlign;
	String name;
	
	public PairedEndAlignment(RefSeqGene left, RefSeqGene right, String name){
		this.leftAlign=left;
		this.rightAlign=right;
		this.name=name;
	}
	
	public PairedEndAlignment(Alignment left, Alignment right){
		this.leftAlign=new RefSeqGene(left);
		this.rightAlign=new RefSeqGene(right);
		this.name=left.getReadName();
	}
	
	public boolean onSameChromosome(){
		if(leftAlign.getChr().equalsIgnoreCase(rightAlign.getChr())){return true;}
		return false;
	}
	
	public boolean oppositeStrands(){
		if(!leftAlign.getOrientation().equalsIgnoreCase(rightAlign.getOrientation())){return true;}
		if(leftAlign.getOrientation().equalsIgnoreCase("*")){return true;}
		return false;
	}
	
	public int getDistance(){
		if(!onSameChromosome()){return -999;}
		int end=Math.min(leftAlign.getEnd(), rightAlign.getEnd());
		int start=Math.max(leftAlign.getStart(), rightAlign.getStart());
		return end-start;
	}
	
	public Alignments getEncompassingAlignment(){
		if(!onSameChromosome()){return null;}
		int start=Math.min(leftAlign.getStart(), rightAlign.getStart());
		int end=Math.max(leftAlign.getEnd(), rightAlign.getEnd());
		Alignments rtrn=new Alignments(leftAlign.getChr(), start, end);
		return rtrn;
	}
	
	public boolean pairOverlaps(){
		return (this.leftAlign.getAlignment().overlaps(rightAlign.getAlignment()));
	}
	
	public Alignments getInsertAlignment(){
		return new Alignments(this.leftAlign.getChr(), Math.min(leftAlign.getEnd(), rightAlign.getEnd()), Math.max(leftAlign.getStart(), rightAlign.getStart()));
	}
	
	private Set getExonPairs(){
		Set rtrn=new TreeSet();
		rtrn.addAll(rightAlign.getExonSet());
		rtrn.addAll(leftAlign.getExonSet());
		return rtrn;
	}
	
	private RefSeqGene getGeneStructure(){
		//TO DO: Add score!!
		
		Alignments align=this.getEncompassingAlignment();
		RefSeqGene gene=new RefSeqGene(align.getChr(), align.getStart(), align.getEnd(), name, "*", getExonPairs());
		return gene;
	}
	
	public String toString(){
		//return this.getGeneStructure().toBED()+"\t"+leftAlign.getOrientation()+"\t"+rightAlign.getOrientation()+"\t"+this.name;
		Alignments insert=this.getInsertAlignment();
		return insert.getChr()+"\t"+insert.getStart()+"\t"+insert.getEnd()+"\t"+"+";
	}
	
	//lower bound
	public int getCDNAInsertSize(){
		//sum all exons contained within the pair
		return sum(this.leftAlign.getExonSizes())+sum(this.rightAlign.getExonSizes());
	}
	
	private int sum(int[] vals){
		int sum=0;
		for(int i=0; i<vals.length; i++){sum+=vals[i];}
		return sum;
	}
	
	//best estimate
	//public int getInferredCDNAInsertSize(){}
	
	//write each record as a SAM
	//link to the paired end read
	/*public String[] toSAM(){
		String left=this.leftAlign.toSAM(this.rightAlign);
		String right=this.rightAlign.toSAM(this.leftAlign);
		String[] rtrn={left, right};
		return rtrn;
	}*/
	
	/*public String toSAM(){
		RefSeqGene gene=this.getGeneStructure();
		gene.setSequence(this.leftAlign.getSequence());
		//return this.leftAlign.toSAM(rightAlign)+"\n"+this.rightAlign.toSAM(leftAlign);
		return this.leftAlign.toSAM(rightAlign);
	}*/
	
	public String getName(){return this.name;}
	
	public int compareTo(Object b){
		PairedEndAlignment al=(PairedEndAlignment)b;
		return this.getEncompassingAlignment().compareTo(al.getEncompassingAlignment());
	}
}
