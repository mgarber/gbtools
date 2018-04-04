package broad.pda.alignment;

import java.util.Set;
import java.util.TreeSet;

import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.alignment.Pair;
import broad.pda.seq.alignment.sam.SAMRecord;

public class PairedEndAlignment {

	private RefSeqGene leftAlign;
	private RefSeqGene rightAlign;
	String name;
	
	public PairedEndAlignment(RefSeqGene left, RefSeqGene right, String name){
		this.leftAlign=left;
		this.rightAlign=right;
		this.name=name;
	}
	
	public PairedEndAlignment(Pair<SAMRecord> pair) {
		this(pair.getValue1().getGene(), pair.getValue2().getGene(), pair.getValue1().getName());
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
	
	public Alignments getInsertAlignment(){
		if(!onSameChromosome()){return null;}
		int start=Math.min(leftAlign.getEnd(), rightAlign.getEnd());
		int end=Math.max(leftAlign.getStart(), rightAlign.getStart());
		Alignments rtrn=new Alignments(leftAlign.getChr(), start, end);
		return rtrn;
	}
	
	public Set getExonPairs(){
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
		return this.getGeneStructure().toBED()+"\t"+leftAlign.getOrientation()+"\t"+rightAlign.getOrientation();
	}
	
	public RefSeqGene getRightMate(){return this.rightAlign;}
	public RefSeqGene getLeftMate(){return this.leftAlign;}
}
