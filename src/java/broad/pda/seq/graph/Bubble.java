package broad.pda.seq.graph;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;

import broad.core.util.CollapseByIntersection;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;

public class Bubble {

	Stack<Alignments> exons;
	//Stack<Alignments> collapsedExons;
	Collection<Alignments> finalizedExons;
	Collection<Alignments> introns;
	int min;
	int max;
	boolean finalized;
	String save="/seq/rinnscratch/mguttman/RNASeq/JoshuaData/TopHat/test";
	
	public Bubble(){
		//exons=new TreeSet<Alignments>();
		//collapsedExons=new Stack<Alignments>();
		exons=new Stack<Alignments>();
		introns=new TreeSet<Alignments>();
		min=Integer.MAX_VALUE;
		max=-Integer.MAX_VALUE;
		finalized=false;
	}
	
	//TODO have it collapse and extract appropriately
	public boolean finalizeBubble(){
		//Go through exons and collapse all
		//Collection<Alignments> collapsedExons=BasicLightweightAnnotation.stitchList(exons, 0);
		//Collection<Alignments> collapsedExons=CollapseByIntersection.CollapseByIntersection(exons, false);
		//collapseExons(exons);
		//Go through each intron and add truncated exons to the set
		//Make sure it leaves all the originals as well
		
		Collection<Alignments> collapsedExons=CollapseByIntersection.collapseByIntersection(exons, false);
		//this.collapseExons(exons);
		finalizedExons=CollapseByIntersection.DecollapseByIntronLocation(collapsedExons, introns);
		
		try {
			//write(save+".uncollapsed", this.exons);
			write(save+".collapsed", this.exons);
			write(save+".finalized", this.finalizedExons);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		this.finalized=true;
		return false;
	}
	
	private void write(String string, Collection<Alignments> exons2) throws IOException {
		FileWriter writer=new FileWriter(string, true);
		
		for(Alignments align: exons2){
			writer.write(align+"\n");
		}
		
		writer.close();
	}

	public void addWindow(RefSeqGene gene){
		//exons.addAll(gene.getExonSet());
		addExons(gene.getExonSet());
		introns.addAll(gene.getIntronSet());
	}

	/*private void addExons(Set<Alignments> exonSet) {
		for(Alignments exon: exonSet){
			min=Math.min(min, exon.getStart());
			max=Math.max(max, exon.getEnd());
			exons.add(exon);
		}
	}*/
	
	private void addExons(Set<Alignments> exonSet) {
		Stack<Alignments> interim = new Stack<Alignments>(); //holds exons that are less than the one on top of the stack
		for(Alignments exon: exonSet){
			if(exons.isEmpty()){
				exons.push(exon);
			}
			else{
				Alignments previous=this.exons.pop();
				Alignments workingExon = new Alignments(exon); //Clone exon so we can modified it later
				if(!interim.isEmpty()) { // First check whether there is an interim exon
					Alignments interimExon = interim.pop();
					if(exon.overlaps(interimExon)) { 
						workingExon = new Alignments(exon.getChr(), Math.min(exon.getStart(), interimExon.getStart()), Math.max(exon.getEnd(), interimExon.getEnd())); //If overlaps interim then merge
					} else { //Otherwise interim is complete put it before the top of the stack.
						Alignments greaterExon = exons.pop();
						exons.push(interimExon);
						exons.push(greaterExon);
						workingExon = exon;
					}
				}
				//Now go as if exon is good
				if(previous.overlaps(workingExon)){
					Alignments newAlign=new Alignments(workingExon.getChr(), Math.min(workingExon.getStart(), previous.getStart()), Math.max(workingExon.getEnd(), previous.getEnd()));
					this.exons.push(newAlign);
				}
				else{
					this.exons.push(previous);
					if(workingExon.compareTo(previous)< 0) {
						interim.push(workingExon); //Check whether exon is less than top of the stack, if it is put it in the interim stack.
					} else {
						this.exons.push(workingExon);
					}
				}
				if(interim.size() > 1) {
					throw new IllegalArgumentException("BUG: interim stack must have at most 1 element");
				}
			}
			min=Math.min(min, exon.getStart());
			max=Math.max(max, exon.getEnd());
		}
	}

	/*private void collapseExons(Collection<Alignments> exons2) {
		for(Alignments exon: exons2){
			if(collapsedExons.isEmpty()){
				collapsedExons.push(exon);
			}
			else{
				Alignments previous=this.collapsedExons.pop();
				//Now go as if exon is good
				if(previous.overlaps(exon)){
					Alignments newAlign=new Alignments(exon.getChr(), Math.min(exon.getStart(), previous.getStart()), Math.max(exon.getEnd(), previous.getEnd()));
					this.collapsedExons.push(newAlign);
				}
				else{
					this.collapsedExons.push(previous);
					this.collapsedExons.push(exon);
					
				}
			}
			min=Math.min(min, exon.getStart());
			max=Math.max(max, exon.getEnd());
		}
	}*/
	
	public boolean overlaps(RefSeqGene gene) {
		if(min==Integer.MAX_VALUE || max==-Integer.MAX_VALUE){return false;}
		Alignments temp=new Alignments(gene.getChr(), min, max);
		return temp.overlaps(gene.getAlignment());
	}

	public boolean isFinalized() {
		return this.finalized;
	}

	public Collection<Alignments> getExons() {
		return this.finalizedExons;
	}
	
	public Collection<Alignments> getIntrons(){
		return this.introns;
	}
	
	
}
