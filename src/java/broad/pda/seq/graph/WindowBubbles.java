package broad.pda.seq.graph;

import java.util.Collection;
import java.util.Stack;
import java.util.TreeSet;

import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;

public class WindowBubbles {

	//Stack<Bubble> plusBubbles;
	//Stack<Bubble> minusBubbles;
	Stack<Bubble> bubbles;
	String chr;
	
	public WindowBubbles(String chr){
		this.chr=chr;
		//plusBubbles=new Stack<Bubble>();
		//minusBubbles=new Stack<Bubble>();
		bubbles=new Stack<Bubble>();
	}
	
	public void addWindow(RefSeqGene gene){
		addWindow(gene, bubbles);
		/*if(gene.getOrientation().equalsIgnoreCase("+")){addWindow(gene, plusBubbles);}
		else if(gene.getOrientation().equalsIgnoreCase("-")){addWindow(gene, minusBubbles);}
		else{addWindow(gene, plusBubbles); addWindow(gene, minusBubbles);}*/
	}
	
	public void addWindow(RefSeqGene gene, Stack<Bubble> bubbles){
		//First pop the stack
		Bubble currentBubble=new Bubble();
		if(!bubbles.isEmpty()){currentBubble=bubbles.pop();}
		//see if overlapping bubble exists
		boolean overlaps=currentBubble.overlaps(gene);
		//If so, add to the current bubble and push
		if(overlaps){currentBubble.addWindow(gene); bubbles.push(currentBubble);}
		//Else make new bubble and add
		else{
			currentBubble.finalizeBubble();
			bubbles.push(currentBubble);
			Bubble newBubble=new Bubble();
			newBubble.addWindow(gene);
			bubbles.push(newBubble);
		}
	}
	
	/*public ChromosomeWithBubbles2[] getGraphs(){
		Collection<Alignments> exons=new TreeSet();
		Collection<Alignments> introns=new TreeSet();
		for(Bubble bubble: minusBubbles){
			if(!bubble.isFinalized()){bubble.finalizeBubble();}
			exons.addAll(bubble.getExons());
			introns.addAll(bubble.getIntrons());
		}
				
		ChromosomeWithBubbles2 minusGraph=new ChromosomeWithBubbles2(chr, exons, introns);
		
		exons=new TreeSet();
		introns=new TreeSet();
		for(Bubble bubble: plusBubbles){
			if(!bubble.isFinalized()){bubble.finalizeBubble();}
			exons.addAll(bubble.getExons());
			introns.addAll(bubble.getIntrons());
		}
				
		ChromosomeWithBubbles2 plusGraph=new ChromosomeWithBubbles2(chr, exons, introns);
		
		ChromosomeWithBubbles2[] rtrn={plusGraph, minusGraph};
		
		return rtrn;
	}*/
	
	public Collection<Alignments> getExons(){
		Collection<Alignments> exons=new TreeSet();
		for(Bubble bubble: bubbles){
			if(!bubble.isFinalized()){bubble.finalizeBubble();}
			exons.addAll(bubble.getExons());
		}
				
		return exons;
	}
	
	/*public Collection<RefSeqGene> getAllPaths(){
		ChromosomeWithBubbles2[] graphs=this.getGraphs();
		Collection<RefSeqGene> rtrn=new TreeSet();
		
		Collection<RefSeqGene> plus=graphs[0].getPaths(0);
		for(RefSeqGene gene: graphs[0].getPaths(0)){
			gene.setOrientation("+");
			rtrn.add(gene);
		}
		
		Collection<RefSeqGene> minus=graphs[1].getPaths(0);
		for(RefSeqGene gene: minus){
			gene.setOrientation("-");
			rtrn.add(gene);
		}
		
		rtrn.addAll(graphs[1].getPaths(0));
		
		
		Collection<RefSeqGene> orphansPlus=graphs[0].getOrphanNodes();
		Collection<RefSeqGene> orphansMinus=graphs[1].getOrphanNodes();
		
		for(RefSeqGene orphan: orphansPlus){
			if(orphansMinus.contains(orphan)){orphan.setOrientation("*"); rtrn.add(orphan);}
		}
		
		return rtrn;
	}*/

	public Collection<Alignments> getIntrons() {
		Collection<Alignments> introns=new TreeSet();
		for(Bubble bubble: bubbles){
			if(!bubble.isFinalized()){bubble.finalizeBubble();}
			introns.addAll(bubble.getIntrons());
		}
		return introns;
	}
	
}
