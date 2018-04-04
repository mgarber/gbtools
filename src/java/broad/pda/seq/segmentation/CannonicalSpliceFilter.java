package broad.pda.seq.segmentation;

import broad.core.sequence.Sequence;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.GeneTools;

public class CannonicalSpliceFilter implements ReadFilter {
	Sequence chromosome;
	
	public CannonicalSpliceFilter(Sequence chromosome) {
		this.chromosome = chromosome;
	}
	
	public boolean passes(Alignments read){
		String orientation=GeneTools.orientationFromSpliceSites(read, chromosome);
		read.setOrientation(orientation);
		if(orientation.equalsIgnoreCase("*")){return false;}
		return true;
	}

}
