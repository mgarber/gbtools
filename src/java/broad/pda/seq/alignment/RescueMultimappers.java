package broad.pda.seq.alignment;

//Goal of this class is to take all reads that are multimappers and try to resolve a unique alignment from them
//The way we'll do this is by using the paired end information to get the best placement of reads along with their pairs
//If a pair is uniquely mapped and their is ambiguoty in the other read then take the alignment that is likely given the mate alignment
//If both pairs are not uniquely mapped and their exist a mapping for the 2 that are consistent (and unique) take this

public class RescueMultimappers {

	public RescueMultimappers(){}
	
}
