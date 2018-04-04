package broad.pda.seq.utils;

import broad.pda.datastructures.Alignments;

public class AlignmentsPair {

	private final Alignments alignment1;
	private final Alignments alignment2;
	
	public AlignmentsPair(Alignments alignment1, Alignments alignment2) {
		this.alignment1 = alignment1;
		this.alignment2 = alignment2;
	}
	
	public Alignments getAlignment1() {
		return this.alignment1;
	}
	
	public Alignments getAlignment2() {
		return this.alignment2;
	}
	
	public boolean hasAlignment(Alignments a) {
		return getAlignment1().equals(a) || getAlignment2().equals(a);
	}
	
	public int hashCode() {
		return alignment1.hashCode() + alignment2.hashCode();
	}
	
	public boolean equals(Object o) {
		AlignmentsPair ap = (AlignmentsPair) o;
		return ((getAlignment1().equals(ap.getAlignment1()) && getAlignment2().equals(ap.getAlignment2()) || (getAlignment1().equals(ap.getAlignment2()) && getAlignment2().equals(ap.getAlignment1()))));
	}
	
}
