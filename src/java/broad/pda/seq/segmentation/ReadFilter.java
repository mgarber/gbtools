package broad.pda.seq.segmentation;
import broad.pda.datastructures.Alignments;

public interface ReadFilter {
	/**
	 * Evaluates whether or not a given read alignment passes the filtering criteria.
	 * Filters can also modify the read itself by, for example setting orientation
	 * @param read, the read alignment
	 * @return true if it passes the filtering criteria.
	 */
	boolean passes(Alignments read);
}
