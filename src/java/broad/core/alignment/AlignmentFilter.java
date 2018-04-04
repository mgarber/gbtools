package broad.core.alignment;

public  interface  AlignmentFilter <T extends AlignmentSummary> {
	/**
	 * Evaluates whether the alignment passes this filter criteria
	 * @param alignment
	 * @return true if the alignment passes filer criteria.
	 */
	boolean accept(T alignment);
}
