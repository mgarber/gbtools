package broad.core.annotation;

import java.util.List;

import broad.core.sequence.Sequence;

public interface GenomicAnnotation extends Cloneable ,Feature, LightweightGenomicAnnotation {
	//String getDisplayName();
	void setReversedOrientation(boolean isInReversedOrientation);
	
	/**
	 * @deprecated Use length() instead
	 * @return
	 */
	int getLength();
	
	Sequence getSequence();
	void setSequence(Sequence seq);
	public void setOrientation(String orientationString);
	public String getOrientation();
	
	/**
	 * Genomic annotation integer identification if such exists.
	 */
	public String getId();
	public void setId(String id);
	
	/**
	 * @return int - the five prime buffer (i.e. bases not belonging to the annotation proper).
	 */
	public int getFivePrimeBases();

	/**
	 * @return int - the three prime buffer (i.e. bases not belonging to the annotation proper).
	 */
	public int getThreePrimeBases();
	
	/**
	 * Checks whether a twoSubjectAnnotation flanks the instance
	 */
	boolean isFlankedBy(TwoSubjectAnnotation twoSubjectAnnotation, int buffer);
	
	/**
	 * Returns the difference (all regions in this genomic annotation not in the given list)
	 * between this genomic annotation and the given list
	 * 
	 * @param others - the annotations to take the difference against
	 * @return A List of genomic annotations of the resulting difference.
	 */
	public List<GenomicAnnotation> minus(List<? extends GenomicAnnotation> others); 
	
	/**
	 * Returns the difference (all regions in this genomic annotation not in the given other annotation)
	 * between this genomic annotation and the one
	 * 
	 * @param other - the annotations to take the difference against
	 * @return A List of genomic annotations of the resulting difference.
	 */
	public List<GenomicAnnotation> minus(GenomicAnnotation other);
	
	/**
	 * Fragments the current annotation if it overlaps the provided one.
	 * it returns an empty list if this annotation does not ovelap the
	 * one passed in.
	 * 
	 * @param GenomicAnnotation annotation to disect the current one
	 * @return List<GenomicAnnotation> of the disected annotations, an empty list is returned if 
	 * 			the annotations do not overlap.
	 */
	List<GenomicAnnotation> disect(GenomicAnnotation a);
	
	
	/**
	 * Fragments the current annotation if it overlaps the provided ones.
	 * It returns a list with one component (this annotation) if no annotation 
	 * in the provided list overlaps the discted annotaion.
	 * 
	 * @param List<GenomicAnnotation> <b>sorted</b> annotations with which to disect the current one.
	 * @return List<GenomicAnnotation> of the disected annotations, a list just this annotation is returned if 
	 * 			the annotations do not overlap.
	 */
	List<GenomicAnnotation> disect(List<? extends GenomicAnnotation> disectors);
	
	
	/**
	 * Gets the start of the annotation considering its orientation 
	 * @return getStart() if the feature is in direct orientation or getEnd() otherwise
	 */
	int getOrientedStart();
	
	/**
	 * Gets the end of the annotation considering its orientation 
	 * @return getEnd() if the feature is in direct orientation or getStart() otherwise
	 */
	int getOrientedEnd();
	
	/**
	 * Returns a string of the form chr<chromosome>:start-end
	 * @return
	 */
	public String getLocationString();
	
	/**
	 * Returns true if the annotation (like a full BED or a RefSeq) has sub annotations like exons or blocks
	 */
	public boolean mayHaveBlocks();
	
	/**
	 * Returns a list of blocks if the annotations has any (@see containsBlocks)
	 */
	public List<? extends GenomicAnnotation> getBlocks();
	
	/**
	 * Adds a block to the annotation if it supports blocks
	 */
	public void addBlock(String name, int start, int end);
	
}
