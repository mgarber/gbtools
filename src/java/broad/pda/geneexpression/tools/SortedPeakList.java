package broad.pda.geneexpression.tools;

import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;

/**
 * This class is a List implementation which sorts the peaks by their genomic location using the
 * comparator specified when constructing a new instance.
 * 
 * @param <T>
 */
public class SortedPeakList<Gene> extends LinkedList<Gene> {
    /**
     * Needed for serialization.
     */
    private static final long serialVersionUID = 1L;
    /**
     * Comparator used to sort the list.
     */
    private Comparator<Gene> comparator = null;
    /**
     * Construct a new instance with the list elements sorted in their
     * {@link java.lang.Comparable} natural ordering.
     */
    public SortedPeakList() {
    	super();
    }
    /**
     * Construct a new instance using the given comparator.
     * 
     * @param comparator
     */
    public SortedPeakList(Comparator<Gene> comparator) {
    	super();
        this.comparator = comparator;
    }
    
    /**
     * Check, if this list contains the given Element. This is faster than the
     * {@link #contains(Object)} method, since it is based on binary search.
     * 
     * @param paramT
     * @return <code>true</code>, if the element is contained in this list;
     * <code>false</code>, otherwise.
     */
    public boolean containsElement(Gene gene) {
        return (Collections.binarySearch(this, gene, comparator) > -1);
    }
    
    /**
     * Add a new entry to the list. The insertion point is calculated using the
     * comparator.
     * 
     * @param paramT
     */
    @Override
    public boolean add(Gene gene) {
        int insertionPoint = Collections.binarySearch(this, gene, comparator);
        super.add((insertionPoint > -1) ? insertionPoint : (-insertionPoint) - 1, gene);
        return true;
    }
    
    /**
     * Adds all elements in the specified collection to the list. Each element
     * will be inserted at the correct position to keep the list sorted.
     * 
     * @param paramCollection
     */
    /*   @Override
   public boolean addAll(Collection<Gene> paramCollection) {
        boolean result = false;
        for (Gene gene:paramCollection) {
            result |= add(gene);
        }
        return result;
    }*/

}