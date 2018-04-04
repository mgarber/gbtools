package broad.core.annotation;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

public class DiscontinuousAnnotation<G extends GenomicAnnotation> extends BasicGenomicAnnotation implements List<G>{
	List<G> parts;
	int maxDistToStitch;
	
	public DiscontinuousAnnotation(String name) {
		super(name);
		parts = new ArrayList<G>();
	}

	public DiscontinuousAnnotation(String name, List<G> list, int maxDistanceToStitch) {
		this(name);
		this.maxDistToStitch = maxDistanceToStitch;
		addAll(list);
	}

	
	public boolean contains(GenomicAnnotation other) {
		int pos = Collections.binarySearch(parts, other, new Comparator<GenomicAnnotation>() {

			public int compare(GenomicAnnotation arg0, GenomicAnnotation arg1) {
				return arg0.contains(arg1) ? 0 : arg0.compareTo(arg1);
			}
			
		});
		return pos >= 0 ;
	}

	public int getDistanceTo(GenomicAnnotation other) {
		GenomicAnnotation closest = getClosestPart(other);
		return closest.getDistanceTo(other);
	}

	public int getEnd() {
		return parts.size() > 0 ? parts.get(parts.size() - 1).getEnd() : 0;
	}


	public int getStart() {
		return parts.size() > 0 ? parts.get(0).getStart() : 0;
	}


	public List<GenomicAnnotation> minus(List<? extends GenomicAnnotation> others) {
		List<GenomicAnnotation> diff = new ArrayList<GenomicAnnotation>(parts.size());
		Iterator<G> partIt = parts.iterator();
		while(partIt.hasNext()) {
			diff.addAll(partIt.next().minus(others));
		}
		return diff;
	}

	public List<GenomicAnnotation> minus(GenomicAnnotation other) {
		List<GenomicAnnotation> oneList = new ArrayList<GenomicAnnotation>(1);
		oneList.add(other);
		return minus(oneList);
	}

	public boolean overlaps(GenomicAnnotation other, int buffer) {
		int pos = Collections.binarySearch(parts, other, new Comparator<GenomicAnnotation>() {

			public int compare(GenomicAnnotation arg0, GenomicAnnotation arg1) {
				return arg0.overlaps(arg1) ? 0 : arg0.compareTo(arg1);
			}
			
		});
		return pos >= 0 ;
	}

	public boolean overlaps(GenomicAnnotation other) {
		return overlaps(other, 0);
	}

	public void setEnd(int end) {
		//do nothing
	}
	
	public void setStart(int start) {
		//do nothing
	}



	public int getEffectiveLength() {
		int length = 0;
		Iterator<G> it = parts.iterator();
		while(it.hasNext()) {
			length += it.next().length();
		}
		
		return length;
	}
	public void takeIntersection(GenomicAnnotation other) {
		GenomicAnnotation closest = getClosestPart(other);
		if(closest.overlaps(other)) {
			closest.takeIntersection(other);
		}
	}

	public void takeUnion(GenomicAnnotation other) {
		GenomicAnnotation closest = getClosestPart(other);
		closest.takeUnion(other);
	}
	
	public  GenomicAnnotation getClosestPart(GenomicAnnotation annot) {
		GenomicAnnotation closest = null;
		if (size() > 0) {
			closest = get(getClosestPartIdx(annot));
		}
		return closest;
	}
	
	private int getClosestPartIdx(GenomicAnnotation annot) {
		int pos = Collections.binarySearch(parts, annot);
		int closestPos = pos;

		if(parts.size() > 0 && pos < 0) {
			int nextPos = Math.min(parts.size() - 1, -pos - 1);
			int priorPos  = Math.max(0, -pos - 2);
			
			GenomicAnnotation prior = parts.get(priorPos);
			GenomicAnnotation next  = parts.get(nextPos);
			
			closestPos = annot.getDistanceTo(prior) > annot.getDistanceTo(next) ? nextPos : priorPos;
		}
		
		return closestPos;
	}
	
	public List<G> getOverlappingParts(GenomicAnnotation annot) {
		ArrayList<G> overlaps = new ArrayList<G>(1);
		if(size() > 0) {
			int closestIdx = getClosestPartIdx(annot);
			ListIterator<G> listIt = listIterator(closestIdx);
			boolean overlapped = true;
			while(listIt.hasNext() && overlapped) {
				G a = listIt.next();
				if(a.overlaps(annot)) {
					overlaps.add(a);
				} else {
					overlapped = false;
				}
			}
			
			if(closestIdx - 1 >= 0) {
				listIt = listIterator(closestIdx - 1);
				overlapped = true;
				while(listIt.hasPrevious() && overlapped) {
					G a = listIt.previous();
					if(a.overlaps(annot)) {
						overlaps.add(a);
					} else {
						overlapped = false;
					}
				}
			}
		}
		
		Collections.sort(overlaps);
		return overlaps;
	}
	/*
	public boolean add(G annotation) {
		if(parts.size() == 0) {
			parts.add(annotation);
			return true;
		}
		
		
		int closestIdx = getClosestPartIdx(annotation);
		GenomicAnnotation closest = get(closestIdx);

		if(!closest.contains(annotation)) {
			parts.add(closest.compareTo(annotation) < 0 ? closestIdx + 1 : closestIdx, annotation);
		
			if(closest.overlaps(annotation,maxDistToStitch) ) {
				parts = BasicGenomicAnnotation.stitchList(parts, maxDistToStitch);
			}
		}
		
		return true;
	}
	 */
	public boolean add(G annotation) {
		return parts.add(annotation);
	}
	
	public void consolidate() {
		Collections.sort(parts);
		parts = BasicGenomicAnnotation.stitchList(parts, maxDistToStitch);
	}
	
	public void add(int arg0, GenomicAnnotation arg1) {
		//Even though it is mandatory mandatory to support, it does not make sense in this context to support it.
		throw new UnsupportedOperationException("The position of added genomic annotations are determined based on their natural ordering, they can't be set");	
	}

	public boolean addAll(Collection<? extends G> arg0) {
		Iterator<? extends G> it = arg0.iterator();
		while(it.hasNext()) {
			G a = it.next();
			add(a);
		}
		return true;
	}

	public boolean addAll(int arg0, Collection<? extends G> arg1) {
		//Even though it is mandatory mandatory to support, it does not make sense in this context to support it.
		throw new UnsupportedOperationException("The position of added genomic annotations are determined based on their natural ordering, they can't be set");	
	}

	public void clear() {
		parts.clear();
	}

	public boolean contains(Object arg0) {
		return parts.contains(arg0);
	}

	public boolean containsAll(Collection<?> arg0) {
		return parts.containsAll(arg0);
	}

	public G get(int i) {
		return parts.get(i);
	}

	public int indexOf(Object arg0) {
		return parts.indexOf(arg0);
	}

	public boolean isEmpty() {
		return parts.isEmpty();
	}

	public Iterator<G> iterator() {
		return parts.iterator();
	}

	public int lastIndexOf(Object arg0) {
		return parts.lastIndexOf(arg0);
	}

	public ListIterator<G> listIterator() {
		return parts.listIterator();
	}

	public ListIterator<G> listIterator(int i) {
		return parts.listIterator(i);
	}

	public boolean remove(Object arg0) {
		return parts.remove(arg0);
	}

	public G remove(int i) {
		return parts.remove(i);
	}

	public boolean removeAll(Collection<?> arg0) {
		return parts.removeAll(arg0);
	}

	public boolean retainAll(Collection<?> arg0) {
		return parts.retainAll(arg0);
	}

	public G set(int i, G arg1) {
		//This is not mandatory to support and it does not make sense in this context to support it.
		throw new UnsupportedOperationException("The position of added genomic annotations are determined based on their natural ordering, they can't be set");	
	}

	public int size() {
		return parts.size();
	}

	public List<G> subList(int i, int j) {
		return parts.subList(i, j);
	}

	public Object[] toArray() {
		return parts.toArray();
	}

	public <T> T[] toArray(T[] arg0) {
		return parts.toArray(arg0);
	}

}
