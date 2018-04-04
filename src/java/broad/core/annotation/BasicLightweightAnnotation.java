package broad.core.annotation;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Stack;
import java.util.TreeSet;


public class BasicLightweightAnnotation implements LightweightGenomicAnnotation{

	private int start;
	private int end;
	private String name;
	private String orientation;
	private String chromosome;
	private double score;
	private List<Double> extraScores;

	public BasicLightweightAnnotation() {
		super();
		extraScores = new ArrayList<Double>();
	}
	
	public BasicLightweightAnnotation(String chr, int start, int end) {
		this(chr, start, end, "*");
	}
	
	public BasicLightweightAnnotation(String chr, int start, int end, String orientation) {
		this.start = start;
		this.end   = end;
		setChromosome(chr);
		this.orientation = orientation;
		extraScores = new ArrayList<Double>();
	}
	
	public BasicLightweightAnnotation(String chr, int start, int end, String orientation, double Scr) {
		this(chr,start,end,orientation);
		score=Scr;
		
	}
	public BasicLightweightAnnotation(String chr, String start, String end) {
		setChromosome(chr);
		this.start=new Integer(start);
		this.end=new Integer(end);
		setOrientation(".");
		extraScores = new ArrayList<Double>();
	}

	public BasicLightweightAnnotation(LightweightGenomicAnnotation annotation) {
		start = annotation.getStart();
		end   = annotation.getEnd();
		name  = annotation.getName();
		setOrientation(annotation.getOrientation());
		setChromosome(annotation.getChromosome());
		score = annotation.getScore();
		extraScores = new ArrayList<Double>();
		for(double score : annotation.getExtraScores()) {
			extraScores.add(score);
		}
	}

	public String getChromosome() {
		return chromosome;
	}

	public String getChromosomeString() {
		if(getChromosome() == null) {
			return null;
		} if(getChromosome().length() < 3) {
			return "chr" + getChromosome();
		} else {
			return getChromosome();
		}

	}

	/**
	 * public void 
	 */
	
	public int getDistanceTo(LightweightGenomicAnnotation other) {
		int dist = 0;
		if(!getChromosome().equals(other.getChromosome())) {
			dist = 1000000000;
		}else if(!overlaps(other)) {
			if(getStart() >= other.getEnd()) {
				dist = getStart() - other.getEnd() ;
			} else {
				dist = other.getStart() - getEnd() ;
			}
		}
		return dist;
	}

	public int getEnd() {
		return end;
	}

	public double getExtraScore(int i) {
		return extraScores.get(i);
	}

	public List<Double> getExtraScores() {
		return extraScores;
	}

	public String getLocationString() {
		return getChromosomeString()+":"+getStart()+"-"+getEnd() ;
	}

	public long getMiddle() {
		return Math.round( (getStart() + getEnd())/(float)2);
	}

	public String getName() {
		return name == null ? getLocationString() : name;
	}

	public String getOrientation() { return orientation == null ? "*" : orientation;}

	public int getOverlap(LightweightGenomicAnnotation other) {
		int overlap = 0;
		if(overlaps(other)) {
			overlap = Math.min(getEnd(), other.getEnd()) - Math.max(getStart(), other.getStart());
		}
		return overlap;
	}


	public double getScore() { return score;}

	public int getStart() {
		return start;
	}

	public void setChromosome(String chromosome) {
		if(chromosome != null) {
			this.chromosome = chromosome.intern();
		}
	}

	public void setEnd(int end) {
		this.end = end;
	}

	protected void setEnd(String data) {
		this.end = Integer.parseInt(data);
	}

	public void setName(String name) {
		this.name = name;
	}

	public void setOrientation(String orientationString) {
		this.orientation = orientationString.intern();
	}

	public void setReversedOrientation(boolean isInReversedOrientation) {
		this.orientation = isInReversedOrientation ? "-" : "=";
	}

	public void setScore(double score) { this.score = score;}

	protected void setStart(String data) {
		this.start = Integer.parseInt(data);
	}

	public void setStart(int start) {
		this.start = start;
	}

	/**
	 * Checks whether two annotations differ by a small (fudge) factor
	 * @param fudge Maximum difference at either end or start to consider similar
	 * @return
	 */
	public boolean almostEqual(LightweightGenomicAnnotation other, int fudge) {
		return (Math.abs(getStart() - other.getStart() ) < fudge) &&  (Math.abs(getEnd() - other.getEnd()) < fudge);
	}

	public boolean contains(LightweightGenomicAnnotation other) {
		//System.out.println("Do this annotation " + toString() + " contains " + other );
		return (other.getChromosome() == null || other.getChromosome().equals(getChromosome())) && getStart() <= other.getStart() && getEnd() >= other.getEnd();
	}

	public boolean equals(Object o) {
		if(o == null || !this.getClass().equals(o.getClass())) {
			return false;
		}
		
		BasicLightweightAnnotation other = (BasicLightweightAnnotation) o;
		
		if(chromosome != null && !chromosome.equals(other.getChromosome())) {
			return false;
		}
		
		if(name != null && !name.equals(other.getName())) {
			return false;
		}
		
		return start == other.getStart() && end == other.getEnd(); 
	}

	public boolean inReversedOrientation() {
		return "-".equals(orientation);
	}

	/**
	 * 
	 * @param other GenomicAnnotation to check of overlap
	 * @return true if the current instance overlaps with the other one.
	 */
	public boolean overlaps(LightweightGenomicAnnotation other) {
		//System.out.print("Does this<"+getName()+"_"+getStart()+"-"+getEnd()+"> overlap <"+other.getName()+"_"+other.getStart()+"-"+other.getEnd()+">");
		//System.out.println(" !<"+(getStart() < other.getEnd() && getEnd() > other.getStart())+">!");
		return overlaps(other, 0);
	}

	/**
	 * 
	 * @param other - other genomic annotation
	 * @param buffer if the overlap is within buffer they will be considered overlapping 
	 *        even if they do not overlap within their original boundaries.
	 * @return true if they overlap in this extended definition
	 */
	public boolean overlaps(LightweightGenomicAnnotation other, int buffer) {
	
		return ((other.getChromosome() == null) || other.getChromosome().equals(getChromosome())) && (getStart() < other.getEnd() + buffer) && (getEnd() > other.getStart() - buffer);
	}

	public boolean overlaps(Collection<LightweightGenomicAnnotation> others, int buffer) {
		
		boolean rtrn = false;
		for(LightweightGenomicAnnotation o : others ){
			if(overlaps(o, buffer)) {
				rtrn = true;
				break;
			}
		}
		return rtrn;
		
	}

	/**
	 * Merge overlapping regions into a single region and leave singletons alone
	 * @param regions The regions to merge
	 * @return The set of merged regions
	 */
	public static Collection<BasicLightweightAnnotation> mergeAllOverlappers(TreeSet<BasicLightweightAnnotation> regions) {
		
		Collection<BasicLightweightAnnotation> rtrn = new TreeSet<BasicLightweightAnnotation>();
		
		// If provided set of regions contains zero or one element, return a copy of the set itself
		if(regions.size() == 0) return rtrn;
		if(regions.size() == 1) {
			rtrn.addAll(regions);
			return rtrn;
		}
		
		Iterator<BasicLightweightAnnotation> iter = regions.iterator();
		// Begin with the first interval in the sorted set
		BasicLightweightAnnotation firstInterval = iter.next();
		// A changing interval that expands when the next interval overlaps it
		// Reset when the next interval does not overlap
		BasicLightweightAnnotation growingInterval = new BasicLightweightAnnotation(firstInterval.getChromosome(), firstInterval.getStart(), firstInterval.getEnd());
		
		while(iter.hasNext()) {
			
			BasicLightweightAnnotation nextInterval = iter.next();
			
			if(nextInterval.overlaps(growingInterval)) {
				// Next interval overlaps growingInterval
				// Update end point of growingInterval
				// Do not need to change start point because set is sorted by start point first
				growingInterval.setEnd(Math.max(growingInterval.getEnd(),nextInterval.getEnd()));
				if(!iter.hasNext()) {
					// This is the last interval in the set
					// Add the last element and leave loop
					BasicLightweightAnnotation lastAdd = new BasicLightweightAnnotation(growingInterval.getChromosome(), growingInterval.getStart(), growingInterval.getEnd());
					rtrn.add(lastAdd);
					continue;					
				}
			} else {
				// Next interval does not overlap growingInterval
				// Add the latest version of growingInterval to the set
				BasicLightweightAnnotation toAdd = new BasicLightweightAnnotation(growingInterval.getChromosome(), growingInterval.getStart(), growingInterval.getEnd());
				rtrn.add(toAdd);
				if(!iter.hasNext()) {
					// This is the last interval in the set
					// Add it and leave loop
					BasicLightweightAnnotation lastAdd = new BasicLightweightAnnotation(nextInterval.getChromosome(), nextInterval.getStart(), nextInterval.getEnd());
					rtrn.add(lastAdd);
					continue;
				}
				// Reset growingInterval to this new interval
				growingInterval.setChromosome(nextInterval.getChromosome());
				growingInterval.setStart(nextInterval.getStart());
				growingInterval.setEnd(nextInterval.getEnd());
				continue;
			}
		}
		
		return rtrn;
	}

	
	public boolean overlaps(Collection<LightweightGenomicAnnotation> others) {
		return overlaps(others, 0);
	}

	public boolean isOrientation() {
		return "+".equals(orientation);
	}

	public void addExtraScore (double score) {
		extraScores.add(score);
	}

	public void removeExtraScores() {
		extraScores = new ArrayList<Double>();
	}

	public void setBoundariesFromAnnoations(List<? extends GenomicAnnotation> annotations) {
		Iterator<? extends GenomicAnnotation> it = annotations.iterator();
		while(it.hasNext()) {
			LightweightGenomicAnnotation annot = it.next();
			start = start > 0 ? Math.min(start, annot.getStart()) : annot.getStart();
			end   = Math.max(end, annot.getEnd());
		}
	}

	public void stitchTo(LightweightGenomicAnnotation other) {
		setStart(Math.min(getStart(), other.getStart()));
		setEnd(Math.max(getEnd(), other.getEnd()));
	}

	/**
	 * Change the current instance to represent its intersection
	 * with the provided annotation if they overlap
	 * @param other
	 */
	public void takeIntersection(LightweightGenomicAnnotation other) {
		if(getStart() < other.getEnd() && getEnd() > other.getStart()) {
			setStart(Math.max(getStart(),other.getStart()));
			setEnd(Math.min(getEnd(), other.getEnd()));
		}
	}

	/**
	 * Change the current instance to represent its union
	 * with the provided annotation 
	 * @param GenomicAnnotation other
	 */
	public void takeUnion(LightweightGenomicAnnotation other) {
		if(overlaps(other)) {
			setStart(Math.min(getStart(),other.getStart()));
			setEnd(Math.max(getEnd(), other.getEnd()));
		}
	}

	public static LightweightGenomicAnnotation createFromUCSC(String ucsc) {
		String [] firstSplit = ucsc.split(":");
		String [] secondSplit = firstSplit[1].split("-");
		
		return new BasicLightweightAnnotation(firstSplit[0],Integer.parseInt(secondSplit[0]), Integer.parseInt(secondSplit[1]));
	}

	/**
	 * Change the current instance to represent its intersection
	 * with the provided annotation if they overlap
	 * @param other
	 */
	public LightweightGenomicAnnotation intersect(LightweightGenomicAnnotation other) {
		LightweightGenomicAnnotation result = null;
		if(chromosome.equals(other.getChromosome()) && getStart() < other.getEnd() && getEnd() > other.getStart()) {
			result  = new BasicLightweightAnnotation(chromosome, Math.max(getStart(),other.getStart()), Math.min(getEnd(), other.getEnd()));
		}
		return result;
	}

	public List<LightweightGenomicAnnotation> intersect(List<? extends LightweightGenomicAnnotation> annotations) {
		List<LightweightGenomicAnnotation> result = new ArrayList<LightweightGenomicAnnotation>() ;
		for(LightweightGenomicAnnotation a : annotations) {
			LightweightGenomicAnnotation aint = intersect(a);
			if(aint != null) {
				result.add(aint);
			}
		}
		return stitchList(result, 0);
	}

	public int hashCode() {
		int hash = 0;
	
		if(chromosome != null) {
			hash = chromosome.hashCode();
		}
		
		if (name != null) {
			hash += name.hashCode();
		}
		
		hash += (start + end)/2;
		return hash;
	}

	public int length() { return getEnd() - getStart(); }

	public String toUCSC() {
		return getLocationString();
	}

	public String toString() {
		return toUCSC();
	}
	/**
	 * 
	 * @param sortedList
	 * @param maxDistanceToStitch
	 * @return
	 */	
	public static <T extends LightweightGenomicAnnotation> List<T> stitchList(Collection<? extends T> sortedList, int maxDistanceToStitch) {

		Stack<T> result = new Stack<T>();

		if(sortedList.size() == 0 ) {
			return result;
		}
		
		Iterator<? extends T> it = sortedList.iterator();
		result.push(it.next());
		while(it.hasNext()) {
			T next = it.next();
			T curr = result.pop();
			if(curr.overlaps(next,maxDistanceToStitch)) {
				curr.stitchTo(next);
				curr.setName(curr.getName()+"-"+next.getName());
				result.push(curr);
			} else {
				result.push(curr);
				result.push(next);
			}
		}
		
		return result;
	}
	
	public int compareTo(LightweightGenomicAnnotation b){
		;
		if(this.equals(b)){return 0;}
		if(!b.getChromosome().equalsIgnoreCase(getChromosome())){return b.getChromosome().compareTo(getChromosome());}
		if(b.getStart() == getStart()) {
			return getEnd() - b.getEnd();
		} else {
			return getStart() - b.getStart();
		}

		
	}
}