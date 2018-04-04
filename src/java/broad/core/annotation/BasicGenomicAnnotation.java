package broad.core.annotation;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import broad.core.sequence.Sequence;

public class BasicGenomicAnnotation extends BasicLightweightAnnotation implements GenomicAnnotation {
	private Sequence sequence;
	private int fivePrimerBuffer;
	private int threePrimerBuffer;
	private String id;
	//This is a hack to speed up sequencial search and operations that involve a list
	int lastBreakInIndex = 0;
	
	protected BasicGenomicAnnotation() {
		super();
	}
	public BasicGenomicAnnotation(String name) {
		super();
		setName(name);
		//this.displayName = name;
	}
	
	public BasicGenomicAnnotation(String name, String chr, int start, int end) {
		this(name);
		setChromosome(chr);
		setStart(start);
		setEnd(end);
	}
	
	public BasicGenomicAnnotation(LightweightGenomicAnnotation ga) {
		super(ga);
	}
	
	public BasicGenomicAnnotation(GenomicAnnotation ga) {
		super(ga);
		setSequence(ga.getSequence());
		fivePrimerBuffer = ga.getFivePrimeBases();
		threePrimerBuffer = ga.getThreePrimeBases();
		id = ga.getId();
	}
	
	/**
	 * A basic genomic annotation may be built from a raw string array
	 * where the first entry is the chromosome, the second the start and
	 * the third the end.
	 */
	public BasicGenomicAnnotation(String [] info) {
		super();
		//? info[0] : info[0] + "_" + info[1] + "-" + info[2])
		if(info[0].contains(":") ) {
			String [] firstSplit = info[0].split(":");
			setChromosome(firstSplit[0]);
			String [] secondSplit = firstSplit[1].split("-");
			setStart(Integer.parseInt(secondSplit[0]));
			setEnd(Integer.parseInt(secondSplit[1]));
			if(info.length > 1) {
				setScore(Double.parseDouble(info[1]));
			}
		} else {
			setChromosome(info[0].substring(3));
			setStart(Integer.parseInt(info[1]));
			setEnd(Integer.parseInt(info[2]));
			if(info.length > 3) {
				setId(info[3]);
				setName(info[3]);
			} else {
				setName(info[0] + "_" + info[1] + "-" + info[2]);
			}
		}
	}

	public int getLength() { return length(); }
	public String toString() {
		return getLocationString() + "\t" + getScore();
	}


	
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id.intern();
	}
	
	
	/**
	 * Change the current instance to represent the difference
	 * of the instance with the other annotation
	 * @param other
	 */
	public void reduceToDifference(GenomicAnnotation other) {
		if(overlaps(other)) {
			//System.out.println("Reduced annotation! " + this);
			if(getStart() < other.getStart()) {
				setEnd(Math.min(getEnd(), other.getStart() - 1));
			} else {
				setStart(Math.max(getStart(), other.getEnd() + 1));
			}
		}
	}
	
	/**
	 * Fragments the current annotation if it overlaps the provided one.
	 * it returns an empty list if this annotation does not ovelap the
	 * one passed in.
	 * 
	 * @param GenomicAnnotation annotation to disect the current one
	 * @return List<GenomicAnnotation> of the disected annotations, an empty list is returned if 
	 * 			the annotations do not overlap.
	 */
	public List<GenomicAnnotation> disect(GenomicAnnotation a) {
		List<GenomicAnnotation> disection = new ArrayList<GenomicAnnotation>();
		if(overlaps(a)) {
			disection.addAll(this.minus(a));
			//System.out.println("minus resulted in " + disection);
			GenomicAnnotation copy = new BasicGenomicAnnotation(this);
			copy.takeIntersection(a);
			//System.out.println("and intersection in " + copy);
			disection.add(copy);
			
			Collections.sort(disection);
		}
		//System.out.println("and disction is " + disection);
		return disection;
	}
	
	
	/**
	 * Fragments the current annotation if it overlaps the provided ones.
	 * It returns a list with one component (this annotation) if no annotation 
	 * in the provided list overlaps the discted annotaion.
	 * 
	 * @param List<GenomicAnnotation> <b>sorted</b> annotations with which to disect the current one.
	 * @return List<GenomicAnnotation> of the disected annotations, a list just this annotation is returned if 
	 * 			the annotations do not overlap.
	 */
	public List<GenomicAnnotation> disect(List<? extends GenomicAnnotation> disectors) {
		List<GenomicAnnotation> disection = new ArrayList<GenomicAnnotation>();
		disection.add(this);
		Iterator<? extends GenomicAnnotation> annotIt = disectors.iterator();
		while(annotIt.hasNext()) {
			GenomicAnnotation a= annotIt.next();
			if(a.getStart() > getEnd()) {
				break;
			}else if(!overlaps(a)) {
				continue;
			}
			List<GenomicAnnotation> newDisection = new ArrayList<GenomicAnnotation>(disection.size());
			Iterator<GenomicAnnotation> oldDisectionIt = disection.iterator();
			while(oldDisectionIt.hasNext()) {
				GenomicAnnotation disectionComponent = oldDisectionIt.next();
				if(disectionComponent.overlaps(a)) {
					newDisection.addAll(disectionComponent.disect(a));
				} else {
					newDisection.add(disectionComponent);
				}
			} 
			disection = newDisection;
		}
		return disection;
	}
	
	
	public boolean isFlankedBy(TwoSubjectAnnotation twoSubjectAnnotation, int buffer) {
		GenomicAnnotation left = new BasicGenomicAnnotation("left");
		left.setStart(getStart() - buffer);
		left.setEnd(getStart() + buffer);
		
		GenomicAnnotation right = new BasicGenomicAnnotation("right");
		right.setStart(getEnd() - buffer);
		right.setEnd(getEnd() + buffer);
		
		
		GenomicAnnotation A = new BasicGenomicAnnotation("A");
		A.setStart(twoSubjectAnnotation.getAStart());
		A.setEnd(twoSubjectAnnotation.getAEnd());
		
		GenomicAnnotation B = new BasicGenomicAnnotation("B");
		B.setStart(twoSubjectAnnotation.getBStart());
		B.setEnd(twoSubjectAnnotation.getBEnd());
		
		return (A.overlaps(left) && B.overlaps(right)) || (A.overlaps(right) && B.overlaps(left));
	}
	
	public List<GenomicAnnotation> minus(List<? extends GenomicAnnotation> others, boolean listIsDisjoint, int startAtIdx) {
		ArrayList<GenomicAnnotation> dif = new ArrayList<GenomicAnnotation>();
		dif.add(this);
		Collections.sort(others);

		Iterator<? extends GenomicAnnotation> it = null;
		if(!listIsDisjoint) {
			List<GenomicAnnotation> stitchedOthers = stitchList(others, 0);
			it = stitchedOthers.iterator();
		} else {
			it = others.listIterator(startAtIdx);
		}

		lastBreakInIndex = startAtIdx;
		while(it.hasNext()) {
			GenomicAnnotation a = it.next();
			if(a.getStart() > getEnd()) {
				break;
			}
			lastBreakInIndex++;
			ArrayList<GenomicAnnotation> newList = new ArrayList<GenomicAnnotation>(dif.size());
			Iterator<GenomicAnnotation> oldListIt = dif.iterator();
			while(oldListIt.hasNext()) {
				GenomicAnnotation oldRegion = oldListIt.next();
				newList.addAll(oldRegion.minus(a));
			}
			dif = newList;
		}
		
		//Collections.sort(dif, new PositionComparator());
		
		return dif;
	}
	
	public List<GenomicAnnotation> minus(List<? extends GenomicAnnotation> others) {
		return minus(others, false, 0);
	}
	
	public List<GenomicAnnotation> minus(GenomicAnnotation other) {
		List<GenomicAnnotation> dif = new ArrayList<GenomicAnnotation>();
		if(!overlaps(other)) {
			dif.add(this);
		} else {
			GenomicAnnotation intersection = new BasicGenomicAnnotation(this);
			intersection.takeIntersection(other);
			GenomicAnnotation first = new BasicGenomicAnnotation(this);
			first.setEnd(intersection.getStart() - 1);
			if(first.length() > 0) {
				dif.add(first);
			}
			GenomicAnnotation second = new BasicGenomicAnnotation(this);
			second.setStart(intersection.getEnd() + 1);
			if(second.length() > 0) {
				dif.add(second);
			}
		}
		return dif;
	}
	
	public int lastIdxWhereListSearchBroke() { return lastBreakInIndex;}
	

	public Sequence getSequence() {
		return sequence;
	}
	public void setSequence(Sequence sequence) {
		this.sequence = sequence;
	}
	/*
	public String getDisplayName() {
		return displayName;
	}
	public void setDisplayName(String displayName) {
		this.displayName = displayName;
	}
	*/

	public void setOrientation(boolean orientation) {
		setOrientation(orientation ? "+" : "-");
	}
	public void setThreePrimeBuffer(int bufferSize) {
		setEnd(getEnd() + bufferSize);
		this.threePrimerBuffer = bufferSize;
	}
	
	public void setFivePrimeBuffer(int bufferSize) {
		setStart(getStart() - bufferSize);
		this.fivePrimerBuffer = bufferSize;
	}
	
	public int getFivePrimeBases() {
		return fivePrimerBuffer;
	}


	public int getThreePrimeBases() {
		return threePrimerBuffer;
	}
	
	/** 
	 * Default compareTo method: If annotations are in the same chromosome
	 * their start/end locations are compared otherwise their chromosomes are.
	 */
	public int compareTo(GenomicAnnotation arg0) {
		if(getChromosome() != null && arg0.getChromosome() != null && !getChromosome().equals(arg0.getChromosome())) {
			return getChromosome().compareTo(arg0.getChromosome());
		}
		return getStart() != arg0.getStart() ? getStart() - arg0.getStart() : getEnd() - arg0.getEnd();
	}
	public int getOrientedEnd() {
		return inReversedOrientation() ? getStart() : getEnd();
	}
	public int getOrientedStart() {
		// TODO Auto-generated method stub
		return inReversedOrientation() ? getEnd() : getStart();
	}
	public void addBlock(String name, int start, int end) {
		//Do nothing if does not support blocks, override if desired.
	}
	public List<? extends GenomicAnnotation> getBlocks() {
		return null;
	}
	public boolean mayHaveBlocks() {
		return false;
	}

}
