package broad.core.alignment;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Stack;

import broad.core.annotation.BED;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.sequence.SequenceRegion;
import broad.core.sequence.WindowSlider;

public class AlignmentList <T extends AlignmentSummary>{
	private List<T> alignmentList;
	private Map<String, List<T>> ambiguousTopHitMap;
	private AlignmentStats stats;
	private List<Short> coverageDistribution; 
	private short[] baseCoverage;

	protected AlignmentList() {
		super();
		alignmentList = new ArrayList<T>();
		stats = new AlignmentStats();
	}
	
	protected void addAlignment(T summary) {
		stats.add(summary);
		alignmentList.add(summary);
	}
	

	public java.util.List<T> getAlignmentList() {
		return alignmentList;
	}
	
	protected void addAlignmentSummary(T aln) {
		alignmentList.add(aln);
		stats.add(aln);
	}
	
	public AlignmentStats getStatistics() {
		return stats;
	}

	public void sortByQueryPosition() {
		Collections.sort(alignmentList, new Comparator<AlignmentSummary>() {

			public int compare(AlignmentSummary arg0, AlignmentSummary arg1) {
				return arg0.getA().getStart() - arg1.getA().getStart();
			}									
		});
	}

	public short[] getCoverageOfQuery() {
		Collections.sort(alignmentList, new Comparator<AlignmentSummary>() {

			public int compare(AlignmentSummary arg0, AlignmentSummary arg1) {
				return arg0.getA().getEnd() - arg1.getA().getEnd() == 0 
					? arg0.getA().getStart() - arg1.getA().getStart()
					: arg0.getA().getEnd() - arg1.getA().getEnd();
			}									
		});
		
		SequenceRegion alignedQuery = getQueryRegion();
		int last = alignedQuery.getEnd();
		//System.out.println("Aligned region " + alignedQuery);
		baseCoverage = new short[last + 1];
		coverageDistribution = new ArrayList<Short>();
		//System.out.println("baseCoverage size " + baseCoverage.length + " last was " + last + " alignmentList last alingment " + alignmentList.get(alignmentList.size() - 1).getA());
		Iterator<T> it = alignmentList.iterator();
		AlignmentSummary aln = null;
		int k = 1;
		while(it.hasNext()) {
			aln = it.next();

			for(int i = aln.getA().getStart() - 1; i < aln.getA().getEnd() ; i++) {
				baseCoverage[i]++;
			}
			k++;
		}
		
		// Now compute coverage distribution stats.
		for(int i = 0; i < last + 1; i ++) {
			if(baseCoverage[i] > 0) {
				coverageDistribution.add(baseCoverage[i]);
			}
		}
		
		Collections.sort(coverageDistribution);
		
		return baseCoverage;
	}
	
	public void filterSubjectBestHits() {
		if(alignmentList.isEmpty()) {
			return;
		}
		
		Collections.sort(alignmentList, new Comparator<T>() {

			public int compare(T arg0, T arg1) {
				int comparison = arg0.getSubject().compareTo(arg1.getSubject());
				if(comparison == 0) {
					comparison = Math.round(100 * (arg1.getScore() - arg0.getScore()));
				}
				return comparison;
			}
			
		}) ;
		
		Stack<T> filteredAlignments = new Stack<T>();
		ambiguousTopHitMap = new HashMap<String, List<T>>();
		
		Iterator<T> it = alignmentList.iterator();
		filteredAlignments.push(it.next());
		while(it.hasNext()) {
			T aln = it.next();
			
			if(ambiguousTopHitMap.get(aln.getSubject()) != null) {
				//Already detected as multihit, lets see if this is another undistinguishable top hit
				List<T> otherHits = ambiguousTopHitMap.get(aln.getSubject());
				if(otherHits.get(0).getScore() == aln.getScore()) {
					otherHits.add(aln);
				}
				
				continue;
			}
			
			T last = filteredAlignments.pop();						
			if(!aln.getSubject().equals(last.getSubject())) {
				filteredAlignments.push(last);
				filteredAlignments.push(aln);
			} else if(last.score == aln.score){
				List<T> multihit = new ArrayList<T>();
				ambiguousTopHitMap.put(aln.getSubject(), multihit);
				multihit.add(last);
				multihit.add(aln);
			} 
			
		}
		
		System.out.println("Filtered alignment # " + filteredAlignments.size());
		alignmentList = filteredAlignments;
	}
	
	public List<BED> getCoverageOfQueryInWindows(int windowSize, int windowOverlap) {
		Collections.sort(alignmentList, new Comparator<AlignmentSummary>() {

			public int compare(AlignmentSummary arg0, AlignmentSummary arg1) {
				return arg0.getA().getEnd() - arg1.getA().getEnd() == 0 
					? arg0.getA().getStart() - arg1.getA().getStart()
					: arg0.getA().getEnd() - arg1.getA().getEnd();
			}									
		});
		
		SequenceRegion alignedQuery = getQueryRegion();
		WindowSlider slider = alignedQuery.getSlider(windowSize, windowOverlap);
		List<BED> windows = new ArrayList<BED>();
		short[] coverage = getCoverageOfQuery();
		
		while(slider.hasNext()) {
			SequenceRegion windowReg = slider.next();
			BED window = new BED(windowReg);
			int windowCoverage = 0;
			for(int i = window.getStart(); i <= window.getEnd(); i++) {
				windowCoverage += coverage[i];
			}
			
			window.setScore(windowCoverage);
			windows.add(window);
		}

		return windows;
	}

	private SequenceRegion getQueryRegion() {
		LightweightGenomicAnnotation firstQuery = alignmentList.get(0).getA();
		LightweightGenomicAnnotation lastQuery = alignmentList.get(alignmentList.size() - 1).getA();
		//System.out.println("Alignmentlist size " + alignmentList.size() + " first query " + firstQuery + ", last " + lastQuery );
		/*Iterator<T> it = alignmentList.iterator();
		int max = 0;
		int min = 1000000000;
		while(it.hasNext()) {
			T a = it.next();
			max = a.getA().getEnd() > max ? a.getA().getEnd() : max;
			min = a.getA().getStart() < min ? a.getA().getStart() : min;
			System.out.println("a A " + a.getAEnd() + " max " + max);
		}
		*/
		SequenceRegion alignedQuery = new SequenceRegion(firstQuery.getName(), firstQuery);
		alignedQuery.setEnd(lastQuery.getEnd());
		//alignedQuery.setStart(min);
		//System.out.println("alignedQuery region " + alignedQuery);
		return alignedQuery;
	}
	


	public short[] getBaseCoverage() {
		return baseCoverage;
	}

	public List<Short> getCoverageDistribution() {
		return coverageDistribution;
	}

	public Map<String, List<T>> getAmbiguousTopHitMap() {
		return ambiguousTopHitMap;
	}
}
