package broad.pda.assembly;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import broad.core.sequence.SequenceRegion;


public class ReadCoverageReader {
	ArrayList<CoveredSequenceRegion> coverageRegions;
	ArrayList<CoveredSequenceRegion> coverageRegionsOrderedByCoverage;
	
	public ReadCoverageReader(String sequenceId, File file) throws IOException {
		super();
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line = null;
		int lastCov = 0;
		coverageRegions = new ArrayList<CoveredSequenceRegion>();

		while((line = br.readLine()) != null) {
			String [] info = line.split("\t");
			int pos = Integer.parseInt(info[0]);
			int cov = Integer.parseInt(info[1]);
			if(pos == 1) {
				CoveredSequenceRegion first = new CoveredSequenceRegion(sequenceId);
				first.setRegionStart(1);
				first.setCoverage(cov);
				coverageRegions.add(first);
			}else if(cov != lastCov) {
				CoveredSequenceRegion priorRegion = coverageRegions.get(coverageRegions.size() - 1);
				priorRegion.setRegionEnd(pos - 1);
				CoveredSequenceRegion newRegion   = new CoveredSequenceRegion(sequenceId);
				newRegion.setRegionStart(pos);
				newRegion.setCoverage(cov);
				coverageRegions.add(newRegion);
			}
			lastCov = cov;
		}
	}	
	
	/**
	 * This gets the coverage, c such that percentOfBasesCovered covered intervals are covered by at least c reads. 
	 * @param percentOfintervals the percent of coverage intervals to include (for 0.5 the returned value is the Median)
	 * @return
	 */
	public int getCoverageForPercentOfBases(float percentOfintervals) {
		if (coverageRegionsOrderedByCoverage == null) {
			coverageRegionsOrderedByCoverage = new ArrayList<CoveredSequenceRegion>(coverageRegions.size());
			for(int i = 0; i < coverageRegions.size(); i++) {
				coverageRegionsOrderedByCoverage.add(coverageRegions.get(i));
			}
			Collections.sort(coverageRegionsOrderedByCoverage, new Comparator<CoveredSequenceRegion>() {

				public int compare(CoveredSequenceRegion arg0, CoveredSequenceRegion arg1) {
					return arg0.coverage - arg1.coverage;
				}
				
			});
		}
		
		System.out.print("#regions <"+coverageRegions.size()+"> percentOfIntervals<"+percentOfintervals+"> multiplication <" +(coverageRegions.size() * percentOfintervals) +"> ");
		int numberOfIntervals = (int) Math.floor((float)coverageRegions.size() * percentOfintervals);
		System.out.println("Number of intervals to use " + numberOfIntervals + " total size " + coverageRegions.size());
		return coverageRegionsOrderedByCoverage.get(coverageRegions.size() - numberOfIntervals).coverage;

	}
	
	public List<CoveredSequenceRegion>  getRegions(int minCoverage) {
		ArrayList<CoveredSequenceRegion> regions = new ArrayList<CoveredSequenceRegion>();
		Iterator<CoveredSequenceRegion> it = coverageRegions.iterator();
		while(it.hasNext()) {
			CoveredSequenceRegion reg = it.next();
			if(reg.getCoverage() >= minCoverage) {
				CoveredSequenceRegion priorReg = regions.size() > 0 
					? regions.get(regions.size() -1 )
					: null;
				if(priorReg != null && (priorReg.getEnd() + 1 == reg.getStart())) {
					//System.out.print("Stitching " + reg + " to " + priorReg);
					priorReg.stitchTo(reg);
					//System.out.println(" obtained " + priorReg);
				} else {
					//System.out.println("Adding " + reg);
					regions.add(reg);
				}
			}
		}
		
		return regions;
	}
	
	public static class CoveredSequenceRegion extends SequenceRegion {
		private int coverage;
		
		public CoveredSequenceRegion(String containingSequenceId) {
			super(containingSequenceId);
		}
		
		public void setCoverage(int coverage) { this.coverage = coverage;}
		public int getCoverage() {return coverage;}

	}
}
