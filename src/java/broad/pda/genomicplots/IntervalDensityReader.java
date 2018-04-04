package broad.pda.genomicplots;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import broad.core.annotation.BasicGenomicAnnotation;

public class IntervalDensityReader {
	private List<DensityInterval> densityIntervals;
	private int maxCoverage;
	private float maxDensity;
	
	public IntervalDensityReader(String fileName) throws IOException {
		maxCoverage = 0;
		maxDensity  = 0;
		densityIntervals = new ArrayList<DensityInterval>();
		File source = new File(fileName);
		BufferedReader br = new BufferedReader(new FileReader(source));
		String line;
		try {
			while((line = br.readLine()) != null) {
				//System.out.println(line);
				if(line.startsWith("#")){
					continue;
				}
				String [] lineArr = line.split("\t");
				DensityInterval interval = new DensityInterval(lineArr);
				maxCoverage = Math.max(maxCoverage, interval.getCoveredBases());
				maxDensity  = Math.max(maxDensity, interval.getDensity());
				
				densityIntervals.add(interval);
			}
		} finally {
			try {
				System.out.print("Closing "+fileName);
				br.close();
				System.out.print(" ..... Closed\n");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	public int getMaximumCoverage() { return maxCoverage; }
	public float getMaximumDensity() { return maxDensity; }
	
	public List<DensityInterval> getDensityIntervals() {
		return densityIntervals;
	}

	
	public static class DensityInterval extends BasicGenomicAnnotation {
		private int  coveredBases;
		
		public DensityInterval(int start, int end) {
			setStart(start);
			setEnd(end);
		}
		
		public DensityInterval(String [] rawData) {
			setStart(Integer.parseInt(rawData[0]));
			setEnd(Integer.parseInt(rawData[1]));
			coveredBases = Integer.parseInt(rawData[2]);
		}
		

		public int  getCoveredBases() { return coveredBases;}
		public void setCoveredBases(int bases) {this.coveredBases = bases;}
		public float getDensity() {return coveredBases/getLength();}
		public boolean inReversedOrientation() {return false;}
		public String getName() {return getStart()+"-" + getEnd();}
	}

}
