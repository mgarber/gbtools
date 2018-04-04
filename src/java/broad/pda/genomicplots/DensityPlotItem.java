package broad.pda.genomicplots;

import java.awt.Graphics2D;

import broad.pda.genomicplots.IntervalDensityReader.DensityInterval;


public class DensityPlotItem extends PlotItem {

	private static final long serialVersionUID = 10003242L;
	DensityInterval interval;

	public DensityPlotItem(DensityInterval annotation) {
		this.interval = annotation;
	}
	
	public void plot(Graphics2D g, BaseToPlotCoordinateMapper plot) {
		if(interval.getCoveredBases() > 0) {
			super.plot(g, plot);
		}
	}
	
	public int getAEnd() {
		return interval.getEnd();
	}

	public String getAName() {
		return interval.getName();
	}

	public int getAStart() {
		return interval.getStart();
	}

	public int getBEnd() {
		return interval.getCoveredBases();
	}

	public String getBName() {
		return interval.getName();
	}

	public int getBStart() {
		return interval.getCoveredBases();
	}

}
