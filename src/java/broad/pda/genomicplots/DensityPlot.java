package broad.pda.genomicplots;

import java.awt.Color;
import java.awt.Graphics2D;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import broad.pda.genomicplots.IntervalDensityReader.DensityInterval;


public class DensityPlot extends Plot {
	private static final long serialVersionUID = 8239324L;
	ArrayList<List<DensityInterval>> densities = new ArrayList<List<DensityInterval>>();
	ArrayList<Color> densityColors = new ArrayList<Color>();

	public DensityPlot(int xSize, int ySize, int maxX, int maxY, int minX,
			int minY, int xTickStop, int yTickStop) {
		super(xSize, ySize, maxX, maxY, minX, minY, xTickStop, yTickStop);
	}

	public DensityPlot(int xSize, int ySize, int maxX, int maxY, int minX,
			int minY, int xTickStop, int yTickStop, int yPosition) {
		super(xSize, ySize, maxX, maxY, minX, minY, xTickStop, yTickStop,
				yPosition);
	}
	
	public void addDensityIntervals(List<DensityInterval> intervals, Color color) {
		densities.add(intervals);
		densityColors.add(color);
	}
	
	public void plot(Graphics2D g) {

		Iterator<DensityInterval> densityIt = null;
		Color color = null;
		DensityInterval interval = null;
		PlotItem priorItem = null;
		for(int i = 0; i < densities.size(); i++) {
			densityIt = densities.get(i).iterator();
			color   = densityColors.get(i);
			g.setColor(color);
			priorItem = null;
			while(densityIt.hasNext()) {
				interval = densityIt.next();
				PlotItem intervalPlotItem = new DensityPlotItem(interval);
				intervalPlotItem.setColor(color);
				intervalPlotItem.plot((Graphics2D)g, this);
				if (priorItem != null) {
					g.drawLine(getPlotXCoordinate(priorItem.getAEnd()), 
						getPlotYCoordinate(priorItem.getBEnd()), 
						getPlotXCoordinate(intervalPlotItem.getAStart()), 
						getPlotYCoordinate(intervalPlotItem.getBStart()));
				}
				priorItem = intervalPlotItem;
			}
		}
	}

}
