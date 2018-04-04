package broad.pda.genomicplots;

import java.awt.Graphics2D;

import broad.core.alignment.AlignmentSummary;

public class RectanglePlotItem extends PlotItem {

	private static final long serialVersionUID = 19071123832L;

	public RectanglePlotItem() {
		super();
	}

	public RectanglePlotItem(AlignmentSummary summary) {
		super(summary);
	}
	
	public void plot(Graphics2D g, BaseToPlotCoordinateMapper plot, boolean fill) {
		g.setPaint(getColor());
		if (fill) {
			g.fillRect(plot.getPlotXCoordinate(getAStart()), plot.getPlotYCoordinate(getBStart()),
				plot.getPlotXDistance(getAEnd(),getAStart()), plot.getPlotYDistance(getBEnd(),getBStart()));
		} else {
			g.drawRect(plot.getPlotXCoordinate(getAStart()), plot.getPlotYCoordinate(getBStart()),
					plot.getPlotXDistance(getAEnd(),getAStart()), plot.getPlotYDistance(getBEnd(),getBStart()));
		}
		if(getCaption() != null) {
			//g.drawString(getCaption(),plot.getPlotXCoordinate(getAStart()), plot.getPlotYCoordinate(getBStart())+1);
		}
	}

}
