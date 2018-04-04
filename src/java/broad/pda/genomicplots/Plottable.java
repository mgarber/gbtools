package broad.pda.genomicplots;

import java.awt.Graphics2D;

public interface Plottable extends Drawable {

	public void plot(Graphics2D g, BaseToPlotCoordinateMapper plot);

}