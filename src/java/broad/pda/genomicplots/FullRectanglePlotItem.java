package broad.pda.genomicplots;

import java.awt.Color;
import java.awt.Graphics2D;

import broad.core.annotation.LightweightGenomicAnnotation;

public class FullRectanglePlotItem {

	private static final long serialVersionUID = 10001221L;
	private boolean vertical;
	private LightweightGenomicAnnotation annotation;
	private Color color;
	
	public FullRectanglePlotItem(LightweightGenomicAnnotation annotation, boolean vertical){
		this.annotation = annotation;
		this.vertical = vertical; 
	}
	
	public void plot(Graphics2D g, BaseToPlotCoordinateMapper plot, boolean fill) {
		g.setPaint(getColor());
		int x = vertical ? plot.getPlotXCoordinate(annotation.getStart()) : plot.getPlotStartX();
		int y = vertical ? plot.getPlotStartY() : plot.getPlotYCoordinate(annotation.getEnd());

		if (x == plot.getPlotStartX()) {
			x = x + 1;
		}
		if (x == plot.getPlotEndX()) {
			x = x - 1;
		}
		if (y == plot.getPlotStartY()) {
			y = y + 1;
		}
		if (y == plot.getPlotEndY()) {
			y = y - 1;
		}
		int xLength = vertical ?
				plot.getPlotXDistance(annotation.getStart(),annotation.getEnd()) : 
				plot.getPlotEndX()-plot.getPlotStartX() - 2;
		
		int yLength = vertical ? 
				plot.getPlotEndY()-plot.getPlotStartY() - 2 :
					plot.getPlotYDistance(annotation.getStart(),annotation.getEnd()) ;
				
		if(fill) {
			g.fillRect(x, y,xLength, yLength);
		}else {
			g.drawRect(x, y,xLength, yLength);
		}
	}
	
	public void setColor(Color color) {
		this.color = color;
	}
	
	public Color getColor() {
		return color;
	}

}
