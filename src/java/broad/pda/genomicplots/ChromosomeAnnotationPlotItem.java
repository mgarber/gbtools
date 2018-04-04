package broad.pda.genomicplots;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Rectangle;

import broad.core.annotation.GenomicAnnotation;

public class ChromosomeAnnotationPlotItem extends UnorientedAnnotationPlotItem implements AnnotationPlotDrawable {
	public ChromosomeAnnotationPlotItem(GenomicAnnotation annot, Color color) {
		super(annot, color);
	}
	
	void drawOrientation(Graphics2D g, Color originalColor, Point orientationLineStart, Point orientationLineEnd) {
		g.setColor(Color.BLACK);
		g.drawLine(orientationLineStart.x, orientationLineStart.y, orientationLineEnd.x, orientationLineEnd.y);
		g.setColor(originalColor);
	}
	
	public void plot(Graphics2D g, GenomicAnnotationPlotter plot, boolean displayName) {
		super.plot(g, plot, displayName);
		Rectangle r = plot.createRectangle(annotation);
		Point orientationLineStart = r.getLocation();
		Point orientationLineEnd   = r.getLocation();
		if(annotation.inReversedOrientation()) {
			orientationLineStart.translate(r.width, 0); //top right corner
			orientationLineEnd.translate(0,r.height);   //bottom left corner
		} else {
			orientationLineEnd.translate(r.width,r.height);	//bottom right corner		
		}
		Color originalColor = g.getColor();		
		drawOrientation(g, originalColor, orientationLineStart, orientationLineEnd);
	}

}
