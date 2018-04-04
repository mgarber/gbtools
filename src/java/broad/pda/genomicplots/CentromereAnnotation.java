package broad.pda.genomicplots;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Polygon;
import java.awt.Rectangle;

import broad.core.annotation.GenomicAnnotation;


public class CentromereAnnotation extends UnorientedAnnotationPlotItem {

	public CentromereAnnotation(GenomicAnnotation annot, Color color) {
		super(annot, color);
	}
	
	public CentromereAnnotation(GenomicAnnotation annot) {
		this(annot, Color.RED);
	}


	public void plot(Graphics2D g, GenomicAnnotationPlotter plot, boolean displayName) {
		Rectangle r = plot.createRectangle(annotation);
		Polygon p = new Polygon();
		Point upperLeftCorner = r.getLocation();
		p.addPoint(upperLeftCorner.x, upperLeftCorner.y);
		p.addPoint(upperLeftCorner.x + r.width, upperLeftCorner.y);
		p.addPoint(upperLeftCorner.x, upperLeftCorner.y + r.height);
		p.addPoint(upperLeftCorner.x + r.width, upperLeftCorner.y + r.height);
		Color originalColor = g.getColor();
		g.setColor(Color.WHITE);
		g.fill(r);
		g.setColor(getColor());
		g.fill(p);
		g.setColor(originalColor);
	}
}
