package broad.pda.genomicplots;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Rectangle;

import broad.core.annotation.GenomicAnnotation;


public class UnorientedAnnotationPlotItem extends AbstractGenomicAnnotationPlotItem {

	protected GenomicAnnotation annotation;
	protected int MAX_SUB_ANNOT_HEIGHT = 4;
	protected int MIN_ANNOT_HIGHT_FOR_CAP = 6;

	public UnorientedAnnotationPlotItem(GenomicAnnotation annot, Color color) {
		super();
		setColor(color);
		this.annotation = annot;
	}

	public int getStart() {
		return annotation.getStart();
	}

	public int getEnd() {
		return annotation.getEnd();
	}

	public void plot(Graphics2D g, GenomicAnnotationPlotter plot, boolean displayName) {
		//System.out.print("Ploting " + annotation.getName());
		Rectangle r = plot.createRectangle(annotation);
		//System.out.println(" as a Rectangle at " + r.getLocation());
		Color originalColor = g.getColor();		
		g.setColor(getColor());
		g.fill(r);
		if(displayName ) {
			Point startOfName = new Point((int)(r.getCenterX() - Math.floor(r.getWidth()/2)), (int)r.getCenterY());
			g.setColor(Color.white);
			java.awt.Font originalFont = g.getFont();
			g.setFont(new java.awt.Font("arial", 0, 2));
			g.drawString(annotation.getName(),(int) startOfName.getX(),(int) startOfName.getY());
			g.setFont(originalFont);
		}

		g.setColor(originalColor);
	}
	

	public void plot(Graphics2D g, GenomicAnnotationPlotter plot) {
		plot(g, plot, false);
		
	}

	@Override
	public GenomicAnnotation getAnnotation() {
		return annotation;
	}

}