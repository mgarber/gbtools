package broad.pda.genomicplots;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;

import broad.core.annotation.GenomicAnnotation;

public class FloatingAnnotation extends AbstractGenomicAnnotationPlotItem {
	static final int SEPARATION = 8;
	GenomicAnnotation annotation;
	
	
	public FloatingAnnotation(GenomicAnnotation annot, Color color) {
		super();
		setColor(color);
		this.annotation = annot;
	}
	public void plot(Graphics2D g, GenomicAnnotationPlotter plot, boolean displayName) {
		int start = plot.getPlotYCoordinate(annotation.getStart());
		int end   = plot.getPlotYCoordinate(annotation.getEnd());
		int x     = plot.getPlotStartX() - SEPARATION;
		Color originalColor = g.getColor();		
		g.setColor(getColor());
		g.drawLine(x, start, x, end);
		System.out.println("Drawing line x" + x + " start " + start + " end " + end);
		//if(displayName ) {
		g.setColor(Color.BLACK);
		java.awt.Font originalFont = g.getFont();
		java.awt.Font textFont     = new java.awt.Font("arial", 0, 8).deriveFont(AffineTransform.getRotateInstance(- Math.PI/2));
		g.setFont(textFont);
		g.drawString(getName(),x - 5,end);
		g.setFont(originalFont);
		//}

		g.setColor(originalColor);
	}
	
	public String getName() {
		return annotation.getName();
	}

	public void plot(Graphics2D g, GenomicAnnotationPlotter plot) {
		plot(g, plot, true);
		
	}
	@Override
	public GenomicAnnotation getAnnotation() {
		return annotation;
	}

}
