package broad.pda.genomicplots;

import java.awt.Color;
import java.awt.Graphics2D;

import broad.core.annotation.GenomicAnnotation;

public abstract class AbstractGenomicAnnotationPlotItem extends AbstractPlotItem{
	Color color;
	
	public AbstractGenomicAnnotationPlotItem() {
		super();
	}

	public abstract void plot(Graphics2D g, GenomicAnnotationPlotter plot, boolean displayAnnotationNames);

	public abstract GenomicAnnotation getAnnotation();
	
	public Color getColor() {return color;}
	public void setColor(Color color) {this.color = color;}

}