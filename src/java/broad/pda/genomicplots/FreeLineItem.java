package broad.pda.genomicplots;

import java.awt.Color;
import java.awt.Graphics2D;

import broad.core.annotation.LightweightGenomicAnnotation;

public class FreeLineItem implements Plottable {
	private LightweightGenomicAnnotation annot;
	private Color color;
	
	public FreeLineItem(LightweightGenomicAnnotation annotation) {
		this.annot = annotation;
	}

	public void plot(Graphics2D g, BaseToPlotCoordinateMapper plot) {
		//Do nothing until we rethink this a bit better....
	}

	public void setCaption(String caption) {
		//Does nothing, the caption is the annotation's name.
	}

	public void setColor(Color color) {
		this.color = color;
	}
	
	public Color getColor() {
		return color;
	}
	
	public long getStart() {
		return annot.getStart();
	}
	
	public long getEnd() {
		return annot.getEnd();
	}

}
