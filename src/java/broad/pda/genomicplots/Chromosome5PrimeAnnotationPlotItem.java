package broad.pda.genomicplots;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Rectangle;

import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.GenomicAnnotation;

public class Chromosome5PrimeAnnotationPlotItem extends
		ChromosomeAnnotationPlotItem {

	public Chromosome5PrimeAnnotationPlotItem(GenomicAnnotation annot, Color color) {
		super(annot, color);
	}

	public void plot(Graphics2D g, GenomicAnnotationPlotter plot, boolean displayName) {
		super.plot(g, plot, displayName);
		BasicGenomicAnnotation telomereCap = new BasicGenomicAnnotation(annotation.getName() + "_telomere");
		telomereCap.setStart(annotation.inReversedOrientation() 
				? annotation.getEnd() - Math.round(annotation.getLength()/3) 
				: annotation.getStart());
		
		telomereCap.setEnd(annotation.inReversedOrientation() 
				? annotation.getEnd() 
				: annotation.getStart() + Math.round(annotation.getLength()/3));
		
		Rectangle r = plot.createRectangle(telomereCap);
		if(r.height > MAX_SUB_ANNOT_HEIGHT) {
			if(annotation.inReversedOrientation()) {
				r.translate(0, r.height - MAX_SUB_ANNOT_HEIGHT);
			}
			r.height = MAX_SUB_ANNOT_HEIGHT;
		}
		g.fill(r);

	}
}
