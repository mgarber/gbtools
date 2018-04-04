package broad.pda.genomicplots;

import java.awt.Color;
import java.awt.Graphics2D;
import java.text.DecimalFormat;

import broad.core.annotation.GenomicAnnotation;

public class FloatingPositionalAnnotation extends FloatingAnnotation {
	private static final DecimalFormat numFormat = new DecimalFormat("##0.0");
	private static final int M = 1000000;
	
	public FloatingPositionalAnnotation(GenomicAnnotation annot, Color color) {
		super(annot, color);
		
	}

	public void plot(Graphics2D g, GenomicAnnotationPlotter plot, boolean displayName)  {
		super.plot(g, plot, true);
	}
	
	public String getName() {
		return numFormat.format(annotation.getChromosome() + ":"+annotation.getStart()/(double)M) + "-" + numFormat.format(annotation.getEnd()/(double)M);
	}
}
