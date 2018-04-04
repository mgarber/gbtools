package broad.pda.genomicplots;

import java.awt.Graphics2D;

public interface AnnotationPlotDrawable extends Drawable {
	public void plot(Graphics2D g, GenomicAnnotationPlotter plot);
	public void plot(Graphics2D g, GenomicAnnotationPlotter plot, boolean displayAnnotationNames);

}
