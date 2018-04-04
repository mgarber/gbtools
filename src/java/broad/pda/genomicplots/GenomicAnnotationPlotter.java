package broad.pda.genomicplots;

import java.awt.Rectangle;

import broad.core.annotation.LightweightGenomicAnnotation;

public interface GenomicAnnotationPlotter  extends BaseToPlotCoordinateMapper{
	Rectangle createRectangle(LightweightGenomicAnnotation a);
}
