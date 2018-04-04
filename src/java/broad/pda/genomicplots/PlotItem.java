package broad.pda.genomicplots;

import java.awt.Color;
import java.awt.Graphics2D;
import java.io.Serializable;

import broad.core.annotation.TwoSubjectAnnotation;

public class PlotItem extends AbstractPlotItem implements Serializable, Plottable {

	TwoSubjectAnnotation annotation;

	private static final long serialVersionUID = 1000451L;

	public PlotItem() {
		super();
		// TODO Auto-generated constructor stub
	}
	
	public PlotItem(TwoSubjectAnnotation annotation){
		this.annotation = annotation;		
	}
	
	public PlotItem(TwoSubjectAnnotation annotation, Color color){
		this.annotation = annotation;		
		setColor(color);
	}

	/* (non-Javadoc)
	 * @see edu.mit.broad.prodinfo.genomicplot.Plottable#plot(java.awt.Graphics2D, edu.mit.broad.prodinfo.genomicplot.BaseToPlotCoordinateMapper)
	 */
	public void plot(Graphics2D g, BaseToPlotCoordinateMapper plot) {
		Color original = g.getColor();
		g.setPaint(getColor());
		//System.out.println(annotation);
		int startx = plot.getPlotXCoordinate(getAStart());
		int starty = plot.getPlotYCoordinate(annotation.isDirect() ? getBStart() : getBEnd());
		int endx   = plot.getPlotXCoordinate(getAEnd());
		int endy   = plot.getPlotYCoordinate(annotation.isDirect() ? getBEnd() : getBStart());
		if((endx == startx) && (endy == starty)) {
			endx++;
			endy++;
		}
		
		//g.drawLine(startx, starty, endx, endy);
		
		g.fillRect(Math.min(startx, endx), Math.min(starty,endy), Math.abs(endx-startx), Math.abs(endy-starty));
		
		//System.err.println("fill rectangle");

		/*
		System.out.println("plotted line from ("+getAStart()+","+getBStart()+") ->("+
				plot.getPlotXCoordinate(getAStart())+","+plot.getPlotYCoordinate(getBStart())+") to " +
				"("+getAEnd()+","+getBEnd()+ ") -> ("+
				plot.getPlotXCoordinate(getAEnd())+","+ plot.getPlotYCoordinate(getBEnd())+")");
		*/	
		g.setColor(original);
	}

	public int getAEnd() {
		return annotation.getAEnd();
	}

	public String getAName() {
		return annotation.getAName();
	}

	public int getAStart() {
		return annotation.getAStart();
	}

	public int getBEnd() {
		return annotation.getBEnd();
	}

	public String getBName() {
		return annotation.getBName();
	}

	public int getBStart() {
		return annotation.getBStart();
	}

	public void plot(Graphics2D g, GenomicAnnotationPlotter plot, boolean displayAnnotationNames) {
		plot(g, plot); 
		//TODO: display name
	}

}
