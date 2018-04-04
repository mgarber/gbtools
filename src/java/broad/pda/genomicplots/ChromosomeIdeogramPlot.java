package broad.pda.genomicplots;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.geom.AffineTransform;
import java.util.ArrayList;
import java.util.Iterator;

import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;

public class ChromosomeIdeogramPlot implements GenomicAnnotationPlotter{

	Dimension dim;
	private Point topLeftCorner;
	BasicGenomicAnnotation ideogramSpan;
	private long plotTickStop;
	AffineTransform transform;
	
	private ArrayList<AbstractGenomicAnnotationPlotItem> children;
	private String title;

	public ChromosomeIdeogramPlot(String title, Point topLeftCornerPosition, int length, int width, long minPosition, long maxPosition, long tickStop) {
		super();
		children = new ArrayList<AbstractGenomicAnnotationPlotItem>();
		this.title = title;
		dim = new Dimension(width, length);
		this.topLeftCorner = topLeftCornerPosition;
		topLeftCorner.x += 50;
		topLeftCorner.y += 50;
		ideogramSpan = new BasicGenomicAnnotation("ideogramSpan");
		ideogramSpan.setStart((int)minPosition);
		ideogramSpan.setEnd((int)maxPosition);

		this.plotTickStop = Math.round(tickStop*length/(maxPosition - minPosition));
		System.out.println("Created ideogram " + title + " at " + topLeftCorner + " with hight " + dim.height);
	}
	
	protected long getMinPosition() { return ideogramSpan.getStart();}
	protected long getMaxPosition() { return ideogramSpan.getEnd();}
	protected Point getTopLeftCorner() { return topLeftCorner;}
	protected Dimension getDimension() { return dim;}
	
	public void translateBy(Point point) {
		topLeftCorner.translate((int)point.getX(),(int) point.getY());
	}
	
	public void addAnnotation(GenomicAnnotation annot, Color color) {
		addItem(new ChromosomeAnnotationPlotItem(annot, color));
	}
	
	public void addUnorientedAnnotation(GenomicAnnotation annot, Color color) {
		addItem(new UnorientedAnnotationPlotItem(annot, color));
	}
	
	public void addCentromere(GenomicAnnotation centromere) {
		addItem(new CentromereAnnotation(centromere));	}
	
	public void addItem(AbstractGenomicAnnotationPlotItem item) {
		children.add(item);
	}
	
	public void plotHorizontal() {
		transform = new AffineTransform();
		// calculate the center of the rectangle this should be
		// the center of rotation.
		double sx = (topLeftCorner.getX() + dim.width)/2;
		double sy = (topLeftCorner.getY() + dim.height)/2; 
		transform.rotate(Math.PI/4f, sx, sy);
		// Transform the four vertices of the rsulting rectangle and translate 
		// so that the vertex closest to the top left corner ends there.
		
	}
	
	public void paintChromosome(Graphics2D g, boolean displayAnnotationNames) {
		// Give a 3% margine to coordinate legends and Plot title
		//System.out.println("\tPlotting at " + topLeftCorner);
		Iterator<AbstractGenomicAnnotationPlotItem> it = children.iterator();
		while(it.hasNext()) {
			AbstractGenomicAnnotationPlotItem item = it.next();
			item.plot(g, this, displayAnnotationNames);
		}
		
		g.setPaint(Color.black);
		Rectangle r = new Rectangle(topLeftCorner, new Dimension(dim.width, getPlotEndY() - getPlotStartY()));		
		g.draw(r);
		if(title != null) {
			g.drawString(title, topLeftCorner.x, (int)  topLeftCorner.y - 25);
		}
		//System.out.println("Drawing rectangle: top left corner "+topLeftCorner+", dim "+dim);
		//Draw ticks
		/*
		float tickNum = Math.round((plotEnd - plotStart)/plotTickStop);		
		for (int i = 1; i <= tickNum; i++) {
			g.drawLine(Math.round(plotStart+i*plotTickStop),plotEndY,Math.round(plotStart+i*plotTickStop),plotEndY+10);
			for (int j = 0 ; j <= 10; j++) {
				g.drawLine(Math.round(plotStart+(i-1)*plotTickStop +j*plotTickStop/10),plotEndY,Math.round(plotStart+(i-1)*plotTickStop +j*plotTickStop/10),plotEndY+5);
			}
		}	
		*/	
	}

	public int getPlotXCoordinate(long sequenceCoord) {
		return 0;
		//return (int) Math.round((sequenceCoord - minPosition)*(getPlotEndX() - getPlotStartX())/(maxPosition - minPosition)) + getPlotStartX();
	}
	
	/* (non-Javadoc)
	 * @see edu.mit.broad.prodinfo.genomicplot.BaseToPlotCoordinateMapper#getPlotYCoordinate(long)
	 */
	public int getPlotYCoordinate(long sequenceCoord) {
		int plotYCoord = getPlotStartY();
		
		if(sequenceCoord > getMinPosition()) {
			plotYCoord += (int) Math.round((Math.min(getMaxPosition(), sequenceCoord) - getMinPosition())*(getPlotEndY() - getPlotStartY())/(getMaxPosition() - getMinPosition()));	
			//plotYCoord += (int) Math.round((Math.min(getMaxPosition(), sequenceCoord) - getMinPosition())*(dim.height)/(getMaxPosition() - getMinPosition()));	
		}
		//System.out.println("minposition " + getMinPosition() + " minY: " + getPlotStartY()+" maxY " + getPlotEndY() +" Y Coordinate "+sequenceCoord+" is translated to "+(getPlotStartY() + (int) Math.round((sequenceCoord - getMinPosition())*(getPlotEndY() - getPlotStartY())/(getMaxPosition() - getMinPosition()))));
		return plotYCoord; 
	}
	
	/* (non-Javadoc)
	 * @see edu.mit.broad.prodinfo.genomicplot.BaseToPlotCoordinateMapper#getPlotYDistance(long, long)
	 */
	public int getPlotYDistance(long start, long end) {		
		return Math.abs(getPlotYCoordinate(start) - getPlotYCoordinate(end)); 
	}
	
	/* (non-Javadoc)
	 * @see edu.mit.broad.prodinfo.genomicplot.BaseToPlotCoordinateMapper#getPlotXDistance(long, long)
	 */
	public int getPlotXDistance(long start, long end) {		
		return Math.abs(getPlotXCoordinate(start) - getPlotXCoordinate(end)); 
	}

	public int getPlotEndX() {
		return (int) (topLeftCorner.getX() + dim.getWidth());
	}

	public int getPlotEndY() {
		return (int) (topLeftCorner.getY() + dim.getHeight());
	}

	public int getPlotStartX() {
		// TODO Auto-generated method stub
		return (int) topLeftCorner.getX();
	}

	public int getPlotStartY() {
		// TODO Auto-generated method stub
		return (int) topLeftCorner.getY();
	}

	public Rectangle createRectangle(LightweightGenomicAnnotation a) {
		// TODO Eventually will need to apply affine transformations to these, for now just create vertical decorations.
		Point topLeftCorner = new Point(getPlotStartX(), getPlotYCoordinate(a.getStart()));
		Dimension rectDim       = new Dimension((int) dim.getWidth(), getPlotYDistance(a.getStart(), a.getEnd()));
		return new Rectangle(topLeftCorner, rectDim);
	}

	public boolean contains(LightweightGenomicAnnotation b) {
		return ideogramSpan.contains(b);
	}
	
	public boolean overlaps(LightweightGenomicAnnotation b) {
		return ideogramSpan.overlaps(b);
	}
}
