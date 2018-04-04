package broad.pda.genomicplots;

import java.awt.Color;
import java.awt.Graphics2D;
import java.io.Serializable;

import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
;

public class Plot implements Serializable, BaseToPlotCoordinateMapper {
	private static final long serialVersionUID = 1000432L;
	public static final long DEFAULT_SCALE = 1000;
	public static final int AXIS_SCALE_SIZE = 350;
	public static final int MARGIN = 100;
	
	int xSize, ySize;
	LightweightGenomicAnnotation xRange;
	LightweightGenomicAnnotation yRange;

	String title;
	int xTickStop, yTickStop;
	// Calculated Coordinates
	int plotStartX;
	int plotStartY;
	int plotEndY;
	int plotEndX;

	long xScale = DEFAULT_SCALE;
	long yScale = DEFAULT_SCALE;
	
	public Plot(int xSize, int ySize, int maxX, int maxY, int minX, int minY, int xTickStop, int yTickStop) {
		super();
		init(xSize, ySize, maxX, maxY, minX, minY, xTickStop, yTickStop, 0, 5);
	}
	
	public Plot(int xSize, int ySize, int maxX, int maxY, int minX, int minY, int xTickStop, int yTickStop, int yPosition) {
		super();
		init(xSize, ySize, maxX, maxY, minX, minY, xTickStop, yTickStop, yPosition, 5);
	}
	
	public Plot (int xSize, int ySize, int maxX, int maxY, int minX, int minY, int xTickStop, int yTickStop, int yPosition, int xPosition) {
		init(xSize, ySize, maxX, maxY, minX, minY, xTickStop, yTickStop, yPosition, xPosition);
	}
	
	private void init(int xSize, int ySize, int maxX, int maxY, int minX, int minY, int xTickStop, int yTickStop, int yPosition, int xPosition) {
		int xMargin = MARGIN;//(int) Math.round(xSize * 0.03);
		int yMargin = MARGIN;//(int) Math.round(ySize * 0.03);
		
		xRange = new BasicGenomicAnnotation("xRange");
		xRange.setStart(minX);
		xRange.setEnd(maxX);
		
		yRange = new BasicGenomicAnnotation("yRange");
		yRange.setStart(minY);
		yRange.setEnd(maxY);
		
		this.xSize = xSize;
		this.ySize = ySize;
		this.xTickStop =  xTickStop;
		this.yTickStop = yTickStop;
		//this.plotStartX = xMargin + 5 + AXIS_SCALE_SIZE;
		this.plotStartX = xMargin + xPosition;
		this.plotStartY = yMargin + yPosition;
		//this.plotEndY = this.ySize - 2*yMargin + yPosition;
		this.plotEndY = this.ySize  + yPosition;
		System.out.println("margin "+ yMargin + " ySize " + ySize + " xSize " + xSize);
		//this.plotEndX = this.xSize - plotStartX + AXIS_SCALE_SIZE;
		this.plotEndX = this.xSize + plotStartX ;
		System.out.println("plotStartX "+plotStartX + " plotStartY "+plotStartY +" plotEndX "+plotEndX+" plotEndY "+plotEndY +" maxY " + maxY + " minY " + minY + " maxX " + maxX + " minX " + minX);
		title = "Plot";
	}
		
	public void paintCoordinateSystem(Graphics2D g) {
		paintXCoordinateSystem(g);
		paintYCoordinateSystem(g);
	}
	
	public  void paintXCoordinateSystem(Graphics2D g) {
		Color originalColor = g.getColor();
		g.setPaint(Color.black);
		g.drawLine(plotStartX,plotEndY+1,plotEndX,plotEndY+1);

		int maxX = getXRange().getEnd();
		int minX = getXRange().getStart();
		int plotXTickStop = (int) Math.round((DEFAULT_SCALE*xTickStop*(plotEndX - plotStartX))/(maxX - minX));
		System.out.println("Plot X tick stop " + plotXTickStop);
		//Draw ticks
		float xTickNum = plotXTickStop > 0 ? Math.round((plotEndX - plotStartX)/plotXTickStop) + 1 : 0;		
		for (int i = 1; i <= xTickNum; i++) {
			long pos = Math.round(xRange.getStart() /xScale) + i*xTickStop;
			int x = Math.round(plotStartX+i*plotXTickStop);
			if(x <= plotEndX) {
				g.drawLine(x,plotEndY+1,x,plotEndY+30);			
				g.drawString(String.valueOf(pos),Math.round(plotStartX+i*plotXTickStop),plotEndY+50);
			}
			/**for (int j = 0 ; j < xTickStop; j++) {
				x = Math.round(plotStartX+(i-1)*plotXTickStop +j*plotXTickStop/xTickStop);
				if(x <= plotEndX) {
					g.drawLine(x,plotEndY+1,x,plotEndY+40);
				}
			}*/

		}
		g.setColor(originalColor);
	}


	public void paintYCoordinateSystem(Graphics2D g) {
		Color originalColor = g.getColor();
		g.setPaint(Color.black);
		g.drawLine(plotStartX-1,plotStartY,plotStartX-1,plotEndY);
		int maxY = getYRange().getEnd();
		int minY = getYRange().getStart();
		
		int plotYTickStop = (int) Math.round((DEFAULT_SCALE*yTickStop*(plotEndY - plotStartY))/(maxY - minY));
		
		int yTickNum = plotYTickStop > 0 ? Math.round((plotEndY - plotStartY)/plotYTickStop) : 0;
		for (int i = yTickNum; i >= 1 ; i--) {
			long pos = Math.round(yRange.getStart()/yScale) + i*yTickStop;
			int y = Math.round(plotEndY - i * plotYTickStop );
			//System.out.println("Ploting "+plotStartX+","+Math.round(plotEndY-i*yTickInterval)+","+(plotStartX-4)+","+Math.round(plotEndY-i*yTickInterval));
			g.drawLine(plotStartX,y,plotStartX-30,y);
			g.drawString(String.valueOf(pos),plotStartX-100,Math.round(plotEndY-i*plotYTickStop));
			/**for (int j = 0 ; j < yTickStop; j++) {
				y = Math.round(plotEndY-(i)*plotYTickStop +j*plotYTickStop/yTickStop);
				if(y >= plotStartY) {
					g.drawLine(plotStartX+1, y,plotStartX-40, y);
				}
			}*/
		}
		g.setColor(originalColor);
	}

	public LightweightGenomicAnnotation getXRange() {
		return xRange;
	}

	public LightweightGenomicAnnotation getYRange() {
		return yRange;
	}
	
	public void reScaleX(long scale) {
		this.xScale = scale;
	}
	
	public void reScaleY(long scale) {
		this.yScale = scale;
	}
	
	/* (non-Javadoc)
	 * @see edu.mit.broad.prodinfo.genomicplot.BaseToPlotCoordinateMapper#getPlotXCoordinate(long)
	 */
	public int getPlotXCoordinate(long sequenceCoord) {
		
		int x = plotStartX;
		if(sequenceCoord > xRange.getEnd()) {
			x = plotEndX;
		} else if (sequenceCoord >= xRange.getStart()) {
			x = (int) Math.round((sequenceCoord - xRange.getStart())*(plotEndX - plotStartX)/(xRange.length())) + plotStartX;
		}
		
		return x;
		
		//return (int) Math.round((sequenceCoord - xRange.getStart())*(plotEndX - plotStartX)/(xRange.getLength())) + plotStartX;
	}
	
	/* (non-Javadoc)
	 * @see edu.mit.broad.prodinfo.genomicplot.BaseToPlotCoordinateMapper#getPlotYCoordinate(long)
	 */
	public int getPlotYCoordinate(long sequenceCoord) {
		
		int y = plotEndY;
		
		if(sequenceCoord > yRange.getEnd()) {
			y = plotStartY;
		} else if (sequenceCoord >= yRange.getStart()) {
		//System.out.println("minY: " + minY+" maxY " + maxY +" Y Coordinate "+sequenceCoord+" is translated to "+(plotEndY - Math.round((sequenceCoord - minY)*(plotEndY - plotStartY)/(maxY - minY))));
			y = plotEndY - Math.round((sequenceCoord - yRange.getStart())*(plotEndY - plotStartY)/(yRange.length()));
		}
		
		return y;
		
		//return plotEndY - Math.round((sequenceCoord - yRange.getStart())*(plotEndY - plotStartY)/(yRange.getLength()));
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

	public void plotItem(Graphics2D g, Plottable item) {
		item.plot(g, this);
	}



	public int getPlotEndX() {
		return plotEndX;
	}

	public int getPlotEndY() {
		return plotEndY;
	}

	public int getPlotStartX() {
		return plotStartX;
	}

	public void setXSeqId(String seqId) {xRange.setChromosome(seqId);}
	public String getXSeqIt() { return xRange.getChromosome();}
	
	public void setYSeqId(String seqId) {yRange.setChromosome(seqId);}
	public String getYSeqIt() { return yRange.getChromosome();}
	
	public int getXPlottedLength() { return xRange.length(); }

	public int getPlotStartY() {
		return plotStartY;
	}

	public String getTitle() {
		return title;
	}

	public void setTitle(String title) {
		this.title = title;
	}

	public long getXScale() {
		return xScale;
	}

	public int getXSize() {
		return xSize;
	}

	public long getYScale() {
		return yScale;
	}

	public int getYSize() {
		return ySize;
	}

	public int getXTickStop() {
		return xTickStop;
	}

	public int getYTickStop() {
		return yTickStop;
	}


}
