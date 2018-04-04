package broad.pda.genomicplots;

import java.awt.Color;
import java.awt.Graphics2D;

public class TriangularAnnotation  {

	private long position;
	private Color color;
	int hight;
	int base;
	boolean inverted = false;
	private String caption;
	private int offset;
	
	public TriangularAnnotation(long position, int hight, int base) {
		this.position = position;
		this.hight = hight;
		this.base  = base;
		this.offset = 0;
	}

	public TriangularAnnotation(long position, int hight, int base, int offset) {
		this.position = position;
		this.hight = hight;
		this.base  = base;
		this.offset = offset;
	}

	public void plot(Graphics2D g, BaseToPlotCoordinateMapper plot, boolean bottom) {
		plot(g, plot, 0, bottom);
	}

	
	public void plot(Graphics2D g, BaseToPlotCoordinateMapper plot, int yOffset, boolean bottom) {
		plot(g, plot, yOffset, !bottom, true);
	}
	
	public void plot (Graphics2D g, BaseToPlotCoordinateMapper plot, int offset, boolean beginning, boolean vertical) {
		g.setPaint(color);
		int yStart = 0;
		int xPos = 0;
		if (vertical) {
			yStart = !beginning ? plot.getPlotEndY() + offset : plot.getPlotStartY() - offset - hight;
			xPos = plot.getPlotXCoordinate(position);
		} else {
			yStart = plot.getPlotYCoordinate(position);
			xPos = beginning ? plot.getPlotStartX() - offset : plot.getPlotEndX() + offset + hight;
			
		}
		int [] verticeXCoords = new int[3];
		int [] verticeYCoords = new int[3];
		
		if (vertical) {
			verticeXCoords[0] =  xPos;
			verticeYCoords[0] = !inverted ? yStart : yStart + hight + offset;
		
			verticeXCoords[1] = Math.round(xPos - base/2);
			verticeYCoords[1] = !inverted ? yStart + hight : yStart + offset;
			verticeXCoords[2] = Math.round(xPos + base/2);
			verticeYCoords[2] = !inverted ? yStart + hight : yStart + offset;
		} else {
			verticeXCoords[0] =  !inverted ? xPos : xPos - hight - offset;
			verticeYCoords[0] = yStart;
		
			verticeYCoords[1] = Math.round(yStart - base/2);
			verticeXCoords[1] = !inverted ? xPos - hight : xPos - offset;
			verticeYCoords[2] = Math.round(yStart + base/2);
			verticeXCoords[2] = !inverted ? xPos - hight : xPos - offset;
			
		}
		
		g.fillPolygon(verticeXCoords, verticeYCoords, 3);
		
		if (caption != null) {
			for(int i = 0; i < caption.length(); i++) {
				g.drawString(String.valueOf(caption.charAt(i)),xPos, yStart + hight + (i+1)*80);
			}
		}
		
	}

	public void setColor(Color color) {
		this.color = color;
	}

	public void setCaption(String caption) {
		this.caption = caption;
	}
	
	public void setInverted(boolean invert) {
		inverted = invert;
	}

}
