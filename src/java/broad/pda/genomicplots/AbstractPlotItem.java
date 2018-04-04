package broad.pda.genomicplots;

import java.awt.Color;

public abstract class AbstractPlotItem implements Drawable{
	private Color color;
	private String caption;
	
	
	/* (non-Javadoc)
	 * @see edu.mit.broad.prodinfo.genomicplot.Plottable#getColor()
	 */
	public Color getColor() {
		return color;
	}

	public void setColor(Color color) {
		this.color = color;
	}

	public String getCaption() {
		return caption;
	}

	/* (non-Javadoc)
	 * @see edu.mit.broad.prodinfo.genomicplot.Plottable#setCaption(java.lang.String)
	 */
	public void setCaption(String caption) {
		this.caption = caption;
	}

}