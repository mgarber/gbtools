package broad.pda.genomicplots;

import java.awt.Color;

public interface Drawable {

	public abstract void setColor(Color color);

	public abstract Color getColor();

	public abstract void setCaption(String caption);

}