package broad.pda.genomicplots;

public interface BaseToPlotCoordinateMapper {

	public abstract int getPlotXCoordinate(long sequenceCoord);

	public abstract int getPlotYCoordinate(long sequenceCoord);

	public abstract int getPlotYDistance(long start, long end);

	public abstract int getPlotXDistance(long start, long end);
	

	public int getPlotEndX();

	public int getPlotEndY();

	public int getPlotStartX();

	public int getPlotStartY();
}