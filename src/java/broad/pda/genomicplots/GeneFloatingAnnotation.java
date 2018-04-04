package broad.pda.genomicplots;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.util.Iterator;

import broad.core.gene.GeneAnnotation;
import broad.core.gene.GeneAnnotation.Exon;

public class GeneFloatingAnnotation extends FloatingAnnotation {
	GeneAnnotation gene;
	
	public GeneFloatingAnnotation(GeneAnnotation gene, Color color) {
		super(gene, color);
		this.gene = gene;		
	}
	
	public void plot(Graphics2D g, GenomicAnnotationPlotter plot, boolean displayName) {
		//int start = plot.getPlotYCoordinate(gene.getStart());
		int end   = plot.getPlotYCoordinate(gene.getEnd());
		int x     = plot.getPlotStartX() - SEPARATION ;
		Color originalColor = g.getColor();		
		g.setColor(getColor());

		//System.out.println("Plotting " + gene + " exon num " + gene.getExonCount() + " but it has " + gene.getExons().size());
		Iterator<Exon> exonIt = gene.getExons().iterator();
		Exon currentExon = exonIt.next();
		Exon nextExon    = null;
		while(exonIt.hasNext()) {
			//System.out.println("\tPlotting " + currentExon);
			int currExonPlotStart = plot.getPlotYCoordinate(currentExon.getStart());
			int currExonPlotEnd   = plot.getPlotYCoordinate(currentExon.getEnd());
			g.fillRect(x - 2, 
					Math.min(currExonPlotStart, currExonPlotEnd), 
					4, 
					Math.abs(currExonPlotEnd - currExonPlotStart) + 1);
			
			nextExon = exonIt.next();
			
			int midIntron = Math.round((nextExon.getStart() + currentExon.getEnd())/2);
			int midIntronY = plot.getPlotYCoordinate(midIntron);
			int midIntronX = gene.inReversedOrientation() ? x + 2 : x - 2;
			
			int nextExonStartY = plot.getPlotYCoordinate(nextExon.getStart());
			int nextExonEndY   = plot.getPlotYCoordinate(nextExon.getEnd());
			g.drawLine(x,
					currExonPlotEnd, 
					midIntronX,
					midIntronY);
			g.drawLine(midIntronX,
					midIntronY, 
					x,
					nextExonStartY);
			
			if(!exonIt.hasNext()) {
				g.fillRect(x - 2, 
						Math.min(nextExonStartY, nextExonEndY), 
						4, 
						Math.abs(nextExonEndY - nextExonStartY) + 1);
			}
			currentExon = nextExon;
		}

		java.awt.Font originalFont = g.getFont();
		java.awt.Font textFont     = new java.awt.Font("arial", 0, 8).deriveFont(AffineTransform.getRotateInstance(- Math.PI/2));
		g.setFont(textFont);
		g.drawString(getName(), x - 6, end);
		g.setFont(originalFont);

		g.setColor(originalColor);
	}
	
	public String getName() {
		return gene.getAltName();
	}
}
