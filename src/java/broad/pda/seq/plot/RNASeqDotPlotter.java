package broad.pda.seq.plot;

import java.awt.BasicStroke;
import java.awt.Color;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jibble.epsgraphics.EpsGraphics2D;

import broad.core.annotation.AnnotationFactoryFactory;
import broad.core.annotation.BED;
import broad.core.annotation.BEDReader;
import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.BasicTwoSubjectAnnotation;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.GenomicAnnotationFilter;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.error.ParseException;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.alignment.PairedEndAlignment;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.genomicplots.Plot;
import broad.pda.genomicplots.PlotItem;


public class RNASeqDotPlotter {
	public static final String USAGE = "Usage: RNASeqDotPlotter TASK=<task_num> <task_args>\n" +
	"\tTasks:\n" +
	"\n\t1. Plot reads in two dimensions -out <Output file> -chr <Chromosome> -start <Region start> -end <Region end> -aln1 <Alignment of first set> -aln2 <Alignment of paireds if available> [-alnformat <Alignment read format, currently BED is supported>] " +
	"\n";
	private static final int RESOLUTION     = 1000;
	private static final int  PLOT_SIZE = (int)(0.9 * RESOLUTION);

	private LightweightGenomicAnnotation region;
	private RegionFilter filter;
	private Plot plot;
	
	public static void main(String [] args) throws IOException, ParseException {
		//ArgumentMap argMap = CLUtil.getParameters(args, USAGE, "1");
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE);
		
		if("1".equals(argMap.getTask()) ) {
			GenomicAnnotation region = new BasicGenomicAnnotation("region",argMap.getMandatory("chr"), argMap.getInteger("start"), argMap.getInteger("end"));
			RNASeqDotPlotter plot = new RNASeqDotPlotter(region);
			String out = argMap.getOutput();
			EpsGraphics2D epsG2D = plot.startPlot(out);
			
			BEDReader reader1 = new BEDReader();
			File aln1File = new File(argMap.getMandatory("aln1"));
			reader1.load(aln1File, AnnotationFactoryFactory.bedFactory, plot.filter);
			Map<String, GenomicAnnotation> aln1Map = new HashMap<String, GenomicAnnotation>();
			List<BED> annotations = reader1.getAnnotationList();
			for (BED a : annotations) {
				aln1Map.put(a.getName(), a);
				for(GenomicAnnotation block : a.getBlocks()) {
					plot.plot(block, epsG2D);
				}
			}
			
			plot.endPlot(epsG2D);
		} else {
			System.err.println("No task Specified " + USAGE);
		}
	}


	public RNASeqDotPlotter (LightweightGenomicAnnotation region) {
		this.region = region;
		System.err.println("REGION: "+ region);
		filter = new RegionFilter(region);
		plot = new Plot(PLOT_SIZE, PLOT_SIZE, region.getEnd(), region.getEnd(), region.getStart(), region.getStart(),5, 5);
	}
	
	protected void plot(GenomicAnnotation read, EpsGraphics2D epsG2D) {
		BasicTwoSubjectAnnotation selfRead = new BasicTwoSubjectAnnotation();
		selfRead.setA(read);
		selfRead.setB(read);
		PlotItem item = new PlotItem(selfRead);
		item.setColor(Color.BLACK);
		plot.plotItem(epsG2D, item);	
	}
	
	protected void plot(RefSeqGene read, EpsGraphics2D epsG2D) {
		//loop through all exon connections and plot self and pairs
		
		Alignments[] exons=read.getExons();
		for(int i=0; i<exons.length; i++){
			Alignments exon=exons[i];
			/*Self plot*/
			BasicTwoSubjectAnnotation selfRead = new BasicTwoSubjectAnnotation();
			selfRead.setA(exon);
			selfRead.setB(exon);
			PlotItem item = new PlotItem(selfRead);
			item.setColor(Color.BLACK);
			plot.plotItem(epsG2D, item);
			
			for(int j=i+1; j<exons.length; j++){
				Alignments nextExon=exons[j];
				selfRead = new BasicTwoSubjectAnnotation();
				selfRead.setA(exon);
				selfRead.setB(nextExon);
				item = new PlotItem(selfRead);
				item.setColor(Color.BLUE);
				plot.plotItem(epsG2D, item);
				
				selfRead = new BasicTwoSubjectAnnotation();
				selfRead.setA(nextExon);
				selfRead.setB(exon);
				item = new PlotItem(selfRead);
				item.setColor(Color.BLUE);
				plot.plotItem(epsG2D, item);
			}
			
		}
		
		
	}
	
	protected void plot(PairedEndAlignment pair, EpsGraphics2D epsG2D) {
		//loop through all exon connections and plot self and pairs
		plot(pair.getLeftMate(), epsG2D);
		plot(pair.getRightMate(), epsG2D);
		
		//connect mates with a gray line
		//left exons to right exons
		Alignments[] left =pair.getLeftMate().getExons();
		Alignments[] right=pair.getRightMate().getExons();
		
		for(int i=0; i<left.length; i++){
			for(int j=0; j<right.length; j++){
				BasicTwoSubjectAnnotation selfRead = new BasicTwoSubjectAnnotation();
				selfRead = new BasicTwoSubjectAnnotation();
				selfRead.setA(left[i]);
				selfRead.setB(right[j]);
				PlotItem item = new PlotItem(selfRead);
				item = new PlotItem(selfRead);
				item.setColor(Color.GRAY);
				plot.plotItem(epsG2D, item);
				
				selfRead = new BasicTwoSubjectAnnotation();
				selfRead = new BasicTwoSubjectAnnotation();
				selfRead.setB(left[i]);
				selfRead.setA(right[j]);
				item = new PlotItem(selfRead);
				item = new PlotItem(selfRead);
				item.setColor(Color.GRAY);
				plot.plotItem(epsG2D, item);
			}
		}
	}
	
	public EpsGraphics2D startPlot(String outputFile) throws IOException {
		FileOutputStream outputStream = new FileOutputStream(outputFile);
		EpsGraphics2D epsG2D = new EpsGraphics2D(region.toUCSC(), outputStream, 0, 0,(int)(1.2*RESOLUTION),(int) (1.2*RESOLUTION));
		java.awt.Font psFont = new java.awt.Font("arial", 0, 12);
		epsG2D.setFont(psFont);
		epsG2D.setStroke(new BasicStroke(2.0f));
		return epsG2D;
	}
	
	public void endPlot(EpsGraphics2D epsG2D) throws IOException {
		plot.paintXCoordinateSystem(epsG2D);
		plot.paintYCoordinateSystem(epsG2D);
		epsG2D.flush();
		epsG2D.close();
	}
	
	private static class RegionFilter implements GenomicAnnotationFilter<BED> {
		LightweightGenomicAnnotation region;
		RegionFilter(LightweightGenomicAnnotation region) {
			this.region = region;
		}
		public boolean accept(BED annotation) {
			return region.overlaps(annotation);
		}

		public boolean isEnough(BED annotation) {
			return false;
		}
		
	}

}
