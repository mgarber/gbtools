/**
 * 
 */
package broad.pda.rap;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import broad.core.error.ParseException;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.slide.Slide;
import broad.pda.seq.slide.SlideConsumer;

/**
 * @author engreitz
 *
 */
public class PlotAggregateRegions {
	Rapture rap;
	int innerLength, outerLength, middleLength, window;
	
	public PlotAggregateRegions(Rapture rap, int inner, int outer, int middle, int window) {
		this.rap = rap;
		innerLength = inner;
		outerLength = outer;
		middleLength = middle;
		this.window = window;
	}
	
	
	public void calculateScores(String in, String out) throws Exception, IOException, ParseException {
		Set<Alignments> regions = BEDFileParser.loadAlignmentData(new File(in));
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		
		Slide slider = new Slide(window, 0);
		for (Alignments region : regions) {
			
			PlotRegions subregions = generateSubregions(region);
			GeneConsumer consumer = new GeneConsumer(bw, region.getName());
			
			consumer.setSubname("beginOuter");
			slider.slide(subregions.beginOuter, consumer);
			
			consumer.setSubname("beginInner");
			slider.slide(subregions.beginInner, consumer);
			
			consumer.setSubname("middle");
			slider.slide(subregions.middle, consumer);
			
			consumer.setSubname("endInner");
			slider.slide(subregions.endInner, consumer);
			
			consumer.setSubname("endOuter");
			slider.slide(subregions.endOuter, consumer);
		}
	}
	
		
	public PlotRegions generateSubregions(Alignments region) throws IOException {
		PlotRegions subregions = new PlotRegions();
		
		Alignments outerLeft = new Alignments(region.getChr(), region.getStart() - outerLength, region.getStart());
		Alignments outerRight = new Alignments(region.getChr(), region.getEnd(), region.getEnd() + outerLength);
		Alignments innerLeft = new Alignments(region.getChr(), region.getStart(), region.getStart() + innerLength);
		Alignments innerRight = new Alignments(region.getChr(), region.getEnd() - innerLength, region.getEnd());
		// TODO:  Adjust so that these regions don't extend beyond the length of a short gene
		
		if (region.getStrand().equals("+")) {
			subregions.beginOuter = outerLeft;
			subregions.beginInner = innerLeft;
			subregions.endOuter = outerRight;
			subregions.endInner = innerRight;
		} else if (region.getStrand().equals("-")) {
			subregions.beginOuter = outerRight;
			subregions.beginInner = innerRight;
			subregions.endOuter = outerLeft;
			subregions.endInner = innerLeft;
		} else {
			throw new IOException("Invalid strand " + region.getStrand());
		}
		
		int midway = (int) Math.floor((region.getEnd() + region.getStart())/2.0);
		int midAdjust = (int) Math.floor(middleLength/2.0);
		subregions.middle = new Alignments(region.getChr(), midway - midAdjust, midway - midAdjust + middleLength);
		
		return subregions;
	}
	
	
	public class GeneConsumer extends SlideConsumer.AbstractSlideConsumer {
		String name, subname;
		BufferedWriter bw;
		int counter = 0;
		public GeneConsumer(BufferedWriter bw, String name) { this.name = name; this.bw = bw; }
		public void setSubname(String name) { subname = name; }
		
		public int getWindowBatch() { return 1000000; }
		public void consume(List<Alignments> windows) throws Exception {
			rap.scoreSegments(windows);
			for (Alignments window : windows) {
				printScore(window);
			}
		}
		
		public void printScore(Alignments align) throws IOException {
			bw.write(name + "\t" + subname + "\t" + counter + "\t");
			bw.write(align.toString());
			for (Double score : align.getScores()) {
				bw.write("\t" + score);
			}
			bw.write("\n");
			counter = counter + window;
		}
	}
	
	
	
	private class PlotRegions {
		public Alignments beginOuter, beginInner, middle, endInner, endOuter;
	}
	
	
	private static String USAGE = "java -cp Rapture.jar broad.pda.rap.PlotAggregateRegion [args] TODO"; // TODO
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception, IOException, ParseException {
		// TODO Auto-generated method stub
		ArgumentMap argmap = CLUtil.getParameters(args, USAGE);
		Rapture rap = Rapture.newFromArgs(argmap);

		String out = argmap.getOutput();
		String in = argmap.getInput();
		
		int innerLength = argmap.getInteger("innerLength", 1000);
		int outerLength = argmap.getInteger("outerLength", 2000);
		int middleLength = argmap.getInteger("middleLength", 1000);
		int window = argmap.getInteger("window", 10);
	
		PlotAggregateRegions p = new PlotAggregateRegions(rap, innerLength, outerLength, middleLength, window);
		p.calculateScores(in, out);
	}

}
