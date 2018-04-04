package broad.pda.seq.slide;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.collections15.CollectionUtils;
import org.apache.commons.collections15.Predicate;

import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.alignment.AlignmentCollection;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;


public class SlideAndCount extends broad.pda.rap.RaptureCommandLineProgram {
    private static final Log log = Log.getInstance(SlideAndCount.class);
	
    @Usage
    public String USAGE = "Slides across the genome and counts reads (or calculates ratios for two Alignment models).\nTesting on Oct 10, 2012 required less than 4G Mem (+9G virtual) for 150M reads.";
   
	@Option(doc="Window size")
	public int WINDOW;
	
	@Option(doc="Overlap between windows")
	public int OVERLAP;
	
	@Option(doc="Control SAM or BAM file for normalization.", optional=true)
	public File CONTROL = null;
	
	@Option(doc="Windows with scores less than this will not be printed", optional=true)
	public Double LOWER_LIMIT = 0.0;
	
	@Option(doc="Windows with scores greater than this will not be printed", optional=true)
	public Double UPPER_LIMIT = Double.MAX_VALUE;
	
	@Option(doc="Invert the filter", optional=true)
	public boolean INVERT_FILTER = false;
	
	@Option(doc="Calculate and print extra statistics")
	public boolean FULL_OUTPUT = false;
	
	@Option(doc="Percent masked allowable per sliding window", optional=true)
	public double PCT_MASKED_ALLOWED = 20.0;
	
	
	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new SlideAndCount().instanceMain(args));
	}
	

	@Override
	protected int doWork() {
		
		try {
			IoUtil.assertFileIsReadable(INPUT);
			IoUtil.assertFileIsWritable(OUTPUT);
		
			if (OUTPUT.exists()) {
				OUTPUT.delete();
			}
			
			List<GenomicAnnotation> annot = getRegions();
			List<LightweightGenomicAnnotation> regions = new ArrayList<LightweightGenomicAnnotation>();
			for (GenomicAnnotation g : annot) regions.add(new BasicLightweightAnnotation(g.getChromosome(), g.getStart(), g.getEnd()));
			
			ContinuousDataAlignmentModel model = getContinuousDataModel(INPUT);
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(OUTPUT,true));
			Slide sac = new Slide(WINDOW, OVERLAP, getMaskedRegions(), PCT_MASKED_ALLOWED);
			
			SlideAndCountConsumer consumer;
			if (CONTROL != null) {
				ContinuousDataAlignmentModel control = getContinuousDataModel(CONTROL);
				consumer = new PrintRatioConsumer(model, control, bw, new CountFilter());
			} else {
				consumer = new PrintScoresConsumer(model, bw, FULL_OUTPUT, new CountFilter());
			}
			
			sac.slide(regions, consumer);
			bw.close();
			
		} catch (Exception e) {
			log.error(e);
		}
		
		return 0;
	}
	
	
	/**
	 * @author engreitz
	 * Consumer class to print out the count scores
	 */
	public static class PrintScoresConsumer extends SlideAndCountConsumer {
		private boolean fullOutput;
		private BufferedWriter bw;
		private Predicate<Alignments> windowFilter = null;
		
		private PrintScoresConsumer(AlignmentCollection data, BufferedWriter bw, boolean fullOutput) throws IOException {
			super(data);
			this.bw = bw;
			this.fullOutput = fullOutput;
		}
		
		private PrintScoresConsumer(AlignmentCollection data, BufferedWriter bw, boolean fullOutput, Predicate<Alignments> p) throws IOException {
			this(data, bw, fullOutput);
			this.windowFilter=p;
		}
		
		
		@Override
		public void initRegion(LightweightGenomicAnnotation region) {
			if (fullOutput) super.initRegion(region);
		}
		
		@Override
		public void consume(List<Alignments> windows) throws Exception {
			super.consume(windows);   // VERY IMPORTANT TO INCLUDE THIS

			if (windowFilter != null) {
				CollectionUtils.filter(windows, windowFilter);
			}
			
			for(Alignments w : windows) {
				bw.write(w.toUCSC());
				bw.write("\t");
				bw.write(String.valueOf(w.getCountScore()));
				
				if (fullOutput) {
					bw.write("\t" + w.getChr() + "\t" + w.getStart() + "\t" + w.getEnd() + "\t"+
							 getRegionReads() + "\t" + getTotalReads());
				}
				
				bw.newLine();
			}
		}
	}
	
	
	public static class PrintRatioConsumer extends CountRatioConsumer {
		private BufferedWriter bw;
		private Predicate<Alignments> p = null;
		
		public PrintRatioConsumer(AlignmentCollection numerator, AlignmentCollection denominator, BufferedWriter bw) throws IOException {
			super(numerator, denominator);
			this.bw = bw;
		}
		
		public PrintRatioConsumer(AlignmentCollection numerator, AlignmentCollection denominator, BufferedWriter bw, Predicate<Alignments> p) throws IOException {
			this(numerator, denominator, bw);
			this.p = p;
		}
		
		
		@Override
		public void consume(List<Alignments> windows) throws Exception {
			super.consume(windows);   // VERY IMPORTANT TO INCLUDE THIS

			if (p != null) {
				CollectionUtils.filter(windows, p);
			}
			
			Alignments.write(bw, windows, false);
		}
	}


	
	/**
	 * Filters Alignments based on their score
	 * Usage:  Predicate<Alignments> filter = new PValueFilter(cutoff, fold);  CollectionUtils.filter(alignments, filter);
	 * @author engreitz
	 *
	 */
	public class CountFilter implements Predicate<Alignments> {
		public boolean evaluate(Alignments obj) {
			boolean result = (obj.getCountScore() >= LOWER_LIMIT) && (obj.getCountScore() <= UPPER_LIMIT);
			if (INVERT_FILTER) result = !result;
			return result;
		}
	}
}
