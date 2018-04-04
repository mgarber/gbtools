package broad.pda.rap;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.ShortBEDReader;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.slide.Slide;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.util.Log;

public abstract class GenomeCommandLineProgram extends CommandLineProgram{
    private static final Log log = Log.getInstance(GenomeCommandLineProgram.class);
    
	@Option(doc="File specifying chromosome sizes.")
	public File SIZES = new File("/seq/lincRNA/data/mm9/sizes");
	
	@Option(doc="File containing masked regions.", optional=true)
	public File MASK_FILE = new File("/seq/mguttman/ChIPData/MaskFiles/MM9Segments/all.mask.mouse.n36.d2.bin");
	
	@Option(doc="Region to process (e.g. chr1, chr1:5000-50230).  Default is the entire genome", optional=true)
	public String REGION = null;
	

	protected Map<String, Integer> sizes;
	protected ShortBEDReader maskedRegions = null;	
	
	
	@Override
	protected String[] customCommandLineValidation() {
		loadSizes();
		loadMaskedRegions();
		return super.customCommandLineValidation();
	}
	

	public List<GenomicAnnotation> getRegions() {
		List<GenomicAnnotation> regions = new ArrayList<GenomicAnnotation>();

		if (REGION != null) {
			if (sizes.containsKey(REGION)) {
				regions.add(new Alignments(REGION, 1, sizes.get(REGION)));
			} else {
				try {
					regions.add(new Alignments(REGION));
				} catch (RuntimeException e) {
					throw new IllegalArgumentException("REGION is improperly formatted");
				}
			}
		} else {
			for (Map.Entry<String,Integer> entry : sizes.entrySet()) {
				regions.add(new Alignments(entry.getKey(), 1, entry.getValue()));
			}
		}
		return regions;
	}
	
	
	
	protected void loadSizes() {
		sizes = BEDFileParser.loadChrSizes(SIZES.getAbsolutePath());
	}
	
	protected void loadMaskedRegions() { 
		try {
			if (MASK_FILE != null) maskedRegions = Slide.loadMaskedRegions(MASK_FILE);
		} catch (Exception e) {
			log.error(e);
		}
	}
	
	protected ShortBEDReader getMaskedRegions() { return maskedRegions; }
	protected Map<String, Integer> getSizes() { return sizes; }
}
