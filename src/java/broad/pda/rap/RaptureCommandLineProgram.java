package broad.pda.rap;

import java.io.File;
import java.io.IOException;

import broad.core.sequence.SequenceUtils;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.util.Log;

public abstract class RaptureCommandLineProgram extends GenomeCommandLineProgram {
    private static final Log log = Log.getInstance(RaptureCommandLineProgram.class);
	
	@Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input SAM or BAM file.")
	public File INPUT;

	@Option(shortName="O", doc="File to write the output to.")
	public File OUTPUT;

	@Option(shortName="E", doc="Extension factor for reads.", optional=true)
	public int EXTENSION_FACTOR = 0;

	@Option(doc="Load read pairs as fragments.", optional=true)
	public boolean PAIRED_END = false;
	
	@Option(doc="Chunk size to load for cached interval tree")
	public int CHUNK_SIZE = 1000000;

	//@Option(doc="Exclude reads with mapping quality less than this.", optional=true)
	//public int MIN_MAPPING_QUALITY = 0;
	
	
	@Override
	protected String[] customCommandLineValidation() {
		if (OUTPUT != null && INPUT.equals(OUTPUT)) {
			return new String[]{"INPUT file and OUTPUT file must differ!"};
		}
		if (PAIRED_END && EXTENSION_FACTOR > 0) {
			return new String[]{"Paired-end and extension factor should not both be specified."};
		}
		return super.customCommandLineValidation();
	}
	
	public ContinuousDataAlignmentModel getContinuousDataModel(File file) throws IOException {
		ContinuousDataAlignmentModel model = SequenceUtils.getDataModel(file.getAbsolutePath(), SIZES.getAbsolutePath(), 0, true, false, null, PAIRED_END);
		model.setExtensionFactor(EXTENSION_FACTOR);
		model.getAlignmentDataModelStats().setChunkSize(CHUNK_SIZE);
		return model;
	}

}
