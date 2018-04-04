package broad.pda.seq.alignment;

import java.io.IOException;

import broad.pda.datastructures.Alignments;

/**
 * @author engreitz
 * Interface to allow for counting over both AnnotationReader and ContinuousDataAlignmentModel.  
 * Should eventually switch to Generic instead of Continuous model
 */
public interface AlignmentCollection {
	public int size() throws IOException;
	public int getCount(Alignments region) throws IOException;
	public int getBasesCovered(Alignments region, int extensionFactor) throws IOException;
	public int getBasesCovered(Alignments region) throws IOException;
}
