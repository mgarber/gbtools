package broad.pda.seq.slide;

import java.io.IOException;

import broad.pda.datastructures.Alignments;

public interface AlignmentModel {
	public double getScore(Alignments region) throws IOException;
}
