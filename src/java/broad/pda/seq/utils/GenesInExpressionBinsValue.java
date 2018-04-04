package broad.pda.seq.utils;

import java.io.IOException;
import java.util.Collection;

import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;

public interface GenesInExpressionBinsValue {
	public void computeValue(Collection<RefSeqGene> genes, ContinuousDataAlignmentModel data, String chr) throws IOException;
	public boolean couldComputeValue();
}
