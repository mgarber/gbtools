/**
 * 
 */
package broad.pda.seq.segmentation;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import broad.pda.datastructures.Alignments;

/**
 * @author engreitz
 * Wrapper for a list of ContinousDataAlignmentModel's
 */
public class MultipleDataAlignmentModel {
	
	// These are two alternative ways to store the data.  models is always used
	// modelMap is only used if names are provided in the constructor
	private List<ContinuousDataAlignmentModel> models;
	private Map<String, ContinuousDataAlignmentModel> modelMap;
	
	public MultipleDataAlignmentModel(List<ContinuousDataAlignmentModel> models) {
		this.models = models;
		modelMap = null;
	}
	
	public MultipleDataAlignmentModel(List<ContinuousDataAlignmentModel> models, List<String> modelNames) {
		this(models);
		if (models.size() != modelNames.size()) {
			throw new IllegalArgumentException();
		}
		modelMap = new HashMap<String, ContinuousDataAlignmentModel>();
		
		for (int i = 0; i < modelNames.size(); i++) {
			modelMap.put(modelNames.get(i), models.get(i));
		}
	}
	
	
	public ContinuousDataAlignmentModel getModel(int index) {
		return models.get(index);
	}
	
	public ContinuousDataAlignmentModel getModel(String key) {
		return modelMap.get(key);
	}
	
	public int size() { return models.size(); }
	public boolean hasNames() { return (modelMap != null); }
	
	/**
	 * @param window
	 * @return				List of counts in the window in each model (in order)
	 * @throws IOException
	 */
	public List<Integer> getCount(Alignments window) throws IOException {
		List<Integer> counts = new ArrayList<Integer>();
		for (ContinuousDataAlignmentModel model : models) {
			counts.add((int) model.count(window));
		}
		return counts;
	}
	
	
	public List<Integer> getTotalReads() throws IOException {
		List<Integer> reads = new ArrayList<Integer>();
		for (ContinuousDataAlignmentModel model : models) {
			reads.add(model.getAlignmentDataModelStats().getTotalReads());
		}
		return reads;
	}
}
