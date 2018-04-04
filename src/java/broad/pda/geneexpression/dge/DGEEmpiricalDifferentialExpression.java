package broad.pda.geneexpression.dge;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math.util.MathUtils;
import org.apache.log4j.Logger;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.differentialExpression.DifferentialExpression;

public class DGEEmpiricalDifferentialExpression {
	static Logger logger = Logger.getLogger(DGEEmpiricalDifferentialExpression.class.getName());

	static final String usage = "Usage: DGEEmpiricalDifferentialExpression "+
			"\n\t-in <Normalized table with expression values of genes by conditions>" +
			"\n\t-sampleInfo <Three column file describing the different samples, first column is the sample name, third column must contain the group to which it belings>" +
			"\n\t-alpha <minimum ignificance level, default is 0.01>" +
			"\n\t-outdir <outputdirectory>" +
			"\n";
	
	public static void main(String [] args) throws IOException, ParseException {
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"DE");
		String outdir = argMap.getOutputDir();
		File outdirFile = new File(outdir);
		if(! outdirFile.exists()) {
			outdirFile.mkdir();
		}
		MatrixWithHeaders data = new MatrixWithHeaders(argMap.getInput()); 
		Map<String, Collection<String>> groups = loadSampleInfo(argMap.getMandatory("sampleInfo"));
		double alpha = argMap.containsKey("alpha") ? argMap.getDouble("alpha") : 0.01;
		logger.debug("Groups: " + groups);
		int maxPermuations = 1000;//computeMaxPermutations(groups);
		logger.debug("Maximum permutations possible: " + maxPermuations);
		List<String> groupNames = new ArrayList<String>(groups.keySet());
		for (int i = 0; i < groupNames.size() - 1; i++) {
			for (int j = i+1; j < groupNames.size() ; j++) {
				String group1 = groupNames.get(i);
				String group2 = groupNames.get(j);
				logger.debug("Comparing " + group1 + ": " + groups.get(group1) + " with " + group2 + ": " + groups.get(group2));
				DifferentialExpression de = new DifferentialExpression(data, groups.get(group1), groups.get(group2), maxPermuations);
				MatrixWithHeaders fdrMatrix = de.getFDRMatrix();
				fdrMatrix.write(outdir + "/" + group1+".vs."+group2+".fdr.txt");
			}
		}
	}

	/**
	 * Take the two largest groups and compute the maximum number of ways you can swap 1 or N-1 elements
	 * @param groups
	 * @return
	 */
	private static int computeMaxPermutations(Map<String, Collection<String>> groups) {
		int [] sizes = new int [groups.size()];
		int i = 0;
		for(Collection<String> groupSamples : groups.values()) {
			sizes[i++] = groupSamples.size();
		}
		Arrays.sort(sizes);
		
		int permNumber = 0;
		if(groups.size() > 1) {
			int largest = sizes[sizes.length - 1];
			int secondLargest = sizes[sizes.length - 2];
			for (int k = 0; k < secondLargest; k++) {
				permNumber += MathUtils.binomialCoefficient(secondLargest, k) * MathUtils.binomialCoefficient(largest, k);
			}
		}
		return permNumber;
	}

	private static Map<String, Collection<String>> loadSampleInfo(String expInfoFile) throws IOException {
		Map<String, Collection<String>> sampleInfo = new HashMap<String, Collection<String>>();
		BufferedReader br = new BufferedReader(new FileReader(expInfoFile));
		String line = null;
		while( (line = br.readLine()) != null) {
			if (line.startsWith("#")) { continue; }
			logger.debug("line: " + line);
			String [] info = line.split("\t");
			String group = info[2];
			String sample = info[0];
			logger.debug("Adding " + sample + " to " + group);
			if(!  sampleInfo.containsKey(group) ) {
				sampleInfo.put(group, new ArrayList<String>());
			}
			Collection<String> groupSamples = sampleInfo.get(group);
			groupSamples.add(sample);
		}
		br.close();
		return sampleInfo;
	}
}
