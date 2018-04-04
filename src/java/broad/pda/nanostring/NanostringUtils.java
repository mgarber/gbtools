package broad.pda.nanostring;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.math.MathUtil;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.differentialExpression.DifferentialExpression;
import broad.pda.differentialExpression.DifferentialScoring;

public class NanostringUtils {
	public static final String USAGE = "Usage: NanostringUtils TASK=<task_num> <task_args>\n" +
	"\tTasks:\n" +
	"\t\t1. Process raw file, that is normalize by spikein controls, then normalize by given experimental controls, and finally compute scores and FDRs for differentially expressed genes " +
	"\n\t\t-in <File with nanostring readdings or standard input>  "+
	"\n\t\t-controlExperiments <comma separated experiment column numbers  to be used as controls (e.g. empty shRNA, fake stimulii) NOTE that first experiment is 0> " +
	"\n\t\t[-controlTranscripts <Comma sperated list of genes to normalize observations accross experiments>] "+
	"\n\t\t-useAsReference <Experiment (Column) to use as reference for normalization purposes if other than the first NOTE: count is zero based, first column is 0>]" +
	"\n\t\t-minConfidence <Only report values for observations with a z-scores extreme enough the ansure at the given confidence level they are not random experimental fluctuationts. All other values are set to 0>]" +
	"\n\t\t-gp <use gene pattern friendly format>" +
	"\n\t\t -groupFile <Points to a file containing groupings. this is a two column, tab separated file, where the first column is the name of the group and the second column is a comma separated list of sample names. Sample names must correspond to the names in the Attributes line of the input nanostring file>" +
	"\n\t\t-cytoscape <Outputs all a .sif, .noa and .eda with data extracted from the putative interactions> -networkName <If specified all files will be prefixed by the provided network name> -ignoreExperiment <Experiments not to output as nodes, add as many names as desired> -ignoreTranscript <Transcripts to avoid output as nodes>" +
	"\n\t\t-fdrType <0: by transcript, 1: by column, 2: by both, and 3: by transcript data (if using this option need to proved -fdrTranscripts <A comman separeted list of transcripts to use as null distribution> The default is to compute an FDR by transcript>" +
	"\n\t\t-scoreToWrite <0: Transcript count 1: confidence value 2: z-score 3: fold change to control average 4: FWER value. It defaults to z-score if control experiments are included otherwise it defaults to raw counts>" +
	"\n\t\t-minFoldChange <Minimum fold change to control SH average transcript exception required to pass confidence test> -max <If using the FWER which is low bound, set up the maximum FWER cutoff>" +
	"\n\t2. Similar to 1 but used to process several nanostring output files together, and can be usef for defining groups when there are replicates the input for this task is a tab separated file describing all nanostring results. Column 1 points to the result file, and column 2 a comma separated list of column numbers indicating control experiments. All other parameters in task 1 are supported" +
	"\n\t\tIt is assumed that all files have the same list of transcripts and thus that all use the same control transcripts" +
	"\n\t3. From a matrix of z-scores (or any scores that can be compared accross transcripts) compute new z-scores using a set of given transcripts -in <Matrix with scores, it must have a header with column names and the first column must contain row names, standard input is assume if this parameter is not specified> -scoreToWrite <0: z-scores 1: confidence> -controlTranscripts <Set of transcripts to use in order to compute new z-scores> -gp <If output should be written as a gene pattern gct format>" +
	"\n\t4. Join confidence matrices: Given to similarly sized confidence matrices combine the into a single one by taking logical operations. AND takes the minimum score, OR the maximum and MAX the average -matrix1 <A confidence matrix> -matrix2 <Another confidence matrix> [-op <[AND],OR,AVG>] " +
	"\n\t5. Filter a score matrix by a compatible confidence score matrix -in <Matrix to filter (standard input by default)> -confidenceMatrix <Matrix with confidnece scores> -minConfidence <Minimum confidence score to pass filter> [-gp <If output should be in gct format> -max <If not confidence but other measure that limits the lower end use this flag> -boolean <Convert passing scores to 1s and -1s> -biotapestry  <Outputs all passing scores a .csv that specifies interactions for knockdowns which caused a transcript to be significantly differentially epressed>]"+
	"\n\t6. Join to matrices sharing same rows by appending the sample data of the second to the first -matrix1 <First matrix to combine> -matrix2 <Second matrix to use> -gp <If output should be in gct format>"+
	"\n\t7. Write a matrix with header in BioTapestry importable format -in <Value matrix, standard input is supported too> -r <add this flag if negative values should be interpreted as inducers and positive values as repressors>"+
	"\n";
	
	public static void main (String [] args) throws IOException, ParseException {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE);
		
		if ("1".equals(argMap.getTask())) {	
			InputStream is = argMap.getInputStream();
			NanostringReader nr = new NanostringReader();
			nr.load(is);
			is.close();
			
			if(argMap.containsKey("controlExperiments") ) {
				List<String> controlExperiments = CLUtil.listFromArray(argMap.getMandatory("controlExperiments").split(","));
				int [] controlCols = new int [controlExperiments.size()];
				for (int i = 0; i < controlCols.length; i++){
					controlCols[i] = Integer.valueOf(Integer.valueOf(controlExperiments.get(i)));
				}
				nr.setControls(controlCols);
			}

			processData(argMap, nr);
			
		}	else if ("2".equals(argMap.getTask())) {
			NanostringReader nr = new NanostringReader();
			BufferedReader br = argMap.getInputReader();
			String line = null;
			while((line = br.readLine()) != null) {
				line.trim();
				int experimentNum = nr.getExperiments().size();
				String [] info = line.split("\t");
				String file = info[0];
				System.err.println("Processing " + file + " controls:  " + info[1]) ;
				String [] controls = info[1].split(",");
				int [] controlCols = new int [controls.length];
				for (int i = 0; i < controlCols.length; i++){
					int controlNum = Integer.valueOf(Integer.valueOf(controls[i])) + experimentNum;
					System.err.println("Adding control " + controlNum + " controlStr " + controls[i]);
					controlCols[i] = controlNum;
				}
				nr.load(file);
				System.err.println("Loaded.");
				nr.addControls(controlCols);
				System.err.println("Controls set");
			}
			
			br.close();

			processData(argMap, nr);
			
		} else if ("3".equals(argMap.getTask())) {
			List<String> controlTranscripts = CLUtil.listFromArray(argMap.getMandatory("controlTranscripts").split(","));
			BufferedReader br = argMap.getInputReader();
			MatrixWithHeaders mwh = null;
			boolean inGCTFormat = argMap.containsKey("gp");
			try {
				mwh = new MatrixWithHeaders(br);
			} finally {
				br.close();
			}
			
			if("1".equals(argMap.get("scoreToWrite"))) {
				System.err.println("computing confidence");
				mwh = mwh.confidenceByColumn(controlTranscripts);
			} else {
				mwh = mwh.zScoresByColumn(controlTranscripts);
			}
			
			BufferedWriter bw = argMap.getOutputWriter();
			if(inGCTFormat) {
				mwh.writeGCT(bw);
			} else {
				mwh.write(bw);
			}
			bw.close();

		}  else if ("4".equals(argMap.getTask())) {
			String mat1File = argMap.getMandatory("matrix1");
			String mat2File = argMap.getMandatory("matrix2");
			String op   = argMap.containsKey("op") ? argMap.getMandatory("op") : "AND";
			
			MatrixWithHeaders mat1 = new MatrixWithHeaders(mat1File);
			MatrixWithHeaders mat2 = new MatrixWithHeaders(mat2File);
			
			if(mat1.columnDimension() != mat2.columnDimension() ) {
				System.err.println("Matrices do not have the same column dimension ");
			}
			if(mat1.rowDimension() != mat2.rowDimension() ) {
				System.err.println("Matrices do not have the same row dimension ");
			}
			
			for (int i = 0; i < mat1.rowDimension(); i++) {
				for(int j =0; j < mat1.columnDimension(); j++) {
					double value = 0;
					if("OR".equals(op)) {
						value = Math.max(mat1.get(i,j), mat2.get(i, j));
					} else if("AVG".equals(op)) {
						value = (mat1.get(i, j) + mat2.get(i, j))/(double)2;
					} else {
						value = Math.min(mat1.get(i, j), mat2.get(i,j));
					}
					mat1.set(i, j, value);
				}
			}
			
			BufferedWriter bw = argMap.getOutputWriter();
			mat1.write(bw);
			bw.close();

		} else if ("5".equals(argMap.getTask())) {
			
			BufferedReader br = argMap.getInputReader();
			MatrixWithHeaders scoreMatrix = null;
			double minConfidence = (!argMap.containsKey("max") || argMap.containsKey("minConfidence") )? argMap.getDouble("minConfidence"): Double.MIN_VALUE;
			double max =  argMap.containsKey("max") ? argMap.getDouble("max"): Double.MAX_VALUE;
			boolean roundUp  = argMap.containsKey("boolean");
			boolean inGCTFormat = argMap.containsKey("gp");
			try {
				scoreMatrix = new MatrixWithHeaders(br);
			} finally {
				br.close();
			}
			String filterFile = argMap.getMandatory("confidenceMatrix");
			MatrixWithHeaders filter = new MatrixWithHeaders(filterFile);
			
			if(scoreMatrix.columnDimension() != filter.columnDimension() ) {
				System.err.println("Matrices do not have the same column dimension ");
			}
			if(scoreMatrix.rowDimension() != filter.rowDimension() ) {
				System.err.println("Matrices do not have the same row dimension ");
			}
			
			for (int i = 0; i < scoreMatrix.rowDimension(); i++) {
				for(int j =0; j < scoreMatrix.columnDimension(); j++) {
					if(filter.get(i, j) > minConfidence && filter.get(i,j) <= max) {
						if(roundUp) {
							scoreMatrix.set(i, j, scoreMatrix.get(i, j) > 0 ? 1 : -1);
						} else {
							scoreMatrix.set(i, j, scoreMatrix.get(i, j));
						}
					} else {
						scoreMatrix.set(i, j, 0);
					}
				}
			}
			
			BufferedWriter bw = argMap.getOutputWriter();
			if(inGCTFormat) {
				scoreMatrix.writeGCT(bw);
			} else {
				scoreMatrix.write(bw);
			}
			bw.close();

		}else if ("6".equals(argMap.getTask())) {			
			String mat1File = argMap.getMandatory("matrix1");
			String mat2File = argMap.getMandatory("matrix2");
			MatrixWithHeaders mat1 = new MatrixWithHeaders(mat1File);
			MatrixWithHeaders mat2 = new MatrixWithHeaders(mat2File);
			boolean inGCTFormat = argMap.containsKey("gp");

			if(mat1.rowDimension() != mat2.rowDimension() ) {
				System.err.println("Matrices do not have the same row dimension ");
			}
			mat1.append(mat2);
			
			BufferedWriter bw = argMap.getOutputWriter();
			if(inGCTFormat) {
				mat1.writeGCT(bw);
			} else {
				mat1.write(bw);
			}
			bw.close();

		} else if ("7".equals(argMap.getTask())){
			BufferedReader br = argMap.getInputReader();
			MatrixWithHeaders scoreMatrix = null;
			try {
				scoreMatrix = new MatrixWithHeaders(br);
			} finally {
				br.close();
			}
			
			BufferedWriter bw = argMap.getOutputWriter();
			try {
				writeInBiotapestryCSV(scoreMatrix, bw, argMap.containsKey("r"));
			} finally {
				br.close();
			}
		}	else {
			System.err.println(USAGE);
		}
	
	}

	private static void processData(ArgumentMap argMap, NanostringReader nr)throws IOException {
		
		float minFoldChange = argMap.containsKey("minFoldChange") ? argMap.getFloat("minFoldChange") : 1;
		double logMinFoldChange = MathUtil.log2(minFoldChange);
		if(argMap.containsKey("useAsReference")) {
			nr.setFirstExperiment(argMap.getInteger("useAsReference"));
		}
		nr.normalizeBySpikeIns();

		if(argMap.containsKey("controlTranscripts")) {
			List<String> controlTranscripts = CLUtil.listFromArray(argMap.getMandatory("controlTranscripts").split(","));
			nr.setControlTranscripts(controlTranscripts);
			nr.normalizeByNormalizationTranscripts();
		}
		
		MatrixWithHeaders normalizedData = nr.getCountMatrix();
		//BufferedWriter tbw = new BufferedWriter(new FileWriter(("normalized.txt")));
		//normalizedData.write(tbw);
		//tbw.close();
		
		normalizedData.log();
		double upRegulatedConfidence =  argMap.containsKey("minConfidence") ? argMap.getDouble("minConfidence") : 0;
		double downRegulatedConfidence = argMap.containsKey("minConfidence") ? argMap.getDouble("minConfidence") : 0;
		
		if(argMap.containsKey("upRegConfidence")) {
			upRegulatedConfidence = argMap.getDouble("upRegConfidence") ;
		}
		if(argMap.containsKey("downRegConfidence")) {
			downRegulatedConfidence = argMap.getDouble("downRegConfidence") ;
		}			
		int scoreToWrite = argMap.containsKey("scoreToWrite") ? argMap.getInteger("scoreToWrite"): 0;
		if(0 == scoreToWrite || nr.getControlExperimentNames().size() == 0) {
			BufferedWriter bw = argMap.getOutputWriter();
			normalizedData.pow();
			normalizedData.round();
			normalizedData.write(bw);
			bw.close();
		}else  {
			Map<String, Collection<String>> groups = setGroups(normalizedData, argMap);
			MatrixWithHeaders result = new MatrixWithHeaders(normalizedData.getRowNames(), new ArrayList<String>(groups.keySet()));
			for(String groupName : groups.keySet()) {

				Collection<String> controls =  new ArrayList<String>(nr.getControlExperimentNames());
				Collection<String> group    = groups.get(groupName);
				
				for(String groupMember : group) {
					if(controls.contains(groupMember) ) {
						controls.remove(groupMember);
					}
				}
				//System.err.println("Controls: " + controls + " exp " + exp);
				DifferentialExpression de = new DifferentialExpression(normalizedData, group, controls) ;
				//MatrixWithHeaders [] permMatrices = de.getPermutationMatrix();
				//for(int f = 0; f < permMatrices.length; f++) {
				//	permMatrices[f].write("tmp/"+exp+"_perm"+f+".txt");
				//}
				MatrixWithHeaders foldMatrix = DifferentialScoring.computeFoldMatrix(normalizedData, group, controls);
				for(String row : result.getRowNames()) {
					
					double fold = foldMatrix.get(row, 0);
					switch (scoreToWrite){
					case 1:
						result.set(row, groupName, Math.abs(fold) > logMinFoldChange ?  1 - de.getFDRMatrix().get(row, 0) : 0);
						break;
					case 2:
						result.set(row, groupName, Math.abs(fold) > logMinFoldChange ? de.getTestStatisticMatrix().get(row, 0) : 0);
						break;
					case 3:
						result.set(row, groupName, MathUtil.log2(fold));
						break;
					case 4:
						result.set(row, groupName,  Math.abs(fold) > logMinFoldChange ?  de.getFWERMatrix().get(row, 0) : 0);
						break;
					}
						
				}
				
			}
			BufferedWriter bw = argMap.getOutputWriter();
			result.write(bw);
			bw.close();
		}
		
		if(argMap.containsKey("cytoscape")) {
			nr.writeCytoscapeNetwork(argMap.containsKey("networkName") ? argMap.get("networkName") : "network", argMap.getAll("ignoreExperiment"), argMap.getAll("ignoreTranscript"), downRegulatedConfidence, upRegulatedConfidence);
		} else {
			OutputStream os = argMap.getOutputStream();
			int scoreToUse = 2;
			if(argMap.containsKey("scoreToWrite")){
				scoreToUse = argMap.getInteger("scoreToWrite");
			} else if (!argMap.containsKey("controlExperiments")) {
				scoreToUse = 0;
				
			}
			nr.writeSpikedInNormalized(os, scoreToUse, argMap.containsKey("gp"),downRegulatedConfidence, upRegulatedConfidence);
			os.close();
		}
	}

	private static Map<String, Collection<String>> setGroups(MatrixWithHeaders data, ArgumentMap argMap) throws IllegalArgumentException, IOException {
		LinkedHashMap<String, Collection<String>> groups = new LinkedHashMap<String, Collection<String>>();
		if(!argMap.containsKey("groupFile")) {
			
			for(String exp : data.getColumnNames()) {
				Collection<String> expGroup = new ArrayList<String>(1);
				expGroup.add(exp);
				groups.put(exp, expGroup);
			}
		} else {
			BufferedReader br = new BufferedReader(new FileReader(argMap.getMandatory("groupFile")));
			String line = null;
			while((line = br.readLine()) != null ) {
				String [] info = line.split("\t");
				String [] group = info[1].split(",\\s*");
				ArrayList<String> groupList = new ArrayList<String>(group.length);
				for(String g : group) {
					groupList.add(g);
				}
				groups.put(info[0], groupList);
			}
			br.close();
		}
		return groups;
	}

	private static void writeInBiotapestryCSV(MatrixWithHeaders scoreMatrix,BufferedWriter bw, boolean reversedSignConversion) throws IOException {
		bw.write("#Command Type,Model Name,Model Parent,,,,,,");
		bw.newLine();
		bw.write("model,root,,,,,,,");
		bw.newLine();
		bw.newLine();
		bw.write("#Standard Interactions,,,,,,,,");
		bw.newLine();
		bw.write("#Command Type,Model Name,Model Parent,,,,,,");
		bw.newLine(); 
		bw.write("#Command Type,Model Name,Source Type,Source Name,Target Type,Target Name,Sign,Source Region Abbrev,Target Region Abbrev");
		bw.newLine();
		
		List<String> columns = scoreMatrix.getColumnNames();
		List<String> rows    = scoreMatrix.getRowNames();
		for(String col : columns) {
			for(String row : rows) {
				double val = scoreMatrix.get(row, col);
				if(val != 0) {
					bw.write("general,root,gene,");
					bw.write(col);
					bw.write(",gene,");
					bw.write(row);
					bw.write(",");
					bw.write((val > 0 && !reversedSignConversion) ||  (val < 0 && reversedSignConversion) ? "positive" : "negative");
					bw.write(",,");
					bw.newLine();
				}
			}
		}
		bw.flush();
		
	}
	
}
