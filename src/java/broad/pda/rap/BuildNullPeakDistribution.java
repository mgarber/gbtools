/**
 * 
 */
package broad.pda.rap;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.broad.igv.Globals;

import broad.core.error.ParseException;
import broad.core.math.EmpiricalDistribution;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.core.util.PipelineUtils;
import broad.pda.annotation.BEDFileParser;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;

/**
 * @author engreitz
 *
 */
public class BuildNullPeakDistribution {
	
	private Runtime run;
	private String jobID;
	private Collection<String> chromosomes;
	
	BuildNullPeakDistribution(String sizeFile) {
		run = Runtime.getRuntime();
		jobID = PipelineUtils.getJobID();
		chromosomes = BEDFileParser.loadChrSizes(sizeFile).keySet();
	}
	
	private static String USAGE = "java -jar RaptureBuildNullDistribution.jar [args] TODO";
		
	private static List<String> loadBamList(String file) throws IOException {
		Collection<String> lines=BEDFileParser.loadList(file, false);
		List<String> list = new ArrayList<String>();
		list.addAll(lines);
		return list;
	}
	
	private void submitJobs(String outdir, List<String> bamList, ArgumentMap argmap, String queue, String rapture) throws IOException, InterruptedException {		
		for (int i = 0; i < bamList.size(); i++) {
			for (int j = 0; j < bamList.size(); j++) {
				if (i == j) continue;
				//for (String chr : chromosomes) {
				//	if (chr.equals("chrM") || chr.equals("chrY")) continue;
					
					String bam1 = bamList.get(i);
					String bam2 = bamList.get(j);

					argmap.put("out", getResultFileName(outdir, bam1, bam2));
					argmap.put("target", bam1);
					argmap.put("control", bam2);
					//argmap.put("chr", chr);
					String command = "-M 4 -P RAP java -Xmx4g -jar " + rapture + " -task distribution " + argmap.toArgString();
					String output = argmap.getOutput() + ".bsub";
					PipelineUtils.bsubProcess(run, jobID, command, output, queue);

					argmap.remove("target");		// not sure why this is necessary but it is
					argmap.remove("control");		// not sure why this is necessary but it is
					argmap.remove("out");
					argmap.remove("chr");
				//}
			}
		}
		
		PipelineUtils.waitForJobs(jobID, run, false);
	}
	
	
	private static String getResultFileName(String outdir, String file1, String file2) {
		String name1 = new File(file1).getName();
		String name2 = new File(file2).getName();
		return outdir + "/" + name1 + "_" + name2;
	}
	
	
	private void combineEmpiricalDistributions(String outdir, List<String> bamList, int[] windows) throws IOException {
		for (int window : windows) {
			EmpiricalDistribution ed = Rapture.getEmptyEmpiricalDistribution();
			for (int i = 0; i < bamList.size(); i++) {
				for (int j = 0; j < bamList.size(); j++) {
					if (i == j) continue;
					//for (String chr : chromosomes) {
					//	if (chr.equals("chrM") || chr.equals("chrY")) continue;
					String filename = Rapture.getEdfFileName(getResultFileName(outdir, bamList.get(i), bamList.get(j)), window);
					EmpiricalDistribution currEd = new EmpiricalDistribution(new File(filename));
					ed.addDistribution(currEd);
					//}
				}
			}
			ed.write(outdir + "/NullPeakDistribution." + window + "_window.edf");
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException, ParseException, InterruptedException {
		Globals.setHeadless(true);
		ArgumentMap argmap = CLUtil.getParameters(args, USAGE);
		
		String outdir = argmap.getOutputDir();
		argmap.remove("outdir");
		String bamlistfile = argmap.getMandatory("bamlist");
		argmap.remove("bamlist");
		List<String> bamList = loadBamList(bamlistfile);
		
		String queue = argmap.getMandatory("queue");
		argmap.remove("queue");
		String rapture = argmap.getMandatory("rapture");
		argmap.remove("rapture");
		
		int[] windows = ContinuousDataAlignmentModel.getWidths(argmap.getMandatory("windows"));
		
		String sizeFile = argmap.getMandatory("sizeFile");
		BuildNullPeakDistribution builder = new BuildNullPeakDistribution(sizeFile);
		if (!argmap.containsKey("combineOnly")) builder.submitJobs(outdir, bamList, argmap, queue, rapture);
		builder.combineEmpiricalDistributions(outdir, bamList, windows);
		
	}
	
}
