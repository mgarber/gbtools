/**
 * 
 */
package broad.pda.rap;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Set;

import broad.core.annotation.AnnotationReader;
import broad.core.annotation.AnnotationReaderFactory;
import broad.core.annotation.GenomicAnnotation;
import broad.core.error.ParseException;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.core.util.PipelineUtils;
import broad.pda.annotation.BEDFileParser;

/**
 * @author engreitz
 * Wrapper class to submit a Scripture ChIP-Seq job across multiple bsub jobs,
 * then aggregate them after all have finished
 */
public class CallPeaks {

	/**
	 * 
	 */
	public CallPeaks(String scripture, String queue, ArgumentMap scriptureParams, Set<String> chromosomes) throws ParseException, IOException, InterruptedException {
		Runtime run=Runtime.getRuntime();
		String jobID=PipelineUtils.getJobID();
		
		String baseOut = scriptureParams.getOutput();
		for (String chr : chromosomes) {
			// Submit a Scripture job
			scriptureParams.put("out", getChrOutput(baseOut, chr));
			scriptureParams.put("chr", chr);
			String command = "-P RAP java -Xmx10g -jar " + scripture + " -task chip " + scriptureParams.toArgString();
			String output = scriptureParams.getOutput() + ".bsub";
			//System.out.println(command);
			//System.out.println(output);
			PipelineUtils.bsubProcess(run, jobID, command, output, queue);
			scriptureParams.remove("chr");  // not sure why this is necessary but it is
		}
		
		PipelineUtils.waitForJobs(jobID, run, false);
		
		// Collect and merge results
		//System.out.println("Merging results ...");
		String combinedResults = baseOut + ".bed";
		BufferedWriter writer = new BufferedWriter(new FileWriter(combinedResults));
		for (String chr : chromosomes) {
			AnnotationReader<? extends GenomicAnnotation> set = AnnotationReaderFactory.create(getChrOutput(baseOut, chr), "BED");
			set.merge();
			writer.write(set.toStringForWriting());
		}
		writer.close();
		
		// Score final segments
		scriptureParams.put("in", combinedResults);
		scriptureParams.put("out", combinedResults + ".scores");
		String command = "-P RAP java -jar " + scripture + " -task scoreSegments " + scriptureParams.toArgString();
		String output = scriptureParams.getOutput() + ".bsub";
		//System.out.println(command);
		//System.out.println(output);
		PipelineUtils.bsubProcess(run, jobID, command, output, queue);
		PipelineUtils.waitForJobs(jobID, run, false);
	}
	
	private String getChrOutput(String base, String chr) {
		return base + "." + chr + ".bed";
	}

	private static String USAGE = "java -jar CallPeaks.jar -queue <queue> -scripture <scripture.jar> -sizeFile <sizes> [...scripture arguments...]";
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws ParseException, IOException, InterruptedException {
		ArgumentMap argmap = CLUtil.getParameters(args, USAGE);
		
		String sizeFile = argmap.getMandatory("sizeFile");
		String queue = argmap.getMandatory("queue");
		argmap.remove("queue");
		String scripture = argmap.getMandatory("scripture");
		argmap.remove("scripture");
		
		Set<String> chromosomes=BEDFileParser.loadChrSizes(sizeFile).keySet();
		new CallPeaks(scripture, queue, argmap, chromosomes);
	}

}
