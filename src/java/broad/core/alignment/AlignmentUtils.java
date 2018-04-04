package broad.core.alignment;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import broad.core.annotation.AnnotationReader;
import broad.core.annotation.AnnotationReaderFactory;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.error.ParseException;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

public class AlignmentUtils {
	public static final String USAGE = "Usage: AlignmentUtils TASK=<task_num> <task_args>\n" +
	"\tTasks:\n" +
	"\t\t1. Clean alignment list from self alignments in a genome wide alignment, given an annotation file with the query locations and a alignment result file -alignment <Alignment file, only PH one line files are supported> -queryAnnotations <Annotation file with locations of aligned sequences in reference> [-queryFormat <format of queryAnntations: [BED],GFF,generic>" +
	"\n";
	
	public static void main (String [] args) throws IOException, ParseException {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE);
		
		if ("1".equals(argMap.getTask())) {	
			String alignment = argMap.getMandatory("alignment");
			String queryAnnotations = argMap.getMandatory("queryAnnotations");
			String queryFormat = argMap.containsKey("queryFormat") ? argMap.get("queryFormat") : "BED";
			
			AnnotationReader<? extends GenomicAnnotation> query = AnnotationReaderFactory.create(queryAnnotations, queryFormat);
			PHTabularReader alnReader = new PHTabularReader(alignment);
			
			Iterator<PHOneLineAlignment> alignmentIt = alnReader.getAlignmentList().iterator();
			BufferedWriter bw = argMap.getOutputWriter();
			while(alignmentIt.hasNext()) {
				AlignmentSummary aln = alignmentIt.next();
				LightweightGenomicAnnotation hit = aln.getA();
				//System.err.println("Hit: " + hit.getName() + " in chr" + hit.getChromosome());
				List<? extends GenomicAnnotation> overlappers = query.getOverlappers(hit);
				boolean print = true;
				if(overlappers.size() > 0) {
					Iterator<? extends GenomicAnnotation> overlappersIt = overlappers.iterator();
					while(overlappersIt.hasNext() && print) {
						GenomicAnnotation overlapper = overlappersIt.next();
						print = !aln.getB().getName().equals(overlapper.getName());
						//System.err.println("Found overlapper for " + aln.getBName() + ": " + overlapper.getName() + " skipping? " + !print);
					}
				}
				if(print) {
					bw.write(aln.toString());
					bw.newLine();
				}
			}
			bw.close();
			
		}  else {
			System.err.println("Task " + argMap.getTask() + " is invalid or no task was specified");
		}
	}
}
