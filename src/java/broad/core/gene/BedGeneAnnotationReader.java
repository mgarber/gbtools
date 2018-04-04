package broad.core.gene;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;

import broad.core.annotation.AnnotationHandler;
import broad.core.annotation.AnnotationReader;
import broad.core.annotation.BED;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.GenomicAnnotationFilter;
import broad.core.error.ParseException;
import broad.core.gene.GeneAnnotation.Exon;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

public class BedGeneAnnotationReader extends AnnotationReader<GeneAnnotation>{
	static final BEDGeneFactory factory = new BEDGeneFactory();
	public BedGeneAnnotationReader() {
		super();
	}
	
	public BedGeneAnnotationReader(BufferedReader br) throws ParseException, IOException {
		super();
		load(br, factory);
	}
	
	public BedGeneAnnotationReader(String input) throws ParseException, IOException {
		super();
		load(new File(input), factory);
	}
	
	public void load(BufferedReader br) throws IOException, ParseException {
		super.load(br, factory);
	}

	@Override	public GeneAnnotation createAnnotation(GenomicAnnotation a) {
		return null;
	}
	
	public static final String USAGE = "Usage: TreeScaler TASK=<task_num> <task_args>\n" +
	"\tTasks:\n" +
	"\t\t1. Get Exon bed file from original annotation file -in <Annotation file (in full BED format) standard input is supported> -out <Output file, standard out is also supported>" +
	"\n\t\t2. Get Gene promters from original annotation file -in <Annotation file (in full BED format) standard input is supported> -out <Output file, standard out is also supported> -upstreamTSS <Bases upstream of TSS> -downstreamTSS <Bases downstream of the TSS>" +
	"\n";
	public static void main (String [] args) throws IOException, ParseException {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE);
		
		if ("1".equals(argMap.getTask())) {	
			BufferedReader br = argMap.getInputReader();
			BedGeneAnnotationReader bgar = new BedGeneAnnotationReader(br);
			br.close();
			
			List<GeneAnnotation> genes = bgar.getAnnotationList();
			BufferedWriter bw = argMap.getOutputWriter();
			for (GeneAnnotation gene : genes) {
				List<Exon> exons = gene.getExons();
				for(Exon exon : exons) {
					bw.write(exon.getChromosomeString());
					bw.write("\t");
					bw.write(String.valueOf(exon.getStart()));
					bw.write("\t");
					bw.write(String.valueOf(exon.getEnd()));
					bw.write("\t");
					bw.write(exon.getName());
					bw.newLine();
				}
			}
			bw.close();
		} else if ("2".equals(argMap.getTask())) {
			int upTSS = argMap.getInteger("upstreamTSS");
			int downTSS = argMap.getInteger("downstreamTSS");
			BufferedReader br = argMap.getInputReader();
			BedGeneAnnotationReader bgar = new BedGeneAnnotationReader(br);
			br.close();
								
			List<GeneAnnotation> genes = bgar.getAnnotationList();
			BufferedWriter bw = argMap.getOutputWriter();
			for (GeneAnnotation gene : genes) {
				int start = 0;
				int end = 0;
				if(gene.inReversedOrientation()) {
					int tss = gene.getEnd();
					start = tss - downTSS;
					end = tss + upTSS;
					
				} else {
					int tss = gene.getStart(); 
					start = tss - upTSS;
					end = tss - downTSS;
				}
				BED bed = new BED(gene);
				bed.setStart(start);
				bed.setEnd(end);
				bw.write(bed.toShortString());
				bw.newLine();
			}
			bw.close();
		}
		
	} 

	@Override
	public int parse(String file,
			GenomicAnnotationFilter<GeneAnnotation> filter,
			AnnotationHandler handler) throws ParseException, IOException {
		// TODO Auto-generated method stub
		return super.parse(file,factory, filter, handler);
	}

	@Override
	public int parse(BufferedReader br,
			GenomicAnnotationFilter<GeneAnnotation> filter,
			AnnotationHandler handler) throws ParseException, IOException {
		// TODO Auto-generated method stub
		return super.parse(br,factory, filter, handler);
	}

}
