package broad.pda.chromosome;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import broad.core.annotation.AnnotationReader;
import broad.core.annotation.AnnotationReaderFactory;
import broad.core.annotation.BED;
import broad.core.annotation.BEDReader;
import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

public class ChromosomeUtils {
	public static String USAGE = "Usage: ChromosomeUtils TASK=<task_num> <task_args>\n" +
	"\tTasks:\n" +
	"\n\t\textract. Extract sequence for annotations: -in <AnnotatIon file> -seqdir <Top level sequence directory for genome sequence> [-format <[BED], GFF, genereic>] -maskSoftMask <Include this flag if softmasked (lower case) sequence should be hard masked (converted to N)> ] "+
	"\n\t\textractThickRegion Extracts the thick region of a BED gene annotation list \n\t\t\t-in <AnnotatIon file> \n\t\t\t-seqdir <Top level sequence directory for genome sequence> \n\t\t\t[-maskSoftMask <Include this flag if softmasked (lower case) sequence should be hard masked (converted to N) -minScore <Only extract peaks that pass minimum score> -pad <extract extra flanking sequence arond the motif>]"+
	"\n";
	public static void main(String [] args) throws Exception {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE, "extract");	
		
		if("extract".equals(argMap.getTask()) || argMap.getTask() == null) {
			String annotationFile = argMap.getInput();
			String annotationFormat = argMap.containsKey("format") ? argMap.get("format") : "BED";
			String seqDir = argMap.getMandatory("seqdir");
			boolean hardMask = argMap.containsKey("maskSoftMask");
			AnnotationReader<? extends GenomicAnnotation> annotations = AnnotationReaderFactory.create(annotationFile, annotationFormat);
			List<Sequence> extractedSeqs = extractAnnotationListSequence(seqDir, hardMask, annotations);
			
			BufferedWriter bw = argMap.getOutputWriter();
			FastaSequenceIO fsio = new FastaSequenceIO();
			fsio.write(extractedSeqs, bw);
			bw.close();
		} if("extractThickRegion".equals(argMap.getTask()) || argMap.getTask() == null) { 
			String annotationFile = argMap.getInput();
			String seqDir = argMap.getMandatory("seqdir");
			boolean hardMask = argMap.containsKey("maskSoftMask");
			BEDReader annotations = new BEDReader(annotationFile);
			Mammal m = new GenericOrganism(new File(seqDir));
			double minScore = argMap.containsKey("minScore") ? argMap.getDouble("minScore") : 0;
			int padding = argMap.containsKey("pad") ? argMap.getInteger("pad") : 0;
			
			List<Sequence> extractedSeqs = new ArrayList<Sequence>();
			Iterator<String> chrIt = annotations.getChromosomeIterator();
			while(chrIt.hasNext()) {
				String chr = chrIt.next();
				System.err.println("Chromosome " + chr);
				String chrSymbol = chr.replace("chr", "");
				Chromosome c = m.getChromosome(chrSymbol);
				System.err.println("Chromosome sequence loaded");
				c.loadSequence();
				List<BED> chrAnnotations = annotations.getChromosomeBEDs(chr);
				for(BED g : chrAnnotations) {
					if(g.getScore() < minScore) {
						continue;
					}
					System.err.println("Score passing annotation " + g);
					LightweightGenomicAnnotation tg = new BasicLightweightAnnotation(g);
					tg.setStart(Math.max(0, g.getThickStart() - padding));
					tg.setEnd(Math.min(c.getSize(),g.getThickEnd()+padding));
					SequenceRegion gReg = new SequenceRegion(chr, tg);
					c.extractRegion(gReg);
					if (hardMask){
						System.err.println("Hard masking");
						gReg.maskSoftmaskedRegions();
						System.err.println("Done hardmasking");
					}
					
					extractedSeqs.add(gReg);
					System.err.println("Done with annotaton");
				}
				c.unloadSequence();
				System.err.println("Chromosome sequence unloaded");
			}
			
			BufferedWriter bw = argMap.getOutputWriter();
			FastaSequenceIO fsio = new FastaSequenceIO();
			fsio.write(extractedSeqs, bw);
			bw.close();
		}else {
			System.err.println(USAGE);
		}

	}
	public static List<Sequence> extractAnnotationListSequence(String seqDir,boolean hardMask, AnnotationReader<? extends GenomicAnnotation> annotations) throws Exception, IOException {
		Mammal m = new GenericOrganism(new File(seqDir));
		
		List<Sequence> extractedSeqs = new ArrayList<Sequence>();
		Iterator<String> chrIt = annotations.getChromosomeIterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			System.err.println("Chromosome " + chr);
			String chrSymbol = chr.replace("chr", "");
			Chromosome c = m.getChromosome(chrSymbol);
			System.err.println("Chromosome sequence loaded");
			c.loadSequence();
			List<? extends GenomicAnnotation> chrAnnotations = annotations.getChromosomeBEDs(chr);
			for(GenomicAnnotation g : chrAnnotations) {
				System.err.println("Annotation " + g);
				SequenceRegion gReg = new SequenceRegion(chr, g);
				c.extractRegion(gReg);
				if (hardMask){
					System.err.println("Hard masking");
					gReg.maskSoftmaskedRegions();
					System.err.println("Done hardmasking");
				}
				
				extractedSeqs.add(gReg);
				System.err.println("Done with annotaton");
			}
			c.unloadSequence();
			System.err.println("Chromosome sequence unloaded");
		}
		return extractedSeqs;
	}

}
