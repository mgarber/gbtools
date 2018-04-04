package broad.core.gene;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import broad.core.annotation.AnnotationFactory;
import broad.core.annotation.AnnotationFactoryFactory;
import broad.core.annotation.AnnotationReader;
import broad.core.annotation.AnnotationReaderFactory;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.error.ParseException;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

public class GeneUtils {
	
	public static final String USAGE = "Usage: GeneUtils TASK=<task_num> <task_args>\n" +
	"\tTasks:\n" +
	"\n\t1. Take the intersection of an annotation set with a gene set: -annotations <File containing the annotation set>  -geneSet <File gene annotation set> [-annotationFormat <[BED],GFF,GENERIC> -geneSetFormat <[REFSEQ], (GTF comming soon)> -use <[CODING_EXONS],EXONS,INTRONS,PROMOTER,FULL_LOCI> -promoterBases <if using promoters you may specify its size, default is 500 bases from start of first exon> -hg17< if using older hg17 refseq annotations add this flag>]" +
	"\n\t2. Take the difference (all nonoverlapping elements) of an annotation set with a gene set: -annotations <File containing the annotation set>  -geneSet <File gene annotation set> [-annotationFormat <[BED],GFF,GENERIC> -geneSetFormat <[REFSEQ], (GTF comming soon)> -use <[CODING_EXONS],EXONS,INTRONS,PROMOTER,FULL_LOCI> -promoterBases <if using promoters you may specify its size, default is 500 bases from start of first exon> -hg17< if using older hg17 refseq annotations add this flag>]" +
	"\n\t3. Find top scoring hits in gene regions [-topHits <# of hits to retun per region default is 1> -annotationFormat <[BED],GFF,GENERIC> -geneSetFormat <[REFSEQ], (GTF comming soon) -use <[CODING_EXONS],EXONS,INTRONS,PROMOTER,FULL_LOCI> -promoterBases <if using promoters you may specify its size, default is 500 bases from start of first exon> -hg17< if using older hg17 refseq annotations add this flag>]" +
	"\n\t4. Get genes overlapping annotations -annotations <File containing the annotation set>  -geneSet <File gene annotation set> [-annotationFormat <[BED],GFF,GENERIC> -geneSetFormat <[REFSEQ], (GTF comming soon)>]" +
	"\n\t5. Convert GFF to full BED -in <Gene list in GFF format or standard input> -out <were to write list. Defaults to standard out>" +
	"\n";
	
	private static final Map<String, String> codonToSymbolMap = new LinkedHashMap<String, String>();
	
	

	public static void main (String [] args) throws IOException, ParseException {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE);
		String geneAnnotationVersion = argMap.containsKey("hg17") ? "HG17" : "HG18";
		String annotFormat = argMap.containsKey("annotationFormat") ? argMap.get("annotationFormat") : "BED";
		String use = argMap.containsKey("use")? argMap.get("use") : "EXONS";
		int promoterBases = argMap.containsKey("promoterBases") ? argMap.getInteger("promoterBases") : 500;
		
		//System.err.println("Using " + geneAnnotationsToUse.size() + " gene annotations");
		
		if ("1".equals(argMap.getTask())) {
			String annotations = argMap.getMandatory("annotations");
			AnnotationReader<? extends GenomicAnnotation> ar = AnnotationReaderFactory.create(annotations, annotFormat);
			String genes = argMap.getMandatory("geneSet");
			RefSeqReader rsr = new RefSeqReader(genes, "CODING_EXONS".equals(use), geneAnnotationVersion);
			if(!use.contains("PROMOTER")) {
				rsr.filterByOverlap(ar.getAnnotationList()); // Avoid dealing with full gene tree, specially in the case were we work with exons.
			}
			List<? extends GenomicAnnotation> geneAnnotationsToUse = 
				getAnnotationsToUse(use, promoterBases, rsr);
			ar.intersect(geneAnnotationsToUse);
			printResults(argMap, ar);
		} 	else if ("2".equals(argMap.getTask())) {	
			String annotations = argMap.getMandatory("annotations");
			AnnotationReader<? extends GenomicAnnotation> ar = AnnotationReaderFactory.create(annotations, annotFormat);
			String genes = argMap.getMandatory("geneSet");
			RefSeqReader rsr = new RefSeqReader(genes, "CODING_EXONS".equals(use), geneAnnotationVersion);
			if(!use.contains("PROMOTER")) {
				rsr.filterByOverlap(ar.getAnnotationList()); // Avoid dealing with full gene tree, specially in the case were we work with exons.
			}
			List<? extends GenomicAnnotation> geneAnnotationsToUse = 
				getAnnotationsToUse(use, promoterBases, rsr);
			ar.minus(geneAnnotationsToUse);
			printResults(argMap, ar);
		}  	else if ("3".equals(argMap.getTask())) {	
			String annotations = argMap.getMandatory("annotations");
			AnnotationReader<? extends GenomicAnnotation> ar = AnnotationReaderFactory.create(annotations, annotFormat);
			String genes = argMap.getMandatory("geneSet");
			RefSeqReader rsr = new RefSeqReader(genes, "CODING_EXONS".equals(use), geneAnnotationVersion);
			List<? extends GenomicAnnotation> geneAnnotationsToUse = 
				getAnnotationsToUse(use, promoterBases, rsr);
			int hits = argMap.containsKey("topHits") ? argMap.getInteger("topHits") : 1;
			Iterator<? extends GenomicAnnotation> geneRegIt = geneAnnotationsToUse.iterator();
			BufferedWriter bw = argMap.getOutputWriter();
			while(geneRegIt.hasNext()) {
				GenomicAnnotation geneRegion = geneRegIt.next();
				List<? extends GenomicAnnotation> overlappers = ar.getOverlappers(geneRegion);
				if(overlappers != null) {
					Collections.sort(overlappers, new Comparator<GenomicAnnotation>() {
	
						public int compare(GenomicAnnotation arg0, GenomicAnnotation arg1){
							return (int) ((arg1.getScore() - arg0.getScore()) * 1000);
						}
						
					});
					for(int i = 0; i < hits; i++) {
						if(overlappers.size() > i) {
							AnnotationFactory<? extends GenomicAnnotation> annotFact = AnnotationFactoryFactory.getFactory(annotFormat);
							
							LightweightGenomicAnnotation hit = annotFact.create(overlappers.get(i));
							hit.setName(geneRegion.getName() + "_"+ i);
							bw.write(hit.toString());
							bw.newLine();
						}
					}
				}
			}
			bw.close();
		} 	else if ("4".equals(argMap.getTask())) {	
			String annotations = argMap.getMandatory("annotations");
			AnnotationReader<? extends GenomicAnnotation> ar = AnnotationReaderFactory.create(annotations, annotFormat);
			String genes = argMap.getMandatory("geneSet");
			RefSeqReader rsr = new RefSeqReader(genes, "CODING_EXONS".equals(use), geneAnnotationVersion);
			Iterator<GeneAnnotation> geneAnnotationIt = rsr.getAnnotationList().iterator();
			BufferedWriter bw = argMap.getOutputWriter();
			while(geneAnnotationIt.hasNext()) {
				GenomicAnnotation gene = geneAnnotationIt.next();
				List<? extends GenomicAnnotation> overlappers = ar.getOverlappers(gene);
				if(overlappers.size() > 0) {
					Iterator<? extends GenomicAnnotation> overlapperIt = overlappers.iterator();
					while(overlapperIt.hasNext()) {
						GenomicAnnotation overlapper = overlapperIt.next();
						gene.setId(overlapper.getName());
						bw.write(gene.toString());
						bw.newLine();
					}
				}
			}
			bw.close();
				
		} 	else if ("5".equals(argMap.getTask())) {	
			BufferedReader br = argMap.getInputReader();
			GFFGeneReader reader = new GFFGeneReader();
			reader.load(br);
			br.close();
			
			List<GeneAnnotation> geneList = reader.getAnnotationList();
			BufferedWriter bw = argMap.getOutputWriter();
			if(reader.getTrackInfo(0) != null) {
				bw.write(reader.getTrackInfo(0));
				bw.newLine();
			}
			for(GeneAnnotation gene : geneList) {
				bw.write(gene.writeAsFullBed());
				bw.newLine(); 
			}
			bw.close();
				
		}else {
			System.err.println("Task " + argMap.getTask() + " is invalid or no task was specified");
		}
}

	private static void printResults(ArgumentMap argMap,
			AnnotationReader<? extends GenomicAnnotation> ar)
			throws IOException {
		Iterator<String> chrIt = ar.getChromosomeAnnotationMap().keySet().iterator();
		BufferedWriter bw = argMap.getOutputWriter();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			IntervalTree<? extends GenomicAnnotation> tree = ar.getChromosomeTree(chr);
			Iterator<? extends GenomicAnnotation> annotIt = tree.valueIterator();
			//System.out.println("chromsome: chr" + chr + " tree size " + tree.size() + " does iterator has a value " + annotIt.hasNext());
			while(annotIt.hasNext()) {
				LightweightGenomicAnnotation annot = annotIt.next();
				bw.write(annot.toString());
				bw.newLine();
			}
		}
		bw.close();
	}

	private static List<? extends GenomicAnnotation> getAnnotationsToUse(
			String use, int promoterBases, RefSeqReader rsr) {
		List<? extends GenomicAnnotation> geneAnnotationsToUse = new ArrayList<GenomicAnnotation>();
		if("EXONS".equals(use) || "CODING_EXONS".equals(use)) {
			geneAnnotationsToUse = rsr.getAllExons(); 
		} else if("INTRONS".equals(use)) { 
			geneAnnotationsToUse = rsr.getAllIntrons(); 
		} else if("PROMOTERS".equals(use) || "PROMOTER".equals(use)) {
			geneAnnotationsToUse = rsr.getAllPromoters(promoterBases);
		} else if("FULL_LOCI".equals(use)) {
			geneAnnotationsToUse = rsr.getAnnotationList();
		} 
		return geneAnnotationsToUse;
	}
	
	public static Map<String, String> getCodonSymbolMap() { return codonToSymbolMap;};
	
	static {
		// *T*
		codonToSymbolMap.put("TTT", "F");
		codonToSymbolMap.put("TTC", "F");
		codonToSymbolMap.put("TTA", "L");
		codonToSymbolMap.put("TTG", "L");
		
		codonToSymbolMap.put("CTT", "L");
		codonToSymbolMap.put("CTC", "L");
		codonToSymbolMap.put("CTA", "L");
		codonToSymbolMap.put("CTG", "L");
		
		codonToSymbolMap.put("ATT", "I");
		codonToSymbolMap.put("ATC", "I");
		codonToSymbolMap.put("ATA", "I");
		codonToSymbolMap.put("ATG", "M");
		
		codonToSymbolMap.put("GTT", "V");
		codonToSymbolMap.put("GTC", "V");
		codonToSymbolMap.put("GTA", "V");
		codonToSymbolMap.put("GTG", "V");

		// *C*
		codonToSymbolMap.put("TCT", "S");
		codonToSymbolMap.put("TCC", "S");
		codonToSymbolMap.put("TCA", "S");
		codonToSymbolMap.put("TCG", "S");
		
		codonToSymbolMap.put("CCT", "P");
		codonToSymbolMap.put("CCC", "P");
		codonToSymbolMap.put("CCA", "P");
		codonToSymbolMap.put("CCG", "P");
		
		codonToSymbolMap.put("ACT", "T");
		codonToSymbolMap.put("ACC", "T");
		codonToSymbolMap.put("ACA", "T");
		codonToSymbolMap.put("ACG", "T");
		
		codonToSymbolMap.put("GCT", "A");
		codonToSymbolMap.put("GCC", "A");
		codonToSymbolMap.put("GCA", "A");
		codonToSymbolMap.put("GCG", "A");
		
		// *A*
		codonToSymbolMap.put("TAT", "Y");
		codonToSymbolMap.put("TAC", "Y");
		codonToSymbolMap.put("TAA", "St");
		codonToSymbolMap.put("TAG", "St");
		
		codonToSymbolMap.put("CAT", "H");
		codonToSymbolMap.put("CAC", "H");
		codonToSymbolMap.put("CAA", "K");
		codonToSymbolMap.put("CAG", "K");
		
		codonToSymbolMap.put("AAT", "N");
		codonToSymbolMap.put("AAC", "N");
		codonToSymbolMap.put("AAA", "K");
		codonToSymbolMap.put("AAG", "K");
		
		codonToSymbolMap.put("GAT", "D");
		codonToSymbolMap.put("GAC", "D");
		codonToSymbolMap.put("GAA", "E");
		codonToSymbolMap.put("GAG", "E");
		
		// *G*
		codonToSymbolMap.put("TGT", "C");
		codonToSymbolMap.put("TGC", "C");
		codonToSymbolMap.put("TGA", "St");
		codonToSymbolMap.put("TGG", "W");
		
		codonToSymbolMap.put("CGT", "R");
		codonToSymbolMap.put("CGC", "R");
		codonToSymbolMap.put("CGA", "R");
		codonToSymbolMap.put("CGG", "R");
		
		codonToSymbolMap.put("AGT", "S");
		codonToSymbolMap.put("AGC", "S");
		codonToSymbolMap.put("AGA", "R");
		codonToSymbolMap.put("AGG", "R");
		
		codonToSymbolMap.put("GGT", "G");
		codonToSymbolMap.put("GGC", "G");
		codonToSymbolMap.put("GGA", "G");
		codonToSymbolMap.put("GGG", "G");
	}
}
