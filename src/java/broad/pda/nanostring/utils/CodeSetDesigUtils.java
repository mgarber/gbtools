package broad.pda.nanostring.utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import broad.core.annotation.AnnotationReader;
import broad.core.annotation.AnnotationReaderFactory;
import broad.core.annotation.BED;
import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.sequence.Sequence;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.chromosome.Chromosome;
import broad.pda.chromosome.GenericOrganism;
import broad.pda.chromosome.Mammal;

public class CodeSetDesigUtils {
	
	public static final String USAGE = "Usage: CodeSetDesigUtils TASK=<task name> <task args>\n" +
	"\tTasks:\n" +
	"\t\taround_peaks Takes a list of sharp peaks, a list of K4 regions and list of gene locii and will design a codeset and an alignment file reports back the highest peak within each region. \n\t\t\t-peaks <peak file> \n\t\t\t-K4 <K4 region file> \n\t\t\t-genes <gene locii file> \n\t\t\t-spacing <Probe spacing>" +
	"\n\tadd_sequence Adds sequence to a bed file with probe design -in <probes set in BED file format> -seqdir <Sequence directory> " +
	"\n";
	
	static final int PROMOTER_RADIUS = 500;
	public static void main (String [] args) throws Exception {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE, "around_peaks");
		
		if ("around_peaks".equals(argMap.getTask())) {
			String geneFile = argMap.getMandatory("genes");
			String peakFile = argMap.getMandatory("peaks");
			String k4File   = argMap.getMandatory("K4");
			int spacing     = argMap.getInteger("spacing");
			String geneInfoFile =  argMap.getMandatory("infoFile") ;
			int maxExtraPeaksToConsider = 2;
			int minUpstreamScore = 40;
			int upstreamDistanceToLook = 30000;
			
			Map<String, String> geneInfo = loadGeneInfo(geneInfoFile);
			

			AnnotationReader<? extends GenomicAnnotation> genes = AnnotationReaderFactory.create(geneFile, "BED");
			AnnotationReader<? extends GenomicAnnotation> peaks = AnnotationReaderFactory.create(peakFile, "BED");
			AnnotationReader<? extends GenomicAnnotation> k4s = AnnotationReaderFactory.create(k4File, "BED");
			
			Map<String, Collection<LightweightGenomicAnnotation>> probes = new LinkedHashMap<String, Collection<LightweightGenomicAnnotation>>();
			Iterator<String> chrIt = genes.getChromosomeIterator();
			
			while(chrIt.hasNext()) {
				String chr = chrIt.next();
				System.err.println("Processing "+chr);
				List<? extends GenomicAnnotation> chrGenes = genes.getChromosomeBEDs(chr);
				Collection<LightweightGenomicAnnotation> chrProbes = new TreeSet<LightweightGenomicAnnotation>();
				probes.put(chr, chrProbes);
				for(GenomicAnnotation gene : chrGenes) {
					if("i".equals(geneInfo.get(gene.getName()))) {
						continue;
					}
					final int tss =  gene.inReversedOrientation() ? gene.getEnd() : gene.getStart();
					List<LightweightGenomicAnnotation> geneProbes =  new ArrayList<LightweightGenomicAnnotation>();
					LightweightGenomicAnnotation promoter = new BasicLightweightAnnotation(chr,tss - PROMOTER_RADIUS, tss + PROMOTER_RADIUS);
					List<? extends LightweightGenomicAnnotation> overlappingPeaks = peaks.getOverlappers(promoter);
					//First design 4 probes around TSS or if available best peak near TSS.
					int probeCenter = tss;
					double highScoringProbe = 0;
					if(!overlappingPeaks.isEmpty()) {
						//System.err.println(gene.getName() + " HAD OVERLAPPING PEAKS AROUND TSS");
						Collections.sort(overlappingPeaks, new Comparator<LightweightGenomicAnnotation>(){

							public int compare(LightweightGenomicAnnotation o1,LightweightGenomicAnnotation o2) {
								//return (int) ( Math.abs(tss - (o1.getMiddle())) - Math.abs(tss - (o2.getMiddle()))); //sort by proximity to TSS
								return (int) (o2.getScore() - o1.getScore()); // sort by score 
							}
							
						});
						//for(int i = 0; i < overlappingPeaks.size(); i++) {
						//	System.err.println("PEAK"+new BED(overlappingPeaks.get(i)).toShortString());
						//}
						probeCenter = (int)overlappingPeaks.get(0).getMiddle();
						highScoringProbe = overlappingPeaks.get(0).getScore();
					} 
					
					int probesToDesign = 4;
					if(!"c".equals(geneInfo.get(gene.getName()))) {
						probeCenter = gene.inReversedOrientation() ? probeCenter - spacing : probeCenter + spacing;
					} else {
						probeCenter = gene.inReversedOrientation() ? probeCenter - spacing/2 : probeCenter + spacing/2;
						probesToDesign = 2;
					}
					int probeNum = 1;
					while(probeNum <= probesToDesign) {
						LightweightGenomicAnnotation probe = new BasicLightweightAnnotation(chr, probeCenter - (spacing/2) , probeCenter + (spacing/2));
						probe.setOrientation(gene.getOrientation());
						probe.setName(gene.getName()+"_probe_"+ probeNum);
						probe.setScore(highScoringProbe);
						geneProbes.add(probe);
						probeCenter = gene.inReversedOrientation() ? probeCenter + spacing :  probeCenter - spacing;
						if((gene.inReversedOrientation() && probeCenter < gene.getStart()) || (!gene.inReversedOrientation() && probeCenter > gene.getEnd())) {
							break;
						}
						probeNum++;
					}
					
					// Now look for 2 probes for next K4 if any.
					int downStreamMostTile = gene.inReversedOrientation() ? geneProbes.get(0).getStart()  - spacing/2 : geneProbes.get(0).getEnd() + spacing/2 ;
					int start =  gene.inReversedOrientation() ? gene.getStart() : downStreamMostTile;
					int end   =  gene.inReversedOrientation() ? downStreamMostTile : gene.getEnd();
					
					if(end <= start) {
						continue;
					}
					
					LightweightGenomicAnnotation restOfGene = new BasicLightweightAnnotation(chr, start, end);
					restOfGene.setOrientation(gene.getOrientation());
					
					List<? extends LightweightGenomicAnnotation> otherCandidatePeaks = peaks.getOverlappers(restOfGene);
					List<LightweightGenomicAnnotation> highestScoringNonPromterK4PolIIpeaks = getRegionCandidatePeaks(k4s, Math.min(minUpstreamScore, highScoringProbe), otherCandidatePeaks);
				
					int peakNum = 0;
					while(peakNum < maxExtraPeaksToConsider && peakNum < highestScoringNonPromterK4PolIIpeaks.size() && highestScoringNonPromterK4PolIIpeaks.get(peakNum).getScore()> highScoringProbe/2) {
						int midPoint = (int)highestScoringNonPromterK4PolIIpeaks.get(peakNum).getMiddle();
						highScoringProbe = Math.max(maxExtraPeaksToConsider, highestScoringNonPromterK4PolIIpeaks.get(peakNum).getScore());
						peakNum++;
						LightweightGenomicAnnotation probe = new BasicLightweightAnnotation(chr, midPoint - spacing , midPoint);
						probe.setOrientation(gene.getOrientation());
						probe.setName(gene.getName()+"_probe_"+ probeNum);
						probe.setScore(highScoringProbe);
						if(!probe.overlaps(geneProbes)) {
							geneProbes.add(probe);
							probeNum++;
						}
						
						LightweightGenomicAnnotation probe2 = new BasicLightweightAnnotation(chr,midPoint , midPoint + spacing);
						probe2.setOrientation(gene.getOrientation());
						probe2.setName(gene.getName()+"_probe_"+ probeNum);
						probe2.setScore(highScoringProbe);
						if(!probe2.overlaps(geneProbes)) {
							geneProbes.add(probe2);
							probeNum++;
						}
						
					}

					//Now we go 30K upstream of the gene to find other good regions to probe.
					int upstreamMostCovered = gene.inReversedOrientation() ? geneProbes.get(0).getStart() + probesToDesign*spacing : geneProbes.get(0).getEnd() - probesToDesign*spacing;
					int upStart = gene.inReversedOrientation() ? upstreamMostCovered : upstreamMostCovered - upstreamDistanceToLook;
					int upEnd   = gene.inReversedOrientation() ? upstreamMostCovered + upstreamDistanceToLook: upstreamMostCovered;
					LightweightGenomicAnnotation upReg = new BasicLightweightAnnotation(gene.getChromosome(), upStart, upEnd);
					upReg.setOrientation(gene.getOrientation());
					
					List<? extends LightweightGenomicAnnotation> upCandidatePeaks = peaks.getOverlappers(upReg);
					List<LightweightGenomicAnnotation> upCandidatePolIIs = getRegionCandidatePeaks(k4s, minUpstreamScore, upCandidatePeaks);
					peakNum = 0;
					probeNum = 1;
					while(peakNum < maxExtraPeaksToConsider && peakNum < upCandidatePolIIs.size() && upCandidatePolIIs.get(peakNum).getScore()> minUpstreamScore) {
						int midPoint = (int)upCandidatePolIIs.get(peakNum).getMiddle();
						highScoringProbe = Math.max(maxExtraPeaksToConsider, upCandidatePolIIs.get(peakNum).getScore());
						peakNum++;
						LightweightGenomicAnnotation probe = new BasicLightweightAnnotation(chr, midPoint - spacing , midPoint);
						probe.setOrientation(gene.getOrientation());
						probe.setName(gene.getName()+"_enhc_probe_"+ probeNum);
						probe.setScore(highScoringProbe);
						if(!probe.overlaps(geneProbes)) {
							geneProbes.add(probe);
							probeNum++;
						}
						
						LightweightGenomicAnnotation probe2 = new BasicLightweightAnnotation(chr,midPoint , midPoint + spacing);
						probe2.setOrientation(gene.getOrientation());
						probe2.setName(gene.getName()+"_enhc_probe_"+ probeNum);
						probe2.setScore(highScoringProbe);
						if(!probe2.overlaps(geneProbes)) {
							geneProbes.add(probe2);
							probeNum++;
						}
						
					}
					
					chrProbes.addAll(geneProbes);
				}
			}
			BufferedWriter bw = argMap.getOutputWriter();
			write(bw, probes);
			bw.close();
		} else if ("add_sequence".equals(argMap.getTask())){
			String bedFile = argMap.getInput();
			String seqdir  = argMap.getMandatory("seqdir");
			
			AnnotationReader<? extends GenomicAnnotation> probeReader = AnnotationReaderFactory.create(bedFile, "BED");
			Mammal m = new GenericOrganism(new File(seqdir));
			Iterator<String> chrIt = probeReader.getChromosomeIterator();
			
			BufferedWriter bw = argMap.getOutputWriter();
			while(chrIt.hasNext()) {
				String chr = chrIt.next();
				Chromosome c = m.getChromosome(chr.replace("chr",""));
				c.loadSequence();
				List<? extends LightweightGenomicAnnotation> chrProbes = probeReader.getChromosomeBEDs(chr);
				for(LightweightGenomicAnnotation probe : chrProbes) {
					Sequence probeSeq = c.getSequence().getSubSequence(probe.getName(), probe.getStart(), probe.getEnd());
					bw.write(new BED(probe).toShortString());
					bw.write("\t");
					bw.write(probeSeq.getSequenceBases());
					bw.newLine();
				}
				c.unloadSequence();
			}
			bw.close();
		}
	}

	private static List<LightweightGenomicAnnotation> getRegionCandidatePeaks(
			AnnotationReader<? extends GenomicAnnotation> k4s,
			double highScoringProbe,
			List<? extends LightweightGenomicAnnotation> otherCandidatePeaks) {
		List<LightweightGenomicAnnotation> highestScoringNonPromterK4PolIIpeaks = new ArrayList<LightweightGenomicAnnotation>();
		for(LightweightGenomicAnnotation polIIPeak : otherCandidatePeaks) {
			List<? extends LightweightGenomicAnnotation> peakK4s = k4s.getOverlappers(polIIPeak);
			if(peakK4s.size() > 0 && polIIPeak.getScore() > highScoringProbe) {
				highestScoringNonPromterK4PolIIpeaks.add(polIIPeak);
			}
		}
		
		Collections.sort(highestScoringNonPromterK4PolIIpeaks, new Comparator<LightweightGenomicAnnotation>(){

			public int compare(LightweightGenomicAnnotation arg0, LightweightGenomicAnnotation arg1) {
				return (int) (arg1.getScore()-arg0.getScore());
			}
			
		});
		return highestScoringNonPromterK4PolIIpeaks;
	}
	
	private static Map<String, String> loadGeneInfo(String geneInfoFile) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(geneInfoFile)); 
		Map<String, String> rtrn = new HashMap<String, String>();
		try {
			br.readLine();
			String line = null;
			while((line = br.readLine())!=null) {
				String [] info = line.split("\t");
				rtrn.put(info[0], info[1]);
			}
		}finally {
			br.close();
		}
		return rtrn;
	}

	public static void write(BufferedWriter bw, Map<String, Collection<LightweightGenomicAnnotation>> peaks) throws IOException {
		for(String chr : peaks.keySet()) {
			for(LightweightGenomicAnnotation peak : peaks.get(chr)) {
				bw.write(new BED(peak).toShortString());
				bw.newLine();
			}
		}
	}

}
