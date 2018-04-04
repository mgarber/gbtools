package broad.pda.seq.pairedend;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import org.apache.commons.math.MathException;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.ReadMate;
import org.broad.igv.sam.SamAlignment;

import broad.core.datastructures.IntervalTree;
import broad.core.math.EmpiricalDistribution;
import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.annotation.GTFFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.gene.RefSeqGeneWithIsoforms;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;
import net.sf.samtools.util.CloseableIterator;

public class PairedEndDistribution {
	static Logger logger = Logger.getLogger(PairedEndDistribution.class.getName());
	private EmpiricalDistribution accurateDistribution;
	private EmpiricalDistribution smallBinDistribution;
	private AlignmentDataModel alignment;
	private double fractionOfPaired;
	private static int DEF_MIN_INS_SIZE = 0;
	private static int DEF_MAX_INS_SIZE = 1500;
	private static int DEF_NUM_BINS = 500;
	private static int DEF_SMALL_NUM_BINS = 30;
	private static double MAX_SD_FROM_INS_SIZE_MEAN = 2;
	
	public PairedEndDistribution(AlignmentDataModel data, BEDFileParser genes) throws IOException{
		//get all simple paths
		alignment = data;
		accurateDistribution = new EmpiricalDistribution(DEF_NUM_BINS, DEF_MIN_INS_SIZE, DEF_MAX_INS_SIZE);
		//smallBinDistribution = new EmpiricalDistribution(DEF_NUM_BINS, DEF_MIN_INS_SIZE, DEF_SMALL_NUM_BINS);
		smallBinDistribution = new EmpiricalDistribution(DEF_SMALL_NUM_BINS, DEF_MIN_INS_SIZE, DEF_MAX_INS_SIZE);
		estimatePairedDistribution(genes);
	}
	
	public PairedEndDistribution(AlignmentDataModel data, RefSeqGene gene, boolean isSecondRead) throws IOException{
		//get all simple paths
		alignment = data;
		accurateDistribution = new EmpiricalDistribution(DEF_NUM_BINS, DEF_MIN_INS_SIZE, DEF_MAX_INS_SIZE);
		smallBinDistribution = new EmpiricalDistribution(DEF_SMALL_NUM_BINS, DEF_MIN_INS_SIZE, DEF_MAX_INS_SIZE);
		estimateStrandedReadSpecificPairedDistribution(gene, isSecondRead);
	}
	
	private void estimatePairedDistribution(BEDFileParser geneReader) throws IOException {
		List<String> shuffledChrList = new ArrayList<String>(alignment.getChromosomeLengths().keySet());
		logger.debug("gene reader chromosome list "+shuffledChrList.toString());
		Collections.shuffle(shuffledChrList);
		Iterator<String> chrIt = shuffledChrList.iterator();
		geneReader.makeGenes(0.1);
		BEDFileParser constituentIsoReader = new BEDFileParser(geneReader.toConstituentIsoformMap(true));

		while(  chrIt.hasNext()) {
			String chr = chrIt.next();
			if(! constituentIsoReader.containChr(chr)) {
				logger.warn("WARN: no constituent isoforms found for chromosome " + chr );
				continue;
			}
			logger.debug("Processing " + chr);
			int numReads = 0;
			int numPairedReads = 0;
			Iterator<RefSeqGeneWithIsoforms> constituentIsoIt = constituentIsoReader.getChrTree(chr).valueIterator();

			while(constituentIsoIt.hasNext()) {
				RefSeqGene constituentIso = constituentIsoIt.next();
				CloseableIterator<Alignment>  readIt = alignment.getAlignmentsOverlappingRegion(constituentIso.getAlignment());
				while(readIt.hasNext()) {
					SamAlignment aln = (SamAlignment) readIt.next();
					numReads++;
					ReadMate mate = aln.getMate();
					if(aln.isFirstOfPair() && mate!=null && mate.getStart() >0) {
						numPairedReads++;
						
						int pairStart = mate.getStart();
						int thisStart = aln.getAlignmentStart();

						int [] pairedAlnExonStarts = {Math.min(pairStart, thisStart), Math.max(pairStart, thisStart)};
						int [] pairedAlnExonEnds = {pairedAlnExonStarts[0] + aln.getReadLength(), pairedAlnExonStarts[1] + aln.getReadLength()}; //TODO: handle case where pairs have different lengths.
						
						RefSeqGene pairedAln = new RefSeqGene(chr,  pairedAlnExonStarts[0], pairedAlnExonEnds[1], 
								aln.getReadName(), aln.getScore(), "*",//aln.isNegativeStrand() ? "-" : "+", 
										pairedAlnExonStarts, pairedAlnExonEnds);
						if(constituentIso.overlaps(pairedAln)) {
							int tStart = constituentIso.genomicToTranscriptPosition(pairedAln.getStart());
							int tEnd   = constituentIso.genomicToTranscriptPosition(pairedAln.getEnd());
							if(tStart > -1 && tEnd >-1) {
								double insertLength = (double) Math.abs(tEnd - tStart);
								accurateDistribution.add(insertLength);
								smallBinDistribution.add(insertLength);
							}
						}
					}
				}
				
				readIt.close();
				if(numReads % 500000 == 0) {
					logger.info("\tProcessed " + numReads + " reads ("+numPairedReads+" paired)");
				}
			}
			this.fractionOfPaired = numPairedReads/(double)numReads;
		}
	}

	
	private void estimateStrandedReadSpecificPairedDistribution(RefSeqGene gene, boolean isSecondRead) throws IOException {

		CloseableIterator<Alignment>  readIt = alignment.getAlignmentsOverlappingRegion(gene.getAlignment());
		while(readIt.hasNext()) {
			SamAlignment aln = (SamAlignment) readIt.next();
			//ONLY PROCESS THIS READ IF ITS ORIENTATION AND FIRST/SECOND IN PAIR AGREE WITH THE REQUIREMENTS
			if((gene.isNegativeStrand() == aln.isNegativeStrand())&&(isSecondRead == aln.isSecondOfPair())){
				ReadMate mate = aln.getMate();
				if(mate!=null && mate.getStart() >0) {
					
					int pairStart = mate.getStart();
					int thisStart = aln.getAlignmentStart();

					int [] pairedAlnExonStarts = {Math.min(pairStart, thisStart), Math.max(pairStart, thisStart)};
					int [] pairedAlnExonEnds = {pairedAlnExonStarts[0] + aln.getReadLength(), pairedAlnExonStarts[1] + aln.getReadLength()}; //TODO: handle case where pairs have different lengths.
						
					RefSeqGene pairedAln = new RefSeqGene(gene.getChr(),  pairedAlnExonStarts[0], pairedAlnExonEnds[1], 
								aln.getReadName(), aln.getScore(), "*",//aln.isNegativeStrand() ? "-" : "+", 
										pairedAlnExonStarts, pairedAlnExonEnds);
					
					int tStart = gene.genomicToTranscriptPosition(pairedAln.getStart());
					int tEnd   = gene.genomicToTranscriptPosition(pairedAln.getEnd());
					if(tStart > -1 && tEnd >-1) {
						double insertLength = (double) Math.abs(tEnd - tStart);
						accurateDistribution.add(insertLength);
						smallBinDistribution.add(insertLength);
					}
				}
			}
		}
		readIt.close();
	}
	
	public EmpiricalDistribution getAccurateSizeDistribution(){return this.accurateDistribution;}
	public EmpiricalDistribution getLowBinSizeDistribution(){return this.smallBinDistribution;}
	
	public double observedPValue(double obs) {
		return accurateDistribution.getPValue(obs, 0);
	}
	
	public PairedEndAnalysisResult doPairedEndCompatibilityAnalysis(RefSeqGene transcript, double pvalToRejectPairedInsertSize) throws IOException {
		CloseableIterator<Alignment> readIt = alignment.getReadIterator(transcript.getAlignment());
		PairedEndAnalysisResult pear = new PairedEndAnalysisResult(transcript,  this);
		while(readIt.hasNext()) {
			SamAlignment read = (SamAlignment)readIt.next();
			ReadMate mate = read.getMate();
			if(read.isFirstOfPair() && mate!=null && mate.getStart() >0) {
				int readStart = Math.min(read.getAlignmentStart(), mate.getStart());
				int readEnd   = Math.max(read.getAlignmentEnd(), mate.getStart() + read.getReadLength());//ASSUMED BOTH PAIRED READS HAVE SAME SIZE
				int transcriptInsertStart = transcript.genomicToTranscriptPosition(readStart);
				int transcriptInsertEnd = transcript.genomicToTranscriptPosition(readEnd);
				int s = Math.min(transcriptInsertStart, transcriptInsertEnd);
				int e = Math.max(transcriptInsertStart, transcriptInsertEnd);
				if(transcriptInsertStart >= 0 && transcriptInsertEnd >= 0 ) {
					int insertSize = Math.abs(transcriptInsertStart - transcriptInsertEnd);
					if(observedPValue(insertSize) < pvalToRejectPairedInsertSize) {
						pear.addInvalidPair(insertSize, s, e);
						logger.trace("discordant paired end insert("+insertSize+"): " + read.getReadName() +" "+read.getChr() +":"+read.getAlignmentStart() +"-" + read.getAlignmentEnd() + " mate: " + mate.getStart());
						//System.err.println("discordant paired end: " + read.getReadName() +" "+read.getChr() +":"+read.getAlignmentStart() +"-" + read.getAlignmentEnd() + " mate: " + mate.getStart());
					} else {
						pear.addValidPair(s, e);
					} 
				}else {
					logger.trace("read " + read.getReadName() +" "+read.getChr() +":"+read.getAlignmentStart() +"-" + read.getAlignmentEnd() +" with mate starting at " + mate.getStart() + " had negative overlapping start or end with transcript " + transcript.toBED());
				}
			}
		}
		readIt.close();
		return pear;
	}
	
	public BEDFileParser mergeReconstructionWithPairs(	Collection<RefSeqGene> genes) throws IOException {
		BEDFileParser geneReader = new BEDFileParser(genes);
		return mergeReconstructionWithPairs(geneReader);
	}

	public BEDFileParser mergeReconstructionWithPairs(BEDFileParser reconstructions) throws IOException {
		//Step2: Try join disjoint paths using paired ends
		Iterator<String>chrIt = reconstructions.getChromosomeIterator();
		BEDFileParser updatedRslt = new  BEDFileParser(); 
		int meanInsertPlusSD = (int) Math.round(accurateDistribution.getMean() + MAX_SD_FROM_INS_SIZE_MEAN * accurateDistribution.getStandardDeviation());
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			IntervalTree<RefSeqGeneWithIsoforms> chrGeneTree = reconstructions.getChrTree(chr);
			Iterator<RefSeqGeneWithIsoforms> geneIt = chrGeneTree.valueIterator();
			SortedMap<RefSeqGene, List<RefSeqGene>> isoformExtensionMap = new TreeMap<RefSeqGene, List<RefSeqGene>>();
			while(geneIt.hasNext()) {
				RefSeqGeneWithIsoforms g = geneIt.next();
				Collection<RefSeqGene> isoforms = g.getAllIsoforms();
				for(RefSeqGene iso : isoforms) {		
					Alignments lastExon = iso.getNumExons() > 0 ? iso.getLastExon() : iso.getAlignment();
					CloseableIterator<Alignment> lastExonAlignmentIt = alignment.getAlignmentsOverlappingRegion(lastExon);
					int coverage = 0;
					int overspanningReads = 0;
					//SortedMap<Alignments, List<ReadMate>> partnerPairMap = new TreeMap<Alignments, List<ReadMate>>();
					SortedMap<RefSeqGene, List<Integer>> partnerPairMap = new TreeMap<RefSeqGene, List<Integer>>();
					while(lastExonAlignmentIt.hasNext()) {
						Alignment aln = lastExonAlignmentIt.next();
						ReadMate mate = aln.getMate();
						boolean singleRead = mate == null || mate.getStart()==0;
						Alignments mateAnnotation = singleRead? null : new Alignments(iso.getChr(), mate.getStart(), mate.getStart()+1);

						if(singleRead ||  iso.overlapsExon(mateAnnotation) && aln.isFirstOfPair() ) {
							coverage++;
						}
						if(mate!=null && mate.getStart() >0 && mate.getStart() > lastExon.getEnd()) {
							overspanningReads ++;
							Iterator<RefSeqGene> partnerIt = partnerPairMap.keySet().iterator();
							int estimatedInsertSize = mate.getStart() + aln.getReadSequence().length() - aln.getStart() ;
							boolean alreadySeen = false;
							while(partnerIt.hasNext()) {
								RefSeqGene partner = partnerIt.next();
								if(partner.overlapsExon(mateAnnotation)){
									List<Integer> partnerMates = partnerPairMap.get(partner);
									partnerMates.add(estimatedInsertSize);
									alreadySeen = true;
								} else if (partner.getEnd() < mate.getStart()){
									break; // That is why we are using a sorted map. 
								}
							}
							if(!alreadySeen) {
								Iterator<RefSeqGeneWithIsoforms> overlapIt = reconstructions.getOverlappers(mateAnnotation).valueIterator();
								while(overlapIt.hasNext()) {
									RefSeqGeneWithIsoforms overlapGene = overlapIt.next();
									Iterator<RefSeqGene> overlapIsoIt = overlapGene.getAllIsoforms().iterator();
									while(overlapIsoIt.hasNext()) {
										RefSeqGene overlapperIso = overlapIsoIt.next();
										if(!overlapperIso.overlaps(iso)) { // This looks confusing but is simple: We only want NON-overlapping transcripts that overlap read mates
											Alignments overlappingExon = overlapperIso.getSingleOverlappingExon(mateAnnotation);
											if(overlappingExon != null) {
												List<Integer> insertSizes = new ArrayList<Integer>();
												insertSizes.add(estimatedInsertSize);
												partnerPairMap.put(overlapperIso, insertSizes);
											}
										}
									}
								}
							}
						}
					}
					lastExonAlignmentIt.close();
					List<RefSeqGene> partnersToStitch = new ArrayList<RefSeqGene>();
					if(!partnerPairMap.isEmpty()) {
						double minExpectedAnchoringPairs =  Math.floor((coverage/(double)lastExon.length())*(accurateDistribution.getMean() - accurateDistribution.getStandardDeviation())/(double)(2) );//This is an understimation!
						if(overspanningReads > minExpectedAnchoringPairs) {
							for(RefSeqGene org : partnerPairMap.keySet()) {
								if(partnerPairMap.get(org).size() > minExpectedAnchoringPairs/(double)2/*At least half of the overspanning reads of iso should land in the partner*/ && 
										Statistics.mean(partnerPairMap.get(org))<(meanInsertPlusSD)){
									partnersToStitch.add(org);
									logger.info("\tpartner exon: " + org.toBED() + " overspanning coverage " + partnerPairMap.get(org).size() + " mean insert size "+Statistics.mean(partnerPairMap.get(org)));
								}
							}
						} else {
							logger.debug("isoform " + iso.toBED() + " had coverage " + coverage + " and  "+ overspanningReads + " reads, which is less than the min expected " + minExpectedAnchoringPairs);
						}
					}
					logger.debug("NOW: Putting partners for, but does it exists already? iso: "+ isoformExtensionMap.containsKey(iso) + " iso is: "+ iso.toBED() +" partners tostitch: " + partnersToStitch + " if I get it back I get: " + isoformExtensionMap.get(iso));
					isoformExtensionMap.put(iso, partnersToStitch);
				}
			}
			updatedRslt.addRefSeqSet(stitch(isoformExtensionMap));
		}
		return updatedRslt;
	}



	private Collection<RefSeqGene> stitch(SortedMap<RefSeqGene, List<RefSeqGene>> isoformExtensionMap) {
		Collection<RefSeqGene> stitchedGenes = new ArrayList<RefSeqGene>();
		logger.debug("First call to stitch. Size of map: "+isoformExtensionMap.keySet().size());
		for(RefSeqGene gene : isoformExtensionMap.keySet()) {
			List<RefSeqGene> partners = isoformExtensionMap.get(gene);
			if(partners == null) {
				logger.error("NULL partners for " + gene.toBED() );
			}
			stitchedGenes.addAll(stitch(gene, isoformExtensionMap, partners));
		}
		return stitchedGenes;
	}

	private Collection<? extends RefSeqGene> stitch(RefSeqGene gene, SortedMap<RefSeqGene, List<RefSeqGene>> isoformExtensionMap, List<RefSeqGene> partners) {
		List<RefSeqGene> stitchedGenes = new ArrayList<RefSeqGene>();
		logger.debug("Stitching " + gene.toBED());
		if( gene.containsAttribute("toKill") ) {
			logger.debug("To kill skipping");
			return stitchedGenes;
		}
		if(partners == null) {
			logger.error("partners list was null for " + gene.toBED() );
		}
		if(partners.isEmpty()) {
			stitchedGenes.add(gene);
		} else {
			boolean wasStitched  = false;
			for(RefSeqGene partner : partners) {
				logger.debug("Partner is:  " + (partner != null ? partner.toBED() : "null!") );
				RefSeqGene stitchedGene = stitch(gene, partner, isoformExtensionMap);
				logger.debug("Stitched gene: " + (stitchedGene != null ? stitchedGene.toBED(): "is null"));
				if(stitchedGene != null) {
					if(!isoformExtensionMap.containsKey(partner)) {
						logger.error("Partner with null isoform extension map: " + partner.toBED() +" going to kill? " + partner.getAttribute("toKill"));
					} else {
						logger.debug("Key exists in isoformExtensionMap. partner is: " + isoformExtensionMap.get(partner));
						stitchedGenes.addAll(stitch(stitchedGene, isoformExtensionMap, isoformExtensionMap.get(partner)));
						wasStitched = true;
					}
				} else {
					logger.debug("Stitched gene was null");
				}
			} 
			if(!wasStitched) { //If no merge worked, just add the gene unchanged
				logger.debug("Was not stitched. adding as is");
				stitchedGenes.add(gene);
			}
		}

		return stitchedGenes;
	}


	/**
	 * Assumes partner is downstream of gene!
	 * @param gene
	 * @param partner
	 * @return
	 */
	private RefSeqGene stitch(RefSeqGene gene, RefSeqGene partner,  SortedMap<RefSeqGene, List<RefSeqGene>> isoformExtensionMap) {
		boolean canStitch = (gene.getOrientation().equals(partner.getOrientation()) || "*".equals(gene.getOrientation()) || "*".equals(partner.getOrientation()) ) & !gene.overlaps(partner);
		//For now lets not stitch together when both genes are multiexonice. 
		canStitch = canStitch && (gene.getNumExons() == 1 || partner.getNumExons() == 1); //TODO: Remove this condition and handle this case
		//For now lots only stitch genes that are within 2 SD of the insert size mean.
		int dist = gene.getAlignment().getDistanceTo(partner.getAlignment());
		canStitch = canStitch && dist>0 &&   dist < accurateDistribution.getMean() + MAX_SD_FROM_INS_SIZE_MEAN * accurateDistribution.getMean(); //TODO: Join such cases via spliced reads if such exists even at low frequencies
		
		RefSeqGene toReturn = null;
		//It is key to modify only the gene and not the partner as the intact partner is required later
		if(canStitch) {
			if(gene.getNumExons() == 1) { //Gene is single exon, we need to update partner's UTR
				RefSeqGene newGene = new RefSeqGene(partner);
				newGene.setUnorientedStart(gene.getStart());
				partner.addAttribute("toKill", "T");
				toReturn = newGene;
			} else if(partner.getNumExons() ==1) {
				toReturn =  new RefSeqGene(gene);
				toReturn.setUnorientedEnd(partner.getEnd());
				partner.addAttribute("toKill", "T");
				//toReturn = gene;
			}
		} 
		return toReturn;
	}

	
	private Map<String, RefSeqGene> findAllSimplePaths(BEDFileParser geneReader) {
		
		Map<String, RefSeqGene> rtrn=new HashMap<String, RefSeqGene>();
		Iterator<String> chrIt = geneReader.getChromosomeIterator();
		
		
		while(  chrIt.hasNext()) {
			String chr = chrIt.next();
			Iterator<RefSeqGeneWithIsoforms> chrGeneIt = geneReader.getChrTree(chr).valueIterator();
			while(chrGeneIt.hasNext()) {
				RefSeqGeneWithIsoforms gene = chrGeneIt.next();
				int overlapNum = 0;
				Iterator<RefSeqGeneWithIsoforms> geneOverlapIt = geneReader.getOverlappers(gene).valueIterator();
				while(geneOverlapIt.hasNext()) {
					overlapNum++;
				}
				
				if(overlapNum == 1) {
					rtrn.put(chr, gene);
				}
						
			}
		}

		return rtrn;
	}
	
	static String usage="Computes paired end derived insert-sized distribution \n -alignment <Paired alignment file in BAM, SAM or Alignemnt format> \n -annotations <transcript annotations set> -minMappingQuality <Minimum mapping quality to use>"+
			"\n";
	
	public static void main(String [] args) throws IOException {
		Globals.setHeadless(true);
		ArgumentMap argmap = CLUtil.getParameters(args, usage, "computeDist");
		if("computeDist".equalsIgnoreCase(argmap.getTask())) {
			String alnFile = argmap.getMandatory("alignment");
			String annotationFile = argmap.getMandatory("annotations");
			int minMapQual = argmap.containsKey("minMappingQuality") ? argmap.getInteger("minMappingQuality") : 0;
			
			AlignmentDataModel dataModel  = new  GenericAlignmentDataModel(alnFile, null, false, minMapQual);
			BEDFileParser annotReader = annotationFile.endsWith(".gtf") || annotationFile.endsWith(".GTF")? new GTFFileParser(annotationFile) : new BEDFileParser(annotationFile);
			PairedEndDistribution ped = new PairedEndDistribution(dataModel, annotReader);
			EmpiricalDistribution ed = ped.getAccurateSizeDistribution();
			BufferedWriter bw = argmap.getOutputWriter();
			ed.write(bw);
			bw.close();
			System.err.println("Mean insert size: " + ed.getMean() + ", sd: " + ed.getStandardDeviation());
		} else {
			System.out.println("Not supported task. Usage is " + usage);
		}

	}

	public static class PairedEndAnalysisResult {
		private RefSeqGene transcript;
		private double meanInsertSize;
		private int [] goodReadCov; 
		private int [] readsWithInsertTooLarge;
		private int [] readsWithInsertTooShort; 
		private EmpiricalDistribution [] baseCoveredInsertSizeDist;
		private boolean deviantIsoform = false;
		private EmpiricalDistribution transcriptInsertSizeDistribution; 
		private PairedEndDistribution parent;
		
		public PairedEndAnalysisResult(RefSeqGene transcript,  PairedEndDistribution parent) {
			this.transcript = transcript;
			this.parent = parent;
			this.meanInsertSize = parent.getAccurateSizeDistribution().getMean();
			
			int transcriptLength = transcript.getSize();
			goodReadCov = new int[transcriptLength];
			readsWithInsertTooLarge = new int[transcriptLength];
			readsWithInsertTooShort = new int[transcriptLength];
			baseCoveredInsertSizeDist = new EmpiricalDistribution[transcriptLength];
			transcriptInsertSizeDistribution = new EmpiricalDistribution(parent.smallBinDistribution.getBinNumber(), parent.smallBinDistribution.getMin(), parent.smallBinDistribution.getMax());
		}

		private void addValidPair(int s, int e) {
			for(int i = s ; i < e; i++) {
				goodReadCov[i]++;
				addToPositionDistribution(i, e-s);
			}
			transcriptInsertSizeDistribution.add(e-s);
		}

		private void addInvalidPair(int insertSize, int s, int e) {
			for(int i = s ; i < e; i++) {
				if(insertSize < meanInsertSize) {
					readsWithInsertTooShort[i]++;
				} else {
					readsWithInsertTooLarge[i]++;
				}
				addToPositionDistribution(i, insertSize);
			}
			transcriptInsertSizeDistribution.add(Math.min(insertSize, parent.accurateDistribution.getMax()-1));
		}
		
		private void addToPositionDistribution(int position, int size) {
			EmpiricalDistribution posDist = baseCoveredInsertSizeDist[position];
			if(posDist == null) {
				posDist = new EmpiricalDistribution(parent.smallBinDistribution.getBinNumber(), parent.smallBinDistribution.getMin(), parent.smallBinDistribution.getMax());
				baseCoveredInsertSizeDist[position] = posDist;
			}
			posDist.add(Math.min(size, parent.smallBinDistribution.getMax()-1));
			
		}

		public boolean isDeviant() {return deviantIsoform;}
		public EmpiricalDistribution getTranscriptInsertSizeDistribution(){return transcriptInsertSizeDistribution;}
		
		public RefSeqGene endFixedTranscript(double pvalToRejectPairedInsertSize) throws IllegalArgumentException, MathException {
		//1) Scan and see if there are regions that make paired ends non-compatible
			RefSeqGene fixedTranscript = new RefSeqGene(transcript);
			boolean isMultiExon = transcript.getNumExons() >1; 
			Alignments fivePExon = isMultiExon ?  transcript.get5PrimeExon() : transcript.getAlignment();
			if(isMultiExon && "*".equals(transcript.getOrientation())) {
				fivePExon = transcript.getFirstExon();
				logger.warn("transcript " + transcript.toBED() + " is multiexonic but un-oriented");
			}
			//System.err.println("transcript # exons " + transcript.getNumExons() + " what is the 5' exon? " + transcript.get5PrimeExon());
			if(fivePExon ==  null ) {
				logger.error("transcript " + transcript.toBED() + " has null 5' exon");
			}
			int endOf5PUtrExon = "-".equals(transcript.getOrientation()) ? transcript.genomicToTranscriptPosition(fivePExon.getStart()) :  transcript.genomicToTranscriptPosition(fivePExon.getEnd());
			Alignments threePExon = isMultiExon ? transcript.get3PrimeExon() : transcript.getAlignment();
			if(isMultiExon && "*".equals(transcript.getOrientation())) {
				threePExon = transcript.getLastExon();
			}
			int startOf3PUtrExon =  "-".equals(transcript.getOrientation()) ? transcript.genomicToTranscriptPosition(threePExon.getEnd()) : transcript.genomicToTranscriptPosition(threePExon.getStart());
			int transcriptLength = transcript.getSize();
			int newStart = 0;
			int newEnd = transcriptLength;
			for(int i = 0; i < transcriptLength; i++) {
				double chiTest =1.0;
				double insertMeanInsSizeAtBase = 0;
				double insertInsSizeAtBaseSD = 0;
				if (baseCoveredInsertSizeDist[i]!= null && baseCoveredInsertSizeDist[i].getMean()>0) {
					chiTest =  parent.smallBinDistribution.testGoodnessOfFit(baseCoveredInsertSizeDist[i]);
					insertMeanInsSizeAtBase = baseCoveredInsertSizeDist[i].getMean();
					insertInsSizeAtBaseSD = baseCoveredInsertSizeDist[i].getStandardDeviation();
				}

				//double expectedDeviantReads = Math.ceil((readsWithInsertTooShort[i] + readsWithInsertTooLarge[i] + goodReadCov[i])* pvalToRejectPairedInsertSize) + 1;
				//if(readsWithInsertTooLarge[i]> expectedDeviantReads || readsWithInsertTooShort[i] > expectedDeviantReads) { //If number of observed bad pairs is higher than expected

				if((goodReadCov[i] == 0 || chiTest < 0.01 && Math.abs(baseCoveredInsertSizeDist[i].getMean() - parent.accurateDistribution.getMean())> 2*parent.accurateDistribution.getStandardDeviation() ) ) { //NO paired ends this typically means that a segmented region was segmented to include an overlapping portion of a different transcript
					if(i < endOf5PUtrExon ){
						newStart = i;
					}else if ( i > startOf3PUtrExon ) {
						newEnd = i;
					} else if(goodReadCov[i] > 0 ) {
						deviantIsoform = true;
						logger.info("Position "+ transcript.transcriptToGenomicPosition(i) + " is deviant ("+goodReadCov[i]+"," + readsWithInsertTooLarge[i] + ","+readsWithInsertTooShort[i]+"), (mean,sd)=("+insertMeanInsSizeAtBase+","+insertInsSizeAtBaseSD+"), chiSquared test: " +chiTest + " background (mean,sd)=("+parent.getAccurateSizeDistribution().getMean() +","+parent.getAccurateSizeDistribution().getStandardDeviation()+")");
						//break;
					}
				}
			}
			if(newStart > 0 ) {
				int newGenomicStart = transcript.transcriptToGenomicPosition(newStart);
				fixedTranscript.set5PrimeEnd(newGenomicStart);
			}
			if(newEnd < transcriptLength) {
				int newGenomicEnd = transcript.transcriptToGenomicPosition(newEnd);
				fixedTranscript.set3PrimeEnd(newGenomicEnd);
			}
			return fixedTranscript;
		}
		
		public RefSeqGene endFixedTranscriptNew(double pvalToRejectPairedInsertSize) throws IllegalArgumentException, MathException {
		//1) Scan and see if there are regions that make paired ends non-compatible
			int window = (int) parent.getAccurateSizeDistribution().getMean()*2;
			int step   = window/2;
			int transcriptLength = transcript.getTranscriptLength();
			int length = window;
			RefSeqGene relativeWindow = transcript.trim(0, window);
			double adjustedPvalue = pvalToRejectPairedInsertSize / (double)(transcriptLength/(window - step));
			//while(length < transcriptLength) 
			RefSeqGene fixedTranscript = new RefSeqGene(transcript);
			boolean isMultiExon = transcript.getNumExons() >1; 
			Alignments fivePExon = isMultiExon ?  transcript.get5PrimeExon() : transcript.getAlignment();
			if(isMultiExon && "*".equals(transcript.getOrientation())) {
				fivePExon = transcript.getFirstExon();
				logger.warn("transcript " + transcript.toBED() + " is multiexonic but un-oriented");
			}
			//System.err.println("transcript # exons " + transcript.getNumExons() + " what is the 5' exon? " + transcript.get5PrimeExon());
			if(fivePExon ==  null ) {
				logger.error("transcript " + transcript.toBED() + " has null 5' exon");
			}
			int endOf5PUtrExon = "-".equals(transcript.getOrientation()) ? transcript.genomicToTranscriptPosition(fivePExon.getStart()) :  transcript.genomicToTranscriptPosition(fivePExon.getEnd());
			Alignments threePExon = isMultiExon ? transcript.get3PrimeExon() : transcript.getAlignment();
			if(isMultiExon && "*".equals(transcript.getOrientation())) {
				threePExon = transcript.getLastExon();
			}
			int startOf3PUtrExon =  "-".equals(transcript.getOrientation()) ? transcript.genomicToTranscriptPosition(threePExon.getEnd()) : transcript.genomicToTranscriptPosition(threePExon.getStart());

			int newStart = 0;
			int newEnd = transcriptLength;
			for(int i = 0; i < transcriptLength; i++) {
				double chiTest =1.0;
				double insertMeanInsSizeAtBase = 0;
				double insertInsSizeAtBaseSD = 0;
				if (baseCoveredInsertSizeDist[i]!= null && baseCoveredInsertSizeDist[i].getMean()>0) {
					chiTest =  parent.smallBinDistribution.testGoodnessOfFit(baseCoveredInsertSizeDist[i]);
					insertMeanInsSizeAtBase = baseCoveredInsertSizeDist[i].getMean();
					insertInsSizeAtBaseSD = baseCoveredInsertSizeDist[i].getStandardDeviation();
				}

				//double expectedDeviantReads = Math.ceil((readsWithInsertTooShort[i] + readsWithInsertTooLarge[i] + goodReadCov[i])* pvalToRejectPairedInsertSize) + 1;
				//if(readsWithInsertTooLarge[i]> expectedDeviantReads || readsWithInsertTooShort[i] > expectedDeviantReads) { //If number of observed bad pairs is higher than expected

				if((goodReadCov[i] == 0 || chiTest < 0.01 && Math.abs(baseCoveredInsertSizeDist[i].getMean() - parent.accurateDistribution.getMean())> 2*parent.accurateDistribution.getStandardDeviation() ) ) { //NO paired ends this typically means that a segmented region was segmented to include an overlapping portion of a different transcript
					if(i < endOf5PUtrExon ){
						newStart = i;
					}else if ( i > startOf3PUtrExon ) {
						newEnd = i;
					} else if(goodReadCov[i] > 0 ) {
						deviantIsoform = true;
						logger.info("Position "+ transcript.transcriptToGenomicPosition(i) + " is deviant ("+goodReadCov[i]+"," + readsWithInsertTooLarge[i] + ","+readsWithInsertTooShort[i]+"), (mean,sd)=("+insertMeanInsSizeAtBase+","+insertInsSizeAtBaseSD+"), chiSquared test: " +chiTest + " background (mean,sd)=("+parent.getAccurateSizeDistribution().getMean() +","+parent.getAccurateSizeDistribution().getStandardDeviation()+")");
						//break;
					}
				}
			}
			if(newStart > 0 ) {
				int newGenomicStart = transcript.transcriptToGenomicPosition(newStart);
				fixedTranscript.set5PrimeEnd(newGenomicStart);
			}
			if(newEnd < transcriptLength) {
				int newGenomicEnd = transcript.transcriptToGenomicPosition(newEnd);
				fixedTranscript.set3PrimeEnd(newGenomicEnd);
			}
			return fixedTranscript;
		}
		
		
	

		
		
	}
	
	public void ensureNonZeroCounts() {
		accurateDistribution.addPsuedocounts(1);
		smallBinDistribution.addPsuedocounts(1);
		
	}
}
