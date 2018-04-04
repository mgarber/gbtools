package broad.pda.chromosome;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Stack;

import broad.core.alignment.Repeat;
import broad.core.alignment.RepeatMaskerReader;
import broad.core.alignment.RepeatMaskerReader.RepeatStatistic;
import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.GFF;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.core.sequence.WindowSlider;
import broad.pda.assembly.AgpEntry;
import broad.pda.assembly.AgpEntryFactory;
import broad.pda.datastructures.Alignments;
import broad.pda.snp.DBSNPReader;
import broad.pda.snp.DBSNPReader.DBSNP;

public class Chromosome {
	GenomicAnnotation centromere;
	List<AgpEntry> gaps;
	List<AgpEntry> clones;
	RepeatMaskerReader repeatReader;
	DBSNPReader dbSNPReader;
	String chromosomeNumber;
	Sequence sequence;
	int size;
	File sequenceFile;
	private AgpEntry shortArm;
	private List<GFF> bands;
	
	protected Chromosome() {
		super();
		gaps = new ArrayList<AgpEntry>();
		clones = new ArrayList<AgpEntry>();
		bands = new ArrayList<GFF>();
		repeatReader = new RepeatMaskerReader();
	}
	
	/**
	 * @param The name of the chromosome    e.g. "chr12"
	 */
	public Chromosome(final String name) {
		this.chromosomeNumber = name;
		gaps = new ArrayList<AgpEntry>();
		clones = new ArrayList<AgpEntry>();
	}
	

	public void loadAGP(String agpFile) throws Exception {

		size = 0;
		File source = new File(agpFile);
		BufferedReader br = new BufferedReader(new FileReader(source));
		String agpWithoutExt = agpFile.substring(0,agpFile.lastIndexOf("."));
		sequenceFile = new File(agpWithoutExt + ".fa");
		String line;

		try {
			while((line = br.readLine()) != null) {
				if(line.startsWith("#") || line.trim().length() ==0){
					continue;
				}
				String[] lineSplit = line.split("\t");
				//System.out.println("Doing chromosome " + lineSplit[0]);
				if(chromosomeNumber == null) {
					chromosomeNumber = lineSplit[0].substring(3);
					if(chromosomeNumber.startsWith("0")) {
						chromosomeNumber = chromosomeNumber.substring(1);
					}
				}
				AgpEntry entry = AgpEntryFactory.createEntry(lineSplit);
				switch (entry.getType()) {
				case AgpEntry.CLONE_TYPE :
					clones.add(entry);
					break;
				case AgpEntry.GAP_TYPE :
					gaps.add(entry);
					break;
				case AgpEntry.CENTROMERE_TYPE :
					centromere = entry;
					break;
				case AgpEntry.SHORT_ARM_TYPE : 
					shortArm = entry;
					break;
				}
				size = Math.max(size, entry.getEnd());
			}
		}  finally {
			try {
				//System.out.print("Closing "+agpFile);
				br.close();
				//System.out.print(" ..... Closed\n");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	public String toString() {
		StringBuffer buf = new StringBuffer("chr");
		buf.append(getSymbol());
		
		return buf.toString();
	}
	
	public int getSize() {
		return size;
	}
	
	public long getUngappedSize() {
		long totalGaps = shortArm != null ? shortArm.getLength() : 0;
		
		Iterator<AgpEntry> gapIt = gaps.iterator();
		while(gapIt.hasNext()) {
			totalGaps += gapIt.next().getLength();
		}
		
		return getSize() - totalGaps;
	}
	
	public void loadSNPs() throws FileNotFoundException {
		if(dbSNPReader == null) {
			System.out.print("loading SNPs for chromosome " + getSymbol()+".....");
			dbSNPReader = new DBSNPReader(getSymbol());
			System.out.println("Done");
		}
	}
	
	public DBSNP getSNP(String dbSNPId) {
		if(dbSNPReader == null) {
			throw new IllegalStateException("You cannot get SNP locations until you call loadSNPs()");
		}
		return dbSNPReader.getSNP(dbSNPId);
	}
	 
	public String getSequenceFileName() { return sequenceFile.getAbsolutePath();}
	
	public void loadSequence() throws IOException {
		FastaSequenceIO fsio = new FastaSequenceIO(sequenceFile);
		//System.err.println("  Extracting ... expted size " + size);
		fsio.extractRecordsWithIDLike(getSymbol(), false, size);
		sequence = fsio.loadAll().get(0);
		//size = sequence.getSequenceBases().length();
	}
	
	public void unloadSequence() { sequence.unloadAllSequences(); System.gc();} 
	
	public Sequence getSequence () { return sequence; }
	
	public WindowSlider getSlider(int windowSize, int overlap) {
		WindowSlider slider = null;
		if(sequence != null) {
			slider = sequence.getSlider(windowSize, overlap);
		} else {
			BasicGenomicAnnotation fullChr = new BasicGenomicAnnotation("chr" + getSymbol(), getSymbol(), 1, (int) getSize());
			SequenceRegion region = new SequenceRegion(getSymbol(), fullChr);
			slider = region.getSlider(windowSize, overlap);
		}
		return slider;
	}
	
	public void getRegion(SequenceRegion region) {
		sequence.getRegion(region);
	}
	
	public void getRegion(SequenceRegion region, boolean softmask) {
		sequence.getRegion(region, softmask);
	}
	
	public void getRegion(SequenceRegion region, boolean softmask, Map<String, IntervalTree<Alignments>> okRepeats) {
		sequence.getRegion(region, softmask, okRepeats);
	}
	
	public List<SequenceRegion> getRegions(List<? extends GenomicAnnotation> annotations, int padding) throws IOException {
		FastaSequenceIO fsIO = new FastaSequenceIO(sequenceFile.getAbsolutePath());
		ArrayList<SequenceRegion> regs = new ArrayList<SequenceRegion>(annotations.size());
		Iterator<? extends GenomicAnnotation> it = annotations.iterator();
		while(it.hasNext()) {
			GenomicAnnotation annot = it.next();
			SequenceRegion reg = new SequenceRegion("chr"+getSymbol());
			reg.setRegionStart(annot.getStart() - padding);
			reg.setRegionEnd(annot.getEnd() + padding);
			reg.setId(annot.getName());
			regs.add(reg);
		}
		fsIO.extractRegions(regs);
		return regs;
	}
	
	public void extractRegions(List<? extends SequenceRegion> regions) throws IOException {
		if(sequence == null ){//|| sequence.getSequenceBases() == null || sequence.getSequenceBases().length() == 0) {
			loadSequence();
		}
		sequence.getRegions(regions);
	}
	
	public void extractRegion(SequenceRegion region) throws IOException {
		ArrayList<SequenceRegion> regionList = new ArrayList<SequenceRegion>(1);
		regionList.add(region);
		extractRegions(regionList);
	}
	
	
	public String getSymbol() {return chromosomeNumber;}
	
	public void loadRepeatInfo(File repeatFile, boolean clean) throws FileNotFoundException {
		if(repeatReader.getRepeatList() == null || repeatReader.getRepeatList().size() == 0) {
			repeatReader.loadRepeats(repeatFile, clean);
		}
	}
	
	public void loadRepeatInfo(File repeatFile, boolean clean, LightweightGenomicAnnotation inRegion) throws FileNotFoundException {
		repeatReader.loadRepeats(repeatFile, clean, inRegion);
	}
	
	/**
	 * Shuffle list in clusters. This method clusters genomic annotations in list then shuffles the clusters as though 
	 * the list provided is not random but you can define clusters that behave randomly.
	 * @param toShuffle The list of annotations to shuffle
	 * @param toAvoid regions to avoid placing shuffled annotations
	 * @param clusterMaxDistance maximum distance from a genomic annotation to a cluster for it to be included in it.
	 * @return cluster preserving shuffle.
	 */
	public List<? extends GenomicAnnotation> shuffleInClusters(List<? extends GenomicAnnotation> toShuffle, 
			List<? extends GenomicAnnotation> toAvoid, 
			int clusterMaxDistance) {
		
		List<GenomicAnnotation> valid = getValidRegions(toAvoid);
		
		
		//Sort toShuffle list, then define clusters
		List<GenomicAnnotation> chrAnnotations = new ArrayList<GenomicAnnotation>(toShuffle.size());
		Iterator<? extends GenomicAnnotation> toShuffleIt = toShuffle.iterator();
		while(toShuffleIt.hasNext()) {
			GenomicAnnotation annot = toShuffleIt.next();
			if(getSymbol().equals(annot.getChromosome())) {
				chrAnnotations.add(annot);
			}
		}
		Collections.sort(chrAnnotations);
		Stack<Stack<GenomicAnnotation>> clusterList = new Stack<Stack<GenomicAnnotation>>();
		Stack<GenomicAnnotation> clusterReps       = new Stack<GenomicAnnotation>();
		
		toShuffleIt = chrAnnotations.iterator();
		if(toShuffleIt.hasNext()) {
			GenomicAnnotation annot = toShuffleIt.next();
			Stack<GenomicAnnotation> firstCluster = new Stack<GenomicAnnotation>();
			clusterList.push(firstCluster);
			firstCluster.push(annot);
			clusterReps.push(annot);		
			//System.out.println("First cluster rep: " + annot.getName());
		}
		while(toShuffleIt.hasNext()) {
			GenomicAnnotation annot = toShuffleIt.next();
			GenomicAnnotation last = clusterList.peek().peek();
			if(annot.getDistanceTo(last) < clusterMaxDistance) {
				clusterList.peek().push(annot);
				//System.out.println("\tAdding " + annot.getName() + " to cluster of rep " + clusterList.peek().firstElement().getName());
			} else {
				clusterReps.push(annot);
				Stack<GenomicAnnotation> newCluster = new Stack<GenomicAnnotation>();
				newCluster.push(annot);
				clusterList.push(newCluster);
				//System.out.println("New cluster: rep " + annot.getName());
			}
		}
		
		//Suffle reps, but do it safely, one at time, we can optimize this later.
		Iterator<GenomicAnnotation> repIt = clusterReps.iterator();
		Iterator<Stack<GenomicAnnotation>> clusterIt = clusterList.iterator();
		ArrayList<GenomicAnnotation> shuffledAnnotations = new ArrayList<GenomicAnnotation>(toShuffle.size());
		while(repIt.hasNext()) {
			GenomicAnnotation rep = repIt.next();
			List<? extends GenomicAnnotation> shuffled = shuffleWithinValid(rep, valid);
			//For simplicity, insist that the shuffle is all whithin one region.
			while(shuffled.size() != 1) { //There should be no risk of infinite look here for annotations that are not ridiculously long.
				shuffled = shuffle(rep, toAvoid);
			}
			GenomicAnnotation shuffledRep = shuffled.get(0);
			
			GenomicAnnotation containingRegion = null;
			int containingRegionNum = 0;
			
			Iterator<GenomicAnnotation> validIt = valid.iterator();
			while(validIt.hasNext()) {
				GenomicAnnotation a = validIt.next();
				//System.out.println("Iterating through valid regions " + a + " does it contain shuffled rep " + shuffledRep);
				if(a.contains(shuffledRep)) {
					//System.out.println ("\t\tYES");
					containingRegion = a;
					break;
				}
				containingRegionNum++;
			}
			shuffledAnnotations.add(shuffledRep); 
			int repShift = shuffledRep.getStart() - rep.getStart();
			Stack<GenomicAnnotation> cluster = clusterIt.next();
			while(cluster.size() > 1) {
				GenomicAnnotation shuffledAnnot = new BasicGenomicAnnotation(cluster.pop());
				shuffledAnnot.setName(shuffledAnnot.getName() + "_suffled");
				shuffledAnnot.setStart(shuffledAnnot.getStart() + repShift);
				shuffledAnnot.setEnd(shuffledAnnot.getEnd() + repShift);
				shuffledAnnotations.add(shuffledAnnot);
				if(!containingRegion.contains(shuffledAnnot)) {
					int overflow = shuffledAnnot.getEnd() - containingRegion.getEnd();
					shuffledAnnot.setEnd(containingRegion.getEnd());
					LightweightGenomicAnnotation nextValidRegion = valid.get((containingRegionNum + 1) % valid.size());
					GenomicAnnotation shuffledAnnotPart2 = new BasicGenomicAnnotation(shuffledAnnot.getName() + "_II");
					shuffledAnnotPart2.setStart(nextValidRegion.getStart());
					shuffledAnnotPart2.setEnd(nextValidRegion.getStart() + overflow - 1);
					shuffledAnnotations.add(shuffledAnnotPart2);
				}

			}
			
		}
		
		return shuffledAnnotations;
	}
	
	public List<? extends GenomicAnnotation> shuffle(GenomicAnnotation toShuffle, List<? extends GenomicAnnotation> toAvoid) {
		ArrayList<GenomicAnnotation> oneMemberList = new ArrayList<GenomicAnnotation>(1);
		oneMemberList.add(toShuffle);
		return shuffle(oneMemberList, toAvoid);
	}
	
	public List<? extends GenomicAnnotation> shuffleWithinValid(GenomicAnnotation toShuffle, List<? extends GenomicAnnotation> valid) {
		ArrayList<GenomicAnnotation> oneMemberList = new ArrayList<GenomicAnnotation>(1);
		oneMemberList.add(toShuffle);
		return shuffleWithinValid(oneMemberList, valid);
	}
	
	/**
	 * Randomizes (using a quasi uniform distribution) a list of genomic annotations
	 * @param toShuffle List of genomic annotations to shuffle
	 * @param toAvoid regions to avoid placing a shuffled annotation.
	 * @return
	 */
	public List<? extends GenomicAnnotation> shuffle(List<? extends GenomicAnnotation> toShuffle, List<? extends GenomicAnnotation> toAvoid) {
		List<GenomicAnnotation> valid = getValidRegions(toAvoid);
		//System.out.println("Regions to avoid: " + (toAvoid != null ? toAvoid.size() : 0) + " so we got valid regions " + (valid != null ? valid.size() : 0));
		return shuffleWithinValid(toShuffle,  valid);
	}

	public List<? extends GenomicAnnotation> shuffleWithinValid(List<? extends GenomicAnnotation> toShuffle, List<? extends GenomicAnnotation> valid) {
		ArrayList<GenomicAnnotation> randomized = new ArrayList<GenomicAnnotation>();
		Random randomizer = new Random(Math.round(Math.random()*1000000) );
		// The next hack attempts to simulate the trivial fact that large regions should be more likely
		// to be drawn by adding as many copies of it as its percent of the total length covered by
		// valid regions. The caveat is that we are rounding the percent thus we used 1000 rather than 100 
		// to account for small but not tiny regions... This may slow things down significantly when randomizing
		// al large list.
		//System.out.println("Total valid regions length " + totalValidLength);
		List<GenomicAnnotation> normalizedValid = normalizeRegionsListBySize(valid);
		
		Iterator<? extends GenomicAnnotation> toRandomizeIt = toShuffle.iterator();
		
		int normalizedValidNumber = normalizedValid.size();
		//System.out.println("Valid num: " + validNumber + " validNumber order of magnitude " + orderOfMagnitudeValid + " normlizedValidSize " + normalizingConstant + " total normalized valid " + normalizedValidNumber);
		while(toRandomizeIt.hasNext()) {
			GenomicAnnotation original = toRandomizeIt.next();
			int randomValidRegion = randomizer.nextInt(normalizedValidNumber);
			GenomicAnnotation validReg = normalizedValid.get(randomValidRegion);
			int randomizedStart = randomizer.nextInt(validReg.getLength()) + validReg.getStart();
			GenomicAnnotation randomizedOriginal = new BasicGenomicAnnotation(original);
			randomizedOriginal.setName(randomizedOriginal.getName() + "_randomized");
			randomized.add(randomizedOriginal);
			randomizedOriginal.setStart(randomizedStart);
			randomizedOriginal.setEnd(Math.min(randomizedStart + original.getLength(), validReg.getEnd()));
			GenomicAnnotation prior = randomizedOriginal;
			
			// Initially we got a random region from the NORMALIZED vector, from now on we work on the small valid region vector
			randomValidRegion = valid.indexOf(validReg);
			int i = 2;
			int left = original.getLength() - prior.getLength();
			while(left > 0) {
				//System.out.print("\tPrior valid region " + randomValidRegion);
				randomValidRegion = (randomValidRegion + 1) % valid.size();
				//System.out.print(" new valid region " + randomValidRegion);
				LightweightGenomicAnnotation nextValidRegion = valid.get(randomValidRegion);
				//System.out.println(" region: " + nextValidRegion);
				GenomicAnnotation nextRandomizedChunk = new BasicGenomicAnnotation(original);
				nextRandomizedChunk.setName(original.getName() + "_randomized_" + i++);
				nextRandomizedChunk.setStart(nextValidRegion.getStart());
				nextRandomizedChunk.setEnd(Math.min(nextValidRegion.getStart() + left, nextValidRegion.getEnd()));
				randomized.add(nextRandomizedChunk);
				left = left - nextRandomizedChunk.getLength();
			}
		}
		return randomized;
	}

	
	public SequenceRegion drawRandomRegion(int size) {
		GenomicAnnotation region = new BasicGenomicAnnotation("initial",getSymbol(), 1, size + 1);
		List<? extends GenomicAnnotation> randomized = shuffleWithinValid(region, getValidRegions(null));
		//TODO: Handle the case one more than one regions are returned!
		if(randomized.size() > 1) {
			System.out.println("WARNING: random region was actually randomized to 2 regions");
		}
		
		SequenceRegion toExtract = new SequenceRegion(getSymbol(), randomized.get(0));
		if(sequence != null && sequence.getSequenceBases().length() > 0) {
			getRegion(toExtract);
		}
		return toExtract;
	}
	
	public void insertAtRandomPoint(SequenceRegion region) {
		GenomicAnnotation insertPoint = new BasicGenomicAnnotation("insertionPoint",getSymbol(), 1, 2);
		List<? extends GenomicAnnotation> randomized = shuffleWithinValid(insertPoint, getValidRegions(null));
		insertPoint = randomized.get(0);
		
		insert(region, insertPoint);

		System.out.println("Inserted " + region + " in " + insertPoint);
	}
	
	public void insert(SequenceRegion region, LightweightGenomicAnnotation insertPoint) {
		StringBuilder newSequence = new StringBuilder(sequence.getSequenceBases().substring(0,insertPoint.getStart()));		
		newSequence.append(region.getSequence().getSequenceBases());
		newSequence.append(sequence.getSequenceBases().substring(insertPoint.getStart()));
		sequence.setSequenceBases(newSequence.toString());
		
		Iterator<AgpEntry> gapIt = gaps.iterator();
		while(gapIt.hasNext()) {
			AgpEntry gap = gapIt.next();
			if(gap.getStart() >= insertPoint.getStart()) {
				gap.setStart(gap.getStart() + region.getLength());
				gap.setEnd(gap.getEnd() + region.getLength());
			}
		}
		
		if(centromere.getStart() >= insertPoint.getStart()) {
			centromere.setStart(centromere.getStart() + region.getLength());
			centromere.setEnd(centromere.getEnd() + region.getLength());
		}
	}
	
	public void delete(SequenceRegion region) {
		StringBuilder newSequence = new StringBuilder(sequence.getSequenceBases().substring(0, region.getStart()));
		newSequence.append(sequence.getSequenceBases().substring(region.getEnd() + 1));
		sequence.setSequenceBases(newSequence.toString());
		
		Iterator<AgpEntry> gapIt = gaps.iterator();
		while(gapIt.hasNext()) {
			AgpEntry gap = gapIt.next();
			if(gap.getStart() >= region.getStart()) {
				gap.setStart(gap.getStart() - region.getLength());
				gap.setEnd(gap.getEnd() - region.getLength());
			}
		}
		
		if(centromere.getStart() >= region.getStart()) {
			centromere.setStart(centromere.getStart() - region.getLength());
			centromere.setEnd(centromere.getEnd() - region.getLength());
		}
	}
	
	public void invert(SequenceRegion region) {
		StringBuilder newSequence = new StringBuilder(sequence.getSequenceBases().substring(0, region.getStart()));
		region.reverse();
		newSequence.append(region.getSequenceBases());
		newSequence.append(sequence.getSequenceBases().substring(region.getEnd() + 1));
		sequence.setSequenceBases(newSequence.toString());
	}
	
	public void computeRepeatStatistics(File repeatFile, boolean clean) 
	throws FileNotFoundException {
		repeatReader = new RepeatMaskerReader(repeatFile, clean);
		System.out.println("Loaded repeats for chr" + getSymbol() + " obtained " + getRepeats().size() + " and did cleanup? " + clean);
		repeatReader.computeStatistics();
	}
	
	public List<Repeat> getRepeats() { return repeatReader.getRepeatList(); }
	
	public List<Repeat> getRepeats(int start, int end) {
		GenomicAnnotation target = new BasicGenomicAnnotation("target");
		
		target.setStart(start);
		target.setEnd(end);
		return repeatReader.getRepeatList(target);
	}
	
	public List<RepeatStatistic> getRepeatStats() {
		return repeatReader.getRepeatStatistics();
	}

	public List<RepeatStatistic> getRepeatFamilyStats() {
		return repeatReader.getRepeatFamilyStatistics();
	}

	public Iterator<? extends GenomicAnnotation> Gaps() {
		return gaps.iterator();
	}

	public GenomicAnnotation getCentromere() {
		if(centromere == null) {
			Iterator<AgpEntry> gapIt = gaps.iterator();
			AgpEntry curGap = null;
			while(gapIt.hasNext()) {
				curGap = gapIt.next();
				if(centromere == null || (curGap.getLength() > centromere.getLength())) {
					centromere = curGap;
				}
			}
		}
		return centromere;
	}
	
	public double getPercentOfGappedSequence(GenomicAnnotation region) {
		List<GenomicAnnotation> overlappingGaps = getGappedSubregions(region);
		
		Iterator<GenomicAnnotation> overlappingGapIt = overlappingGaps.iterator();
		int totalGappedLength = 0;
		while(overlappingGapIt.hasNext()) {
			GenomicAnnotation gap = overlappingGapIt.next();
			//System.out.println("Overlapping " + gap);
			totalGappedLength += gap.getLength();
		}

		return totalGappedLength/(double)region.getLength();
	}

	private List<GenomicAnnotation> getGappedSubregions(GenomicAnnotation region) {
		List<GenomicAnnotation> overlappingGaps = new ArrayList<GenomicAnnotation>();
		
		if(centromere != null && centromere.overlaps(region)) {
			GenomicAnnotation centromereCopy = new BasicGenomicAnnotation(centromere);
			centromereCopy.takeIntersection(region);
			overlappingGaps.add(centromereCopy);
		}if(shortArm != null && shortArm.overlaps(region)) {
			GenomicAnnotation shortArmCopy = new BasicGenomicAnnotation(centromere);
			shortArmCopy.takeIntersection(region);
			overlappingGaps.add(shortArmCopy);			
		}
		
		Iterator<AgpEntry> gapIt = gaps.iterator();
		while(gapIt.hasNext()) {
			AgpEntry gap = gapIt.next();
			if(gap.getStart() > region.getEnd()) {
				break;
			}
			if(gap.overlaps(region)) {
				GenomicAnnotation gapCopy = new BasicGenomicAnnotation(gap);
				gap.takeIntersection(region);
				overlappingGaps.add(gapCopy);
			}
		}
		
		//System.out.println(overlappingGaps);
		return overlappingGaps;
	}

	public void unloadRepeats() {
		repeatReader.clearRepeats();
	}
	
	public AgpEntry getShortArm() { return shortArm; }
	
	public List<AgpEntry> gaps() { return gaps;}

	public DBSNPReader getDbSNPReader() {
		return dbSNPReader;
	}

	public void setDbSNPReader(DBSNPReader dbSNPReader) {
		this.dbSNPReader = dbSNPReader;
	}

	public RepeatMaskerReader getRepeatReader() {
		return repeatReader;
	}

	public void setBanding(List<GFF> bands) {
		this.bands = bands;
	}
	
	public List<GFF> getBanding() { return bands;}

	private List<GenomicAnnotation> getValidRegions(List<? extends GenomicAnnotation> toAvoid) {
		if(toAvoid == null) {
			toAvoid = new ArrayList<GenomicAnnotation>(); //initialize if null to avoid null pointer exceptions....
		}
		ArrayList<GenomicAnnotation> toAvoidInChr = new ArrayList<GenomicAnnotation>();
			
		if(centromere != null)
			toAvoidInChr.add(centromere);
		if(shortArm != null) 
			toAvoidInChr.add(shortArm);
		
		toAvoidInChr.addAll(gaps());
		
		Iterator<? extends GenomicAnnotation> avoidIt = toAvoid.iterator();
		while(avoidIt.hasNext()) {
			GenomicAnnotation avoid = avoidIt.next();
			if(getSymbol().equals(avoid.getChromosome())) {
				toAvoidInChr.add(avoid);
				//System.out.println("\tAdded " + avoid + " to avoid in chr " + getSymbol());
			}
		}
		
		GenomicAnnotation allOfMe = new BasicGenomicAnnotation("chr" + getSymbol(),getSymbol(), 1, (int)getSize());
		
		List<GenomicAnnotation> valid = allOfMe.minus(toAvoidInChr);
		return valid;
	}
	
	public List<SequenceRegion> chunk(int chunkSize, int chunkOverlap) {
		return getSequence() == null 
			? chunkWithNoSequence(chunkSize, chunkOverlap)
			: getSequence().chunk(chunkSize, chunkOverlap);
	}

	public List<SequenceRegion> chunkWithNoSequence(int chunkSize, int chunkOverlap) {
		int start = 0;
		int end = chunkSize;
		List<SequenceRegion> chunks = new ArrayList<SequenceRegion>(getSize()/(chunkSize - chunkOverlap+1));
		while(start < getSize()) {
			SequenceRegion chunk = new SequenceRegion(getSymbol());
			chunk.setChromosome(getSymbol());
			chunk.setStart(start);
			chunk.setEnd(end);
			chunks.add(chunk);
			
			start = start + chunkSize - chunkOverlap;
			end = Math.min(getSize() , start + chunkSize);
		}
		return chunks;
	}

	private List<GenomicAnnotation> normalizeRegionsListBySize(List<? extends GenomicAnnotation> valid) {
		Iterator<? extends GenomicAnnotation> validIt = valid.iterator();
		int totalValidLength = 0;
		while(validIt.hasNext()) {
			totalValidLength += validIt.next().length();
		}
		
		int valids = valid.size();	
		validIt = valid.iterator();
		double orderOfMagnitudeValid = Math.ceil(Math.log10(valids));
		int normalizingConstant =(int) Math.pow(10, (orderOfMagnitudeValid + 2));

		List<GenomicAnnotation> normalizedValid = new ArrayList<GenomicAnnotation>(normalizingConstant);
		while(validIt.hasNext()) { 
			GenomicAnnotation reg = validIt.next();
			double roundedPctCoveratge = Math.round(reg.length()/(float)totalValidLength * normalizingConstant);
			//System.out.println("\tRegion " + reg + " is " + roundedPctCoveratge + " of total ");
			for(int i = 0; i < roundedPctCoveratge; i++) {
				//System.out.println("\t\tAdding region " + reg.getName() + ", total newValid size: " + newNormalizedValid.size());
				normalizedValid.add(reg);
			}
		}
		return normalizedValid;
	}



}
