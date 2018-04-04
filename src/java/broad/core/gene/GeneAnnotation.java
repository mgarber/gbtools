package broad.core.gene;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Stack;

import broad.core.annotation.BED;
import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.error.ParseException;

public class GeneAnnotation extends BasicGenomicAnnotation {

	private int cdsStart;
	private int cdsEnd;
	private int exonCount;
	private List<Exon> exons;
	
	// HG18 new fields
	public enum CDSStat {none, unk, incmpl, cmpl}
	private String  altName;
	private boolean completeCDSStart;
	private boolean completeCDSEnd;
	
	
	
	public GeneAnnotation(String[] data, boolean cdsOnly) throws ParseException {
		setFromRawData(data, cdsOnly, "HG17");
	}
	
	public GeneAnnotation(String [] strings, boolean cdsOnly, String version) throws ParseException {
		setFromRawData(strings, cdsOnly, version);
	}
	
	public GeneAnnotation(BED fullBed, boolean cdsOnly) {
		super(fullBed);
		List<GenomicAnnotation> blocks = fullBed.getBlocks();
		cdsStart = fullBed.getThickStart();
		cdsEnd   = fullBed.getThickEnd();
		exonCount = blocks.size();
		exons = new ArrayList<Exon>(blocks.size());
		Iterator<GenomicAnnotation> blockIt = blocks.iterator();
		while(blockIt.hasNext()) {
			GenomicAnnotation block = blockIt.next();
			if(!cdsOnly){
				Exon exon = new Exon(this, block);
				this.addExon(exon);
			} else if(block.overlaps(getCDS()) ) {
				block.takeIntersection(getCDS());
				Exon exon = new Exon(this, block);
				this.addExon(exon);
			}
		}
		
		//shift(-1);
		
	}
	

	public GeneAnnotation(BED bed) {
		this(bed, false);
	}
	
	public void setFromRawData(String [] data, boolean cdsOnly, String version) throws ParseException {
		int i = "HG18".equalsIgnoreCase(version) ? 1 : 0;
		
		setName(data[i++]);
		setChromosomeFromString(data[i++]);
		setOrientation(data[i++]);
		setStartFromString(data[i++]);
		setEndFromString(data[i++]);
		setCdsStartFromString(data[i++]);
		setCdsEndFromString(data[i++]);
		setExonCountFromString(data[i++]);
		String exonStarts = data[i++];
		String exonEnds   = data[i++];
		String exonFrames = null;
		if("HG18".equals(version)) {
			setId(data[i++]);
			altName = data[i++];
			completeCDSStart = "cmpl".equals(data[i++]);
			completeCDSEnd   = "cmpl".equals(data[i++]);
			exonFrames       = data[i++];
		}
		setExonsFromString(exonStarts, exonEnds, cdsOnly, exonFrames);
	}
	
	public String writeAsExonBed() {
		StringBuffer buf = new StringBuffer();
		int i = 1;
		Iterator<Exon> exonIt = exons.iterator();
		while(exonIt.hasNext()) {
			Exon exon = exonIt.next();
			buf.append(exon.getChromosomeString()).append("\t")
				.append(exon.getStart()).append("\t")
				.append(exon.getEnd()).append("\t")
				.append(getName()).append("_").append(i);
			if(exonIt.hasNext()) {buf.append("\n");}
			i++;
		}
		return buf.toString();
	}
	
	public String writeAsFullBed() {
		StringBuffer buf = new StringBuffer("chr");
		buf.append(getChromosome())
			.append("\t")
			.append(getStart())
			.append("\t")
			.append(getEnd())
			.append("\t")
			.append(getName())
			.append("\t")
			.append(getScore())
			.append("\t")
			.append(getOrientation())
			.append("\t")
			.append(getCdsStart())
			.append("\t")
			.append(getCdsEnd())
			.append("\t")
			.append("0,0,0")
			.append("\t");
		
		
			buf.append(exons.size()).append("\t");
			for(int i = 0; i < exons.size() - 1; i++) {
				buf.append(exons.get(i).getLength()).append(",");
			}
			
			if(exons.size() > 0){buf.append(exons.get(exons.size() - 1).getLength());}
			buf.append("\t");
			
			for(int i = 0; i < exons.size() - 1; i ++) {
				buf.append(exons.get(i).getStart() - getStart()).append(",");
			}
			if(exons.size() > 0){buf.append(exons.get(exons.size() - 1).getStart() - getStart());}
					
		return buf.toString();
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder("0"); //damn bin, should not include in refseq downloads.
		sb.append("\t").append(getName())
			.append("\t").append(getChromosomeString())
			.append("\t").append(getOrientation())
			.append("\t").append(getStart())
			.append("\t").append(getEnd())
			.append("\t").append(getCdsStart())
			.append("\t").append(getCdsEnd())
			.append("\t").append(exons.size()).append("\t");
		
		StringBuilder exonsStartSB = new StringBuilder();
		StringBuilder exonsEndSB = new StringBuilder();
		StringBuilder exonsFrameSB = new StringBuilder();
		Iterator<Exon> exonIt = exons.iterator();
		while(exonIt.hasNext()) {
			Exon exon = exonIt.next();
			exonsStartSB.append(exon.getStart());
			exonsEndSB.append(exon.getEnd());
			exonsFrameSB.append(exon.getFrame());
			if(exonIt.hasNext() ) {
				exonsStartSB.append(",");
				exonsEndSB.append(",");
				exonsFrameSB.append(",");
			}
		}
		sb.append(exonsStartSB.toString()).append("\t");
		sb.append(exonsEndSB.toString()).append("\t");
		sb.append(getId()).append("\t")
			.append(altName).append("\t")
			.append(completeCDSStart ? "cmpl" : "no").append("\t")
			.append(completeCDSEnd ? "cmpl" : "no").append("\t");
		sb.append(exonsFrameSB.toString());
		
		return sb.toString();
	}

	public GeneAnnotation(String geneName) {
		super();
		setName(geneName);
		exons = new ArrayList<Exon>();
	}

	public boolean containsExonById(String exonId) {
		Exon exon = getExonById(exonId);
		return exon != null;
	}
	
	public Exon getExonById(String exonId) {
		Iterator<Exon> it = exons.iterator();
		Exon exon = null;
		while(it.hasNext()) {
			exon = it.next();
			if(exonId.equals(exon.getName())) {
				return exon;
			}
		}
		return null;
		
	}
	
	public void shift(int shift) {
		setStart(getStart()  + shift);
		setEnd(getEnd() + shift);
		cdsEnd = cdsEnd + shift;
		cdsStart = cdsStart + shift;
		Iterator<Exon> it = exons.iterator();
		while(it.hasNext()) {
			Exon exon = it.next();
			exon.setStart(exon.getStart() + shift);
			exon.setEnd(exon.getEnd() + shift);
		}
	}
	
	public List<Integer[]> getCodons() {
		ArrayList<Integer[]> codons = new ArrayList<Integer[]>(getExonicLength()/3); 
		List<Exon> codingExons = getCodingExons();
		if(codingExons.size() == 0) {
			System.err.println(getName() + " with annotation " + toString() + " has no coding exons");
			return new ArrayList<Integer[]>();
		}
		Iterator<Exon> exonIt = codingExons.iterator();
		LightweightGenomicAnnotation cds = getCDS();
		Exon exon = exonIt.next();
		int genomicLocation = exon.getStart();
		//System.out.println("Gene " + toString() + "("+getOrientation()+") starting at " + genomicLocation);
		while(genomicLocation < cds.getEnd()) {	
			Integer [] codon = new Integer[3];
			for(int i = 0; i < 3; i++) {
				codon [i] = genomicLocation++;				
				if( genomicLocation > exon.getEnd() && genomicLocation < cds.getEnd() ) {
					exon = exonIt.next(); 
					genomicLocation = exon.getStart(); //+ 1; //TODO: This was a hack to test BED loads -  commented out afterwards.
				} 
			}
			codons.add(codon);
		}
		return codons;
	}
	
	public List<Integer[]> getOverlappingCodons(LightweightGenomicAnnotation annotation) {
		//TODO: this is wasteful. This should just build a new list of only overlapping exons rather than building the full list then extracting desired exons.
		List<Integer[]> allCodons = getCodons();
		List<Integer[]> overlapping = new ArrayList<Integer[]>();
		Iterator<Integer[]> codonIt = allCodons.iterator();
		while(codonIt.hasNext()) {
			Integer [] codon = codonIt.next();
			if(codon[0] >= annotation.getStart() && codon[0] < annotation.getEnd() || 
					codon[1] >= annotation.getStart() && codon[1] < annotation.getEnd() ||
					codon[2] > annotation.getStart() && codon[2] < annotation.getEnd() ) {
				overlapping.add(codon);
			}
			if(codon[0] > annotation.getEnd()) {
				break;
			}
		}
		return overlapping;
	}

	public int getCdsEnd() {
		return cdsEnd;
	}
	

	public int getCdsStart() {
		return cdsStart;
	}
	
	protected void setCdsStart(int start) { cdsStart = start;}
	protected void setCdsEnd(int end) { cdsEnd = end;}

	public int getExonCount() {
		return exonCount == 0 ? getExons().size() : exonCount;
	}
	
	/**
	 * Gets promoter bases, where promoter is taken as the bases preceding the start of the first exon.
	 * @param bases - Number of promoter bases to include
	 * @return GenomicAnnotation for the promoter
	 */
	public GenomicAnnotation getPromoter(int bases) {
		BasicGenomicAnnotation promoter = new BasicGenomicAnnotation(this);
		
		promoter.setStart(inReversedOrientation() ? getCdsEnd() : getStart() - bases);
		promoter.setEnd(inReversedOrientation() ? getEnd() + bases : getCdsStart());
		
		return promoter;
	}
	
	protected void setChromosomeFromString(String chrString) {
		setChromosome(chrString.substring(3));
	}
	protected void setExonCountFromString(String number) {
		exonCount = Integer.parseInt(number);
	}

	protected void setStartFromString(String string) {
		setStart(Integer.parseInt(string));
	}

	protected void setEndFromString(String string) {
		setEnd(Integer.parseInt(string));
	}
	
	protected void setCdsStartFromString(String string) {
		cdsStart = Integer.parseInt(string);
	}
	
	protected void setCdsEndFromString(String string) {
		cdsEnd = Integer.parseInt(string);
	}

	public void setExonsFromString(String startCoordinates, String endCoordinates, boolean cdsOnly) throws ParseException{	
		setExonsFromString(startCoordinates, endCoordinates, cdsOnly, null);
	}
	
	public void setExonsFromString(String startCoordinates, String endCoordinates, boolean cdsOnly, String exonFrames) throws ParseException{		
		int [] startCoords = parseCoordinates(startCoordinates, exonCount);
		int [] endCoords   = parseCoordinates(endCoordinates, exonCount);
		int [] frames      = null;
		if(exonFrames != null) {
			frames = parseCoordinates(exonFrames, exonCount);
		}
		
		this.exons = new ArrayList<Exon>(exonCount);
		for(int i = 0; i < exonCount; i++) {
			Exon exon = new Exon(this, getName()+"_"+i);
			exon.setStart(startCoords[i]);
			exon.setEnd(endCoords[i]);
			exon.setFrame(frames != null ? frames[i] : -1);
			exon.setOrientation(getOrientation());
			if(cdsOnly) {
				intersectCDS(exon);
			}
			if(exon.length() > 0) {
				exons.add(exon);
			}
		}
		Collections.sort(this.exons);
	}
	
	
	protected int[] parseCoordinates(String commaSeparatedCoords, int expectedNumber) throws ParseException{
		int [] coords = new int[expectedNumber];
		String [] rawCoords = commaSeparatedCoords.split(",");
		if (rawCoords.length != expectedNumber) {
			throw new ParseException("Expected number or coordinates <"+expectedNumber+"> does not match the number in data<"+rawCoords.length+">");
		}
		
		for (int i = 0; i < expectedNumber; i++ ) {
			coords[i] = Integer.parseInt(rawCoords[i]);
		}
		return coords;
	}
	
	public List<Exon> getExons() {
		return exons;
	}
	
	public List<GenomicAnnotation> getIntrons() {
		if(getExons() == null) {
			return null;
		}
		List<GenomicAnnotation> introns = new ArrayList<GenomicAnnotation>(getExons().size() - 1);
		Iterator<Exon> exonIt = getExons().iterator();
		
		Exon lastExon = exonIt.hasNext() ? exonIt.next() : null;
		int intronNum = 0;
		while(exonIt.hasNext()) {
			Exon exon = exonIt.next();
			GenomicAnnotation intron = new BasicGenomicAnnotation(this);
			intron.setName(getName() + "_intron_" + intronNum);
			intron.setStart(lastExon.getEnd());
			intron.setEnd(exon.getStart() + 1);
			lastExon = exon;
			introns.add(intron);
		}
		
		return introns;
		
	}
	
	public List<Exon> getExonsInStrandOrder() {
		List<Exon> orientedExons = exons;
		if(inReversedOrientation()) {
			orientedExons = new ArrayList<Exon>(exons.size());
			for(int i = exons.size() - 1; i >= 0; i--) {
				orientedExons.add(exons.get(i));
			}
		}
		
		return orientedExons;
	}
	
	
	
	public void intersectCDS(Exon exon) {
		exon.setStart((int)Math.max(getCdsStart(),exon.getStart()));
		exon.setEnd((int)Math.min(getCdsEnd(), exon.getEnd()));
	}
	
	public List<Exon> getCodingExons() {
		List<Exon> codingExons = new ArrayList<Exon>();
		GenomicAnnotation codingRegion = getCDS();
		Iterator<Exon> exonIt = getExons().iterator();
		while(exonIt.hasNext()) {
			Exon exon = exonIt.next();
			if(codingRegion.contains(exon)) {
				codingExons.add(exon);
			} else if(exon.overlaps(codingRegion)) {
				Exon codingExon = new Exon(this, exon);
				codingExon.takeIntersection(codingRegion);
				codingExons.add(codingExon);
			}
		}
		return codingExons;
	}
	
	
	public void resetExons() {
		this.exons = new ArrayList<Exon>();
	}
	
	public void addExon(GenomicAnnotation exon) {
		addExon(exon, 0);
	}
	
	public void addExon(GenomicAnnotation tmpExon, int frame) {
		Exon exon = new Exon(this, tmpExon);
		exon.frame = frame;
		exons.add(exon);
		Collections.sort(exons, new Comparator<Exon>() {

			public int compare(Exon arg0, Exon arg1) {
				return (int) (arg0.getStart() - arg1.getStart());
			}
			
		});
		
	}
	
	public void takeExonsUnion(List<Exon> otherExons) {
		Iterator<Exon> annotIt = otherExons.iterator();
		List<String> exonsThatOvelaped = new ArrayList<String>();
		System.out.println("First Pass at finding overlapping exons");
		while(annotIt.hasNext()) {
			Exon otherExon = annotIt.next();
			if(exonsThatOvelaped.contains(otherExon.getName())) {
				System.out.println("\texon "+otherExon+" is contained in a processed exon, skipping....");
				continue;
			}
			boolean overlapped = false;
			for(int i = 0; i < exons.size(); i++) {
				//exons are sorted, so it is safe to break if we are too far out
				Exon exon = exons.get(i);
				if (otherExon.getEnd() < exon.getStart()) {
					break;
				}
				
				if (exon.overlaps(otherExon)) {
					overlapped = true;
					exon.takeUnion(otherExon);
					exonsThatOvelaped.add(exon.getName());
					System.out.println("\t"+otherExon.toString() + " and " + exon.toString() + " overlapped, not adding.");
				}
			}
			if (!overlapped) {
				System.out.println("\tExon "+ otherExon.toString()+" did not overlap anything. Adding");
				exons.add(otherExon);
			}
		}
		
		if(exons.size() == 0) {
			return;
		}
		Collections.sort(exons);
		// Newly extended exons may actually overlap, this will merge those that do
		System.out.println("Second pass at overlapping exons");
		Stack<Exon> newList = new Stack<Exon>();
		newList.push(exons.get(0));
		exons.get(0).setName(getName() + "_" + 0); // Rename first one as it is not renamed later on.
		int j = 1;
		for(int i = 1; i < exons.size(); i++) {
			Exon lastModifiedExon = newList.pop();
			if(lastModifiedExon.overlaps(exons.get(i))) {
				lastModifiedExon.takeUnion(exons.get(i));
				System.out.println("\tExon " + exons.get(i) + " overlapped " + lastModifiedExon);
				newList.push(lastModifiedExon);
			} else {
				newList.push(lastModifiedExon);
				newList.push(exons.get(i));
				System.out.print("\tExon " + exons.get(i) + " did not overlap last exon " + lastModifiedExon);
				exons.get(i).setName(getName()+"_"+j++);
				System.out.println("... new name " + exons.get(i).getName());
				
			}
		}
		exons = newList;
		exonCount = newList.size();
	}
	
	/**
	 * Stitches together any overlapping exons.
	 */
	public void collapseExons() {
		List<Exon> collapesedExons = stitchList(exons, 0);
		exons = collapesedExons;
		exonCount = exons.size();
	}
	
	public static class Exon extends BasicGenomicAnnotation {
		private GeneAnnotation gene;
		private int frame;
		
		public Exon(GeneAnnotation gene, String name) {
			super(name);
			this.gene = gene;
		}
		
		public void setFrame(int frame) {
			this.frame = frame;
			
		}
		
		public int getFrame() { return frame;}

		public Exon(GeneAnnotation gene, BasicGenomicAnnotation exonInfo) {
			super(exonInfo.getName());
			if(!gene.getChromosome().equals(exonInfo.getChromosome())) {
				throw new IllegalArgumentException("Traying to add an exon " + exonInfo + " to gene " + gene +" which are on different chromosomes");
			} else if (!gene.contains(exonInfo)) {
				throw new IllegalArgumentException("Traying to add an exon " + exonInfo + " to gene " + gene +" but gene does not contain it");
			}
						
			this.gene = gene;
			setStart(exonInfo.getStart());
			setEnd(exonInfo.getEnd());
			setChromosome(exonInfo.getChromosome());
			setReversedOrientation(exonInfo.inReversedOrientation());
		}

		
		public Exon(GeneAnnotation gene, GenomicAnnotation block) {
			super(block);
			this.gene = gene;
		}

		public int length() {
			return getEnd() - getStart() ;
		}

		public int compareTo(Exon arg0) {
			long diff = getStart() == arg0.getStart() ? arg0.getEnd() - getEnd()  : getStart() - arg0.getStart();
			int result = 0;

			// Lets avoid weird cast errors 
			if (diff < Integer.MIN_VALUE) {
				result = Integer.MIN_VALUE;
			} else if(diff > Integer.MAX_VALUE) {
				result = Integer.MAX_VALUE;
			} else {
				result = (int) diff;
			}
			return result;
		}

		public String getChromosome() {
			return gene.getChromosome();
		}

		public void extendToIntron(int includeIntronicRegionSize) {
			super.setThreePrimeBuffer(includeIntronicRegionSize);
			super.setFivePrimeBuffer(includeIntronicRegionSize);
		}
		
		public GeneAnnotation getGene() {
			return gene;
		}

		public void setGene(GeneAnnotation gene) {
			this.gene = gene;
			
		}
	}

	public void extendExons(int includeIntronicRegionSize) {
		Iterator<Exon> it = exons.iterator();
		Exon exon = null;
		while(it.hasNext()) {
			exon = it.next();
			exon.extendToIntron(includeIntronicRegionSize);
		}
		
	}

	public String getAltName() {
		return altName == null || altName.length() == 0 ? getName() : altName;
	}

	public boolean isCompleteCDSEnd() {
		return completeCDSEnd;
	}

	public boolean isCompleteCDSStart() {
		return completeCDSStart;
	}
	
	public int getExonicLength() {
		int length = 0;
		Iterator<Exon> it = exons.iterator();
		while(it.hasNext()) {
			length += it.next().length();
		}
		return length;
	}

	public GenomicAnnotation getCDS() {
		BasicGenomicAnnotation cds = new BasicGenomicAnnotation(this);
		cds.setStart(getCdsStart());
		cds.setEnd(getCdsEnd());
		return cds;
	}
	
	public void addBlock(String name, int start, int end) {
		Exon exon = new Exon(this, name);
		exon.setStart(start);
		exon.setEnd(end);
		addExon(exon);
	}
	
	public List<Exon>  getBlocks() {return getExons();}

	public boolean mayHaveBlocks() {
		return true;
	}

}
