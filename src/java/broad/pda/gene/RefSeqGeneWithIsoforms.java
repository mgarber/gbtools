package broad.pda.gene;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import org.broad.igv.sam.Alignment;

import broad.core.annotation.GFF;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.util.CollapseByIntersection;
import broad.pda.datastructures.Alignments;


/* This class was generated to resolve the problem of treating isoforms with the same exact
stat-end coordinates while storing genes in a sorted manner in an intervalTree.
The class extends RefSeqGene.
The fields inherited by RefSeqGene will represent the isoform with the longest transcript.
Isoforms may be with different orientation
 */
public class RefSeqGeneWithIsoforms extends RefSeqGene{

	Collection <RefSeqGene> isoforms;
	int numOfIsoforms;
	
	
	//Ctor:
	public RefSeqGeneWithIsoforms(Alignment align) {
		super(align);
		this.numOfIsoforms=1;  this.isoforms=new ArrayList <RefSeqGene>();
	}
	public RefSeqGeneWithIsoforms(Collection<? extends LightweightGenomicAnnotation> exons){
		super(exons);
		this.numOfIsoforms=1; this.isoforms=new ArrayList <RefSeqGene>();
	}
	public RefSeqGeneWithIsoforms(LightweightGenomicAnnotation align){
		super(align);
		this.numOfIsoforms=1; this.isoforms=new ArrayList <RefSeqGene>();
	}
	public RefSeqGeneWithIsoforms(String pslString){
		super(pslString);
		this.numOfIsoforms=1; this.isoforms=new ArrayList <RefSeqGene>();
	}
	public RefSeqGeneWithIsoforms(String chrString, int start, int end){
		super(chrString, start,end);
		this.numOfIsoforms=1; this.isoforms=new ArrayList <RefSeqGene>();
	}
	public RefSeqGeneWithIsoforms(String chr, int start, int end, String name, String orientation, Collection<Alignments> exons){
		super(chr,start,end,name,orientation,exons);
		this.numOfIsoforms=1; this.isoforms=new ArrayList <RefSeqGene>();
	}
	public RefSeqGeneWithIsoforms(String chr, int start, int end, String name, String orientation, int[] exonsStart, int[] exonsEnd){
		super(chr,start,end,name,orientation,exonsStart,exonsEnd);
		this.numOfIsoforms=1; this.isoforms=new ArrayList <RefSeqGene>();
	}
	public RefSeqGeneWithIsoforms(String chr, int start, int end, String name, String orientation, List<Integer> exonsStart, List<Integer> exonsEnd){
		super(chr,start,end,name,orientation,exonsStart,exonsEnd);
		this.numOfIsoforms=1; this.isoforms=new ArrayList <RefSeqGene>();
	}
	public RefSeqGeneWithIsoforms(String chr, int start, int end, String name, String orientation, List<Integer> exonsStart, List<Integer> exonsEnd, int blockStart, int blockEnd){
		super(chr,start,end,name,orientation,exonsStart,exonsEnd,blockStart,blockEnd);
		this.numOfIsoforms=1; this.isoforms=new ArrayList <RefSeqGene>();
	}
	
	public RefSeqGeneWithIsoforms(String chr, int start, int end, String name, String orientation, List<Integer> exonsStart, List<Integer> exonsEnd, int blockStart, int blockEnd,String [] extraData){
		super(chr,start,end,name,orientation,exonsStart,exonsEnd,blockStart,blockEnd, extraData);
		this.numOfIsoforms=1; this.isoforms=new ArrayList <RefSeqGene>();
	}
	
	public RefSeqGeneWithIsoforms(String chr, int start, int end, String name, String orientation, List<Integer> exonsStart, List<Integer> exonsEnd, String [] extraData){
		super(chr,start,end,name,orientation,exonsStart,exonsEnd, extraData);
		this.numOfIsoforms=1; this.isoforms=new ArrayList <RefSeqGene>();
	}
		
	
	public RefSeqGeneWithIsoforms(RefSeqGene gene) {
		//super (gene.toBED(),false);
		super(gene);
		this.numOfIsoforms=1; 
		isoforms=new ArrayList <RefSeqGene>();
		setSequence(gene.getSequence());
		setCountScore(gene.getCountScore());
		setBedScore(gene.getBedScore());
		
		if (gene.getExons()!=null) {this.setExons(gene.getExons(), gene.getExonsScores());}
		
		this.numOfIsoforms=1; this.isoforms=new ArrayList <RefSeqGene>();
	}
	public boolean IsExactOverlappingIsoform(RefSeqGene gene){
		return gene.getChr().equalsIgnoreCase(this.getChr()) && gene.getStart()== this.getStart() && gene.getEnd()==this.getEnd();
	}
	
	//Add an isoform
	public boolean AddIsoform (RefSeqGene gene) {
		return AddIsoform(gene, true);
	}
	
	public boolean AddIsoform (RefSeqGene gene, boolean updateMainIsoform){
		
		if (this.IsExactOverlappingIsoform(gene)){
			if (!updateMainIsoform || gene.getTranscriptLength() <= this.getTranscriptLength()){
				this.isoforms.add(gene);
				this.numOfIsoforms++;
			}
			else{
				RefSeqGene currGene= this.makeRefSeqGeneInstance();
				this.isoforms.add(currGene);
				this.numOfIsoforms++;
				this.updateMainIsoform(gene);
				
			}
			return true;
		} else {
			System.err.println("Could not add isoform " + gene.toUCSC() + " to iso set with main iso: " + toUCSC());
		}
		return false;
	}
	
	
	public boolean addContainedIsoform (RefSeqGene isoform){
		boolean okToAdd = contains(isoform) && overlaps(isoform);
		if(okToAdd) {

			isoforms.add(isoform);
			numOfIsoforms++;
		}
		
		return okToAdd;
	}
	
	public boolean AddAllIsoforms(RefSeqGeneWithIsoforms overlapper) {
		return AddAllIsoforms(overlapper, true);
	}
	
	public boolean AddAllIsoforms(RefSeqGeneWithIsoforms overlapper, boolean updateMainIsoform) {
		boolean res=true;
		Collection<RefSeqGene> set=overlapper.getAllIsoforms();
		if (set.isEmpty()) return false;
		for (RefSeqGene G : set){
			res=!res ? false : this.AddIsoform(G, updateMainIsoform);
			//if(!res) return false;
		}
		return res;
		
	}
	
	public boolean AddAllIsoNotIncludedforms(RefSeqGeneWithIsoforms other) {
		boolean res=true;
		Collection<RefSeqGene> otherIsos=other.getAllIsoforms();
		Collection<RefSeqGene> thisIsos = getAllIsoforms();
		if (otherIsos.isEmpty()) return false;
		for (RefSeqGene G : otherIsos){
			boolean contained = false;
			Iterator<RefSeqGene> thisIsoIt = thisIsos.iterator();
			while(!contained && thisIsoIt.hasNext()) {
				RefSeqGene iso = thisIsoIt.next();
				contained = iso.compareTo(G) == 0;
			}
			if(!contained) {
				res=!res ? false : this.AddIsoform(G, false);
			}
		}
		return res;
	}
	
	public Collection<RefSeqGene> getAllIsoforms() {

		return getAllIsoforms(true);
	}
	
	/**
	 * In cases were the containg gene is artificial, like when building gene loci and adding
	 * all isoforms contain in the loci. The containing gene is not desired. This methods allows
	 * the caller to exclude the containing isoform from the list.
	 * @param includeContainingIsoform - boolean flag indicating whether or not to include the containing isoform
	 * @return Collection<RefSeqGene> of all the isoforms for this gene
	 */
	public Collection<RefSeqGene> getAllIsoforms(boolean includeContainingIsoform) {
		
		//RefSeqGene currGene= makeRefSeqGeneInstance();
		Collection <RefSeqGene> rtrn= new ArrayList<RefSeqGene>();
		if (! this.isoforms.isEmpty())
			rtrn.addAll(this.isoforms);
		if(includeContainingIsoform) {
			rtrn.add(this);
		}
		return rtrn;
	}
	
	//This function assigns the same bed score to all isoforms
	public void setBedScore(double scr){
		this.bedScore=scr;
		if (! this.isoforms.isEmpty()){
			for (RefSeqGene iso:this.isoforms)
				iso.setBedScore(scr);
		}
	}
	
	public void setExtraFields(double [] scores) {
		super.setExtraFields(scores);
		/**
		 * Commented out by skadri on 04/02/12
		 * Replaces extra fields for all children of a gene. 
		 */
		/*if (! this.isoforms.isEmpty()){
			for (RefSeqGene iso: isoforms)
				iso.setExtraFields(scores);
		}*/
	}
	
	//This function assigns the update the count score to equal the bed score in all isoforms
	public void updateScrToBedScore(){
		super.updateScrToBedScore();
		for (RefSeqGene Iso:this.isoforms)
			Iso.updateScrToBedScore();
	}
	
	public List<LightweightGenomicAnnotation> getScoredExons(){
		List<LightweightGenomicAnnotation> allExons= new ArrayList<LightweightGenomicAnnotation>();
		for (RefSeqGene Iso:this.getAllIsoforms()){
			allExons.addAll(Iso.getScoredExons());
		}
		return allExons;
	}
	//Get Merged Transcript
	
	
	private void updateMainIsoform(RefSeqGene gene) {

		chr=gene.getChr();
		start=gene.getStart();
		stop=gene.getEnd();
		name=gene.getName();
		orientation=gene.getOrientation();
		sequence=gene.getSequence();
		bedScore=gene.getBedScore(); //the transcript score as it appears in a bed file
		extraFields=gene.getExtraFields();
		setExons(gene.getExons(),gene.getExonsScores());
		
		
	}
	
	private RefSeqGene makeRefSeqGeneInstance() {
		
		//return new RefSeqGene(this.getChr(), this.getStart(), this.getEnd(), this.getName(),this.getCountScore(),this.getBedScore(), this.getOrientation(), this.exonStarts, this.exonEnds,this.exonScores,this.sequence);
		//RefSeqGene g=new RefSeqGene(toBED(),false);
		RefSeqGene g=new RefSeqGene(this);
		g.setSequence(this.sequence);
		g.setBedScore(this.getBedScore());
		g.setCountScore(this.getCountScore());
		return  g;
		
	}
	
	public boolean overlapsByGenomicRegion(
			IntervalTree<RefSeqGeneWithIsoforms> otherTree, boolean considerOrientation) {
			
			Iterator<Node< RefSeqGeneWithIsoforms>> geneIt=otherTree.overlappers(this.getStart(),this.getEnd());
			while (geneIt.hasNext()){
				RefSeqGene otherGene=geneIt.next().getValue();
				if (considerOrientation && this.overlapsGene(otherGene))
					return true;
				if (!considerOrientation && this.overlapsGeneInAnyOrientation(otherGene))
					return true;
			}
			return false;
	}

	public void dedup() {
		isoforms = new TreeSet<RefSeqGene>(isoforms);
	}

	public RefSeqGene findCompatibleGenes(IntervalTree<RefSeqGeneWithIsoforms> overlapTree, int[] numIntrons) {
		
		RefSeqGene rtrn=null;
		int bestIntronNum=0;
		Iterator <RefSeqGeneWithIsoforms> gIt=overlapTree.valueIterator();
		while(gIt.hasNext()){
			Collection<RefSeqGene> genes=gIt.next().getAllIsoforms();
			for (RefSeqGene g: genes){
				Collection<RefSeqGene> myIso=this.getAllIsoforms();
				for (RefSeqGene iso:myIso){
					int intronNum=iso.numOfCompatibleIntrons(g);
					if (intronNum > bestIntronNum){
						bestIntronNum=intronNum;
						rtrn=g;
					}
				}
			}
		}
		
		numIntrons[0]=bestIntronNum;
		return rtrn;
	}
	public void expandUtrs(Integer utr1, Integer utr2) {
		super.expandUtrs(utr1,utr2);
		for(RefSeqGene g: this.isoforms){
			g.expandUtrs(utr1,utr2);
		}
		
	}
	
	public void set3PrimeEnd(int updated3Prime) {
		super.set3PrimeEnd(updated3Prime);
		for(RefSeqGene g: this.isoforms){
			g.set3PrimeEnd(updated3Prime);
		}
	}
	
	public void set5PrimeEnd(int updated5Prime) {
		super.set5PrimeEnd(updated5Prime);
		for(RefSeqGene g: this.isoforms){
			g.set5PrimeEnd(updated5Prime);
		}
	}
	
	
	//adds a suffix to all the isoforms
	public void addSuffixToName(String refName) {

		this.setName(this.getName()+refName);
		for (RefSeqGene g:this.isoforms){
			g.setName(g.getName()+refName);
		}
		
	}
	
	
	
	public boolean overlapsExon(LightweightGenomicAnnotation exon){
		boolean overlaps = super.overlapsExon(exon);
		
		Iterator<RefSeqGene> isoformIt = isoforms.iterator();
		while(!overlaps && isoformIt.hasNext()) {
			overlaps = isoformIt.next().overlapsExon(exon);
		}
		
		return overlaps;
	}
	
	
	
	public boolean overlapsExon(RefSeqGene gene) {
		boolean overlaps = super.overlaps(gene);
		
		Iterator<RefSeqGene> isoformIt = isoforms.iterator();
		while(!overlaps && isoformIt.hasNext()) {
			overlaps = isoformIt.next().overlapsExon(gene);
		}
		
		return overlaps;
	}
	
	
	public boolean overlaps(RefSeqGene other) {
		boolean overlaps = super.overlaps(other);
		
		Iterator<RefSeqGene> isoformIt = isoforms.iterator();
		while(!overlaps && isoformIt.hasNext()) {
			RefSeqGene isoform = isoformIt.next(); 
			if(!isoform.equals(this)) {
				overlaps = isoform.overlaps(other);
			}
		}
		
		return overlaps;
	}
	
	public RefSeqGene getMerged() {

		RefSeqGene curr=this;
		for (RefSeqGene iso: this.isoforms)
				curr=curr.takeUnion(iso);
		return curr;
	}

	
public void cleanIsoforms() {
		isoforms = new ArrayList<RefSeqGene>();
		numOfIsoforms = 1;
	}
	
	public Set<Alignments> getExonSet(){
		Set<Alignments> rtrn=new TreeSet<Alignments>();
		if (exonStarts != null){
			for(int i=0; i<exonStarts.length; i++){
				Alignments exon = new Alignments(chr, exonStarts[i], exonEnds[i]);
				exon.setOrientation(getOrientation());
				exon.setName(getName()+"_"+i);
				rtrn.add(exon);
			}
		} else {
			Alignments exon = new Alignments(chr, getStart(), getEnd());
			rtrn.add(exon);
		}
		return rtrn;
	}
	
	public Alignments[] getIntronsBlocks(){
		
		Object[] sortedUniqueExons=this.getSortedAndUniqueExons().toArray();
		if(sortedUniqueExons.length==0){return new Alignments[0];}
		Alignments[] rtrn=new Alignments[sortedUniqueExons.length-1];
		
		
		for(int i=0; i<sortedUniqueExons.length-1; i++){
			Alignments current=(Alignments)sortedUniqueExons[i];
			Alignments next=(Alignments) sortedUniqueExons[i+1];
			Alignments align=new Alignments(chr, current.getEnd(), next.getStart());
			rtrn[i]=align;
		}
		return rtrn;
	}
	
	public Collection<Alignments> getAllIsoformExonSet(){
		Collection<Alignments> rtrn=new TreeSet<Alignments>(); 
		
		for(RefSeqGene isoform : isoforms) {
			rtrn.addAll(isoform.getExonSet());
		}

		return rtrn;
	}
	
	public Collection<Alignments> getAllIsoformIntronSet(){
		Collection<Alignments> rtrn=new TreeSet<Alignments>(); 
		
		for(RefSeqGene isoform : isoforms) {
			rtrn.addAll(isoform.getIntronSet());
		}

		return rtrn;
	}
	
	public RefSeqGene constituentIsoform() {
		Collection<Alignments> exons = getAllIsoformExonSet();
		exons = CollapseByIntersection.collapseByIntersection(exons, true);
		Set<Alignments> constituentExons = new TreeSet<Alignments>();
		Collection<RefSeqGene> isoforms = getAllIsoforms();
		for(Alignments exon : exons) {
			Iterator<RefSeqGene> isoformIt = isoforms.iterator();
			boolean overlaps = true;
			while(overlaps &&  isoformIt.hasNext()) {
				RefSeqGene isoform = isoformIt.next();
				overlaps = isoform.overlapsExon(exon);
				//System.err.println("exon " + exon.toUCSC() +" overlaps " + isoform.toUCSC() + "("+isoform.getName()+")? " + overlaps);
			}
			if(overlaps) {
				//System.err.println("Adding exon " + exon.toUCSC());
				constituentExons.add(new Alignments(exon));
			}
		}
		
		
		RefSeqGene constituentIsoform = null;
		
		if(constituentExons.size() >0) {
			constituentIsoform = new RefSeqGene(constituentExons);
			constituentIsoform.setName(getName());
			constituentIsoform.setBedScore(getBedScore());
			constituentIsoform.setOrientation(getOrientation());
		}
		return constituentIsoform;
	}
	
	/**
	 * Constructs an artificial gene where the constituent introns (regions that are intronic in all isoforms) are the exons of the gene.
	 * @return
	 */
	public RefSeqGene constituentIntrons() {
		Collection<Alignments> introns = getAllIsoformIntronSet();
		introns = CollapseByIntersection.collapseByIntersection(introns, true);
		Set<Alignments> constituentIntrons = new TreeSet<Alignments>();
		Collection<RefSeqGene> isoforms = getAllIsoforms();
		for(Alignments intron : introns) {
			Iterator<RefSeqGene> isoformIt = isoforms.iterator();
			boolean isContained = true;
			while(isContained &&  isoformIt.hasNext()) {
				RefSeqGene isoform = isoformIt.next();
				Iterator<Alignments> isoIntronIt = isoform.getIntronSet().iterator();
				boolean containedInIsoformIntron = false;
				while(!containedInIsoformIntron && isoIntronIt.hasNext()) {
					Alignments isoIntron = isoIntronIt.next();
					containedInIsoformIntron = isoIntron.contains(intron);
				}
				isContained = containedInIsoformIntron;
				//System.err.println("exon " + exon.toUCSC() +" overlaps " + isoform.toUCSC() + "("+isoform.getName()+")? " + overlaps);
			}
			if(isContained) {
				constituentIntrons.add(new Alignments(intron));
			}
		}
		
		RefSeqGene constituentIntroform = new RefSeqGene(constituentIntrons);
		constituentIntroform.setName(getName());
		constituentIntroform.setBedScore(getBedScore());
		constituentIntroform.setOrientation(getOrientation());
		
		return constituentIntroform;
	}
	
	public String toMISO(String source){
		StringBuilder rtrn = new StringBuilder();
		GFF geneGFF = new GFF(getName());
		geneGFF.setSource(source);
		geneGFF.setOrientation(this.getOrientation());
		geneGFF.setFeature("gene");
		
		geneGFF.addAttribute("ID", getName());//As in version 1353 to compile with cufflinks
		geneGFF.addAttribute("Name", getName());//As in version 1353 to compile with cufflinks

		//Add +1 when we print GFF as it is 1 based rather than 0 based
		geneGFF.setStart(getStart());  //BUG FIX : MORAN AUG 17TH, ADD +1; 2nd bug fix - 11.29.10 only during print we add +1 to start pos
		geneGFF.setEnd(getEnd()); //BUG FIX : MORAN AUG 17TH, END IS CORRECT AS GFF IS INCLUSIVE BUT ALSO USES 1 BASE CORRDINATES (WE USE 0 BASED)
		geneGFF.setChromosome(getChr());
		rtrn.append(geneGFF.toString(false, true));
		int i = 0;
		Collection<RefSeqGene> isoforms = getAllIsoforms();
		for (RefSeqGene iso : isoforms) {
			rtrn.append("\n");
			i++;
			rtrn.append(iso.toMISO(getName(), iso.getName()+"-" +i, source));
		}

		return rtrn.toString();
	}
	public void colapse(int exonFudgeFactor) {
		List<RefSeqGene> collapsedIsoforms = new ArrayList<RefSeqGene>();
		Iterator<RefSeqGene> isoIt = isoforms.iterator(); 
		while ( isoIt.hasNext()) {
			RefSeqGene iso = isoIt.next();
			boolean foundAlmostEqual = this.almostEqual(iso, exonFudgeFactor);
			Iterator<RefSeqGene> collapsedIsoIt = collapsedIsoforms.iterator();
			while(!foundAlmostEqual && collapsedIsoIt.hasNext()) {
				foundAlmostEqual = collapsedIsoIt.next().almostEqual(iso, exonFudgeFactor);
			}
			if(!foundAlmostEqual) {
				collapsedIsoforms.add(iso);
			}
		}
		
		isoforms = collapsedIsoforms;
		/*if(numOfIsoforms != isoforms.size()+1) {
			System.err.println("Collapsed isoforms for " + getName() + " from  " + numOfIsoforms + " to " + (isoforms.size()+1));
		}*/
		numOfIsoforms = collapsedIsoforms.size() + 1; //TODO: why do we have this field? it is completely redundant.
		
	}
	
	public boolean almostContains(RefSeqGene other, int exonFudgeFactor) {
		boolean isAlmostContained = false;
		//System.err.println("Checking whether refseqGene " + other.toBED() + " is contained in refseqgenewithisos " + toBED() );
		Iterator<RefSeqGene> isoIt = getAllIsoforms().iterator();
		while(!isAlmostContained && isoIt.hasNext()) {
			RefSeqGene iso = isoIt.next();
			isAlmostContained = iso.almostContainsRefSeq(other, exonFudgeFactor);
		}
		
		return isAlmostContained;
	
	}
	
	public boolean almostContains(RefSeqGeneWithIsoforms other, int exonFudgeFactor) {
		boolean isAlmostContained = false;
		Iterator<RefSeqGene> otherIsoIt = other.getAllIsoforms().iterator();
		while(!isAlmostContained && otherIsoIt.hasNext()) {
			RefSeqGene otherIso = otherIsoIt.next();
			isAlmostContained = almostContains(otherIso, exonFudgeFactor);
		}
		
		return isAlmostContained;
	
	}
	
	public boolean almostContainsStructure(RefSeqGene other, int exonFudgeFactor) {
		boolean isAlmostContained = false;
		//System.err.println("Checking whether refseqGene " + other.toBED() + " is contained in refseqgenewithisos " + toBED() );
		Iterator<RefSeqGene> isoIt = getAllIsoforms().iterator();
		while(!isAlmostContained && isoIt.hasNext()) {
			RefSeqGene iso = isoIt.next();
			isAlmostContained = iso.almostContainsStructureRefSeq(other, exonFudgeFactor);
		}
		
		return isAlmostContained;
	}
	
	
	public boolean almostContainsStructure(RefSeqGeneWithIsoforms other, int exonFudgeFactor) {
		boolean isAlmostContained = false;
		Iterator<RefSeqGene> otherIsoIt = other.getAllIsoforms().iterator();
		while(!isAlmostContained && otherIsoIt.hasNext()) {
			RefSeqGene otherIso = otherIsoIt.next();
			isAlmostContained = almostContainsStructure(otherIso, exonFudgeFactor);
		}
		
		return isAlmostContained;
	}
	
	
	public boolean hasOverlapping3PUTR(RefSeqGeneWithIsoforms other) {
		boolean overlaps = false;
		if("*".equals(getOrientation()) || "*".equals(other.getOrientation())) {
			overlaps = false;
		} else {
			overlaps =  get3PrimeExon().overlaps(other.get3PrimeExon());
		}
		return overlaps;
	}
	
	public boolean hasOverlapping5PUTR(RefSeqGeneWithIsoforms other) {
		return get5PrimeExon().overlaps(other.get5PrimeExon());
	}

	public String allIsoformsToBED() {
		StringBuilder sb = new StringBuilder();
		for(RefSeqGene iso : getAllIsoforms()) {
			sb.append(iso.toBED());
			sb.append("\n");
		}
		return sb.toString();
	}
	public void standardizeOreintation() {
		int plus = 0;
		int minus = 0;
		for(RefSeqGene iso : getAllIsoforms()) {
			if ("+".equals(iso.getOrientation())) {
				plus++;
			} else if ("-".equals(iso.getOrientation())) {
				minus++;
			} 
		}
		
		if(plus > minus) {
			for(RefSeqGene iso : getAllIsoforms()) {
				iso.setOrientation("+");
			}
		} else 	if(minus > plus) {
			for(RefSeqGene iso : getAllIsoforms()) {
				iso.setOrientation("-");
			}
		}else{
			for(RefSeqGene iso : getAllIsoforms()) {
				iso.setOrientation(RefSeqGene.UNORIENTED_STRING);
			}
		}
		
	}
	
 	public String toString(String name){
		String rtrn="";
		Collection<RefSeqGene> genes=this.getAllIsoforms();
		for(RefSeqGene gene: genes){
			gene.setName(name);
			rtrn+=gene.toBED()+"\n";
		}
		return rtrn;
 	}
	public Integer getNumberOfIsoforms() {
		
		return this.isoforms.size();
	}
	
	
	public Collection<? extends RefSeqGene> selectRandIsoSubset(
			Integer maxIsoPerLoci) {
		
		Collection<RefSeqGene> rtrn = new TreeSet <RefSeqGene>();
		int i=0;
		for (RefSeqGene g: this.getAllIsoforms()){
			if (i< maxIsoPerLoci){
				rtrn.add(g); i++;
			}
			else
				break;
		}
			
		return rtrn;
	}


}
