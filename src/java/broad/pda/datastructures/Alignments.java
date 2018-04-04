package broad.pda.datastructures;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.broad.igv.sam.Alignment;

import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.math.Statistics;
import broad.pda.chromosome.Chromosome;
import broad.pda.gene.RefSeqGene;
import broad.pda.rnai.ExtractSequence;


public class Alignments extends BasicGenomicAnnotation {

	int id = -1;
	List<Double> scores;
	String type;
	
	//Added for SHRiMP output
	String editString;
	int shrimpScore;
	int readCount;
	private double countScore;
	int startIndex;
	int endIndex;
	String line;
	
	List<Integer> countScores;
	
		
	public Alignments(String ucsc){
		String chr=ucsc.split(":")[0];
		String start=ucsc.split(":")[1].split("-")[0];
		String end=ucsc.split(":")[1].split("-")[1];
		setChromosome(chr);
		setStart(new Integer(start));
		setEnd(new Integer(end));
	}

	public Alignments(String chr, int start, int end){
		super("", chr, start, end);
	}
	
	public Alignments(String chr, int start, int end, String strand){
		this(chr, start, end);
		super.setOrientation(strand);
	}
	
	public Alignments(String chr, String start, String end){
		this(chr, Integer.parseInt(start), Integer.parseInt(end));
	}

	public Alignments(String chr, int start, int stop, ContinuousData data){
		this.startIndex=start;
		this.endIndex=stop;
		
		setStart(data.getLocations(chr)[startIndex]);
		setEnd(data.getLocations(chr)[endIndex]);
		setChromosome(chr);
		
	}

	public Alignments(String chr, int start, int end, String strand, String type){
		this(chr, start, end, strand);
		this.type=type;
	}
	
	public Alignments(String fileName, boolean fromFileName){
		String chr=fileName.split("_")[0];
		String start=fileName.split("_")[1].split("-")[0];
		String end=fileName.split("_")[1].split("-")[1];
		setChromosome(chr);
		setStart(new Integer(start));
		setEnd(new Integer(end));
	}
	
	public Alignments(LightweightGenomicAnnotation diffElement) {
		super(diffElement);
	}

	
	public Alignments(Alignment record) {
		this(record.getChromosome(), record.getStart(), record.getEnd());
	}
	
	public static Set<Alignments> getAnchorPoints(Collection<Alignments> col1, Collection<Alignments> col2){
		Set<Alignments> anchors=new TreeSet<Alignments>();
		for(Alignments align: col1){
			if(col2.contains(align)){anchors.add(align);}
		}
		return anchors;
	}

	public static Set<Alignments> getAnchorPoints(Collection[] cols){
		Set anchors=new TreeSet();
		Set allAligns=new TreeSet();
		for(int i=0; i<cols.length; i++){
			allAligns.addAll(cols[i]);
		}
		for(Object align: allAligns){
			int count=0;
			for(int i=0;i<cols.length; i++){
				if(cols[i].contains(align)){count++;}
			}
			if(count>=2){anchors.add(align);}
		}
		return anchors;
	}

	/*************Might be bad************/
	public String getChr(){return getChromosome();}

	private Integer getChrNum(String chr){
		int rtrn=25;
		String num=chr.replaceAll("chr", "");
		if(num.equalsIgnoreCase("X")){rtrn=23;}
		else if(num.equalsIgnoreCase("Y")){rtrn=24;}
		else{rtrn=new Integer(num);}
		return rtrn;
	}

	public double getCountScore(){return this.countScore;}

	public String getEditString(){
		return this.editString;
	}

	public int getEndIndex() {
		return endIndex;
	}

	public int getID(){
		return id;
	}

	/**
	 * Intersect this region with a collection of regions
	 * @param exons the collection of regions
	 * @return A collection consisting of the overlaps with each region. Only includes nonempty overlaps.
	 */
	public Collection<Alignments> getIntersections(Collection<Alignments> exons) {
		Collection<Alignments> overlaps = new ArrayList<Alignments>();
		for(Alignments exon : exons) {
			if(this.getChr().equals(exon.getChr()) && this.getStart() < exon.getEnd() && this.getEnd() > exon.getStart()) {
				Alignments b  = new Alignments(this.getChr(), Math.max(getStart(),exon.getStart()), Math.min(getEnd(), exon.getEnd()));
				overlaps.add(b);
			}
		}
		return overlaps;
	}

	public double getMeanScore() {
		return scores != null ? Statistics.mean(scores) : 0;
	}

	public double getMedianScore() {
		return scores != null ? Statistics.median(scores) : 0;
	}

	public double getNumberPositions(int rate){
		return ((this.getEnd()-this.getStart())/rate)+1;
		
	}

	public int getReadCount(int count){
		return this.readCount;
	}

	/*public double getScore(){
		double sum=0;
		if(scores != null) {
			for(int i=0; i<scores.size(); i++){
				sum+=scores.get(i);
			}
		}
		return  sum;
	}*/

	public List<Double> getScores() { return this.scores;}

	public String getSequence(Chromosome chr, boolean repeatMask) throws Exception {
		return ExtractSequence.getSequenceUnoriented(this, chr, repeatMask, new TreeMap());
	}

	public int getShrimpScore(){
		return this.shrimpScore;
	}

	public int getSize(){
		return getEnd()-this.getStart();
	}

	public int getStartIndex() {
		return startIndex;
	}

	public String getStrand(){return getOrientation();}

	public String getString(){
		return toUCSC();
	}

	public String getType(){return this.type;}

	//get unique (non-overlapping) length of both sides
	//Returns a negative number if they dont overlap
	public int getUniqueLength(Alignments exon) {
		if(!overlaps(exon)){return -99;}
		int rightDistance=exon.getEnd()-this.getEnd();
		int leftDistance=this.getStart()-exon.getStart();
		return Math.max(rightDistance, leftDistance);
	}
	
	public List<Integer> getCountScores() { return countScores; }
	public void setCountScores(List<Integer> countScores) {this.countScores = countScores;}

	public void setCountScore(double score){this.countScore=score;}

	public void setEndIndex(int endIndex) {
		this.endIndex = endIndex;
	}

	public void setID(int ID){
		this.id = ID;
	}

	public void setLine(String line){this.line=line;}

	public void setScores(List<Double> values) { this.scores = values;}

	public void setStartIndex(int startIndex) {
		this.startIndex = startIndex;
	}

	public void setStrand(String strand) {
		setOrientation(strand);
		
	}

	public void addScore(double score){
		if(this.scores==null){this.scores=new ArrayList<Double>();}
		this.scores.add(score);
	}
	
	//SHRiMP-specific (LAG 6-17-06)
	public void addShrimpScore(int ss){
		this.shrimpScore=ss;
	}
	
	public void addEditString(String edits){
		//differences between read sequence in colorspace and genomic alignment in DNA space (e.g. "5x11x5x3")
		this.editString = edits;
	}
	
	public void addReadCount(int count){
		this.readCount = count;
	}
	
	public int numberOfScores(){return this.scores.size();}
	
	/*public static Alignments collapseSetOfAlignments(Collection<Alignments> overlappingAlignments){
		Object[] array=overlappingAlignments.toArray();
		Arrays.sort(array);
		Alignments start=(Alignments)array[0];
		Alignments end=(Alignments)array[array.length-1];
		return new Alignments(start.chr, start.start, end.end);
	}*/
	
	public boolean fullyContained(Alignments g){
		if(!g.getChr().equals(getChr())){return false;}
		if(g.getStart()>=getStart() && g.getEnd()<=getEnd()){return true;}
		return false;
	}
	
	public boolean equals(Object o){
		Alignments a=(Alignments) o;
		
		if(a.getChr().equalsIgnoreCase(getChr()) && a.getStart()==getStart() && a.getEnd()==getEnd()){return true;}
		
		if(!a.getChr().equalsIgnoreCase(this.getChr())){return false;}
		if(a.getStart()!=this.getStart()){return false;}
		if(a.getEnd()!=getEnd()){return false;}
		return true;
	}
	
	
	
	/*public int compareTo(Object b){
		Alignments al=(Alignments)b;
		if(!al.chr.equalsIgnoreCase(chr)){return al.chr.compareTo(chr);}
		return -(new Integer(al.getStart()).compareTo(new Integer(start)));
		
	}*/
	
	
	/**
	 * First compare chr, then start, then end
	 */
	public int compareTo(LightweightGenomicAnnotation al){
		/*;
		if(this.equals(b)){return 0;}
		if(!b.getChromosome().equalsIgnoreCase(getChr())){return b.getChromosome().compareTo(getChr());}
		int mid1=(b.getStart()+b.getEnd())/2;
		int mid2=(getStart()+getEnd())/2;
		
		return -(new Integer(mid1).compareTo(new Integer(mid2)));*/
		
		//Alignments al=(Alignments)b;
		if(!al.getChromosomeString().equalsIgnoreCase(getChromosomeString())){return getChromosomeString().compareTo(al.getChromosomeString());}
		//int mid1=(al.start+al.stop)/2;
		//int mid2=(start+stop)/2;
		
		//return mid1-mid2;
		int counter1=al.getStart()-getStart();
		counter1 = counter1==0  ? al.getEnd()- getEnd() : counter1;
		return -counter1;
		
	}
	
	public Alignments extendRHS(int extension, Map dataMap){
		double[] extensionScores=new double[extension];
		int counter=0;
		for(int i=0; i<extension; i++){
			Object o=dataMap.get(getEnd()+(counter*25));
			while(o==null){
				counter++;
				o=dataMap.get(getEnd()+(counter*25));
			}
			extensionScores[i]=(Double)o;
			counter++;
		}
		
		Alignments align= new Alignments(this.getChr(), this.getStart(), getEnd()+(counter*25));
		align.scores=new ArrayList<Double>();
		align.scores.addAll(this.scores);
		for(int i=0; i<extensionScores.length; i++){align.addScore(extensionScores[i]);}
		return align;
	}
	
	public Alignments extendLHS(int extension){
		return new Alignments(this.getChr(), this.getStart()-extension, getEnd());
		
	}
	
	public String toString(){
		if(this.line==null || this.line.isEmpty()){return getChromosomeString()+"\t"+getStart()+"\t"+getEnd()+"\t"+getName()+"\t"+getScore()+"\t"+getOrientation();}
		else{return this.line;}
	}
	
	public String toFileName(){return getChr()+"_"+getStart()+"-"+getEnd();}
	public String toUCSC(){return this.getChromosomeString()+":"+getStart()+"-"+getEnd();}
	
	public String toStringNum(){
		return getChr().replaceAll("chr", "")+"\t"+getStart()+"\t"+getEnd();
		
	}
	
	public String toStringNumIGV(){
		return getChr().replaceAll("chr", "")+"\t"+getStart();
		
	}
	
	public String toSAM() {
		RefSeqGene gene=new RefSeqGene(this);
		return gene.toSAM();
	}

	/*public int compareTo(Object b){
		Alignments al=(Alignments)b;
		if(!al.chr.equalsIgnoreCase(chr)){return al.chr.compareTo(chr);}
		return -(new Integer(al.getStart()).compareTo(new Integer(start)));
		
	}*/
	
	
	public boolean overlaps2(Alignments next){
		if(next==null){return false;}
		if(!getChr().equals(next.getChr())){return false;}
		if(next.getStart()>=getStart() && next.getStart()<=getEnd()){return true;}
		if(next.getEnd()>=getStart() && next.getEnd()<=getEnd()){return true;}
		//else if(track.getStart()<=start &&track.end<=end){return true;}
		return false;
	}

	public boolean overlapsCollection(Collection<Alignments> c){
		for(Alignments align: c){
			if(overlaps(align)){return true;}
		}
		return false;
	}

	//This is closed-closed
	public boolean overlapsAtAll(Alignments align){
		if(!align.getChr().equalsIgnoreCase(getChr())){return false;}
		if(align.getStart()>=getStart() && align.getStart()<=getEnd()){return true;}
		if(getStart()>=align.getStart() && getStart()<=align.getEnd()){return true;}
		return false;
	}
	
	
	public boolean overlapsAtAll(Set<Alignments> set){
		
		for(Alignments align: set){
			if(overlapsAtAll(align)){return true;}
		}
		return false;
	}
	
	public boolean within(Alignments bigger){
		if(!bigger.getChr().equalsIgnoreCase(getChr())){return false;}
		if(getStart()>=bigger.getStart() && getStart()<=bigger.getEnd()){return true;}
		if(getEnd()>=bigger.getStart() && getEnd()<=bigger.getEnd()){return true;}
		return false;
	}
	
	public int length(){
		return getEnd()-this.getStart();
	}
	
	
	public boolean containedWithin(Alignments align){
		if(!align.getChr().equalsIgnoreCase(getChr())){return false;}
		
		if(align.getStart()>getStart() && align.getEnd()<getEnd()){return true;}
		return false;
	}
	
	public boolean containedWithin(Alignments interval, Alignments gene){
		if(!gene.getChr().equalsIgnoreCase(interval.getChr())){return false;}
		
		if(gene.getStart()>=interval.getStart() && gene.getEnd()<=interval.getEnd()){return true;}
		return false;
	}
	
	public static int distanceBetween(Alignments align1, Alignments align2){
		int start=Math.max(align1.getStart(), align2.getStart());
		int end=Math.min(align1.getEnd(), align2.getEnd());
		return start-end;
	}
	
	
	/*public static Alignments collapseSetOfAlignments(Collection<Alignments> overlappingAlignments){
		Object[] array=overlappingAlignments.toArray();
		Arrays.sort(array);
		Alignments start=(Alignments)array[0];
		Alignments end=(Alignments)array[array.length-1];
		return new Alignments(start.chr, start.start, end.end);
	}*/
	
	//Must be overlapping
	public static Alignments collapseSetOfAlignments(Collection<Alignments> overlappingAlignments){
		Set start=new TreeSet();
		Set end=new TreeSet();
		String chr="";
		for(Alignments align: overlappingAlignments){
			chr=align.getChromosome();
			start.add(align.getStart());
			end.add(align.getEnd());
		}
		Integer startAlign=(Integer)start.toArray()[0];
		Integer endAlign=(Integer)end.toArray()[end.size()-1];
		
		Alignments rtrn=new Alignments(chr, startAlign, endAlign);
		rtrn.scores=new ArrayList();
		for(Alignments align: overlappingAlignments){
			rtrn.scores.addAll(align.scores);
		}
		return rtrn;
	}

	public Alignments merge(Alignments a1, Alignments a2){
		Alignments temp=new Alignments(getChr(), Math.min(a1.getStart(), a2.getStart()), Math.max(a1.getEnd(), a2.getEnd()));
		temp=new Alignments(getChr(), Math.min(temp.getStart(), getStart()), Math.max(temp.getEnd(), getEnd()));
		//System.err.println(a1.toUCSC()+" "+a2.toUCSC()+" "+this.toUCSC()+temp.toUCSC());
		
		return temp;
	}
	
	public Alignments merge(Alignments a1){
		Alignments temp=new Alignments(getChr(), Math.min(a1.getStart(), getStart()), Math.max(a1.getEnd(), getEnd()));
		
		return temp;
	}
	
	
	
	public Alignments excludeRegion(Alignments intron){
		Alignments rtrn=null;
		int compareTo=this.compareTo(intron);
				
		//if compareTo >0 then its on the rhs if <0 then LHS
		
		if(compareTo>0){
			rtrn=new Alignments(getChr(), Math.max(intron.getEnd(), getStart()), getEnd());
		}
		else if(compareTo<0){
			rtrn=new Alignments(getChr(), getStart(), Math.min(intron.getStart(), getEnd()));
		}
		
		return rtrn;
	}
	
	
	public Collection<Alignments> excludeRegionFullyContained(Alignments intron){
		Collection<Alignments> rtrn=new TreeSet<Alignments>();
		int compareTo=this.compareTo(intron);
				
		//if compareTo >0 then its on the rhs if <0 then LHS
		
		if(fullyContained(intron)){
			rtrn.add(new Alignments(getChr(), getStart(), intron.getStart()));
			rtrn.add(new Alignments(getChr(), intron.getEnd(), getEnd()));
		}
		
		else if(compareTo>0){
			rtrn.add(new Alignments(getChr(), Math.max(intron.getEnd(), getStart()), getEnd()));
		}
		else if(compareTo<0){
			rtrn.add(new Alignments(getChr(), getStart(), Math.min(intron.getStart(), getEnd())));
		}
		
		return rtrn;
	}
	
	
	public static void write(String save, Collection<Alignments> alignments) throws IOException {
		write(save, alignments, true);
	}
	
	/**
	 * Write alignments 
	 * @param save
	 * @param alignments
	 * @param sort				sort by genomic position
	 * @throws IOException
	 */
	public static void write(String save, Collection<Alignments> alignments, boolean sort) throws IOException{
		BufferedWriter writer=new BufferedWriter(new FileWriter(save));
		write(writer, alignments, sort);
		writer.close();
	}
	
	
	public int getMidPoint() {
		return (getStart()+getEnd())/2;
	}
	
	/**
	 * Write alignments sorted by score i
	 * @param args
	 */
	public static void write(String save, Collection<Alignments> alignments, final int sortOnScore) throws IOException {
		List<Alignments> keys = new LinkedList<Alignments>(alignments);
		Collections.sort(keys, new Comparator<Alignments>() {
	         @Override
		        public int compare(Alignments o1, Alignments o2) {
		            return o1.getScores().get(sortOnScore).compareTo(o2.getScores().get(sortOnScore));
		         }
		 });
		write(save, alignments, false);
	}
	
	public static void write(BufferedWriter writer, Collection<Alignments> alignments) throws IOException {
		write(writer, alignments, false);
	}
	
	public static void write(BufferedWriter writer, Collection<Alignments> alignments, boolean sort) throws IOException {
		List<Alignments> keys = new LinkedList<Alignments>(alignments);
		if (sort) {
			Collections.sort(keys, new Comparator<Alignments>() {
		         @Override
		         public int compare(Alignments o1, Alignments o2) {
		             return o1.compareTo(o2);
		         }
		     });
		}
		
		for(Alignments align: keys) {
			writer.write(align.toStringWithScores());
			writer.write("\n");
		}
	}
	
	
	public String toStringWithScores() {
		StringBuffer sb = new StringBuffer(toString() + "\t" + countScore);
		for (Double score: getScores()) {
			sb.append("\t" + score);
		}
		return sb.toString();
	}
	
	
	public static void main(String[] args){
		Alignments exon=new Alignments("chr21:17895585-17895672");
		Alignments intron=new Alignments("chr21:17893013-17899048");
		
		System.err.println(intron.fullyContained(exon));
		
	}

}
