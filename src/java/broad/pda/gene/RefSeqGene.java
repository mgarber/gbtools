package broad.pda.gene;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.log4j.Logger;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.AlignmentBlock;

import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.annotation.GFF;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.sequence.Sequence;
import broad.core.util.CollapseByIntersection;
import broad.pda.chromosome.Chromosome;
import broad.pda.datastructures.Alignments;
import broad.pda.rnai.ExtractSequence;
import broad.pda.seq.segmentation.GenomeWithGaps2;
import net.sf.samtools.SAMRecord;

public class RefSeqGene implements Comparable<RefSeqGene>{
	static Logger logger = Logger.getLogger(RefSeqGene.class.getName());
	String chr;
	int start;
	int stop;
	ArrayList<Double> scores;
	Map<String, String> attributes;
	String name;
	int[] exonStarts;
	int[] exonEnds;
	String orientation;
	int blockStart; // beginning of CDS
	int blockEnd; // end of CDS
	double[] exonScores;
	String sequence;
	private String samRecord;
	private double countScore=0; //Moran -added this as init value
	private Collection<Alignments> sortedExons;
	double bedScore; //the transcript score as it appears in a bed file
	String [] extraFields;
	boolean isSorted=false;
	
	
	public static final Pattern START_CODON_PATRN= Pattern.compile("ATG",  Pattern.CASE_INSENSITIVE);
	public static final String [] STOP_CODONS =  {"TAG", "TAA", "TGA"};
	public static final String UNORIENTED_STRING = "*";
	
	public RefSeqGene(String pslString){
		this(pslString, true);
	}
	
	public RefSeqGene(Collection<? extends LightweightGenomicAnnotation> exons, String name, String orientation){
		//System.err.println(exons.size());
		
		int[] exonStart=new int[exons.size()];
		int[] exonEnd=new int[exons.size()];

		int i=0;
		this.start=Integer.MAX_VALUE;
		this.stop=0;
		String chr="";
		this.orientation = orientation;
		
		this.blockStart=stop;
		this.blockEnd=stop;
		
		for(LightweightGenomicAnnotation exon: exons){
			chr=exon.getChromosome();
			start=Math.min(start, exon.getStart());
			stop=Math.max(stop, exon.getEnd());
			exonStart[i]=exon.getStart();
			exonEnd[i]=exon.getEnd();
			orientation = exon.getOrientation() == null ? "*" : exon.getOrientation();
			i++;
		}
		this.chr=chr;
		this.name=name == null ? getAlignment().toUCSC() : name;
		
		this.exonStarts=exonStart;
		this.exonEnds=exonEnd;
		this.exonScores=new double[this.exonStarts.length];
	}
	
	public RefSeqGene(String rawData, boolean isPSLFormat) {
		if(isPSLFormat) {
			String[] tokens=rawData.split("\t");
			this.chr=tokens[13];
			this.start=new Integer(tokens[15]);
			this.stop=new Integer(tokens[16]);
			this.name=tokens[9];
			this.orientation=tokens[8];
			if (this.orientation.equals("."))
				this.orientation="*";
			
			int numExons=new Integer(tokens[17]);
			int[] exonStart=new int[numExons];
			int[] exonEnd=new int[numExons];
			
			String[] exonStarts=tokens[20].replaceAll("\"", "").trim().split(",");
			String[] exonSizes=tokens[18].replaceAll("\"", "").trim().split(",");
			
			for(int i=0; i<exonStarts.length; i++){
				exonStart[i]=new Integer(exonStarts[i].toString());
				exonEnd[i]=exonStart[i]+new Integer(exonSizes[i].toString());
			}
			
			this.exonStarts=exonStart;
			this.exonEnds=exonEnd;
			this.exonScores=new double[this.exonStarts.length];
			
			
			Sequence seq=new Sequence(name);
			seq.setSequenceBases(tokens[21].replaceAll(",", "").trim().toUpperCase());
			//System.err.println(tokens[21]);
			this.addSequence(seq);
			
		} else { //Assume regular BED format
			//System.err.println(rawData);
	       	String[] tokens=rawData.split("\t");
			chr=(tokens[0]);
			start=new Integer(tokens[1]);
			stop=new Integer(tokens[2]);
			orientation = "*";
			if(tokens.length > 3) {
				name=tokens[3];
				if(tokens.length > 4) {
					bedScore = new Double(tokens[4]);
					if(tokens.length > 5){
						orientation=tokens[5];
						if (this.orientation.equals("."))
							this.orientation="*";
						
						if(tokens.length > 6) {
							blockStart=new Integer(tokens[6]);
							blockEnd=new Integer(tokens[7]);
							
							String[] blockSizes=tokens[10].split(",");
							String[] blockStarts=tokens[11].split(",");
							List<Integer>[] exonStartEnd=setBlockStartsAndEnds(blockStarts, blockSizes, new Integer(tokens[9]), start);
							exonStarts = new int[blockSizes.length];
							exonEnds   = new int[blockSizes.length];
							for(int i = 0; i < blockSizes.length; i++ ) {
								exonStarts[i] = exonStartEnd[0].get(i);
								exonEnds[i] = exonStartEnd[1].get(i);
							}
							String [] extraColumns = null;
							if(tokens.length > 12) {
								extraColumns = new String[tokens.length - 12];
								for(int j = 12; j < tokens.length; j++) {
									extraColumns[j - 12] = tokens[j];
								}
							}
							extraFields = extraColumns;
						}
					}
				}
			}
		}
	}
	
	public RefSeqGene (String chrString, int start, int end){
		this(new Alignments(chrString, start,end));
	}
	
	public RefSeqGene(String chr, int start, int end, String name, String orientation, int[] exonsStart, int[] exonsEnd){
		this.chr=chr;
		this.start=start;
		this.stop=end;
		this.name=name;
		this.orientation=orientation;
		if (this.orientation.equals("."))
			this.orientation="*";
		
		
		this.exonStarts=exonsStart;
		this.exonEnds=exonsEnd;
		this.exonScores=new double[this.exonStarts.length];
		//this.blockStart = exonsStart[0];
		//this.blockEnd = exonsEnd[exonsEnd.length - 1];
	}
	
	public RefSeqGene(String chr, int start, int end, String name, String orientation, int[] exonsStart, int[] exonsEnd, int blockStart, int blockEnd){
		this(chr, start, end, name, orientation, exonsStart, exonsEnd);
		this.blockStart = blockStart;
		this.blockEnd = blockEnd;
	}
	
	
	
	public RefSeqGene(String chr2, int start2, int end, String name2,
			double bedScore2, String orientation2, String sequence2,String[] extraFields2) {
		this(chr2,start2,end,name2,bedScore2,orientation2,null,null,null,sequence2,extraFields2);
	}
	
	public RefSeqGene(String chr2, int start2, int end, String name2,
			double bedScore2, String orientation2, String sequence2) {
		this(chr2,start2,end,name2,bedScore2,orientation2,null,null,null,sequence2);
	}
	
	public RefSeqGene(String chr, int start, int end, String name,double mybedScore, String orientation, int[] exonsStart, int[] exonsEnd){
		this( chr,  start,  end,  name,  orientation,  exonsStart, exonsEnd);
		this.bedScore=mybedScore;	
	}
	
	public RefSeqGene(String chr, int start, int end, String name, String orientation, List<Integer> exonsStart, List<Integer> exonsEnd){
		this(chr, start, end, name, orientation, exonsStart, exonsEnd, null);	
	}
	
	public RefSeqGene(String chr, int start, int end, String name, String orientation, Collection<Alignments> exons){
		this.chr=chr;
		this.start=start;
		this.stop=end;
		this.name=name;
		this.orientation=orientation;
		if (this.orientation.equals("."))
			this.orientation="*";
		
		int[] exonStart=new int[exons.size()];
		int[] exonEnd=new int[exons.size()];
	
		int i=0;
		for(Alignments exon: exons){
			exonStart[i]=exon.getStart();
			exonEnd[i]=exon.getEnd();
			i++;
		}
		
		this.exonStarts=exonStart;
		this.exonEnds=exonEnd;
		this.exonScores=new double[this.exonStarts.length];
	}
	
	public RefSeqGene(String chr, int start, int end, String name, String orientation, List<Integer> exonsStart, List<Integer> exonsEnd, int blockStart, int blockEnd){
		this(chr,start,end,name, orientation, exonsStart, exonsEnd);
		this.blockStart=blockStart;
		this.blockEnd=blockEnd;
	}
	
	//Bug fix Dec 28th 2010 , was not intializing extraFields
	public RefSeqGene(String chr, int start, int end, String name, String orientation, List<Integer> exonsStart, List<Integer> exonsEnd, String [] extraData){
		this.chr=chr;
		this.start=start;
		this.stop=end;
		this.name=name;
		this.orientation=orientation;
		
		if (this.orientation.equals("."))
			this.orientation="*";
		
		int[] exonStart=new int[exonsStart.size()];
		int[] exonEnd=new int[exonsStart.size()];
	
		for(int i=0; i<exonStart.length; i++){
			exonStart[i]=exonsStart.get(i);
			exonEnd[i]=exonsEnd.get(i);
		}
		
		this.exonStarts=exonStart;
		this.exonEnds=exonEnd;
		this.exonScores=new double[this.exonStarts.length];
		this.blockStart = exonsStart.get(0);
		this.blockEnd = exonsEnd.get(exonsEnd.size() - 1);
		this.bedScore=0; //initialize with zero
		
		if(extraData != null) 
			setExtraFields(extraData);
	}
	
	public RefSeqGene(String chr2, int start2, int end, String name2,
			double bedScore2, String orientation2, int[] exonStarts2,
			int[] exonEnds2, double[] exonScores2,String sequence2) {
			this.chr=chr2;
			this.start=start2;
			this.stop=end;
			this.name=name2;
			this.bedScore=bedScore2;
			this.orientation=orientation2;
			if (this.orientation.equals("."))
				this.orientation="*";
			this.exonStarts=exonStarts2;
			this.exonEnds=exonEnds2;
			this.exonScores=exonScores2;
			this.sequence=sequence2;
	}
	
	public RefSeqGene(String chr2, int start2, int end, String name2,
			double bedScore2, String orientation2, int[] exonStarts2,
			int[] exonEnds2, double[] exonScores2,String sequence2,String[] extraFields2) {
			this.chr=chr2;
			this.start=start2;
			this.stop=end;
			this.name=name2;
			this.bedScore=bedScore2;
			this.orientation=orientation2;
			if (this.orientation.equals("."))
				this.orientation="*";
			this.exonStarts=exonStarts2;
			this.exonEnds=exonEnds2;
			this.exonScores=exonScores2;
			this.sequence=sequence2;
			this.extraFields=extraFields2;
	}
	
	public RefSeqGene(String chr, int start, int end, String name,
			double score, String strand, List<Integer> exonsStart,
			List<Integer> exonsEnd, int blockStart2, int blockEnd2) {
		this( chr,  start,  end,  name, strand,  exonsStart,  exonsEnd,  blockStart2, blockEnd2);
		this.bedScore=score;
	}
	
	public RefSeqGene(String chr, int start, int end, String name, String orientation, List<Integer> exonsStart, List<Integer> exonsEnd, int blockStart, int blockEnd, String [] extraData){
		this(chr,start,end,name, orientation, exonsStart, exonsEnd);
		this.blockStart=blockStart;
		this.blockEnd=blockEnd;
		if(extraData != null) {
			this.extraFields = extraData;
		}
	}
	
	public RefSeqGene(String chr, int start, int end, String name,
			double score, String strand, List<Integer> exonsStart,
			List<Integer> exonsEnd, int blockStart2, int blockEnd2, String [] extraColumns) {
		this( chr,  start,  end,  name, strand,  exonsStart,  exonsEnd,  blockStart2, blockEnd2, extraColumns);
		this.bedScore=score;
	}
	
	public RefSeqGene(Alignment align){
		AlignmentBlock[] blocks=align.getAlignmentBlocks();
		
		Collection<Alignments> exons=new TreeSet();
		
		if(blocks==null){
			Alignments exon=new Alignments(align.getChromosome(), align.getAlignmentStart(), align.getAlignmentEnd());
			exons.add(exon);
		}
		else{
			for(int i=0; i<blocks.length; i++){Alignments exon=new Alignments(align.getChromosome(), blocks[i].getStart(), blocks[i].getEnd()); exons.add(exon);}
		}
		int[] exonStart=new int[exons.size()];
		int[] exonEnd=new int[exons.size()];
	
		int i=0;
		this.start=Integer.MAX_VALUE;
		this.stop=0;
		String chr="";
		
		for(Alignments exon: exons){
			chr=exon.getChr();
			start=Math.min(start, exon.getStart());
			stop=Math.max(stop, exon.getEnd());
			exonStart[i]=exon.getStart();
			exonEnd[i]=exon.getEnd();
			i++;
		}
		
		this.orientation="*";
		this.chr=chr;
		this.name=this.getAlignment().toUCSC();
		
		
		this.exonStarts=exonStart;
		this.exonEnds=exonEnd;
		this.exonScores=new double[this.exonStarts.length];
		
	}
	
	public RefSeqGene(RefSeqGene gene) {
		this(gene.getChr(), gene.getStart(), gene.getEnd(), gene.getName(), gene.getOrientation(), gene.getExonSet());
		this.blockStart=gene.blockStart;
		this.blockEnd=gene.blockEnd;
		
		if (gene.extraFields != null)
			setExtraFields(gene.getExtraFields());
		if (gene.scores !=null)
			setScores(gene.getScores());
		
		if (gene.attributes !=null)
			setAttributes(gene.getAttributes());
		this.isSorted=gene.isSorted;
		if (gene.samRecord!=null)
			this.samRecord=gene.getSAMString();
		if (gene.countScore!=0)  
			this.countScore=gene.getCountScore();
		if (gene.sortedExons!=null)
			setSortedExons(gene.getSortedAndUniqueExons());
		
	}
	
	public RefSeqGene(LightweightGenomicAnnotation align){
		this.chr=align.getChromosome();
		this.start=align.getStart();
		this.stop=align.getEnd();
		this.name=align.getName();
		this.orientation=align.getOrientation();
		if (this.orientation.equals("."))
			this.orientation="*";
		int[] exonStart=new int[1];
		int[] exonEnd=new int[1];
	
		exonStart[0]=align.getStart();
		exonEnd[0]=align.getEnd();
			
		//this.sequence=setN(align.length());
		
		
		this.exonStarts=exonStart;
		this.exonEnds=exonEnd;
		this.exonScores=new double[this.exonStarts.length];
		
	}
	
	public RefSeqGene(Collection<? extends LightweightGenomicAnnotation> exons){
		this(exons, null);
	}
	
	public RefSeqGene(Collection<? extends LightweightGenomicAnnotation> exons, String name){
		//System.err.println(exons.size());
		
		int[] exonStart=new int[exons.size()];
		int[] exonEnd=new int[exons.size()];
	
		int i=0;
		this.start=Integer.MAX_VALUE;
		this.stop=0;
		String chr="";
		orientation = "*";
		
		
		for(LightweightGenomicAnnotation exon: exons){
			chr=exon.getChromosome();
			start=Math.min(start, exon.getStart());
			blockStart = start;
			stop=Math.max(stop, exon.getEnd());
			blockEnd = stop;
			exonStart[i]=exon.getStart();
			exonEnd[i]=exon.getEnd();
			orientation = exon.getOrientation() == null ? "*" : exon.getOrientation();
			i++;
		}
		this.chr=chr;
		this.name=name == null ? getAlignment().toUCSC() : name;
		
		this.exonStarts=exonStart;
		this.exonEnds=exonEnd;
		this.exonScores=new double[this.exonStarts.length];
	}
	
	public void setSequence(String seq){this.sequence=seq;}
	
	public void setSequenceFromChromosome(Sequence chrSequence) {
		Sequence  geneSequence = new Sequence(getName(),getGappedSize());
		Set<Alignments> exons = getExonSet();
		for(Alignments exon : exons) {
			Sequence exonSeq = chrSequence.getSubSequence(getName(), exon.getStart(), exon.getEnd());
			geneSequence.appendToSequence(exonSeq.getSequenceBases());
		}
		
		if("-".equals(orientation)) {
			geneSequence.reverse();
		}
		this.sequence = geneSequence.getSequenceBases();
	}
	
	public void setOrientation(String orientation){this.orientation=orientation;}
	
	public void setChromosome(String chr) {
		this.chr = chr;
	}
	
	public void setName(String name){this.name=name;}
	
	private void setStart(int i) { this.start=i; }
	
	private void setEnd(int i) { this.stop=i; }
	
	public void setUnorientedStart(int newStart) {
		if(!isSorted) {
			
		}
		int oldStart = getStart();
		setStart(newStart);
		if(this.blockStart == oldStart) {
			this.blockStart = newStart;
		}
		if(exonStarts[0] == oldStart) {
			this.exonStarts[0] = newStart;
		}
	}
	
	public void setUnorientedEnd(int newEnd) {
		int oldEnd = getEnd();
		setEnd(newEnd);
		if(blockEnd == oldEnd) {
			this.blockEnd = newEnd;
		} if (exonEnds[this.exonEnds.length - 1]==oldEnd) {
			this.exonEnds[this.exonEnds.length - 1] = newEnd;
		}
	}
	
	public void setSAMString(String sam){this.samRecord=sam;}
	
	public void setExtraFields(double[] d) {
		extraFields = new String[d.length];
		for (int i = 0; i < d.length; i++) {
			String s = (new Double(d[i])).toString();
			extraFields[i] = s;
		}
	}
	
	public void setExtraFields(double val, int i) {
		if (i>0 && i<=this.extraFields.length)
			this.extraFields[i-1]=String.valueOf(val);
		
	}
	
	private void setExtraFields(String[] extraData){
		String[] extraDataCpy= new String[extraData.length];
		for (int i=0; i<extraData.length; i++)
			extraDataCpy[i]=extraData[i];
		this.extraFields = extraDataCpy;
	}
	
	private void setScores(ArrayList<Double> scores2) {
		ArrayList<Double> newscores = new ArrayList<Double>();
		for (int i=0; i< scores2.size(); i++)
			newscores.add(new Double(scores2.get(i)));
		this.scores=newscores;
		
	}
	
	private static List<Integer> [] setBlockStartsAndEnds(String[] blockStarts, String[] blockSizes, int size, int start){
		List<Integer> starts=new ArrayList<Integer> ();
		List<Integer>  end=new ArrayList<Integer> ();
		for(int i=0; i<size; i++){
			starts.add(start+new Integer(blockStarts[i].replaceAll("\"", "").trim()));
			end.add((Integer)starts.get(i)+new Integer(blockSizes[i].replaceAll("\"", "").trim()));
		}
		List [] rtrn={starts, end};
		return rtrn;
	}
	
	private void setAttributes(Map<String, String> attrs) {
		Map<String, String> m= new HashMap<String, String>();
		for (String key: attrs.keySet())
			m.put(key, attrs.get(key));
		this.attributes=m;
	}
	
	public void set5PrimeEnd(int updated5Prime) {
		if("-".equals(orientation)) {
			for(int i = 0 ; i < exonEnds.length; i++){
				if(exonEnds[i] ==  stop) {
					exonEnds[i] = updated5Prime;
				}
			}
			this.stop = updated5Prime;
		} else {
			for(int i = 0 ; i < exonEnds.length; i++){
				if(exonStarts[i] ==  start) {
					exonStarts[i] = updated5Prime;
				}
			}
			this.start = updated5Prime;
			
		}
	}
	
	public void set3PrimeEnd(int updated3Prime) {
		if("-".equals(orientation)) {
			for(int i = 0 ; i < exonEnds.length; i++){
				if(exonStarts[i] ==  start) {
					exonStarts[i] = updated3Prime;
				}
			}
			this.start = updated3Prime;
		} else {
			for(int i = 0 ; i < exonEnds.length; i++){
				if(exonEnds[i] ==  stop) {
					exonEnds[i] = updated3Prime;
				}
			}
			this.stop = updated3Prime;
			
		}
	}
	
	public void set5PrimeEndForGeneOnly(int updated5Prime) {
		if("-".equals(orientation)) {
			for(int i = 0 ; i < exonEnds.length; i++){
				if(exonEnds[i] ==  stop) {
					exonEnds[i] = updated5Prime;
				}
			}
			this.stop = updated5Prime;
		} else {
			for(int i = 0 ; i < exonEnds.length; i++){
				if(exonStarts[i] ==  start) {
					exonStarts[i] = updated5Prime;
				}
			}
			this.start = updated5Prime;
			
		}
	}
	
	public void set3PrimeEndForGeneOnly(int updated3Prime) {
		if("-".equals(orientation)) {
			for(int i = 0 ; i < exonEnds.length; i++){
				if(exonStarts[i] ==  start) {
					exonStarts[i] = updated3Prime;
				}
			}
			this.start = updated3Prime;
		} else {
			for(int i = 0 ; i < exonEnds.length; i++){
				if(exonEnds[i] ==  stop) {
					exonEnds[i] = updated3Prime;
				}
			}
			this.stop = updated3Prime;
			
		}
	}
	
	/**
	 * Set coding region coordinates (includes start codon but not stop codon)
	 * @param start start position
	 * @param end end position
	 */
	public void setCDSRegion(int start, int end){
		if(start > end) {
			throw new IllegalArgumentException("Start position " + start + " is greater than end position " + end);
		}
		if((start < this.start && start < this.stop) || (start > this.start && start > this.stop)) {
			throw new IllegalArgumentException("Start position " + start + " is out of range of gene start and end " + this.start + "-" + this.stop);
		}
		if((end < this.start && end < this.stop) || (end > this.start && end > this.stop)) {
			throw new IllegalArgumentException("End position " + end + " is out of range of gene start and end " + this.start + "-" + this.stop);
		}
		this.blockStart = start;
		this.blockEnd = end;
	}
	
	protected void setExons(Alignments[] exons,	double[] giveExonsScores) {
		
		this.exonStarts=new int[exons.length] ;
		this.exonEnds=new int[exons.length];
		this.exonScores=new double[exons.length];
		for (int i=0; i<exons.length; i++) {
			this.exonStarts[i]=exons[i].getStart();
			this.exonEnds[i]=exons[i].getEnd();
			this.exonScores[i]= giveExonsScores != null ? giveExonsScores [i] : 0;
		}
		
	}
	
	public void setExonsScores(Alignments exon, double score){
		int index=-1;
		for(int i=0; i<this.getNumExons(); i++){
			if(exon.getStart()==exonStarts[i]){index=i;}
		}
		if (index>-1 )
			this.exonScores[index]=score;
	}
	
	private void setSortedExons(Collection<Alignments> sortedAndUniqueExons) {
		Collection<Alignments> set=new LinkedList<Alignments>();
		for(Alignments a : sortedAndUniqueExons )
			set.add(new Alignments (a));
		this.sortedExons=set;
	}
	
	public void setBedScore (double d){bedScore=d;}
	
	public void setCountScore(double score){this.countScore=score;}
	
	private String setN(int len){
		String rtrn="";
		for(int i=0; i<len; i++){rtrn+="N";}
		return rtrn;
	}
	
	public String getName(){return this.name;}
	
	public int getStart(){return this.start;}
	
	public int getEnd(){return this.stop;}
	
	public String getChr(){return this.chr;}
	
	public int getSize(){
		int rtrn=0;
		Collection<Alignments> exons=this.getExonSet();
		
		for(Alignments exon: exons){rtrn+=exon.getSize();}
		
		return rtrn;
	}
	
	public int getGenomicLength() {return stop - start;}
	
	public String getOrientation(){return this.orientation;}
	
	public boolean isNegativeStrand(){return this.orientation.equals("-");}
	
	public String getSequence(){return this.sequence;}
	
	public Sequence getSequence(Sequence chrSequence) {
		Sequence seq=new Sequence(getName());
		Collection<Alignments> exons=this.getSortedAndUniqueExons();
		
		for(Alignments exon: exons){
			Sequence sequence=chrSequence.getSubSequence("", exon.getStart(), exon.getEnd());
			seq.appendToSequence(sequence.getSequenceBases());
		}
		
		return seq;
	}
	
	public String getSequence(Chromosome chr2, boolean repeatMask, boolean stranded) throws Exception {
		return ExtractSequence.getSequenceForGene(this, chr2, repeatMask, new TreeMap(), stranded);
	}
	
	public Sequence getSequenceObject(){
		Sequence seq=new Sequence(this.name);
		seq.setSequenceBases(this.sequence);
		return seq;
	}
	
	public String getString(){
		return getChr()+":"+this.start+"-"+this.stop;
	}
	
	public String getSAMString(){return this.samRecord;}
	
	public int getOrientedStart() {
		if (this.orientation.equals("-"))
			return stop;
		else
			return start;	
	}
	
	public int getOrientedEnd() {
		if (this.orientation.equals("-"))
			return start;
		else
			return stop;	
	}
	
	/**
	 * 
	 * @param name
	 * @return attribute value or 0 if the attribute is not defined.
	 */
	public String getAttribute(String name) {
		return attributes == null || !attributes.containsKey(name) ? null : attributes.get(name);
	}
	
	public Map getAttributes() { return attributes;} //TODO: Return a copy, attributes should not be modified outside of accessor methods.
	
	/**
	public Integer[] getExtraFields() {
		// TODO Auto-generated method stub
		return null;
	}*/
		public String[] getExtraFields() {
			// TODO Auto-generated method stub
			if (this.extraFields!=null)
				return extraFields;
			else 
				return null;
		}
		
	public String getExtraFields(int i) {
		
		if (this.extraFields!=null &&  this.extraFields.length>=(i+1))
			return this.extraFields[i];
		else
		{	
			//System.err.println("no extra fields");
			return "";
		}
		
	}
	
	public double getScore(){
		double score=0;
		if (this.exonScores !=null)
			for(int i=0; i<this.exonScores.length; i++){score+=exonScores[i];}
		return score;
	}
	
	private ArrayList<Double> getScores() {
		return this.scores;
	}
	
	public double getNormalizedScore(){
		double score=this.getScore();
		Set<Alignments> exons=this.getExonSet();
		
		if(exons.isEmpty()) return 0;
		
		return score/getTranscriptLength();
	}
	
	public double getRPKM() {
		if (extraFields != null)
			return Double.parseDouble(extraFields[4]);
		else
			return -1;
	}
	
	public double getPValue() {
		if (extraFields != null)
			return Double.parseDouble(extraFields[0]);
		else
			return -1;
	}
	
	public double getCountScore(){return this.countScore;}
	
	public double[] getExonsScores(){
		return this.exonScores;
	}
	
	/*Returns -1 if score for this exon does not exist*/
	public double getExonsScores(Alignments exon){
		int index=-1;
		if (this.exonStarts==null) return -1;
		for(int i=0; i<this.exonStarts.length; i++){
			if(exon.getStart()==exonStarts[i]){index=i;}
		}
		if (index==-1) return -1;
		return this.exonScores[index];
	}
	
	public List<LightweightGenomicAnnotation> getScoredExons() {
		ArrayList<LightweightGenomicAnnotation> rtrn=new ArrayList<LightweightGenomicAnnotation>();
		for(int i=0; i<exonStarts.length; i++){
			rtrn.add(new BasicLightweightAnnotation(this.chr, this.exonStarts[i], this.exonEnds[i],this.getOrientation(),this.getBedScore()));
		}
		return rtrn;
		
	}
	
	public Alignments getSingleOverlappingExon(Alignments a) {
		for (Alignments exon : getExonSet()) {
			//if (exon.overlaps2(a))
			if(exon.overlaps(a)) 
				return exon;
		}
		
		return null;
	}
	
	public int getTranscriptLength(){
		
		int length=0;
		Set<Alignments> exons=this.getExonSet();
		
		if(exons.size() > 0) {
			for(Alignments exon: exons){
				length = length + exon.length();
			}
		} else {
			length = getSize();
		}
		return length;
	}
	
	/**
	 * Get the set of exons comprising the 5' UTR
	 * @return the set of exons comprising the 5' UTR as a collection of Alignments objects
	 */
	public Collection<Alignments> get5Utr() {
		return this.get5UtrIncludingIntrons().getIntersections(this.getSortedAndUniqueExons());
	}
	
	public RefSeqGene get3UTRGene() {
		if(!hasCDS()){return new RefSeqGene(this);}
		Alignments UTRRegion=get3UTR();
		RefSeqGene rtrn=this.trimAbsolute(UTRRegion.getStart(), UTRRegion.getEnd());		
		return rtrn;
	}
	
	
	
	public RefSeqGene get5UTRGene() {
		if(!hasCDS()){return new RefSeqGene(this);}
		Alignments UTRRegion=get5UTR();
		RefSeqGene rtrn=this.trimAbsolute(UTRRegion.getStart(), UTRRegion.getEnd());		
		return rtrn;
	}
	
	public Alignments get5UTR(){
		if(orientation.equalsIgnoreCase("+")){
			Alignments rtrn=new Alignments(chr, start, this.getCDSRegion().getStart());
			return rtrn;
		}
		Alignments rtrn=new Alignments(chr, this.getCDSRegion().getEnd(), this.stop);
		return rtrn;
	}
	
	
	public Alignments get3UTR(){
		if(orientation.equalsIgnoreCase("-")){
			Alignments rtrn=new Alignments(chr, start, this.getCDSRegion().getStart());
			return rtrn;
		}
		Alignments rtrn=new Alignments(chr, this.getCDSRegion().getEnd(), this.stop);
		return rtrn;
	}
	
	public boolean hasCDS(){
		if(this.blockStart==this.blockEnd){return false;}
		return true;
	}
	
	/**
	 * Get the set of exons comprising the 3' UTR
	 * @return the set of exons comprising the 3' UTR as a collection of Alignments objects
	 */
	public Collection<Alignments> get3Utr() {
		return this.get3UtrIncludingIntrons().getIntersections(this.getSortedAndUniqueExons());
	}

	
	public Collection<RefSeqGene> getWindows(int windowSize) {
		Collection<RefSeqGene> subGenes=new TreeSet<RefSeqGene>();
		for(int i=0; i<getTranscriptLength(); i++){
			RefSeqGene subGene=trim(i, i+windowSize);
			if(subGene!=null){
				subGenes.add(subGene);
			}
		}
		return subGenes;
	}
	
	public Collection<RefSeqGene> getWindows(int windowSize, int stepSize, int start) {
		Collection<RefSeqGene> subGenes=new TreeSet<RefSeqGene>();
		for(int i=start; i<getTranscriptLength(); i=i+stepSize){
			RefSeqGene subGene=trim(i, i+windowSize);
			if(subGene!=null){
				subGenes.add(subGene);
			}
		}
		return subGenes;
	}
	
	public RefSeqGene getStartCodon() {
		if(this.getCDS() == null) return null;
		RefSeqGene gene=this.getCDS();
		//System.err.println(gene);
		if(this.orientation == null) return null;
		if(gene.getOrientation().equalsIgnoreCase("+")){
			return gene.trim(0, 3);
		}
		else if(gene.getOrientation().equalsIgnoreCase("-")){
			return gene.trim(gene.getTranscriptLength()-3, gene.getTranscriptLength());
		}
		else{
			System.err.println("UNKNOWN Orientation... assuming +");
			return gene.trimAbsolute(gene.getStart(), gene.getStart()+3);
		}
	}
	

	
	/**
	 * Get the region of the gene span that lies 5' of start codon, including introns and the start codon itself
	 * @return the region of the gene span that lies 5' of start codon as an Alignments object
	 */
	public Alignments get5UtrIncludingIntrons(){
		if(orientation.equalsIgnoreCase("+")){
			Alignments rtrn=new Alignments(chr, start, this.getCDSRegion().getStart() - 1);
			return rtrn;
		}
		Alignments rtrn=new Alignments(chr, this.getCDSRegion().getEnd() + 1, this.stop);
		return rtrn;
	}
	
	/**
	 * Return the chromosome number without the chr
	 * @return the chromosome number without the chr
	 */
	public String getChrNum() {
		String chr=getChr();
		return chr.replaceAll("chr", "");
	}
	
	/**
	 * Get the region of the gene span that lies 3' of the stop codon, including introns but not the stop codon itself
	 * @return the region of the gene span that lies 3' of the stop codon as an Alignments object
	 */
	public Alignments get3UtrIncludingIntrons(){
		if(orientation.equalsIgnoreCase("-")){
			Alignments rtrn=new Alignments(chr, start, this.getCDSRegion().getStart() - 1);
			return rtrn;
		}
		Alignments rtrn=new Alignments(chr, this.getCDSRegion().getEnd() + 1, this.stop);
		return rtrn;
	}
	
	public int getNumExons(){
	if (this.exonEnds != null)
		return this.exonEnds.length;
	else
		return 0;
	}
	
	public Alignments[] getExons(){
		
		if(this.exonStarts==null) 
			return null;
		Alignments[] rtrn=new Alignments[this.exonStarts.length];
		for(int i=0; i<rtrn.length; i++){
			rtrn[i]=new Alignments(this.chr, this.exonStarts[i], this.exonEnds[i]);
		}
		return rtrn;
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
	
	public Collection<Alignments> getSortedAndUniqueExons(){
		if(isSorted){return this.sortedExons;}
		Collection<Alignments> rtrn=new ArrayList();
		
		Set<Alignments> exons=this.getExonSet();
		//System.err.println("\t\tGetting sorted and unique exons. exons start " + Arrays.toString(exonStarts) + " - " + Arrays.toString(exonEnds));
		rtrn=CollapseByIntersection.collapseByIntersection(exons, false);
	
		this.sortedExons=rtrn;
		this.isSorted=true;
		return rtrn;
	}
	
	public int[] getExonSizes(){
		if(this.exonEnds==null) 
			return null;
		int[] rtrn=new int[this.exonEnds.length];
		for(int i=0; i<rtrn.length; i++){
			rtrn[i]=exonEnds[i]-exonStarts[i];
		}
		return rtrn;
	}
	
	public Alignments getFirstExon(){
		if(!isSorted){sortedExons=this.getSortedAndUniqueExons();}
		return sortedExons.iterator().next();
	}
	
	public Alignments getLastExon(){
		if(!isSorted){sortedExons=this.getSortedAndUniqueExons();}
		return (Alignments)sortedExons.toArray()[sortedExons.size()-1];
	}
	
	public IntervalTree<Alignments> getExonTree() {
		IntervalTree<Alignments> rtrn=new IntervalTree();
		
		for(Alignments exon: this.getExonSet()){rtrn.put(exon.getStart(), exon.getEnd(), exon);}
		
		return rtrn;
	}
	
	public Alignments get5PrimeExon() {
		//if + return the first exon
		if(orientation.equalsIgnoreCase("+")){
			Alignments rtrn=getExons()[0];
			rtrn.setOrientation("-");
			return rtrn;
		}
		//if minus return the last exon
		if(orientation.equalsIgnoreCase("-")){
			Alignments rtrn=getExons()[this.getNumExons()-1];
			rtrn.setOrientation("+");
			return rtrn;
		}
		//else return null;
		return null;
	}
	
	public Alignments get3PrimeExon(int size) {
		//if + return the last exon
		if(orientation.equalsIgnoreCase("+")){
			Alignments rtrn=getExons()[this.getNumExons()-1];
			if(rtrn.getSize()>size){
				rtrn=new Alignments(rtrn.getChr(), rtrn.getEnd()-size, rtrn.getEnd());
			}
			rtrn.setOrientation("+");
			return rtrn;
		}
		//if minus return the first exon
		if(orientation.equalsIgnoreCase("-")){
			Alignments rtrn=getExons()[0];
			if(rtrn.getSize()>size){
				rtrn=new Alignments(rtrn.getChr(), rtrn.getStart(), rtrn.getStart()+size);
			}
			rtrn.setOrientation("-");
			return rtrn;
		}
		//else return null;
		return null;
	}
	
	public Alignments get3PrimeExon() {
		
		//if + return the last exon
		if(orientation.equalsIgnoreCase("+")){
			Alignments rtrn=getExons()[this.getNumExons()-1];
			rtrn.setOrientation("+");
			return rtrn;
		}
		//if minus return the first exon
		if(orientation.equalsIgnoreCase("-")){
			Alignments rtrn=getExons()[0];
			rtrn.setOrientation("-");
			return rtrn;
		}
		//else return null;
		return null;
	}
	
	public Alignments getOrientedLastExon() {
		if (!isSorted){sortedExons=this.getSortedAndUniqueExons();}
		if (this.orientation.equals("-"))
			return (Alignments)sortedExons.toArray()[0];
		else
			return (Alignments)sortedExons.toArray()[sortedExons.size() - 1];
	}
	
	public RefSeqGene getIntrons(){
		return new RefSeqGene(this.getIntronSet());
	}
	
	public Alignments getLastIntron() {
		RefSeqGene introns = getIntrons();
		if(introns.getNumExons() == 0) throw new IllegalArgumentException("Gene has no introns");
		return introns.getOrientedLastExon();
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
	
	public Collection<Alignments> getIntronSet(){
		Collection<Alignments> rtrn=new TreeSet<Alignments>(); 
		
		Object[] sortedUniqueExons=this.getSortedAndUniqueExons().toArray();
		
		for(int i=0; i<sortedUniqueExons.length-1; i++){
			Alignments current=(Alignments)sortedUniqueExons[i];
			Alignments next=(Alignments) sortedUniqueExons[i+1];
			//Alignments align=new Alignments(this.chr, current.getEnd()+1, next.getStart()-1); // If we assume Exons are in [start,end) format there should be no need add one to the end.  
			Alignments align=new Alignments(this.chr, current.getEnd(), next.getStart());   
			//if(align.getSize()<0){System.err.println(align.toUCSC()+" "+current.toUCSC()+" "+next.toUCSC()+" "+toBED());}
			align.setOrientation(this.orientation);
			rtrn.add(align);
		}
		return rtrn;
	}
	
	public RefSeqGene getIntronTranscript() {
	
		if (this.getNumExons()<=1)
			return null;
		Alignments[] iarr=this.getIntronsBlocks();
		int g_st=iarr[0].getStart();
		int g_end=iarr[iarr.length -1].getEnd();
		RefSeqGene g =new RefSeqGene(this.getChr(),g_st,g_end,this.getName(),this.getOrientation(),this.getIntronSet());
		return g;
	}
	
	/*public RefSeqGene getCDS(){
		Alignments cds=getCDSRegion();
		Collection<Alignments> exons=new TreeSet<Alignments>();
		
		
		for(Alignments exon: this.getExonSet()){
			if(exon.overlapsAtAll(cds)){exons.add(exon);}	
		}
		
		Alignments firstExon=exons.iterator().next();
		Alignments lastExon=(Alignments)exons.toArray()[exons.size()-1];
		
		exons.remove(firstExon);
		exons.remove(lastExon);
		
		firstExon=trimFirst(firstExon, cds);
		lastExon=trimLast(lastExon, cds);
		
		exons.add(firstExon);
		exons.add(lastExon);
		
		return new RefSeqGene(exons);
	}*/
	
	public RefSeqGene getCDS(){
		Alignments cds=getCDSRegion();
		RefSeqGene rtrn=this.trimAbsolute(cds.getStart(), cds.getEnd());		
		return rtrn;
	}
	
	public Alignments getCDSRegion(){
		Alignments align=new Alignments(this.chr, this.blockStart, this.blockEnd);
		return align;
	}
	
	public Alignments getAlignment(){
		Alignments align=new Alignments(this.chr, this.start, this.stop, this.orientation);
		return align;
	}
	
	public Alignments getPromoter(int fudgeFactor)
	{ return getPromoter(fudgeFactor,fudgeFactor);}
	
	public Alignments getPromoter(int upstream,int downstream){
		Alignments align=null;
		if(this.orientation.equalsIgnoreCase("+")){
			align=new Alignments(this.chr, this.start-upstream, this.start+downstream);
		}
		else{
			align=new Alignments(this.chr, this.stop-downstream, this.stop+upstream);
		}
		align.setName(getName());
		align.setOrientation(getOrientation());
		return align;
	}
	
	public double getNumberPositions(int rate){
		return ((this.stop-this.start)/rate)+1;
		
	}
	
	public int getGappedSize(){
		int sum=0;
		Alignments[] exons=this.getExons();
		if(exons!=null){
			for(int i=0; i<exons.length; i++){
				sum+=exons[i].getSize();
			}
		}
		return sum;
	}
	
	public double getBedScore(){return bedScore;}
	
	public int get3PrimeGenomicDistance(LightweightGenomicAnnotation annot) {
		
		int res=0;
		if (this.getOrientation().equalsIgnoreCase("-"))
			res=annot.getEnd()-this.getStart();
		else
			res=this.getEnd()-annot.getStart();
			
	
		return res;
	}
	
	public int get3PrimeTranscriptDistance(LightweightGenomicAnnotation annot) {
		
		int res=0;
		//find the overlapping block
		int stBlck=0;
		int enBlck=0;
		for (int i=0;i<this.exonStarts.length; i++){
			if (annot.getStart() <= this.exonEnds[i] && annot.getStart() >= this.exonStarts[i] )
				stBlck=i;
			if (annot.getEnd() <= this.exonEnds[i] && annot.getEnd() >= this.exonStarts[i] )
				enBlck=i;
			
		}
		if (this.getOrientation().equalsIgnoreCase("-")){
			for (int i=0; i<enBlck; i++)
				res+=(this.exonEnds[i]-this.exonStarts[i]);
			res+=annot.getEnd()-this.exonStarts[enBlck];
		}
		else{
			for (int i=this.exonStarts.length-1; i>stBlck; i--)
				res+=(this.exonEnds[i]-this.exonStarts[i]);
			res+=this.exonEnds[stBlck]-annot.getStart();
		}
			
	
		return res;
	}
	
	public RefSeqGene getExtended5primeIsoform(int buffer) {
	
		RefSeqGene newGene=this.copy();
		if (this.getOrientation().equalsIgnoreCase("-"))
			newGene.setEnd(newGene.getEnd()+buffer);
		else
			newGene.setStart(Math.max(newGene.getStart()-buffer,0));
		
			
		return newGene;
	}
	
	public RefSeqGene getExtended3primeIsoform(int buffer) {
	
		RefSeqGene newGene=this.copy();
		if (this.getOrientation().equalsIgnoreCase("-"))
			newGene.setStart(Math.max(newGene.getStart()-buffer,0));
		else
			newGene.setEnd (newGene.getEnd()+buffer);
		
		return newGene;
	}
	
	public Alignments get5primeRegion(int buffer, boolean limitByFirstExon) {
		int pos=0;
		int top=0;
		int st=0;
		Alignments res=null;
		if (this.getOrientation().equalsIgnoreCase("-")){
			 pos=this.getEnd();
			 st=Math.max(pos-buffer,this.exonEnds[0]);
			 top=pos+buffer;
		}
		else{
			pos=this.getStart();
			top=Math.min(pos+buffer,this.exonStarts[0]);
			st=Math.max(pos-buffer,0);
		}
		if(limitByFirstExon){
			res=new Alignments (this.chr,st,top); 
		}
		else{
			res=new Alignments (this.chr,Math.max(pos-buffer,0),pos+buffer);
		}
		return res;
	}
	
	public RefSeqGene getOverlap(RefSeqGene second) {
		Set<Alignments> otherExons = second.getExonSet();
		List<Alignments> thisExons  = new ArrayList<Alignments>(getExonSet());
		Set<LightweightGenomicAnnotation> overlappingExons = new TreeSet<LightweightGenomicAnnotation>();
		for(Alignments exon: otherExons) {
			Alignments exonClone = new Alignments(exon);
			overlappingExons.addAll(exonClone.intersect(thisExons));
		}
		
		RefSeqGene rtrn = null;
		if(!overlappingExons.isEmpty()) {
			rtrn = new RefSeqGene(overlappingExons);
			rtrn.setOrientation(orientation);
		}
		
		return rtrn;
	}
	
	/**
	 * Get all overlaps between this gene and any gene in the collection
	 * @param others The collection of genes to check for overlap
	 * @return A RefSeqGene object whose "exons" consist of overlapping intervals between this gene and the collection
	 */
	public RefSeqGene getOverlap(Collection<RefSeqGene> others) {
		Set<Alignments> otherExons = new TreeSet<Alignments>();
		for(RefSeqGene other : others) otherExons.addAll(other.getExonSet());
		List<Alignments> thisExons  = new ArrayList<Alignments>(getExonSet());
		Set<LightweightGenomicAnnotation> overlappingExons = new TreeSet<LightweightGenomicAnnotation>();
		for(Alignments exon: otherExons) {
			Alignments exonClone = new Alignments(exon);
			overlappingExons.addAll(exonClone.intersect(thisExons));
		}
		
		RefSeqGene rtrn = null;
		if(!overlappingExons.isEmpty()) {
			rtrn = new RefSeqGene(overlappingExons);
			rtrn.setOrientation(orientation);
		}
		
		return rtrn;		
	}
	
	private int getOverlap(Alignments exon, Alignments align) {
		LightweightGenomicAnnotation intersect=exon.intersect(align);
		return intersect.length();
	}
	
	private Map getOverlapping(Collection<Alignments> alignments){
		Map<Alignments, Set> overlapping=new TreeMap();
		for(Alignments peak1: alignments){
			Set overlaps=new TreeSet();
			overlaps.add(peak1);
			Alignments ext=new Alignments(peak1.getChr(), peak1.getStart(), peak1.getEnd());
			for(Alignments peak2: alignments){
				if(ext.overlapsAtAll(peak2) || peak2.overlapsAtAll(ext)){overlaps.add(peak2);}
			}
			overlapping.put(peak1, overlaps);
		}
		
		return overlapping;
	}
	
	//position of splice junction relative to full length sequence (0 index)
	public ArrayList<Integer> getSpliceJunctionCoordinates()throws Exception{
		ArrayList<Integer> coordinates=new ArrayList();
		
		Collection<Alignments> exons=this.getExonSet();
		
		int length=0;
		int counter=0;
		for(Alignments exon: exons){
			length+=exon.getSize();
			counter++;
			if(counter<exons.size()){coordinates.add(length);}
		}
		
		
		return coordinates;
	}
	
	public IntervalTree<Integer> getSpliceJunctionCoordinatesTree() throws Exception {
		ArrayList<Integer> junctions=getSpliceJunctionCoordinates();
		IntervalTree<Integer> tree=new IntervalTree();
		
		for(Integer pos: junctions){
			tree.put(pos, pos+1, pos);
		}
		
		return tree;
	}
	
	/*public String toSAM(){
		String qualityString="255";
		String qualityLetters="";
		String mateReference="*";
		String matePosition="0";
		String insertSize="0";
		
		int strandNum=0;
		if(this.orientation.equalsIgnoreCase("-")){strandNum=16;}
		if(this.sequence==null){this.sequence="";}
		
		for(int i=0; i<sequence.toCharArray().length; i++){qualityLetters=qualityLetters+"I";}
		
		String cigar=makeCigar();
		
		String rtrn=name+"\t"+strandNum+"\t"+chr+"\t"+(start+1)+"\t"+qualityString+"\t"+cigar+"\t"+mateReference+"\t"+matePosition+"\t"+insertSize+"\t"+sequence+"\t"+qualityLetters;
		
		return rtrn;
	}
	
	public String toSAM(RefSeqGene mate){
		String qualityString="255";
		String qualityLetters="";
		String mateReference="=";
		String matePosition=new Integer(mate.getStart()+1).toString();
		String insertSize=new Integer(mate.getStart()-start).toString();
		
		int strandNum=0;
		if(this.orientation.equalsIgnoreCase("-")){strandNum=16;}
		
		for(int i=0; i<sequence.toCharArray().length; i++){qualityLetters=qualityLetters+"I";}
		
		String cigar=makeCigar();
		
		String rtrn=name+"\t"+strandNum+"\t"+chr+"\t"+(start+1)+"\t"+qualityString+"\t"+cigar+"\t"+mateReference+"\t"+matePosition+"\t"+insertSize+"\t"+sequence+"\t"+qualityLetters;
				
		return rtrn;
	}*/
	
	/*private String makeCigar(){
		String rtrn="";
		Alignments[] exons=this.getExons();
		for(int i=0; i<exons.length; i++){
			Alignments exon=(Alignments)exons[i];
			if(i!=0){
				Alignments previous=(Alignments)exons[i-1];
				Alignments intron=new Alignments(exon.getChr(), previous.getEnd(), exon.getStart());
				if(intron.getSize()>0){rtrn=rtrn+(intron.getSize())+"N";}
			}
			rtrn=rtrn+(exon.getSize())+"M";
		}
		return rtrn;
	}*/
	
	/*public String toSAM(){
		String qualityString="255";
		String qualityLetters="";
		String mateReference="*";
		String matePosition="0";
		String insertSize="0";
		
		int strandNum=0;
		if(this.orientation.equalsIgnoreCase("-")){strandNum=16;}
		if(this.sequence==null){this.sequence="";}
		
		for(int i=0; i<sequence.toCharArray().length; i++){qualityLetters=qualityLetters+"I";}
		
		String cigar=makeCigar();
		
		String rtrn=name+"\t"+strandNum+"\t"+chr+"\t"+(start+1)+"\t"+qualityString+"\t"+cigar+"\t"+mateReference+"\t"+matePosition+"\t"+insertSize+"\t"+sequence+"\t"+qualityLetters;
		
		return rtrn;
	}
	
	public String toSAM(RefSeqGene mate){
		String qualityString="255";
		String qualityLetters="";
		String mateReference="=";
		String matePosition=new Integer(mate.getStart()+1).toString();
		String insertSize=new Integer(mate.getStart()-start).toString();
		
		int strandNum=0;
		if(this.orientation.equalsIgnoreCase("-")){strandNum=16;}
		
		for(int i=0; i<sequence.toCharArray().length; i++){qualityLetters=qualityLetters+"I";}
		
		String cigar=makeCigar();
		
		String rtrn=name+"\t"+strandNum+"\t"+chr+"\t"+(start+1)+"\t"+qualityString+"\t"+cigar+"\t"+mateReference+"\t"+matePosition+"\t"+insertSize+"\t"+sequence+"\t"+qualityLetters;
				
		return rtrn;
	}*/
	
	/*private String makeCigar(){
		String rtrn="";
		Alignments[] exons=this.getExons();
		for(int i=0; i<exons.length; i++){
			Alignments exon=(Alignments)exons[i];
			if(i!=0){
				Alignments previous=(Alignments)exons[i-1];
				Alignments intron=new Alignments(exon.getChr(), previous.getEnd(), exon.getStart());
				if(intron.getSize()>0){rtrn=rtrn+(intron.getSize())+"N";}
			}
			rtrn=rtrn+(exon.getSize())+"M";
		}
		return rtrn;
	}*/
	
	public String toString(){
		return this.toBED();
		//return this.name;
	}
	
	public String toRefSeq(){
		Set<Alignments> exons = getExonSet();
		String rtrn=name+"\t"+(name == null ? toUCSC() : name)+"\t"+chr+"\t"+orientation+"\t"+start+"\t"+stop+"\t"+start+"\t"+stop+"\t"+exons.size();
		String starts="";
		String ends="";
		for(Alignments exon : exons) {
			starts+=exon.getStart();
			ends+=exon.getEnd();
		}
		rtrn=rtrn+"\t"+starts+"\t"+ends;
		return rtrn;
	}
	
	/*public String toSAM(){
		String qualityString="255";
		String qualityLetters="";
		String mateReference="*";
		String matePosition="0";
		String insertSize="0";
		
		int strandNum=0;
		if(this.orientation.equalsIgnoreCase("-")){strandNum=16;}
		if(this.sequence==null){this.sequence="";}
		
		for(int i=0; i<sequence.toCharArray().length; i++){qualityLetters=qualityLetters+"I";}
		
		String cigar=makeCigar();
		
		String rtrn=name+"\t"+strandNum+"\t"+chr+"\t"+(start+1)+"\t"+qualityString+"\t"+cigar+"\t"+mateReference+"\t"+matePosition+"\t"+insertSize+"\t"+sequence+"\t"+qualityLetters;
		
		return rtrn;
	}
	
	public String toSAM(RefSeqGene mate){
		String qualityString="255";
		String qualityLetters="";
		String mateReference="=";
		String matePosition=new Integer(mate.getStart()+1).toString();
		String insertSize=new Integer(mate.getStart()-start).toString();
		
		int strandNum=0;
		if(this.orientation.equalsIgnoreCase("-")){strandNum=16;}
		
		for(int i=0; i<sequence.toCharArray().length; i++){qualityLetters=qualityLetters+"I";}
		
		String cigar=makeCigar();
		
		String rtrn=name+"\t"+strandNum+"\t"+chr+"\t"+(start+1)+"\t"+qualityString+"\t"+cigar+"\t"+mateReference+"\t"+matePosition+"\t"+insertSize+"\t"+sequence+"\t"+qualityLetters;
				
		return rtrn;
	}*/
	
	/*private String makeCigar(){
		String rtrn="";
		Alignments[] exons=this.getExons();
		for(int i=0; i<exons.length; i++){
			Alignments exon=(Alignments)exons[i];
			if(i!=0){
				Alignments previous=(Alignments)exons[i-1];
				Alignments intron=new Alignments(exon.getChr(), previous.getEnd(), exon.getStart());
				if(intron.getSize()>0){rtrn=rtrn+(intron.getSize())+"N";}
			}
			rtrn=rtrn+(exon.getSize())+"M";
		}
		return rtrn;
	}*/
	
	public String toUCSC() {
		return getChr()+":" +getStart() + "-" + getEnd();
	}
	
	public String toBED (){
		return toBED(true);
	}
	
	public String toBED(int r, int g, int b) {
		return toBED(true, r, g, b);
	}
	
	/**
	 * Write a set of genes to a bed file
	 * @param genes The set of genes
	 * @param outFile The output file
	 * @throws IOException
	 */
	public static void writeBedFile(Collection<RefSeqGene> genes, String outFile) throws IOException {
		FileWriter w = new FileWriter(outFile);
		for(RefSeqGene gene : genes) w.write(gene.toBED() + "\n");
		w.close();
	}
	
	public String toBED(boolean useExtraFields){
		return toBED(true, 0, 0, 0);
	}
	
	public String toBED(boolean useExtraFields, int r, int g, int b){
		if(r < 0 || r > 255 || g < 0 || g > 255 || b < 0 || b > 255) {
			throw new IllegalArgumentException("RGB values must be between 0 and 255");
		}
		String rgb = r + "," + g + "," + b;
		Set<Alignments> exons = getExonSet();
		//String rtrn=this.getChr()+"\t"+this.getStart()+"\t"+this.getEnd()+"\t"+(name == null ? toUCSC() : this.name)+"\t"+getCountScore()+"\t"+orientation+"\t"+start+"\t"+stop+"\t"+rgb+"\t"+exons.size();
		String rtrn=this.getChr()+"\t"+this.getStart()+"\t"+this.getEnd()+"\t"+(name == null ? toUCSC() : this.name)+"\t"+this.bedScore+"\t"+orientation+"\t"+this.blockStart+"\t"+this.blockEnd+"\t"+rgb+"\t"+exons.size();
		String sizes="";
		String starts="";
		for(Alignments exon : exons){
			sizes=sizes+(exon.length())+",";
			starts=starts+(exon.getStart()-getStart())+",";
		}
		rtrn=rtrn+"\t"+sizes+"\t"+starts;
		
		if(extraFields != null & useExtraFields) {
			for(String field : extraFields) {
				rtrn = rtrn + "\t" + field;
			}
		}
		return rtrn;
	}
	
	public String toBEDwithBedScore(){
		return toBED();
	}
	
	public String[] toBEDArr() {
	
		String[] rtrn; ;
		if(extraFields != null)
			 rtrn=new String[12+extraFields.length];
		else
			 rtrn=new String[12];
		
		rtrn[0]= this.getChr();
		rtrn[1]= String.valueOf(this.getStart());
		rtrn[2]= String.valueOf(this.getEnd());
		rtrn[3]= this.getName();
		rtrn[4]= String.valueOf(this.getBedScore());
		rtrn[5]= this.getOrientation();
		rtrn[6]= String.valueOf(this.getStart());
		rtrn[7]= String.valueOf(this.getEnd());
		rtrn[8]= "0,0,0";
		rtrn[9]= String.valueOf(this.exonStarts.length);
		
		String sizes="";
		String starts="";
		for(int i=0; i<this.exonStarts.length; i++){
			sizes=sizes+(this.exonEnds[i]-this.exonStarts[i])+",";
			starts=starts+(this.exonStarts[i]-this.getStart())+",";
		}
		rtrn[10]=sizes;
		rtrn[11]=starts;
		
		if(extraFields != null) {
			int i=12;
			for(String field : extraFields) {
				rtrn[i] =  field;
				i++;
			}
		}
		return rtrn;
	}
	
	public String toGTF(String source){
		return this.toGTF(source,this.getName(),getName() +".0");
	
	}
	
	public String toGTF(String source,String geneName ,String transcriptName){
		StringBuilder rtrn = new StringBuilder();
		Set<Alignments> exons = getExonSet();
		GFF transcriptGFF = new GFF(getName());
		transcriptGFF.setFeature("mRNA");
		transcriptGFF.setSource(source);
		transcriptGFF.setOrientation(this.getOrientation());
		transcriptGFF.addAttribute("gene_id", transcriptName);
		if(attributes != null && !attributes.isEmpty()) {
			for(String attr : attributes.keySet()) {
				transcriptGFF.addAttribute(attr, String.valueOf(attributes.get(attr)));
			}
		}
		//Add +1 when we print GFF as it is 1 based rather than 0 based
		transcriptGFF.setStart(getStart());  // BUG FIX : MORAN AUG 17TH, ADD +1 ; 2nd bug fix - 11.29.10 only during print we add +1 to start pos
		transcriptGFF.setEnd(getEnd()); //BUG FIX : MORAN AUG 17TH, END IS CORRECT AS GFF IS INCLUSIVE BUT ALSO USES 1 BASE CORRDINATES (WE USE 0 BASED)
		transcriptGFF.setChromosome(getChr());
		rtrn.append(transcriptGFF.toString(true));
		rtrn.append("\n");
		int i = 0;
		for(Alignments exon : exons) {
			GFF exonGFF = new GFF(exon);
			exonGFF.setFeature("exon");
			exonGFF.setSource(source);
			exonGFF.setName(exon.getName());
			exonGFF.addAttribute("gene_id", exon.getName());
			exonGFF.addAttribute("transcript_id", transcriptName);
			exonGFF.addAttribute("Parent",   transcriptName);
			if(attributes != null && !attributes.isEmpty()) {
				for(String attr : attributes.keySet()) {
					exonGFF.addAttribute(attr, String.valueOf(attributes.get(attr)));
				}
			}
			rtrn.append(exonGFF.toString(true));
			rtrn.append("\n");
		}
	
		return rtrn.toString();
	}
	
	public String toMISO(String parentId, String id, String source) {
		StringBuilder rtrn = new StringBuilder();
		Set<Alignments> exons = getExonSet();
		GFF transcriptGFF = new GFF(id);
		transcriptGFF.setFeature("mRNA");
		transcriptGFF.setSource(source);
		transcriptGFF.setOrientation(this.getOrientation());
		transcriptGFF.addAttribute("ID", id);
		transcriptGFF.addAttribute("Parent", parentId);
	
		//Add +1 when we print GFF as it is 1 based rather than 0 based
		transcriptGFF.setStart(getStart());  // BUG FIX : MORAN AUG 17TH, ADD +1 ; 2nd bug fix - 11.29.10 only during print we add +1 to start pos
		transcriptGFF.setEnd(getEnd()); //BUG FIX : MORAN AUG 17TH, END IS CORRECT AS GFF IS INCLUSIVE BUT ALSO USES 1 BASE CORRDINATES (WE USE 0 BASED)
		transcriptGFF.setChromosome(getChr());
		rtrn.append(transcriptGFF.toString(false, true));
		int i = 0;
		for(Alignments exon : exons) {
			rtrn.append("\n");
			GFF exonGFF = new GFF(exon);
			exonGFF.setFeature("exon");
			exonGFF.setSource(source);
			exonGFF.setName(exon.getName());
			exonGFF.addAttribute("ID", id+"-"+i);
			exonGFF.addAttribute("Parent",   id);
			rtrn.append(exonGFF.toString(false, true));
			i++;
		}
	
		return rtrn.toString();
	}
	
	/*public String toSAM(){
		String qualityString="255";
		String qualityLetters="";
		String mateReference="*";
		String matePosition="0";
		String insertSize="0";
		
		int strandNum=0;
		if(this.orientation.equalsIgnoreCase("-")){strandNum=16;}
		if(this.sequence==null){this.sequence="";}
		
		for(int i=0; i<sequence.toCharArray().length; i++){qualityLetters=qualityLetters+"I";}
		
		String cigar=makeCigar();
		
		String rtrn=name+"\t"+strandNum+"\t"+chr+"\t"+(start+1)+"\t"+qualityString+"\t"+cigar+"\t"+mateReference+"\t"+matePosition+"\t"+insertSize+"\t"+sequence+"\t"+qualityLetters;
		
		return rtrn;
	}
	
	public String toSAM(RefSeqGene mate){
		String qualityString="255";
		String qualityLetters="";
		String mateReference="=";
		String matePosition=new Integer(mate.getStart()+1).toString();
		String insertSize=new Integer(mate.getStart()-start).toString();
		
		int strandNum=0;
		if(this.orientation.equalsIgnoreCase("-")){strandNum=16;}
		
		for(int i=0; i<sequence.toCharArray().length; i++){qualityLetters=qualityLetters+"I";}
		
		String cigar=makeCigar();
		
		String rtrn=name+"\t"+strandNum+"\t"+chr+"\t"+(start+1)+"\t"+qualityString+"\t"+cigar+"\t"+mateReference+"\t"+matePosition+"\t"+insertSize+"\t"+sequence+"\t"+qualityLetters;
				
		return rtrn;
	}*/
	
	/*private String makeCigar(){
		String rtrn="";
		Alignments[] exons=this.getExons();
		for(int i=0; i<exons.length; i++){
			Alignments exon=(Alignments)exons[i];
			if(i!=0){
				Alignments previous=(Alignments)exons[i-1];
				Alignments intron=new Alignments(exon.getChr(), previous.getEnd(), exon.getStart());
				if(intron.getSize()>0){rtrn=rtrn+(intron.getSize())+"N";}
			}
			rtrn=rtrn+(exon.getSize())+"M";
		}
		return rtrn;
	}*/
	
	//MORAN: added GTF in cufflinks format without parent
	//Note: To use the attributes transcript_id or gene_id  as is, set geneName and transcriptNme to empy strings
	//else the name will be concatenated to the attributes values
	public String toCufflinksGTF(String source,String geneName ,String transcriptName,String nonParsedAttrs){
		StringBuilder rtrn = new StringBuilder();
		Set<Alignments> exons = getExonSet();
		GFF transcriptGFF = new GFF(getName());
		transcriptGFF.setSource(source);
		transcriptGFF.setOrientation(this.getOrientation());
		transcriptGFF.setFeature("transcript");//As in version 1353 
		
		if ( nonParsedAttrs.isEmpty()){
			transcriptGFF.clearAttribute("gene_id");
			if (! geneName.isEmpty())
				transcriptGFF.addAttribute("gene_id", geneName);//As in version 1353 to compile with cufflinks
			transcriptGFF.clearAttribute("transcript_id");
			if (! transcriptName.isEmpty())
				transcriptGFF.addAttribute("transcript_id", transcriptName);//
			if(attributes != null && !attributes.isEmpty()) {
				for(String attr : attributes.keySet()) {
					transcriptGFF.addAttribute(attr, String.valueOf(attributes.get(attr)));
				}
			}
		}
		else {
			transcriptGFF.setAttributes(nonParsedAttrs);
		}
		//Add +1 when we print GFF as it is 1 based rather than 0 based
		transcriptGFF.setStart(getStart());  //BUG FIX : MORAN AUG 17TH, ADD +1; 2nd bug fix - 11.29.10 only during print we add +1 to start pos
		transcriptGFF.setEnd(getEnd()); //BUG FIX : MORAN AUG 17TH, END IS CORRECT AS GFF IS INCLUSIVE BUT ALSO USES 1 BASE CORRDINATES (WE USE 0 BASED)
		transcriptGFF.setChromosome(getChr());
		rtrn.append(transcriptGFF.toCufflinksString(true));
		rtrn.append("\n");
		
		for(Alignments exon : exons) {
			GFF exonGFF = new GFF(exon);
			exonGFF.setFeature("exon");
			exonGFF.setSource(source);
			
			if ( nonParsedAttrs.isEmpty()){
				exonGFF.setName(exon.getName());
				exonGFF.clearAttribute("gene_id");
				if (! geneName.isEmpty())
					exonGFF.addAttribute("gene_id", geneName); //As in version 1353 to compile with cufflinks
				exonGFF.clearAttribute("transcript_id");
				if (! transcriptName.isEmpty())
					exonGFF.addAttribute("transcript_id", transcriptName);
				//exonGFF.addAttribute("Parent",   transcriptName);
				if(attributes != null && !attributes.isEmpty()) {
					for(String attr : attributes.keySet()) {
						exonGFF.addAttribute(attr, String.valueOf(attributes.get(attr)));
					}
				}
			}
			else
				exonGFF.setAttributes(nonParsedAttrs);
			
			rtrn.append(exonGFF.toCufflinksString(true));
			rtrn.append("\n");
		}
	
		return rtrn.toString();
	}
	
	//Fixed off by one error
	public String toSAM(){
		String qualityString="255";
		String qualityLetters="";
		String mateReference="*";
		String matePosition="0";
		String insertSize="0";
		
		int strandNum=0;
		if(this.orientation.equalsIgnoreCase("-")){strandNum=16;}
		if(this.sequence==null){this.sequence="";}
		
		for(int i=0; i<sequence.toCharArray().length; i++){qualityLetters=qualityLetters+"I";}
		
		String cigar=makeCigar();
		
		String rtrn=name+"\t"+strandNum+"\t"+chr+"\t"+(start+1)+"\t"+qualityString+"\t"+cigar+"\t"+mateReference+"\t"+matePosition+"\t"+insertSize+"\t"+sequence+"\t"+qualityLetters;
		
		return rtrn;
	}
	
	public String toSAM(String originalSAMLine) {
		return toSAM(originalSAMLine, false);
	}
	
	public String toSAM(String originalSAMLine, boolean flip) {
		String[] tokens=originalSAMLine.split("\t");
		//all we want to do is update the relevant columns which are the chromosome position and cigar string
		//orientation is [1]
		//chr is [2]
		//position is [3]
		//cigar is [5]
		//sequence is [9]
		//quality is [10]
		
		String sequence=tokens[9];
		//String quality=tokens[10];
		String orientation=tokens[1];
		if(flip){
			if(orientation.equalsIgnoreCase("0")){orientation="16";}
			else if(orientation.equalsIgnoreCase("16")){orientation="0";}
			else{System.err.println(orientation);}
			sequence=Sequence.reverseSequence(sequence);
			//quality=reverse(quality);
		}
		 
		String rtrn="";
		for(int i=0; i<tokens.length; i++){
			if(i!=1 && i!=2 && i!=3 && i!=5 && i!=9){rtrn+=tokens[i]+"\t";}
			else if(i==1){rtrn+=orientation+"\t";}
			else if(i==2){rtrn+=this.getChr()+"\t";}
			else if(i==3){rtrn+=(this.start+1)+"\t";}
			else if(i==5){String cigar=makeCigar(); rtrn+=cigar+"\t";}
			else if(i==9){rtrn+=sequence+"\t";}
			//else if(i==10){rtrn+=quality+"\t";}
		}
		return rtrn;
	}
	
	public void setCDS(Alignments cdsRegion){
		this.blockStart=cdsRegion.getStart();
		this.blockEnd=cdsRegion.getEnd();
	}
	
	public int absoluteToRelativePosition(int position) {
		Collection<Alignments> exons=this.getExonSet();
		
		int count=0;
		for(Alignments exon: exons){
			if(exon.fullyContained(new Alignments(chr, position, position))){
				count+=(position-exon.getStart()-1);
			}
			else if(exon.getEnd()<position){count+=exon.length();}
		}
		
		return count;
	}
	public boolean isInExon(int snp) {
		Collection<Alignments> exons=this.getExonSet();
		
		
		for(Alignments exon: exons){
			if(exon.getStart()<snp && exon.getEnd()>snp){return true;}
		}
		return false;
	}
	
	/**
	 * Get all ORFs as gene objects
	 * @param trim3UTRs whether to trim 3'UTRs to end at the beginning of next ORF
	 * @return the collection of ORFs
	 */
	public Collection<RefSeqGene> findAllORFs(boolean trim3UTRs) {
		
		if(trim3UTRs) return findAllOrfsConservative3UTRs();
		
		Collection<int[]> orfs = findAllORFs(sequence);
		String cdsOrientation = orientation;
		if("*".equals(orientation)) {
			//TODO If unknown strand get all 6 frames
			/*Collection<int[]> reverseOrfs = findAllORF(Sequence.reverseSequence(sequence));
			cdsOrientation = "+";
			if(reverseOrf[1] - reverseOrf[0] > orf[1]- orf[0]) {
				orf[0] = reverseOrf[0];
				orf[1] = reverseOrf[1];
				cdsOrientation = "-";
			}*/
			System.err.println("Not finding ORFs in gene " + name + " because strand is unknown.");
		}
		
		Collection<RefSeqGene> rtrn=new TreeSet<RefSeqGene>();
		for(int[] orf: orfs){
			int orfGenomicStart = "-".equals(cdsOrientation) ? mapToGenomic(orf[1], cdsOrientation):  mapToGenomic(orf[0], cdsOrientation)  ; 
			int orfGenomicEnd   = "-".equals(cdsOrientation) ? mapToGenomic(orf[0], cdsOrientation):  mapToGenomic(orf[1], cdsOrientation);
			
			//System.err.println("ORF genomic start end " + orfGenomicStart +"-" + orfGenomicEnd);
			
			LightweightGenomicAnnotation geneCDS = new BasicLightweightAnnotation(getChr(), orfGenomicStart, orfGenomicEnd);
			Set<LightweightGenomicAnnotation> cdsExons = new TreeSet<LightweightGenomicAnnotation>();
			Set<Alignments> exonSet = getExonSet();
			for(Alignments e : exonSet) {
				if(geneCDS.overlaps(e)) {
					cdsExons.add(e.intersect(geneCDS));
				}
			}
			
			RefSeqGene cds = null;
			if(cdsExons.size() > 0 ) {
				cds = new RefSeqGene(cdsExons);
				cds.setOrientation(cdsOrientation);
				cds.setName(getName());
				cds.setCountScore((orf[1]- orf[0])/(double)(sequence.length()));
			}
			
			//TODO Add 5' UTR and 3' UTR to each ORF
			//cds.setCDS(cds.getAlignment());
			RefSeqGene testCDS=this.copy();
			testCDS.setCDS(cds.getAlignment());
						
			//rtrn.add(cds);
			rtrn.add(testCDS);
		}
		return rtrn;
	}
	
	/**
	 * Get all ORFs as genes, with 3'UTRs trimmed so they end at the beginning of next downstream ORF
	 * @return The collection of all ORFs
	 */
	private Collection<RefSeqGene> findAllOrfsConservative3UTRs() {
		
		if("*".equals(orientation)) {
			System.err.println("Not finding ORFs in gene " + name + " because strand is unknown.");
		}
		
		//System.err.println("Finding ORFs with shortened 3'UTRs for gene " + name + " " + chr + ":" + start + "-" + stop);
		
		TreeSet<RefSeqGene> orfs = (TreeSet<RefSeqGene>) findAllORFs(false);
		
		Collection<RefSeqGene> rtrn = new TreeSet<RefSeqGene>();
		
		if(orientation.equals("+")) {
			
			Iterator<RefSeqGene> iter = orfs.iterator();
			if(!iter.hasNext()) return rtrn;
			
			while(iter.hasNext()) {
				RefSeqGene orf = iter.next();
				int cdsEnd = orf.blockEnd;
				int endPos = orf.stop;
				Iterator<RefSeqGene> after = orfs.tailSet(orf).iterator();
				while(after.hasNext()) {
					int nextStart = after.next().blockStart;
					if(nextStart > cdsEnd) {
						endPos = nextStart;
						break;
					}
				}
				//System.err.println("Setting end from " + orf.stop + " to " + endPos);
				RefSeqGene trimmed = orf.trimAbsolute(orf.getStart(), endPos);			
				//System.err.println("New end is " + trimmed.stop);
				if(trimmed != null) {
					rtrn.add(trimmed);
					//System.err.println("Added ORF with CDS " + trimmed.blockStart + "-" + trimmed.blockEnd);
				}
			}
		}
		
		if(orientation.equals("-")) {
			Iterator<RefSeqGene> iter = orfs.descendingIterator();
			if(!iter.hasNext()) return rtrn;
			while(iter.hasNext()) {
				RefSeqGene orf = iter.next();
				int cdsStart = orf.blockStart;
				int startPos = orf.start;
				Iterator<RefSeqGene> after = ((TreeSet) orfs.headSet(orf)).descendingIterator();
				while(after.hasNext()) {
					int nextEnd = after.next().blockEnd;
					if(nextEnd < cdsStart) {
						startPos = nextEnd;
						break;
					}
				}
				//System.err.println("Setting start from " + orf.start + " to " + startPos);
				RefSeqGene trimmed = orf.trimAbsolute(startPos, orf.getEnd());			
				//System.err.println("New start is " + trimmed.start);
				if(trimmed != null) {
					rtrn.add(trimmed);
					//System.err.println("Added ORF with CDS " + orf.blockStart + "-" + orf.blockEnd);
				}
			}
		}
		
		return rtrn;
		
	}
	
	/**
	 * Find start codons
	 * @return collection of all start codon genomic coordinates
	 */
	public TreeSet<BasicLightweightAnnotation> findAllStartCodons() {
		Collection<int[]> orfs = findAllORFs(sequence);
		String cdsOrientation = orientation;
		if("*".equals(orientation)) {
			System.err.println("Skipping because we dont know strand");
			return null;
		}
		
		TreeSet<BasicLightweightAnnotation> rtrn = new TreeSet<BasicLightweightAnnotation>();
		for(int[] orf: orfs){
			int startCodonGenomicStart = "-".equals(cdsOrientation) ? mapToGenomic(orf[1], cdsOrientation):  mapToGenomic(orf[0], cdsOrientation)  ; 
			int startCodonGenomicEnd   = "-".equals(cdsOrientation) ? mapToGenomic(orf[1] + 2, cdsOrientation):  mapToGenomic(orf[0] + 2, cdsOrientation);
			
			BasicLightweightAnnotation startCodon = new BasicLightweightAnnotation(getChr(), startCodonGenomicStart, startCodonGenomicEnd);
			rtrn.add(startCodon);
			System.err.println("Found start codon at " + getChr() + " " + startCodonGenomicStart + "-" + startCodonGenomicEnd);
		}
		return rtrn;
		
	}
	
	
	public static Collection<int[]> findAllORFs(String sequence) {
		Collection<int[]> allORFs=new ArrayList<int[]>();
		int lastStart = 0;
		int lastEnd = 0;
		Matcher m = START_CODON_PATRN.matcher(sequence);
		while(m.find()) {
			int startCodonPos = m.start(); //start codon
			int thisORFEnd = startCodonPos;
			boolean foundStopCodon = false;
			//System.err.println("new orf start: " + startCodonPos);
			while(thisORFEnd < sequence.length() - 3  && !foundStopCodon) {
				String currentCodon = sequence.substring(thisORFEnd, thisORFEnd+3);
				//System.err.print(" "+currentCodon+ " ");
				thisORFEnd += 3;
				foundStopCodon = isStopCodon(currentCodon);
				//int[] pos={startCodonPos, thisORFEnd};
				//allORFs.add(pos);
			}
			if(foundStopCodon){
				int[] pos={startCodonPos, thisORFEnd};
				allORFs.add(pos);
			}
			//System.err.println("Orf length: " + (thisORFEnd - startCodonPos) );
			if(lastEnd - lastStart < thisORFEnd - startCodonPos) {
				lastStart = startCodonPos;
				lastEnd   = thisORFEnd;
				//System.err.println("It was a winner");
			}
		}
		int [] startEnd = {lastStart, lastEnd};
		return allORFs;
	}
	
	
	
	public Collection<RefSeqGene> findAllORFs(Sequence chrSeq, boolean trim3UTRs) {
		this.setSequenceFromChromosome(chrSeq);
		return findAllORFs(trim3UTRs);
	}
	
	
	public String toRefFlat (){
		
		Set<Alignments> exons = getExonSet();
		//String rtrn=this.getChr()+"\t"+this.getStart()+"\t"+this.getEnd()+"\t"+(name == null ? toUCSC() : this.name)+"\t"+getCountScore()+"\t"+orientation+"\t"+start+"\t"+stop+"\t"+"0,0,0"+"\t"+exons.size();
		String rtrn=this.getName()+"\t"+this.getName()+"\t"+this.getChr()+"\t"+orientation+"\t"+this.getStart()+"\t"+this.getEnd()+"\t"+start+"\t"+stop+"\t"+exons.size()+"\t";
		
		for(Alignments exon : exons){
			rtrn= rtrn+ exon.getStart()+",";
		}
		
		rtrn=rtrn+"\t";
		
		for(Alignments exon : exons){
			rtrn= rtrn+ exon.getEnd()+",";
		}
		
		
		return rtrn;
		
	}
	
	public void addSequence(Sequence seq){
		this.sequence=seq.getSequenceBases();
	}
	
	public void addSuffixToName(String refName) {
		this.setName(this.getName()+refName);
		
	}
	
	public void addAttribute(String name, String val){
		if(attributes == null) {
			attributes = new HashMap<String, String>();
		}
		attributes.put(name, val);
	}
	
	public void addExtraField (String value) {
		String [] newExtraFields = null;
		if(extraFields != null) {
			newExtraFields = new String[extraFields.length + 1];
			for(int i = 0; i < extraFields.length; i++) {
				newExtraFields[i] = extraFields[i];
			}
		} else {
			newExtraFields = new String [1];
		}
		
		newExtraFields[newExtraFields.length - 1] = value;
		
		this.extraFields = newExtraFields;
	}
	
	public void expandUtrs(Integer utr1, Integer utr2) {
		this.setStart(this.start-utr1);
		this.setEnd(this.stop+utr2);
	}
	
	public void updateLastExonWithNewEnd(int newEnd) {
		if (this.orientation.equals("-")) {
			this.start = newEnd;
			this.blockStart = newEnd;
			this.exonStarts[0] = newEnd;
		} else {
			this.stop = newEnd;
			this.blockEnd = newEnd;
			this.exonEnds[this.exonEnds.length - 1] = newEnd;
		}
		
		this.isSorted = false;
		//this.name = this.chr + ":" + String.valueOf(this.start) + "-" + String.valueOf(this.stop);
	}
	
	public void updateScrToBedScore() {
		this.countScore=this.bedScore;
	}
	
	public void updateFirstExonWithNewStart(int newStart) {
		if (this.orientation.equals("-")) {
			this.stop = newStart;
			this.blockEnd = newStart;
			this.exonEnds[this.exonEnds.length - 1] = newStart;
		} else {
			this.start = newStart;
			this.blockStart = newStart;
			this.exonStarts[0] = newStart;
		}
		
		this.isSorted = false;
		//this.name = this.chr + ":" + String.valueOf(this.start) + "-" + String.valueOf(this.stop);
	}
	
	public void updatePicardSAM (SAMRecord picardSAM, boolean flip, int referenceSequenceIdx) {
		String sequence=picardSAM.getReadString();
		//String quality=tokens[10];
		boolean reversedOrientation=picardSAM.getReadNegativeStrandFlag();
		if(flip){
			picardSAM.setReadString(Sequence.reverseSequence(sequence));
			picardSAM.setReadNegativeStrandFlag(!reversedOrientation);
		}
		
		picardSAM.setAlignmentStart(getStart()+1);
		//picardSAM.setAlignmentEnd(getEnd()+1);
		picardSAM.setCigarString(makeCigar());
		picardSAM.setReferenceName(getChr());
		picardSAM.setReferenceIndex(referenceSequenceIdx);
	}
	
	public static void main(String[] args){
		Collection<Alignments> exons=new TreeSet();
		exons.add(new Alignments("chr19", 500, 1000));
		exons.add(new Alignments("chr19", 1500, 2000));
		exons.add(new Alignments("chr19", 2500,3000));
		
		RefSeqGene gene=new RefSeqGene(exons);
		System.err.println("Trimmed "+gene.trim(5, 5));
	}
	
	/*private String makeCigar(){
		String rtrn="";
		Alignments[] exons=this.getExons();
		for(int i=0; i<exons.length; i++){
			Alignments exon=(Alignments)exons[i];
			if(i!=0){
				Alignments previous=(Alignments)exons[i-1];
				Alignments intron=new Alignments(exon.getChr(), previous.getEnd(), exon.getStart());
				if(intron.getSize()>0){rtrn=rtrn+(intron.getSize())+"N";}
			}
			rtrn=rtrn+(exon.getSize())+"M";
		}
		return rtrn;
	}*/
	
	//if the gene has an exon that overlaps the specified exon
	public boolean hasExon(Alignments exon) {
		if(exon==null || !exon.getChr().equalsIgnoreCase(chr)){return false;}
		IntervalTree<Alignments> exonTree=CollapseByIntersection.makeIntervalTree(getExonSet()).get(exon.getChr());
		return exonTree.overlappers(exon.getStart(), exon.getEnd()).hasNext();
	}
	
	public boolean contains(RefSeqGene gene){
		return gene.getStart() >= getStart() && gene.getEnd() <= getEnd();
	}
	
	public boolean containsAttribute(String key) {
		// TODO Auto-generated method stub
		return attributes != null && attributes.containsKey(key);
	}
	
	public static boolean isStopCodon(String currentCodon) {
		boolean isStopCodon = false;
		for( String stop : STOP_CODONS) {
			if(currentCodon.equalsIgnoreCase(stop)) {
				isStopCodon = true;
				break;
			}
		}
		return isStopCodon;
	}
	
	private boolean isDone(Map<Alignments, Set> overlapping){
		for(Alignments align: overlapping.keySet()){
			int size=overlapping.get(align).size();
			if(size>1){return false;}
		}
		return true;
	}
	
	public boolean isUnoriented() {
		return !"-".equals(orientation) && !"+".equals(orientation) ;
	}
	
	/*private String makeCigar(){
		String rtrn="";
		Alignments[] exons=this.getExons();
		for(int i=0; i<exons.length; i++){
			Alignments exon=(Alignments)exons[i];
			if(i!=0){
				Alignments previous=(Alignments)exons[i-1];
				Alignments intron=new Alignments(exon.getChr(), previous.getEnd(), exon.getStart());
				if(intron.getSize()>0){rtrn=rtrn+(intron.getSize())+"N";}
			}
			rtrn=rtrn+(exon.getSize())+"M";
		}
		return rtrn;
	}*/
	
	public boolean isFullyCompatible(RefSeqGene iso) {
	
		boolean res=false;
		RefSeqGene myIntron = getIntronTranscript();
		RefSeqGene isoIntron=iso.getIntronTranscript();
		if (myIntron==null | isoIntron==null)
			return false;
		if (myIntron.isEqualStructure(isoIntron))
			res=true;
		return res;
	}
	
	private boolean isEqualStructure(RefSeqGene iso) {
	
		boolean res=true;
		if (this.getNumExons()!=iso.getNumExons())
			return false;
		Alignments[] myEx= this.getExons();
		Alignments[] isoEx= iso.getExons();
		for (int i=0; i<this.getNumExons();i++){
			if (!(myEx[i].equals(isoEx[i])))
				res=false;
				
		}
		
		return res;
	}
	
	public boolean overlapsExon(RefSeqGene gen) {
		Collection<Alignments> exonSet=this.getExonSet();
		for(Alignments align: exonSet){
			if(align.overlaps(new BasicGenomicAnnotation(gen.getName(),gen.getChr(),gen.getStart(),gen.getEnd())) ) {return true;}
		}
		return false;
	}
	
	public boolean overlapsExon(LightweightGenomicAnnotation exon){
		Collection<Alignments> exonSet=this.getExonSet();
		
		for(Alignments align: exonSet){
			if(align.overlaps(exon)){return true;}
		}
		
		return false;
	}
	
	//Overlaps in the exon level ; i.e at least one exon of "this" overlap an exon of other 
	public boolean overlaps(RefSeqGene other) {
		//System.err.println("Findig overlap betweein this " + toBED() + "\nand\n"+other.toBED());
		boolean overlaps = false;
		//System.err.println("\t\t this gene's orientation [" + orientation + "] others [" + other.getOrientation()+ "]");
		if(orientation.equals(other.getOrientation()) || "*".equals(orientation) || "*".equals(other.getOrientation())) {
			//System.err.println("\t\tOrientation is compatible");
			Collection<Alignments> exons = getSortedAndUniqueExons();
			Collection<Alignments> otherExons = other.getSortedAndUniqueExons();
			Iterator<Alignments> otherExonsIt  = otherExons.iterator();
			
			while (!overlaps && otherExonsIt.hasNext()) {
				Alignments otherExon = otherExonsIt.next();
				Iterator<Alignments> exonIt  = exons.iterator();
				while(!overlaps && exonIt.hasNext()) {
					Alignments exon = exonIt.next();
					overlaps = exon.overlaps(otherExon);
					if(exon.compareTo(otherExon) > 0) {
						break;
					}
					//System.err.println("exon " + exon.toUCSC() + " otherExon " + otherExon.toUCSC() + " overlaps? " + overlaps + " exons size " + exons.size() + " other exons size " + otherExons.size());
				}
			}
		}
		return overlaps;
	}
	
	/**
	 * Whether this gene overlaps any gene in the collection at the exon level
	 * @param others The genes to check for overlap
	 * @return True iff this gene overlaps at least one gene in the collection at the exon level
	 */
	public boolean overlaps(Collection<RefSeqGene> others) {
		for(RefSeqGene other : others) {
			if(overlaps(other)) return true;
		}
		return false;
	}
	
	/*private String makeCigar(){
		String rtrn="";
		Alignments[] exons=this.getExons();
		for(int i=0; i<exons.length; i++){
			Alignments exon=(Alignments)exons[i];
			if(i!=0){
				Alignments previous=(Alignments)exons[i-1];
				Alignments intron=new Alignments(exon.getChr(), previous.getEnd(), exon.getStart());
				if(intron.getSize()>0){rtrn=rtrn+(intron.getSize())+"N";}
			}
			rtrn=rtrn+(exon.getSize())+"M";
		}
		return rtrn;
	}*/
	
	//asks if the whole genomic regions of the two genes overlap in the same orientation
	//if one gene does have a specified orientation, assume they 2 genes are in the same orientation
	public boolean overlapsGene(RefSeqGene gen) {
		if (  this.getOrientation().equals(gen.getOrientation()) || this.getOrientation().equalsIgnoreCase("*") || gen.getOrientation().equalsIgnoreCase("*")) 
		{	
		   Alignments g1 =new Alignments(this.getChr(),this.getStart(),this.getEnd());
		   Alignments g2 =new Alignments(gen.getChr(),gen.getStart(),gen.getEnd());
		   return g1.overlaps(g2);
		}
	   else
		   return false;
	}
	
	/**
	 * Create a new gene by merging overlapping exons and retaining non overlapping exons, ignoring and erasing orientation and CDS
	 * @param other The other gene
	 * @return This gene merged with the other gene, erasing orientation and CDS
	 */
	public RefSeqGene mergeExonsAnyOrientation(RefSeqGene other) {
		if(!this.getChr().equals(other.getChr())) {
			throw new IllegalArgumentException("Cannot merge genes on different chromosomes.");
		}
		TreeSet<BasicLightweightAnnotation> exons = new TreeSet<BasicLightweightAnnotation>();
		exons.addAll(this.getExonSet());
		exons.addAll(other.getExonSet());
		Collection<BasicLightweightAnnotation> mergedExons = BasicLightweightAnnotation.mergeAllOverlappers(exons);
		return new RefSeqGene(mergedExons);
	}
	
	/**
	 * Merge exon sets of genes whose genomic spans overlap in any orientation, and leave singleton genes the same
	 * @param genes The gene set
	 * @return New set of merged genes without orientation information
	 */
	public static Collection<RefSeqGene> mergeAllExonsAnyOrientation(TreeSet<RefSeqGene> genes) {
		Collection<RefSeqGene> rtrn = new TreeSet<RefSeqGene>();
		
		// If provided set of regions contains zero or one element, return a copy of the set itself
		if(genes.size() == 0) return rtrn;
		if(genes.size() == 1) {
			rtrn.addAll(genes);
			return rtrn;
		}
		
		Iterator<RefSeqGene> iter = genes.iterator();
		// Begin with the first gene in the sorted set
		RefSeqGene firstGene = iter.next();
		// A changing exon set that expands when the next gene overlaps it
		// Reset when the next gene does not overlap
		Collection<Alignments> growingExonSet = firstGene.getExonSet();
		
		while(iter.hasNext()) {
			
			RefSeqGene growingGene = new RefSeqGene(growingExonSet);
			RefSeqGene nextGene = iter.next();
			
			if(nextGene.overlapsGeneInAnyOrientation(growingGene)) {
				// Next gene span overlaps growing exon set genomic span
				// Merge the exon sets
				Collection<Alignments> tmpExons = new TreeSet<Alignments>();
				tmpExons.addAll(growingExonSet);
				tmpExons.addAll(nextGene.getExonSet());
				growingExonSet.clear();
				growingExonSet.addAll(growingGene.mergeExonsAnyOrientation(nextGene).getExonSet());
				
				if(!iter.hasNext()) {
					// This is the last gene in the set
					// Add the last element and leave loop
					RefSeqGene lastAdd = new RefSeqGene(growingExonSet);
					rtrn.add(lastAdd);
					continue;					
				}
			} else {
				// Next gene does not overlap growing exon set
				// Make gene from latest version of growing exon set and add to the return set
				rtrn.add(growingGene);
				if(!iter.hasNext()) {
					// This is the last gene in the set
					// Add it and leave loop
					rtrn.add(nextGene);
					continue;
				}
				// Reset growing exon set to this new gene
				growingExonSet.clear();
				growingExonSet.addAll(nextGene.getExonSet());
				continue;
			}
		}
		
		return rtrn;

	}
	
	
	
	/*private String makeCigar(){
		String rtrn="";
		Alignments[] exons=this.getExons();
		for(int i=0; i<exons.length; i++){
			Alignments exon=(Alignments)exons[i];
			if(i!=0){
				Alignments previous=(Alignments)exons[i-1];
				Alignments intron=new Alignments(exon.getChr(), previous.getEnd(), exon.getStart());
				if(intron.getSize()>0){rtrn=rtrn+(intron.getSize())+"N";}
			}
			rtrn=rtrn+(exon.getSize())+"M";
		}
		return rtrn;
	}*/
	
	public boolean overlapsGeneInAnyOrientation(RefSeqGene otherGene) {
		   Alignments g1 =new Alignments(this.getChr(),this.getStart(),this.getEnd());
		   Alignments g2 =new Alignments(otherGene.getChr(),otherGene.getStart(),otherGene.getEnd());
		   return g1.overlaps(g2);
	}
	
	public boolean overlaps2(LightweightGenomicAnnotation next){
		if(next==null){return false;}
		if(!this.chr.equals(next.getChromosome())){return false;}
		if(next.getStart()>=start && next.getStart()<=stop){return true;}
		if(next.getEnd()>=start && next.getEnd()<=stop){return true;}
		//else if(track.start<=start &&track.end<=end){return true;}
		return false;
	}
	
	/*private String makeCigar(){
		String rtrn="";
		Alignments[] exons=this.getExons();
		for(int i=0; i<exons.length; i++){
			Alignments exon=(Alignments)exons[i];
			if(i!=0){
				Alignments previous=(Alignments)exons[i-1];
				Alignments intron=new Alignments(exon.getChr(), previous.getEnd(), exon.getStart());
				if(intron.getSize()>0){rtrn=rtrn+(intron.getSize())+"N";}
			}
			rtrn=rtrn+(exon.getSize())+"M";
		}
		return rtrn;
	}*/
	
	public boolean equals(Object o){
		broad.pda.gene.RefSeqGene a=(broad.pda.gene.RefSeqGene)o;
		if(!a.chr.equalsIgnoreCase(this.chr)){return false;}
		if(a.start!=this.start){return false;}
		if(a.stop!=this.stop){return false;}
		if(!a.orientation.equals(this.orientation)){return false;}
		if (a.exonStarts!=null & this.exonStarts!=null ){
			if( a.exonStarts.length!=this.exonStarts.length){return false;}
			
			for(int i=0; i<a.exonStarts.length; i++){
				if(a.exonStarts[i]!=this.exonStarts[i]){return false;}
				if(a.exonEnds[i]!=this.exonEnds[i]){return false;}
			}
		}
		
		return true;
	}
	
	
	/*private String makeCigar(){
		String rtrn="";
		Alignments[] exons=this.getExons();
		for(int i=0; i<exons.length; i++){
			Alignments exon=(Alignments)exons[i];
			if(i!=0){
				Alignments previous=(Alignments)exons[i-1];
				Alignments intron=new Alignments(exon.getChr(), previous.getEnd(), exon.getStart());
				if(intron.getSize()>0){rtrn=rtrn+(intron.getSize())+"N";}
			}
			rtrn=rtrn+(exon.getSize())+"M";
		}
		return rtrn;
	}*/
	
	public boolean almostEqual(RefSeqGene other, int exonFudgeFactor) {
		List<Alignments> otherExons = new ArrayList<Alignments>(other.getExonSet());
		List<Alignments> thisExons = new ArrayList<Alignments>(getExonSet());
		boolean areAlmostEqual = otherExons.size() == thisExons.size();
		
		int i = 0;
		while( areAlmostEqual && i < thisExons.size()) {
			Alignments thisExon = thisExons.get(i);
			Alignments otherExon = otherExons.get(i);
			areAlmostEqual = Math.abs(thisExon.getStart() - otherExon.getStart()) < exonFudgeFactor && Math.abs(thisExon.getEnd() - otherExon.getEnd()) < exonFudgeFactor;
			i++;
		}
		
		return areAlmostEqual;
	}
	
	public boolean almostSameStructure(RefSeqGene other, int intronFudgeFactor) {
		List<Alignments> otherIntrons = new ArrayList<Alignments>(other.getIntronSet());
		List<Alignments> thisIntrons = new ArrayList<Alignments>(getIntronSet());
		boolean almostSameStructure = otherIntrons.size() == thisIntrons.size();
		
		int i = 0;
		while( almostSameStructure && i < thisIntrons.size()) {
			Alignments thisIntron = thisIntrons.get(i);
			Alignments otherIntron = otherIntrons.get(i);
			almostSameStructure = Math.abs(thisIntron.getStart() - otherIntron.getStart()) < intronFudgeFactor && Math.abs(thisIntron.getEnd() - otherIntron.getEnd()) < intronFudgeFactor;
			i++;
		}
		
		return almostSameStructure;
	}
	
	public boolean almostContains(RefSeqGene other, int exonFudgeFactor) {
		return almostContainsRefSeq(other, exonFudgeFactor);
	}
	
	protected boolean almostContainsRefSeq(RefSeqGene other, int exonFudgeFactor) {
		if(!other.getOrientation().equals(getOrientation())) {
			return false;
		}
		List<Alignments> otherExons = new ArrayList<Alignments>(other.getExonSet());
		List<Alignments> thisExons = new ArrayList<Alignments>(getExonSet());
		boolean isAlmostContained = otherExons.size() <= thisExons.size();
		
		if( "-".equals(getOrientation())) {
			int i = otherExons.size() - 1;
			int j = thisExons.size() - 1;
			while( isAlmostContained && i >=0) {
				Alignments thisExon = thisExons.get(j);
				Alignments otherExon = otherExons.get(i);
				isAlmostContained = Math.abs(thisExon.getStart() - otherExon.getStart()) < exonFudgeFactor && Math.abs(thisExon.getEnd() - otherExon.getEnd()) < exonFudgeFactor;
				i--;
				j--;
			}
		} else {
			int i = 0;
			while( isAlmostContained && i < otherExons.size()) {
				Alignments thisExon = thisExons.get(i);
				Alignments otherExon = otherExons.get(i);
				isAlmostContained = Math.abs(thisExon.getStart() - otherExon.getStart()) < exonFudgeFactor && Math.abs(thisExon.getEnd() - otherExon.getEnd()) < exonFudgeFactor;
				i++;
			}
		} 
		
		return isAlmostContained;
	}
	
	public boolean almostContainsStructure(RefSeqGene other, int exonFudgeFactor) {
		return almostContainsStructureRefSeq(other, exonFudgeFactor);
	}
	
	protected boolean almostContainsStructureRefSeq(RefSeqGene other, int exonFudgeFactor) {
		if(!other.getOrientation().equals(getOrientation())) {
			return false;
		}
		List<Alignments> otherIntrons = new ArrayList<Alignments>(other.getIntronSet());
		List<Alignments> thisIntrons = new ArrayList<Alignments>(getIntronSet());
		if(otherIntrons.size() == 0 && thisIntrons.size() == 0) {
			return  (Math.abs(getStart() - other.getStart() ) < exonFudgeFactor) &&  (Math.abs(getEnd() - other.getEnd()) < exonFudgeFactor);
		}
		
		boolean isAlmostContained = otherIntrons.size() <= thisIntrons.size();
		
		
		Iterator<Alignments> otherIntronIt = otherIntrons.iterator();
		Iterator<Alignments> thisIntronIt  = thisIntrons.iterator();
		
		Alignments lastIntron = thisIntronIt.next();
		while(isAlmostContained && otherIntronIt.hasNext()) {
			Alignments otherIntron = otherIntronIt.next();
			boolean found = lastIntron.almostEqual(otherIntron, exonFudgeFactor);
			while(!found && lastIntron.compareTo(otherIntron)< 0 && thisIntronIt.hasNext()  ){
				lastIntron = thisIntronIt.next();
				found = lastIntron.almostEqual(otherIntron, exonFudgeFactor); 
			}
			isAlmostContained = found;
		}
		
		return isAlmostContained;
	}
	
	public int hashCode() {
		Set<Alignments> exons = getExonSet();
		
		String rtrn=this.getChr()+"\t"+this.getStart()+"\t"+this.getEnd()+"\t"+orientation+"\t"+start+"\t"+stop+"\t"+exons.size();
		String sizes="";
		String starts="";
		for(Alignments exon : exons){
			sizes=sizes+(exon.length())+",";
			starts=starts+(exon.getStart()-getStart())+",";
		}
		rtrn=rtrn+"\t"+sizes+"\t"+starts;
		
		return rtrn.hashCode();
	}
	
	/*public String toSAM(){
		String qualityString="255";
		String qualityLetters="";
		String mateReference="*";
		String matePosition="0";
		String insertSize="0";
		
		int strandNum=0;
		if(this.orientation.equalsIgnoreCase("-")){strandNum=16;}
		if(this.sequence==null){this.sequence="";}
		
		for(int i=0; i<sequence.toCharArray().length; i++){qualityLetters=qualityLetters+"I";}
		
		String cigar=makeCigar();
		
		String rtrn=name+"\t"+strandNum+"\t"+chr+"\t"+(start+1)+"\t"+qualityString+"\t"+cigar+"\t"+mateReference+"\t"+matePosition+"\t"+insertSize+"\t"+sequence+"\t"+qualityLetters;
		
		return rtrn;
	}
	
	public String toSAM(RefSeqGene mate){
		String qualityString="255";
		String qualityLetters="";
		String mateReference="=";
		String matePosition=new Integer(mate.getStart()+1).toString();
		String insertSize=new Integer(mate.getStart()-start).toString();
		
		int strandNum=0;
		if(this.orientation.equalsIgnoreCase("-")){strandNum=16;}
		
		for(int i=0; i<sequence.toCharArray().length; i++){qualityLetters=qualityLetters+"I";}
		
		String cigar=makeCigar();
		
		String rtrn=name+"\t"+strandNum+"\t"+chr+"\t"+(start+1)+"\t"+qualityString+"\t"+cigar+"\t"+mateReference+"\t"+matePosition+"\t"+insertSize+"\t"+sequence+"\t"+qualityLetters;
				
		return rtrn;
	}*/
	
	/*private String makeCigar(){
		String rtrn="";
		Alignments[] exons=this.getExons();
		for(int i=0; i<exons.length; i++){
			Alignments exon=(Alignments)exons[i];
			if(i!=0){
				Alignments previous=(Alignments)exons[i-1];
				Alignments intron=new Alignments(exon.getChr(), previous.getEnd(), exon.getStart());
				if(intron.getSize()>0){rtrn=rtrn+(intron.getSize())+"N";}
			}
			rtrn=rtrn+(exon.getSize())+"M";
		}
		return rtrn;
	}*/
	
	//TODO: This should look at the Exon starts and Exon Ends NOT making the exons and then looping through
	//sort by midpoint
	public int compareTo(RefSeqGene b){
		int result = getAlignment().compareTo(b.getAlignment());
		if(result == 0) {
			List<Alignments> aExons =new ArrayList<Alignments> (getExonSet());
			List<Alignments> bExons =new ArrayList<Alignments> (b.getExonSet());
			
			int minLength = Math.min(aExons.size(), bExons.size());
			int idx = 0;
			while(idx < minLength && result == 0) {
				result = aExons.get(idx).compareTo(bExons.get(idx));
				idx++;
			}
			if(result == 0) {
				result = aExons.size() - bExons.size();
			}
			//check for location of the cds
			if(result ==0){
				result=getCDSRegion().compareTo(b.getCDSRegion());
			}
		}

		
		return result;
	}
	
	public RefSeqGene copy(){
		return new RefSeqGene(this.getChr(), this.getStart(), this.getEnd(), this.getName(), this.getOrientation(), this.exonStarts, this.exonEnds, this.blockStart, this.blockEnd);
	}
	
	/*public String toSAM(){
		String qualityString="255";
		String qualityLetters="";
		String mateReference="*";
		String matePosition="0";
		String insertSize="0";
		
		int strandNum=0;
		if(this.orientation.equalsIgnoreCase("-")){strandNum=16;}
		if(this.sequence==null){this.sequence="";}
		
		for(int i=0; i<sequence.toCharArray().length; i++){qualityLetters=qualityLetters+"I";}
		
		String cigar=makeCigar();
		
		String rtrn=name+"\t"+strandNum+"\t"+chr+"\t"+(start+1)+"\t"+qualityString+"\t"+cigar+"\t"+mateReference+"\t"+matePosition+"\t"+insertSize+"\t"+sequence+"\t"+qualityLetters;
		
		return rtrn;
	}
	
	public String toSAM(RefSeqGene mate){
		String qualityString="255";
		String qualityLetters="";
		String mateReference="=";
		String matePosition=new Integer(mate.getStart()+1).toString();
		String insertSize=new Integer(mate.getStart()-start).toString();
		
		int strandNum=0;
		if(this.orientation.equalsIgnoreCase("-")){strandNum=16;}
		
		for(int i=0; i<sequence.toCharArray().length; i++){qualityLetters=qualityLetters+"I";}
		
		String cigar=makeCigar();
		
		String rtrn=name+"\t"+strandNum+"\t"+chr+"\t"+(start+1)+"\t"+qualityString+"\t"+cigar+"\t"+mateReference+"\t"+matePosition+"\t"+insertSize+"\t"+sequence+"\t"+qualityLetters;
				
		return rtrn;
	}*/
	
	/*private String makeCigar(){
		String rtrn="";
		Alignments[] exons=this.getExons();
		for(int i=0; i<exons.length; i++){
			Alignments exon=(Alignments)exons[i];
			if(i!=0){
				Alignments previous=(Alignments)exons[i-1];
				Alignments intron=new Alignments(exon.getChr(), previous.getEnd(), exon.getStart());
				if(intron.getSize()>0){rtrn=rtrn+(intron.getSize())+"N";}
			}
			rtrn=rtrn+(exon.getSize())+"M";
		}
		return rtrn;
	}*/
	
	public RefSeqGene trim(int relativeStart, int relativeEnd){
		if(this.exonStarts.length>1){
		
			Collection<Alignments> gaps=new TreeSet<Alignments>(); //introns in relative space
			for(Alignments intron: getIntronSet()){
				Alignments relativeIntron=new Alignments(chr, intron.getStart()-this.getStart(), intron.getEnd()-this.getStart());
				gaps.add(relativeIntron);
			}
				
			//System.err.println(gaps);
			GenomeWithGaps2 gwg=new GenomeWithGaps2(gaps, getChr(), getAlignment().getSize(), getOrientation());
			
			RefSeqGene relativeWindow = gwg.getRelativeWindow(relativeStart, relativeEnd);
			RefSeqGene relativeGene = null; 
			if(relativeWindow != null) {	//TODO: REMOVE THIS CHECK ONLY TEMPORARY	
				relativeGene= relativeWindow.addConstantFactor(this.getStart());
				relativeGene.setOrientation(orientation);
			}
			
			return relativeGene;
		}
		
		else{
			Alignments exon=new Alignments(this.chr, this.getStart()+(relativeStart), this.getStart()+relativeEnd);
			return new RefSeqGene(exon);
		}
	}
	
	public RefSeqGene findLongestORF() {
		int [] orf = findLongestORF(sequence);
		String cdsOrientation = orientation;
		if("*".equals(orientation)) {
			int [] reverseOrf = findLongestORF(Sequence.reverseSequence(sequence));
			cdsOrientation = "+";
			if(reverseOrf[1] - reverseOrf[0] > orf[1]- orf[0]) {
				orf[0] = reverseOrf[0];
				orf[1] = reverseOrf[1];
				cdsOrientation = "-";
			}
		}
		
		int orfGenomicStart = "-".equals(cdsOrientation) ? mapToGenomic(orf[1], cdsOrientation):  mapToGenomic(orf[0], cdsOrientation)  ; 
		int orfGenomicEnd   = "-".equals(cdsOrientation) ? mapToGenomic(orf[0], cdsOrientation):  mapToGenomic(orf[1], cdsOrientation);
		
		//System.err.println("ORF genomic start end " + orfGenomicStart +"-" + orfGenomicEnd);
		
		LightweightGenomicAnnotation geneCDS = new BasicLightweightAnnotation(getChr(), orfGenomicStart, orfGenomicEnd);
		Set<LightweightGenomicAnnotation> cdsExons = new TreeSet<LightweightGenomicAnnotation>();
		Set<Alignments> exonSet = getExonSet();
		for(Alignments e : exonSet) {
			if(geneCDS.overlaps(e)) {
				cdsExons.add(e.intersect(geneCDS));
			}
		}
		
		RefSeqGene cds = null;
		if(cdsExons.size() > 0 ) {
			cds = new RefSeqGene(cdsExons);
			cds.setOrientation(cdsOrientation);
			cds.setName(getName());
			cds.setCountScore((orf[1]- orf[0])/(double)(sequence.length()));
		}
		return cds;
	}
	
	private Alignments trimLast(Alignments lastExon, Alignments cds) {
		Alignments rtrn=new Alignments(cds.getChr(), lastExon.getStart(), Math.min(cds.getEnd(), lastExon.getEnd()));
		return rtrn;
	}
	
	public int numOverlappingExons(RefSeqGene iso) {
		int res=0;
		for (Alignments myEx: this.getExonSet()){
			for (Alignments isoEx: iso.getExonSet()){
				if (myEx.overlaps(isoEx)){
					res++; break;
				}
			}
		}
		return res;
	}
	
	public int numOfCompatibleIntrons(RefSeqGene g) {
		int rtrn=0;
		Alignments[] myIntrons= this.getIntronsBlocks();
		Alignments[] gIntrons= g.getIntronsBlocks();
		int myNan=-10;
		int match=0;
		int[]ix=new int[gIntrons.length];
		for(int i=0; i<gIntrons.length; i++){
			ix[i]=myNan;
			for(int j=0; j<myIntrons.length; j++){
				if (gIntrons[i].equals(myIntrons[j])){
					ix[i]=j;
					match++;
				}
			}
		}
		if (match==0)
			return 0;
		int currLength=0;
		double Prev=-1;
		for (int i=0;i<ix.length;i++){
			if(ix[i]-Prev==1 | (Prev==-1 & ix[i]!= myNan)){ //continue chain, or start chain 
				Prev=ix[i];
				currLength++;
			}
			else
			{
				if (currLength > rtrn)
					rtrn=currLength;
				currLength=0;
				Prev=-1;
			}
		}
		if (currLength > rtrn)
			rtrn=currLength;
		return rtrn;
	}
	
	//@return A RefSeqGene that represent the trimmed piece if the region overlaps any exonic regions
	/*public RefSeqGene trimAbsolute(int alignmentStart, int alignmentEnd) {
		Collection<Alignments> rtrn=new TreeSet<Alignments>();
		//get exons overlapping from start to end
		Alignments region=new Alignments(chr, alignmentStart, alignmentEnd);
		
		//go through each exon
		Collection<Alignments> exons=this.getSortedAndUniqueExons();
		for(Alignments exon: exons){
			if(exon.overlaps(region)){rtrn.add(exon);}
		}
	
		if(rtrn.isEmpty()){return null;}
	
		//now trim by ends
		//Get first and last exons
		Alignments first=(rtrn.iterator().next());
		Alignments last=(Alignments)rtrn.toArray()[rtrn.size()-1];
				
		//TODO Need to check that alignmentStart and alignmentEnd are within exonic regions or else return null as well
		if(!(alignmentStart>=first.getStart() && alignmentStart<=first.getEnd())){return null;}
		if(! (alignmentEnd>=last.getStart()&& alignmentEnd<=last.getEnd())) {return null;} 
	
		//remove them from the current collection
		rtrn.remove(first);
		rtrn.remove(last);
		
		if(first.equals(last)){
			//trim first starting at alignmentStart
			Alignments newFirst=new Alignments(first.getChr(), alignmentStart, alignmentEnd);
			rtrn.add(newFirst);
		}
		else{
			//trim first starting at alignmentStart
			Alignments newFirst=new Alignments(first.getChr(), alignmentStart, first.getEnd());
			//trim last ending alignmentEnd
			Alignments newEnd=new Alignments(last.getChr(), last.getStart(), alignmentEnd);
			
			//add new first and last to the collection 
			rtrn.add(newFirst);
			rtrn.add(newEnd);
		}
		
		//make and return RefSeqGene
		return new RefSeqGene(rtrn);
	}*/
	
	//@return A RefSeqGene that represent the trimmed piece if the region overlaps any exonic regions
	//MG: Replaced
	public RefSeqGene trimAbsolute(int alignmentStart, int alignmentEnd) {
		Collection<Alignments> rtrn=new TreeSet<Alignments>();
		//get exons overlapping from start to end
		Alignments region=new Alignments(chr, alignmentStart, alignmentEnd);
		
		//go through each exon
		Collection<Alignments> exons=this.getSortedAndUniqueExons();
		for(Alignments exon: exons){
			if(exon.overlaps(region)){
				//System.out.println(exon);
				rtrn.add(exon);
			}
		}

		if(rtrn.isEmpty()){return null;}

		//now trim by ends
		//Get first and last exons
		Alignments first=(rtrn.iterator().next());
		Alignments last=(Alignments)rtrn.toArray()[rtrn.size()-1];
				
		//TODO Need to check that alignmentStart and alignmentEnd are within exonic regions or else return null as well
		/*if(!(alignmentStart>=first.getStart() && alignmentStart<=first.getEnd())){
			System.err.println("NULL");
			return null;
		}
		if(! (alignmentEnd>=last.getStart()&& alignmentEnd<=last.getEnd())) {
			System.err.println("NULL");
			return null;
		} */

		//remove them from the current collection
		rtrn.remove(first);
		rtrn.remove(last);
		
		if(first.equals(last)){
			//trim first starting at alignmentStart
			Alignments newFirst=new Alignments(first.getChr(), alignmentStart, alignmentEnd);
			rtrn.add(newFirst);
		}
		else{
			//trim first starting at alignmentStart
			Alignments newFirst=new Alignments(first.getChr(), alignmentStart, first.getEnd());
			//trim last ending alignmentEnd
			Alignments newEnd=new Alignments(last.getChr(), last.getStart(), alignmentEnd);
			
			//add new first and last to the collection 
			rtrn.add(newFirst);
			rtrn.add(newEnd);
		}
		
		//make and return RefSeqGene
		//System.err.println(rtrn);
		RefSeqGene rtrnGene = new RefSeqGene(rtrn, getName(), getOrientation());
		rtrnGene.setCDS(this.getCDSRegion());
		return rtrnGene;
	}
	
	private Alignments trimFirst(Alignments firstExon, Alignments cds) {
		Alignments rtrn=new Alignments(cds.getChr(), Math.max(firstExon.getStart(), cds.getStart()), Math.min(firstExon.getEnd(), cds.getEnd()));
		return rtrn;
	}
	
	//computes the percentage of this that overlaps with gene2
	public double percentOverlapping(RefSeqGene gene2){
		if(!this.chr.equalsIgnoreCase(gene2.getChr())){return 0.0;}
		IntervalTree<Alignments> test=makeExonTree();
		
		Collection<Alignments> exons=gene2.getExonSet();
		int totalOverlap=0;
		for(Alignments exon: exons){
			Iterator<Node<Alignments>> overlappers=test.overlappers(exon.getStart(), exon.getEnd());
			while(overlappers.hasNext()){
				Alignments align=overlappers.next().getValue();
				int overlap=getOverlap(exon, align);
				totalOverlap+=overlap;
			}
		}
		return (double)totalOverlap/getTranscriptLength();
	}
	
	//computes the percentage of the genomic region spanned by this object that overlaps with that of gene2
	public Double percentGenomeOverlapping(RefSeqGene gene2) {
		
		if(!this.chr.equalsIgnoreCase(gene2.getChr())){return 0.0;}
		Alignments g1 = new Alignments(this.getChr(),this.getStart(),this.getEnd());
		Alignments g2 = new Alignments(gene2.getChr(),gene2.getStart(),gene2.getEnd());
		int overlap=getOverlap(g1, g2);
		return (double)overlap/this.getGenomicLength();
		
	}
	
	private IntervalTree<Alignments> makeExonTree() {
		IntervalTree<Alignments> tree=new IntervalTree<Alignments>();
		
		Collection<Alignments> exons= getExonSet();
		for(Alignments exon: exons){
			tree.put(exon.getStart(), exon.getEnd(), exon);
		}
		
		return tree;
	}
	
	//Fixed off-by-one error
	private String makeCigar(){
		String rtrn="";
		Alignments[] exons=this.getExons();
		for(int i=0; i<exons.length; i++){
			Alignments exon=(Alignments)exons[i];
			if(i!=0){
				Alignments previous=(Alignments)exons[i-1];
				Alignments intron=new Alignments(exon.getChr(), previous.getEnd(), exon.getStart());
				if(intron.getSize()>0){rtrn=rtrn+(intron.getSize())+"N";}
			}
			rtrn=rtrn+(exon.getSize())+"M";
		}
		return rtrn;
		
	}
	
	/*public String toSAM(){
		String qualityString="255";
		String qualityLetters="";
		String mateReference="*";
		String matePosition="0";
		String insertSize="0";
		
		int strandNum=0;
		if(this.orientation.equalsIgnoreCase("-")){strandNum=16;}
		if(this.sequence==null){this.sequence="";}
		
		for(int i=0; i<sequence.toCharArray().length; i++){qualityLetters=qualityLetters+"I";}
		
		String cigar=makeCigar();
		
		String rtrn=name+"\t"+strandNum+"\t"+chr+"\t"+(start+1)+"\t"+qualityString+"\t"+cigar+"\t"+mateReference+"\t"+matePosition+"\t"+insertSize+"\t"+sequence+"\t"+qualityLetters;
		
		return rtrn;
	}
	
	public String toSAM(RefSeqGene mate){
		String qualityString="255";
		String qualityLetters="";
		String mateReference="=";
		String matePosition=new Integer(mate.getStart()+1).toString();
		String insertSize=new Integer(mate.getStart()-start).toString();
		
		int strandNum=0;
		if(this.orientation.equalsIgnoreCase("-")){strandNum=16;}
		
		for(int i=0; i<sequence.toCharArray().length; i++){qualityLetters=qualityLetters+"I";}
		
		String cigar=makeCigar();
		
		String rtrn=name+"\t"+strandNum+"\t"+chr+"\t"+(start+1)+"\t"+qualityString+"\t"+cigar+"\t"+mateReference+"\t"+matePosition+"\t"+insertSize+"\t"+sequence+"\t"+qualityLetters;
				
		return rtrn;
	}*/
	
	/*private String makeCigar(){
		String rtrn="";
		Alignments[] exons=this.getExons();
		for(int i=0; i<exons.length; i++){
			Alignments exon=(Alignments)exons[i];
			if(i!=0){
				Alignments previous=(Alignments)exons[i-1];
				Alignments intron=new Alignments(exon.getChr(), previous.getEnd(), exon.getStart());
				if(intron.getSize()>0){rtrn=rtrn+(intron.getSize())+"N";}
			}
			rtrn=rtrn+(exon.getSize())+"M";
		}
		return rtrn;
	}*/
	
	private RefSeqGene addConstantFactor(int factor){
		Collection<LightweightGenomicAnnotation> set=new TreeSet();
		for(Alignments exon: this.getExonSet()){
			Alignments abs=new Alignments(exon.getChr(), exon.getStart()+factor, exon.getEnd()+factor);
			set.add(abs);
		}
		return new RefSeqGene(set);
	}
	
	/**
	 * Maps cDNA coordinates to genomic coordinates
	 * @param start within transcript
	 * @param end end in transcript
	 * @param cdsOrientation
	 * @return A RefSeqGene that that maps to genome space
	 */
	private int mapToGenomic(int position, String cdsOrientation) {
		int coveredGround = 0;
		int exonIdx = -1;
		if("-".equals(cdsOrientation)) {
			position  = sequence.length() - position;
		}
		
		int genomicPosition = getStart();
		if(position > 0) {
		
			while(coveredGround < position) {
				exonIdx++;
				coveredGround += exonEnds[exonIdx] - exonStarts[exonIdx];
			}
	
			coveredGround = coveredGround - (exonEnds[exonIdx] - exonStarts[exonIdx]); //Substract the length of the exon that contains start
			
			genomicPosition = /*getStart() +*/ exonStarts[exonIdx] + (position - coveredGround ); //Genomic Start is the start of the exon + whatever extra sequence.
		}		
		
		return genomicPosition;
	}
	
	public int transcriptToGenomicPosition(int transcriptPosition) {
		int relativePosition= transcriptPosition;
		int genomicPosition = -1;
		if("-".equals(getOrientation())){
			relativePosition=(getSize())-(transcriptPosition) - 1; //Recall that the transcript is [0,L) closed open from 0 to the length L
		}
		RefSeqGene samRecord= trim(relativePosition, relativePosition+1 ); //TODO: Should we do relativePosition-1, relativePosition when transcript in reversed orientation?
		if(samRecord != null) {
			genomicPosition = samRecord.getStart() ;
		} else {
			Exception t = new Exception ("n");
			logger.warn("Could not map back to genome position relative position was: "+ relativePosition+ " this transicript: " +toBED() + " stack trace: ", t );
		}
		return genomicPosition;
	}
	
	/**
	 * Convert an interval in transcriptome space to genome space
	 * @param startPosOnTranscript Start position in transcriptome space
	 * @param endPosOnTranscript End position in transcriptome space
	 * @return The interval (possibly spliced) in genome space
	 */
	public RefSeqGene transcriptToGenomicPosition(int startPosOnTranscript, int endPosOnTranscript) {
		if(startPosOnTranscript > endPosOnTranscript) {
			throw new IllegalArgumentException("Start position in transcriptome space cannot be greater than end position.");
		}
		int genomicStart = isNegativeStrand() ? transcriptToGenomicPosition(endPosOnTranscript) : transcriptToGenomicPosition(startPosOnTranscript);
		int genomicEnd = isNegativeStrand() ? transcriptToGenomicPosition(startPosOnTranscript) : transcriptToGenomicPosition(endPosOnTranscript);
		RefSeqGene genomicInterval = new RefSeqGene(getChr(), genomicStart, genomicEnd);
		return getOverlap(genomicInterval);
	}
	
	/**
	 * Maps genomic coordinates to cDNA coordinates
	 * @param genomicPosition within transcript
	 * @return 0-based position within oriented start of transcript
	 */
	public int genomicToTranscriptPosition(int genomicPosition) {
		if(genomicPosition < getStart() || genomicPosition > getEnd() ) {
			return -1;
		}
		
		List<Alignments> exons = new ArrayList<Alignments>(getExonSet());
				
		int position = 0;
		if("-".equals(getOrientation())) {
			for(int i = exons.size() -1 ; i >=0; i--) {
				Alignments e = exons.get(i);
				if(genomicPosition < e.getStart()) {
					position += e.length();
				} else if( e.getStart() <= genomicPosition && genomicPosition < e.getEnd()) {
					position += e.getEnd() - 1 - genomicPosition; //Recall that ends are open, so the first position (0) in the exon when going backwards is end -1
					break;
				} else {
					return -1;
				}
				/*if(e.getStart() >= genomicPosition) {
					position += e.length();
					if(e.getStart() == genomicPosition) {
						break;
					}
				} else {
					if(e.getEnd() <= genomicPosition){//genomicPosition landed in an intron, return -1
						return -1;
					} else {
						position += e.getEnd() - genomicPosition;
					}
					break;
				}*/
			}
		} else {
			for(int i = 0; i < exons.size() ; i++) {
				Alignments e = exons.get(i);
				
				if(genomicPosition > e.getEnd()) {
					position += e.length();
				} else if (e.getStart() <= genomicPosition && genomicPosition < e.getEnd()) {
					position +=  genomicPosition - e.getStart();
					break;
				} else {
					return -1;
				}
				//System.err.println("e: " + e.toUCSC());
				/*if(e.getEnd() < genomicPosition) {
					position += e.length();
					//System.err.println("\tadded to position");
					if(e.getEnd() == genomicPosition-1) {
						break;
					}
				} else {
					//System.err.println("\tfinishing..");
					if(e.getStart() > genomicPosition){//genomicPosition landed in an intron, return -1
						//System.err.println("\treturning -1");
						return -1;
					} else {
						position +=  genomicPosition - e.getStart();
						//System.err.println("\tadding to position ");
					}
					break;
				}*/
			}
		}
	
		
		return position;
	}
	
	public static int [] findLongestORF(String sequence) {
		int lastStart = 0;
		int lastEnd = 0;
		Matcher m = START_CODON_PATRN.matcher(sequence);
		while(m.find()) {
			int startCodonPos = m.start(); 
			int thisORFEnd = startCodonPos;
			boolean foundStopCodon = false;
			//System.err.println("new orf start: " + startCodonPos);
			while(thisORFEnd < sequence.length() - 3  && !foundStopCodon) {
				String currentCodon = sequence.substring(thisORFEnd, thisORFEnd+3);
				//System.err.print(" "+currentCodon+ " ");
				thisORFEnd += 3;
				foundStopCodon = isStopCodon(currentCodon);
			}
			//System.err.println("Orf length: " + (thisORFEnd - startCodonPos) );
			if(lastEnd - lastStart < thisORFEnd - startCodonPos) {
				lastStart = startCodonPos;
				lastEnd   = thisORFEnd;
				//System.err.println("It was a winner");
			}
		}
		int [] startEnd = {lastStart, lastEnd};
		return startEnd;
	}
	
	/***************************************************************************/
	//Collapse methods
	private Set collapse(Collection alignments){
		Set rtrn=new TreeSet();
		boolean overlappingB=true;
		boolean finished=false;
		
		while(!finished){
		Map<Alignments, Set> overlapping=getOverlapping(alignments);
		boolean done=isDone(overlapping);
		if(!done){overlappingB=false;}
		Set<Alignments> collapsed=collapse(overlapping);
		rtrn.addAll(collapsed);
		if(overlappingB){finished=true;}
		}
		
		return rtrn;
	}
	
	private Set<Alignments> collapse(Map<Alignments, Set> overlappingMap){
		Set rtrn=new TreeSet();
		
		for(Alignments align: overlappingMap.keySet()){
			Set<Alignments> alignments=overlappingMap.get(align);
			int start=align.getStart();
			int end=align.getEnd();
			for(Alignments peak: alignments){
				start=Math.min(start, peak.getStart());
				end=Math.max(end, peak.getEnd());
			}
			Alignments collapsed=new Alignments(align.getChr(), start, end);
			rtrn.add(collapsed);
		}
		return rtrn;
	}
	
	public RefSeqGene takeUnion(RefSeqGene other) {
		TreeSet<LightweightGenomicAnnotation> combinedExons = new TreeSet<LightweightGenomicAnnotation>();
		combinedExons.addAll(getExonSet());
		combinedExons.addAll(other.getExonSet());
		List<LightweightGenomicAnnotation> mergedExons = BasicLightweightAnnotation.stitchList(combinedExons, 0);
		RefSeqGene union = new RefSeqGene(mergedExons);
		union.setName(getName());
		union.setOrientation(orientation);
		return union;
		
	}

	public double getPercentCDS() {
		return (double)this.getCDS().getSize()/(double)this.getSize();
	}

	public double getPercent3UTR() {
		if(this.get3UTRGene()!=null){
			return (double)this.get3UTRGene().getSize()/(double)this.getSize();
		}
		return 0.0;
	}
	
	
	
	/*public String toSAM(){
		String qualityString="255";
		String qualityLetters="";
		String mateReference="*";
		String matePosition="0";
		String insertSize="0";
		
		int strandNum=0;
		if(this.orientation.equalsIgnoreCase("-")){strandNum=16;}
		if(this.sequence==null){this.sequence="";}
		
		for(int i=0; i<sequence.toCharArray().length; i++){qualityLetters=qualityLetters+"I";}
		
		String cigar=makeCigar();
		
		String rtrn=name+"\t"+strandNum+"\t"+chr+"\t"+(start+1)+"\t"+qualityString+"\t"+cigar+"\t"+mateReference+"\t"+matePosition+"\t"+insertSize+"\t"+sequence+"\t"+qualityLetters;
		
		return rtrn;
	}
	
	public String toSAM(RefSeqGene mate){
		String qualityString="255";
		String qualityLetters="";
		String mateReference="=";
		String matePosition=new Integer(mate.getStart()+1).toString();
		String insertSize=new Integer(mate.getStart()-start).toString();
		
		int strandNum=0;
		if(this.orientation.equalsIgnoreCase("-")){strandNum=16;}
		
		for(int i=0; i<sequence.toCharArray().length; i++){qualityLetters=qualityLetters+"I";}
		
		String cigar=makeCigar();
		
		String rtrn=name+"\t"+strandNum+"\t"+chr+"\t"+(start+1)+"\t"+qualityString+"\t"+cigar+"\t"+mateReference+"\t"+matePosition+"\t"+insertSize+"\t"+sequence+"\t"+qualityLetters;
				
		return rtrn;
	}*/
	
	/*private String makeCigar(){
		String rtrn="";
		Alignments[] exons=this.getExons();
		for(int i=0; i<exons.length; i++){
			Alignments exon=(Alignments)exons[i];
			if(i!=0){
				Alignments previous=(Alignments)exons[i-1];
				Alignments intron=new Alignments(exon.getChr(), previous.getEnd(), exon.getStart());
				if(intron.getSize()>0){rtrn=rtrn+(intron.getSize())+"N";}
			}
			rtrn=rtrn+(exon.getSize())+"M";
		}
		return rtrn;
	}*/
	
	
    
	
	
	//Moran TODO: debug . Not working
	/*public boolean addExon(RefSeqGene exon) {
		
		boolean res=false;
		
		
		if (this.getOrientation().equals(exon.getOrientation()) || exon.getOrientation().equals("*") ){
			
			//if this is a first exon added and this is a partial bed- make it a complete first:
			if (this.getNumExons()==0){
			this.blockStart=this.start;
			this.blockEnd=this.stop;
			this.exonStarts[0]=this.start;
			this.exonEnds[0]=this.stop;
			this.exonScores[0]=0;
			}
			
			this.start=Math.min(this.start, exon.getStart());
			this.stop= Math.max(this.stop,exon.getEnd());
			int i= this.getNumExons();
			this.blockStart=this.start;
			this.blockEnd=this.stop;
			this.exonStarts[i]=exon.start;
			this.exonEnds[i]=exon.stop;
			this.exonScores[i]=0;
			
			res=true; //added the exon
		}
		return res;
	
	}
	*/
	public RefSeqGene extendAnnotation(int length) { 
		RefSeqGene subAnnotation = new RefSeqGene(this.copy());

		if("-".equals(orientation)){
			for(int i = 0 ; i < subAnnotation.exonEnds.length; i++){
				if(subAnnotation.exonStarts[i] ==  subAnnotation.getStart()) {
					subAnnotation.exonStarts[i] = Math.max(subAnnotation.getStart()-length,0);
				}
			}
			subAnnotation.setStart(Math.max(subAnnotation.getStart()-length,0));

		}
		else{
			for(int i = 0 ; i < subAnnotation.exonEnds.length; i++){
				if(subAnnotation.exonEnds[i] ==  subAnnotation.getEnd()) {
					subAnnotation.exonEnds[i] = subAnnotation.getEnd()+length;
				}
				else{
				//	System.out.println("Annotation End: "+subAnnotation.getEnd()+" Exon end: "+subAnnotation.exonEnds[i]);
				}
			}
			subAnnotation.setEnd(subAnnotation.getEnd()+length);
			//System.out.println(this.getEnd()+" to "+subAnnotation.getEnd());

		}
		return subAnnotation;
	}
}
