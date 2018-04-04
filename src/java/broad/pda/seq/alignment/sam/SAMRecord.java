package broad.pda.seq.alignment.sam;

import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;

//TODO: I need to implement this!
public class SAMRecord implements Comparable<SAMRecord>{

	private static final int READ_IS_PAIRED_END_SEQ = 0x0001;
	private static final int READ_IS_PAIRED_MAPPED_PROPER_PAIR = 0x0002;
	private static final int READ_UNMAPPED = 0x0004;
	private static final int UNMAPPED_MATE = 0x0008;
    private static final int READ_STRAND_FLAG = 0x10;
    private static final int MATE_STRAND_FLAG = 0x20;	
	private static final int FIRST_OF_PAIR_FLAG = 0x40;
    private static final int SECOND_OF_PAIR_FLAG = 0x80;

    
	
	int flag;
	int mappingQuality;
	String mateReferenceName="*";
	int matePosition;
	int inferredInsertSize;
	String seq="*";
	String quality="*";
	RefSeqGene gene;
	SAMTags tags;
	String samLine;
	boolean hasNoAlignment=true;
	String name;
	String cigar;
	
	
	public SAMRecord(String samLine){
		String[] tokens=samLine.split("\t");
		this.samLine=samLine;
		this.name=tokens[0];
		this.flag=new Integer(tokens[1]);
		this.seq=tokens[9];
		if(tokens.length>10){this.quality=tokens[10];}
		setTags(tokens);
		if(flag!=4){
			hasNoAlignment=false;
			this.cigar=tokens[5];
			this.matePosition=new Integer(tokens[7]);
			this.inferredInsertSize=new Integer(tokens[8]);
			this.mappingQuality=new Integer(tokens[4]);
			gene=SAMUtils.SAMFormatFullBED(samLine);
			String[] tags=samLine.split("\t"); //starting at 11
		}
	}
	
	
	public SAMRecord(String name, RefSeqGene alignment, String seq) {
		this.name=name;
		this.gene=alignment;
		this.seq=seq;
		cigar=makeCigar();
		this.tags=new SAMTags();
		flag = 0;
	}

	


	//Fixed off-by-one error
	private String makeCigar(){
		String rtrn="";
		Alignments[] exons=gene.getExons();
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
	

	public String toString(){
		String rtrn=name+"\t"+flag+"\t"+gene.getChr()+"\t"+(gene.getStart()+1)+"\t"+mappingQuality+"\t"+cigar+"\t"+mateReferenceName+"\t"+matePosition+"\t"+inferredInsertSize+"\t"+seq+"\t"+quality+"\t"+tags.toString();
		return rtrn;
	}
	

	
	private void setTags(String[] tokens) {
		this.tags=new SAMTags();
		//tags start at [11]
		for(int i=11; i<tokens.length; i++){
			this.tags.setTags(tokens[i]);
		}
	}

	public String getSAMLine(){return this.samLine;}

	public boolean equals(Object other){
		if(other==null){return false;}
		SAMRecord s=(SAMRecord)other;
		if(!name.equalsIgnoreCase(s.name)){return false;}
		if(!seq.equalsIgnoreCase(s.seq)){return false;}
		if(!gene.equals(s.gene)){return false;}
		return true;
	}
	

	public String getName() {
		return this.name;
	}

	public RefSeqGene getGene() {
		return this.gene;
	}

	public int compareTo(SAMRecord other) {
		//first sort by position
		int coor=gene.compareTo(other.gene);
		if(coor==0){
			//then sort by name
			coor=name.compareTo(other.name);
			if(coor==0){
				//then sort by sequence
				coor=seq.compareTo(other.seq);
			}
		}
		return coor;
	}
	
	public int getNumMismatches() {
		return tags.getNumberOfMismatches();
	}

	public boolean overlaps(SAMRecord record2) {
		return this.gene.getAlignment().overlaps(record2.gene.getAlignment());
	}
		
	public void setWeight(double weight){this.tags.setWeight(weight);}


	public void setMappingQuality(int i) {
		this.mappingQuality=i;
	}


	public String getProgramName() {return tags.getProgramName();}


	public void setName(String string) {
		this.name=string;
	}


	public void setProgramName(String name) {
		tags.setProgram(name);
	}


	public void setNumMismatches(int mismatches) {
		tags.setMismatches(mismatches);
	}


	public int getLength() {
		return this.gene.getTranscriptLength();
	}


	public double getScore() {
		double score=(getLength()-getNumMismatches())-getNumMismatches();
		return score;
	}


	public String getSequence() {
		return this.seq;
	}


	public int getFlag() {
		return this.flag;
	}


	public Integer getPosition() {
		return this.gene.getStart();
	}


	public Object getMatePosition() {
		return this.matePosition;
	}
	
	public boolean getSecondOfPairFlag() {
		return (this.getFlag() & SECOND_OF_PAIR_FLAG) != 0;
	}
	
	public void setSecondOfPairFlag() {
		flag = flag | SECOND_OF_PAIR_FLAG;
	}
	
	public boolean getFirstOfPairFlag() {
		return (this.getFlag() & FIRST_OF_PAIR_FLAG) != 0;
	}
	
	public void setFirstOfPairFlag() {
		flag = flag | FIRST_OF_PAIR_FLAG;
	}
	
	public boolean getPairedSequenceFlag()  {
		return (flag & READ_IS_PAIRED_END_SEQ) != 0;
	}
	
	public void setPairedSequenceFlag()  {
		flag =  flag | READ_IS_PAIRED_END_SEQ;
	}
	
	public boolean getProperlyPairedMappedFlag() {
		return  (flag & READ_IS_PAIRED_MAPPED_PROPER_PAIR) != 0;
	}
	
	public void setProperlyPairedMappedFlag() {
		flag = flag | READ_IS_PAIRED_MAPPED_PROPER_PAIR;
	}
	
	public boolean getReadUnmappedFlag() {
		return (READ_UNMAPPED & flag) != 0;
	}
	
	public void setReadUnmappedFlag() {
		flag = READ_UNMAPPED | flag;
	}
	
	public boolean getMateUnmappedFlag() {
		return (UNMAPPED_MATE & flag) != 0;
	}
	
	public void setMateUnmappedFlag() {
		flag = UNMAPPED_MATE | flag;
	}
	
    public boolean getMateNegativeStrandFlag() {
        return (flag & MATE_STRAND_FLAG) != 0;
    }
    
    public void setMateNegativeStrandFlag() {
        flag = flag | MATE_STRAND_FLAG;
    }
	
	/**
     * strand of the query (false for forward; true for reverse strand).
     */
    public boolean getReadNegativeStrandFlag() {
        return (this.getFlag() & READ_STRAND_FLAG) != 0;
    }
    
    public void setReadNegativeStrandFlag() {
        flag = flag | READ_STRAND_FLAG;
    }

    public static  int getReadStrandFlag (){
    	return READ_STRAND_FLAG;
    }
    
    public void reverseOrientation() {
    	flag = flag ^ READ_STRAND_FLAG;
    }

	
}
