package broad.core.annotation;

import java.util.Collection;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.alignment.sam.SAMRecord;

public class PSL {
	
	private String name;
	private String [] rawData;
	private String orientation;
	private String chr;
	private String [] blockSizes;
	private String [] blockStarts;
	
	private static final Pattern goodLine = Pattern.compile("^[0-9]");

	public  PSL(String pslLine) {
		rawData = pslLine.split("\t");
		orientation = rawData[8];
		this.name=rawData[9];
		blockSizes=rawData[18].split(",");
		blockStarts = rawData[20].split(",");
		chr = rawData[13];
	}

	public String getName() {
		return this.name;
	}

	public String toPSL() {
		StringBuilder sb = new StringBuilder(rawData[0]);
		for(int i = 1; i < rawData.length; i++) {
			sb.append("\t").append(rawData[i]);
		}
		return sb.toString();
	}

	public Integer getSeqNumber() {
		return new Integer(this.name.replaceAll("seq.", "").replaceAll("a", "").trim());
	}

	public RefSeqGene toGene() {
		
		if(rawData.length<13){System.err.println(toPSL());}
		
		Collection<Alignments> exons=new TreeSet<Alignments>();		
		
		for(int i=0; i<blockStarts.length; i++){
			int start=new Integer(blockStarts[i]);
			int size=new Integer(blockSizes[i]);
			Alignments align=new Alignments(chr, start, start+size);
			exons.add(align);
		}
		
		RefSeqGene gene=new RefSeqGene(exons);
		gene.setName(this.getName());
		return gene;
	}

	public int getMinCoverage() {
		int min=Integer.MAX_VALUE;
		RefSeqGene gene=this.toGene();
		Collection<Alignments> exons=gene.getSortedAndUniqueExons();
		for(Alignments exon: exons){min=Math.min(min, exon.length());}
		if(min==Integer.MAX_VALUE){min=0;}
		return min;
	}

	public double getPercentMapped() {
		int readSize=new Integer(rawData[10]);
		int mappedSize=new Integer(rawData[12])-new Integer(rawData[11]);
		
		return (double)mappedSize/readSize;
	}

	public String getSequence() {
		if(rawData.length>21){return rawData[21].replaceAll(",", "").toUpperCase();}
		return "*";
	}

	public int getMismatches() {
		return new Integer(rawData[1]);
	}

	public SAMRecord toSAM() {
		String name=getName().replaceAll("@", "");
		SAMRecord record=new SAMRecord(name, toGene(), getSequence());
		record.setNumMismatches(getMismatches());
		if("-".equals(orientation)) {
			record.setReadNegativeStrandFlag();
		}
		return record;
	}
	
	public static PSL create(String rawData) {
    	PSL psl = null;
    	if(rawData!=null) {
    		rawData = rawData.trim();
        	Matcher m = goodLine.matcher(rawData);
        	if(m.find()) {
        		psl=new PSL(rawData);
        	}
    	}

    	return psl;
	}

}
