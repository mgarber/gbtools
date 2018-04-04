package broad.pda.snp;

import broad.pda.snp.DBSNPReader.DBSNP;

public class SNPSampleInfo {
	String sampleId;
	boolean sangerCovered;
	boolean four54Covered;
	char[] sangerGenotype = null;
	char[] four54Genotype = null;
	int numOf454ForwardReads;
	int numOf454ReverseReads;
	float four54Score;
	BroadSNP snp;
	
	

	public SNPSampleInfo(String sampleId, BroadSNP snp) {
		super();
		this.snp = snp;
		this.sampleId = sampleId;
	}

	public boolean is454WildtypeHomozygous() {
		
		return has454Genotype() && (snp.getReferenceAlele() == four54Genotype[0] && snp.getReferenceAlele() == four54Genotype[1]); 
	}
	
	public boolean has454Genotype() { return four54Genotype != null; }

	public boolean isSangerWildtypeHomozygous() {
		return hasSangerGenotype() && (snp.getReferenceAlele() == sangerGenotype[0] && snp.getReferenceAlele() == sangerGenotype[1]); 	
	}
	
	public boolean hasSangerGenotype() { return sangerGenotype != null; }
	
	public float get454Call() {
		float result = -1;
		if(has454Genotype()) {
			result = getFour54Score();
		} else if(!isFour54Covered()) {
			result = 0;
		}
		return result;
	}

	public int getSangerCallAsInt() {
		int result = 0;
		if(hasSangerGenotype()) {
			result = isSangerWildtypeHomozygous() ? -1 : 1;
		}
		return result;
	}

	public String getSampleId() {
		return sampleId;
	}


	public DBSNP getSnp() {
		return snp;
	}


	public boolean isFour54Covered() {
		return four54Covered;
	}


	public char[] getFour54Genotype() {
		return four54Genotype;
	}


	public float getFour54Score() {
		return four54Score;
	}


	public boolean isSangerCovered() {
		return sangerCovered;
	}


	public char[] getSangerGenotype() {
		return sangerGenotype;
	}


	public void setFour54Covered(boolean four54Covered) {
		this.four54Covered = four54Covered;
	}


	public void setFour54Genotype(char[] four54Genotype) {
		this.four54Genotype = four54Genotype;
	}


	public void setFour54Score(float score) {
		this.four54Score = score;
	}


	public void setSampleId(String sampleId) {
		this.sampleId = sampleId;
	}


	public void setSangerCovered(boolean sangerCovered) {
		this.sangerCovered = sangerCovered;
	}


	public void setSangerGenotype(char[] sangerGenotype) {
		this.sangerGenotype = sangerGenotype;
	}


	public void setSnp(BroadSNP snp) {
		this.snp = snp;
	}


	public int getNumOf454ForwardReads() {
		return numOf454ForwardReads;
	}


	public int getNumOf454ReverseReads() {
		return numOf454ReverseReads;
	}


	public void setNumOf454ForwardReads(int numOf454ForwardReads) {
		this.numOf454ForwardReads = numOf454ForwardReads;
	}


	public void setNumOf454ReverseReads(int numOf454ReverseReads) {
		this.numOf454ReverseReads = numOf454ReverseReads;
	}



}
