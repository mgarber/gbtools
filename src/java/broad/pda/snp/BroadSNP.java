package broad.pda.snp;

import java.util.Collection;
import java.util.HashMap;

import broad.pda.snp.DBSNPReader.DBSNP;


public class BroadSNP extends DBSNP {
	int sangerCategoryId;
	String dbSNPId;
	char referenceAlele;
	HashMap<String, SNPSampleInfo> samples = new HashMap<String, SNPSampleInfo>(); 

	public BroadSNP(String id) {
		super(id);
		
	}
	
	public void setCategoryIdFromString(String categoryId) {
		sangerCategoryId = Integer.parseInt(categoryId);
	}
	
	public void addSample(SNPSampleInfo sample) {
		//System.out.println("Adding sample " + sample.getSampleId() + " to SNP " + getName());
		samples.put(sample.getSampleId(), sample);
	}
	
	public Collection<SNPSampleInfo> getSamples() { return samples.values(); }

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

	public String getDbSNPId() {
		return dbSNPId;
	}

	public char getReferenceAlele() {
		return referenceAlele;
	}

	public int getSangerCategoryId() {
		return sangerCategoryId;
	}

	public void setDbSNPId(String dbSNPId) {
		this.dbSNPId = dbSNPId;
	}

	public void setReferenceAlele(char referenceAlele) {
		this.referenceAlele = referenceAlele;
	}

	public void setSangerCategoryId(int sangerCategoryId) {
		this.sangerCategoryId = sangerCategoryId;
	}
	
	public SNPSampleInfo getSampleInfo(String sampleId) {
		return samples.get(sampleId);
	}

}
