package broad.pda.agilent;

import java.io.BufferedWriter;
import java.io.IOException;

import broad.pda.gene.RefSeqGene;
import broad.pda.gene.RefSeqGeneWithIsoforms;

public class GeneAnnotationInfo {

	int highConfidence=0;
	int bothAssemblies=0;
	int assemblyAndExternal=0;
	int bwnSets=0;
	int numIsoforms=0;
	String geneName;
	RefSeqGeneWithIsoforms mergedRef;
	RefSeqGene longestChainIso=null;
	int K4K36=0;
	int iK4K36=0;
	int K4K36lincs =0;
	int scanStat=0;
	double maxUniformityScr_cov=0;
	double maxUniformityScr=0;
	double maxExpVal=0;
	double maxJSspecificity=0;
	String maxJSspecificity_centroid="";
	String neighborString=""; //first 12 tabs has bed info
	int NR=0;
	int GENCODE=0;
	int UCSC=0;
	int SEQ=0;
	
	public GeneAnnotationInfo(RefSeqGeneWithIsoforms g){
		this.mergedRef=g;
	}
	
	public RefSeqGeneWithIsoforms getMergedRef() { return this.mergedRef;}
	public void setNeighborsLine(String line) { this.neighborString=line; }

	public void setMaxExp(double d) {this.maxExpVal=d;}
		

	public void setMaxJSspecificity(double d) {this.maxJSspecificity=d;	}

	public void setMaxJSspecificity_centroid(double d){this.maxJSspecificity_centroid=String.valueOf(d);}
	public void setMaxJSspecificity_centroid(String string) { this.maxJSspecificity_centroid=string;	}

	public void setk4k36Lincs(int i) {this.K4K36lincs=i;	}
	
	public void setNR(int i) {this.NR=i;	}
	public void setGENCODE(int i) {this.GENCODE=i;	}
	public void setUCSC(int i) {this.UCSC=i;	}

	
	//bw.write("chr\tstart\tend\tname\tscore\tstrand\ttStart\ttEnd\trgb\tblockCnt\tblockSizes\tblockStarts");
	//bw.write("\tNumIsoforms\tHighConfidenceSet\tCompatibleBetweenAssemblies\tCompatibleBetweenAssemblyAndExternal\tCompatibleBetweenSets");
	//bw.write("MaxUniformityScore\tCoverage\tScanStatExp\tK4K36\tK4K36Independan\tK4K36lincs\");
	//neighbors

	public  void writeInfo(BufferedWriter bw) throws IOException {

		int length=0;
		int exNum=0;
		if (this.longestChainIso!=null){
			exNum=this.longestChainIso.getNumExons();
			length=this.longestChainIso.getTranscriptLength();
		}
		bw.write(this.mergedRef.toBED(false));
		bw.write("\t"+this.numIsoforms);
		bw.write("\t"+this.highConfidence+"\t"+this.bothAssemblies+"\t"+this.assemblyAndExternal+"\t"+this.bwnSets);
		bw.write("\t"+this.maxUniformityScr+"\t"+this.maxUniformityScr_cov+"\t");
		bw.write(this.scanStat+"\t"+this.K4K36+"\t"+this.iK4K36+"\t"+this.K4K36lincs);
		bw.write("\t"+this.maxExpVal+"\t"+this.maxJSspecificity+"\t"+this.maxJSspecificity_centroid);
		bw.write("\t"+this.NR+"\t"+this.UCSC+"\t"+this.GENCODE+"\t"+this.SEQ);
		bw.write("\t"+exNum+"\t"+length);
		String[] tmp= neighborString.split("\t");
		for (int i=12;i<tmp.length;i++)
			bw.write("\t"+tmp[i]);
		bw.write("\n");
	}

	public void incrementNumIso() { this.numIsoforms++;}

	public void updateHighConfidence(Integer i, String[] a) {
		if (Integer.valueOf(a[i])> 0) this.highConfidence=1;	}

	public void updateCompatibleBetweenAssemblies(Integer i, String[] a) {
		if (Integer.valueOf(a[i])> 0) this.bothAssemblies=1;	}
	
	public void updateCompatibleBetweenAssemblyAndExternal(Integer i, String[] a) {
		if (Integer.valueOf(a[i])> 0) this.assemblyAndExternal=1;	}
	
	public void updateCompatibleBetweenSets(Integer i, String[] a) {
		if (Integer.valueOf(a[i])> 0) this.bwnSets=1;	}
	
	public void updateNR(Integer i, String[] a) {
		if (Integer.valueOf(a[i])> 0){
			this.NR=1;
			this.highConfidence=1;
		}
	}
	public void updateGENCODE(Integer i, String[] a) {if (Integer.valueOf(a[i])> 0) this.GENCODE=1;	}
	public void updateUCSC(Integer i, String[] a) {if (Integer.valueOf(a[i])> 0) this.UCSC=1;	}
	public void updateSEQ(boolean a) {if (a) this.SEQ=1;	}
	
	public void updateUniformityScore(Integer i,Integer c, String[] a) {
		Double d=Double.valueOf(a[i]);
		if (d>this.maxUniformityScr){
			this.maxUniformityScr=d;
			this.maxUniformityScr_cov=Double.valueOf(a[c]);
		}
	}
	
	public void updateScanStatExp(Integer i, String[] a) {
		Double d=Double.valueOf(a[i]);
		if (d> 0) this.scanStat=1;	}
	
	public void updateK4K36(Integer i, String[] a) {
		Double d=Double.valueOf(a[i]);
		if (d> 0) this.K4K36=1;	}
	
	public void updateK4K36Independent(Integer i, String[] a) {
		Double d=Double.valueOf(a[i]);
		if (d> 0) this.iK4K36=1;
	}

	public void writeInfo(BufferedWriter bw, String emptyArr) throws IOException {
		//System.err.println(this.mergedRef.getName()+"\t"+this.neighborString);
		if (this.neighborString.equalsIgnoreCase(""))
			this.neighborString=emptyArr;
		writeInfo(bw);
		
	}

	public void setlongestIso(RefSeqGene longestIso) { this.longestChainIso=longestIso;}
		
	

	
	
	
}