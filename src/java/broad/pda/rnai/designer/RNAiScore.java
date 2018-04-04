package broad.pda.rnai.designer;

import broad.pda.datastructures.Alignments;

public class RNAiScore {

	String for5Prime="CCGG";
	String loop="CTCGAG";
	String for3Prime="TTTTTG";
	String rev5Prime="AATTCAAAAA";
	String rev3Prime="";
	
	String kmer;
	Alignments align;
	double rs8Score;
	int numBlastHits;
	String forwardOligo;
	String reverseOligo;
	
	
	public RNAiScore(String kmer, Alignments alignment, double rs8Score, int numBlastHits){
		this.kmer=kmer;
		this.align=alignment;
		this.rs8Score=rs8Score;
		this.numBlastHits=numBlastHits;
		makeCloningOligos(kmer);
	}
	
	
	private void makeCloningOligos(String kmer){
		String reverseComp=ComputeMIRScore.reverseComplement(kmer);
		this.forwardOligo=this.for5Prime+kmer+this.loop+reverseComp+for3Prime;
		this.reverseOligo=this.rev5Prime+kmer+loop+reverseComp+this.rev3Prime;	
	}
	
	public String getKmer(){return this.kmer;}
	public double getScore(){return this.rs8Score;}
	public Alignments getPosition(){return this.align;}
	public int getNumBlastHits(){return this.numBlastHits;}
	
	
	public String toString(){return align.toUCSC()+"\t"+kmer+"\t"+rs8Score+"\t"+forwardOligo+"\t"+reverseOligo+"\t"+numBlastHits;}
	
}
