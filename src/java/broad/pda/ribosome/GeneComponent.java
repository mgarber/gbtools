package broad.pda.ribosome;

import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.GeneScore;

public class GeneComponent implements Comparable<GeneComponent>{

	private RefSeqGene gene;
	private RefSeqGene CDS;
	private RefSeqGene UTR5;
	private RefSeqGene UTR3;
	private RefSeqGene intron;
	
	private GeneScore geneScore;
	private GeneScore UTR5Score;
	private GeneScore UTR3Score;
	private GeneScore intronScore;
	private GeneScore CDSScore;
	
	boolean hasWindowScores;
	
	private Map<RefSeqGene, GeneScore> geneWindowScores;
	private Map<RefSeqGene, GeneScore> predictedORFScores;
	private Map<RefSeqGene, GeneScore> predictedORFRNASeqScores;
	
	private GeneScore maxGeneWindow;
	private GeneScore max5UTRWindow;
	private GeneScore max3UTRWindow;
	private GeneScore maxCDSWindow;
	
	private double geneTE;
	private double CDSTE;
	private double UTR5TE;
	private double UTR3TE;
	private double intronTE;
	
	private GeneScore pairedGeneExpression;
	private GeneScore paired3UTRExpression;
	private GeneScore paired5UTRExpression;
	private GeneScore pairedCDSExpression;
	
	private GeneScore geneRNASeqScore;
	private GeneScore cdsRNASeqScore;
	private GeneScore utr3RNASeqScore;
	private GeneScore utr5RNASeqScore;
	private GeneScore intronRNASeqScore;
	
	
	
	public GeneComponent(RefSeqGene gene){
		this.gene=gene;
		this.CDS=gene.getCDS();
		this.UTR5=gene.get5UTRGene();
		this.UTR3=gene.get3UTRGene();
		this.intron=gene.getIntrons();
		this.hasWindowScores=false;
	}
	
	public RefSeqGene getGene(){return gene;}
	public RefSeqGene get3UTR(){return this.UTR3;}
	public RefSeqGene get5UTR(){return this.UTR5;}
	public RefSeqGene getIntron(){return this.intron;}
	public RefSeqGene getCDS(){return this.CDS;}

	public String getGeneName() {return this.gene.getName();}
	
	public void setGeneScore(double[] scores) {
		this.geneScore = new GeneScore(this.gene, scores);
	}

	public GeneScore getGeneScore() {
		return geneScore;
	}

	public void setUTR5Score(double[] scores) {
		UTR5Score = new GeneScore(this.UTR5, scores);
	}

	public GeneScore getUTR5Score() {
		return UTR5Score;
	}

	public void setUTR3Score(double[] scores) {
		UTR3Score = new GeneScore(this.UTR3, scores);
	}
	
	public GeneScore getUTR3Score() {
		return UTR3Score;
	}

	public void setIntronScore(double[] scores) {
		this.intronScore = new GeneScore(this.intron, scores);
	}

	public GeneScore getIntronScore() {
		return intronScore;
	}

	public void setCDSScore(double[] scores) {
		CDSScore = new GeneScore(this.CDS, scores);
	}

	public GeneScore getCDSScore() {
		return CDSScore;
	}

	public String getName() {
		return this.gene.getName();
	}

	public double getGeneEnrichment() {
		return this.geneScore.getEnrichment();
	}

	public double getCDSEnrichment() {
		return this.CDSScore.getEnrichment();
	}

	public double getUTR5Enrichment() {
		return this.UTR5Score.getEnrichment();
	}

	public double getUTR3Enrichment() {
		return this.UTR3Score.getEnrichment();
	}

	public double getIntronEnrichment() {
		return this.intronScore.getEnrichment();
	}

	public double getSignificance() {
		return this.geneScore.getScanPvalue();
	}

	public void setGeneWindowScores(Map<RefSeqGene, double[]> geneWindowScores, int windowSize) {
		//TODO Use these scores to populate CDS and UTR counts;
		
		this.hasWindowScores=true;
		Map<RefSeqGene, GeneScore> rtrn= new TreeMap<RefSeqGene, GeneScore>();
		for(RefSeqGene gene: geneWindowScores.keySet()){
			GeneScore score=new GeneScore(gene, geneWindowScores.get(gene));
			rtrn.put(gene, score);
		}
		this.geneWindowScores=rtrn;
		
		//Populate variables
		this.maxGeneWindow=determineMax();
		if(this.hasCDS()){
			this.maxCDSWindow=determineMax(getCDS(), windowSize);
			this.max3UTRWindow=determineMax(get3UTR(), windowSize);
			this.max5UTRWindow=determineMax(get5UTR(), windowSize);
		}
	}

	private GeneScore determineMax(RefSeqGene region, int windowSize) {
		if(region==null){return null;}
		Collection<RefSeqGene> subGenes=region.getWindows(windowSize);
		return determineMax(subGenes);
	}

	
	private GeneScore determineMax() {
		Collection<RefSeqGene> subGenes=geneWindowScores.keySet();
		return determineMax(subGenes);
	}
	
	//TODO Get gene max region and value, UTRs, CDS
	private GeneScore determineMax(Collection<RefSeqGene> subGenes) {
		GeneScore max=null;
		// go through each window and get max
		//save region and score
		for (RefSeqGene window: subGenes){
			GeneScore score=geneWindowScores.get(window);
			max=max(score, max);
		}
		return max;
	}

	private GeneScore max(GeneScore score, GeneScore max) {
		if(max==null){return score;}
		if(score==null){return max;}
		
		double v1=score.getFullyContainedNumberOfReads();
		double v2=max.getFullyContainedNumberOfReads();
		
		if(v1>v2){return score;}
		return max;
		
	}

	public Map<RefSeqGene, GeneScore> getGeneWindowScores() {
		return geneWindowScores;
	}

	public boolean hasCDS() {
		if(this.getCDS()==null){return false;}
		return true;
	}

	public boolean hasIntron() {
		if(this.gene.getNumExons()>1){return true;}
		return false;
	}

	//TODO Get max gene score
	public GeneScore getMaxGeneWindow(){
		return this.maxGeneWindow;
	}
	
	public GeneScore getMax5UTRWindow(){return this.max5UTRWindow;}
	public GeneScore getMax3UTRWindow(){return this.max3UTRWindow;}
	public GeneScore getMaxCDSWindow(){return this.maxCDSWindow;}

	public void setGeneTE(double geneTE) {
		this.geneTE = geneTE;
	}

	public double getGeneTE() {
		return geneTE;
	}

	public void setCDSTE(double cDSTE) {
		CDSTE = cDSTE;
	}

	public void setUTR5TE(double uTR5TE) {
		UTR5TE = uTR5TE;
	}

	public double getUTR5TE(double alpha) {
		if(this.UTR5Score!=null && this.utr5RNASeqScore!=null){
			return RibosomeScoring.computeTE(UTR5Score, utr5RNASeqScore, alpha);
		}
		else{
			System.err.println("Variables are not set");
			return -999;
		}
	}

	public void setUTR3TE(double uTR3TE) {
		UTR3TE = uTR3TE;
	}

	public double getUTR3TE(double alpha) {
		if(this.UTR3Score!=null && this.utr3RNASeqScore!=null){
			return RibosomeScoring.computeTE(UTR3Score, utr3RNASeqScore, alpha);
		}
		else{
			System.err.println("Variables are not set");
			return -999;
		}
	}

	
	public void setPairedGeneExpression(GeneScore pairedGeneExpression) {
		this.pairedGeneExpression = pairedGeneExpression;
	}

	public GeneScore getPairedGeneExpression() {
		return pairedGeneExpression;
	}

	public void setPaired3UTRExpression(GeneScore paired3UTRExpression) {
		this.paired3UTRExpression = paired3UTRExpression;
	}

	public GeneScore getPaired3UTRExpression() {
		return paired3UTRExpression;
	}

	public void setPaired5UTRExpression(GeneScore paired5UTRExpression) {
		this.paired5UTRExpression = paired5UTRExpression;
	}

	public GeneScore getPaired5UTRExpression() {
		return paired5UTRExpression;
	}

	public void setPairedCDSExpression(GeneScore pairedCDSExpression) {
		this.pairedCDSExpression = pairedCDSExpression;
	}

	public GeneScore getPairedCDSExpression() {
		return pairedCDSExpression;
	}

	public void setGeneRNASeqScore(GeneScore rnaSeqGeneScore) {
		this.geneRNASeqScore = rnaSeqGeneScore;
	}
	
	public void setGeneRNASeqScore(double[] rnaSeqGeneScore) {
		this.geneRNASeqScore = new GeneScore(this.gene, rnaSeqGeneScore);
	}

	public GeneScore getGeneRNASeqScore() {
		return geneRNASeqScore;
	}

	public void setCDSRNASeqScore(GeneScore cdsRNASeqScore) {
		this.cdsRNASeqScore = cdsRNASeqScore;
	}
	
	public void setCDSRNASeqScore(double[] cdsRNASeqScore) {
		this.cdsRNASeqScore = new GeneScore(this.getCDS(), cdsRNASeqScore);
	}

	public GeneScore getCDSRNASeqScore() {
		return cdsRNASeqScore;
	}

	public void setUTR3RNASeqScore(GeneScore utr3RNASeqScore) {
		this.utr3RNASeqScore = utr3RNASeqScore;
	}
	
	public void setUTR3RNASeqScore(double[] utr3RNASeqScore) {
		this.utr3RNASeqScore = new GeneScore(this.get3UTR(), utr3RNASeqScore);
	}

	
	public GeneScore getUTR3RNASeqScore() {
		return utr3RNASeqScore;
	}

	public void setUTR5RNASeqScore(GeneScore utr5RNASeqScore) {
		this.utr5RNASeqScore = utr5RNASeqScore;
	}
	
	public void setUTR5RNASeqScore(double[] utr5RNASeqScore) {
		this.utr5RNASeqScore = new GeneScore(this.get5UTR(), utr5RNASeqScore);
	}

	public GeneScore getUTR5RNASeqScore() {
		return utr5RNASeqScore;
	}

	public void setIntronRNASeqScore(GeneScore intronRNASeqScore) {
		this.intronRNASeqScore = intronRNASeqScore;
	}
	
	public void setIntronRNASeqScore(double[] intronRNASeqScore) {
		this.intronRNASeqScore = new GeneScore(this.getIntron(), intronRNASeqScore);
	}

	public GeneScore getIntronRNASeqScore() {
		return intronRNASeqScore;
	}

	public void setPredictedORFScores(Map<RefSeqGene, double[]> predictedORFScores) {
		Map<RefSeqGene, GeneScore> tmp=new TreeMap<RefSeqGene, GeneScore>();
		for(RefSeqGene gene: predictedORFScores.keySet()){
			GeneScore score=new GeneScore(gene, predictedORFScores.get(gene));
			tmp.put(gene, score);
		}
		this.predictedORFScores = tmp;
	}
	
	public Map<RefSeqGene, GeneScore> getPredictedORFScores() {
		return predictedORFScores;
	}

	public void setPredictedORFRNASeqScores(Map<RefSeqGene, GeneScore> predictedORFRNASeqScores) {
		this.predictedORFRNASeqScores = predictedORFRNASeqScores;
	}

	public Map<RefSeqGene, GeneScore> getPredictedORFRNASeqScores() {
		return predictedORFRNASeqScores;
	}

	public boolean hasPredictedORFs() {
		if(this.predictedORFScores!=null && !this.predictedORFScores.isEmpty()){return true;}
		return false;
	}

	public void setIntronTE(double computeTE) {
		this.intronTE=computeTE;
	}
	
	public double getIntronTE(){
		return this.intronTE;
	}

	public double getGeneTE(double alpha) {
		if(this.geneScore!=null && this.geneRNASeqScore!=null){
			return RibosomeScoring.computeTE(geneScore, geneRNASeqScore, alpha);
		}
		else{
			System.err.println("Variables are not set");
			return -999;
		}
	}
	
	public double getCDSTE(double alpha) {
		if(this.CDSScore!=null && this.cdsRNASeqScore!=null){
			return RibosomeScoring.computeTE(CDSScore, cdsRNASeqScore, alpha);
		}
		else{
			System.err.println("Variables are not set");
			return -999;
		}
	}
	
	public double getCDSTEAllReads(double alpha) {
		if(this.CDSScore!=null && this.cdsRNASeqScore!=null){
			return RibosomeScoring.computeTEAllReads(CDSScore, cdsRNASeqScore, alpha);
		}
		else{
			System.err.println("Variables are not set");
			return -999;
		}
	}
	
	public double getGeneTEAllReads(double alpha) {
		if(this.geneScore!=null && this.geneRNASeqScore!=null){
			return RibosomeScoring.computeTEAllReads(geneScore, geneRNASeqScore, alpha);
		}
		else{
			System.err.println("Variables are not set");
			return -999;
		}
	}

	@Override
	public String toString() {
		return this.gene.toString();
	}
	
	/**
	 * Returns comparison of the underlying genes as implemented in RefSeqGene
	 */
	@Override
	public int compareTo(GeneComponent g) {
		return this.gene.compareTo(g.gene);
	}

	/**
	 * Get RRS score
	 * @param alpha Scan P value cutoff for gene
	 * @param normalizeByExpression Whether to normalize ribosome ratio by background expression ratio
	 * @param fullyContainedReadsInCds Whether to count only fully contained reads in CDS (alternative is all reads overlapping CDS)
	 * @return RRS score = ribosome read count ratio of CDS to UTR
	 */
	public double getRRS(double alpha, boolean normalizeByExpression, boolean fullyContainedReadsInCds) {
		
		if(gene.hasCDS() && gene.get3UTRGene() != null && getGeneScore().getScanPvalue() < alpha){
			
			double utr3reads = getUTR3Score().getFullyContainedNumberOfReads();
			double utr3readsRnaSeq = getUTR3RNASeqScore().getFullyContainedNumberOfReads();
			
			double cdsReads;
			double cdsReadsRnaSeq;
			if(fullyContainedReadsInCds) {
				cdsReads = getCDSScore().getFullyContainedNumberOfReads();
				cdsReadsRnaSeq = getCDSRNASeqScore().getFullyContainedNumberOfReads();
			}
			else {
				cdsReads = getCDSScore().getNumberOfReads();
				cdsReadsRnaSeq = getCDSRNASeqScore().getNumberOfReads();
			}
			
			int cdsSize = gene.getCDS().getSize();
			int utrSize = gene.get3UTRGene().getSize();
			
			if(cdsReads == 0 && utr3reads == 0) return -99;
			
			double cdsRibosomeRPKM = (cdsReads + 1)/cdsSize;
			double utr3RibosomeRPKM = (utr3reads + 1)/utrSize;
			double ribosomeRatio = (cdsRibosomeRPKM)/(utr3RibosomeRPKM);
						
			if(!normalizeByExpression) return ribosomeRatio;
			
			if(utr3reads == 0 && utr3readsRnaSeq == 0) return -99;
			
			double cdsExpressionRPKM = (cdsReadsRnaSeq + 1)/cdsSize;
			double utr3ExpressionRPKM = (utr3readsRnaSeq + 1)/utrSize;
			double expressionRatio = (cdsExpressionRPKM)/(utr3ExpressionRPKM);
			
			return ribosomeRatio / expressionRatio;
			
		}
		else return -99;
	}





	
	
}
