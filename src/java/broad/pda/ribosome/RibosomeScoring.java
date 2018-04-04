package broad.pda.ribosome;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.TreeMap;

import org.broad.igv.sam.Alignment;

import broad.core.datastructures.IntervalTree;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.AlignmentDataModelStats;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GeneScore;

public class RibosomeScoring {

	
	public static GeneComponent scoreFeature(RefSeqGene gene, ContinuousDataAlignmentModel model, AlignmentDataModelStats data) throws IOException{
		IntervalTree<Alignment> tree=data.getIntervalTreeCached(gene.getChr(), gene.getStart(), gene.getEnd());
		if(model.hasDataForChromosome(gene.getChr())){
			//for each gene get UTRs, CDSes
			GeneComponent component=new GeneComponent(gene);
			component.setGeneScore(model.scoreGene(component.getGene(), tree));
			if(component.hasCDS()){component.setCDSScore(model.scoreGene(component.getCDS(), tree));}
			if(component.get3UTR()!=null){component.setUTR3Score(model.scoreGene(component.get3UTR(), tree));}
			if(component.get5UTR()!=null){component.setUTR5Score(model.scoreGene(component.get5UTR(),tree));}
			if(component.hasIntron()){component.setIntronScore(model.scoreGene(component.getIntron(), tree));}
			return component;
		}
		return null;
	}
	
	
	public static GeneComponent scoreFeature(RefSeqGene gene, ContinuousDataAlignmentModel ribosome, ContinuousDataAlignmentModel expression) throws IOException{
		GeneComponent rtrn=scoreFeature(gene, ribosome, ribosome.getAlignmentDataModelStats());
		addExpression(rtrn, expression, expression.getAlignmentDataModelStats());
		return rtrn;
	}
	
	public static Collection<GeneComponent> scoreFeature(RefSeqGene gene, Collection<RefSeqGene> orfs, ContinuousDataAlignmentModel ribosomeModel, ContinuousDataAlignmentModel expressionModel) throws IOException{
		Collection<GeneComponent> rtrn=new HashSet<GeneComponent>();
		//Make ribosome tree for gene
		IntervalTree<Alignment> ribosomeTree=ribosomeModel.getAlignmentDataModelStats().getIntervalTreeCached(gene.getChr(), gene.getStart(), gene.getEnd());
		//Make expression tree for gene
		IntervalTree<Alignment> expressionTree=expressionModel.getAlignmentDataModelStats().getIntervalTreeCached(gene.getChr(), gene.getStart(), gene.getEnd());
		for(RefSeqGene orf: orfs){
			GeneComponent component=scoreFeature(orf, ribosomeModel, expressionModel);
			
			//TODO Consider using trees
			/*GeneComponent component=new GeneComponent(orf);
			component.setGeneScore(ribosomeModel.scoreGene(component.getGene(), ribosomeTree));
			component.setGeneRNASeqScore(expressionModel.scoreGene(component.getGene(), expressionTree));*/
			
			rtrn.add(component);
		}
		return rtrn;	
	}
	
	
	public static void addExpression(GeneComponent feature, ContinuousDataAlignmentModel model, AlignmentDataModelStats data) throws IOException {
		RefSeqGene gene=feature.getGene();
		IntervalTree<Alignment> tree=data.getIntervalTree(gene.getChr(), gene.getStart(), gene.getEnd());
		if(model.hasDataForChromosome(gene.getChr())){
			feature.setGeneRNASeqScore(model.scoreGene(feature.getGene(), tree));
			//for each gene get UTRs, CDSes
			if(feature.hasCDS()){feature.setCDSRNASeqScore(model.scoreGene(feature.getCDS(), tree));}
			if(feature.get3UTR()!=null){feature.setUTR3RNASeqScore(model.scoreGene(feature.get3UTR(), tree));}
			if(feature.get5UTR()!=null){feature.setUTR5RNASeqScore(model.scoreGene(feature.get5UTR(),tree));}
			if(feature.hasIntron()){feature.setIntronRNASeqScore(model.scoreGene(feature.getIntron(), tree));}
		}	
	}
	
	public static double computeTE(GeneScore ribo, GeneScore rna, double alpha) {
		double TE=-99;
		if(ribo!=null && rna!=null){
			if(rna.getScanPvalue()<alpha){
				if(!(ribo.getFullyContainedNumberOfReads() == 0 && rna.getFullyContainedNumberOfReads() == 0)) {
					TE=(ribo.getFullyContainedNumberOfReads()+1.0)/(rna.getFullyContainedNumberOfReads()+1.0);
				}
			}
		}		
		return TE;
	}
	
	public static double computeTEAllReads(GeneScore ribo, GeneScore rna, double alpha) {
		double TE=-99;
		if(ribo!=null && rna!=null){
			if(rna.getScanPvalue()<alpha){
				if(!(ribo.getNumberOfReads() == 0 && rna.getNumberOfReads() == 0)) {
					TE=(ribo.getNumberOfReads()+1.0)/(rna.getNumberOfReads()+1.0);
				}
			}
		}		
		return TE;
	}
	
	
	public static Map<GeneComponent, Double> getRRS(Collection<GeneComponent> genes, double alpha, boolean normalizeByExpression, boolean fullyContainedReadsInCds) {
		Map<GeneComponent, Double> rtrn = new TreeMap<GeneComponent, Double>();
		for(GeneComponent gene : genes) {
			rtrn.put(gene, gene.getRRS(alpha, normalizeByExpression, fullyContainedReadsInCds));
		}
		return rtrn;
	}
	
	public static double getRRS(GeneComponent gene, double alpha, boolean normalizeByExpression, boolean fullyContainedReadsInCds) {
		return gene.getRRS(alpha, normalizeByExpression, fullyContainedReadsInCds);
	}
	
	
	public static Map<String, GeneComponent> getMaxRRSPerGene(Collection<GeneComponent> genes, double alpha, boolean normalizeByExpression, boolean fullyContainedReadsInCds) {
		
		Map<String, GeneComponent> rtrn = new TreeMap<String, GeneComponent>();
		
		if(genes == null) return rtrn;
		
		for(GeneComponent gene : genes) {
			if(!rtrn.containsKey(gene.getGeneName())) {
				rtrn.put(gene.getGeneName(), gene);
				continue;
			}
			if(gene.getRRS(alpha, normalizeByExpression, fullyContainedReadsInCds) > rtrn.get(gene.getGeneName()).getRRS(alpha, normalizeByExpression, fullyContainedReadsInCds)) {
				rtrn.put(gene.getGeneName(), gene);
			}
		}
		
		return rtrn;
		
	}
	
	public static void writeMaxRRSPerGene(Collection<GeneComponent> genes, double alpha, String save, boolean append, boolean normalizeByExpression, boolean fullyContainedReadsInCds) throws IOException  {
		
		Map<String, GeneComponent> max = getMaxRRSPerGene(genes,alpha, normalizeByExpression, fullyContainedReadsInCds);
		FileWriter w = new FileWriter(save, append);
		
		if(!append) w.write("Chr\tStart\tEnd\tName\tScore\tStrand\tCDS_start\tCDS_end\tRGB\tExonCount\tExonSizes\tExonStarts\tCDS_RPKM\tCDS_TE\tUTR3_RPKM\tRRS\n");
		
		for(GeneComponent gene : max.values()) {
			double cds = -99;
			double cdste = -99;
			if(gene.hasCDS()) {
				cds = gene.getCDSScore().getRPKM();
				cdste = gene.getCDSTE(alpha);
			}
			double utr3 = -99;
			if(gene.hasCDS() && gene.get3UTR() != null) {
				utr3 = gene.getUTR3Score().getRPKM();
			}
			double rrs = getRRS(gene, alpha, normalizeByExpression, fullyContainedReadsInCds);
			w.write(gene.toString() + "\t" + cds + "\t" + cdste + "\t" + utr3 + "\t" + rrs + "\n");
		}
		w.close();
	}
	
	/*public static double getRRS(RefSeqGene gene, ContinuousDataAlignmentModel ribosomeData, ContinuousDataAlignmentModel expressionData, double alpha) throws IOException {
		
		System.err.println("Getting RRS for gene " + gene.getName());
		System.err.println("Scoring gene...");
		GeneComponent c = scoreFeature(gene, ribosomeData, expressionData);
		GeneScore geneScore = c.getGeneScore();
		
		if(gene.hasCDS() && gene.get3UTRGene() != null && geneScore.getScanPvalue() < alpha) {
			
			System.err.println("alpha = " + alpha);
			double cdsRPKM = (c.getCDSScore().getFullyContainedNumberOfReads() + 1) / gene.getCDS().getSize();
			System.err.println("CDS RPKM = " + cdsRPKM);
			RefSeqGene newUtr = gene.get3UtrUpToNextStartCodon();
			System.err.println("Got conservative 3'UTR up to next start codon.");
			double utr3RPKM = (ribosomeData.getData().getData().getCountsPerAlignmentFullyContained(newUtr, ribosomeData.getData().getIntervalTreeCached(newUtr.getChr(), newUtr.getStart(), newUtr.getEnd()), 0) + 1) / newUtr.getSize();
			System.err.println("3'UTR RPKM = " + utr3RPKM);
			return cdsRPKM / utr3RPKM;
			
		}
		return -99;
	}*/
	
	public static void writeRRS(String save, Collection<GeneComponent> genes, boolean append, double alpha, boolean normalizeByExpression, boolean fullyContainedReadsInCds) throws IOException{
		FileWriter writer=new FileWriter(save, append);

		if(!append) writer.write("Chr\tStart\tEnd\tName\tScore\tStrand\tCDS_start\tCDS_end\tRGB\tExonCount\tExonSizes\tExonStarts\tCDS_RPKM\tCDS_TE\tUTR3_RPKM\tRRS\n");
		
		if(genes != null) {
			for(GeneComponent gene: genes){
				double cds = -99;
				double cdste = -99;
				if(gene.hasCDS()) {
					cds = gene.getCDSScore().getRPKM();
					cdste = gene.getCDSTE(alpha);
				}
				double utr = -99;
				if(gene.hasCDS() && gene.get3UTR() != null) {
					utr = gene.getUTR3Score().getRPKM();
				}
				double rrs = getRRS(gene, alpha, normalizeByExpression, fullyContainedReadsInCds);
				writer.write(gene.toString() + "\t" + cds + "\t" + cdste + "\t" + utr + "\t" + rrs + "\n");
			}
		}
		
		writer.close();
		
	}
	
}
