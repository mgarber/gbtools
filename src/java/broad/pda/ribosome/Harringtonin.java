package broad.pda.ribosome;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.broad.igv.sam.Alignment;

import broad.core.datastructures.IntervalTree;
import broad.core.sequence.SequenceUtils;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GeneScore;

public class Harringtonin {

	/**
	 * RPKM of start codon
	 */
	private Map<RefSeqGene, GeneScore> startPeak;
	
	ContinuousDataAlignmentModel harringtoninData;
	ContinuousDataAlignmentModel expressionData;
	Collection<RefSeqGene> genes;
	
	/**
	 * The significant peaks per gene
	 */
	Map<RefSeqGene, Map<RefSeqGene,GeneScore>> peakScores;
	private double alpha=0.01;
	
	//Score the start position of a transcript
	public Harringtonin(ContinuousDataAlignmentModel harringtoninData, ContinuousDataAlignmentModel expressionModel, Collection<RefSeqGene> genes) throws IOException{
		this.harringtoninData=harringtoninData;
		this.expressionData=expressionModel;
		this.genes=genes;
		scoreStart(this.genes);
		scanWindows(this.genes);
	}
	
	/**
	 * Get the part of the gene overlapping the intervals
	 * @param gene the gene
	 * @param alignments the intervals
	 * @return a RefSeq gene object consisting of the parts of exons overlapping the intervals
	 */
	private static RefSeqGene getOverlapAsGene(RefSeqGene gene, Collection<Alignments> alignments) {
		int min=Integer.MAX_VALUE;
		int max=-Integer.MAX_VALUE;
		
		for(Alignments align: alignments){
			min=Math.min(align.getStart(), min);
			max=Math.max(align.getEnd(), max);
		}

		min=Math.max(min, gene.getStart());
		max=Math.min(max, gene.getEnd());
				
		RefSeqGene rtrn=gene.trimAbsolute(min, max);
		return rtrn;
	}

	
	private void scanWindows(Collection<RefSeqGene> genesToScan) throws IOException {
		this.peakScores=new TreeMap<RefSeqGene, Map<RefSeqGene, GeneScore>>();
		
		
		int i=0;
		for(RefSeqGene gene: genesToScan){
			
			System.err.println("Scanning "+gene.getName() +" "+i+" "+genesToScan.size());
			int geneSize = gene.getSize();
			IntervalTree<Alignment> tree = this.harringtoninData.getAlignmentDataModelStats().getIntervalTree(gene.getChr(), gene.getStart(), gene.getEnd());
			GeneComponent fullScore=RibosomeScoring.scoreFeature(gene, this.harringtoninData, this.harringtoninData.getAlignmentDataModelStats());
			double lambda = fullScore.getGeneScore().getNumberOfReads()/geneSize;
			Collection<RefSeqGene> windows=gene.getWindows(3);
			Map<RefSeqGene, GeneScore> list = new TreeMap<RefSeqGene, GeneScore>();
			
			// Temporary container for consecutive overlapping significant 3mers
			TreeSet<RefSeqGene> overlappingSignificantWindows = new TreeSet<RefSeqGene>();
			
			Iterator<RefSeqGene> iter = windows.iterator();
			
			// Go through the windows in order of location
			while(iter.hasNext()) {
			
				RefSeqGene window = iter.next();
				
				GeneScore windowScore = this.harringtoninData.scoreGene(window, tree, lambda);
				
				if(windowScore.getScanPvalue() < this.alpha) {
					
					if(overlappingSignificantWindows.isEmpty()) {
						overlappingSignificantWindows.add(window);
						continue;
					}
					
					if(window.overlapsGeneInAnyOrientation(overlappingSignificantWindows.last())) {
						overlappingSignificantWindows.add(window);
						continue;						
					}
					
					// if this window does not overlap the previous cluster
					// make new gene out of all the previous overlapping windows and add to list
					Collection<Alignments> mergedWindows = new TreeSet<Alignments>();
					Alignments interval = new Alignments(overlappingSignificantWindows.first().getChr(),overlappingSignificantWindows.first().getStart(),overlappingSignificantWindows.last().getEnd());
					mergedWindows.add(interval);
					RefSeqGene peak = getOverlapAsGene(gene, mergedWindows);
					GeneScore score = this.harringtoninData.scoreGene(peak,tree,lambda);
					list.put(peak, score);
					
					// clear the set of overlapping windows
					
					overlappingSignificantWindows.clear();
					
					// start new set of overlapping windows
					
					overlappingSignificantWindows.add(window);
					
				}
				
			}
			
			// add the last peak
			if(!overlappingSignificantWindows.isEmpty()) {
				Collection<Alignments> mergedWindows = new TreeSet<Alignments>();
				Alignments interval = new Alignments(overlappingSignificantWindows.first().getChr(),overlappingSignificantWindows.first().getStart(),overlappingSignificantWindows.last().getEnd());
				mergedWindows.add(interval);
				RefSeqGene peak = getOverlapAsGene(gene, mergedWindows);
				GeneScore score = this.harringtoninData.scoreGene(peak,tree,lambda);
				GeneScore global=new GeneScore(peak, this.harringtoninData.scoreGene(peak, tree));
				if(global.getScanPvalue()>score.getScanPvalue()){score=global;} //TODO: Test this
				list.put(peak, score);
			}
			
			// keep all peaks for the gene
			this.peakScores.put(gene, list);
			i++;
		}
		
	}

	/**
	 * Score the start codon of each gene; score is RPKM
	 * @param genes
	 * @throws IOException
	 */
	private void scoreStart(Collection<RefSeqGene> genes) throws IOException {
		this.startPeak=new TreeMap<RefSeqGene, GeneScore>();
		for(RefSeqGene gene: genes){
			GeneComponent fullScore=RibosomeScoring.scoreFeature(gene, this.harringtoninData, this.harringtoninData.getAlignmentDataModelStats());
			double lambda=fullScore.getGeneScore().getNumberOfReads()/gene.getSize();
			RefSeqGene start=gene.getStartCodon();
			RefSeqGene peak=this.harringtoninData.getPeak(gene, start);
			
			if(peak!=null){
				GeneScore startScore=new GeneScore(peak, this.harringtoninData.scoreGene(peak, lambda)); //TODO Fix gene size
				this.startPeak.put(gene, startScore);
			}
		}
	}

	/**
	 * Write start codong score for each gene
	 * @param save output file
	 * @param append do not overwrite
	 * @throws IOException
	 */
	private void write(String save, boolean append) throws IOException {
		FileWriter writer=new FileWriter(save, append);
		
		if(!append){
			String headers = "Gene\t";
			headers += "GeneLength\t";
			headers += "ExpressionLocalLambda\t";
			headers += "ExpressionPValue\t";
			headers += "ExpressionRPKM\t";
			headers += "ExpressionEnrichment\t";
			headers += "StartCodonLocalLambda\t";
			headers += "StartCodonPvalue\t";
			headers += "StartCodonRPKM\t";
			headers += "StartCodonEnrichment\t";
			writer.write(headers + "\n");
		}
		
		for(RefSeqGene gene: this.genes){
			
			GeneScore start = this.startPeak.get(gene);
			GeneScore expression = new GeneScore(gene, this.expressionData.scoreGene(gene));
			
			String line = gene.getName() + "\t";
			line += expression.getGeneLength() + "\t";
			line += expression.getLocalLambda() + "\t";
			line += expression.getScanPvalue() + "\t";
			line += expression.getRPKM() + "\t";
			// expression enrichment with respect to whole chromosome
			line += expression.getEnrichment() + "\t";
			
			if(start != null){
				line += start.getLocalLambda() + "\t";
				line += start.getScanPvalue() + "\t";
				line += start.getRPKM() + "\t";
				// peak enrichment with respect to gene
				line += start.getEnrichment() + "\t";
			}
			else{
				line += "0.0\t1.0\t0.0\t0.0";
			}
			
			writer.write(line + "\n");
			
		}
		
		writer.close();
	}
	
	private void writePeaks(String save, boolean append) throws IOException {
		FileWriter writer=new FileWriter(save, append);
		
		for(RefSeqGene gene: this.peakScores.keySet()){
			for(RefSeqGene w: this.peakScores.get(gene).keySet()){
				GeneScore window= this.peakScores.get(gene).get(w);
				double p=window.getScanPvalue();
				w.setName("p="+p);
				if(p<this.alpha){writer.write(w+"\n");}
			}
		}
		
		writer.close();
	}
	
	private void writeStartCodonPeaks(String save, boolean append) throws IOException {
		FileWriter writer=new FileWriter(save, append);
		
		for(RefSeqGene gene : this.startPeak.keySet()) {
			
			RefSeqGene start=gene.getStartCodon();
			RefSeqGene peak=this.harringtoninData.getPeak(gene, start);

			GeneScore window = this.startPeak.get(gene);
			double p = window.getScanPvalue();
				
			peak.setName("p="+p);
			if(p < this.alpha) writer.write(peak + "\n");
			
		}
		
		writer.close();
		
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>4){
			String bamFile=args[0];
			Map<String, Collection<RefSeqGene>> genesByChr=BEDFileParser.loadDataByChr(new File(args[1]));
			String sizeFile=args[2];
			String save=args[3];
			String expressionBAM=args[4];
			String chrToUse=null;
			if(args.length>5){chrToUse=args[5];}
			boolean append=false;
			for(String chr: genesByChr.keySet()){
				if(chrToUse==null || !chrToUse.startsWith("chr") || chrToUse.equalsIgnoreCase(chr)){
					ContinuousDataAlignmentModel model=SequenceUtils.getDataModel(bamFile, sizeFile, chr, false);
					ContinuousDataAlignmentModel expressionModel=SequenceUtils.getDataModel(expressionBAM, sizeFile, chr, false);
					Harringtonin h=new Harringtonin(model, expressionModel, genesByChr.get(chr));
					h.write(save, append);
					h.writePeaks(save+".peaks.bed", append);
					h.writeStartCodonPeaks(save + ".peaks.startCodon.bed", append);
					append = true;
				}
			}
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam file \n args[1]=genes \n args[2]=size file \n args[3]=save \n args[4]=expression bam \n args[5]=chr (optional)";
	
}
