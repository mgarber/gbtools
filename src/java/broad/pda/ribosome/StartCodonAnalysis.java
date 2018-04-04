package broad.pda.ribosome;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.broad.igv.Globals;
import org.broad.igv.sam.Alignment;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceUtils;
import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;
import broad.pda.ribosome.misc.FindORFs;
import broad.pda.seq.segmentation.AlignmentDataModelStats;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GeneScore;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

public class StartCodonAnalysis {

	//TODO Add Expression levels and significance calculation
	//TODO Compute significance of peak using local lambda estimated off the gene
	//TODO Score random occurring AUGs for pile-ups
	
	private Map<RefSeqGene, RefSeqGene> geneToStart;
	Map<RefSeqGene, GeneScore> coverage;
	private double alpha=0.01;
	
	public StartCodonAnalysis(String bamFile, String sizeFile, Collection<RefSeqGene> genes, String save, String chr, Sequence chrSequence) throws IOException{
		Globals.setHeadless(true);
		AlignmentDataModelStats data=new AlignmentDataModelStats(new GenericAlignmentDataModel(bamFile, sizeFile), null, chr);
		ContinuousDataAlignmentModel model=new ContinuousDataAlignmentModel(data);
		Map<RefSeqGene, Collection<GeneScore>> scores=scoreStartCodon(genes, data, model, chrSequence);
		Map<RefSeqGene, Collection<RefSeqGene>> peaks=getPeaks(scores, model, data);
		Map<RefSeqGene, Collection<GeneScore>> randomScores=randomScores(genes, data, model, peaks);
		write(scores, randomScores, save);
		write(randomScores, save+".random.bed");
		write(scores, save+".orfs.bed");
	}
	
	private Map<RefSeqGene, Collection<RefSeqGene>> getPeaks(Map<RefSeqGene, Collection<GeneScore>> scores, ContinuousDataAlignmentModel model, AlignmentDataModelStats data) throws IOException {
		Map<RefSeqGene, Collection<RefSeqGene>> rtrn=new TreeMap<RefSeqGene, Collection<RefSeqGene>>();
		
		for(RefSeqGene gene: scores.keySet()){
			IntervalTree<Alignment> tree=data.getIntervalTree(gene.getChr(), gene.getStart(), gene.getEnd());
			Collection<RefSeqGene> peaks=new TreeSet<RefSeqGene>();
			for(GeneScore score: scores.get(gene)){
				GeneScore peak=model.getPeak(gene, score.getGene(), tree);
				peaks.add(peak.getGene());
			}
			rtrn.put(gene, peaks);
		}
		
		return rtrn;
	}

	private Map<RefSeqGene, Collection<GeneScore>> randomScores(Collection<RefSeqGene> genes, AlignmentDataModelStats data, ContinuousDataAlignmentModel model, Map<RefSeqGene, Collection<RefSeqGene>> peaks) throws IOException {
		Map<RefSeqGene, Collection<GeneScore>> rtrn=new TreeMap<RefSeqGene, Collection<GeneScore>>();
		
		int counter=0;
		for(RefSeqGene gene: genes){
			if(counter%100 ==0){System.err.println("Random "+counter+" "+genes.size());}
			IntervalTree<Alignment> tree=data.getIntervalTree(gene.getChr(), gene.getStart(), gene.getEnd());
			GeneComponent gc=new GeneComponent(gene);
			Collection<RefSeqGene> windows=gc.getGene().getWindows(3);
			gc.setGeneWindowScores(model.scoreGenes(windows, gene.getChr(), tree), 3);
			Collection<GeneScore> regions=gc.getGeneWindowScores().values();
			regions=filter(regions, peaks.get(gene));
			rtrn.put(gene, regions);
			counter++;
		}
		
		
		return rtrn;
	}

	private Collection<GeneScore> filter(Collection<GeneScore> regions,	Collection<RefSeqGene> peaks) {
		Collection<GeneScore> rtrn=new ArrayList<GeneScore>();
		
		IntervalTree<RefSeqGene> tree=new IntervalTree<RefSeqGene>();
		for(RefSeqGene peak: peaks){
			tree.put(peak.getStart(), peak.getEnd(), peak);
		}
		
		for(GeneScore window: regions){
			Iterator<Node<RefSeqGene>> iter=tree.overlappers(window.getGene().getStart(), window.getGene().getEnd());
			if(!iter.hasNext()){rtrn.add(window);}
		}
		return rtrn;
	}

	/*private Map<RefSeqGene, Double> score(Map<RefSeqGene, RefSeqGene> genesToStart, AlignmentDataModelStats data, ContinuousDataAlignmentModel model) throws IOException {
		Map<RefSeqGene, Double> rtrn=new TreeMap<RefSeqGene, Double>();
		this.expression=new TreeMap<RefSeqGene, GeneScore>();
		
		for(RefSeqGene gene: genesToStart.keySet()){
			IntervalTree<Alignment> tree=data.getIntervalTree(gene.getChr(), gene.getStart(), gene.getEnd());
			RefSeqGene peak=genesToStart.get(gene);
			
			GeneScore peakScore=new GeneScore(peak, model.scoreGene(peak, tree));
			GeneScore fullScore=new GeneScore(gene, model.scoreGene(gene, tree));
			
			double normScore=normalize(peakScore, fullScore);
			rtrn.put(peak, normScore);
			this.expression.put(gene, fullScore);
		}
		
		return rtrn;
	}*/

	/*private void write(Map<RefSeqGene, Double> scores, Map<RefSeqGene, Double> expressionScores, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(RefSeqGene gene: this.geneToStart.keySet()){
			RefSeqGene peak=this.geneToStart.get(gene);
			double score=scores.get(peak);
			double expressionScore=expressionScores.get(peak);
			GeneScore coverage=this.coverage.get(gene);
			GeneScore expression=this.expression.get(gene);
			writer.write(peak+"\t"+score+"\t"+expressionScore+"\t"+expression.getScanPvalue()+"\t"+coverage.getScanPvalue()+"\n");
		}
		
		writer.close();
	}*/

	private void write(Map<RefSeqGene, Collection<GeneScore>> scores, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(RefSeqGene gene: scores.keySet()){
			System.err.println(gene.getName());
			Collection<GeneScore> peaks=scores.get(gene);
			GeneScore coverage=this.coverage.get(gene);
			for(GeneScore peak: peaks){
				double norm=this.normalize(peak, coverage);
				double p=peak.getScanPvalue();
				if(norm>0 && p<alpha){
					RefSeqGene g=peak.getGene();
					g.setName("s="+norm);
					writer.write(g+"\n");
				}
			}
		}
		
		writer.close();
	}
	

	
	private void write(Map<RefSeqGene, Collection<GeneScore>> scores, Map<RefSeqGene, Collection<GeneScore>> randomScores, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("chr\tstart\tend\tratio\tStart Codon (Norm)\tMax non-AUG (Norm)\tStart Codon (Avg)\tMax non-AUG(Avg)\tSignificance\n");
		
		for(RefSeqGene gene: scores.keySet()){
			System.err.println(gene.getName());
			Collection<GeneScore> peaks=scores.get(gene);
			GeneScore random=max(randomScores.get(gene));
			GeneScore peak=max(peaks);
			if(peak!=null){
				GeneScore coverage=this.coverage.get(gene);
				double norm=this.normalize(peak, coverage);
				double randomNorm=this.normalize(random, coverage);
				double ratio=norm/randomNorm;
				if(norm>0){
					writer.write(peak.getGene().getChr()+"\t"+peak.getGene().getStart()+"\t"+peak.getGene().getEnd()+"\t"+ratio+"\t"+norm+"\t"+randomNorm+"\t"+peak.getAvergeNumberOfReads()+"\t"+random.getAvergeNumberOfReads()+"\t"+coverage.getScanPvalue()+"\n");
				}
			}
		}
		
		writer.close();
	}
	
	private GeneScore max(Collection<GeneScore> scores) {
		GeneScore rtrn=null;
		
		for(GeneScore score: scores){
			if(rtrn==null){rtrn=score;}
			else{
				if(score.getAvergeNumberOfReads()>rtrn.getAvergeNumberOfReads()){rtrn=score;}
			}
		}
		
		return rtrn;
	}

	private Map<RefSeqGene, Collection<GeneScore>> scoreStartCodon(Collection<RefSeqGene> genes, AlignmentDataModelStats data, ContinuousDataAlignmentModel model, Sequence chrSequence) throws IOException {
		this.geneToStart=new TreeMap<RefSeqGene, RefSeqGene>();
		this.coverage=new TreeMap<RefSeqGene, GeneScore>();
		
		Map<RefSeqGene, Collection<GeneScore>> rtrn=new TreeMap<RefSeqGene, Collection<GeneScore>>();
		int counter=0;
		for(RefSeqGene gene: genes){
			if(counter%100 ==0){System.err.println(counter+" "+genes.size());}
			IntervalTree<Alignment> tree=data.getIntervalTree(gene.getChr(), gene.getStart(), gene.getEnd());
			GeneScore fullScore=new GeneScore(gene, model.scoreGene(gene, tree));
			this.coverage.put(gene, fullScore);
			Collection<RefSeqGene> startCodons=new TreeSet<RefSeqGene>();
			Collection<RefSeqGene> orfs=FindORFs.findAllORFs(gene, chrSequence, false);
			for(RefSeqGene orf: orfs){
				RefSeqGene start=orf.getStartCodon();
				startCodons.add(start);
			}
			Collection<GeneScore> peakScores=new ArrayList<GeneScore>();
			for(RefSeqGene startCodon: startCodons){
				//TODO This is where the work needs to be done
				//GeneScore startCodonScore=model.getPeak(gene, startCodon, tree);
				GeneScore startCodonScore=new GeneScore(startCodon, model.scoreGene(startCodon, tree));
				this.geneToStart.put(gene, startCodonScore.getGene());
				peakScores.add(startCodonScore);
			}
			
			rtrn.put(gene, peakScores);
			counter++;
		}
		return rtrn;
	}

	private double normalize(GeneScore startCodonScore, GeneScore fullScore) {
		double avg=startCodonScore.getNumberOfReads()/startCodonScore.getGene().getTranscriptLength();
		return avg/fullScore.getAvergeNumberOfReads();
	}

	/*private double normalize(GeneScore startCodonScore,	Map<RefSeqGene, double[]> randomScores) {
		//for now lets score by median
		double[] scores=new double[randomScores.size()];
		
		int i=0;
		for(RefSeqGene gene: randomScores.keySet()){
			GeneScore r1=new GeneScore(gene, randomScores.get(gene));
			scores[i]=r1.getNumberOfReads();
			i++;
		}
		
		double median=Statistics.median(scores);
		return startCodonScore.getNumberOfReads()/median;
	}*/
	

	private Map<RefSeqGene, double[]> computeRandomScores(RefSeqGene gene,IntervalTree<Alignment> tree, ContinuousDataAlignmentModel model) throws IOException {
		//Basically go through each codon and compute the overlapping scores, we'll then compute the median as a norm factor
		Collection<RefSeqGene> windows=gene.getCDS().getWindows(3);
		Map<RefSeqGene, double[]> scores=model.scoreGenes(windows, gene.getChr(), tree);
		return scores;
	}

	public static void main(String[] args)throws IOException{
		if(args.length>4){
			String bamFile=args[0];
			Map<String, Collection<RefSeqGene>> genesByChr=BEDFileParser.loadDataByChr(new File(args[1]));
			String genomeDir=args[2];
			String save=args[3];
			String chr=args[4];
			String sizeFile=genomeDir+"/sizes";
			
			Sequence chrSeq=SequenceUtils.getChrSequence(genomeDir, chr); 
			
			new StartCodonAnalysis(bamFile, sizeFile, genesByChr.get(chr), save, chr, chrSeq);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam file \n args[1]=genes \n args[2]=genome directory \n args[3]=save \n args[4]=chr";
	
}
