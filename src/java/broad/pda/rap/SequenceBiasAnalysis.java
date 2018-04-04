package broad.pda.rap;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import Jama.Matrix;
import broad.core.sequence.SequenceUtils;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GeneScore;

public class SequenceBiasAnalysis {
	
	double alignmentScoreCutoff;
	int extension=300;

	//Take the regions and permute them to compute the enrichment
	public SequenceBiasAnalysis(ContinuousDataAlignmentModel model, String alignmentScoreFile, int numPerm, String sizeFile, String save, double alignmentScoreCutoff) throws IOException{
		Map<String, Integer> chrSizes=BEDFileParser.loadChrSizes(sizeFile);
		
		this.alignmentScoreCutoff=alignmentScoreCutoff;
		
		//Take BEDGraph of sequence alignments and make features
		Collection<Alignments> features=makeSequenceRegions(alignmentScoreFile);
		BEDFileParser.writeBED("test.bed", features);
		
		//Score features
		List<GeneScore> scores=scoreFeatures(model, features);
		System.err.println(scores.size()+" elements passed alignment score "+alignmentScoreCutoff);
		//write("scores.bed", scores);
		
		//permute and score
		List<GeneScore>[] permutedScores=permute(model, features, numPerm, chrSizes);
		
		//print results
		write(scores, permutedScores, save);
	}
	
	private void write(String save, List<GeneScore> scores) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(GeneScore score: scores){
			writer.write(score.getGene().getChr()+"\t"+score.getGene().getStart()+"\t"+score.getGene().getEnd()+"\t"+score.getNumberOfReads()+"\n");
		}
		
		writer.close();
	}

	private Collection<Alignments> makeSequenceRegions(String alignmentScoreFile) throws IOException {
		Collection<String> list=BEDFileParser.loadList(alignmentScoreFile);
		
		Collection<Alignments> alignmentScores=new TreeSet<Alignments>();
		
		for(String line: list){
			String[] tokens=line.split("\t");
			Alignments region=new Alignments(tokens[0], new Integer(tokens[1]), new Integer(tokens[2]));
			double score=new Double(tokens[3]);
			if(score>this.alignmentScoreCutoff){
				Alignments extended=extend(region, this.extension);
				alignmentScores.add(extended);
			}
		}
		
		//TODO Collapse alignments
		return alignmentScores;
	}

	private Alignments extend(Alignments region, int extension) {
		return new Alignments(region.getChr(), region.getStart()-extension, region.getEnd()+extension);
	}

	//This will just write a row for each
	private void write(Collection<GeneScore> scores, Collection<GeneScore>[] permutedScores, String save) throws IOException {
		Matrix matrix=new Matrix(scores.size(), permutedScores.length+1);
		
		//set the first column to scores
		int row=0;
		for(GeneScore score: scores){
			matrix.set(row, 0, score.getEnrichment());
			row++;
		}
		
		for(int column=0; column<permutedScores.length; column++){
			row=0;
			for(GeneScore score: permutedScores[column]){
				matrix.set(row, column+1, score.getEnrichment());
				row++;
			}
		}
		
		matrix.write(save);
	}

	private List<GeneScore>[] permute(ContinuousDataAlignmentModel model,Collection<Alignments> features, int numPerm, Map<String, Integer> chrSizes) throws IOException {
		List<GeneScore>[] rtrn=new List[numPerm];
		for(int i=0; i<numPerm; i++){
			Collection<Alignments> permutedFeatures=BEDFileParser.permute(features, chrSizes);
			rtrn[i]=scoreFeatures(model, permutedFeatures);
			//write("perm"+i+".bed", rtrn[i]);
		}
		return rtrn;
	}

	private List<GeneScore> scoreFeatures(ContinuousDataAlignmentModel model, Collection<Alignments> features) throws IOException{
		List<GeneScore> scores=new ArrayList<GeneScore>();
		
		for(Alignments align: features){
			RefSeqGene gene=new RefSeqGene(align);
			GeneScore score=new GeneScore(gene, model.scoreGene(gene));
			scores.add(score);
		}
		
		return scores;
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>5){
			String BAMFile=args[0];
			String sizeFile=args[1];
			String bedGraph=args[2];
			int numPerm=new Integer(args[3]);
			String save=args[4];
			double alignmentScoreCutoff=new Double(args[5]);
			ContinuousDataAlignmentModel dataModel=SequenceUtils.getDataModel(BAMFile, sizeFile, false);
			new SequenceBiasAnalysis(dataModel, bedGraph, numPerm, sizeFile, save, alignmentScoreCutoff);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=BAM file \n args[1]=sizeFile \n args[2]=bedGraph \n args[3]=numPerm \n args[4]=save \n args[5]=alignment score cutoff";
	
}
