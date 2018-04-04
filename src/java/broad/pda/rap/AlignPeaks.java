package broad.pda.rap;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import broad.core.math.Statistics;
import broad.core.motif.SearchException;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;
import broad.pda.chromosome.Chromosome;
import broad.pda.datastructures.Alignments;
import jaligner.Alignment;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixGenerator;

public class AlignPeaks {
	private static final double UNALIGNABLE = -99;
	Matrix matrix=MatrixGenerator.generate(1.0f, -1.0f);
	int numPerm=100;
	
	//Take the peaks called and align directly to mRNA sequence
	public AlignPeaks(Collection<Alignments> regions, Collection<Sequence> referenceSeq, String genomeDirectory,String outFile) throws Exception{
		//Read the sizes file
		Map<String, Integer> chrSizes=BEDFileParser.loadChrSizes(genomeDirectory+"/sizes");
		//Align the regions to the ref seq 
		Map<Alignments, Double> regionScores=align(regions, referenceSeq, genomeDirectory,outFile+".log");
		write(outFile, regionScores);
		//permute and normalize
		//Map<Alignments, PermScore> normalized=normalizeScores(regionScores, numPerm, chrSizes, referenceSeq, genomeDirectory);
		
	}

	private Map<Alignments, PermScore> normalizeScores(Map<Alignments, Double> regionScores, int numPerm, Map<String, Integer> chrSizes, Collection<Sequence> referenceSeq, String genomeDirectory) throws Exception {
		Map<Alignments, PermScore> rtrn=new TreeMap<Alignments, PermScore>();
		
		for(Alignments region: regionScores.keySet()){
			PermScore perm=new PermScore(regionScores.get(region));
			for(int i=0; i<numPerm; i++){
				Alignments permRegion=BEDFileParser.permute(region, chrSizes);
				double score=align(permRegion, referenceSeq, genomeDirectory, "log");
				perm.addPermutationScore(score);
			}
			System.err.println(region.toUCSC()+" "+perm.observed+" "+perm.normalizedScore()+" "+perm.zScore());
			rtrn.put(region, perm);
		}
		
		return rtrn;
	}

	private void write(String save, Map<Alignments, Double> regionScores) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Alignments align: regionScores.keySet()){
			writer.write(align.getChr()+"\t"+align.getStart()+"\t"+align.getEnd()+"\t"+regionScores.get(align)+"\n");
		}
		
		writer.close();
	}


	private Map<Alignments, Double> align(Collection<Alignments> regions, Collection<Sequence> referenceSeq, String genomeDirectory,String outFile) throws Exception {
		Map<Alignments, Double> rtrn=new TreeMap<Alignments, Double>();
		int counter=0;
		//For each region in bed file
		for (Alignments region: regions) {
			//Align to the reference
			double score=align(region, referenceSeq, genomeDirectory,outFile);
			rtrn.put(region, score);
			counter++;
			//double score=align(region, referenceSeq, genomeDirectory);
			System.err.println(counter+" "+regions.size());
		}
		return rtrn;
	}
	
	private double align(Alignments region, Collection<Sequence> referenceSeq, String genomeDirectory,String outFile) throws Exception {
		String chr=null;
		Chromosome chrom=null;
		
		if(chrom==null || !chr.equalsIgnoreCase(region.getChr())){
			chr=region.getChr();
			String sequenceFile=genomeDirectory+"/"+region.getChr().replaceAll("chr", "").trim()+"/"+region.getChr()+".agp";
			chrom = new Chromosome(sequenceFile);
			chrom.loadSequence();
		}
		
		//first get the sequence for the window
		String seq=region.getSequence(chrom, false);
		Sequence sequence=new Sequence(region.getName());
		sequence.setSequenceBases(seq);
		
		//align the sequence to the reference
		double score=localAlignment(referenceSeq, sequence,outFile);
		//Map<Alignment,Double> score = localAlignment(referenceSeq, sequence);
		return score;
		
	}
		
	
	private double localAlignment(Collection<Sequence> referenceSeqs, Sequence seq,String outFile) throws SearchException, IOException {
		//if seq is NNNN skip and return null;
		//Map<Alignment,Double> rtrn = new HashMap<Alignment,Double>();
		BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));

		boolean allN=true;
		char[] chars=seq.getSequenceBases().toCharArray();
		for(int i=0; i<chars.length; i++){
			if(chars[i]!='N'){allN=false; break;}
		}
		
		if(allN){
			//rtrn.put(null, UNALIGNABLE);
			return UNALIGNABLE;}
				
		//Lets try to seed again
		/*List<String> seeds=this.makeKmers(seq, kmer);
		boolean hasSeed=this.seed(seeds, referenceSeq);
		if(!hasSeed){return NO_SEED;}*/
		List<Double> scores=new ArrayList<Double>();
		double max = Double.MIN_VALUE;
		for(Sequence referenceSeq: referenceSeqs){
			//align left and rightRev
			jaligner.Sequence refSeq=new jaligner.Sequence(referenceSeq.getSequenceBases().toUpperCase());
			jaligner.Sequence windowSeq=new jaligner.Sequence(seq.getSequenceBases().toUpperCase());
			jaligner.Sequence windowASSeq=new jaligner.Sequence(Sequence.reverseSequence(seq.getSequenceBases()).toUpperCase());
					
			//Because it can hybridize to either DNA strand
			Alignment align=SmithWatermanGotoh.align(refSeq, windowSeq, matrix, 2, 1);
			Alignment alignAS=SmithWatermanGotoh.align(refSeq, windowASSeq, matrix, 2, 1);
					
			double score= Math.max(align.calculateScore(), alignAS.calculateScore());
			if(score>max){
			//	rtrn = new HashMap<Alignment,Double>();
				if(score==align.calculateScore())
					bw.write(align.getSummary()+"\n");
				else
					bw.write(alignAS.getSummary()+"\n");
				max= score;
			}
			//double score= Math.max(align.getNumberOfMatches2(), alignAS.getNumberOfMatches2());
			scores.add(score);
		}
		
		bw.close();
		//logger.info(refSeq.getSequence()+"\t"+windowSeq.getSequence()+"\t"+score);
		//Returning the max scored alignment
		return Statistics.max(scores);
		//return rtrn;
	}
	
	private class PermScore{
		double observed;
		ArrayList<Double> permutedValues;
		
		PermScore(double observed){
			this.observed=observed;
			this.permutedValues=new ArrayList<Double>();
		}
		
		void addPermutationScore(double val){permutedValues.add(val);}
		
		double normalizedScore(){
			return observed-Statistics.median(permutedValues);
		}
		
		double zScore(){
			return Statistics.zScore(observed, permutedValues);
		}
	}
	
	public static void main(String[] args)throws Exception{
		if(args.length>3){
			//Read peaks
			Collection<Alignments> peaks=BEDFileParser.loadAlignmentData(new File(args[0]));
			//Read the fasta sequences
			FastaSequenceIO geneFasta=new FastaSequenceIO(args[1]);
			List<Sequence> geneSequences=geneFasta.loadAll();
			String genomeDir=args[2];
			String outFile = args[3];
			new AlignPeaks(peaks, geneSequences, genomeDir,outFile);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=peaks \n args[1]=fasta file \n args[2]=genome directory \n args[3] = output file name";
	
}
