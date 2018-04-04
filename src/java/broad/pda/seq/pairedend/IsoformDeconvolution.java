package broad.pda.seq.pairedend;

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
import broad.core.datastructures.IntervalTree.Node;
import broad.core.math.EmpiricalDistribution;
import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.graph.Path;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;
import net.sf.samtools.util.CloseableIterator;


//Take all possible isoforms determined by segmentations and score each one based on probability

public class IsoformDeconvolution {
	
	public IsoformDeconvolution(Collection<RefSeqGene> transcripts, AlignmentDataModel pairs, String save) throws IOException{
		EstimatePairedEndDistribution estimateDist=new EstimatePairedEndDistribution(pairs, transcripts);
		EmpiricalDistribution pairedDist=estimateDist.getSizeDistribution();
		
		System.err.println(pairedDist.getProbability(pairedDist.getMean()));
		
		//assign probabilities to each transcript
		Map<Path, Double> probabilities=assignProbabilities(transcripts, pairs, pairedDist);
		
		Map<Path, Double> weights=assignWeights(probabilities);
		
		//TODO: Weight isoforms
		write(save, weights);
	}

	private Map<Path, Double> assignWeights(Map<Path, Double> probabilities) {
		Map<Path, Double> rtrn=new TreeMap();
		
		Map<String, IntervalTree<Path>> pathTree=TreeUtils.makeIntervalTreeByPath(probabilities.keySet());
		
		
		for(Path path: probabilities.keySet()){
			Iterator<Node<Path>> overlappers=pathTree.get(path.getChromosome()).overlappers(path.getStart(), path.getEnd());
			double relative=assignRelativeWeight(path, overlappers, probabilities);
			rtrn.put(path, relative);
		}
		
		return rtrn;
	}

	private double assignRelativeWeight(Path path,Iterator<Node<Path>> overlappers, Map<Path, Double> probabilities) {
		double value=probabilities.get(path);
		
		double total=0;
		while(overlappers.hasNext()){
			Path path2=overlappers.next().getValue();
			total+=probabilities.get(path2);
		}
		
		return value/total;
	}

	private Map<Path, Double> assignProbabilities(Collection<RefSeqGene> transcripts, AlignmentDataModel pairs, EmpiricalDistribution pairedDist) throws IOException {
		Collection<Path> paths=new TreeSet<Path>();
		
		for(RefSeqGene transcript: transcripts){
			if(transcript.getNumExons()>1){
				Path path=new Path(transcript, null); //TODO Replace this with the parent graph
				System.err.println(path.toGene());
				CloseableIterator<Alignment> iter= pairs.getAlignmentsOverlappingRegion(transcript.getAlignment());
				while(iter.hasNext()){
					Alignment align=iter.next();
					path.addPairedEnd(align);
				}
				iter.close();
				paths.add(path);
			}
		}
		
		return assignProbabilities(paths, pairedDist);
	}
	
	private Map<Path, Double> assignProbabilities(Collection<Path> paths, EmpiricalDistribution pairedDist) {
		Map<Path, Double> rtrn=new TreeMap<Path, Double>();
		
		for(Path path: paths){
			Collection<Alignment> pairs=path.getPairedEndEdges();
			double score=scoreEachRead(path, pairs, pairedDist);
			rtrn.put(path, score);
		}
		
		return rtrn;
	}

	

	private double scoreEachRead(Path path, Collection<Alignment> pairs, EmpiricalDistribution pairedDist) {
		double score=0;
		double count=0;
		for(Alignment read: pairs){
			double size=estimateSize(read, path.toGene());
			//double cdf=pairedDist.getCummulativeProbability(size);
			double prob=1-pairedDist.getProbability(size);
			score+=prob;
			count++;
		}
		return score;
	}

	private double estimateSize(Alignment read, RefSeqGene transcript) {
		double distance=-1;
		RefSeqGene insert=transcript.trimAbsolute(read.getAlignmentStart(), read.getAlignmentEnd());
		if(insert==null){distance=-1;}
		else{distance=insert.getSize();}
		return distance;
	}

	private void write(String save, Map<Path, Double> probabilities) throws IOException {
		FileWriter writer=new FileWriter(save);
	
		for(Path gene: probabilities.keySet()){
			writer.write(gene.toGene()+"\t"+probabilities.get(gene)+"\n");
		}
		
		writer.close();
	}

	public static void main(String[] args)throws IOException{
		if(args.length>3){
			Collection<RefSeqGene> transcripts=BEDFileParser.loadData(new File(args[0]));
			AlignmentDataModel pairs=new GenericAlignmentDataModel(args[1], args[2]);
			String save=args[3];
			new IsoformDeconvolution(transcripts, pairs, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=transcripts \n args[1]=pair slignments \n args[2]=sizes \n args[3]=save";
	
}
