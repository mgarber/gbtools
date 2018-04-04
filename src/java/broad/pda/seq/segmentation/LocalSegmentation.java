package broad.pda.seq.segmentation;

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
import broad.core.util.CollapseByIntersection;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;

//Local Segmentation really should be done after updating by paired ends to get as best an estimate as possible of the true start and end of the graph
//TODO Also should really rescan using the local rate since there are clear exon strucutures that are blurred when junctions arent present.
//TODO Steps first filter then segment (Example Cox5b)

public class LocalSegmentation {

	public LocalSegmentation(Collection<RefSeqGene> genes, AlignmentDataModel data, String save, double alpha) throws IOException{
		Map<RefSeqGene, double[]> localScores=localSegment(genes, data);
		write(save, localScores, alpha);
	}
	
	private Map<RefSeqGene, double[]> localSegment(Collection<RefSeqGene> genes, AlignmentDataModel data) throws IOException{
		Map<RefSeqGene, double[]> rtrn=new TreeMap();
		
		ContinuousDataAlignmentModel model=new ContinuousDataAlignmentModel(data);
		Map<String, IntervalTree<RefSeqGene>> geneTree=CollapseByIntersection.makeIntervalTreeForGenes(genes);
		
		int i=0;
		for(RefSeqGene gene: genes){
			//get all overlapping genes
			Iterator<Node<RefSeqGene>> bubble=geneTree.get(gene.getChr()).overlappers(gene.getStart(), gene.getEnd());
			//estimate lambda
			double lambda=estimateLambda(bubble, data, model, gene.getChr());
			bubble=geneTree.get(gene.getChr()).overlappers(gene.getStart(), gene.getEnd());
			Map<RefSeqGene, double[]> temp=localScore(geneTree, gene, lambda, data);
			rtrn.putAll(temp);
			i++;
			if(i%10 ==0){System.err.println(i+" "+gene.getAlignment().toUCSC()+" "+lambda+" "+model.getLambda(gene.getChr()));}
		}
		
		return rtrn;
	}
	
	
	//estimate lambda from introns and single exons in all
	private double estimateLambda(Iterator<Node<RefSeqGene>> bubble, AlignmentDataModel data, ContinuousDataAlignmentModel model, String chr)   throws IOException{
		/*Collection<Alignments> exons=new TreeSet();
		
		int i=0;
		while(bubble.hasNext()){
			RefSeqGene gene=bubble.next().getValue();
			if(gene.getNumExons()>1){
				exons.addAll(gene.getExonSet());
			}
			i++;
		}
		
		//System.err.println("Bubble size "+i);
		
		if(exons!=null && !exons.isEmpty()){
		exons=CollapseByIntersection.CollapseByIntersection(exons, false);
		RefSeqGene gene=new RefSeqGene(exons);
		
		Collection<Alignments> introns=gene.getIntronSet();
		
		double sum=0;
		double count=0;
		IntervalTree<Alignment> tree=data.getIntervalTreeCached(gene.getChr(), gene.getStart(), gene.getEnd());
		
		for(Alignments intron: introns){
			sum+=data.getCountsPerAlignment(intron, tree, 0);
			count+=intron.getSize();
		}
		double lambda=Math.max(sum/count, model.getLambda(chr));
		return lambda;
		}
		else{
			return model.getLambda(chr);
		}*/
		
		//Attempt 2: Just take the whole region and compute lambda
		// compute exon to intron ratio
		
		Collection<Alignments> exons=new TreeSet();
		while(bubble.hasNext()){
			RefSeqGene gene=bubble.next().getValue();
			//if(gene.getNumExons()>1){
				exons.addAll(gene.getExonSet());
			//}
		}
		
		//System.err.println("Bubble size "+i);
		
		if(exons!=null && !exons.isEmpty()){
			exons=CollapseByIntersection.collapseByIntersection(exons, false);
			RefSeqGene gene=new RefSeqGene(exons);
			
			if(gene.getNumExons()==1){return model.getLambda(chr);}
			double ratio=sum(gene.getExonSet())/sum(gene.getIntronSet());
			
			IntervalTree<Alignment> tree=data.getIntervalTreeCached(gene.getChr(), gene.getStart(), gene.getEnd());
			double counts=data.getCountsPerAlignment(gene.getAlignment(), tree, 0);
			double lambda=Math.max(counts/gene.getAlignment().getSize(), model.getLambda(chr));
			
			//System.out.println(gene+"\t"+lambda+"\t"+ratio);
			
			return lambda;
		}
		
		return model.getLambda(chr);
	}

	private double sum(Collection<Alignments> collection) {
		double sum=0;
		for(Alignments align: collection){sum+=align.getSize();}
		return sum;
	}

	private Map<RefSeqGene, double[]> localScore(Map<String, IntervalTree<RefSeqGene>> geneTree, RefSeqGene gene, double lambda, AlignmentDataModel data) throws IOException{
		Map<RefSeqGene, double[]> rtrn=new TreeMap<RefSeqGene, double[]>();
		Iterator<Node<RefSeqGene>> bubble=geneTree.get(gene.getChr()).overlappers(gene.getStart(), gene.getEnd());
		Alignments region=sourceToSink(bubble);
		
		if(region==null){return rtrn;}
		bubble=geneTree.get(gene.getChr()).overlappers(gene.getStart(), gene.getEnd());
		IntervalTree<Alignment> tree=data.getIntervalTreeCached(region.getChr(), region.getStart(), region.getEnd());
		while(bubble.hasNext()){
			RefSeqGene path=bubble.next().getValue();
			//System.err.println(path);
			double[] scanP=scanPRate(path, data, lambda, region.getSize());
			rtrn.put(path, scanP);
		}
		return rtrn;
	}

	private Map<RefSeqGene, double[]> localScore(Collection<RefSeqGene> paths, AlignmentDataModel data, Alignments region) throws IOException{
		Map<RefSeqGene, double[]> rtrn=new TreeMap<RefSeqGene, double[]>();
		IntervalTree<Alignment> tree=data.getIntervalTreeCached(region.getChr(), region.getStart(), region.getEnd());
		
		//TODO 2 ideas about how to compute lambda
		//1) scale lambda by the exon/intron ratio size
		//2)compute lambda as the rate within the introns (ignoring the exonic counts)
		double lambda=data.getCountsPerAlignment(region, tree, 0)/region.getSize();
		
		//System.err.println(region.toUCSC()+" "+lambda+" "+region.getSize());
		
		for(RefSeqGene gene: paths){
			double[] scanP=scanPRate(gene, data, lambda, region.getSize());
			rtrn.put(gene, scanP);
		}
			
		return rtrn;
	}
	
	public double[] scanPRate(RefSeqGene gene, AlignmentDataModel data, double lambda, double numMarkers) throws IOException{
		String chr=gene.getChr();
		int count=gene.getTranscriptLength();
		
		IntervalTree<Alignment> tree=data.getIntervalTreeCached(chr, gene.getStart(), gene.getEnd());
		double sum=data.getCountsPerAlignment(gene, tree, 0);
		//System.err.println(sum+" "+numMarkers);
		
		double enrich=(sum/count)/lambda;
		double[] rtrn={calculatePVal(new Double(sum).intValue(), lambda, count, numMarkers), enrich};
		return rtrn;
	}
	
	private double calculatePVal(int k, double lambda, double w, double T){
		//if(k<=2){return 1;}
		double lambdaW=lambda*w;
		double a=((k-lambdaW)/k)*(lambda*(T-w)*poisson(k-1, lambdaW));
		double result=Fp(k-1, lambdaW)*Math.exp(-a);
		double p=1-result;
		p=Math.abs(p);
		p=Math.min(1, p);
		//p=Math.max(0, p);
		return p;
	}
	
	private double poisson(int k, double lambda){
		cern.jet.random.Poisson poiss=new cern.jet.random.Poisson(lambda, new cern.jet.random.engine.DRand());
		return poiss.pdf(k);
	}
	
	private double Fp(int k,double lambdaW){
		double sum=0;
		for(int i=0; i<=k; i++){
			sum+=poisson(i, lambdaW);
		}
		return sum;
	}
	
	//requires that paths.size() is at least 1
	private Alignments sourceToSink(Collection<RefSeqGene> paths){
		String chr="";
		int start=Integer.MAX_VALUE;
		int end=-Integer.MAX_VALUE;
		
		for(RefSeqGene path: paths){
			chr=path.getChr();
			start=Math.min(path.getStart(), start);
			end=Math.max(path.getEnd(), end);
		}
		
		return new Alignments(chr, start, end);
	}
	
	
	//requires that paths.size() is at least 1
	private Alignments sourceToSink(Iterator<Node<RefSeqGene>> bubble){
		if(bubble==null || !bubble.hasNext()){return null;}
		
		String chr="";
		int start=Integer.MAX_VALUE;
		int end=-Integer.MAX_VALUE;
		
		while(bubble.hasNext()){
			RefSeqGene path=bubble.next().getValue();
			chr=path.getChr();
			start=Math.min(path.getStart(), start);
			end=Math.max(path.getEnd(), end);
		}
		
		return new Alignments(chr, start, end);
	}
	
	/*private static ChromosomeWithBubbles2 makeGraph(Collection<RefSeqGene> genes, String chr){
		//collapse all exons and  keep all intons
		Collection<Alignments> exons=new TreeSet();
		Collection<Alignments> introns=new TreeSet();
		
		for(RefSeqGene gene: genes){
			exons.addAll(gene.getExonSet());
			introns.addAll(gene.getIntronSet());
		}
		
		exons=CollapseByIntersection.CollapseByIntersection(exons, false);
		
		if(!introns.isEmpty()){
			exons=CollapseByIntersection.DecollapseByIntronLocation(exons, introns);
		}
		
		//populate graph
		ChromosomeWithBubbles2 bubbles=new ChromosomeWithBubbles2(chr, exons, introns);
		return bubbles;
	}*/
	
	
	private void write(String save, Map<RefSeqGene, double[]> localScores, double alpha) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(RefSeqGene gene: localScores.keySet()){
			double[] vals=localScores.get(gene);
			if(vals[0]<alpha){writer.write(gene+"\t"+vals[0]+"\t"+vals[1]+"\n");}
		}
		
		writer.close();
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>4){
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[0]));
			AlignmentDataModel data=new GenericAlignmentDataModel(args[1], args[2]);
			String save=args[3];
			double alpha=new Double(args[4]);
			new LocalSegmentation(genes, data, save, alpha);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=Genes (Full BED) \n args[1]=alignment file (SAM) \n args[2]=sizes \n args[3]=save \n args[4]=alpha";
	
}
