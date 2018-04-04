package broad.pda.seq.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.util.CollapseByIntersection;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.graph.Path;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

public class TFsWithinGene {
	double alpha=0.05;
	//int[] windowSizes={25,50,100,150,200,300,400,500,600,700,800,900,1000,1500,2000,2500,3000};
	int extensionLength=10000;
	int[] windowSizes={100};
	
	public TFsWithinGene(AlignmentDataModel alignmentModel, Map<String, RefSeqGene> genes, String save, int extensionFactor, Collection<String> queryGenes) throws IOException{
		Map<String, IntervalTree<RefSeqGene>> trees=CollapseByIntersection.makeIntervalTreeForGenes(genes.values());
		
		this.extensionLength=extensionFactor;
		ContinuousDataAlignmentModel dataModel=new ContinuousDataAlignmentModel(alignmentModel);
		Map<String, Collection<Alignments>> segments=getSegments(alignmentModel, dataModel, this.windowSizes, genes, trees, queryGenes);
		
		//Now that we have all the segments find all that are within the gene +/- extension factor
		//get only those regions queried
		
		
		writePeaks(save, segments);
	}

	private void writePeaks(String save, Map<String, Collection<Alignments>> segments) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String name: segments.keySet()){
			Collection<Alignments> peaks=segments.get(name);
			for(Alignments region: peaks){writer.write(region+"\t"+name+"\n");}
		}
		
		writer.close();
	}

	private Map<String, Collection<Alignments>> getSegments(AlignmentDataModel alignmentModel, ContinuousDataAlignmentModel dataModel, int[] windowSizes, Map<String, RefSeqGene> genes, Map<String, IntervalTree<RefSeqGene>> trees, Collection<String> queryGenes) throws IOException {
		Map<String, Collection<Alignments>> rtrn=new TreeMap<String, Collection<Alignments>>();
		for(String name: queryGenes){
			System.err.println(name);
			RefSeqGene gene=genes.get(name.toUpperCase().trim());
			//System.err.println(name+" "+gene);
			//Alignments extended=getExtension(gene, trees);
			Alignments extended=extend5Primer(gene);
			for(int i=0; i<windowSizes.length; i++){
				Map<Path, double[]> segments=dataModel.scanFromGraph(windowSizes[i], alpha, true, extended.getChr(), extended.getStart(), extended.getEnd(), null, false, 0.1);
				rtrn.put(gene.getName(), refToAlign(segments.keySet()));
			}
		}
		
		return rtrn;
		
	}
	
	private Alignments extend5Primer(RefSeqGene gene) {
		//System.err.println(gene);
		if(gene.getOrientation().equalsIgnoreCase("+")){
			Alignments align=new Alignments(gene.getChr(), gene.getStart()-this.extensionLength, gene.getStart()+this.extensionLength);
			return align;
		}
		else{
			Alignments align=new Alignments(gene.getChr(), gene.getEnd()-this.extensionLength, gene.getEnd()+this.extensionLength);
			return align;
		}
	}

	private Alignments getExtension(RefSeqGene gene, Map<String, IntervalTree<RefSeqGene>> trees) {
		Alignments extended=new Alignments(gene.getChr(), gene.getStart()-this.extensionLength, gene.getEnd()+this.extensionLength);
		
		//check if has overlapping gene that isnt self
		Iterator<Node<RefSeqGene>> overlappers=trees.get(extended.getChr()).overlappers(extended.getStart(), extended.getEnd());
		
		int start=Integer.MAX_VALUE;
		int end=-Integer.MAX_VALUE;
		int counter=0;
		while(overlappers.hasNext()){
			RefSeqGene overlapping=overlappers.next().getValue();
			counter++;
			start=Math.min(overlapping.getEnd(), start);
			end=Math.max(overlapping.getStart(), end);
		}
		
		if(counter<=1){return extended;}
		else{return new Alignments(extended.getChr(), start, end);}
		//if not then return +/- extension
		//else return from end of previous gene to start of next gene
	}

	private Collection<Alignments> refToAlign(Collection<Path> genes) {
		Collection<Alignments> rtrn=new TreeSet<Alignments>();
		
		for(Path gene: genes){rtrn.add(gene.toGene().getAlignment());}
		
		return rtrn;
	}

	public static void main(String[] args)throws IOException{
		if(args.length>4){
			AlignmentDataModel alignments=new GenericAlignmentDataModel(args[0], args[1]);
			Map<String, RefSeqGene> genes=BEDFileParser.loadDataByName(new File(args[2]));
			String save=args[4];
			int extensionFactor=50000;
			Collection<String> queryGenes=parseNames(args[3]);
			if(args.length>5){extensionFactor=new Integer(args[5]);}
			new TFsWithinGene(alignments, genes, save, extensionFactor, queryGenes);
		}
		else{System.err.println(usage);}
	}
	
	private static Collection<String> parseNames(String file) throws IOException {
		Collection<String> rtrn=new TreeSet<String>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine; 
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {rtrn.add(nextLine.split("\t")[0]);}
		reader.close();
		return rtrn;
	}

	static String usage=" args[0]=alignment file (SAM, BAM, aligned) \n args[1]=genome sizes \n args[2]=genes \n args[3]=query genes \n args[4]=save \n args[5]=max distance around gene (default 50,000)";
	
}
