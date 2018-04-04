package broad.pda.seq.segmentation;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.AlignmentBlock;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;

public class TranscriptUniformityScore {

	private static int DEFAULT_MIN_MAPPING_QUALITY = 5;
	private static int extensionFactor=0;
	AlignmentDataModel data;
	
	
	public TranscriptUniformityScore(GenericAlignmentDataModel newdata){
		this.data=newdata;
	}
	
	public TranscriptUniformityScore(AlignmentDataModel newdata){
		this.data=newdata;
	}
	
	public static  void calcUniformityScore (String alignmentFile,boolean useConstituentExons,
			double minMappingQuality ,BEDFileParser annotationParser,String save,
			String sizes,File[] maskFiles , Map<String, Integer> maskFileData, boolean isStranded) throws IOException {
		
		Map<String, Collection<RefSeqGene>> annotations = null;
		if(useConstituentExons) {
			annotationParser.makeGenes();
			//annotationParser.writeFullBed("Madegenes.bed");
			annotations = annotationParser.toConstituentIsoformMap();
		} else {
			annotations = annotationParser.toMap();
		}
		Collection<RefSeqGene> annotationCollection = new ArrayList<RefSeqGene>();
		for(Collection<RefSeqGene> chrAnnotations : annotations.values()) {
			annotationCollection.addAll(chrAnnotations);
		}
		//BEDFileParser.writeFullBED("constituentIsoforms.bed", annotationCollection);
			
		
		if(!isStranded) {
			System.err.println("Scoring using all reads ");
			//Here we initilaize the reads
			AlignmentDataModel alignments=new GenericAlignmentDataModel(alignmentFile, sizes, false, minMappingQuality);
			Map<RefSeqGene, double[]> scores=new TreeMap<RefSeqGene, double[]>();				
			runScore(annotations, save, maskFileData, alignments, scores,extensionFactor);
			writeFullBED(save, scores);
		} else {
			System.err.println("Scoring minus reads");
			AlignmentDataModel alignments=new GenericAlignmentDataModel(alignmentFile, sizes, false, minMappingQuality);
			alignments.setNegativeStranded();
			Map<RefSeqGene, double[]> scores=new TreeMap<RefSeqGene, double[]>();				
			runScore(annotations, save, maskFileData, alignments, scores,extensionFactor);
			writeFullBED(save+".minus", scores);

			System.err.println("Scoring plus reads");
			alignments.setPositiveStranded();
			scores=new TreeMap<RefSeqGene, double[]>();				
			runScore(annotations, save, maskFileData, alignments, scores,extensionFactor);
			writeFullBED(save+".plus", scores);
		}
	}
	
	
	private static void writeFullBED(String save,
			Map<RefSeqGene, double[]> scores) throws IOException {
		
		FileWriter writer=new FileWriter(save);
		for(RefSeqGene align: scores.keySet()){
			double[] scr=scores.get(align);
			writer.write(align.toBED());
			for(double p : scr) 
				writer.write("\t"+p);
			writer.write("\n");
		}
		writer.close();
	}


	public static  void runScore(Map<String, Collection<RefSeqGene>> annotations, String save,
			Map<String, Integer> maskFileData, AlignmentDataModel alignments,	Map<RefSeqGene, double[]> scores,int extensionFactor) throws IOException {
		AlignmentDataModelStats alignmentData = new AlignmentDataModelStats(alignments, maskFileData);
		ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(alignmentData, maskFileData, 0, 1, null, null);
		TranscriptUniformityScore uniforScr= new TranscriptUniformityScore(alignmentData.getData());
		for(String chr : annotations.keySet()) {
			System.err.println("processig " + chr);
			Collection<RefSeqGene> chrAnnotations = annotations.get(chr);
			scores.putAll(uniforScr.scoreGenes( data,chrAnnotations, chr,extensionFactor));
			//System.err.print("done. saving data  to " + save+"."+chr);
		}
	}		
	
	
	public  Map<RefSeqGene, double[]> scoreGenes(ContinuousDataAlignmentModel model,Collection<RefSeqGene> set, String chrToUse,int extensionFactor) throws IOException{
		
		Map<RefSeqGene, double[]> rtrn=new TreeMap<RefSeqGene, double[]>();
		int i=0;
		for(RefSeqGene align: set){
			if(chrToUse==null || chrToUse.equalsIgnoreCase(align.getChr())){
				if(model.getChromosomeLength(align.getChr())!=0){
					IntervalTree<Alignment> tree=model.data.getIntervalTreeCached(align.getChr(), align.getStart()-1, align.getEnd()+1);
					double[] pRate=calcUniformScore(align, tree,extensionFactor);
					align.setBedScore(pRate[2]);//setting uniformity ratio as gene scores
					rtrn.put(align, pRate);
					i++;
				}
			}
		}

		model.data.resetTreeCache();

		return rtrn;
	}

	/**
	 * computes the mean coverage of the gene and the mean exact coverage
	 * of all the introns 
	 * @param model 
	 * @param gene , gene to calculate score for
	 * @param tree, a tree with all the SAM records spanning the gene
	 * @return [a=global coverage, b=mean intron covergae, the ratio b/a ]
	 */
	public double[] calcUniformScore( RefSeqGene gene,
			IntervalTree<Alignment> alignTree,int extensionFactor) {
		
		double[] res = new double[3];
		//System.err.println("\n"+"Gene: " +gene.getName());
		//Row 636 in GenericAligmentDataModel
		double sum=data.getCountsPerAlignment(gene, alignTree, extensionFactor);
		double sum2=data.getBasesCoveredPerAlignment(gene, alignTree, extensionFactor);
		int geneLength=gene.getTranscriptLength();
		double avgCoverage = sum2/(double) geneLength;
		
		double tmp=sum*60;
		//System.err.println("countedPerAlign: " + sum + "after multiplying with read length =60 "+ tmp + "Covered bases: "+sum2 );
		
		//get the coordinates that flank the intron
		Collection<Alignments> introns= gene.getIntronSet();
		IntervalTree <Alignments> intronTree= new IntervalTree<Alignments>();
		HashMap<Alignments,Double> cntMap= new HashMap<Alignments,Double>();
		for (Alignments a: introns){
			intronTree.put(a.getStart(),a.getEnd(),a);
			cntMap.put(a,new Double(0));
			//System.err.println("Intron:"+a.getStart()+"\t"+a.getEnd()+"\n");
		}
		
		
		
		Iterator<Node<Alignment>> it= alignTree.iterator();
		double minReadLength=76;
		int numReads=0;
		while(it.hasNext()){
			
			Node<Alignment> alnNode=it.next();
			int replicates=alnNode.getNumReplicates();
			numReads+=alnNode.getNumReplicates();
			Alignment read=alnNode.getValue();
			//System.err.println("GapRead: "+read.getReadName());
			double coveredBase=0;
			//Get the coordinates of the gapps= introns 
			//LEARN IF THESE START IN THE FIRST BASE OF THE INTRON , OR THE LAST OF THE EXON
			AlignmentBlock[] blocks=read.getAlignmentBlocks();
			
			for (int i=0;i < blocks.length-1 ;i++){
				int iStrt=blocks[i].getEnd();
				int iEnd=blocks[i+1].getStart();
				//System.err.println("GapRead: "+read.getReadName()+"\t"+blocks[i].toString() +iStrt+"\t"+iEnd+"\n");
				Node<Alignments> suportedIntronNode= intronTree.find(iStrt, iEnd);
				if (suportedIntronNode != null){
					Alignments suportedIntron=suportedIntronNode.getValue();
					cntMap.put(suportedIntron,cntMap.get(suportedIntron)+replicates);
				}
				coveredBase += (blocks[i].getEnd()-blocks[i].getStart());
			}
			//if (coveredBase > 0)
				//minReadLength=Math.min(minReadLength, coveredBase);
		}
		
		//System.err.println("Num reads: " + numReads);
	
		
		//Take the max intron count 
		double maxCnt=0.0;
		for (Alignments a: introns){
			//System.err.println("Intron count: "+ cntMap.get(a));
			maxCnt=Math.max(cntMap.get(a),maxCnt);
			//System.err.println("CurrMax: "+ maxCnt);
		}
		
		double val=avgCoverage;
		res[0]=val;
		res[1]=maxCnt;
		if (val >0)
			res[2]=maxCnt/val;
		else
			res[2]=0;
		
		//System.err.println(avgCoverage+"\t"+res[0]+"\t"+res[1]+"\t"+res[2]+"\n");
		
		return res;
	}
	
	
	//get the number of alignments overlapping a given region with a cached interval tree
	/*
	public double getCoveragePerAlignment(Alignments align, IntervalTree<Alignment> tree, int EF){
		double counter=0;
		Iterator<Node<Alignment>> iter=tree.overlappers(align.getStart(), align.getEnd()+1); //TODO Consider getting rid of the plus 1
		
		while(iter.hasNext()){
			Node<Alignment> node=iter.next();
			int num=node.getNumReplicates();
			Alignment alignment=node.getValue();
			int coverdBases=getBaseCoverdPerAlignment(alignment, align, EF);
			
			counter+=coverdBases*num;
		}
		return counter;
	}
	*/
	
	
	
	static String usage="Usage: GeneTools -task <task name> "+
	"\n\nTask: UniformityScore -  Computes 3 expression related scores for a set of annotations [a=coverage,b=max compatible intron reads,b/a ] -in <Full BED file with annotations to score> \n\t-alignment <Alignment file in BAM, SAM or Alignemnt format> \n\t-sizeFile <Chromosome size file> \n\t-out <Output file name> \n\t -maskFileDir <Mask File directory>\n\t -useConstituentExons <For each gene a set of contituent exons will be chosen and scored>"+
	"\n" ;

	public static void main(String [] args) throws Exception  {
		ArgumentMap argmap = CLUtil.getParameters(args, usage, "full");
		
		if("uniformityScore".equalsIgnoreCase(argmap.getTask())) {
			
			String alignmentFile = argmap.getMandatory("alignment");
			boolean useConstituentExons = argmap.containsKey("useConstituentExons");
			double minMappingQuality = argmap.containsKey("minMappingQuality") ? argmap.getDouble("minMappingQuality") : DEFAULT_MIN_MAPPING_QUALITY;
			//Map<String, Collection<RefSeqGene>> annotations = BEDFileParser.loadDataByChr(new File(argmap.getInput()));
			BEDFileParser annotationParser =  new BEDFileParser(argmap.getInput());
			String save = argmap.getOutput();
			String sizes = argmap.getMandatory("sizeFile");
			File[] maskFiles  = argmap.isPresent("maskFileDir") ?  new File(argmap.get("maskFileDir")).listFiles(): null ;
			Map<String, Integer> maskFileData = ContinuousDataAlignmentModel.parseMaskFiles(maskFiles);
			System.err.println("Minimum mapping quality to count reads: " + minMappingQuality);
			boolean isStranded = argmap.containsKey("stranded");
		
			calcUniformityScore (alignmentFile, useConstituentExons, minMappingQuality , annotationParser,
					save,sizes,maskFiles,maskFileData,isStranded );
		}
		else{
			System.err.println(usage);
		}
	}
		
}

