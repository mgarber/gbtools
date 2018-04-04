package broad.pda.ribosome;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.broad.igv.Globals;
import org.broad.igv.sam.Alignment;

import broad.core.datastructures.IntervalTree;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceUtils;
import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;
import broad.pda.ribosome.misc.FindORFs;
import broad.pda.seq.segmentation.AlignmentDataModelStats;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GeneScore;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

public class ComputeRibosomeOccupancyByFeature {

	int windowSize=90;
	public static double ALPHA=0.01;
	
	private ContinuousDataAlignmentModel ribosomeData;
	private ContinuousDataAlignmentModel expressionData;
	
	public ComputeRibosomeOccupancyByFeature(ContinuousDataAlignmentModel ribosomeData, ContinuousDataAlignmentModel expressionData) throws IOException{
			this.ribosomeData=ribosomeData;
			this.expressionData=expressionData;
	}
	
	public ComputeRibosomeOccupancyByFeature(ContinuousDataAlignmentModel ribosomeData) throws IOException{
		this.ribosomeData=ribosomeData;
}
	
	public Collection<GeneComponent> scoreFeatures(Collection<RefSeqGene> genes) throws IOException{
		//Step 1: Score ribosome
		Collection<GeneComponent> scoresByFeature=scoreFeatures(genes, this.ribosomeData, this.ribosomeData.getAlignmentDataModelStats());
		System.err.println("Scored features");
		
		if(expressionData!=null){
			//Step 2: Add expression
			addExpression(scoresByFeature, this.expressionData, this.expressionData.getAlignmentDataModelStats());
			System.err.println("Added Expression");
		}
		
		return scoresByFeature;
	}
	
	
	public ComputeRibosomeOccupancyByFeature(String bamFile, String expressionFile, String sizeFile, Map<String, Collection<RefSeqGene>> genesByChr, String save, String chrToUse) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		Map<String, Integer> chromosomeSizes=BEDFileParser.loadChrSizes(sizeFile);
		
		for(String chr: chromosomeSizes.keySet()){
			if((chrToUse==null || chrToUse.equalsIgnoreCase(chr)) && genesByChr.get(chr)!=null){
			System.err.println(chr+"\t"+genesByChr.get(chr).size());
			Collection<RefSeqGene> genes=genesByChr.get(chr);
		
			//Step 1: Set up data model for ribosome
			Globals.setHeadless(true);
			AlignmentDataModelStats data=new AlignmentDataModelStats(new GenericAlignmentDataModel(bamFile, sizeFile), null, chr);
			ContinuousDataAlignmentModel model=new ContinuousDataAlignmentModel(data);
			
			//Step 2: Score features
			Collection<GeneComponent> scoresByFeature=scoreFeatures(genes, model, data);
			System.err.println("Scored features");
			
			//Step 3: Score ORFs
			/*scoreEachORF(scoresByFeature, data, model, chrSequence);
			System.err.println("Got and scored predicted ORFs");*/
			
			//Step 4: Setup data model for expression
			data=new AlignmentDataModelStats(new GenericAlignmentDataModel(expressionFile, sizeFile), null, chr);
			model=new ContinuousDataAlignmentModel(data);
			
			//Step 5: Add expression to features
			addExpression(scoresByFeature, model, data);
			System.err.println("Added Expression");
			
			//Step 6: Add expression to ORFs
			/*addExpressionToPredictedORFs(scoresByFeature, data, model);
			System.err.println("Added Expression to ORFs");*/
			
			//Step 7: Write features
			writeFeatures(scoresByFeature, writer); //write all components of the gene as part of the same line
			
			//Step 8: Write ORFs
			//writeORFScores(scoresByFeature, save+".maxORF.bed");
			
			//TODO Write coding, noncoding, and intron CDS
			//writePredictedORFs(scoresByFeature, save+".introns");
			
			System.err.println("Finished scoring features");
			}
		}
		
		writer.close();
	}
	
	
	public static void writeRRS(String save, Collection<GeneComponent> genes, boolean started, boolean normalizeByExpression, boolean fullyContainedReadsInCds) throws IOException{
		RibosomeScoring.writeRRS(save, genes, started, ALPHA, normalizeByExpression, fullyContainedReadsInCds);
	}
	
	
	private void writeScatter(String save, Collection<GeneComponent> genes, boolean started) throws IOException {
		FileWriter writer=new FileWriter(save, started);
		
		for(GeneComponent gene: genes){
			double expression=gene.getGeneRNASeqScore().getRPKM();
			double ribosome=gene.getGeneScore().getRPKM();
			writer.write(gene.getName()+"\t"+expression+"\t"+ribosome+"\n");
		}
		
		writer.close();
	}
	
	private void writePredictedORFs(Collection<GeneComponent> scoresByFeature,	String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(GeneComponent gene: scoresByFeature){
			if(gene.hasCDS()){
				writer.write(gene.getName()+"\tCDS\t"+gene.getCDSRNASeqScore().getScanPvalue()+"\t"+this.computeTE(gene.getCDSScore(), gene.getCDSRNASeqScore())+"\n");
			}
			if(gene.hasPredictedORFs()){
				GeneScore[] maxORF=getMaxORF(gene);
				if(gene.hasCDS()){
					writer.write(gene.getName()+"_intron"+"\tIntron\t"+gene.getIntronRNASeqScore().getScanPvalue()+"\t"+this.computeTE(maxORF[0], maxORF[1])+"\n");
				}
				else{
					writer.write(gene.getName()+"\tNonCoding\t"+gene.getGeneRNASeqScore().getScanPvalue()+"\t"+this.computeTE(maxORF[0], maxORF[1])+"\n");
				}
			}
		}
		
		writer.close();
	}


	private void writeORFScores(Collection<GeneComponent> scoresByFeature, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(GeneComponent gc: scoresByFeature){
			if(gc.hasPredictedORFs()){
				GeneScore[] maxORF=getMaxORF(gc);
				GeneScore ribo=maxORF[0];
				GeneScore rna=maxORF[1];	
				if(ribo!=null && rna!=null){
					double TE=this.computeTE(ribo, rna);
					RefSeqGene gene=ribo.getGene();
					gene.setBedScore(TE);
					writer.write(gene+"\n");
				}
				
			}
		}
		
		writer.close();
	}


	private void addExpressionToPredictedORFs(Collection<GeneComponent> features, AlignmentDataModelStats data, ContinuousDataAlignmentModel model) throws IOException {
		for(GeneComponent feature: features){
			if(feature.hasPredictedORFs()){
				RefSeqGene gene=feature.getGene();
				Map<RefSeqGene, GeneScore>  predictedORFs=feature.getPredictedORFScores();
				Map<RefSeqGene, GeneScore> predictedORFsRNASeq=new TreeMap<RefSeqGene, GeneScore>();
				if(predictedORFs!=null && !predictedORFs.isEmpty()){
					IntervalTree<Alignment> tree=data.getIntervalTree(gene.getChr(), gene.getStart(), gene.getEnd());
					for(RefSeqGene orf: predictedORFs.keySet()){
						GeneScore score=new GeneScore(orf, model.scoreGene(orf, tree));
						predictedORFsRNASeq.put(orf, score);
					}
				}
				feature.setPredictedORFRNASeqScores(predictedORFsRNASeq);
			}		
		}
	}


	private void scoreEachORF(Collection<GeneComponent> scoresByFeature, AlignmentDataModelStats data, ContinuousDataAlignmentModel model, Sequence chrSequence, boolean trim3UTRs) throws IOException {
		for(GeneComponent gc: scoresByFeature){
			RefSeqGene gene=null;
			if(!gc.hasCDS()){gene=gc.getGene();}
			else{gene=gc.getIntron();}
			if(gene!=null){
				IntervalTree<Alignment> tree=data.getIntervalTree(gene.getChr(), gene.getStart(), gene.getEnd());
				Collection<RefSeqGene> orfs=FindORFs.findAllORFs(gene, chrSequence, trim3UTRs);
				Map<RefSeqGene, double[]> scores=model.scoreGenes(orfs, gene.getChr(), tree);
				gc.setPredictedORFScores(scores);
			}
		}
	}


	private void addExpression(Collection<GeneComponent> features, ContinuousDataAlignmentModel model, AlignmentDataModelStats data) throws IOException {
		for(GeneComponent feature: features){
			RibosomeScoring.addExpression(feature, model, data);
		}		
	}
	
	private Collection<RefSeqGene> getSubGenes(RefSeqGene gene, int windowSize) {
		Collection<RefSeqGene> subGenes=new TreeSet<RefSeqGene>();
		for(int i=0; i<gene.getTranscriptLength(); i++){
			RefSeqGene subGene=gene.trim(i, i+windowSize);
			subGenes.add(subGene);
		}
		return subGenes;
	}

	/*private void writeFeatures(Collection<GeneComponent> scoresByFeature, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("Name\tGene Type \tNumber of Exons \tRNA-Seq Significance\tSignificance\tGene Counts\tCDS Counts\t5'UTR Counts\t3'UTR Counts\tIntron Counts\tGene Expression\tCDS Expression\t5'UTR Expression\t3'UTR Expression\tIntron Expression\tCDS TE\n");
		for(GeneComponent gene: scoresByFeature){
			GeneScore[] maxORF=getMaxORF(gene);
			String geneType="noncoding";
			if(gene.hasCDS()){geneType="coding";}
			writer.write(gene.getName()+"\t"+geneType+"\t"+gene.getGene().getNumExons()+"\t"+gene.getGeneRNASeqScore().getScanPvalue()+"\t"+gene.getSignificance()+"\t"+gene.getGeneScore().getFullyContainedNumberOfReads());
			if(gene.hasCDS()){
				writer.write("\t"+gene.getCDSScore().getFullyContainedNumberOfReads()+"\t"+gene.getUTR5Score().getFullyContainedNumberOfReads()+"\t"+gene.getUTR3Score().getFullyContainedNumberOfReads());
			}
			else{
				if(maxORF!=null && maxORF[0]!=null){
					writer.write("\t"+maxORF[0].getNumberOfReads()+"\tNA\tNA");
				}
				else{
					writer.write("\tNA\tNA\tNA");
				}
			}
			
			if(gene.hasIntron()){
				writer.write("\t"+gene.getIntronScore().getFullyContainedNumberOfReads());
			}
			else{
				writer.write("\tNA");
			}
			
			writer.write("\t"+gene.getGeneRNASeqScore().getFullyContainedNumberOfReads());
			if(gene.hasCDS()){
				writer.write("\t"+gene.getCDSRNASeqScore().getFullyContainedNumberOfReads()+"\t"+gene.getUTR5RNASeqScore().getFullyContainedNumberOfReads()+"\t"+gene.getUTR3RNASeqScore().getFullyContainedNumberOfReads());
			}
			else{
				if(maxORF!=null && maxORF[1]!=null){
					writer.write("\t"+maxORF[1].getNumberOfReads()+"\tNA\tNA");
				}
				else{
					writer.write("\tNA\tNA\tNA");
				}
			}
			
			if(gene.hasIntron()){
				writer.write("\t"+gene.getIntronRNASeqScore().getFullyContainedNumberOfReads());
			}
			else{
				writer.write("\tNA");
			}
			
			if(gene.hasCDS()){
				System.err.println(gene.getName()+" "+gene.getCDSScore().getNumberOfReads()+" "+gene.getCDSRNASeqScore().getNumberOfReads()+" "+this.computeTE(gene.getCDSScore(), gene.getCDSRNASeqScore()));
				writer.write("\t"+this.computeTE(gene.getCDSScore(), gene.getCDSRNASeqScore()));
			}
			else{
				if(maxORF!=null && maxORF[0]!=null && maxORF[1]!=null){
					writer.write("\t"+this.computeTE(maxORF[0], maxORF[1]));
				}
				else{writer.write("\t-999");}
			}
			
			writer.write("\n");
		}
		
		writer.close();
	}*/
	
	private void writeFeatures(Collection<GeneComponent> scoresByFeature, FileWriter writer) throws IOException {
		//writer.write("Name\tRibosome p-value\tExpression p-value\tRibosome counts\tExpression counts\tTranslational Efficiency\n");
		System.err.println("writing... "+scoresByFeature.size());
		
		for(GeneComponent gene: scoresByFeature){
			//try{
			//System.err.println(gene.getName()+" "+gene.getGene().getAlignment());
			
			RefSeqGene g=gene.getGene();
			g.setName("s="+gene.getGeneTE());
			
			//if(gene.hasCDS()){
				writer.write(g+"\n");
				//writer.write("\t"+gene.getGeneScore().getScanPvalue());
				//writer.write("\t"+gene.getGeneRNASeqScore().getScanPvalue());
				//writer.write("\t"+gene.getGeneScore().getFullyContainedNumberOfReads());
				//writer.write("\t"+gene.getGeneRNASeqScore().getFullyContainedNumberOfReads());
				//writer.write("\t"+gene.getGeneTE()+"\n");
				
				
				/*double intronP=1.0;
				if(gene.getIntron()!=null){
					intronP=gene.getIntronRNASeqScore().getScanPvalue();
				}
				
				gene.getGeneRNASeqScore().getScanPvalue();
				gene.getCDSRNASeqScore().getScanPvalue();
				gene.getUTR5RNASeqScore().getScanPvalue();
				gene.getUTR3RNASeqScore().getScanPvalue();
				
				writer.write("\t"+gene.getGeneRNASeqScore().getScanPvalue()+"\t"+gene.getCDSRNASeqScore().getScanPvalue()+"\t"+intronP+"\t"+gene.getUTR5RNASeqScore().getScanPvalue()+"\t"+gene.getUTR3RNASeqScore().getScanPvalue());
				writer.write("\t"+gene.getGeneTE()+"\t"+gene.getCDSTE()+"\t"+gene.getIntronTE()+"\t"+gene.getUTR5TE()+"\t"+gene.getUTR3TE()+"\n");*/
			//}
			//}catch(NullPointerException ex){ex.printStackTrace();}
		}

	}

	/*private Map<String, Map<RefSeqGene, double[]>> scoreFeatures(Collection<RefSeqGene> genes) throws IOException{
		Collection<RefSeqGene> UTR3=new TreeSet<RefSeqGene>();
		Collection<RefSeqGene> UTR5=new TreeSet<RefSeqGene>();
		Collection<RefSeqGene> CDS=new TreeSet<RefSeqGene>();
		Collection<RefSeqGene> introns=new TreeSet<RefSeqGene>();
		
		for(RefSeqGene gene: genes){
			//for each gene get UTRs, CDSes
			UTR3.add(gene.get3UTRGene());
			UTR5.add(gene.get5UTRGene());
			CDS.add(gene.getCDS());
			introns.add(gene.getIntrons());
		}
		
		Map<String, Map<RefSeqGene, double[]>> rtrn=new TreeMap<String, Map<RefSeqGene, double[]>>();
		rtrn.put("UTR3", model.scoreGenes(UTR3));
		rtrn.put("UTR5", model.scoreGenes(UTR5));
		rtrn.put("CDS", model.scoreGenes(CDS));
		rtrn.put("intron", model.scoreGenes(introns));
		rtrn.put("gene", model.scoreGenes(genes));
		
		return rtrn;
	}*/
	
	private GeneScore[] getMaxORF(GeneComponent gene) {
		//Get the max TE
		if(gene.hasPredictedORFs()){
			double maxTE=0;
			GeneScore ribo1=null;
			GeneScore rna1=null;
			Map<RefSeqGene, GeneScore> orfs=gene.getPredictedORFScores();
			for(RefSeqGene orf: orfs.keySet()){
				GeneScore ribo=orfs.get(orf);
				GeneScore rna=gene.getPredictedORFRNASeqScores().get(orf);
				//double TE=computeTE(ribo, rna);
				double TE=ribo.getAvergeNumberOfReads();
				if(TE>maxTE){
					maxTE=TE;
					ribo1=ribo;
					rna1=rna;
				}
			}
			//System.err.println(gene.getName()+" "+gene.getPredictedORFScores().size()+" "+maxTE);
			GeneScore[] rtrn={ribo1, rna1};
			return rtrn;
		}
		return null;
	}


	private double computeTE(GeneScore ribo, GeneScore rna) {
		return RibosomeScoring.computeTE(ribo, rna, ALPHA);
	}
	
	

	private Collection<GeneComponent> scoreFeatures(Collection<RefSeqGene> genes, ContinuousDataAlignmentModel model, AlignmentDataModelStats data) throws IOException{
		Collection<GeneComponent> rtrn=new ArrayList<GeneComponent>();
		
		for(RefSeqGene gene: genes){
			GeneComponent component=RibosomeScoring.scoreFeature(gene, model, data);
			rtrn.add(component);
		}
		
		return rtrn;
	}

	private void writeTEByExpression(String save,Collection<GeneComponent> geneComponents, boolean started) throws IOException {
		FileWriter writer=new FileWriter(save, started);
		
		for(GeneComponent gene: geneComponents){
			writer.write(gene.getGeneName()+"\t"+gene.getGeneTE(ALPHA)+"\t"+gene.getGeneTEAllReads(ALPHA)+"\t"+gene.getGeneRNASeqScore().getAvergeNumberOfReads()+"\n");
		}
		
		writer.close();
	}

	
	public static void main(String[] args)throws IOException{
		if(args.length>3){
			String BAMFile=args[0];
			String sizeFile=args[1];
			Map<String, Collection<RefSeqGene>> genesByChr=BEDFileParser.loadDataByChr(new File(args[2]));
			String save=args[3];
			String chrToUse=null;
			String expressionFile=null;
			
			if(args.length>4){
				chrToUse=args[4];
			}
			if(args.length>5){
				expressionFile=args[5];
			}
			
			boolean started=false;
			for(String chr: genesByChr.keySet()){
				if(chrToUse==null || !chrToUse.startsWith("chr") || chrToUse.equalsIgnoreCase(chr) ){
					try{
					System.err.println(chr);
					ContinuousDataAlignmentModel ribosomeModel=SequenceUtils.getDataModel(BAMFile, sizeFile, chr, false);
					ContinuousDataAlignmentModel expModel=SequenceUtils.getDataModel(expressionFile, sizeFile, chr, false);
					ComputeRibosomeOccupancyByFeature crof=	new ComputeRibosomeOccupancyByFeature(ribosomeModel, expModel);
					Collection<GeneComponent> geneComponents=crof.scoreFeatures(genesByChr.get(chr));
					
					crof.writeTEByExpression(save, geneComponents, started);
					
					//crof.writeORF3UTRRatio(save, geneComponents, started);
					//crof.writeScatter(save+".scatter", geneComponents, started);
					// normalize RRS by expression ratio
					ComputeRibosomeOccupancyByFeature.writeRRS(save, geneComponents, started, true, true);
					crof.writeScatter(save+".scatter", geneComponents, started);
					started=true;
					}catch(NullPointerException ex){System.err.println("skipping "+chr);}
				}
			}
		}
		else{
			System.err.println(usage);
		}
	}
	

	static String usage=" args[0]=BAM file \n args[1]=size file \n args[2]=genes (full BED) \n args[3]=save \n args[4]=chr (optional)";
	
}
