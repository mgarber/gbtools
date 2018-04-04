package broad.pda.ribosome;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import org.broad.igv.Globals;
import org.broad.igv.sam.Alignment;

import broad.core.datastructures.IntervalTree;
import broad.core.sequence.SequenceUtils;
import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.AlignmentDataModelStats;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GeneScore;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

public class ComputeRibosomeOccupancyInWindows {
	
	int windowSize=90;
	double alpha=0.01;
	Collection<GeneComponent> scoresByFeature;
	Map<String, GeneComponent> maxExpression;
	
	public ComputeRibosomeOccupancyInWindows(String bamFile, String expressionFile, String sizeFile, Collection<RefSeqGene> genes, int windowSize, String chr) throws IOException{
		Globals.setHeadless(true);
		
		this.windowSize=windowSize;
		
		AlignmentDataModelStats data=new AlignmentDataModelStats(new GenericAlignmentDataModel(bamFile, sizeFile), null, chr);
		ContinuousDataAlignmentModel model=new ContinuousDataAlignmentModel(data);
		this.scoresByFeature=scanWindow(data, model, genes);
		
		data=new AlignmentDataModelStats(new GenericAlignmentDataModel(expressionFile, sizeFile), null, chr);
		model=new ContinuousDataAlignmentModel(data);
		Collection<GeneComponent> expressionByFeature=scanWindow(data, model, genes); //max of the expression as well
		
		this.maxExpression=getMap(expressionByFeature);
		//computeTranslationEfficiency(scoresByFeature, model, maxExpression);
		
		//writeWindowWiggle(scoresByFeature, maxExpression, save);
		
		//TODO We should also classify each window by whether it is in an "ORF" or not
		
		//writeTE(scoresByFeature, maxExpression, save+".wig");
		//writeBED(scoresByFeature, maxExpression, save+".bed");
		//writeMaxTE(scoresByFeature, maxExpression, save+".max.bed");
	}
	
	
	
	private void writeTE(Collection<GeneComponent> scoresByFeature, Map<String, GeneComponent> maxExpression, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(GeneComponent ribosome: scoresByFeature){
			Map<RefSeqGene, GeneScore> ribosomeScores=ribosome.getGeneWindowScores();
			Map<RefSeqGene, GeneScore> expressionScores=maxExpression.get(ribosome.getName()).getGeneWindowScores();
			
			for(RefSeqGene gene: ribosomeScores.keySet()){
				GeneScore ribo=ribosomeScores.get(gene);
				GeneScore expr=expressionScores.get(gene);
				double score=RibosomeScoring.computeTEAllReads(ribo, expr, alpha);
				if(score>0){writer.write(gene.getChr()+"\t"+gene.getStart()+"\t"+(gene.getStart()+1)+"\t"+score+"\n");}
				//writer.write(gene.getChr()+"\t"+gene.getStart()+"\t"+gene.getEnd()+"\t"+ribo.getAvergeNumberOfReads()+"\n");
			}
			
		}
		
		writer.close();
	}
	
	private void writeMaxTE(Collection<GeneComponent> scoresByFeature, Map<String, GeneComponent> maxExpression, String save, boolean started) throws IOException {
		FileWriter writer=new FileWriter(save, started);
		
		for(GeneComponent ribosome: scoresByFeature){
			Map<RefSeqGene, GeneScore> ribosomeScores=ribosome.getGeneWindowScores();
			Map<RefSeqGene, GeneScore> expressionScores=maxExpression.get(ribosome.getName()).getGeneWindowScores();
			
			RefSeqGene max=null;
			double maxTE=-99;
			
			for(RefSeqGene gene: ribosomeScores.keySet()){
				GeneScore ribo=ribosomeScores.get(gene);
				GeneScore expr=expressionScores.get(gene);
				double score=RibosomeScoring.computeTEAllReads(ribo, expr, alpha);
				if(score>maxTE){
					maxTE=score;
					max=gene;
				}
				//if(score>0){writer.write(gene.getChr()+"\t"+gene.getStart()+"\t"+(gene.getStart()+1)+"\t"+score+"\n");}
				//writer.write(gene.getChr()+"\t"+gene.getStart()+"\t"+gene.getEnd()+"\t"+ribo.getAvergeNumberOfReads()+"\n");
			}
			
			if(max!=null){
				GeneScore ribo=ribosomeScores.get(max);
				GeneScore expr=expressionScores.get(max);
				writer.write(ribosome.getName()+"\t"+maxTE+"\t"+expr.getAvergeNumberOfReads()+"\t"+expr.getScanPvalue()+"\n");
				//max.setName("s="+maxTE);
				//writer.write(max+"\n");
			}
			
		}
		
		writer.close();
	}

	
	private void writeBED(Collection<GeneComponent> scoresByFeature, Map<String, GeneComponent> maxExpression, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(GeneComponent ribosome: scoresByFeature){
			Map<RefSeqGene, GeneScore> ribosomeScores=ribosome.getGeneWindowScores();
			Map<RefSeqGene, GeneScore> expressionScores=maxExpression.get(ribosome.getName()).getGeneWindowScores();
			
			for(RefSeqGene gene: ribosomeScores.keySet()){
				GeneScore ribo=ribosomeScores.get(gene);
				GeneScore expr=expressionScores.get(gene);
				double score=RibosomeScoring.computeTEAllReads(ribo, expr, alpha);
				if(score>0){
					gene.setName("s="+score);
					writer.write(gene+"\n");
					//writer.write(gene.getChr()+"\t"+gene.getStart()+"\t"+(gene.getStart()+1)+"\t"+score+"\n");
				}
				//writer.write(gene.getChr()+"\t"+gene.getStart()+"\t"+gene.getEnd()+"\t"+ribo.getAvergeNumberOfReads()+"\n");
			}
			
		}
		
		writer.close();
	}


	private Map<String, GeneComponent> getMap(Collection<GeneComponent> expressionByFeature) {
		Map<String, GeneComponent> rtrn=new TreeMap<String, GeneComponent>();
		
		for(GeneComponent gc: expressionByFeature){
			rtrn.put(gc.getName(), gc);
		}
		
		return rtrn;
	}

	private static void writeWindowWiggle(Collection<GeneComponent> scoresByFeature, Map<String, GeneComponent> maxExpression, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		writer.write("Name\tGene Type\tRNASeq Significance\tGeneReads\t3'UTR\t5'UTR\tCDS\tGene Max Expression\t3' UTR Max Expression\t5' UTR Max Expression\tCDS Max Expression\tGene Paired Expression\t3'UTR Paired Expression\t5' UTR Paired Expression\tCDS Paired Expression\n");
		
		for(GeneComponent gc: scoresByFeature){
			GeneComponent gcMaxExpression=maxExpression.get(gc.getName());
			
			String max="NA";
			String UTR3="NA";
			String UTR5="NA";
			String CDS="NA";
			
			String geneTE="NA";
			String UTR3TE="NA";
			String UTR5TE="NA";
			String CDSTE="NA";
			
			String geneMaxExpression="NA";
			String UTR3MaxExpression="NA";
			String UTR5MaxExpression="NA";
			String CDSMaxExpression="NA";
			
			String geneType="Non-coding";
			
			if(gc.hasCDS()){geneType="Coding";}
			
			try{
				max=""+gc.getMaxGeneWindow().getFullyContainedNumberOfReads();
				geneTE=""+gc.getPairedGeneExpression().getFullyContainedNumberOfReads();
				geneMaxExpression=""+gcMaxExpression.getMaxGeneWindow().getFullyContainedNumberOfReads();
				}catch(NullPointerException ex){}
			try{
				UTR3=""+gc.getMax3UTRWindow().getFullyContainedNumberOfReads(); 
				UTR3MaxExpression=""+gcMaxExpression.getMax3UTRWindow().getFullyContainedNumberOfReads();
				UTR3TE=""+gc.getPaired3UTRExpression().getFullyContainedNumberOfReads();
				}catch(NullPointerException ex){}
			try{
				UTR5=""+gc.getMax5UTRWindow().getFullyContainedNumberOfReads(); 
				UTR5MaxExpression=""+gcMaxExpression.getMax5UTRWindow().getFullyContainedNumberOfReads();
				UTR5TE=""+gc.getPaired5UTRExpression().getFullyContainedNumberOfReads();
				}catch(NullPointerException ex){}
			try{
				CDS=""+gc.getMaxCDSWindow().getFullyContainedNumberOfReads(); 
				CDSMaxExpression=""+gcMaxExpression.getMaxCDSWindow().getFullyContainedNumberOfReads();
				CDSTE=""+gc.getPairedCDSExpression().getFullyContainedNumberOfReads();
				}catch(NullPointerException ex){}
			writer.write(gc.getName()+"\t"+geneType+"\t"+gc.getGeneRNASeqScore().getScanPvalue()+"\t"+max+"\t"+UTR3+"\t"+UTR5+"\t"+CDS+"\t"+geneMaxExpression+"\t"+UTR3MaxExpression+"\t"+UTR5MaxExpression+"\t"+CDSMaxExpression+"\t"+geneTE+"\t"+UTR3TE+"\t"+UTR5TE+"\t"+CDSTE+"\n");
		}
		writer.close();
	}

	//Compute the score for max window and return the ratio
	private void computeTranslationEfficiency(Collection<GeneComponent> scoresByFeature, ContinuousDataAlignmentModel expressionModel, Map<String, GeneComponent> maxExpression) throws IOException {
		
		for(GeneComponent gc: scoresByFeature){
			GeneScore geneWindow=gc.getMaxGeneWindow();
			GeneScore UTR5Window=gc.getMax5UTRWindow();
			GeneScore UTR3Window=gc.getMax3UTRWindow();
			GeneScore CDSWindow=gc.getMaxCDSWindow();
			
			if(geneWindow!=null){
				GeneScore geneWindowExpression=new GeneScore(geneWindow.getGene(), expressionModel.scoreGene(geneWindow.getGene()));
				gc.setGeneTE(geneWindow.getFullyContainedNumberOfReads()/geneWindowExpression.getFullyContainedNumberOfReads());
				gc.setPairedGeneExpression(geneWindowExpression);
				gc.setGeneRNASeqScore(maxExpression.get(gc.getName()).getGeneScore());
			}
			if(UTR5Window!=null){
				GeneScore UTR5WindowExpression=new GeneScore(UTR5Window.getGene(), expressionModel.scoreGene(UTR5Window.getGene()));
				gc.setUTR5TE(UTR5Window.getFullyContainedNumberOfReads()/UTR5WindowExpression.getFullyContainedNumberOfReads());	
				gc.setPaired5UTRExpression(UTR5WindowExpression);
			}
			if(UTR3Window!=null){
				GeneScore UTR3WindowExpression=new GeneScore(UTR3Window.getGene(), expressionModel.scoreGene(UTR3Window.getGene()));
				gc.setUTR3TE(UTR3Window.getFullyContainedNumberOfReads()/UTR3WindowExpression.getFullyContainedNumberOfReads());
				gc.setPaired3UTRExpression(UTR3WindowExpression);
			}
			
			if(CDSWindow!=null){
				GeneScore CDSWindowExpression=new GeneScore(CDSWindow.getGene(), expressionModel.scoreGene(CDSWindow.getGene()));
				gc.setCDSTE(CDSWindow.getFullyContainedNumberOfReads()/CDSWindowExpression.getFullyContainedNumberOfReads());
				gc.setPairedCDSExpression(CDSWindowExpression);
			}
						
		}
		
	}

	private Collection<GeneComponent> scanWindow(AlignmentDataModelStats data, ContinuousDataAlignmentModel model, Collection<RefSeqGene> genes) throws IOException{
		Collection<GeneComponent> rtrn=new ArrayList<GeneComponent>();
				
		int counter=0;
		for(RefSeqGene gene: genes){
			IntervalTree<Alignment> tree=data.getIntervalTreeCached(gene.getChr(), gene.getStart(), gene.getEnd());
			GeneComponent gc=new GeneComponent(gene);
			long start=System.currentTimeMillis();
			Collection<RefSeqGene> windows=gc.getGene().getWindows(windowSize);
			long getWindows=(System.currentTimeMillis()-start);
			start=System.currentTimeMillis();
			gc.setGeneWindowScores(model.scoreGenes(windows, gene.getChr(), tree), windowSize);
			long scoreWindows=System.currentTimeMillis()-start;
			gc.setGeneScore(model.scoreGene(gc.getGene(), tree));
			rtrn.add(gc);
			counter++;
			System.err.println(counter+" "+genes.size()+" "+(getWindows/1000.0)+" "+(scoreWindows/1000.0));
		}
		
		return rtrn;
	}
	
	
	public static void main(String[] args)throws IOException{
		if(args.length>6){
			String bamFile=args[0];
			String expressionFile=args[1];
			String sizeFile=args[2];
			Map<String, Collection<RefSeqGene>> genesByChr=BEDFileParser.loadDataByChr(new File(args[3]));
			String save=args[4];
			int windowSize=new Integer(args[5]);
			String chrToUse=args[6];
			
			
			boolean started=false;
			for(String chr: genesByChr.keySet()){
				if(chrToUse==null || !chrToUse.startsWith("chr") || chrToUse.equalsIgnoreCase(chr) ){
					try{
						ContinuousDataAlignmentModel expressionData=SequenceUtils.getDataModel(bamFile, sizeFile, chr, false);
						System.err.println("Starting for "+chr);
						
						ComputeRibosomeOccupancyInWindows windows=new ComputeRibosomeOccupancyInWindows(bamFile, expressionFile, sizeFile,  genesByChr.get(chr), windowSize, chr);
						windows.writeMaxTE(windows.scoresByFeature, windows.maxExpression, save, started);
						
						started=true;
					}catch(Exception ex){
						System.err.println("Skipping "+chr);
						ex.printStackTrace();
					}
				}
			}
			
			
		}
		else{
			System.err.println(usage);
		}
	}
	
	static String usage=" args[0]=BAM file \n argsp[1]=expression BAM \n args[2]=size file \n args[3]=genes (full BED) \n args[4]=save \n args[5]=window size \n args[6]=chr";
	
}
