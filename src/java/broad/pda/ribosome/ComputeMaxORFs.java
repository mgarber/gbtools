package broad.pda.ribosome;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.broad.igv.sam.Alignment;

import broad.core.datastructures.IntervalTree;
import broad.core.parser.CommandLineParser;
import broad.core.sequence.SequenceUtils;
import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;
import broad.pda.ribosome.misc.FindORFs;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GeneScore;

public class ComputeMaxORFs {

	private static double alpha=0.01;
	private ContinuousDataAlignmentModel ribosomeData;
	private ContinuousDataAlignmentModel expressionData;
	Map<RefSeqGene, Collection<RefSeqGene>> orfs;
	Map<RefSeqGene, Collection<RefSeqGene>> nonORFs;
		
	Collection<RefSeqGene> genes;
	boolean hasCDS;
	
	
	//score the max ORF and non-ORF
	public ComputeMaxORFs(ContinuousDataAlignmentModel ribosomeData, ContinuousDataAlignmentModel expressionData, Collection<RefSeqGene> genes, String genomeDir, boolean trim3UTRs) throws IOException{
		this.ribosomeData=ribosomeData;
		this.expressionData=expressionData;
		this.genes=genes;
		
		//To speed this up and utilize the known coding sequences we will check if the gene hasCDS() and use it
		this.hasCDS=checkIfHasCDS(genes);
		
		if(!hasCDS){
			FindORFs orfFinder=new FindORFs(genes, genomeDir, trim3UTRs);
			this.orfs=orfFinder.getAllORFs();
			this.nonORFs=orfFinder.getAllNonORFs();
			System.err.println("Found ORFs");
		}
		else{
			this.orfs=getCDS(genes);
			//this.nonORFs=getUTRs(genes);
		}
		
		
	}
	
	/*private Map<RefSeqGene, Collection<RefSeqGene>> getUTRs(Collection<RefSeqGene> genes2) {
		Map<RefSeqGene, Collection<RefSeqGene>> rtrn=new TreeMap<RefSeqGene, Collection<RefSeqGene>>();
		
		for(RefSeqGene gene: genes2){
			Collection<RefSeqGene> list=new TreeSet<RefSeqGene>();
			if(gene.get3UTRGene()!=null){
				gene.get3UTRGene()
				list.add(gene.get3UTRGene());}
			if(gene.get5UTRGene()!=null){list.add(gene.get5UTRGene());}
			rtrn.put(gene, list);
		}
		
		return rtrn;
	}*/

	private Map<RefSeqGene, Collection<RefSeqGene>> getCDS(Collection<RefSeqGene> genes2) {
		Map<RefSeqGene, Collection<RefSeqGene>> rtrn=new TreeMap<RefSeqGene, Collection<RefSeqGene>>();
		
		for(RefSeqGene gene: genes2){
			Collection<RefSeqGene> list=new TreeSet<RefSeqGene>();
			list.add(gene);
			rtrn.put(gene, list);
		}
		
		return rtrn;
	}

	private boolean checkIfHasCDS(Collection<RefSeqGene> genes2) {
		boolean hasCDS=true;
		
		for(RefSeqGene gene: genes2){
			if(!gene.hasCDS()){return false;}
		}
		
		return hasCDS;
	}

	private Map<RefSeqGene, GeneScore> scoreLocal(Map<RefSeqGene, Collection<RefSeqGene>> orfs) throws IOException {
		Map<RefSeqGene, GeneScore> rtrn=new TreeMap<RefSeqGene, GeneScore>();
		
		int i=0;
		for(RefSeqGene gene: orfs.keySet()){
			System.err.println(gene.getName()+" "+i+" "+orfs.size());
			GeneComponent fullScore=RibosomeScoring.scoreFeature(gene, ribosomeData, ribosomeData.getAlignmentDataModelStats());
			double lambda=fullScore.getGeneScore().getNumberOfReads()/gene.getSize();
			IntervalTree<Alignment> tree=ribosomeData.getAlignmentDataModelStats().getIntervalTree(gene.getChr(), gene.getStart(), gene.getEnd());
			
			for(RefSeqGene orf: orfs.get(gene)){
				GeneScore windowScore=ribosomeData.scoreGene(orf.getCDS(), tree, lambda, gene.getSize());
				GeneScore windowScoreGlobal=new GeneScore(orf.getCDS(), ribosomeData.scoreGene(orf.getCDS(), tree));
				
				if(windowScoreGlobal.getScanPvalue()>windowScore.getScanPvalue()){
					windowScore=windowScoreGlobal;
				}
				
				//System.err.println(windowScore.getScanPvalue()+" "+windowScoreGlobal.getScanPvalue());
				
				rtrn.put(orf, windowScore);
			}
			i++;
		}
		return rtrn;
	}

	private Map<RefSeqGene, Collection<GeneComponent>> score(Map<RefSeqGene, Collection<RefSeqGene>> orfs) throws IOException {
		Map<RefSeqGene, Collection<GeneComponent>> orfScores=new TreeMap<RefSeqGene, Collection<GeneComponent>>();
		
		int i=0;
		//iterate over each gene and get ORFs
		for(RefSeqGene gene: orfs.keySet()){
			System.err.println("Scoring " + gene.getName());
			Collection<GeneComponent> scores=RibosomeScoring.scoreFeature(gene, orfs.get(gene), ribosomeData, expressionData);
			orfScores.put(gene, scores);
			i++;
			//System.err.println("scoring "+i+" "+orfs.size());
		}
		
		return orfScores;
	}

	
	
	private static void writeRatio(FileWriter writer, Map<RefSeqGene, GeneComponent> maxORF, Map<RefSeqGene, GeneComponent> maxNonORF) throws IOException {
		
		for(RefSeqGene gene: maxORF.keySet()){
			GeneComponent orf=maxORF.get(gene);
			GeneComponent nonORF=maxNonORF.get(gene);
			if(orf!=null && nonORF!=null){
				double orfTE=orf.getGeneTE(alpha);
				double nonORFTE=nonORF.getGeneTE(alpha);
				if(orfTE!=-99 && nonORFTE!=-99){
					writer.write(gene.getName()+"\t"+(orfTE/nonORFTE)+"\t"+(orf.getGeneTEAllReads(alpha)/nonORF.getGeneTEAllReads(alpha))+"\t"+(orf.getGeneScore().getFullyContainedNumberOfReads()/nonORF.getGeneScore().getFullyContainedNumberOfReads())+"\t"+(orf.getGeneScore().getNumberOfReads()/nonORF.getGeneScore().getNumberOfReads())+"\n");
				}
			}
		}
		
	}

	
	
	
	private static void writeAll(String save, Map<RefSeqGene, Collection<GeneComponent>> orfScores, boolean started) throws IOException {
		FileWriter writer=new FileWriter(save, started);
		for(RefSeqGene gene: orfScores.keySet()){
			int counter=0;
			for(GeneComponent orf: orfScores.get(gene)){
				//String name=gene.getName()+"_"+counter;
				double te=orf.getCDSTE(alpha);
				double teAll=orf.getCDSTEAllReads(alpha);
				RefSeqGene g=orf.getCDS();
				g.setName("s="+teAll);
				//if(te!=-99 && teAll!=-99){
					writer.write(g+"\n");
				
				//}
				counter++;
			}
		}
		writer.close();
	}

	private static void writeMax(String save, Map<RefSeqGene, GeneComponent> max, boolean started) throws IOException {
		FileWriter writer=new FileWriter(save, started);
		
		writeMax(writer, max);
		
		writer.close();
	}
	
	private static void writeMax(FileWriter writer, Map<RefSeqGene, GeneComponent> max) throws IOException {
		
		int i=0;
		for(RefSeqGene gene: max.keySet()){
			GeneComponent gc=max.get(gene);
			double te=gc.getCDSTE(alpha);
			double teAll=gc.getCDSTEAllReads(alpha);
			if(te!=-99 && teAll!=-99){
				writer.write(gc.getGene()+"\n");
			}
			i++;
			System.err.println("writing "+i+" "+max.size());
		}
		
	}

	
	private void scoreAndWriteEachORF(String save, boolean started) throws IOException {
		FileWriter writer=new FileWriter(save, started);
		
		//got through each ORF and score the significance
		for(RefSeqGene gene: this.orfs.keySet()){
			GeneComponent fullScore=RibosomeScoring.scoreFeature(gene, ribosomeData, ribosomeData.getAlignmentDataModelStats());
			double lambda=fullScore.getGeneScore().getNumberOfReads()/gene.getSize();
			IntervalTree<Alignment> tree=ribosomeData.getAlignmentDataModelStats().getIntervalTree(gene.getChr(), gene.getStart(), gene.getEnd());
			
			for(RefSeqGene orf: this.orfs.get(gene)){
				GeneScore windowScore=ribosomeData.scoreGene(orf.getCDS(), tree, lambda);
				RefSeqGene w=windowScore.getGene();
				//TODO Write other format
				if(windowScore.getScanPvalue()<alpha){
					double te=RibosomeScoring.scoreFeature(w, ribosomeData, expressionData).getGeneTEAllReads(alpha);
					w.setName("te="+te);
					writer.write(w+"\n");
				}
			}
		}
		
		writer.close();
	}
	
	private void scoreAndWriteEachORF(String save, Map<RefSeqGene, GeneScore> orfScores, boolean started) throws IOException {
		FileWriter writer=new FileWriter(save, started);
		
		//got through each ORF and score the significance
		
			
			for(RefSeqGene orf: orfScores.keySet()){
				GeneScore windowScore=orfScores.get(orf);
				double p=windowScore.getScanPvalue();
				if(p<alpha){
					RefSeqGene w=windowScore.getGene();
					w.setName("plocal="+p);
					writer.write(w+"\n");
				}
			}
		
		
		writer.close();
	}

	
	private Map<RefSeqGene, GeneComponent> getMax(Map<RefSeqGene, Collection<GeneComponent>> orfScores) {
		Map<RefSeqGene, GeneComponent> rtrn=new TreeMap<RefSeqGene, GeneComponent>();
		
		for(RefSeqGene gene: orfScores.keySet()){
			
			System.err.println("Getting max TE ORF for gene " + gene.getName());
			
			GeneComponent max=null;
			for(GeneComponent orf: orfScores.get(gene)){
				if(max!=null){
					if(orf.getCDSTEAllReads(alpha)>max.getCDSTEAllReads(alpha)){max=orf;}
				}
				else{max=orf;}
			}
			
			if(max!=null){rtrn.put(gene, max);}
		}
		
		return rtrn;
	}

	private static void writeBED(FileWriter writer, Map<RefSeqGene, Collection<GeneComponent>> orfScores) throws IOException {
		for(RefSeqGene gene: orfScores.keySet()){
			for(GeneComponent orf: orfScores.get(gene)){
				RefSeqGene orfGene=orf.getGene();
				orfGene.setName("TE="+orf.getGeneTE(alpha));
				writer.write(orfGene+"\n");
			}
		}
	}
	
	private static void scoreKnownCDS(String save, ContinuousDataAlignmentModel ribosomeData, ContinuousDataAlignmentModel expressionData, Collection<RefSeqGene> genes, boolean started) throws IOException {
		FileWriter writer=new FileWriter(save, started);
		
		if(!started){
			writer.write("Name\tRibosome p-value\tExpression p-value\tORF global p-value\tORF local p-value\t%CDS\t%3'UTR\n");
		}
		
		//got through each ORF and score the significance
		int counter=0;
		for(RefSeqGene gene: genes){
			System.err.println(gene.getName()+" "+counter+" "+genes.size());
			GeneComponent fullScore=RibosomeScoring.scoreFeature(gene, ribosomeData, ribosomeData.getAlignmentDataModelStats());
			double lambda=fullScore.getGeneScore().getNumberOfReads()/gene.getSize();
			IntervalTree<Alignment> tree=ribosomeData.getAlignmentDataModelStats().getIntervalTree(gene.getChr(), gene.getStart(), gene.getEnd());
			
			GeneScore orfScore=ribosomeData.scoreGene(gene.getCDS(), tree, lambda, gene.getSize());
			GeneComponent geneScore=RibosomeScoring.scoreFeature(gene, ribosomeData, expressionData);
			
			double percentCDS=gene.getPercentCDS();
			double percent3UTR=gene.getPercent3UTR();
			
			writer.write(gene.getName()+"\t"+geneScore.getGeneScore().getScanPvalue()+"\t"+geneScore.getGeneRNASeqScore().getScanPvalue()+"\t"+geneScore.getCDSScore().getScanPvalue()+"\t"+orfScore.getScanPvalue()+"\t"+percentCDS+"\t"+percent3UTR+"\n");
			counter++;
		}
		
		writer.close();
		
	}

	/*private static void scoreKnownCDS(String save, ContinuousDataAlignmentModel ribosomeData, ContinuousDataAlignmentModel expressionData, Map<RefSeqGene, Collection<RefSeqGene>> orfs, boolean started) throws IOException {
		FileWriter writer=new FileWriter(save, started);
		
		if(!started){
			writer.write("Name\tRibosome p-value\tExpression p-value\tORF global p-value\tORF local p-value\t%CDS\t%3'UTR\n");
		}
		
		//got through each ORF and score the significance
		int counter=0;
		for(RefSeqGene gene: orfs.keySet()){
			
			System.err.println(gene.getName()+" "+counter+" "+genes.size());
			GeneComponent fullScore=RibosomeScoring.scoreFeature(gene, ribosomeData, ribosomeData.getAlignmentDataModelStats());
			double lambda=fullScore.getGeneScore().getNumberOfReads()/gene.getSize();
			IntervalTree<Alignment> tree=ribosomeData.getAlignmentDataModelStats().getIntervalTree(gene.getChr(), gene.getStart(), gene.getEnd());
			
			GeneScore orfScore=ribosomeData.scoreGene(gene.getCDS(), tree, lambda, gene.getSize());
			GeneComponent geneScore=RibosomeScoring.scoreFeature(gene, ribosomeData, expressionData);
			
			double percentCDS=gene.getPercentCDS();
			double percent3UTR=gene.getPercent3UTR();
			
			writer.write(gene.getName()+"\t"+geneScore.getGeneScore().getScanPvalue()+"\t"+geneScore.getGeneRNASeqScore().getScanPvalue()+"\t"+geneScore.getCDSScore().getScanPvalue()+"\t"+orfScore.getScanPvalue()+"\t"+percentCDS+"\t"+percent3UTR+"\n");
			counter++;
		}
		
		writer.close();
		
	}*/
	
	private static void writeTable(String save, Map<RefSeqGene, GeneComponent> maxScores,	Map<RefSeqGene, GeneComponent> maxNonORFScores, boolean started, boolean normalizeByExpression, boolean fullyContainedReadsInCds) throws IOException {
		FileWriter writer=new FileWriter(save, started);
		
		if(!started){
			writer.write("Gene\tExpression_significance\tMax_ORF_TE\tMax_NON_ORF_TE\tMax_ORF_Ratio\tMax_Non_ORF_Ratio\n");
		}
		
		for(RefSeqGene gene: maxScores.keySet()){
			GeneComponent orf=maxScores.get(gene);
			if(maxNonORFScores==null || maxNonORFScores.get(gene)==null){
				writer.write(gene.getName()+"\t"+orf.getGeneRNASeqScore().getScanPvalue()+"\t"+orf.getCDSTEAllReads(alpha)+"\t"+"NA"+"\t"+orf.getRRS(alpha, normalizeByExpression, fullyContainedReadsInCds)+"\t"+"NA"+"\n");
			}
			else{
				GeneComponent nonORF=maxNonORFScores.get(gene);
				writer.write(gene.getName()+"\t"+orf.getGeneRNASeqScore().getScanPvalue()+"\t"+orf.getCDSTEAllReads(alpha)+"\t"+nonORF.getCDSTEAllReads(alpha)+"\t"+orf.getRRS(alpha, normalizeByExpression, fullyContainedReadsInCds)+"\t"+nonORF.getRRS(alpha, normalizeByExpression, fullyContainedReadsInCds)+"\n");
			}
		}
		
		writer.close();
	}
	

	private static void writeTable(String save, Map<RefSeqGene, Collection<GeneComponent>> orfScores, boolean started, boolean normalizeByExpression, boolean fullyContainedReadsInCds) throws IOException {
		FileWriter writer=new FileWriter(save, started);
		
		if(!started){
			writer.write("Gene_name\tORF_region\tORF_TE\tcdsVsUtrRatio\n");
		}
		
		for(RefSeqGene gene: orfScores.keySet()){
			for(GeneComponent orf: orfScores.get(gene)){
				//System.err.println(orf.getGene());
				writer.write(gene.getName()+"\t"+orf.getCDS().getAlignment().toUCSC()+"\t"+orf.getCDSTEAllReads(alpha)+"\t"+orf.getRRS(alpha, normalizeByExpression, fullyContainedReadsInCds)+"\n");
			}
		}
		
		writer.close();
	}

	

	public static void main(String[] args) throws IOException{
			
			CommandLineParser p = new CommandLineParser();
			p.addStringArg("-g", "Gene bed file", true);
			p.addStringArg("-s", "Sequence directory", true);
			p.addStringArg("-o", "Output file prefix", true);
			p.addStringArg("-r", "Ribosome bam", true);
			p.addStringArg("-e", "Expression bam", true);
			p.addStringArg("-c", "Chr to use", false, null);
			p.addBooleanArg("-u", "Trim 3'UTRs to beginning of next ORF", false, Boolean.valueOf(false));
			p.addBooleanArg("-n", "Normalize RRS by expression", true);
			p.addBooleanArg("-f", "Only count fully contained reads in CDS", false, Boolean.valueOf(true));
			p.parse(args);

			String ribosomeBAM = p.getStringArg("-r");
			String expressionBAM=(p.getStringArg("-e"));
			Map<String, Collection<RefSeqGene>> genesByChr=BEDFileParser.loadDataByChr(new File(p.getStringArg("-g")));
			String genomeDir=p.getStringArg("-s");
			String save=p.getStringArg("-o");
			String chrToUse=p.getStringArg("-c");			
			boolean trimUTRs = p.getBooleanArg("-u").booleanValue();
			boolean normalizeByExpression = p.getBooleanArg("-n");
			boolean fullyContainedReadsInCds = p.getBooleanArg("-f");
			
			String sizes=SequenceUtils.getChrSizesFile(genomeDir);
			
			
			ContinuousDataAlignmentModel ribosomeData=SequenceUtils.getDataModel(ribosomeBAM, sizes, false);
			ContinuousDataAlignmentModel expressionData=SequenceUtils.getDataModel(expressionBAM, sizes, false);
			
			boolean started=false;
			for(String chr: genesByChr.keySet()){
				System.err.println(chr);
				if(chrToUse==null || chrToUse.equalsIgnoreCase(chr)){
					//Sequence chrSeq=SequenceUtils.getChrSequence(genomeDir, chr);
					//if(chrSeq!=null){
						ComputeMaxORFs maxORF=new ComputeMaxORFs(ribosomeData, expressionData, genesByChr.get(chr), genomeDir, trimUTRs);
						Map<RefSeqGene, Collection<GeneComponent>> orfScores=maxORF.score(maxORF.orfs);
						Map<RefSeqGene, GeneComponent> maxScores=maxORF.getMax(orfScores);
						writeMax(save+".maxORF.bed", maxScores, started);
						
						Map<RefSeqGene, GeneComponent> maxNonORFScores=null;
						if(!maxORF.hasCDS){
							Map<RefSeqGene, Collection<GeneComponent>> nonORFScores=maxORF.score(maxORF.nonORFs);
							maxNonORFScores=maxORF.getMax(nonORFScores);
							writeMax(save+".maxNonORF.bed", maxNonORFScores, started);
							writeTable(save+".allNonORFRatio", nonORFScores, started, normalizeByExpression, fullyContainedReadsInCds);
						}
						
						String outmax = save + ".maxORFRatio";
						String outall = save + ".allORFRatio";
						
						if(trimUTRs) {
							outmax += "_3UTRsTrimmed";
							outall += "_3UTRsTrimmed";
						}
						
						writeTable(outmax, maxScores, maxNonORFScores, started, normalizeByExpression, fullyContainedReadsInCds);
						writeTable(outall, orfScores, started, normalizeByExpression, fullyContainedReadsInCds);
						
						
						//Map<RefSeqGene, GeneScore> localScores=maxORF.scoreLocal(maxORF.orfs);
						//maxORF.scoreAndWriteEachORF(save+".orf.scores.bed", localScores, started);
						started=true;
					//}
				}
			}
						
	
		
	}
	

	






	
}
