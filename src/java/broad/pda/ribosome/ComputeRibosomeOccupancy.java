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

import broad.core.parser.CommandLineParser;
import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

public class ComputeRibosomeOccupancy {

	ContinuousDataAlignmentModel model;
	int windowSize=90;
	double alpha=0.05;
	
	public ComputeRibosomeOccupancy(AlignmentDataModel data, Collection<RefSeqGene> genes, String save) throws IOException{
		Globals.setHeadless(true);
		model=new ContinuousDataAlignmentModel(data);
				
		System.err.println("Finished building the model");		
		//lets compute counts for each genomic feature within a gene model
		Collection<GeneComponent> scoresByFeature=scoreFeatures(genes);
		
		writeFeatures(scoresByFeature, save); //write all components of the gene as part of the same line
		
		System.err.println("Finished scoring features");
		
		//scan fixed window and write max
		/*Map<RefSeqGene, Map<RefSeqGene, double[]>> scanScores=scanWindow(genes);
		writeMax(scanScores, save+".max");		
		System.err.println("Finished scanning windows");*/
	}
	
	
	private void writeMax(Map<RefSeqGene, Map<RefSeqGene, double[]>> scanScores, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(RefSeqGene gene: scanScores.keySet()){
			double max=0;
			Map<RefSeqGene, double[]> scores=scanScores.get(gene);
			for(RefSeqGene subScore: scores.keySet()){
				double[] vals=scores.get(subScore);
				max=Math.max(max, vals[1]);
			}
			writer.write(gene.getName()+"\t"+gene.getTranscriptLength()+"\t"+max+"\n");
		}
		
		writer.close();
	}

	private Map<RefSeqGene, Map<RefSeqGene, double[]>> scanWindow(Collection<RefSeqGene> genes) throws IOException{
		Map<RefSeqGene, Map<RefSeqGene, double[]>> rtrn=new TreeMap<RefSeqGene, Map<RefSeqGene, double[]>>();
		
		for(RefSeqGene gene: genes){
			Collection<RefSeqGene> subGenes=getSubGenes(gene, windowSize);
			rtrn.put(gene, model.scoreGenes(subGenes));
		}
		
		return rtrn;
	}
	
	
	private Collection<RefSeqGene> getSubGenes(RefSeqGene gene, int windowSize) {
		Collection<RefSeqGene> subGenes=new TreeSet<RefSeqGene>();
		for(int i=0; i<gene.getTranscriptLength(); i++){
			RefSeqGene subGene=gene.trim(i, i+windowSize);
			subGenes.add(subGene);
		}
		return subGenes;
	}

	private void writeFeatures(Collection<GeneComponent> scoresByFeature, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		/*writer.write("Name\t Significance\tGene Enrich\tCDS Enrich\t5'UTR Enrich\t3'UTR Enrich\tIntron Enrich\n");
		for(GeneComponent gene: scoresByFeature){
			writer.write(gene.getName()+"\t"+gene.getSignificance()+"\t"+gene.getGeneEnrichment()+"\t"+gene.getCDSEnrichment()+"\t"+gene.getUTR5Enrichment()+"\t"+gene.getUTR3Enrichment()+"\t"+gene.getIntronEnrichment()+"\n");
		}*/
		
		
		writer.write("Name\tGeneScore\tCdsScore\tIntronScore\tUtr3Score\tUtr5Score\n");
		
		for(GeneComponent gene : scoresByFeature) {
			writer.write(gene.getName() + "\t" + gene.getGeneScore() + "\t "+ gene.getCDSScore() + "\t" + gene.getIntronScore() + "\t" + gene.getUTR3Score() + "\t" + gene.getUTR5Score() + "\n");
		}
		
		
		writer.close();
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
	
	private Collection<GeneComponent> scoreFeatures(Collection<RefSeqGene> genes) throws IOException{
		Collection<GeneComponent> rtrn=new ArrayList<GeneComponent>();
		
		for(RefSeqGene gene: genes){
			try{
			//for each gene get UTRs, CDSes
			GeneComponent component=new GeneComponent(gene);
			component.setCDSScore(model.scoreGene(component.getCDS()));
			component.setGeneScore(model.scoreGene(component.getGene()));
			component.setIntronScore(model.scoreGene(component.getIntron()));
			component.setUTR3Score(model.scoreGene(component.get3UTR()));
			component.setUTR5Score(model.scoreGene(component.get5UTR()));
			rtrn.add(component);	
			}catch (NullPointerException ex){gene.getName();}
		}
		
		return rtrn;
	}

	public static void main(String[] args)throws IOException{

		CommandLineParser p = new CommandLineParser();
		
		p.addStringArg("-b", "bam file", true);
		p.addStringArg("-s", "size file", true);
		p.addStringArg("-g", "genes", true);
		p.addStringArg("-o", "output file", true);
		
		p.parse(args);
		
		String BAMFile = p.getStringArg("-b");
		String sizeFile = p.getStringArg("-s");
		AlignmentDataModel data=new GenericAlignmentDataModel(BAMFile, sizeFile);
		Collection<RefSeqGene> genes = BEDFileParser.loadData(new File(p.getStringArg("-g")));
		String save = p.getStringArg("-o");
		new ComputeRibosomeOccupancy(data, genes, save);
			
	}
	
	static String usage=" args[0]=BAM file \n args[1]=size file \n args[2]=genes (full BED) \n args[3]=save";
	
}
