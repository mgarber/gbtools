package broad.pda.seq.joshuaAnalyses;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import broad.core.datastructures.MatrixWithHeaders;
import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.AlignmentDataModelStats;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

public class ExpressionEstimation {
	
	//String chr="chr19";

	public ExpressionEstimation(File samFile, File geneFile, String sizes, String save)throws IOException{
		AlignmentDataModel alignments=new GenericAlignmentDataModel(samFile.getAbsolutePath(), sizes);
		Map<String, Collection<RefSeqGene>> genes=BEDFileParser.loadDataByChr(geneFile);
		Map<RefSeqGene, double[]> geneExpression=scoreGenes(alignments, genes, 0);
		write(geneExpression, save);
	}
	
	public ExpressionEstimation(File[] samFiles, File geneFile, String sizes, String save, int numberOfReads)throws IOException{
		double lambda=computeLambda(numberOfReads, sizes);
		Map<String, Collection<RefSeqGene>> genes=BEDFileParser.loadDataByChr(geneFile);
		
		List<String> geneNames=new ArrayList();
		List<String> fileNames=new ArrayList();
		for(int i=0; i<samFiles.length; i++){fileNames.add(samFiles[i].getName());}
		for(String chr: genes.keySet()){for(RefSeqGene gene: genes.get(chr)){geneNames.add(gene.getName());}}
		
		MatrixWithHeaders pvalues=new MatrixWithHeaders(geneNames, fileNames);
		MatrixWithHeaders RPKM=new MatrixWithHeaders(geneNames, fileNames);
		
		for(int i=0; i<samFiles.length; i++){
			AlignmentDataModel alignments=new GenericAlignmentDataModel(samFiles[i].getAbsolutePath(), sizes);
			Map<RefSeqGene, double[]> geneExpression=scoreGenes(alignments, genes, lambda);
			for(RefSeqGene gene: geneExpression.keySet()){
				double[] vals=geneExpression.get(gene);
				pvalues.set(gene.getName(), samFiles[i].getName(), vals[0]);
				RPKM.set(gene.getName(), samFiles[i].getName(), vals[4]);
			}
		}
		
		//TODO Convert RPKMs to ranks
		
		RPKM.writeGCT(save+".RPKM.gct");
		pvalues.writeGCT(save+".pvalues.gct");
	}
	
	private double computeLambda(int numberOfReads, String sizes) {
		Map<String, Integer> chromosomeSizes=BEDFileParser.loadChrSizes(sizes);
		
		double sum=0;
		for(String chr: chromosomeSizes.keySet()){sum=sum+chromosomeSizes.get(chr);}
		
		double lambda=(double)numberOfReads/sum;
		System.err.println("NumReads "+numberOfReads+" total length "+sum+" lambda "+lambda);
		
		return lambda;
	}

	private void write(Map<RefSeqGene, double[]> geneExpression, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(RefSeqGene gene: geneExpression.keySet()){
			double[] vals=geneExpression.get(gene);
			writer.write(gene.getName());
			for(int i=0; i<vals.length; i++){writer.write("\t"+vals[i]);}
			writer.write("\n");
		}
		
		writer.close();
	}

	public Map<RefSeqGene, double[]> scoreGenes(AlignmentDataModel alignments, Map<String, Collection<RefSeqGene>> genes, double lambda) throws IOException{
		Map<RefSeqGene, double[]> rtrn=new TreeMap<RefSeqGene, double[]>();
		
		for(String chr : genes.keySet()) {
			System.err.println(chr);
			System.err.println("Lambda: "+lambda);
			try{
			AlignmentDataModelStats alignmentData = new AlignmentDataModelStats(alignments, null, chr, lambda);
			ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(alignmentData);
			
			System.err.println("Made continuous data");
			
			Collection<RefSeqGene> chrAnnotations = genes.get(chr);
			Map<RefSeqGene, double[]> scores=data.scoreGenes(chrAnnotations, chr);
			rtrn.putAll(scores);
			}catch (NullPointerException ex){System.err.println("Skipping "+chr);}
		}
		
		return rtrn;
	}
	
	private Map<Integer, List<String>> binData(MatrixWithHeaders data, int numBins){
		//sort by expression
		MatrixWithHeaders sorted=data.sortList(0);
		
		int valsPerBin=sorted.getNumberRows()/numBins;
		
		Map<Integer, List<String>> valsByBin=new TreeMap<Integer, List<String>>();
		
		
		
		int count=0;
		for(String gene: sorted.getRowNames()){
			int binNumber=count%valsPerBin;
			List<String> vals=valsByBin.get(binNumber);
			if(vals==null){vals=new ArrayList();}
			vals.add(gene);
			valsByBin.put(binNumber, vals);
			count++;
		}
		
		return valsByBin;
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>4){
			File samFile=new File(args[0]);
			File geneFile=new File(args[1]);
			String sizes=args[2];
			String save=args[3];
			int numReads=new Integer(args[4]);
			File[] samFiles={samFile};
			new ExpressionEstimation(samFiles, geneFile, sizes, save, numReads);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=samFile \n args[1]=geneFile \n args[2]=sizes \n args[3]=save \n args[4]=numReads";
}
