package broad.pda.seq.alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Map;

import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

public class TopHatAlignmentPipeline {

	private static int runTopHat(String save, String referenceLocation, String fastq, Runtime run) throws IOException, InterruptedException{
		String command="/seq/mguttman/scripts/TopHat/bin/tophat -o "+save +" "+referenceLocation+" "+fastq;
		Process p=run.exec(command);
		int completed=p.waitFor();
		int exitVal=p.exitValue();
		writeError(p.getErrorStream());
				
		return exitVal;
	}
	
	private static void writeError(InputStream errorStream) throws IOException {
		BufferedReader reader=	new BufferedReader(new InputStreamReader(errorStream));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			System.err.println(nextLine);
		}
		System.err.println();
	}

	private static int computeExpression(String alignmentFile, String sizes, String geneFile, String save, Runtime run) throws IOException, InterruptedException{
		AlignmentDataModel data=new GenericAlignmentDataModel(alignmentFile, sizes);
		ContinuousDataAlignmentModel model=new ContinuousDataAlignmentModel(data);
		Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(geneFile));
		Map<RefSeqGene, double[]> expression=model.scoreGenes(genes);
		write(save, expression);
		
		//String chr="chr1";
		//String command="/seq/mguttman/scripts/ForAlon/EstimateExpressionNew.jar "+alignmentFile+" "+sizes+" "+genes+" "+save+" "+chr;
		//Process p=run.exec(command);
		//int completed=p.waitFor();
		//int exitVal=p.exitValue();
		//writeError(p.getErrorStream());
		//return exitVal;
		return 0;
	}
	
	private static void write(String save, Map<RefSeqGene, double[]> expression) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("Gene Name \t P-Value \t NormalizedEnrichment \t RawEnrichment \n");
		for(RefSeqGene gene: expression.keySet()){
			double[] scanP=expression.get(gene);
			writer.write(gene.getName()+"\t"+scanP[0]+"\t"+scanP[1]+"\t"+scanP[2]+"\n");
		}
		
		writer.close();
	}

	private static int sortAlignmentFile(String alignmentFile, String save, Runtime run) throws IOException, InterruptedException{
		String command="/xchip/igv/tools/igvtools sort "+alignmentFile+" "+save;
		Process p=run.exec(command);
		int completed=p.waitFor();
		int exitVal=p.exitValue();
		writeError(p.getErrorStream());
		return exitVal;
	}
	
	private static int indexAlignmentFile(String alignmentFile, Runtime run) throws IOException, InterruptedException{
		String command="/xchip/igv/tools/igvtools index "+alignmentFile;
		Process p=run.exec(command);
		int completed=p.waitFor();
		int exitVal=p.exitValue();
		writeError(p.getErrorStream());
		return exitVal;
	}
	
	public static void main(String[] args)throws IOException, InterruptedException{
		if(args.length>2){
			String fastq=args[0];
			String bowtie=args[1];
			String alignmentDir=args[2];
						
			Runtime run=java.lang.Runtime.getRuntime();
			
			String alignmentFile=alignmentDir+"/accepted_hits.sam";
			
			int tophatCompleted=runTopHat(alignmentDir, bowtie, fastq, run);
			int sortCompleted= sortAlignmentFile(alignmentFile, alignmentFile+".sorted.sam", run);
			int indexCompleted=indexAlignmentFile(alignmentFile+".sorted.sam", run);
			//int expressionCompleted=computeExpression(alignmentFile+".sorted.sam", sizes, genes, save, run);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=fastq file \n args[1]=Bowtie reference \n args[2]=alignment output directory";
	
}
