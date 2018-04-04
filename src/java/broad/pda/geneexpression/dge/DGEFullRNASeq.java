/**
 * 
 */
package broad.pda.geneexpression.dge;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.math.MathException;
import org.apache.log4j.Logger;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.annotation.GTFFileParser;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.AlignmentDataModelStats;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

/**
 * @author skadri
 *
 */
public class DGEFullRNASeq {

	
	static Logger logger = Logger.getLogger(DGEFullRNASeq.class.getName());
	static int DEFAULT_MIN_MAPPING_QUALITY = -1;
	static final int COUNT_SCORE = 2;
	static final int RPKM_SCORE = 4;
	static final int RPK_SCORE = 3;
	
	public DGEFullRNASeq(ArgumentMap argmap) throws IOException, MathException{
		
		/*
		 * Read the names of the alignment file into an array
		 */
		BufferedReader br = new BufferedReader(new FileReader(argmap.getMandatory("alignments")));
		List<String> alignmentFiles = new ArrayList<String>();
		String s;
		while((s = br.readLine())!= null){
			alignmentFiles.add(s);
		}
		
		boolean useConstituentExons = argmap.containsKey("useConstituentExons");
		boolean useConstituentIntrons = argmap.containsKey("useConstituentIntrons");

		/*
		 * Read the annotation file
		 */
		String annotationFile = argmap.getMandatory("annotations");
		/*
		 * Check the format of the annotation file and call the GTF or BED parser accordingly
		 */
		BEDFileParser annotationParser =  annotationFile.endsWith(".gtf") || annotationFile.endsWith(".GTF")? new GTFFileParser(annotationFile) : new BEDFileParser(annotationFile);
		Map<String, Collection<RefSeqGene>> annotations = null;
		List<String> rows = new ArrayList<String>();
		
		if(useConstituentExons) {
			annotationParser.makeGenes(0.1);
			annotations = annotationParser.toConstituentIsoformMap();
			annotationParser.writeFullBed("MadeGenes.bed");
			BufferedWriter ciw = new BufferedWriter(new FileWriter("constituentExons.bed"));
			for(String chr : annotations.keySet()) {
				Collection<RefSeqGene> constituentIsoforms = annotations.get(chr);
				for (RefSeqGene g : constituentIsoforms) {
					ciw.write(g.toBED());
					ciw.newLine();
				}
			}
			ciw.close();
		} else if (useConstituentIntrons){
			annotationParser.makeGenes(0.1);
			annotationParser.writeFullBed("MadeGenes.bed");
			annotations = annotationParser.toConstituentIntroformMap();
			
		}else {
			annotations = annotationParser.toMap();
			
		}
		Collection<RefSeqGene> annotationCollection = new ArrayList<RefSeqGene>();
		for(Collection<RefSeqGene> chrAnnotations : annotations.values()) {
			annotationCollection.addAll(chrAnnotations);
		}
		
		/*
		 * HashMap of gene to rowName
		 */
		HashMap<RefSeqGene, String> duplicateNameMap = new HashMap<RefSeqGene, String>();
		/*
		 * HashMap of gene name to the number of duplicates
		 */
		HashMap<String, Integer> duplicateMap = new HashMap<String, Integer>();
		
		int duplicates=0;
		for(RefSeqGene gene: annotationCollection){
			if(!rows.contains(gene.getName())) {
				String name = gene.getName();
				if(duplicateNameMap.containsKey(gene)){
					logger.info("Entry for "+name+" already exists");
				}
				else{
					rows.add(name);
					duplicateMap.put(name, 1);
					duplicateNameMap.put(gene, name);
				}
			} 
			// If the gene name has another annotation
			else {
				
				if(duplicateNameMap.containsKey(gene)){
					logger.info("Entry for "+gene.getName()+" already exists in "+duplicateNameMap.get(gene));
				}
				else{
					//Row name is now the geneName appended with the duplicate number
					duplicateMap.put(gene.getName(), (duplicateMap.get(gene.getName())+1));
					String name = (gene.getName()+"_"+duplicateMap.get(gene.getName()));
					rows.add(name);
					duplicateNameMap.put(gene, name);
					//logger.warn("Duplicated annotation : "+ gene.toBED());
					duplicates++;
				}
			}
		}
		
		for(RefSeqGene g : annotationCollection){
			logger.info(g.toUCSC()+" has "+duplicateNameMap.get(g));
		}

		/*
		 * Parameters if not provided are set to defaults:
		 * maxPrematureEnd : default 0
		 * max3Pextension : default 0
		 * window: default 500
		 * minMappingQuality : default -1
		 */
		int minimumMappingQuality = argmap.isPresent("minMappingQuality")? argmap.getInteger("minMappingQuality") : DEFAULT_MIN_MAPPING_QUALITY;

		/*
		 * FLAG for WEIGHING READS BY NH FLAG
		 * TRUE by default
		 * Convert string to boolean
		 */
		boolean weighReadsFlag = argmap.isPresent("weighReadsFlag")? (Boolean.parseBoolean(argmap.get("weighReadsFlag"))): true;
		
		/*
		 * FLAG for REMOVING PCR DUPLICATES 
		 * FALSE by default
		 * Convert string to boolean
		 */
		boolean removePCRDuplicatesFlag = argmap.isPresent("removePCRDuplicatesFlag")? (Boolean.parseBoolean(argmap.get("removePCRDuplicatesFlag"))): false;

		/*
		 * The columns of the matrix will be each alignment file name
		 */
		List<String> cols = new ArrayList<String>();
		for(String name: alignmentFiles){
			cols.add(name);
		}
		
		/* 
		 * Initialize the matrix with row and column names
		 */
		MatrixWithHeaders countsMatrix = new MatrixWithHeaders(rows,cols);
		MatrixWithHeaders rpkmMatrix = new MatrixWithHeaders(rows,cols);
		MatrixWithHeaders rpkMatrix = new MatrixWithHeaders(rows,cols);
		/*
		 * Initialize the data models for all alignment files
		 * @param: <alignment flieName> <load_chromosome_stats> <minMappingQuality> <remove_PCR_duplicates> <weigh reads by NH flag>
		 */
/*		ContinuousDataAlignmentModel[] libDataModels = new ContinuousDataAlignmentModel[alignmentFiles.size()];
		for(int i=0;i<alignmentFiles.size();i++){
			libDataModels[i] = AlignmentUtils.loadAlignmentData(alignmentFiles.get(i),true,minimumMappingQuality,removePCRDuplicatesFlag,weighReadsFlag);
		}*/
		
		logger.info("Minimum mapping quality to count reads: " + minimumMappingQuality);
		BEDFileParser.writeFullBED("constituentIsoforms.bed", annotationCollection);
		String save = argmap.getOutput();
		boolean isStranded = argmap.containsKey("stranded");
		
		File[] maskFiles  = argmap.isPresent("maskFileDir") ?  new File(argmap.get("maskFileDir")).listFiles(): null ;

		String sizes = argmap.get("sizeFile");
		for(int i=0;i<alignmentFiles.size();i++){
			String outputname = alignmentFiles.get(i)+"."+save;
			Map<String, Integer> maskFileData = ContinuousDataAlignmentModel.parseMaskFiles(maskFiles);
			if(!isStranded) {
				logger.info("Scoring using all reads for "+alignmentFiles.get(i));
				AlignmentDataModel alignments=new GenericAlignmentDataModel(alignmentFiles.get(i), sizes, false, minimumMappingQuality,removePCRDuplicatesFlag, weighReadsFlag);
				Map<RefSeqGene, double[]> scores=new TreeMap<RefSeqGene, double[]>();				
				runScore(annotations, outputname, maskFileData, alignments, scores);
				ContinuousDataAlignmentModel.writeFullBED(outputname, scores, 1.1);
				for(RefSeqGene g : annotationCollection){
			//		logger.info(g.toUCSC()+" has duplicate name map: "+duplicateNameMap.get(g));
					if(!scores.containsKey(g))
						scores.put(g,new double[9]);
			//		logger.info(g.toUCSC()+" has score:"+ scores.get(g)[COUNT_SCORE]);
					
					countsMatrix.set(duplicateNameMap.get(g), alignmentFiles.get(i), scores.get(g)[COUNT_SCORE]);
					rpkmMatrix.set(duplicateNameMap.get(g), alignmentFiles.get(i), scores.get(g)[RPKM_SCORE]);
					rpkMatrix.set(duplicateNameMap.get(g), alignmentFiles.get(i), scores.get(g)[RPK_SCORE]);
				}
				
			} else {
				logger.info("Scoring minus reads");
				AlignmentDataModel alignments=new GenericAlignmentDataModel(alignmentFiles.get(i), sizes, false, minimumMappingQuality);
				alignments.setNegativeStranded();
				Map<RefSeqGene, double[]> scoresN=new TreeMap<RefSeqGene, double[]>();				
				scoresN = runScore(annotations, outputname, maskFileData, alignments, scoresN);
				ContinuousDataAlignmentModel.writeFullBED(outputname+".minus", scoresN, 1.1);

				logger.info("Scoring plus reads");
				alignments.setPositiveStranded();
				Map<RefSeqGene, double[]> scoresP=new TreeMap<RefSeqGene, double[]>();				
				scoresP = runScore(annotations, outputname, maskFileData, alignments, scoresP);
				ContinuousDataAlignmentModel.writeFullBED(outputname+".plus", scoresP, 1.1);
				
				for(RefSeqGene g : annotationCollection){
					if(!scoresP.containsKey(g))
						scoresP.put(g,new double[9]);
					if(!scoresN.containsKey(g))
						scoresN.put(g,new double[9]);
					countsMatrix.set(duplicateNameMap.get(g), alignmentFiles.get(i), scoresN.get(g)[COUNT_SCORE]+scoresP.get(g)[COUNT_SCORE]);
					rpkmMatrix.set(duplicateNameMap.get(g), alignmentFiles.get(i), scoresN.get(g)[RPKM_SCORE]+scoresP.get(g)[RPKM_SCORE]);
					rpkMatrix.set(duplicateNameMap.get(g), alignmentFiles.get(i), scoresN.get(g)[RPK_SCORE]+scoresP.get(g)[RPK_SCORE]);
				}		
			}	
		}
		
		BufferedWriter outBw = new BufferedWriter(new FileWriter(argmap.getOutput()+".counts.txt"));
		countsMatrix.write(outBw);
		outBw.close();
		BufferedWriter bww = new BufferedWriter(new FileWriter(argmap.getOutput()+".rpkm.txt"));
		rpkmMatrix.write(bww);
		bww.close();
		outBw = new BufferedWriter(new FileWriter(argmap.getOutput()+".rpk.txt"));
		rpkMatrix.write(outBw);
		outBw.close();
		
		BufferedWriter bedBw = new BufferedWriter(new FileWriter(argmap.getOutput()+".geneNameMap.bed"));
		for(RefSeqGene gene: annotationCollection){
			bedBw.write(gene.toBED());
			bedBw.append("\t"+duplicateNameMap.get(gene));
			bedBw.newLine();
		}
		bedBw.close();
		/*
		 * FLAG to return a normalized matrix 
		 * FALSE by default
		 * Convert string to boolean
		 */
		boolean normalizedOutput = argmap.isPresent("normalizedOutput")? (Boolean.parseBoolean(argmap.get("normalizedOutput"))): false;
		if(normalizedOutput){
			writeNormalizedMatrix(argmap.getOutput(), rpkMatrix);
		}
		
	}

	public static void writeNormalizedMatrix(String fileName, MatrixWithHeaders resultMatrix) throws IOException{
		
		Map<String,Double> factors = DGE.calculateNormalizationFactors(resultMatrix);
		MatrixWithHeaders normalizedMatrix = resultMatrix.multiplyColumnsWithConstants(factors);
		String normalizedFileName = fileName+".normalized";
		String normalizedFactorFile = fileName+".coverage";
		BufferedWriter outBw = new BufferedWriter(new FileWriter(normalizedFileName));
		normalizedMatrix.write(outBw);
		outBw.close();
		
		BufferedWriter fBw = new BufferedWriter(new FileWriter(normalizedFactorFile));
		for(String ss:factors.keySet()){
			fBw.write(ss+"\t"+((double)1.0/factors.get(ss))+"\n");
		}
		fBw.close();
	}


	private Map<RefSeqGene, double[]> runScore(Map<String, Collection<RefSeqGene>> annotations, String save,
			Map<String, Integer> maskFileData, AlignmentDataModel alignments,	Map<RefSeqGene, double[]> scores) throws IOException {
		
		AlignmentDataModelStats alignmentData = new AlignmentDataModelStats(alignments, maskFileData);
		ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(alignmentData, maskFileData, 0, 1, null, null);
		for(String chr : annotations.keySet()) {
			logger.info("processing " + chr);
			Collection<RefSeqGene> chrAnnotations = annotations.get(chr);
			scores.putAll(data.scoreGenes(chrAnnotations, chr));
			//System.err.print("done. saving data  to " + save+"."+chr);
		}
		return scores;
	}

}
