package broad.pda.geneexpression.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import broad.core.datastructures.IntervalTree;
import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.gene.RefSeqGeneWithIsoforms;
import broad.pda.geneexpression.dge.DGE;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;

/**
 * Helper class to ReadPileupScanner to perform downstream analyses
 * @author skadri
 *
 */
public class PileupScannerAnalyzer {

	private int ZSCORE_COLUMN = 1;
	private int COUNT_COLUMN = 0;
	private double ZSCORE_THRESHOLD = 3.0;
	private double COUNT_THRESHOLD = 2.0;
	private int DIST_THRESHOLD = 300;
	
	Map<String,RefSeqGene> primaryPeaks;
	Map<String,RefSeqGene> secondaryPeaks;
	
	static final String usage = "Usage: PileupScannerAnalyzer -task <task name> "+
			"\n\tanalyze: Performs a bunch of downsstream analytical tasks " + 
			"\n\t\t-in <Table from ReadPileupScanner with GeneName	Window	Count	ZScore> "+
			"\n\t\t-out <Output name> "+
			"\n\t\t-annotation <Annotation bed file>"+
			"\n\t\t-dge <DGE result table>"+
			"\n\t\t-cols <column indices separated by commas>"+
			"\n\t\t-k4file <File with K4me3peaks>";
	
	
	public PileupScannerAnalyzer(String args[]) throws IllegalArgumentException, IOException{
		
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"analyze");
		/*
		 * 1. WRITE A SUMMARY RESULT TABLE
		 */
		writeTable(argMap.getInput(),argMap.getOutput(),argMap.getMandatory("annotation"));
		
		int[] colIndex = ContinuousDataAlignmentModel.getWidths(argMap.getMandatory("cols"));
		makeResultTable(argMap.get("dge"), colIndex, argMap.getOutput()+".results",argMap.getMandatory("k4File"),argMap.getMandatory("annotation"));
	}
	
	/**
	 * OUTPUT is a table with
	 * GENE_NAME	PRIMARY_PEAK	Z-SCORE		DIST_FROM_ANNOTATED_END		
	 * 						SECONDARY_PEAK		Z-SCORE		DIST_FROM_ANNOTATED_END		DIST_BET_TWO_PEAKS
	 * @param inFile
	 * @param outFile
	 * @throws IOException
	 */
	public void writeTable(String inFile,String outFile,String annotationFile) throws IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
		BufferedWriter bwP = new BufferedWriter(new FileWriter(outFile+".primary.wig"));
		BufferedWriter bwS = new BufferedWriter(new FileWriter(outFile+".secondary.wig"));
		
		int windowSize = 5;
		BEDFileParser geneParser =  new BEDFileParser(annotationFile);
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(new File(inFile))));
		String nextLine;
				
		primaryPeaks = new HashMap<String,RefSeqGene>();
		secondaryPeaks = new HashMap<String,RefSeqGene>();

		Map<RefSeqGene,Double[]> peakToScoreMap = new HashMap<RefSeqGene,Double[]>();
		
		//For each line in the input file
		while ((nextLine = reader.readLine()) != null ) {
	
			if(looksLikeData(nextLine)){
				
				String[] tokens = nextLine.split(DGE.whitespaceDelimiter);
				
				//First column is gene name
				String geneName = tokens[0];
				RefSeqGene peak = new RefSeqGene(tokens[1].split(":")[0],new Integer(tokens[1].split(":")[1].split("-")[0])
									,new Integer(tokens[1].split(":")[1].split("-")[1]));
				Double[] scores = getScores(tokens);
				if(scores.length==2)
					peakToScoreMap = addPeak(peakToScoreMap,peak,scores,geneName);	
			}
		}
		
		bw.write("GENE_NAME\tPRIMARY_PEAK\tZ-SCORE\tDIST_FROM_ANNOTATED_END\tSECONDARY_PEAK\t"+
					"Z-SCORE\tDIST_FROM_ANNOTATED_END\tDIST_BET_TWO_PEAKS\n");
		//WRITING TO TABLE
		for(String geneName:primaryPeaks.keySet()){
			//First column geneName
			bw.write(geneName+"\t");
			
			RefSeqGene primary = primaryPeaks.get(geneName);
			//Second column primary peak
			//This exists since the keyset is used for the loop
			bw.write(primary.toUCSC()+"\t");
			
			//Third column is z-score
			if(!peakToScoreMap.containsKey(primary)){
				System.err.println(geneName+" does not have a score stored for primary "+primary.toUCSC());
			}
			bw.write(peakToScoreMap.get(primary)[ZSCORE_COLUMN]+"\t");
			RefSeqGene gene = geneParser.get(geneName);
			//System.out.println(geneName);
			//Fourth column is distance from annotated end
			//Negative number = upstream of gene
			if(gene.isNegativeStrand()){
				bw.write(gene.getEnd()-primary.getEnd()+"\t");
			}
			else{
				bw.write(primary.getStart()-gene.getStart()+"\t");
			}
			
			if(secondaryPeaks.containsKey(geneName)){
				RefSeqGene secondary = secondaryPeaks.get(geneName);
				//bw.write(secondary.toUCSC()+"\t");
				if(!peakToScoreMap.containsKey(secondary))
					System.out.println(gene.getName());
				else{
				bw.write(peakToScoreMap.get(secondary)[ZSCORE_COLUMN]+"\t");
				if(gene.isNegativeStrand()){
					bw.write(secondary.getEnd()-gene.getEnd()+"\t");
					bw.write(Math.abs(primary.getEnd()-secondary.getEnd())+"\n");
				}
				else{
					bw.write(gene.getStart()-secondary.getStart()+"\t");
					bw.write(Math.abs(primary.getStart()-secondary.getStart())+"\n");
				}}
			}
			else{
				bw.write("NA\tNA\tNA\tNA\n");
			}
		}
		
		//WRITING PRIMARY and secondary WIG FILE
		Iterator<String> chrom_iter = geneParser.getChromosomeIterator();
		while(chrom_iter.hasNext()) {
			String chr = chrom_iter.next();
			bwP.write("variableStep chrom="+chr+" span="+windowSize+"\n");
			bwS.write("variableStep chrom="+chr+" span="+windowSize+"\n");
			Iterator<RefSeqGeneWithIsoforms> annotation_iter = geneParser.getChrTree(chr).valueIterator();
			while(annotation_iter.hasNext()){
				/*
				 * annotation = current node of RefSeqGeneWithIsoforms in interval tree
				 */
				RefSeqGeneWithIsoforms gene = annotation_iter.next();
				
				if(primaryPeaks.containsKey(gene.getName())){
					RefSeqGene peak = primaryPeaks.get(gene.getName());
					bwP.write(peak.getStart()+"\t"+peakToScoreMap.get(peak)[ZSCORE_COLUMN]+"\n");
				}
				if(secondaryPeaks.containsKey(gene.getName())){
					RefSeqGene peak = secondaryPeaks.get(gene.getName());
					if(!peakToScoreMap.containsKey(peak))
						;//System.out.println(gene.getName());
					else
						bwS.write(peak.getStart()+"\t"+peakToScoreMap.get(peak)[ZSCORE_COLUMN]+"\n");
				}
			}
		}
		bwP.close();
		bwS.close();
		bw.close();
		
	}
	
	public void makeResultTable(String dgeTable,int[] colIndex,String outFile,String k4me3File,String annotationFile) throws IOException{
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));

		// GET THE DGE SCORE MAP
		Map<String,Double> geneToScoreMap = readDGETable(dgeTable,colIndex);
		/*
		 * 1. FOR EACH PERCENTILE interval, get the list of genes
		 */
		List<Double> vals = new ArrayList<Double>();
		
		Map<String,Set<RefSeqGene>> geneToK4Map = assignK4me3Peaks(k4me3File,annotationFile);
		//Only add expressed to list
		for(Double d:geneToScoreMap.values()){
			if(d>0.0)
				vals.add(d);
		}
		
		bw.write("Interval"+"\t"+"Thresholds"+"\t"+"#Genes"+"\t"+"#GenesWithTSS"+"\t"+"#GenesWithMultipleTSS"+"\t"+"#GeneWithPrimaryTSSOverlapsK4me3"+"\t"+"#GeneWithSecondaryTSSOverlapsK4me3"+"\t"+"#GenesWithK4Peaks"+"\n");
		double increment = 10;
		for(double i=0.0;i<100;i+=increment){
			double j=i+increment;
			bw.write((i)+"%-"+(j)+"%"+"\t");
			//GET THE THRESHOLDS
			double minThreshold = getPercentileThreshold(vals,i);
			double maxThreshold = getPercentileThreshold(vals,j);
			
			bw.write((minThreshold)+"-"+(maxThreshold)+"\t");
			//GET THE LIST OF GENES IN THAT INTERVAL
			List<String> selectedGenes = new ArrayList<String>();
			for(String g:geneToScoreMap.keySet()){
				if(geneToScoreMap.get(g)>=minThreshold && geneToScoreMap.get(g)<=maxThreshold){
					selectedGenes.add(g);
				}
			}
			
			//Total number of genes
			bw.write(selectedGenes.size()+"\t");
			
			//#Genes with TSSs
			int tssCounter = 0;
			//#Genes with multiple TSSs
			int multTSSCounter = 0;
			//#Genes with K4 peak
			int k4Counter = 0;
			for(String geneName:selectedGenes){
				if(primaryPeaks.containsKey(geneName)){
					tssCounter++;
				}
				if(secondaryPeaks.containsKey(geneName)){
					multTSSCounter++;
				}
				if(geneToK4Map.containsKey(geneName)){
					k4Counter++;
				}
			}
			
			bw.write(tssCounter+"\t");
			bw.write(multTSSCounter+"\t");
			
			int cnt[] = getOverlapWithK4me3(selectedGenes, geneToK4Map);
			bw.write(cnt[0]+"\t");
			bw.write(cnt[1]+"\t");
			bw.write(k4Counter+"\t");
			
			bw.write("\n");
		}
		bw.close();
	}
	
	/**
	 * Assigns K4 peaks from a bed file to the genes in the annotation file
	 * @param k4me3PeakFile
	 * @param annotationFile
	 * @return
	 * @throws IOException
	 */
	private Map<String,Set<RefSeqGene>> assignK4me3Peaks(String k4me3PeakFile,String annotationFile) throws IOException{
		
		BEDFileParser k4me3peaks = new BEDFileParser(k4me3PeakFile);
		BEDFileParser geneParser =  new BEDFileParser(annotationFile);
		
		Map<String,Set<RefSeqGene>> geneTok4PeakMap = new HashMap<String,Set<RefSeqGene>>();
		
		int extension = 1000;
		
		/*
		 * ASSIGN THE K4ME3 PEAKS TO THE GENES FIRST
		 */
		Iterator<String> chrIt = geneParser.getChromosomeIterator();
		/*
		 * For each chromosome
		 */
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			/*
			 * If the alignment data has data from that chromosome
			 */
			if(k4me3peaks.containChr(chr)){
				/*
				 * Obtain an iterator over the interval tree values built from the annotations (genes) on that particular chromosome.
				 * IntervalTree of RefSeqGeneWithIsoforms
				 */
				Iterator<RefSeqGeneWithIsoforms> annotation_iter = geneParser.getChrTree(chr).valueIterator();
				/*
				 * Parse annotation tree
				 * While there is an annotated RefSeqGeneWithIsoforms in Interval tree to analyze
				 * Thus, for each gene
				 */
				while(annotation_iter.hasNext()){
					/*
					 * annotation = current node of RefSeqGeneWithIsoforms in interval tree
					 */
					RefSeqGeneWithIsoforms gene = annotation_iter.next();
					Set<RefSeqGene> associatedPeaks = new HashSet<RefSeqGene>();
					RefSeqGene geneStart;
					
					double minDistance = Integer.MAX_VALUE;
					//Get all peaks intoGene bases inside this gene
					if(gene.getTranscriptLength()<extension){
						geneStart = gene;
					}
					else{
						/*
						 * get annotation for region of length intoGeneExtension, 0 length from start of transcript
						 */
						geneStart = DGE.getSubAnnotationFromStart(gene, extension, 0);
					}
					if(geneStart!=null){
						IntervalTree<RefSeqGeneWithIsoforms> overlappersTree = k4me3peaks.getOverlappers(geneStart);
								
						if(!overlappersTree.isEmpty()){
							Iterator<RefSeqGeneWithIsoforms> overlappersIter = overlappersTree.valueIterator();
							while(overlappersIter.hasNext()){
								RefSeqGene overlapper = overlappersIter.next();
								associatedPeaks.add(overlapper);
							}
						}	
					}
					//All peaks downstream of the gene
					Alignments start = null;
					if(gene.getOrientation().equals("+"))
						start = new Alignments(gene.getChr(), gene.getStart() - extension, gene.getStart());
					else
						start = new Alignments(gene.getChr(), gene.getEnd(), gene.getEnd() + extension);
					//Get all peaks inside the alignment
					IntervalTree<RefSeqGeneWithIsoforms> overlappersTree = k4me3peaks.getOverlappers(start);
					if(!overlappersTree.isEmpty()){
						Iterator<RefSeqGeneWithIsoforms> overlappersIter = overlappersTree.valueIterator();
						while(overlappersIter.hasNext()){
							RefSeqGene overlapper = overlappersIter.next();
							associatedPeaks.add(overlapper);
						}
					}
					
					if(associatedPeaks.size()>0)
						geneTok4PeakMap.put(gene.getName(),associatedPeaks);
				}
			}
		}
		
		return geneTok4PeakMap;

	}
	/**
	 * Get the expression in that percentile
	 * @param vals
	 * @param percent
	 * @return
	 */
	private double getPercentileThreshold(List<Double> vals,double percent){
		
		if(percent<0.0 || percent>100.00){
			throw new IllegalArgumentException("Percentile should be between 0.0 and 100.0. Supplied: "+percent);
		}
		double p = percent/100.00;			
		Collections.sort(vals);
		if(percent==0.0)
			return vals.get(0).doubleValue();
		if(percent==100.00)
			return vals.get(vals.size()-1).doubleValue();
		return Statistics.quantile(vals, p);
	}
	/**
	 * Reads the DGE table and returns a map of gene to DGE score (averaged across the reps indices specified)
	 * @param dgeTable
	 * @param colIndex
	 * @return
	 * @throws IOException 
	 */
	private Map<String,Double> readDGETable(String dgeTable,int[] colIndex) throws IOException{
		
		Map<String,Double> geneToScoreMap = new HashMap<String,Double>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(new File(dgeTable))));
				
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			
			if(looksLikeDGEData(nextLine)){
				String[] tokens=nextLine.split(DGE.whitespaceDelimiter);
				String geneName = tokens[0];
				
				//Average the score
				double score = 0.0;
				for(int c:colIndex){
					score +=(new Double(tokens[c]));
				}
				score = score/colIndex.length;
				geneToScoreMap.put(geneName, score);
			}
		}
		
		return geneToScoreMap;
	}
	
	/**
	 * 
	 * @param genes
	 * @param geneTok4PeakMap
	 * @return 	[0] overlap with primary peaks
	 * 			[1] overlap with secondary peaks
	 */
	private int[] getOverlapWithK4me3(List<String> genes,Map<String,Set<RefSeqGene>> geneTok4PeakMap){
		
		int[] tssWithK4 = new int[2];
		tssWithK4[0] = 0;
		tssWithK4[1] = 0;
		int upstream = 0;
		int downstream =0;
		for(String gene:genes){
			//If there is a TSS site for this gene
			if(primaryPeaks.containsKey(gene)){
				//If there is a k4 peak for the specified gene
				if(geneTok4PeakMap.containsKey(gene)){
					
					//Check for overlap
					//this peak already passes thresholds
					RefSeqGene tss = primaryPeaks.get(gene);
					if(overlaps(tss,geneTok4PeakMap.get(gene),upstream,downstream)){
						tssWithK4[0]++;
					}
					
				}
				else{ 
					//TSS does not have supporting K4me3 peak
				}
			}
			if(secondaryPeaks.containsKey(gene)){
				//If there is a k4 peak for the specified gene
				if(geneTok4PeakMap.containsKey(gene)){
					
					//Check for overlap
					//this peak already passes thresholds
					RefSeqGene tss = secondaryPeaks.get(gene);
					if(overlaps(tss,geneTok4PeakMap.get(gene),upstream,downstream)){
						tssWithK4[1]++;
					}
					
				}
				else{ 
					//TSS does not have supporting K4me3 peak
				}
			}
		}
		return tssWithK4;
	}
	
	private boolean overlaps(RefSeqGene x,RefSeqGene y,int th1,int th2){
		
		return((y.getStart()<=x.getEnd() && x.getStart()<=y.getEnd())||((y.getStart()-x.getEnd())<(th2))||((y.getEnd()-x.getStart())>(th1)));
	}
	
	private boolean overlaps(RefSeqGene x,Set<RefSeqGene> y,int th1,int th2){
		for(RefSeqGene g:y){
			if(overlaps(x,g,th1,th2))
				return true;
		}
		return false;
	}
	
	/**
	 * USES THE OUTPUT OF WRITE TABLE TO COMPARE TWO CONDITIONS
	 * @param inFile1
	 * @param inFile2
	 * @param outFile
	 * @throws IOException 
	 * @throws NumberFormatException 
	 */
	public void getAlternativeTSSs(String inFile1,String inFile2,String outFile) throws NumberFormatException, IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));

		Map<String,RefSeqGene> file1PeakMap = new HashMap<String,RefSeqGene>();
		Map<String,RefSeqGene> file2PeakMap = new HashMap<String,RefSeqGene>();

		//Read file 1
		//For each line in the input file
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(new File(inFile1))));
		String nextLine;
		while ((nextLine = reader.readLine()) != null ) {
	
			if(looksLikeTableData(nextLine)){
				
				String[] tokens = nextLine.split(DGE.whitespaceDelimiter);
				
				//First column is gene name
				String geneName = tokens[0];
				RefSeqGene peak = new RefSeqGene(tokens[1].split(":")[0],new Integer(tokens[1].split(":")[1].split("-")[0])
									,new Integer(tokens[1].split(":")[1].split("-")[1]));				
				file1PeakMap.put(geneName, peak);
				Double[] scores = getScores(tokens);
			}
		}
		
		//Read file2
		reader=new BufferedReader(new InputStreamReader(new FileInputStream(new File(inFile2))));
		while ((nextLine = reader.readLine()) != null ) {
	
			if(looksLikeTableData(nextLine)){
				
				String[] tokens = nextLine.split(DGE.whitespaceDelimiter);
				
				//First column is gene name
				String geneName = tokens[0];
				RefSeqGene peak = new RefSeqGene(tokens[1].split(":")[0],new Integer(tokens[1].split(":")[1].split("-")[0])
									,new Integer(tokens[1].split(":")[1].split("-")[1]));				
				file2PeakMap.put(geneName, peak);	
			}
		}
		
		//For genes in both files
		for(String geneName:file1PeakMap.keySet()){
			if(file2PeakMap.containsKey(geneName)){
				
				//Compare the peaks
				if((distanceBetPeaks(file1PeakMap.get(geneName),file2PeakMap.get(geneName))>DIST_THRESHOLD)){
					bw.write(geneName+"\t"+file1PeakMap.get(geneName).toUCSC()+"\t"+file2PeakMap.get(geneName).toUCSC()+
							"\n");
				}
			}
		}
		
		bw.close();
	}
	/**
	 * This function will add the peak to the map 
	 * 				if it makes the zscore threshold
	 * 					if there is less than or just one other peak in the map
	 * 					if there are already 2 peaks in the map, keep the top highest peaks in the map
	 * @param peakToScoreMap
	 * @param peak
	 * @param scores
	 * @return
	 */
	private Map<RefSeqGene,Double[]> addPeak(Map<RefSeqGene,Double[]> peakToScoreMap,RefSeqGene peak,Double[] scores,String geneName){
		
		//if it makes the zscore threshold
		if(scores[ZSCORE_COLUMN]>ZSCORE_THRESHOLD && scores[COUNT_COLUMN]>COUNT_THRESHOLD){
			//if there is no peak assigned
			if(!primaryPeaks.containsKey(geneName)){
				primaryPeaks.put(geneName,peak);
				peakToScoreMap.put(peak, scores);
			}
			else{
				//if no secondary peak only primary
				if(!secondaryPeaks.containsKey(geneName)){
					//If the distance between the two peaks is more than threshold
					if((distanceBetPeaks(peak,primaryPeaks.get(geneName)))>DIST_THRESHOLD){
						//if this peak is better than primary 
						if((scores[ZSCORE_COLUMN]>peakToScoreMap.get(primaryPeaks.get(geneName))[ZSCORE_COLUMN])){
							//replace
							secondaryPeaks.put(geneName, primaryPeaks.get(geneName));
							primaryPeaks.put(geneName, peak);
							peakToScoreMap.put(peak, scores);
						}
						//make this peak secondary. primary remains same
						else{
							//peakToScoreMap.remove(secondaryPeaks.get(geneName));
							secondaryPeaks.put(geneName, peak);
							peakToScoreMap.put(peak, scores);
						}
					}
					//Check if the other peak is better, then replace
					else{
						//if this peak is better than primary 
						if((scores[ZSCORE_COLUMN]>peakToScoreMap.get(primaryPeaks.get(geneName))[ZSCORE_COLUMN])){
							//replace and remove primary
							//peakToScoreMap.remove(primaryPeaks.get(geneName));
							primaryPeaks.put(geneName, peak);
							peakToScoreMap.put(peak, scores);
						}
					}
				}
				//If there is a secondary peak
				else{
					//if this peak is better than primary
					if(scores[ZSCORE_COLUMN]>peakToScoreMap.get(primaryPeaks.get(geneName))[ZSCORE_COLUMN]){
						if((distanceBetPeaks(peak,primaryPeaks.get(geneName)))>DIST_THRESHOLD){
							//replace secondary with previous primary
							//peakToScoreMap.remove(secondaryPeaks.get(geneName));
							secondaryPeaks.put(geneName, primaryPeaks.get(geneName));
							primaryPeaks.put(geneName, peak);
							peakToScoreMap.put(peak, scores);
						}
						//this peak is better than the primary and too close
						//Replace primary with this peak
						else{
							//replace and remove primary
							//peakToScoreMap.remove(primaryPeaks.get(geneName));
							primaryPeaks.put(geneName, peak);
							peakToScoreMap.put(peak, scores);
						}
					}
					else{
						//if it is better than secondary
						if(scores[ZSCORE_COLUMN]>peakToScoreMap.get(secondaryPeaks.get(geneName))[ZSCORE_COLUMN]){
							//If it isnt too close to the primary, replace the current secondary
							if((distanceBetPeaks(peak,primaryPeaks.get(geneName)))>DIST_THRESHOLD){
								//peakToScoreMap.remove(secondaryPeaks.get(geneName));
								secondaryPeaks.put(geneName, peak);
								peakToScoreMap.put(peak, scores);
							}
							//if the distance to primary is too close do not replace the secondary
							else{
								
							}
						}
					}
					
				}
			}
			
		}
		return peakToScoreMap;
	}
	
	/**
	 * Returns the distance between the two specified peaks in genomic space
	 * @param peak1
	 * @param peak2
	 * @return
	 */
	private int distanceBetPeaks(RefSeqGene peak1,RefSeqGene peak2){
		
		return (Math.max(peak1.getStart()-peak2.getEnd(), peak2.getStart()-peak1.getEnd()));
	}
	
	/**
	 * This function returns the scores for the peak
	 * @param tokens
	 * @return
	 */
	private Double[] getScores(String[] tokens){
		
		Double[] scores = new Double[2];
		//COUNTS
		scores[0] = new Double(tokens[2]);
		//Z-SCORE
		scores[1] = new Double(tokens[3]);
		return scores;
	}
	/**
	 * Returns true if the line has some data
	 * @param line
	 * @return
	 */
	private boolean looksLikeData(String line){
		return line.split(DGE.whitespaceDelimiter).length==4 && !line.contains("Infinity");
	}
	
	private boolean looksLikeTableData(String line){
		return line.split(DGE.whitespaceDelimiter).length==8 && !line.contains("GENE_NAME");
	}
	
	private boolean looksLikeDGEData(String line){
		return line.split(DGE.whitespaceDelimiter).length>0 && line.contains("gene_");
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {

		new PileupScannerAnalyzer(args);
	}
	
			
}
