package broad.pda.geneexpression.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;
import broad.pda.gene.RefSeqGeneWithIsoforms;
import broad.pda.geneexpression.dge.DGE;

/**
 * Helper class to ReadPileupScanner to write the result table
 * @author skadri
 *
 */
public class PileupScannerResultWriter {

	private int ZSCORE_COLUMN = 1;
	private int COUNT_COLUMN = 0;
	private double ZSCORE_THRESHOLD = 6.0;
	private double COUNT_THRESHOLD = 5.0;
	private int DIST_THRESHOLD = 300;
	
	Map<String,RefSeqGene> primaryPeaks;
	Map<String,RefSeqGene> secondaryPeaks;
	
	public PileupScannerResultWriter(){
		
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
			bw.write(peakToScoreMap.get(primary)[ZSCORE_COLUMN]+"\t");
			
			RefSeqGene gene = geneParser.get(geneName);
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
				bw.write(secondary.toUCSC()+"\t");
				bw.write(peakToScoreMap.get(secondary)[ZSCORE_COLUMN]+"\t");
				if(gene.isNegativeStrand()){
					bw.write(secondary.getEnd()-gene.getEnd()+"\t");
					bw.write(Math.abs(primary.getEnd()-secondary.getEnd())+"\n");
				}
				else{
					bw.write(gene.getStart()-secondary.getStart()+"\t");
					bw.write(Math.abs(primary.getStart()-secondary.getStart())+"\n");
				}
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
					bwS.write(peak.getStart()+"\t"+peakToScoreMap.get(peak)[ZSCORE_COLUMN]+"\n");
				}
			}
		}
		bwP.close();
		bwS.close();
		bw.close();
		
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
		Map<RefSeqGene,Double[]> peak1ToScoreMap = new HashMap<RefSeqGene,Double[]>();
		Map<RefSeqGene,Double[]> peak2ToScoreMap = new HashMap<RefSeqGene,Double[]>();

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
							peakToScoreMap.remove(secondaryPeaks.get(geneName));
							secondaryPeaks.put(geneName, peak);
							peakToScoreMap.put(peak, scores);
						}
					}
				}
				//If there is a secondary peak
				else{
					//if this peak is better than primary
					if(scores[ZSCORE_COLUMN]>peakToScoreMap.get(primaryPeaks.get(geneName))[ZSCORE_COLUMN]){
						if((distanceBetPeaks(peak,primaryPeaks.get(geneName)))>DIST_THRESHOLD){
							//replace
							peakToScoreMap.remove(secondaryPeaks.get(geneName));
							secondaryPeaks.put(geneName, primaryPeaks.get(geneName));
							primaryPeaks.put(geneName, peak);
							peakToScoreMap.put(peak, scores);
						}
					}
					else{
						//if it is better than secondary
						if((distanceBetPeaks(peak,secondaryPeaks.get(geneName)))>DIST_THRESHOLD){
							if(scores[ZSCORE_COLUMN]>peakToScoreMap.get(secondaryPeaks.get(geneName))[ZSCORE_COLUMN]){
								peakToScoreMap.remove(secondaryPeaks.get(geneName));
								secondaryPeaks.put(geneName, peak);
								peakToScoreMap.put(peak, scores);
							}
							else{
								//NOTHING
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
		return line.split(DGE.whitespaceDelimiter).length==4;
	}
	
	private boolean looksLikeTableData(String line){
		return line.split(DGE.whitespaceDelimiter).length==8 && !line.contains("GENE_NAME");
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {

		PileupScannerResultWriter dummy = new PileupScannerResultWriter();
		//args[0] = Input file
		//args[1] = Output file
		//args[2] = RefSeqGene annotation file
//		dummy.writeTable(args[0], args[1], args[2]);
		//args[0] = table1
		//args[1] = table2
		//args[2] = output
		dummy.getAlternativeTSSs(args[0], args[1], args[2]);
	}
	
	static final String usage = "Usage: PileupScannerWriter -task <task name> "+
			"\n\tscore3P: Computes expression of a given annotation set " + 
			"\n\t\t-alignment <Alignment (mapped to genome)> ";
	
			
}
