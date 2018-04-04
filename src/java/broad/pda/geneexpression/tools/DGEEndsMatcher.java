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

import org.apache.log4j.Logger;

import broad.core.datastructures.IntervalTree;
import broad.core.math.Statistics;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.gene.RefSeqGeneWithIsoforms;
import broad.pda.geneexpression.dge.DGE;

/**
 * Given the peaks assigned to 5' end of a gene and 3' end of a gene, this class functions to match the 3' and 5' ends of the genes
 * @author skadri
 *
 */
public class DGEEndsMatcher {

	Map<String,List<RefSeqGene>> dge3pMap;
	Map<RefSeqGene,Double> dge3PeakToScoreMap;
	Map<String,List<RefSeqGene>> dge5pMap;
	Map<RefSeqGene,Double> dge5PeakToScoreMap;
	Set<String> geneNames;
	
	private int MAX_EXTENSION = 1500;
	private static int MAX_5P_DISTANCE = 5000;
	private static double THRESHOLD = 10.0;
	private static double COV_THRESHOLD = 0.5;
	
	static Logger logger = Logger.getLogger(DGEEndsMatcher.class.getName());
	
	public DGEEndsMatcher(String dge3pFile,String dge5pFile,int[] colIndex) throws IOException{
		
		geneNames = new HashSet<String>();
		dge3pMap = readAnnotationFile(dge3pFile);
		dge5pMap = readAnnotationFile(dge5pFile);
		
		dge3PeakToScoreMap = readAnnotationScores(dge3pFile,colIndex);
		dge5PeakToScoreMap = readAnnotationScores(dge5pFile,colIndex);
	}
	
	public void matchGeneEnds(String outFile) throws IOException{
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
		BufferedWriter bw3one5mult = new BufferedWriter(new FileWriter(outFile+"3pOne.5pMultiple"));
		BufferedWriter bw3mult5one = new BufferedWriter(new FileWriter(outFile+"3pMultiple.5pOne"));
		BufferedWriter bw35mult = new BufferedWriter(new FileWriter(outFile+"3pMultiple.5pMultiple"));
		BufferedWriter bw3Only = new BufferedWriter(new FileWriter(outFile+"3pOnly"));
		BufferedWriter bw5Only = new BufferedWriter(new FileWriter(outFile+"5pOnly"));
		
		for(String gene:geneNames){
			//Gene is expressed in 3p
			if(dge3pMap.containsKey(gene)){
				//BEST: Gene is expressed in 3p and 5p
				if(dge5pMap.containsKey(gene)){
					//ONE 3P PEAK
					if(dge3pMap.get(gene).size()==1){
						//ONE 3P AND ONE 5P PEAK
						if(dge5pMap.get(gene).size()==1){							
							bw.write(gene+"\t"+dge5pMap.get(gene).get(0).toUCSC()+"\t"+dge3pMap.get(gene).get(0).toUCSC()+"\n");
						}
						//ONE 3P AND MULTIPLE 5P PEAKS
						else{
							for(RefSeqGene p:dge5pMap.get(gene)){
								bw3one5mult.write(gene+"\t"+p.toUCSC()+"\t"+dge3pMap.get(gene).get(0).toUCSC()+"\n");
							}
						}
					}
					//MULTIPLE 3P PEAKS
					else{
						//MULTIPLE 3P AND ONE 5P PEAK
						if(dge5pMap.get(gene).size()==1){
							for(RefSeqGene p:dge3pMap.get(gene)){
								bw3mult5one.write(gene+"\t"+dge5pMap.get(gene).get(0).toUCSC()+"\t"+p.toUCSC()+"\n");
							}
						}
						//MULTIPLE 3P AND 5P PEAKS
						else{
							for(RefSeqGene p5:dge5pMap.get(gene)){
								for(RefSeqGene p3:dge3pMap.get(gene)){
									bw35mult.write(gene+"\t"+p5.toUCSC()+"\t"+p3.toUCSC()+"\n");
								}
							}
						}
					}
				}
				//Gene expressed in 3p but not in 5p
				else{
					for(RefSeqGene p:dge3pMap.get(gene)){
						if(dge3PeakToScoreMap.get(p)>THRESHOLD)
							bw3Only.write(gene+"\t"+p.toUCSC()+"\n");
					}
				}
			}
			//Gene not expressed n 3p
			else{
				//Gene is expressed in 5p
				if(dge5pMap.containsKey(gene)){
					for(RefSeqGene p:dge5pMap.get(gene)){
						if(dge5PeakToScoreMap.get(p)>THRESHOLD)
							bw5Only.write(gene+"\t"+p.toUCSC()+"\n");
					}
				}
			}
		}
		bw.close();
		bw3one5mult.close();
		bw3mult5one.close();
		bw35mult.close();
		bw3Only.close();
		bw5Only.close();
	}
	
	/**
	 * This function takes as input a 5'DGE peak file and a 3' DGE peak file, and gives approximate locations of transcripts based on the peak matching.
	 * @param dge5pPeakFile
	 * @param dge3pPeakFile
	 * @param outFile
	 * @throws IOException 
	 */
	public static void matchDGEPeaksToFormTranscripts(String dge5pPeakFile,String dge3pPeakFile,String outFile,String annotationFile) throws IOException{
		
		//Set<String> geneNames = new HashSet<String>();
		List<Map<String, SortedPeakList<RefSeqGene>>> peaks3p = readPeakFile(dge3pPeakFile);
		List<Map<String, SortedPeakList<RefSeqGene>>> peaks5p = readPeakFile(dge5pPeakFile);
		
		BEDFileParser annotationParser = new BEDFileParser(annotationFile);		

		List<RefSeqGene> transcripts = new ArrayList<RefSeqGene>();
		BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
		int counter=0;
		//For POSITIVE STRAND
		int k=0;
		//For each chromosome
		for(String chr:peaks5p.get(k).keySet()){
			if(peaks3p.get(k).containsKey(chr)){
			//5' index
			int i=0;
			RefSeqGene start = getNextPeak(peaks5p.get(k).get(chr),i);
			i = peaks5p.get(k).get(chr).indexOf(start)+1;
			//3' index
			int j=0;
			boolean inTxt = false;
			List<RefSeqGene> curr5pList = new ArrayList<RefSeqGene>();
			List<RefSeqGene> curr3pList = new ArrayList<RefSeqGene>();
			curr5pList.add(start);
			//System.err.println(peaks5p.get(k).get(chr).size()+ " "+peaks3p.get(k).get(chr).size());
			//WHILE THERE IS PEAK DATA IN 5' AND 3' PEAKS
			while(i<peaks5p.get(k).get(chr).size() && j<peaks3p.get(k).get(chr).size()){
				//Get the next 5' peak
				RefSeqGene curr5pGene = getNextPeak(peaks5p.get(k).get(chr),i);
				//Get the next 3' peak
				RefSeqGene curr3pGene = getNextPeak(peaks3p.get(k).get(chr),j);
				if(curr3pGene==null){
					//System.err.println(peaks5p.get(k).get(chr).size()+ " "+peaks3p.get(k).get(chr).size());
					//System.err.println("i = "+i+"j = "+j);
					break;
				}
				if(!inTxt){
					if(curr5pGene==null){
						inTxt = true;
						curr3pList.add(curr3pGene);
						j = peaks3p.get(k).get(chr).indexOf(curr3pGene)+1;
					}
					else{
						//IF 5' PEAK IS BEFORE 3' PEAK
						if(compareTo(curr5pGene, curr3pGene)<0){
							//IF DISTANCE BETWEEN CURR 5' AND THIS 5' PEAK IS LESS THAN THRESHOLD
							if(curr5pGene.getStart()-start.getStart() < MAX_5P_DISTANCE){
								//ADD ITH 5P TO LIST
								curr5pList = addPeakToList(curr5pList,curr5pGene);
							}
							else{
								//NEW LIST AND START
								curr5pList = new ArrayList<RefSeqGene>();
								curr5pList.add(curr5pGene);
							}
							start = curr5pGene;
							i = peaks5p.get(k).get(chr).indexOf(start)+1;;
						}
						//3' PEAK BEFORE THE NEXT 5' PEAK
						else{
							inTxt = true;
							curr3pList.add(curr3pGene);
							j = peaks3p.get(k).get(chr).indexOf(curr3pGene)+1;
						}
					}
				}
				//IN TRANSCRIPT
				else{
					if(curr5pGene==null){
						curr3pList.add(curr3pGene);
						j = peaks3p.get(k).get(chr).indexOf(curr3pGene)+1;
					}
					else{
						//IF 5' PEAK IS BEFORE 3' PEAK
						if(compareTo(curr5pGene, curr3pGene)<0){
							//FINISH TXT
							inTxt = false;
							//OUTPUT
							for(RefSeqGene p5:curr5pList){
								for(RefSeqGene p3:curr3pList){
									if(p5.getStart()<p3.getEnd()){
										RefSeqGene g = new RefSeqGene(chr,p5.getStart(),p3.getEnd());
										g.setOrientation("+");
										g.setName("gene"+counter);
										counter++;
										g.setBedScore(p5.getBedScore()+p3.getBedScore());
										bw.write(g.toBED()+"\n");
										transcripts.add(g);
									}
									else{
										System.out.println(p5.getStart()+" is greater than "+p3.getEnd());
									}
								}
							}
							//RE-INITIALIZE
							curr5pList = new ArrayList<RefSeqGene>();
							curr3pList = new ArrayList<RefSeqGene>();
							start = curr5pGene;
							curr5pList.add(start);
							i = peaks5p.get(k).get(chr).indexOf(start)+1;;
						}
						//IF 3' PEAK BEFORE 5'
						else{
							curr3pList.add(curr3pGene);
							j = peaks3p.get(k).get(chr).indexOf(curr3pGene)+1;
						}
					}
				}
			}
		}
		}
		//FOR NEGATIVE STRAND
		k=1;
		//For each chromosome
		for(String chr:peaks5p.get(k).keySet()){
			if(peaks3p.get(k).containsKey(chr)){
			//5' index
			int i=0;
			//3' index
			int j=0;
			RefSeqGene end = getNextPeak(peaks3p.get(k).get(chr),j);
			j = peaks3p.get(k).get(chr).indexOf(end)+1;
			boolean inTxt = false;
			List<RefSeqGene> curr5pList = new ArrayList<RefSeqGene>();
			List<RefSeqGene> curr3pList = new ArrayList<RefSeqGene>();
			curr3pList.add(end);
			//WHILE THERE IS PEAK DATA IN 5' AND 3' PEAKS
			while(i<peaks5p.get(k).get(chr).size() && j<peaks3p.get(k).get(chr).size()){
				//Get the next 5' peak
				RefSeqGene curr5pGene = getNextPeak(peaks5p.get(k).get(chr),(i));
				//logger.info(chr+"\t"+peaks5p.get(k).get(chr).indexOf(curr5pGene));
				//Get the next 3' peak
				RefSeqGene curr3pGene = getNextPeak(peaks3p.get(k).get(chr),(j));
				if(curr5pGene==null)
					break;
				if(!inTxt){
					if(curr3pGene==null){
						inTxt = true;
						curr5pList = addPeakToList(curr5pList,curr5pGene);
						i = peaks5p.get(k).get(chr).indexOf(curr5pGene)+1;
					}
					else{
						//IF 3' PEAK IS BEFORE 5' PEAK
						if(compareTo(curr3pGene, curr5pGene)<0){
							//IF DISTANCE BETWEEN CURR 3' AND THIS 3' PEAK IS LESS THAN THRESHOLD
							if(curr3pGene.getStart()-end.getStart() < MAX_5P_DISTANCE){
								//ADD ITH 5P TO LIST
								curr3pList.add(curr3pGene);
							}
							else{
								//NEW LIST AND START
								curr3pList = new ArrayList<RefSeqGene>();
								curr3pList.add(curr3pGene);
							}
							end = curr3pGene;
							j = peaks3p.get(k).get(chr).indexOf(curr3pGene)+1;
						}
						//5' PEAK BEFORE THE NEXT 3' PEAK
						else{
							inTxt = true;
							curr5pList = addPeakToList(curr5pList,curr5pGene);
							i = peaks5p.get(k).get(chr).indexOf(curr5pGene)+1;
						}
					}
				}
				//IN TRANSCRIPT
				else{
					if(curr3pGene==null){
						curr5pList = addPeakToList(curr5pList,curr5pGene);
						i = peaks5p.get(k).get(chr).indexOf(curr5pGene)+1;
					}
					else{
						//IF 3' PEAK IS BEFORE 5' PEAK
						if(compareTo(curr3pGene, curr5pGene)<0){
							//FINISH TXT
							inTxt = false;
							//OUTPUT
							for(RefSeqGene p3:curr3pList){
								for(RefSeqGene p5:curr5pList){
									if(p3.getStart()<p5.getEnd()){
										RefSeqGene g = new RefSeqGene(chr,p3.getStart(),p5.getEnd());
										g.setOrientation("-");
										g.setName("gene"+counter);
										counter++;
										g.setBedScore(p5.getBedScore()+p3.getBedScore());
										bw.write(g.toBED()+"\n");
										transcripts.add(g);
									}
									else{
										System.out.println(p3.getStart()+" is greater than "+p5.getEnd());
									}
								}
							}
							//RE-INITIALIZE
							curr5pList = new ArrayList<RefSeqGene>();
							curr3pList = new ArrayList<RefSeqGene>();
							end = curr3pGene;
							curr3pList.add(end);
							j = peaks3p.get(k).get(chr).indexOf(curr3pGene)+1;
						}
						//IF 3' PEAK BEFORE 5'
						else{
							curr5pList = addPeakToList(curr5pList,curr5pGene);
							i = peaks5p.get(k).get(chr).indexOf(curr5pGene)+1;
						}
					}
				}
			}
		}		
		}
		bw.close();
		
		BEDFileParser annotations = new BEDFileParser(transcripts);
		annotations.merge();
		bw = new BufferedWriter(new FileWriter(outFile+".merged.bed"));
		annotations.writeFullBed(bw);
		bw.close();
		
		bw = new BufferedWriter(new FileWriter(outFile+".novel.bed"));
		for(RefSeqGene g:annotations.GetGenes()){
			if(annotationParser.getOverlappers(g).isEmpty())
				bw.write(g.toBED()+"\n");
		}
		bw.close();
		
	}
	
	private static RefSeqGene getNextPeak(SortedPeakList<RefSeqGene> peakList, int i){
		
		//System.err.println("Starting i:"+i);
		while(i<peakList.size()){
			//System.err.println("Peak score"+peakList.get(i).getBedScore());
			if(peakList.get(i).getBedScore()>COV_THRESHOLD){
				return peakList.get(i);
			}
			else{
				i++;
			}
		}
		//System.err.println("Returning null coz no more peaks above threshold");
		return null;
		
	}
	
	private static List<RefSeqGene> addPeakToList(List<RefSeqGene> peakList,RefSeqGene peak){
		
		double max = Double.MIN_VALUE;
		//Get the most enriched peak
		for(RefSeqGene g:peakList){
			if(g.getBedScore()>max){
				max = g.getBedScore();
			}
		}		
		//Check that the added peak is at least 50% abundance of the 
		if(peak.getBedScore()>=(0.5*max)){
			peakList.add(peak);
		}
		
		if(peak.getBedScore()>max){
			List<RefSeqGene> newList = new ArrayList<RefSeqGene>();
			for(RefSeqGene g:peakList){
				if(g.getBedScore()>=(0.5*peak.getBedScore())){
					newList.add(g);
				}
			}
			return newList;
		}
		else
			return peakList;
		
/*		if(peakList.size()>0){
			if(peak.getBedScore()>peakList.get(0).getBedScore()){
				peakList = new ArrayList<RefSeqGene>();
				peakList.add(peak);
			}
		}
		else{
			peakList = new ArrayList<RefSeqGene>();
			peakList.add(peak);
		}
		return peakList;*/
	}
	
	private static List<Map<String, SortedPeakList<RefSeqGene>>> readPeakFile(String fileName) throws IOException{
		
		logger.info("For: "+fileName);
		List<Map<String, SortedPeakList<RefSeqGene>>> peaks = new ArrayList<Map<String, SortedPeakList<RefSeqGene>>>();
		for(int i=0;i<2;i++)
			peaks.add(new HashMap<String, SortedPeakList<RefSeqGene>>());
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(new File(fileName))));
		String nextLine;
		int i=0;
		while ((nextLine = reader.readLine()) != null ) {
	
			if(looksLikeBEDData(nextLine) ){
				RefSeqGene gene = new RefSeqGene(nextLine, false);
				//System.err.println("Gene: " + gene.toBED());
	
				SortedPeakList<RefSeqGene> data=new SortedPeakList<RefSeqGene>();
				
				if(gene.getOrientation().equalsIgnoreCase("+")){
					if(peaks.get(0).containsKey(gene.getChr())){
						data=peaks.get(0).get(gene.getChr());
					}
					data.add(gene);
					peaks.get(0).put(gene.getChr(), data);
				}
				else{
					if(peaks.get(1).containsKey(gene.getChr())){
						data=peaks.get(1).get(gene.getChr());
					}
					data.add(gene);
					peaks.get(1).put(gene.getChr(), data);
				}
			}
			i++;
			if(i%10000==0){logger.info(i);}
		}
	
		logger.info("# peaks on positive strand");
		for(String chr:peaks.get(0).keySet()){
			logger.info(chr+"\t"+peaks.get(0).get(chr).size());
		}
		logger.info("# peaks on negative strand");
		for(String chr:peaks.get(1).keySet()){
			logger.info(chr+"\t"+peaks.get(1).get(chr).size());
		}
		reader.close();
		return peaks;
		
	}
	
	private Map<String,List<RefSeqGene>> readAnnotationFile(String fileName/*,String peaksFile*/) throws IOException{
		
		Map<String,List<RefSeqGene>> geneToPeaksMap = new HashMap<String,List<RefSeqGene>>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(new File(fileName))));
		
		//BufferedWriter bw = new BufferedWriter(new FileWriter(peaksFile));
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			
			if(looksLikeData(nextLine)){
				String[] tokens=nextLine.split(this.whitespaceDelimiter);
				String geneName = tokens[0].split("chr")[0].substring(0, tokens[0].split("chr")[0].length()-1);
				String peakStr = tokens[0].split("_")[tokens[0].split("_").length-1];
				//System.out.println(peakStr.split(":")[1].split("-")[0]);
				RefSeqGene peak = new RefSeqGene(peakStr.split(":")[0],new Integer(peakStr.split(":")[1].split("-")[0]),new Integer(peakStr.split(":")[1].split("-")[1]));
				geneNames.add(geneName);
				
				//bw.write(peak.getChr()+"\t"+peak.getStart()+"\t"+peak.getEnd()+"\n");
				if(geneToPeaksMap.containsKey(geneName)){
					geneToPeaksMap.get(geneName).add(peak);
				}
				else{
					geneToPeaksMap.put(geneName, new ArrayList<RefSeqGene>());
					geneToPeaksMap.get(geneName).add(peak);
				}
			}
		}
		//bw.close();
		return geneToPeaksMap;
	}
	
	private Map<RefSeqGene,Double> readAnnotationScores(String fileName,int[] colIndex) throws IOException{
		
		Map<RefSeqGene,Double> peakToScoreMap = new HashMap<RefSeqGene,Double>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(new File(fileName))));
				
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			
			if(looksLikeData(nextLine)){
				String[] tokens=nextLine.split(this.whitespaceDelimiter);
				String peakStr = tokens[0].split("_")[tokens[0].split("_").length-1];
				RefSeqGene peak = new RefSeqGene(peakStr.split(":")[0],new Integer(peakStr.split(":")[1].split("-")[0]),new Integer(peakStr.split(":")[1].split("-")[1]));
				
				//Average the score
				double score = 0.0;
				for(int c:colIndex){
					score +=(new Double(tokens[c]));
				}
				score = score/colIndex.length;
				peakToScoreMap.put(peak, score);
			}
		}
		
		return peakToScoreMap;
	}
	
	private void assignK4me3PeaksAndCompare(String k4me3PeakFile,String annotationFile,int extension,String outFile) throws IOException{
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
		BufferedWriter bwNot = new BufferedWriter(new FileWriter(outFile+".noOverlap"));
		BEDFileParser k4me3peaks = new BEDFileParser(k4me3PeakFile);
		BEDFileParser geneParser =  new BEDFileParser(annotationFile);
	//	BEDFileParser dge5pPeaks = new BEDFileParser(dge5PeakToScoreMap.keySet());
		
		Map<String,Set<RefSeqGene>> geneTok4PeakMap = new HashMap<String,Set<RefSeqGene>>();
		
		int intoGeneExtension = extension;
		
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
				logger	.info("Processing " + chr);
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
					if(gene.getTranscriptLength()<intoGeneExtension){
						geneStart = gene;
					}
					else{
						/*
						 * get annotation for region of length intoGeneExtension, 0 length from start of transcript
						 */
						geneStart = DGE.getSubAnnotationFromStart(gene, intoGeneExtension, 0);
					}
					if(geneStart!=null){
						IntervalTree<RefSeqGeneWithIsoforms> overlappersTree = k4me3peaks.getOverlappers(geneStart);
								
						if(!overlappersTree.isEmpty()){
							Iterator<RefSeqGeneWithIsoforms> overlappersIter = overlappersTree.valueIterator();
							while(overlappersIter.hasNext()){
								RefSeqGene overlapper = overlappersIter.next();
						/*		int distance;
								if(!gene.isNegativeStrand()){
									distance = Math.abs(gene.getStart()-overlapper.getStart());
								}
								else{
									distance = Math.abs(gene.getEnd()-overlapper.getEnd());
								}
								if(distance<minDistance){
									minDistance = distance;
									geneTok4PeakMap.put(gene.getName(), overlapper);
								}*/
								associatedPeaks.add(overlapper);
							}
						}	
					}
					//All peaks downstream of the gene
					Alignments start = null;
					if(gene.getOrientation().equals("+"))
						start = new Alignments(gene.getChr(), gene.getStart() - MAX_EXTENSION, gene.getStart());
					else
						start = new Alignments(gene.getChr(), gene.getEnd(), gene.getEnd() + MAX_EXTENSION);
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
//						System.out.println(gene.getName()+" has 0 peaks assigned.");
						geneTok4PeakMap.put(gene.getName(),associatedPeaks);
				}
			}
		}
		
		BufferedWriter bw1 = new BufferedWriter(new FileWriter(outFile+".K4peaks.assigned.summary.txt"));
		BufferedWriter bw2 = new BufferedWriter(new FileWriter(outFile+".K4peaks.assigned.txt"));
		int d = 0;
		for(String g:geneTok4PeakMap.keySet()){
			bw1.write(g+"\t"+geneTok4PeakMap.get(g).size()+"\n");
			for(RefSeqGene x:geneTok4PeakMap.get(g)){
				bw2.write(g+"\t"+x.toUCSC()+"\t"+x.getBedScore()+"\n");
			}
			d += geneTok4PeakMap.get(g).size();
		}
		bw1.close();
		bw2.close();
		System.out.println("Total Genes with K4 peaks:"+geneTok4PeakMap.keySet().size()+" with "+d+" peaks.");
		
		int c = 0;
		for(String g:dge5pMap.keySet()){
			c += dge5pMap.get(g).size();
		}
		System.out.println("Total Genes with 5'DGE peaks:"+dge5pMap.keySet().size()+" with "+c+" peaks.");

		int thresholdPeakCtr=0;
		int thresholdGeneCtr = 0;
		int overlapPeakCtr = 0;
		int overlapGeneCtr = 0;
		int noK4data = 0;
		int nooverlap = 0;
		for(String gene:dge5pMap.keySet()){
			if(geneTok4PeakMap.containsKey(gene)){
				boolean notOverlaps = true;
				boolean geneHas5pData = false;
				for(RefSeqGene peak:dge5pMap.get(gene)){
					//Only allow peaks above a certain threshold
					if(peakPassesThreshold(dge5PeakToScoreMap,peak)){
						thresholdPeakCtr++;
						geneHas5pData = true;
						if(overlaps(peak,geneTok4PeakMap.get(gene),-1000,100)){
							bw.write(gene+"\t"+/*geneTok4PeakMap.get(gene).toUCSC()+*/"\n");
							notOverlaps = false;
							overlapPeakCtr++;
						}
					}
					else{
						//logger.info("Gene: "+gene+" has peak "+peak.toUCSC()+" but score "+dge5PeakToScoreMap.get(peak)+" is under threshold of 10");
					}
				}
				if(geneHas5pData){
					thresholdGeneCtr++;
					if(!notOverlaps){
						overlapGeneCtr++;
					}
					else{
						bwNot.write(gene+"\n");
						nooverlap++;
					}
				}
			}			
			else{
				boolean geneHas5pData = false;
				for(RefSeqGene peak:dge5pMap.get(gene)){
					//Only allow peaks above a certain threshold
					if(peakPassesThreshold(dge5PeakToScoreMap,peak)){
						//System.out.println(gene);
						geneHas5pData = true;
						thresholdPeakCtr++;
					}
				}
				if(geneHas5pData){
					thresholdGeneCtr++;
					noK4data++;
					bwNot.write(gene+"\n");
				}
			}
		}
		
		logger.info(thresholdPeakCtr+" peaks for "+thresholdGeneCtr+" genes pass the threshold.");
		logger.info(overlapPeakCtr+" peaks for "+overlapGeneCtr+" genes have a K4me3 domain.");
		logger.info("Total number of genes without K4 data: "+noK4data+" no overlap: "+nooverlap);
		bw.close();
		bwNot.close();
	}
	
	private boolean peakPassesThreshold(Map<RefSeqGene,Double> peakToScoreMap,RefSeqGene peak){
		
		//return (peakToScoreMap.get(peak)>=THRESHOLD);
		
		//75th quantile
		List<Double> vals = new ArrayList<Double>();
		
		//Only add expressed peaks to list
		for(Double d:peakToScoreMap.values()){
			if(d>0)
				vals.add(d);
		}
		Collections.sort(vals);
		double thresh = Statistics.quantile(vals, 0.25);
		//System.out.println(thresh);*/
		
		//List<Double> vals = new ArrayList<Double>(peakToScoreMap.values());
		//Collections.sort(vals);
	//	double thresh = Statistics.median(vals);
		//System.out.println(thresh);
		return (peakToScoreMap.get(peak)>thresh);
	}

	
	/**
	 * Compares two genes and returns a negative number if gene x is before y in linear space. If they have the same start then, if x ends before y
	 * @param x
	 * @param y
	 * @return
	 */
	public static int compareTo(RefSeqGene x, RefSeqGene y) {
		if(x.getChr() != null && y.getChr() != null && !x.getChr().equals(y.getChr())) {
			return x.getChr().compareTo(y.getChr());
		}
		return x.getStart() != y.getStart() ? x.getStart() - y.getStart() : x.getEnd() - y.getEnd();
	}
	
	/**
	 * Returns true if gene x and gene y overlap each other
	 * @param x
	 * @param y
	 * @param threshold
	 * @return
	 */
	private boolean overlaps(RefSeqGene x,RefSeqGene y,int threshold){
		
		return((y.getStart()<=x.getEnd() && x.getStart()<=y.getEnd())||(Math.abs(x.getEnd()-y.getStart())<threshold)||(Math.abs(y.getEnd()-x.getStart())<threshold));
	}
	
	private boolean overlaps(RefSeqGene x,RefSeqGene y,int th1,int th2){
		
		return((y.getStart()<=x.getEnd() && x.getStart()<=y.getEnd())||((y.getStart()-x.getEnd())<(th2))||((y.getEnd()-x.getStart())>(th1)));
	}
	
	/**
	 * Returns true if the gene overlaps any gene in the given set
	 * @param x
	 * @param y
	 * @param threshold
	 * @return
	 */
	private boolean overlaps(RefSeqGene x,Set<RefSeqGene> y,int threshold){
		for(RefSeqGene g:y){
			if(overlaps(x,g,threshold))
				return true;
		}
		return false;
	}
	
	private boolean overlaps(RefSeqGene x,Set<RefSeqGene> y,int th1,int th2){
		for(RefSeqGene g:y){
			if(overlaps(x,g,th1,th2))
				return true;
		}
		return false;
	}
	
	/**
	 * Returns true if the line does not start with GeneName
	 * @param line
	 * @return
	 */
	private boolean looksLikeData(String line){
		return line.trim().length() > 0 && ! line.startsWith("GeneName");
	}
	
	/**
	 * Returns true if the line looks like BED data
	 * @param nextLine
	 * @return
	 */
	public static boolean looksLikeBEDData(String nextLine) {
		return nextLine.trim().length() > 0 && ! nextLine.startsWith("#") && !nextLine.startsWith("track") && !nextLine.startsWith("browser");
	}
	
	public static void compareNovelTranscipts(String trans1,String trans2,String outFile) throws IOException{
		
		BEDFileParser genes1 = new BEDFileParser(trans1);	
		BEDFileParser genes2 = new BEDFileParser(trans2);
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
		for(RefSeqGene g:genes2.GetGenes()){
			if(genes1.getOverlappers(g).isEmpty()){
				bw.write(g+"\n");
			}
		}
		bw.close();
	}
	
	private String whitespaceDelimiter = "\\s++"; //$NON-NLS-1$
	
	static final String usage = "Usage: DGEEndsMatcher -task <task name> "+
			"\n\t\t"
			;
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		//ArgumentMap argmap = CLUtil.getParameters(args,usage,"score");
		int[] colIndex = new int[2];
		colIndex[0] = 5;
		colIndex[1] = 6;
		DGEEndsMatcher dummy = new DGEEndsMatcher("dge5p.genes.annotations.secondread.final.merge.strict.0.1","dge5p.genes.annotations.secondread.final.merge.strict.all.table.normalized",colIndex);
		//dummy.matchGeneEnds("dc.4h.endmatching");
		//dummy.assignK4me3PeaksAndCompare("K4me3_0.all.peaks.edited.bed","refseq.mm9.04_12.uniq.noNRs.nonrandom.genes.bed",1000,"dc.0h.K4.5p.1kBwiggle.final.merge.strict.overlaps");
		dummy.assignK4me3PeaksAndCompare("K4me3.compressed.peaks.edited.bed","refseq.mm9.04_12.uniq.noNRs.nonrandom.genes.bed",1000,"dc.4h.k4me3.5p.final.merge.strict.all.75percentile.-1000.100.overlaps");
		
		//dummy.matchDGEPeaksToFormTranscripts("dge5p.genes.annotations.secondread.5endDC_4h_rep1.peaks.bed",
//				"dge3p.genes.annotations.secondread.3endDC_4h_rep1.peaks.bed","dc.4h.rep1.single.transcripts.bed","refseq.mm9.04_12.uniq.noNRs.nonrandom.genes.bed");
//		dummy.compareNovelTranscipts("dc.4h.rep1.single.transcripts.bed.novel.bed", "dc.0h.single.transcripts.bed.novel.bed", "dc.0hWRT4h.novelTranscripts.bed");
	}

}
