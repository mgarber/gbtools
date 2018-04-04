package broad.pda.geneexpression.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.annotation.GTFFileParser;
import broad.pda.gene.RefSeqGeneWithIsoforms;
/**
 * @author skadri
 *
 */
public class UniqueOrthologs{

	
	static final String usage = "Usage: UniqueOrthologs -task <task name> "+
			"\n\tmake: Computes the unique orthologs of a gene list" +
			"\n\t\t-annotations1 <RefSeq bed file for species 1 [BED by default]> "+
			"\n\t\t-annotations2 <RefSeq bed file for species 2 [BED by default]> "+
			"\n\t\t-genelist <Orthlog list> "+
			"\n\t\t-minOverlap MinimumOverlap required for isoforms. [Default is 40%] "+
			"\n\t\t-out <Output file [Defaults to stdout]> "
			;
	
	static Logger logger = Logger.getLogger(UniqueOrthologs.class.getName());
	double minOverlap;
	
	public UniqueOrthologs (String[] args)throws IOException{
	
		Globals.setHeadless(true);
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"make");
		
		/*
		 * Read the annotation file for first species
		 */
		String annotations1File = argMap.getMandatory("annotations1");
		/*
		 * Read the annotation file for second species
		 */
		String annotations2File = argMap.getMandatory("annotations2");
		/*
		 * Check the format of the annotation files and call the GTF or BED parser accordingly
		 */
		BEDFileParser annotations1 =  annotations1File.endsWith(".gtf") || annotations1File.endsWith(".GTF")? new GTFFileParser(annotations1File) : new BEDFileParser(annotations1File);
		BEDFileParser annotations2 =  annotations2File.endsWith(".gtf") || annotations2File.endsWith(".GTF")? new GTFFileParser(annotations2File) : new BEDFileParser(annotations2File);

		/*  
		 * This is actually a bidirectional map but not implemented as such.
		 */
		Map<String,String> orthologs = Collections.synchronizedMap(new HashMap<String,String>());
		
		minOverlap = argMap.isPresent("minOverlap")? argMap.getDouble("minOverlap") : 0.4;
		
		/*
		 * Output Filename
		 */
		String outputFileName = argMap.get("out");
		
		/*
		 * Read the gene list of orthologs
		 */
		String genelistFile = argMap.getMandatory("genelist");
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(genelistFile))); 
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			
			//FOR EACH ENTRY
			String[] tokens=nextLine.split("\t");
			//There should be two tokens 
			if(tokens.length!=2){
				if(tokens.length>2)
					logger.info("Error in line "+tokens.toString()+" : Too many columns");
				else
					logger.info("Error in line "+tokens.toString()+" : Too few columns");
			}
			else{
				String gene1 = tokens[0];
				String gene2 = tokens[1];
				//IF GENE1 HAS 1 ID
				if(gene1.split(",").length==1&annotations1.get(gene1)!=null){
					
					//IF GENE1 HAS NOT BEEN CONSIDERED ALREADY
					if(!orthologs.containsKey(gene1)){
						//IF GENE2 HAS SINGLE ID
						if(gene2.split(",").length==1&&annotations2.get(gene2)!=null){
							//IF GENE2 NOT ALREADY CONSIDERED
							if(!orthologs.containsValue(gene2)){
								//ADD TO LIST
								orthologs.put(gene1, gene2);
							}
							else{//GENE2 is already in list
								//ORTH = ORTHOLOG OF GENE 2 IN LIST
								String orth = findKeyForValue(orthologs,gene2);
								//IF PREVIOUS ORTHOLOG OF GENE2 is ISOFORM OF GENE1
								if (passTwoGeneIsoformTest(annotations1,gene1,orth)){
									//IGNORE GENE1 ENTRY
									logger.info(gene2+" has an ortholog "+orth+" which is an isoform of "+gene1+".Ignore "+gene1+"\t"+gene2);
								}
								//PREVIOUS ORTHOLOG OF GENE2 IS A DIFFERENT GENE
								else{
									//REMOVE ENTRY FROM HASHMAP
									logger.info(gene2+" has an ortholog "+orth+" and "+gene1+" which are not overlapping isoforms. Ignore both "+gene1+"\t"+gene2+" and "+orth+"\t"+gene2);
									orthologs.remove(orth);
								}
							}
						}
						//GENE2 HAS MULTIPLE IDS
						else{
							//IF ALL GENES IN GENE2 ARE NOT OVERLAPPING ISOFORMS
							String[] genes = gene2.split(",");
							if(!passIsoformTest(annotations2,genes)){
								//IGNORE ENTRY
								logger.info("All genes in "+gene2+" are not overlapping isoforms. Ignore "+gene1+"\t"+gene2);
							}
							//IF ALL GENES IN GENE2 ARE OVERLAPPING ISOFORMS
							else{
								boolean gene2AlreadyConsidered=false;
								// FOR EACH GENE IN GENE2
								for(int i=0;i<genes.length;i++){
									String g2 = genes[i];
									//IF G2 IS ALREADY IN LIST
									if(orthologs.containsValue(g2)){
										gene2AlreadyConsidered=true;
										//ORTH = ORTHOLOG OF G2 IN LIST
										System.out.println(g2);
										String orth2 = findKeyForValue(orthologs,g2);
										//IF ORTH2 AND GENE1 ARE ISOFORMS
										if (passTwoGeneIsoformTest(annotations1,gene1,orth2)){
											//IGNORE GENE1 ENTRY
											logger.info(g2+" in "+gene2+" has an ortholog "+orth2+" which is an isoform of "+gene1+".Ignore "+gene1+"\t"+g2);
										}
										//PREVIOUS ORTHOLOG OF one gene in GENE2 IS A DIFFERENT GENE
										else{
											//REMOVE ENTRY FROM HASHMAP
											logger.info(gene2+" has an ortholog "+orth2+" and "+gene1+" which are not overlapping isoforms. Ignoring entry for "+orth2+"\t"+g2);
											orthologs.remove(g2);
										}
									}
								}
								//IF SOME GENE in GENE2 is ALREADY IN LIST
								if(gene2AlreadyConsidered){
									//IGNORE ENTRY FOR GENE1
									logger.info(gene2+" orthologs were already processed. Ignore "+gene1+"\t"+gene2);
								}
								//NO GENE in GENE2 was ALREADY IN LIST
								else{
									//ADD ONE ISOFORM TO LIST
									orthologs.put(gene1,getClosestInLengthIsoform(annotations2,genes,annotations1.get(gene1)));
								}
							}
						}
					}
					//GENE 1 IS ALREADY IN LIST
					else{
						//ORTH = ORTHOLOG OF GENE 1 IN LIST
						String orth = orthologs.get(gene1);
						//SINGLE OR MORE GENE2 IDS
						String[] genes = gene2.split(",");
						//FOR EACH GENE IN GENE2
						for(int i=0;i<genes.length;i++){
							String g2 = genes[i];
							//IF PREVIOUS ORTHOLOG OF GENE1 is ISOFORM OF GENE2
							if (passTwoGeneIsoformTest(annotations2,gene2,orth)){
								//IGNORE GENE2 ENTRY
								logger.info(gene1+" has an ortholog "+orth+" which is an isoform of "+gene2+".Ignore "+gene1+"\t"+gene2);
							}
							//PREVIOUS ORTHOLOG OF GENE1 IS A DIFFERENT GENE
							else{
								//REMOVE ENTRY FROM HASHMAP
								logger.info(gene1+" has an ortholog "+orth+" and "+gene2+" which are not overlapping isoforms. Ignore both "+gene1+"\t"+gene2+" and "+gene1+"\t"+orth);
								orthologs.remove(gene1);
							}
						}
					}
				}
				//GENE1 HAS MULTIPLE IDS
				else{
					//IF ALL GENES IN GENE1 ARE NOT OVERLAPPING ISOFORMS
					String[] g1s = gene1.split(",");
					if(!passIsoformTest(annotations1,g1s)){
						//IGNORE ENTRY
						logger.info("All genes in "+gene1+" are not overlapping isoforms. Ignore "+gene1+"\t"+gene2);
					}
					//IF ALL GENES IN GENE1 ARE OVERLAPPING ISOFORMS
					else{
						//IF GENE2 HAS SINGLE ID
						if(gene2.split(",").length==1&&annotations2.get(gene2)!=null){
							//IF GENE 2 NOT ALREADY IN LIST
							if(!orthologs.containsValue(gene2)){
								boolean gene1AlreadyConsidered=false;
								// FOR EACH GENE IN GENE1
								for(int i=0;i<g1s.length;i++){
									String g1 = g1s[i];
									//IF G1 IS ALREADY IN LIST
									if(orthologs.containsKey(g1)){
										gene1AlreadyConsidered=true;
										//ORTH1 = ORTHOLOG OF G1 IN LIST
										String orth1 = orthologs.get(g1);
										//IF ORTH1 AND GENE2 ARE ISOFORMS
										if (passTwoGeneIsoformTest(annotations2,gene2,orth1)){
											//IGNORE GENE1 ENTRY
											logger.info(g1+" in "+gene1+" has an ortholog "+orth1+" which is an isoform of "+gene2+".Ignore "+g1+"\t"+gene2);
										}
										//PREVIOUS ORTHOLOG OF one gene in GENE1 IS A DIFFERENT GENE
										else{
											//REMOVE ENTRY FROM HASHMAP
											logger.info(gene1+" has an ortholog "+orth1+" and "+gene2+" which are not overlapping isoforms. Ignoring entry for "+g1+"\t"+orth1);
											orthologs.remove(g1);
										}
									}
								}
								//IF SOME GENE in GENE1 is ALREADY IN LIST
								if(gene1AlreadyConsidered){
									//IGNORE ENTRY FOR GENE1
									logger.info(gene1+" orthologs were already processed. Ignore "+gene1+"\t"+gene2);
								}
								//NO GENE in GENE1 was ALREADY IN LIST
								else{
									//ADD ONE ISOFORM TO LIST
									orthologs.put(getClosestInLengthIsoform(annotations1,g1s,annotations2.get(gene2)),gene2);
								}
							}
							//GENE2 ALREADY IN LIST
							else{
								//ORTH = ORTHOLOG OF GENE2 IN LIST
								String orth2 = findKeyForValue(orthologs,gene2);
								//FOR EACH GENE IN GENES1
								for(int i=0;i<g1s.length;i++){
									String g1 = g1s[i];
									//IF G1 IS ALREADY IN LIST
									if(orthologs.containsKey(g1)){
										//ORTH1 = ORTHOLOG OF G1 IN LIST
										String orth1 = orthologs.get(g1);
										//IF ORTH1 AND GENE2 ARE NOT ISOFORMS
										if (!passTwoGeneIsoformTest(annotations2,gene2,orth1)){
											//REMOVE ENTRY FOR ORTH1 FROM HASHMAP
											logger.info(g1+" has an ortholog "+orth1+" and "+gene2+" which are not overlapping isoforms. Ignoring entry for "+g1+"\t"+orth1);
											orthologs.remove(g1);
										}
										//PREVIOUS ORTHOLOG OF one gene in GENE1 AND GENE2 are ISOFORMS
										else{
											//IGNORE GENE1 ENTRY
											logger.info(g1+" in "+gene1+" has an ortholog "+orth1+" which is an isoform of "+gene2+".Ignore "+g1+"\t"+gene2);
										}
									}
									//IF ORTH2 AND G1 ARE NOT ISOFORMS
									if (!passTwoGeneIsoformTest(annotations1,g1,orth2)){
										//REMOVE ENTRY FOR ORTH1 FROM HASHMAP
										logger.info(gene2+" has an ortholog "+orth2+" and "+g1+" which are not overlapping isoforms. Ignoring entry for "+orth2+"\t"+gene2);
										orthologs.remove(orth2);
									}
									//PREVIOUS ORTHOLOG OF GENE2 AND ONE gene in GENE1 are ISOFORMS
									else{
										//IGNORE GENE1 ENTRY
										logger.info(gene2+" has an ortholog "+orth2+" which is an isoform of "+g1+".Ignore "+g1+"\t"+gene2);
									}
								}
								//IGNORE ENTRY FOR GENE2
								logger.info(gene2+" orthologs were already processed. Ignore "+gene1+"\t"+gene2);
							}
						}
						// IF GENE2 AND GENE1 BOTH HAVE MULTIPLE ISOFORMS
						else{
							//IF ALL GENES IN GENE2 ARE NOT OVERLAPPING ISOFORMS
							String[] g2s = gene2.split(",");
							if(!passIsoformTest(annotations2,g2s)){
								//IGNORE ENTRY
								logger.info("All genes in "+gene2+" are not overlapping isoforms. Ignore "+gene1+"\t"+gene2);
							}
							//GENE1 AND GENE2 HAVE OVERLAPPING ISOFORMS
							else{
								boolean geneInList=false;
								//FOR EACH GENE G1 IN GENE1
								for(int i=0;i<g1s.length;i++){
									String g1 = g1s[i];
									//IF G1 IS ALREADY IN LIST
									if(orthologs.containsKey(g1)){
										geneInList=true;
										//ORTH1 = ORTHOLOG OF G1 IN LIST
										String orth1 = orthologs.get(g1);
										//FOR EACH GENE G2 IN GENE2
										for(int j=0;j<g2s.length;j++){
											String g2 = g2s[j];
											//IF ORTH1 AND G2 ARE NOT ISOFORMS
											if (!passTwoGeneIsoformTest(annotations2,g2,orth1)){
												//REMOVE ENTRY FOR ORTH1 FROM HASHMAP
												logger.info(g1+" has an ortholog "+orth1+" and "+g2+"in "+gene2+" are not overlapping isoforms. Ignoring entry for "+g1+"\t"+orth1);
												orthologs.remove(g1);
												break;
											}
											else{
												//IGNORE g2
												logger.info(g1+" in "+gene1+" has an ortholog "+orth1+" which is an isoform of "+g2+" in "+gene2+ ".Ignore "+g1+"\t"+g2);
											}
										}
									}
								}
								//FOR EACH GENE G2 in GENE2
								for(int j=0;j<g2s.length;j++){
									String g2 = g2s[j];
									//IF G2 IS ALREADY IN LIST
									if(orthologs.containsValue(g2)){
										geneInList=true;
										//ORTH2 = ORTHOLOG OF G2 IN LIST
										String orth2 = findKeyForValue(orthologs,g2);
										//FOR EACH GENE G1 IN GENE1
										for(int i=0;i<g1s.length;i++){
											String g1 = g1s[i];
											//IF ORTH2 AND G1 ARE NOT ISOFORMS
											if (!passTwoGeneIsoformTest(annotations1,g1,orth2)){
												//REMOVE ENTRY FOR ORTH2 FROM HASHMAP
												logger.info(g2+"in "+gene2+" has an ortholog "+orth2+" and "+g1+"in "+gene1+" are not overlapping isoforms. Ignoring entry for "+orth2+"\t"+g2);
												orthologs.remove(orth2);
												break;
											}
											else{
												//IGNORE g1
												logger.info(g2+" in "+gene2+" has an ortholog "+orth2+" which is an isoform of "+g1+" in "+gene1+ ".Ignore "+g1+"\t"+g2);
											}
										}
									}
								}
								//IF NEITHER LISTS HAVE GENES ALREADY CONSIDERED
								if(!geneInList){
									//Find longest isoform in gene1 and find the closest isoform in gene2
									String g1 = getLongestIsoform(annotations1,gene1.split(","));
									orthologs.put(g1,getClosestInLengthIsoform(annotations2,gene2.split(","),annotations1.get(g1)));
								}
								else{
									//NOTHING
								}
							}
						}
					}
				}
			}
		}
		writeToFile(orthologs,argMap.getOutput());
	}
	
	//FINDS
	private String findKeyForValue(Map<String,String> orthologs,String value){ 
		
		Collection<String> keys = orthologs.keySet();
		//System.out.println(keys.toString());
		Iterator<String> iter = keys.iterator();
		String k;
		while(iter.hasNext()){
			k = iter.next();
			if(value.equals(orthologs.get(k))){
				return k;  
			}
		}
		return null;
	}
	
	private boolean passTwoGeneIsoformTest(BEDFileParser annotations,String gene1, String gene2){
		
		if(annotations.get(gene2)==null || annotations.get(gene1)==null)
			return false;
		else
			return (annotations.isOverlapCompatible(annotations.get(gene1), annotations.get(gene2), minOverlap));
	}
	
	private boolean passIsoformTest(BEDFileParser annotations,String[] genes){
		
		boolean passTest=false;
		for(int i=0;i<genes.length-1;i++){
			for(int j=i+1;j<genes.length;j++){
				passTest = passTwoGeneIsoformTest(annotations,genes[i],genes[j]);
				if(!passTest)
					break;
			}
		}
		return passTest;
	}
	
	private String getLongestIsoform(BEDFileParser annotations,String[] genes){
		
		int maxlen = 0;
		int index=0;
		for(int i=0;i<genes.length-1;i++){
			if(annotations.get(genes[i])!=null){
				if (annotations.get(genes[i]).getTranscriptLength()>maxlen){
					maxlen = annotations.get(genes[i]).getTranscriptLength();
					index = i;
				}
			}
		}
		return genes[index];
	}
	
	private String getClosestInLengthIsoform(BEDFileParser annotations,String[] genes,RefSeqGeneWithIsoforms targetGene){
		
		int targetLen = targetGene.getTranscriptLength();
		int index=0;
		
		for(int i=0;i<genes.length-1;i++){
			if(annotations.get(genes[i])!=null){
				if (Math.abs(annotations.get(genes[i]).getTranscriptLength()-targetGene.getTranscriptLength())<targetLen){
					targetLen = Math.abs(annotations.get(genes[i]).getTranscriptLength()-targetGene.getTranscriptLength());
					index = i;
				}
			}
		}
		return genes[index];
	}
	
	
	private void writeToFile(Map<String,String> orthologs,String outputFileName) throws IOException{
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputFileName));
		
		for(String gene:orthologs.keySet()){
			bw.write(gene+"\t"+orthologs.get(gene)+"\n");
		}
		bw.close();
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException{

		UniqueOrthologs dummy = new UniqueOrthologs(args);
	}

}
