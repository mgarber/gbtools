package broad.pda.agilent;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;

import Jama.Matrix;
import broad.core.annotation.AnnotationReader;
import broad.core.annotation.AnnotationReaderFactory;
import broad.core.annotation.BED;
import broad.core.annotation.BEDReader;
import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.annotation.GFF;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.annotation.GTFFileParser;
import broad.pda.annotation.Locus;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.GeneTools;
import broad.pda.gene.RefSeqGene;
import broad.pda.gene.RefSeqGeneWithIsoforms;
import broad.pda.gene.transcriptsNullModel;


//Created by Moran 2/18/10
//The class supports annotation tasks of probes on the array


public class ArrAnnotUtils {
	
	public static final String USAGE = "Usage: ArrAnnotationUtils TASK=<task_num> <task_args>\n" +
	"\tTasks:\n" +
	"\t1. Take the intersection of a BED annotation set with probes on the array. Require that the EXONS will intersect the array regions: \n\t\t -set1 <File containing transcripts set>  \n\t\t -set2 <File containing array probe set> " + "\n"+
	"\t2. Calc the # of intersecting regions of all subsets: \n\t\t-lstOfSets <File with set's file names> " +"\n"+
	"\t3. Calc the # of intersecting regions of all subsets while considering orientation, sets should be full refseq genes!!: \n\t\t-lstOfSets <File with set's file names> " +"\n"+
	"\t4. Calc for each genomic region (merge all transcripts) in how many sets it is covered: \n\t\t-lstOfSets <File with set's file names> "+ "set2format <[BED], GFF or generic> \n\t\t-outformat <desired output format BED, GFF or generic, default is se1format>" +"\n"+
	"\t5. Calc for each transcript in a bed file in how many sets it is covered: \n\t\t-set1 <inputBedFile> \n\t\t-lstOfSets <File with set's file names> "+" \n\t\t-outformat <desired output format BED, GFF or generic, default is se1format>" +"\n"+
	"\tPickExon. Generate an extended bed for probe design for each transcript in the input bed based on its overlap among different sets: \n\t\t-set1 <inputBedFile> \n\t\t-lstOfSets <File with set's file names> \n\t\t-outPrefix <outFilePrefix> \n\t\t-dirOutput <outDir>" +  "\n"+
	"\tTranscriptSetHeatmap. Generate a transcript vs. set heatmap: \n\t\t-set1 <inputBedFile> \n\t\t-lstOfSets <File with set's file names> " +"\n"+
	"\tTrimBy3Prime. Load a bedfile, print how many isoforms where lost: \n\t\t-set1 <fname> \n\t\t-trim <3primeTrimSize,defualt is not to trim the 5 prime> \n\t\t-out <output File, defualt stdout>"+ "\n"+
	"\tMapName. Map names in set1 to the name of the corrosponding overlapping gene in set2: \n\t\t-set1 <fname> \n\t\t-set2 <fname> \n\t\t-out <output File, defualt stdout>"+"\n"+
	"\tClassifyLincs. Classify lincs to by relation to protein\n\t\t-coding genes: \n\t\t-set <fname> \n\t\t-refSet <fname> \n\t\t-out <output File, defualt stdout>"+"\n"+
	"\tMedianScrPerQuantile. Calc the median score per quntile: \n\t\t-set1 <Bed file. optional, if analyze one file> \n\t\t-setLst <optional file with a list of bed files. optional> \n\t\t-quantile <number of quntiles> \n\t\t-out <outFile> \n\t\t-rpkm <optional field, extra field col (number of extra feild e.g 1-5 for scripture's output); if input is rpkm file>  -filterBy <optional field,extra field col;number of extra feild e.g 1-5 for scripture's output, if input is rpkm file> -lowFilterThreshold <optional field,lower Threshold value to filter by the filterBy col; if input is rpkm file> " +  "\n"+
	"\tIsoformStatistics. For each transcript in the ref file provide statistics of isoforms in different sets and isoform compatability: \n\t\t -ref <inputBedFile> \n\t\t-lstOfSets <File with set's file names> "+" \n\t\t-out <outFile> "+" \n\t\t-knownGenes <BED of known genes to calc distance from> -set1Names <subset1names for comaprison> -set2Names<subset2names for comaprison> \n\t\t -betweenAllCompatibleSet<sets for all vs all compatability comaprison> \n\t\t -TabbedLsts <if subsetLsts are tabbed>"+ "-uniqIsoBed <BED file of unique reference isoforms> -multiExonMode"+ "\n"+
	"\tIsoformSelection. Given a reference bed extract all isoforms that overlap it, as well as isoforms that are comaptible between sets: \n\t\t -ref <inputBedFile> \n\t\t-lstOfSets <File with set's file names> "+" \n\t\t-outPrefix <outFile prefix> "+"  -set1Names <subset1names for comaprison> -set2Names<subset2names for comaprison> \n\t\t -betweenAllCompatibleSet<sets for all vs all compatability comaprison> \n\t\t -TabbedLsts <if subsetLsts are tabbed> \n\t\t -cuffCompareSet <name of CuffCompareOutSet>"+"\n"+
	"\tScoreByMaxCompIsoform. Given a bed of choosen isoforms (such as in IsoformSelection) score them by the comaptible set isoform with the maximal score: \n\t\t -ref <inputBedFile> \n\t\t-lstOfSets <File with set's file names> "+" \n\t\t-outFile <outFile name > "+"\n"+
	"\tSelectComplexIsoPerLoci. Given  a reference bed extract the isoform that has the longest intron chain and overlaps that loci: \n\t\t -ref <inputBedFile> \n\t\t-lstOfSets <File with set's file names> "+" \n\t\t-outFile <outFile name > \n\t\t -cuffCompareSet <name of CuffCompareOutSet>"+"\n"+
	"\tConservedExonPerGene. Given a reference set of a SiPhy output file of a set of exons report for the gene in the input set the best scoring exon (lowest omega): \n\t\t -set1 <BED> \n\t\t-siphy <siphy out> \n\t\t-blThreshold <Branch length threshold> \n"+
	"\tRemoveSingleExons. Given a BED/GTF, remove the single exon trancript form the file. \n\t\t -in <BED/GTF file name> \n\t\t-informat <BED/GTF> \n\t\t-out <BED/GTF file name> \n\t\t-outformat<BED/GTF> \n\t\t sizeT <single exon larger than this value will not be filtered> \n"+ 
	"\tCompareAssembleyAndRefLength. Given a BED/GTF assembly and a reference bed, for every ref transcript that has a comaptible assembly transcript print the transcript lengths and (len(asm)/len(ref)) for 3' and 5' exons. \n\t\t -in <BED/GTF file name> \n\t\t-informat <BED/GTF> \n\t\t -ref <BED> \n\t\t chr <chrZ ; optional, load one chr> \n\t\t-out <outfile name> \n"+ 
	"\tGetK4K36Overlaps. Given a input bed and a k4 bed and k36 bed report all transcripts in the input that have K4K36 support. \n\t\t -in <bed> \n\t\t-k4<bed> \n\t\t-k36<bed> \n\t\t -refGenes <bed of known annotations to filter from .independent out file>\n\t\t-out<outname>  \n"+
	"\tUniqRefIsoformStatistics. For each transcript in the ref file provide statistics of isoforms in different sets and isoform compatability: \n\t\t -ref <inputBedFile> \n\t\t-lstOfSets <File with set's file names> "+" \n\t\t-out <outFile> "+" \n\t\t-knownGenes <BED of known genes to calc distance from> -set1Names <subset1names for comaprison> -set2Names<subset2names for comaprison> \n\t\t -betweenAllCompatibleSet<sets for all vs all compatability comaprison> \n\t\t -TabbedLsts <if subsetLsts are tabbed>"+ "-uniqIsoBed <BED file of unique reference isoforms> -multiExonMode"+ "-Ik4k36 <optional bed> -k4k36 <-optional bed> -uniformityScr <optionalGCT> \n"+
	"\tMakeDetectionByGene. convert a detection heatmap to a by gene detection map. \n\t\t-gct <> \n\t\t-gtf <> \n\t\t -out <> \n"+
	"\tAddToUniqRefIsoformStatistics. Add columns to uniqIsoStat: \n\t\t -inTab<out table from UniqRefIsoformStatistic> \n\t\t -out <outFile> "+" \n\t\t -updateCol <int> \n\t\t -bed <bed with Isoforms with new trait> \n\t\t -title <string>  \n\t\t -byIsoform <optional , for data given by iso > \n\t\t -isoGeneMap <if byIsoform, map isoToGene> \n\t\t -byMax <if byIsoform, select representative iso by max val,otherwise by min> \n\t\t -columnCriteria <if byIsoform, select representative isoform by this crieria> \n\t\t -mode(extendMode/updateMode ; for by isoform run) \n\t\t -outPrefix \n\t\t -extendByTable(true if extenting by a gct rather than a bed)\n"+
	"\tParseTransmap. Given a bed file of annotations find transmap transcripts that overlap it. \n\t\t-in <BED> -refgene <transmapReggene.bed mapping> \n\t\t-est <transmapEST.bed mapping> \n\t\t-mRNA <transmapmRNA.bed mapping> \n\t\t -ucsc <transmapmUCSC.bed mapping> \n\t\t -ucscCodingMap <kgTxInfo.tab>\n\t\t -transcriptSpecieMapF <transmapId\tspecies> \n\t\t -chainAlnTab <netChainFullTable> \n\t\t -chainSpecie <e.g mm9> \n\t\t -chainSpecieTranscripts \n\t\t -specieNameFile <sp build\tname>\n\t\t-out \n"+
	"\tmakeSpeciesHeatmap. Given a bed file of annotations generate a heatmap between transcripts and species marking a mapping between the 2. \n\t\t-in <BED> -refgene <transmapReggene.bed mapping> \n\t\t-est <transmapEST.bed mapping> \n\t\t-mRNA <transmapmRNA.bed mapping> \n\t\t -ucsc <transmapmUCSC.bed mapping> \n\t\t  -transcriptSpecieMapF <transmapId\tspecies> \n\t\t -out \n"+
	"\tparseFSAAlignemnt. Given an FSA batch outfile and a pairfile give statistics on alignments.\n\t\t -aln <fsa batch out f> \n\t\t -pair <> \n\t\t -stockholm <if stockholm format>\n\t\t -out <> \n"+
	"\tParseTFbind. given a bed file of lincRNAs and a list of chipseq peaks, match peaks to lincRNAs promoters. \n\t\t -in<BED> \n\t\t -outTable <name> \n\t\t -peakLst <tab file with peak file apth and set name> -upstream <int, expend utr> -downstream <int, expend utr> -tconsXlocMap <mapping between xlocs and tcons>  \n "+
	"\tLevel2Analysis. Parse CSF,PFAM and Pseudogene attributes from a Level 2 lincRNA catalog \n"+
	"\tselectRandomNonAnnotatedBlocks . Select B blocks of size S fro the unannoated portion of the genome (walk through genes and select a random block between them). \n\t\t -setSize <int> \n\t\t -blockSize <int> \n\t\t -annotList <file with a list of annotation bed files> -out<fname> \n"+
	"\tMakeGWASTable. extract meta table and expression table for all lincs in GWAS- no gene LD regions \n\t\t -metaBed <metaBed> \n\t\t -expMat <gct> \n\t\t -gwasTab <> \n\t\t -detectionMat <> \n\t\t -hg18Bed \n\t\t -outprefix \n\t\t-EnhancerBed \n"+
	"\t level1Level2BugFix. given a level1 and level2 GTF files removes from level1 the genes that overlap level2 and are not in refseq. Generate the 2 complementry output files.  \n\t\t -l1gtf <fname> \n\t\t -l2gtf <fname> \n\t\t -refseqSrcStr <hg19_refGene>\n"+ 
	"\t selectMostComplexTranscriptPerGeneId . given a gtf select the longest/ max isoform transcript per gene id . \n\t\t -gtf \n\t\t -maxExNum <optional, else use longest transcript> \n\t\t -out \n"+
	"\t enhancerEnrichment. calc enrichment with enhancer regions; count how many of ref set overlap in their exons with enhancer segments. \n\t\t -enh <BED> \n\t\t-set <ref set BED> \n\t\t -annotTofilter <annotations to filter> \n\t\t -centromers <centromer location> \n\t\t -permPval <num of rand perm>\n"+
	"\t classifyEnhStates. Given Encode enh regions and all HMM calles, find enh regions that are classified as enh by the majority. \n\t\t -enhRegion \n\t\t -hmmList \n\t\t -out \n\t\t -enhStateFile \n\t\t -null state " +
	"\t GenerateIsoformsTissueMat. \n\t\t -refSetGTF <GTF> \n\t\t -setLstFile<file name>  \n\t\t  -isTab  \n\t\t -multiExonMode  \n\t\t -out \n" +
	"\t divergentK4Overlap. \n\t\t -neighborTab  \n\t\t -lincBed  \n\t\t -refGeneBed \n\t\t -k4bed \n\t\t -out \n " +
	"\t mapNames. Concat the name field of transcripts of set 1 with the names of the first overlapping transcript in set2 \n\t\t -set1 <bed> \n\t\t -set2 <bed> \n\t\t -defOrder <boolean; 1 if set1 name is first> \n\t\t -out <bed>  \n" +
	"\t extractLongestIsoFromGTF. Given a gtf and a list of gene Ids , extract the longest transcript for this gene and output as GTF \n\t\t -set1 <gtf> \n\t\t -geneNames <lst> \n\t\t -geneMap<tabTxt>  \n\t\t -out <gtf>  \n" +
	"\t BestMatchCuffdiffSegmentPerGene. \n\t\t -gtf \n\t\t -cuffdiff <genes.diff out file> \n\t\t -segmentSizeThreshold <minimal fragment size to consider> \n\t\t  -out <outf>. Given a gtf a segmentSizeThreshold and a cuffdiff out file, report for each gene name the highest fold change segment (that has a significant enrichment, and over certain size), report an additional file , with all segments per gene\n"+
	"\t ReannotateGTF.  Given a GTF(e.g. cuffmerge output) add GTF attributes from new GTF files. To annotate by several ref GTFs , run this command sequantialy. \n\t\t -inGTF (gtf) \n\t\t -annotGTF (gtf) \n\t\t -annotPrefix (str)\n\t\t -out (str) -src(str) \n" +
	"\t ReannotateGTFfromTable.  Given a GTF(e.g. cuffmerge output) add GTF attributes from a table using the same ids. To annotate by several refs , run this command sequantialy. \n\t\t -inGTF (gtf) \n\t\t -annotTAB (gtf) \n\t\t -annotPrefixList (comme Seperated prefix list of fields to add)\n\t\t -annotHeaderList (comme Seperated header list of fields to add)\n\t\t -headerIdenitifier (header name to idenitfy the transcript by) \n\t\t  -attrIdenitifier (gene attr name to idenitfy the transcript by) \n\t\t -out (str)  \n" +
	"\t CapNumIsoformsPerLoci. given a GTF select the first  'maxIsoPerLoci' isofrms for genes with a higher isoform number. \n\t\t inGTF \n\t\t outf \n\t\t maxIsoPerLoci \n\t\t src"+
	"\n";
	
	
	public static void main (String [] args) throws IOException, ParseException {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE);
		
		if ("1".equals(argMap.getTask())) {	//calc intersection of set with array probes
			String set1In = argMap.getMandatory("set1");
			String set2In = argMap.getMandatory("set2");
			String set2Format = argMap.containsKey("set2format") ? argMap.get("set2format") : "BED";
			
			
			AnnotationReader<? extends GenomicAnnotation> set2 = AnnotationReaderFactory.create(set2In, set2Format);
			BEDReader set1= new BEDReader(set1In);
			
			//change this : write a function that takes annotation readers and writes
			// the first annotation an the name 
			BufferedWriter bw = argMap.getOutputWriter();
			writeIntersectionToExtendedBed(set2,set1.getAnnotationList(),bw);
			bw.close();
			}
		else if ("2".equals(argMap.getTask()))//calc overlap between sets
		{
			String setLst = argMap.getMandatory("lstOfSets");
			String setFormat = argMap.containsKey("setsformat") ? argMap.get("setsformat") : "BED";
			BufferedWriter bw = argMap.getOutputWriter();
			bw.write("SetName\t"+"SetOverlapWithSuperSet"+"\t"+"SuperSetOverlapWithSet"+"\t"+"SetSize"+"\t"+"SuperSetSize"+"\t"+"AdditionOfNewSet"+"\n");
			ArrayList <String> nameLst= new ArrayList <String> ();
			//read to set of sets
			Map <String,AnnotationReader<? extends GenomicAnnotation>> setOfSets= uploadFiles(setLst,setFormat,nameLst);
			//merged set= Union set
			AnnotationReader<? extends GenomicAnnotation> superSet = AnnotationReaderFactory.create(setFormat);
			//Iterator<String> setsIt = setOfSets.keySet().iterator();
			Iterator<String> setsIt = nameLst.iterator();
			int newSetSize=0;
			while(setsIt.hasNext()) {
				String currName=setsIt.next();
				AnnotationReader<? extends GenomicAnnotation> currSet= setOfSets.get(currName);
				//increase subset by one set on calc overlap with Union set, write to file
				currSet.merge(); //TODO - MAKE MERGE ORIENTAION SENSATIVE
				int l1= currSet.getOverlappers(superSet.getAnnotationList()).size();
				int l2= superSet.getOverlappers(currSet.getAnnotationList()).size();
				int s1=currSet.getAnnotationList().size();
				int s2= superSet.getAnnotationList().size();
				superSet.concatAnnotationSet(0, currSet.getAnnotationList());
				superSet.merge();
				newSetSize=superSet.getAnnotationList().size();;
				bw.write(currName+"\t"+l1+"\t"+l2+"\t"+s1+"\t"+s2+"\t"+(newSetSize-s2)+"\n");
				
			}
			
			bw.write("final super set size : "+ newSetSize+ "\n");
			bw.close();
					
		}
		
		else if ("2".equals(argMap.getTask()))//calc overlap between sets
		{
			String setLst = argMap.getMandatory("lstOfSets");
			String setFormat = argMap.containsKey("setsformat") ? argMap.get("setsformat") : "BED";
			BufferedWriter bw = argMap.getOutputWriter();
			calcOverlapOfMultipleSets(setLst,setFormat,bw);	
			bw.close();
		}
		
		else if ("3".equals(argMap.getTask()))//calc overlap between sets while considering orientation
		{
			String setLst = argMap.getMandatory("lstOfSets");
			BufferedWriter bw = argMap.getOutputWriter();
			calcOrientedOverlapOfMultipleSets(setLst,bw);	
			bw.close();
					
		}
		
		
		else if ("4".equals(argMap.getTask()))//Number Of occurrences In Sets
		{
			String setLst = argMap.getMandatory("lstOfSets");
			String setFormat = argMap.containsKey("setsformat") ? argMap.get("setsformat") : "BED";
			BufferedWriter bw = argMap.getOutputWriter();
			
			//read to set of sets
			ArrayList <String> nameLst= new ArrayList <String> ();
			Map <String,AnnotationReader<? extends GenomicAnnotation>> setOfSets= uploadFiles(setLst,setFormat,nameLst);
			//merged set= Union set
			AnnotationReader<? extends GenomicAnnotation> superSet = AnnotationReaderFactory.create(setFormat);
			Iterator<String> setsIt = setOfSets.keySet().iterator();
			while(setsIt.hasNext()) {
				String currName=setsIt.next();
				AnnotationReader<? extends GenomicAnnotation> currSet= setOfSets.get(currName);
				currSet.merge();
				superSet.concatAnnotationSet(0, currSet.getAnnotationList());
				superSet.merge();
			}
			setsIt = setOfSets.keySet().iterator();
			Boolean setToZero=true;
			while(setsIt.hasNext()) {
				String currName=setsIt.next();
				AnnotationReader<? extends GenomicAnnotation> currSet= setOfSets.get(currName);
				superSet.IncrementScoreIfOverlap(currSet, 0, setToZero);
				setToZero=false;
			}
						
			//write superset
			writeAnnotationReader(superSet,bw);
			bw.close();		
		}
		
		else if ("5".equals(argMap.getTask()))//Number Of occurrences In Sets , using input bed
		{
			String setLst = argMap.getMandatory("lstOfSets");
			String refSet = argMap.getMandatory("set1");
			String setFormat = argMap.containsKey("setsformat") ? argMap.get("setsformat") : "BED";
			BufferedWriter bw = argMap.getOutputWriter();
			
			//read to set of sets
			ArrayList <String> nameLst= new ArrayList <String> ();
			Map <String,AnnotationReader<? extends GenomicAnnotation>> setOfSets= uploadFiles(setLst,setFormat,nameLst);
			AnnotationReader<? extends GenomicAnnotation> superSet = AnnotationReaderFactory.create(refSet, setFormat); 
			
			superSet =CountOverlaps1 (superSet, setOfSets);
				
			//write superset
			writeAnnotationReader(superSet,bw);
			bw.close();		
		}
		else if ("TranscriptSetHeatmap".equals(argMap.getTask()))//Number Of occurrences In Sets , using input bed
		{
			String setLst = argMap.getMandatory("lstOfSets");
			String refSet = argMap.getMandatory("set1");
			BufferedWriter bw = argMap.getOutputWriter();
			
			//read to set of sets
			ArrayList <String> nameLst= new ArrayList <String> ();
			Map <String,AnnotationReader<? extends GenomicAnnotation>> setOfSets= uploadFiles(setLst,"BED",nameLst);
			AnnotationReader<? extends GenomicAnnotation> superSet = AnnotationReaderFactory.create(refSet, "BED"); 
			
			MatrixWithHeaders heatmap =GenerateHeatmap (superSet, setOfSets, nameLst);
				
			//write superset
			heatmap.writeGCT(bw);
			bw.close();		
		}
		else if ("PickExon".equals(argMap.getTask())) //choose exon to design a probe set to
		{
			String setLst = argMap.getMandatory("lstOfSets");
			String outPrefix = argMap.getMandatory("outPrefix");
			String outDir = argMap.getMandatory("dirOutput");
			String refSetName = argMap.containsKey("set1") ? argMap.get("set1") : "";
			BufferedWriter bw = argMap.getOutputWriter();
			String trim = argMap.containsKey("trim") ? argMap.get("trim") : "";
			int probeSize=60;  /// can change to be as a parameter, but set of Agilent
			//read ref set and all dataset
			BEDFileParser refSet;
			if (! refSetName.equalsIgnoreCase(""))
				refSet= new BEDFileParser(refSetName);
			else
				refSet=null;
			ExtendedBedForBrobeDesign (refSet,setLst,outPrefix,outDir,bw,trim,probeSize);
			bw.close();
			
		}
		else if ("TrimBy3Prime".equals(argMap.getTask())) 
		{
			String inf = argMap.getMandatory("set1");
			String trimSize = argMap.containsKey("trim") ? argMap.get("trim") : "";
			BEDFileParser bed;
			if (trimSize.equalsIgnoreCase(""))
				bed= new BEDFileParser(inf);
			else
				bed=new BEDFileParser(inf,Integer.valueOf(trimSize));
			
			BufferedWriter bw = argMap.getOutputWriter();
			bed.writeFullBed(bw);
			bw.close();
			
		}
		else if ("MapName".equalsIgnoreCase(argMap.getTask())) {
			String set1In = argMap.getMandatory("set1");
			String set2In = argMap.getMandatory("set2");
			BufferedWriter bw = argMap.getOutputWriter();
			BEDFileParser bed= MapNames(set1In,set2In);
			bed.writeFullBed(bw);
			bw.close();
			
			
		}
		else if ("ClassifyLincs".equalsIgnoreCase(argMap.getTask())) {
			String lincsFile = argMap.getMandatory("set");
			String refGenesFile = argMap.getMandatory("refSet");
			String bufferSize = argMap.getMandatory("bufferSize");
			String outName = argMap.getMandatory("output");
		
			//BEDFileParser bed= ClassifyLincs(lincsFile,refGenesFile,Integer.valueOf(bufferSize));
			Map<RefSeqGene, String> lincs= ClassifyLincs(lincsFile,refGenesFile,Integer.valueOf(bufferSize));
			BEDFileParser.writeFullBED(outName, lincs) ;
		
		}
		
		else if ("MedianScrPerQuantile".equalsIgnoreCase(argMap.getTask())){
			int q=(int)(Double.valueOf(argMap.getMandatory("quantile"))).doubleValue();
			String file = argMap.containsKey("set1") ? argMap.get("set1") : "";
			String lstfile = argMap.containsKey("setLst") ? argMap.get("setLst") : "";
			int rpkmField = argMap.containsKey("rpkm") ? (int)(Double.valueOf(argMap.get("rpkm")).doubleValue()) : 0;
			int filterField = argMap.containsKey("filterBy") ? (int)(Double.valueOf(argMap.get("filterBy")).doubleValue()) : 0;
			double lowT = argMap.containsKey("UpFilterThreshold") ? (Double.valueOf(argMap.get("UpFilterThreshold")).doubleValue()) : 0;
			BufferedWriter bw = argMap.getOutputWriter();
			getMedianScrPerQauntile(file,lstfile,q,bw,rpkmField,filterField,lowT);
			bw.close();
		}
		
		//
		else if ("IsoformStatistics".equals(argMap.getTask()))
		{
			String setLst = argMap.getMandatory("lstOfSets");
			String refSet = argMap.getMandatory("ref");
			String geneFile = argMap.getMandatory("knownGenes");
			String lst1= argMap.getMandatory("set1Names");
			String lst2= argMap.getMandatory("set2Names");
			String lst3= argMap.getMandatory("betweenAllCompatibleSet");
			Boolean isTab=argMap.containsKey("TabbedLsts") ? true: false;
			Boolean isMultiExon=argMap.containsKey("multiExonMode") ? true: false;
			//String excludeFromCompatible=argMap.containsKey("excludeFromCompatible") ?argMap.get("excludeFromCompatible"): "";
			String uniqIsoBed= argMap.containsKey("uniqIsoBed")? argMap.get("uniqIsoBed") : "";
			BufferedWriter bw = argMap.getOutputWriter();
			GenerateIsoformsStatistics (refSet,setLst,geneFile,bw,lst1,lst2,isTab,lst3,uniqIsoBed,isMultiExon);
			
			bw.close();		
		}
		else if ("UniqRefIsoformStatistics".equals(argMap.getTask()))
		{
			String setLst = argMap.getMandatory("lstOfSets");
			String refSet = argMap.getMandatory("ref");
			String geneFile = argMap.getMandatory("knownGenes");
			String lst1= argMap.getMandatory("set1Names");
			String lst2= argMap.getMandatory("set2Names");
			String lstref= argMap.getMandatory("setExternalNames");
			String lst3= argMap.getMandatory("betweenAllCompatibleSet");
			Boolean isTab=argMap.containsKey("TabbedLsts") ? true: false;
			Boolean isMultiExon=argMap.containsKey("multiExonMode") ? true: false;
			String uniqIsoBed= argMap.containsKey("uniqIsoBed")? argMap.get("uniqIsoBed") : "";
			String Ik4k36= argMap.containsKey("ik4k36")? argMap.get("ik4k36") : "";
			String k4k36= argMap.containsKey("k4k36")? argMap.get("k4k36") : "";
			String uniformityScr= argMap.containsKey("uniformityScr")? argMap.get("uniformityScr") : "";
			String cov= argMap.containsKey("cov")? argMap.get("cov") : "";
			String scanStatF= argMap.containsKey("scanStat")? argMap.get("scanStat") : "";
			
			BufferedWriter bw = argMap.getOutputWriter();
			GenerateIsoformsStatisticsByUniqRef (refSet,setLst,geneFile,bw,lst1,lst2,isTab,lst3,lstref,uniqIsoBed,isMultiExon,Ik4k36,k4k36,uniformityScr,cov,scanStatF);
			
			bw.close();		
		}
		else if ("AddToUniqRefIsoformStatistics".equals(argMap.getTask())){
			 String inTab =  argMap.getMandatory("inTab");
			 BufferedWriter bw = argMap.getOutputWriter();
			 String updateCol = argMap.getMandatory("updateCol");
			 String bed = argMap.getMandatory("bed");
			 String title = argMap.getMandatory("title");
			 //byIsoform,isoGeneMap,byMax,columnCriteria
			 if (argMap.containsKey("byIsoform")){ //adds to the table several columns by choosing the isoforms some max/min criteria on a specific extra column
				 String isoGeneMap= argMap.getMandatory("isoGeneMap");
				 boolean byMaxCriteria= argMap.containsKey("byMax")? true : false;
				 boolean sumCriteria= argMap.containsKey("sumCriteria")? true : false;
				 String columnCriteria= argMap.getMandatory("columnCriteria");
				 String mode= argMap.containsKey("mode")? argMap.getMandatory("mode"): "extendMode";
				 boolean extendByBED = argMap.containsKey("extendByTable")? false: true;
				 AddToUniqRefIsoformStatistics(isoGeneMap,byMaxCriteria,columnCriteria, inTab,title,bed,bw,mode,extendByBED,sumCriteria);
			 }
			 else
				 AddToUniqRefIsoformStatistics(inTab,updateCol,title,bed,bw);
			 bw.close();
		}
		else if ("MakeDetectionByGene".equals(argMap.getTask())){
			
			String expMat = argMap.getMandatory("gct");  
			String refset = argMap.getMandatory("refGTF");
			BufferedWriter bw = argMap.getOutputWriter();
			makeDetectionByGene (expMat,refset,bw);
			bw.close();
		}
		else if("GeneStatistics".equals(argMap.getTask()))
		{
			String neighborMap = argMap.getMandatory("neighborMap");
			String expMat = argMap.getMandatory("expMat");
			String isoformStatMat = argMap.getMandatory("isoformStatMat");
			String refSet = argMap.getMandatory("refGTF");
			String JSgct = argMap.getMandatory("JSgct");
			String k4k36lincs = argMap.getMandatory("k4k36lincs");
			BufferedWriter bw = argMap.getOutputWriter();
			GenerateGeneStatMat(refSet,neighborMap, expMat, isoformStatMat,JSgct,k4k36lincs,bw );
			bw.close();
			
		}
		else if ("IsoformSelection".equals(argMap.getTask())){

			String setLst = argMap.getMandatory("lstOfSets");
			String refSet = argMap.getMandatory("ref");
			String lst1= argMap.getMandatory("set1Names");
			String lst2= argMap.getMandatory("set2Names");
			String lst3= argMap.getMandatory("betweenAllCompatibleSet");
			Boolean isTab=argMap.containsKey("TabbedLsts") ? true: false;
			String cuffcmprSet= argMap.getMandatory("cuffCompareSet");
			String prefix= argMap.getMandatory("outPrefix");
			ExtractIsoforms (refSet,setLst,lst1,lst2,isTab,lst3,cuffcmprSet,prefix);
			
		}
		
		else if ("ScoreByMaxCompIsoform".equals(argMap.getTask())){
			String setLst = argMap.getMandatory("lstOfSets");
			String refSet = argMap.getMandatory("ref");
			String outFile= argMap.getMandatory("outFile");
			ScoreByMaxScoredCompatibleIsoform (refSet,setLst,outFile,true);
		}
		else if ("SelectComplexIsoPerLoci".equals(argMap.getTask())){
			String setLst = argMap.getMandatory("lstOfSets");
			String refSet = argMap.getMandatory("ref");
			String outFile= argMap.getMandatory("outFile");
			String cuffcmprSet= argMap.getMandatory("cuffCompareSet");
			SelectLongestIntornChainIsoPerLoci (refSet,setLst,outFile,cuffcmprSet);
		}
		
		else if ("ConservedExonPerGene".equals(argMap.getTask())){
			String set = argMap.getMandatory("set1");
			String siphy = argMap.getMandatory("siphy");
			double bl= argMap.getDouble("blThreshold");
			BufferedWriter bw = argMap.getOutputWriter();
			GetConservedExonScrPerGene (set,siphy,bl,bw);
			bw.close();
		}
		
		else if ("RemoveSingleExons".equals(argMap.getTask())){
			String in = argMap.getInput();
			String informat= argMap.getMandatory("informat");
			String outformat= argMap.getMandatory("outformat");
			String src=  argMap.containsKey("src")? argMap.get("src"): "unknown";
			double sizeT=argMap.getDouble("sizeT");
			String out= argMap.getOutput();
			BufferedWriter bw = new BufferedWriter(new FileWriter(out));
			removeSingleExons (in,informat,outformat,bw,src,(int)sizeT);
			bw.close();
		}
		
		else if ("CompareAssembleyAndRefLength".equals(argMap.getTask())){
			String in = argMap.getInput();
			String informat= argMap.getMandatory("informat");
			String ref=argMap.getMandatory("ref");
			BufferedWriter bw  = argMap.getOutputWriter();
			String chr=argMap.containsKey("chr")? argMap.get("chr"): "";
			compareAssembleyAndRefLength(in,informat,ref,bw,chr);
			bw.close();
		}
		else if("GetK4K36Overlaps".equals(argMap.getTask())){
			String in = argMap.getInput();
			String k4= argMap.getMandatory("k4");
			String k36=argMap.getMandatory("k36");
			String genesF=argMap.getMandatory("refGenes");
			Integer promoterInterval = argMap.containsKey("promoterInterval")? argMap.getInteger("promoterInterval"): 2000;
			Integer randPerm = argMap.containsKey("permPval")? argMap.getInteger("permPval"):0; 
			String outPre=argMap.getOutput();
			String[] infiles=new String[3];
			if (randPerm > 0){
					infiles[2] = argMap.getMandatory("chrSizes");
					infiles[0]  = argMap.getMandatory("annotTofilter"); 
					infiles[1] = argMap.getMandatory("centromers"); 
		    }
			GetK4K36Overlaps(in,k4,k36,outPre,promoterInterval,genesF,randPerm,infiles);
		}
		else if ("ParseTransmap".equals(argMap.getTask())){
			String in = argMap.getInput();
			String refgeneF = argMap.getMandatory("refgene");
			String estF = argMap.getMandatory("est");
			String mrnaF = argMap.getMandatory("mRNA");
			String ucscF = argMap.getMandatory("ucsc");
			String ucscCodingMap  = argMap.getMandatory("ucscCodingMap");
			String transcriptSpecieMapF= argMap.getMandatory("transcriptSpecieMap");
			String ChainTabF= argMap.getMandatory("chainAlnTab");
			String ChainSpecie =  argMap.getMandatory("chainSpecie");
			String outPrefix = argMap.getMandatory ("outPrefix");
			String chainSpecieTranscripts = argMap.getMandatory ("chainSpecieTranscripts");
			BufferedWriter bw  = argMap.getOutputWriter();
			ParseTransmap(in,refgeneF,estF,mrnaF,ucscF,ucscCodingMap,transcriptSpecieMapF, ChainTabF, ChainSpecie,bw,outPrefix,chainSpecieTranscripts);
			bw.close();
			
		}
		else if ("makeSpeciesHeatmap".equals(argMap.getTask())){
			String in = argMap.getInput();
			String refgeneF = argMap.getMandatory("refgene");
			String estF = argMap.getMandatory("est");
			String mrnaF = argMap.getMandatory("mRNA");
			String ucscF = argMap.getMandatory("ucsc");
			String transcriptSpecieMapF= argMap.getMandatory("transcriptSpecieMap");
			String specieNameFile= argMap.getMandatory("specieNameFile");
			BufferedWriter bw  = argMap.getOutputWriter();
			makeSpHeatMap(in,refgeneF,estF,	mrnaF,ucscF,transcriptSpecieMapF,specieNameFile,bw);
			bw.close();
			
		}
		else if ("Level2Analysis".equals(argMap.getTask())){
		    String in =argMap.getInput();
		    BufferedWriter bw  = argMap.getOutputWriter();
		    Level2Analysis(in,bw);
		    bw.close();
		}
		else if ("ParseTFbind".equalsIgnoreCase(argMap.getTask())){
			String in = argMap.getInput();
			String outStr  = argMap.getOutput();
			String listF = argMap.getMandatory("peakLst");
			Integer upstream = argMap.getInteger("upstream");
			Integer downstream = argMap.getInteger("downstream");
			String tconsXlocMap = argMap.getMandatory("tconsXlocMap");
			boolean multimap =  argMap.containsKey("singleFile")? false :true;
			PraseTFbind (in,outStr,listF,upstream,downstream,tconsXlocMap,multimap);
		}
		else if ("parseFSAAlignemnt".equalsIgnoreCase(argMap.getTask())){
			String aln=argMap.getMandatory("aln");
			String pair =argMap.getMandatory("pair");
			String isoGeneMapF =argMap.getMandatory("isoGeneMap");
			boolean format = argMap.containsKey("stockholm")? true:false;
			String outprefix  = argMap.getOutput();
			parseFSAAlignment(aln,pair,outprefix,isoGeneMapF,format);
			
			
		}
		else if ("selectRandomNonAnnotatedBlocks".equalsIgnoreCase(argMap.getTask())){
			Integer setSize=argMap.getInteger("setSize");
			Integer blockSize =argMap.getInteger("blockSize");
			String fileList=argMap.getMandatory("annotList");
			String outprefix  = argMap.getOutput();
			selectRandomNonAnnotatedBlocks(setSize,blockSize,fileList,outprefix);
			
		}
		else if ("MakeGWASTable".equalsIgnoreCase(argMap.getTask())){
			String metaBed = argMap.getMandatory("metaBed");
			String expMat  = argMap.getMandatory("expMat");
			String gwasTab = argMap.getMandatory("gwasTab");
			String detectionMat = argMap.getMandatory ("detectionMat"); 
			String buildMapFile = argMap.getMandatory("buildMapBed");
			String outprefix = argMap.getMandatory("outprefix");
			String filterF = argMap.containsKey("calcEnrich")? argMap.getMandatory("filterF"):"";
			String EnhancerBed =argMap.containsKey("EnhancerBed")? argMap.getMandatory("EnhancerBed"):"";
			String snpPos = argMap.containsKey("snpPos")? argMap.getMandatory("snpPos"):"";
 			String snpsInLD =  argMap.containsKey("snpsInLD")? argMap.getMandatory("snpsInLD"):"";
			MakeGWASTable (metaBed,expMat,gwasTab,detectionMat,outprefix,buildMapFile,filterF,EnhancerBed,snpPos,snpsInLD);
		}
		else if("level1Level2BugFix".equalsIgnoreCase(argMap.getTask())){
			String l1gtf = argMap.getMandatory("l1gtf");
			String l2gtf = argMap.getMandatory("l2gtf");
			String refseqSrcStr = argMap.getMandatory("refseqSrcStr");
			String outpre =argMap.getMandatory("outpre");
			level1Level2BugFix (l1gtf,l2gtf,refseqSrcStr,outpre+"L1_",outpre+"L2_");
			
		}
		else if ("selectMostComplexTranscriptPerGeneId".equalsIgnoreCase(argMap.getTask())){
			
			String GTF=  argMap.getMandatory("gtf");
			String out = argMap.getOutput();
			selectMostComplexTranscriptPerGeneI(GTF,out);
		}
		else if ("enhancerEnrichment".equalsIgnoreCase(argMap.getTask())){
			String set= argMap.getMandatory("set");
			String enh = argMap.getMandatory("enh");
			String informat =  argMap.containsKey("informat")? argMap.get("informat"):"gtf";
			if (argMap.containsKey("permPval")){
				Integer permNum = argMap.getInteger("permPval");
				String centromers =argMap.getMandatory("centromers");
				String chrSizes =argMap.getMandatory("chrSizes");
				String annotFilter =argMap.getMandatory("annotTofilter");
				String outPrefix = argMap.getMandatory("outPrefix");
				
				overlapEnrichment (set,enh,permNum,centromers,annotFilter,outPrefix,chrSizes,informat);
			}
			else
				overlapEnrichment (set,enh,0,null,null,null,null,informat);
				
		}
		
		else if ("classifyEnhStates".equalsIgnoreCase(argMap.getTask())) {
				String enhRegions = argMap.getMandatory("enhRegions");
				String hmmList = argMap.getMandatory("hmmList");
				String outPrefix = argMap.getMandatory("outprefix");
				String enhStateFile = argMap.getMandatory("enhStateFile");
				String nullstate = argMap.getMandatory("nullstate");
				Integer binSize = argMap.containsKey("binSize")? argMap.getInteger("binSize") : 200;
				classifyEnhStates( enhRegions, hmmList, outPrefix, enhStateFile, nullstate,binSize);
		}
		else if ("GenerateIsoformsTissueMat".equalsIgnoreCase(argMap.getTask())){
			String refSetGTF =  argMap.getMandatory("refSetGTF");
			String setLstFile =  argMap.getMandatory("setLstFile");
			 Boolean multiExonMode = argMap.containsKey("multiExonMode")? true : false;
			String out  = argMap.getOutput();
			GenerateIsoformsTissueMat(refSetGTF ,setLstFile, true,multiExonMode,out,true) ; 

		}
		else if ("divergentK4Overlap".equalsIgnoreCase(argMap.getTask())){
			String neighborTab =  argMap.getMandatory("neighborTab");
			String lincBed =  argMap.getMandatory("lincBed");
			String refGeneBed =  argMap.getMandatory("refGeneBed");
			String k4bed =  argMap.getMandatory("k4bed");
			String out =  argMap.getMandatory("outf");
			
			analyzeDivergentK4Overlap (neighborTab, lincBed, refGeneBed, k4bed,out);
		}
		
		else if ("mapNames".equalsIgnoreCase(argMap.getTask())) {
			BEDFileParser set1 = new BEDFileParser( argMap.getMandatory("set1"));
			BEDFileParser set2In = new BEDFileParser (argMap.getMandatory("set2"));
			boolean defOrder =  argMap.containsKey("defOrder")? new Boolean (argMap.getMandatory("defOrder") ): true;
			String out = argMap.getOutput();
			BEDFileParser outbed = new BEDFileParser ();
			
			List<RefSeqGene> set1genes = set1.GetGenes();
			for (RefSeqGene g: set1genes){
				IntervalTree<RefSeqGeneWithIsoforms> overlaps= set2In.getOverlappers(g);
				if (!overlaps.isEmpty()){
					Iterator <RefSeqGeneWithIsoforms> it= overlaps.valueIterator();
					RefSeqGeneWithIsoforms og = it.next();
					String name1 = g.getName();
					String name2 = og.getName();
					String newName = name1+"_"+name2;
					if (!defOrder)
						newName = name2+"_"+name1;
					g.setName(newName);
					outbed.addRefSeq(g);
				}
			}
			outbed.writeFullBed(out);
			
			
			
			
		}
			
		else if ("extractLongestIsoFromGTF".equalsIgnoreCase(argMap.getTask())){
			String gtf =  argMap.getMandatory("gtf");
			String geneLst=  argMap.getMandatory("geneLst");
			String geneMap =  argMap.getMandatory("GeneTconMap");
			String out =  argMap.getMandatory("outf");
			
			extractLongestIsoFromGTF (gtf,geneLst,geneMap, out);
		}

		else if ("BestMatchCuffdiffSegmentPerGene".equalsIgnoreCase(argMap.getTask())){
			String gtf =  argMap.getMandatory("gtf");
			String diff=  argMap.getMandatory("cuffdiff");
			String out =  argMap.getMandatory("outf");
			Integer segmentSizeThreshold = argMap.getInteger("segmentSizeThreshold");
			
			BestMatchCuffdiffSegmentPerGene (gtf,diff, out,segmentSizeThreshold);
		}
	
		else if ("ReannotateGTF".equalsIgnoreCase(argMap.getTask())){
			String GTF = argMap.getMandatory("inGTF");
			String annotGTF = argMap.getMandatory("annotGTF");
			String annotPrefix = argMap.getMandatory("annotPrefix");
			String src = argMap.getMandatory("src");
			String outf =  argMap.getMandatory("outf");
			
			ReannotateGTF(GTF,annotGTF,annotPrefix,outf,src);
			
			
		}
		else if ("ReannotateGTFfromTable".equalsIgnoreCase(argMap.getTask())){
			String GTF = argMap.getMandatory("inGTF");
			String annotTab = argMap.getMandatory("annotTAB");
			String annotPrefix = argMap.getMandatory("annotPrefixList");
			String annotHeaderList = argMap.getMandatory("annotHeaderList");
			String  headerIdenitifier = argMap.getMandatory("headerIdenitifier");
			String attrIdenitifier = argMap.getMandatory("attrIdenitifier");
			String outf =  argMap.getMandatory("outf");
			ReannotateGTFfromTable (GTF,annotTab,annotPrefix,annotHeaderList,headerIdenitifier, attrIdenitifier,outf);
		}	
		
		else if ("CapNumIsoformsPerLoci".equalsIgnoreCase(argMap.getTask())){
			String GTF = argMap.getMandatory("inGTF");
			String outf =  argMap.getMandatory("outf");
			Integer maxIsoPerLoci = argMap.getInteger("maxIsoPerLoci");
			String src = argMap.getMandatory("src");
			CapNumIsoformsPerLoci (GTF, outf,maxIsoPerLoci,src);
		}
		
		else {
			System.err.println("Task " + argMap.getTask() + " is invalid or no task was specified");
		}
 }


		





			//Input: Transcript-gene name map.
	//		Uniq reference transcripts (in GTF , so we can generate the mergedByGene Ref).
	//		List of tissue beds.  Scripture ref, Cufflinks ref, External ref, union ref 
	//		K4K36 BED, K4K36 "independent" BED, Uniformity Score, ScanStat BED
	/*Function Description: 1) go over all refs , create Locus set.
	 * 	2) go over every tissue file, place transcript with compatible
	 *  3) Go over the merged Genes , mark the ones that have a compatible iso between 2 sets or in the NRs
	 *  4) report table for each transcript : name of gene, num uniq iso, High Conf, and specific detection info
	
	*/
			private static void GenerateIsoformsStatisticsByUniqRef(String refSetGTF,String setLstFile,String annotatedRefFile, BufferedWriter bw,String set1LstF,String set2LstF,Boolean isTab, String betweenAllCompatibleSet,String externalLsts,String uniqIsoF, Boolean MultiExonMode, String ik4k36F, String k4k36F, String uniformityScrF,String covF, String scanStatBed) throws IOException, ParseException{
			
			BEDFileParser annotatedRef= new BEDFileParser(annotatedRefFile);
			ArrayList<String> set1Lst=loadNameLst(set1LstF,isTab);
			ArrayList<String> set2Lst=loadNameLst(set2LstF,isTab);
			ArrayList<String> set3Lst=loadNameLst(betweenAllCompatibleSet,isTab);
			ArrayList<String> setRefLst=loadNameLst(externalLsts,isTab);
			
			Map<String,Locus> refLociMap=new  HashMap<String,Locus>();
			
			BEDFileParser uniqTconRef=new BEDFileParser();
			uniqTconRef.loadGTF(new File(refSetGTF),"");
			//HashMap<RefSeqGeneWithIsoforms, List<RefSeqGeneWithIsoforms>> geneTconMap=new HashMap<RefSeqGeneWithIsoforms, List<RefSeqGeneWithIsoforms>>();
			//Map<RefSeqGeneWithIsoforms, RefSeqGeneWithIsoforms> TconGeneMap=new HashMap<RefSeqGeneWithIsoforms,RefSeqGeneWithIsoforms>();
			HashMap<String, List<String>> geneTconMap=new HashMap<String, List<String>>();
			Map<String, String> TconGeneMap=new HashMap<String,String>();
			Map<String , RefSeqGeneWithIsoforms> geneNameGeneMap= new HashMap<String , RefSeqGeneWithIsoforms>();
			MatrixWithHeaders uniScr=null;
			MatrixWithHeaders cov=null;
			HashMap<String,Double> k4k36=null;
			HashMap<String,Double> ik4k36=null;
			HashMap<String,Double> scanStat=null;
			
			if (!k4k36F.equals(""))
				k4k36=getNamesFromBed(k4k36F);
			if (!ik4k36F.equals(""))
				ik4k36=getNamesFromBed(ik4k36F);
			if (!scanStatBed.equals(""))
				scanStat=getNamesFromBed(scanStatBed);
			if (! uniformityScrF.equals(""))
				uniScr=new MatrixWithHeaders(uniformityScrF);
			if (! covF.equals(""))
				cov=new MatrixWithHeaders(covF);
				
			
			
			//1.
			BEDFileParser uniqGeneRef= AssociateGeneTranscript(uniqTconRef,geneTconMap,TconGeneMap,geneNameGeneMap);
			
			//2.
			LinkedList<String> setsLst =loadCompatibleToRefLociIsoformMap (uniqTconRef,setLstFile,refLociMap,MultiExonMode,true) ;
			
			//3.
			HashMap<String,ArrayList<Boolean>> comaptibleIso=markComaptibleIso (refLociMap,set1Lst,set2Lst,setRefLst,set3Lst);
			HashMap<String,Boolean> comaptibleGene=markComaptibleGene (uniqGeneRef,comaptibleIso,geneTconMap);
			
			//4. Print:			
			bw.write("chr\tstart\tend\tname\tscore\tstrand\ttStart\ttEnd\trgb\tblockCnt\tblockSizes\tblockStarts");
			bw.write("\tis5'EqOrient\tt5'NeighborDist\tis3'EqOrient\t3'NeighborDist\t");
			bw.write("CompatibleBetweenAssemblies\tCompatibleBetweenAssemblyAndExternal\tCompatibleBetweenSets");
			//number of isoforms in set
			for (String s:setsLst)
				bw.write("\t"+s);
			bw.write("\tGeneName\tHighConfSet\t");
			bw.write("MaxUniformityScore\tCoverage\tScanStatExp\tK4K36\tK4K36Independan\n");
			//bw.write("\tK4K36\tK4K36Independant\tMaxSpecificity\tMaxExpVal\tMaxOmega\tRepeat\tGwas\tEnhancer\n");
			
			for (String isoName: refLociMap.keySet()){
				
				Locus locus=refLociMap.get(isoName);
				RefSeqGene iso=locus.getReferenceTranscript();
				bw.write(iso.toBED(false));
				
				String geneName=TconGeneMap.get(isoName);
				RefSeqGeneWithIsoforms gene= geneNameGeneMap.get(geneName) ;
				double[] res= locus.getDistanceTo5primeNeighbor(annotatedRef);
				bw.write("\t"+res[0]+"\t"+res[1]);
				res= locus.getDistanceTo3primeNeighbor(annotatedRef);
				bw.write("\t"+res[0]+"\t"+res[1]);
				
				for (int i=0;i<3;i++){
					Integer x=comaptibleIso.get(isoName).get(i)? 1:0 ;
					bw.write("\t"+x);
				}
				
					
				HashMap<String,Integer> setNumIsoMap=locus.getNumIsoformsPerSet();
				for (String s:setsLst ){
					if (setNumIsoMap.containsKey(s))
						bw.write("\t"+setNumIsoMap.get(s));
					else
						bw.write("\t0");
				}
				Integer x=comaptibleGene.get(geneName)? 1:0 ;
				bw.write("\t"+gene.getName()+"\t"+x);
				
				double[] tmp= new double[5];
				//Uniformity, scan k4
				if(uniScr != null  && uniScr.containsRow(isoName) && cov != null  && cov.containsRow(isoName)){
					double[]a1=uniScr.getRow(isoName);
					double[]c1=cov.getRow(isoName);
					for (int j=0;j<a1.length;j++){
						if (a1[j]>tmp[0]){
							tmp[0]=a1[j]; 
							tmp[1]=c1[j];
						}
					}
				}
						
				if (scanStat != null && scanStat.containsKey(isoName)) 	tmp[2]=1;
				if (k4k36 != null && k4k36.containsKey(isoName)) 	tmp[3]=1;
				if (ik4k36 != null && ik4k36.containsKey(isoName)) 	tmp[4]=1;
				
				bw.write("\t"+tmp[0]+"\t"+tmp[1]+"\t"+tmp[2]+"\t"+tmp[3]+"\t"+tmp[4]+"\n");
					
			}
		
			
			
			
	}

	
	

	private static void AddToUniqRefIsoformStatistics(String inTab,
			String updateCol, String title, String bedF, BufferedWriter bw) throws IOException {
		
		HashMap<String,Integer> markedGenes = new HashMap<String,Integer>();
		BEDFileParser bed= new BEDFileParser(bedF);
		for (RefSeqGene g: bed.GetGenes())
			markedGenes.put(g.getName(),1 );
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(inTab)));
		String nextLine;
		int col=-1;
		String[] Header=reader.readLine().split("\t");
		ArrayList <String> newHeader = new ArrayList <String>();
		for(int i=0; i<Header.length; i++){
			newHeader.add(Header[i]);
			if (updateCol.equalsIgnoreCase(Header[i]))
				col=i;
		}
		
		if (col==-1)
			newHeader.add(title);
		else
			newHeader.set(col,title);
		writeTabLine( newHeader,bw);
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			String[] a=nextLine.split("\t");
			int val=0;
			if (markedGenes.containsKey(a[3])){
				//System.err.println("Found: "+a[3]);
				val=1;
			}
			
			if(col==-1)
				bw.write(nextLine+"\t"+val+"\n");
			else{
				a[col]=String.valueOf(val);
				writeTabLine(a,bw);
			}
			
		}
	}

	
 // by using the inBED option one can provide an extended bed (inBED=true) or a GCT (inBED false)
	private static void AddToUniqRefIsoformStatistics(String isoGeneMapF,
			boolean byMaxCriteria, String columnCriteriaString, String inTab,
			 String title, String bedF, BufferedWriter bw, String mode, boolean inBED, boolean sumCriteria) throws IOException, ParseException {

		System.err.println("Mode is:  " +mode);
		String nextLine;
		//Read isoform bed file
		BEDFileParser trait_bed=null;
		MatrixWithHeaders trait_mat=null;
		//TODO- read the first line of the file , find the criteria column and write the header
		String[] header2=null;
		int shiftVal=0;
		
		if (inBED){
			trait_bed= new BEDFileParser(bedF);
			BufferedReader r=new BufferedReader(new InputStreamReader(new FileInputStream(bedF)));
			header2 =r.readLine().split("\t");
			r.close();
			shiftVal=12;
		}
		else{
			trait_mat = new MatrixWithHeaders(bedF);
			List<String> d_arr =  trait_mat.getColumnNames();
			header2 = new String[d_arr.size()];
			for (int i=0; i<d_arr.size(); i++) header2[i]=d_arr.get(i);
		}
		
		
		
		
		//Read mapping of isoforms to genes
		HashMap<String, LinkedList<String>> isoGeneMap = loadIsoGeneMap(isoGeneMapF);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(inTab)));
		
		
		int col=-1;
		int updateColumn=4;
		String[] Header=reader.readLine().split("\t");
		ArrayList <String> newHeader = new ArrayList <String>();
		for(int i=0; i<Header.length; i++){
			newHeader.add(Header[i]);
			if (Header[i].equalsIgnoreCase(columnCriteriaString))
				updateColumn=i;
		}
		System.err.println("Update col is :  " + updateColumn);
		int columnCriteria =0;
		String currCriteria = columnCriteriaString;
		
		String[] emptyExtraFields = new String[header2.length-shiftVal];
		for(int i=shiftVal; i<header2.length; i++){
			if (mode.equalsIgnoreCase("extendMode")){ //add headers only in extend mode/not update
				newHeader.add(header2[i]);
				emptyExtraFields[i-shiftVal]="0";
			}
			if (currCriteria.equalsIgnoreCase(header2[i]))
					columnCriteria=i-shiftVal;
		}
		
		
		System.err.println("Column criteria is : " + columnCriteria);
		writeTabLine( newHeader,bw);
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			String[] a=nextLine.split("\t");
			String[] addFields = emptyExtraFields;
			if (inBED){
				RefSeqGene refIso= selectRepresentativeIsoform(isoGeneMap,a[3],trait_bed.getNameGeneMap(),byMaxCriteria,columnCriteria);	
				if (mode.equalsIgnoreCase("extendMode")){
					
					if (refIso !=null && refIso.getExtraFields() != null)
						addFields=refIso.getExtraFields();
					bw.write(nextLine+"\t");
					writeTabLine(addFields,bw);
				}
				else { //update mode
					a[updateColumn]=refIso.getExtraFields(columnCriteria);
					//if (a[3].equalsIgnoreCase("XLOC_004399")){
					//	System.err.println("updatedColumn is : "+updateColumn + "value is " +  a[updateColumn]);
					//}
					writeTabLine(a,bw);
					}
				}
			
			else { // input was in  table and not bed format
				
				String iso_name = selectRepresentativeIsoform(isoGeneMap,a[3],trait_mat,byMaxCriteria,columnCriteriaString);
				if (mode.equalsIgnoreCase("extendMode")){
					if (! iso_name.isEmpty()){
						 double [] numArr= trait_mat.getRow(iso_name);
						 for (int i=0;i<numArr.length;i++) addFields[i]=String.valueOf(numArr[i]);
					}
					bw.write(nextLine+"\t");
					writeTabLine(addFields,bw);
					
				}
				else{
					if (! iso_name.isEmpty() & trait_mat.containsColumn(columnCriteriaString)){
						a[updateColumn]=String.valueOf(trait_mat.get(iso_name,columnCriteriaString));
					}
					writeTabLine(a,bw);
				}
			}
			
			
			
		}
	}

	
	
	


	private static void GenerateIsoformsTissueMat(String refSetGTF,String setLstFile,Boolean isTab ,boolean MultiExonMode,String outf,boolean fullCompatible) throws IOException, ParseException{
		
		GTFFileParser uniqTconRef= new GTFFileParser(refSetGTF);
		ArrayList<String[]> setLstPaths=getNameLst(setLstFile);
		HashMap<String, List<String>> TconSetsMap=new HashMap<String, List<String>>();
		
		List <String> setNames = new LinkedList <String>();
		List <String> refIsoNames = new LinkedList <String>();
		List <String> refGeneNames = new LinkedList <String>();
		HashMap<String,String> isoGenMap = new  HashMap<String,String>();
		//Go over every set and mark all the isoforms that are compatible with it
		
		
		for (int i=0; i < setLstPaths.size();i++){
			String setName =setLstPaths.get(i)[1];
			String setPath  =setLstPaths.get(i)[0];
			 setNames.add(setName);
			BEDFileParser bed= new BEDFileParser(setPath);
			if (MultiExonMode)
				bed=bed.removeSingeleExon();
			for (RefSeqGene refiso: uniqTconRef.GetGenes()){ //go over the lincs ref file 
				if (i==0){
					refIsoNames.add(refiso.getName());
					String g =refiso.getAttribute("gene_id");
					refGeneNames.add(g);
					isoGenMap.put(refiso.getName(),g);
				}
				Iterator<RefSeqGeneWithIsoforms> tconIt=bed.getOverlappers(refiso).valueIterator();
					while(tconIt.hasNext()){
						RefSeqGeneWithIsoforms MainGene =tconIt.next();
						for(RefSeqGene g_iso:MainGene.getAllIsoforms()){
							boolean hit=false;
								if (fullCompatible)
									hit= Locus.isFullyComaptible(g_iso,refiso,false) ;
								else
									hit = refiso.overlaps(g_iso);

								if (hit){
									if (! TconSetsMap.containsKey(refiso.getName()))
										TconSetsMap.put(refiso.getName(),new LinkedList <String>() );
									TconSetsMap.get(refiso.getName()).add(setName);
								}
						}	
					}
				}
		}
								
		
		//Make a matrix
		MatrixWithHeaders mat = new MatrixWithHeaders (refIsoNames,setNames,refGeneNames);
		
		for (String isoname: TconSetsMap.keySet()){
			double[] vals = new double[setNames.size()];
			for (int i=0; i<vals.length;i++) {vals[i]=0;} 
			//mat.addRow(isoname,isoGenMap.get(isoname), vals);
			for(String set: TconSetsMap.get(isoname)){
				mat.set(isoname,set,1);
			}
		}
		
		//write
		mat.writeGCT(outf);
		
		
		
		
	}
	
	
	//selects the isoform with max value of a give trait given an input matrix
	//matric description equals the gene name, and name equals the isoform name
	private static String selectRepresentativeIsoform(
			HashMap<String, LinkedList<String>> isoGeneMap, String geneName,
			MatrixWithHeaders traitMat, boolean byMaxCriteria,
			String columnCriteriaString) {

		String res_iso_name="";
		double resVal =0;
		LinkedList<String> isoList = isoGeneMap.get(geneName);
		boolean first=true;
		if (! traitMat.containsColumn(columnCriteriaString))
			return res_iso_name;
		for (String currIso : isoList){
			if (first){
				resVal = traitMat.get(currIso,columnCriteriaString);
				res_iso_name= currIso;
				first=false;
			}
			else {
				double v = traitMat.get(currIso,columnCriteriaString);
				if ((byMaxCriteria & (v > resVal)) | (byMaxCriteria==false & (v < resVal))){
					resVal=v;
					res_iso_name=currIso;
				}
			}
		}
		
		return res_iso_name;
	}	
		


	
	private static RefSeqGene selectRepresentativeIsoform(
			HashMap<String, LinkedList<String>> isoGeneMap, String geneName,
			HashMap<String,RefSeqGene> nameIsoMap, boolean byMaxCriteria, int columnCriteria) {

		LinkedList<String> isoNameLst = isoGeneMap.get(geneName);
		if (isoNameLst==null)
			System.err.println(geneName);
		boolean first = true;
		double currVal=0;
		RefSeqGene selectedIso = null;
		for (String name: isoNameLst){
			if (geneName.equalsIgnoreCase("XLOC_004399")){
				System.err.println(name);
			}
			RefSeqGene iso=nameIsoMap.get(name);
			if (iso==null || iso.getExtraFields(columnCriteria).equalsIgnoreCase("")) //in case the isoform was not annotated with the updated value
				continue;
			if (first){
				selectedIso=iso;
				currVal= Double.valueOf(selectedIso.getExtraFields(columnCriteria));
				first=false;
				if (geneName.equalsIgnoreCase("XLOC_004399")){
					System.err.println(name + " curr val : " + currVal);
				}
			}
			else{
				
				if (geneName.equalsIgnoreCase("XLOC_004399")){
					System.err.println("Not First" + name + " curr val : " + currVal);
				}
				
				if (byMaxCriteria & currVal < Double.valueOf(iso.getExtraFields(columnCriteria))){
					selectedIso=iso;
					currVal= Double.valueOf(selectedIso.getExtraFields(columnCriteria));
					if (geneName.equalsIgnoreCase("XLOC_004399")){
						System.err.println(name + " curr val : " + currVal);
					}
				}
				if (byMaxCriteria==false & currVal > Double.valueOf(iso.getExtraFields(columnCriteria))){
					selectedIso=iso;
					currVal= Double.valueOf(selectedIso.getExtraFields(columnCriteria));
				}
				
			}
			
		}
		return selectedIso;
	}





	private static HashMap<String, LinkedList<String>> loadIsoGeneMap(
			String isoGeneMapF) throws IOException {
		
		String nextLine;
		HashMap<String, LinkedList<String>> isoGeneMap= new HashMap<String, LinkedList<String>>();
		BufferedReader mapreader=new BufferedReader(new InputStreamReader(new FileInputStream(isoGeneMapF)));
		
		while ((nextLine = mapreader.readLine()) != null && (nextLine.trim().length() > 0)) {
			String[] a=nextLine.split("\t");
			if (! isoGeneMap.containsKey(a[1]))
				isoGeneMap.put(a[1],new LinkedList<String>());
			isoGeneMap.get(a[1]).add(a[0]);
			
		}
		return isoGeneMap;
	}





	private static void writeTabLine(String[] a, BufferedWriter bw) throws IOException {
		int i=0;
		for( i=0; i<a.length-1; i++)
			bw.write(a[i]+"\t");
		bw.write(a[i]+"\n");
	}
	
	private static void writeTabLine(ArrayList<String> a, BufferedWriter bw) throws IOException {
		int i=0;
		for( i=0; i<a.size()-1; i++)
			bw.write(a.get(i)+"\t");
		bw.write(a.get(i)+"\n");
	}




	private static HashMap<String, Double> getNamesFromBed(String file) throws IOException {
		
		HashMap<String, Double> res= new HashMap<String, Double>();
		BEDFileParser bed= new BEDFileParser(file);
		for (RefSeqGene g: bed.GetGenes()){
			res.put(g.getName(), new Double (g.getBedScore()));
		}

		return res;
	}





	private static HashMap<String, Boolean> markComaptibleGene(
					BEDFileParser uniqGeneRef, HashMap<String, ArrayList<Boolean>> comaptibleIso,
					HashMap<String, List<String>> geneTconMap) {
		
		HashMap<String, Boolean> res=new HashMap<String, Boolean>();
		
		for(String geneName : geneTconMap.keySet()){
			Boolean b=false;
			for (String isoName: geneTconMap.get(geneName)){
				if (comaptibleIso.get(isoName).get(0)==true || comaptibleIso.get(isoName).get(1)==true || comaptibleIso.get(isoName).get(2)==true){
					b=true; 
					break;
				}
			}
			res.put(geneName,b);
		}
		
		return res;
	}




//Returns: 3 overlaps : 1. by two assemblers, 2. by assembler and external reference , 3. by two samples(no ref) 
	private static HashMap<String, ArrayList<Boolean>> markComaptibleIso(
					Map<String, Locus> refLociMap,
					ArrayList<String> BwnSetLst1, ArrayList<String> BwnSetLst2,ArrayList<String> BwnSetLst3,
					ArrayList<String> WithinSetsLst) {
			
		HashMap<String, ArrayList<Boolean>> res= new HashMap<String, ArrayList<Boolean>>();
		
		for (String g: refLociMap.keySet()){
			ArrayList<Boolean> resArr=new ArrayList<Boolean>();
			
			Locus l = refLociMap.get(g);
			int a=l.numOfOverlapSets(BwnSetLst1);
			int b=l.numOfOverlapSets(BwnSetLst2);
			int ref=l.numOfOverlapSets(BwnSetLst3);
			int c=l.numOfOverlapSets(WithinSetsLst);
			//Scripture and Cufflinks
			if (a>0 & b>0)
				resArr.add(true);
			else
				resArr.add(false);
			// assembler and external reference
			if ((a>0 || b> 0) & ref >0)
				resArr.add(true);
			else
				resArr.add(false);
			// by two samples (not ref)
			if (c>1) 
				resArr.add(true);
			else
				resArr.add(false);
			res.put(g, resArr);
		}

		return res;
	}



	//uniqTconRef has to be loaded from a gtf
	private static BEDFileParser AssociateGeneTranscript(
					BEDFileParser uniqTconRef,
					HashMap<String, List<String>> geneTconMap,
					Map<String, String> tconGeneMap, Map<String, RefSeqGeneWithIsoforms> geneNameGeneMap) {
			
			BEDFileParser uniqGenes = new BEDFileParser();
			HashMap <String,Locus> geneMap = new HashMap <String,Locus>();
			
			for (RefSeqGene iso: uniqTconRef.GetGenes()){
				String gene_id=GFF.getGTFAttr(iso.getExtraFields()[0], "gene_id");
				gene_id=gene_id.replace("\"","");
				if (!geneMap.containsKey(gene_id))
					geneMap.put(gene_id,new Locus(new RefSeqGeneWithIsoforms(iso)) );
				geneMap.get(gene_id).addIsoform(new RefSeqGeneWithIsoforms(iso), "def");
			}
			System.err.println("Error.Locus was merged to more than 1 transcript");
			
			for (String geneName: geneMap.keySet()){
				RefSeqGeneWithIsoforms gene= new RefSeqGeneWithIsoforms( geneMap.get(geneName).getMerged());
				gene.setName(geneName);
				uniqGenes.addRefSeq(gene);
				LinkedList<String> lst=new LinkedList <String> ();
				geneTconMap.put(gene.getName(),lst);
				geneNameGeneMap.put(gene.getName(),gene);
				
				Collection<RefSeqGene> allIso=(geneMap.get(geneName).getAllIsoforms(""));
				for (RefSeqGene i:allIso){
					tconGeneMap.put(i.getName(),gene.getName());
					geneTconMap.get(gene.getName()).add(i.getName());
				}
				
			}

			return uniqGenes;
	}


	private static LinkedList<String> loadCompatibleToRefLociIsoformMap(BEDFileParser uniqTconRef,
			String setLstFile,Map<String, Locus> refLociMap,boolean multiExonMode,boolean fullCompatible) throws IOException {

		Map<String , RefSeqGene> nameLociMap=new HashMap<String , RefSeqGene>();
		for(RefSeqGene g:uniqTconRef.GetGenes()){
			refLociMap.put(g.getName(),new Locus(new RefSeqGeneWithIsoforms(g)));
			nameLociMap.put(g.getName(),g);
		}
			
		String line;
		//Go over the different sets and load to the COMPATIBLE Locus
		BufferedReader br=new BufferedReader(new FileReader(setLstFile));
		Map<String,String> setsPath= new HashMap<String,String>();
		 while((line = br.readLine()) != null) {
				line = line.trim();
				String [] lineSplit = line.split("\t");
				setsPath.put(lineSplit[1], lineSplit[0]);	
		 }
		 return ( loadCompatibleToRefLociIsoformMap(uniqTconRef, setsPath, refLociMap, multiExonMode, fullCompatible));
	}
		 
	private static LinkedList<String> loadCompatibleToRefLociIsoformMap(BEDFileParser uniqTconRef,
			Map<String,String> setsPath,Map<String, Locus> refLociMap,boolean multiExonMode,boolean fullCompatible) throws IOException {
 
		LinkedList<String> setsLst=new LinkedList(setsPath.keySet());
		Collections.sort(setsLst);
		for (String setName: setsLst){
			BEDFileParser bed= new BEDFileParser(setsPath.get(setName));
			if (multiExonMode)
				bed=bed.removeSingeleExon();
			for (RefSeqGene iso: bed.GetGenes()){
				Iterator<RefSeqGeneWithIsoforms> lociIt=uniqTconRef.getOverlappers(iso).valueIterator();
					while(lociIt.hasNext()){
						RefSeqGeneWithIsoforms MainGene =lociIt.next();
						for(RefSeqGene g_iso:MainGene.getAllIsoforms()){
							if (refLociMap.containsKey(g_iso.getName())){
								//RefSeqGene lookupInstance= nameLociMap.get(g_iso.getName());
								Locus g_locus = refLociMap.get(g_iso.getName()) ;
								if (fullCompatible){ // Case where you want all exons of one
									//of the compared transcripts of one to overlap with the other
									if (Locus.isFullyComaptible(g_locus.getReferenceTranscript(),iso,false))
										g_locus.addIsoform(new RefSeqGeneWithIsoforms(iso), setName);
								}
								else{ //ask for one exon overlap
									if (  g_locus.getReferenceTranscript().overlaps(iso))
										g_locus.addIsoform(new RefSeqGeneWithIsoforms(iso), setName);
				
								}
							}
						}
					}
			}
		}
		return setsLst;
	}

		
	



	//This function takes a refset of merged loci,and finds in which sets it was detected.
	//It then compares if the loci had an isoform that s compatible between 2 sets lists (set1 or set2) or between any 2 samples of sets (betweenAllCompatibleSet)
		private static void GenerateIsoformsStatistics(String refSet,String setLstFile,String annotatedRefFile, BufferedWriter bw,String set1LstF,String set2LstF,Boolean isTab, String betweenAllCompatibleSet,String uniqIsoF, Boolean MultiExonMode) throws IOException{

			
			BEDFileParser annotatedRef= new BEDFileParser(annotatedRefFile);
			ArrayList<String> set1Lst=loadNameLst(set1LstF,isTab);
			ArrayList<String> set2Lst=loadNameLst(set2LstF,isTab);
			ArrayList<String> set3Lst=loadNameLst(betweenAllCompatibleSet,isTab);
			Map<RefSeqGene,Locus> refLociMap=new  HashMap<RefSeqGene,Locus>();
			
			LinkedList<String> setsLst =loadRefLociIsoformMap (refSet,setLstFile,refLociMap,MultiExonMode) ;
			
			
			//If a cuffcompare reference is available - compare to that reference
			BEDFileParser uniqIsoBed=null;
			if (! uniqIsoF.equals(""))
			 uniqIsoBed=new BEDFileParser(uniqIsoF);
			
			/*
			//Load the reference loci of lincRNAs
			BEDFileParser ref= new BEDFileParser(refSet);
			Iterator<String> chrIt = ref.getChromosomeIterator();
			while (chrIt.hasNext()){
				Iterator<RefSeqGeneWithIsoforms> geneIt= ref.getChrTree(chrIt.next()).valueIterator();
				while(geneIt.hasNext()){
					RefSeqGeneWithIsoforms g=geneIt.next();
					refLociMap.put(g,new Locus(g));
				}
			}
			String line;
			//Go over the different sets and load to the overlapping Locus
			BufferedReader br=new BufferedReader(new FileReader(setLst));
			Map<String,String> setsPath= new HashMap<String,String>();
			 while((line = br.readLine()) != null) {
					line = line.trim();
					String [] lineSplit = line.split("\t");
					setsPath.put(lineSplit[1], lineSplit[0]);	
			 }		
			LinkedList<String> setsLst=new LinkedList(setsPath.keySet());
			Collections.sort(setsLst);
			for (String setName: setsLst){
				BEDFileParser bed= new BEDFileParser(setsPath.get(setName));
				chrIt=bed.getChromosomeIterator();
				while (chrIt.hasNext()){
					Iterator<RefSeqGeneWithIsoforms> geneIt= bed.getChrTree(chrIt.next()).valueIterator();
					while (geneIt.hasNext()){
						RefSeqGeneWithIsoforms iso =geneIt.next();
						Iterator<RefSeqGeneWithIsoforms> lociIt=ref.getOverlappers(iso).valueIterator();
						while(lociIt.hasNext()){
							RefSeqGeneWithIsoforms locus =lociIt.next();
							if (locus.getOrientation().equalsIgnoreCase(iso.getOrientation()))
								refLociMap.get(locus).addIsoform(iso,setName);
						}
					}
				}
			}
			*/
			
			bw.write("chr\tstart\tend\tname\tscore\tstrand\ttStart\ttEnd\trgb\tblockCnt\tblockSizes\tblockStarts");
			bw.write("\tis5'EqOrient\tt5'NeighborDist\tis3'EqOrient\t3'NeighborDist\tCompatibleBetweenAssemblies\tCompatibleBetweenSets");
			//number of isoforms in set
			for (String s:setsLst)
				bw.write("\t"+s);
			//
			bw.write("\tmaxScore\n");
			for (RefSeqGene g: refLociMap.keySet()){
				Locus locus=refLociMap.get(g);
				bw.write(g.toBED());
				double[] res= locus.getDistanceTo5primeNeighbor(annotatedRef);
				bw.write("\t"+res[0]+"\t"+res[1]);
				res= locus.getDistanceTo3primeNeighbor(annotatedRef);
				bw.write("\t"+res[0]+"\t"+res[1]);
				
				if (uniqIsoBed==null){
				bw.write("\t"+locus.getIsComaptibleBetweenSubsets(set1Lst, set2Lst, true));
				bw.write("\t"+locus.getIsComaptibleBetweenAny2Sets(set3Lst, true));
				}
				else{
					int x= locus.getIsComaptibleBetweenSubsets(set1Lst, set2Lst, uniqIsoBed);
					bw.write("\t"+x);
					 x= locus.getIsComaptibleBetweenAny2Sets(set3Lst, uniqIsoBed);
					bw.write("\t"+x);	
				}	
				HashMap<String,Integer> setNumIsoMap=locus.getNumIsoformsPerSet();
				HashMap<String, Double> setMaxScoreMap=locus.getMaxScorePerSet();
				double maxScr=0.0;
				for (String s:setsLst ){
					if (setNumIsoMap.containsKey(s)){
						bw.write("\t"+setNumIsoMap.get(s));
						maxScr=Math.max(maxScr, setMaxScoreMap.get(s));
					}
					else
						bw.write("\t0");
				}
				bw.write("\t"+maxScr+"\n");
			}
		
	}




		//Given a list of Loci and a list of samples, the function will extract for every loci:
		//1. all cuffCmprIso = uniq set of isoform structures
		//2. all the isofroms in the loci from all samples (not including the cuffCmpr)
		//3. all compatible between subsets
		//4. all compatible between any 2 sets
		//5. Select 1 measuring isoform: by the next order- longest Intron Chain isoform from line 3, line 4 and finally line 1
		private static void ExtractIsoforms(String refSet, String setLstFile,String set1LstF, String set2LstF,
				Boolean isTab, String betweenAllCompatibleSet,	String cuffcmprSet, String prefix) throws IOException {
		
			ArrayList<String> set1Lst=loadNameLst(set1LstF,isTab);
			ArrayList<String> set2Lst=loadNameLst(set2LstF,isTab);
			ArrayList<String> set3Lst=loadNameLst(betweenAllCompatibleSet,isTab);
			Map<RefSeqGene,Locus> refLociMap=new  HashMap<RefSeqGene,Locus>();
			LinkedList<String> setsLst =loadRefLociIsoformMap (refSet,setLstFile,refLociMap,false) ;
			
			BEDFileParser unq=new BEDFileParser();
			BEDFileParser all=new BEDFileParser();
			BEDFileParser compBetweenSubsets=new BEDFileParser();
			BEDFileParser compBetweenAny2=new BEDFileParser();
			BEDFileParser probeIso=new BEDFileParser();
			
			for (RefSeqGene g: refLociMap.keySet()){
				Locus locus=refLociMap.get(g);
				//report unique isoforms
				locus.updateIsoformsNameWithLocusName();
				Collection<RefSeqGene> unqIsoSet=locus.getSetIsoforms(cuffcmprSet);
				unq.addRefSeqSet(unqIsoSet);
				//report all isoforms that overlap the loci
				all.addRefSeqSet(locus.getAllIsoforms(cuffcmprSet));
				//report isoforms that are compatible between at least 2 samples
				Collection<RefSeqGene> between2IsoSet=locus.getUnqComaptibleBetweenAny2Sets(set3Lst,unqIsoSet,true);
				compBetweenAny2.addRefSeqSet(between2IsoSet);
				Collection< RefSeqGene> betweenSubsetsIsoSet =locus.getUnqComaptibleBetweenSubsets(set1Lst,set2Lst,unqIsoSet,true);
				compBetweenSubsets.addRefSeqSet(betweenSubsetsIsoSet);
				RefSeqGene prb=selectBestProbePerLoci(betweenSubsetsIsoSet,between2IsoSet,unqIsoSet);
				if (prb != null)
					probeIso.addRefSeq(prb);
				
			}
			
			unq.writeFullBed(prefix+"UniqIso.bed");
			all.writeFullBed(prefix+"AllIso.bed");
			compBetweenSubsets.writeFullBed(prefix+"CompatibleBwnSubsetIso.bed");
			compBetweenAny2.writeFullBed(prefix+"CompatibleBwnAny2Iso.bed");
			probeIso.writeFullBed(prefix+"MeasuresIso.bed");
	}

		
	private static RefSeqGene selectBestProbePerLoci(Collection<RefSeqGene> betweenSubsetsIsoSet,
				Collection<RefSeqGene> between2IsoSet,Collection<RefSeqGene> unqIsoSet) {
		
		Collection<RefSeqGene> res=new LinkedList<RefSeqGene>();
		if (!betweenSubsetsIsoSet.isEmpty())
			res=Locus.SelectLongestIntronChainCandidate(betweenSubsetsIsoSet);
		else if (!between2IsoSet.isEmpty())
			res=Locus.SelectLongestIntronChainCandidate(between2IsoSet);
		else
			res=Locus.SelectLongestIntronChainCandidate(unqIsoSet);
		return selectMaxScrIso(res);
	}




//At The moment this is arbitrary since the score value is not set for the reference isoform set: i.e all =0
	private static RefSeqGene selectMaxScrIso(Collection<RefSeqGene> set) {
			double scr=-1;
			for (RefSeqGene g: set )
				scr=Math.max(scr,g.getBedScore());
			for (RefSeqGene g: set ){
				if (g.getBedScore()==scr)
					return g;
			}
			return null;
	}

//
	private static void ScoreByMaxScoredCompatibleIsoform(String refSet,String setLstFile, String savefile, boolean isFully) throws IOException {
		
		Map<RefSeqGene,Locus> refLociMap=new  HashMap<RefSeqGene,Locus>();
		LinkedList<String> setsLst =loadRefLociIsoformMap (refSet,setLstFile,refLociMap,false) ;
		BEDFileParser resBed=new BEDFileParser();
		double maxScr=0.0;
		
		for (RefSeqGene g: refLociMap.keySet()){
			Locus locus=refLociMap.get(g);
			Collection <RefSeqGene> set=locus.getComaptibleIsoforms(g,isFully);
			RefSeqGene maxIso= selectMaxScrIso(set)	;
			if (maxIso!= null)
				maxScr=maxIso.getBedScore();
			g.setBedScore(maxScr);
			resBed.addRefSeq(g);
		}
		
		resBed.writeFullBed(savefile);
	}


	private static void SelectLongestIntornChainIsoPerLoci(String refBED,String setLstFile, String savefile,String refSet) throws IOException {
		Map<RefSeqGene,Locus> refLociMap=new  HashMap<RefSeqGene,Locus>();
		LinkedList<String> setsLst =loadRefLociIsoformMap (refBED,setLstFile,refLociMap,false) ;
		BEDFileParser resBed=new BEDFileParser();
		System.err.println("Ref set is: "+refSet);
		for (RefSeqGene g: refLociMap.keySet()){
			Locus locus=refLociMap.get(g);
			RefSeqGene iso=locus.getLongestIntornChainIso( refSet);
			if (iso!= null)
				resBed.addRefSeq(iso);
			else
				System.err.println("No isoforms found for " + g.toBED());
		}
		resBed.writeFullBed(savefile);
	}

		private static LinkedList<String>  loadRefLociIsoformMap (String refSet,String setLstFile,Map<RefSeqGene,Locus> refLociMap, Boolean multiExonMode) throws IOException{
			
			BEDFileParser ref= new BEDFileParser(refSet);
			Iterator<String> chrIt = ref.getChromosomeIterator();
			while (chrIt.hasNext()){
				Iterator<RefSeqGeneWithIsoforms> geneIt= ref.getChrTree(chrIt.next()).valueIterator();
				while(geneIt.hasNext()){
					RefSeqGeneWithIsoforms g=geneIt.next();
					refLociMap.put(g,new Locus(g));
				}
			}
			String line;
			//Go over the different sets and load to the overlapping Locus
			BufferedReader br=new BufferedReader(new FileReader(setLstFile));
			Map<String,String> setsPath= new HashMap<String,String>();
			 while((line = br.readLine()) != null) {
					line = line.trim();
					String [] lineSplit = line.split("\t");
					setsPath.put(lineSplit[1], lineSplit[0]);	
			 }		
			LinkedList<String> setsLst=new LinkedList(setsPath.keySet());
			Collections.sort(setsLst);
			for (String setName: setsLst){
				BEDFileParser bed= new BEDFileParser(setsPath.get(setName));
				if (multiExonMode)
					bed=bed.removeSingeleExon();
				
				chrIt=bed.getChromosomeIterator();
				while (chrIt.hasNext()){
					Iterator<RefSeqGeneWithIsoforms> geneIt= bed.getChrTree(chrIt.next()).valueIterator();
					while (geneIt.hasNext()){
						RefSeqGeneWithIsoforms iso =geneIt.next();
						Iterator<RefSeqGeneWithIsoforms> lociIt=ref.getOverlappers(iso).valueIterator();
						while(lociIt.hasNext()){
							RefSeqGeneWithIsoforms locus =lociIt.next();
							if (locus.getOrientation().equalsIgnoreCase(iso.getOrientation()))
								refLociMap.get(locus).addIsoform(iso,setName);
						}
					}
				}
			}
			
			return setsLst;
			
		}

		
		
		




		private static void writeAnnotationReader(
			AnnotationReader<? extends GenomicAnnotation> set,
			BufferedWriter bw) throws IOException {
		   
			Iterator<String> chrIt = set.getChromosomeAnnotationMap().keySet().iterator();
			while(chrIt.hasNext()) {
				String chr = chrIt.next();
				IntervalTree<? extends GenomicAnnotation> tree = set.getChromosomeTree(chr);
				Iterator<? extends GenomicAnnotation> annotIt = tree.valueIterator();
				while(annotIt.hasNext()) {
					LightweightGenomicAnnotation annot = annotIt.next();
					bw.write(annot.toString());
					bw.newLine();
				}
			}
		
	}



		/**
		 * ADDED BY MORAN Feb 19th 2010
		 * Filters this reader by removing all elements that do not (overlap AND HAVE THE SAME ORIENTATION) as the given list of annotations 
		 * @param annotations a list of Annotations to filter
		 * @throws IOException 
		 * 
		 */
		@SuppressWarnings("unchecked")
		public static void writeIntersectionToExtendedBed(AnnotationReader<? extends GenomicAnnotation> arrayAnnot , List<BED> annotations, BufferedWriter bw ) throws IOException {
			
			Iterator<BED> it = annotations.iterator(); 
			
			//loop over annotation list
			while(it.hasNext()) {
				BED annot = it.next();
				String lst="\t";
				 
				IntervalTree<? extends GenomicAnnotation> tree = arrayAnnot.getChromosomeTree(annot.getChromosome());
				if(tree == null) {
					continue;
				}
				//Find all array probes that overlap with the annotation 
				Iterator<?> overlaperIt = tree.overlappers(annot.getStart(), annot.getEnd());
				while(overlaperIt.hasNext()) {
					Node<GenomicAnnotation> overlapperNode = (Node<GenomicAnnotation>) overlaperIt.next();
					GenomicAnnotation overlapper = overlapperNode.getValue();
					//if an exon (block) of this annotation intersects a probe/region on the array
					//write to extended bed file the name of this probe
					if (annot.IntersectBlocks(overlapper) )
					{
						if ( overlapper.getOrientation().equalsIgnoreCase("*") ||  overlapper.getOrientation().equalsIgnoreCase(annot.getOrientation()))
						{
							String n=overlapper.getName();
							
							lst=lst.concat(n);
							lst=lst.concat(",");
							
						}
					}
					
				}
				
				//write to the file
				bw.write(annot.toString());
				bw.write(lst);
				bw.newLine();
			}
					
		
		}
		
		
		
		public static HashMap<String,AnnotationReader<? extends GenomicAnnotation>> uploadFiles( String setLst, String setFormat, ArrayList<String> nameLst) throws IOException, ParseException {
			
			HashMap <String,AnnotationReader< ? extends GenomicAnnotation>> annotSet= new HashMap<String,AnnotationReader<? extends GenomicAnnotation>>(); 
			BufferedReader br = new BufferedReader(new FileReader(setLst));
			String line;
			
			 while((line = br.readLine()) != null) {
					line = line.trim();
					String [] lineSplit = line.split("\t");
					AnnotationReader<? extends GenomicAnnotation> set = AnnotationReaderFactory.create(lineSplit[0], setFormat);
					annotSet.put(lineSplit[1], set);
					nameLst.add(lineSplit[1]);
			 }		
			
			return annotSet;	
		}
		
   
		private static void calcOverlapOfMultipleSets(String setLst,String setFormat, BufferedWriter bw) throws IOException, ParseException {
			
			bw.write("SetName\t"+"SetOverlapWithSuperSet"+"\t"+"SuperSetOverlapWithSet"+"\t"+"SetSize"+"\t"+"SuperSetSize"+"\t"+"AdditionOfNewSet"+"\n");
			ArrayList <String> nameLst= new ArrayList <String> ();
			//read to set of sets
			Map <String,AnnotationReader<? extends GenomicAnnotation>> setOfSets= uploadFiles(setLst,setFormat,nameLst);
			//merged set= Union set
			AnnotationReader<? extends GenomicAnnotation> superSet = AnnotationReaderFactory.create(setFormat);
			//Iterator<String> setsIt = setOfSets.keySet().iterator();
			Iterator<String> setsIt = nameLst.iterator();
			int newSetSize=0;
			while(setsIt.hasNext()) {
				String currName=setsIt.next();
				AnnotationReader<? extends GenomicAnnotation> currSet= setOfSets.get(currName);
				//increase subset by one set on calc overlap with Union set, write to file
				currSet.merge(); //TODO - MAKE MERGE ORIENTAION SENSATIVE
				int l1= currSet.getOverlappers(superSet.getAnnotationList()).size();
				int l2= superSet.getOverlappers(currSet.getAnnotationList()).size();
				int s1=currSet.getAnnotationList().size();
				int s2= superSet.getAnnotationList().size();
				superSet.concatAnnotationSet(0, currSet.getAnnotationList());
				superSet.merge();
				newSetSize=superSet.getAnnotationList().size();;
				bw.write(currName+"\t"+l1+"\t"+l2+"\t"+s1+"\t"+s2+"\t"+Math.max((newSetSize-s2),0)+"\n");
				
			}
			
			bw.write("final super set size : "+ newSetSize+ "\n");
			
			
		}
		

		
     private static void calcOrientedOverlapOfMultipleSets(String setLst,BufferedWriter bw) throws IOException {
			
			
			bw.write("SetName\t"+"SetOverlapWithSuperSet"+"\t"+"SuperSetOverlapWithSet"+"\t"+"SetSize"+"\t"+"SuperSetSize"+"\t"+"AdditionOfNewSet"+"\n");
			ArrayList <String[]> nameLst= getNameLst(setLst);
			
			Map<String, IntervalTree<RefSeqGene>> superSet= new LinkedHashMap<String, IntervalTree<RefSeqGene>>();
			int newSetSize=0;
			//iterate Over Sets Name
			for (String[] names:nameLst)
			{
			    //read set to a refseq gene set
				Map<String, Collection<RefSeqGene>> set=BEDFileParser.loadDataByChr(new File(names[0]));
				//merge with itself
				Map<String, IntervalTree<RefSeqGene>> mergedSet= GeneTools.merge(set);
				int s1= numOfGenes(mergedSet);
				int s2= numOfGenes(superSet);
				//overlap with superset
				int l1 = numOfGenes(GeneTools.overlap(mergedSet,superSet,false));
				int l2 = numOfGenes(GeneTools.overlap(superSet,mergedSet,false));
				//merge with superset
				superSet=GeneTools.mergeSets(mergedSet,superSet);
				//bw.write("the size of concatenated superset is :"+ numOfGenes(superSet) +"\n");
				superSet=GeneTools.merge( makeCollectionMap(superSet));
				newSetSize=numOfGenes(superSet);
				//bw.write("the size of concatenated and merged superset is :"+ newSetSize +"\n");
				bw.write(names[1]+"\t"+l1+"\t"+l2+"\t"+s1+"\t"+s2+"\t"+Math.max((newSetSize-s2),0)+"\n");
			
			}
			bw.write("final super set size : "+ newSetSize+ "\n");
			
			
		
	}


     private static Map<String, Collection<RefSeqGene>> makeCollectionMap(
			Map<String, IntervalTree<RefSeqGene>> set) {
    	 Map<String, Collection<RefSeqGene>> res= new LinkedHashMap<String, Collection<RefSeqGene>>();
    	 for (String chr: set.keySet()){
    		 res.put(chr,set.get(chr).toCollection());
    	 }
		
		return res;
	}


	private static int numOfGenes( Map<String, IntervalTree<RefSeqGene>> set) {
 		int res=0;
 		for (String s:set.keySet() ) {res+= set.get(s).size();}
 		return res;
 	}


	private static ArrayList<String[]> getNameLst(String setLst) throws IOException {
		
		BufferedReader br = new BufferedReader(new FileReader(setLst));
		String line;
		ArrayList<String[]> nameLst= new ArrayList<String[]>();
		
		while((line = br.readLine()) != null) {
				line = line.trim();
				String [] lineSplit = line.split("\t");
				nameLst.add(lineSplit);
				
		 }
		
		return nameLst;
	}
	
	private static ArrayList<String> loadNameLst(String setLst,Boolean isTab) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(setLst));
		String line;
		ArrayList<String> nameLst= new ArrayList<String>();
		while((line = br.readLine()) != null) {
				line = line.trim();
				if (isTab){//(line.matches("\t")){
					String [] lineSplit = line.split("\t");
					nameLst.add(lineSplit[1]);
				}
				else
					nameLst.add(line);
		 }
		return nameLst;
	}

	
	//support orientation 
	private static void ExtendedBedForBrobeDesign(BEDFileParser refSet , String setLst, String outPrefix,String outDir, BufferedWriter stdout, String trim, int probeSize) throws IOException
	{
		
		//extendedBed, intervalBed
		BufferedWriter bw1 = new BufferedWriter(new FileWriter(outDir+outPrefix + "_extended.bed") );
		BufferedWriter bw2 = new BufferedWriter(new FileWriter(outDir+outPrefix + "_interval.bed"));
		BufferedWriter bwFasta = new BufferedWriter(new FileWriter(outDir+outPrefix + "_probeCandidates.bed"));
		
		BufferedWriter bwKey = new BufferedWriter(new FileWriter(outDir+outPrefix + "_KeyProbeCandidates.bed"));
		
		System.err.println("Opened probe design output file: " + outDir+outPrefix + "_extended.bed");
		
		ArrayList <String> nameLst= new ArrayList <String> ();
		   
		Map <String,BEDFileParser> setOfSets= uploadBedFiles(setLst,nameLst,trim);
		
		
		
		Iterator<RefSeqGene> itConGene=refSet.GetGenes().iterator();
		//For every consensus
		while(itConGene.hasNext()){
			RefSeqGene conGene=itConGene.next();
			int setCntr=0;
			int isoformCntr=0;
			//Generate exon-Isoform-Tree 
			IntervalTree<LightweightGenomicAnnotation> exonIsoformTree =new IntervalTree<LightweightGenomicAnnotation>();
			//For every Set : add overlappers to the isoform tree with a weighted score
			for (String setName: nameLst){
				
				BEDFileParser set=setOfSets.get(setName);
				//set.writeFullBed(stdout);
				IntervalTree<RefSeqGeneWithIsoforms> overlapTree=set.getOverlappers(conGene);
				if (overlapTree.isEmpty())
					continue;
				setCntr++;
				double treesize=BEDFileParser.getNumberOfIsoforms(overlapTree) ;
				double scr= 1.0/treesize;
				isoformCntr+=treesize;
				//Convert the  overlapTree to an exonTree - each exon as a separate gene (with the score of its original gene)
				IntervalTree<LightweightGenomicAnnotation> overlapExonTree = toExonTree(overlapTree,scr);
			
				Iterator<LightweightGenomicAnnotation> itG=overlapExonTree.valueIterator();
				while(itG.hasNext()){
					LightweightGenomicAnnotation g=itG.next();
					
					//(if 2 exons span the same exact region- just update the weighted score)
					Node<LightweightGenomicAnnotation> N=exonIsoformTree.find(g.getStart(),g.getEnd());
					
					if( N== null)
						exonIsoformTree.put(g.getStart(),g.getEnd(), g);
					else{
						N.getValue().setScore(N.getValue().getScore()+g.getScore());
					}
				}
			}
			
			if (exonIsoformTree.isEmpty()) //there are no isoforms overlapping this consensus seq.
				continue;
						
			//Extract all start and end points of exons in this region ,(Sort start and end points)
			ArrayList<Integer> posArr=getIntervalBounderies(exonIsoformTree);
			
		
			IntervalTree<LightweightGenomicAnnotation> IntervalIsoformTree= new IntervalTree<LightweightGenomicAnnotation> ();
			//For each interval, count num of overlapping exons in the exonIsoformTree
			//Insert the interval to the IntervalIsoformTree
			for (int i=0; i<posArr.size()-1;i++){
				int st=posArr.get(i);
				int end=posArr.get(i+1);
				LightweightGenomicAnnotation interval=new BasicLightweightAnnotation(conGene.getChr(),st,end,conGene.getOrientation(),0);
				Iterator<Node<LightweightGenomicAnnotation>> exonsIt= exonIsoformTree.overlappers(st, end);
				while (exonsIt.hasNext()){
					LightweightGenomicAnnotation currExon=exonsIt.next().getValue();
					double scr=currExon.getScore();
					interval.setScore(interval.getScore()+currExon.getScore());
				}
				if (interval.getScore()>0)
					IntervalIsoformTree.put(st, end, interval);
			}
		
			//write extended bed line
			//System.err.println(conGene.toBED());
			//System.err.println("length of position array:" + posArr.size() +" length of interval Tree " + IntervalIsoformTree.size() );
			printExtendedBedLine(conGene,IntervalIsoformTree,bw1,bw2,setCntr,isoformCntr);
			printBedLinePerExon(conGene,IntervalIsoformTree,bwFasta,bwKey,setCntr,probeSize,isoformCntr);
			
		}	
		
		bw1.close();
		bw2.close();
		bwFasta.close();
		bwKey.close();
		System.err.println("Closed probe design output file");
		return;
	}
	
	
	
	//ASSUMES hat the intervals in the intervalIsoformTree are not overlapping
	private static ArrayList<String> printExtendedBedLine(RefSeqGene conGene,
			IntervalTree<LightweightGenomicAnnotation> intervalIsoformTree,
			BufferedWriter bw_extendedBed, BufferedWriter bw_intervalBed, int setCntr,int IsoformCntr) throws IOException {
			
			
			
			//System.err.println(intervalIsoformTree.size()+"\n");
			int[] intervalStarts=new int[intervalIsoformTree.size()];
			int[] intervalEnds=new int[intervalIsoformTree.size()];
			double[] intervalScores=new double[intervalIsoformTree.size()];
			int intervalCount=intervalIsoformTree.size();
			
			int i=0;
			Iterator <LightweightGenomicAnnotation> it=intervalIsoformTree.valueIterator();
			while(it.hasNext()){
				LightweightGenomicAnnotation annot=it.next();
				intervalStarts[i]=annot.getStart();
				intervalEnds[i]=annot.getEnd();
				intervalScores[i]=annot.getScore();
				i++;
			}
			
			String sizes="";
			String starts="";
			String scores="";
			
			int chromStart =intervalStarts[0];
			int chromEnd=intervalEnds[intervalCount-1];
			
			NumberFormat nf = NumberFormat.getInstance();
			nf.setMaximumFractionDigits(3);

			
			for( i=0; i<intervalCount; i++){
				sizes=sizes+(intervalEnds[i]-intervalStarts[i])+",";
				starts=starts+(intervalStarts[i]-chromStart)+",";
				//scores=scores+(intervalScores[i])+",";
				scores=scores+(nf.format(intervalScores[i]))+",";
			}
			
			double[] overlapStat=intervalOverlapStatistics(setCntr,intervalStarts,intervalEnds,intervalScores);
			
			String  name=conGene.getChr()+":"+String.valueOf(conGene.getStart())+"-"+String.valueOf(conGene.getEnd())+"_"+conGene.getOrientation()+"("+ conGene.getName() +")";
			String[] bl =conGene.toBEDArr();
			
			//print line to bed files
			String str=bl[0]+"\t"+bl[1]+"\t"+bl[2]+"\t"+name+
			"\t"+setCntr+"\t"+bl[5]+"\t"+bl[6]+"\t"+bl[7]+"\t"+bl[8]+"\t"+bl[9]+"\t"+bl[10]+
			"\t"+bl[11]+"\t"+intervalCount+"\t"+sizes+"\t"+starts+"\t"+scores+"\t"+overlapStat[0]+"\t"+overlapStat[1]+"\t"+IsoformCntr+"\n";
			
			bw_extendedBed.write(str);
			
			str=bl[0]+"\t"+chromStart+"\t"+chromEnd+"\t"+name+
			"\t"+setCntr+"\t"+bl[5]+"\t"+chromStart+"\t"+chromEnd+"\t"+bl[8]+
			"\t"+intervalCount+"\t"+sizes+"\t"+starts+"\n";
			
			bw_intervalBed.write(str);
			
			
		return null;
	}

	
	private static void printBedLinePerExon(RefSeqGene conGene,
			IntervalTree<LightweightGenomicAnnotation> intervalIsoformTree,
			BufferedWriter bwFasta, BufferedWriter bwKey, int setCntr, int probeSize,int IsoformCntr) throws IOException {
		
		NumberFormat nf = NumberFormat.getInstance();
		nf.setMaximumFractionDigits(3);
		int intervalNum=0;

		String name=conGene.getChr()+":"+String.valueOf(conGene.getStart())+"-"+String.valueOf(conGene.getEnd())+"_"+conGene.getOrientation()+"("+ conGene.getName() +")";
		
		Iterator <LightweightGenomicAnnotation> it=intervalIsoformTree.valueIterator();
		while(it.hasNext()){
			LightweightGenomicAnnotation annot=it.next();
			int Start=annot.getStart();
			int End=annot.getEnd();
			if ((End-Start) >= probeSize){
				intervalNum++;
				double isoScore=annot.getScore();
				int dist=conGene.get3PrimeTranscriptDistance(annot);
				String intervalName=name+"_Interval"+intervalNum;
				//String header=name+"_"+"start="+Start+"_end="+End+"_setNum="+setCntr+"_isoformNum="+IsoformCntr+"_isoformScr="+isoScore+"_3primeDist="+dist;
				int[] blockStarts=new int[1];
				int[] blockEnd=new int[1];
				blockStarts[0]=Start;
				blockEnd[0]=End;
				double newScr=Double.valueOf(nf.format(isoScore))/setCntr;
				RefSeqGene newGene=new RefSeqGene (conGene.getChr(),Start,End,intervalName,newScr,conGene.getOrientation(),blockStarts,blockEnd);
				bwFasta.write(newGene.toBEDwithBedScore()+"\n");
				
				String str= conGene.getChr()+ "\t" + Start +"\t" + End + "\t" + name + "\t" + 
				conGene.getOrientation()+ "\t" + intervalName +"\t" + isoScore + "\t" + IsoformCntr + "\t" + 
				setCntr + "\t" + dist + "\n";
			
				bwKey.write(str);
			}
		}
			
		
		
	}


	//res: <normalized Score in 20% top-overlapping region> <% length overlap of the 20% top-overlapping region >
	private static double[] intervalOverlapStatistics(int setCntr,
			int[] intervalStarts, int[] intervalEnds, double[] intervalScores) {

		double[] res= new double[2];
		
		ArrayList<Double> intervalScoreList=Statistics.toDoubleArrayList(intervalScores);
		ArrayList<Double> unq=Statistics.toUniqueSortedList(intervalScoreList);
		int index=(int) Math.ceil(0.8*unq.size())-1;
		double topScr=unq.get(index);
		ArrayList<Double> topScrs=new ArrayList<Double>();
		double totalLen=0;
		double totalTopLen=0;
		for (int i=0; i<intervalScores.length; i++){
			totalLen+=intervalEnds[i]-intervalStarts[i];
			if(intervalScores[i]>=topScr){
				topScrs.add(intervalScores[i]/setCntr);
				totalTopLen+=intervalEnds[i]-intervalStarts[i];
			}
		}
		res[0]=Statistics.median(Statistics.toDoubleArray(topScrs));
		res[1]=totalTopLen/totalLen;
		
		return res;
	}











	private static ArrayList<Integer> getIntervalBounderies(
			IntervalTree<LightweightGenomicAnnotation> exonTree) {
		
		ArrayList<Integer> rtrn=new ArrayList<Integer>();
		Iterator<LightweightGenomicAnnotation> itAnnot=exonTree.valueIterator();
		while(itAnnot.hasNext()){
			LightweightGenomicAnnotation annot=itAnnot.next();
			if (! rtrn.contains(annot.getStart())) {rtrn.add(annot.getStart());}
			if (! rtrn.contains(annot.getEnd())) {rtrn.add(annot.getEnd());}
		}
		Collections.sort(rtrn);
		 
		return rtrn;
	}


	private static IntervalTree<LightweightGenomicAnnotation> toExonTree(
			IntervalTree<RefSeqGeneWithIsoforms> overlapTree, double scr) {
		
		IntervalTree<LightweightGenomicAnnotation> exonTree=new IntervalTree<LightweightGenomicAnnotation>();
		
		Iterator<RefSeqGeneWithIsoforms> itGene=overlapTree.valueIterator();
		while(itGene.hasNext()){
			RefSeqGeneWithIsoforms gene=itGene.next();
			gene.setBedScore(scr);
			List<LightweightGenomicAnnotation> exons=gene.getScoredExons();
			for (int i=0; i<exons.size(); i++){
				Node<LightweightGenomicAnnotation> curr=exonTree.find(exons.get(i).getStart(), exons.get(i).getEnd());
				if (curr==null)
					exonTree.put(exons.get(i).getStart(), exons.get(i).getEnd(),exons.get(i));
				else
					curr.getValue().setScore(curr.getValue().getScore()+exons.get(i).getScore());
			}
		}
			
		return exonTree;
	}


	private static Map<String, BEDFileParser> uploadBedFiles(String setLst,
			ArrayList<String> nameLst, String trimGene) throws IOException {
		
		HashMap <String,BEDFileParser> annotSet= new HashMap<String,BEDFileParser>(); 
		BufferedReader br = new BufferedReader(new FileReader(setLst));
		String line;
		int bound=-1;
		
		if (! trimGene.equalsIgnoreCase("")) 
			 bound=Integer.valueOf(trimGene);
		
		 while((line = br.readLine()) != null) {
				line = line.trim();
				String [] lineSplit = line.split("\t");
				BEDFileParser set;
				if (bound == -1)
					set = new BEDFileParser(lineSplit[0]);
				else
					set= new BEDFileParser(lineSplit[0],bound);
				annotSet.put(lineSplit[1], set);
				nameLst.add(lineSplit[1]);
		 }		
		
		return annotSet;	
		
	}



	private static AnnotationReader<? extends GenomicAnnotation> CountOverlaps1(AnnotationReader<? extends GenomicAnnotation> superSet, Map<String, AnnotationReader<? extends GenomicAnnotation>> setOfSets) {
		
		Iterator<String> setsIt = setOfSets.keySet().iterator();
		setsIt = setOfSets.keySet().iterator();
		Boolean setToZero=true;
		AnnotationReader<? extends GenomicAnnotation> newSuperSet=superSet;
		while(setsIt.hasNext()) {
			String currName=setsIt.next();
			AnnotationReader<? extends GenomicAnnotation> currSet= setOfSets.get(currName);
			newSuperSet.IncrementScoreIfOverlap(currSet, 0, setToZero);
			setToZero=false;
		}
				
		return newSuperSet;
		
	}
	
	
	private static MatrixWithHeaders GenerateHeatmap(
			AnnotationReader<? extends GenomicAnnotation> superSet,
			Map<String, AnnotationReader<? extends GenomicAnnotation>> setOfSets,
			ArrayList<String> nameLst) {
		
		ArrayList<String> l= new ArrayList<String>();
		
		List<? extends GenomicAnnotation> lst=superSet.getAnnotationList();
		Iterator<? extends GenomicAnnotation> annotIt=lst.iterator();
		while (annotIt.hasNext()){ l.add(annotIt.next().getName()); }
		
		MatrixWithHeaders heatMap=new MatrixWithHeaders(l,	nameLst,0);
		
		Iterator<String> setsIt = setOfSets.keySet().iterator();
		setsIt = setOfSets.keySet().iterator();
		AnnotationReader<? extends GenomicAnnotation> newSuperSet=superSet;
		Boolean setToZero=true;
		while(setsIt.hasNext()) {
			String currName=setsIt.next();
			AnnotationReader<? extends GenomicAnnotation> currSet= setOfSets.get(currName);
			newSuperSet.IncrementScoreIfOverlap(currSet, 0, setToZero);
			updateMatrixByScore(newSuperSet,currName,heatMap);
		}
		
		
		return heatMap;
	}



	private static void updateMatrixByScore(
			AnnotationReader<? extends GenomicAnnotation> Annotations,
			String SetName, MatrixWithHeaders mat) {
			
		List<? extends GenomicAnnotation> lst=Annotations.getAnnotationList();
		Iterator<? extends GenomicAnnotation> annotIt=lst.iterator();
		while (annotIt.hasNext()){
			GenomicAnnotation annot=annotIt.next();
			double scr=annot.getScore();
			if (scr!=0){
				mat.set(annot.getName(), SetName, scr);
			}
		}
		
	}
	
	

	private static BEDFileParser MapNames(String set1In, String set2In) throws IOException {
	
		BEDFileParser bed= new BEDFileParser();
		BEDFileParser set1= new BEDFileParser(set1In);
		BEDFileParser set2= new BEDFileParser(set2In);
		
		int multiNameCntr=0;
		List<RefSeqGene> set1genes=set1.GetGenes();
		String newName;
		
		for(RefSeqGene gene : set1genes)
		{
			IntervalTree<RefSeqGeneWithIsoforms> tree=set2.getOverlappers(gene);
			if (tree.size()>1){
				System.err.println("Has "+tree.size() +"  overlapping genes :" +gene.toBED());
				multiNameCntr++;
			}
			if (! tree.isEmpty()){
			newName=tree.findByIndex(0).getValue().getName();
			gene.setName(newName);
			}
			bed.addRefSeq(gene);
		}
	System.err.println( "\nmultiple mapping was identified for  "+ multiNameCntr + " genes");	
	return bed;
}


	private static Map<RefSeqGene, String> ClassifyLincs(String lincsFile,
		String refGenesFile,int buffer) throws IOException {
		
		Map<RefSeqGene, String> res= new TreeMap<RefSeqGene, String>();
		BEDFileParser lincsBed= new BEDFileParser(lincsFile);
		BEDFileParser refGenesBed= new BEDFileParser(refGenesFile);
		
		
		List<RefSeqGene> lincs=lincsBed.GetGenes();
		String newName;
		
		for(RefSeqGene linc : lincs)
		{
			RefSeqGene ext5primeLinc= linc.getExtended5primeIsoform(buffer);
			RefSeqGene ext3primeLinc= linc.getExtended3primeIsoform(buffer);
			
			//Is this linc an alternative 5' utr of its 3' neighbor?
			int[] dist3prime= distanceTo3primeNeighbor(linc,ext3primeLinc,refGenesBed);
			
			//Is this linc an alternative 3' utr of its 5' neighbor?
			int[] dist5prime= distanceTo5primeNeighbor(linc,ext5primeLinc,refGenesBed);
			
			//Is this a linc transcribed in the opposite direction of a gene's promotor?
			int[] distAntiTranscription = distanceToOppositeStrand5primeNeighbor(linc,ext5primeLinc,refGenesBed);
			
			String[] extField= new String[6];
			
			extField[0]=String.valueOf(dist3prime[0]);
			extField[1]=String.valueOf(dist3prime[1]);
			extField[2]=String.valueOf(dist5prime[0]);
			extField[3]=String.valueOf(dist5prime[1]);
			extField[4]=String.valueOf(distAntiTranscription[0]);
			extField[5]=String.valueOf(distAntiTranscription[1]);
			/*
			linc.addExtraField(String.valueOf(dist3prime[0]));
			linc.addExtraField(String.valueOf(dist3prime[1]));
			linc.addExtraField(String.valueOf(dist5prime[0]));
			linc.addExtraField(String.valueOf(dist5prime[1]));
			linc.addExtraField(String.valueOf(distAntiTranscription[0]));
			linc.addExtraField(String.valueOf(distAntiTranscription[1]));
			*/
			String str=extField[0]+"\t"+extField[1]+"\t"+extField[2]+"\t"+extField[3]+"\t"+extField[4]+"\t"+extField[5];
			res.put(linc, str) ;
		}
	
	return res;
		
	
}


private static  void MapNamesToTable(String set1In, String set2In, String out) throws IOException {
		
		BEDFileParser bed= new BEDFileParser();
		BEDFileParser set1= new BEDFileParser(set1In);
		BEDFileParser set2= new BEDFileParser(set2In);
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		int multiNameCntr=0;
		List<RefSeqGene> set1genes=set1.GetGenes();
		String newName;
		
		for(RefSeqGene gene : set1genes)
		{
			String name1 = gene.getName();
			Iterator<RefSeqGeneWithIsoforms> treeIt= set2.getOverlappers(gene).valueIterator();
			if (treeIt.hasNext()){
				RefSeqGeneWithIsoforms ogene  = treeIt.next();
				for (RefSeqGene g: ogene.getAllIsoforms())
				{
					String name2 = g.getName();
					bw.write(name1+"\t"+name2+"\n");
				}
			}
			
		}
		bw.close();
	
	}

	


	private static int[] distanceTo3primeNeighbor(RefSeqGene linc, RefSeqGene ext3primeLinc, BEDFileParser refGenesBed) {
		
		int[] res=new int[2];
		IntervalTree<RefSeqGeneWithIsoforms> overlapers=refGenesBed.getOverlappers(ext3primeLinc);
		if (overlapers.isEmpty())
			{res[0]=0; res[1]=0;}
		else{
			Iterator<RefSeqGeneWithIsoforms> itG=overlapers.valueIterator();
			res[0]=1; 
			boolean first=true;
			while(itG.hasNext()){
				RefSeqGeneWithIsoforms gene=itG.next();
				if(first)
				{	
					if(linc.getOrientation().equalsIgnoreCase("-"))
						res[1]=linc.getStart()-gene.getEnd();
					else
						res[1]=gene.getStart()-linc.getEnd();
					first=false;
				}
				else
				{ 
					if(linc.getOrientation().equalsIgnoreCase("-"))
						res[1]=Math.min(linc.getStart()-gene.getEnd(), res[1]);
					else
						res[1]=Math.min(gene.getStart()-linc.getEnd(), res[1]);
					
				}
				
			}
		}
		return res;		
	}

	private static int[] distanceTo5primeNeighbor(RefSeqGene linc,RefSeqGene ext5primeLinc, BEDFileParser refGenesBed) {
		
		int[] res=new int[2];
		IntervalTree<RefSeqGeneWithIsoforms> overlapers=refGenesBed.getOverlappers(ext5primeLinc);
		if (overlapers.isEmpty())
			{res[0]=0; res[1]=0;}
		else{
			Iterator<RefSeqGeneWithIsoforms> itG=overlapers.valueIterator();
			res[0]=1; 
			boolean first=true;
			while(itG.hasNext()){
				RefSeqGeneWithIsoforms gene=itG.next();
				if(first)
				{	
					if(linc.getOrientation().equalsIgnoreCase("-"))
						res[1]=gene.getStart()-linc.getEnd();
					else
						res[1]=linc.getStart()-gene.getEnd();
					first=false;
				}
				else
				{ 
					if(linc.getOrientation().equalsIgnoreCase("-"))
						res[1]=Math.min(gene.getStart()-linc.getEnd(), res[1]);
					else
						res[1]=Math.min(linc.getStart()-gene.getEnd(), res[1]);
					
				}
				
			}
		}
		return res;	
		
	}



public static int[] distanceToOppositeStrand5primeNeighbor(
			RefSeqGene linc, RefSeqGene ext5primeLinc, BEDFileParser refGenesBed) {
		
		int[] res=new int[2];
		BEDFileParser lincBED=new BEDFileParser();
		lincBED.addRefSeq(ext5primeLinc);
		
		//get otherSet genes that overlap ext5primeLinc
		BEDFileParser overlapBed=refGenesBed.overlapByGenomicRegion(lincBED,"false");
		
		
		if (overlapBed.getChrTree(linc.getChr())==null)
			{res[0]=0; res[1]=0;}
		else{
			IntervalTree<RefSeqGeneWithIsoforms> overlapers=overlapBed.getChrTree(linc.getChr());
			Iterator<RefSeqGeneWithIsoforms> itG=overlapers.valueIterator();
			
			boolean first=true;
			while(itG.hasNext()){
				RefSeqGeneWithIsoforms gene=itG.next();
				if (! linc.getOrientation().equalsIgnoreCase(gene.getOrientation())) 
				{
					if(first)
					{	
						if(linc.getOrientation().equalsIgnoreCase("-") && linc.getEnd()<gene.getStart())
							{res[1]=gene.getStart()-linc.getEnd();
							res[0]=1; first=false;}
						if(linc.getOrientation().equalsIgnoreCase("+") && linc.getStart()>gene.getEnd())
							{res[1]=linc.getStart()-gene.getEnd();
							res[0]=1; first=false;}
						
					}
					else
					{ 
						if(linc.getOrientation().equalsIgnoreCase("-") && linc.getEnd()<gene.getStart())
							res[1]=Math.min(gene.getStart()-linc.getEnd(), res[1]);
						if(linc.getOrientation().equalsIgnoreCase("+") && linc.getStart()>gene.getEnd())
							res[1]=Math.min(linc.getStart()-gene.getEnd(), res[1]);
						
					}
				}
				
			}
		}
		return res;	
		
		
		
	}

		

	private static void getMedianScrPerQauntile(String file, String setLst,int q, BufferedWriter bw, int rpkmField, int filterField, double lowT ) throws IOException {
		
		Map <String,BEDFileParser> setOfSets=new HashMap <String,BEDFileParser>() ;
		if (!setLst.equals("")){
			
			ArrayList <String> nameLst= new ArrayList <String> ();
			setOfSets= uploadBedFiles(setLst,nameLst, "");
		}
		else if (! file.equals("")){
			BEDFileParser b= new BEDFileParser(file);
			setOfSets.put(file,b);
		}
		else {return;}
		
		LinkedList<String> setsLst=new LinkedList(setOfSets.keySet());
		Collections.sort(setsLst);
		
		for (String f: setsLst){
			BEDFileParser bed=setOfSets.get(f);
			Vector<Double> arr= new Vector<Double>();
			for (RefSeqGene g : bed.GetGenes()){
				if (filterField !=0){
					if (Double.valueOf(g.getExtraFields(filterField)) > lowT)
						continue;
				}
				
				if (rpkmField !=0){
					String str= g.getExtraFields(rpkmField);
					
					//if (str.matches("\\d"))
						arr.add(Double.valueOf(str));
				}
				else
					arr.add(g.getBedScore());
			}
			bw.write(f);
			//System.err.println(arr.size());
			Collections.sort(arr);
			double [] medArr=new double[q];
			int s=(int) Math.floor(arr.size()/q);
			int st=0;
			int en=st+s/2;
			//System.err.println(st+"\t"+en+"\t"+q+"\t"+arr.size());
			for (int i=0; i<(q); i++){
			    medArr[i]=arr.get(en);
			    st=st+s;
			    en=st+s/2;
			    //System.err.println(st+"\t"+en+"\t"+q+"\t"+arr.size());
			    bw.write("\t"+medArr[i]);
			}
			bw.write("\n");
			
			
		}
	
		
	}





	private static void GetConservedExonScrPerGene(String setf, String siphyf,double bl,BufferedWriter bw) throws IOException {

		BEDFileParser set=new BEDFileParser(setf);
		BEDFileParser siphy=readSiphyFile(siphyf);
		BEDFileParser out=new BEDFileParser();
		//Collection <RefSeqGene> c= new LinkedList <RefSeqGene>();
		set.merge();
		int flag=0;
		int FirstExon=0;
		int lastExon=0;
		for (RefSeqGene g: set.GetGenes()){
			Collection<RefSeqGeneWithIsoforms> overlappers=siphy.getOverlappers(g).toCollection();
			double minScr=5;
			RefSeqGeneWithIsoforms currExon=null;
			int exNum=0;
			int selectedExNum=0;
			
			for (RefSeqGeneWithIsoforms r:overlappers){
				exNum++;
				String field=r.getExtraFields(3);
				if ((!field.equals("")) && Double.valueOf(field)>= bl ){
					
					if (Double.valueOf(r.getExtraFields(0)) <= minScr){
						minScr=Math.min(minScr, Double.valueOf(r.getExtraFields(0)));
						currExon=r;
						selectedExNum=exNum;
					}
				}
			}
			if (currExon != null)
			{
				currExon.setBedScore(minScr);
				out.addRefSeq(currExon);
				if (selectedExNum != 1 && selectedExNum != exNum )
					flag++;
				if (selectedExNum == 1)
					FirstExon++;
				if (selectedExNum == exNum )
					lastExon++;	
				//c.add(currExon);
			}
		}
		//out.addRefSeqSet(c);
		out.writeFullBed(bw);
		System.err.println("Number of genes for which a middle exons was chosen is: " + flag);
		System.err.println("Number of genes for which the First exon was chosen is: " + FirstExon);
		System.err.println("Number of genes for which the last exon was chosen is: " + lastExon);

	}



	private static BEDFileParser readSiphyFile(String siphyf) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(siphyf)));
    	String nextLine;
    	BEDFileParser bed=new BEDFileParser();
    	while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
    		String[] tokens=nextLine.split("\t");
			String chr=(tokens[0]);
			int start=new Integer(tokens[1]);
			int stop=new Integer(tokens[2]);
			RefSeqGene g=new RefSeqGene(chr,start,stop);
			for(int j = 3; j < tokens.length; j++) {
				g.addExtraField( tokens[j]);
			}
			bed.addRefSeq(g);
			
    	}
		return bed;
	}

	
	
	private static void removeSingleExons(String in, String informat,
			String outformat, BufferedWriter bw, String src, int sizeT) throws IOException {
		
		BEDFileParser bed= new BEDFileParser();
		if (informat.equalsIgnoreCase("bed"))
			bed=new BEDFileParser(in);
		else
			bed.loadGTF(new File(in),"");
		
		BEDFileParser outbed=bed.removeSmallTranscripts(2,sizeT);
		
		if (outformat.equalsIgnoreCase("bed"))
			outbed.writeFullBed(bw);
		else
			outbed.bed2CuffGtf(bw, "", true);
		
		
	}
	
	

	private static void compareAssembleyAndRefLength(String in,
		String informat, String ref_bed, BufferedWriter bw, String chr ) throws IOException {
	
		BEDFileParser bed= new BEDFileParser();
		if (informat.equalsIgnoreCase("bed"))
			bed=new BEDFileParser(in);
		else
			bed.loadGTF(new File(in),chr);
		
		BEDFileParser ref=new BEDFileParser(ref_bed);
		bw.write("refGeneName\tIsoformName\tgeneLen+\tisoLen\tleftEndFracCovered\trightEndFracCovered\tOrientation\tIsFullyCompatible\tleftNoUTRFracCovered\trightNoUTRFracCovered\t\n");
		
		for (RefSeqGene g: ref.GetGenes()){
			Iterator<RefSeqGeneWithIsoforms> it=bed.getOverlappers(g).valueIterator();
			while(it.hasNext()){
				Collection<RefSeqGene> iso_set=it.next().getAllIsoforms();
				for (RefSeqGene iso:iso_set){
					if (iso.getNumExons()<=1 && g.getNumExons() >1) continue; 
					Alignments[] g_ex =g.getExons();
					Alignments[] iso_ex =iso.getExons();
					int orient =1;
					if (g.getOrientation().equals("-")) orient=-1;
					if (g.isFullyCompatible(iso)){
						int isoLen=iso.getTranscriptLength();
						int genLen=g.getTranscriptLength();
						int[] isoEx= iso.getExonSizes();
						int[] gEx=g.getExonSizes();
						Double leftEnd=new Double(isoEx[0])/new Double(gEx[0]);
						Double rightEnd=new Double(isoEx[g.getNumExons()-1])/new Double(gEx[g.getNumExons()-1]);
						bw.write(g.getName()+"\t"+iso.getName()+"\t"+genLen+"\t"+isoLen+"\t"+leftEnd+"\t"+rightEnd+"\t"+orient+"\t"+1+"\t");
					}
					else { //NOT COMPATIBLE

						int isoLen=iso.getTranscriptLength();
						int genLen=g.getTranscriptLength();
						
						double leftEnd =0.0;
						double rightEnd =0.0;
						if (iso_ex[0].getStart() < g_ex[0].getEnd())
							leftEnd=new Double(Math.abs(iso_ex[0].getStart() - g_ex[0].getEnd()))/new Double(Math.abs(g_ex[0].getStart() - g_ex[0].getEnd()));
						if (iso_ex[iso.getNumExons()-1].getEnd() > g_ex[g.getNumExons()-1].getStart())
							rightEnd=new Double(iso_ex[iso.getNumExons()-1].getEnd() - g_ex[g.getNumExons()-1].getStart())/new Double(g_ex[g.getNumExons()-1].getEnd() - g_ex[g.getNumExons()-1].getStart());
						bw.write(g.getName()+"\t"+iso.getName()+"\t"+genLen+"\t"+isoLen+"\t"+leftEnd+"\t"+rightEnd+"\t"+orient+"\t"+0+"\t");
					}	
					
					double leftEnd =0.0;
					double rightEnd =0.0;
					//UTR calc:
					if (iso_ex[0].getStart() < g_ex[0].getEnd())
						leftEnd=new Double(Math.abs(iso_ex[0].getStart() - g_ex[0].getEnd()))/new Double(Math.abs(g.getCDSRegion().getStart() - g_ex[0].getEnd()));
					if (iso_ex[iso.getNumExons()-1].getEnd() > g_ex[g.getNumExons()-1].getStart()){
						rightEnd=new Double(iso_ex[iso.getNumExons()-1].getEnd() - g_ex[g.getNumExons()-1].getStart())/new Double(  g.getCDSRegion().getEnd()- g_ex[g.getNumExons()-1].getStart());
					
					Alignments a= g.getCDSRegion();	
					a.getEnd();	
					}
					bw.write(leftEnd+"\t"+rightEnd+"\n");
				}
			
			
			}
		}
	}
	
	private static void GetK4K36Overlaps(String in, String k4, String k36,String outprefix ,int promoterInterval,String genesF, int numPerm, String[] infiles) throws IOException {

		BEDFileParser bed=new BEDFileParser(in);
		BEDFileParser genesBed= new BEDFileParser(genesF);
		Map<String, IntervalTree<Alignments>>  k4TreeMap= BEDFileParser.loadAlignmentDataToTree(new File(k4));
		Map<String, IntervalTree<Alignments>>  k36TreeMap= BEDFileParser.loadAlignmentDataToTree(new File(k36));
		
		BEDFileParser[] trueRes=GetK4K36Overlaps (bed,genesBed,k4TreeMap,k36TreeMap,outprefix,promoterInterval , true);
		
		int cnt1=0;
		int cnt2=0;
		BEDFileParser b1=trueRes[0];
		BEDFileParser b2=trueRes[1];
		 b1.merge(); b2.merge();
		int t1=b1.GetGenes().size();
		int t2=b2.GetGenes().size();
		
		
		if (numPerm>0){
			int [] val1 = new int[numPerm];
			int [] val2 = new int[numPerm];
			int i=0;
			String annotTofilter=infiles[0]; String centromers=infiles[1]; String chrSizes=infiles[2];
			transcriptsNullModel nullmodel = new transcriptsNullModel (new BEDFileParser (annotTofilter),new BEDFileParser(centromers),chrSizes); 
	    	nullmodel.makeRandTranscriptList(numPerm, bed);
	    	for (BEDFileParser  randbed: nullmodel.getRandomSets()){
	    		BEDFileParser[] randRes=GetK4K36Overlaps (randbed,genesBed,k4TreeMap,k36TreeMap,outprefix,promoterInterval , true);
	    		BEDFileParser rb1=randRes[0];
	    		BEDFileParser rb2=randRes[1];
	    		 rb1.merge(); rb2.merge();
	    		int rt1=rb1.GetGenes().size();
	    		int rt2=rb2.GetGenes().size();
	    		if (rt1>=t1)
	    			cnt1++;
	    		if (rt2>=t2)
	    			cnt2++;
	    		val1[i]=rt1;
	    		val2[i]=rt2;
	    		i++;
	    	}
	    	
	    	BufferedWriter bw = new BufferedWriter(new FileWriter(outprefix + "_countsOfRandSets.txt"));
			for (int j=0; j<i;j++) {bw.write(val1[j]+"\t");} 
			bw.write("\n");
			for (int j=0; j<i;j++) {bw.write(val2[j]+"\t");} 
			bw.write("\n");
			bw.close();
			
		}
		
		double p1= new Double(cnt1)/new Double(numPerm);
		double p2= new Double(cnt2)/new Double(numPerm);
		
		System.err.println("Input set size : " + bed.GetGenes().size());
		System.err.println("Number K4K36 in real set : " + t1 + " independent : " + t2);
		System.err.println("Significance by perm pvalue  : " + p1 + " independent : " + p2);
		
		
		
	}
	
    private static BEDFileParser[] GetK4K36Overlaps(BEDFileParser bed, BEDFileParser genesBed, Map<String, IntervalTree<Alignments>>  k4TreeMap, Map<String, IntervalTree<Alignments>> k36TreeMap , String outprefix ,int promoterInterval, boolean toWrite) throws IOException {
    	
    	BEDFileParser res=new BEDFileParser();
		BEDFileParser res2=new BEDFileParser(); //will maintain genes that their k4-k36 overlap does not overlap a gene
		
		Iterator <String > chrIt =bed.getChromosomeIterator();
		while(chrIt.hasNext()){
			String chr = chrIt.next();
			IntervalTree<Alignments>  k4tree=new IntervalTree<Alignments> ();
			IntervalTree<Alignments>  k36tree=new IntervalTree<Alignments> ();
			if (k4TreeMap.containsKey(chr))
				  k4tree=k4TreeMap.get(chr);
			if (k36TreeMap.containsKey(chr))
				 k36tree=k36TreeMap.get(chr);
			Iterator<RefSeqGeneWithIsoforms> gIt= bed.getChrTree(chr).valueIterator();
			while(gIt.hasNext()){
				RefSeqGeneWithIsoforms g= gIt.next();
				for (RefSeqGene iso:g.getAllIsoforms()){
					boolean hask4=false;
					boolean hask36=false;
					boolean hasGene=false;
					Alignments promoter= iso.get5primeRegion(promoterInterval,true);
					if (! k4tree.isEmpty()){
						Iterator<Node<Alignments>> k4Overlap=k4tree.overlappers(promoter.getStart(), promoter.getEnd());
						if(k4Overlap.hasNext()){
							hask4=true;
							while(k4Overlap.hasNext()){
								Alignments region=k4Overlap.next().getValue();
								if (! genesBed.getOverlappers(new RefSeqGene(region)).isEmpty())
								{
									hasGene=true; break;
								}
							}
							
						}
					}
					if (! k36tree.isEmpty()){
						Iterator<Node<Alignments>> k36Overlap=k36tree.overlappers(iso.getStart(), iso.getEnd());
						if(k36Overlap.hasNext()){
							hask36=true;
							while(k36Overlap.hasNext()){
								Alignments region=k36Overlap.next().getValue();
								if (! genesBed.getOverlappers(new RefSeqGene(region)).isEmpty())
								{
									hasGene=true; break;
								}
							}
						}
					}
					if (hask4 & hask36){
						res.addRefSeq(iso);
						if (hasGene==false)
							res2.addRefSeq(iso);
					}
				}
				
			}
		}
		if (toWrite) {
			String out = outprefix  + ".bed";
			res.writeFullBed(out);
			out = outprefix  + ".independent.bed";
			res2.writeFullBed(out);
		}
		BEDFileParser[] resArr = new BEDFileParser[2];
		resArr[0]=res; resArr[1]=res2;
		return  resArr;
	}
	
	
	
	

	private static void GenerateGeneStatMat(String refSetGTF, String neighborMap, String expMatF,
	String isoformStatMat,String JSspGct,String k4k36LincsBed, BufferedWriter bw) throws IOException, ParseException {
		
		BEDFileParser uniqTconRef=new BEDFileParser();
		uniqTconRef.loadGTF(new File(refSetGTF),"");
		HashMap<String, List<String>> geneTconMap=new HashMap<String, List<String>>();
		Map<String, String> TconGeneMap=new HashMap<String,String>();
		Map<String , RefSeqGeneWithIsoforms> geneNameGeneMap= new HashMap<String , RefSeqGeneWithIsoforms>();
		Map<String , GeneAnnotationInfo> geneSummaryMap= new HashMap<String , GeneAnnotationInfo>();
		MatrixWithHeaders JSmat=new MatrixWithHeaders(JSspGct);
		MatrixWithHeaders expMat=new MatrixWithHeaders(expMatF);
		BEDFileParser k4k36Lincs= new BEDFileParser(k4k36LincsBed);
		
		//1.
		BEDFileParser uniqGeneRef= AssociateGeneTranscript(uniqTconRef,geneTconMap,TconGeneMap,geneNameGeneMap);
		System.err.println("Total genes: "+ uniqGeneRef.GetGenes().size());
		
		//2.init geneSummaryMap
		for (String name : geneNameGeneMap.keySet())
			geneSummaryMap.put(name,new GeneAnnotationInfo(geneNameGeneMap.get(name)));
		
		//3.parse isoformStat
		loadIsoformStatTabToGeneSummaryMap(isoformStatMat,geneSummaryMap,TconGeneMap);
		
		//4.load to neighbors
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(neighborMap)));
		String nextLine;
		boolean first=true;
		String[] NeighborHeader=null;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			String[] a=nextLine.split("\t");
			String name= a[3];
			
			if (first){
				NeighborHeader=a;
				first=false;
			}
			else{
				if(geneSummaryMap.containsKey(name)){
					//System.err.println(name+"\t"+nextLine);
					GeneAnnotationInfo info=geneSummaryMap.get(name);
					info.setNeighborsLine(nextLine);
				}
			}
		}
		
		//5.parse exp,//6.parse JS
		for (String name : geneNameGeneMap.keySet()){
			String[] tmpCenter= new String[1];
			double[] e= getMaxGctRow(expMat,name,tmpCenter);
			double[] j= getMaxGctRow(JSmat,name,tmpCenter);
			GeneAnnotationInfo info=geneSummaryMap.get(name);
			info.setMaxExp(e[0]);
			info.setMaxJSspecificity(j[0]);
			info.setMaxJSspecificity_centroid(tmpCenter[0]);
		}
		
		//7. k4k36Lincs
		for (RefSeqGene iso: k4k36Lincs.GetGenes()){
			if (TconGeneMap.containsKey(iso.getName())){
				String geneName= TconGeneMap.get(iso.getName());
				if (geneSummaryMap.containsKey(geneName))
					(geneSummaryMap.get(geneName)).setk4k36Lincs(1);
			}
		}
		
		//8. longest exon-chain isoform
		SetLongestIsoInGeneSummaryMap(uniqTconRef,geneTconMap,TconGeneMap,geneSummaryMap);
		
		//Write:
		bw.write("chr\tstart\tend\tname\tscore\tstrand\ttStart\ttEnd\trgb\tblockCnt\tblockSizes\tblockStarts");
		bw.write("\tNumIsoforms\tHighConfidenceSet\tCompatibleBetweenAssemblies\tCompatibleBetweenAssemblyAndExternal\tCompatibleBetweenSets");
		bw.write("\tMaxUniformityScore\tCoverage\tScanStatExp\tK4K36\tK4K36Independan\tK4K36lincs");
		bw.write("\tMaxExp\tMaxJSspecificity\tMaxJSspecificity_cov\tNR\tUCSC\tGENCODE\tSEQ");
		bw.write("\tLongChainIsoExonNum\tLongChainIsoTranscriptLength");
		for (int i=12;i<NeighborHeader.length;i++)
			bw.write("\t"+NeighborHeader[i]);
		bw.write("\n");
		
		String emptyArr="";
		for (int i=0;i<NeighborHeader.length-1 ; i++) emptyArr=emptyArr+"NaN"+"\t";
		emptyArr=emptyArr+"NaN";
		for (String name : geneNameGeneMap.keySet()){
			geneSummaryMap.get(name).writeInfo(bw,emptyArr);	
		}
	}




	private static void makeDetectionByGene (String expMatF, String refgtf, BufferedWriter bw) throws IOException, ParseException {
	 
		BEDFileParser uniqTconRef=new BEDFileParser();
		uniqTconRef.loadGTF(new File(refgtf),"");
		HashMap<String, List<String>> geneTconMap=new HashMap<String, List<String>>();
		Map<String, String> TconGeneMap=new HashMap<String,String>();
		Map<String , RefSeqGeneWithIsoforms> geneNameGeneMap= new HashMap<String , RefSeqGeneWithIsoforms>();
		Map<String ,ArrayList<Double>> geneSummaryMap= new HashMap<String , ArrayList<Double>>();
		MatrixWithHeaders expMat=new MatrixWithHeaders(expMatF);
		
		//1.
		BEDFileParser uniqGeneRef= AssociateGeneTranscript(uniqTconRef,geneTconMap,TconGeneMap,geneNameGeneMap);
		System.err.println("Total genes: "+ uniqGeneRef.GetGenes().size());
				
		//2.init geneSummaryMap
		for (String name : geneNameGeneMap.keySet()){
			ArrayList<Double>  d=new ArrayList<Double>();
			for (int i=0; i<expMat.columnDimension() ;i++)
				d.add(new Double(0));
			geneSummaryMap.put(name,d);
		}
		
		//3. go over all transcripts and update the corresponding gene
		for (String tcon: TconGeneMap.keySet()){
			if (expMat.containsRow(tcon )&&  geneSummaryMap.containsKey( TconGeneMap.get(tcon))){
			double []tconExp= expMat.getRow(tcon);
			ArrayList<Double> a= geneSummaryMap.get( TconGeneMap.get(tcon));
			
			for (int i=0; i<tconExp.length ; i++){
				a.add(i, a.get(i)+tconExp[i]);
			}
			geneSummaryMap.put( TconGeneMap.get(tcon),a);
			}
		}
	 
		//4. write
		bw.write ("name\tdescription");
		for (int i=0 ; i<expMat.getColumnNames().size(); i++)
			bw.write("\t"+expMat.getColumnNames().get(i));
		bw.write("\n");
		for (String gene : geneSummaryMap.keySet()){
			bw.write(gene+"\t"+gene);
			ArrayList<Double> a = geneSummaryMap.get(gene);
			for (int i=0 ; i<a.size(); i++)
				bw.write("\t"+a.get(i));
			bw.write("\n");
		}
		return;
	}





	private static void SetLongestIsoInGeneSummaryMap(
			BEDFileParser uniqTconRef, 	HashMap<String, List<String>> geneTconMap,
			Map<String, String> tconGeneMap, Map<String, GeneAnnotationInfo> geneSummaryMap) {

		for(String geneName: geneSummaryMap.keySet()){
			GeneAnnotationInfo info=geneSummaryMap.get(geneName);
			HashMap<String,RefSeqGene> isoMap=uniqTconRef.getNameGeneMap();
			if (geneTconMap.containsKey(geneName)){
				Locus locus=new Locus(info.getMergedRef());
				for (String isoName:geneTconMap.get(geneName)){
					if (isoMap.containsKey(isoName))
						locus.addIsoform(new RefSeqGeneWithIsoforms(isoMap.get(isoName)), "def");
						
				}
				RefSeqGene longestIso= locus.getLongestIntornChainIso("def");
				info.setlongestIso(longestIso);
			}
		}
		
	}





	private static void loadIsoformStatTabToGeneSummaryMap(
			String isoformStatMat,Map<String, GeneAnnotationInfo> geneSummaryMap,
			Map<String, String> tconGeneMap) throws IOException {

		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(isoformStatMat)));
		String nextLine;
		Map<String,Integer> headerIndex= new HashMap<String,Integer>();
		
		String[] NeighborHeader=reader.readLine().split("\t");
		for (int i=0; i< NeighborHeader.length ;i++){
			headerIndex.put(NeighborHeader[i],new Integer(i));
		}
		
		
		//("chr\tstart\tend\tname\tscore\tstrand\ttStart\ttEnd\trgb\tblockCnt\tblockSizes\tblockStarts");
		//Isoforms\tHighConfidenceSet\tCompatibleBetweenAssemblies\tCompatibleBetweenAssemblyAndExternal\tCompatibleBetweenSets");
		//uniformityScore\tCoverage\tScanStatExp\tK4K36\tK4K36Independent\tK4K36lincs");
	
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			String[] a=nextLine.split("\t");
			String isoName= a[3];
			
			if (tconGeneMap.containsKey(isoName)){
				String geneName= tconGeneMap.get(isoName);
				GeneAnnotationInfo info=geneSummaryMap.get(geneName);
				info.incrementNumIso();
				info.updateHighConfidence(headerIndex.get("HighConfSet"),a); 
				info.updateCompatibleBetweenAssemblies(headerIndex.get("CompatibleBetweenAssemblies"),a);
				info.updateCompatibleBetweenAssemblyAndExternal(headerIndex.get("CompatibleBetweenAssemblyAndExternal"),a); 
				info.updateCompatibleBetweenSets(headerIndex.get("CompatibleBetweenSets"),a); 
				info.updateUniformityScore(headerIndex.get("MaxUniformityScore"),headerIndex.get("Coverage"),a); 
				info.updateScanStatExp(headerIndex.get("ScanStatExp"),a); 
				info.updateK4K36(headerIndex.get("K4K36"),a); 
				info.updateK4K36Independent(headerIndex.get("K4K36Independent"),a);
				info.updateNR(headerIndex.get("Refseq_NR_lincs"),a);
				info.updateGENCODE(headerIndex.get("gencode_lincs"),a); 
				info.updateUCSC(headerIndex.get("ucsc_lincs"),a); 
				info.updateSEQ(detectedBySeq(a,headerIndex));
				
				
				
			}
		}
		
		
	}



	private static boolean detectedBySeq(String[] a,Map<String, Integer> headerIndex) {

		for(String head: headerIndex.keySet()){
			if (head.contains("cufflinks") || head.contains("scripture")){
				int  ix=headerIndex.get(head);
				if(Double.valueOf(a[ix]) > 0)
					return true;
			}
		}
		return false;
	}





	private static double[] getMaxGctRow(MatrixWithHeaders gct, String name, String[] tmpname) {

		double[] res=new double[2];
		
		res[0]=Double.NaN; res[1]=0;
		if(gct.containsRow(name)){
			double[] r= gct.getRow(name) ;
			double minval= Statistics.min(r);
			res[0]=minval;
			for (int i=0; i<r.length; i++){
				if(r[i]>=res[0]){
					res[0]=r[i]; res[1]=i;
				}
			}
		}
		tmpname[0]=gct.getColoumnName((int) res[1]);
		return res;
	}
	
	
	//output: transcriptBed\tnumMapped\tpctMapped\tCompatible\tMappedToCoding
	private static void ParseTransmap(String in, String refgeneF, String estF,
			String mrnaF, String ucscF, String ucscCodingMap,String transcriptSpecieMapF, String ChainTabF,String ChainSpecie, BufferedWriter bw, String outPrefix, String chainSpecieTranscriptsFile) throws IOException {

		BEDFileParser refbed = new BEDFileParser(in);
		HashMap<String, Locus> geneVsMapped = new HashMap<String,Locus>();
		HashMap<String, Integer> isCoding = new HashMap<String,Integer>();
		for (RefSeqGene g: refbed.GetGenes()){
			geneVsMapped.put(g.getName(),new Locus(new RefSeqGeneWithIsoforms(g))) ;
		}
		HashMap<String , String> setLstFile = new HashMap<String , String>(); 
		setLstFile.put("refgene",refgeneF);
		setLstFile.put("est",estF);
		setLstFile.put("mRNA",mrnaF);
		setLstFile.put("ucsc",ucscF);
		
		HashMap<String,String> plainNames = new HashMap<String,String> ();
		
		LinkedList<String> setLst= loadCompatibleToRefLociIsoformMap(refbed,setLstFile,geneVsMapped,false,false);
		HashMap<String,String> transmapTranscriptNames = getMappedTranscriptIds(geneVsMapped,plainNames,true);
		HashMap <String,String> ucscCategory= loadCodingMap(ucscCodingMap,transmapTranscriptNames,plainNames); 
		HashMap <String,String> transcriptSpecieMap = loadLociSpecieMap(transcriptSpecieMapF,transmapTranscriptNames,true); 
		System.err.println("loaded coding map");
		BEDFileParser NetChainTab = loadNetChainTab(ChainTabF,refbed);
		System.err.println("loaded net chain");
		BEDFileParser otherSpecieBed = new BEDFileParser();
		BEDFileParser chainSpecieTranscripts = loadSpTranscripts(chainSpecieTranscriptsFile,plainNames);
		System.err.println("loaded transcripts from other specie");
		HashMap<String,RefSeqGene> cpTconMap=chainSpecieTranscripts.getNameGeneMap();
		HashMap <String , Set<RefSeqGene>> refSpMap = new  HashMap <String , Set<RefSeqGene>>();
				
		bw.write("#chr\tstart\tend\tname\tscore\tstrand\ttStart\ttEnd\trgb\tblockCnt\tblockSizes\tblockStarts\t");
		bw.write("transmap_numMapped\ttransmap_pctMapped\ttransmap_pctGenomeMapped\ttransmap_maxExonsCovered\ttransmap_FullyCompatible\ttransmap_PartiallyComaptible" +
				"\ttransmap_MappedToCoding\ttransmap_MappedToNearCoding\tPositionBias\tMappedToMouse\tMm9NetLevel\tMm9NetType\tMm9Unq2OrthoNetType\tNetInfo\n");
		for (String g: geneVsMapped.keySet()){
			Locus loci=geneVsMapped.get(g);
			int numIso=loci.getAllIsoforms("").size();
			Set <RefSeqGene> orthologs = new TreeSet <RefSeqGene>();
			int foundInRefSpecie = isInRefSpecie (loci.getAllIsoforms(""),transcriptSpecieMap,ChainSpecie, orthologs);
			System.err.println("Number of ortholgs " + orthologs.size());
			updateMapSpTranscripts(orthologs, cpTconMap, chainSpecieTranscripts, refSpMap,otherSpecieBed,g,loci.getReferenceTranscript() );
			
			HashMap<String,Integer> numPartialCompMap=loci.getIsPartiallyComaptibleWithRefPerSet();
			HashMap<String,Integer> numFullCompMap=loci.getIsComaptibleWithRefPerSet();
			int tnumFullCompMap=0; 
			int tnumPartialCompMap=0;
			for (String s:numPartialCompMap.keySet() ){
				 tnumFullCompMap+=numFullCompMap.get(s); 
				tnumPartialCompMap+=numPartialCompMap.get(s);
			}
			double maxExonsCovered = loci.getMaxExonsCovered();
			double maxPctCovered = loci.getMaxPctCovered();
			double maxPctGenomeCovered = loci.getMaxPctGenomeCovered();
			int posBias =  loci.isTerminiCoverageBias(); //-1,0,1 to a 5'/none/3' bias 
			int coding=0;
			int nearCoding=0;
			for(RefSeqGene iso : loci.getSetIsoforms("ucsc")){
				String isoname=iso.getName();
				boolean sameOrientation= (iso.getOrientation().equalsIgnoreCase(loci.getReferenceTranscript().getOrientation()));
				String[] clnName =isoname.split("-");
				String ctg = ucscCategory.get(clnName[0]);
				if (ctg != null && sameOrientation){
					if (ctg.equalsIgnoreCase("coding"))
						coding=1;
					if (ctg.equalsIgnoreCase("nearCoding"))
						nearCoding=1;
				}
				
			}
			Integer[] netLevel= new Integer[2];
			Integer[]  netType = new Integer [2];
			String str = findAlignmentType( loci.getReferenceTranscript(),NetChainTab,netLevel, netType); //TODO - WHAT IS THE TYPE OF ALIGNEMENT ON THIS LEVEL- IF GAP ON TOP WHAT IS NEXT
			String str2 = findNeighborAlignmentType( loci.getReferenceTranscript(), orthologs,NetChainTab,netLevel, netType); //TODO - WHAT IS THE TYPE OF ALIGNEMENT ON THIS LEVEL- IF GAP ON TOP WHAT IS NEXT
			
			/*//Do not use refseq NM as representing Coding as it also includes some predicated genes, UCSC is suppose to cover all refseq 
			 * in its coding/noncoding annotations
			for(RefSeqGene iso : loci.getSetIsoforms("refgene")){
				String isoname=iso.getName();
				String[] clnName =isoname.split("-");
				if (clnName[0].contains("NM"))
					coding=true;
			}
			*/
			String outStr=  loci.getReferenceTranscript().toBED() +"\t" + numIso +"\t"+ maxPctCovered +"\t"+maxPctGenomeCovered+"\t"+maxExonsCovered+"\t" ;
			outStr += tnumFullCompMap+"\t"+tnumPartialCompMap+"\t"+coding+"\t"+nearCoding+"\t"+posBias+"\t"+foundInRefSpecie+"\t"+netLevel[0]+"\t"+netType[0]+"\t"+netType[1]+"\t"+str+"\n";
			
			
			bw.write( outStr);
			//System.err.println( outStr);
		}
		
		otherSpecieBed.writeFullBed(outPrefix+ChainSpecie+"_mappedTranscripts.bed");
		BufferedWriter bw2 = new BufferedWriter(new FileWriter(outPrefix+ChainSpecie+"_mappedTranscriptPairs.tab"));
		BEDFileParser bed2 = new BEDFileParser();
		HashMap<String ,RefSeqGene> refGeneNames = refbed.getNameGeneMap();
		for (String g: refSpMap.keySet()){
			bed2.addRefSeq(refGeneNames.get(g));
			for (RefSeqGene g2: refSpMap.get(g)){
				String[] newname=g2.getName().split("\\.");
				bw2.write(g+"\t"+newname[0]+"\t"+g2.getBedScore()+"\t"+g2.getAttribute("spB_pctAln")+"\n");
			}
		}
		bw2.close();
		bed2.writeFullBed(outPrefix+"Ref"+"_mappedTranscripts.bed");
	}

	//Given transmapped data , and a file with species identifier output a gct of the mapping between 
	//species to transcripts
	private static void makeSpHeatMap(String in, String refgeneF, String estF,
			String mrnaF, String ucscF, String transcriptSpecieMapF,String speciesNamesF,BufferedWriter bw) throws IOException{
		
		BEDFileParser refbed = new BEDFileParser(in);
		HashMap<String, Locus> geneVsMapped = new HashMap<String,Locus>();
		HashMap<String, Integer> isCoding = new HashMap<String,Integer>();
		for (RefSeqGene g: refbed.GetGenes()){
			geneVsMapped.put(g.getName(),new Locus(new RefSeqGeneWithIsoforms(g))) ;
		}
		HashMap<String , String> setLstFile = new HashMap<String , String>(); 
		setLstFile.put("refgene",refgeneF);
		setLstFile.put("est",estF);
		setLstFile.put("mRNA",mrnaF);
		setLstFile.put("ucsc",ucscF);
		
		HashMap<String,String> plainNames = new HashMap<String,String> ();
		
		LinkedList<String> setLst= loadCompatibleToRefLociIsoformMap(refbed,setLstFile,geneVsMapped,false,false);
		HashMap<String,String> transmapTranscriptNames = getMappedTranscriptIds(geneVsMapped,plainNames,true);
		HashMap <String,String> transcriptSpecieMap = loadLociSpecieMap(transcriptSpecieMapF,transmapTranscriptNames,true); 
		HashMap<String,String> foundGenes= new HashMap<String,String> ();
		//Hash all species-gene interactions
		HashMap <String,Set<String>> speciesGeneMap = new HashMap <String,Set<String>>();
		for (String g: geneVsMapped.keySet()){
			Locus loci=geneVsMapped.get(g);
			updateSpecieGeneMap(speciesGeneMap,loci.getAllIsoforms(""),transcriptSpecieMap,g,foundGenes);
		}
		
		//Read speciesNames 
		LinkedList<String> spNames=new LinkedList();
		BufferedReader reader= new BufferedReader( new FileReader(speciesNamesF));
		String nextLine;
		HashMap<String, String> spNameMap = new HashMap<String, String>();
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			String[] a = nextLine.split("\t");
			if(! spNameMap.containsValue(a[1]))
				spNames.add(a[1]);
			spNameMap.put(a[0],a[1]);
		}
		
		//Make specie gene GCT
		LinkedList<String> species=new LinkedList(speciesGeneMap.keySet());
		LinkedList<String> genes=new LinkedList(foundGenes.keySet());
		Matrix mat=new Matrix( genes.size(),spNames.size(),0);
		MatrixWithHeaders gct=new MatrixWithHeaders(mat,genes,spNames);
		
		for (String sp: species){
			String spn=spNameMap.get(sp);
			if (spn  != null){
				for (String g: speciesGeneMap.get(sp) ){
					gct.set(g, spn, 1);
				}
			}
		}
	    gct.writeGCT(bw);
		
	}


	private static BEDFileParser loadSpTranscripts(
			String chainSpecieTranscriptsFile,
			HashMap<String, String> transmapTranscriptNames) throws IOException {
		
		BEDFileParser res = new BEDFileParser();
		BufferedReader reader= new BufferedReader( new FileReader(chainSpecieTranscriptsFile));
		String nextLine;
		nextLine = reader.readLine();
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			String[] a = nextLine.split("\t");
			if(transmapTranscriptNames.containsKey(a[3]))
				res.addRefSeq(new RefSeqGene(nextLine,false));
		}
		reader.close();
		
		return res;
	}



	private static void updateMapSpTranscripts(Set<RefSeqGene> orthologs,
			HashMap<String,RefSeqGene> spTconMap, BEDFileParser chainSpecieTranscripts,
			HashMap<String, Set<RefSeqGene>> refSpMap,BEDFileParser otherSpecieBed, String refGeneName, RefSeqGene refGene) {
		
		Set<RefSeqGene> spSet = new TreeSet<RefSeqGene>();
		for (RefSeqGene g : orthologs) {
			String[] arr = g.getName().split("\\.");
			String name = arr[0];
			
			if (spTconMap.containsKey(name)){
				RefSeqGene gene = new RefSeqGene( spTconMap.get(name));
				chainSpecieTranscripts.addRefSeq(gene);
				double d=refGene.percentOverlapping(g);
				double spb_d=g.percentOverlapping(refGene);
				
				gene.setBedScore(d);
				gene.addAttribute("spB_pctAln", String.valueOf(spb_d));
				spSet.add(gene);
			}
			else
				System.err.println ("Did not find " +name + " " + g.getName());
			
		}
		refSpMap.put(refGeneName,spSet);
		otherSpecieBed.addRefSeqSet(spSet);
		
		
	}



	//Find the type of net alignments that overlap the transcript : top/inv/nonSyn/Syn
	private static String findAlignmentType(
			RefSeqGene referenceTranscript,
			BEDFileParser netChainTab, Integer[] netLevel, Integer[] netType) {
		
		String all="";
		netLevel[0]=0;
		netType[0]=0;
		Iterator<RefSeqGeneWithIsoforms> gIt= netChainTab.getOverlappers(referenceTranscript).valueIterator();
		while(gIt.hasNext()){
			for(RefSeqGene g: gIt.next().getAllIsoforms()){
				String v = g.getAttribute("netLevel");
				String t =  g.getAttribute("netType");
				if (v != null & t != null){
					all = all + v+"_"+t+";"; 
					if ( t.equalsIgnoreCase("NonSyn") || t.equalsIgnoreCase("Inv") ){
						netLevel[0]= Integer.valueOf(v);
						netType[0]= 2;
						if (t.equalsIgnoreCase("NonSyn"))
							netType[0]= 1;
						
					}
				}
			}
		}
		return all;
	}


	//Find the type of net alignments that overlap the regions spanning the mouse transmapped  transcript
	//but not the human reference transcript : top/inv/nonSyn/Syn
	//update netLevel/netType on index 1
	private static String findNeighborAlignmentType(
			RefSeqGeneWithIsoforms referenceTranscript, Collection <RefSeqGene> orthologs,
			BEDFileParser netChainTab, Integer[] netLevel, Integer[] netType) {
		
		String all="";
		netLevel[1]=0;
		netType[1]=0;
		for ( RefSeqGene or : orthologs) {
			Iterator<RefSeqGeneWithIsoforms> gIt= netChainTab.getOverlappers(or).valueIterator();
			while(gIt.hasNext()){
				for(RefSeqGene g: gIt.next().getAllIsoforms()){
					String v = g.getAttribute("netLevel");
					String t =  g.getAttribute("netType");
					if (v != null & t != null){
						all = all + v+"_"+t+";"; 
						if ( (t.equalsIgnoreCase("NonSyn") || t.equalsIgnoreCase("Inv")) && !g.overlaps(referenceTranscript) ){
							netLevel[0]= Integer.valueOf(v);
							netType[0]= 2;
							if (t.equalsIgnoreCase("NonSyn"))
								netType[0]= 1;
							
						}
					}
				}
			}
		}
		return all;
	}

	private static int isInRefSpecie(Collection<RefSeqGene> allIsoforms,
			HashMap<String, String> transcriptSpecieMap, String chainSpecie, Set<RefSeqGene> orthologs) {
		
		int res=0;
		for(RefSeqGene g: allIsoforms){
			String name=g.getName();
			String[] clnName =name.split("-");
			String specie = transcriptSpecieMap.get(clnName[0]);
			if ( specie!=null && specie.equalsIgnoreCase(chainSpecie) ) {
				res=1;
				orthologs.add(g);
			}
		}
		return res;
	}
	
	//Update the specie map with the genes that map to it
	private static void updateSpecieGeneMap(HashMap<String, Set<String>> speciesGeneMap,
			Collection<RefSeqGene> allIsoforms,HashMap<String, String> transcriptSpecieMap, String refgene, HashMap<String, String> foundGenes) {

		for(RefSeqGene g: allIsoforms){
			String name=g.getName();
			String[] clnName =name.split("-");
			String species = transcriptSpecieMap.get(clnName[0]);
			if ( species!=null  ) {
				if (! speciesGeneMap.containsKey(species))
					speciesGeneMap.put(species, new TreeSet<String>());
				speciesGeneMap.get(species).add(refgene);
				if (! foundGenes.containsKey(refgene))
					foundGenes.put(refgene,"");
			}
			else {
				System.err.println("Not found in species file "+ clnName +" mapped to "+refgene);
			}
		}
		
	}


	// Receive a net/chain alignment table (UCSC) and load it to a BEDfile, while maintaining info 
	// in the extra fields variables. Load only the regions that overlap transcripts in the reference bed
	private static BEDFileParser loadNetChainTab(String chainTabF, BEDFileParser refbed) throws IOException {
		BEDFileParser bed = new BEDFileParser();
		BufferedReader reader= new BufferedReader( new FileReader(chainTabF));
		String nextLine;
		nextLine = reader.readLine();
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			String[] a = nextLine.split("\t");
			String chr = a[2];
			int st= Integer.valueOf(a[3]);
			int end= Integer.valueOf(a[4]);
			String level = a[1];
			String type = a[15];
			RefSeqGene g= new RefSeqGene(chr,st,end);
			if (! refbed.getOverlappers(g).isEmpty()){
				g.addAttribute("netLevel", level);
				g.addAttribute("netType", type);
				bed.addRefSeq(g);
			}
		}
		reader.close();
		return bed;
	}

	private static HashMap<String, String> getMappedTranscriptIds(HashMap<String, Locus> geneVsMapped,HashMap<String, String> plainNames, boolean splitByUnderscr) {
		
		HashMap<String, String> res = new HashMap<String, String>(); 
		for (String g: geneVsMapped.keySet()){
			for (RefSeqGene iso: geneVsMapped.get(g).getAllIsoforms("")){
				String name=iso.getName();
				String newname=name;
				if (splitByUnderscr){
					String[] cln = name.split("-");
					newname=cln[0];
					String [] cln2 = name.split("\\.");
					plainNames.put(cln2[0],"");
				}
				res.put(newname,"");
			}
		}
		
		return res;
	}
	
	private static HashMap<String, String> loadLociSpecieMap(
			String transcriptSpecieMapF, HashMap<String, String> transmapTranscriptNames, boolean splitByUnderscr) throws IOException {

		BufferedReader reader= new BufferedReader( new FileReader(transcriptSpecieMapF));
		String nextLine;
		HashMap<String, String> res = new HashMap<String, String>();
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			String[] a = nextLine.split("\t");
			String sp=a[0];
			String t=a[1];
			if (splitByUnderscr)
			{
				String[] cln = t.split("-");
				t=cln[0];
			}
			if( transmapTranscriptNames.containsKey(t))	
				res.put(t,sp);
		}
		reader.close();
		return res ;
	}


	private static HashMap<String, String> loadCodingMap(String ucscCodingMap, HashMap<String, String> transmapTranscriptNames,HashMap<String, String>  plainTransmapNames) throws IOException {
		
		HashMap<String, String> res= new HashMap<String, String>();
		BufferedReader reader= new BufferedReader(new FileReader(ucscCodingMap));
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			String[] a=nextLine.split("\t");
			if (transmapTranscriptNames.containsKey(a[0]) || plainTransmapNames.containsKey(a[0]))
				res.put(a[0],a[1]);
		}
		reader.close();
		return res;
	}
	
	
	//Input:  a list of TFs binding peaks and a transcript bedfile, multilist states if we have 1 bed file chr/start/end/tf_name or multple chr/st/end files per each TF 
	//Output: 1. a table that marks a binding upstream (-1), downstream (1), or none (0) for each transcript and TF pair
	//		  2. a table with stats on each data set - how many peaks in total, how many crossed the transcript BED
	private static void PraseTFbind(String in, String outStr, String listF,Integer upstream, Integer downstream,String isoGeneMapF,boolean multiList) throws IOException {

		BEDFileParser bed= new BEDFileParser(in);
			
		BufferedWriter tab1= new BufferedWriter (new FileWriter(outStr+"_bindingMatrix"));
		BufferedWriter tab2= new BufferedWriter (new FileWriter(outStr+"_setStatistics"));
		BufferedWriter tab3= new BufferedWriter (new FileWriter(outStr+"_bindingMatrixByGene"));
		tab2.write("TF_set\tTF\ttotalPeakNum\toverlapPeakNum\n");
		
		HashMap<String, Locus> geneVsUpsPeaks = new HashMap<String,Locus>();
		HashMap<String, Locus> geneVsDownsPeaks = new HashMap<String,Locus>();
		for (RefSeqGene g: bed.GetGenes()){
			geneVsUpsPeaks.put(g.getName(),new Locus(new RefSeqGeneWithIsoforms(g))) ;
			geneVsDownsPeaks.put(g.getName(),new Locus(new RefSeqGeneWithIsoforms(g))) ;
		}
		ArrayList<String> TFsetsList= new ArrayList<String>();
		
		if (multiList )
			loadBindingInfoFromMultipleFiles (listF,TFsetsList,bed,geneVsUpsPeaks,geneVsDownsPeaks,upstream,downstream ,tab2);
		else
			loadBindingInfoFromSingleFiles (listF,TFsetsList,bed,geneVsUpsPeaks,geneVsDownsPeaks,upstream,downstream ,tab2);

			
		
		
		/*
		BufferedReader br = new BufferedReader (new FileReader(listF));
		String nextLine=br.readLine();
		while((nextLine = br.readLine()) != null && (nextLine.trim().length() > 0)){
			String[] a =nextLine.split("\t");
			String path=a[0];
			String TF_set=a[1];
			String TF=a[2];
			TFsetsList.add(TF_set);
			BEDFileParser peaks = parsePeakFile(path);
			int totalPeakNum= 0;
			int overlapPeakNum= 0;
			for (RefSeqGene bs: peaks.GetGenes()){
				totalPeakNum++;
				overlapPeakNum+=updateOverlappingPeaks (bed,bs,TF_set,geneVsUpsPeaks,upstream,0); 
				overlapPeakNum+=updateOverlappingPeaks (bed,bs,TF_set,geneVsDownsPeaks,0,downstream); 
			}
			tab2.write(TF_set+"\t"+TF+"\t"+totalPeakNum+"\t"+overlapPeakNum+"\n");	
		}
		*/
		
		
		tab1.write("TranscriptName");
		tab3.write("LociName");
		int totalTF=0;
		for (String tf: TFsetsList){
			tab1.write("\t"+tf);
			tab3.write("\t"+tf);
			 totalTF++;
		}
		tab1.write("\n");
		tab3.write("\n");
		
		HashMap<String , LinkedList<String>> geneIsoMap  = ArrAnnotUtils.loadIsoGeneMap(isoGeneMapF);
		HashMap<String , ArrayList<Integer>> isoBindingTab = new HashMap<String , ArrayList<Integer>>();
		
		for (RefSeqGene g: bed.GetGenes()){
			ArrayList<Integer> bindRow= new ArrayList<Integer>();
			tab1.write(g.getName());
			HashMap<String,Integer> hits = new HashMap<String,Integer>();
			updateTFHits(g,geneVsDownsPeaks,1,hits,2);
			updateTFHits(g,geneVsUpsPeaks,-1,hits,2);
			for (String tf: TFsetsList){
				int val =0;
				if (hits.containsKey(tf)){
					val=hits.get(tf);
				}
				tab1.write("\t"+val);
				bindRow.add(val);
			}
			tab1.write("\n");
			isoBindingTab.put(g.getName(), bindRow);
		}
		tab1.close();
		tab2.close();
		
		for ( String gene:geneIsoMap.keySet() ){
			int[] geneHits=new int[totalTF];
			tab3.write(gene);
			for(String iso: geneIsoMap.get(gene)){
				if (isoBindingTab.containsKey(iso)){
					ArrayList<Integer> hits = isoBindingTab.get(iso);
					for(int i=0;i<hits.size();i++){
						if (geneHits[i]==0)	
							geneHits[i]=hits.get(i);
						else{
							if (geneHits[i]*hits.get(i) < 0) //if they are from 2 diferent dire ctions mark by 2, else, leave as is
								geneHits[i]=2;
						}
					}
				}
			}
			for(int i=0;i<geneHits.length;i++)
				tab3.write("\t"+geneHits[i]);
			tab3.write("\n");
		}
		
		tab3.close();
		
		
	}
	

	
	

	private static void loadBindingInfoFromMultipleFiles(String listF,
			ArrayList<String> TFsetsList, BEDFileParser bed,
			HashMap<String, Locus> geneVsUpsPeaks,
			HashMap<String, Locus> geneVsDownsPeaks, Integer upstream,
			Integer downstream, BufferedWriter tab2 ) throws IOException {


		BufferedReader br = new BufferedReader (new FileReader(listF));
		String nextLine=br.readLine();
		while((nextLine = br.readLine()) != null && (nextLine.trim().length() > 0)){
			String[] a =nextLine.split("\t");
			String path=a[0];
			String TF_set=a[1];
			String TF=a[2];
			TFsetsList.add(TF_set);
			BEDFileParser peaks = parsePeakFile(path);
			int totalPeakNum= 0;
			int overlapPeakNum= 0;
			for (RefSeqGene bs: peaks.GetGenes()){
				totalPeakNum++;
				overlapPeakNum+=updateOverlappingPeaks (bed,bs,TF_set,geneVsUpsPeaks,upstream,0); 
				overlapPeakNum+=updateOverlappingPeaks (bed,bs,TF_set,geneVsDownsPeaks,0,downstream); 
			}
			tab2.write(TF_set+"\t"+TF+"\t"+totalPeakNum+"\t"+overlapPeakNum+"\n");	
		}
		
	}



//each line in a 4 col bed file with the name describing the TF
	private static void loadBindingInfoFromSingleFiles(String infile,
			ArrayList<String> TFsetsList, BEDFileParser bed,
			HashMap<String, Locus> geneVsUpsPeaks,
			HashMap<String, Locus> geneVsDownsPeaks, Integer upstream,
			Integer downstream, BufferedWriter tab2 ) throws IOException {

		HashMap<String ,Integer> totalPeakNum= new HashMap<String ,Integer>();
		HashMap<String ,Integer> overlapPeakNum= new HashMap<String ,Integer>();
		BufferedReader br = new BufferedReader (new FileReader(infile));
		HashMap<String ,String> TFS = new HashMap<String,String>();
		BEDFileParser peaks = new BEDFileParser();
		String nextLine;
		
		
		while((nextLine = br.readLine()) != null && (nextLine.trim().length() > 0)){
			String[] a =nextLine.split("\t");
			String chr=a[0];
			int start= Integer.valueOf(a[1]);
			int end= Integer.valueOf(a[2]);
			String TF=a[3];
			TFS.put(TF,"");
			 LinkedList<Integer> b1 = new LinkedList<Integer>();
			 LinkedList<Integer> b2 = new LinkedList<Integer>();
			 
			 b1.add(start);
			 b2.add(start);
			RefSeqGene bs = new RefSeqGene (chr,start,end,TF,"*",b1,b2);
			if (! totalPeakNum.containsKey(TF) ){
				totalPeakNum.put(TF, 0);
				overlapPeakNum.put(TF, 0);
			}
			int tmp =0;
			
			tmp +=updateOverlappingPeaks (bed,bs,TF,geneVsUpsPeaks,upstream,0); 
			tmp += updateOverlappingPeaks (bed,bs,TF,geneVsDownsPeaks,0,downstream); 
			
			totalPeakNum.put(TF, totalPeakNum.get(TF)+1);
			overlapPeakNum.put(TF, overlapPeakNum.get(TF)+tmp);
			
			
			
			
		}
	
		for (String t : totalPeakNum.keySet() ){
			tab2.write(t+"\t"+t+"\t"+totalPeakNum.get(t)+"\t"+overlapPeakNum.get(t)+"\n");	
			TFsetsList.add(t);
		}
	}






	private static BEDFileParser parsePeakFile(String path) throws IOException {
		
		BEDFileParser bed = new BEDFileParser();
		BufferedReader br = new BufferedReader (new FileReader(path));
		System.err.println(path);
		String nextLine;
		nextLine = br.readLine();
		while((nextLine = br.readLine()) != null && (nextLine.trim().length() > 0)){
			String a[] = nextLine.split("\t");
			LinkedList<Integer> exSt= new LinkedList<Integer>(); 
			LinkedList<Integer> exEnd= new LinkedList<Integer>(); 
			exSt.add(Integer.valueOf(a[3]));
			exEnd.add(Integer.valueOf(a[4]));
			RefSeqGene g= new RefSeqGene (a[2],Integer.valueOf(a[3]),Integer.valueOf(a[4]),"","*",exSt,exEnd);
			
			
			bed.addRefSeq(g);
			
		}
		br.close();
		return bed;
	}


	private static void updateTFHits(RefSeqGene gene,HashMap<String, Locus> geneVsPeaks, int updateVal,
			HashMap<String, Integer> hits, int overlapUpdateVal) {

		Locus locus = geneVsPeaks.get(gene.getName());
		HashMap<String,Integer> setSize=locus.getNumIsoformsPerSet();
		for(String set:setSize.keySet()){
			if (setSize.get(set)==0)
				continue;
			if (hits.containsKey(set)   ){
				if ((hits.get(set)!=updateVal))
					hits.put(set,overlapUpdateVal);
			}
			else
				hits.put(set,updateVal);
		}
		
	}



	private static int updateOverlappingPeaks(BEDFileParser expendedBed,
			RefSeqGene bs,String TF_set, HashMap<String, Locus> geneVsPeaks,int upstream, int downstream) {
	
		int res=0;
		IntervalTree<RefSeqGeneWithIsoforms> chrTree=expendedBed.getChrTree(bs.getChr());
		if (chrTree==null) return res;
		int addVal=Math.max(upstream,  downstream);
		Iterator<Node<RefSeqGeneWithIsoforms>> overlapIt = chrTree.overlappers(Math.max(bs.getStart(),2)-addVal, bs.getEnd()+addVal);
		//Iterator<RefSeqGeneWithIsoforms> overlapIt = expendedBed.getOverlappers(bs).valueIterator();
		if (overlapIt==null) return res;
		while(overlapIt.hasNext()){
			RefSeqGeneWithIsoforms multiIso=overlapIt.next().getValue();
			for (RefSeqGene gene: multiIso.getAllIsoforms()){
				Alignments promoter= gene.getPromoter(upstream, downstream);
				Alignments bsAlgn= new Alignments(bs.getChr(),bs.getStart(),bs.getEnd());
				if (bsAlgn.overlaps(promoter)){
					String g_name=gene.getName();
					Locus locus= geneVsPeaks.get(g_name);
					locus.addIsoform(new RefSeqGeneWithIsoforms(bs), TF_set);
					res=1;
				}
			}
		}
		
		return res;

	}



	private static void Level2Analysis(String inGTF, BufferedWriter bw) throws IOException {

	    BEDFileParser gtf = new BEDFileParser();
	    gtf.loadGTF(new File(inGTF), "");
	    int row=gtf.GetGenes().size();
	    bw.write("#1.2\n"+row+"\t9\n");
	    String str="Transcript"+"\t"+"Gene"+"\t"+"antisense"+"\t"+"TranscriptLength"+"\t"+"NumExons\tPseudo\tPfam\tCSF\tORFSize\tFracOfOrf\tSrOrfRank\n";
	    bw.write(str);
	    for (RefSeqGene g : gtf.GetGenes()){
		String  gtfAttr = g.getExtraFields()[0];
		double transcriptSize =g.getTranscriptLength();
		int exNum= g.getNumExons();
		String tname = GFF.getGTFAttr(gtfAttr, "transcript_id");
		String gname = GFF.getGTFAttr(gtfAttr,"gene_id");
		int pfam = (GFF.getGTFAttr(gtfAttr,"PFAM_Domains").replace("\"", "")).equals("None")? 0:1;
		int psuedo =( GFF.getGTFAttr(gtfAttr,"pseudogene").replace("\"", "")).equals("1")? 1:0;
		int antisense = gtfAttr.contains("antisense")? 1:0 ;
		double csfSt = Double.valueOf(GFF.getGTFAttr(gtfAttr,"csf_start").replace("\"", "") ).doubleValue();
		double csfEnd = Double.valueOf( GFF.getGTFAttr(gtfAttr,"csf_end").replace("\"", "") ).doubleValue();
		double csf = Double.valueOf( GFF.getGTFAttr(gtfAttr,"csf").replace("\"", "") ).doubleValue();
		double csflength= csfEnd-csfSt;
		double csfFrac= csflength/transcriptSize;
		double stFrac = csfSt/transcriptSize;
		bw.write(tname+"\t"+gname+"\t"+antisense+"\t"+transcriptSize+"\t"+exNum+"\t"+psuedo+"\t"+pfam+"\t"+csf+"\t"+csflength+"\t"+csfFrac+"\t"+stFrac+"\n");
		
	    }

	}
	
	
	//Given an FSA batch output alignemnt file (> \n >transcript1 \n aln1 \n >transcript2 \n aln2) 
	//provide stats on aln length, seq length, amount identity 
	private static void parseFSAAlignment(String aln, String pair,String outprefix,String  isoGeneMapF, boolean stockholmformat) throws IOException {
		
		BufferedWriter bwTcon= new BufferedWriter (new FileWriter(outprefix+"_byTranscript.tab"));
		BufferedWriter bwXloc= new BufferedWriter (new FileWriter(outprefix+"_byGene.tab"));
		BufferedWriter bwXlocBias= new BufferedWriter (new FileWriter(outprefix+"_alnBias_byGene.gct"));
		
		BufferedReader br = new BufferedReader (new FileReader(aln));
		String nextLine;
		HashMap < String, ArrayList<String> > refVsStatStr = new  HashMap < String, ArrayList<String> >();
		HashMap < String, ArrayList<Double> > refVsStatVal = new  HashMap < String, ArrayList<Double> >();
		
		HashMap < String, HashMap <String,ArrayList<Double>> > refVsBiasStats = new  HashMap < String, HashMap <String,ArrayList<Double>> >();
		
		
		int i=0;
		double accuracy =0;
		String tcon1,seq1,tcon2 ,seq2,tmp;
		String arr[];
		bwTcon.write("tcon1\ttcon2\tlen1\tlen2\talignmentLen\tnumIdentity\tpctIdentityFromAlign\tpctIdentityFromRef1\tpctAlignementFromRef1\tGF_Acuuracy\n");
		while((nextLine = br.readLine()) != null && (nextLine.trim().length() > 0)){
			if (stockholmformat){
				br.readLine();
				tmp=br.readLine();
				//System.err.println(tmp);
				tmp=tmp.replace("#=GF", "");
				tmp=tmp.replace("Accuracy", "");
				tmp=tmp.replace(" ", "");
				//arr=tmp.split("\s");
				accuracy = Double.valueOf(tmp);
				tmp=br.readLine();
				arr=tmp.split("\\s+");
				tcon1=arr[0];
				seq1=arr[1];
				tmp=br.readLine();
				arr=tmp.split("\\s+");
				tcon2=arr[0];
				seq2=arr[1];
				tmp=br.readLine();
				tmp=br.readLine();
				i+=7;
				//System.err.println(accuracy + "\n" + tcon1 + "\n"+ seq1+ "\n" + tcon2 + "\n"+ seq2+ "\n" );
			}
			else{
				 tcon1=br.readLine();
				 seq1=br.readLine();
				 tcon2=br.readLine();
				 seq2=br.readLine();
				 i+=5;
			}
			tcon1=tcon1.replace(">", "");
			tcon2=tcon2.replace(">", "");
			
			if (seq1==null || seq2==null || tcon1==null || tcon2==null){
				System.err.println("error in line " +i);
			}
			else{
				ArrayList<Double> stat= computeAlnStat(seq1,seq2);
				ArrayList<Double> statBias = computeAlnBiasStat(seq1,seq2);
				stat.add(accuracy);
				String str="";
				str += tcon1+"\t"+tcon2;
				for (Double d: stat)
					str += "\t"+d;
				bwTcon.write(str+"\n");
				if (! refVsStatStr.containsKey(tcon1)){
					refVsStatStr.put(tcon1, new ArrayList<String>() );
					refVsStatVal.put(tcon1, new ArrayList<Double>() );
					refVsBiasStats.put(tcon1, new HashMap< String,ArrayList<Double>>());
				}
				refVsStatStr.get(tcon1).add(str);
				refVsStatVal.get(tcon1).add(stat.get(6));
				refVsBiasStats.get(tcon1).put(tcon2,statBias);
						
			}
			
		}
		bwXloc.write("xloc\ttcon1\ttcon2\tlen1\tlen2\talignmentLen\tnumIdentity\tpctIdentityFromAlign\tpctIdentityFromRef1\tpctAlignementFromRef1\tGF_Acuuracy\n");
		HashMap<String, LinkedList<String>> geneIsoMap = loadIsoGeneMap(isoGeneMapF);
		for (String gene: geneIsoMap.keySet()){
			List<String> isos= geneIsoMap.get(gene);
			String bestStr="";
			double bestScr=0;
			for (String iso : isos ){
				if (refVsStatVal.containsKey(iso)){
					ArrayList<Double> arr1 = refVsStatVal.get(iso);
					for (int j=0; j<arr1.size(); j++)
					{
						if (arr1.get(j) > bestScr){
							bestScr=arr1.get(j);
							bestStr=refVsStatStr.get(iso).get(j);
						}
					}
				}
			}
			if (! bestStr.equals("")){
				bwXloc.write(gene+"\t"+bestStr+"\n");
				String[] tmp1= bestStr.split("\t");
				ArrayList<Double> statBias=refVsBiasStats.get(tmp1[0]).get(tmp1[1]);
				bwXlocBias.write(gene+"\t"+tmp1[1]);
				
				for (Double d: statBias)
					bwXlocBias.write("\t"+d);
				bwXlocBias.write("\n");
			}
		}
		bwXloc.close();
		bwTcon.close();
		bwXlocBias.close();
	}



	


// out : len1\tlen2\talignmentLen\tnumIdentity\tpctIdentityFromAlign\tpctIdentityFromRef1\tpctAlignementFromRef1\n");


	private static ArrayList<Double> computeAlnStat(String seq1, String seq2) {
		
		double len1=0.0;
		double len2=0.0;
		double identity=0.0;
		double aln=0.0;
		char a;
		char b;
		char[] s1=seq1.toCharArray();
		char[] s2=seq2.toCharArray();
		for (int i=0 ; i< s1.length; i++){
			a=s1[i];
			b=s2[i];
			if (a != '-')
				len1++;
			if (b != '-')
				len2++;
			if (a != '-' && b !='-')
				aln++;
			if (a==b)
				identity++;
			
		}
		
		ArrayList<Double> res = new ArrayList<Double>();
		res.add(len1);
		res.add(len2);
		res.add(aln);
		res.add(identity);
		res.add(identity/aln);
		res.add(identity/len1);
		res.add(aln/len1);
		return res;
	}
	
	
	//runs over the refernce transcript calc accumulation across bins
	private static ArrayList<Double>  computeAlnBiasStat(String seq1, String seq2){
		
		ArrayList<Double> statBias = new ArrayList<Double>();
		
		
		double identity=0.0;
		char a;
		char b;
		char[] s1=seq1.toCharArray();
		char[] s2=seq2.toCharArray();
		String clnSeq= seq1.replace("-", "");
		int size1 = clnSeq.length();
		int binSize = size1/10;
		double cnt=0.0;
		for (int i=0 ; i< s1.length; i++){
			a=s1[i];
			b=s2[i];
			if (a != '-') {
				cnt++;
				if (a==b)
					identity++;
				if (cnt == binSize){
					statBias.add( identity / cnt);
					cnt=0;
					identity=0;
				}
			}
		}
		return statBias;
	}
	
	private static BEDFileParser selectRandomNonAnnotatedBlocks(Integer setSize,
			Integer blockSize, String fileList, String outprefix) throws IOException {
		
		//read bed files to one bed
		BufferedReader br = new BufferedReader (new FileReader(fileList));
		String nextLine;
		BEDFileParser annot = new BEDFileParser();
		while((nextLine = br.readLine()) != null ){
			 String [] arr=nextLine.split("\t");
			 System.err.println(arr[0]);
			 String file = arr[0];
			 BEDFileParser bed = new BEDFileParser(file);
			 annot.addRefSeqSet(bed.GetGenes());
		}
		br.close();
		System.err.println("size of annot file " + annot.GetGenes().size());
		BEDFileParser res=new BEDFileParser();
		//make random annots
		Iterator<String> chrIt = annot.getChromosomeIterator();
		int i=0;
		while(chrIt.hasNext()){
			String chr = chrIt.next();
			Iterator<RefSeqGeneWithIsoforms> geneIt = annot.getChrTree(chr).valueIterator();
			if (geneIt.hasNext()){
				RefSeqGeneWithIsoforms prev = geneIt.next();
				while (geneIt.hasNext()) {
					RefSeqGeneWithIsoforms curr = geneIt.next();
					int space = curr.getStart()-prev.getEnd();
					if (space - blockSize > blockSize){
						int randAdd= (int)Math.floor(Math.random()*(space - blockSize-1)); 
						Alignments nBlock= new Alignments (prev.getExons()[0]);
						nBlock.setStart(prev.getEnd()+randAdd);
						nBlock.setEnd(prev.getEnd()+randAdd+blockSize-1);
						nBlock.setName("random_"+i);
						i++;
						res.addRefSeq(new RefSeqGene(nBlock));
						prev=curr;
					}
				}
			}
			System.err.println("size of res file " + res.GetGenes().size() + " " + chr);
		}
		res.writeFullBed(outprefix+".bed");
		return res;
	}

	
	

	private static void MakeGWASTable(String metaBedf, String expMatf,
			String gwasTabf, String detectionMatf, String outprefix,String buildMapBedf
			,String filterf,String EnhancerBedF,String snpPos,String snpAssociation) throws IOException, ParseException {
		
		//store meta bed line by gene ID
		BEDFileParser enhancerBed =new BEDFileParser();
		if (!EnhancerBedF.equals(""))
			enhancerBed= new BEDFileParser(EnhancerBedF);
		HashMap<String , String[]> metaMap= new HashMap<String , String[]>(); 
		BufferedReader br = new BufferedReader (new FileReader(metaBedf));
		String nextLine ;
		String MetaHeaderLine = br.readLine();
		MetaHeaderLine.trim();
		String[] MetaHeader = MetaHeaderLine.split("\t");
		int conservativeIx=0;
		int enhancerIx=0;
		for (int i=0; i<MetaHeader.length;i++){
			if (MetaHeader[i].equalsIgnoreCase("HighConfidenceSet")) conservativeIx=i;
			if (MetaHeader[i].contains("Enhancer")) enhancerIx =i;
		}
		while((nextLine = br.readLine()) != null && (nextLine.trim().length() > 0) ){
				 String [] arr=nextLine.split("\t");
				 metaMap.put(arr[3],arr);
		}
		br.close();
		
		//GCT expMat by gene id
		MatrixWithHeaders expmat=  new MatrixWithHeaders(expMatf);
		//GCT detection mat
		MatrixWithHeaders detectmat=  new MatrixWithHeaders(detectionMatf);
		
		//gwas Tab by bed and Map (rs_pubmed) ; list of all associated rs_pubmed, and all LD_trait
		Map <String, String[]>  rsPubmed_data_map =new HashMap <String, String[]>();
		Map <String, ArrayList<String>>  ldTrait_rsPubmed_map =new HashMap <String, ArrayList<String>>();
		Map <String, String[]>  ldTrait_data_map =new HashMap <String, String[]>();
		BEDFileParser LD_bed= new BEDFileParser();
		br = new BufferedReader (new FileReader(gwasTabf));
		
		String  gwasHeaderLine = br.readLine();
		gwasHeaderLine.trim();
		String[] gwasHeader = gwasHeaderLine.split("\t");
		
		//Process SNP data 
		HashMap<String,LinkedList<String>> snp_ldSnps = new HashMap<String,LinkedList<String>>();
		BEDFileParser snpLoc = new  BEDFileParser ();
		LoadLDSnpData(snp_ldSnps,snpLoc,snpPos,snpAssociation);
		Map <String, RefSeqGene>snpNameLocMap = snpLoc.getNameGeneMap();
		
		int snpIx=0;
		int traitIx=0;
		int pubmedIx=3;
		int pIx=0;
		int posSnpIx=0;
		int nogeneI=0;
		int reportedGenI=0;
		for (int i=0; i<gwasHeader.length; i++){
			if (gwasHeader[i].equalsIgnoreCase("SNPs"))
				snpIx=i;
			if (gwasHeader[i].equalsIgnoreCase("disease/trait"))
				traitIx=i;
			if (gwasHeader[i].equalsIgnoreCase("p-Value"))
				pIx=i;
			if (gwasHeader[i].equalsIgnoreCase("genes"))
				nogeneI=i;
			if (gwasHeader[i].equalsIgnoreCase("posSt"))
				posSnpIx=i;
			if (gwasHeader[i].equalsIgnoreCase("Reported Gene(s)"))
				reportedGenI=i;
				
		}
		while((nextLine = br.readLine()) != null && (nextLine.trim().length() > 0)){
				 String [] arr=nextLine.split("\t");
				 if (arr[pIx].matches("NS") || arr[pIx].matches("E"))
					arr[pIx] ="1";
				 rsPubmed_data_map.put(arr[3],arr);
				 String trait=arr[traitIx];
				 trait=trait.replace(" ","_");
				 String ldTrait=trait+"_"+arr[0]+":"+arr[1]+"-"+arr[2];
				 RefSeqGene LD_region = new RefSeqGene(new Alignments(arr[0],arr[1],arr[2]));
				 LD_region.setOrientation("*");
				 IntervalTree<RefSeqGeneWithIsoforms> overlaps= LD_bed.getOverlappers(LD_region);
				 boolean update=false;
				 if (!overlaps.isEmpty()) {
					 Iterator <RefSeqGeneWithIsoforms> it= overlaps.valueIterator();
				 	 while(it.hasNext()){
				 		RefSeqGeneWithIsoforms g = it.next();
				 		if (g.getStart()==(new Double(arr[1])) & g.getEnd()==(new Double(arr[2]))){
				 			g.addAttribute("trait",g.getAttribute("trait")+";"+ ldTrait);
				 			update=true;
				 		}
				 	 }
				 }
				 if (update == false) { 
					 LD_region.addAttribute("trait", ldTrait );
					 LD_bed.addRefSeq(LD_region);
				}
				 
				 if  (! ldTrait_rsPubmed_map.containsKey(ldTrait)){
					 ldTrait_rsPubmed_map.put(ldTrait, new ArrayList<String>());
					 ldTrait_data_map.put(ldTrait, arr);
					 
				 }
				 ldTrait_rsPubmed_map.get(ldTrait).add(arr[3]);
				 
				 if (new Double (ldTrait_data_map.get(ldTrait)[pIx]) < new Double(arr[pIx]) )
					 ldTrait_data_map.put(ldTrait,arr);
				 
		}
		br.close();
		
		//build map BED 
		BEDFileParser mapbed = new BEDFileParser(buildMapBedf); 
		Map <String,RefSeqGene> name_gene_map= mapbed.getNameGeneMap();
		
		//Count total overlaps (get pval later)
		double totLDtraitHits=0.0;
		double noGeneHit=0.0;
		double intergenicHit=0;
		HashMap <String , int[]> unqLDcnt = new HashMap <String , int[]>();
		double noGenLD=0.0;
		double intergenicLD=0.0;
		double conservativeHit=0.0;
		int ldEnhCnt=0;
		Set<RefSeqGene> uniqXlocSet = new HashSet<RefSeqGene>();
		//write out put : line for every (xloc : rs_pubmed) or (xloc : LD trait)
		//for (xloc: LD trait ) - LD,trait,rs_pumed list, best rs_pebmed, best pval, xloc , xloc exp 
		
		BufferedWriter bw_all = new BufferedWriter(new FileWriter(outprefix+"_allSnpsTranscripts.tab"));
		BufferedWriter bw_LDtrait = new BufferedWriter(new FileWriter(outprefix+"_LDRTranscripts.tab"));
		BufferedWriter bw_LDtrait_Exp = new BufferedWriter(new FileWriter(outprefix+"LDR_exp.tab"));
		BufferedWriter bw_LDtrait_dExp = new BufferedWriter(new FileWriter(outprefix+"LDR_detectMat.tab"));
		
		bw_all.write(gwasHeaderLine+"\tLD_enh\tXloc\tintron-exon\tOthorOverlappingSnps\tNoGene\t"+MetaHeaderLine+"\n");
		bw_LDtrait.write("LD_trait\t"+gwasHeaderLine+"\tLD_enh\tXloc\tintron-exon\tOthorOverlappingSnps\tNoGene\t"+MetaHeaderLine+"\n");
		bw_LDtrait_Exp.write("Xloc\tLD_trait\tNoGene\tLD_enh");
		for (int i=0; i<expmat.columnDimension();i++){
			bw_LDtrait_Exp.write("\t"+expmat.getColoumnName(i));
		}
		bw_LDtrait_Exp.write("\n");
		bw_LDtrait_dExp.write("Xloc\tLD_trait\tNoGene\tLD_enh");
		for (int i=0; i<detectmat.columnDimension();i++){
			bw_LDtrait_dExp.write("\t"+detectmat.getColoumnName(i));
		}
		bw_LDtrait_dExp.write("\n");
		
		String[] emptyS= new String[1]; emptyS[0]="";
		double[] emptyD= new double[0]; 
		
		//cross build map with gwas bed
		for (RefSeqGene ldr : LD_bed.GetGenes()){
			
			String ldr_pos = ldr.getChr()+":"+ldr.getStart()+"-"+ldr.getEnd();
			IntervalTree<RefSeqGeneWithIsoforms> xlocOverlap = mapbed.getOverlappers(ldr);
				if (!xlocOverlap.isEmpty()){
				Iterator <RefSeqGeneWithIsoforms> xlocOverlapIt=  xlocOverlap.valueIterator();
				String [] traitsArr = ldr.getAttribute("trait").split(";");
				
				for ( int i=0; i<traitsArr.length; i++){
					String ldTrait=traitsArr[i];
					//Get LD trait region infor
					 String[]ld_trait_info =  ldTrait_data_map.get(ldTrait);
					 boolean noGen= ld_trait_info.length > nogeneI? ld_trait_info[nogeneI].equalsIgnoreCase("nogene"):false;
					 boolean intergenic =  ld_trait_info.length > reportedGenI? (ld_trait_info[reportedGenI].equalsIgnoreCase("intergenic")|| ld_trait_info[reportedGenI].equalsIgnoreCase("NR")):false;
					//Ld region overlap with enh
					 int ld_enh = (enhancerBed.getOverlappers(ldr).isEmpty())? 0:1;
					 //count LD regions and the number of genes in each region , and howmany of which have no gene
					 if (! unqLDcnt.containsKey(ldr_pos)){
						 int[] a= new int[2];
						 a[0]=xlocOverlap.size();
						 a[1]= noGen? 1:0;
						 unqLDcnt.put(ldr_pos,a);
						 if (noGen){
							 noGenLD++;
							 if (intergenic)
								 intergenicLD++;
							 if (ld_enh == 1)
								 ldEnhCnt++;
						 }
					 }
					 while (xlocOverlapIt.hasNext()){
						 totLDtraitHits++;
						 if (noGen)
							 noGeneHit++;
						 if (noGen && intergenic)
							 intergenicHit++;
						 
						 for (RefSeqGene g : xlocOverlapIt.next().getAllIsoforms()){
							 String [] xloc_info= metaMap.containsKey(g.getName())? metaMap.get(g.getName()): emptyS;
							double[] xloc_exp = expmat.containsRow(g.getName())? expmat.getRow(g.getName()): emptyD;
							double[] xloc_detect = detectmat.containsRow(g.getName())? detectmat.getRow(g.getName()): emptyD;
							if (noGen) uniqXlocSet.add(g);
							if (xloc_info.length > conservativeIx && xloc_info[conservativeIx].equals("1") && noGen )
								 conservativeHit++;
							int intExOverlap = 0;
							int overlapOtherSnps=0;
							
							
							for (String rsPubmed: ldTrait_rsPubmed_map.get(ldTrait))
							 {
								
								String[] rs_info= rsPubmed_data_map.get(rsPubmed);
								int curr =0 ;
								int curr2=0;
								//overlap an intron or exon
								if (g.overlapsGene(new RefSeqGene(new Alignments(g.getChr(),new Integer(rs_info[posSnpIx]),new Integer(rs_info[posSnpIx])+1)) ))
									curr =1;
								//overlap exon
								if (g.overlaps(new RefSeqGene(new Alignments(g.getChr(),new Integer(rs_info[posSnpIx]),new Integer(rs_info[posSnpIx])+1)) ))
									curr =2;
								if (curr> intExOverlap)
									intExOverlap=curr;
								//Does it overlap other SNP in the region
								if (snp_ldSnps.containsKey(rs_info[snpIx])){
									
									LinkedList<String> linkedSnps= snp_ldSnps.get(rs_info[snpIx]);
									for (String snp: linkedSnps){
										System.err.println(snp+"\t"+"hasNeighborSnps");
										if (snpNameLocMap.containsKey(snp) && snpNameLocMap.get(snp).overlapsGene(g) ){
											curr2=1;
											if (snpNameLocMap.get(snp).overlaps(g) )
												curr2=2;
										}
										if (curr2 > overlapOtherSnps ){
											overlapOtherSnps=curr2;
											System.err.println(g.getName()+"\t"+snp+"\t"+curr2);
										}
											
										}
								}
								
								bw_all.write(toTabSepString(rs_info)+ 
										"\t"+ld_enh+"\t"+g.getName()+"\t"+curr+"\t"+curr2+"\t"+noGen+"\t"+toTabSepString(xloc_info)+"\n");;
							 
							 }
							
							bw_LDtrait.write(ldTrait+"\t"+toTabSepString(ld_trait_info)+ 
									"\t"+ld_enh+"\t"+g.getName()+"\t"+intExOverlap+"\t"+overlapOtherSnps+"\t"+noGen+"\t"+toTabSepString(xloc_info)+"\n");
							String exp=toTabSepString(xloc_exp);
							if (!exp.isEmpty())
								bw_LDtrait_Exp.write(g.getName()+"\t"+ldTrait+"\t"+noGen+"\t"+ld_enh+"\t"+exp+"\n");
							
							
							String dexp=toTabSepString(xloc_detect);
							if (!dexp.isEmpty())
								bw_LDtrait_dExp.write(g.getName()+"\t"+ldTrait+"\t"+noGen+"\t"+ld_enh+"\t"+dexp+"\n");
							
							
						 }
						 
					 }
					
					
						
				}
			}
			
		}
			
		 bw_all.close() ;
		 bw_LDtrait.close() ;
		 bw_LDtrait_Exp.close();
		 bw_LDtrait_dExp.close();
		
		 int consGen=0;
		 int enhGen=0;
		 for (RefSeqGene g :uniqXlocSet){
			 String [] xloc_info= metaMap.containsKey(g.getName())? metaMap.get(g.getName()): emptyS;
			 if (xloc_info.length>enhancerIx){
				 if (xloc_info[conservativeIx].equals("1")) consGen++;
				 if (xloc_info[enhancerIx].equals("1")) enhGen++;
			 }
		 }
		 System.err.println("Total Number of unique LD regions: " + unqLDcnt.size() + " of which "+ noGenLD + " have no gene in LD");
		 System.err.println("Total Number of LDTrait with no gen: " + noGenLD + " of which "+ intergenicLD + " are also reported as intergenic , and " +ldEnhCnt + " overlap Enhancers");
		 System.err.println("Total Number of LDTrait_Transcript  pairs: " + totLDtraitHits + " of which "+ noGeneHit + " have no gene in LD (" + intergenicHit +"; reported as intergenic ; "+ conservativeHit+" are conservative lincs) " );
		 if (conservativeIx != 0 &  enhancerIx != 0)
			 System.err.println("In total there are  " +uniqXlocSet.size() + " unq transcripts associated with no genes of which  " + consGen +" are conservative set and " + enhGen+ " overlap enhancers  ");
		//TODO: Calc enrichment according to a permutation test
		 if (! filterf.equals("")){
			 
		 }
	}





	private static void LoadLDSnpData(
			HashMap<String, LinkedList<String>> snpLdSnps, BEDFileParser snpLoc, String snpPos, String snpAssociationF) throws IOException {
		
		if (snpPos.equals("") || snpAssociationF.equals(""))
			return;
		BufferedReader br = new BufferedReader (new FileReader(snpAssociationF));
		String nextLine;
		HashMap<String,String> allSnp= new HashMap<String,String>();
		while((nextLine = br.readLine()) != null && (nextLine.trim().length() > 0)){
		    String [] arr= nextLine.split("\t");
		    LinkedList<String> l = new LinkedList<String>();
		    for (int i=0; i<arr.length; i++){ l.add(arr[i]); allSnp.put(arr[i],"");}
		    for (int i=0; i<arr.length; i++) {snpLdSnps.put(arr[i],l);  allSnp.put(arr[i],"");}
		   
		    
		}
		br.close();
		
		br = new BufferedReader (new FileReader(snpPos));
		while((nextLine = br.readLine()) != null && (nextLine.trim().length() > 0)){
		    String [] arr= nextLine.split("\t");
		    if( allSnp.containsKey(arr[3])){
		    	Alignments a = new Alignments (arr[0] ,new Integer(arr[1]), new Integer( arr[1])+1);
		    	a.setName(arr[3]);
		    	snpLoc.addRefSeq(new RefSeqGene(a));
		    }
		
		}
		br.close();
		
		
	}




	private static String toTabSepString(String[] vec) {
		String res="";
		
		for (int i=0; i<vec.length-1;i++)
			res+=vec[i]+"\t";
		res+=vec[vec.length-1];
		return res;
	}

	private static String toTabSepString(double[] vec) {
		String res="";
		if (vec.length==0)
			return res;
		for (int i=0; i<vec.length-1;i++)
			res+=vec[i]+"\t";
		res+=vec[vec.length-1];
		return res;
	}
	
	
//Reports the trascripts ids that remain or not
	private static void level1Level2BugFix(String l1gtf, String l2gtf,String refseqSrcStr,String l1OutF,String l2OutF) throws IOException {

		
		
		GTFFileParser l1= new GTFFileParser(l1gtf);
		GTFFileParser l2= new GTFFileParser(l2gtf);
		BEDFileParser l2merge= l2.mergeByGtfGeneId();
		
		HashMap<String,ArrayList<String>> l1XlocTcon = l1.getGTFGeneTranscriptMap () ; 
		HashMap<String,ArrayList<String>> l2XlocTcon = l2.getGTFGeneTranscriptMap ();
		
		HashMap<String,ArrayList<String>> l1Res = l1.getGTFGeneTranscriptMap ();
		HashMap<String,ArrayList<String>> l2Res = l2.getGTFGeneTranscriptMap ();
		
		
		//If a gene was already removed from level2 or decided to be retain in level 1,  remove the other
		HashMap<String,String> l2Removed = new HashMap<String,String>();
		HashMap<String,String> l1Stayed = new HashMap<String,String>();
		int l2removeCnt=0;
		
		for (RefSeqGene g : l2merge.GetGenes() ){
			String gid=g.getAttribute("gene_id");
			IntervalTree <RefSeqGeneWithIsoforms> l1_overlap = l1.getOverlappers(g);
			if (! l1_overlap.isEmpty()){
				Iterator<RefSeqGeneWithIsoforms>l1it = l1_overlap.valueIterator();
				while (l1it.hasNext()){
					RefSeqGeneWithIsoforms curr =l1it.next();
					for (RefSeqGene iso: curr.getAllIsoforms()){
						if (iso.getOrientation().equalsIgnoreCase(g.getOrientation())){
							String gid2=iso.getAttribute("gene_id");
							String src=iso.getAttribute("source");
							if (src.equalsIgnoreCase(refseqSrcStr)){
								if(l2Res.containsKey(gid)) { 
									l2Res.remove(gid);
									l2removeCnt++;
									l2Removed.put(gid,"");
								}
							}
						}
					}
				}
			}
		}
		
		//Go over level2, with anything that crosses: check if level 1 has refseq
			// if does cross out level 2
			// else cross out level 1
		
		for (RefSeqGene g : l2merge.GetGenes() ){
			
			String gid=g.getAttribute("gene_id");
	
			if (l2Removed.containsKey(gid))
				continue;
			IntervalTree <RefSeqGeneWithIsoforms> l1_overlap = l1.getOverlappers(g);
			if (! l1_overlap.isEmpty()){
				Set<String> inRef = new HashSet<String>();
				Set<String> notRef = new HashSet<String>();
				Iterator<RefSeqGeneWithIsoforms>l1it = l1_overlap.valueIterator();
				while (l1it.hasNext()){
					RefSeqGeneWithIsoforms curr =l1it.next();
					for (RefSeqGene iso: curr.getAllIsoforms()){
						if (iso.getOrientation().equalsIgnoreCase(g.getOrientation())){
							String gid2=iso.getAttribute("gene_id");
							notRef.add(gid2);
						}
						
					}
				}
				//If in ref is empty remove all annotaion from level 1 
				if (! notRef.isEmpty()){
					for (String gid3:notRef){
						if (l1Res.containsKey(gid3)) l1Res.remove(gid3);
					}
				}
			
			}
		}
		
		//Write output xloc and tcon ids for each
		writeXlocTconIds ( l1Res,l1OutF);
		writeXlocTconIds ( l2Res,l2OutF);
		
		System.err.println("l2removeCnt : "+l2removeCnt);

	}




private static void writeXlocTconIds(HashMap<String, ArrayList<String>> gene_tcon_map,	String outf) throws IOException {

	BufferedWriter xlocs = new BufferedWriter(new FileWriter (outf+"gene.ids"));
	BufferedWriter tcons = new BufferedWriter(new FileWriter (outf+"transcript.ids"));
	for (String g:gene_tcon_map.keySet()){
		xlocs.write(g+"\n");
		for (String t: gene_tcon_map.get(g)) 
			tcons.write(t+"\n");
	}
	xlocs.close();
	tcons.close();
}
	
	


	//for every gene ID : report one file with the max exon transcript and afile with the max transcript
	private static void selectMostComplexTranscriptPerGeneI(String gtf,	String out) throws IOException {
		
		BEDFileParser inbed= new BEDFileParser();
		inbed.loadGTF(new File(gtf),"");
		
		HashMap <String,Locus> geneMap = new HashMap <String,Locus>();
		BEDFileParser outBed_1=new BEDFileParser();
		BEDFileParser outBed_2=new BEDFileParser();
	
		for (RefSeqGene iso: inbed.GetGenes()){
			String gene_id="";
			if (iso.getExtraFields() != null)
				gene_id=GFF.getGTFAttr(iso.getExtraFields()[0], "gene_id");
			else
				gene_id=iso.getAttribute("gene_id");
			//System.err.println(iso.getName()+"\n");
			if (!geneMap.containsKey(gene_id))
				geneMap.put(gene_id,new Locus(new RefSeqGeneWithIsoforms(iso)) );
			geneMap.get(gene_id).addIsoform(new RefSeqGeneWithIsoforms(iso), "def");
		}
		
		for (String geneName: geneMap.keySet()){
			
			Locus loc = geneMap.get(geneName);
			RefSeqGene g= loc.getLongestIntornChainIso("def");
			g.setName(geneName);
			outBed_1.addRefSeq(g);
			
			
		}
		outBed_1.writeFullBed(out+"_maxIntronChainIso.bed");
		
	}

	
	
	private static void overlapEnrichment(String inSet, String element_bed,
			Integer permNum, String centromers, String annotFilter,	String outPrefix,String chrSizes,String informat) throws IOException {

		BEDFileParser set= new BEDFileParser();
		BEDFileParser geneMerge_set;
		if (informat.equalsIgnoreCase("gtf")){
			set.loadGTF(new File(inSet),"");
			geneMerge_set = set.mergeByGtfGeneId();
		}
		else
			geneMerge_set = new BEDFileParser(inSet);
		
		BEDFileParser merged_elements= new BEDFileParser(element_bed);  
		merged_elements.merge();
		
		int cnt=0;
		int cnt2=0;
		int OverlapSize=countExonLevelOverlap(geneMerge_set,merged_elements);
		
		BEDFileParser overlap2 = geneMerge_set.overlap_GenomicLevel(merged_elements,"stay","false");
		int overlapGenomicSize = overlap2.GetGenes().size();
		
		if (permNum>0){
			String [] val1 = new String[permNum];
			String [] val2 = new String[permNum];
			
			int i=0;
			transcriptsNullModel nullmodel = new transcriptsNullModel (new BEDFileParser (annotFilter),new BEDFileParser(centromers),chrSizes); 
	    	nullmodel.makeRandTranscriptList(permNum, geneMerge_set);
	    	
	    	for (BEDFileParser b : nullmodel.getRandomSets()){
	    		int cntRand=countExonLevelOverlap(b,merged_elements);
	    		if (cntRand >= OverlapSize & cntRand >0)
	    			cnt++;
	    		val1[i]=String.valueOf(cntRand);
	    		
	    		
	    		
	    		BEDFileParser overlapRand = b.overlap_GenomicLevel(merged_elements,"stay","false");
	    		int overlapGenomicSizeRand = overlapRand.GetGenes().size();
	    		if (overlapGenomicSizeRand >= overlapGenomicSize &  overlapGenomicSizeRand > 0)
	    			cnt2++;
	    		val2 [i] = String.valueOf(overlapGenomicSizeRand);
	    		
	    		i++;	
	    	}
	    	BufferedWriter bw = new BufferedWriter (new FileWriter(outPrefix+"randomOverlapCnts.txt"));
	    	writeTabLine(val1,bw);
	    	bw.write("\n");
	    	writeTabLine(val2,bw);
	    	bw.write("\n");
	    	bw.close();
		}
		
		int sizeGene = geneMerge_set.GetGenes().size();
		double pct = ((double)OverlapSize)/((double) sizeGene);
		double pct2 = ((double)overlapGenomicSize)/((double) sizeGene);
		System.err.println(" Number of elements: "+ merged_elements.GetGenes().size()+   "Number of unique genes:  " + sizeGene + " out of which " + OverlapSize +"("+pct+") overlap the element at the exon level"   );
		System.err.println(" permPval :" + ((double)cnt)/((double) permNum) + " " + cnt);
		System.err.println( overlapGenomicSize + " overlaped in the genomic level ;" + pct2 + " permPval :" + ((double)cnt2)/((double) permNum) + " " + cnt2 );
		
		
		
	}




	private static int countExonLevelOverlap(BEDFileParser geneMergeSet,
			BEDFileParser elements) {
		
		int cnt=0;
	
		//count the number of genes that overlap the element across their exons
		for (RefSeqGene g: geneMergeSet.GetGenes()){
			
			IntervalTree<RefSeqGeneWithIsoforms> overlap = elements.getOverlappers(g);
			if (! overlap.isEmpty())
				cnt++;
			
			
		}
		
				
		
		return cnt;
	}
	
	
	
	private static void classifyEnhStates(String enhRegions, String hmmList,
			String outPrefix, String enhStateFile, String nullstate, int segmentSize) throws IOException {


		BEDFileParser refRegions= new BEDFileParser(enhRegions);
		Map<String,String> setsPath =  loadSetsPath(hmmList);
		Map<String, Locus> refLociMap = new  HashMap<String, Locus>  ();
		List<String> refStatesLst= loadNameLst(enhStateFile,false); 
		BEDFileParser refSegments = segmentRegions(refRegions,segmentSize);
		
		for(RefSeqGene g:refSegments.GetGenes()){
			refLociMap.put(g.getName(),new Locus(new RefSeqGeneWithIsoforms(g)));
			
		}
		
		LinkedList<String> cellNames= loadCompatibleToRefLociIsoformMap(refSegments,
				 setsPath,refLociMap,false,false) ;
		
		BEDFileParser resBed = new BEDFileParser();
		for (RefSeqGene segmentGen : refSegments.GetGenes()){
			
			String segment = segmentGen.getName();
			Locus loc = refLociMap.get(segment);
			Collection<RefSeqGene> isoSet=loc.getAllIsoforms("");
			double enhCnt=0.0;
			double nonCnt=0.0;
			double elseCnt=0.0;
			for (RefSeqGene iso : isoSet){
				String state = iso.getName();
				if (state.equalsIgnoreCase(nullstate))
					nonCnt++;
				else if (strmatch(state, refStatesLst)){ 
					enhCnt++;
				}
				else{
					elseCnt++;
				}
				
			}
			
			if ((enhCnt / (enhCnt+elseCnt) ) >= 0.5)
				resBed.addRefSeq(segmentGen);
				
		
		}
		
		resBed.writeFullBed(outPrefix + "majoritySelectedRegions.bed");
	 
		
	}

	public static boolean strmatch(String state, List<String> refStatesLst) {

	
		for (int i=0; i<refStatesLst.size();i++)
		{
			if (state.equalsIgnoreCase(refStatesLst.get(i)))
				return true;
		}
		return false;
	}


	//Given a bed file, segment it start stop regions to X bases regions
	private static BEDFileParser segmentRegions(BEDFileParser refRegions,
			int segmentSize) {
		
		BEDFileParser b = new BEDFileParser();
		for (RefSeqGene g: refRegions.GetGenes()){
			
			int start = g.getStart();
			int end= g.getEnd();
			int numBins= (int) (end-start)/segmentSize;
			if (numBins <2){
				b.addRefSeq(g);
			}
			else{
				int cstart=start;
				int cend= start+segmentSize-1;
				for (int i=0; i<numBins; i++){
					List<Integer> es = new LinkedList<Integer>(); 
					es.add(cstart);
					List<Integer> ee = new LinkedList<Integer>(); 
					ee.add(cend);
					
					String name = g.getChr()+":"+cstart+"-"+cend;
					RefSeqGene gcurr= new RefSeqGene(g.getChr(),cstart,cend,name,g.getOrientation(),es,ee);
					b.addRefSeq(gcurr);
					cstart= cend+1;
					cend= cstart+segmentSize-1;
					if (i==numBins-2)
						cend=end;
				}
			}
			
		}
		return b;
	}


	private static Map<String, String> loadSetsPath(String hmmList) throws IOException {

		BufferedReader bw = new BufferedReader(new FileReader(hmmList)) ;
		String line;
		Map<String, String> res = new HashMap<String, String>();
		while((line=bw.readLine())!=null){
			String[] a= line.split ("\t");
			res.put(a[0],a[1]);
		}
		
		return res;
	}


	
	
	

	private static void analyzeDivergentK4Overlap(String neighborTab,
	String lincBedF, String refGeneBedF, String k4bedF, String out) throws IOException {

		BEDFileParser lincsBed = new BEDFileParser (lincBedF);
		BEDFileParser genesBed = new BEDFileParser (refGeneBedF);
		BEDFileParser k4bed = new BEDFileParser (k4bedF);
		
		HashMap <String, RefSeqGene> lincMap = lincsBed.getNameGeneMap();
		HashMap <String, RefSeqGene> geneMap = genesBed.getNameGeneMap();
		
		//
		BufferedReader bw = new BufferedReader(new FileReader( neighborTab)) ;
		String line;
		line = bw.readLine();
		String[] header = line.split("\t");
		HashMap <String, Integer> headerHash = new HashMap();
		for (int i=0; i<header.length; i++) {headerHash.put( header[i],i);}
		
		int isDivColR = headerHash.get("isDivR");
		int isDivColL = headerHash.get("isDivL");
		int gR = headerHash.get("rightNeighborName");
		int gL = headerHash.get("leftNeighborName");
		int dR = headerHash.get("distRight");
		int dL = headerHash.get("distLeft");
		int dist = 0;
		
		BufferedWriter bw2 = new BufferedWriter(new FileWriter( out)) ;
		bw2.write("geneName\tlincName\tStartGene\tEndGene\tStartlnc\tEndlnc\tStartK4\tEndK4\tDistance\n");
		
		while((line=bw.readLine())!=null){
			String[] a= line.split ("\t");
			String gene ="";
			if ( a[isDivColR].equals("1")){
				gene = a[gR];
				dist = Integer.valueOf(a[dR]);
			}
			if ( a[isDivColL].equals("1")){
				gene = a[gL];
				dist = Integer.valueOf(a[dL]);
			}
			String l=a[headerHash.get("name")];
			if (! gene.equals("") & geneMap.containsKey(gene) & lincMap.containsKey(l) ){
				RefSeqGene g = new RefSeqGene(geneMap.get(gene));
				RefSeqGene lnc = new RefSeqGene (lincMap.get(l));
				g.expandUtrs(1000,1000);
				lnc.expandUtrs(1000,1000);
				
				if (k4bed.containChr(g.getChr())) {
				Iterator<Node<RefSeqGeneWithIsoforms>> intervals =k4bed.getChrTree(g.getChr()).iterator(g.getStart(), g.getEnd());
				
					while (intervals.hasNext()){
						RefSeqGeneWithIsoforms k4 = intervals.next().getValue();
						if (k4.overlapsGeneInAnyOrientation(lnc) & k4.overlapsGeneInAnyOrientation(g) ){
							bw2.write(g.getName()+"\t"+lnc.getName()+"\t"+g.getStart()+"\t"+g.getEnd()+"\t"+ lnc.getStart()+"\t"+lnc.getEnd()+"\t"+k4.getStart()+"\t"+k4.getEnd()+"\t"+dist+"\n");
							
						}
					}
				}
			}
		}
		
		bw2.close();
		bw.close();
	}
	
	
	
	
	
private static void extractLongestIsoFromGTF (String gtf_name, String geneLst_name, String geneMap_name, String out) throws IOException{
		
		
		GTFFileParser gtf= new GTFFileParser(gtf_name);
		HashMap <String, RefSeqGene> lincMap = gtf.getNameGeneMap();
		HashMap<String, LinkedList<String>>  geneIsoMap = loadIsoGeneMap(geneMap_name);
		ArrayList<String> genelst=loadNameLst(geneLst_name,false);
		BEDFileParser outbed = new  BEDFileParser();
		
		for (String geneId : genelst){
			if (geneIsoMap.containsKey(geneId)) {
				ArrayList <RefSeqGene> isoSet = new ArrayList <RefSeqGene>();
				LinkedList<String> isoIds = geneIsoMap.get(geneId);
				for (String iso : isoIds){
					if (lincMap.containsKey(iso)){
						RefSeqGene tcon = lincMap.get(iso);
						isoSet.add(tcon);
					}
				}
				if   (! isoSet.isEmpty()) {
					int maxlen = 0; 
					RefSeqGene selectTcon =isoSet.get(0) ;
					for (RefSeqGene r : isoSet){
						int len = r.getTranscriptLength();
						if (len > maxlen){
							selectTcon = r;
							maxlen= len;
						}
					}
					outbed.addRefSeq(selectTcon);
				}
			}
		}
		BufferedWriter bw2 = new BufferedWriter(new FileWriter( out)) ;
		outbed.bed2CuffGtf(bw2, "longTcon", true);
		bw2.close();
			
		
		
 }
	
	
	
	private static void BestMatchCuffdiffSegmentPerGene(String gtf,	String diff, String out , int segmentSizeThreshold) throws IOException {


		//merge GTF by gene id
		BEDFileParser inbed= new BEDFileParser();
		inbed.loadGTF(new File(gtf),"");
		BEDFileParser mergedBed =  inbed.mergeByGtfGeneId();
	
		HashMap<String,RefSeqGene> nameGeneMap = mergedBed.getNameGeneMap();
		HashMap<String,HashMap<String,String>> geneDiffMap = new HashMap<String,HashMap<String,String>>();
		HashMap<String,HashMap<String,Double>> geneBestDiffMap = new HashMap<String,HashMap<String,Double>>();
		HashMap<String,String> geneBestEnrichMap_name = new HashMap<String,String>();
		HashMap<String,Double> geneBestEnrichMap_val = new HashMap<String,Double>();
			
		//read diff - store in a hash of gene ids , maintain best scoring
		String line;
		BufferedReader br = new BufferedReader (new FileReader (diff));
		String headerline=br.readLine();
		while((line=br.readLine())!=null){
			String[] lineArr= line.split ("\t");
			String segment = lineArr[3];
			Double expFPKM = Double.valueOf(lineArr[8]);
			Double log2enrich = Double.valueOf(lineArr[9]);
			String pvalSig = lineArr[13]; 
			
			String[] a = segment.split(":");
			String[] b = a[1].split("-");
			Alignments segAln = new Alignments (a[0],b[0],b[1]);
			IntervalTree<RefSeqGeneWithIsoforms> overlappers = mergedBed.getOverlappers(segAln);
			
			Iterator<RefSeqGeneWithIsoforms> overlapperIt = overlappers.valueIterator();
			while (overlapperIt.hasNext()) {
				RefSeqGeneWithIsoforms gene = overlapperIt.next();
				String geneName = gene.getName();
				if (! geneDiffMap.containsKey(geneName))
				{
					geneDiffMap.put(geneName, new HashMap<String,String>());
					geneBestDiffMap.put(geneName, new HashMap<String,Double>());
					geneBestEnrichMap_name.put(geneName,"");
					geneBestEnrichMap_val.put(geneName,new Double(0));
					
				}
				geneDiffMap.get(geneName).put(segment, line);
				geneBestDiffMap.get(geneName).put(segment, expFPKM);
				
				if ( (geneBestEnrichMap_val.get(geneName) <= expFPKM & segAln.length() > segmentSizeThreshold ) || geneBestEnrichMap_name.get(geneName).equals("") ){
					geneBestEnrichMap_val.put(geneName,expFPKM);
					geneBestEnrichMap_name.put(geneName,segment) ; 
				}
				
				
			}
		}
			
				
		//Write out file
		BufferedWriter bw2 = new BufferedWriter(new FileWriter( out + "bestHit.diff")) ;
		BufferedWriter bw3 = new BufferedWriter(new FileWriter( out + "allHit.diff")) ;
		
		bw3.write("GeneName\t"+headerline+"\n");
		bw2.write("GeneName\t"+headerline+"\n");
		
		for (String geneName: nameGeneMap.keySet()) {
			
			if (geneBestEnrichMap_name.containsKey(geneName)) {
				
				String bestSegment = geneBestEnrichMap_name.get(geneName); 
				String bestLine = geneDiffMap.get(geneName).get(bestSegment);
				bw2.write(geneName+"\t"+bestLine+"\n");
				
				for (String segs : geneDiffMap.get(geneName).keySet() ){
					String cline = geneDiffMap.get(geneName).get(segs);
					bw3.write(geneName+"\t"+cline+"\n");
					
				}
			}
			
			
		}
		
		
		bw2.close();
		bw3.close();
		
	}
	

	/**
	 * add attributes to isoforms in the GTF that are represented in the annotGTF
	 * @param gtf input GTF that would be annotated with additional attributes based on annot GTF
	 * @param annotGTF has the aditional attributes
	 * @param annotPrefix  prefix to be added for each attribute added to the GTF
	 * @param outf  output file name
	 * 
	 * @throws IOException 
	 */
	
 
	private static void ReannotateGTF(String gtf_name, String annotGTF_name,
	String annotPrefix, String outf,String src) throws IOException  {

		GTFFileParser gtf= new GTFFileParser(gtf_name);
		GTFFileParser annotGtf= new GTFFileParser(annotGTF_name);
		
		for (RefSeqGene annotG : annotGtf.GetGenes()){
			
			IntervalTree<RefSeqGeneWithIsoforms> overlappers = gtf.getOverlappers(annotG);
			
			Iterator<RefSeqGeneWithIsoforms> it = overlappers.valueIterator();
			while (it.hasNext()){
				RefSeqGeneWithIsoforms g = it.next();
				for (RefSeqGene iso: g.getAllIsoforms()){
					 boolean isComp= annotG.isFullyCompatible(iso);
					 if (isComp ){
						 Map<String,String > attrs= annotG.getAttributes();
						 for (String key : attrs.keySet() )
						 {
							 if (! key.equals("exon_number"))
							 iso.addAttribute(annotPrefix+key, attrs.get(key));
						 }
					 }
					
				}
				
			}
			
		}
		BufferedWriter bw = new BufferedWriter(new FileWriter( outf)) ;
		gtf.GTFWrite(bw, src,true);
		bw.close();

}


	/** TO DO!
	 * add attributes to isoforms in the GTF that are represented in the annotGTF
	 * @param gtf input GTF that would be annotated with additional attributes based on annot GTF
	 * @param annotTabf has the aditional attributes mapped to  headerIdenitifier (which is key in the GTF)
	 * @param annotPrefix  prefix to be added for each attribute added to the GTF
	 * @param annotHeaderList  List of attr from the tab to be added as attrs to gtf
	 * @param headerIdenitifier  the tab identifer that matches the gtf identifier to map accordingly 
	 * @param attrIdenitifier  the gtf attribute identifier to map accordingly 
	 * @param outf  output file name
	 * 
	 * @throws IOException 
	 */
	
	
	private static void ReannotateGTFfromTable (String gtf_name,
			String annotTabf,String annotPrefix,String annotHeaderList, String headerIdenitifier,
			String attrIdenitifier, String outf)throws IOException  {
		/*
				GTFFileParser gtf= new GTFFileParser(gtf_name);
				GTFFileParser annotTAB= new GTFFileParser(annotTabf);
				
				for (RefSeqGene annotG : annotGtf.GetGenes()){
					
					IntervalTree<RefSeqGeneWithIsoforms> overlappers = gtf.getOverlappers(annotG);
					
					Iterator<RefSeqGeneWithIsoforms> it = overlappers.valueIterator();
					while (it.hasNext()){
						RefSeqGeneWithIsoforms g = it.next();
						for (RefSeqGene iso: g.getAllIsoforms()){
							 boolean isComp= annotG.isFullyCompatible(iso);
							 if (isComp ){
								 Map<String,String > attrs= annotG.getAttributes();
								 for (String key : attrs.keySet() )
								 {
									 iso.addAttribute(annotPrefix+key, attrs.get(key));
								 }
							 }
							
						}
						
					}
					
				}
				BufferedWriter bw = new BufferedWriter(new FileWriter( outf)) ;
				gtf.bed2CuffGtf(bw, annotPrefix,true);
				bw.close();
*/
		}

	
	
	private static void CapNumIsoformsPerLoci(String gtf_name, String outf,
			Integer maxIsoPerLoci, String src) throws IOException {
		
		GTFFileParser gtf= new GTFFileParser(gtf_name);	
		GTFFileParser gtfout = new GTFFileParser();
		HashMap <String,Locus> geneIDIsoMap =  gtf.getGeneIDIsoMap();
		
		BufferedWriter bw2 = new BufferedWriter(new FileWriter( outf+"FilteredGeneIDs.txt")) ;
		
		for (String geneId: geneIDIsoMap.keySet())
		{
			
			Locus loc = geneIDIsoMap.get(geneId);
			int numIso= loc.getNumberOfIsoforms();
			if ( numIso> maxIsoPerLoci) 
				{
					gtfout.addRefSeqSet(loc.selectRandIsoSubset(maxIsoPerLoci));
					bw2.write(geneId+"\t"+numIso+"\n");
				}
			
			else
				gtfout.addRefSeqSet(loc.getAllIsoforms());
			
			
		}
		
		
		BufferedWriter bw = new BufferedWriter(new FileWriter( outf)) ;
		gtfout.GTFWrite(bw, src, true);
		bw.close();
		bw2.close();
	}

}