package broad.pda.annotation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.core.annotation.AnnotationReader;
import broad.core.annotation.BED;
import broad.core.annotation.BEDReader;
import broad.core.annotation.BasicAnnotationReader;
import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.GFF;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.annotation.PSL;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.gene.RefSeqGeneWithIsoforms;


public class BEDFileParser {

	//////	
	//Addition by Moran May 18th 2010

	private Map<String, IntervalTree<RefSeqGeneWithIsoforms>>  annotationSetMap;
	private Map<String, List<RefSeqGeneWithIsoforms>> nameAnnotationMap;
	static Logger logger = Logger.getLogger(BEDFileParser.class.getName());

	public BEDFileParser() {
		this.annotationSetMap = new TreeMap<String, IntervalTree<RefSeqGeneWithIsoforms>>();
		this.nameAnnotationMap = new HashMap<String, List<RefSeqGeneWithIsoforms>>();
	}

	public BEDFileParser(String fileName) throws IOException {
		this.annotationSetMap = new TreeMap<String, IntervalTree<RefSeqGeneWithIsoforms>>();
		this.nameAnnotationMap = new HashMap<String, List<RefSeqGeneWithIsoforms>>();
		loadIsoformDataByChrToTree(new File(fileName));
	}	
	
	public BEDFileParser(File file) throws IOException {
		this.annotationSetMap = new TreeMap<String, IntervalTree<RefSeqGeneWithIsoforms>>();
		this.nameAnnotationMap = new HashMap<String, List<RefSeqGeneWithIsoforms>>();
		loadIsoformDataByChrToTree(file);
	}


	//This Ctor upload only the part of a transcript that is 3PrimeTrimSize in length from the 3 prime end
	public BEDFileParser(String fileName, int ThreePrimeTrimSize ) throws IOException {
		this.annotationSetMap = new TreeMap<String, IntervalTree<RefSeqGeneWithIsoforms>>();
		this.nameAnnotationMap = new HashMap<String, List<RefSeqGeneWithIsoforms>>();
		trimAndLoadIsoformDataByChrToTree(new File(fileName),ThreePrimeTrimSize);
	}

	public BEDFileParser(String fileName, String chr) throws IOException {
		this.annotationSetMap = new TreeMap<String, IntervalTree<RefSeqGeneWithIsoforms>>();
		this.nameAnnotationMap = new HashMap<String, List<RefSeqGeneWithIsoforms>>();
		loadIsoformDataByChrToTree(new File(fileName),chr);
	}
	
	public List<RefSeqGeneWithIsoforms> getAnnotations(String name) {
		return nameAnnotationMap.get(name);
	}

	/*public BEDFileParser(String fileName , boolean GFF) throws IOException {
		this.annotationSetMap = new TreeMap<String, IntervalTree<RefSeqGeneWithIsoforms>>();
		parseGFF(new File(fileName));
	}	*/


	public BEDFileParser(Collection<RefSeqGene> genes) {
		this();
		for(RefSeqGene gene : genes) {
			addRefSeq(gene);
		}
	}

	public BEDFileParser(Map<String, Collection<RefSeqGene>> transcripts) {
		this();
		for (String chr : transcripts.keySet()) {
			Collection<RefSeqGene> chrTranscripts = transcripts.get(chr);
			for(RefSeqGene transcript : chrTranscripts) {
				addRefSeq(transcript);
			}
		}
	}

	public RefSeqGeneWithIsoforms get(String geneName) {
		if(nameAnnotationMap.containsKey(geneName)) {
			return nameAnnotationMap.get(geneName).get(0);
		} else {
			return null;
		}
	}

	/////////////////////////////////////////////////////////////////////
	
	// TODO: GFF TO BED, implement mitch's version
	
	///////////////////////////////////////////////////////////////////////////////////////////////////
	
	public List<RefSeqGene> GetGenes(){
	
		List<RefSeqGene> rtrn=new ArrayList<RefSeqGene>();
		Iterator<String> chrIt=this.getChromosomeIterator();
		while(chrIt.hasNext()){
			String chr=chrIt.next();
			rtrn.addAll(getChrRefSeqGenes(chr));
		}
		return rtrn;
	}

	
	public List<RefSeqGeneWithIsoforms> GetGenesWithIsoforms(){
		
		List<RefSeqGeneWithIsoforms> rtrn=new ArrayList<RefSeqGeneWithIsoforms>();
		Iterator<String> chrIt=this.getChromosomeIterator();
		while(chrIt.hasNext()){
			String chr=chrIt.next();
			rtrn.addAll(getChrRefSeqGenesWithIsoforms(chr));
		}
		return rtrn;
	}
	
	public IntervalTree<RefSeqGeneWithIsoforms> getChrTree(String chr){
		if (this.annotationSetMap.containsKey(chr))
			return this.annotationSetMap.get(chr);
		else
			return null;
	}

	private List<RefSeqGeneWithIsoforms> getChrRefSeqGenesWithIsoforms(String chr) {
		return (List<RefSeqGeneWithIsoforms>) getChrRefSeqGenes(chr,false);
	}
	
	private List<? extends RefSeqGene> getChrRefSeqGenes(String chr) {
		return getChrRefSeqGenes(chr,true);
	}
	
	private List<? extends RefSeqGene> getChrRefSeqGenes(String chr, boolean seperateIsoformes) {
	
		List<RefSeqGene> rtrn=new ArrayList<RefSeqGene>();
		List<RefSeqGeneWithIsoforms> rtrn2 = new ArrayList<RefSeqGeneWithIsoforms>();
	
			
		
		if (this.annotationSetMap.containsKey(chr))
		{
			IntervalTree<RefSeqGeneWithIsoforms> tree=this.getChrTree(chr);
			Iterator<RefSeqGeneWithIsoforms> geneIt=tree.valueIterator();
			while(geneIt.hasNext())
			{
				RefSeqGeneWithIsoforms gene=geneIt.next();
				rtrn.addAll(gene.getAllIsoforms());
				rtrn2.add(gene);
			}
		}
		if (seperateIsoformes)
			return rtrn;
		return rtrn2;
	}

	public Map<String, IntervalTree<RefSeqGeneWithIsoforms>> getIntervalTreeWithIsoforoms(){return this.annotationSetMap;}
	
	public Iterator<String> getChromosomeIterator() {
	return this.annotationSetMap.keySet().iterator();}

	//Return only things that overlap an exon!
	public IntervalTree<RefSeqGeneWithIsoforms> getOverlappers(RefSeqGene gene) {
	
		IntervalTree<RefSeqGeneWithIsoforms> overlapTree= new IntervalTree<RefSeqGeneWithIsoforms>();
	
		if (gene!= null & this.containChr(gene.getChr())){
			IntervalTree<RefSeqGeneWithIsoforms> chrTree= this.getChrTree(gene.getChr());
	
			//System.err.println("gene " + gene.toBED() );
	
			if (chrTree!=null /*&& ! chrTree.doesOverlap(gene.getStart(), gene.getEnd())*/) {
				//System.err.println("\tchrTree overlaps gene");
				Iterator<RefSeqGeneWithIsoforms> overlapperIt = new IntervalTree.ValuesIterator<RefSeqGeneWithIsoforms>(chrTree.overlappers(gene.getStart(), gene.getEnd()));
				//System.err.println("Gene " + gene.toBED() + " has overlappers? " + overlapperIt.hasNext());
				while(overlapperIt.hasNext()) {
					RefSeqGeneWithIsoforms overlapper = overlapperIt.next();
					//System.err.println("\tGoing through loci overlaper, testing " + overlapper.toBED() );
					//if(overlapper == null) {throw new IllegalStateException("Ovelpper was null when looking for overlappers for gene " + gene.toBED());}
					if ( overlapper.overlaps(gene)  ) {
						//System.err.println("\t\tYES. Overlapper " + overlapper.toUCSC() + " overlaps with gene, going to add");
						Node<RefSeqGeneWithIsoforms> alreadyFound = overlapTree.find(overlapper.getStart(), overlapper.getEnd());
						if (alreadyFound == null ) {
							overlapTree.put(overlapper.getStart(), overlapper.getEnd(), overlapper);
						}
						else //there is already one isoform that spans the same interval choose the larger isoform
						{
							alreadyFound.getValue().AddAllIsoforms(overlapper);
	
						}
					} 
				}
			}	
		}
		return overlapTree;
	}

	public IntervalTree<RefSeqGeneWithIsoforms> getOverlappers(LightweightGenomicAnnotation region) {
		RefSeqGene artificialGene = new RefSeqGene(region);
		return getOverlappers(artificialGene);
	}

	public Iterator<RefSeqGeneWithIsoforms> getGenomicRegionOverlappers(LightweightGenomicAnnotation region) {
	
		Iterator<RefSeqGeneWithIsoforms> overlapperIt= new ArrayList<RefSeqGeneWithIsoforms>().iterator();
	
		if (region!= null & this.containChr(region.getChromosome())){
			IntervalTree<RefSeqGeneWithIsoforms> chrTree= this.getChrTree(region.getChromosome());
	
			//System.err.println("gene " + gene.toBED() );
	
			if (chrTree!=null /*&& ! chrTree.doesOverlap(gene.getStart(), gene.getEnd())*/) {
				//System.err.println("\tchrTree overlaps gene");
				overlapperIt = new IntervalTree.ValuesIterator<RefSeqGeneWithIsoforms>(chrTree.overlappers(region.getStart(), region.getEnd()));
	
			}	
		}
		return overlapperIt;
	}

	public BEDFileParser getMergedCopy() {
	
		BEDFileParser bed = new BEDFileParser();
		bed.addRefSeqSet(this.GetGenes());
		bed.merge();
		return bed;
	}

	public Map<String, IntervalTree<RefSeqGeneWithIsoforms>> getMergedAnnotationMap(double minPctOverlapToMerge) {
			Map<String, IntervalTree<RefSeqGeneWithIsoforms>> mergedTreeMap = new LinkedHashMap<String, IntervalTree<RefSeqGeneWithIsoforms>>();
	
		Iterator<String> chrIt = this.annotationSetMap.keySet().iterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			IntervalTree<RefSeqGeneWithIsoforms> newTree = new IntervalTree<RefSeqGeneWithIsoforms>();
			mergedTreeMap.put(chr, newTree);
			//Collection<? extends RefSeqGene> oldTree = this.getChrRefSeqGenes(chr);
			IntervalTree<RefSeqGeneWithIsoforms> oldTree = getChrTree(chr);
			Iterator<RefSeqGeneWithIsoforms> oldGeneIt=oldTree.valueIterator();
			while(oldGeneIt.hasNext()) {
				RefSeqGeneWithIsoforms currentElement =oldGeneIt.next();
				//System.err.println("Isoforms " + currentElement.getAllIsoforms().size() + " CurrentElement " + currentElement.toBED());
				Iterator<RefSeqGeneWithIsoforms> overlapperIt = new IntervalTree.ValuesIterator<RefSeqGeneWithIsoforms>(newTree.overlappers(currentElement.getStart(), currentElement.getEnd()));
				RefSeqGeneWithIsoforms overlapper = null;
				int overlapperNum = 0;
				while(overlapperIt.hasNext() ) {
					RefSeqGeneWithIsoforms overlapperCandidate= overlapperIt.next();
					if( isOverlapCompatible(currentElement,overlapperCandidate, minPctOverlapToMerge) ) {
						overlapper = overlapperCandidate;
					}
				
					if(overlapper != null) {
						overlapperNum++;
						newTree.remove(overlapper.getStart(),overlapper.getEnd());
						//RefSeqGene mergedElement1=  overlapper.takeUnion(currentElement);//Bug fix : Dec 27 2010
						RefSeqGene mergedElement1=  overlapper.takeUnion(currentElement.getMerged());
						RefSeqGeneWithIsoforms mergedElement=new RefSeqGeneWithIsoforms(mergedElement1);
						newTree.put(mergedElement.getStart(), mergedElement.getEnd(), mergedElement);
					} 
					if(overlapperNum > 1 ) {
						//System.err.println("WARN: Merging annotations error, more than one different gene overlaps current element " + currentElement.toBED() );
					}
					overlapper = null;
				}
				if(overlapperNum == 0) {
					//Bug fix : Dec 27 2010 - the function was not merging a refseq gene with isoforms that had several isoforms
					RefSeqGeneWithIsoforms mergedElement=new RefSeqGeneWithIsoforms(currentElement.getMerged());
					newTree.put(currentElement.getStart(), currentElement.getEnd(),  mergedElement);
					//newTree.put(currentElement.getStart(), currentElement.getEnd(),  currentElement);
				}
			}
		}
		return mergedTreeMap;
	}

	public Map<String, IntervalTree<RefSeqGeneWithIsoforms>> getMergedAnnotationMap(boolean mergall) {
		return getMergedAnnotationMap(0);
	}

	public Map<String, IntervalTree<RefSeqGeneWithIsoforms>> getAnnotationSetMap() { return this.annotationSetMap;}

	public int getNumberOfIsoforms(String chr) {
	
		int res=0;
		if (this.annotationSetMap.containsKey(chr)){
			res=this.getChrRefSeqGenes(chr).size();
		}
		return res;
	}

	/**
	 * 
	 * @return Number of transcripts in the reader
	 */
	public int getNumberOfIsoforms() {
	
		int res=0;
		Iterator<String> chrIt = getChromosomeIterator();
		while(chrIt.hasNext()) {
			res += this.getChrRefSeqGenes(chrIt.next()).size();
		}
	
		return res;
	}

	/////////////////////////////////////////////////////////////////////
	
	// TODO: GFF TO BED, implement mitch's version
	
	
	
	
	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////
	
	public static double getNumberOfIsoforms(
			IntervalTree<RefSeqGeneWithIsoforms> overlapTree) {
		double res=0;
		Iterator<RefSeqGeneWithIsoforms> itG=overlapTree.valueIterator();
		while(itG.hasNext())
		{
			RefSeqGeneWithIsoforms g=itG.next();
			res+=g.getAllIsoforms().size();
		}
		return res;
	}

	public RefSeqGene getExactIsoform(RefSeqGeneWithIsoforms loci) {
		RefSeqGene res=null;
		Iterator<RefSeqGeneWithIsoforms> itIso=this.getOverlappers(loci).valueIterator();
		while (itIso.hasNext()){
			Collection<RefSeqGene> gSet = itIso.next().getAllIsoforms();
			for (RefSeqGene g: gSet){
				if (g.equals(loci)){
					res=g;
					return res;
				}
			}
		}
		return res;
	}

	public HashMap<String, RefSeqGene> getNameGeneMap() {
		HashMap<String, RefSeqGene> nameGeneMap = new HashMap<String, RefSeqGene>();
		for (RefSeqGene g : this.GetGenes())
			nameGeneMap.put(g.getName(), g);
		return nameGeneMap;
	}

	public BEDReader get3PRegion(int toIn, int toOut) {
		BEDReader utrs = new BEDReader();
		Iterator <String> chrIt = getChromosomeIterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			Iterator<RefSeqGeneWithIsoforms> geneIt = getChrTree(chr).valueIterator();
			while(geneIt.hasNext()) {
				RefSeqGeneWithIsoforms gene = geneIt.next();
				int end = gene.getOrientedEnd();
				int threePStart = gene.getOrientation().equals("-") ? end - toOut : end - toIn;
				int threePEnd = gene.getOrientation().equals("-") ? end + toIn : end + toOut;
				BED utr = new BED(gene.getName(), chr, threePStart, threePEnd);
				utr.setOrientation(gene.getOrientation());
				utrs.addAnnotation(utr);
			}
		}
		return utrs;
	}

	public BEDReader getPromoters() {
		int defaultIn = 250;
		int defaultOut = 500;
		
		return getPromoters(250, 500);
	}

	public BEDReader getPromoters(int toIn, int toOut) {
		BEDReader promoters = new BEDReader();
		Iterator <String> chrIt = getChromosomeIterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			Iterator<RefSeqGeneWithIsoforms> geneIt = getChrTree(chr).valueIterator();
			while(geneIt.hasNext()) {
				RefSeqGeneWithIsoforms gene = geneIt.next();
				int start = gene.getOrientedStart();
				int promoterStart = gene.getOrientation().equals("-") ? start - toIn : start - toOut;
				int promoterEnd = gene.getOrientation().equals("-") ? start + toOut : start + toIn;
				BED promoter = new BED(gene.getName(), chr, promoterStart, promoterEnd);
				promoter.setOrientation(gene.getOrientation());
				promoters.addAnnotation(promoter);
			}
		}
		return promoters;
	}

	

	/////////////////////////////////////////////////////////////////////
	
	// TODO: GFF TO BED, implement mitch's version
	
	
	
	
	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////
	
	private static List<Integer> [] setBlockStartsAndEnds(String[] blockStarts, String[] blockSizes, int size, int start){
		List<Integer> starts=new ArrayList<Integer> ();
		List<Integer>  end=new ArrayList<Integer> ();
		for(int i=0; i<size; i++){
			starts.add(start+new Integer(blockStarts[i].replaceAll("\"", "").trim()));
			end.add((Integer)starts.get(i)+new Integer(blockSizes[i].replaceAll("\"", "").trim()));
		}
		List [] rtrn={starts, end};
		return rtrn;
	}

	private void setAnnotationSetMap(Map<String, IntervalTree<RefSeqGeneWithIsoforms>> mergedTreeMap) {
		this.annotationSetMap=mergedTreeMap;
		this.nameAnnotationMap.clear();
		Iterator<String> chrIt = this.annotationSetMap.keySet().iterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			Iterator<RefSeqGeneWithIsoforms> it = annotationSetMap.get(chr).valueIterator();
			while(it.hasNext()) {
				RefSeqGeneWithIsoforms geneWithIsos = it.next();
				Iterator<RefSeqGene> isos = geneWithIsos.getAllIsoforms().iterator();
				while(isos.hasNext()) {
					RefSeqGene g = isos.next();
					updateAnnotationNameMap(g.getName(), geneWithIsos);
				}
			}
		}
	
	
	}
	
	private void updateAnnotationNameMap(String name, RefSeqGeneWithIsoforms geneWithIsos) {
		if(nameAnnotationMap.containsKey(name)) {
			nameAnnotationMap.get(name).add(geneWithIsos);
		} else {
			List<RefSeqGeneWithIsoforms> nameList = new ArrayList<RefSeqGeneWithIsoforms>();
			nameList.add(geneWithIsos);
			nameAnnotationMap.put(name, nameList);	
		}
	}

	public static Set<RefSeqGene> loadData(File file) throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	
		Set<RefSeqGene> data=new TreeSet<RefSeqGene>();
		String nextLine;
		int i=0;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			if(looksLikeData(nextLine)){
				RefSeqGene track=new RefSeqGene(nextLine, false);
				data.add(track);
				i++;
				if(i%10000 ==0){System.err.println(i);}
			}
	
		}
	
	
		reader.close();
		return data;
	
	}

	public static Set<RefSeqGene> loadData(File file, Alignments region) throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	
		Set<RefSeqGene> data=new TreeSet<RefSeqGene>();
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			String[] tokens=nextLine.split("\t");
			if(looksLikeData(nextLine) && tokens.length>11){          	
				String chr=(tokens[0]);
				int start=new Integer(tokens[1]);
				int end=new Integer(tokens[2]);
	
				String name=tokens[3];
				String strand=tokens[5];
	
				String[] blockSizes=tokens[10].split(",");
				String[] blockStarts=tokens[11].split(",");
				List<Integer>[] exonStartEnd=setBlockStartsAndEnds(blockStarts, blockSizes, new Integer(tokens[9]), start);
	
				String [] extraColumns = null;
				if(tokens.length > 12) {
					extraColumns = new String[tokens.length - 12];
					for(int j = 12; j < tokens.length; j++) {
						extraColumns[j - 12] = tokens[j];
					}
				}
				RefSeqGene track=new RefSeqGene(chr, start, end, name, strand, exonStartEnd[0], exonStartEnd[1], extraColumns);
				if(track.getAlignment().overlapsAtAll(region)){data.add(track);}
			}
	
		}
	
	
		reader.close();
		return data;
	
	}

	public static Map<String, RefSeqGene> loadDataByName(File file) throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	
		Map<String, RefSeqGene> data=new TreeMap<String, RefSeqGene>();
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	
			if(looksLikeData(nextLine)){
	
				String[] tokens=nextLine.split("\t");
				String chr=(tokens[0]);
				int start=new Integer(tokens[1]);
				int end=new Integer(tokens[2]);
	
				String name=tokens[3].toUpperCase().trim();
				String strand=tokens[5];
	
				int blockStart=new Integer(tokens[6]);
				int blockEnd=new Integer(tokens[7]);
	
				String[] blockSizes=tokens[10].split(",");
				String[] blockStarts=tokens[11].split(",");
				List<Integer>[] exonStartEnd=setBlockStartsAndEnds(blockStarts, blockSizes, new Integer(tokens[9]), start);
	
				String [] extraColumns = null;
				if(tokens.length > 12) {
					extraColumns = new String[tokens.length - 12];
					for(int j = 12; j < tokens.length; j++) {
						extraColumns[j - 12] = tokens[j];
					}
				}
	
				RefSeqGene track=new RefSeqGene(chr, start, end, name, strand, exonStartEnd[0], exonStartEnd[1], blockStart, blockEnd, extraColumns);
				data.put(name, track);
			}
	
		}
	
	
		reader.close();
		return data;
	
	}

	public static Map<String, Collection<RefSeqGene>> loadDataByChr(File file) throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	
		Map<String, Collection<RefSeqGene>> rtrn=new TreeMap<String, Collection<RefSeqGene>>();
		String nextLine;
		int i=0;
		while ((nextLine = reader.readLine()) != null ) {
	
			if(looksLikeData(nextLine) ){
				RefSeqGene gene = new RefSeqGene(nextLine, false);
				//System.err.println("Gene: " + gene.toBED());
	
				Collection<RefSeqGene> data=new TreeSet<RefSeqGene>();
				if(rtrn.containsKey(gene.getChr())){
					data=rtrn.get(gene.getChr());
				}
				data.add(gene);
				rtrn.put(gene.getChr(), data);
			}
			i++;
			if(i%10000==0){System.err.println(i);}
		}
	
	
		reader.close();
		return rtrn;
	
	}

	public static Map<String, IntervalTree<RefSeqGene>> loadDataByChrToTree(File file) throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	
		Map<String, IntervalTree<RefSeqGene>> rtrn=new TreeMap<String, IntervalTree<RefSeqGene>>();
		String nextLine;
		int i=0;
		int IsoformMissCntr=0;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	
			if(looksLikeData(nextLine)){
	
				RefSeqGene gene = new RefSeqGene(nextLine, false);
	
	
				IntervalTree<RefSeqGene> data = rtrn.get(gene.getChr());
				if(data == null){
					data = new IntervalTree<RefSeqGene>();
					rtrn.put(gene.getChr(), data);
				}
				if (data.find(gene.getStart(), gene.getEnd()) != null)
					IsoformMissCntr++;
				data.put(gene.getStart(), gene.getEnd(), gene);
			}
			i++;
			if(i%10000==0){System.err.println(i);}
		}
	
	
		reader.close();
		System.err.println("While loading "+file.getAbsolutePath()+ " missed " + IsoformMissCntr + " isoforms"  );
		return rtrn;		 
	}

	public static Map<String, Integer> loadChrSizes(String file){
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
	
		try{	
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	
	
			String nextLine;
			while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
				String[] tokens=nextLine.split("\t| +");
				//System.err.println(nextLine+" "+tokens.length);
	
				rtrn.put(tokens[0], new Integer(tokens[1]));
			}
	
	
			reader.close();
		}catch(IOException ex){ex.printStackTrace();}
		return rtrn;
	
	}

	public static List<String> loadList(String file) throws IOException {
	    BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	
	    List<String> data=new ArrayList<String>();
	        String nextLine;
	        int i=0;
	        while ((nextLine = reader.readLine()) != null) {
	          data.add(nextLine);
	
	        }
	
	
	        reader.close();
	        return data;
	}

	public static List<String> loadList(String file, boolean skipFirstLine) throws IOException {
	    BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	
	    List<String> data=new ArrayList<String>();
	        String nextLine;
	        int i=0;
	        while ((nextLine = reader.readLine()) != null) {
	            if(i>0){
	                data.add(nextLine);
	            }
	            else{
	                if(!skipFirstLine){data.add(nextLine);}
	            }
	
	            i++;
	        }
	
	
	        reader.close();
	        return data;
	}

	protected void loadIsoformDataByChrToTree(File file) throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		int i = 0;
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			if(looksLikeData(nextLine)){
				
				RefSeqGene gene=new RefSeqGene(nextLine, false);	    			
				this.addRefSeq(gene);
			}
			i++;
			if(i%10000==0){logger.debug(i);}
		}
		reader.close();
	}

	protected void loadIsoformDataByChrToTree(File file, String chr) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			if(looksLikeData(nextLine) && isChrEqual(nextLine,chr)){
				RefSeqGene gene=new RefSeqGene(nextLine, false);	    			
				this.addRefSeq(gene);
			}
		}
		reader.close();
	}

	// This will be moving to a new class that will implement an abstract class called RefSeqGeneSet
	//This code loads cufflinks GTF to a bed file parser: the saved name is the transcript_id, the score is the FPKM
	//The gene_id is the same as the transcripts id longest prefix (remove the last .x)
	//if chr is empty - load all the data, else load only the specific chr
	public void loadGTF(File file, String chr) throws IOException {
	
		Map<String, Collection> exonMap=new TreeMap();
		Map<String, String> signMap=new HashMap();
		Map<String, String> srcMap=new HashMap();
		Map<String, String> attrMap=new HashMap();
		boolean getOneChr= !(chr.equals(""));
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		
		String nextLine;
	        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	           if(nextLine.startsWith("chr")){
					try{
	        	   String[] tokens=nextLine.split("\t");
	        	   if (getOneChr && !(tokens[0].equals(chr)))
	        		   continue;
	        	   // -1 on start point because gff is 1 shifted while BED is 0 shifted.
	        	   //On the other hand, gff is inclusive while bed is not-> do not subtract 1 from the end coordinate
					Alignments align=new Alignments(tokens[0], new Integer(tokens[3])-1, new Integer(tokens[4]));
					//String name1=tokens[8].split(";")[1]; //transcripts_id - this assumes that transcript_id is in the second column - buggy
					//String name=name1.split(" ")[2];
					String name=  GFF.getGTFAttr(tokens[8], "transcript_id");
					name=name.replaceAll("\\s", "_");
					name=name.replaceAll("\"", "");
					
					String sign= tokens[6];
					String attributes =tokens[8];
					
					Collection c=new TreeSet();
					if(exonMap.containsKey(name)){c=exonMap.get(name);}
					if (tokens[2].equalsIgnoreCase("exon"))
						c.add(align);
					if (!signMap.containsKey(name)){
						signMap.put(name, sign);
						attrMap.put(name,attributes);
						srcMap.put(name, tokens[1]);
					}
					if (tokens[2].equalsIgnoreCase("transcript")){
						signMap.put(name, tokens[6]);
						attrMap.put(name,attributes);
						srcMap.put(name, tokens[1]);
					}
					exonMap.put(name, c);
					}catch(Exception ex){System.err.println(nextLine);}
	           }
	        }
	                 
	        reader.close();
	        
	        for(String transcript: exonMap.keySet()){
				Collection<Alignments> exons=exonMap.get(transcript);
				RefSeqGene gene=new RefSeqGene(exons);
				gene.setName(transcript);
				gene.setOrientation(signMap.get(transcript));
				gene.addExtraField(attrMap.get(transcript));
				gene.addAttribute("source",srcMap.get(transcript));
				this.addRefSeq(gene);
			}
	        
	
	}

	public static Map<Alignments, String> loadDataLine(File file) throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		Map<Alignments, String> alignments =  loadDataLine(reader);
		//reader.close();
		return alignments;
	}

	public static Set<Alignments> loadAlignmentData(File file) throws IOException{
		if(file==null){return null;}
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		Set<Alignments> alignments =  loadAlignmentData(reader);
		//reader.close();
		return alignments;
	}

	public static Set<Alignments> loadAlignmentData(BufferedReader reader) throws IOException{
		Set<Alignments> data=new TreeSet<Alignments>();
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	
			if(looksLikeData(nextLine)){
				String[] tokens=nextLine.split("\t");
				String chr=(tokens[0]);
				int start=new Integer(tokens[1]);
				int end=new Integer(tokens[2]);
				Alignments track=new Alignments(chr, start, end);
				if(tokens.length>3) track.setName(tokens[3]);
				if(tokens.length>5) track.setStrand(tokens[5]);
				track.setLine(nextLine);
				data.add(track);
			}
	
		}
	
	
		reader.close();
		return data;
	
	}

	public static Set<Alignments> loadAlignmentData(File file, boolean ucscDownload) throws IOException{
		if(!ucscDownload){return loadAlignmentData(file);}
	
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		Set<Alignments> data=new TreeSet<Alignments>();
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	
	
			String[] tokens=nextLine.split("\t");
			String chr=(tokens[1]);
			int start=new Integer(tokens[2]);
			int end=new Integer(tokens[3]);
	
			Alignments track=new Alignments(chr, start, end);
			track.setLine(nextLine);
			data.add(track);
	
	
		}
	
	
		reader.close();
		return data;
	
	}

	public static Map<String, Collection<Alignments>> loadAlignmentDataByChr(File file) throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	
		Map<String, Collection<Alignments>> data=new TreeMap<String, Collection<Alignments>>();
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	
			if(looksLikeData(nextLine)){
	
				String[] tokens=nextLine.split("\t");
				String chr=(tokens[0]);
				int start=new Integer(tokens[1]);
				int end=new Integer(tokens[2]);
	
				String name=null;
				if(tokens.length>3){name=tokens[3];}
				/*String strand=tokens[5];
	
				String[] blockSizes=tokens[10].split(",");
				String[] blockStarts=tokens[11].split(",");
				ArrayList[] exonStartEnd=setBlockStartsAndEnds(blockStarts, blockSizes, new Integer(tokens[9]), start);
				 */
				Alignments track=new Alignments(chr, start, end);
				track.setName(name);
				Collection<Alignments> set=new ArrayList<Alignments>();
				if(data.containsKey(chr)){set=data.get(chr);}
				track.setLine(nextLine);
				set.add(track);
				data.put(chr, set);
			}
	
		}
	
	
		reader.close();
		return data;
	
	}

	public static Map<String, IntervalTree<Alignments>> loadAlignmentDataToTree(File file) throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	
		Map<String, IntervalTree<Alignments>> data=new TreeMap<String, IntervalTree<Alignments>>();
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	
			if(looksLikeData(nextLine)){
	
				String[] tokens=nextLine.split("\t");
				String chr=(tokens[0]);
				int start=new Integer(tokens[1]);
				int end=new Integer(tokens[2]);
	
				String name="";
				if(tokens.length>3){name=tokens[3];}
				Alignments track=new Alignments(chr, start, end);
	
	
				IntervalTree<Alignments> tree=new IntervalTree();
	
				if(data.containsKey(chr)){tree=(IntervalTree)data.get(chr);}
				tree.put(start, end, track);
				data.put(chr, tree);
			}
	
		}
	
	
		reader.close();
		return data;
	
	}

	public static Map<String, List<Alignments>> loadAlignmentDataByChrArrayList(File file) throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	
		Map<String, List<Alignments>>  data=new TreeMap<String, List<Alignments>> ();
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	
			if(looksLikeData(nextLine)){
	
				String[] tokens=nextLine.split("\t");
				String chr=(tokens[0]);
				int start=new Integer(tokens[1]);
				int end=new Integer(tokens[2]);
	
				Alignments track=new Alignments(chr, start, end);
				if(tokens.length>3){track=new Alignments(chr, start, end, tokens[3]);}
	
	
				List<Alignments> list= data.get(chr);
				if(data.containsKey(chr)){
					list=new ArrayList<Alignments>();
					data.put(chr, list);
				}
				list.add(track);
			}
	
		}
	
	
		reader.close();
		return data;
	
	}

	public static Map<String, Set<Alignments>> loadAlignmentDataByChr3Prime(File file, int extensionFactor) throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	
		Map<String, Set<Alignments>> data=new TreeMap<String, Set<Alignments>>();
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	
	
			String[] tokens=nextLine.split("\t");
			String chr=(tokens[0]);
			int startOriginal=new Integer(tokens[1]);
			int endOriginal=new Integer(tokens[2]);
	
	
			int start=Math.max(0, startOriginal-extensionFactor);
			int end=endOriginal+extensionFactor;
	
			if(tokens.length>5){
				String strand=tokens[5];
				if(strand.equalsIgnoreCase("+")){start=startOriginal;}
				else if(strand.equalsIgnoreCase("-")){end=endOriginal;}
			}
	
	
			Alignments track=new Alignments(chr, start, end);
			if(tokens.length>5){track=new Alignments(chr, start, end, tokens[5]);}
			//System.err.println(track);
			Set<Alignments> set = data.get(chr);
			if(set == null){
				set= new TreeSet<Alignments>();
				data.put(chr, set);
			}
			set.add(track);
	
	
		}
	
	
		reader.close();
		return data;
	
	}

	public static Map<String, Set<Alignments>> loadAlignmentDataByChr5Prime(File file, int extensionFactor) throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	
		Map<String, Set<Alignments>> data=new TreeMap<String, Set<Alignments>>();
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	
			if(looksLikeData(nextLine)){
	
				String[] tokens=nextLine.split("\t");
				String chr=(tokens[0]);
				int startOriginal=new Integer(tokens[1]);
				int endOriginal=new Integer(tokens[2]);
	
	
				int start=startOriginal-extensionFactor;
				int end=endOriginal+extensionFactor;
	
				if(tokens.length>5){
					String strand=tokens[5];
					if(strand.equalsIgnoreCase("+")){end=endOriginal;}
					else if(strand.equalsIgnoreCase("-")){start=startOriginal;}
				}
	
	
				Alignments track=new Alignments(chr, start, end);
				if(tokens.length>5){track=new Alignments(chr, start, end, tokens[5]);}
				Set<Alignments> set = data.get(chr);
				if(set == null){
					set= new TreeSet<Alignments>();
					data.put(chr, set);
				}
				set.add(track);
			}
	
		}
	
	
		reader.close();
		return data;
	
	}

	public static Map<Alignments, String> loadDataLine(BufferedReader reader) throws IOException{
		Map<Alignments, String> data=new TreeMap<Alignments, String>();
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	
			if(nextLine.startsWith("chr")){
	
				String[] tokens=nextLine.split("\t");
				String chr=(tokens[0]);
				int start=new Integer(tokens[1]);
				int end=new Integer(tokens[2]);
	
				Alignments track=new Alignments(chr, start, end);
				if(tokens.length>5){track=new Alignments(chr, start, end, tokens[5]);}
				data.put(track, nextLine);
			}
	
		}
	
	
		reader.close();
		return data;
	
	}

	

	/////////////////////////////////////////////////////////////////////
	
	// TODO: GFF TO BED, implement mitch's version
	
	/////////////////////////////////////////////////////////////////////
	
	// TODO: GFF TO BED, implement mitch's version
	
	
	
	
	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////
	
	public static Map<Integer, Collection<PSL>> loadPSLBySeqName(File[] pslFiles) throws IOException {
	
		Map<Integer, Collection<PSL>> rtrn=new TreeMap<Integer, Collection<PSL>>();
	
		for(int i=0; i<pslFiles.length; i++){
			System.err.println(pslFiles[i]);
	
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(pslFiles[i])));
	
	
			String nextLine;
			while ((nextLine = reader.readLine()) != null) {
				try{
					PSL psl=new PSL(nextLine);              	
					Collection<PSL> list=new ArrayList();
					if(rtrn.containsKey(psl.getSeqNumber())){list=rtrn.get(psl.getSeqNumber());}
					list.add(psl);
					rtrn.put(psl.getSeqNumber(), list);
				}catch(Exception ex){System.err.println("Skipping: "+nextLine);}
	
			}
	
	
			reader.close();
	
		}
		return rtrn;
	
	
	}

	public static Map<Integer, Collection<PSL>> loadPSLBySeqName(File file) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	
		Map<Integer, Collection<PSL>> rtrn=new TreeMap<Integer, Collection<PSL>>();
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			try{
				PSL psl=new PSL(nextLine);              	
				Collection<PSL> list=new ArrayList();
				if(rtrn.containsKey(psl.getSeqNumber())){list=rtrn.get(psl.getSeqNumber());}
				list.add(psl);
				rtrn.put(psl.getSeqNumber(), list);
			}catch(Exception ex){System.err.println("Skipping: "+nextLine);}
	
		}
	
	
		reader.close();
		return rtrn;
	}

	public void addRefSeq (RefSeqGene gene) {
		String chr=gene.getChr();
		if (  this.annotationSetMap.containsKey(chr) == false) {
			IntervalTree<RefSeqGeneWithIsoforms> data = new IntervalTree<RefSeqGeneWithIsoforms>();
			this.annotationSetMap.put(gene.getChr(), data);
		}
		IntervalTree<RefSeqGeneWithIsoforms> chrTree= this.annotationSetMap.get(gene.getChr());
		if (chrTree.find(gene.getStart(), gene.getEnd())==null){
			RefSeqGeneWithIsoforms geneWithIsos = new RefSeqGeneWithIsoforms(gene);
			try{
				chrTree.put(gene.getStart(), gene.getEnd(), geneWithIsos);
			} catch(IllegalArgumentException e) {
				System.err.println(e.getMessage());
				throw new IllegalArgumentException("Gene name: " + gene.getName());
			}
			updateAnnotationNameMap(gene.getName(), geneWithIsos);
		}
		else{
			RefSeqGeneWithIsoforms geneWithIsos = chrTree.find(gene.getStart(), gene.getEnd()).getValue();
			geneWithIsos.AddIsoform(gene);
			updateAnnotationNameMap(gene.getName(), geneWithIsos);
		}
	}

	public void addRefSeq (RefSeqGeneWithIsoforms gene) {
		String chr=gene.getChr();
		if (  this.annotationSetMap.containsKey(chr) == false) {
			IntervalTree<RefSeqGeneWithIsoforms> data = new IntervalTree<RefSeqGeneWithIsoforms>();
			this.annotationSetMap.put(gene.getChr(), data);
		}
		IntervalTree<RefSeqGeneWithIsoforms> chrTree= this.annotationSetMap.get(gene.getChr());
		if (chrTree.find(gene.getStart(), gene.getEnd())==null){
			chrTree.put(gene.getStart(), gene.getEnd(), gene);
			for(RefSeqGene iso : gene.getAllIsoforms()) {
				updateAnnotationNameMap(iso.getName(), gene);
			}
		}
		else{
			RefSeqGeneWithIsoforms existingGene = chrTree.find(gene.getStart(), gene.getEnd()).getValue();
			existingGene.AddAllIsoforms(gene);
			for(RefSeqGene iso : gene.getAllIsoforms()) {
				updateAnnotationNameMap(iso.getName(), existingGene);
			}
			
		}
	}

	public void addRefSeqSet (Collection<? extends RefSeqGene> collection) {
		for (RefSeqGene g : collection ) {this.addRefSeq(g);}
	}

	public static void writeBED(String save, Collection<Alignments> regions)throws IOException{
		FileWriter writer=new FileWriter(save);
	
		for(Alignments align: regions){writer.write(align+"\n");}
	
		writer.close();
	}

	public static void writeBED(String save, Map<Alignments, double[]> scores)throws IOException{
		writeBED(save, scores, false, false);
	}
	
	public static void writeSortedBED(String save, Map<Alignments, double[]> scores)throws IOException {
		writeBED(save, scores, false, true);
	}
	
	public static void writeSortedBED(String save, Collection<Alignments> alignments) throws IOException {
		writeBED(save, alignments, true);
	}
	
	public static void writeSortedBEDWithScores(String save, Map<Alignments, double[]> scores) throws IOException {
		writeBED(save, scores, true, true);
	}

	public static void writeBEDWithScores(String save, Map<Alignments, double[]> scores)throws IOException{
		writeBED(save, scores, true, false);
	}

	public static void writeBED(String save, Collection<Alignments> alignments, boolean sort) throws IOException{
		Map<Alignments, double[]> scores = new HashMap<Alignments, double[]>();
		for (Alignments align : alignments) {
			scores.put(align, null);
		}
		writeBED(save, scores, false, sort);
	}
	
	public static void writeBED(String save, Map<Alignments, double[]> scores, boolean withScores, boolean sort) throws IOException{
		FileWriter writer=new FileWriter(save);
	
		List<Alignments> keys = new LinkedList<Alignments>(scores.keySet());
		if (sort) {
			Collections.sort(keys, new Comparator<Alignments>() {
		         @Override
		         public int compare(Alignments o1, Alignments o2) {
		             return o1.compareTo(o2);
		         }
		     });
		}
		
		for(Alignments align: keys){
			BED bed = new BED(align);
			double[] vals = scores.get(align);
			if (vals != null) {
				vals=scores.get(align);
				bed.setScore(vals[1]);
			}
			
			writer.write(bed.toShortString());
			
			if (withScores && vals != null) {
				for(double score : vals) {
					writer.write("\t");
					writer.write(String.valueOf(score));
				}
			}
			writer.write("\n");
		}
		writer.close();
	}

	public void writeFullBed(String save) throws IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(save));
		writeFullBed(bw);
		bw.close();
	}

	public void writeFullBed(String save, boolean useExtraFields) throws IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(save));
		writeFullBed(bw, useExtraFields);
		bw.close();
	}

	public void writeFullBed(BufferedWriter bw) throws IOException{
		writeFullBed(bw, true);
	}

	public void writeFullBed(BufferedWriter bw , boolean useExtraFields) throws IOException{
		Iterator<String> chrIt=this.getChromosomeIterator();
		while(chrIt.hasNext()){
			String chr=chrIt.next();
			Iterator<RefSeqGeneWithIsoforms> geneIt = getChrTree(chr).valueIterator();
			while(geneIt.hasNext()) {
				RefSeqGeneWithIsoforms gene = geneIt.next();
				Collection<RefSeqGene> isoforms = gene.getAllIsoforms();
				for(RefSeqGene isoform: isoforms){
					bw.write(isoform.toBED(useExtraFields));
					bw.newLine();
				}
			}
		}
	}

	public static void writeFullBED(String save, Collection<RefSeqGene> genes)throws IOException{
		FileWriter writer=new FileWriter(save);
	
		for(RefSeqGene gene: genes){writer.write(gene+"\n");}
	
		writer.close();
	}

	public static void writeFullBED(String save, Map<RefSeqGene, ? extends Object> genes)throws IOException{
		FileWriter writer=new FileWriter(save);
	
		for(RefSeqGene gene: genes.keySet()){writer.write(gene+"\t"+genes.get(gene)+"\n");}
	
		writer.close();
	}

	public void makeGenes() {
		makeGenes(0);
	}

	public void makeGenes(double minPctOverlapToMerge) {
		//System.err.println("Make genes");
		Map<String, IntervalTree<RefSeqGeneWithIsoforms>> mergedGenes = getMergedAnnotationMap(minPctOverlapToMerge);
		Iterator<String> chrIt=mergedGenes.keySet().iterator();
		while(chrIt.hasNext()){
			String chr=chrIt.next();
			IntervalTree<RefSeqGeneWithIsoforms> tree=mergedGenes.get(chr);
			Iterator<RefSeqGeneWithIsoforms> geneIt=tree.valueIterator();
			while(geneIt.hasNext()) {
				RefSeqGeneWithIsoforms gene = geneIt.next();
				gene.setName("gene");
				//System.err.println("Merged gene " + gene.toBED());
				gene.cleanIsoforms();
				int isoformsAdded = 0;
				Iterator<RefSeqGeneWithIsoforms> overlappers = getOverlappers(gene).valueIterator();
				while(overlappers.hasNext()) {
					RefSeqGeneWithIsoforms overlapper = overlappers.next();
					if(isOverlapCompatible(gene, overlapper, minPctOverlapToMerge)) {
					//System.err.println("Adding overlapper " + overlapper.toBED());
						Collection<RefSeqGene> isoforms = overlapper.getAllIsoforms();
						for(RefSeqGene iso: isoforms) {
							gene.setName(gene.getName()+"_"+iso.getName());
							boolean couldAdd = gene.addContainedIsoform(iso);
							if(!couldAdd) {
								//System.err.println("WARN: Could not add isoform " + overlapper.toBED() + " to " + gene.toBED());
							} else {
								isoformsAdded++;
								//System.err.println("Added isoform " + overlapper.getName() + " to " + gene.getName());
							}
						}
					} 
					else{
						//System.err.println("is NOT overlap compatible.");
					}
				}
				if(isoformsAdded ==0) {
					System.err.println("ERROR: Gene " + gene.getName() + " " + gene.toBED() + "  had no overlapping isoforms");
				}
			}
		}		
		
		annotationSetMap = mergedGenes;
	}

	public static Map<String, IntervalTree<RefSeqGene>> makeIntervalTreeFor3Prime(File geneFile, int size) throws IOException {
		Map<String,Collection<RefSeqGene>> genesByChr=loadDataByChr(geneFile);
		Map<String, IntervalTree<RefSeqGene>> rtrn=new TreeMap();
	
		for(String chr: genesByChr.keySet()){
			Collection<RefSeqGene> genes=genesByChr.get(chr);
			IntervalTree<RefSeqGene> tree=new IntervalTree();
			for(RefSeqGene gene: genes){
				Alignments end=gene.get3PrimeExon(size); //TODO should do this right later
				tree.put(end.getStart(), end.getEnd(), gene);
			}
			rtrn.put(chr, tree);
		}
	
		return rtrn;
	}

	public boolean isEmpty()  {
		boolean isEmpty = true;
		Iterator<IntervalTree<RefSeqGeneWithIsoforms>> annotationTreeIt = annotationSetMap.values().iterator();
		while(isEmpty && annotationTreeIt.hasNext()) {
			isEmpty = annotationTreeIt.next().isEmpty();
		}
		
		return isEmpty;
	}



	private static boolean looksLikeData(String nextLine) {
		return nextLine.trim().length() > 0 && ! nextLine.startsWith("#") && !nextLine.startsWith("track") && !nextLine.startsWith("browser");
	}

	public boolean contains(String geneName) {
		return nameAnnotationMap.containsKey(geneName);
	}

	public boolean containChr(String chr) {
		return this.annotationSetMap.containsKey(chr);
	}

	private boolean isChrEqual(String data, String chr) {
		String[] tokens=data.split("\t");
		return tokens[0].equalsIgnoreCase(chr);
	}

	public static boolean isOverlapCompatible(RefSeqGene first,  RefSeqGene second, double minPctOverlap) {
		// This is a bit of a hack. We want to discard antisense overlappers but also sense overlappers
		// which just overlap their ends.
		boolean compatible = first.overlaps(second);
		//System.err.println("Value of compatible flag: "+compatible);
		if(compatible) {
			//RefSeqGene overlap = first.getOverlap(second);
			double pctOverlap = Math.max(first.percentOverlapping(second), second.percentOverlapping(first));
			//System.err.println("first.percentOverlapping(second) " + first.percentOverlapping(second) + " second.percentOverlapping(first) " + second.percentOverlapping(first));
			compatible = pctOverlap > minPctOverlap;
			//Alignments overlapRegion = new Alignments(overlap.getChr(), overlap.getStart(), overlap.getEnd());
			
			//compatible = overlap.getExons().length>0 || (first.getExons().length <= 2  || (!first.get3PrimeExon().contains(overlapRegion) && !first.get5PrimeExon().contains(overlapRegion)) ) &&
			//(second.getExons().length <= 2  || (!second.get3PrimeExon().contains(overlapRegion) && !second.get5PrimeExon().contains(overlapRegion)) );
			/*if(!compatible) {
				//System.err.println("Even thouth they overlap " + first.toBED() + " and " + second.toBED() + " are not compatible their only overlap by " + pctOverlap);
			}*/
		}
		return compatible;
	}

	public void clear() {
		this.annotationSetMap.clear();
		this.nameAnnotationMap.clear();
	}

	// merge the current annotation set by combining annotations that overlap in the same orientation into one.
	public void merge(){
		this.setAnnotationSetMap( getMergedAnnotationMap(true)); 
	}

	/**
	 * Merges transcripts that either 
	 * a) Are just truncated version of a transcript in the set
	 * b) Have similar exonic structure but some exons differ by at most the given exonFudgeFactor
	 * c) Are truncated versions of a transcript in the set module exons that differ by the given exonFudgeFactor
	 * d) Merge transcripts that differ only in their 5' and/or 3' end
	 * It also removes single exon genes that overlap multiexonic transcripts.
	 * @param exonFudgeFactor
	 * 
	 */
	public void collapse(int exonFudgeFactor) {
		Iterator<String> chrIt = getChromosomeIterator();
		collapseUTRs();
		System.err.println("Collapsed UTRs");
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			Iterator<RefSeqGeneWithIsoforms> chrGeneIt = getChrTree(chr).valueIterator();
			while(chrGeneIt.hasNext()) {
				RefSeqGeneWithIsoforms gene = chrGeneIt.next();
				gene.colapse(exonFudgeFactor);
			}
		}
		System.err.println("Collapsed transcript with same start and end");
		
		filterSingleExons();
		
		
		//BEDFileParser tmpReader = new BEDFileParser();
		List<RefSeqGeneWithIsoforms> toRemove = new ArrayList<RefSeqGeneWithIsoforms>();
		chrIt = getChromosomeIterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			Iterator<RefSeqGeneWithIsoforms> chrGeneIt = getChrTree(chr).valueIterator();
			while(chrGeneIt.hasNext()) {
				RefSeqGeneWithIsoforms gene = chrGeneIt.next();
				if(toRemove.contains(gene)) {
					continue;
				}
				String orientation = gene.getOrientation();
				gene.setOrientation("*");
				Iterator<RefSeqGeneWithIsoforms> overlapIt = getOverlappers(gene).valueIterator();
				gene.setOrientation(orientation);
				while(overlapIt.hasNext()) {
					RefSeqGeneWithIsoforms overlapper = overlapIt.next();
					if(!overlapper.equals(gene) && (overlapper.almostEqual(gene, exonFudgeFactor) || gene.almostContainsStructure(overlapper, exonFudgeFactor))) {
						toRemove.add(overlapper);
					}
				}
			}
		}
		System.err.println("Removed shorter transcript with that have an exon structure already covered; going to remove  " + toRemove.size() + " transcripts");
		for(RefSeqGeneWithIsoforms g : toRemove) {
			remove(g);
		}
		
	}

	/**
	 * Homogenizes UTRs for overlapping transcripts but lengthening UTRs until all overlapping (from 5' to 3') transcripts same UTRs
	 * and transcripts that only differ in their UTR are removed. 
	 * 
	 */
	public void collapseUTRs() {
		BEDFileParser tmpReader = new BEDFileParser();
		int numMergedGenes=0;
		 Iterator<String> chrIt = getChromosomeIterator();
		 while(chrIt.hasNext()) {
			 String chr = chrIt.next();
			 Iterator<RefSeqGeneWithIsoforms> chrGeneIt = getChrTree(chr).valueIterator();
			 while(chrGeneIt.hasNext()) {
				 RefSeqGeneWithIsoforms gene = chrGeneIt.next();
				 gene.standardizeOreintation();
				 if(gene.getNumExons() == 1 || gene.isUnoriented()) {
					 tmpReader.addRefSeq(gene);
					 continue;
				 }
				 Iterator<RefSeqGeneWithIsoforms> overlapIt = tmpReader.getOverlappers(gene).valueIterator();
				 if(!overlapIt.hasNext()) {
					 tmpReader.addRefSeq(gene);
				 } else {
					 boolean mergedToAnotherGene = false;
					 while(overlapIt.hasNext() && !mergedToAnotherGene) {
						 RefSeqGeneWithIsoforms overlapper = overlapIt.next();
						 if(overlapper.getNumExons() == 1 || !overlapper.getOrientation().equals(gene.getOrientation()) ) {
							 continue;
						 }
						 //temporaryOrientTranscriptPair(gene, overlapper);
						 //System.err.println("Gene orientation " + gene.getOrientation() + " overlapper oreintation " + overlapper.getOrientation());
						 if(overlapper.get3PrimeExon().overlaps(gene.get3PrimeExon()) && overlapper.get5PrimeExon().overlaps(gene.get5PrimeExon())) {
							 mergedToAnotherGene = true;
							 //System.err.println("GENE before updating: \n" + gene.allIsoformsToBED());
							 //System.err.println("OVERLPR before updating: \n" + overlapper.allIsoformsToBED());
							 numMergedGenes++;
							 tmpReader.remove(overlapper);
							 //System.err.println("Merging\n\t" + gene.toBED() + "\n\t with\n\t " + overlapper.toBED());
							 //System.err.println("\tgene's 5' and 3' exons: " + gene.get5PrimeExon().toUCSC() + " --- " + gene.get3PrimeExon().toUCSC());
							 //System.err.println("\toverlapper's 5' and 3' exons: " + overlapper.get5PrimeExon().toUCSC() + " --- " + overlapper.get3PrimeExon().toUCSC());
							 if("+".equals(gene.getOrientation())) {
								 if(gene.getEnd() > overlapper.getEnd()) {
									 overlapper.set3PrimeEnd(gene.getEnd());
								 } else {
									 gene.set3PrimeEnd(overlapper.getEnd());
								 }
								 if(gene.getStart() < overlapper.getStart()) {
									 overlapper.set5PrimeEnd(gene.getStart());
								 } else {
									 gene.set5PrimeEnd(overlapper.getStart());
								 }
								 
							 } else if ("-".equals(gene.getOrientation())) {
								 if(gene.getEnd() > overlapper.getEnd()) {
									 overlapper.set5PrimeEnd(gene.getEnd());
								 }else {
									 gene.set5PrimeEnd(overlapper.getEnd());
								 }
								 if(gene.getStart() < overlapper.getStart()) {
									 overlapper.set3PrimeEnd(gene.getStart());
								 } else {
									 gene.set3PrimeEnd(overlapper.getStart());
								 }
								 
							 }
							 //System.err.println("GENE after updating: \n" + gene.allIsoformsToBED());
							 //System.err.println("OVERLPR after updating: \n" + overlapper.allIsoformsToBED());
							 //System.err.println("\n\n");
							 boolean result = overlapper.AddAllIsoNotIncludedforms(gene);
	
							 if(!result){								 
								 System.err.println("WARN: Could not add "+gene.toUCSC() + " to  " + overlapper.toUCSC() );
							 }
							 
							 tmpReader.addRefSeq(overlapper);
						 } 
					 }
					 if(!mergedToAnotherGene) {
						 tmpReader.addRefSeq(gene);
					 }
				 }
			 }
		 }
		 this.annotationSetMap = tmpReader.annotationSetMap;
		 this.nameAnnotationMap = tmpReader.nameAnnotationMap;
		 System.err.println("3' UTR merged "+ numMergedGenes + " genes");
	}

	/*	private void temporaryOrientTranscriptPair(RefSeqGeneWithIsoforms gene1,
			RefSeqGeneWithIsoforms gene2) {
		if(gene2.isUnoriented() && !gene1.isUnoriented()) {
			 gene2.setOrientation(gene1.getOrientation());
		 } else if (gene1.isUnoriented() && !gene2.isUnoriented()) {
			 gene1.setOrientation(gene2.getOrientation());
		 } else if (gene1.isUnoriented() && gene2.isUnoriented()) {
			 gene1.setOrientation("+");
			 gene2.setOrientation("+");
		 }
	}*/
	
	/**
	 * Merges all single exon transcripts
	 */
	public void collapseSingleExonTranscripts () {
		AnnotationReader<BasicGenomicAnnotation> singleExonSet = new BasicAnnotationReader();
		 Iterator<String> chrIt = getChromosomeIterator();
		 while(chrIt.hasNext()) {
			 String chr = chrIt.next();
			 Iterator<RefSeqGeneWithIsoforms> chrGeneIt = getChrTree(chr).valueIterator();
			 while(chrGeneIt.hasNext()) {
				 RefSeqGeneWithIsoforms gene = chrGeneIt.next();
				 Alignments [] geneExons = gene.getExons();
				 if(geneExons.length == 1) {
					 singleExonSet.addAnnotation(new BasicGenomicAnnotation(geneExons[0]));
					 remove(gene);
				 }
			 
			 }
		 }
		 
		 singleExonSet.merge();
		 chrIt = singleExonSet.getChromosomeIterator();
		 while(chrIt.hasNext()) {
			 String chr = chrIt.next();
			 Iterator<BasicGenomicAnnotation> mergedTranscriptIt = singleExonSet.getChromosomeTree(chr).valueIterator();
			 while(mergedTranscriptIt.hasNext()) {
				 RefSeqGene singleExonGene = new RefSeqGene(new Alignments(mergedTranscriptIt.next()));
				 addRefSeq(singleExonGene);
			 }
		 }
		 
	}

	public void collapseBy3pUTR() {
		BEDFileParser tmpParser = new BEDFileParser();
		//List<RefSeqGeneWithIsoforms> removed = new ArrayList<RefSeqGeneWithIsoforms>();
		int numMergedGenes = 0;
		 Iterator<String> chrIt = getChromosomeIterator();
		 while(chrIt.hasNext()) {
			 String chr = chrIt.next();
			 Iterator<RefSeqGeneWithIsoforms> chrGeneIt = getChrTree(chr).valueIterator();
			 while(chrGeneIt.hasNext()) {
				 RefSeqGeneWithIsoforms gene = chrGeneIt.next();
				 String orientation = gene.getOrientation();
				 gene.setOrientation("*");
				 Iterator<RefSeqGeneWithIsoforms> overlapIt = tmpParser.getOverlappers(gene).valueIterator();
				 gene.setOrientation(orientation);
				 if(!overlapIt.hasNext()) {
					 tmpParser.addRefSeq(gene);
				 } else {
					 boolean mergedToAnotherGene = false;
					 while(overlapIt.hasNext()) {
						 RefSeqGeneWithIsoforms overlapper = overlapIt.next();
						 if(!overlapper.equals(gene) && (overlapper.hasOverlapping3PUTR(gene) )) {
							 tmpParser.remove(overlapper);
							 if(("-".equals(gene.getOrientation()) && gene.getStart() < overlapper.getStart()) || ("+".equals(gene.getOrientation()) && gene.getEnd() > overlapper.getEnd())) {
								 overlapper.set3PrimeEnd(gene.getOrientedEnd());
							 } 
							 tmpParser.addRefSeq(overlapper);
							 mergedToAnotherGene = true;
							 //System.err.println("Merged genes " + gene.toBED() + " with " + overlapper.toBED() );
							 numMergedGenes++;
							 break;
						 } 
					 }
					 if(!mergedToAnotherGene) {
						 tmpParser.addRefSeq(gene);
					 }
				 }
			 }
		 }
		 
		 this.annotationSetMap = tmpParser.annotationSetMap;
		 this.nameAnnotationMap = tmpParser.nameAnnotationMap;
		 System.err.println("3' UTR merged "+ numMergedGenes + " genes");	
	
	}

	public void IncrementScoreIfOverlap(BEDFileParser bed,int setNumber, Boolean setToZero) {

		Iterator<String> chrIt = this.annotationSetMap.keySet().iterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			IntervalTree<RefSeqGeneWithIsoforms> oldTree = this.getChrTree(chr);
			Iterator<RefSeqGeneWithIsoforms> valueIt = new IntervalTree.ValuesIterator<RefSeqGeneWithIsoforms>(oldTree.iterator());
			while(valueIt.hasNext()) {
				RefSeqGene currentElement = valueIt.next();
				if (setToZero)
				{
					currentElement.setBedScore(0);
				}
				if (! bed.getOverlappers(currentElement).isEmpty())
				{
					currentElement.setBedScore(currentElement.getBedScore()+1);
				}
			}
		}
	}


    public void updateScrToBedScore() {

		Iterator<String> chrIt=this.getChromosomeIterator();
		while(chrIt.hasNext()){
			String chr=chrIt.next();
			IntervalTree<RefSeqGeneWithIsoforms> currTree= this.getChrTree(chr);
			Iterator <RefSeqGeneWithIsoforms> geneIt=currTree.valueIterator();
			while(geneIt.hasNext()){
				RefSeqGeneWithIsoforms gene=geneIt.next();
				gene.updateScrToBedScore();
			}

		}	

	}

	public void filterByOverlap (Collection<RefSeqGene> genes){
		//TODO: filter out all refseq genes that overlap "genes", from the current instance of this.
	}

	private void filterSingleExons() {
		List<RefSeqGeneWithIsoforms> toRemove = new ArrayList<RefSeqGeneWithIsoforms>();
		 Iterator<String> chrIt = getChromosomeIterator();
		 while(chrIt.hasNext()) {
			 String chr = chrIt.next();
			 Iterator<RefSeqGeneWithIsoforms> chrGeneIt = getChrTree(chr).valueIterator();
			 while(chrGeneIt.hasNext()) {
				 RefSeqGeneWithIsoforms gene = chrGeneIt.next();
				 gene.standardizeOreintation();
				 if(gene.getNumExons() == 1 ) {
					 Iterator<RefSeqGeneWithIsoforms> overlapIt = getOverlappers(gene).valueIterator();
					 boolean overlapsMultiexonic = false;
					 while(overlapIt.hasNext() && !overlapsMultiexonic) {
						 RefSeqGeneWithIsoforms overlapper = overlapIt.next();
						 overlapsMultiexonic = overlapper.getNumExons() > 1;
					 }
					 if(overlapsMultiexonic) {
						 toRemove.add(gene);
					 }
				 }
			 }
		 }
		 
		 for(RefSeqGeneWithIsoforms geneToKill : toRemove) {
			 remove(geneToKill);
		 }
		 
		 System.err.println("Removed " + toRemove.size() + " single exon genes");
		
	}

	

	/*	private void temporaryOrientTranscriptPair(RefSeqGeneWithIsoforms gene1,
			RefSeqGeneWithIsoforms gene2) {
		if(gene2.isUnoriented() && !gene1.isUnoriented()) {
			 gene2.setOrientation(gene1.getOrientation());
		 } else if (gene1.isUnoriented() && !gene2.isUnoriented()) {
			 gene1.setOrientation(gene2.getOrientation());
		 } else if (gene1.isUnoriented() && gene2.isUnoriented()) {
			 gene1.setOrientation("+");
			 gene2.setOrientation("+");
		 }
	}*/
	
	public void intersect(BEDFileParser bed) {
		//TODO: replace the current state with the intersection with this set.
	}

	public  void bed2gtf(BufferedWriter bw,String src)  throws IOException{
		bed2gtf(bw,src,false,false);
	}

	//If extra fields==true, all the GTF attributes will be set based on the extra fields
	public  void bed2CuffGtf(BufferedWriter bw,String src,boolean useExtraFileds) throws IOException {
		bed2gtf(bw,src,true,useExtraFileds);
		
	}

	//If extra fields==true, all the GTF attributes will be set based on the extra fields
	private  void bed2gtf(BufferedWriter bw,String src,boolean cuffForm,boolean CuffExtraFields) throws IOException {
	
		Map <String, Integer> lociNames=new HashMap <String,Integer>();
	
		Iterator<String> chrIt =this.getChromosomeIterator();
		while(chrIt.hasNext()){
			IntervalTree<RefSeqGeneWithIsoforms> chrTree=	this.getChrTree(chrIt.next());
			Iterator<RefSeqGeneWithIsoforms> gIt=chrTree.valueIterator();
			while (gIt.hasNext()){
				RefSeqGeneWithIsoforms g=gIt.next();
				Collection <RefSeqGene> geneSet= g.getAllIsoforms();
				int i=0;
				for (RefSeqGene gi : geneSet){
					String geneName =gi.getName(); 
					geneName=geneName.replaceAll("\\s", "_");
					geneName=geneName.replaceAll("\"", "");
	
					if (lociNames.containsKey(geneName)){
						int x=lociNames.get(geneName);
						lociNames.put(geneName, x+1);
						geneName=geneName+"_"+x;
					}
					else
						lociNames.put(geneName, 1);
	
					String transcriptName= geneName +"."+i;
					String str;
					if (cuffForm){
						String attr=gi.getExtraFields(0);
						if (! CuffExtraFields) 
							attr="";
						str=gi.toCufflinksGTF(src,geneName,transcriptName,attr);
					}
					else
						str=gi.toGTF(src,geneName,transcriptName);
					bw.write(str);
					i++;	
				}
			}
		}
	
	
	}

	//Supports cufflinks
	public void bed2gtf_bundleOverlappers(BufferedWriter bw, String src,BEDFileParser mergedData) throws IOException {
	
		Iterator<String> chrIt = mergedData.getChromosomeIterator();
		Map <String, Integer> lociNames=new HashMap <String,Integer>();
	
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			IntervalTree<RefSeqGeneWithIsoforms> mergedTree=mergedData.getChrTree(chr);
			Iterator <RefSeqGeneWithIsoforms> geneIt = mergedTree.valueIterator();
			//for each merged transcript get all of its overlapping isoforms and update their name 
			while(geneIt.hasNext()) { 
				RefSeqGene mergedTranscript=geneIt.next();
				String lociName= mergedTranscript.getName();
				String tmp="\"";
	
				lociName=lociName.replaceAll(" ", "_");
				lociName=lociName.replaceAll("\"", "");
				if (lociNames.containsKey(lociName)){
					int x=lociNames.get(lociName);
					lociNames.put(lociName, x+1);
					lociName=lociName+"_"+x;
				}
				else
					lociNames.put(lociName, 1);
				IntervalTree<RefSeqGeneWithIsoforms> overlappers=this.getOverlappers(mergedTranscript);
				Iterator<RefSeqGeneWithIsoforms> it = overlappers.valueIterator();
				int i=0;
				while (it.hasNext()){
					RefSeqGeneWithIsoforms g=it.next();
					//g.setName(lociName);
					for(RefSeqGene gi : g.getAllIsoforms()){
						String geneName =lociName;//gi.getName();
						geneName=geneName.replaceAll("\\s", "_");
						geneName=geneName.replaceAll("\"", "");
						String transcriptName= lociName +"."+i; //geneName +"."+i;
						bw.write(gi.toCufflinksGTF(src,geneName,transcriptName,""));
						i++;	
	
					}
				}
			}
		}
	
	
	}

	/////////////////////////////////////////////////////////////////////
	
	// TODO: GFF TO BED, implement mitch's version
	
	
	
	
	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////
	
	public static void convertPSLToBED(String file, String save) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		FileWriter writer=new FileWriter(save);
	
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			try{
				PSL psl=new PSL(nextLine);              	
				//Collection<PSL> list=new ArrayList();
				//if(rtrn.containsKey(psl.getSeqNumber())){list=rtrn.get(psl.getSeqNumber());}
				//list.add(psl);
				//rtrn.put(psl.getSeqNumber(), list);
				writer.write(psl.toGene()+"\n");
			}catch(Exception ex){System.err.println("Skipping: "+nextLine);}
	
		}
	
	
		reader.close();
		writer.close();
	}

	public void dedup() {
		for(String chr : annotationSetMap.keySet()) {
			Iterator<RefSeqGeneWithIsoforms> geneIt = annotationSetMap.get(chr).valueIterator();
			while(geneIt.hasNext()) {
				RefSeqGeneWithIsoforms gene = geneIt.next();
				gene.dedup();
			}
		}
	
	}

	/**
	 * Removes from this reader all elements that overlap the given list.
	 * @param recursive - if true not only transcripts overlapping the given set but also transcripts that overlap the overlappers will be removed
	 */
	public void minus(Collection<RefSeqGene> genes, boolean recursive) { 
		for(RefSeqGene gene : genes) {
			IntervalTree<RefSeqGeneWithIsoforms> overlapTree = getOverlappers(gene);
			Iterator<RefSeqGeneWithIsoforms> overlapIt = overlapTree.valueIterator();
			while(overlapIt.hasNext()) {
				RefSeqGeneWithIsoforms overlapper = overlapIt.next();
				remove(overlapper);
				if(recursive) {
					IntervalTree<RefSeqGeneWithIsoforms> overlapOverlapTree = getOverlappers(overlapper);
					Iterator<RefSeqGeneWithIsoforms> overlapOverlapperIt = overlapOverlapTree.valueIterator();
					while(overlapOverlapperIt.hasNext()) {
						RefSeqGeneWithIsoforms overlapOverlapper = overlapOverlapperIt.next();
						remove(overlapOverlapper);
					}
				}
			}
		}
	}


	public static void main(String[] args)throws IOException{
	
		Collection<RefSeqGene> genes=loadData(new File(args[0]));
		String save=args[1];
		FileWriter writer=new FileWriter(save+".plus");
		for(RefSeqGene gene: genes){
			if(gene.getOrientation().equalsIgnoreCase("+")){writer.write(gene+"\n");}
		}
		writer.close();
	
		writer=new FileWriter(save+".minus");
		for(RefSeqGene gene: genes){
			if(gene.getOrientation().equalsIgnoreCase("-")){writer.write(gene+"\n");}
		}
		writer.close();
	
	}

	//This function enables one to upload only exons that are in ThreePrimeTrimSize from the 3 prime
	private void trimAndLoadIsoformDataByChrToTree(File file,int ThreePrimeTrimSize) throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	
		String nextLine;
		int i=0;
		int trimCntr=0;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	
			if(looksLikeData(nextLine)){
	
				String[] tokens=nextLine.split("\t");
				String chr=(tokens[0]);
				int start=new Integer(tokens[1]);
				int end=new Integer(tokens[2]);
				RefSeqGene gene;
				gene=new RefSeqGene(new Alignments(chr, start, end));
				if ((end-start) >  ThreePrimeTrimSize && tokens.length>=6 ){
					if(tokens[5].equalsIgnoreCase("-"))
						gene=new RefSeqGene(new Alignments(chr, start, start+ThreePrimeTrimSize));
					else
						gene=new RefSeqGene(new Alignments(chr, end-ThreePrimeTrimSize, end));
				}	
	
				if(tokens.length>9){	
					String name=tokens[3];
					double score=new Double(tokens[4]);
					String strand=tokens[5];
	
					int blockStart=new Integer(tokens[6]);
					int blockEnd=new Integer(tokens[7]);
	
					String[] blockSizes=tokens[10].split(",");
					String[] blockStarts=tokens[11].split(",");
					List<Integer>[] exonStartEnd=setBlockStartsAndEnds(blockStarts, blockSizes, new Integer(tokens[9]), start);
	
					String [] extraColumns = null;
					if(tokens.length > 12) {
						extraColumns = new String[tokens.length - 12];
						for(int j = 12; j < tokens.length; j++) {
							extraColumns[j - 12] = tokens[j];
						}
					}
					if ((end-start) <=  ThreePrimeTrimSize)
						gene=new RefSeqGene(chr, start, end, name,score, strand, exonStartEnd[0], exonStartEnd[1], blockStart, blockEnd, extraColumns);
	
	
					else{
						trimCntr++;
						if (strand.equalsIgnoreCase("-"))
							gene=trimMinusStrand(chr,start,end,name,score,strand,exonStartEnd[0],exonStartEnd[1],blockSizes,ThreePrimeTrimSize);
						else
							gene=trimPlusStrand(chr,start,end,name,score,strand,exonStartEnd[0],exonStartEnd[1],blockSizes,ThreePrimeTrimSize);
					}
				}
	
				this.addRefSeq(gene);
	
			}
			i++;
			if(i%10000==0){logger.debug(i);}
		}
	
		logger.debug("Trimmed by 3 prime :" + trimCntr);
		reader.close();
	
	}

	//trim gene such that the transcript only includes "trimSize" from the 3' end. will break an exon in the middle
	private RefSeqGene trimPlusStrand(String chr, int start, int end, String name, double score, 
			String strand, List<Integer> blockStarts, List<Integer> blockEnds, String[] blockSizes,int trimSize) {
		int sum=0;
		int lastIndex=blockSizes.length -1;
		int prvSum;
		for (int i=blockSizes.length -1; i >= 0; i--){
			lastIndex=i;
			prvSum=sum;
			sum+=Integer.valueOf(blockSizes[i]);
			if (sum > trimSize ){
	
				blockStarts.set(i,blockEnds.get(i)-(trimSize-prvSum));
				blockSizes[i]=String.valueOf(blockEnds.get(i)-blockStarts.get(i));
	
				break;
			}
	
		}
		List<Integer> newBlockStarts=new ArrayList<Integer>();
		List<Integer> newBlockEnds=new ArrayList<Integer>();
	
		for (int i=lastIndex; i<blockSizes.length; i++){
			newBlockStarts.add(blockStarts.get(i));
			newBlockEnds.add(blockEnds.get(i));
		}
		int newStart=newBlockStarts.get(0);
		int newEnd=newBlockEnds.get(newBlockEnds.size()-1);
	
		RefSeqGene gene=new RefSeqGene(chr, newStart, newEnd, name,score, strand, newBlockStarts, newBlockEnds, newStart, newEnd);
	
		return gene;
	}

	//trim gene such that the transcript only includes "trimSize" from the 3' end. will break an exon in the middle		
	private RefSeqGene trimMinusStrand(String chr, int start, int end, String name, double score, 
			String strand, List<Integer> blockStarts, List<Integer> blockEnds, String[] blockSizes,int trimSize) {
		int sum=0;
		int lastIndex=0;
		int prvSum;
		for (int i=0; i < blockSizes.length; i++){
			lastIndex=i;
			prvSum=sum;
			sum+=Integer.valueOf(blockSizes[i]);
			if (sum > trimSize ){
				blockEnds.set(i,blockStarts.get(i)+(trimSize-prvSum));
				blockSizes[i]=String.valueOf(blockEnds.get(i)-blockStarts.get(i));
				break;
			}
		}
		List<Integer> newBlockStarts=new ArrayList<Integer>();
		List<Integer> newBlockEnds=new ArrayList<Integer>();
	
		for (int i=0; i<=lastIndex; i++){
			newBlockStarts.add(blockStarts.get(i));
			newBlockEnds.add(blockEnds.get(i));
		}
		int newStart=newBlockStarts.get(0);
		int newEnd=newBlockEnds.get(newBlockEnds.size()-1);
	
		RefSeqGene gene=new RefSeqGene(chr, newStart, newEnd, name,score, strand, newBlockStarts, newBlockEnds, newStart, newEnd);
	
		return gene;	
	}

	public RefSeqGeneWithIsoforms remove(RefSeqGeneWithIsoforms record) {
		IntervalTree<RefSeqGeneWithIsoforms> tree = getChrTree(record.getChr());
		RefSeqGeneWithIsoforms removed = null;
		if(tree != null) {
			removed = tree.remove(record.getStart(), record.getEnd());
			if(removed != null) {
				Iterator<RefSeqGene> overlapperIsos = removed.getAllIsoforms().iterator();
				while(overlapperIsos.hasNext()) {
					RefSeqGene g = overlapperIsos.next();
					this.nameAnnotationMap.remove(g.getName());
				}
			}
		}
		return removed;
	}

	public Map<String, Collection<RefSeqGene>> toMap() {
		Map<String, Collection<RefSeqGene>> rtrn=new LinkedHashMap<String, Collection<RefSeqGene>>();
		Iterator<String> chrIt=this.getChromosomeIterator();
		while(chrIt.hasNext()){
			String chr=chrIt.next();
			IntervalTree<RefSeqGeneWithIsoforms> tree=this.getChrTree(chr);
			Iterator<RefSeqGeneWithIsoforms> geneIt=tree.valueIterator();
			Collection<RefSeqGene> chrGenes = new ArrayList<RefSeqGene>(tree.size());
			rtrn.put(chr, chrGenes);
			while(geneIt.hasNext())
			{
				chrGenes.addAll(geneIt.next().getAllIsoforms());
			}
		}
		return rtrn;
	}

	public Map<String, Collection<RefSeqGene>> toConstituentIsoformMap() {
		return toConstituentIsoformMap(false);
	}

	public Map<String, Collection<RefSeqGene>> toConstituentIsoformMap(boolean removeGenesWithoutConstituentIsoforms) {
		Map<String, Collection<RefSeqGene>> rtrn=new LinkedHashMap<String, Collection<RefSeqGene>>();
		Iterator<String> chrIt=this.getChromosomeIterator();
		while(chrIt.hasNext()){
			String chr=chrIt.next();
			IntervalTree<RefSeqGeneWithIsoforms> tree=this.getChrTree(chr);
			Iterator<RefSeqGeneWithIsoforms> geneIt=tree.valueIterator();
			Collection<RefSeqGene> chrGenes = new ArrayList<RefSeqGene>(tree.size());
			rtrn.put(chr, chrGenes);
			while(geneIt.hasNext())
			{
				RefSeqGeneWithIsoforms gene = geneIt.next();
				RefSeqGene  constituentIsoform = gene.constituentIsoform();
				if(constituentIsoform != null) {
					chrGenes.add(constituentIsoform);
				} else if(!removeGenesWithoutConstituentIsoforms){
					logger.warn("Could not find constituent isoform for gene " + gene.getName() + " adding original gene ");
					chrGenes.add(gene);
				}
			}
			logger.debug("Processed " + chr + " added " + chrGenes.size());
		}
		return rtrn;
	}

	/** 
	 * Builds an artificial map of isoforms where "constituent" introns are the exons.
	 * @return
	 */
	public Map<String, Collection<RefSeqGene>> toConstituentIntroformMap() {
		Map<String, Collection<RefSeqGene>> rtrn=new LinkedHashMap<String, Collection<RefSeqGene>>();
		Iterator<String> chrIt=this.getChromosomeIterator();
		while(chrIt.hasNext()){
			String chr=chrIt.next();
			IntervalTree<RefSeqGeneWithIsoforms> tree=this.getChrTree(chr);
			Iterator<RefSeqGeneWithIsoforms> geneIt=tree.valueIterator();
			Collection<RefSeqGene> chrGenes = new ArrayList<RefSeqGene>(tree.size());
			rtrn.put(chr, chrGenes);
			while(geneIt.hasNext())
			{
				chrGenes.add(geneIt.next().constituentIntrons());
			}
		}
		return rtrn;
	}

	/*	private void temporaryOrientTranscriptPair(RefSeqGeneWithIsoforms gene1,
			RefSeqGeneWithIsoforms gene2) {
		if(gene2.isUnoriented() && !gene1.isUnoriented()) {
			 gene2.setOrientation(gene1.getOrientation());
		 } else if (gene1.isUnoriented() && !gene2.isUnoriented()) {
			 gene1.setOrientation(gene2.getOrientation());
		 } else if (gene1.isUnoriented() && gene2.isUnoriented()) {
			 gene1.setOrientation("+");
			 gene2.setOrientation("+");
		 }
	}*/
	
	public BEDFileParser copy() {
		BEDFileParser copy = new BEDFileParser();
		Iterator<String> chrIt = getChromosomeIterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			Iterator<RefSeqGeneWithIsoforms> geneIt = getChrTree(chr).valueIterator();
			while(geneIt.hasNext()) {
				copy.addRefSeq(geneIt.next());
			}
		}
		
		return copy;
	}

	// Fulfills 2 tasks:
	//Filter: filter out all the overlapping genes. result will include all
	//gene in "this" that DID NOT overlap the other set
	//Stay: result will include all genes in "this" that DID overlap the other set
	public BEDFileParser overlap_GenomicLevel(BEDFileParser other, String task,String considerOrientaion) {
	
		BEDFileParser resBed=new BEDFileParser();
		Iterator<String> chrIt=this.getChromosomeIterator();
		boolean doesOverlap=false;
		while(chrIt.hasNext()) {
			String chr=chrIt.next();
			IntervalTree<RefSeqGeneWithIsoforms> otherTree=other.getChrTree(chr);
			if(otherTree==null  ){
				if (task.equalsIgnoreCase("Filter"))
					resBed.addRefSeqSet(this.getChrRefSeqGenes(chr));
			}
			else{
				for (RefSeqGene gene: this.getChrRefSeqGenes(chr) ){
					RefSeqGeneWithIsoforms gene2=new RefSeqGeneWithIsoforms(gene);
	
					if (considerOrientaion.equalsIgnoreCase("true"))
						doesOverlap= gene2.overlapsByGenomicRegion(otherTree,true);
					else
						doesOverlap= gene2.overlapsByGenomicRegion(otherTree,false);
	
					if(doesOverlap==false && task.equalsIgnoreCase("Filter") ){
						resBed.addRefSeq(gene);
					}
					if(doesOverlap==true && task.equalsIgnoreCase("Stay") ){
						resBed.addRefSeq(gene);
					}
				}
	
			}
		}
		return resBed;
	
	}

	public BEDFileParser overlapByGenomicRegion(BEDFileParser other,String considerOrientaion) {

		return this.overlap_GenomicLevel(other,"Stay",considerOrientaion);
	}



	public  BEDFileParser mergeByGtfGeneId (){
		HashMap <String,Locus> geneMap = new HashMap <String,Locus>();
		BEDFileParser outBed=new BEDFileParser();
	
		for (RefSeqGene iso: GetGenes()){
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
		//Bug fix: if a gene id is merged across several sperated fragments al will be included with the same gene_id
		for (String geneName: geneMap.keySet()){
			List<RefSeqGene> genes= geneMap.get(geneName).getAllMerged();
			for (RefSeqGene gene : genes){
				gene.setName(geneName);
				gene.addAttribute("gene_id",geneName);
				outBed.addRefSeq(gene);	
			}
		}
	return outBed;
	}

	/*	private void temporaryOrientTranscriptPair(RefSeqGeneWithIsoforms gene1,
			RefSeqGeneWithIsoforms gene2) {
		if(gene2.isUnoriented() && !gene1.isUnoriented()) {
			 gene2.setOrientation(gene1.getOrientation());
		 } else if (gene1.isUnoriented() && !gene2.isUnoriented()) {
			 gene1.setOrientation(gene2.getOrientation());
		 } else if (gene1.isUnoriented() && gene2.isUnoriented()) {
			 gene1.setOrientation("+");
			 gene2.setOrientation("+");
		 }
	}*/
	
	
	public BEDFileParser equalizeTranscriptEnds(BEDFileParser other) {
		BEDFileParser result = new BEDFileParser();
		Iterator<String> chrIt = getChromosomeIterator();
		
		//FOR EACH CHROMOSOME
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			Iterator<RefSeqGeneWithIsoforms> geneIt = getChrTree(chr).valueIterator();
			List<RefSeqGene> addedNewIsoforms = new ArrayList<RefSeqGene>();
			
			//FOR EACH ISOFORM ON THIS CHROMOSOME
			while(geneIt.hasNext()) {
				//GENE UNDER CONSIDERATION
				Iterator<RefSeqGene> isoform_iter = geneIt.next().getAllIsoforms().iterator();
				while(isoform_iter.hasNext()){
						
					RefSeqGene gene = isoform_iter.next();
					int numExons = gene.getNumExons();
					
					boolean added = false;
					
					//ITERATOR OVER ALL OVERLAPPING GENES WITH CURRENT GENE IN THE "OTHER" BED FILE 
					Iterator<RefSeqGeneWithIsoforms> overlapIt = other.getOverlappers(gene).valueIterator();
					
					//FOR EVERY OVERLAPPING GENE IN THE OTHER BED FILE
					while(overlapIt.hasNext()) {
						RefSeqGeneWithIsoforms overlapper = overlapIt.next();
						// IF GENE AND  OVERLAPPING ISOFORM HAVE <=1 EXON OR IF GENE IS UNORIENTED(?)
							//ADD OVERLAPPING GENE TO THE RESULT LIST
						if(numExons <= 1 && overlapper.getNumExons() <= 1 || gene.isUnoriented()) {
							result.addRefSeq(overlapper);
							added = true;
						} 
						
						// ELSE IF GENE HAS MORE THAN 1 EXON
						else if(numExons > 1) {
							// GET ALL ISOFORMS OF OVERLAPPER
							Collection<RefSeqGene> isoforms = overlapper.getAllIsoforms();
							
							//GET THE LAST EXON OF THE GENE 
							Alignments gene3pUTR = gene.get3PrimeExon();
							Alignments gene5pUTR = gene.get5PrimeExon();
							
							//FOR EVERY OVERLAPPING ISOFORM
							for(RefSeqGene isoform: isoforms){
								
								// IF ISOFORM HAS MORE THAN 1 EXON
								if(isoform.getNumExons() > 1) {
									//???? OVERLAPPING INTRONS BETWEEN GENE AND OVERLAPPING ISOFORM???
									int numCompatibleIntrons = gene.numOfCompatibleIntrons(overlapper);
									//IF GENE AND OVERLAPPER HAVE SAME ORIENTATION && #COMPATIBLE INTRONS>HALF OF #INTRONS IN GENE
									if(gene.getOrientation().equals(isoform.getOrientation()) && numCompatibleIntrons > Math.floor((numExons - 1)/2.0)) {//TODO: can this be a bit more principled?
										//GET LAST EXON OF OVERLAPPER
										Alignments threePUTR = isoform.get3PrimeExon();
										Alignments fivePUTR = isoform.get5PrimeExon();
										
										//IF NO ISOFORMS ARE ADDED TO RESULT && GENE AND OVERLAPPER HAVE SAME LAST EXON
										if(!added && (threePUTR.overlaps(gene3pUTR)&& fivePUTR.overlaps(gene5pUTR))) {	
											if(threePUTR.overlaps(gene3pUTR)){
												//3'
												//TODO IS THE ASSUMPTION CORRECT	
												System.out.println(gene.getName()+" 3'End is being Set to "+isoform.toUCSC());
												gene.set3PrimeEndForGeneOnly(( gene.getOrientation().equals("+") ? threePUTR.getEnd() : threePUTR.getStart())); //This can happen many times, but we assume all isoform with overlapping 3' UTRS have same ends.
											}
											if(fivePUTR.overlaps(gene5pUTR)){
												//5'
												System.out.println(gene.getName()+" 5'End is being Set to "+isoform.toUCSC());
												gene.set5PrimeEndForGeneOnly(gene.getOrientation().equals("+") ? fivePUTR.getStart() : fivePUTR.getEnd());
											}
											
											result.addRefSeq(gene);
											added = true;
										} 
										//IF 
										else if((!gene.overlapsExon(threePUTR) || !gene.overlapsExon(fivePUTR))  && !addedNewIsoforms.contains(isoform)) {
											System.err.println("Entered for "+gene.getName());
											RefSeqGene altUTRIso = new RefSeqGene(isoform);
											if(!gene.overlapsExon(threePUTR))
												if(!gene.overlapsExon(fivePUTR))
													altUTRIso.setName(gene.getName()+"_alt3p5p_"+isoform.getOrientedEnd());
												else
													altUTRIso.setName(gene.getName()+"_alt3p_"+isoform.getOrientedEnd());
											else
												altUTRIso.setName(gene.getName()+"_alt5p_"+isoform.getOrientedEnd());
											result.addRefSeq(altUTRIso);
											addedNewIsoforms.add(isoform);
										}
									}
								}
							}
		
						}
					}
					//IF NO OVERLAPPING ISOFORMS WERE ADDED,
						//ADD THE GENE TO THE RESULT LIST
					if(!added) {
						System.out.println("Nothing added for "+gene.getName());
						result.addRefSeq(gene);
					}
				}
			}
			
			Iterator<RefSeqGeneWithIsoforms> geneIt2 = getChrTree(chr).valueIterator();
			while(geneIt2.hasNext()) {
				Iterator<RefSeqGene> isoform_iter2 = geneIt2.next().getAllIsoforms().iterator();
				while(isoform_iter2.hasNext())
					System.out.println(isoform_iter2.next().toBED());
			}

		}
		
		return result;
	}

	public BEDFileParser differenceByGenomicRegion(BEDFileParser other,String considerOrientaion) {

		return this.overlap_GenomicLevel(other,"Filter",considerOrientaion);
	}

	

	
	/////////////////////////////////////////////////////////////////////

	// TODO: GFF TO BED, implement mitch's version

	


	

	///////////////////////////////////////////////////////////////////////////////////////////////////

	public static Map<String, BEDFileParser> uploadSetOfFiles(String bedLst) throws IOException {

		Map<String, BEDFileParser> res= new HashMap<String, BEDFileParser>();
		BufferedReader br = new BufferedReader(new FileReader(bedLst));
		String line;

		while((line = br.readLine()) != null) {
			line = line.trim();
			String [] lineSplit = line.split("\t");
			BEDFileParser set = new BEDFileParser(lineSplit[0]);
			res.put(lineSplit[1], set);

		}		


		return res;
	}


	public static BEDFileParser expanedUtrs(BEDFileParser bed,Integer utr1, Integer utr2) {

		BEDFileParser newBed= new BEDFileParser();
		for (RefSeqGene g:bed.GetGenes()){
			RefSeqGene newG=g.copy();
			newG.expandUtrs(utr1, utr2);
			newBed.addRefSeq(newG);
		}

		return newBed;

	}

	public static BEDFileParser readSiphyOutToBed(String siphyFile) throws NumberFormatException, IOException {
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(siphyFile)));
		BEDFileParser out= new BEDFileParser();
		int i = 0;
    	String nextLine;
    	while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
    		    String[] tokens=nextLine.split("\t");
    			RefSeqGene gene=new RefSeqGene(tokens[0],Integer.valueOf(tokens[1]),Integer.valueOf(tokens[2]));
    			gene.addExtraField(tokens[3]);
    			gene.addExtraField(tokens[4]);
    			gene.addExtraField(tokens[5]);
    			gene.addExtraField(tokens[6]);
    			out.addRefSeq(gene);
    			i++;
        		if(i%10000==0){System.err.println(i);}
    	}
    	reader.close();
		
		return out;
	}

	public BEDFileParser removeSmallTranscripts(int minExnNum, int minTranscriptLength) {

		BEDFileParser res=new BEDFileParser();
		
		for (RefSeqGene g: this.GetGenes()){
			if (g.getNumExons()>=minExnNum || g.getTranscriptLength()>=minTranscriptLength ){
				res.addRefSeq(g);
			}
		}
		
		return res;
	}

	public  BEDFileParser removeSingeleExon() {
		BEDFileParser res=new BEDFileParser();
		for (RefSeqGene g: this.GetGenes()){
			if (g.getNumExons()>=2 )
				res.addRefSeq(g);
		}
		return res;
	}

	public RefSeqGeneWithIsoforms getClosestUpstream(RefSeqGene annot) {
		IntervalTree<RefSeqGeneWithIsoforms> chrTree = getChrTree(annot.getChr());
		RefSeqGeneWithIsoforms closestUpstream = null;
		if(chrTree != null) {
			Node<RefSeqGeneWithIsoforms> max  = chrTree.max(annot.getStart()-1, annot.getStart() );
			if(max != null) {
				closestUpstream = max.getValue();
			}
		}
		return closestUpstream;
		
	}
	
	public RefSeqGeneWithIsoforms getClosestDownstream(RefSeqGene annot) {
		IntervalTree<RefSeqGeneWithIsoforms> chrTree = getChrTree(annot.getChr());
		RefSeqGeneWithIsoforms closestDownstream = null;
		if(chrTree != null) {
			Node<RefSeqGeneWithIsoforms> min  = chrTree.min(annot.getEnd() + 1, annot.getEnd() + 2);
			if(min != null) {
				closestDownstream = min.getValue();
			}
		}
		return closestDownstream;
		
	}

	public RefSeqGeneWithIsoforms getClosestUpstream(LightweightGenomicAnnotation annot) {
		IntervalTree<RefSeqGeneWithIsoforms> chrTree = getChrTree(annot.getChromosome());
		RefSeqGeneWithIsoforms closestUpstream = null;
		if(chrTree != null) {
			Node<RefSeqGeneWithIsoforms> max  = chrTree.max(annot.getStart()-1, annot.getStart() );
			if(max != null) {
				closestUpstream = max.getValue();
			}
		}
		return closestUpstream;
		
	}
	


	public RefSeqGeneWithIsoforms getClosestDownstream(LightweightGenomicAnnotation annot) {
		IntervalTree<RefSeqGeneWithIsoforms> chrTree = getChrTree(annot.getChromosome());
		RefSeqGeneWithIsoforms closestDownstream = null;
		if(chrTree != null) {
			Node<RefSeqGeneWithIsoforms> min  = chrTree.min(annot.getEnd() + 1, annot.getEnd() + 2);
			if(min != null) {
				closestDownstream = min.getValue();
			}
		}
		return closestDownstream;
		
	}
	
	/**
	 * This method will create a new set with permuted features of the same length
	 * @param features
	 * @param chrSizes 
	 * @return
	 */
	public static Collection<Alignments> permute(Collection<Alignments> features, Map<String, Integer> chrSizes) {
		Collection<Alignments> rtrn=new ArrayList<Alignments>();
		//iterate through each feature and make a new feature of the same size placed on the same chromosome in random locations
		for(Alignments feature: features){
			Alignments permuted=permute(feature, chrSizes);
			rtrn.add(permuted);
		}
		return rtrn;
	}

	public static Alignments permute(Alignments feature, Map<String, Integer> chrSizes) {
		//create a new feature with the same chromosome and length but random start
		String chr=feature.getChr();
		int chrLength=chrSizes.get(chr);
		int width=feature.getSize();
		
		//pick a random location between 0 and chrLength-width
		int start = 0 + new Double(Math.random()*(chrLength-width)).intValue();
		int end = start + width;
		
		Alignments rtrn=new Alignments(chr, start, end);
		return rtrn;
	}

	

/*	private void temporaryOrientTranscriptPair(RefSeqGeneWithIsoforms gene1,
			RefSeqGeneWithIsoforms gene2) {
		if(gene2.isUnoriented() && !gene1.isUnoriented()) {
			 gene2.setOrientation(gene1.getOrientation());
		 } else if (gene1.isUnoriented() && !gene2.isUnoriented()) {
			 gene1.setOrientation(gene2.getOrientation());
		 } else if (gene1.isUnoriented() && gene2.isUnoriented()) {
			 gene1.setOrientation("+");
			 gene2.setOrientation("+");
		 }
	}*/
	

	
    

}
