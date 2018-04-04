package broad.core.gene;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import broad.core.annotation.AnnotationHandler;
import broad.core.annotation.AnnotationReader;
import broad.core.annotation.GFF;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.GenomicAnnotationFilter;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.error.ParseException;
import broad.core.gene.GeneAnnotation.Exon;

public abstract class GeneAnnotationReader extends AnnotationReader<GeneAnnotation>{
	private HashMap<String,GeneAnnotation> geneList;
	//private HashMap<String, List<GeneAnnotation>> chromosomeGeneMap = new HashMap<String, List<GeneAnnotation>>();
	//private HashMap<String, IntervalTree<GeneAnnotation>> chromosomeGeneMap = new HashMap<String, IntervalTree<GeneAnnotation>>();
	
	public GeneAnnotationReader() {
		super();
		AnnotationSet<GeneAnnotation> trivialSet = new AnnotationSet<GeneAnnotation>();
		addAnnotationSet(trivialSet);
	}
	
	/**
	 * Loads gene models from exon gff file
	 * @param fileName input gff file
	 * @param chr Only loads those that lay in given chr or all if chr is null.
	 * @throws IOException
	 */
	public void loadExonGffFile(String fileName, String chr) throws IOException {
		geneList = new HashMap<String,GeneAnnotation>();
		File source = new File(fileName);
		BufferedReader br = new BufferedReader(new FileReader(source));
		String line;
		while((line = br.readLine()) != null) {
			if(line.startsWith("#")){
				continue;
			}
			//System.out.println(line);
			GFF exon = new GFF(line.split("\t"));
			//System.out.println("Got exon " + exon.getName() + " " + exon.getLocationString());
			
			GeneAnnotation gene = geneList.get(exon.getName());
			if(gene == null) {
				
				gene = new GeneAnnotation(exon.getName());
				gene.setChromosome(exon.getChromosome());
				gene.setStart(exon.getStart());
				gene.setCdsStart(exon.getStart());
				gene.setEnd(exon.getEnd());
				gene.setCdsEnd(exon.getEnd());
				gene.setReversedOrientation(exon.inReversedOrientation());
				geneList.put(exon.getName(), gene);
				//System.out.println("Adding new gene " + gene.getName() + " " + gene.getLocationString());
				if(chr == null || chr.equals(gene.getChromosome())) {
					addToChromsomeMap(gene);
				} 
				
			}
			if(exon.getStart() < gene.getStart()) {
				gene.setStart(exon.getStart());
				gene.setCdsStart(exon.getStart());
			}
			
			if(exon.getEnd() > gene.getEnd()) {
				gene.setEnd(exon.getEnd());
				gene.setCdsEnd(exon.getEnd());
			}
			gene.addExon(exon);
			
		}

		//System.out.print("loaded " + geneList.values().size() + " Closing "+fileName );
		/*Iterator<String> chrIt = chromosomeGeneMap.keySet().iterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			System.out.println("\tchr" + chr + " loaded " + chromosomeGeneMap.get(chr).size() + " genes");
		}
		*/
		br.close();
		//System.out.print(" ..... Closed\n");
	}
	
	public void  loadChromosomeGeneAnnotations(String fileName, String chr, boolean loadCDSOnly, String version) 
	throws FileNotFoundException, ParseException, IOException{
		geneList = new HashMap<String,GeneAnnotation>();
		File source = new File(fileName);
		BufferedReader br = new BufferedReader(new FileReader(source));
		String line;
		while((line = br.readLine()) != null) {
			if(line.startsWith("#")){
				continue;
			}
			//System.out.println(line);
			GeneAnnotation entry = new GeneAnnotation(line.split("\t"), loadCDSOnly, version);
			if (chr.equals(entry.getChromosome())){
				geneList.put(entry.getName(),entry);
				addToChromsomeMap(entry);
			}
		}

		//System.out.print("Closing "+fileName);
		br.close();
		//System.out.print(" ..... Closed\n");
	}
	
	public void  loadChromosomeGeneAnnotations(String fileName, String chr, boolean loadCDSOnly) 
	throws FileNotFoundException, ParseException, IOException{
		loadChromosomeGeneAnnotations(fileName, chr, loadCDSOnly, null);
	}
	
	public  void load(String fileName, boolean loadCDSOnly, String  version, GenomicAnnotationFilter<GeneAnnotation> filter) 
		throws FileNotFoundException, ParseException, IOException{
			geneList = new HashMap<String,GeneAnnotation>();
			File source = new File(fileName);
			BufferedReader br = new BufferedReader(new FileReader(source));
			String line;
			while((line = br.readLine()) != null) {
				if(line.startsWith("#")){
					continue;
				}
				//System.out.println(line);
				GeneAnnotation entry = new GeneAnnotation(line.split("\t"), loadCDSOnly, version);
				if(filter == null || filter.accept(entry)) {
					geneList.put(entry.getName(),entry);
					addToChromsomeMap(entry);
				}
				
				if(filter != null && filter.isEnough(entry)) {
					break;
				}

			}

			//System.out.print("Closing "+fileName);
			br.close();
			//System.out.print(" ..... Closed\n");
	}
	
	public void load(String fileName, boolean loadCDSOnly, String version) throws ParseException, FileNotFoundException, IOException {
		geneList = new HashMap<String,GeneAnnotation>();
		File source = new File(fileName);
		BufferedReader br = new BufferedReader(new FileReader(source));
		String line;
		while((line = br.readLine()) != null) {
			if(line.startsWith("#")){
				continue;
			}
			GeneAnnotation entry = new GeneAnnotation(line.split("\t"), loadCDSOnly, version);
			geneList.put(entry.getName(),entry);
			addToChromsomeMap(entry);
		}

		//System.out.print("Closing "+fileName);
		br.close();
		//System.out.print(" ..... Closed\n");
	}
	public void load(String fileName, boolean loadCDSOnly) throws ParseException, FileNotFoundException, IOException {
		load(fileName, loadCDSOnly, null);
	}
	
	public List<Exon> getAllExons() {
		List<GeneAnnotation> genes = getAnnotationList();
		Iterator<GeneAnnotation> it = genes.iterator();
		List<Exon> exons = new ArrayList<Exon>(genes.size() * 4);
		while(it.hasNext()) {
			GeneAnnotation gene = it.next();
			exons.addAll(gene.getExons());
		}
		
		return exons;
	}
	
	public List<GenomicAnnotation> getAllIntrons() {
		List<GeneAnnotation> genes = getAnnotationList();
		Iterator<GeneAnnotation> it = genes.iterator();
		List<GenomicAnnotation> introns = new ArrayList<GenomicAnnotation>(genes.size() * (4 - 1));
		while(it.hasNext()) {
			GeneAnnotation gene = it.next();
			introns.addAll(gene.getIntrons());
		}
		
		return introns;		
		
	}
	
	public List<GenomicAnnotation> getAllPromoters(int promoterBases) {
		List<GeneAnnotation> genes = getAnnotationList();
		Iterator<GeneAnnotation> it = genes.iterator();
		List<GenomicAnnotation> promoters = new ArrayList<GenomicAnnotation>(genes.size() );
		while(it.hasNext()) {
			GeneAnnotation gene = it.next();
			promoters.add(gene.getPromoter(promoterBases));
		}
		
		return promoters;
	}
	

	public List<GeneAnnotation> getAllGenesBetween(GeneAnnotation start, GeneAnnotation end) 
	throws Exception{
		ArrayList<GeneAnnotation> result = new ArrayList<GeneAnnotation>();
		if(!start.getChromosome().equals(end.getChromosome())) {
			throw new Exception(start.getName() + " is in chromosome "+ start.getChromosome() + " while " + 
					end.getName() + " is in chromosome " + end.getChromosome() + " they must be in the same chromosome");
		}
		
		Map<String, IntervalTree<GeneAnnotation>> chromosomeGeneMap = super.getChromosomeTreeMap();
		Iterator<GeneAnnotation> geneIt = new IntervalTree.ValuesIterator<GeneAnnotation>(chromosomeGeneMap.get(start.getChromosome()).iterator());
		boolean pastEnd = false;
		while(geneIt.hasNext() && !pastEnd) {
			GeneAnnotation gene = geneIt.next();
			if (gene.getStart() > end.getEnd()) {
				pastEnd = true;
			}
			if(gene.getEnd() > start.getStart() && gene.getStart() < end.getEnd() ) {
				System.out.println("Adding new gene: "+gene + " with symbols " + getSymbolsForGene(gene.getName()));				
				result.add(gene);
			}
		}
		return result;
	}	
	
	/**
	 * This method intends to remove splice variants and redundant annotations.
	 * Since it inspects the loaded annotation tree, the annotation kept will be that that starts (or ends, if in opposite strant) the earliest.
	 * Unfortunately, overlapping transcripts that are not redundant nor variants will also be eliminated.
	 */
	public void eliminateOverlappingGenes() {
		for(String chr : getChromosomeTreeMap().keySet()) {
			IntervalTree<GeneAnnotation> originalTree = getChromosomeTree(chr);
			IntervalTree<GeneAnnotation> newTree      = new IntervalTree<GeneAnnotation>();
			GeneAnnotation last = null;
			Iterator<GeneAnnotation> annotationIt = originalTree.valueIterator();
			while(annotationIt.hasNext()) {
				GeneAnnotation current = annotationIt.next();
				if(last == null || ! current.overlaps(last)) {
					newTree.put(current.getStart(), current.getEnd(), current);
				}
				last = current;
			}
			System.err.println("Eliminated " + (originalTree.size() - newTree.size()));
			getChromosomeTreeMap().put(chr, newTree);
		}
	}
	
	public List<GeneAnnotation> getAllOverlappingGenes(LightweightGenomicAnnotation ga) {
		ArrayList<GeneAnnotation> result = new ArrayList<GeneAnnotation>();
		Map<String, IntervalTree<GeneAnnotation>> chromosomeGeneMap = super.getChromosomeTreeMap();
		IntervalTree<GeneAnnotation> geneTree = chromosomeGeneMap.get(ga.getChromosome());
		if(geneTree!= null) {
			Iterator<GeneAnnotation> geneIt = new IntervalTree.ValuesIterator<GeneAnnotation>(geneTree.overlappers(ga.getStart(), ga.getEnd()));
			while(geneIt.hasNext()) {
				GeneAnnotation gene = geneIt.next();
				result.add(gene);
			}
		}
		return result;
	}
	
	public GeneAnnotation getGene(String geneId) {
		return geneList.get(geneId);
	}
	
	public abstract void loadAliasTable(String fileName) throws IOException;
	public abstract List<GeneAnnotation> getAnnotationsForSymbol(String symbol);
	public abstract List<String> getSymbolsForGene(String id);
	
	public IntervalTree<Exon> getExonTree(String chromosome) {
		IntervalTree<Exon> exonTree = new IntervalTree<Exon>();
		IntervalTree<GeneAnnotation> chrGeneTree = getChromosomeTreeMap().get(chromosome);
		if(chrGeneTree != null) {
			Iterator<GeneAnnotation> geneIt = new IntervalTree.ValuesIterator<GeneAnnotation>(chrGeneTree.iterator());
			while(geneIt.hasNext()) {
				Iterator<Exon> exonIt = geneIt.next().getExons().iterator();
				while(exonIt.hasNext()) {
					Exon exon = exonIt.next();
					exonTree.put(exon.getStart(), exon.getEnd(), exon);
				}
			}
		}
		return exonTree;
	}
	
	private void addToChromsomeMap(GeneAnnotation entry) {
		super.addAnnotation(entry);
	}
	
	@Override
	public int parse(String file,
			GenomicAnnotationFilter<GeneAnnotation> filter,
			AnnotationHandler handler) throws ParseException, IOException {
		// TODO Auto-generated method stub
		return 0;
	}
	@Override
	public int parse(BufferedReader br,
			GenomicAnnotationFilter<GeneAnnotation> filter,
			AnnotationHandler handler) throws ParseException, IOException {
		// TODO Auto-generated method stub
		return 0;
	}
}
