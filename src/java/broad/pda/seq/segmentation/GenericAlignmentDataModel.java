package broad.pda.seq.segmentation;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.AlignmentBlock;
import org.broad.igv.sam.reader.AlignmentQueryReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.util.ResourceLocator;

import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.sequence.Sequence;
import broad.core.util.CollapseByIntersection;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.GeneTools;
import broad.pda.gene.RefSeqGene;
import broad.pda.samtools.SAMPairedEndFileReader;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CloseableIterator;


public class GenericAlignmentDataModel implements AlignmentDataModel{
	Logger logger =  Logger.getLogger(GenericAlignmentDataModel.class.getName());
	AlignmentQueryReader  reader;
	Map<String, Integer> chromosomeSizes;
	//Map<String, Integer> firstReads;
	Map<String, IntervalTree<Alignments>> cachingScoringTrees;
	Map<String,IntervalTree<Alignments>> gapTree;
	int chunkSize=1000000;
	int chunkEnd;
	int chunkStart;
	IntervalTree<Alignment> chunkAlignmentTree;
	String chunkChr;
	Node<Alignment> currentAlignment;
	boolean startedIteration;
	static String usage=" args[0]=file \n args[1]=sizes \n args[2]=region \n args[3]=save";
	CloseableIterator<Alignment> readIterator;   // JE - this looks like a hack to me
	double spliceWeightFactor;
	Map<String, IntervalTree<Alignments>> fullAlignmentTreesAsAlignmentsCached;
	double minMappingQuality;
	boolean isStranded = false;
	boolean isNegativeStrand;
	boolean isPaired = false;
	boolean isSecondRead;
	private ReadCountNormalizer normalizer;
	private String alignmentFileName;
	boolean pairsAsFragments;  // tells if it's reading a paired Alignment file.  Will only read paired reads, and they will be output as a single fragment (count).
	private int extensionFactor;
	
	/*
	 * Remove PCR duplicates
	 */
	private boolean removeDuplicatesFlag=true;
	/*
	 * Flag to weigh the read counts by the NH flag value
	 */
	private boolean weighReadCounts=false;
	
	/*
	 * OVERLOADED CONSTRUCTORS:
	 * String, String: 									filename, lengthFile
	 * String, String, boolean:							filename, lengthFile, upweightSplices
	 * String, String, boolean, boolean					filename, lengthFile, removeDuplicatesFlag, weighReadCounts
	 * String, String, boolean, boolean, boolean		filename, lengthFile, upweightSplices, removeDuplicatesFlag, weighReadCounts
	 * String, String, boolean, int						filename, lengthFile, upweightSplices, minMappingQuality
	 * String, String, boolean, double					filename, lengthFile, upweightSplices, minMappingQuality
	 * String, String, boolean, int, boolean			filename, lengthFile, upweightSplices, minMappingQuality, removeDuplicatesFlag
 	 * String, String, int, boolean						filename, lengthFile, minMappingQuality, removeDuplicatesFlag
	 *
	 // NOTE: Minimum mapping quality will be set to zero if weighReadCounts is set to true
	 * String, String, boolean, int, boolean, boolean	filename, lengthFile, upweightSplices, minMappingQuality, removeDuplicatesFlag, weighReadCounts
	 * 
	 */
	public GenericAlignmentDataModel(String filename, String lengthFile, boolean upweightSplices, double minMappingQuality, boolean removeDuplicatesFlag, boolean weighReadCounts, String strand, boolean loadPairsAsFragments) throws IOException{
		this.minMappingQuality = minMappingQuality;
		this.alignmentFileName = filename;
		this.pairsAsFragments = loadPairsAsFragments;
		
		if (pairsAsFragments) {
			this.reader = new SAMPairedEndFileReader(new File(filename));
		} else {
			this.reader = AlignmentReaderFactory.getReader(new ResourceLocator(filename), false);
		}
		
		if(lengthFile != null && !lengthFile.isEmpty()) {
			this.chromosomeSizes=BEDFileParser.loadChrSizes(lengthFile);
		} else {
			this.chromosomeSizes = new TreeMap<String, Integer>();
			SAMSequenceDictionary dictionary = reader.getHeader().getSequenceDictionary();
			if(dictionary != null && !dictionary.getSequences().isEmpty() ) {
				List<SAMSequenceRecord> refSequences = dictionary.getSequences();
				for(SAMSequenceRecord record : refSequences) {
					chromosomeSizes.put(record.getSequenceName(), record.getSequenceLength());
				}
			} else {
				throw new IllegalArgumentException("The size file describing reference sequence length was empty and the alignment header did not contain a data dictionary");
			}
		}
		
		//this.firstReads=firstReads();
		this.cachingScoringTrees=initializeCachingTree();
		this.chunkAlignmentTree=new IntervalTree<Alignment>();
		this.chunkEnd=0;
		this.chunkStart=0;
		this.chunkChr="";
		this.currentAlignment=null;
		this.startedIteration=false;
		spliceWeightFactor=setSpliceWeight(upweightSplices);
		this.alignmentFileName = filename;
		
		Globals.setHeadless(true);
		
		
		setRemoveDuplicatesFlag(removeDuplicatesFlag);
		
		// NOTE: Minimum mapping quality will be set to zero if weighReadCounts is set to true
		setWeighReadCountsFlag(weighReadCounts);
		
		if (strand != null) {
			if("+".equals(strand)){
				this.setPositiveStranded();
			}
			else{
				this.setNegativeStranded();
			}
		}
	}
	
	public GenericAlignmentDataModel(String filename, String lengthFile, boolean upweightSplices, double minMappingQuality) throws IOException {
		this(filename, lengthFile, upweightSplices, minMappingQuality, false, false, null, false);
	}
	
	public GenericAlignmentDataModel(String filename, String lengthFile) throws IOException{
		this(filename, lengthFile, false,0);
	}

	public GenericAlignmentDataModel(String filename, String lengthFile, int minMappingQuality) throws IOException{
		this(filename, lengthFile, false, minMappingQuality);
	}

	public GenericAlignmentDataModel(String filename, String lengthFile, boolean upweightSplices) throws IOException{
		this(filename, lengthFile, false, 0);
	}

	//public int getFirstReads(String chr){return this.firstReads.get(chr);}
	
	public GenericAlignmentDataModel(String filename, String lengthFile, boolean upweightSplices, int minMappingQuality, boolean removeDuplicatesFlag) throws IOException{
		this(filename, lengthFile, upweightSplices, minMappingQuality, removeDuplicatesFlag, false);
	}
	
	public GenericAlignmentDataModel(String filename, String lengthFile, int minMappingQuality, boolean removeDuplicatesFlag) throws IOException{
		this(filename, lengthFile, false, minMappingQuality, removeDuplicatesFlag);
	}
	
	public GenericAlignmentDataModel(String filename, String lengthFile, boolean removeDuplicatesFlag, boolean weighReadCounts) throws IOException{
		this(filename, lengthFile, false, removeDuplicatesFlag, weighReadCounts);
	}
	
	public GenericAlignmentDataModel(String filename, String lengthFile, boolean upweightSplices, boolean removeDuplicatesFlag, boolean weighReadCounts) throws IOException{
		this(filename, lengthFile, upweightSplices, 0, removeDuplicatesFlag, weighReadCounts);
	}
	
	public GenericAlignmentDataModel(String filename, String lengthFile, boolean upweightSplices, int minMappingQuality, boolean removeDuplicatesFlag, boolean weighReadCounts) throws IOException{
		this(filename, lengthFile, upweightSplices, minMappingQuality, removeDuplicatesFlag, weighReadCounts, null);
	}

	public GenericAlignmentDataModel(String filename, String lengthFile, boolean upweightSplices, int minMappingQuality, boolean removeDuplicatesFlag, boolean weighReadCounts,String strand) throws IOException{
		this(filename, lengthFile, upweightSplices, minMappingQuality, removeDuplicatesFlag, weighReadCounts, strand, false);
	}
	
	
	//Get all alignments overlapping a specified region
	@Override
	public CloseableIterator<Alignment> getAlignmentsOverlappingRegion(Alignments align) throws IOException{
		CloseableIterator<Alignment> iter =  reader.query(align.getChromosome(), align.getStart(), align.getEnd(), false);
		return iter;
	}

	@Override
	public double getBasesCoveredPerAlignment(RefSeqGene gene, IntervalTree<Alignment> tree, int EF){
		double sum=0;
		double counter=0;
		//get Alignments over the whole region
		Iterator<Node<Alignment>> iter=tree.overlappers(gene.getStart(), gene.getEnd()+1); //TODO Consider getting rid of the plus 1
		
		while(iter.hasNext()){
			Node<Alignment> node=iter.next();
			int num=node.getNumReplicates();
			Alignment alignment=node.getValue();
			//only count those overlapping one of the exons
			counter+=(getBasesCoveredPerAlignment(alignment, gene, EF)*num);
		}
		return counter;
	}

	
	
	@Override
	public double getCountsPerAlignmentFullyContained(RefSeqGene gene, IntervalTree<Alignment> tree, int EF){
		double counter=0;
		//get Alignments over the whole region
		Iterator<Node<Alignment>> iter=tree.overlappers(gene.getStart(), gene.getEnd()+1);
		
		while(iter.hasNext()){
			Node<Alignment> node=iter.next();
			int num=node.getNumReplicates();
			Alignment alignment=node.getValue();
			//only count those overlapping one of the exons
			counter+=(this.getCountsPerAlignmentFullyContained(alignment, gene, EF)*num);
		}
		return counter;
	}
	
	private double getCountsPerAlignmentFullyContained(Alignment record, RefSeqGene gene, int EF){
		double counter=0;
		if(!passesQualAndStrandness(record)) {
    		return 0;
    	}
		Collection<Alignments> exons=gene.getExonSet();
		
		boolean count=false;
	    if(record.getAlignmentBlocks()!=null){
	    	for(AlignmentBlock b : record.getAlignmentBlocks()) {
	    		Alignments exon=new Alignments(record.getChromosome(), b.getStart(), b.getEnd());
	        	Alignments read=this.extendRead(exon, record.isNegativeStrand(), EF);
	        	for(Alignments align: exons){
	        		if(align.fullyContained(read)){
	        			count=true; 
	        			continue;
	        		}
	        	}
	        }
	    }
	    else{
	    	Alignments read=this.extendRead(record, record.isNegativeStrand(), EF);
	        for(Alignments align: exons){
	    	if(align.fullyContained(read)){count=true;}
	        }
	   }
	   if(count){counter++;}
	    
	   if(record.getAlignmentBlocks()!=null && record.getAlignmentBlocks().length>1){
		   counter=counter*this.spliceWeightFactor;
	   }
	   
		return normalizer == null ? counter : normalizer.normalize(counter);
	}
	
		//Moran:Counts the number of bases that were actually covered  (for exact coverage calculation)
		private double getBasesCoveredPerAlignment(Alignment record, RefSeqGene gene, int EF){
		double counter=0;
		if(!passesQualAndStrandness(record)) {
			return 0;
		}
		Collection<Alignments> exons=gene.getExonSet();
	    if(record.getAlignmentBlocks()!=null){
	    	for(AlignmentBlock b : record.getAlignmentBlocks()) {
	    		Alignments exon=new Alignments(record.getChromosome(), b.getStart(), b.getEnd());
	        	Alignments extended=this.extendRead(exon, record.isNegativeStrand(), EF);
	        	for(Alignments align: exons){
	        		if(extended.overlaps(align))
	        			counter+= extended.getOverlap(align);
	        		}
	        }
	    }
	    else{
	    	Alignments extended=this.extendRead(record, record.isNegativeStrand(), EF);
	        for(Alignments align: exons){
	    	if(extended.overlaps(align)){counter+= extended.getOverlap(align);}
	        }
	   }
		return counter;
	}
	
		
		/**
		 * @author Jesse Engreitz
		 * @param region	region of interest
		 * @param EF		extension factor
		 * @return Number of bases covered by at least one read
		 * Warning: This is not a very efficient function with lots of reads  TODO rewrite using merge
		 */
		@Override
		public int getBasesCovered(Alignments region, int EF) throws IOException {
			Alignments extended = extendRead(region, false, EF);
			boolean covered[] = new boolean[region.getSize()];
			CloseableIterator<Alignment> iter =  reader.query(extended.getChromosome(), extended.getStart(), extended.getEnd(), false);
			while (iter.hasNext()) {
				Alignment record = iter.next();
				Alignments extendedRecord = extendRead(new Alignments(record), record.isNegativeStrand(), EF);
				if(passesQualAndStrandness(record)) {
					int startIndex = Math.max(0, extendedRecord.getStart() - region.getStart());
					int endIndex = Math.min(region.getSize(), extendedRecord.getEnd());
					for (int i = startIndex; i < endIndex; i++) {
						covered[i] = true;
					}
				}
			}
			
			int counter = 0;
			for (int i = 0; i < covered.length; i++) counter += covered[i] ? 1 : 0;
			return counter;
		}

		
		@Override
		public RefSeqGene getPeak(RefSeqGene gene, RefSeqGene startCodon, IntervalTree<Alignment> tree, int EF){
			//get Alignments over the whole region
			Iterator<Node<Alignment>> iter=tree.overlappers(gene.getStart(), gene.getEnd()+1); //TODO Consider getting rid of the plus 1
			
			Collection<Alignments> alignments=new TreeSet<Alignments>();
			
			int min=Integer.MAX_VALUE;
			int max=-Integer.MAX_VALUE;
			boolean updated=false;
			
			while(iter.hasNext()){
				Node<Alignment> node=iter.next();
				Alignment alignment=node.getValue();
				//only count those overlapping one of the exons
				Collection<Alignments> exons=startCodon.getExonSet();
				
				 if(alignment.getAlignmentBlocks()!=null){
					boolean overlapsExon=false;
					Collection<Alignments> tmp=new TreeSet<Alignments>();
			    	for(AlignmentBlock b : alignment.getAlignmentBlocks()) {
			    		Alignments exon=new Alignments(alignment.getChromosome(), b.getStart(), b.getEnd());
			        	Alignments extended=this.extendRead(exon, alignment.isNegativeStrand(), EF);
			        	tmp.add(extended);
			        	for(Alignments align: exons){
			        		if(extended.overlaps(align)){
			        			overlapsExon=true;
			        		}
			        	}
			        }
			    	if(overlapsExon){
			    		alignments.addAll(tmp);
			    		updated=true;
			    	}
			    }
			    else{
			    	Alignments extended=this.extendRead(alignment, alignment.isNegativeStrand(), EF);
			        for(Alignments align: exons){
			    	if(extended.overlaps(align)){alignments.add(extended); updated=true;}
			        }
			   }
			}
			   
				
			//System.err.println("UPDATED "+updated);
			if(!updated){return null;}
			return getPeak(gene, alignments);
		}
		
		
		private RefSeqGene getPeak(RefSeqGene gene, Collection<Alignments> alignments) {
			int min=Integer.MAX_VALUE;
			int max=-Integer.MAX_VALUE;
			
			for(Alignments align: alignments){
				min=Math.min(align.getStart(), min);
				max=Math.max(align.getEnd(), max);
			}
			
			//System.err.println(min+" "+max);
			
			min=Math.max(min, gene.getStart());
			max=Math.min(max, gene.getEnd());
			
			//System.err.println(min+" "+max);
			//System.out.println(gene.getChr()+"\t"+min+"\t"+max);
			//System.out.println(gene.getChr()+"\t"+gene.getStart()+"\t"+gene.getEnd());
			
			//TODO Check if doesnt overlap exon if not first exon
			
			RefSeqGene rtrn=gene.trimAbsolute(min, max);
			return rtrn;
		}

	//Get all represented chromosomes and their lengths in the genome
		@Override
		public Map<String, Integer> getChromosomeLengths(){
			return this.chromosomeSizes;
		}

		@Override
		public int getChromosomeLength(String chr) { return chromosomeSizes.get(chr);}

	@Override
	public CloseableIterator<Alignment> getChromosomeReadIterator(String chromosome) throws IOException {
		CloseableIterator<Alignment> rtrn = new EmptyCloseableIterator();
		if(chromosomeSizes.containsKey(chromosome)) {
			Alignments chromosomeRegion = new Alignments(chromosome, 0, chromosomeSizes.get(chromosome));
			rtrn = getAlignmentsOverlappingRegion(chromosomeRegion);
		}
		
		return rtrn;
	}

	
	@Override
	public int getChunkEnd(){return this.chunkEnd;}

	@Override
	public int getChunkStart(){return this.chunkStart;}

	//get the number of alignments overlapping a given region
	@Override
	public double getCounts(String chr) throws IOException{
		CloseableIterator<Alignment> iter =  reader.query(chr, 0, getChromosomeLength(chr), false);
		double counter=getCountsPerAlignmentUncached(iter);
	   	iter.close();
	    return counter;
	}

	
	@Override
	public double getCountsOfUniqueOverlappers(Alignments target, Alignments exclude, IntervalTree<Alignment> tree, int EF) {
		double counter=0;
		Iterator<Node<Alignment>> iter=tree.overlappers(target.getStart(), target.getEnd()); //It was target.getEnd()+1 which will affect the right end.
		//System.err.println("\t\ttarget: " + target.toUCSC() + " exclude: " + exclude.toUCSC() + " found target overlapping reads? " + iter.hasNext());
		while(iter.hasNext()){
			Node<Alignment> node=iter.next();
			int num=node.getNumReplicates();
			Alignment read =node.getValue();
			//System.err.println("Target " + target.toUCSC() + " exclude " + exclude.toUCSC() + " min qual "+minMappingQuality+" read quality " + read.getMappingQuality() + " passes filter? " + passesQualAndStrandness(read) + " numPlaces " + num);
			if(passesQualAndStrandness(read)) {
				Alignments readAnnotation = new Alignments(read.getChr(), read.getStart(), read.getEnd());
				if(target.overlaps(readAnnotation)) { // Interval tree seems to return overlapping reads to the closed-closed interval rather than closed-open interval, so recheck we got what we wanted.
					counter += getCountsThatDontOverlapAlignment(node.getValue(),  exclude, EF) * num;
					
				}
			}
		}
		return counter;	
		
	}

	@Override
	public double getCountsOfUniqueOverlappers(Alignments target, RefSeqGene exclude, IntervalTree<Alignment> tree, int EF) {
		double counter=0;
		Iterator<Node<Alignment>> iter=tree.overlappers(target.getStart(), target.getEnd());
		while(iter.hasNext()){
			Node<Alignment> node=iter.next();
			Alignment read =node.getValue(); 
			if(passesQualAndStrandness(read)) {
				Alignments readAnnotation = new Alignments(read.getChr(), read.getStart(), read.getEnd());
				if(target.overlaps(readAnnotation)) { // Interval tree seems to return overlapping reads to the closed-closed interval rather than closed-open interval, so recheck we got what we wanted.
					int num=node.getNumReplicates();
					counter+=(getCountsThatDontOverlapAlignment(read,  exclude.getExonSet(), EF)*num);
				}
			}
		}
		return counter;	
		
	}

	@Override
	public double getCountsOnBase(String chr, int index) throws IOException {
		CloseableIterator<Alignment> iter = reader.query(chr, index, index + 1, false);
		double counter = getCountsPerAlignmentUncached(iter);
		iter.close();
		return counter;
	}

	
	public double getCountsOnWindow(String chr, int startIndex, int endIndex) throws IOException {
		return getCountsPerAlignmentUncached(chr, startIndex, endIndex);
	}

	
	public double getCountsPerAlignmentUncached(String chr, int startIndex, int endIndex) throws IOException {
		CloseableIterator<Alignment> iter = reader.query(chr, startIndex, endIndex, false);
		double counter = getCountsPerAlignmentUncached(iter);
		iter.close();
		return counter;
	}

	private double getCountsPerAlignmentUncached(Iterator<Alignment> iter){
		double counter=0;
		while (iter.hasNext()) {
			Alignment record = iter.next();
			if(passesQualAndStrandness(record)) {
				counter = counter + countReads(record);
			}
		}
		return normalizer == null ? counter : normalizer.normalize(counter);
	}
	
	private double getCountsPerAlignmentUncached(Alignments target, Iterator<Alignment> iter, int EF){
		double counter=0;
		while (iter.hasNext()) {
			Alignment align = iter.next();
			counter += getCountsPerAlignment(align, target, EF);
		}
		return counter;
	}
	
	private double getCountsPerAlignmentUncached(Alignments target, int EF) throws IOException {
		CloseableIterator<Alignment> iter = reader.query(target.getChromosome(), target.getStart()-EF, target.getEnd()+EF, false);
		double count = getCountsPerAlignmentUncached(target, iter, EF);
		iter.close();
		return count;
	}
	
	
	//get the number of alignments overlapping a given region
	public double getCountsPerAlignment(Alignments align, int EF) throws IOException {
		// JE 9/27/12 - Use the cached interval tree instead of accessing the reader directly every time
		// Old version:
		//CloseableIterator<Alignment> iter =  reader.query(extended.getChromosome(), extended.getStart(), extended.getEnd(), false);
		//double counter=getCountsPerAlignment(iter, align, EF);
	   	//iter.close();
		
		double count = 0;
		
		// If the requested target is larger than the chunkSize, do an uncached iteration
		if (align.getEnd() - align.getStart() > chunkSize) {
			count = getCountsPerAlignmentUncached(align, EF);
		} else {  // otherwise, load a cached interval tree
			IntervalTree<Alignment> tree = getIntervalTreeCached(align.getChromosome(), align.getStart()-EF, align.getEnd()+EF);
			count = getCountsPerAlignment(align, tree, EF);
		}
		
	    return count;
	}
	
	//get the number of alignments overlapping a given region
	public double getCountsPerAlignment(Alignments align) throws IOException{
		int EF = extensionFactor;
		return getCountsPerAlignment(align, EF);
	}
	


	private double getCountsPerAlignment(Alignment record, Alignments align, int EF){
		// TODO:  implement an abstract "AlignmentFilter" class and pass a list of these 
		// to this function, so that other functions like "getCountsPerAlignmentForPolyA" can just
		// customize the filter class.  JE
		
		double counter=0;
		if(!passesQualAndStrandness(record)) {
			return 0;
		}
		boolean count=false;
	    if(record.getAlignmentBlocks()!=null){
	    	
	    	for(AlignmentBlock b : record.getAlignmentBlocks()) {
	    		Alignments exon=new Alignments(record.getChromosome(), b.getStart(), b.getEnd());
	        	Alignments extended=extendRead(exon, record.isNegativeStrand(), EF);
	        	if(extended.overlaps(align)){count=true;}
	        }
	    }
	    else{
	    	Alignments extended=extendRead(record, record.isNegativeStrand(), EF);
	        if(extended.overlaps(align)){count=true;}
	   }
	   if(count ){counter = counter + countReads(record);}
	    
		return normalizer == null ? counter : normalizer.normalize(counter);
	}
	
	
	private double getCountsPerAlignment(Iterator<Node<Alignment>> iter, Alignments align, int EF) {
		double counter = 0;
		
		while(iter.hasNext()){
			Node<Alignment> node=iter.next();
			int num=node.getNumReplicates();
			Alignment alignment=node.getValue();
			counter+=(getCountsPerAlignment(alignment, align, EF)*num);
		}
		
		return counter;
	}

	//get the number of alignments overlapping a given region with a cached interval tree
	//Note the node in the tree has a counter of how many replicates spanned that start-end region
	@Override
	public double getCountsPerAlignment(Alignments align, IntervalTree<Alignment> tree, int EF) {
		Iterator<Node<Alignment>> iter = tree.overlappers(align.getStart(), align.getEnd()+1); //TODO Consider getting rid of the plus 1
		return getCountsPerAlignment(iter, align, EF);
	}
	
	
	public boolean loadedPairsAsFragments() {
		return pairsAsFragments;
	}
	
	public double getCountsPerAlignmentForPolyA(Alignments align, int EF) throws IOException{
		Alignments extended=new Alignments(align.getChromosome(), align.getStart()-EF, align.getEnd()+EF);
		CloseableIterator<Alignment> iter =  reader.query(extended.getChromosome(), extended.getStart(), extended.getEnd(), false);
		double counter=getCountsPerAlignmentForPolyA(iter, align, EF);
	   	iter.close();
	    return counter;
	}

	// TODO: merge with getCountsPerAlignment JE
	private double getCountsPerAlignmentForPolyA(Iterator<Alignment> iter, Alignments align, int EF){
		double counter=0;
		while (iter.hasNext()) {
			boolean count=false;
	        Alignment record = iter.next();
	        if(passesQualAndStrandness(record)) {
		        if(record.getAlignmentBlocks()!=null){
		        	for(AlignmentBlock b : record.getAlignmentBlocks()) {
		        		Alignments exonExtended=new Alignments(record.getChromosome(), b.getStart()-EF, b.getEnd()+EF);
		        		if(exonExtended.overlaps(align)){count=true;}
		        	}
		        }
		        else{
		        	Alignments extended=new Alignments(record.getChromosome(), record.getStart()-EF, record.getEnd()+EF);
		        	if(extended.overlaps(align) ){count=true;}
		        }
		        if(count){counter = counter + countReads(record);}
	        }
	    }
		return normalizer == null ? counter : normalizer.normalize(counter);
	}

	
	
	@Override
	public double getCountsPerAlignment(RefSeqGene gene, IntervalTree<Alignment> tree, int EF){
		double counter=0;
		//get Alignments over the whole region
		Iterator<Node<Alignment>> iter=tree.overlappers(gene.getStart(), gene.getEnd()+1); //TODO Consider getting rid of the plus 1
		
		while(iter.hasNext()){
			Node<Alignment> node=iter.next();
			int num=node.getNumReplicates();
			Alignment alignment=node.getValue();
			//only count those overlapping one of the exons
			counter+=(getCountsPerAlignment(alignment, gene, EF)*num);
		}
		return counter;
	}

	@Override
	public double getCountsPerAlignmentStranded(RefSeqGene gene, IntervalTree<Alignment> tree, int EF, String orientation){
		double sum=0;
		//Alignments[] exons=gene.getExons();
		
		double counter=0;
		//get Alignments over the whole region
		Iterator<Node<Alignment>> iter=tree.overlappers(gene.getStart(), gene.getEnd()+1); //TODO Consider getting rid of the plus 1
		
		while(iter.hasNext()){
			Node<Alignment> node=iter.next();
			int num=node.getNumReplicates();
			Alignment alignment=node.getValue();
			//only count those overlapping one of the exons
			boolean neg=alignment.isNegativeStrand();
			String strand="+";
			if(neg){strand="-";}
			if(strand.equalsIgnoreCase(orientation)){
				counter+=(getCountsPerAlignment(alignment, gene, EF)*num);
			}
		}
		return counter;
	}

	//get counts for an array of alignments
	@Override
	public double getCountsPerAlignment(Alignments[] alignments, IntervalTree<Alignment> tree, int EF){
		double counter=0;
		for(int i=0; i<alignments.length; i++){
			counter+=this.getCountsPerAlignment(alignments[i], tree, EF);
		}
		return counter;
	}

	//get counts for an array of alignments
	@Override
	public double getCountsPerAlignment(Alignments[] alignments,  int EF) throws IOException{
		double counter=0;
		for(int i=0; i<alignments.length; i++){
			//IntervalTree<Alignment> tree = getIntervalTreeCached(alignments[i].getChr(), alignments[i].getStart(), alignments[i].getEnd());
			//counter+=this.getCountsPerAlignment(alignments[i], tree, EF);
			counter += this.getCountsPerAlignment(alignments[i], EF);  // JE 9/26/12 this function was changed to use the intervalTree already
		}
		return counter;
	}

	private double getCountsPerAlignment(Alignment record, RefSeqGene gene, int EF){
		double counter=0;
		if(!passesQualAndStrandness(record)) {
			return 0;
		}
		Collection<Alignments> exons=gene.getExonSet();
		
		boolean count=false;
	    if(record.getAlignmentBlocks()!=null){
	    	for(AlignmentBlock b : record.getAlignmentBlocks()) {
	    		Alignments exon=new Alignments(record.getChromosome(), b.getStart(), b.getEnd());
	        	Alignments extended=this.extendRead(exon, record.isNegativeStrand(), EF);
	        	for(Alignments align: exons){
	        		if(extended.overlaps(align)){
	        			count=true; 
	        			continue;
	        		}
	        	}
	        }
	    }
	    else{
	    	Alignments extended=this.extendRead(record, record.isNegativeStrand(), EF);
	        for(Alignments align: exons){
	    	if(extended.overlaps(align)){count=true;}
	        }
	   }
	   if(count){counter = counter + countReads(record);}
	    
	   if(record.getAlignmentBlocks()!=null && record.getAlignmentBlocks().length>1){
		   counter=counter*this.spliceWeightFactor;
	   }
	
		return normalizer == null ? counter : normalizer.normalize(counter);
	}
	
	@Override
	public double getCountsPerAlignment(Alignments align,
			Map<String, IntervalTree<Alignments>> goodExonTree,
			IntervalTree<Alignment> tree, int extensionFactor) {
	
		//Alignments[] exons=gene.getExons();
		
		double counter=0;
		//get Alignments over the whole region
		Iterator<Node<Alignment>> iter=tree.overlappers(align.getStart(), align.getEnd()+1); //TODO Consider getting rid of the plus 1
		
		while(iter.hasNext()){
			Node<Alignment> node=iter.next();
			int num=node.getNumReplicates();
			Alignment alignment=node.getValue();
			//only count those overlapping one of the exons
			counter+=(getCountsPerAlignment(alignment, align, goodExonTree, extensionFactor)*num);
		}
		return  counter;
	}

	private double getCountsPerAlignment(Alignment record, Alignments align,
			Map<String, IntervalTree<Alignments>> goodExonTree,
			int extensionFactor) {
		
		if(!passesQualAndStrandness(record)) {
			return 0;
		}
		
		double counter=0;
		
		boolean count=false;
	    if(record.getAlignmentBlocks()!=null){
	    	for(AlignmentBlock b : record.getAlignmentBlocks()) {
	    		Alignments exon=new Alignments(record.getChromosome(), b.getStart(), b.getEnd());
	        	Alignments extended=this.extendRead(exon, record.isNegativeStrand(), extensionFactor);
	        	if(extended.overlaps(align) && !goodExonTree.get(extended.getChr()).overlappers(extended.getStart(), extended.getEnd()).hasNext()){count=true;}
	        }
	    }
	    else{
	    	Alignments extended=this.extendRead(record, record.isNegativeStrand(), extensionFactor);
	        if(extended.overlaps(align) && !goodExonTree.get(extended.getChr()).overlappers(extended.getStart(), extended.getEnd()).hasNext()){count=true;}
	   }
	   if(count){counter = counter + countReads(record);}
	    counter = normalizer != null ? (int)normalizer.normalize(counter) : counter;
		return counter;
		
	}

	@Override
	public double getCountsPerAlignment(LightweightGenomicAnnotation region,int EF) throws IOException {
		Alignments align=new Alignments(region.getChromosome(), region.getStart(), region.getEnd());
		return getCountsPerAlignment(align, EF);
	}


	private double getCountsPerAlignmentWithSameEndpointForPolyA(Iterator<Alignment> iter, Alignments align, int EF, boolean polyAIsOnRight, int typicalOriginalAlignmentLength, int maxDistanceFromOriginalAlignmentEndpoint) {
		double counter=0;
		while (iter.hasNext()) {
			boolean count=false;
	        Alignment record = iter.next();
	        if(record.getAlignmentBlocks()!=null){
	        	if(!passesQualAndStrandness(record)) {
	        		continue;
	        	}
	        	for(AlignmentBlock b : record.getAlignmentBlocks()) {	        			
	        			if (b.getEnd() - b.getStart() < typicalOriginalAlignmentLength - 2)
	        				continue;
	        			if (polyAIsOnRight) {
	        				if (Math.abs(b.getStart() - align.getStart()) <= maxDistanceFromOriginalAlignmentEndpoint)
	        					count = true;
	        			} else {
	        				if (Math.abs(b.getEnd() - align.getEnd()) <= maxDistanceFromOriginalAlignmentEndpoint)
	        					count = true;
	        			}
	        	}
	        }
	        else{
	        		if (record.getEnd() - record.getStart() < typicalOriginalAlignmentLength - 2)
	        			continue;
	    			if (polyAIsOnRight) {
	    				if (Math.abs(record.getStart() - align.getStart()) <= maxDistanceFromOriginalAlignmentEndpoint)
	    					count = true;
	    			} else {
	    				if (Math.abs(record.getEnd() - align.getEnd()) <= maxDistanceFromOriginalAlignmentEndpoint)
	    					count = true;
	    			}
	        }
	        if(count){counter = counter + countReads(record);}
	    }
		return normalizer == null ? counter : normalizer.normalize(counter);
	}


	public double getCountsPerAlignmentWithSameEndpointForPolyA(Alignments align, int EF, boolean polyAIsOnRight, int typicalOriginalAlignmentLength, int maxDistanceFromOriginalAlignmentEndpoint) throws IOException{
		Alignments extended=new Alignments(align.getChromosome(), align.getStart()-EF, align.getEnd()+EF);
		CloseableIterator<Alignment> iter =  reader.query(extended.getChromosome(), extended.getStart(), extended.getEnd(), false);
		double counter=getCountsPerAlignmentWithSameEndpointForPolyA(iter, align, EF, polyAIsOnRight, typicalOriginalAlignmentLength, maxDistanceFromOriginalAlignmentEndpoint);
	   	iter.close();
	    return counter;
	}


	/**
	 * Returns the number of reads that do not overlap the given Alignment object
	 * @param iter Iterator of reads to test
	 * @param align Alignment object
	 * @param EF extension factor
	 * @return Number of reads that do not overlap the Alignment object
	 */
	//Cant overlap align at all
	//also cant overlap an intron at all
	private double getCountsThatDontOverlapAlignment(Iterator<Alignment> iter,  Alignments align, int EF){
		double counter=0;
		while (iter.hasNext()) {
			boolean count=false;
			boolean hasExon=false;
	        Alignment record = iter.next();
			if(!passesQualAndStrandness(record)) {
	    		continue;
	    	}
	        if(record.getAlignmentBlocks()!=null){
	        	//cant overlap ANY of the blocks
	        	for(AlignmentBlock b : record.getAlignmentBlocks()) {
	        		Alignments exon=new Alignments(record.getChromosome(), b.getStart(), b.getEnd());
	        		Alignments extended=extendRead(exon, record.isNegativeStrand(), EF);
	        		if(extended.overlaps(align)){
	        			count=true;
	        			break;
	        		}
	        		hasExon=true;
	        	}
	        }
	        else{
	        	Alignments extended=this.extendRead(record, record.isNegativeStrand(), EF);
	        	if(extended.overlaps(align)){count=true;}
	        	hasExon=true;
	        	//System.err.println(position+" "+extended.toUCSC());
	        }
	        if(!count && hasExon){counter++;}
	    }
		return normalizer == null ? counter : normalizer.normalize(counter);
	}

	private double getCountsThatDontOverlapAlignment(Alignment record,  Alignments align, int EF){
		if(!passesQualAndStrandness(record)) {
			return 0;
		}
		boolean overlaps = false;
		AlignmentBlock[] blocks = record.getAlignmentBlocks();
		if(blocks == null || blocks.length == 0) {
			Alignments extendedRecord  = extendRead(record, record.isNegativeStrand(), EF);
			overlaps = extendedRecord.overlaps(align);
			//System.err.println("\t\tExtendedread: " + extendedRecord.toUCSC() + " overlaps align " + align.toUCSC() + "? " + overlaps);
		} else {
			for (AlignmentBlock b : blocks) {
				Alignments block=new Alignments(record.getChromosome(), b.getStart(), b.getEnd());
				Alignments extended=extendRead(block, record.isNegativeStrand(), EF);
				if(extended.overlaps(align)){
					overlaps = true;
					//System.err.println("\t\tExtended block: " + extended.toUCSC() + " overlaps align " + align.toUCSC() + "? " + overlaps);
					break;
				}			
	
			}
		}
		//System.err.println("\t\trecord: " + record.getChr()+":"+record.getStart()+"-"+record.getEnd() + " will be excluded? " +!overlaps);
		/*double counter=0;
		boolean count=false;
		boolean hasExon=false;
		if(record.getAlignmentBlocks()!=null){
		   //cant overlap ANY of the blocks
	       	for(AlignmentBlock b : record.getAlignmentBlocks()) {
	       		Alignments exon=new Alignments(record.getChromosome(), b.getStart(), b.getEnd());
	       		Alignments extended=extendRead(exon, record.isNegativeStrand(), EF);
	       		if(extended.overlaps(align)){
	       			count=true;
	       			break;
	       		}
	       		hasExon=true;
	       	}
	    }
	    else{
	    	Alignments extended=this.extendRead(record, record.isNegativeStrand(), EF);
	    	count = extended.overlaps(align);
	        hasExon=true;
	        //System.err.println(position+" "+extended.toUCSC());
	    }
	    if(!count && hasExon){counter++;}
		 */
		return overlaps ? 0 :
			(normalizer == null ? 1 : normalizer.normalize(1)) ;
	
	}

	private double getCountsThatDontOverlapAlignment(Iterator<Alignment> iter,  Collection<Alignments> alignments, int EF){
		
		Map<String, IntervalTree<Alignments>> trees=makeTree(alignments);
		
		double counter=0;
		while (iter.hasNext()) {
			boolean count=false;
	        Alignment record = iter.next();
			if(!passesQualAndStrandness(record)) {
	    		continue;
	    	}
	        if(record.getAlignmentBlocks()!=null){
	        	//cant overlap ANY of the blocks
	        	for(AlignmentBlock b : record.getAlignmentBlocks()) {
	        		Alignments exon=new Alignments(record.getChromosome(), b.getStart(), b.getEnd());
	        		Alignments extended=extendRead(exon, record.isNegativeStrand(), EF);
	        		int numOverlaps=trees.get(extended.getChr()).numOverlappers(extended.getStart(), extended.getEnd());
	        		if(numOverlaps>0){
	        			count=true;
	        			break;
	        		}
	        	}
	        }
	        else{
	        	Alignments extended=this.extendRead(record, record.isNegativeStrand(), EF);
	        	int numOverlaps=trees.get(extended.getChr()).numOverlappers(extended.getStart(), extended.getEnd());
	        	if(numOverlaps>0){count=true;}
	        	//System.err.println(position+" "+extended.toUCSC());
	        }
	        if(!count){counter = counter + countReads(record);}
	    }
		return normalizer == null ? counter : normalizer.normalize(counter);
	}

	private double getCountsThatDontOverlapAlignment(Alignment record,  Collection<Alignments> alignments, int EF){
		if(!passesQualAndStrandness(record)) {
			return 0;
		}
		Map<String, IntervalTree<Alignments>> trees=makeTree(alignments);
		
		double counter=0;
		boolean count=false;
	    if(record.getAlignmentBlocks()!=null){
	    	//cant overlap ANY of the blocks
	        for(AlignmentBlock b : record.getAlignmentBlocks()) {
	        	Alignments exon=new Alignments(record.getChromosome(), b.getStart(), b.getEnd());
	        	Alignments extended=extendRead(exon, record.isNegativeStrand(), EF);
	        	int numOverlaps=trees.get(extended.getChr()).numOverlappers(extended.getStart(), extended.getEnd());
	        	if(numOverlaps>0){
	        		count=true;
	        		break;
	        	}
	        }
	    }
	    else{
	       	Alignments extended=this.extendRead(record, record.isNegativeStrand(), EF);
	       	int numOverlaps=trees.get(extended.getChr()).numOverlappers(extended.getStart(), extended.getEnd());
	       	if(numOverlaps>0){count=true;}
	       	//System.err.println(position+" "+extended.toUCSC());
	    }
	    if(!count){counter = counter + countReads(record);}
	    
	    return normalizer == null ? counter : normalizer.normalize(counter);
	}

	@Override
	public int getCountsWithinExons(Alignments align, IntervalTree<Alignment> tree, int EF){
		int counter=0;
		Iterator<Node<Alignment>> iter=tree.overlappers(align.getStart(), align.getEnd()+1); //TODO Consider getting rid of the plus 1
		
		while(iter.hasNext()){
			Node<Alignment> node=iter.next();
			int num=node.getNumReplicates();
			Alignment alignment=node.getValue();
			 if(alignment.getAlignmentBlocks()!=null && alignment.getAlignmentBlocks().length==1){
				 counter+=(getCountsPerAlignment(alignment, align, EF)*num);
			 }
		}
		return counter;
	}

	@Override
	public int getCountsWithinExons(Alignments align, Iterator<Alignment> iter, int EF){
		int counter=0;
		//Iterator<Node<Alignment>> iter=tree.overlappers(align.getStart(), align.getEnd()+1); //TODO Consider getting rid of the plus 1
		
		int i=0;
		while(iter.hasNext()){
			//System.err.println("Read "+i);
			Alignment alignment=iter.next();
			//int num=node.getNumReplicates();
			//Alignment alignment=node.getValue();
			 if(alignment.getAlignmentBlocks()!=null && alignment.getAlignmentBlocks().length==1){
				 counter+=(getCountsPerAlignment(alignment, align, EF));
			 }
			 i++;
		}
		return counter;
	}

	
	//get the exonic regions overlapping a given region
	@Override
	public Iterator<Alignments> getExonAlignmentsOverlappingRegion(Alignments align) throws IOException{
		Set<Alignments> s=new TreeSet<Alignments>();
		CloseableIterator<Alignment> iter=getAlignmentsOverlappingRegion(align);
		while(iter.hasNext()){
			Alignment record=iter.next();
			if(passesQualAndStrandness(record)) {
				Collection<Alignments> alignments=toCollection(record);
				s.addAll(alignments);
			}
		}
		iter.close();
		
		return s.iterator();
		
	}

	public Iterator<Alignments> getExonAlignmentsOverlappingRegion(IntervalTree<Alignment> tree){
		Set<Alignments> s=new TreeSet<Alignments>();
		Iterator<Node<Alignment>> iter=tree.iterator();
		while(iter.hasNext()){
			Alignment record=iter.next().getValue();
			if(passesQualAndStrandness(record)) {
				Collection<Alignments> alignments=toCollection(record);
				s.addAll(alignments);
			}
		}
		
		return s.iterator();
	}


	@Override
	public IntervalTree<Alignments> getFullIntervalTreeAsAlignments(String chr) throws IOException {
		if (this.fullAlignmentTreesAsAlignmentsCached == null || !this.fullAlignmentTreesAsAlignmentsCached.containsKey(chr) || this.fullAlignmentTreesAsAlignmentsCached.get(chr) == null || this.fullAlignmentTreesAsAlignmentsCached.get(chr).size() == 0) {
		 	IntervalTree<Alignments> tree = new IntervalTree<Alignments>();
		 	int start = 0;
		 	int end = this.getChromosomeLengths().get(chr);
			Alignments entireChrAsAlignments = new Alignments(chr, start, end);
		 	CloseableIterator<Alignment> iter = this.getAlignmentsOverlappingRegion(new Alignments(chr, start, end));
			while (iter.hasNext()) {
				Alignment record = iter.next();
				Alignments t1 = new Alignments(record.getChromosome(), record.getAlignmentStart(), record.getAlignmentEnd());
				if (t1.overlaps(entireChrAsAlignments)  ) {
					Node<Alignments> node = tree.find(record.getAlignmentStart(), record.getAlignmentEnd()+1);
					if (passesQualAndStrandness(record)) {
						if (node != null) {
							node.incrementCount();
						} else {
							tree.put(record.getAlignmentStart(), record.getAlignmentEnd()+1, t1);
						}
					}
				}
			}
			iter.close();
			
			if (this.fullAlignmentTreesAsAlignmentsCached == null)
				this.fullAlignmentTreesAsAlignmentsCached = new TreeMap<String, IntervalTree<Alignments>>();
			
			this.fullAlignmentTreesAsAlignmentsCached.put(chr, tree);
		}
		
		return this.fullAlignmentTreesAsAlignmentsCached.get(chr);
	}

	
	@Override
	public GeneCounts getGeneCounts(RefSeqGene gene, int extensionFactor) throws IOException {
		//IntervalTree<Alignment> alnTree = getIntervalTreeCached(gene.getChr(), gene.getStart(), gene.getEnd());
		
		double exon = getCountsPerAlignment(gene.getExons(),  extensionFactor);
		double intron = getCountsPerAlignment(gene.getIntronsBlocks(),  extensionFactor);
		double UTR3 = getCountsPerAlignment(gene.get3UtrIncludingIntrons(),  extensionFactor);
		double UTR5 = getCountsPerAlignment(gene.get5UtrIncludingIntrons(), extensionFactor);
		
		return new GeneCounts(exon, intron, UTR3, UTR5);
	}
	
	public String getModelFilePath() { return alignmentFileName;}

	/**
	public IntervalTree<Alignment> getIntervalTree(String chr, int start, int end){
	 	//System.err.print("Building tree from "+chr+" "+start+"-"+end + " - " );
		int counter=0;
	 	IntervalTree<Alignment> tree=new IntervalTree<Alignment>();
		CloseableIterator<Alignment> iter=getAlignmentsOverlappingRegion(new Alignments(chr, start, end));
		while(iter.hasNext()){
			Alignment record=iter.next();
			Alignments t1=new Alignments(record.getChromosome(), record.getAlignmentStart(), record.getAlignmentEnd());
			Alignments t2=new Alignments(chr, start, end);
			if(t1.overlaps(t2)){
				//System.err.println("Getting interval overlappers "+start+"-"+end+" "+record.getChr()+":"+record.getAlignmentStart()+"-"+record.getAlignmentEnd());
				
				Node<Alignment> node=tree.find(record.getAlignmentStart(), record.getAlignmentEnd()+1);
				if(record.getMappingQuality()>1){
					if(node!=null){node.incrementCount();}
					else{tree.put(record.getAlignmentStart(), record.getAlignmentEnd()+1, record);}
					counter++;
				}
			}
			//else{System.err.println("WARN: Doesnt overlap queried: "+chr+":"+start+"-"+end+" got "+record.getChromosome()+":"+record.getAlignmentStart()+"-"+record.getAlignmentEnd());}
		}
		iter.close();
		//System.err.println(tree.size()+" "+counter);
		return tree;
	}
	 * @throws IOException 
	*/
	
	@Override
	public IntervalTree<Alignment> getIntervalTree(String chr, int start, int end) throws IOException{
	 	//System.err.print("Building tree from "+chr+" "+start+"-"+end + " - " );
		int counter=0;
	 	IntervalTree<Alignment> tree=new IntervalTree<Alignment>();
		CloseableIterator<Alignment> iterReadsOverlappingRegion=getAlignmentsOverlappingRegion(new Alignments(chr, start, end));
		//System.err.println("interval: " + chr+":"+start+"-"+end+" free mem: " + Runtime.getRuntime().freeMemory());
		Alignments regionAsAlignments=new Alignments(chr, start, end);
		while(iterReadsOverlappingRegion.hasNext()){
			//if (counter % 1000000 == 1) System.out.print("getting " + counter);
			Alignment record=iterReadsOverlappingRegion.next();
			Alignments recordAsAlignments=new Alignments(record.getChromosome(), record.getAlignmentStart(), record.getAlignmentEnd());
			if(recordAsAlignments.overlaps(regionAsAlignments)){
				//System.err.println("Getting interval overlappers "+start+"-"+end+" "+record.getChr()+":"+record.getAlignmentStart()+"-"+record.getAlignmentEnd());
				
				Node<Alignment> node=tree.find(record.getAlignmentStart(), record.getAlignmentEnd()+1);
				if(passesQualAndStrandness(record) ){
					if(node!=null){node.incrementCount();}
					else{tree.put(record.getAlignmentStart(), record.getAlignmentEnd()+1, record);}
					counter++;
				}
			}
			//else{System.err.println("WARN: Doesnt overlap queried: "+chr+":"+start+"-"+end+" got "+record.getChromosome()+":"+record.getAlignmentStart()+"-"+record.getAlignmentEnd());}
		}
		iterReadsOverlappingRegion.close();
		//System.err.println("end interval free mem: " + Runtime.getRuntime().freeMemory() + " tree size " + tree.size()+" counter "+counter);
		return tree;
	}
	
	
	@Override
	public IntervalTree<Alignment> getIntervalTreeCached(String chr, int start, int end) throws IOException{
		//System.err.println("Cached interval tree called for: " + chr +":" + start + "-" + end + " --- chunks start-end " + chunkStart + "-" + chunkEnd);
		if(!chr.equalsIgnoreCase(this.chunkChr)){resetTreeCache(chr);}
		//if end>chunkEnd (includes havent started)
		int newEnd = 0;
		int newStart = 0;
		if(end>this.chunkEnd){
			if(start > this.chunkEnd) {
				newStart = start;
				newEnd = Math.max(end, start + this.chunkSize);
			} else {
				newEnd=Math.max(end, this.chunkEnd+this.chunkSize);
				newStart=Math.min(start, this.chunkEnd);
			}
			//System.err.println("Getting new 1 "+newStart+"-"+newEnd+" "+start+"-"+end);
			this.chunkAlignmentTree=getIntervalTree(chr, newStart, newEnd);
			this.chunkEnd=newEnd;
			this.chunkStart=newStart;
		}
		if(start<this.chunkStart){//Should this be else if?
			if(end > this.chunkStart) {
				newEnd = end;
				newStart = Math.min(start, end - chunkSize);
			} else {
				newStart=Math.min(start, chunkStart-this.chunkSize);
				newEnd=Math.max(end, chunkEnd);
			}
			//System.err.println("Getting new 2 "+newStart+"-"+newEnd+" "+start+"-"+end);
			this.chunkAlignmentTree=getIntervalTree(chr, newStart, newEnd);
			this.chunkEnd=newEnd;
			this.chunkStart=newStart;
		}
		
		//System.err.println(chunkStart+"-"+chunkEnd);
		
		return this.chunkAlignmentTree;
	}
	

	// What does truncated mean? JE
	public IntervalTree<Alignment> getIntervalTreeTruncated(String chr, int start, int end) throws IOException{
	 	//System.err.print("Building tree from "+chr+" "+start+"-"+end + " - " );
		Alignments region=new Alignments (chr, start, end);
		int counter=0;
	 	IntervalTree<Alignment> tree=new IntervalTree<Alignment>();
		CloseableIterator<Alignment> iter=getAlignmentsOverlappingRegion(new Alignments(chr, start, end));
		while(iter.hasNext()){
			Alignment record=iter.next();
			
			Alignments readAlign=new Alignments(record.getChromosome(), record.getAlignmentStart(), record.getAlignmentEnd());
			//Node<Alignment> node=tree.find(record.getAlignmentStart(), record.getAlignmentStart()+1);
			if(passesQualAndStrandness(record) ){
				if(readAlign.overlaps(region)){
					//System.err.println("Skipping low quality reads "+readAlign.toUCSC());
					//if(node!=null){node.incrementCount();}
					tree.put(record.getAlignmentStart(), record.getAlignmentStart()+1, record);
					tree.put(record.getAlignmentEnd()-1, record.getAlignmentEnd(), record);
					counter++;
				}
			}
		}
		iter.close();
		//System.err.println(tree.size()+" "+counter);
		return tree;
	}

	
	//TODO Fix cache should only do min/max if overlaping current interval otherwise just jump
	@Override
	public IntervalTree<Alignment> getIntervalTreeTruncatedCached(String chromosome, int start, int end) throws IOException {
				
		if(!chromosome.equalsIgnoreCase(this.chunkChr)){resetTreeCache(chromosome);}
		//if end>chunkEnd (includes havent started)
		if(end>this.chunkEnd){
			int newStart=start;			
			//if overlaps do this
			Alignments region=new Alignments(chromosome, start, end);
			if(region.overlaps(new Alignments(chromosome, chunkStart, chunkEnd))){
				newStart=Math.min(start, this.chunkEnd);
			}
			
			int newEnd=Math.max(end, this.chunkEnd+this.chunkSize);
			System.err.println("Updating cache 1 "+chromosome+":"+start+"-"+end+" new cache "+newStart+"-"+newEnd+" old cache "+chunkStart+"-"+chunkEnd);
			this.chunkAlignmentTree=getIntervalTreeTruncated(chromosome, newStart, newEnd);
			this.chunkEnd=newEnd;
			this.chunkStart=newStart;
		}
		if(start<this.chunkStart){//Should this be else if?
			int newStart=Math.min(start, chunkStart-this.chunkSize);
			int newEnd=Math.max(end, chunkEnd);
			System.err.println("Updating cache 2 "+chromosome+":"+start+"-"+end+" new cache "+newStart+"-"+newEnd+" old cache "+chunkStart+"-"+chunkEnd);
			this.chunkAlignmentTree=getIntervalTreeTruncated(chromosome, newStart, newEnd);
			this.chunkEnd=newEnd;
			this.chunkStart=newStart;
		}
		
		return this.chunkAlignmentTree;
	}

	private Collection<Alignments> getIntrons(Alignment record){
		Collection<Alignments> exons=this.toCollection(record);
		RefSeqGene gene=new RefSeqGene(exons);
		Collection<Alignments> rtrn=gene.getIntronSet();
		return rtrn;
	}

	@Override
	public double getMinimumMappingQuality() { return minMappingQuality;}


	@Override
	public Alignments getNextExon(String chr) throws IOException{
		Alignment align = null;
		boolean isSpliced = true;
		while(hasNextExon(chr) && isSpliced) {
			align = readIterator.next();
			if( align.getAlignmentBlocks()==null || align.getAlignmentBlocks().length==1){
				isSpliced = false;
			}			
		}
	
		 return align == null ? null : new Alignments(align.getChromosome(), align.getAlignmentStart(), align.getAlignmentEnd());
	}	
	
	//get the number of spliced reads overlapping a given region
	@Override
	public double getNumberOfSplicedReads(Alignments align) throws IOException{
		CloseableIterator<Alignment> iter =  reader.query(align.getChromosome(), align.getStart(), align.getEnd(), false);
		
		double counter=0;
		while (iter.hasNext()) {
			boolean count=false;
			int blocks=0;
	        Alignment record = iter.next();
			if(!passesQualAndStrandness(record)) {
	    		continue;
	    	}
	        for(AlignmentBlock b : record.getAlignmentBlocks()) {
	           	Alignments exon=new Alignments(record.getChromosome(), b.getStart(), b.getEnd());
	        	if(exon.overlapsAtAll(align)){count=true;}
	        	blocks++;
	        }
	        if(count && blocks>1){counter = counter + countReads(record);}
	    }
	    iter.close();
	    return normalizer == null ? counter : normalizer.normalize(counter);
	}

	public double getNumberOfSplicedReads(){
		CloseableIterator<Alignment> iter =  reader.iterator();
		
		double splicedReads=0;
		while (iter.hasNext()) {
			Alignment record = iter.next();
	        if(passesQualAndStrandness(record)){
	        	splicedReads++;
	        }
	    }
	    iter.close();
	    return splicedReads;
	}

	//TODO Get this working
	@Override
	public Collection<Alignments> getOverlappingRegions(String chr) throws IOException{
		Collection<Alignments> rtrn=new TreeSet<Alignments>();
		Collection<Alignments> splicedExons=new TreeSet<Alignments>();
		Alignments chrRegion=new Alignments(chr, 0, this.chromosomeSizes.get(chr));
		CloseableIterator<Alignment> reads =  reader.query(chrRegion.getChromosome(), chrRegion.getStart(), chrRegion.getEnd(), false);
				
		Alignments current=null;
		while(reads.hasNext()){
			Collection<Alignments> alignments=toCollection(reads.next());
			while(alignments.size()>1){splicedExons.addAll(alignments); alignments=toCollection(reads.next());}
			
				if(current==null){current=alignments.iterator().next();}
				Alignments exonic=current;
				while(reads.hasNext()){
					alignments=toCollection(reads.next());
					while(alignments.size()>1){splicedExons.addAll(alignments); alignments=toCollection(reads.next());}
					
						Alignments next=alignments.iterator().next();;
						if(current.overlapsAtAll(next)){
							exonic=exonic.merge(next);
							current=next;
						}
						else{current=next; rtrn.add(exonic); rtrn.add(current); break;}
					
				}
				if(exonic!=null){rtrn.add(exonic);}
			
		}
		reads.close();
		rtrn.addAll(splicedExons);
		rtrn=CollapseByIntersection.collapseByIntersection(rtrn, false);
		return rtrn;
	}

	@Override
	public CloseableIterator<Alignment> getReadIterator(){
		return this.reader.iterator();
	}
	
	@Override
	public CloseableIterator<Alignment> getReadIterator(Alignments region) throws IOException{
		return this.reader.query(region.getChr(), region.getStart(), region.getEnd(), false);
	}
	
	public SAMFileHeader getSamHeader () throws IOException{
		return this.reader.getHeader();
	}
	
	
	//optimized for iterative scoring
	//keep track of previous scores
	@Override
	public double getScorePerAlignmentFromCache(RefSeqGene window, IntervalTree chunkAlignmentTree, int extensionFactor){
		//keep a global tree of scores
		//for a new gene first check if overlaps anything previously seen, if so pull the previous score and compute the unique pieces
		
		IntervalTree<Alignments> cachingTree=this.cachingScoringTrees.get(window.getChr());
		Set<Alignments> exons=window.getExonSet();
		double score=0;
		for(Alignments exon: exons){
			//get closest to perfect overlapper with this one
			Node<Alignments> maxNode=cachingTree.max(exon.getStart(), exon.getEnd()); //if null then needs to compute from scratch and store
			//compute score of exon given that we know maxNode
			double exonScore=computeNumberOfAlignmentsBasedOnCache(exon, maxNode, chunkAlignmentTree, extensionFactor);
			score+=exonScore;
			
			//remove the previous
			if(maxNode!=null){cachingTree.remove(maxNode.getStart(), maxNode.getEnd());}
			
			//add the new score to the exon
			exon.setCountScore(exonScore);
			
			//add the exon and its score to the cache
			cachingTree.put(exon.getStart(), exon.getEnd(), exon);	
		}
		
		this.cachingScoringTrees.put(window.getChr(), cachingTree);
				
		//TODO Optimize memory efficiency by purging cache
		
		return score;
	}

	@Override
	public Collection<Alignments> getSplicedReadExons(String chr) throws IOException{
		Collection<Alignments> regions=new TreeSet();
		Alignments chrRegion=new Alignments(chr, 0, this.getChromosomeLengths().get(chr));
		this.readIterator=reader.query(chrRegion.getChromosome(), chrRegion.getStart(), chrRegion.getEnd(), false);
		while(this.readIterator.hasNext()){
			Collection<Alignments> exons=toCollection(this.readIterator.next());
			if(exons!=null && exons.size()>1){regions.addAll(exons);}
		}
		this.readIterator.close();
		return regions;
	}
	
	public  void setExtensionFactor(int extensionFactor) { this.extensionFactor = extensionFactor;}

	@Override
	public Map<Alignments, Integer> getSplicedReads(Alignments align) throws IOException{
		Collection<ReadFilter> list=new ArrayList<ReadFilter>();
		return getSplicedReads(align, list);
	}

	@Override
	public Map<Alignments, Integer> getSplicedReads(Alignments align, final int minIntronSize, final int maxIntronSize , int minNumIntrons, Sequence chrSeq) throws IOException{
		Collection<ReadFilter> list=new ArrayList<ReadFilter>();
		list.add(new ReadFilter(){
			@Override
			public boolean passes(Alignments read) {
				return read.length() > minIntronSize && read.length()<maxIntronSize;
			}
	    	
	    });
		return getSplicedReads(align, list, minNumIntrons, chrSeq);
	}

	@Override
	public Map<Alignments, Integer> getSplicedReads(Alignments align, final int minIntronSize, final int maxIntronSize) throws IOException{
	
	    return getSplicedReads( align, new ReadFilter(){
			@Override
			public boolean passes(Alignments read) {
				return read.length() > minIntronSize && read.length()<maxIntronSize;
			}
	    	
	    });
	}

	@Override
	public Map<Alignments, Integer> getSplicedReads(Alignments region, ReadFilter filter) throws IOException {
	    List<ReadFilter> filters = new ArrayList<ReadFilter>(1);
	    filters.add(filter);
	    return getSplicedReads(region, filters);
	}

	@Override
	public Map<Alignments, Integer> getSplicedReads(Alignments region, Collection<ReadFilter> filters, int minNumIntrons) throws IOException{
		return getSplicedReads(region, filters, minNumIntrons, null);
	}

	@Override
	public Map<Alignments, Integer> getSplicedReads(Alignments region, Collection<ReadFilter> filters, int minNumIntrons, Sequence chrSeq) throws IOException{
		Map<Alignments, Integer> rtrn=new TreeMap<Alignments, Integer>(); 
		
		
		CloseableIterator<Alignment> iter =  getAlignmentsOverlappingRegion(region);
		while (iter.hasNext()) {
	       	Alignment record = iter.next();
	       	if(passesQualAndStrandness(record) ){
		       	Collection<Alignments> introns=getIntrons(record);
		       	for(Alignments intron: introns){
		       		boolean passes=true;
		       		if(chrSeq!=null){
		       			String orientation=GeneTools.orientationFromSpliceSites(intron, chrSeq);
		       			intron.setOrientation(orientation);
		       			if(orientation.equalsIgnoreCase("*")){passes=false;}
		       		}
		       		Iterator<ReadFilter> filterIt = filters.iterator();
		       		while(passes && filterIt.hasNext()) {
		       			passes = filterIt.next().passes(intron);
		       		}
		       		
		       		if(passes) {
			       		int num=0;
			       		if(rtrn.containsKey(intron)){num=rtrn.get(intron);}
			       		num++;
			       		rtrn.put(intron, num);
		       		}
		       	}	       	
	       	}
	    }
		
		Map<Alignments, Integer> newRtrn=new TreeMap<Alignments, Integer>();
		for(Alignments intron: rtrn.keySet()){
		if(rtrn.get(intron)>minNumIntrons){newRtrn.put(intron, rtrn.get(intron));}
		}
		
	    iter.close();
	    return newRtrn;
	}

	@Override
	public Map<Alignments, Integer> getSplicedReads(Alignments region, Collection<ReadFilter> filters) throws IOException{
		return getSplicedReads(region, filters, 0);
	}

	public Map<Alignments, Integer> getSplicedReads(Alignments region, Collection<ReadFilter> filters, Sequence chrSeq) throws IOException{
		return getSplicedReads(region, filters, 0, chrSeq);
	}

	@Override
	public Map<Alignments, Integer> getSplicedReads(int minIntronSize, int maxIntronSize) throws IOException{
		Map<Alignments, Integer> rtrn=new TreeMap<Alignments, Integer>(); 
		for(String chr: this.chromosomeSizes.keySet()){
			//System.err.println(chr);
			Alignments align=new Alignments(chr, 0, this.chromosomeSizes.get(chr));
			rtrn.putAll(getSplicedReads(align,  minIntronSize, maxIntronSize));
		}
	    return rtrn;
	}

	@Override
	public double getSpliceWeightFactor(){return this.spliceWeightFactor;}
	
	
	
	//public int getFirstReads(String chr){return this.firstReads.get(chr);}
	
	public double[] getTotalNumberOfMappedReads(){
		CloseableIterator<Alignment> iter=reader.iterator();
		double numReads=0;
		double numMapped=0;
		while(iter.hasNext()){
			Alignment a=iter.next();
			if(a.getMappingQuality() > minMappingQuality) {
				numReads++;
				if (a.isMapped()) {numMapped++;}
			}
		}
		
		iter.close();
		double [] res= new double[2];
		res[0]=numReads; res[1]=numMapped;
		return res;
	}
	
	public double[] getTotalNumberOfMappedStrandedReads(){
		CloseableIterator<Alignment> iter=reader.iterator();
		double numReads=0;
		double numMapped=0;
		while(iter.hasNext()){
			Alignment a=iter.next();
			numReads++;
			if (passesQualAndStrandness(a) && a.isMapped() ) {numMapped++;}
		}
		
		iter.close();
		double [] res= new double[2];
		res[0]=numReads; res[1]=numMapped;
		return res;
	}
	
	
	
	@Override
	public double getTotalNumberOfStrandedReads() {
		CloseableIterator<Alignment> iter=reader.iterator();
		double numReads=0;
		while(iter.hasNext()){
			Alignment record = iter.next();
			if(passesQualAndStrandness(record)){
				numReads++;
			}
		}
		iter.close();
		return numReads;
	}

	//public int getFirstReads(String chr){return this.firstReads.get(chr);}
	
	public double getTotalNumberOfReads(){
		CloseableIterator<Alignment> iter=reader.iterator();
		double numReads=0;
		while(iter.hasNext()){
			Alignment record = iter.next();
			if(record.getMappingQuality() > minMappingQuality) {
				numReads++;
			}
		}
		
		iter.close();
		return numReads;
	}

	public double getWeight() throws IOException{
		CloseableIterator<Alignment> iter=findLargestChromsomeAlignmentIterator();
		double numReads=0;
		double splicedReads=0;
		while(iter.hasNext()){
			Alignment record=iter.next();
			if(record.getMappingQuality() > minMappingQuality) {
				if(passesQualAndStrandness(record)){splicedReads++;}
				numReads++;
			}
		}
		
		iter.close();
		return numReads/splicedReads;
	}

	@Override
	public void setChunkSize(int chunkSize){this.chunkSize=chunkSize;}

	@Override
	public void setMinimumMappingQuality(double minimumMappingQual) { this.minMappingQuality = minimumMappingQual;}

	@Override
	public void setNegativeStranded() {
		this.isStranded = true;
		this.isNegativeStrand = true;
	}

	public void setSecondRead() {
		this.isPaired = true;
		this.isSecondRead = true;
	}

	public void setFirstRead() {
		this.isPaired = true;
		this.isSecondRead = false;
	}

	
		
		@Override
		public void setNormalizationStrategy(ReadCountNormalizer normalizer) {
			this.normalizer = normalizer;
		}

	@Override
	public void setPositiveStranded() {
		this.isStranded = true;
		this.isNegativeStrand = false;
	}

	private double setSpliceWeight(boolean upweightSplices) throws IOException {
		logger.debug("Computing weights..... upweighting? " + upweightSplices + " weight: ");
		double weight= upweightSplices ? getWeight() : 1;
		logger.debug(weight);
		return weight;
	}

		
		//Given a specific region, checks whether all spliced reads apn the same exact cooridinates and sequence (Moran 5/4/2010)	
		public boolean AreSplicedReadsIdentical(Alignments region, Collection<ReadFilter> filters,  Sequence chrSeq) throws IOException{
			
			
			boolean res=true;
			ArrayList<Alignment> splicedReads= new ArrayList<Alignment>();
			
			CloseableIterator<Alignment> iter =  getAlignmentsOverlappingRegion(region);
			while (iter.hasNext()) {
		       	Alignment record = iter.next();
		       	if(record.getMappingQuality()>minMappingQuality && (!isNegativeStrand || record.isNegativeStrand() == isNegativeStrand)){
			       	Collection<Alignments> introns=getIntrons(record);
			       	for(Alignments intron: introns){
			       		boolean passes=true;
			       		if(chrSeq!=null){
			       			String orientation=GeneTools.orientationFromSpliceSites(intron, chrSeq);
			       			intron.setOrientation(orientation);
			       			if(orientation.equalsIgnoreCase("*")){passes=false;}
			       		}
			       		if (! (filters ==null)){
				       		Iterator<ReadFilter> filterIt = filters.iterator();
				       		while(passes && filterIt.hasNext()) {
				       			passes = filterIt.next().passes(intron);
				       		}
			       		}
			       		
			       		if(passes) {
				       		
				       		if( ! splicedReads.contains(record)){
				       			splicedReads.add(record);
				       		}
			       		}
			       	}	       	
		       	}
		    }
			iter.close();
			
			if (splicedReads.size()<2) return false;
			
			for(int i=0; i <splicedReads.size(); i++){
				for (int j=i; j <splicedReads.size(); j++) {
					
					Alignment a =splicedReads.get(i);
					Alignment b =splicedReads.get(j);
					if (a.getEnd()!= b.getEnd() || a.getStart() != b.getStart() || ! a.getCigarString().equals(b.getCigarString()) || ! a.getReadSequence().equals(b.getReadSequence()))
					{
						res=false; return res;
					}
				}
			}
			
		  
		    return res;
		}

	@Override
	public void clearFullIntervalTreeAsAlignmentsCached(String chr) {
		if (this.fullAlignmentTreesAsAlignmentsCached == null)
			return;
		
		this.fullAlignmentTreesAsAlignmentsCached.remove(chr);
	}

	private double computeNumberOfAlignmentsBasedOnCache(Alignments exon, Node<Alignments> node, IntervalTree chunkAlignmentTree, int extensionFactor){
		//if doesnt exist then just score from scratch
		if(node==null){return getCountsPerAlignment(exon, chunkAlignmentTree, extensionFactor);}
		Alignments previousAlignment=node.getValue();
		//if doesnt overlap at all then just do the stupid thing and compute the score directly from the alignments
		if(!previousAlignment.overlapsAtAll(exon) && exon.overlapsAtAll(previousAlignment)){return getCountsPerAlignment(exon, chunkAlignmentTree, extensionFactor);}
		if(previousAlignment.equals(exon)){return previousAlignment.getCountScore();}
		
		//if cached interval starts before exon and ends after exon then compute from scratch and remove from tree
		//TODO Make sure this is efficient perhaps try keeping it and computing off of it
		if(previousAlignment.getStart()<exon.getStart() && previousAlignment.getEnd()>exon.getEnd()){
			double score= getCountsPerAlignment(exon, chunkAlignmentTree, extensionFactor);
			this.cachingScoringTrees.get(previousAlignment.getChr()).remove(previousAlignment.getStart(), previousAlignment.getEnd());
			return score;
		}
		
		double previousScore=previousAlignment.getCountScore();
		Alignments startPosition = new Alignments(exon.getChr(), previousAlignment.getStart(), exon.getStart()); //is positions not in exon but in node
		double subtractVal=getCountsOfUniqueOverlappers(startPosition, exon, chunkAlignmentTree, extensionFactor); //this computes the values to subtract from the score
		Alignments endPosition = new Alignments(exon.getChr(), previousAlignment.getEnd()+1, exon.getEnd()+1); //this is the positions in exon but not in node
		double addVal=getCountsOfUniqueOverlappers(endPosition, previousAlignment, chunkAlignmentTree, extensionFactor); //this computes the values to add to the score
		
		double score=((previousScore-subtractVal)+addVal);
		return score;
	}

	private Alignments extendRead(Alignments exon, boolean negStrand, int EF){
		Alignments rtrn=new Alignments(exon.getChromosome(), exon.getStart(), exon.getEnd()+EF);
		if(negStrand){rtrn=new Alignments(exon.getChromosome(), exon.getStart()-EF, exon.getEnd());}
		return rtrn;
	}

	private Alignments extendRead(Alignment exon, boolean negStrand, int EF){
		Alignments rtrn=new Alignments(exon.getChromosome(), exon.getAlignmentStart(), exon.getAlignmentEnd()+EF);
		if(negStrand){rtrn=new Alignments(exon.getChromosome(), exon.getAlignmentStart()-EF, exon.getAlignmentEnd());}
		return rtrn;
	}

	/**
	 * Returns the largest chromosomes THAT HAS ALIGNMENTS (changed by manuel 01/24/10
	 * @return
	 * @throws IOException 
	 */
	private CloseableIterator<Alignment> findLargestChromsomeAlignmentIterator() throws IOException {
		List<String> chrList = new ArrayList<String>(chromosomeSizes.keySet());
		Collections.sort(chrList, new Comparator<String>(){
	
			@Override
			public int compare(String o1, String o2) {
				return chromosomeSizes.get(o2) - chromosomeSizes.get(o1);
			}
			
		});
		Iterator<String> chrIt = chrList.iterator();
		CloseableIterator<Alignment> largest = null;
		while(largest == null && chrIt.hasNext()) {
			String candidateChr = chrIt.next();
			CloseableIterator<Alignment> iter=reader.query(candidateChr, 0, chromosomeSizes.get(candidateChr), true);
			if(iter.hasNext()) {
				largest = iter;
			} else {
				iter.close();
			}
		}
		
		return largest;
	}

	private Map<String, Integer> firstReads() throws IOException{
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		for(String chr: this.chromosomeSizes.keySet()){
			logger.info(chr+" "+this.chromosomeSizes.get(chr));
			CloseableIterator<Alignment> iter=getAlignmentsOverlappingRegion(new Alignments(chr, 0, this.chromosomeSizes.get(chr)));
			if(iter.hasNext()){Alignment align=iter.next(); rtrn.put(chr, align.getAlignmentStart());}
			iter.close();
		}
		
		return rtrn;
	}

	@Override
	public boolean hasNextExon(String chr) throws IOException{
		if(!this.startedIteration || !this.chunkChr.equalsIgnoreCase(chr)){
			this.chunkChr=chr;
			Alignments chrRegion=new Alignments(chr, 0, this.getChromosomeLengths().get(chr));
			this.readIterator=reader.query(chrRegion.getChromosome(), chrRegion.getStart(), chrRegion.getEnd(), false);
		}
		this.startedIteration=true;
		return this.readIterator.hasNext();
	}

	private Map<String, IntervalTree<Alignments>> initializeCachingTree(){
		Map<String, IntervalTree<Alignments>> rtrn=new TreeMap<String, IntervalTree<Alignments>>();
		
		for(String chr: chromosomeSizes.keySet()){
			rtrn.put(chr, new IntervalTree<Alignments>());
		}
		
		return rtrn;
	}

	@Override
	public boolean isPositiveStranded() {return isStranded && !isNegativeStrand;}

	@Override
	public boolean isNegativeStranded() {return isStranded && isNegativeStrand;}

	@Override
	public boolean isStranded() { return isStranded;}

	private Map<String, IntervalTree<Alignments>> makeTree(Collection<Alignments> alignments){
		Map<String, IntervalTree<Alignments>> rtrn=new TreeMap<String, IntervalTree<Alignments>>();
		
		for(Alignments align: alignments){
			IntervalTree<Alignments> tree=new IntervalTree<Alignments>();
			String chr=align.getChr();
			if(rtrn.containsKey(chr)){tree=rtrn.get(chr);}
			//if(align.getStart()>=align.getEnd()){System.err.println("ERROR: " +align.toUCSC());}
			tree.put(align.getStart(), align.getEnd()+1, align);
			rtrn.put(align.getChr(), tree);
		}
		
		return rtrn;
	}

	@Override
	public boolean passes(RefSeqGene gene, IntervalTree<Alignment> tree, int EF, double cutoff){
		//double sum=0;
		//Alignments[] exons=gene.getExons();
		long start=System.currentTimeMillis();
		double counter=0;
		
		Collection<Alignments> exons=gene.getExonSet();
		
		//TODO We can do a fast pass through where we grab the exons and sum if less than critVal then return false else sum correctly
		
		for(Alignments exon: exons){
			//get Alignments over the whole region
			Iterator<Node<Alignment>> iter=tree.overlappers(exon.getStart(), exon.getEnd()+1); 
			while(iter.hasNext()){
				Node<Alignment> node=iter.next();
				int num=node.getNumReplicates();
				Alignment alignment=node.getValue();
				if(cutoff==0 && getCountsPerAlignment(alignment, exon, EF)>0){return true;}
			}
		}
		if(cutoff==0){return false;}
		
		Iterator<Node<Alignment>> iter=tree.overlappers(gene.getStart(), gene.getEnd()+1); 
		int i=0;
		while(iter.hasNext()){
			Node<Alignment> node=iter.next();
			int num=node.getNumReplicates();
			Alignment alignment=node.getValue();
			//only count those overlapping one of the exons
			counter+=(getCountsPerAlignment(alignment, gene, EF)*num);
			if(counter>cutoff){return true;}
			i++;
		}
		return false;
	}

	//get the number of alignments overlapping a given region with a cached interval tree
	@Override
	public long passesCutoff(Alignments align, IntervalTree<Alignment> tree, int EF, double threshold){
		double counter=0;
		
		long start=System.currentTimeMillis();
		Iterator<Node<Alignment>> iter=tree.overlappers(align.getStart(), align.getEnd());
		long end=System.currentTimeMillis();
		
		while(iter.hasNext()){
			Node<Alignment> node=iter.next();
			int num=node.getNumReplicates();
			counter+=(getCountsPerAlignment(node.getValue(), align, EF)*num);
		}
		return end-start;
		
	}

	private boolean passesQualAndStrandness(Alignment record) {
		return record.getMappingQuality() > minMappingQuality  && (!isStranded || record.isNegativeStrand() == isNegativeStrand ) && record.isMapped() && !record.isDuplicate() && (!isPaired || record.isSecondOfPair() == isSecondRead);
	}

	@Override
	public void resetTreeCache(String chr){
		//System.err.println(chr);
		this.chunkChr=chr;
		this.chunkAlignmentTree=new IntervalTree<Alignment>();
		this.chunkStart=0;
		this.chunkEnd=0;
	}
	
	@Override
	public void resetTreeCache(){
		this.chunkChr="";
		this.chunkStart=0;
		this.chunkEnd=0;
	}
	
	
	@Override
	public void restartIterator(){this.startedIteration=false;}

	private Collection<Alignments> toCollection(Alignment record){
		Collection<Alignments> rtrn=new ArrayList<Alignments>();
		AlignmentBlock[] blocks=record.getAlignmentBlocks();
		if(blocks==null){rtrn.add(new Alignments(record.getChromosome(), record.getAlignmentStart(), record.getAlignmentEnd())); return rtrn;}
		for(int i=0; i<blocks.length; i++){
			Alignments exon=new Alignments(record.getChromosome(), blocks[i].getStart(), blocks[i].getEnd());
			rtrn.add(exon);
		}
		
		return rtrn;
	}

	@Override
	public void unsetStranded()  { isStranded = false;}

	private Alignments updateEnd(Alignments exon, IntervalTree tree){
		Iterator<Node<Collection<Alignment>>> iter=tree.overlappers(exon.getStart(), exon.getEnd());
		boolean updated=false;
		int locationEnd=-99;
		
		while(iter.hasNext()){
			Alignment record=iter.next().getValue().iterator().next();
			
			if(record.getAlignmentBlocks()!=null){
				for(AlignmentBlock b : record.getAlignmentBlocks()) {
	        		Alignments read=new Alignments(record.getChromosome(), b.getStart(), b.getEnd());
	        		if(read.overlapsAtAll(exon)){
	        			if(locationEnd>=0){locationEnd=Math.max(locationEnd, read.getEnd());}
	    				else{locationEnd=read.getEnd();}
	        			updated=true;
	        		}
	        	}
			}
			else{
				if(locationEnd>=0){locationEnd=Math.max(locationEnd, record.getAlignmentEnd());}
				else{locationEnd=record.getAlignmentEnd();}
				updated=true;
			}
		}
		
		
		if(updated){
			Alignments updatedExon=new Alignments(exon.getChr(), exon.getStart(), Math.min(exon.getEnd(), locationEnd));
			//System.err.println("End "+exon.toUCSC()+" "+updatedExon.toUCSC()+" "+locationEnd);
			return updatedExon;
		}
		return exon;
	}

	//TODO Currently only works if there is some read within the first and last exons. Need to return to the case of where the first read occurs in some other exon
	//TODO idea would be to first count numbers for each exon and then make new gene starting from first non-0 exon and ending at last non-0 exon
	@Override
	public RefSeqGene updateGeneByFirstCounts(RefSeqGene gene, IntervalTree<Alignment> tree, int EF){
		double sum=0;
		Collection<Alignments> exons=gene.getSortedAndUniqueExons();
		Collection<Alignments> truncated=new TreeSet();
		
		int i=0;
		for(Alignments exon: exons){
			if(i==0){
				Alignments updated=updateStart(exon, tree);
				truncated.add(updated);
			}
			else if(i==(exons.size()-1)){
				Alignments updated=this.updateEnd(exon, tree);
				truncated.add(updated);
			}
			else{truncated.add(exon);}
			i++;
		}
		
		RefSeqGene updated=new RefSeqGene(truncated);
		
		return updated;
	}

	private Alignments updateStart(Alignments exon, IntervalTree tree){
		Iterator<Node<Collection<Alignment>>> iter=tree.overlappers(exon.getStart(), exon.getEnd());
		boolean updated=false;
		int locationStart=-99;
		
		while(iter.hasNext()){
			Alignment record=iter.next().getValue().iterator().next();
			
			if(record.getAlignmentBlocks()!=null){
				for(AlignmentBlock b : record.getAlignmentBlocks()) {
	        		Alignments read=new Alignments(record.getChromosome(), b.getStart(), b.getEnd());
	        		if(read.overlapsAtAll(exon)){
	        			if(locationStart>=0){locationStart=Math.min(locationStart, read.getStart());}
	    				else{locationStart=read.getStart();}
	        			updated=true;
	        		}
	        	}
			}
			else{
				if(locationStart>=0){locationStart=Math.min(locationStart, record.getAlignmentStart()); }
				else{locationStart=record.getAlignmentStart();}
				updated=true;
			}
		}
		
		
		if(updated){
			Alignments updatedExon=new Alignments(exon.getChr(), Math.max(exon.getStart(), locationStart), exon.getEnd());
			//System.err.println("Start "+exon.toUCSC()+" "+updatedExon.toUCSC()+" "+locationStart);
			return updatedExon;
		}
		return exon;
		
	}

	
	public static void main(String[] args) throws IOException{
		if(args.length>3){
		String file=args[0];
		String sizes=args[1];
		Alignments region=new Alignments(args[2]);
		String save=args[3];
		
		AlignmentDataModel data=new GenericAlignmentDataModel(file, sizes);
		CloseableIterator<Alignment> iter=data.getAlignmentsOverlappingRegion(region);
		
		FileWriter writer=new FileWriter(save);
		while(iter.hasNext()){
			Alignment align=iter.next();
			writer.write(align.getChromosome()+"\t"+align.getStart()+"\t"+align.getEnd()+"\t"+align.getReadName()+"\n");
		}
		iter.close();
		writer.close();
		}
		else{System.err.println(usage);}
	}
	

	public void setRemoveDuplicatesFlag(boolean flag){
		
		this.removeDuplicatesFlag = flag;
	}
	
	/**
	 * @author: skadri
	 * This function will set the weighReadCounts flag which decides how the reads are counted for each alignment
	 * @param flag: true if the read counts should be weighed by the value of the NH flag 
	 */
	public void setWeighReadCountsFlag(boolean flag){
		
		if (flag){
			// Set minimum mapping quality to -1 since we will be weighing reads
			this.minMappingQuality=-1.0;
		}
		this.weighReadCounts = flag;
	}
	
	/**
	 * @author: skadri
	 * This function will calculate the count value for an alignment read record based on the alue of the weigh ReadCounts flag
	 * @param record
	 * @return Count value for alignment
	 */
	public double countReads(Alignment record){
		
		if(weighReadCounts && record.getAttribute("NH") != null){
			return ((double)1/(Double.valueOf((record.getAttribute("NH")).toString())));
		}
		else{
			return (1.0);
		}
	}
	
	/**
	 * An empty iterator to avoid returning null 
	 * @author mgarber
	 *
	 */
	public static class EmptyCloseableIterator implements CloseableIterator<Alignment> {

		@Override
		public void remove() {
			// Nothing to do

		}

		@Override
		public Alignment next() {
			throw new NoSuchElementException("Empty iterator, no elements for this iterator");
		}

		@Override
		public boolean hasNext() {
			return false;
		}

		@Override
		public void close() {
			// Nothing to do

		}
	}
	
	
	
}
