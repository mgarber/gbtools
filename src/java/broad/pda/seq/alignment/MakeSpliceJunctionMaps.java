package broad.pda.seq.alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.AlignmentBlock;
import org.broad.igv.sam.reader.AlignmentQueryReader;
import org.broad.igv.sam.reader.BAMQueryReader;

import broad.core.annotation.PSL;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.GeneTools;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.alignment.sam.SAMUtils;
import net.sf.samtools.util.CloseableIterator;

public class MakeSpliceJunctionMaps {

	int readLength=76;
	boolean filterBySplice=true;
	private boolean repeatMask=false;
	int minIntronSize=10;
	int maxIntronSize=1000000;
	double minPercentMapped=.9;
	int minCoverageOnEachSide=10;
	
	//TODO: Try filtering BLAT results first
	//TODO: Run mdust
	//TODO: Only get unique BLAT hits
	
	public MakeSpliceJunctionMaps(Map<String, Collection<RefSeqGene>> genesByChr, String save, String genomeDir, int minIntron, int maxIntron) throws IOException{
		//Make pairwise junctions
		Map<String, Collection<RefSeqGene>> junctionsByChr=new TreeMap<String, Collection<RefSeqGene>>();
		
		for(String chr: genesByChr.keySet()){
			System.err.println(chr);
			Collection<RefSeqGene> genes=genesByChr.get(chr);
			Collection<RefSeqGene> junctions=new TreeSet();
			for(RefSeqGene gene: genes){
				Collection<RefSeqGene> junc=getJunctions(gene, readLength, minIntron, maxIntron);
				junctions.addAll(junc);
			}
			junctionsByChr.put(chr, junctions);
		}
		
		FileWriter bedWriter=new FileWriter(save+".bed");
		FileWriter seqWriter=new FileWriter(save+".fa");
		
		//I want to get the concatanated sequence and see how many more I get with just junctions and the full thing
		writeAll(bedWriter, seqWriter, genesByChr, genomeDir, "fragment");
		writeAll(bedWriter, seqWriter, junctionsByChr, genomeDir, "junction");
		
		bedWriter.close();
		seqWriter.close();
	}
	
	public MakeSpliceJunctionMaps(Map<String, Collection<RefSeqGene>> genesByChr, Map<String, Collection<RefSeqGene>> junctionsByChr, String save, String genomeDir) throws IOException{
		FileWriter bedWriter=new FileWriter(save+".bed");
		FileWriter seqWriter=new FileWriter(save+".fa");
		
		//I want to get the concatanated sequence and see how many more I get with just junctions and the full thing
		writeAll(bedWriter, seqWriter, genesByChr, genomeDir, "fragment");
		writeAll(bedWriter, seqWriter, junctionsByChr, genomeDir, "junction");
		
		bedWriter.close();
		seqWriter.close();
	}
	
	public MakeSpliceJunctionMaps(File samFile, String save, String genomeDir, int readLength) throws IOException{
		this.readLength=readLength;
		//FileWriter writer=new FileWriter(save);
		Map<String, Collection<RefSeqGene>> junctionsByChr=new TreeMap<String, Collection<RefSeqGene>>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(samFile)));
		String nextLine;
		int counter=0;
		int total=0;
		while ((nextLine = reader.readLine()) != null) {
			RefSeqGene gene=SAMUtils.SAMFormatFullBED(nextLine);
			if(getPercentMapped(gene, readLength)>this.minPercentMapped && getMinCoverage(gene)>this.minCoverageOnEachSide){
			Collection<RefSeqGene> junc=getJunctions(gene, readLength);
			Collection<RefSeqGene> junctions=new TreeSet();
			if(junctionsByChr.containsKey(gene.getChr())){
				junctions=junctionsByChr.get(gene.getChr());
			}
			junctions.addAll(junc);
			junctionsByChr.put(gene.getChr(), junctions);
			counter++;
			//writer.write(nextLine+"\n");
			}
			total++;
			if(total %100000 ==0){System.err.println(counter+" "+total);}
		}
		reader.close();
		//writer.close();
		write(save, junctionsByChr, genomeDir);
	}
	
	
	private int getMinCoverage(RefSeqGene gene) {
		int min=Integer.MAX_VALUE;
		Collection<Alignments> exons=gene.getSortedAndUniqueExons();
		for(Alignments exon: exons){min=Math.min(min, exon.length());}
		if(min==Integer.MAX_VALUE){min=0;}
		return min;
	}

	private double getPercentMapped(RefSeqGene gene, int readLength) {
		return (double)gene.getTranscriptLength()/readLength;
	}

	//TODO Make this take a SAM file
	public MakeSpliceJunctionMaps(File[] files, String save, String genomeDir) throws IOException {
		Map<String, Collection<RefSeqGene>> junctionsByChr=new TreeMap<String, Collection<RefSeqGene>>();
		
		for(int i=0; i<files.length; i++){
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(files[i])));
		
    		String nextLine;
    		while ((nextLine = reader.readLine()) != null) {
        		try{
        			//TODO Filter by minimum number of bases on each side of the intron
    				//TODO Require a minimum percent of read covered (if 76bp read then at least 70bp needs to be mapped)
        			PSL psl=new PSL(nextLine);              		
	        			if(psl.getPercentMapped()>this.minPercentMapped && psl.getMinCoverage()>this.minCoverageOnEachSide){
	        			RefSeqGene gene=psl.toGene();
	        			Collection<RefSeqGene> junc=getJunctions(gene, readLength);
	        			Collection<RefSeqGene> junctions=new TreeSet();
	        			if(junctionsByChr.containsKey(gene.getChr())){
	        				junctions=junctionsByChr.get(gene.getChr());
	        			}
	        			junctions.addAll(junc);
	        			junctionsByChr.put(gene.getChr(), junctions);
        			}
        		}catch(Exception ex){System.err.println("Skipping: "+nextLine);}	
        	}
        	reader.close();
		}
		
		write(save, junctionsByChr, genomeDir);
		
	}


	


	public MakeSpliceJunctionMaps(File samFile, String save, String genomeDir, int minIntron, int maxIntron) throws IOException {
		AlignmentQueryReader reader = new BAMQueryReader(samFile);//SamQueryReaderFactory.getReader(new ResourceLocator(samFile.getAbsolutePath()), false);
		CloseableIterator<Alignment> iter=reader.iterator();
		
		Map<String, Collection<RefSeqGene>> junctionsByChr=new TreeMap<String, Collection<RefSeqGene>>();
		Collection<Alignments> introns=new TreeSet<Alignments>();
		
		int counter=0;
		while(iter.hasNext()){
			Alignment read=iter.next();
			introns.addAll(getIntrons(read));
			/*Collection<RefSeqGene> junc=getJunctions(read, readLength, minIntron, maxIntron);
			if(junc!=null && !junc.isEmpty()){
				Collection<RefSeqGene> junctions=new TreeSet<RefSeqGene>();
				String chr=junc.iterator().next().getChr();
				if(junctionsByChr.containsKey(chr)){junctions=junctionsByChr.get(chr);}
				junctions.addAll(junc);
				junctionsByChr.put(chr, junctions);
			}*/
			if(counter% 1000000 ==0){System.err.println(counter+" "+introns.size());}
			counter++;
		}
		
		Map<String, Collection<Alignments>> intronsByChr=splitByChr(introns);
		
		FileWriter bedWriter=new FileWriter(save+".bed");
		FileWriter seqWriter=new FileWriter(save+".fa");
		
		//I want to get the concatanated sequence and see how many more I get with just junctions and the full thing
		writeAll(bedWriter, seqWriter, intronsByChr, genomeDir);
		
		bedWriter.close();
		seqWriter.close();
		
		BEDFileParser.writeBED(save, introns);
		
	}
	
	private Map<String, Collection<Alignments>> splitByChr(Collection<Alignments> exons) {
		Map<String, Collection<Alignments>> rtrn=new TreeMap<String, Collection<Alignments>>();
		
		for(Alignments exon: exons){
			Collection<Alignments> list=new TreeSet<Alignments>();
			if(rtrn.containsKey(exon.getChr())){
				list=rtrn.get(exon.getChr());
			}
			list.add(exon);
			rtrn.put(exon.getChr(), list);
		}
		
		return rtrn;
	}

	private Collection<RefSeqGene> getJunctions(RefSeqGene gene, int readLength2, int minIntron, int maxIntron) {
		Collection<RefSeqGene> junctions=new TreeSet<RefSeqGene>();
		Collection<Alignments> introns=gene.getIntronSet();
		
		for(Alignments intron: introns){
			Alignments exon1=new Alignments(intron.getChr(), intron.getStart()-readLength2, intron.getStart());
			Alignments exon2=new Alignments(intron.getChr(), intron.getEnd(), intron.getEnd()+readLength2);
			Collection<Alignments> exons=new TreeSet<Alignments>();
			exons.add(exon1);
			exons.add(exon2);
			RefSeqGene jun=new RefSeqGene(exons);
			junctions.add(jun);
		}
		
		return junctions;
	}
	
	private Collection<RefSeqGene> getJunctions(Alignment read, int readLength2, int minIntron, int maxIntron) {
		Collection<RefSeqGene> junctions=new TreeSet();
		Collection<Alignments> introns=getIntrons(read);
		
		for(Alignments intron: introns){
			Alignments exon1=new Alignments(intron.getChr(), intron.getStart()-readLength2, intron.getStart());
			Alignments exon2=new Alignments(intron.getChr(), intron.getEnd(), intron.getEnd()+readLength2);
			Collection<Alignments> exons=new TreeSet<Alignments>();
			exons.add(exon1);
			exons.add(exon2);
			RefSeqGene jun=new RefSeqGene(exons);
			junctions.add(jun);
		}
		
		return junctions;
	}
	
	private Collection<Alignments> toExons(Alignment record){
		Collection<Alignments> rtrn=new ArrayList<Alignments>();
		AlignmentBlock[] blocks=record.getAlignmentBlocks();
		if(blocks==null){rtrn.add(new Alignments(record.getChromosome(), record.getAlignmentStart(), record.getAlignmentEnd())); return rtrn;}
		for(int i=0; i<blocks.length; i++){
			Alignments exon=new Alignments(record.getChromosome(), blocks[i].getStart(), blocks[i].getEnd());
			rtrn.add(exon);
		}
		
		return rtrn;
	}
	
	private Collection<Alignments> getIntrons(Alignment record){
		Collection<Alignments> exons=this.toExons(record);
		RefSeqGene gene=new RefSeqGene(exons);
		Collection<Alignments> rtrn=gene.getIntronSet();
		return rtrn;
	}
	
	private Collection<RefSeqGene> getJunctions(RefSeqGene gene, int readLength2) {
		return getJunctions(gene, readLength2, -Integer.MAX_VALUE, Integer.MAX_VALUE);
	}

	
	private void write(String save, Collection<RefSeqGene> junctions)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(RefSeqGene junction: junctions){writer.write(junction.toBED()+"\n");}
		
		writer.close();
	}

	private void write(String save, Map<String, Collection<RefSeqGene>> junctionsByChr, String genomeDir) throws IOException {
		FileWriter bedWriter=new FileWriter(save+".bed");
		FileWriter seqWriter=new FileWriter(save+".fa");
		
		int total=0;
		int filtered=0;
		
		for(String chr: junctionsByChr.keySet()){
			System.err.println(chr);
			Collection<RefSeqGene> junctions=junctionsByChr.get(chr);
			String chrSequenceFile=genomeDir+"/"+chr.replaceAll("chr", "")+"/"+chr+".fa";
			FastaSequenceIO fsio = new FastaSequenceIO(chrSequenceFile);
			List<Sequence> seqs = fsio.loadAll();
			Sequence chrSequence = seqs.get(0);
			
			for(RefSeqGene junction: junctions){
				Alignments intron=junction.getIntronSet().iterator().next();
				String orientation=GeneTools.orientationFromSpliceSites(intron, chrSequence);
				
				if(!orientation.equalsIgnoreCase("*") && intron.getSize()>this.minIntronSize && intron.getSize()<this.maxIntronSize){
					filtered++;
					bedWriter.write(junction.toBED()+"\n");
					seqWriter.write(">"+junction.getName()+"\n"+junction.getSequence(chrSequence).getSequenceBases()+"\n");
				}
				total++;
			}
			System.err.println("Total number of reads "+total+" filtered number "+filtered);
		}
			
		bedWriter.close();
		seqWriter.close();
	}
	
	private void writeAll(FileWriter bedWriter, FileWriter seqWriter, Map<String, Collection<RefSeqGene>> junctionsByChr, String genomeDir, String string) throws IOException {
		//FileWriter bedWriter=new FileWriter(save+".bed");
		//FileWriter seqWriter=new FileWriter(save+".fa");
		
		int total=0;
		int filtered=0;
		int i=0;
		for(String chr: junctionsByChr.keySet()){
			try{
			System.err.println(chr);
			Collection<RefSeqGene> junctions=junctionsByChr.get(chr);
			String chrSequenceFile=genomeDir+"/"+chr.replaceAll("chr", "")+"/"+chr+".fa";
			FastaSequenceIO fsio = new FastaSequenceIO(chrSequenceFile);
			List<Sequence> seqs = fsio.loadAll();
			Sequence chrSequence = seqs.get(0);
			
			for(RefSeqGene junction: junctions){
				junction.setName(string+""+i);
				bedWriter.write(junction.toBED()+"\n");
				seqWriter.write(">"+junction.getName()+"\n"+junction.getSequence(chrSequence).getSequenceBases()+"\n");
				i++;
				total++;
			}
			System.err.println("Total number of reads "+total+" filtered number "+filtered);
			}catch(FileNotFoundException ex){System.err.println("Skipping "+chr);}
		}
			
		//bedWriter.close();
		//seqWriter.close();
	}
	
	private void writeAll(FileWriter bedWriter, FileWriter seqWriter, Map<String, Collection<Alignments>> junctionsByChr, String genomeDir) throws IOException {
		//FileWriter bedWriter=new FileWriter(save+".bed");
		//FileWriter seqWriter=new FileWriter(save+".fa");
		
		int total=0;
		int filtered=0;
		int i=0;
		for(String chr: junctionsByChr.keySet()){
			try{
			System.err.println(chr);
			Collection<Alignments> junctions=junctionsByChr.get(chr);
			String chrSequenceFile=genomeDir+"/"+chr.replaceAll("chr", "")+"/"+chr+".fa";
			FastaSequenceIO fsio = new FastaSequenceIO(chrSequenceFile);
			List<Sequence> seqs = fsio.loadAll();
			Sequence chrSequence = seqs.get(0);
			
			for(Alignments intron: junctions){
				RefSeqGene junction=makeJunction(intron, this.readLength);
				junction.setName("junction"+i);
				bedWriter.write(junction.toBED()+"\n");
				seqWriter.write(">"+junction.getName()+"\n"+junction.getSequence(chrSequence).getSequenceBases()+"\n");
				i++;
				total++;
			}
			System.err.println("Total number of reads "+total+" filtered number "+filtered);
			}catch(FileNotFoundException ex){System.err.println("Skipping "+chr);}
		}
			
		//bedWriter.close();
		//seqWriter.close();
	}


	private RefSeqGene makeJunction(Alignments intron, int readLength2) {
		Alignments exon1=new Alignments(intron.getChr(), intron.getStart()-readLength2, intron.getStart());
		Alignments exon2=new Alignments(intron.getChr(), intron.getEnd(), intron.getEnd()+readLength2);
		Collection<Alignments> exons=new TreeSet<Alignments>();
		exons.add(exon1);
		exons.add(exon2);
		RefSeqGene jun=new RefSeqGene(exons);
		return jun;
	}

	public static void main(String[] args)throws IOException{
		if(args.length>4){
			File samFile=new File(args[0]);
			String save=args[1];
			String genomeDir=args[2];
			int minIntron=new Integer(args[3]);
			int maxIntron=new Integer(args[4]);
			new MakeSpliceJunctionMaps(samFile, save, genomeDir, minIntron, maxIntron);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=SAM file \n args[1]=save \n args[2]=genome directory \n args[3]=min intron \n args[4]=max intron";
	
}
