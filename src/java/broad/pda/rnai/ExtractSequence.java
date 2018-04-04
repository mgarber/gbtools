package broad.pda.rnai;

import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.annotation.GenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.pda.annotation.BEDFileParser;
import broad.pda.chromosome.Chromosome;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;


public class ExtractSequence {
	
	static int numNs=0;

	public static Map getSequenceForGene(Collection<RefSeqGene> alignments, String genomeDirectory, boolean repeatMask, Map<String, IntervalTree<Alignments>> okRepeats)throws Exception{
		Map rtrn=new TreeMap();
		Map<String, Set> alignmentsByChr=splitByChr(alignments);
		
		for(String chr: alignmentsByChr.keySet()){
			Set<RefSeqGene> set=alignmentsByChr.get(chr);
			String sequenceFile=genomeDirectory+"/"+chr.replaceAll("chr", "").trim()+"/"+chr+".agp";
			Chromosome chrom = new Chromosome(sequenceFile);
			chrom.loadSequence();
			for(RefSeqGene align: set){
				//System.err.println(align.toBED());
				rtrn.put(align, getSequenceForGene(align, chrom, repeatMask, okRepeats));
			}
		}
		return rtrn;
	}
	
	public static Map splitByChr(Collection<RefSeqGene> genes){
		Map<String, Set> rtrn=new TreeMap();
		
		for(RefSeqGene gene: genes){
			Set set=new TreeSet();
			String chr=gene.getChr();
			if(rtrn.containsKey(chr)){
				set=rtrn.get(chr);
			}
			set.add(gene);
			rtrn.put(chr, set);
		}
		
		return rtrn;
	}
	
	public static String getSequenceForGene(RefSeqGene gene, Chromosome chrom, boolean repeatMask, Map<String, IntervalTree<Alignments>> okRepeats, boolean stranded)throws Exception{
		ArrayList seq=new ArrayList();
		String sequenceString="";
		int counter=0;
		for(int i=0; i<gene.getExons().length; i++){
			Alignments exon=gene.getExons()[i];
			seq.add(getSequenceUnoriented(exon, chrom, repeatMask, okRepeats));
		}
		String Ns="";
		for(int i=0; i<numNs; i++){Ns+="N";}
		
		for(int i=0; i<seq.size(); i++){sequenceString+=(Ns+seq.get(i));}
		if(gene.getOrientation().equalsIgnoreCase("-") && stranded){sequenceString=Sequence.reverseSequence(sequenceString);}
		
		return sequenceString;
	}
	
	
	public static String getSequenceForGene(RefSeqGene gene, Chromosome chrom, boolean repeatMask, Map<String, IntervalTree<Alignments>> okRepeats)throws Exception{
		return getSequenceForGene(gene, chrom, repeatMask, okRepeats, true);
	}
	
	public static String getSequenceUnoriented(Alignments align, Chromosome chrom, boolean repeatMask, Map<String, IntervalTree<Alignments>> okRepeats)throws Exception{
		//System.err.println(align);
		SequenceRegion target = new SequenceRegion("chr"+chrom.getSymbol());
		target.setRegionStart(align.getStart());
		target.setRegionEnd(align.getEnd());
		target.setChromosome(align.getChr());
		chrom.getRegion(target, repeatMask, okRepeats);
		Sequence seq=target.getSequence();
		
		//System.err.println(align+" "+seq+" "+seq.getSequenceBases());
		
		if(align.getOrientation().equalsIgnoreCase("+")){return seq.getSequenceBases();}
		else if(align.getOrientation().equalsIgnoreCase("-")){return reverseComplement(seq.getSequenceBases());}
		//else{System.err.println("NO STRAND INFO ");}
		return seq.getSequenceBases();
		
	}
	
	public static String getSequenceUnoriented(Alignments align, Chromosome chrom, boolean repeatMask)throws Exception{
		//System.err.println(align);
		SequenceRegion target = new SequenceRegion("chr"+chrom.getSymbol());
		target.setRegionStart(align.getStart());
		target.setRegionEnd(align.getEnd());
		target.setChromosome(align.getChr());
		chrom.getRegion(target, repeatMask);
		Sequence seq=target.getSequence();
		
		//System.err.println(align+" "+seq+" "+seq.getSequenceBases());
		
		if(align.getOrientation().equalsIgnoreCase("+")){return seq.getSequenceBases();}
		else if(align.getOrientation().equalsIgnoreCase("-")){return reverseComplement(seq.getSequenceBases());}
		//else{System.err.println("NO STRAND INFO ");}
		return seq.getSequenceBases();
		
	}
	
	public static String getSequenceForGene(RefSeqGene gene, String genomeDirectory, boolean repeatMask)throws Exception{
		/*String[] seq=new String[gene.getExons().length];
		String sequenceString="";
		for(int i=0; i<seq.length; i++){
			seq[i]=getSequenceUnoriented(gene.getExons()[i], genomeDirectory, repeatMask);
		}
		String Ns="";
		for(int i=0; i<numNs; i++){Ns+="N";}
		
		for(int i=0; i<seq.length; i++){sequenceString+=(Ns+seq[i]);}
		if(gene.getOrientation().equalsIgnoreCase("-")){sequenceString=reverseComplement(sequenceString);}
		
		Sequence rtrn=new Sequence(gene.getName());
		rtrn.setSequenceBases(sequenceString);
		return rtrn;*/
		String sequenceFile=genomeDirectory+"/"+gene.getChr().replaceAll("chr", "").trim()+"/"+gene.getChr()+".agp";
		Chromosome chrom = new Chromosome(sequenceFile);
		chrom.loadSequence();
		return getSequenceForGene(gene, chrom, repeatMask, new TreeMap());
	}
	
	
	/**
	 * Suggestion:  Use classes in net.sf.picard.reference instead, which are WAY faster because they use an indexed FASTA file
	 * @param annot
	 * @param genomeDirectory
	 * @param repeatMask
	 * @return
	 * @throws Exception
	 */
	public static String getGenomeSequence(GenomicAnnotation annot, String genomeDirectory, boolean repeatMask) throws Exception {
		String sequenceFile=genomeDirectory+"/"+annot.getChromosome().replaceAll("chr", "").trim()+"/"+annot.getChromosome()+".agp";
		Chromosome chrom = new Chromosome(sequenceFile);
		chrom.loadSequence();
		SequenceRegion target = new SequenceRegion("chr" + chrom.getSymbol(), annot);
		chrom.getRegion(target, repeatMask);
		return target.getSequenceBases();
	}
	
	
	private static String getSequenceUnoriented(Alignments align, String genomeDirectory, boolean repeatMask, Map<String, IntervalTree<Alignments>> okRepeats)throws Exception{
		String sequenceFile=genomeDirectory+"/"+align.getChr().replaceAll("chr", "").trim()+"/"+align.getChr()+".agp";
		Chromosome chrom = new Chromosome(sequenceFile);
		chrom.loadSequence();
		SequenceRegion target = new SequenceRegion("chr"+chrom.getSymbol());
		target.setChromosome(align.getChr());
		target.setRegionStart(align.getStart());
		target.setRegionEnd(align.getEnd());
		chrom.getRegion(target, repeatMask, okRepeats);
		Sequence seq=target.getSequence();
		
		if(align.getOrientation().equalsIgnoreCase("+")){return seq.getSequenceBases();}
		else if(align.getOrientation().equalsIgnoreCase("-")){return reverseComplement(seq.getSequenceBases());}
		//else{System.err.println("NO STRAND INFO ");}
		return seq.getSequenceBases();
		
	}
	
	private static String getSequenceUnoriented(Alignments align, String genomeDirectory, boolean repeatMask)throws Exception{
		return getSequenceUnoriented(align, genomeDirectory, repeatMask, null);
	}
	
	public static String reverseComplement(String seq){
		char[] seqChar=seq.toCharArray();
		char[] reverse=new char[seqChar.length];
		
		int j=0;
		for(int i=seqChar.length-1; i>=0; i--){
			reverse[j]=seqChar[i];
			j++;
		}
		
		String rtrn="";
		for(int i=0; i<reverse.length; i++){
			if(reverse[i]=='A' || reverse[i]=='a'){rtrn+='T';}
			if(reverse[i]=='T' || reverse[i]=='t'){rtrn+='A';}
			if(reverse[i]=='C' || reverse[i]=='c'){rtrn+='G';}
			if(reverse[i]=='G' || reverse[i]=='g'){rtrn+='C';}
			if(reverse[i]=='N' || reverse[i]=='n'){rtrn+='N';}
		}
		
		return rtrn;
	}
	
	public static void writeFasta(String save, Map<String, Collection<RefSeqGene>> genesByChr, String genomeDir, boolean repeatMask)throws Exception{
		FileWriter writer=new FileWriter(save);
		
		int i=0;
		for(String chr: genesByChr.keySet()){
			try{
			System.err.println(chr);
			String sequenceFile=genomeDir+"/"+chr.replaceAll("chr", "").trim()+"/"+chr+".agp";
			Chromosome chrom = new Chromosome(sequenceFile);
			chrom.loadSequence();
			Collection<RefSeqGene> genes=genesByChr.get(chr);
			for(RefSeqGene gene: genes){
				String seq=getSequenceForGene(gene, chrom, repeatMask, new TreeMap());
				writer.write(">"+gene.getName()+"\n");
				writer.write(seq+"\n");
				i++;
				if(i% 1000 ==0){System.err.println(i);}
			}
			}catch(Exception ex){System.err.println("Skipping "+chr);}
		}
				
		writer.close();
	}

	
	
	public static void main(String[] args)throws Exception{
		if(args.length>3){
			Map<String, Collection<RefSeqGene>> genes=BEDFileParser.loadDataByChr(new File(args[0]));
			String genomeDir=args[1];
			boolean repeatMask=new Boolean(args[2]);
			String save=args[3];
			writeFasta(save, genes, genomeDir, repeatMask);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=genes \n args[1]=genomeDir \n args[2]=repeatMask \n args[3]=save";
	
	
}
