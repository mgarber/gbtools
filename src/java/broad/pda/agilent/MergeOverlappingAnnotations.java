package broad.pda.agilent;

import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.core.util.CollapseByIntersection;
import broad.pda.annotation.BEDFileParser;
import broad.pda.chromosome.Chromosome;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;

public class MergeOverlappingAnnotations {

	public MergeOverlappingAnnotations(File[] files, String save, String genomeDirectory)throws Exception{
		Set<RefSeqGene> annotations=new TreeSet();
		for(int i=0; i<files.length; i++){System.err.println(files[i]);annotations.addAll(BEDFileParser.loadData(files[i]));}
		Map<String, IntervalTree<RefSeqGene>> trees=CollapseByIntersection.makeIntervalTreeForGenes(annotations);
		//Map<RefSeqGene, ArrayList> overlappingAnnotations=overlappingAnnotations(annotations, trees);
		write(save, annotations, trees, genomeDirectory);
	}
	
	
	private void write(String save, Collection<RefSeqGene> annotations, Map<String, IntervalTree<RefSeqGene>> trees, String genomeDirectory)throws Exception{
		FileWriter writer=new FileWriter(save);
		
		String chr="";
		Chromosome chrom=null;
		
		for(RefSeqGene gene: annotations){
			if(!gene.getChr().equalsIgnoreCase(chr)){
				System.err.println(gene.getChr());
				chr=gene.getChr();
				String sequenceFile=genomeDirectory+"/"+chr.replaceAll("chr", "").trim()+"/"+chr+".agp";
				if(new File(sequenceFile).exists()){chrom = new Chromosome(sequenceFile);
					chrom.loadSequence();
				}
				else{chrom=null;}
			}
			//System.err.println(gene.getName());
			Iterator<Node<RefSeqGene>> overlappers=trees.get(gene.getChr()).overlappers(gene.getStart(), gene.getEnd());
			int i=0;
			if(chrom!=null){
			writer.write(">");
			while(overlappers.hasNext()){
				if(i>0){writer.write(",");} 
				RefSeqGene next=overlappers.next().getValue();
				writer.write(next.getName()+","+next.getName()+"_F"+","+next.getName()+"_R"); i++;
			}
			writer.write("\n");
			writer.write(this.getSequenceForGene(gene, chrom));
			writer.write("\n");
			}
		}
		
		writer.close();
	}
	
	private String getSequenceForGene(RefSeqGene gene, Chromosome chrom)throws Exception{
		ArrayList seq=new ArrayList();
		String sequenceString="";
		int counter=0;
		for(int i=0; i<gene.getExons().length; i++){
			Alignments exon=gene.getExons()[i];
			if(exon.getSize()>10){seq.add(getSequenceUnoriented(exon, chrom));}
		}
				
		for(int i=0; i<seq.size(); i++){sequenceString+=(seq.get(i));}
		if(gene.getOrientation().equalsIgnoreCase("-")){sequenceString=reverseComplement(sequenceString);}
		
		return sequenceString;
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
	
	private String getSequenceUnoriented(Alignments align, Chromosome chrom)throws Exception{
		//System.err.println(align);
		SequenceRegion target = new SequenceRegion("chr"+chrom.getSymbol());
		target.setRegionStart(align.getStart());
		target.setRegionEnd(align.getEnd());
		chrom.getRegion(target);
		Sequence seq=target.getSequence();
		
		//System.err.println(align+" "+seq+" "+seq.getSequenceBases());
		
		if(align.getOrientation().equalsIgnoreCase("+")){return seq.getSequenceBases();}
		else if(align.getOrientation().equalsIgnoreCase("-")){return reverseComplement(seq.getSequenceBases());}
		//else{System.err.println("NO STRAND INFO ");}
		return seq.getSequenceBases();
		
	}
	
	private Map<RefSeqGene, ArrayList> overlappingAnnotations(Set<RefSeqGene> annotations, Map<String, IntervalTree> trees){
		Map<RefSeqGene, ArrayList> rtrn=new TreeMap();
		for(RefSeqGene gene: annotations){
			ArrayList<RefSeqGene> list=new ArrayList();
			Iterator<Node<RefSeqGene>> overlappers=trees.get(gene.getChr()).overlappers(gene.getStart(), gene.getEnd());
			while(overlappers.hasNext()){
				RefSeqGene over=overlappers.next().getValue();
				list.add(over);
			}
			rtrn.put(gene, list);
		}
		return rtrn;
	}
	
	public static void main(String[] args)throws Exception{
		if(args.length>2){
			File[] files=new File(args[0]).listFiles();
			String save=args[1];
			String genomeDirectory=args[2];
			new MergeOverlappingAnnotations(files, save, genomeDirectory);
		}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=annotation directory \n args[1]=save file \n args[2]=genome directory";
}
