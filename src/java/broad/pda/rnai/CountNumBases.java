package broad.pda.rnai;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.pda.annotation.BEDFileParser;
import broad.pda.chromosome.Chromosome;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;


public class CountNumBases {

	public CountNumBases(Collection<Alignments> regions, Map<String, Collection<RefSeqGene>> genes, String genomeDirectory, File okRepeatFile, String save)throws Exception{
		Map<String, IntervalTree<Alignments>> okRepeats=BEDFileParser.loadAlignmentDataToTree(okRepeatFile);
		Collection<RefSeqGene> rtrn=new TreeSet();
		for(String chr: genes.keySet()){
			System.err.println(chr);
			Collection<RefSeqGene> set=genes.get(chr);
			String sequenceFile=genomeDirectory+"/"+chr.replaceAll("chr", "").trim()+"/"+chr+".agp";
			Chromosome chrom = new Chromosome(sequenceFile);
			chrom.loadSequence();
			for(RefSeqGene align: set){
				String seq=ExtractSequence.getSequenceForGene(align, chrom, true, okRepeats);
				align.setSequence(seq);
				rtrn.add(align);
			}
		}
		
		Collection<RefSeqGene> collapsed=getLargestUnmaskedRegion(regions, rtrn);
		
		write(save, collapsed);
	}
	
	private Collection<RefSeqGene> getLargestUnmaskedRegion(Collection<Alignments> regions, Collection<RefSeqGene> isoforms){
		Collection<RefSeqGene> rtrn=new TreeSet();
		
		Map<String, IntervalTree<RefSeqGene>> isoformTree=makeTree(isoforms);
		
		for(Alignments align: regions){
			Iterator<Node<RefSeqGene>> iter=isoformTree.get(align.getChr()).overlappers(align.getStart(), align.getEnd());
			RefSeqGene best=getBest(iter);
			if(best!=null){rtrn.add(best);}
			else{System.err.println(align);}
		}
		
		return rtrn;
		
	}
	
	private RefSeqGene getBest(Iterator<Node<RefSeqGene>> genes){
		int max=-10;
		RefSeqGene rtrn=null;
		
		while(genes.hasNext()){
			RefSeqGene gene=genes.next().getValue();
			String seq=gene.getSequence();
			int[] counts=getCounts(seq);
			if(counts[1]>=max){max=counts[1]; rtrn=gene;}
		}
		
		return rtrn;
	}
	
	private Map<String, IntervalTree<RefSeqGene>> makeTree(Collection<RefSeqGene> genes){
		Map<String, IntervalTree<RefSeqGene>> rtrn=new TreeMap<String, IntervalTree<RefSeqGene>>();
		
		for(RefSeqGene align: genes){
			IntervalTree<RefSeqGene> tree=new IntervalTree<RefSeqGene>();
			String chr=align.getChr();
			if(rtrn.containsKey(chr)){tree=rtrn.get(chr);}
			tree.put(align.getAlignment().getStart(), align.getAlignment().getEnd(), align);
			rtrn.put(align.getChr(), tree);
		}
		
		return rtrn;
		
	}
	
	private void write(String save, Collection<RefSeqGene> genes)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(RefSeqGene gene: genes){
			String seq=gene.getSequence();
			int[] counts=getCounts(seq);
			writer.write(gene+"\t"+counts[0]+"\t"+counts[1]+"\t"+counts[2]+"\t"+counts[3]+"\n");
		}
		
		writer.close();
	}
	
	private int[] getCounts(String seq){
		int numNs=0;
		int numLower=0;
		int numUpper=0;
		int other=0;
		int total=0;
		
		char[] chars=seq.toCharArray();
		for(int i=0; i<chars.length; i++){
			if(chars[i]=='N' || chars[i]=='n'){numNs++;}
			else if(chars[i]=='a' || chars[i]=='c'|| chars[i]=='g' || chars[i]=='t'){numLower++;}
			else if(chars[i]=='A' || chars[i]=='C' || chars[i]=='G' || chars[i]=='T'){numUpper++;}
			else{other++;}
			total++;
		}
		
		int[] rtrn={total, numUpper, numLower, numNs};
		return rtrn;
	}
	
	static String usage=" args[0]=regions (K4-K36) \n args[1]=isoforms(Full bed) \n args[2]=genomeDirectory \n args[3]=ok repeats \n args[4]=save file";
	
	public static void main(String[] args)throws Exception{
		if(args.length>4){
			Collection<Alignments> regions=BEDFileParser.loadAlignmentData(new File(args[0]));
			Map<String, Collection<RefSeqGene>> isoforms=BEDFileParser.loadDataByChr((new File(args[1])));
			String genomeDirectory=args[2];
			File okRepeats=new File(args[3]);
			String save=args[4];
			new CountNumBases(regions, isoforms, genomeDirectory, okRepeats, save);
		}
		else{System.err.println(usage);}
	}
	
}
