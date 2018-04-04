package broad.pda.rnai.designer;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.util.CollapseByIntersection;
import broad.pda.annotation.BEDFileParser;
import broad.pda.chromosome.Chromosome;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.rnai.ExtractSequence;
import broad.pda.rnai.RNAiGeneAnnotation;

public class MatchReportAndStructure {

	File okRepeatFile;
	
	public MatchReportAndStructure(Map<String, Collection<RefSeqGene>> genesByChr, String genomeDirectory, File rnaiReport, String save, File okRepeatFile, Map<String, IntervalTree<RefSeqGene>> tree)throws Exception{
		this.okRepeatFile=okRepeatFile;
		Map<String, RNAiGeneAnnotation> rnai=RNAiFileFormatUtils.parseRNAiReportFileBySequence(rnaiReport);
		Map<String, RefSeqGene> sequences=getSequences(genesByChr, genomeDirectory);
		Map<RNAiGeneAnnotation, RefSeqGene> match=match(rnai, sequences);
		
		//for all non-matching
		for(RNAiGeneAnnotation annotation: match.keySet()){
			RefSeqGene gene=match.get(annotation);
			if(gene==null){
				Iterator<Node<RefSeqGene>> overlappers=tree.get(annotation.getRegion().getChr()).overlappers(annotation.getRegion().getStart(), annotation.getRegion().getEnd());
				while(overlappers.hasNext()){
					RefSeqGene temp=overlappers.next().getValue();
					System.out.println(annotation+"\t"+temp);
				}
			}
		}
		
		write(save, match);
	}
	
	
	private Map<String, RefSeqGene> getSequences(Map<String, Collection<RefSeqGene>> genesByChr, String genomeDirectory) throws Exception {
		Map<String, RefSeqGene> rtrn=new TreeMap();
		Map<String, IntervalTree<Alignments>> okRepeats=BEDFileParser.loadAlignmentDataToTree(okRepeatFile);
		
		
		for(String chr: genesByChr.keySet()){
			System.err.println(chr);
			String sequenceFile=genomeDirectory+"/"+chr.replaceAll("chr", "").trim()+"/"+chr+".agp";
			Chromosome chrom = new Chromosome(sequenceFile);
			chrom.loadSequence();
			for(RefSeqGene align: genesByChr.get(chr)){
				String seq=ExtractSequence.getSequenceForGene(align, chrom, true, okRepeats).toUpperCase();
				rtrn.put(seq, align);
			}
		}
		
		return rtrn;
	}


	private Map<RNAiGeneAnnotation, RefSeqGene> match(Map<String, RNAiGeneAnnotation> rnai, Map<String, RefSeqGene> sequences) {
		Map<RNAiGeneAnnotation, RefSeqGene> rtrn=new HashMap<RNAiGeneAnnotation, RefSeqGene>();
		
		for(String seq: rnai.keySet()){
			RNAiGeneAnnotation annotation=rnai.get(seq);
			RefSeqGene gene=sequences.get(seq);
			rtrn.put(annotation, gene);
		}
		
		return rtrn;
	}


	private void write(String save, Map<RNAiGeneAnnotation, RefSeqGene> match) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(RNAiGeneAnnotation rnai: match.keySet()){
			writer.write(rnai.toString()+"\t"+match.get(rnai)+"\n");
		}
		
		writer.close();
	}


	public static void main(String[] args)throws Exception{
		if(args.length>4){
			Map<String, Collection<RefSeqGene>> genesByChr=BEDFileParser.loadDataByChr(new File(args[0]));
			String genomeDir=args[1];
			File rnaiReport=new File(args[2]);
			String save=args[3];
			File okRepeats=new File(args[4]);
			
			Map<String, IntervalTree<RefSeqGene>> tree=CollapseByIntersection.makeIntervalTreeForGenes(BEDFileParser.loadData(new File(args[0])));
			new MatchReportAndStructure(genesByChr, genomeDir, rnaiReport, save, okRepeats, tree);
		}
		else{
			System.err.println(usage);
		}
	}
	
	static String usage=" args[0]=Full BED \n args[1]=genome Directory \n args[2]=RNAi report \n args[3]=save \n args[4]=okRepeats";
}
