package broad.pda.rnaseq.expression;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import broad.core.primer3.Primer3SequenceInputTags.SequenceRegionCoordinates;
import broad.core.primer3.PrimerPair;
import broad.core.primer3.qPCRPrimerDesigner;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.chromosome.Chromosome;
import broad.pda.gene.RefSeqGene;

public class DigitalExpressionPrimerDesigner {

	static String usage=" -in <alignment file (Full BED)> \n -seqdir <sequence directory> \n -out <save file or stdout if none specified> \n -target  <target region (distance from end of gene), include as many as desired> \n -repeatMask <Add this flag to repeat mask> \n -numDesigns <options, default=5> \n";
	private String seqDir;
	private List<Integer> targets;
	private int numDesigns;
	private boolean repeatMask;

	public DigitalExpressionPrimerDesigner(String seqDir, List<Integer> targets, int numDesigns, boolean repeatMask) {
		this.seqDir = seqDir;
		this.targets = targets;
		this.numDesigns = numDesigns;
		this.repeatMask = repeatMask;
	}
	

	private Map<RefSeqGene, Collection<PrimerPair>> design(Map<String, Collection<RefSeqGene>> geneMap) throws Exception {
		Map<RefSeqGene, Collection<PrimerPair>> rtrn=new TreeMap<RefSeqGene, Collection<PrimerPair>>();
		for(String chr: geneMap.keySet()){
			//System.err.println("working on " +chr+" ...");
			String sequenceFile=seqDir+"/"+chr.replaceAll("chr", "").trim()+"/"+chr+".agp";
			Chromosome chrom = new Chromosome(sequenceFile);
			chrom.loadSequence();
			Map<RefSeqGene, Collection<PrimerPair>> map=design(chrom, geneMap.get(chr));
			rtrn.putAll(map);
		}
		
		return rtrn;
	}
	

	
	private Map<RefSeqGene, Collection<PrimerPair>> design(Chromosome chrom, Collection<RefSeqGene> genes) throws Exception {
		Map<RefSeqGene, Collection<PrimerPair>> rtrn=new TreeMap<RefSeqGene, Collection<PrimerPair>>();
		
		for(RefSeqGene gene: genes){
			Collection<PrimerPair> primers=design(gene, chrom);
			rtrn.put(gene, primers);
		}
		
		return rtrn;
	}


	private Collection<PrimerPair> design(RefSeqGene gene, Chromosome chrom) throws Exception {
		Collection<PrimerPair> primerList = new ArrayList<PrimerPair> (numDesigns * targets.size());
		for(int i = 0; i < targets.size(); i++) {
			SequenceRegionCoordinates src = null;
			System.err.println("Designing primers for " + gene.getName() + " orientation ["+gene.getOrientation()+"] target: " + targets.get(i));
			if("-".equals(gene.getOrientation())) {
				src = new SequenceRegionCoordinates(targets.get(i), targets.get(i)+ 1) ;
			} else {
				src = new SequenceRegionCoordinates(gene.getTranscriptLength() - targets.get(i), gene.getTranscriptLength() - targets.get(i)+ 1) ;
			}
			try{
			Collection<PrimerPair> targetPrimers = qPCRPrimerDesigner.designPCRPrimers(chrom, gene, repeatMask, numDesigns, src);
			for(PrimerPair pp: targetPrimers) {
				pp.setPrimerPairId(pp.getPrimerPairId() + "_t" + i);
				primerList.add(pp);
			}
			}catch (Exception ex){}
		}
		return primerList;
	}

	private static void write(Map<RefSeqGene, Collection<PrimerPair>> design,BufferedWriter bw) throws IOException {
		bw.write("WellPosition\tRow\tColumn\tName\t Sequence_left\tPriduct Size\tPrimer pair penalty\n");
		int col = 0;
		char row = 65;
		for(RefSeqGene g : design.keySet()) {
			Collection<PrimerPair> gPairs = design.get(g);
			for(PrimerPair pp: gPairs) {
				bw.write(row);
				bw.write(String.valueOf((col+1)));
				bw.write("\t");
				bw.write(row);
				bw.write("\t");
				bw.write(String.valueOf((col+1)));
				bw.write("\t");
				bw.write(pp.getPrimerPairId());
				bw.write("\t");
				bw.write(pp.getLeftPrimer());
				bw.write("\t");
				bw.write(String.valueOf(pp.getProductSize()));
				bw.write("\t");
				bw.write(String.valueOf(pp.getPrimerPairPenalty()));
				bw.newLine();
				col =  (col+1) % 8 ;
				if(col ==0 ) {
					row++;
				}
			}

		}
		bw.write("WellPosition\tRow\tColumn\tName\t Sequence_Right\tPriduct Size\tPrimer pair penalty\n");
		col = 0;
		row = 65;
		for(RefSeqGene g : design.keySet()) {
			Collection<PrimerPair> gPairs = design.get(g);
			for(PrimerPair pp: gPairs) {
				bw.write(row);
				bw.write(String.valueOf((col+1)));
				bw.write("\t");
				bw.write(row);
				bw.write("\t");
				bw.write(String.valueOf((col+1)));
				bw.write("\t");
				bw.write(pp.getPrimerPairId());
				bw.write("\t");
				bw.write(pp.getRightPrimer());
				bw.write("\t");
				bw.write(String.valueOf(pp.getProductSize()));
				bw.write("\t");
				bw.write(String.valueOf(pp.getPrimerPairPenalty()));
				bw.newLine();
				col =  (col+1) % 8 ;
				if(col ==0 ) {
					row++;
				}
			}
		}
		
	}
	
	private static void writeISPCR(Map<RefSeqGene, Collection<PrimerPair>> design, BufferedWriter bw) throws IOException {
		bw.write("ID\tleft_primer\tright_primer\tproduct_size\tleft_primer_position\tright_primer_position\tleft_primer_penalty\tright_primer_penalty\tprimer_pair_penalty\n");
		for(RefSeqGene g : design.keySet()) {
			Collection<PrimerPair> gPairs = design.get(g);
			for(PrimerPair pp: gPairs) {
				bw.write(pp.getPrimerPairId());
				bw.write("\t");
				bw.write(pp.getLeftPrimer());
				bw.write("\t");
				bw.write(pp.getRightPrimer());
				bw.write("\t");
				bw.write(String.valueOf(pp.getProductSize()));
				bw.write("\t");
				bw.write(String.valueOf(pp.getLeftPrimerPosition()));
				bw.write("\t");
				bw.write(String.valueOf(pp.getRightPrimerPosition()));
				bw.write("\t");
				bw.write(String.valueOf(pp.getLeftPrimerPenalty()));
				bw.write("\t");
				bw.write(String.valueOf(pp.getRightPrimerPenalty()));
				bw.write("\t");
				bw.write(String.valueOf(pp.getPrimerPairPenalty()));				
				bw.newLine();
			}
		}		
	}

	public static void main(String[] args)throws Exception{
		ArgumentMap argMap = CLUtil.getParameters(args, usage);
		
		String seqdir = argMap.getMandatory("seqdir");
		List<String> targets = argMap.getAllMandatory("target");
		List<Integer> targetDistanceFromEnd = new ArrayList<Integer>(targets.size());
		for(String t : targets) {
			targetDistanceFromEnd.add(Integer.parseInt(t));
		}
		boolean repeatMask = argMap.containsKey("repeatMask");
		int numDesigns = argMap.containsKey("numDesigns") ? argMap.getInteger("numDesigns") : 5;
		
		DigitalExpressionPrimerDesigner depd = new DigitalExpressionPrimerDesigner(seqdir, targetDistanceFromEnd, numDesigns, repeatMask);
		String in = argMap.getInput();
		Map<String, Collection<RefSeqGene>> genes=BEDFileParser.loadDataByChr(new File(in));
		Map<RefSeqGene, Collection<PrimerPair>> design = depd.design(genes);
		
		
		BufferedWriter bw = argMap.getOutputWriter(); 
		write(design, bw);
		bw.close();
		
		bw = new BufferedWriter(new FileWriter(in+".ispcr"));
		writeISPCR(design,bw);
		bw.close();
	}







}
