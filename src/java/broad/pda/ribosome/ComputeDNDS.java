package broad.pda.ribosome;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.broad.tribble.FeatureSource;
import org.broad.tribble.source.BasicFeatureSource;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;
import broad.pda.ribosome.misc.FindORFs;

public class ComputeDNDS {
	
	static String geneticCodeFile="/seq/mguttman/Annotations/GeneticCode.txt";
	Map<String, Collection<RefSeqGene>> genesByChr;
	String genomeDir;
	File vcfFile;
	
	public ComputeDNDS(File bedFile, File vcfFile, String genomeDir) throws IOException{
		this.genesByChr=BEDFileParser.loadDataByChr(bedFile);
		this.genomeDir=genomeDir;
		this.vcfFile=vcfFile;
	}
	
	private void scoreCDS(String save, int minNumSNPs) throws IOException{
		FileWriter writer=new FileWriter(save);
		FeatureSource<VariantContext> source = BasicFeatureSource.getFeatureSource(vcfFile.getAbsolutePath(), new VCFCodec(), true);
		
		for(String chr: genesByChr.keySet()){
			System.err.println(chr);
			String chrSeqFile=genomeDir+"/"+chr.replaceAll("chr", "")+"/"+chr+".fa";
			if(new File(chrSeqFile).exists()){
				FastaSequenceIO fsio = new FastaSequenceIO(chrSeqFile);
				Sequence chrSeq=fsio.loadAll().get(0);
				Collection<RefSeqGene> genes=genesByChr.get(chr);
				for(RefSeqGene gene: genes){
					RefSeqGene orf=gene.getCDS();
					Iterator<VariantContext> iter = source.query(orf.getChrNum(), orf.getStart(), orf.getEnd());// chr17:35,647,006-35,647,114
					double[] dnds=dNdS(orf, iter, chrSeq,0);
					double total=dnds[0]+dnds[1];
					if(total>=minNumSNPs){writer.write(gene.getName()+"\t"+gene.getCDS().getSize()+"\tCDS\t"+dnds[0]+"\t"+dnds[1]+"\n");}
				}
			}
		}
		source.close();
		writer.close();
	}
	
	
	private void scoreIntrons(String save, int minNumSNPs) throws IOException{
		FileWriter writer=new FileWriter(save);
		FeatureSource<VariantContext> source = BasicFeatureSource.getFeatureSource(vcfFile.getAbsolutePath(), new VCFCodec(), true);
		
		for(String chr: genesByChr.keySet()){
			System.err.println(chr);
			String chrSeqFile=genomeDir+"/"+chr.replaceAll("chr", "")+"/"+chr+".fa";
			if(new File(chrSeqFile).exists()){
				FastaSequenceIO fsio = new FastaSequenceIO(chrSeqFile);
				Sequence chrSeq=fsio.loadAll().get(0);
				Collection<RefSeqGene> genes=genesByChr.get(chr);
				for(RefSeqGene gene: genes){
					if(gene.getIntrons()!=null){
						Collection<RefSeqGene> intronORFs=FindORFs.findAllORFs(gene.getIntrons(), chrSeq, false);
						int i=0;
						for(RefSeqGene intronORF: intronORFs){
							Iterator<VariantContext> iter = source.query(intronORF.getChrNum(), intronORF.getStart(), intronORF.getEnd());// chr17:35,647,006-35,647,114
							double[] dnds=dNdS(intronORF, iter, chrSeq, 0);
							double total=dnds[0]+dnds[1];
							if(total>=minNumSNPs){writer.write(gene.getName()+"\t"+intronORF.getSize()+"\tIntron_"+i+"\t"+dnds[0]+"\t"+dnds[1]+"\n");}
							i++;
						}
					}
				}
			}
		}
		source.close();
		writer.close();
	}
	
	private void scoreAllORFs(String save, int minNumSNPs) throws IOException{
		FileWriter writer=new FileWriter(save);
		//FileWriter writer2=new FileWriter(save+".bed");
		FeatureSource<VariantContext> source = BasicFeatureSource.getFeatureSource(vcfFile.getAbsolutePath(), new VCFCodec(), true);
		
		for(String chr: genesByChr.keySet()){
			System.err.println(chr);
			String chrSeqFile=genomeDir+"/"+chr.replaceAll("chr", "")+"/"+chr+".fa";
			if(new File(chrSeqFile).exists()){
				FastaSequenceIO fsio = new FastaSequenceIO(chrSeqFile);
				Sequence chrSeq=fsio.loadAll().get(0);
				Collection<RefSeqGene> genes=genesByChr.get(chr);
				for(RefSeqGene gene: genes){
					Collection<RefSeqGene> intronORFs=FindORFs.findAllORFs(gene, chrSeq, false);
						int i=0;
						for(RefSeqGene intronORF: intronORFs){
							Iterator<VariantContext> iter = source.query(intronORF.getChrNum(), intronORF.getStart(), intronORF.getEnd());// chr17:35,647,006-35,647,114
							double[] dnds=dNdS(intronORF, iter, chrSeq, 0);
							double ratio=(dnds[0]+1)/(dnds[1]+1);
							double total=dnds[0]+dnds[1];
							if(total>=minNumSNPs){
								writer.write(gene.getName()+"\t"+intronORF.getSize()+"\tORF_"+i+"\t"+dnds[0]+"\t"+dnds[1]+"\n");
								intronORF.setBedScore(ratio);
								//writer2.write(intronORF+"\n");
							}
							i++;
					}
				}
			}
		}
		source.close();
		writer.close();
		//writer2.close();
	}
	
	
	
	private void scoreDefinedORF(String save, int minNumSNPs) throws IOException{
		FileWriter writer=new FileWriter(save);
		FileWriter writer2=new FileWriter(save+".bed");
		FeatureSource<VariantContext> source = BasicFeatureSource.getFeatureSource(vcfFile.getAbsolutePath(), new VCFCodec(), true);
		
		for(String chr: genesByChr.keySet()){
			System.err.println(chr);
			String chrSeqFile=genomeDir+"/"+chr.replaceAll("chr", "")+"/"+chr+".fa";
			if(new File(chrSeqFile).exists()){
				FastaSequenceIO fsio = new FastaSequenceIO(chrSeqFile);
				Sequence chrSeq=fsio.loadAll().get(0);
				Collection<RefSeqGene> genes=genesByChr.get(chr);
				for(RefSeqGene gene: genes){
					RefSeqGene orf=gene.getCDS();
						Iterator<VariantContext> iter = source.query(orf.getChrNum(), orf.getStart(), orf.getEnd());// chr17:35,647,006-35,647,114
						double[] dnds=dNdS(orf, iter, chrSeq, 0);
						double ratio=(dnds[0]+1)/(dnds[1]+1);
						double total=dnds[0]+dnds[1];
						if(total>=minNumSNPs){
							writer.write(gene.getName()+"\t"+orf.getSize()+"\t"+dnds[0]+"\t"+dnds[1]+"\t"+ratio+"\n");
						}
						orf.setName("dN/dS="+ratio);
						writer2.write(orf+"\n");
							
					
				}
			}
		}
		source.close();
		writer.close();
		writer2.close();
	}
	
	private void scanWindows(String save, int minNumSNPs, int windowSize) throws IOException{
		FileWriter writer=new FileWriter(save);
		FeatureSource<VariantContext> source = BasicFeatureSource.getFeatureSource(vcfFile.getAbsolutePath(), new VCFCodec(), true);
		
		for(String chr: genesByChr.keySet()){
			System.err.println(chr);
			String chrSeqFile=genomeDir+"/"+chr.replaceAll("chr", "")+"/"+chr+".fa";
			if(new File(chrSeqFile).exists()){
				FastaSequenceIO fsio = new FastaSequenceIO(chrSeqFile);
				Sequence chrSeq=fsio.loadAll().get(0);
				Collection<RefSeqGene> genes=genesByChr.get(chr);
				for(RefSeqGene gene: genes){
					GeneComponent gc=new GeneComponent(gene);
					Collection<RefSeqGene> windowsFrame0=gc.getGene().getWindows(windowSize, 3, 0);
					double minFrame0=Double.MAX_VALUE;
					RefSeqGene minFrame0Window=null;
					for(RefSeqGene window: windowsFrame0){
						Iterator<VariantContext> iter = source.query(window.getChrNum(), window.getStart(), window.getEnd());// chr17:35,647,006-35,647,114
						double[] dnds=dNdS(window, iter, chrSeq, 0);//NEED to iterate frame
						double total=dnds[0]+dnds[1];
						//if(total>minNumSNPs){
							double ratio=(dnds[0]+1)/(dnds[1]+1);
							if(ratio<minFrame0){
								minFrame0=ratio;
								minFrame0Window=window;
							}
						//}
					}
					writer.write(minFrame0Window+"\t"+minFrame0+"\n");
				}
			}
		}
		source.close();
		writer.close();
	}
	
	
	private static double[] dNdS(RefSeqGene orf, String vcfFile, Sequence chrSequence, int frameOffset) throws IOException{
		double countSyn=0;
		double countNonSyn=0;
		FeatureSource<VariantContext> source = BasicFeatureSource.getFeatureSource(vcfFile, new VCFCodec(), true);
		Iterator<VariantContext> iter = source.query(orf.getChrNum(), orf.getStart(), orf.getEnd());// chr17:35,647,006-35,647,114
	    while (iter.hasNext()) {
	    	VariantContext variant = iter.next();
	    	variant.getStart(); //determine the frame of the variant
	    	boolean isCDS=isInCDS(orf, variant);//test if its in the exon
	    	if(isCDS){
	    		boolean isSyn=isSyn(orf, variant, chrSequence, frameOffset);//if so, determine the frame
	    		if(isSyn){countSyn++;}
	    		else{countNonSyn++;}
	    	}
	    }
	    	
       source.close();
	
       double[] scores={countNonSyn, countSyn};
       return scores;
	}
	
	private static double[] dNdS(RefSeqGene orf, Iterator<VariantContext> iter, Sequence chrSequence, int frame) throws IOException{
		double countSyn=0;
		double countNonSyn=0;
		while (iter.hasNext()) {
	    	VariantContext variant = iter.next();
	    	variant.getStart(); //determine the frame of the variant
	    	boolean isCDS=isInCDS(orf, variant);//test if its in the exon
	    	if(isCDS){
	    		boolean isSyn=isSyn(orf, variant, chrSequence, frame);//if so, determine the frame
	    		if(isSyn){countSyn++;}
	    		else{countNonSyn++;}
	    	}
	    }
	    	
	
       double[] scores={countNonSyn, countSyn};
       return scores;
	}
	
	private static boolean isSyn(RefSeqGene orf, VariantContext variant, Sequence chrSequence, int frame) throws IOException {
		Collection<String> codons=getCodon(orf, variant, chrSequence, frame);
		
		Collection<String> translation=translate(codons);
		
		if(translation.isEmpty()){throw new IllegalArgumentException("Empty translation");}
		else if(translation.size()==1){return true;}
		return false;
	}
	
	private static Collection<String> translate(Collection<String> codons) throws IOException {
		Map<String, String> map=parseGeneticCode(geneticCodeFile);
		
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String codon: codons){
			String trans=map.get(codon);
			rtrn.add(trans);
		}
		
		return rtrn;
	}

	private static Map<String, String> parseGeneticCode(String geneticCodeFile2) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		Collection<String> list=BEDFileParser.loadList(geneticCodeFile2);
		
		for(String line: list){
			String[] tokens=line.split("\t");
			rtrn.put(tokens[0].trim(), tokens[1].trim());
		}
		
		return rtrn;
	}

	private static Collection<String> getCodon(RefSeqGene orf, VariantContext variant, Sequence chrSequence, int frame) throws IOException{
		//determine the frame
		int relativePosition=orf.absoluteToRelativePosition(variant.getStart());
		int codonPosition=(relativePosition+frame)%3;//make sure this works
		
		RefSeqGene codon=orf.trim(relativePosition, relativePosition+3);
		if(codonPosition==1){
			codon=orf.trim(relativePosition-1, relativePosition+2);
		}
		if(codonPosition==2){
			codon=orf.trim(relativePosition-2, relativePosition+1);
		}
		
		Sequence codonSeq=codon.getSequence(chrSequence);
		
		Collection<String> alternativeCodons=getAlternativeCodons(codonSeq, codonPosition, variant, orf);
		
		Map<String, String> code=parseGeneticCode(geneticCodeFile);
		
		//System.out.println(codon.getChr()+"\t"+codon.getStart()+"\t"+codon.getEnd()+"\tREF_"+codonSeq.getSequenceBases()+"_"+code.get(codonSeq.getSequenceBases()));
				
		for(String alternative: alternativeCodons){
			System.out.println(codon.getChr()+"\t"+codon.getStart()+"\t"+codon.getEnd()+"\t"+alternative+"_"+code.get(alternative));
		}
		
		return alternativeCodons;
	}

	private static Collection<String> getAlternativeCodons(Sequence codonSeq, int codonPosition, VariantContext variant, RefSeqGene orf) {
		Collection<String> rtrn=new TreeSet<String>(); 
		Map<String, Genotype> genotypes= variant.getGenotypes();
		 for(String sample: genotypes.keySet()){
			 Genotype genotype=genotypes.get(sample);
	        	List<Allele> alleles=genotype.getAlleles();
	        	for(Allele allele: alleles){
	        		if(allele.isNonReference() && allele.isCalled()){
	        			String seq=getVariantSeq(codonSeq.getSequenceBases(), codonPosition, allele.getBaseString());
	        			if(orf.getOrientation().equalsIgnoreCase("-")){seq=Sequence.reverseSequence(seq);}
	        			rtrn.add(seq.toUpperCase());
	        		}
	        	}
		 }
		 
		 String codon=codonSeq.getSequenceBases();
		 if(orf.getOrientation().equalsIgnoreCase("-")){codon=Sequence.reverseSequence(codon);}
		 rtrn.add(codon.toUpperCase());
		 return rtrn;
	}

	private static String getVariantSeq(String sequenceBases, int codonPosition, String baseString) {
		char[] seq=sequenceBases.toCharArray();
		
		String variant="";
		
		for(int i=0; i<seq.length; i++){
			if(i==codonPosition){
				variant+=baseString;
			}
			else{variant+=seq[i];}
		}
		
		return variant;
	}

	private static boolean isInCDS(RefSeqGene orf, VariantContext variant) {
		return orf.isInExon(variant.getStart());
	}

	public static void main(String[] args)throws IOException{
		if(args.length>5){
			File bedFile=new File(args[0]);
			File vcfFile=new File(args[1]);
			String genomedir=args[2];
			String save=args[3];
			int minNumSNPs=new Integer(args[4]);
			boolean allORFs=new Boolean(args[5]);
			ComputeDNDS dnds=new ComputeDNDS(bedFile, vcfFile, genomedir);
			if(allORFs){
				dnds.scoreAllORFs(save, minNumSNPs);
			}
			else{dnds.scoreDefinedORF(save, minNumSNPs);}
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=bed file \n args[1]=vcf file \n args[2]=genome directory \n args[3]=save \n args[4]=min Num SNPs \n args[5]=all ORFs?";
	
}
