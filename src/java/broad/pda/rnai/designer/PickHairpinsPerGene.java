package broad.pda.rnai.designer;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.motif.SearchException;
import broad.core.motif.SequenceMotif;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.pda.chromosome.Chromosome;
import broad.pda.datastructures.Alignments;
import broad.pda.rnai.RNAiGeneAnnotation;


//Take all hairpins for all isoforms and rank how often it was picked and where it lands in each isoform
public class PickHairpinsPerGene {

	int minNumber=5;
	int maxNumber=10;
	String genomeDir="/seq/genome/mouse/mouse_Mm9/";
	
	//Merge hairpins across all isoforms
	//Prioritize selection of ones that hit multiple isoforms
	public PickHairpinsPerGene(File rnaiAnnotationFile, File[] hairpinFiles, String save)throws Exception{
		Collection<RNAiGeneAnnotation> rnaiAnnotations=RNAiFileFormatUtils.parseRNAiReportFile(rnaiAnnotationFile);
		Map<String, Collection<HairpinKmer>> hairpins=RNAiFileFormatUtils.parseHairpinDesignFiles(hairpinFiles);
		
		Map<String, Collection<RNAiGeneAnnotation>> geneVariants=getVariantsPerGene(rnaiAnnotations);
		Map<String, Collection<HairpinKmer>> hairpinsPerGene=collapseVariantsToGenes(rnaiAnnotations, hairpins, geneVariants);
		
		Map<String, Collection<HairpinKmer>> geneHairpins=getBestHairpins(rnaiAnnotations, hairpinsPerGene, geneVariants);
		
		/**For debugging**/
		//Map each kmer to the exact location on DNA that it maps
		//Print the genome coordinates for each kmer
		//Map<String, Collection<Alignments>> mapped=mapHairpinCoordinates(geneVariants, geneHairpins, genomeDir);
		//write(save+".mapped.bed", mapped);
		/*****************/
		
		write(save, geneHairpins, geneVariants);
	}

	

	private Map<String, Collection<Alignments>> mapHairpinCoordinates(Map<String, Collection<RNAiGeneAnnotation>> geneVariants, Map<String, Collection<HairpinKmer>> geneHairpins, String genomeDirectory) throws Exception {
		Map<String, Alignments> geneLocations=getGeneCoordinates(geneVariants); //gene by region
		Map<String, Collection<Alignments>> hairpinLocations=new TreeMap();
		
		Chromosome chrom = null;
		String chr="";
		
		for(String gene: geneHairpins.keySet()){
			Alignments geneRegion=geneLocations.get(gene);
			System.err.println(gene+" "+geneRegion);
			Collection<HairpinKmer> hairpins=geneHairpins.get(gene);
			if(!geneRegion.getChr().equalsIgnoreCase(chr)){
				chr=geneRegion.getChr();
				System.err.println("updating chromosome to "+chr);
				String sequenceFile=genomeDirectory+"/"+chr.replaceAll("chr", "").trim()+"/"+chr+".agp";
				chrom=new Chromosome(sequenceFile);
				chrom.loadSequence();
			}
			SequenceRegion target = new SequenceRegion("chr"+chrom.getSymbol());
			target.setRegionStart(geneRegion.getStart());
			target.setRegionEnd(geneRegion.getEnd());
			target.setChromosome(geneRegion.getChr());
			chrom.getRegion(target, false);
			Sequence seq=target.getSequence();
			
			for(HairpinKmer kmer: hairpins){
				String kmerSeq=kmer.getKmerSequence();
				if(geneRegion.getOrientation().equalsIgnoreCase("-")){kmerSeq=ComputeMIRScore.reverseComplement(kmerSeq); System.err.println("Reverse Comp");}
				SequenceMotif motif=new SequenceMotif(kmerSeq, 1);
				List<SequenceRegion> regions=motif.match(seq);
				hairpinLocations.put(kmer.getKmerSequence(), add(regions, geneRegion.getStart()));
			}
		}
		return hairpinLocations;
	}

	private Collection<Alignments> add(List<SequenceRegion> regions, int start) {
		Collection<Alignments> rtrn=new TreeSet();
		
		for(SequenceRegion region: regions){
			Alignments align=new Alignments(region.getChromosome(), region.getAbsoluteStart(start), region.getAbsoluteEnd(start));
			rtrn.add(align);
		}
		
		return rtrn;
	}



	private Map<String, Alignments> getGeneCoordinates(Map<String, Collection<RNAiGeneAnnotation>> geneVariants){
		Map<String, Alignments> geneLocations=new TreeMap();
		for(String gene: geneVariants.keySet()){
			String chr="";
			int start=Integer.MAX_VALUE;
			int end=-Integer.MAX_VALUE;
			String strand="+";
			Collection<RNAiGeneAnnotation> transcripts=geneVariants.get(gene);
			for(RNAiGeneAnnotation transcript: transcripts){
				Alignments transcriptRegion=transcript.getRegion();
				chr=transcriptRegion.getChr();
				start=Math.min(start, transcriptRegion.getStart());
				end=Math.max(end, transcriptRegion.getEnd());
				strand=transcriptRegion.getOrientation();
			}
			Alignments align=new Alignments(chr, start, end);
			align.setOrientation(strand);
			geneLocations.put(gene, align);
		}
		return geneLocations;
	}

	private void write(String string, Map<String, Collection<Alignments>> mapped) throws IOException{
		FileWriter writer=new FileWriter(string);
		
		for(String str: mapped.keySet()){
			Collection<Alignments> regions=mapped.get(str);
			for(Alignments region: regions){
				writer.write(region+"\t"+str+"\n");
			}
		}
		
		writer.close();
	}



	//For each gene find the best hairpins
	//best is defined as targetting as many isoforms as possible
	private Map<String, Collection<HairpinKmer>> getBestHairpins(Collection<RNAiGeneAnnotation> rnaiAnnotations, Map<String, Collection<HairpinKmer>> hairpinsPerGene, Map<String, Collection<RNAiGeneAnnotation>> geneVariants) throws SearchException {
		Map<String, Collection<HairpinKmer>> rtrn=new TreeMap();

		for(String gene: hairpinsPerGene.keySet()){
			Collection<HairpinKmer> hairpins=hairpinsPerGene.get(gene);
			//try to get min number covering all
			hairpins=computeVariantCoverage(hairpins, geneVariants.get(gene));
			Collection<HairpinKmer> topKmers=pickTopHairpins(hairpins, geneVariants.get(gene));
			rtrn.put(gene, topKmers);
		}
		
		return rtrn;
	}

	
	private Collection<HairpinKmer> computeVariantCoverage(Collection<HairpinKmer> hairpins, Collection<RNAiGeneAnnotation> transcriptVariants) throws SearchException {
		Collection<HairpinKmer> rtrn=new TreeSet();
		//for each hairpin score how many variants would be covered
		for(HairpinKmer hairpin: hairpins){
			double sum=0;
			double count=0;
			for(RNAiGeneAnnotation transcript: transcriptVariants){
				boolean matches=match(hairpin.getKmerSequence(), transcript.getSequence());
				if(matches){sum++; hairpin.addTranscript(transcript);}
				//System.out.println(hairpin.getKmerSequence()+"\t"+transcript.getTranscriptName()+"\t"+matches);
				count++;
			}
			hairpin.setPercentVariantsCovered(sum/count);
			rtrn.add(hairpin);
		}
		return rtrn;
	}

	private boolean match(String kmerSequence, Sequence sequence) throws SearchException {
		SequenceMotif motif=new SequenceMotif(kmerSequence,1);
		List<SequenceRegion> matches=motif.match(sequence);
		return matches.size()>0;
	}

	private Collection<HairpinKmer> pickTopHairpins(Collection<HairpinKmer> list, Collection<RNAiGeneAnnotation> transcripts) throws SearchException {
		//go through a sorted list and pick until 1) the list runs out or 2) we get to min number for all variants
		
		Collection<HairpinKmer> hairpins=new TreeSet(list);
		
		Map<RNAiGeneAnnotation, IntervalTree<HairpinKmer>> best=new HashMap();
		for(RNAiGeneAnnotation gene: transcripts){
			IntervalTree<HairpinKmer> tree=new IntervalTree();
			best.put(gene, tree);
		}
		
		//The hairpin list is sorted and so we'll go through it in order
		for(HairpinKmer hp: hairpins){
			//we got the "best" hairpin lets see how many are covered per isoform
			best=getIsoformCoverage(hp, best);
			boolean done=isDone(best);
			if(done){break;}
			//if it doesnt add any new info to the list (bec of overlaps throw it away) or because the transcript is fully covered
		}
		
		Collection<HairpinKmer> rtrn=convertToCollection(best);
		
		return rtrn;
	}
	
	
	private Collection<HairpinKmer> convertToCollection(Map<RNAiGeneAnnotation, IntervalTree<HairpinKmer>> best) {
		Collection<HairpinKmer> rtrn=new TreeSet();
		
		for(RNAiGeneAnnotation gene: best.keySet()){
			IntervalTree<HairpinKmer> tree=best.get(gene);
			Iterator<Node<HairpinKmer>> iter=tree.iterator();
			while(iter.hasNext()){
				HairpinKmer kmer=iter.next().getValue();
				rtrn.add(kmer);
			}
		}
		
		return rtrn;
	}

	private boolean isDone(Map<RNAiGeneAnnotation, IntervalTree<HairpinKmer>> best) {
		boolean rtrn=true;
		
		int totalCounts=countAllKmers(best);
		if(totalCounts>=this.maxNumber){return true;}
		
		for(RNAiGeneAnnotation transcript: best.keySet()){
			IntervalTree<HairpinKmer> kmers=best.get(transcript);
			int num=count(kmers);
			if(num<this.minNumber){return false;}
		}
		
		return rtrn;
	}

	private int countAllKmers(Map<RNAiGeneAnnotation, IntervalTree<HairpinKmer>> best) {
		Collection set=convertToCollection(best);
		
		return set.size();
	}

	private int count(IntervalTree<HairpinKmer> kmers) {
		if(kmers==null){return 0;}
		Collection<String> set=new TreeSet();
		
		Iterator<Node<HairpinKmer>> iter=kmers.iterator();
		while(iter.hasNext()){
			HairpinKmer kmer=iter.next().getValue();
			set.add(kmer.getKmerSequence());
		}
		
		return set.size();
	}

	//TODO If all we are doing is adding to ones that are already covered then dont add them
	private Map<RNAiGeneAnnotation, IntervalTree<HairpinKmer>> getIsoformCoverage(HairpinKmer hp, Map<RNAiGeneAnnotation, IntervalTree<HairpinKmer>> covered) throws SearchException {
		Collection<HairpinKmer> rtrn=new TreeSet();
		//get all variants covered by new hairpin
		
		for(RNAiGeneAnnotation transcript: hp.getTranscripts()){
			System.err.println(transcript.getTranscriptName());
			boolean transcriptIsDone=false;
			int count=count(covered.get(transcript));
			if(count>=this.minNumber){transcriptIsDone=true;}
			List<SequenceRegion> regions=hp.getRelativePositions(transcript);
			boolean overlaps=false;
			if(covered.containsKey(transcript) && !transcriptIsDone){
				for(SequenceRegion region: regions){
					if(covered.get(transcript).overlappers(region.getStart(), region.getEnd()).hasNext()){overlaps=true;}
				}
			}
			if(!overlaps && !transcriptIsDone){
				for(SequenceRegion region: regions){
					IntervalTree<HairpinKmer> tree=new IntervalTree();
					if(covered.containsKey(transcript)){tree=covered.get(transcript);}
					tree.put(region.getStart(), region.getEnd(), hp);
					covered.put(transcript, tree);
				}
			}
		}
	
		return covered;
	}

	//As a stupid first pass just put them all together
	private Map<String, Collection<HairpinKmer>> collapseVariantsToGenes(Collection<RNAiGeneAnnotation> rnaiAnnotations, Map<String, Collection<HairpinKmer>> hairpins, Map<String, Collection<RNAiGeneAnnotation>> geneVariants) {
		Map<String, Collection<HairpinKmer>> rtrn=new TreeMap();
		
		
		
		for(String gene: geneVariants.keySet()){
			Collection<RNAiGeneAnnotation> variants=geneVariants.get(gene);
			Collection<HairpinKmer> kmers=getHairpinsForGenes(hairpins, variants);
			rtrn.put(gene, kmers);
		}
		
		return rtrn;
	}

	private Map<String, Collection<RNAiGeneAnnotation>> getVariantsPerGene(Collection<RNAiGeneAnnotation> rnaiAnnotations) {
		Map<String, Collection<RNAiGeneAnnotation>> rtrn=new TreeMap();
		
		for(RNAiGeneAnnotation rnai: rnaiAnnotations){
			String gene=rnai.getGeneName();
			String transcript=rnai.getTranscriptName();
			Collection<RNAiGeneAnnotation> list=new HashSet();
			if(rtrn.containsKey(gene)){
				list=rtrn.get(gene);
			}
			list.add(rnai);
			rtrn.put(gene, list);
		}
		
		return rtrn;
	}

	private Collection<HairpinKmer> getHairpinsForGenes(Map<String, Collection<HairpinKmer>> hairpins, Collection<RNAiGeneAnnotation> variants) {
		Collection<HairpinKmer> rtrn=new TreeSet();
		for(RNAiGeneAnnotation variant: variants){
			rtrn.addAll(hairpins.get(variant.getTranscriptName()));
		}
		return rtrn;
	}

	private void write(String save,	Map<String, Collection<HairpinKmer>> geneHairpins, Map<String, Collection<RNAiGeneAnnotation>> geneVariants) throws IOException, SearchException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("GENE.SOURCEID\tTRANS.SOURCEID\tRelative Position\tRS8 Score\tHairpin Sequence\n");
		//write for each transcript version the hairpin sequence and relative position in the sequence
		for(String gene: geneHairpins.keySet()){
			Collection<HairpinKmer> hairpins=geneHairpins.get(gene);
			/*for(HairpinKmer hairpin: hairpins){
				Collection<RNAiGeneAnnotation> transcripts=hairpin.getTranscripts();
				for(RNAiGeneAnnotation transcript: transcripts){
					List<SequenceRegion> relative=hairpin.getRelativePositions(transcript);
					writer.write(transcript.getGeneName()+"\t"+transcript.getTranscriptName()+"\t"+relative.get(0).getStart()+"\t"+hairpin.getRS8Score()+"\t"+hairpin.getKmerSequence()+"\n");
				}
			}*/
			writer.write(gene+"\t"+hairpins.size()+"\n");
		}
		writer.close();
	}


	private Collection<String> makeCollection(Collection<HairpinKmer> hairpins) {
		Collection<String> rtrn=new TreeSet();
		
		for(HairpinKmer kmer: hairpins){
			rtrn.add(kmer.getKmerSequence());
		}
		
		return rtrn;
	}

	public static void main(String[] args)throws Exception{
		if(args.length>2){
		File rnaiFile=new File(args[0]);
		File[] haiprinFiles=new File(args[1]).listFiles();
		String save=args[2];
		new PickHairpinsPerGene(rnaiFile, haiprinFiles, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=RNAi Report File \n args[1]=hairpin files \n args[2]=save";
}
