package broad.pda.rnai.designer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.motif.SearchException;
import broad.core.motif.SequenceMotif;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.rnai.RNAiGeneAnnotation;

public class DesignRNAiHairpinOligos {

	private int k=21;
	int num=5;
	
	Map<String, Double> mirLookup;
	Set<String> lookup2;
	Set<String> rcLookup2;
	List<Sequence> geneSequences;
	private int seed=7;
	Map<String, RefSeqGene> geneMap;
	double minScore=1;
	int numMismatches=3;
	boolean searchHarderForN=true;
	boolean doBlast=true;
	
	public DesignRNAiHairpinOligos(File rnaiReport, File geneFastaFile, File geneCoordinateFile, String saveDir, File mirSeed, File mirLookup)throws Exception{
		initializeMirLookups(mirSeed, mirLookup);
		initializeGeneSequence(geneFastaFile, geneCoordinateFile);
		
		//Parse input file
		Collection<RNAiGeneAnnotation> seqMap=parseRNAiReportFile(rnaiReport);
		
		//enumerate and score kmers
		Map<RNAiGeneAnnotation, List<HairpinKmer>> kmers=scoreKmersPerTranscript(seqMap);
		
		//write kmers for testing with SMATCH
		//writeKmers(saveDir+"/kmers.fa", kmers);
		
		//BLAST Kmers
		//blastKmers(kmers, geneFastaFile);
		
		//Pick optimal oligos
		Map<RNAiGeneAnnotation, Collection<HairpinKmer>> topOligos=pickBestScoringOligos(kmers, num, new ArrayList());
		
		//Check if the number of hairpins per transcript is less than num
		//If so fill in with loser criteria
		Map<RNAiGeneAnnotation, Collection<HairpinKmer>> findRemaining=searchHarder(kmers, topOligos);
		
		write(saveDir, topOligos, "top");
		write(saveDir, findRemaining, "remaining");
		write(saveDir, kmers, "kmers");
	}
	
	
	/*public DesignRNAiHairpinOligos(File rnaiReport, File geneFastaFile, File geneCoordinateFile, String saveDir, File hairpinFile)throws Exception{
		initializeMirLookups();
		initializeGeneSequence(geneFastaFile, geneCoordinateFile);
		
		//Parse input file
		Collection<RNAiGeneAnnotation> seqMap=parseRNAiReportFile(rnaiReport);
		
		//enumerate and score kmers
		Map<RNAiGeneAnnotation, List<HairpinKmer>> kmers=scoreKmersPerTranscript(seqMap);
		
		Collection<String> alreadyChosen=parseHairpinFile(hairpinFile);
		
		//write kmers for testing with SMATCH
		//writeKmers(saveDir+"/kmers.fa", kmers);
		
		//BLAST Kmers
		//blastKmers(kmers, geneFastaFile);
		
		//Pick optimal oligos
		Map<RNAiGeneAnnotation, Collection<HairpinKmer>> topOligos=pickBestScoringOligos(kmers, num, alreadyChosen);
		
		//Check if the number of hairpins per transcript is less than num
		//If so fill in with loser criteria
		Map<RNAiGeneAnnotation, Collection<HairpinKmer>> findRemaining=searchHarder(kmers, topOligos);
		
		write(saveDir, topOligos, "top");
		write(saveDir, findRemaining, "remaining");
		write(saveDir, kmers, "kmers");
	}*/	
	
	
	private Collection<String> parseHairpinFile(File hairpinFile) throws IOException {
		Collection<String> rtrn=new TreeSet();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(hairpinFile)));
		String nextLine;
		int i=0;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			if(i>0){
				String[] tokens=nextLine.split("\t");
				String kmer=tokens[0];
				rtrn.add(kmer);
			}
			i++;
		}
		reader.close();
		
		return rtrn;
	}


	private Map<RNAiGeneAnnotation, Collection<HairpinKmer>> searchHarder(Map<RNAiGeneAnnotation, List<HairpinKmer>> kmers,	Map<RNAiGeneAnnotation, Collection<HairpinKmer>> topOligos) throws SearchException{
		if(!this.searchHarderForN){return topOligos;}
		//Want to find at least n number but didnt so we'll search harder loosening the rules until we do
		Map<RNAiGeneAnnotation, Collection<HairpinKmer>> rtrn=new HashMap();
		
		for(RNAiGeneAnnotation rnai: kmers.keySet()){
			Collection<HairpinKmer> temp=topOligos.get(rnai);
			Collection<HairpinKmer> hairpins=new TreeSet();
			hairpins.addAll(temp);
			System.err.println(hairpins.size()+" "+this.num);
			if(hairpins.size()<this.num){
				Collection<HairpinKmer> additionalHairpins=findAdditional(kmers.get(rnai), hairpins, rnai, 11);
				hairpins.addAll(additionalHairpins);
				additionalHairpins=findAdditional(kmers.get(rnai), hairpins, rnai, 14);
				hairpins.addAll(additionalHairpins);
			}
			rtrn.put(rnai, hairpins);
			
		}
		
		return rtrn;
	}



	private Collection<HairpinKmer> findAdditional(List<HairpinKmer> list, Collection<HairpinKmer> hairpins2, RNAiGeneAnnotation rnai, int overlap) throws SearchException {
		//sort the list
		TreeSet<HairpinKmer> hairpins=new TreeSet<HairpinKmer>(list);
		
		IntervalTree<HairpinKmer> tree=makeTree(hairpins2);
		int numGotten=hairpins2.size();
		System.err.println(numGotten+" "+tree.size());
		for(HairpinKmer hp: hairpins){
			if(numGotten<this.num){
				System.err.println(numGotten+" "+num);
				Iterator<Node<HairpinKmer>> overlappers=tree.overlappers(hp.getKmerStartPosition(), hp.getKmerEndPosition());
				int maxOverlap=maxOverlap(overlappers, hp);
				if(maxOverlap<=overlap){
					if(hp.getRS8Score()>=this.minScore){
						boolean blast=this.blastKmer(hp, rnai);
						if(!blast){
							numGotten++;
							tree.put(hp.getKmerStartPosition(), hp.getKmerEndPosition(), hp);
						}
						else{System.out.println("REJECTED: "+hp.getKmerSequence()+" "+hp.getRS8Score());}
						}
				}
			}
			else{break;}
		}
		Collection<HairpinKmer> kmers= tree.toCollection();
		return kmers;
	}

	private int maxOverlap(Iterator<Node<HairpinKmer>> iter, HairpinKmer kmer){
		int maxOverlap=0;
		while(iter.hasNext()){
			HairpinKmer hp=iter.next().getValue();
			int overlap=kmer.getOverlap(hp);
			maxOverlap=Math.max(overlap, maxOverlap);
		}
		return maxOverlap;
	}

	private IntervalTree<HairpinKmer> makeTree(Collection<HairpinKmer> hairpins2) {
		IntervalTree<HairpinKmer> tree=new IntervalTree();
		if(hairpins2==null){return tree;}
		
		for(HairpinKmer kmer: hairpins2){
			tree.put(kmer.getKmerStartPosition(), kmer.getKmerEndPosition(), kmer);
		}
		
		return tree;
	}
	
	
	


	private void writeKmers(String save, Map<RNAiGeneAnnotation, List<HairpinKmer>> kmers) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(RNAiGeneAnnotation gene: kmers.keySet()){
			List<HairpinKmer> hps=kmers.get(gene);
			for(HairpinKmer kmer: hps){writer.write(">"+kmer.toString()+"\n"+kmer.getKmerSequence()+"\n");}
		}
		
		writer.close();
	}



	private void initializeGeneSequence(File geneSequence, File geneCoordinateFile) throws IOException {
		if(geneSequence==null || !geneSequence.exists()){this.geneSequences=new ArrayList<Sequence>(); System.err.println("Warning: No gene sequences were found, proceeding without specificity test");}
		else{
			FastaSequenceIO fsio = new FastaSequenceIO(geneSequence);
			this.geneSequences = fsio.loadAll();
		}
		
		if(geneCoordinateFile==null || !geneCoordinateFile.exists()){this.geneMap=new TreeMap(); System.err.println("Warning: No gene coordinates were found, no way to ensure not targetting self (might be too conservative)");}
		else{
			this.geneMap=BEDFileParser.loadDataByName(geneCoordinateFile);
		}
	}



	private void blastKmers(Map<RNAiGeneAnnotation, List<HairpinKmer>> kmers, File geneSequence) throws IOException, SearchException {
		FastaSequenceIO fsio = new FastaSequenceIO(geneSequence);
		List<Sequence> seqIt = fsio.loadAll();
		
		for(RNAiGeneAnnotation rnai: kmers.keySet()){
			List<HairpinKmer> hps=kmers.get(rnai);
			for(HairpinKmer kmer: hps){
				//First pass go through and see if there are perfect matches to a gene
				SequenceMotif perfectMotif=new SequenceMotif(kmer.getKmerSequence(), 1);
				if(matches(perfectMotif, seqIt)){kmer.setBlastHit(true);}
				//If none then get all kmers that are n-1 matches
		
				//If none then get all kmers that are n-2 matches
			}
		}
	}
	
	/*private boolean blastKmer(HairpinKmer kmer) throws SearchException{
		//First pass go through and see if there are perfect matches to a gene
		SequenceMotif perfectMotif=new SequenceMotif(kmer.getKmerSequence(), 1);
		if(matches(perfectMotif, this.geneSequences)){return true;}
		//If none then get all kmers that are n-1 matches
		//If none then get all kmers that are n-2 matches
		Collection<SequenceMotif> mismatchedMotif=makeSequenceMotifWithMismatches(kmer.getKmerSequence(), 2);
		if(matches(mismatchedMotif, this.geneSequences)){return true;}
		
		return false;
	}*/
	
	
	private boolean blastKmer(HairpinKmer kmer, RNAiGeneAnnotation rnai)throws SearchException{
		if(!doBlast){return false;}
		System.err.println(kmer.getKmerSequence()+" Forward Scan");
		SmatchLike smatch=new SmatchLike(kmer.getKmerSequence(), this.geneSequences, seed);
		Collection<String> possibleForwardMatches=smatch.getForwardTargets(numMismatches);
		if(possibleForwardMatches.size()>0){System.err.println(kmer.getKmerSequence()+" "+possibleForwardMatches.size()+" "+possibleForwardMatches);}
		//check if possible matches are actually itself
		boolean passes=isNonSelfHit(possibleForwardMatches, rnai, geneMap);
		
		if(passes){System.err.println("Skipping reverse scan"); return passes;}
	
		System.err.println(kmer.getKmerSequence()+" Reverse Scan");
		Collection<String> possibleReverseMatches=smatch.getReverseTargets(numMismatches-1);
		if(possibleReverseMatches.size()>0){System.err.println(kmer.getKmerSequence()+" "+possibleReverseMatches.size()+" "+possibleReverseMatches);}
		//check if possible matches are actually itself
		passes=isNonSelfHit(possibleReverseMatches, rnai, geneMap);
		
		return passes;
	}

	private boolean isNonSelfHit(Collection<String> possibleMatches, RNAiGeneAnnotation rnai, Map<String, RefSeqGene> geneMap2) {
		for(String matches: possibleMatches){
			Alignments region=rnai.getRegion();
			RefSeqGene gene=geneMap2.get(matches);
			if(gene==null){
				gene=geneMap2.get(matches.split("\\.")[0]);
				if(gene==null){System.out.println("Couldnt find coordinates for this RefSeq Gene --> Assuming cross hyb"); return true;}
			}
			if(!region.overlaps(gene.getAlignment())){return true;}
			else{System.out.println(matches+" and "+rnai.getTranscriptName()+" are the same gene");}
		}
		return false;
	}



	//TODO: Replace with SMATCH which might be much faster
	private static Collection<SequenceMotif> makeSequenceMotifWithMismatches(String kmerSequence, int num) throws SearchException {
		Collection<SequenceMotif> rtrn=new ArrayList();
		
		char[] chars=kmerSequence.toCharArray();
		for(int i=0; i<chars.length; i++){
			for(int j=i; j<chars.length; j++){
				if(i!=j){
				String degenerate=string(chars, i, j);
				rtrn.add(new SequenceMotif(degenerate, 1));
				}
			}
		}
		
		return rtrn;
	}

	
	private static String string(char[] chars, int p1, int p2) {
		String rtrn="";
		
		for(int i=0; i<chars.length; i++){
			if(i==p1 || i==p2){rtrn+="N";}
			else{rtrn+=chars[i];}
		}
		
		return rtrn;
	}



	private boolean matches(Collection<SequenceMotif> mismatchedMotif, List<Sequence> geneSequences) {
		int i=0;
		for(SequenceMotif motif: mismatchedMotif){
			System.err.println(i+" "+motif.getMotif());
			if(matches(motif, geneSequences)){return true;}
			i++;
		}
		return false;
	}


	private boolean matches(SequenceMotif motif, List<Sequence> geneSequence){
		for(Sequence seq: geneSequence) {
			if(!seq.getId().startsWith("NR_") || seq.getId().startsWith("nr_")){//TODO: Currently only looking against protein coding ref seq genes
			List<SequenceRegion> matches = motif.match(seq);
			if(matches.size()>0){return true;}
			}
		}
		return false;
	}

	private void initializeMirLookups(File mirSeedLookup, File miRNALookup) throws IOException {
		mirLookup=ComputeMIRScore.parseMirScores(mirSeedLookup);
		lookup2=ComputeMIRScore.parse(miRNALookup, false);
		rcLookup2=ComputeMIRScore.parse(miRNALookup, true);
	}



	private Map<RNAiGeneAnnotation, Collection<HairpinKmer>> pickBestScoringOligos(Map<RNAiGeneAnnotation, List<HairpinKmer>> kmers, int num2, Collection<String> alreadyChosen) throws SearchException {
		Map<RNAiGeneAnnotation, Collection<HairpinKmer>> rtrn=new HashMap();
				
		for(RNAiGeneAnnotation rnai: kmers.keySet()){
			System.err.println(rnai.getTranscriptName());
			List<HairpinKmer> hairpins=kmers.get(rnai);
			Collection<HairpinKmer> topHairpins=pickTopN(hairpins, rnai, num2, alreadyChosen);
			rtrn.put(rnai, topHairpins);
		}
		
		return rtrn;
	}



	private Collection<HairpinKmer> pickTopN(List<HairpinKmer> list, RNAiGeneAnnotation rnai, int num2, Collection<String> alreadyChosen) throws SearchException {
		TreeSet<HairpinKmer> hairpins=new TreeSet(list);
		
		IntervalTree<HairpinKmer> tree=makeTree(alreadyChosen, list);
		int numGotten=0;
		System.err.println("Starting at "+numGotten);
		for(HairpinKmer hp: hairpins){
			if(numGotten<num2){
				if(!tree.overlappers(hp.getKmerStartPosition(), hp.getKmerEndPosition()).hasNext()){
					if(hp.getRS8Score()>=this.minScore){
					boolean blast=this.blastKmer(hp, rnai);
					if(!blast){
						numGotten++;
						tree.put(hp.getKmerStartPosition(), hp.getKmerEndPosition(), hp);
					}
					else{System.out.println("REJECTED: "+hp.getKmerSequence()+" "+hp.getRS8Score());}
					}
				}
			}
			else{break;}
		}
		Collection<HairpinKmer> kmers= tree.toCollection();
		return kmers;
	}



	private IntervalTree<HairpinKmer> makeTree(Collection<String> alreadyChosen, Collection<HairpinKmer> allKmers) {
		IntervalTree<HairpinKmer> rtrn=new IntervalTree();
		
		for(HairpinKmer kmers: allKmers){
			if(alreadyChosen.contains(kmers.getKmerSequence())){
				rtrn.put(kmers.getKmerStartPosition(), kmers.getKmerEndPosition(), kmers);
			}
		}
		
		return rtrn;
	}


	private Map<RNAiGeneAnnotation, List<HairpinKmer>> scoreKmersPerTranscript(Collection<RNAiGeneAnnotation> seqMap){
		Map<RNAiGeneAnnotation, List<HairpinKmer>> scores=new HashMap<RNAiGeneAnnotation, List<HairpinKmer>>();
		for(RNAiGeneAnnotation id: seqMap){
			System.err.println(id.getTranscriptName());
			String seq=id.getTranscriptSequence();;
			List<HairpinKmer> kmerScore=scoreKmers(seq);
			scores.put(id, kmerScore);
		}
		return scores;
	}
	
	private List<HairpinKmer> scoreKmers(String seq){
		return ComputeOriginalScore.enumerateAllKMers(seq, k, mirLookup, lookup2, rcLookup2);
	}
	
	private Collection<RNAiGeneAnnotation> parseRNAiReportFile(File rnaiReport) throws IOException {
		return RNAiFileFormatUtils.parseRNAiReportFile(rnaiReport);
	}
	
	private void write(String saveDir, Map<RNAiGeneAnnotation, ? extends Collection<HairpinKmer>> scores) throws IOException {
		write(saveDir, scores, "");
	}
	
	private void write(String saveDir, Map<RNAiGeneAnnotation, ? extends Collection<HairpinKmer>> scores, String addition) throws IOException {
		for(RNAiGeneAnnotation id: scores.keySet()){
			String save=saveDir+"/"+id.getTranscriptName()+"."+addition+".txt";
			FileWriter writer=new FileWriter(save);
			Collection<HairpinKmer> kmerScores=scores.get(id);
			for(HairpinKmer kmer: kmerScores){
				writer.write(kmer.getKmerSequence()+"\t"+kmer.getKmerStartPosition()+"\t"+kmer.getRS8Score()+"\t"+kmer.getOriginalScore()+"\t"+kmer.getMirScore()+"\n");
			}
			writer.close();
		}
	}
	
	public static void main(String[] args)throws Exception{
		if(args.length>5){
			File file=new File(args[0]);
			File geneFastaFile=new File(args[1]);
			File geneCoordinateFile=new File(args[2]);
			String save=args[3];
			File mirSeed=new File(args[4]);
			File mirLookup=new File(args[5]);
			new DesignRNAiHairpinOligos(file, geneFastaFile, geneCoordinateFile, save, mirSeed, mirLookup);
		}
		else{System.err.println(usage);}		
	}	
	
	/*public static void main(String[] args) throws SearchException{
		String str="ACTGACTG";
		Collection<SequenceMotif> motifs=makeSequenceMotifWithMismatches(str, 2);
		for(SequenceMotif motif: motifs){
			System.err.println(motif.getMotif());
		}
		System.err.println(motifs.size());
	}*/
	
	static String usage=" args[0]=File (RNAi Report) \n args[1]=gene fasta file \n args[2]=gene coordinate file (Full BED) \n args[3]=save directory \n args[4]=mir seed \n args[5]=mir lookup";
	
}
