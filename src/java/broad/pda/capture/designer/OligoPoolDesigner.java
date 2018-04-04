package broad.pda.capture.designer;


import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.motif.SearchException;
import broad.core.primer3.PrimerPair;
import broad.core.primer3.PrimerUtils;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.core.sequence.WindowSlider;
import broad.pda.annotation.BEDFileParser;
import broad.pda.rnai.designer.SmatchLike;
import broad.pda.seq.alignment.Pair;
import jaligner.Alignment;
import jaligner.NeedlemanWunschGotoh;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixGenerator;

public class OligoPoolDesigner {

	int oligoSize=120;
	int overlapDensity=15; //consider making even denser
	
	int numStockNeg=20;
	int matchCutoff=30;
	
	int numberOfProbesPerArray=55000; //because smallest pool is 6,000 and we'll do sense and antisense
	private double TMCutoff=50;
	boolean maskLower=true; //Whether to count lower case as repeat sequence
	private int repeatLimit=50; //Number of repeat bases tolerated	
	
	//Make sure it doesnt pull anything else in transcriptome (maybe be conservative and use genome and transcriptome)
	
	//Things to do:
	//1) BLAST verify uniqueness
	//2) Make sure primers are in correct orientation when we add the complement probe
	//3) Exclude Nanostring probes (if we have them)
	//4) Make sure spacing for multiple isoforms isnt way off--done
	//5) Have option to build primer tails by fold x coverage rather than overlap (will want overlap later but for now)
	//6) Shouldnt be about repeat or not but rather corss-hyb or not.
	//7) Random genes should come from selected genes
	//8) Optimally place on arrays such that no array has more than 60 primers
	
	public OligoPoolDesigner(Collection<Sequence> sequences, List<Sequence> allGenes, Collection<String> negatives, String save) throws Exception{
		//Because some genes have multiple isoforms we'll first classify each into a group by gene name
		Map<String, Collection<Sequence>> sequencesByName=getSequencesByName(sequences);		
		
		//Take regions and design oligos from tiling path
		Map<Sequence, Map<Integer, Collection<Sequence>>> initialPools=designTilingPaths(sequences); 
		
		//Filter repetitive sequences (for now just kill any probe with more than 30 Ns in it)
		initialPools=filterRepetitive(initialPools);
		
		//writeSequences(save+".temp.fa", initialPools);
		
		for(String neg: negatives){
			if(!sequencesByName.containsKey(neg)){System.err.println("ERROR "+neg);}
		}
		
		
		//Check if probes are redundant within a gene
		Map<String, Map<Integer, Collection<Sequence>>> pools=filterRedundantProbes(sequencesByName, initialPools);
		
		//Design artificial negative controls off of the same sequence properties as random picked genes
		pools.putAll(this.designScrambledPools(sequencesByName, negatives));
		
		System.err.println("Finished making pools");
		
		//Test if probes cross-hyb in genome or transcriptome
		//pools=filterTilingPaths(pools, allGenes);//TODO Took this out and will BLAT verify all at the end against transcriptome and genome
		
		
		System.err.println("Filtered probes");
		
		//Design PCR tails with no cross-priming or complementarity to genes
		Map<PrimerPair, Collection<String>> primers=designPrimers(pools, allGenes);

		FileWriter writer=new FileWriter("primers.txt");
		for(PrimerPair primer: primers.keySet()){
			Collection<String> sub=primers.get(primer);
			writer.write(primer.getLeftPrimer()+"\t"+primer.getRightPrimer());
			for(String s: sub){writer.write("\t"+s);}
			writer.write("\n");
		}
		writer.close();
		
		System.err.println("Designed primers");
		
		//Assign sequences to arrays
		Map<Integer, Collection<String>> genesPerArray=assignSequencesToArrays(pools);
		
		System.err.println("Assigned genes to array total of "+genesPerArray.size());
		
		//Assign primers to genes
		Map<String, PrimerPair> primersPerSequence=assignPrimersToGenes(pools, primers, genesPerArray);
		
		System.err.println("Assigned primers to gene");
		
		//Write output file
		writeReport(save, pools, genesPerArray, primersPerSequence, primers);
		
		System.err.println("wrote output");
				
		//TODO For each primer in an array pool, check what possible amplicons exist
		
		
		//Add new subtail to probes in nanostring set
		
		//TODO Add a BLAST (or BLAT) step
	}

	
	private void writeSequences(String save, Map<Sequence, Map<Integer, Collection<Sequence>>> initialPools) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		int counter=0;
		for(Sequence seq: initialPools.keySet()){
			Map<Integer, Collection<Sequence>> path=initialPools.get(seq);
			for(Integer num: path.keySet()){
				int counter2=0;
				for(Sequence probe: path.get(num)){
					writer.write(">"+seq.getId()+"_isoform"+counter+"_"+num+"_"+counter2+"\n"+probe.getSequenceBases()+"\n");
					counter2++;
				}
			}
			
			counter++;
		}
		
		writer.close();
	}


	private Map<String, Collection<Sequence>> getSequencesByName(Collection<Sequence> sequences) {
		Map<String, Collection<Sequence>> rtrn=new TreeMap<String, Collection<Sequence>>();
		
		
		for(Sequence seq: sequences){
			Collection<Sequence> list=new ArrayList<Sequence>();
			String name=seq.getId();
			if(rtrn.containsKey(name)){list=rtrn.get(name);}
			list.add(seq);
			rtrn.put(name, list);
		}
		
		
		return rtrn;
	}


	private Map<String, Map<Integer, Collection<Sequence>>> filterRedundantProbes(Map<String, Collection<Sequence>> sequencesByName, Map<Sequence, Map<Integer, Collection<Sequence>>> pools) {
		Map<String, Map<Integer, Collection<Sequence>>> rtrn=new TreeMap<String, Map<Integer, Collection<Sequence>>>(); 
		//Returned map will correspond to the gene not isoforms now
		//Also, this is where we'll sync up the tiling path numbers so that there are independent non-overlapping tiling paths
		
		for(String gene: sequencesByName.keySet()){
			Collection<Sequence> sequences=sequencesByName.get(gene);
			System.err.println("Reducing gene "+gene+" "+sequences.size());
			//Here we will collapse the tiling paths into unique probes and synchronize the tiling paths so that they are non-overlapping
			Map<Integer, Collection<Sequence>> geneProbes=getUniqueAndOptimalTilingPath(sequences, pools);
			rtrn.put(gene, geneProbes);
		}
		
		return rtrn;
	} 


	//get all tiling paths for a single gene
	//filter non unique probes 
	//synchronize tiling paths
	private Map<Integer, Collection<Sequence>> getUniqueAndOptimalTilingPath(Collection<Sequence> sequences, Map<Sequence, Map<Integer, Collection<Sequence>>> pools) {
		Map<Integer, Collection<Sequence>> rtrn=new TreeMap<Integer, Collection<Sequence>>();
		
		Map<String, Integer> allProbesSoFar=new TreeMap<String, Integer>();
		//Pick a reference isoform (first?)
		//Make this the baseline tiling path
		//For each other get all probes, if completely overlaps a probe already designed then ignore
		//If has overlap less than cutoff, keep but put in same tiling path as the overlapper
		
		for(Sequence isoform: sequences){
			Map<Integer, Collection<Sequence>> tilingPath=pools.get(isoform);
			if(rtrn.isEmpty()){rtrn=tilingPath;}
			else{
				for(Integer poolNum: tilingPath.keySet()){
					Collection<Sequence> allProbes=tilingPath.get(poolNum);
					Map<Sequence, Integer> probeAssignment=assignProbes(allProbes, allProbesSoFar, poolNum);
					for(Sequence probe: probeAssignment.keySet()){
						int assignedPoolNum=probeAssignment.get(probe);
						Collection<Sequence> probes=rtrn.get(assignedPoolNum);
						
						if(probes!=null){
							probes.add(probe);
						}
						else{
							probes=new ArrayList<Sequence>();
							probes.add(probe);
						}
						rtrn.put(assignedPoolNum, probes);
					}
				}
			}
			updateAllProbes(rtrn, allProbesSoFar);
		}
		
		return rtrn;
	}


	private void updateAllProbes(Map<Integer, Collection<Sequence>> rtrn, Map<String, Integer> allProbesSoFar) {
		for(Integer poolNum: rtrn.keySet()){
			Collection<Sequence> seqs=rtrn.get(poolNum);
			for(Sequence seq: seqs){
				allProbesSoFar.put(seq.getSequenceBases(), poolNum);
			}
		}
	}


	private Map<Sequence, Integer> assignProbes(Collection<Sequence> allTilingProbes, Map<String, Integer> allProbesAlreadyUsed, Integer poolNum) {
		Map<Sequence, Integer> rtrn=new HashMap<Sequence, Integer>();
		
		//For each probe, if completely overlaps a probe already designed then ignore
		//If has overlap less than cutoff, keep but put in same tiling path as the overlapper
		
		for(Sequence probe: allTilingProbes){
			if(!allProbesAlreadyUsed.containsKey(probe.getSequenceBases())){
				//Step 1: Align this probe to all already used
				boolean hasPerfect=false;
				boolean hasGood=false;
				int goodIndexNumber=-1;
				double bestPercent=-1;
				for(String other: allProbesAlreadyUsed.keySet()){
					//local alignment
					Alignment align=localAlign(probe.getSequenceBases(), other);
					double percentMatch=(double)align.getSequence1().length/(double)this.oligoSize;
					//TODO Think about what todo with probes that have good overlap match with 2 tiling paths
					if(percentMatch>.30 && percentMatch<.9){//Has to be greater than the defined overlap or else its meaningless
						//System.err.println("Probes have decent match "+probe.getSequenceBases()+" "+other+" "+percentMatch+" "+align.getSequence1().length);
						int num=allProbesAlreadyUsed.get(other);
						hasGood=true;
						if(percentMatch>bestPercent){
							bestPercent=percentMatch;
							goodIndexNumber=num;
						}
					
						//rtrn.put(probe, num);
					}
					else if(percentMatch>=.9){
					//	System.err.println("Probes have great match "+probe.getSequenceBases()+" "+other+" "+percentMatch +" "+align.getSequence1().length);
						hasPerfect=true;
					}
				}
				if(!hasPerfect){
					if(hasGood){
						//System.err.println("Reassigned "+probe.getSequenceBases()+" from pool "+poolNum+" to "+goodIndexNumber);
						rtrn.put(probe, goodIndexNumber);
					}
					else{
						rtrn.put(probe, poolNum);
					}
				}
			}
		}
		
		
		return rtrn;
	}


	private Map<Sequence, Map<Integer, Collection<Sequence>>> filterRepetitive(Map<Sequence, Map<Integer, Collection<Sequence>>> pools) {
		Map<Sequence, Map<Integer, Collection<Sequence>>> rtrn=new HashMap<Sequence, Map<Integer, Collection<Sequence>>>();
		
		for(Sequence gene: pools.keySet()){
			Map<Integer, Collection<Sequence>> map=pools.get(gene);
			Map<Integer, Collection<Sequence>> filteredMap=new TreeMap<Integer, Collection<Sequence>>();
			for(Integer poolNum: map.keySet()){
				Collection<Sequence> list=map.get(poolNum);
				Collection<Sequence> filtered=new ArrayList<Sequence>();
				for(Sequence seq: list){
					int numRepeat=countRepeats(seq);
					if(numRepeat<repeatLimit){filtered.add(seq);}
					//else{System.out.println(seq.getSequenceBases()+"\t"+numRepeat);}
				}
				//System.err.println(list.size()+" "+filtered.size());
				filteredMap.put(poolNum, filtered);
				//System.err.println("Gene "+gene.getId()+" pool num "+poolNum+" has "+filtered.size()+" reduced from "+list.size());
			}
			rtrn.put(gene, filteredMap);
		}
		
		
		return rtrn;
	}


	private int countRepeats(Sequence seq) {
		char[] bases=seq.getSequenceBases().toCharArray();
		int good=0;
		
		for(int i=0; i<bases.length; i++){
			if(bases[i]=='A' || bases[i]=='C' || bases[i]=='G' ||bases[i]=='T'){good++;}
			else if(!this.maskLower &&	(bases[i]=='a' || bases[i]=='c' || bases[i]=='g' ||bases[i]=='t')){good++;}
		}
		
		int count=bases.length-good;
		
		return count;
	}


	private void writeReport(String save, Map<String, Map<Integer, Collection<Sequence>>> pools, Map<Integer, Collection<String>> genesPerArray, Map<String, PrimerPair> primersPerSequence, Map<PrimerPair, Collection<String>> primers) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("Gene\tArray number\tTiling path number\tGene left primer \tGene right primer\tTiling Path left primer\tFull Probe Sequence\tNumber repeat bases\n");
		
		FileWriter writer2=new FileWriter(save+".fa");
		
		for(Integer arrayNum: genesPerArray.keySet()){
			Collection<String> genes=genesPerArray.get(arrayNum);
			for(String gene: genes){
				PrimerPair primer=primersPerSequence.get(gene);
				Iterator<String> subprimers=primers.get(primer).iterator();
				Map<Integer, Collection<Sequence>> tilingPath=pools.get(gene);
				int numProbes=0;
				for(Integer tilingPathNumber: tilingPath.keySet()){
					int counter=0;
					Collection<Sequence> probes=tilingPath.get(tilingPathNumber);
					String tilingPrimer=subprimers.next();
					for(Sequence probe: probes){
						String fullProbeSeq=(primer.getLeftPrimer()+Sequence.get3Prime(tilingPrimer, 5)+probe.getSequenceBases()+Sequence.reverseSequence(primer.getRightPrimer())).toUpperCase();
						writer.write(gene+"\t"+arrayNum+"\t"+tilingPathNumber+"\t"+primer.getLeftPrimer()+"\t"+primer.getRightPrimer()+"\t"+tilingPrimer+"\t"+fullProbeSeq+"\t"+countRepeats(probe)+"\n");
						writer2.write(">"+gene+"_"+tilingPathNumber+"_"+(counter++)+"\n"+probe.getSequenceBases()+"\n");
						numProbes++;
					}
				}
				System.err.println(gene+" had "+numProbes+" probes");
			}
		}
		
		writer2.close();
		writer.close();
	}


	private Map<String, PrimerPair> assignPrimersToGenes(Map<String, Map<Integer, Collection<Sequence>>> pools,	Map<PrimerPair, Collection<String>> originalPrimers, Map<Integer, Collection<String>> genesPerArray) {
		Map<String, PrimerPair> rtrn=new HashMap<String, PrimerPair>();
		
		for(Integer arrayNum: genesPerArray.keySet()){
			Map<PrimerPair, Collection<String>> primers=copy(originalPrimers);
			Collection<String> genes=genesPerArray.get(arrayNum);
			for(String gene: genes){
				Map<Integer, Collection<Sequence>> pool=pools.get(gene);
				//System.err.println(gene.getId()+" "+primers.size()+" "+pool);
				PrimerPair primer=getPrimers(primers, pool.size()); //remove used primers from map
				if(primer==null){throw new IllegalArgumentException("There arent enough primers to support this");}
				rtrn.put(gene, primer);
			}
		}
		
		
		
		return rtrn;
	}


	private Map<PrimerPair, Collection<String>> copy(Map<PrimerPair, Collection<String>> originalPrimers) {
		Map<PrimerPair, Collection<String>> rtrn=new TreeMap<PrimerPair, Collection<String>>();
		
		for(PrimerPair primer: originalPrimers.keySet()){
			Collection<String> primers=originalPrimers.get(primer);
			Collection<String> list=new TreeSet<String>();
			for(String p: primers){list.add(p);}
			rtrn.put(primer, list);
		}
		
		return rtrn;
	}


	private PrimerPair getPrimers(Map<PrimerPair, Collection<String>> primers, int size) {
		
		for(PrimerPair primer: primers.keySet()){
			Collection<String> sub=primers.get(primer);
			if(sub.size()>size){
				primers.remove(primer);
				return primer;
			}
		}
		
		//There arent enough primers
		//TODO Be smart and make more
		return null;
	}


	private Map<Integer, Collection<String>> assignSequencesToArrays(Map<String, Map<Integer, Collection<Sequence>>> pools) {
		Map<Integer, Collection<String>> genesPerArray=new TreeMap<Integer, Collection<String>>();
		
		int totalProbes=0;
		for(String gene: pools.keySet()){
			Map<Integer, Collection<Sequence>> pool=pools.get(gene);
			int probesPerGene=0;
			for(Integer num: pool.keySet()){probesPerGene+=pool.get(num).size();}
			totalProbes+=probesPerGene;
			int arrayNum=totalProbes/ this.numberOfProbesPerArray;
			Collection<String> list=new ArrayList<String>();
			if(genesPerArray.containsKey(arrayNum)){list=genesPerArray.get(arrayNum);}
			list.add(gene);
			genesPerArray.put(arrayNum, list);
		}
		
		return genesPerArray;
	}


	//Test if probes cross-hyb in genome or transcriptome
	private Map<Sequence, Map<Integer, Collection<Sequence>>> filterTilingPaths(Map<Sequence, Map<Integer, Collection<Sequence>>> pools, List<Sequence> allGenes) throws SearchException {
		Map<Sequence, Map<Integer, Collection<Sequence>>> rtrn=new HashMap<Sequence, Map<Integer, Collection<Sequence>>>();
		
		for(Sequence gene: pools.keySet()){
			System.err.println("Searching for cross-hybs for gene "+gene.getId());
			Map<Integer, Collection<Sequence>> map=new TreeMap<Integer, Collection<Sequence>>();
			for(Integer num: pools.get(gene).keySet()){
				Collection<Sequence> filtered=new HashSet<Sequence>();
				for(Sequence seq: pools.get(gene).get(num)){
					boolean crossHybs=crossHybs(seq, allGenes, gene);
					if(!crossHybs){filtered.add(seq);}
				}
				map.put(num, filtered);
				System.err.println(gene.getId()+" For pool number "+num+" started with "+pools.get(gene).get(num).size()+" and ended with "+filtered.size());
			}
			rtrn.put(gene, map);
		}
		return pools;
	}


	//Does this sequence match anywhere in all genes (whe
	private boolean crossHybs(Sequence seq, List<Sequence> allGenes, Sequence gene) throws SearchException {
		//TODO Test if there is a perfect cross-hyb of minMatch
		
		String kmer1=seq.getSequenceBases();
		SmatchLike smatch=new SmatchLike(kmer1, allGenes, this.matchCutoff);
		Collection<String> crossHybs=smatch.getAllPossibleTargets(0); //TODO Can speed this up by ending as soon as non-self found
		
		if(crossHybs.size()>1){
			System.err.println(crossHybs);
			return true;
		}
		
		return false;
		
			/*else{
				for(String subprimer: primers.get(primer)){
					String kmer3=subprimer;
					smatch=new SmatchLike(kmer3, allGenes, kmer3.length());
					int size3=smatch.getAllPossibleTargets(0).size();
					if(size3>0){excludedSubprimers.add(subprimer);}
				}
			}*/
				
		
		
		/*for(Sequence other: allGenes){
			if(!other.equals(gene)){
				boolean crossHyb=crossHyb(seq, other);
				if(crossHyb){return true;}
			}
		}
		return false;*/
	}

	private boolean crossHyb(Sequence seq, Sequence other) {
		//For now, align using a local alignment and filter if has more than N matches
		Pair<Alignment> alignments=align(seq, other);
		
		//TO just perfect seed match
		
		Alignment sense=alignments.getValue1();
		Alignment antisense=alignments.getValue2();
		
		if(sense.getNumberOfMatches()>this.matchCutoff || antisense.getNumberOfMatches()>this.matchCutoff){
			double senseTM=PrimerUtils.computeMaxTM(sense);
			double antisenseTM=PrimerUtils.computeMaxTM(antisense);
			if(senseTM>this.TMCutoff || antisenseTM>this.TMCutoff){
				System.err.println(seq.getId()+" and "+other.getId()+" crosshybed with "+sense.getNumberOfMatches()+" "+antisense.getNumberOfMatches() +" "+senseTM+" "+antisenseTM);
				return true;
			}
		}
		return false;
	}


	private Pair<Alignment> align(Sequence seq, Sequence other) {
		jaligner.Sequence s1=new jaligner.Sequence(Sequence.reverseSequence(seq.getSequenceBases()));
		jaligner.Sequence s3=new jaligner.Sequence((seq.getSequenceBases()));
		jaligner.Sequence s2=new jaligner.Sequence(other.getSequenceBases());
		Matrix matrix=MatrixGenerator.generate(1.0f, -1.0f);
		jaligner.Alignment alignment=SmithWatermanGotoh.align(s1, s2, matrix, 1000000, 1000000);
		jaligner.Alignment alignment2=SmithWatermanGotoh.align(s3, s2, matrix, 1000000, 1000000);
		return new Pair<Alignment>(alignment, alignment2);
	}
	
	private Alignment localAlign(String seq1, String seq2) {
		jaligner.Sequence s1=new jaligner.Sequence(seq1);
		jaligner.Sequence s2=new jaligner.Sequence(seq2);
		Matrix matrix=MatrixGenerator.generate(1.0f, -100000.0f); //To get the perfect match which we expect here
		jaligner.Alignment alignment=SmithWatermanGotoh.align(s1, s2, matrix, 1000000, 1000000);
		return alignment;
	}
	
	private Pair<Alignment> globalAlign(Sequence seq, Sequence other) {
		jaligner.Sequence s1=new jaligner.Sequence(Sequence.reverseSequence(seq.getSequenceBases()));
		jaligner.Sequence s3=new jaligner.Sequence((seq.getSequenceBases()));
		jaligner.Sequence s2=new jaligner.Sequence(other.getSequenceBases());
		Matrix matrix=MatrixGenerator.generate(1.0f, -1.0f);
		jaligner.Alignment alignment=NeedlemanWunschGotoh.align(s1, s2, matrix, 1000000, 1000000);
		jaligner.Alignment alignment2=NeedlemanWunschGotoh.align(s3, s2, matrix, 1000000, 1000000);
		return new Pair<Alignment>(alignment, alignment2);
	}


	private Map<PrimerPair, Collection<String>> designPrimers(Map<String, Map<Integer, Collection<Sequence>>> pools, List<Sequence> allGenes) throws Exception {
		int numberOfGenesPerArray=getMaxNumberOfGenesPerArray(pools);
		
		//System.err.println(pools.size()+" "+numberOfGenesPerArray);
		
		DesignPCRTails pcr=new DesignPCRTails(new Double((numberOfGenesPerArray+1)*2).intValue());
		
		
		Map<PrimerPair, Collection<String>> primers=pcr.getPrimers();
		primers=filterPerfectCrossHyb(primers, allGenes); //TODO Add back to speed up
		return primers;
	}

	
	private Map<PrimerPair, Collection<String>> filterPerfectCrossHyb(Map<PrimerPair, Collection<String>> primers, List<Sequence> allGenes) throws SearchException {
		Map<PrimerPair, Collection<String>> rtrn=new TreeMap<PrimerPair, Collection<String>>();
		
		
		Collection<PrimerPair> excludedPrimers=new TreeSet<PrimerPair>();
		Collection<String> excludedSubprimers=new TreeSet<String>();
		
		for(PrimerPair primer: primers.keySet()){
			String kmer1=Sequence.get3Prime(primer.getLeftPrimer(),15);
			String kmer2=primer.getRightPrimer();
			SmatchLike smatch=new SmatchLike(kmer1, allGenes, kmer1.length());
			Collection<String> set1=smatch.getAllPossibleTargets(0);
			smatch=new SmatchLike(kmer2, allGenes, kmer2.length());
			Collection<String> set2=smatch.getAllPossibleTargets(0);
			if(set1.size()>0 || set2.size()>0){excludedPrimers.add(primer);}
			/*else{
				for(String subprimer: primers.get(primer)){
					String kmer3=subprimer;
					smatch=new SmatchLike(kmer3, allGenes, kmer3.length());
					int size3=smatch.getAllPossibleTargets(0).size();
					if(size3>0){excludedSubprimers.add(subprimer);}
				}
			}*/
		}
		
		int excludedPrimerCounter=0;
		for(PrimerPair primer: primers.keySet()){
			if(!excludedPrimers.contains(primer)){
				Collection<String> subprimers=primers.get(primer);
				Collection<String> list=new TreeSet<String>();
				for(String sub: subprimers){
					if(!excludedSubprimers.contains(sub)){
						list.add(sub);
					}
				}
				rtrn.put(primer, list);
			}
			else{excludedPrimerCounter++;}
		}
		
		System.err.println("Excluded "+excludedPrimerCounter+" primers due to high cross-hybs with genes");
		
		return rtrn;
	}



	private int getMaxNumberOfGenesPerArray(Map<String, Map<Integer, Collection<Sequence>>> pools){
		Map<Integer, Collection<String>> genesPerArray=getNumberOfGenesPerArray(pools);
		int max=-1;
		for(Integer num: genesPerArray.keySet()){
			max=Math.max(max, genesPerArray.get(num).size());
		}
		return max;
	}


	private Map<Integer, Collection<String>> getNumberOfGenesPerArray(Map<String, Map<Integer, Collection<Sequence>>> pools) {
		Map<Integer, Collection<String>> genesPerArray=new TreeMap<Integer, Collection<String>>();
		
		int totalProbes=0;
		for(String gene: pools.keySet()){
			Map<Integer, Collection<Sequence>> pool=pools.get(gene);
			int probesPerGene=0;
			for(Integer num: pool.keySet()){probesPerGene+=pool.get(num).size();}
			totalProbes+=probesPerGene;
			int arrayNum=totalProbes/this.numberOfProbesPerArray;
			Collection<String> list=new HashSet<String>();
			if(genesPerArray.containsKey(arrayNum)){list=genesPerArray.get(arrayNum);}
			list.add(gene);
			genesPerArray.put(arrayNum, list);
		}
		
		return genesPerArray;
	}
	
	

	private Map<Sequence, Map<Integer, Collection<Sequence>>> filterCrossHybs(Map<Sequence, Map<Integer, Collection<Sequence>>> pools,Collection<Sequence> allGenes) {
		Map<String, Alignment> alignmentScores=align(pools, allGenes);
		Map<Sequence, Map<Integer, Collection<Sequence>>> filteredPools=filterSequences(alignmentScores, pools);
		return filteredPools;
	}


	//TODO
	/*private Map<Pair<String>, Collection<String>> filterPrimers(Map<Pair<String>, Collection<String>> primers,	Map<Sequence, Map<Integer, Collection<Sequence>>> filteredPools, Collection<Sequence> allGenes) {
		//For each primer make sure that it wont amplify internally to a probe
		
		//ensure tails dont crosshyb
		Collection<String> tails=new TreeSet<String>();
		Collection<String> rightPrimer=new TreeSet<String>();
		for(Pair<String> primer: primers.keySet()){
			Collection<String> barcodes=primers.get(primer);
			for(String barcode: barcodes){
				String seq1=primer.getValue1()+barcode;
				tails.add(seq1);
			}
			rightPrimer.add(primer.getValue2());
		}
		
		tails=filterCrossHybs(tails, allGenes);
		rightPrimer=filterCrossHybs(rightPrimer, allGenes);
	}*/



	private void assignTails(String save,Map<Sequence, Map<Integer, Collection<Sequence>>> filteredPools,Map<Pair<String>, Collection<String>> filteredPrimers) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<Pair<String>> usedPrimers=new HashSet<Pair<String>>();
		
		for(Sequence geneSequence: filteredPools.keySet()){
			Map<Integer, Collection<Sequence>> probes=filteredPools.get(geneSequence);
			
			int poolSize=probes.keySet().size();
			Pair<String> primers=getPrimers(filteredPrimers, poolSize, usedPrimers);
			usedPrimers.add(primers);
			Iterator<String> barcodes=filteredPrimers.get(primers).iterator();
			for(Integer poolNum: probes.keySet()){
				String barcode=barcodes.next();
				Collection<Sequence> probeSeqs=filteredPools.get(geneSequence).get(poolNum);
				for(Sequence probeSeq: probeSeqs){
					String fullProbeSeq=primers.getValue1()+barcode+probeSeq.getSequenceBases()+Sequence.reverseSequence(primers.getValue2());
					writer.write(geneSequence.getId()+"\t"+poolNum+"\t"+primers.getValue1()+"\t"+primers.getValue2()+"\t"+barcode+"\t"+probeSeq.getSequenceBases()+"\t"+fullProbeSeq+"\n");
				}
			}
		}
		
		writer.close();
	}



	private Pair<String> getPrimers(Map<Pair<String>, Collection<String>> filteredPrimers, int poolSize, Collection<Pair<String>> usedPrimers) {
		for(Pair<String> primers: filteredPrimers.keySet()){
			Collection<String> barcodes=filteredPrimers.get(primers);
			if(barcodes.size()>=poolSize && !usedPrimers.contains(primers)){return primers;}
		}
		return null;
	}



	private Map<Sequence, Map<Integer, Collection<Sequence>>> filterSequences(Map<String, Alignment> alignmentScores, Map<Sequence, Map<Integer, Collection<Sequence>>> pools) {
		Map<Sequence, Map<Integer, Collection<Sequence>>> rtrn=new HashMap();
		
		for(Sequence seq: pools.keySet()){
			Map<Integer, Collection<Sequence>> sub=new TreeMap();
			for(Integer num: pools.get(seq).keySet()){
				Collection<Sequence> seqs=pools.get(seq).get(num);
				Collection<Sequence> filtered=filterProbes(seqs, alignmentScores);
				sub.put(num, filtered);
			}
			rtrn.put(seq, sub);
		}
		
		
		
		
		return rtrn;
	}

	private Collection<Sequence> filterProbes(Collection<Sequence> probes,	Map<String, Alignment> alignmentScores) {
		Collection<Sequence> rtrn=new ArrayList<Sequence>();
		for(Sequence probe: probes){
			//System.err.println(probe.getSequenceBases());
			//System.err.println(alignmentScores.keySet());
			Alignment alignment=alignmentScores.get(probe.getSequenceBases());
			if(alignment.getNumberOfMatches()<matchCutoff){rtrn.add(probe);}
			else{System.err.println("Filtered "+probe.getId()+" "+probe.getSequenceBases()+" "+alignment.getNumberOfMatches());}
		}
		return rtrn;
	}

	private void write(String save, Map<Sequence, Map<Integer, Collection<Sequence>>> filteredPools) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Sequence seq: filteredPools.keySet()){
			for(Integer poolNum: filteredPools.get(seq).keySet()){
				for(Sequence oligo: filteredPools.get(seq).get(poolNum)){
					writer.write(seq.getId()+"\t"+poolNum+"\t"+oligo.getSequenceBases()+"\t"+Sequence.reverseSequence(oligo.getSequenceBases())+"\n");
				}
			}
		}
		
		writer.close();
	}

	private Map<String, Alignment> align(Map<Sequence, Map<Integer, Collection<Sequence>>> pools, Collection<Sequence> sequences) {
		Map<String, Alignment> scores=new HashMap<String, Alignment>();
		for(Sequence s: pools.keySet()){
			for(Integer n: pools.get(s).keySet()){
				for(Sequence random: pools.get(s).get(n)){
					Alignment align=alignBestNonSelf(random, sequences, s.getId());
					scores.put(random.getSequenceBases(), align);
				}
			}
		}
		return scores;
	}
	
	private Map<String, Alignment> align(Collection<String> primers, Collection<Sequence> sequences) {
		Map<String, Alignment> scores=new HashMap<String, Alignment>();
		
		for(String random: primers){
			Alignment align=alignBestNonSelf(random, sequences);
			scores.put(random, align);
		}
		return scores;
	}
	
	

	private Alignment alignBestNonSelf(Sequence random, Collection<Sequence> sequences, String name) {
		double maxScore=-1;
		jaligner.Alignment best=null;
		for(Sequence seq: sequences){
			//Perform pairwise alignment
			jaligner.Sequence s1=new jaligner.Sequence(Sequence.reverseSequence(random.getSequenceBases()));
			jaligner.Sequence s3=new jaligner.Sequence((random.getSequenceBases()));
			jaligner.Sequence s2=new jaligner.Sequence(seq.getSequenceBases());
			Matrix matrix=MatrixGenerator.generate(1.0f, -1.0f);
			jaligner.Alignment alignment=SmithWatermanGotoh.align(s1, s2, matrix, 1000000, 1000000);
			jaligner.Alignment alignment2=SmithWatermanGotoh.align(s3, s2, matrix, 1000000, 1000000);
			float score=alignment.getScore();
			if(score>maxScore && !(name.equalsIgnoreCase(seq.getId()))){best=alignment; maxScore=score;}
			if(alignment2.getScore()>maxScore){best=alignment2; maxScore=alignment2.getScore();}
		}
		return best;
	}
	
	private Alignment alignBestNonSelf(String random, Collection<Sequence> sequences) {
		double maxScore=-1;
		jaligner.Alignment best=null;
		for(Sequence seq: sequences){
			//Perform pairwise alignment
			jaligner.Sequence s1=new jaligner.Sequence(Sequence.reverseSequence(random));
			jaligner.Sequence s3=new jaligner.Sequence((random));
			jaligner.Sequence s2=new jaligner.Sequence(seq.getSequenceBases());
			Matrix matrix=MatrixGenerator.generate(1.0f, -1.0f);
			jaligner.Alignment alignment=SmithWatermanGotoh.align(s1, s2, matrix, 1000000, 1000000);
			jaligner.Alignment alignment2=SmithWatermanGotoh.align(s3, s2, matrix, 1000000, 1000000);
			float score=alignment.getScore();
			if(score>maxScore){best=alignment; maxScore=score;}
			if(alignment2.getScore()>maxScore){best=alignment2; maxScore=alignment2.getScore();}
		}
		return best;
	}
	
	private Map<Integer, Collection<Sequence>> getPool(Sequence seq){
		Map<Integer, Collection<Sequence>> pool=new HashMap<Integer, Collection<Sequence>>();
		int counter=0;
		WindowSlider slider=WindowSlider.getSlider(seq, oligoSize, oligoSize-this.overlapDensity);
		while(slider.hasNext()){
			SequenceRegion region=slider.next();
			int poolNumber=(counter%oligoSize);
			Collection<Sequence> set=new ArrayList<Sequence>();
			if(pool.containsKey(poolNumber)){set=pool.get(poolNumber);}
			Sequence probe=region.getSequence();
			Sequence sense=new Sequence(seq.getId());
			sense.setSequenceBases(Sequence.reverseSequence(probe.getSequenceBases()));
			set.add(sense);
			//System.err.println(poolNumber+" "+counter+" "+sense.getSequenceBases());
			pool.put(poolNumber, set);
			counter+=this.overlapDensity;
		}
		return pool;
	}

	private Collection<Sequence> designStockNegatives(int numStockNeg) {
		Collection<Sequence> rtrn=new ArrayList<Sequence>();
		
		for(int i=0; i<numStockNeg; i++){
			String string=Sequence.generateRandomSequence(this.oligoSize);
			Sequence seq=new Sequence("random_"+i);
			seq.setSequenceBases(string);
			rtrn.add(seq);
		}
		
		return rtrn;
	}

	private Map<String, Map<Integer, Collection<Sequence>>> designScrambledPools(Map<String, Collection<Sequence>> sequences, Collection<String> negatives) {
		Map<String, Map<Integer, Collection<Sequence>>> pools=new HashMap<String, Map<Integer, Collection<Sequence>>>();
		//We'll scramble each probe to equalize GC content per probe
		
		for(String seq: negatives){
			Collection<Sequence> seqs=sequences.get(seq);
			Sequence scrambled=seqs.iterator().next().scramble();
			Map<Integer, Collection<Sequence>> pool=getPool(scrambled);
			pools.put(seq+"_Random", pool);
		}
		
		System.err.println("Random pool size "+pools.size());
		
		return pools;
	}

	private Map<Sequence, Map<Integer, Collection<Sequence>>> designTilingPaths(Collection<Sequence> sequences) {
		Map<Sequence, Map<Integer, Collection<Sequence>>> pools=new HashMap<Sequence, Map<Integer, Collection<Sequence>>>();
		for(Sequence seq: sequences){
			Map<Integer, Collection<Sequence>> pool=getPool(seq);
			pools.put(seq, pool);
		}
		return pools;
	}
	
	public static void main(String[] args) throws Exception{
		if(args.length>2){
			FastaSequenceIO fsio = new FastaSequenceIO(new File(args[0]));
			Collection<Sequence> seq= fsio.loadAll();
			fsio= new FastaSequenceIO(new File(args[1]));
			List<Sequence> all= fsio.loadAll();
			all.addAll(seq);
			Collection<String> negatives=parseList(args[2]);
			String save=args[3];
			new OligoPoolDesigner(seq, all,negatives, save);
		}
		else{
			System.err.println(usage);
		}
	}
	
	private static Collection<String> parseList(String string) throws IOException {
		Collection<String> rtrn=BEDFileParser.loadList(string);
		return rtrn;
	}

	static String usage=" args[0]=seqeunces \n args[1]=all genes \n args[2]=negatives \n args[3]=save";
	
}
