package broad.pda.seq.pairedend;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.math.EmpiricalDistribution;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.AlignmentDataModelStats;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.utils.GenesInExpressionBins;

public class ConnectLastExonToSingleExonsWithPairedEnds {

	private String chr;
	private int maxDistanceForSingleExonsAfterLastExon;
	private Collection<RefSeqGene> genes;
	private ContinuousDataAlignmentModel data;
	private AlignmentDataModelStats pairedEndData;
	private IntervalTree<Alignments> pairedEndTree;
	private Comparator<RefSeqGene> sortSingleExonsByDistanceFromLastExonComparatorPlusOrientation;
	private Comparator<RefSeqGene> sortSingleExonsByDistanceFromLastExonComparatorMinusOrientation;
	private IntervalTree<RefSeqGene> singleExonsTree;
	private EmpiricalDistribution pairedInsertDistribution;
	private double minProbabilityInPairedInsertDistributionToUseWhenConnectingToSingleExons;
	private double maxProbabilityInPairedInsertDistributionToUseWhenConnectingToSingleExons;
	private GenesInExpressionBins<PairedEndSignificanceInExpressionBin> genesInBins;
	private double minPoissonForConnectingExons;
	
	public ConnectLastExonToSingleExonsWithPairedEnds(String chr, Collection<RefSeqGene> genes, ContinuousDataAlignmentModel data, int maxDistanceForSingleExonsAfterLastExon, double minProbabilityInPairedInsertDistributionToUseWhenConnectingToSingleExons, double maxProbabilityInPairedInsertDistributionToUseWhenConnectingToSingleExons, int numOfBinsForGenes, double minPoissonForConnectingExons, double middlePercentOfExpressionForExtractingMiddleGenes) throws IOException {
		this.chr = chr;
		this.genes = genes;
		this.data = data;
		this.pairedEndData = data.getPairedData();
		this.maxDistanceForSingleExonsAfterLastExon = maxDistanceForSingleExonsAfterLastExon;
		this.minProbabilityInPairedInsertDistributionToUseWhenConnectingToSingleExons = minProbabilityInPairedInsertDistributionToUseWhenConnectingToSingleExons;
		this.maxProbabilityInPairedInsertDistributionToUseWhenConnectingToSingleExons = maxProbabilityInPairedInsertDistributionToUseWhenConnectingToSingleExons;
		this.minPoissonForConnectingExons = minPoissonForConnectingExons;
		
		initiateSortSingleExonsByDistanceFromLastExonComparators();
		System.out.println("Generating single exons tree");
		generateSingleExonsTree();
		System.out.println("Generating paired end tree");
		generatePairedEndTree();		
		System.out.println("Generating paired insert distribution");
		generatePairedInsertDistribution();
		System.out.println("Placing genes in bins based on expression");
		generateGenesInExpressionBins(numOfBinsForGenes, middlePercentOfExpressionForExtractingMiddleGenes);
	}

	// chr1:60389312-60389784 -- too many paired ends to connect
	
	public void initiateSortSingleExonsByDistanceFromLastExonComparators() {
		this.sortSingleExonsByDistanceFromLastExonComparatorPlusOrientation = new Comparator<RefSeqGene>() {
			public int compare(RefSeqGene gene1, RefSeqGene gene2) {
				int d = gene1.getEnd() - gene2.getEnd();
				if (d == 0)
					d = gene1.getOrientedLastExon().length() - gene2.getOrientedLastExon().length();
				return d;
			}
		};
		
		this.sortSingleExonsByDistanceFromLastExonComparatorMinusOrientation = new Comparator<RefSeqGene>() {
			public int compare(RefSeqGene gene1, RefSeqGene gene2) {
				int d = gene2.getStart() - gene1.getStart();
				if (d == 0)
					d = gene1.getOrientedLastExon().length() - gene2.getOrientedLastExon().length();
				return d;
			}
		};
	}
	
	private void generateSingleExonsTree() {
		this.singleExonsTree = new IntervalTree<RefSeqGene>();
		for (RefSeqGene gene : this.genes) {
			if (!gene.getChr().equals(this.chr))
				continue;
			
			if (gene.getNumExons() == 1)
				this.singleExonsTree.put(gene.getStart(), gene.getEnd(), gene);
		}
	}
	
	private void generatePairedInsertDistribution() {
		Collection<RefSeqGene> simplePaths = EstimatePairedEndDistribution.findSimplePathsRefSeqGene(this.genes);
		EmpiricalDistribution ed = EstimatePairedEndDistribution.estimatePairedInsertDistributionRefSeqGene(simplePaths, this.pairedEndTree);
		this.pairedInsertDistribution = ed;
	}
	
	private void generateGenesInExpressionBins(int numOfBinsForGenes, double middlePercent) {
		this.genesInBins = new GenesInExpressionBins<PairedEndSignificanceInExpressionBin>(genes, numOfBinsForGenes, PairedEndSignificanceInExpressionBin.class, this.data, this.chr, middlePercent);
	}
	
	/**
	private void generatePairedEndTree() {
		this.pairedEndTree = new IntervalTree<Alignments>();
		CloseableIterator<Alignment> readIt = pairedEndData.getReadIterator(new Alignments(chr, 0, pairedEndData.getChromosomeLengths().get(chr)));
		while (readIt.hasNext()) {
			Alignment a = readIt.next();
			
			if (!a.getChr().equals(chr))
				continue;
			
			Alignments a1 = new Alignments(a.getChr(), a.getAlignmentStart(), a.getAlignmentEnd());
			this.pairedEndTree.put(a.getAlignmentStart(), a.getAlignmentEnd(), a1);
		}
	}
	 * @throws IOException 
	*/
	
	private void generatePairedEndTree() throws IOException {
		this.pairedEndData.getData().clearFullIntervalTreeAsAlignmentsCached(this.chr);
		
		boolean isStranded = this.pairedEndData.getData().isStranded();
		boolean isNegativeStranded = this.pairedEndData.getData().isNegativeStranded();
		this.pairedEndData.getData().unsetStranded();
		
	 	this.pairedEndTree = this.pairedEndData.getData().getFullIntervalTreeAsAlignments(this.chr);

	 	if (isStranded) {
	 		if (isNegativeStranded)
	 			this.pairedEndData.getData().setNegativeStranded();
	 		else
	 			this.pairedEndData.getData().setPositiveStranded();
	 	}
	}
	
	public void writePairedInsertDistributionToFile(String saveToFile) throws IOException {
		this.pairedInsertDistribution.write(saveToFile + ".pairedInsertSizeDistribution.txt");
	}
	
	private boolean pairedEndIsValidForSingleExon(Alignments pairedEnd, RefSeqGene singleExon, RefSeqGene gene, Alignments orientedLastExon) {
		int pairedEndSingleExonEnd;
		int pairedEndLastExonEnd;
		if (gene.getOrientation().equals("+")) {
			pairedEndSingleExonEnd = pairedEnd.getEnd() - 1;
			pairedEndLastExonEnd = pairedEnd.getStart();
		} else {
			pairedEndSingleExonEnd = pairedEnd.getStart();
			pairedEndLastExonEnd = pairedEnd.getEnd() - 1;
		}

		if (!(pairedEndSingleExonEnd >= singleExon.getStart() && pairedEndSingleExonEnd < singleExon.getEnd()))
			return false;

		Alignments a = new Alignments(pairedEnd.getChr(), pairedEndLastExonEnd, pairedEndLastExonEnd + 1);
		if (!gene.hasExon(a))
			return false;

		return true;
	}
	
	private boolean singleExonIsValidForGene(RefSeqGene singleExon, RefSeqGene gene) {
		if (!singleExon.getOrientation().equals("*") && !singleExon.getOrientation().equals(gene.getOrientation()))
			return false;
		
		int geneEnd = gene.getOrientedEnd();
		if (gene.getOrientation().equals("+")) {
			if (singleExon.getStart() < geneEnd)
				return false;
		} else {
			if (singleExon.getEnd() > geneEnd)
				return false;
		}
		
		return true;
	}
		
	public Collection<RefSeqGene> connectLastExonToSingleExonsWithPairedEnds() throws IOException{	
		Set<RefSeqGene> singleExonsLinkedTo = new TreeSet<RefSeqGene>();
		Map<RefSeqGene, RefSeqGene> genesWithUpdatedEnds = new TreeMap<RefSeqGene, RefSeqGene>();

		for (RefSeqGene gene : this.genes) {
			if (!gene.getChr().equals(this.chr))
				continue;
			if (gene.getOrientation().equals("*"))
				continue;
			
			Alignments orientedLastExon = gene.getOrientedLastExon();
			Alignments spaceToSearchForSingleExons;
			if (gene.getOrientation().equals("+"))
				spaceToSearchForSingleExons = new Alignments(gene.getChr(), orientedLastExon.getEnd(), orientedLastExon.getEnd() + this.maxDistanceForSingleExonsAfterLastExon);
			else
				spaceToSearchForSingleExons = new Alignments(gene.getChr(), orientedLastExon.getStart() - this.maxDistanceForSingleExonsAfterLastExon, orientedLastExon.getStart());

			Set<RefSeqGene> singleExonsFollowingGene = new HashSet<RefSeqGene>();
			Map<RefSeqGene, Set<Alignments>> singleExonsFollowingGeneWithPairedEndAlignment = new HashMap<RefSeqGene, Set<Alignments>>();

			Iterator<Node<RefSeqGene>> singleExonsFollowing = this.singleExonsTree.overlappers(spaceToSearchForSingleExons.getStart(), spaceToSearchForSingleExons.getEnd());
			while (singleExonsFollowing.hasNext()) {
				Node<RefSeqGene> singleExonFollowingNode = singleExonsFollowing.next();
				RefSeqGene singleExon = singleExonFollowingNode.getValue();
				if (!singleExonIsValidForGene(singleExon, gene))
					continue;
				singleExonsFollowingGene.add(singleExon);

				Alignments spaceToSearchForPairedEnds;
				if (gene.getOrientation().equals("+"))
					spaceToSearchForPairedEnds = new Alignments(gene.getChr(), orientedLastExon.getStart(), singleExon.getEnd());
				else
					spaceToSearchForPairedEnds = new Alignments(gene.getChr(), singleExon.getStart(), orientedLastExon.getEnd());
				
				Set<Alignments> pairedEndsForThisSingleExon = new HashSet<Alignments>();
				Iterator<Node<Alignments>> pairedEndsNearby = this.pairedEndTree.overlappers(spaceToSearchForPairedEnds.getStart(), spaceToSearchForPairedEnds.getEnd());
				while (pairedEndsNearby.hasNext()) {
					Node<Alignments> pairedEndNode = pairedEndsNearby.next();
					Alignments pairedEnd = pairedEndNode.getValue();
					if (pairedEndIsValidForSingleExon(pairedEnd, singleExon, gene, orientedLastExon)) {
						pairedEndsForThisSingleExon.add(pairedEnd);
					}
				}
				
				if (pairedEndsForThisSingleExon.size() > 0){
					singleExonsFollowingGeneWithPairedEndAlignment.put(singleExon, pairedEndsForThisSingleExon);
				}
			}

			if (singleExonsFollowingGeneWithPairedEndAlignment.size() == 0)
				continue;

			RefSeqGene[] singleExonsFollowingGeneWithPairedEndSortedByDistanceFromLastExon = (RefSeqGene[]) singleExonsFollowingGeneWithPairedEndAlignment.keySet().toArray(new RefSeqGene[singleExonsFollowingGeneWithPairedEndAlignment.keySet().size()]);
			if (gene.getOrientation().equals("+"))
				Arrays.sort(singleExonsFollowingGeneWithPairedEndSortedByDistanceFromLastExon, sortSingleExonsByDistanceFromLastExonComparatorPlusOrientation);
			else
				Arrays.sort(singleExonsFollowingGeneWithPairedEndSortedByDistanceFromLastExon, sortSingleExonsByDistanceFromLastExonComparatorMinusOrientation);
			
			RefSeqGene singleExonToConnectLastExonTo = null;
			RefSeqGene geneWithUpdatedEnd = null;
			int numOfPairedEndsSupportingConnectionFromLastExon = -1;
			for (int i = singleExonsFollowingGeneWithPairedEndSortedByDistanceFromLastExon.length - 1; i >= 0; i--) {
				RefSeqGene singleExonFollowingGeneWithPairedEnd = singleExonsFollowingGeneWithPairedEndSortedByDistanceFromLastExon[i];
				Set<Alignments> pairedEndsForThisSingleExon = singleExonsFollowingGeneWithPairedEndAlignment.get(singleExonFollowingGeneWithPairedEnd);
				int numOfPairedEndsSupportingThisSingleExon = pairedEndsForThisSingleExon.size();
				if (!thereAreEnoughPairedEndsForConnection(pairedEndsForThisSingleExon, gene, singleExonFollowingGeneWithPairedEnd)) {
					continue;
				}
				
				RefSeqGene geneClone = new RefSeqGene(gene.toBED(), false);
				int newEnd;
				if (gene.getOrientation().equals("+"))
					newEnd = singleExonFollowingGeneWithPairedEnd.getEnd();
				else
					newEnd = singleExonFollowingGeneWithPairedEnd.getStart();
				geneClone.updateLastExonWithNewEnd(newEnd);
				
				boolean singleExonShouldBeConnected = false;
				for (Alignments pairedEnd : pairedEndsForThisSingleExon) {
					double insertSize = estimateInsertSize(pairedEnd, geneClone);
					double cumulativeProbability = pairedInsertDistribution.getCummulativeProbability(insertSize);
					if (cumulativeProbability >= this.minProbabilityInPairedInsertDistributionToUseWhenConnectingToSingleExons && cumulativeProbability <= this.maxProbabilityInPairedInsertDistributionToUseWhenConnectingToSingleExons) {
						singleExonShouldBeConnected = true;
						break;
					}
				}

				if (singleExonShouldBeConnected) {
					singleExonToConnectLastExonTo = singleExonFollowingGeneWithPairedEnd;
					geneWithUpdatedEnd = geneClone;
					numOfPairedEndsSupportingConnectionFromLastExon = numOfPairedEndsSupportingThisSingleExon;
					break;
				}
			}
			
			if (singleExonToConnectLastExonTo == null)
				continue;
			
			for (RefSeqGene singleExon : singleExonsFollowingGene) {
				if (gene.getOrientation().equals("+")) {
					if (singleExon.getEnd() <= singleExonToConnectLastExonTo.getEnd())
						singleExonsLinkedTo.add(singleExon);
				} else {
					if (singleExon.getStart() >= singleExonToConnectLastExonTo.getStart())
						singleExonsLinkedTo.add(singleExon);
				}
			}
			
			Set<RefSeqGene> possibleSingleExonsToConnectThisSingleExonTo = new TreeSet<RefSeqGene>();
			this.extendSingleExonConnectedToLastExonToOtherSingleExonsUsingRecursivePairedEnds(singleExonToConnectLastExonTo, gene, possibleSingleExonsToConnectThisSingleExonTo);
			
			int newEnd = -1;
			for (RefSeqGene singleExon : possibleSingleExonsToConnectThisSingleExonTo) {
				if (gene.getOrientation().equals("+")) {
					if (newEnd == -1 || singleExon.getEnd() > newEnd)
						newEnd = singleExon.getEnd();
				} else {
					if (newEnd == -1 || singleExon.getStart() < newEnd)
						newEnd = singleExon.getStart();
				}
				
				singleExonsLinkedTo.add(singleExon);
			}
			
			if (newEnd > -1)
				geneWithUpdatedEnd.updateLastExonWithNewEnd(newEnd);
			
			genesWithUpdatedEnds.put(gene, geneWithUpdatedEnd);
			
			System.out.println("Modifying " + gene.getName() + " ::: to ::: " + geneWithUpdatedEnd.getName());
		}
		
		Collection<RefSeqGene> modifiedGenes = new TreeSet<RefSeqGene>();
		for (RefSeqGene gene : this.genes) {
			if (singleExonsLinkedTo.contains(gene))
				continue;
			
			if (genesWithUpdatedEnds.containsKey(gene)) {
				RefSeqGene modifiedGene = genesWithUpdatedEnds.get(gene);
				modifiedGenes.add(modifiedGene);
				continue;
			}
			
			modifiedGenes.add(gene);
		}
		
		return modifiedGenes;
	}
	
	private boolean thereAreEnoughPairedEndsForConnection(Set<Alignments> pairedEndsSupportingConnection, RefSeqGene connectFrom, RefSeqGene connectTo) throws IOException{
		int numOfPairedEndsAtConnectionSite = pairedEndsSupportingConnection.size();
		
		int connectionSiteStart = Integer.MAX_VALUE;
		int connectionSiteEnd = -Integer.MAX_VALUE;
		for (Alignments pairedEnd : pairedEndsSupportingConnection) {
			connectionSiteStart = Math.min(connectionSiteStart, pairedEnd.getStart());
			connectionSiteEnd = Math.max(connectionSiteEnd, pairedEnd.getEnd());
		}
		
		int lengthOfConnectionSite;
		if (connectFrom.hasExon(new Alignments(this.chr, connectionSiteStart, connectionSiteStart + 1))) {
			int connectFromLength = connectFrom.trimAbsolute(connectionSiteStart, connectFrom.getEnd()).getSize();
			int connectToLength = connectTo.trimAbsolute(connectTo.getStart(), connectionSiteEnd).getSize();
			lengthOfConnectionSite = connectFromLength + connectToLength;
		} else {
			int connectFromLength = connectFrom.trimAbsolute(connectFrom.getStart(), connectionSiteEnd).getSize();
			int connectToLength = connectTo.trimAbsolute(connectionSiteStart, connectTo.getEnd()).getSize();
			lengthOfConnectionSite = connectFromLength + connectToLength;
		}

		double rpkmOfConnectFrom = connectFrom.getRPKM();
		double rpkmOfConnectTo = connectTo.getRPKM();
		
		double rpkmOfBoth = Math.min(rpkmOfConnectFrom, rpkmOfConnectTo); // THIS IS ONE WAY OF TAKING THE RPKM OF BOTH TRANSCRIPTS IN THE CONNECTION; TO USE THE MEAN OF THE TWO RPKM'S, COMMENT OUT THIS LINE AND USE THE BELOW LINE
		//double rpkmOfBoth = (rpkmOfConnectFrom + rpkmOfConnectTo) / 2;
		
		PairedEndSignificanceInExpressionBin pesieb = this.genesInBins.getValueForExpressionUsingNearbyExpressionBinsIfNeeded(rpkmOfBoth);
		double gammaPForExpressionBin = pesieb.getGammaP();
		double expectedNumOfPairedEnds = gammaPForExpressionBin * lengthOfConnectionSite;
		
		if (numOfPairedEndsAtConnectionSite >= expectedNumOfPairedEnds) {
			return true;
		} else {
			double poisson = AlignmentDataModelStats.poisson(numOfPairedEndsAtConnectionSite, expectedNumOfPairedEnds);
			return poisson >= this.minPoissonForConnectingExons;
		}
	}
	
	private void extendSingleExonConnectedToLastExonToOtherSingleExonsUsingRecursivePairedEnds(RefSeqGene singleExonConnectedToLastExon, RefSeqGene gene, Set<RefSeqGene> extendedSingleExons) throws IOException{
		if (extendedSingleExons.contains(singleExonConnectedToLastExon))
			return;

		Map<RefSeqGene, Set<Alignments>> pairedEndsSupportingConnectionToSingleExon = new TreeMap<RefSeqGene, Set<Alignments>>();
		
		Iterator<Node<Alignments>> pairedEndsOverlappingSingleExonConnectedToLastExon = this.pairedEndTree.overlappers(singleExonConnectedToLastExon.getStart(), singleExonConnectedToLastExon.getEnd());
		while (pairedEndsOverlappingSingleExonConnectedToLastExon.hasNext()) {
			Node<Alignments> pairedEndNode = pairedEndsOverlappingSingleExonConnectedToLastExon.next();
			Alignments pairedEnd = pairedEndNode.getValue();

			double insertSize = pairedEnd.getEnd() - pairedEnd.getStart(); // estimating the insert size this way works because we're estimating the insert size of the paired read when singleExonConnectedToLastExon is connected to another single exon -- pairedEnd therefore overlaps only exons
			double cumulativeProbability = pairedInsertDistribution.getCummulativeProbability(insertSize);

			if (cumulativeProbability < this.minProbabilityInPairedInsertDistributionToUseWhenConnectingToSingleExons || cumulativeProbability > this.maxProbabilityInPairedInsertDistributionToUseWhenConnectingToSingleExons)
				continue;

			Alignments spaceToSearchForSingleExons;

			if (gene.getOrientation().equals("+")) {
				if (!(pairedEnd.getStart() >= singleExonConnectedToLastExon.getStart() && pairedEnd.getStart() <= singleExonConnectedToLastExon.getEnd()))
					continue;
				if (!(pairedEnd.getEnd() > singleExonConnectedToLastExon.getEnd()))
					continue;
				
				spaceToSearchForSingleExons = new Alignments(gene.getChr(), pairedEnd.getEnd(), pairedEnd.getEnd() + 1);
			} else {
				if (!(pairedEnd.getEnd() >= singleExonConnectedToLastExon.getStart() && pairedEnd.getEnd() <= singleExonConnectedToLastExon.getEnd()))
					continue;
				if (!(pairedEnd.getStart() < singleExonConnectedToLastExon.getStart()))
					continue;
				
				spaceToSearchForSingleExons = new Alignments(gene.getChr(), pairedEnd.getStart(), pairedEnd.getStart() + 1);
			}

			Iterator<Node<RefSeqGene>> singleExonsAttachedToPairedEnd = this.singleExonsTree.overlappers(spaceToSearchForSingleExons.getStart(), spaceToSearchForSingleExons.getEnd());
			while (singleExonsAttachedToPairedEnd.hasNext()) {
				Node<RefSeqGene> singleExonFollowingNode = singleExonsAttachedToPairedEnd.next();
				RefSeqGene singleExon = singleExonFollowingNode.getValue();

				if (gene.getOrientation().equals("+")) {
					if (!(singleExon.getEnd() > singleExonConnectedToLastExon.getEnd()))
						continue;
				} else {
					if (!(singleExon.getStart() < singleExonConnectedToLastExon.getStart()))
						continue;
				}

				if (!singleExon.getOrientation().equals("*") && !singleExon.getOrientation().equals(gene.getOrientation()))
					continue;

				if (pairedEndsSupportingConnectionToSingleExon.containsKey(singleExon)) {
					pairedEndsSupportingConnectionToSingleExon.get(singleExon).add(pairedEnd);
				} else {
					TreeSet<Alignments> pairedEndsSupportingConnection = new TreeSet<Alignments>();
					pairedEndsSupportingConnection.add(pairedEnd);
					pairedEndsSupportingConnectionToSingleExon.put(singleExon, pairedEndsSupportingConnection);
				}
			}
		}
		
		for (RefSeqGene singleExon : pairedEndsSupportingConnectionToSingleExon.keySet()) {
			Set<Alignments> pairedEndsSupportingConnection = pairedEndsSupportingConnectionToSingleExon.get(singleExon);

			if (!thereAreEnoughPairedEndsForConnection(pairedEndsSupportingConnection, singleExonConnectedToLastExon, singleExon))
				continue;

			extendedSingleExons.add(singleExon);
			
			this.extendSingleExonConnectedToLastExonToOtherSingleExonsUsingRecursivePairedEnds(singleExon, gene, extendedSingleExons);
		}
	}
	
	private static double estimateInsertSize(Alignments pairedEnd, RefSeqGene gene) {
		return gene.trimAbsolute(pairedEnd.getStart(), pairedEnd.getEnd()).getSize();
	}
	
}
