package broad.pda.seq.pairedend;

import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.AlignmentDataModelStats;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.utils.GenesInExpressionBinsValue;

public class PairedEndSignificanceInExpressionBin implements GenesInExpressionBinsValue {

	private Collection<RefSeqGene> genesInBin;
	private ContinuousDataAlignmentModel data;
	private AlignmentDataModelStats pairedEndData;
	private String chr;
	
	private int totalNumOfPairedEndsAtAllConnectionSites;
	private int totalLengthOfAllConnectionSites;
	private double totalAcrossAllConnectionSitesOfNumOfPairedEndsAtConnectionSiteDividedByLengthOfConnectionSite;
	private int numOfConnectionSites;
	
	private double gammaP_1;
	private double gammaP_2;
	
	private boolean couldCalculateGammaP;
	
	public PairedEndSignificanceInExpressionBin() { }

	public void computeValue(Collection<RefSeqGene> genes, ContinuousDataAlignmentModel data, String chr) throws IOException {
		this.couldCalculateGammaP = false;
		
		this.genesInBin = genes;
		this.data = data;
		this.pairedEndData = data.getPairedData();
		this.chr = chr;
		
		this.totalNumOfPairedEndsAtAllConnectionSites = 0;
		this.totalLengthOfAllConnectionSites = 0;
		
		this.totalAcrossAllConnectionSitesOfNumOfPairedEndsAtConnectionSiteDividedByLengthOfConnectionSite = 0;
		this.numOfConnectionSites = 0;
		
		IntervalTree<Alignments> pairedEndTree = this.pairedEndData.getData().getFullIntervalTreeAsAlignments(this.chr);

		TreeMap<Alignments, Set<Alignments>> connectionSitesAlreadyHandled = new TreeMap<Alignments, Set<Alignments>>();
		for (RefSeqGene gene : this.genesInBin) {
			Object[] sortedUniqueExons = gene.getSortedAndUniqueExons().toArray();
			Alignments[] introns = gene.getIntronsBlocks();
			for (int i = 0; i < introns.length; i++) {
				Alignments intron = introns[i];
				
				Alignments prevExon = (Alignments)sortedUniqueExons[i];
				Alignments nextExon = (Alignments)sortedUniqueExons[i + 1];
				
				Collection<Alignments> validPairedEndsForConnectionSite = new TreeSet<Alignments>();

				Iterator<Node<Alignments>> pairedEndsOverlapping = pairedEndTree.overlappers(intron.getStart(), intron.getEnd());				
				while (pairedEndsOverlapping.hasNext()) {
					Node<Alignments> pairedEndNode = pairedEndsOverlapping.next();
					Alignments pairedEnd = pairedEndNode.getValue();

					if ((pairedEnd.getStart() >= prevExon.getStart() && pairedEnd.getStart() < prevExon.getEnd()) && (pairedEnd.getEnd() > nextExon.getStart() && pairedEnd.getEnd() <= nextExon.getEnd()))
						validPairedEndsForConnectionSite.add(pairedEnd);
				}

				if (validPairedEndsForConnectionSite.size() == 0)
					continue;

				int connectionSiteStart = Integer.MAX_VALUE;
				int connectionSiteEnd = -Integer.MAX_VALUE;
				for (Alignments validPairedEnd : validPairedEndsForConnectionSite) {
					connectionSiteStart = Math.min(connectionSiteStart, validPairedEnd.getStart());
					connectionSiteEnd = Math.max(connectionSiteEnd, validPairedEnd.getEnd());
				}
				Alignments connectionSite = new Alignments(this.chr, connectionSiteStart, connectionSiteEnd);
				
				Set<Alignments> intronsHandledForConnectionSite = connectionSitesAlreadyHandled.get(connectionSite);
				if (intronsHandledForConnectionSite != null && intronsHandledForConnectionSite.contains(intron))
					continue;
				
				int numOfPairedEndsAtConnectionSite = validPairedEndsForConnectionSite.size();
				int lengthOfConnectionSite = gene.trimAbsolute(connectionSiteStart, connectionSiteEnd).getSize();

				this.totalNumOfPairedEndsAtAllConnectionSites += numOfPairedEndsAtConnectionSite;
				this.totalLengthOfAllConnectionSites += lengthOfConnectionSite;
				
				double numOfPairedEndsDividedByLength = (double)numOfPairedEndsAtConnectionSite / (double)lengthOfConnectionSite;
				this.totalAcrossAllConnectionSitesOfNumOfPairedEndsAtConnectionSiteDividedByLengthOfConnectionSite += numOfPairedEndsDividedByLength;
				this.numOfConnectionSites++;
			}
		}
		
		if (this.totalLengthOfAllConnectionSites > 0 && this.numOfConnectionSites > 0) {
			this.couldCalculateGammaP = true;	
			
			this.gammaP_1 = (double)totalNumOfPairedEndsAtAllConnectionSites / (double)totalLengthOfAllConnectionSites;
			this.gammaP_2 = totalAcrossAllConnectionSitesOfNumOfPairedEndsAtConnectionSiteDividedByLengthOfConnectionSite / (double)numOfConnectionSites;
		}
	}
	
	public double getGammaP_1() {
		return this.gammaP_1;
	}
	
	public double getGammaP_2() {
		return this.gammaP_2;
	}
	
	public double getGammaP() {
		return getGammaP_1(); // RIGHT NOW THIS USES THE FIRST METHOD TO CALCULATE GAMMA; CHANGE THIS TO getGammaP_2() TO USE THE SECOND METHOD
	}
	
	public boolean couldComputeValue() {
		return this.couldCalculateGammaP;
	}
	
}
