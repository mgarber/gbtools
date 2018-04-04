package broad.pda.seq.segmentation;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.math.MathException;
import org.broad.igv.Globals;
import org.junit.Test;

import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import junit.framework.TestCase;

public class TestPolyASeqAnalyzer extends TestCase{

	@Test
	/*public void testStrandedNess() throws IOException {
		
		Globals.setHeadless(true);
		System.out.println("Using Version R4.4");
		
		AlignmentDataModel alignments=new GenericAlignmentDataModel("temp.bam", "dm3.sizes", false, 5,true,false);
		
		System.out.println("AlignmentDataModel loaded, initializing model stats");
		AlignmentDataModelStats alignmentData = new AlignmentDataModelStats(alignments);
		ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(alignmentData);
		
		Alignments al = new Alignments("chrX","2599157","2600737");
		double XX = data.count(al);
		System.out.println("All: "+XX);
		
		alignments.setPositiveStranded();
		XX = data.count(al);
		System.out.println("Positive strand: "+XX);
		
		
		alignments.setSecondRead();
		XX = data.count(al);
		System.out.println("Positive second read strand: "+XX);
		BEDFileParser ann = new BEDFileParser("dm3.refseq.bed");
		
		RefSeqGene rf = new RefSeqGene(al);
		//Get the PairedEndDistribution for the entire alignment model
		PairedEndDistribution ped = new PairedEndDistribution(alignments,rf,true);
		ped.ensureNonZeroCounts();
		//Get the iterator only over the refseq gene. Gets all reads. Positive/Negative strand and both mates in that region
		//Will Have to check for strandedness and first/second read in the ietrator loop
		CloseableIterator<Alignment> readIt = alignments.getReadIterator(rf.getAlignment());
		while(readIt.hasNext()) {
			SamAlignment aln = (SamAlignment) readIt.next();
			ReadMate mate = aln.getMate();
			if(aln.isNegativeStrand()){
				System.out.println("Negative Strand");
			}
			else{
			if(aln.isSecondOfPair()){
				System.out.println("read is second of pair: Mate is first");
				//System.out.println(mate.toString());
			}
			else{
				System.out.println("Uhhhh OOOHhhhh");
			//	System.out.println(mate.toString());
			}
			}
		}
		readIt.close();
//		XX =alignments.getTotalNumberOfStrandedReads();
		alignments.setNegativeStranded();
		XX = data.count(al);
		System.out.println("Negative strand: "+XX);
		
			
	}*/
	
	public void testGenericStrandedNess() throws IOException, IllegalArgumentException, MathException{
		
		Globals.setHeadless(true);
		System.out.println("Using Version R4.4");
		System.out.println("Processing " + "temp2.bam");
		int [] defaultWindows = {0};
		AlignmentDataModelStats dataModelStats = null;
		GenericAlignmentDataModel dataModel = new GenericAlignmentDataModel("temp2.bam", null, false, 5, true,true);
		//dataModel.setPositiveStranded();
		Alignments al = new Alignments("chr2L","3479630","3479813");

		ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(dataModel);
		System.out.println(data.getCount(al));
		
		BEDFileParser segmentation = data.segmentChromosome(0.1, defaultWindows, "chr2L", false, false, "chr2L.fa",0.05, 0, Integer.MAX_VALUE, 0);
		BufferedWriter bw = new BufferedWriter(new FileWriter("temp2.bed"));
		segmentation.writeFullBed(bw);
		bw.close();
		
		//dataModel.setNegativeStranded();
/*		GenericAlignmentDataModel dataModelNeg = new GenericAlignmentDataModel("temp2.bam", null, false, 5, true,true,"-");
		ContinuousDataAlignmentModel dataNeg = new ContinuousDataAlignmentModel(dataModel);
		System.out.println(dataNeg.getCount(al));
		
		segmentation = dataNeg.segmentChromosome(0.1, defaultWindows, "chr2L", false, false, "chr2L.fa",0.05, 0, Integer.MAX_VALUE, 0.1);
		bw = new BufferedWriter(new FileWriter("temp2.negative.bed"));
		segmentation.writeFullBed(bw);
		bw.close();*/
/*		CloseableIterator<Alignment> readIt = dataModel.getReadIterator(al);
		while(readIt.hasNext()) {
			SamAlignment aln = (SamAlignment) readIt.next();
			if(aln.isNegativeStrand()){
				System.out.println("Negative Strand");
			}
			else{
				System.out.println("Postive Strand");
			}
		}*/
	}

}
