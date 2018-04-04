package broad.pda.gene;

public class TestRefSeqPositionMapping  extends junit.framework.TestCase{
	public void testPlusMapping() throws Exception {
		int [] exonStarts = {1000000, 1001250,1002000, 1005000};
		int [] exonEnds = {1000100, 1001350,1002200, 1005200};
		RefSeqGene g = new RefSeqGene("chr1", 1000000, 1005200, "test_gene", 0, "+", exonStarts, exonEnds);
		
		int testPositionShort = 900000;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1000001;
		assertEquals(1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1000050;
		assertEquals(50, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1001250;
		assertEquals(100, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1001249;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1001350;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1001351;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1002000;
		assertEquals(100+100, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1002100;
		assertEquals(100+100+100, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1002200;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1005000;
		assertEquals(400, g.genomicToTranscriptPosition(testPositionShort));		

		testPositionShort = 1005199;
		assertEquals(400+199, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1005200;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1000000;
		assertEquals(0, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1005200;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
	}
	
	public void testMinusMapping() throws Exception {
		int [] exonStarts = {1000000, 1001250,1002000, 1005000};
		int [] exonEnds = {1000100, 1001350,1002200, 1005200};
		RefSeqGene g = new RefSeqGene("chr1", 1000000, 1005200, "test_gene", 0, "-", exonStarts, exonEnds);
		
		int testPositionShort = 900000;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1000001;
		assertEquals(200+200+100+98, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1000050;
		assertEquals(200+200+100+50-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1001250;
		assertEquals(200+200+100-1, g.genomicToTranscriptPosition(testPositionShort));
		
		
		testPositionShort = 1001249;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1001350;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1001351;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1002000;
		assertEquals(200+200-1, g.genomicToTranscriptPosition(testPositionShort));

		testPositionShort = 1002001;
		assertEquals(200+199-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1002100;
		assertEquals(200+100-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1002200;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1005001;
		assertEquals(199-1, g.genomicToTranscriptPosition(testPositionShort));		

		testPositionShort = 1005199;
		assertEquals(0, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1000000;
		assertEquals(200+200+100 + 100-1, g.genomicToTranscriptPosition(testPositionShort));
		
	}
	public void testPositivePositionMapping() throws Exception {
		int [] exonStarts = {1000000, 1001250,1002000, 1005000};
		int [] exonEnds = {1000100, 1001350,1002200, 1005200};
		RefSeqGene g = new RefSeqGene("chr1", 1000000, 1005200, "test_gene", 0, "+", exonStarts, exonEnds);
		
		assertEquals(1000000, g.transcriptToGenomicPosition(0));
		
		int testPositionShort = 1000001;
		assertEquals(1, g.genomicToTranscriptPosition(testPositionShort));
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));
		
		testPositionShort = 1000050;
		assertEquals(50, g.genomicToTranscriptPosition(testPositionShort));
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));


		testPositionShort = 1002001;
		assertEquals(100+100+1, g.genomicToTranscriptPosition(testPositionShort));
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));
		
		testPositionShort = 1002100;
		assertEquals(100+100+100, g.genomicToTranscriptPosition(testPositionShort));
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));
		
		testPositionShort = 1005001;
		assertEquals(100+100+200+1, g.genomicToTranscriptPosition(testPositionShort));		
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));
		
		testPositionShort = 1005199;
		assertEquals(400+199, g.genomicToTranscriptPosition(testPositionShort));
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));
		
		
	}	
	
	public void testCHD7CoordinateMapping() throws Exception {
		int [] exonStarts = {0,62494,102235,116221,121623,122763,129452,137622,141242,143025,143259,143738,145075,149898,151557,157308,158052,158903,159311,162879,163082,166099,166485,169750,170286,171728,172267,172497,173254,173733,174064,175598,177210,177680,182139,183431,183783,186251};
		int [] exonLengths = { 318,1839,431,142,138,66,56,115,84,138,122,244,177,144,256,211,196,168,180,111,206,200,160,90,104,130,73,58,229,209,672,161,228,444,222,141,105,3012};
		
		int [] exonEnds = new int[exonStarts.length];
		int geneStart = 61591323;
		int totalLength = 0;
		for (int i = 0; i < exonStarts.length; i++) {
			exonEnds[i] = exonStarts[i] + exonLengths[i];
			totalLength += exonLengths[i];
			System.out.println("Exon " + i + " at " + (geneStart + exonStarts[i]) + " takes total transcript size to "+ totalLength);
		}

		RefSeqGene g = new RefSeqGene("chr8", 61591323, 61780586, "CHD7", 0, "+", exonStarts, exonEnds);
		
		int snpCoord = 61735063;
		
		System.out.println("Premature stop codon transcript position: " +g.genomicToTranscriptPosition(snpCoord) );
	}

	
	public void testNegativePositionMapping() throws Exception {
		int [] exonStarts = {1000000, 1001250,1002000, 1005000};
		int [] exonEnds = {1000100, 1001350,1002200, 1005200};
		RefSeqGene g = new RefSeqGene("chr1", 1000000, 1005200, "test_gene", 0, "-", exonStarts, exonEnds);
		
		assertEquals(1005200-1, g.transcriptToGenomicPosition(0));
		assertEquals(0, g.genomicToTranscriptPosition(1005200-1));
		
		assertEquals(200, g.genomicToTranscriptPosition(1002199));
		assertEquals(1002199, g.transcriptToGenomicPosition(200));

		int testPositionShort = 1000001;
		assertEquals(200+200+100+99-1, g.genomicToTranscriptPosition(testPositionShort));
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));
		
		testPositionShort = 1000050;
		assertEquals(200+200+100+50-1, g.genomicToTranscriptPosition(testPositionShort));
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));


		testPositionShort = 1002001;
		assertEquals(200+199-1, g.genomicToTranscriptPosition(testPositionShort));
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));
		
		testPositionShort = 1002100;
		assertEquals(200+100-1, g.genomicToTranscriptPosition(testPositionShort));
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));
		
		testPositionShort = 1005001;
		assertEquals(199-1, g.genomicToTranscriptPosition(testPositionShort));		
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));
		
		testPositionShort = 1005199;
		assertEquals(0, g.genomicToTranscriptPosition(testPositionShort));
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));
		
		
	}
	
	
}
