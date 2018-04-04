package broad.pda.seq.graph;

import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.graph.ChromosomeWithBubbles2.WindowIterator;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;
import junit.framework.TestCase;

public class TestWindowIterator extends TestCase {
	public void testSkipGapSpanningWindows() {
		List<LightweightGenomicAnnotation> vertices = new ArrayList<LightweightGenomicAnnotation>();
		List<LightweightGenomicAnnotation> edges = new ArrayList<LightweightGenomicAnnotation>();
		
		vertices.add(new BasicLightweightAnnotation("chr6", 50,51));

		vertices.add(new BasicLightweightAnnotation("chr6", 100,101));
		vertices.add(new BasicLightweightAnnotation("chr6", 150,151));

		vertices.add(new BasicLightweightAnnotation("chr6", 200,201));
		vertices.add(new BasicLightweightAnnotation("chr6", 250,251));

		vertices.add(new BasicLightweightAnnotation("chr6", 300,301));
		vertices.add(new BasicLightweightAnnotation("chr6", 350,351));
		
		//All exon splice 
		edges.add(new BasicLightweightAnnotation("chr6", 51,100));
		edges.add(new BasicLightweightAnnotation("chr6", 151,200));
		edges.add(new BasicLightweightAnnotation("chr6", 251,300));

		//alternative splice form (skip econ 200-250)
		edges.add(new BasicLightweightAnnotation("chr6", 151,300));		
		
		ChromosomeWithBubbles2 cwb = new ChromosomeWithBubbles2("chr6", vertices, edges, null, null, 0,0,0,0);

		assertTrue("should be flagged", cwb.isSpanningGap(90, 103));
		
		WindowIterator it = cwb.iterator(200, 90);
		Collection<RefSeqGene> refseqs = it.next();
		
		
		
	}
	public void testGraphWalking() throws IOException {
		ChromosomeWithBubbles2 cwb = loadData();
		int windowSize = 2500;
		int regionStart = 29694786; //29690634
		WindowIterator it = cwb.iterator(windowSize, regionStart);
		System.out.println("ITERATOR TEST");
		long startTime = System.nanoTime();
		while(it.hasNext()) {
 			List<RefSeqGene> paths = new ArrayList<RefSeqGene>(it.next());
 			int start = paths.get(0).getStart();
 			List<RefSeqGene> unCachedPaths = new ArrayList<RefSeqGene>(cwb.getPaths(start, start+windowSize));
 			
 			assertEquals("Cached and uncached returned different number of paths",unCachedPaths.size(), paths.size());
 			
 			for(int i = 0; i < paths.size(); i++) {
 				RefSeqGene cached = paths.get(i);
 				RefSeqGene unCached = unCachedPaths.get(i);
 				
 				assertEquals("Path "+i+ " start did not agree for window " + start,unCached.getStart(), cached.getStart());
 				assertEquals("Path "+i+ " end did not agree for window " + start,unCached.getEnd(), cached.getEnd());
 				Alignments [] cachedExons = cached.getExons();
 				Alignments [] unCachedExons = unCached.getExons();
 				assertEquals("Path "+i+ " exon # was different for window " + start,unCachedExons.length, cachedExons.length);
 				for(int j = 0; j < unCachedExons.length; j++) {
 					assertTrue("Exon "+ j + " differ for Path "+i+ " for window " + start + ", cached exon: " + cachedExons[j].toUCSC() + " uncached exon " + unCachedExons[j].toUCSC(), unCachedExons[j].equals(cachedExons[j]));
 				}
 				
 			}
 			//double start = System.nanoTime();
 			for(RefSeqGene p : paths) {
 				System.out.println(p.toBED());
 			}
 			//System.err.println("Print time: " + (System.nanoTime() - start));
 			if(paths.iterator().next().getStart() > 29709380) {
 				break;
 			}
		}
		System.out.println("ITERATOR TEST END total time: " + (System.nanoTime() - startTime));	
		
		System.out.println("total windows " + it.getNumWindowsDone() + " data access " + it.getDataAccess() + " trivial windows " + it.getWinWithNoJumps());

	}
	
 
	
	private ChromosomeWithBubbles2 loadData() throws IOException {
		URL testAlnURL = getClass().getResource("chr22.alignments.sam");
		URL sizeFileURL = getClass().getResource("sizes");
		AlignmentDataModel data = new GenericAlignmentDataModel(testAlnURL.getPath(), sizeFileURL.getPath());
		ChromosomeWithBubbles2 cwb = new ChromosomeWithBubbles2("chr22");
		cwb.loadBubbles(data, null, 50, true);
		return cwb;
	}
}
