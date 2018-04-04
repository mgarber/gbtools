package broad.pda.seq.graph;

import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;
import junit.framework.TestCase;


public class TestBubbleBuilding2 extends TestCase{
	
	
	public void testIntronWithMutltipleNodes3() {
		List<LightweightGenomicAnnotation> vertices = new ArrayList<LightweightGenomicAnnotation>(6);
		vertices.add(new BasicLightweightAnnotation("chr2", 10,20));
		vertices.add(new BasicLightweightAnnotation("chr2", 50,60));
		vertices.add(new BasicLightweightAnnotation("chr2", 100,110));
		vertices.add(new BasicLightweightAnnotation("chr2", 200,210));
		
		List<LightweightGenomicAnnotation> edges = new ArrayList<LightweightGenomicAnnotation>(6);
		edges.add(new BasicLightweightAnnotation("chr2", 20, 200,"+"));
		edges.add(new BasicLightweightAnnotation("chr2", 20, 50,"+"));
		edges.add(new BasicLightweightAnnotation("chr2", 60, 100, "+"));
		edges.add(new BasicLightweightAnnotation("chr2", 110, 200,"+"));
		


		ChromosomeWithBubbles2 cwb = new ChromosomeWithBubbles2("chr2", vertices, edges, null, null, 0,0,0,0);
		Collection<Path> paths = cwb.getPaths(1);
		assertEquals("Wrong number of paths", 2, paths.size());
		
		for(Path p : paths) {
			RefSeqGene r = p.toGene();
			System.err.println(r.toBED());
		}

				
	}
	
	
	public void testAddNonAbuttingEdge() {
		List<LightweightGenomicAnnotation> vertices = new ArrayList<LightweightGenomicAnnotation>(6);
		vertices.add(new BasicLightweightAnnotation("chr2", 10,20));
		vertices.add(new BasicLightweightAnnotation("chr2", 10,21));
		vertices.add(new BasicLightweightAnnotation("chr2", 50,60));
		vertices.add(new BasicLightweightAnnotation("chr2", 50,61));
		vertices.add(new BasicLightweightAnnotation("chr2", 49,60));
		vertices.add(new BasicLightweightAnnotation("chr2", 100,110));
		vertices.add(new BasicLightweightAnnotation("chr2", 99,111));
		
		List<LightweightGenomicAnnotation> edges = new ArrayList<LightweightGenomicAnnotation>(6);
		edges.add(new BasicLightweightAnnotation("chr2", 20, 50));
		edges.add(new BasicLightweightAnnotation("chr2", 20, 100));
		edges.add(new BasicLightweightAnnotation("chr2", 20, 49));
		edges.add(new BasicLightweightAnnotation("chr2", 20, 99));
		edges.add(new BasicLightweightAnnotation("chr2", 20, 50));
		
		ChromosomeWithBubbles2 cwb = new ChromosomeWithBubbles2("chr2", vertices, edges, null, null, 0,0,0,0);
		
		//boolean cooldAdd = cwb.addPairedEndConnections(new BasicLightweightAnnotation("chr2", 15, 55), 2);
		//assertTrue ("Error adding edge", cooldAdd);
		//cooldAdd = cwb.addPairedEndConnections(new BasicLightweightAnnotation("chr2", 55, 105), 2);
		//assertTrue ("Error adding edge", cooldAdd);
		
	}
	
	
	//TODO: Add Jira ticket #
	public void testCloseByJunctions() {
		LightweightGenomicAnnotation lwg = new BasicLightweightAnnotation("chr2", 50, 100);
		LightweightGenomicAnnotation lwg1 = new BasicLightweightAnnotation("chr2", 150, 300);
		LightweightGenomicAnnotation lwg2 = new BasicLightweightAnnotation("chr2", 155, 300);
		ChromosomeWithBubbles2 cwb = new ChromosomeWithBubbles2("chr2");
		
		cwb.addEdgeAddingMissingVertices(lwg, 5);
		cwb.addEdgeAddingMissingVertices(lwg1, 3);
		cwb.addEdgeAddingMissingVertices(lwg2, 8);
		
		//Sanity check, all systems are GO
		Collection<RefSeqGene> paths = cwb.getPaths(40,60);
		
		assertEquals(1, paths.size());
		RefSeqGene path = paths.iterator().next();
		assertEquals(110,path.getEnd());
		
		paths = cwb.getPaths(140,180);
		
		assertEquals(2, paths.size());
		
		
		
	}
	public void testIntronWithMutltipleNodes() {
		List<LightweightGenomicAnnotation> vertices = new ArrayList<LightweightGenomicAnnotation>(6);
		vertices.add(new BasicLightweightAnnotation("chr2", 10,20));
		vertices.add(new BasicLightweightAnnotation("chr2", 10,21));
		vertices.add(new BasicLightweightAnnotation("chr2", 50,60));
		vertices.add(new BasicLightweightAnnotation("chr2", 50,61));
		vertices.add(new BasicLightweightAnnotation("chr2", 49,60));
		vertices.add(new BasicLightweightAnnotation("chr2", 100,110));
		vertices.add(new BasicLightweightAnnotation("chr2", 99,111));
		
		List<LightweightGenomicAnnotation> edges = new ArrayList<LightweightGenomicAnnotation>(6);
		edges.add(new BasicLightweightAnnotation("chr2", 20, 50, "+"));
		edges.add(new BasicLightweightAnnotation("chr2", 20, 100,"+"));
		edges.add(new BasicLightweightAnnotation("chr2", 20, 49,"+"));
		edges.add(new BasicLightweightAnnotation("chr2", 20, 99,"+"));
		edges.add(new BasicLightweightAnnotation("chr2", 20, 50,"+"));


		ChromosomeWithBubbles2 cwb = new ChromosomeWithBubbles2("chr2", vertices, edges, null, null, 0,0,0,0);
		
		List<RefSeqGene> paths = new ArrayList<RefSeqGene>(cwb.getGenePaths(0));
		Collections.sort(paths);
		assertEquals("Wrong # paths", 5, paths.size());
		Iterator<RefSeqGene> pathIt = paths.iterator();
		RefSeqGene path = pathIt.next();

		assertEquals("Wrong start ", 10, path.getStart());
		assertEquals("Wrong end #1", 60, path.getEnd());
		Alignments [] exons = path.getExons();
		assertEquals("Wrong exon #", 2, exons.length);
		assertEquals("Wrong exon start",10, exons[0].getStart());
		assertEquals("Wrong exon start",20, exons[0].getEnd());
		assertEquals("Wrong exon start",49, exons[1].getStart());
		assertEquals("Wrong exon start",60, exons[1].getEnd());
		
		path = pathIt.next();
		assertEquals("Wrong start ", 10, path.getStart());
		assertEquals("Wrong end #2", 60, path.getEnd());
		exons = path.getExons();
		assertEquals("Wrong exon #", 2, exons.length);
		assertEquals("Wrong exon start",10, exons[0].getStart());
		assertEquals("Wrong exon start",20, exons[0].getEnd());
		assertEquals("Wrong exon start",50, exons[1].getStart());
		assertEquals("Wrong exon start",60, exons[1].getEnd());
		
		path = pathIt.next();
		assertEquals("Wrong start ", 10, path.getStart());
		assertEquals("Wrong end #3", 61, path.getEnd());
		exons = path.getExons();
		assertEquals("Wrong exon #", 2, exons.length);
		assertEquals("Wrong exon start",10, exons[0].getStart());
		assertEquals("Wrong exon start",20, exons[0].getEnd());
		assertEquals("Wrong exon start",50, exons[1].getStart());
		assertEquals("Wrong exon start",61, exons[1].getEnd());

		path = pathIt.next();
		assertEquals("Wrong start ", 10, path.getStart());
		assertEquals("Wrong end #4", 110, path.getEnd());
		exons = path.getExons();
		assertEquals("Wrong exon #", 2, exons.length);
		assertEquals("Wrong exon start",10, exons[0].getStart());
		assertEquals("Wrong exon start",20, exons[0].getEnd());
		assertEquals("Wrong exon start",100, exons[1].getStart());
		assertEquals("Wrong exon start",110, exons[1].getEnd());
		
		path = pathIt.next();
		assertEquals("Wrong start ", 10, path.getStart());
		assertEquals("Wrong end #5", 111, path.getEnd());
		exons = path.getExons();
		assertEquals("Wrong exon #", 2, exons.length);
		assertEquals("Wrong exon start",10, exons[0].getStart());
		assertEquals("Wrong exon start",20, exons[0].getEnd());
		assertEquals("Wrong exon start",99, exons[1].getStart());
		assertEquals("Wrong exon start",111, exons[1].getEnd());
		

		
		edges.add(new BasicLightweightAnnotation("chr2", 21, 50,"+"));
		cwb = new ChromosomeWithBubbles2("chr2", vertices, edges, null, null, 0,0,0,0);
		
		paths = new ArrayList<RefSeqGene>(cwb.getGenePaths(0));
		Collections.sort(paths);
		assertEquals("Wrong # paths", 7, paths.size());
		pathIt = paths.iterator();
		
		
		path = pathIt.next();
		assertEquals("Wrong start ", 10, path.getStart());
		assertEquals("Wrong end #1", 60, path.getEnd());
		exons = path.getExons();
		assertEquals("Wrong exon #", 2, exons.length);
		assertEquals("Wrong exon start",10, exons[0].getStart());
		assertEquals("Wrong exon start",20, exons[0].getEnd());
		assertEquals("Wrong exon start",49, exons[1].getStart());
		assertEquals("Wrong exon start",60, exons[1].getEnd());
		
		path = pathIt.next();
		assertEquals("Wrong start ", 10, path.getStart());
		assertEquals("Wrong end #2", 60, path.getEnd());
		exons = path.getExons();
		assertEquals("Wrong exon #", 2, exons.length);
		assertEquals("Wrong exon start",10, exons[0].getStart());
		assertEquals("Wrong exon start",20, exons[0].getEnd());
		assertEquals("Wrong exon start",50, exons[1].getStart());
		assertEquals("Wrong exon start",60, exons[1].getEnd());
		
		path = pathIt.next();
		assertEquals("Wrong start ", 10, path.getStart());
		assertEquals("Wrong end #3", 60, path.getEnd());
		exons = path.getExons();
		assertEquals("Wrong exon #", 2, exons.length);
		assertEquals("Wrong exon start",10, exons[0].getStart());
		assertEquals("Wrong exon start",21, exons[0].getEnd());
		assertEquals("Wrong exon start",50, exons[1].getStart());
		assertEquals("Wrong exon start",60, exons[1].getEnd());
		
		
		path = pathIt.next();
		assertEquals("Wrong start ", 10, path.getStart());
		assertEquals("Wrong end #4", 61, path.getEnd());
		exons = path.getExons();
		assertEquals("Wrong exon #", 2, exons.length);
		assertEquals("Wrong exon start",10, exons[0].getStart());
		assertEquals("Wrong exon start",20, exons[0].getEnd());
		assertEquals("Wrong exon start",50, exons[1].getStart());
		assertEquals("Wrong exon start",61, exons[1].getEnd());
		
		path = pathIt.next();
		assertEquals("Wrong start ", 10, path.getStart());
		assertEquals("Wrong end #5", 61, path.getEnd());
		exons = path.getExons();
		assertEquals("Wrong exon #", 2, exons.length);
		assertEquals("Wrong exon start",10, exons[0].getStart());
		assertEquals("Wrong exon start",21, exons[0].getEnd());
		assertEquals("Wrong exon start",50, exons[1].getStart());
		assertEquals("Wrong exon start",61, exons[1].getEnd());
	
		path = pathIt.next();
		assertEquals("Wrong start ", 10, path.getStart());
		assertEquals("Wrong end #6", 110, path.getEnd());
		exons = path.getExons();
		assertEquals("Wrong exon #", 2, exons.length);
		assertEquals("Wrong exon start",10, exons[0].getStart());
		assertEquals("Wrong exon start",20, exons[0].getEnd());
		assertEquals("Wrong exon start",100, exons[1].getStart());
		assertEquals("Wrong exon start",110, exons[1].getEnd());

		path = pathIt.next();
		assertEquals("Wrong start ", 10, path.getStart());
		assertEquals("Wrong end #7", 111, path.getEnd());
		exons = path.getExons();
		assertEquals("Wrong exon #", 2, exons.length);
		assertEquals("Wrong exon start",10, exons[0].getStart());
		assertEquals("Wrong exon start",20, exons[0].getEnd());
		assertEquals("Wrong exon start",99, exons[1].getStart());
		assertEquals("Wrong exon start",111, exons[1].getEnd());
				
	}
	
	public void testIntronWithMutltipleNodes2() {
		List<LightweightGenomicAnnotation> vertices = new ArrayList<LightweightGenomicAnnotation>(6);
		vertices.add(new BasicLightweightAnnotation("chr2", 10,20));
		vertices.add(new BasicLightweightAnnotation("chr2", 9,20));
		vertices.add(new BasicLightweightAnnotation("chr2", 50,60));
		vertices.add(new BasicLightweightAnnotation("chr2", 50,61));
		vertices.add(new BasicLightweightAnnotation("chr2", 100,110));
		vertices.add(new BasicLightweightAnnotation("chr2", 100,111));
		
		List<LightweightGenomicAnnotation> edges = new ArrayList<LightweightGenomicAnnotation>(6);
		edges.add(new BasicLightweightAnnotation("chr2", 20, 50,"+"));

		ChromosomeWithBubbles2 cwb = new ChromosomeWithBubbles2("chr2", vertices, edges, null, null, 0,0,0,0);
		
		Collection<RefSeqGene> paths = new TreeSet<RefSeqGene>(cwb.getGenePaths(0));
		assertEquals("Wrong # paths", 4, paths.size());
		Iterator<RefSeqGene> pathIt = paths.iterator();
		
		RefSeqGene path = pathIt.next();
		assertEquals("Wrong start ", 9, path.getStart());
		assertEquals("Wrong end #1", 60, path.getEnd());
		Alignments [] exons = path.getExons();
		assertEquals("Wrong exon #", 2, exons.length);
		assertEquals("Wrong exon start",9, exons[0].getStart());
		assertEquals("Wrong exon start",20, exons[0].getEnd());
		assertEquals("Wrong exon start",50, exons[1].getStart());
		assertEquals("Wrong exon start",60, exons[1].getEnd());
		
		path = pathIt.next();
		assertEquals("Wrong start ", 9, path.getStart());
		assertEquals("Wrong end #2", 61, path.getEnd());
		exons = path.getExons();
		assertEquals("Wrong exon #", 2, exons.length);
		assertEquals("Wrong exon start",9, exons[0].getStart());
		assertEquals("Wrong exon start",20, exons[0].getEnd());
		assertEquals("Wrong exon start",50, exons[1].getStart());
		assertEquals("Wrong exon start",61, exons[1].getEnd());
		
		path = pathIt.next();
		assertEquals("Wrong start ", 10, path.getStart());
		assertEquals("Wrong end #3", 60, path.getEnd());
		exons = path.getExons();
		assertEquals("Wrong exon #", 2, exons.length);
		assertEquals("Wrong exon start",10, exons[0].getStart());
		assertEquals("Wrong exon start",20, exons[0].getEnd());
		assertEquals("Wrong exon start",50, exons[1].getStart());
		assertEquals("Wrong exon start",60, exons[1].getEnd());
		
		path = pathIt.next();
		assertEquals("Wrong start ", 10, path.getStart());
		assertEquals("Wrong end #4", 61, path.getEnd());
		exons = path.getExons();
		assertEquals("Wrong exon #", 2, exons.length);
		assertEquals("Wrong exon start",10, exons[0].getStart());
		assertEquals("Wrong exon start",20, exons[0].getEnd());
		assertEquals("Wrong exon start",50, exons[1].getStart());
		assertEquals("Wrong exon start",61, exons[1].getEnd());

		
		edges.add(new BasicLightweightAnnotation("chr2", 60, 100,"+"));
		edges.add(new BasicLightweightAnnotation("chr2", 61, 100,"+"));
		
		//Now add pahts to the third exon.
		cwb = new ChromosomeWithBubbles2("chr2", vertices, edges, null, null, 0,0,0,0);
		paths = cwb.getGenePaths(1);
		assertEquals("Wrong # paths", 8, paths.size());
		pathIt = paths.iterator();
		
		//Finally add a splice form that skips exon 2
		edges.add(new BasicLightweightAnnotation("chr2", 20, 100,"+"));
		cwb = new ChromosomeWithBubbles2("chr2", vertices, edges, null, null, 0,0,0,0);
		paths = cwb.getGenePaths(0);
		assertEquals("Wrong # paths", 12, paths.size());
		pathIt = paths.iterator();
	}
	

	
	public void testSimple() {
		ChromosomeWithBubbles2 cwb = new ChromosomeWithBubbles2("chr2");

		LightweightGenomicAnnotation lwg = new BasicLightweightAnnotation("chr2", 250, 300);
		LightweightGenomicAnnotation lwg2 = new BasicLightweightAnnotation("chr2", 375, 600);
		LightweightGenomicAnnotation lwg3 = new BasicLightweightAnnotation("chr2", 375, 400);
		LightweightGenomicAnnotation lwg4 = new BasicLightweightAnnotation("chr2", 500, 600);
		cwb.addEdgeAddingMissingVertices(lwg2, 6);

		//Window does not intersec any subgraph, should not do nothing
		Collection<RefSeqGene> windows = cwb.getPaths(200, 260);
		assertEquals("Only one window expected", 1, windows.size());
		
		RefSeqGene window = windows.iterator().next();
		assertEquals("Bad size", 60, window.getGappedSize());
		assertEquals("Bad size", 60, window.getGenomicLength()); 
		
		//Only one edge, should yield one path with two exons.
		windows = cwb.getPaths(200, 650);
		assertEquals("Only one window expected", 1, windows.size());
		RefSeqGene path = windows.iterator().next();
		assertEquals("Bad Start",200, path.getStart());
		assertEquals("Bad End", 875, path.getEnd());
		assertEquals("Bad Gapped size",450, path.getGappedSize());
		Alignments [] exons = path.getExons();
		assertEquals("Bad exon #", 2, exons.length);
		assertEquals ("Bad First Exon End", 375, exons[0].getEnd() );
		assertEquals ("Bad second Exon Start", 600, exons[1].getStart());
		
		//Added a second edge within window, should obtain now 3 exons
		cwb.addEdgeAddingMissingVertices(lwg, 5);
		windows = cwb.getPaths(200, 650);
		assertEquals("Only one window expected", 1, windows.size());
		path = windows.iterator().next();
		assertEquals("Bad Start",200, path.getStart());
		assertEquals("Bad End", 925, path.getEnd());
		assertEquals("Bad Gapped size",450, path.getGappedSize());
		exons = path.getExons();
		assertEquals("Bad exon #", 3, exons.length);
		
		assertEquals ("Bad First Exon End", 250, exons[0].getEnd() );
		assertEquals ("Bad second Exon Start", 300, exons[1].getStart());
		assertEquals ("Bad second Exon End", 375, exons[1].getEnd());
		assertEquals ("Bad third Exon Start", 600, exons[2].getStart());
		assertEquals ("Bad third Exon End", 925, exons[2].getEnd());
		
		//Added an alternative splice form, should now get two paths, the first should be
		//the same as last, the second should jave 4 exons.
		cwb.addEdgeAddingMissingVertices(lwg3, 5);
		cwb.addEdgeAddingMissingVertices(lwg4, 5);
		
		windows = cwb.getPaths(200, 650);
		assertEquals("Only one window expected", 2, windows.size());
		
		Iterator<RefSeqGene> windowIt = windows.iterator();
		path = windowIt.next();
		assertEquals("Bad Start",200, path.getStart());
		assertEquals("Bad End", 925, path.getEnd());
		assertEquals("Bad Gapped size",450, path.getGappedSize());
		exons = path.getExons();
		assertEquals("Bad exon #", 3, exons.length);
		
		assertEquals ("Bad First Exon End", 250, exons[0].getEnd() );
		assertEquals ("Bad second Exon Start", 300, exons[1].getStart());
		assertEquals ("Bad second Exon End", 375, exons[1].getEnd());
		assertEquals ("Bad third Exon Start", 600, exons[2].getStart());
		assertEquals ("Bad third Exon End", 925, exons[2].getEnd());
		
		path = windowIt.next();
		assertEquals("Bad Start",200, path.getStart());
		assertEquals("Bad End", 825, path.getEnd());
		assertEquals("Bad Gapped size",450, path.getGappedSize());
		exons = path.getExons();
		assertEquals("Bad exon #", 4, exons.length);
		
		assertEquals ("Bad First Exon End", 250, exons[0].getEnd() );
		assertEquals ("Bad second Exon Start", 300, exons[1].getStart());
		assertEquals ("Bad second Exon End", 375, exons[1].getEnd());
		assertEquals ("Bad third Exon Start", 400, exons[2].getStart());
		assertEquals ("Bad third Exon End", 500, exons[2].getEnd());
		assertEquals ("Bad third Exon Start", 600, exons[3].getStart());
		
	}
	

	
	
	public void testGetPathsEdgeCondition() {
		ChromosomeWithBubbles2 cwm = new ChromosomeWithBubbles2("chr1");
		cwm.addEdgeAddingMissingVertices(new BasicLightweightAnnotation("chr2", 10,200), 5);
		Collection<RefSeqGene> paths = cwm.getPaths(8,108);
		assertEquals(1, paths.size());
		RefSeqGene path = paths.iterator().next();
		assertEquals(8, path.getStart());
		assertEquals(298, path.getEnd());
		Alignments [] exons = path.getExons();
		assertEquals(2, exons.length);
		assertEquals(2, exons[0].length());
		
		paths = cwm.getPaths(9,109);
		path = paths.iterator().next();
		exons = path.getExons();
		assertEquals(9, path.getStart());
		assertEquals(299, path.getEnd());
		assertEquals(2, exons.length);
		assertEquals(1, exons[0].length());
		
		
		paths = cwm.getPaths(10,110);
		path = paths.iterator().next();
		assertEquals(10, path.getStart());
		assertEquals(110, path.getEnd());
		exons = path.getExons();
		assertEquals(1,exons.length);
		
		paths = cwm.getPaths(0,9);
		path = paths.iterator().next();
		exons = path.getExons();
		assertEquals(0, path.getStart());
		assertEquals(9, path.getEnd());
		assertEquals(1, exons.length);
		assertEquals(9, exons[0].length());
		
		paths = cwm.getPaths(0,10);
		path = paths.iterator().next();
		exons = path.getExons();
		assertEquals(0, path.getStart());
		assertEquals(10, path.getEnd());
		assertEquals(1, exons.length);
		assertEquals(10, exons[0].length());
		
		paths = cwm.getPaths(0,11);
		path = paths.iterator().next();
		exons = path.getExons();
		assertEquals(0, path.getStart());
		assertEquals(201, path.getEnd());
		assertEquals(2, exons.length);
		assertEquals(1, exons[1].length());
		
	}


	public void testGetPaths() {
		List<LightweightGenomicAnnotation> vertices = new ArrayList<LightweightGenomicAnnotation>();
		List<LightweightGenomicAnnotation> edges = new ArrayList<LightweightGenomicAnnotation>();
		
		vertices.add(new BasicLightweightAnnotation("chr6", 49,50));

		vertices.add(new BasicLightweightAnnotation("chr6", 100,101));
		vertices.add(new BasicLightweightAnnotation("chr6", 149,150));

		vertices.add(new BasicLightweightAnnotation("chr6", 200,201));
		vertices.add(new BasicLightweightAnnotation("chr6", 249,250));

		vertices.add(new BasicLightweightAnnotation("chr6", 300,301));
		vertices.add(new BasicLightweightAnnotation("chr6", 349,350));
		
		//All exon splice 
		edges.add(new BasicLightweightAnnotation("chr6", 50,100));
		edges.add(new BasicLightweightAnnotation("chr6", 150,200));
		edges.add(new BasicLightweightAnnotation("chr6", 250,300));

		//alternative splice form (skip econ 200-250)
		edges.add(new BasicLightweightAnnotation("chr6", 150,300));		
		
		ChromosomeWithBubbles2 cwb = new ChromosomeWithBubbles2("chr6", vertices, edges, null, null, 0,0,0,0);
		
		Collection<RefSeqGene> paths = cwb.getPaths(0, 300);
		assertTrue("Paths were empty",!paths.isEmpty());
		assertEquals("There should only be 2 paths", 2 , paths.size());

		Iterator<RefSeqGene> pathIt = paths.iterator();
		
		assertTrue("No First path?",pathIt.hasNext());
		RefSeqGene path1 = pathIt.next();
		
		assertEquals("path1 end should be 450", 450, path1.getEnd());
		assertEquals("Should have 4 exons",4, path1.getExons().length);
		int [] path1ExpectedExonsStart = {0,100,200,300};
		int [] path1ExpectedExonsEnd   = {50,150,250,450};
		Alignments [] path1Exons = path1.getExons();
		for(int i = 0; i < 4; i++) {
			assertEquals("Bad Exon " + i + " start ",path1ExpectedExonsStart[i], path1Exons[i].getStart());
			assertEquals("Bad Exon " + i + " end ",path1ExpectedExonsEnd[i], path1Exons[i].getEnd());
		}
		
		RefSeqGene path2 = pathIt.next();
		
		assertEquals("path2 start should be 2", 0, path2.getStart());
		assertEquals("path2 end should be 500", 500, path2.getEnd());
		assertEquals("Should have 3 exons",3, path2.getExons().length);
		int [] path2ExpectedExonsStart = {0,100,300};
		int [] path2ExpectedExonsEnd   = {50,150,500};
		Alignments [] path2Exons = path2.getExons();
		for(int i = 0; i < 3; i++) {
			assertEquals("Bad Exon " + i + " start ",path2ExpectedExonsStart[i], path2Exons[i].getStart());
			assertEquals("Bad Exon " + i + " end ",path2ExpectedExonsEnd[i], path2Exons[i].getEnd());
		}
		
		
	}
	
	public void testGraphWalking() throws IOException {
		ChromosomeWithBubbles2 cwb = loadData();
		
		Collection<RefSeqGene> paths = cwb.getPaths(200000, 200020);
		assertEquals("There should be only one path", 1, paths.size());
		RefSeqGene window = paths.iterator().next();
		assertEquals(200000, window.getStart());
		assertEquals(1, window.getExonSet().size());
		
		paths = cwb.getPaths(17543740, 17543760);
		
		for(RefSeqGene p : paths) {
			//System.out.println(p.toBED());
		}
		
		System.out.println("=====IMPORTANT=====");
		for(int p = 29695634; p < 29705380; p++) {
			paths = cwb.getPaths(p, p+10);
			int i = 0;
			for(RefSeqGene pth : paths) {
				pth.setName("w"+p +"_"+i);
				System.out.println(pth.toBED());
			}
		}
		System.out.println("=====IMPORTANT END=====");
	}
	
	public void testJump() throws IOException {
		ChromosomeWithBubbles2 cwb = loadData();
		
		Collection<RefSeqGene> paths = cwb.getPaths(29695465,29700728);
		//assertEquals("There should be only one path", 1, paths.size());
		RefSeqGene window = paths.iterator().next();
		//assertEquals(200000, window.getStart());
		//assertEquals(1, window.getExonSet().size());
		
		//paths = cwb.getPaths(17543740, 17543760);
		
		for(RefSeqGene p : paths) {
			System.out.println(p.toBED());
		}
		
	}
	
	private ChromosomeWithBubbles2 loadData() throws IOException {
		URL testAlnURL = getClass().getResource("chr22.alignments.sam");
		URL sizeFileURL = getClass().getResource("sizes");
		AlignmentDataModel data = new GenericAlignmentDataModel(testAlnURL.getPath(), sizeFileURL.getPath());
		ChromosomeWithBubbles2 cwb = new ChromosomeWithBubbles2("chr22");
		cwb.loadBubbles(data, null, 20, true);
		return cwb;
	}
	
}
