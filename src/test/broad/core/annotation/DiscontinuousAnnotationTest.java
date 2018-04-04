package broad.core.annotation;

import junit.framework.TestCase;

public class DiscontinuousAnnotationTest extends TestCase {
	public void testAdd1() {
		DiscontinuousAnnotation<GenomicAnnotation> da = new DiscontinuousAnnotation<GenomicAnnotation>("test");
		
		GenomicAnnotation a1 = new BasicGenomicAnnotation("a1", "1", 1000, 1500);

		GenomicAnnotation a2 = new BasicGenomicAnnotation("a2", "1", 2000, 2500);
		
		da.add(a1);
		da.add(a2);
		
		assertEquals("Oh oh, added to non overlapping items but do not have 2 in the list",2, da.size());
		assertEquals("Oh oh the largest item was " + a2.getEnd(),a2.getEnd(),da.getEnd());
		assertEquals("Oh oh the largest item was " + a1.getStart(),a1.getStart(),da.getStart());
		assertEquals("Effective length is not correct", a1.getLength() + a2.getLength(), da.getEffectiveLength());

		GenomicAnnotation a3 = new BasicGenomicAnnotation("a3", "1", 1900, 2200);
		assertTrue("a3 " + a3 + " overlaps!", da.overlaps(a3));
		da.add(a3);
		da.consolidate();

		assertEquals("Added overlapping annotation total size should not have changed", 2, da.size());
		assertEquals("Oh oh the largest item was " + a2.getEnd(),a2.getEnd(),da.getEnd());
		assertEquals("Oh oh the largest item was " + a1.getStart(),a1.getStart(),da.getStart());
		assertEquals("Effective length is not correct", a1.getLength() + 600, da.getEffectiveLength());
		
		GenomicAnnotation a4 = new BasicGenomicAnnotation("a4", "1", 1400, 1950);
		assertEquals("a4 " + a4 + " overlaps!", 2, da.getOverlappingParts(a4).size());
		da.add(a4);
		da.consolidate();
		assertEquals("Added overlapping annotation bridging the existing two, there should be only one", 1, da.size());
		assertEquals("Oh oh the largest item was " + a2.getEnd(),a2.getEnd(),da.getEnd());
		assertEquals("Oh oh the largest item was " + a1.getStart(),a1.getStart(),da.getStart());
		assertEquals("Effective length is not correct", 1500, da.getEffectiveLength());
	}
	
	
}
