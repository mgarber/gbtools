package broad.core.util;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import broad.pda.datastructures.Alignments;
import junit.framework.TestCase;

public class TestCollapseByIntersection extends TestCase {
	public void testDecolapse() {
		Alignments exon = new Alignments("chr22", 100,200);
		Alignments intron = new Alignments("chr22", 140,150);
		
		List<Alignments> exons = new ArrayList<Alignments>();
		exons.add(exon);
		List<Alignments> introns = new ArrayList<Alignments>();
		introns.add(intron);
		
		Collection<Alignments> decolapse = CollapseByIntersection.DecollapseByIntronLocation(exons, introns);

		assertEquals(2, decolapse.size());
		Iterator<Alignments> it = decolapse.iterator();
		Alignments first = it.next();
		Alignments second = it.next();
		assertEquals(100, first.getStart());
		assertEquals(140, first.getEnd());
		assertEquals(150, second.getStart());
		assertEquals(200, second.getEnd());
		
		
		introns.add(new Alignments("chr22", 130,150));
		decolapse = CollapseByIntersection.DecollapseByIntronLocation(exons, introns);
		assertEquals(3, decolapse.size());
		it = decolapse.iterator();
		first = it.next();
		second = it.next();
		Alignments third = it.next();
		assertEquals(100, first.getStart());
		assertEquals(130, first.getEnd());
		assertEquals(100, second.getStart());
		assertEquals(140, second.getEnd());
		assertEquals(150, third.getStart());
		assertEquals(200, third.getEnd());
		
		introns.add(new Alignments("chr22", 170,180));
		
		decolapse = CollapseByIntersection.DecollapseByIntronLocation(exons, introns);
		assertEquals(4, decolapse.size());
		it = decolapse.iterator();
		first = it.next();
		second = it.next();
		third = it.next();
		Alignments fourth = it.next();
		assertEquals(100, first.getStart());
		assertEquals(130, first.getEnd());
		assertEquals(100, second.getStart());
		assertEquals(140, second.getEnd());
		assertEquals(150, third.getStart());
		assertEquals(170, third.getEnd());
		assertEquals(180, fourth.getStart());
		assertEquals(200, fourth.getEnd());
		
		exons.remove(0);
		exons.add(new Alignments("chr22", 175,210));
		decolapse = CollapseByIntersection.DecollapseByIntronLocation(exons, introns);
		
		assertEquals(1, exons.size());
		it = decolapse.iterator();
		first = it.next();
		
		assertEquals(180, first.getStart());
		assertEquals(210, first.getEnd());
		
	}
}
