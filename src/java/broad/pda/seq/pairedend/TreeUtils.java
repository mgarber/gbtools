package broad.pda.seq.pairedend;

import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.pda.seq.graph.Path;

public class TreeUtils {

	public static Map<String, IntervalTree<Path>> makeIntervalTreeByPath( Collection <Path> paths) {
		Map<String, IntervalTree<Path>> rtrn=new TreeMap();
		
		for(Path align: paths){
			//System.err.println(align.toGene());
			IntervalTree<Path> tree=new IntervalTree();
			String chr=align.getChromosome();
			//System.err.println(chr);
			if(rtrn.containsKey(chr)){tree=rtrn.get(chr);}
			Node<Path> node = tree.find(align.getStart(), align.getEnd()+1);
			if (node != null)
				node.incrementCount();
			else
				tree.put(align.getStart(), align.getEnd()+1, align);
			rtrn.put(chr, tree);
		}
		
		return rtrn;
	}

}
