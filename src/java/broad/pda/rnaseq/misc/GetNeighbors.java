package broad.pda.rnaseq.misc;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;


public class GetNeighbors {

	Map<Alignments, Alignments[]> neighbors;
	
	public GetNeighbors(Set<Alignments> alignments, Set<Alignments> genes){
		Map<String, IntervalTree> trees=makeIntervalTree(genes);
		neighbors=getNeighbors(alignments, trees);
	}
	
	
	private Map<Alignments, Alignments[]> getNeighbors(Set<Alignments> alignments, Map<String, IntervalTree> trees){
		Map rtrn=new TreeMap();
		
		for(Alignments align: alignments){
			Alignments left=this.getLeftNeighbor(align, trees);
			Alignments right=this.getRightNeighbor(align, trees);
			Alignments[] array={left, right};
			int leftDistance=distance(align, left);
			int rightDistance=distance(align, right);
			//System.out.println(align+"\t"+leftDistance+"\t"+rightDistance);
			rtrn.put(align, array);
		}
		
		return rtrn;
	}
	
	public Map<Alignments, Alignments[]> getNeighbors() {
		
		return this.neighbors;
	}
	
	private int distance(Alignments align, Alignments left) {
		return Math.max(align.getStart(), left.getStart())-Math.min(align.getEnd(), left.getEnd());
	}


	public static Map<String, IntervalTree> makeIntervalTree(Set<Alignments> alignments){
		Map<String, IntervalTree> rtrn=new TreeMap();
		
		for(Alignments align: alignments){
			IntervalTree tree=new IntervalTree();
			String chr=align.getChr();
			if(rtrn.containsKey(chr)){tree=rtrn.get(chr);}
			tree.put(align.getStart(), align.getEnd(), align);
			rtrn.put(align.getChr(), tree);
		}
		
		return rtrn;
	}
	
	
	private static Alignments getRightNeighbor(Alignments align, Set<Alignments> genes){
		Map<String, IntervalTree> trees=makeIntervalTree(genes);
		
		IntervalTree tree=trees.get(align.getChr());
		
		Iterator iter=tree.iterator(align.getStart(), align.getEnd());
		
		Alignments rtrn=null;
		boolean overlap=true;
		while(iter.hasNext() && overlap){
			Node<Alignments> node=(Node<Alignments>)iter.next();
			rtrn=node.getValue();
			overlap=rtrn.overlapsAtAll(align);
		}
		return rtrn;
		
	}
	
	
	private Alignments getRightNeighbor(Alignments align, Map<String, IntervalTree> trees){
		IntervalTree tree=trees.get(align.getChr());
		
		Iterator iter=tree.iterator(align.getStart(), align.getEnd());
		
		Alignments rtrn=null;
		boolean overlap=true;
		while(iter.hasNext() && overlap){
			Node<Alignments> node=(Node<Alignments>)iter.next();
			rtrn=node.getValue();
			overlap=rtrn.overlapsAtAll(align);
		}
		return rtrn;
		
	}
	
	private Alignments getLeftNeighbor(Alignments align, Map<String, IntervalTree> trees){
		IntervalTree tree=trees.get(align.getChr());
		
		Iterator iter=tree.reverseIterator(align.getStart(), align.getEnd());
		
		Alignments rtrn=null;
		boolean overlap=true;
		while(iter.hasNext() && overlap){
			Node<Alignments> node=(Node<Alignments>)iter.next();
			rtrn=node.getValue();
			overlap=rtrn.overlapsAtAll(align);
		}
		return rtrn;
	}
	
	
	
	private static Alignments getLeftNeighbor(Alignments align, Set<Alignments> genes){
		Map<String, IntervalTree> trees=makeIntervalTree(genes);
		
		IntervalTree tree=trees.get(align.getChr());
		
		Iterator iter=tree.reverseIterator(align.getStart(), align.getEnd());
		
		Alignments rtrn=null;
		boolean overlap=true;
		while(iter.hasNext() && overlap){
			Node<Alignments> node=(Node<Alignments>)iter.next();
			rtrn=node.getValue();
			overlap=rtrn.overlapsAtAll(align);
		}
		return rtrn;
	}
	
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
		Set<Alignments> alignments=BEDFileParser.loadAlignmentData(new File(args[0])); //genes
		Set<Alignments> regions=BEDFileParser.loadAlignmentData(new File(args[1]));
		
		String save=args[2];
		FileWriter writer=new FileWriter(save);
		
		for(Alignments region: regions){
			Alignments neighbor=getRightNeighbor(region, alignments);
			Alignments leftNeighbor=getLeftNeighbor(region, alignments);
			if(leftNeighbor!=null && neighbor!=null){
			//System.err.println(leftNeighbor.getName()+" "+neighbor.getName());
			writer.write(region+"\t"+leftNeighbor.toUCSC()+"\t"+neighbor.toUCSC()+"\n");
			}
			else{System.err.println(region);}
		}
		writer.close();
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=genes \n args[1]=regions \n args[2]=save file";

	

	
}
