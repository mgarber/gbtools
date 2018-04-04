package broad.pda.seq.pairedend;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.broad.igv.sam.Alignment;

import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.seq.segmentation.AlignmentDataModel;
import net.sf.samtools.util.CloseableIterator;


//Try and add paired end junctions and get best possible junctions for each paired end read from seperate alignments of left and right
public class PairedEndMapping {

	static int maxDistance=1000000;
	
	public PairedEndMapping(Map<String, Collection> leftAlignments, Map<String, Collection> rightAlignments, String save)throws IOException{
		System.err.println("V2 started");
		getFilterWritePairs(leftAlignments, rightAlignments, save);
		
		//Map<String, Collection> allPairs=getAllPossiblePairs(leftAlignments, rightAlignments);
		//Map<String, Collection> filteredPairs=filterPairs(allPairs);
		//write(save, filteredPairs);
	}
	
	
	
	
	private void getFilterWritePairs(Map<String, Collection> leftAlignments, Map<String, Collection> rightAlignments, String save)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(String readName: leftAlignments.keySet()){
			Collection<RefSeqGene> left=leftAlignments.get(readName);
			Collection<RefSeqGene> right=rightAlignments.get(readName);
			if(left!=null && right!=null){
			Collection<PairedEndAlignment> pairs=getAllPossiblePairs(left, right, readName);
			for(PairedEndAlignment pair: pairs){if(goodPair(pair)){writer.write(pair.getInsertAlignment()+"\n");}}
			}
		}
		writer.close();
	}
	
	
	/*private void getFilterWritePairsSAM(Map<String, Collection> leftAlignments, Map<String, Collection> rightAlignments, String save)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(String readName: leftAlignments.keySet()){
			Collection<RefSeqGene> left=leftAlignments.get(readName);
			Collection<RefSeqGene> right=rightAlignments.get(readName);
			if(left!=null && right!=null){
			Collection<PairedEndAlignment> pairs=getAllPossiblePairs(left, right, readName);
			for(PairedEndAlignment pair: pairs){
				if(goodPair(pair)){
				String[] records=pair.toSAM();
				for(int i=0; i<records.length; i++){
					writer.write(records[i]+"\n");
				}
				}
			}	
			
			}
		}
		writer.close();
	}*/
	
	private static boolean goodPair(PairedEndAlignment pair){
		if(pair.onSameChromosome() && pair.getDistance()<maxDistance && pair.oppositeStrands()){return true;}
		return false;
	}
	
	private Map filterPairs(Map<String, Collection> allPairs){
		Map rtrn=new TreeMap();
		
		for(String readName: allPairs.keySet()){
			Collection<PairedEndAlignment> pairs=allPairs.get(readName);
			Collection<PairedEndAlignment> filteredPairs=filterPairs(pairs);
			rtrn.put(readName, filteredPairs);
		}
		
		return rtrn;
	}
	
	private Collection<PairedEndAlignment> filterPairs(Collection<PairedEndAlignment> pairs){
		Collection<PairedEndAlignment> rtrn=new HashSet();
		
		for(PairedEndAlignment pair: pairs){
			if(pair.onSameChromosome()){
				rtrn.add(pair);
			}
		}
		
		return rtrn;
	}
	
	private Map getAllPossiblePairs(Map<String, Collection> leftAlignments, Map<String, Collection> rightAlignments){
		Map rtrn=new TreeMap();
		
		for(String readName: leftAlignments.keySet()){
			Collection<RefSeqGene> left=leftAlignments.get(readName);
			Collection<RefSeqGene> right=rightAlignments.get(readName);
			Collection<PairedEndMapping> pairs=getAllPossiblePairs(left, right, readName);
			rtrn.put(readName, pairs);
		}
		
		return rtrn;
	}
	
	private static Collection getAllPossiblePairs(Collection<RefSeqGene> left, Collection<RefSeqGene> right, String readName){
		Collection<PairedEndAlignment> rtrn=new HashSet();
		for(RefSeqGene leftAlign: left){
			for(RefSeqGene rightAlign: right){
				PairedEndAlignment paired=new PairedEndAlignment(leftAlign, rightAlign, readName);
				rtrn.add(paired);
			}
		}
		return rtrn;
	}
	
	private static Collection getAllPossiblePairs(RefSeqGene leftAlign, Collection<RefSeqGene> right, String readName){
		Collection<PairedEndAlignment> rtrn=new HashSet();
		for(RefSeqGene rightAlign: right){
			PairedEndAlignment paired=new PairedEndAlignment(leftAlign, rightAlign, readName);
			rtrn.add(paired);
		}
		return rtrn;
	}
	
	
	private void write(String save, Map<String, Collection> filteredPairs)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(String name: filteredPairs.keySet()){
			Collection<PairedEndAlignment> list=filteredPairs.get(name);
			for(PairedEndAlignment paired: list){
				writer.write(paired.toString()+"\n");
			}
		}
		
		writer.close();
	}
	
	
	public static Collection<PairedEndAlignment> getAllPairs(File leftFile, File rightFile, String chr)throws IOException{
		Collection<PairedEndAlignment> rtrn=new TreeSet();
		
		Map<String, Collection> leftAlignments=SAMToAligned.parseSAMByName(leftFile, chr);
		Map<String, Collection> rightAlignments=SAMToAligned.parseSAMByName(rightFile, chr);
		
		for(String readName: leftAlignments.keySet()){
			Collection<RefSeqGene> left=leftAlignments.get(readName);
			Collection<RefSeqGene> right=rightAlignments.get(readName);
			if(left!=null && right!=null && left.size()==1 && right.size()==1){//TODO MG: Added support for unique only mapping
				Collection<PairedEndAlignment> pairs=getAllPossiblePairs(left, right, readName);
				for(PairedEndAlignment pair: pairs){if(goodPair(pair)){rtrn.add(pair);}}
			}
		}
		return rtrn;
	}
	
	public static Collection<PairedEndAlignment> getAllPairs(AlignmentDataModel data1, AlignmentDataModel data2, String chr)throws IOException{
		Collection<PairedEndAlignment> rtrn=new TreeSet();
		
		Alignments chrAlign=new Alignments(chr, 0, data1.getChromosomeLengths().get(chr.replaceAll("chr", "")));
		CloseableIterator<Alignment> iter=data1.getAlignmentsOverlappingRegion(chrAlign);
		
		Map<String, Collection> rightAlignments=makeCollection(data2.getAlignmentsOverlappingRegion(chrAlign));
		
		
		while(iter.hasNext()){
			Alignment left=iter.next();
			String readName=left.getReadName();
			Collection<RefSeqGene> rightCollection=rightAlignments.get(readName);
			
			
			
			if(left!=null && rightCollection!=null){
				Collection<PairedEndAlignment> pairs=getAllPossiblePairs(new RefSeqGene(left), rightCollection, readName);
				for(PairedEndAlignment pair: pairs){if(goodPair(pair)){rtrn.add(pair);}}
			}
			
		}
		
		return rtrn;
	}
	
	
	private static Map<String, Collection> makeCollection(CloseableIterator<Alignment> iter){
		Map<String, Collection> rtrn=new TreeMap();
		
		while(iter.hasNext()){
			Alignment align=iter.next();
			Collection list=new TreeSet();
			if(rtrn.containsValue(align.getReadName())){list=rtrn.get(align.getReadName());}
			list.add(new RefSeqGene(align));
			rtrn.put(align.getReadName(), list);
		}
		
		return rtrn;
	}
	
	
	/*public static Collection<PairedEndAlignment> getAllPairs(File leftFile, File rightFile, String chr){
		Collection<PairedEndAlignment> rtrn=new TreeSet();
		
		Map<String, Collection> leftAlignments=SAMToAligned.parseSAMByName(leftFile, chr);
		Map<String, Collection> rightAlignments=SAMToAligned.parseSAMByName(rightFile, chr);
		
		for(String readName: leftAlignments.keySet()){
			Collection<RefSeqGene> left=leftAlignments.get(readName);
			Collection<RefSeqGene> right=rightAlignments.get(readName);
			if(left!=null && right!=null){
			Collection<PairedEndAlignment> pairs=getAllPossiblePairs(left, right, readName);
			rtrn.addAll(pairs);
			//for(PairedEndAlignment pair: pairs){if(goodPair(pair)){writer.write(pair+"\n");}}
			}
		}
		return rtrn;
	}*/
	
	public static void main(String[] args)throws IOException{
		if(args.length>3){
			File left=new File(args[0]);
			File right=new File(args[1]);
			String save=args[2];
			String chr=args[3];
			new PairedEndMapping(SAMToAligned.parseSAMByName(left, chr), SAMToAligned.parseSAMByName(right, chr), save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=left read alignments (SAM format) \n args[1]=right read alignments (SAM Format) \n args[2]=save file \n args[3]=chr";
}
