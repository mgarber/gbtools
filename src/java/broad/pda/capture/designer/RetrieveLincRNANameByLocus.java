package broad.pda.capture.designer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;

public class RetrieveLincRNANameByLocus {

	
	public RetrieveLincRNANameByLocus(Collection<Alignments> loci, File reportFile, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		Map<String, IntervalTree<Alignments>> report=parseReportToTree(reportFile); //Make sure the alignments has the id of the lincRNA
		
		for(Alignments region: loci){
			Iterator<Node<Alignments>> over=report.get(region.getChr()).overlappers(region.getStart(), region.getEnd());
			writer.write(region.toString());
			while(over.hasNext()){
				Alignments gene=over.next().getValue();
				writer.write("\t"+gene.getName());
			}
			writer.write("\n");
		}
		writer.close();
	}

	private String getIds(Alignments region, Map<String, IntervalTree<Alignments>> report, String strand) {
		String rtrn="";
		Iterator<Node<Alignments>> over=report.get(region.getChr()).overlappers(region.getStart(), region.getEnd());
		Collection<String> names=new TreeSet<String>();
		while(over.hasNext()){
			Alignments align=over.next().getValue();
			if(align.getOrientation().equals(strand)){
				names.add(align.getName());
			}
		}
		
		int counter=0;
		for(String name: names){
			rtrn+=name;
			if(names.size()>1 && (counter<names.size()-1)){rtrn+="-";}
			counter++;
		}
		
		return rtrn;
	}

	private Collection<String> parseList(File listFile) throws IOException {
	
		Collection<String> rtrn=new TreeSet<String>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(listFile)));
		String nextLine;
		 while ((nextLine = reader.readLine()) != null) {
        	rtrn.add(nextLine);
        }
		
        return rtrn;
	}

	private Map<String, IntervalTree<Alignments>> parseReportToTree(File reportFile) throws IOException {
		Map<String, IntervalTree<Alignments>> rtrn=new TreeMap<String, IntervalTree<Alignments>>();
		
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(reportFile)));
		String nextLine;
		int i=0;
        while ((nextLine = reader.readLine()) != null) {
        	if(i!=0){
	        	String[] tokens=nextLine.split("\t");
	        	String id=tokens[0];
	        	Alignments region=new Alignments(tokens[3], new Integer(tokens[4]), new Integer(tokens[5]));
	        	region.setOrientation(tokens[6]);
	        	region.setName(id);
	        	IntervalTree<Alignments> tree=new IntervalTree<Alignments>();
	        	if(rtrn.containsKey(region.getChr())){tree=rtrn.get(region.getChr());}
	        	tree.put(region.getStart(), region.getEnd(), region);
	        	rtrn.put(region.getChr(), tree);
        	}
        	i++;
        }
		
		return rtrn;
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
			Collection<Alignments> loci=BEDFileParser.loadAlignmentData(new File(args[0]));
			File reportFile=new File(args[1]);
			String save=args[2];
			new RetrieveLincRNANameByLocus(loci, reportFile, save);
		}
		else{
			System.err.println(usage);
		}
	}
	
	static String usage=" args[0]=list of regions (BED) \n args[1]=report file \n args[2]=save";
	
}
