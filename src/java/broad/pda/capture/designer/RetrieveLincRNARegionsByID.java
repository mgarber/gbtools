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
import broad.pda.gene.RefSeqGeneWithIsoforms;

public class RetrieveLincRNARegionsByID {

	
	public RetrieveLincRNARegionsByID(File listFile, File reportFile, Map<String, IntervalTree<RefSeqGeneWithIsoforms>> genes, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		Collection<String> list=parseList(listFile);
		Map<String, Alignments> report=parse(reportFile);
		
		for(String id: list){
			Alignments region=report.get(id);
			if(region!=null){
				Iterator<Node<RefSeqGeneWithIsoforms>> over=genes.get(region.getChr()).overlappers(region.getStart(), region.getEnd());
				while(over.hasNext()){
					RefSeqGeneWithIsoforms gene=over.next().getValue();
					writer.write(gene.toString(id));
				}
			}
			else{System.err.println("Skipping "+id+" because wasnt in report");}
		}
		writer.close();
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

	private Map<String, Alignments> parse(File reportFile) throws IOException {
		Map<String, Alignments> rtrn=new TreeMap<String, Alignments>();
		
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(reportFile)));
		String nextLine;
		int i=0;
        while ((nextLine = reader.readLine()) != null) {
        	if(i!=0){
	        	String[] tokens=nextLine.split("\t");
	        	String id=tokens[0];
	        	Alignments region=new Alignments(tokens[3], new Integer(tokens[4]), new Integer(tokens[5]));
	        	region.setOrientation(tokens[6]);
	        	rtrn.put(id, region);
        	}
        	i++;
        }
		
		return rtrn;
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>3){
			File list=new File(args[0]);
			File reportFile=new File(args[1]);
			BEDFileParser parser = new BEDFileParser(args[2]);
			Map<String, IntervalTree<RefSeqGeneWithIsoforms>> genes= parser.getAnnotationSetMap();
			String save=args[3];
			new RetrieveLincRNARegionsByID(list, reportFile, genes, save);
		}
		else{
			System.err.println(usage);
		}
	}
	
	static String usage=" args[0]=list of IDs \n args[1]=report file \n args[2]=BED file \n args[3]=save";
	
}
