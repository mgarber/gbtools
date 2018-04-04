package broad.pda.rnai;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Map;
import java.util.TreeMap;

import broad.core.datastructures.IntervalTree;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;

public class AnnotateRNAiReport {

	String header;
	
	public AnnotateRNAiReport(File rnaiFile, File annotationFile, String save) throws IOException{
		Map<Alignments, String> rnaiInfo=parse(rnaiFile);
		Map<String, IntervalTree<Alignments>> annotationTree=BEDFileParser.loadAlignmentDataToTree(annotationFile);
		write(save, rnaiInfo, annotationTree);
	}
	
	private void write(String save, Map<Alignments, String> rnaiInfo, Map<String, IntervalTree<Alignments>> annotationTree) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write(header+"\n");
		
		for(Alignments region: rnaiInfo.keySet()){
			boolean overlaps=annotationTree.get(region.getChr()).overlappers(region.getStart(), region.getEnd()).hasNext();
			String line=rnaiInfo.get(region);
			writer.write(line+"\t"+overlaps+"\n");
		}
		
		writer.close();
	}

	private Map<Alignments, String> parse(File rnaiFile) throws IOException {
		Map<Alignments, String> rtrn=new TreeMap<Alignments, String>();
		
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(rnaiFile)));
		String nextLine;
		int i=0;
        while ((nextLine = reader.readLine()) != null) {
        	if(i!=0){
        	String[] tokens=nextLine.split("\t");
        	Alignments align=new Alignments("chr"+tokens[16], new Integer(tokens[19]), new Integer(tokens[20]));
        	rtrn.put(align, nextLine);
        	}
        	else{header=nextLine;}
        	i++;
        }
		
		return rtrn;
	}

	public static void main(String[] args)throws IOException{
		if(args.length>2){
			File rnaiFile=new File(args[0]);
			File files=new File(args[1]);
			String save=args[2];
			new AnnotateRNAiReport(rnaiFile, files, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=rnai file \n args[1]=annotation File \n args[2]=save";
	
}
