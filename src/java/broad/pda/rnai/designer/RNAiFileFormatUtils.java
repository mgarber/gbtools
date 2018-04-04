package broad.pda.rnai.designer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.pda.rnai.RNAiGeneAnnotation;

public class RNAiFileFormatUtils {

	
	public static Collection<RNAiGeneAnnotation> parseRNAiReportFile(File rnaiReport) throws IOException {
		Collection<RNAiGeneAnnotation> rtrn=new ArrayList<RNAiGeneAnnotation>();	
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(rnaiReport)));
		String nextLine;
		int i=0;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			if(i>0){rtrn.add(new RNAiGeneAnnotation(nextLine));}
			i++;
		}
		return rtrn;
	}

	//Return map of isoform name and a collection of the hairpins designed
	public static Map<String, Collection<HairpinKmer>> parseHairpinDesignFiles(File[] hairpinFiles) throws IOException {
		Map<String, Collection<HairpinKmer>> rtrn=new TreeMap();
		
		for(int i=0; i<hairpinFiles.length; i++){
			File file=hairpinFiles[i];
			String isoformName=file.getName().split("\\.")[0];
			Collection<HairpinKmer> hps=parseHairpinFile(file);
			rtrn.put(isoformName, hps);
		}
		
		return rtrn;
	}

	private static Collection<HairpinKmer> parseHairpinFile(File file) throws IOException {
		Collection<HairpinKmer> rtrn=new TreeSet();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			String[] tokens=nextLine.split("\t");
			HairpinKmer hp=new HairpinKmer(tokens[0], new Integer(tokens[1]), new Integer(tokens[1])+21);
			hp.setRS8Score(new Double(tokens[2]));
			hp.setOriginalScore(new Double(tokens[3]));
			hp.setMirScore(new Double(tokens[4]));
			rtrn.add(hp);
		}
		
		return rtrn;
	}

	public static Map<String, RNAiGeneAnnotation> parseRNAiReportFileBySequence(File rnaiReport) throws IOException {
		Collection<RNAiGeneAnnotation> rnai=parseRNAiReportFile(rnaiReport);
		Map<String, RNAiGeneAnnotation> rtrn=new TreeMap();
		for(RNAiGeneAnnotation annotation: rnai){
			rtrn.put(annotation.getSequence().getSequenceBases().toUpperCase(), annotation);
		}
		
		return rtrn;
	}
	
}
