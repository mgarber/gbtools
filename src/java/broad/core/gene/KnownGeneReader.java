package broad.core.gene;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import broad.core.annotation.GenomicAnnotation;
import broad.core.error.ParseException;

public class KnownGeneReader extends GeneAnnotationReader {
	private HashMap<String, List<String>> symbolToKnownGeneIds = new HashMap<String, List<String>>();
	private HashMap<String, List<String>> knownGeneIdSymbols   = new HashMap<String, List<String>>();
	
	
	public KnownGeneReader(String fileName, String chr,  boolean loadCDSOnly, String version) 
	throws ParseException, IOException {
		super.loadChromosomeGeneAnnotations(fileName, chr, loadCDSOnly);
	}
	
	
	public KnownGeneReader(String fileName, String chr,  boolean loadCDSOnly) 
	throws ParseException, IOException {
		this.loadChromosomeGeneAnnotations(fileName, chr, loadCDSOnly, "HG17");
	}

	public KnownGeneReader(String fileName,  boolean loadCDSOnly, String version) throws ParseException, IOException {
		super.load(fileName, loadCDSOnly); //known gene format is the same for both HG18 and HG17
	}

	public KnownGeneReader(String fileName,  boolean loadCDSOnly) throws ParseException, IOException {
		this.load(fileName, loadCDSOnly, "HG17");
	}

	@Override
	public void loadAliasTable(String fileName) throws IOException {
		File source = new File(fileName);
		BufferedReader br = new BufferedReader(new FileReader(source));
		String line;
		while((line = br.readLine()) != null) {
			if(line.startsWith("#")){
				continue;
			}
			String[] splitLine = line.split("\t");
			String symbol = splitLine[1];
			List<String> symbolListForId = knownGeneIdSymbols.get(splitLine[0]);
			if(symbolListForId == null) {
				symbolListForId = new ArrayList<String>(0);
				knownGeneIdSymbols.put(splitLine[0], symbolListForId);
			}
			symbolListForId.add(symbol);
			
			List<String> knownGeneListForSymbol = symbolToKnownGeneIds.get(symbol);
			if(knownGeneListForSymbol == null) {
				knownGeneListForSymbol = new ArrayList<String>();
				symbolToKnownGeneIds.put(symbol, knownGeneListForSymbol);
			}
			knownGeneListForSymbol.add(splitLine[0]);
		}
		System.out.print("Closing "+fileName);
		br.close();
		System.out.print(" ..... Closed\n");	
	}

	@Override
	public List<GeneAnnotation> getAnnotationsForSymbol(String symbol) {
		List<String> knwonGeneForSymbol = symbolToKnownGeneIds.get(symbol);
		if(knwonGeneForSymbol == null) {
			return new ArrayList<GeneAnnotation>();
		}
		//System.out.println("knwonGeneForSymbol " + knwonGeneForSymbol);
		ArrayList<GeneAnnotation> annotations = new ArrayList<GeneAnnotation>(knwonGeneForSymbol.size());
		Iterator<String> it = knwonGeneForSymbol.iterator();
		while (it.hasNext()) {
			String id = it.next();
			GeneAnnotation gene = super.getGene(id);
			if(gene == null) {
				System.out.println("ERROR: expected knwon gene id<" + id +"> for symbol <"+ symbol + "> was not found in loaded refSeq table");
			} else {
				annotations.add(gene);
			}
		}
		return annotations;
	}

	public List<String> getSymbolsForGene(String id) {
		if (knownGeneIdSymbols == null) {
			return null;
		}
		
		return knownGeneIdSymbols.get(id);
	}


	@Override
	public GeneAnnotation createAnnotation(GenomicAnnotation a) {
		// TODO Auto-generated method stub modify parent class to delegate gene annotation creation here.
		return null;
	}
}
