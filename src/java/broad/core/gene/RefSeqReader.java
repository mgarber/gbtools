package broad.core.gene;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.GenomicAnnotationFilter;
import broad.core.error.ParseException;

public class RefSeqReader extends GeneAnnotationReader {
	private RefLinkReader refLink;
	
	public RefSeqReader(String fileName, String chr, boolean loadCDSOnly) 
	throws ParseException, IOException {
		super.loadChromosomeGeneAnnotations(fileName, chr, loadCDSOnly);
	}
	
	public RefSeqReader(String fileName, String chr, boolean loadCDSOnly, String version) 
	throws ParseException, IOException {
		super.loadChromosomeGeneAnnotations(fileName, chr, loadCDSOnly, version);
	}
	
	public RefSeqReader(String fileName, boolean loadCDSOnly) throws ParseException, IOException {
		super.load(fileName, loadCDSOnly, "HG17");
	}
	
	public RefSeqReader(String fileName, boolean loadCDSOnly, String version) throws ParseException, IOException {
		super.load(fileName, loadCDSOnly, version);
	}
	
	public RefSeqReader(String fileName, boolean loadCDSOnly, String version, GenomicAnnotationFilter<GeneAnnotation> filter) throws ParseException, IOException {
		super.load(fileName, loadCDSOnly, version, filter);
	}

	@Override
	public void loadAliasTable(String fileName) throws IOException {
		refLink = new RefLinkReader(fileName);
		
	}
	
	public String getRefSeqForEntrez(int entrezId) {
		return refLink.getRefSeqForEntrezId(entrezId);
	}
	
	public List<String> getSymbolsForGene(String id) {
		ArrayList<String> result = null;
		if (refLink == null) {
			return null;
		}
		
		String name = refLink.getRefSeqName(id);
		if(name != null) {
			result = new ArrayList<String>();
			result.add(name);
		}
		return result;
	}

	@Override
	public List<GeneAnnotation> getAnnotationsForSymbol(String symbol) {
		List<String> refSeqsForSymbol = refLink.getRefSeqsForSymbol(symbol);
		ArrayList<GeneAnnotation> annotations = new ArrayList<GeneAnnotation>(refSeqsForSymbol.size());
		Iterator<String> it = refSeqsForSymbol.iterator();
		while (it.hasNext()) {
			String id = it.next();
			GeneAnnotation gene = super.getGene(id);
			if(gene == null) {
				System.out.println("ERROR: expected refseq id<" + id +"> for symbol <"+ symbol.toUpperCase() + "> was not found in loaded refSeq table");
			} else {
				annotations.add(gene);
			}
		}
		return annotations;
	}

	@Override
	public GeneAnnotation createAnnotation(GenomicAnnotation a) {
		// TODO Auto-generated method stub - modify parent class to delegate annotation creation logic here.
		return null;// This is not used yet, 
	}

}
