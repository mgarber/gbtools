package broad.core.gene;

import java.io.IOException;
import java.util.List;

import broad.core.annotation.GenomicAnnotation;

public class GenericGeneAnnotationReader extends GeneAnnotationReader {

	public GenericGeneAnnotationReader(String gffExonFile, String chr) throws IOException {
		super.loadExonGffFile(gffExonFile, chr);
	}
	@Override
	public List<GeneAnnotation> getAnnotationsForSymbol(String symbol) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getSymbolsForGene(String id) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void loadAliasTable(String fileName) throws IOException {
		// TODO Auto-generated method stub

	}
	@Override
	public GeneAnnotation createAnnotation(GenomicAnnotation a) {
		// TODO Auto-generated method stub
		return null;
	}

}
