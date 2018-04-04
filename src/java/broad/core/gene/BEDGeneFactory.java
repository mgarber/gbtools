package broad.core.gene;

import broad.core.annotation.AnnotationFactory;
import broad.core.annotation.AnnotationFactoryFactory;
import broad.core.annotation.BED;
import broad.core.annotation.GenomicAnnotation;
import broad.core.error.ParseException;

public class BEDGeneFactory implements AnnotationFactory<GeneAnnotation> {
	
	public GeneAnnotation create(String[] rawFields) throws ParseException {
		BED bed = AnnotationFactoryFactory.bedFactory.create(rawFields);
		return new GeneAnnotation(bed);
	}

	public GeneAnnotation create(GenomicAnnotation baseAnnotation) {
		return null;
	}

	public GeneAnnotation create(String name) {
		// TODO Auto-generated method stub
		return null;
	}
	

}
