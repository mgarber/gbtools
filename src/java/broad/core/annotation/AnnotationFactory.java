package broad.core.annotation;

import broad.core.error.ParseException;

public interface AnnotationFactory<T extends LightweightGenomicAnnotation> {
	T create(String [] rawFields) throws ParseException;
	T create(GenomicAnnotation baseAnnotaiton);
	T create(String name);
}
