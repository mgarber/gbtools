package broad.core.alignment;

import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.GenomicAnnotation;

public class SmatchAlignmentSummary extends AlignmentSummary{

	public SmatchAlignmentSummary(String[] primerData, String[] queryData) {

		super();

		String primerId = primerData[1].trim();
		String primerSeq = primerData[2].trim();

		GenomicAnnotation primer = new BasicGenomicAnnotation(primerId);
		primer.setChromosome(primer.getName().length() > 3 ? primer.getName().substring(3) : primer.getName());

		GenomicAnnotation query = new BasicGenomicAnnotation(queryData[1].trim());
		query.setChromosome(query.getName().length() > 3 ? query.getName().substring(3) : query.getName());

		int numMismatches = Integer.parseInt(queryData[4].trim());

		String strand = queryData[6].trim();
		query.setOrientation(strand);
		if (strand.equals("-"))
			setReversedOrientation(true);

		setAlignmentLength(queryData[12].trim().length());

		int start = Integer.parseInt(queryData[11].trim());
		query.setStart(start);
		query.setEnd(start + getAlignmentLength() - 1);

		String querySeq = queryData[12].trim();

		primer.setStart(primerSeq.indexOf(querySeq.toLowerCase()) + 1);
		primer.setEnd(primer.getStart() + getAlignmentLength() - 1);

		setPid((getAlignmentLength() - numMismatches) / getAlignmentLength());
		setGapOpenings(0);
		setScore(getAlignmentLength() - numMismatches);

		setA(query);
		setB(primer);
	}
}

