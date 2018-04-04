package broad.core.alignment;

import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.GenomicAnnotation;

public  class PHOneLineAlignment extends AlignmentSummary {
	public PHOneLineAlignment(String [] rawData) {
		super();
		String queryName = rawData[0].trim();
		if(queryName.contains("chr")) {
			int chrStart = queryName.indexOf("chr");
			queryName = queryName.substring(chrStart);
		}
		GenomicAnnotation a = new BasicGenomicAnnotation(queryName);
		a.setChromosome(a.getName().startsWith("chr") && a.getName().length() < 6 ? a.getName().substring(3) : a.getName());
		GenomicAnnotation b = new BasicGenomicAnnotation(rawData[1].trim());
		b.setChromosome(a.getName().startsWith("chr") && a.getName().length() < 6 ? b.getName().substring(3) : b.getName());
		int i = 2;
		setPid(Float.parseFloat(rawData[i++].trim()));
		setAlignmentLength(Integer.parseInt(rawData[i++].trim()));
		setGapOpenings(Integer.parseInt(rawData[i++].trim()));
		setMissmatches(getAlignmentLength() - Integer.parseInt(rawData[i++].trim()));

		int start = Integer.parseInt(rawData[i++].trim());
		int end   = Integer.parseInt(rawData[i++].trim());
		if(start < end) {
			a.setStart(start);
			a.setEnd(end);
		} else {
			a.setStart(end);
			a.setEnd(start);
			a.setOrientation("-");
			setReversedOrientation(true);
		}
		
		start = Integer.parseInt(rawData[i++].trim());
		end   = Integer.parseInt(rawData[i++].trim());
		if(start < end) {
			b.setStart(start);
			b.setEnd(end);
		} else {
			b.setStart(end);
			b.setEnd(start);
			b.setOrientation("-");
			setReversedOrientation(true);
		}
		
		setEValue(Float.parseFloat(rawData[i++].trim()));
		setScore(Float.parseFloat(rawData[i++].trim()));

		setA(a);
		setB(b);
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder(getA().getName());
		sb.append("\t").append(getB().getName())
			.append("\t").append(getPid())
			.append("\t").append(getAlignmentLength())
			.append("\t").append(getGapOpenings())
			.append("\t").append(getMissmatches())
			.append("\t");
			if(!getA().inReversedOrientation()) {
				sb.append(getA().getStart()).append("\t")
					.append(getA().getEnd());
			} else {
				sb.append(getA().getEnd()).append("\t")
					.append(getA().getStart());
			}
			sb.append("\t");
			sb.append(getB().getStart()).append("\t").append(getB().getEnd());
			
			sb.append("\t").append(getEValue()).append("\t").append(getScore());
			
		return sb.toString();
	}
}
