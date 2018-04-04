package broad.core.alignment;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class Arachne454AlignmentReader {
	ArrayList<AlignmentSummary> alignments;

	public Arachne454AlignmentReader(String inputFileName) throws IOException {
		super();
		alignments = new ArrayList<AlignmentSummary>();
		BufferedReader br = new BufferedReader(new FileReader(inputFileName));
		String line = null;
		while((line = br.readLine()) != null) {
			String[] lineInfo = line.split(", ");
			AlignmentSummary aln = new AlignmentSummary();
			String [] seqInfo = lineInfo[0].split(" vs ");
			String subject = seqInfo[0].substring(0, seqInfo[0].length() - 2);
			String orientation = seqInfo[0].substring(seqInfo[0].length() - 3);
			aln.setSubject(subject);
			aln.setQuery(seqInfo[1]);
			aln.setReversedOrientation("rc".equals(orientation));

			String [] startEndInfo = lineInfo[2].split(" to ");
			String [] subjStartEndInfo = startEndInfo[0].substring(5).split("-");
			String [] queryStartEndInfo = startEndInfo[1].split("-");

			aln.setSubjectStart(Integer.parseInt(subjStartEndInfo[0]));
			aln.setSubjectEnd(Integer.parseInt(subjStartEndInfo[1]));
			aln.setQueryStart(Integer.parseInt(queryStartEndInfo[0]));
			aln.setQueryEnd(Integer.parseInt(queryStartEndInfo[1]));
			
			alignments.add(aln);
		}
		br.close();
	}
	
	public List<AlignmentSummary> getAlignments() { return alignments;}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
