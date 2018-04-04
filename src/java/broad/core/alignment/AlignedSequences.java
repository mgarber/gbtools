package broad.core.alignment;

import broad.core.sequence.Sequence; 

public class AlignedSequences extends AlignmentSummary {
	Sequence querySequence;
	Sequence subjectSequence;
	
	public AlignedSequences() {
		super();
		
		querySequence = new Sequence("query");
		subjectSequence = new Sequence("subject");
	}
	
	public void setQuerySequenceString(String queryAlignedString) {
		querySequence.setSequenceBases(queryAlignedString);
	}
	
	public void setSubjectSequenceString(String subjectAlignedString) {
		subjectSequence.setSequenceBases(subjectAlignedString);
	}
	
	public String getQuerySequenceString() { return querySequence.getSequenceBases();}
	public String getSubjectSequenceString() { return subjectSequence.getSequenceBases();}
	
	public boolean isConsistent() {
		return querySequence.getSequenceBases() != null && subjectSequence.getSequenceBases() != null && 
		 querySequence.getSequenceBases().length() == subjectSequence.getSequenceBases().length();
	}
	
}
