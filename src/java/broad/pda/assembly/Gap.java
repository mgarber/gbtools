package broad.pda.assembly;

public class Gap extends AgpEntry {
	private static final String GAP = "gap";
	public Gap(AgpEntry entry) {
		super(entry.getContainingSequenceId());
		setName(GAP);
		setType(AgpEntry.GAP_TYPE);
		setStart(entry.getStart());
		setEnd(entry.getEnd());
		setChromosome(entry.getChromosome());
	}
	
	public Gap(String containingSeqId) {
		super(containingSeqId);
	}

	public void setSequenceBases(String bases) {
		// do nothing.
	}
	
	public String getSequenceBases() {
		StringBuffer buf = new StringBuffer((int)(getEnd() - getStart()));
		for (int i = 1; i <= (int)(getEnd() - getStart()) + 1; i++) {
			buf.append("N");
		}
		return buf.toString();
	}
	
	public int getLength() { 
		return (int) (getEnd() - getStart());
	}
}
