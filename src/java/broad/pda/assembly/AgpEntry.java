package broad.pda.assembly;

import broad.core.sequence.SequenceRegion;

public class AgpEntry extends SequenceRegion{

	static public final int CLONE_TYPE       = 0;
	static public final int CENTROMERE_TYPE  = 1;
	static public final int GAP_TYPE         = 2;
	static public final int CONTIG_TYPE      = 4;
	static public final int OTHER_TYPE		 = 3;
	public static final int SHORT_ARM_TYPE   = 5;
	private boolean reversedOrientation;
	private int type;
	private int number;
	
	public AgpEntry(String parentSequence) {
		super(parentSequence);
	}
	

	protected void setInReverseOrientation(boolean isInreverseOrientation) {
		this.reversedOrientation = isInreverseOrientation;
	}
	public boolean inReversedOrientation() {
		return reversedOrientation;
	}
		
	public void setType(int type) {
		this.type = type;
	}
	
	public int getType() {
		return type;
	}
	
	public int getLength() { return super.getLength() + 1;}
	
	public String toString() {
		StringBuffer buf = new StringBuffer("chr"+getContainingSequenceId());
		
		buf.append("\t")
			.append(getStart())
			.append("\t")
			.append(getEnd())
			.append("\t")
			.append(getNumber())
			.append("\t");
		
		switch (getType()) {
		case AgpEntry.CLONE_TYPE:
			buf.append("F\t");
			break;
		case AgpEntry.CENTROMERE_TYPE:
		case AgpEntry.SHORT_ARM_TYPE:
		case AgpEntry.GAP_TYPE:
			buf.append("N\t");
			break;
		case AgpEntry.CONTIG_TYPE:
			buf.append("W\t");
			break;
		default:
			buf.append("O\t");
			break;
		}
		
		buf.append(getName())
			.append("\t")
			.append("0\t")
			.append(getLength() + 1)
			.append("\t+");
		
		return buf.toString();
	}


	public int getNumber() {
		return number;
	}


	public void setNumber(int number) {
		this.number = number;
	}



}


