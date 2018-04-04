package broad.pda.seq.alignment.sam;

public class SAMTags {

	private int numberOfMismatches;
	private double weight;
	private int strata;
	private String mismatchString;
	private int numOtherAlignments;
	private String programName;
	
	public SAMTags(){
		numberOfMismatches=0;
		weight=1.0;
		strata=0;
		mismatchString="76";
		numOtherAlignments=0;
	}

	public void setTags(String vals){
		try{
		String[] tokens=vals.split(":");
		String tag=tokens[0];
		String type=tokens[1];
		String value=tokens[2];
		if(tag.equalsIgnoreCase("NM")){this.numberOfMismatches=new Integer(value);}
		if(tag.equalsIgnoreCase("XW")){this.weight=new Double(value);}
		if(tag.equalsIgnoreCase("XA")){this.strata=new Integer(value);}
		if(tag.equalsIgnoreCase("MD")){this.mismatchString=value;}
		if(tag.equalsIgnoreCase("XM")){this.numOtherAlignments=new Integer(value);}
		if(tag.equalsIgnoreCase("PG")){this.programName=value;}
		}catch(Exception ex){ex.printStackTrace();} //TODO Need to fix this
	}
	
	public String toString(){
		String rtrn="NM:i:"+numberOfMismatches+"\tXW:f:"+new Double(weight).floatValue()+"\tXA:i:"+strata+"\tMD:Z:"+mismatchString+"\tXM:i:"+numOtherAlignments;
		if(this.programName!=null && !programName.isEmpty()){rtrn+="\tPG:Z:"+this.programName;}
		return rtrn;
	}
	
	public void setWeight(double weight){this.weight=weight;}
	
	public int getNumberOfMismatches() {return numberOfMismatches;}
	
	public int getStrata(){return strata;}
	
	public String getProgramName(){return this.programName;}

	public void setProgram(String name) {
		this.programName=name;
	}

	public void setMismatches(int mismatches) {
		this.numberOfMismatches=mismatches;
	}
}
