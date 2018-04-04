package broad.core.alignment;

import java.text.DecimalFormat;

import broad.core.annotation.BasicGenomicAnnotation;

public class Repeat extends BasicGenomicAnnotation {
	static final int MAX_SCORE_SIZE = 5;
	static final int MAX_SEQ_NAME_SIZE = 25;
	static final int MAX_GENOMIC_POS_SIZE = 11;
	static final int MAX_REPEAT_POS_SIZE  = 7;
	
	float percentDivergence;
	private float percentDeletions;
	private float percentInsertions;
	private String query;
	private int queryBasesLeft;
	private String repeatName;
	private String repeatClass;
	private String repeatFamily;
	private int repeatStart;
	private int repeatEnd;
	private boolean isRemarkable;
	private int repeatBasesLeft;
	private String id;
	
	Repeat(String [] data) {
		int i = 0;//data[0].length() == 0 ? 1 : 0;
		setScore(Integer.parseInt(data[i++]));
		this.percentDivergence = Float.parseFloat(data[i++]);
		this.percentDeletions  = Float.parseFloat(data[i++]);
		this.percentInsertions = Float.parseFloat(data[i++]);
		this.query = data[i++].intern();
		setChromosome(query.length() > 3 ? query.substring(3) : query);
		setStart(Integer.parseInt(data[i++]) - 1); // -1 because Repeat masker output is 1 based
		setEnd(Integer.parseInt(data[i++]));   // no -1 because Repeat masker uses closed rather than semi-closed intervals
		this.queryBasesLeft = Integer.parseInt(data[i].substring(1,data[i++].length() - 1));
		setOrientation("C".equals(data[i++]) ? "-" : "+");
		this.repeatName = data[i++].intern();
		this.repeatFamily = data[i++].intern();
		// Repeat start and end are tricky, depends on whether hit is direct or not
		if(inReversedOrientation()) {
			this.repeatStart = Integer.parseInt(data[13]);
			this.repeatEnd   = Integer.parseInt(data[12]);
			this.repeatBasesLeft  = Integer.parseInt(data[11].substring(1,data[11].length() - 1));
		} else {
			this.repeatStart = Integer.parseInt(data[11] );
			this.repeatEnd = Integer.parseInt(data[12]);
			this.repeatBasesLeft  = Integer.parseInt(data[13].substring(1,data[13].length() - 1));
		}
		this.id = data[14];
		this.isRemarkable = data.length == 16 && "*".equals(data[15]);
	}
	
	public Repeat() {
		super();
	}
	
	public void initFromUCSCTableData(String [] data) {
		int i = 0;//data[0].length() == 0 ? 1 : 0;
		setScore(Integer.parseInt(data[i++]));
		this.percentDivergence = Float.parseFloat(data[i++])/10f;
		this.percentDeletions  = Float.parseFloat(data[i++])/10f;
		this.percentInsertions = Float.parseFloat(data[i++])/10f;
		this.query = data[i++].intern();
		setChromosome(query.length() > 3 ? query.substring(3) : query);
		setStart(Integer.parseInt(data[i++]) - 1); // -1 because Repeat masker output is 1 based
		setEnd(Integer.parseInt(data[i++]));   // no -1 because Repeat masker uses closed rather than semi-closed intervals
		this.queryBasesLeft = -Integer.parseInt(data[i++]);
		setOrientation(data[i++]);
		this.repeatName = data[i++].intern();
		this.repeatClass = data[i++].intern();
		this.repeatFamily = data[i++].intern();
		// Repeat start and end are tricky, depends on whether hit is direct or not
		if(inReversedOrientation()) {
			this.repeatStart = Integer.parseInt(data[14]);
			this.repeatEnd   = Integer.parseInt(data[13]);
			this.repeatBasesLeft  = -Integer.parseInt(data[12]);
		} else {
			this.repeatStart = Integer.parseInt(data[12] );
			this.repeatEnd = Integer.parseInt(data[13]);
			this.repeatBasesLeft  = -Integer.parseInt(data[14]);
		}
		this.id = data[15];
	}
	
	public void setQueryBasesLeft(int basesLeft) {
		queryBasesLeft = basesLeft;
	}

	public String format() {
		DecimalFormat df = new DecimalFormat();
		df.setMaximumFractionDigits(1);
		
		StringBuffer buf = new StringBuffer(formatNumber(getScore(), MAX_SCORE_SIZE, false));
		
		buf.append(" ")
			.append(rightJustify(df.format(percentDivergence),4))
			.append(" ")
			.append(rightJustify(df.format(percentDeletions),4))
			.append(" ")
			.append(rightJustify(df.format(percentInsertions),4))
			.append(" ")
			.append(leftJustify(query, MAX_SEQ_NAME_SIZE))
			.append(" ");

			
		buf.append(formatNumber(getStart() + 1, MAX_GENOMIC_POS_SIZE, false)) // + 1 because Repeat masker output is 1 based
			.append(" ")
			.append(formatNumber(getEnd() , MAX_GENOMIC_POS_SIZE, false))  // no + 1 because Repeat masker uses closed rather than semi-closed intervals
			.append(" ")
			.append(formatNumber(getQueryBasesLeft(), MAX_GENOMIC_POS_SIZE, true))
			.append(" ");

		buf.append(inReversedOrientation() ? "C" : "+")
			.append(" ")
			.append(leftJustify(repeatName, MAX_SEQ_NAME_SIZE))
			.append(" ")
			.append(leftJustify(repeatFamily, MAX_SEQ_NAME_SIZE))
			.append(" ");
		
		if(inReversedOrientation()) {
			buf.append(formatNumber(repeatBasesLeft, MAX_REPEAT_POS_SIZE, true))
				.append(" ")
				.append(formatNumber(repeatEnd, MAX_REPEAT_POS_SIZE, false))
				.append(" ")
				.append(formatNumber(repeatStart, MAX_REPEAT_POS_SIZE, false));
		} else {
			buf.append(formatNumber(repeatStart, MAX_REPEAT_POS_SIZE, false))
				.append(" ")	
				.append(formatNumber(repeatEnd, MAX_REPEAT_POS_SIZE, false))
				.append(" ")
				.append(formatNumber(repeatBasesLeft, MAX_REPEAT_POS_SIZE, true));
		}
		buf.append(" ")
			.append(rightJustify(id, 4))
			.append(" ");
		if(isRemarkable) {
			buf.append("*");
		}
		
		return buf.toString();
	}
	
	public static String formatNumber(long position, int cellSize, boolean inParen) {
		String positionStr = String.valueOf(position);
		if(inParen) {
			positionStr = "(" + positionStr + ")";
		}

		return rightJustify(positionStr, cellSize);
	}
	
	public static String formatNumber(double number, int cellSize, boolean inParen) {
		return formatNumber((int) number, cellSize, inParen);
	}
	
	public static String leftJustify(String name, int maxSize) {
		String newName = name.length() > maxSize ? name.substring(0,maxSize - 1) : name;
		StringBuffer buf = new StringBuffer(newName);
		for(int i = 0; i < maxSize - newName.length(); i++) {
			buf.append(" ");
		}
		return buf.toString();
	}
	
	public static String rightJustify(String name, int maxSize) {
		StringBuffer buf = new StringBuffer();
		String newName = name.length() > maxSize ? name.substring(0, maxSize - 1) : name;
		for(int i = 0; i < maxSize - newName.length(); i++) {
			buf.append(" ");
		}
		buf.append(newName);
		return buf.toString();
	}
	

	public float getPercentDeletions() {
		return percentDeletions;
	}

	public float getPercentDivergence() {
		return percentDivergence;
	}

	public float getPercentInsertions() {
		return percentInsertions;
	}

	public String getQuery() {
		return query;
	}
	
	protected void setQuery(String query) {
		this.query = query.intern();
	}

	public int getQueryBasesLeft() {
		return queryBasesLeft;
	}

	public int getRepeatEnd() {
		return repeatEnd;
	}

	public String getRepeatFamily() {
		return repeatFamily;
	}
	
	public String getGeneralFamily() {
		return repeatFamily.replaceAll("/.+","");
	}

	public String getRepeatName() {
		return repeatName;
	}

	public int getRepeatStart() {
		return repeatStart;
	}

	public String getName() {
		return getRepeatName();
	}
	
	public String getId() {
		return id;
	}

}
