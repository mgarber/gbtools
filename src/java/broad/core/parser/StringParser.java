/**
 * 
 */
package broad.core.parser;

/**
 * @author prussell
 * Parse a string around whitespace or another specified delimiter
 * Get full array of parsed strings or individual values by index, parsed to various types
 */
public class StringParser {

	public StringParser() {this.tokens = null;}
	
	/**
	 * Parses the string around whitespace and stores tokens
	 * @param s The string to parse
	 */
	public void parse(String s) {
		this.str = s;
		if(this.str == null) return;
		if(this.str == "") return;
		this.tokens = s.split(whitespaceDelimiter);
	}
	
	/**
	 * Parses the string around whitespace and stores tokens
	 * @param s The string to parse
	 */
	public static String[] getTokens(String s) {
		return s.split(whitespaceDelimiter);
	}

	
	/**
	 * Parses the string around the specified delimiter and stores tokens
	 * @param s The string to parse
	 * @param delimiter The delimiter
	 */
	public void parse(String s, String delimiter) {
		this.tokens = s.split(delimiter);
	}
	
	/**
	 * Get number of fields
	 * @return number of fields
	 */
	public int getFieldCount() {
		if(this.tokens == null) return 0;
		return this.tokens.length;
	}
	
	/**
	 * Gets the token at specified position
	 * @param index The position
	 * @return the desired token as a String
	 */
	public String asString(int index) {
		return this.tokens[index];
	}
	
	/**
	 * Gets the token at specified position and parses to an int
	 * @param index The position
	 * @return the desired token as an int
	 */
	public int asInt(int index) {
		try {
			return Integer.parseInt(this.tokens[index]);
		} catch (NumberFormatException e) {
			throw new NumberFormatException("Field " + index + " cannot be parsed to int in string: " + this.str);
		}
	}
	
	/**
	 * Gets the token at specified position and parses to a double
	 * @param index The position
	 * @return the desired token as a double
	 */
	public double asDouble(int index) {
		try {
			return Double.parseDouble(this.tokens[index]);
		} catch (NumberFormatException e) {
			throw new NumberFormatException("Field " + index + " cannot be parsed to double in string: " + this.str);
		}
	}
	
	/**
	 * Gets the token at specified position and parses to a float
	 * @param index The position
	 * @return the desired token as a float
	 */
	public float asFloat(int index) {
		try {
			return Float.parseFloat(this.tokens[index]);
		} catch (NumberFormatException e) {
			throw new NumberFormatException("Field " + index + " cannot be parsed to float in string: " + this.str);
		}
	}
	
	/**
	 * Gets the array of tokens
	 * @return all the parsed tokens as a String[]
	 */
	public String[] getStringArray() {
		return this.tokens;
	}
	
	private static String whitespaceDelimiter = "\\s++"; //$NON-NLS-1$
	private String[] tokens;
	private String str;
	
}
