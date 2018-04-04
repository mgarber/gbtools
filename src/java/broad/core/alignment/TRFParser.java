package broad.core.alignment;

import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import broad.core.error.ParseException;

public class TRFParser {
	private static final Pattern hitLineStart = Pattern.compile("^[0-9]+\\s[0-9]+\\s");
	
	
	public void parse(java.io.BufferedReader reader, TRFHandler handler) throws  IOException, ParseException {
		String line = null;
		String currentSequence = null;
		while((line = reader.readLine()) != null) {
			if(line.startsWith("Sequence: ")) {
				currentSequence = line.replaceAll("Sequence: ", "");
				handler.sequence(currentSequence);
			}
			Matcher matcher = hitLineStart.matcher(line);
			if(matcher.find()) {
				handler.tandemRepeat(new TRFHit(currentSequence, line.split(" ")));
			}
		}
	}

}
