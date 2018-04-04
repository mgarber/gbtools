package broad.core.alignment;

import broad.core.error.ParseException;

public interface TRFHandler {
	void sequence(String sequenceName);
	void tandemRepeat(TRFHit tandemRepeat) throws ParseException;

}
