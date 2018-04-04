package broad.pda.seq.alignment;

import java.io.IOException;

import broad.pda.annotation.BEDFileParser;

public class PSLToBED {

	public static void main(String[] args) throws IOException{
		BEDFileParser.convertPSLToBED(args[0], args[1]);
	}
	
}
