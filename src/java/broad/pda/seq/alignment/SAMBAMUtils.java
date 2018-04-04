package broad.pda.seq.alignment;

import java.util.Map;

import broad.pda.annotation.BEDFileParser;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

public class SAMBAMUtils {
	public static void createSequenceDictionary(SAMFileReader reader, String sizeFile) {
		Map<String, Integer> sizes = BEDFileParser.loadChrSizes(sizeFile);
		SAMFileHeader header = reader.getFileHeader();
		SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
		for(String key : sizes.keySet()) {
			dictionary.addSequence(new SAMSequenceRecord(key, sizes.get(key)));
		}
		header.setSequenceDictionary(dictionary);
	}
}
