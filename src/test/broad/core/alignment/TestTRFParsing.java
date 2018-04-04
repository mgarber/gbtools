package broad.core.alignment;

import java.io.FileInputStream;
import java.net.URL;
import java.util.List;

import junit.framework.TestCase;

public class TestTRFParsing extends TestCase {
	static final String XIST_TRF_OUT = "XIST_hg18.fa.2.7.7.80.10.50.2000.dat";
	
	public void testParse() throws Exception {
		URL trfDataUrl = getClass().getResource(XIST_TRF_OUT);
		TRFReader reader = new TRFReader();
		FileInputStream fis = new FileInputStream(trfDataUrl.getFile());
		reader.load(fis);
		fis.close();
		List<TRFHit> TSIXTandemRepeats = reader.getSequenceTandemRepeats("NR_003255");
		assertEquals("Wrong number of TSIX repeats", 8, TSIXTandemRepeats.size());
		
		List<TRFHit> XISTTandemRepeats = reader.getSequenceTandemRepeats("NR_001564");
		assertEquals("Wrong number of distinct XIST repeats",15, XISTTandemRepeats.size());
	}
}
