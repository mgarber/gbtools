package umms.glab.nio;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.net.URL;
import java.nio.channels.FileChannel;

import org.junit.Test;

import junit.framework.TestCase;

public class FileChannelBufferedReaderTest extends TestCase {
	//static final String FQ_TEST = "test.fq"; //AARG it is not working. lets hack in the meantime
	static final String FQ_TEST = "/Users/mgarber/Dropbox/eclipsews/gbtools/bin/umms/glab/nio/test.fq";

	@Test
	public void testBasicFunction() throws IOException {
		URL fqURL = getClass().getResource(FQ_TEST);
		//RandomAccessFile fqRAF = new RandomAccessFile(fqURL.getFile(), "r");
		RandomAccessFile fqRAF = new RandomAccessFile(FQ_TEST, "r");
		FileChannel fqChannel = fqRAF.getChannel();
		FileChannelBufferedReader fqChannelBR = new FileChannelBufferedReader(fqChannel);
		fqChannelBR.init();
		
		//BufferedReader br = new BufferedReader(new FileReader(fqURL.getFile()));
		BufferedReader br = new BufferedReader(new FileReader(FQ_TEST));
		
		String brLine;
		String fcbrLine;
		while( (brLine = br.readLine() ) != null) {
			fcbrLine = fqChannelBR.readLine();
			assertEquals("Lines did not match: br: "+brLine + " fcbr: " + fcbrLine, brLine, fcbrLine);
		}
		
		fqRAF.close();
		br.close();
	}
	
	public void testBasicFunctionSmallBuffer() throws IOException {
		URL fqURL = getClass().getResource(FQ_TEST);
		//RandomAccessFile fqRAF = new RandomAccessFile(fqURL.getFile(), "r");
		RandomAccessFile fqRAF = new RandomAccessFile(FQ_TEST, "r");
		FileChannel fqChannel = fqRAF.getChannel();
		FileChannelBufferedReader fqChannelBR = new FileChannelBufferedReader(fqChannel, 2); //This is the key of this test. Only 5 bytes in the buffer
		fqChannelBR.init();
		
		//BufferedReader br = new BufferedReader(new FileReader(fqURL.getFile()));
		BufferedReader br = new BufferedReader(new FileReader(FQ_TEST));
		
		String brLine;
		String fcbrLine;
		while( (brLine = br.readLine() ) != null) {
			fcbrLine = fqChannelBR.readLine();
			assertEquals("Lines did not match: br: "+brLine + " fcbr: " + fcbrLine, brLine, fcbrLine);
		}
		
		fqRAF.close();
		br.close();
	}

}
