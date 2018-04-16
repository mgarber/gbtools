package umms.glab.nio;

import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;

import org.apache.log4j.Logger;

public class FileChannelBufferedReader {
	static Logger logger = Logger.getLogger(FileChannelBufferedReader.class.getName());
	
	public static final int DEFALUT_BUFFER_SIZE = 65536;
	private int bufferSize;
	private FileChannel fc;
	private ByteBuffer buffer;
	private int currentBufferLine;
	private String [] bufferLines;
	private long previousLinePosition;
	private long position;

	
	public FileChannelBufferedReader(FileChannel fc) {
		this.fc = fc;
		this.bufferSize = DEFALUT_BUFFER_SIZE;
	}
	
	public FileChannelBufferedReader(FileChannel fc, int bufferSize) {
		this.fc = fc;
		this.bufferSize = bufferSize;
	}
	
	public void init() throws IOException {
		position = fc.position();
		previousLinePosition = 0;
		this.buffer = ByteBuffer.allocate(bufferSize);
		readBuffer();

		
	}
	
	public String readLine() throws UnsupportedEncodingException, IOException {
		String line = null;
		boolean atEOF = false;
		if(currentBufferLine >= bufferLines.length - 1) { //If we are at the last line in the buffer
			atEOF = readBuffer() == -1;
		} 
		
		if(!atEOF & bufferLines != null && bufferLines.length > 0) {
			line = bufferLines[currentBufferLine];
			currentBufferLine++;
			previousLinePosition = position;
			position = position + 1 + line.getBytes().length; // I am still not sure why the new line is just worth 1 and not 2 bytes.
		}
		
		return line;
	}

	private int readBuffer() throws IOException, UnsupportedEncodingException {
		buffer.clear();
		fc.position(position);
		int retVal = 0;
		StringBuffer lineStringBuffer = new StringBuffer(bufferSize);
		currentBufferLine = 0;
		while (retVal != -1 && !lineStringBuffer.toString().contains("\n") )  {//Again, this should be \R, but it is unclear why \R does not compile
			retVal = fc.read(buffer);
			buffer.flip();
			lineStringBuffer.append(new String(buffer.array(), "ASCII"));
		}
			//System.out.println("Called readBuffer, leastLine " + (bufferLines != null && bufferLines.length > 0 ? bufferLines[bufferLines.length - 1] : "NONE"));
		bufferLines = lineStringBuffer.toString().split("\n"); //It should be \R but somehow it is not compiling. \R is supposed to be available in Java 8

		//System.out.println("Line after reloading: " +  (bufferLines.length > 0 ? bufferLines[0] : "NONE"));
		if(bufferLines.length > 0) {
			logger.debug("FileChannelBufferedReader initialized. Buffer loaded, Number of lines: " + bufferLines.length + " first: " + bufferLines[0]);
		} else {
			logger.info("Tried to initialize FileChannelBufferedReader, but nothing was read in while initializing buffer");
		}
		
		return retVal;
	}
	
	public long getPosition () { return position;}
	public long getPreviousLinePosition() { return previousLinePosition;}

}
