/**
 * 
 */
package broad.pda.rap;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

//import broad.core.annotation.BED;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

/**
 * @author engreitz
 *
 */
public class SamToBed {

	public SamToBed(File in, BufferedWriter bw, int extend) throws IOException {
		final SAMFileReader reader = new SAMFileReader(in);
        final CloseableIterator<SAMRecord> iterator = reader.iterator();
        while (iterator.hasNext()) {
            final SAMRecord record = iterator.next();
            
         // BED format is 0-based [), whereas SAM is 1-based []
            int start = record.getAlignmentStart() - 1;  
            int end = record.getAlignmentEnd();
            
			if (record.getReadNegativeStrandFlag()) { 
				start = Math.max(end - extend, 0); 
			} else { 
				end = Math.min(start + extend, record.getHeader().getSequence(record.getReferenceName()).getSequenceLength()); 
			}
			bw.write(record.getReferenceName() + "\t" + start + "\t" + end + "\n");   
        }
        iterator.close();
        bw.close();
	}
	
	private static String USAGE = "java -jar SamToBed -in <in.sam> -out <out.bed> [-extendTo <# of bases to extend each read to]\n" +
								  "\tTODO: support paired-end read format\n";
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE);
		
		// TODO support paired-end read format
		File in = new File(argMap.getInput());
		BufferedWriter bw = argMap.getOutputWriter();
		int extend = (argMap.containsKey("extendTo")) ? argMap.getInteger("extendTo") : -1;
		new SamToBed(in, bw, extend);
	}

}
