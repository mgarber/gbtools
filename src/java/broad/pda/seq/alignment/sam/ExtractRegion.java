package broad.pda.seq.alignment.sam;

import java.io.File;

import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

public class ExtractRegion {
	public static String USAGE = "Extracts all alignment records within specified region, assumes the alignment is sorted and indexed. \n\t-in <Alignment file in bam or sam format, extension must be .bam or .sam> \n\t-chr <Chromosome> \n\t-start<Start of region>\n\t-end<End of region> -out <Output extension needs to be .sam or .bam> \n";
	
	public static void main(String [] args) throws Exception {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE, "extract");
		
		File inputFile = new File(argMap.getInput());
		File outFile   = new File(argMap.getOutput());
		
		String chr = argMap.getMandatory("chr");
		int start = argMap.getInteger("start");
		int end = argMap.getInteger("end");
		final SAMFileReader inputSam = new SAMFileReader(inputFile);
		inputSam.setValidationStringency(ValidationStringency.SILENT);
		CloseableIterator<SAMRecord> readIt =  inputSam.query(chr, start, end, false);
		
		final SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(inputSam.getFileHeader(), true, outFile);

		while(readIt.hasNext()) {
			SAMRecord read = readIt.next();
			outputSam.addAlignment(read);
		}

		readIt.close();
        inputSam.close();
        outputSam.close();
		
	}
}
