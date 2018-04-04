package broad.pda.seq.alignment.sam;

import java.io.BufferedWriter;
import java.io.File;

import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;

public class ExtractHeader {
	public static String USAGE = "Extracts the header of a SAM or BAM file. \n\t-in <Alignment file in bam or sam format, extension must be .bam or .sam> -out <If specified header will be written to this file otherwise it is written to standard out> \n";
	
	public static void main(String [] args) throws Exception {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE, "extract");
		
		File inputFile = new File(argMap.getInput());

		final SAMFileReader inputSam = new SAMFileReader(inputFile);
		inputSam.setValidationStringency(ValidationStringency.SILENT);
		SAMFileHeader header = inputSam.getFileHeader();
		BufferedWriter bw = argMap.getOutputWriter();
		bw.write("Sort order: " + header.getSortOrder().toString());
		bw.close();		
	}
}
