package broad.pda.seq.alignment;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMRecord;

public class RecordWriter {
	BufferedWriter writer;	
	SAMFileWriter samWriter;
	boolean writeFullBed;
	int minimumMappingScore = 0;
	public RecordWriter(String outputFile, SAMFileHeader header, String applicationName) throws IOException {

		if(outputFile.endsWith(".sam")|| outputFile.endsWith(".bam")) {
			header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
			header.addProgramRecord(new SAMProgramRecord(applicationName));
			SAMFileWriterFactory factory = new SAMFileWriterFactory();
			samWriter = factory.makeSAMOrBAMWriter(header, false, new File(outputFile));
		} else if(outputFile.endsWith(".aligned") || outputFile.endsWith(".bed")){
			writer = new BufferedWriter(new FileWriter(outputFile));
			writeFullBed =  outputFile.endsWith(".bed");
		} else if (outputFile.endsWith(".bed")){

		}else {
			throw new IllegalArgumentException("Output file type not supported, the extension of paired end output file can be only one of .aligned, .sam, .bam or .bed");
		}
	}

	public void setMinimumMappingScore(int minimumMappingScore) {
		this.minimumMappingScore = minimumMappingScore;

	}

	public void addAlignment(SAMRecord record) throws IOException {
		if(samWriter != null) {
			samWriter.addAlignment(record);
		} else if (writer != null && record.getMappingQuality() > minimumMappingScore) {

			writer.write(record.getReferenceName());
			writer.write("\t");
			writer.write(String.valueOf(record.getAlignmentStart()));
			writer.write("\t");
			writer.write(String.valueOf(record.getAlignmentEnd()));
			if(writeFullBed) {
				writer.write("\t");
				writer.write(record.getReadName());
				writer.write("\t");
				writer.write(String.valueOf(record.getMappingQuality()));
				writer.write("\t");
				writer.write(record.getReadNegativeStrandFlag() ? "-" : "+");
			}
			writer.newLine();
		}

	}

	public void close() throws IOException {
		if(samWriter != null ) {
			samWriter.close();
		}

		if(writer != null) {
			writer.close();
		}
	}



}
