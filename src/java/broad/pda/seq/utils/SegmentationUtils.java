package broad.pda.seq.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import org.broad.igv.sam.Alignment;

import broad.core.annotation.AnnotationReader;
import broad.core.annotation.AnnotationReaderFactory;
import broad.core.annotation.BED;
import broad.core.annotation.GenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.error.ParseException;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

public class SegmentationUtils {
	static final int chunkSize=50000000;
	
	public static final String USAGE = "Usage: SegmentationUtils TASK=<task name> <task args>\n" +
	"\tTasks:\n" +
	"\t\tcall_high_peaks Takes a list of regions and an alignment file reports back the highest peak within each region. \n\t\t\t-regions <Region file> \n\t\t\t-alignment <alignment file>  \n\t\t\t-maskFileDir <Mask File directory> \n\t\t\t-sizeFile <Chromosome size file>" +
	"\t\ttotalNumberOfReads Takes an alignment file and reports back the number of aligned reads in the file. \n\t\t\t-alignment <alignment file>  \n\t\t\t-name <name of set> \n\t\t\t-sizeFile <Chromosome size file>" +
	"\n";
	public static void main (String [] args) throws IOException, ParseException {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE, "call_high_peaks");
		
		if ("call_high_peaks".equals(argMap.getTask())) {
			String regionFormat = argMap.containsKey("regionFormat") ? argMap.get("regionFormat") : "BED";
			String regionFile = argMap.getMandatory("regions");
			String alignmentFile = argMap.getMandatory("alignment");
			File[] maskFiles  = new File(argMap.getMandatory("maskFileDir")).listFiles();
			String sizes = argMap.getMandatory("sizeFile");

			AlignmentDataModel data =new GenericAlignmentDataModel(alignmentFile, sizes);
			System.err.println("Loaded data");
			AnnotationReader<? extends GenomicAnnotation> regions = AnnotationReaderFactory.create(regionFile, regionFormat);
			System.err.println("Loaded regions");
			
			Iterator<String> chrIt = regions.getChromosomeIterator();
			Map<String, Collection<BED>> peaks = new LinkedHashMap<String, Collection<BED>>();
			while(chrIt.hasNext()) {
				String chr = chrIt.next();
				System.err.println("Processing " + chr);
				int chunkNumber=1;
				IntervalTree<Alignment> chunkAlignmentTree=data.getIntervalTree(chr, 0, chunkNumber*chunkSize);
				TreeSet<BED> chrPeaks = new TreeSet<BED>();
				peaks.put(chr, chrPeaks);
				List<? extends GenomicAnnotation> chrRegions = regions.getChromosomeBEDs(chr);
				
				for(GenomicAnnotation chrRegion : chrRegions) {
					//System.err.println("Processing region " + chrRegion);
					if(chrRegion.length() == 0) {
						System.err.println("WARNING: found zero length segment: " + chrRegion + " IGNORING");
						continue;
					}
					double st = System.nanoTime();
					Alignments highestPeak = null;
					for(int p  = chrRegion.getStart(); p < chrRegion.getEnd(); p++) {
						Alignments peak = new Alignments(chr, p, p+1);
						peak.setOrientation(chrRegion.getOrientation());
						if(p+1>=(chunkNumber*chunkSize-1)){
							chunkNumber++; 
							chunkAlignmentTree=data.getIntervalTree(chr, p, Math.max(chunkNumber*chunkSize, p+1));
						}
						double score=data.getCountsPerAlignment(peak, chunkAlignmentTree, 0);
						

						if(highestPeak == null || highestPeak.getScore() < score) {
							highestPeak = peak;
							highestPeak.addScore(score);
						} else if (highestPeak.getScore() == score && highestPeak.getEnd() == peak.getStart()) {
							highestPeak.setEnd(peak.getEnd());
						}
					}

					chrPeaks.add(new BED(highestPeak));
					//System.err.println("finished region, took " + (System.nanoTime() - st));
				}
			
			}
			BufferedWriter bw = argMap.getOutputWriter();
			write(bw, peaks);
			bw.close();
		}
		else if ("totalNumberOfReads".equals(argMap.getTask())){
			String alignmentFile = argMap.getMandatory("alignment");
			String sizes = argMap.getMandatory("sizeFile");
			String name=argMap.getMandatory("name");
			//This Only CALCULATES THE NUMBER OF MAPPED READS- WE WANT ALSO THE NUMBER OF READS:
			//GenericAlignmentDataModel data =new GenericAlignmentDataModel(alignmentFile, sizes);
			//double[]  total = data.getTotalNumberOfMappedReads();
			
			double[]  total = getTotalNumberOfMappedReads(alignmentFile);
			
			System.out.println(name+"\t"+total[0]+"\t"+total[1]+"\t"+total[2]);
			
		}
		else{
				System.err.println("Invalid task " + argMap.getTask() + "\n" + USAGE);
		}
	}
	
	
	public static void write(BufferedWriter bw, Map<String, Collection<BED>> peaks) throws IOException {
		for(String chr : peaks.keySet()) {
			for(BED peak : peaks.get(chr)) {
				bw.write(peak.toShortString());
				bw.newLine();
			}
		}
	}
	
	//Returns an array with: [total num of reads] [total mappings] [total of unique reads that were mapped]
	private static double[] getTotalNumberOfMappedReads(String alignmentFile) {
		
		final SAMFileReader inputSam = new SAMFileReader (new File(alignmentFile));
    	Iterator <SAMRecord> readIter =inputSam.iterator();
    	
    	double numReads=0;
		double numMapped=0;
		
		Map <String,Integer> fragMap=new HashMap<String,Integer>();
		
    	while(readIter.hasNext()){
    		SAMRecord sam=readIter.next();
    		numReads++;
			if (!sam.getReadUnmappedFlag()) {
				numMapped++;
				if (! fragMap.containsKey(sam.getReadName()))
						fragMap.put(sam.getReadName(), 0);
				fragMap.put(sam.getReadName(), fragMap.get(sam.getReadName())+1);
				
				}
    		
    	}
		
    	inputSam.close();
		double [] res= new double[3];
		res[0]=numReads;
		res[1]=numMapped; 
		res[2]=fragMap.size(); //num of unique map
		return res;
		
		
	}


}
