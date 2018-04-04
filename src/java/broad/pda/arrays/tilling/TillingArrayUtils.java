package broad.pda.arrays.tilling;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.annotation.AnnotationFactoryFactory;
import broad.core.annotation.AnnotationHandler;
import broad.core.annotation.AnnotationReader;
import broad.core.annotation.AnnotationReaderFactory;
import broad.core.annotation.BED;
import broad.core.annotation.BEDReader;
import broad.core.annotation.BasicAnnotationReader;
import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.GenomicAnnotationFilter;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.datastructures.ContinuousData;
import broad.pda.datastructures.LocationAwareMatrix;
import broad.pda.gene.RefSeqGene;
import broad.pda.gene.RefSeqGeneWithIsoforms;

public class TillingArrayUtils {
	public static final double DEFAULT_FLOOR = 200.0;
	public static final String USAGE = "Usage: TillingArrayUtils TASK=<task_num> <task_args>\n" +
	"\tTasks:\n" +
	"\n\t1. Normalize data (in GFF format or other annotation formats) \n\t\t-in <Input file, standard input is not supported> \n\t\t-out <output file or standard output if none is specified> \n\t\tOptional: \n\t\t\t-indir <If instead of a single file a set of files should be normalized together (need to provide -ext)> \n\t\t\t-ext <comma separated list of extensions to load> ] \n\t\t\t-gp <for GCT output instead of IGV> [\n\t\t\t-format <Input data format [GFF], BED, Generic> \n\t\t\t-floor <Minimum observation value, all observation less than floor are set to floor> -qn <Include this flag for quantile normalization rather than sample mean centering>  \n\t\t\t -noNorm  <do not apply any normalization on the data>] -noMeanNorm  <processes the data without row mean normalization>] " +
	"\n\t2. Segment a previously normalized set (or just one) of expression experiments \n\t\t-in <Normalized experiment table> \n\t\t-outdir <If instead of a single sample a set of samples are segmented at once> \n\t\tOptional:\n\t\t\t-colNames <Comma separated list of column names to segment if a subset rather than all columns are to be segmented> \n\t\t\t-chr <If only a data for a specific chromosome should be segmented> \n\t\t\t-windows <Comma separated list of window sizes (in consecutive probes) to run default (5,7,8,9,10,15,20)> \n\t\t\t-twoTailed <Default is one tailed test> \n\t\t\t-mergeDistance <Maximum base distance to merge probes into contiguous segment [1000]> \n\t\t\t-randomizations <How many randomizations to perform, the more the more accurate pvalues become [1000] \n\t\t\t-alpha <max FWR corrected pvalue [0.1]> " +
	"\n\t3. Find differentially expressed regions between two groups of samples \n\t\t-in <Normalized experiment table, standard input is assumed if none is provided> \n\t\t-out <Output file with segments exceeding FWER corrected pvalue output goes to standard output if no out is provided> \n\t\t-group1 <Comma separated list of column names corresponding to the first group> -group2 <Comma separated list of column names corresponding to the second group> \n\t\tOptional:\n\t\t\t-chr <If only a data for a specific chromosome should be segmented>  \n\t\t\t-mergeDistance <Maximum base distance to merge probes into contiguous segment [1000]> \n\t\t\t-randomizations <How many randomizations to perform, the more the more accurate pvalues become [1000] \n\t\t\t-alpha <max FWER corrected pvalue [0.1]> -filterByFDR <use alpha to control FDR rather than FWER> \n\t\t\t-regions <File containing regions to restric the analysis too, if none is suppied, then all segments are used> -regionFileFormat <[BED], GFF, GENERIC> " +
	"\n\t4. Score contiguous regions from normalized data \n\t\t-in <Normalized experiment table, standard input is assumed if none is provided> \n\t\t-out <Mean scored contiguous regions> \n\t\t\t-mergeDistance <Maximum base distance to merge probes into contiguous segment [1000]>\n\t\t\t-regions <File containing regions to restric the analysis too, if none is suppied, then all segments are used> -regionFileFormat <[BED], GFF, GENERIC> " +
	"\n\t5. Substract columns \n\t\t-in <Normalized experiment table, standard input is assumed if none is provided> \n\t\t-out <Matrix with column differences> \n\t\t\tcolumns1 <First set of columns (comma separated list of column names)> \n\t\t\tcolumns2 <comma separated list of column names, must be of same size of column1>" +
	"\n\t6. Shuffle a previously normalized set \n\t\t-in <Normalized experiment table> \n\t\t-out <output file or standard output if none specified>"+
	"\n\t7. Generate shuffle value and extreme value distributions \n\t\t-in <Normalized experiment table> \n\t\t-outdir <output directory to write data> -\n\t\t-colName <Column to use> -\n\t\t-permutations <How many permutations to do for extreme value distribution> \n\t\t-window <Window to scan (in consecutve probes)> \n\t\t-chr <Chromosome to load> \n\t\t-mergedDistance <Maximum base distance to merge probes into contiguous segment [100]>"+
	"\n\t8. Filter a location aware matrix with a given set of annotations (e.g. Repeats) \n\t\t-regions <Region file (Only BED format supported right now)> \n\t\t-filterRegions <Regions to filter> \n\t\t-maxPctInFilterRegions <Minimum % of locations in regions>" +
	"\n\t9. Load raw data from different files to one matrix (in GFF format or other annotation formats) \n\t\t-in <Input file, standard input is not supported> \n\t\t-out <output file or standard output if none is specified> \n\t\tOptional: \n\t\t\t-indir <If instead of a single file a set of files should be normalized together (need to provide -ext)> \n\t\t\t-ext <comma separated list of extensions to load> ] \n\t\t\t-gp <for GCT output instead of IGV> [\n\t\t\t-format <Input data format [GFF], BED, Generic> ]  " +
	"\n\t10. Map lincs to their overlapping segments/regions. report 2 files: 1. the list of segments 2. mapping between linc to segment  \n\t\t-in <Input file, bed format> \n\t\t-out <output Prefix> \n\t\t -regions <regions file>" +
	"\n";

	public static void main(String [] args) throws ParseException, IOException {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE, "1");
		TillingArrayUtils tau = new TillingArrayUtils();
		if("1".equals(argMap.getTask()) || argMap.getTask() == null) {
			List<File> dataFiles = tau.getWorkFiles(argMap);
			double floor = argMap.containsKey("floor") ? argMap.getDouble("floor") : DEFAULT_FLOOR;
			
			if(dataFiles.isEmpty() ) {
				return;
			}
			String format = argMap.containsKey("format") ? argMap.get("format") : "GFF";
			File first = dataFiles.get(0);
			AnnotationReader<? extends GenomicAnnotation> reader = AnnotationReaderFactory.create(format); 
			LocationAwareMatrix matrix = new LocationAwareMatrix( );
			int rowsRead = reader.parse(first.getAbsolutePath(),  null, matrix);
		
			reader = null; //Help the GC by explicitly telling it now that the reader will no longer be referenced.
			List<String> colNames = new ArrayList<String>(dataFiles.size());
			for(File file : dataFiles) {
				colNames.add(CLUtil.replaceFileExtension(file.getName(),""));
			}
			matrix.setColumns(colNames);
			matrix.setDataDimensions(rowsRead, colNames.size());

			for(int i = 0; i < dataFiles.size(); i++) {
				tau.addData(dataFiles.get(i), matrix, format, colNames.get(i), floor);
			}

			if (argMap.containsKey("noNorm"))
			{
				//do not do anything- ADDED BY MORAN
				System.err.println("Data is not normalized \n");
				
			}
			else {
				if(argMap.containsKey("qn")) {
			
				System.err.println("Quanlile normalizing, may take a while");
				matrix.quantileNormalizeColumns();	
				}
				if (! argMap.containsKey("noMeanNorm")) matrix.normalizeColumnsByMean();
				matrix.log();
		   
				
			}
			
			
			BufferedWriter bw = argMap.getOutputWriter();
			if(argMap.containsKey("gp")) {
				matrix.writeGCT(bw);
			} else {
				matrix.write(bw);
			}
			bw.close();
			
		} else if ("2".equals(argMap.getTask())) {
			String outdir = argMap.getOutputDir();
			String chr = null;
			if(argMap.containsKey("chr")) {
				chr = argMap.getMandatory("chr");
			}
			boolean twoTailed = argMap.containsKey("twoTailed");
			double alpha = argMap.containsKey("alpha") ? argMap.getDouble("alpha") : 0.1;
			int randomizations = argMap.containsKey("randomizations") ? argMap.getInteger("randomizations") : 1000;
			int mergeDistance = argMap.containsKey("mergeDistance") ? argMap.getInteger("mergeDistance") : 1000;
			int[] windowSizes = getWindowSizes(argMap);
			BufferedReader br = argMap.getInputReader();
			BasicAnnotationReader combinedSegments = new BasicAnnotationReader();
			try {
				LocationAwareMatrix matrix = new LocationAwareMatrix(br, chr);
				List<String> cols = new ArrayList<String>();
				if(argMap.containsKey("colNames")) {
					String [] colNames = argMap.getMandatory("colNames").split(",");
					for(String colName : colNames) {
						if (!matrix.hasColumn(colName)) {
							throw new IllegalArgumentException("Column " + colName + " is not in input data");
						}
						cols.add(colName);
					}
				} else {
					cols = matrix.getColumnNames();
				}
				
				for(String col : cols) {
					System.err.println("Processing sample " + col);
					ContinuousData data = new ContinuousData(matrix, col);
					WindowPermuteSegmenter wps = new WindowPermuteSegmenter(data, mergeDistance);
					Map<Alignments,Double> allSegments = new TreeMap<Alignments, Double>();
					
					if(chr != null) {
						allSegments = wps.segment(twoTailed ? 2 : 1, windowSizes, randomizations, 1 - alpha, chr);						
					} else {
						Map<String, Map<Alignments, Double>> segments = wps.segment(twoTailed ? 2 : 1, windowSizes, randomizations, 1 - alpha);
						allSegments = new TreeMap<Alignments, Double>();
						for(String chromosome : segments.keySet()) {							
							allSegments.putAll(segments.get(chromosome));
						}
					}
					String outfile = outdir + "/" + col  + ".segments";
					System.err.println("Printing to " + outfile);
					wps.writeMerged(allSegments, outfile);
					for(Alignments a: allSegments.keySet()) {
						
						combinedSegments.addAnnotation(new BasicGenomicAnnotation(a));
					}
					
				}
			} finally {
				if(br != null) {
					br.close();
				}
			}
			combinedSegments.merge();
			List<BasicGenomicAnnotation> combinedSegmentList = combinedSegments.getAnnotationList();
			String outfile = outdir + "/all.segments";
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
			for(BasicGenomicAnnotation bga : combinedSegmentList) {
				bw.write(new BED(bga).toString());
				bw.newLine();
			} 
			bw.close();
		} else if ("3".equals(argMap.getTask())) {
			String outdir = argMap.getOutputDir();
			String chr = argMap.get("chr");
			if(chr!= null) {
				chr = chr.replace("chr","");
			}
			double alpha = argMap.containsKey("alpha") ? argMap.getDouble("alpha") : 0.1;
			int randomizations = argMap.containsKey("randomizations") ? argMap.getInteger("randomizations") : 1000;
			int mergeDistance = argMap.containsKey("mergeDistance") ? argMap.getInteger("mergeDistance") : 100;
			int[] windowSizes = getWindowSizes(argMap);
			String [] group1Columns = argMap.getMandatory("group1").split(",");
			String [] group2Columns = argMap.getMandatory("group2").split(",");
			
			List<? extends LightweightGenomicAnnotation> regions = null;
			if(argMap.containsKey("regions")) {
				String regionFile = argMap.getMandatory("regions");
				String format = argMap.containsKey("regionFormat")? argMap.get("regionFormat") : "BED";
				AnnotationReader<? extends LightweightGenomicAnnotation> reader = AnnotationReaderFactory.create(regionFile, format);
				regions = (List<? extends LightweightGenomicAnnotation>) reader.getAnnotationList();
				
			}
			
			BufferedReader br = argMap.getInputReader();
			try {
				LocationAwareMatrix matrix = new LocationAwareMatrix(br, chr);

				WindowPermuteDifferentialExpression wpde = 
					new WindowPermuteDifferentialExpression(matrix, group1Columns, group2Columns, chr, mergeDistance, regions);
				MatrixWithHeaders result = wpde.ttest(alpha, randomizations, argMap.containsKey("filterByFDR"));
				BufferedWriter bw = argMap.getOutputWriter();
				if(argMap.containsKey("gp")) {
					result.writeGCT(bw);
				} else {
					result.write(bw);
				}
				bw.close();

			} finally {
				if(br != null) {
					br.close();
				}
			}
		} else if ("4".equals(argMap.getTask())) {
			int mergeDistance = argMap.containsKey("mergeDistance") ? argMap.getInteger("mergeDistance") : 100;
			
			List<? extends LightweightGenomicAnnotation> regions = null;
			if(argMap.containsKey("regions")) {
				String regionFile = argMap.getMandatory("regions");
				String format = argMap.containsKey("regionFormat")? argMap.get("regionFormat") : "BED";
				AnnotationReader<? extends LightweightGenomicAnnotation> reader = AnnotationReaderFactory.create(regionFile, format);
				regions = (List<? extends LightweightGenomicAnnotation>) reader.getAnnotationList();
				
			}
			
			BufferedReader br = argMap.getInputReader();
			
			try {
				LocationAwareMatrix matrix = new LocationAwareMatrix(br, null);
				List<String> cols = matrix.getColumnNames();
				if(cols.size() == 0) {
					System.err.println("Matrix had no columns");
					return;
				}
				List<Alignments> contiguousRegions = new ArrayList<Alignments>();
				ContinuousData data = new ContinuousData(matrix, cols.get(0));

				for(String chromosome : data.getChromosomes()) {
					Set<Alignments> chrAlignments = regions != null ? 
							data.makeContiguousRegions(chromosome, regions) : 
								data.makeContiguousRegions(chromosome, mergeDistance);
					contiguousRegions.addAll(chrAlignments);
				}
				LocationAwareMatrix result = new LocationAwareMatrix(contiguousRegions, cols);
				
				for(int i = 0 ; i < cols.size(); i++) {
					if(i > 0) {
						data = new ContinuousData(matrix, cols.get(i));
					}
					System.err.println("Processing sample " + cols.get(i));
					for(Alignments region : contiguousRegions) {
						result.set(region.toUCSC(), i, Statistics.mean(data.getDataForAlignments(region)));
					}
					
				}
				
				BufferedWriter bw = argMap.getOutputWriter();
				if(argMap.containsKey("gp")) {
					result.writeGCT(bw);
				} else {
					result.write(bw);
				}
				bw.close();
			} finally {
				if(br != null) {
					br.close();
				}
			}
			

		} else if ("5".equals(argMap.getTask())) {
			String [] columns1 = argMap.getMandatory("columns1").split(",");
			String [] columns2 = argMap.getMandatory("columns2").split(",");
			
			if(columns1.length != columns2.length) {
				System.err.println("The sizes of both columns groups does not match: columns1 has " + columns1.length + " while columns2 " + columns2.length);
			}
			
				
			BufferedReader br = argMap.getInputReader();
			
			try {
				LocationAwareMatrix matrix = new LocationAwareMatrix(br, null);
				LocationAwareMatrix matrix1 = (LocationAwareMatrix) matrix.submatrixByColumnNames(columns1);
				LocationAwareMatrix matrix2 = (LocationAwareMatrix) matrix.submatrixByColumnNames(columns2);
				
				matrix1.minus(matrix2);
				
				BufferedWriter bw = argMap.getOutputWriter();
				if(argMap.containsKey("gp")) {
					matrix1.writeGCT(bw);
				} else {
					matrix1.write(bw);
				}
				bw.close();
			} finally {
				if(br != null) {
					br.close();
				}
			}
			

		} else if ("6".equals(argMap.getTask())) {
			BufferedReader br = argMap.getInputReader();
			try {
				LocationAwareMatrix matrix = new LocationAwareMatrix(br, null);
				List<String> cols = matrix.getColumnNames();
				LocationAwareMatrix shuffledMatrix = new LocationAwareMatrix(new ArrayList<LightweightGenomicAnnotation>(matrix.getRowAnnotations()), matrix.getColumnNames());
				for(String col : cols) {
					ContinuousData data = new ContinuousData(matrix, col);
					
					double [] shuffled = data.getShuffledData();
					for(int i = 0; i < shuffled.length; i++) {
						shuffledMatrix.set(i, col, shuffled[i]);
					}
					
				}
				BufferedWriter bw = argMap.getOutputWriter();
				shuffledMatrix.write(bw);
				bw.close();
			} finally {
				if(br != null) {
					br.close();
				}
			}
		} else if ("7".equals(argMap.getTask())) {
			String outdir = argMap.getOutputDir();
			String chr = argMap.get("chr");
			int randomizations = argMap.getInteger("permutations") ;
			int mergeDistance = argMap.containsKey("mergeDistance") ? argMap.getInteger("mergeDistance") : 1000;
			int window = argMap.getInteger("window");
			String col = argMap.getMandatory("colName");

			List<Double> dataDistribution = new ArrayList<Double>();
			List<Double> minDistribution  = new ArrayList<Double>();
			List<Double> maxDistribution  = new ArrayList<Double>();
			BufferedReader br = argMap.getInputReader();
			try {
				LocationAwareMatrix matrix = new LocationAwareMatrix(br, chr);
				ContinuousData data = new ContinuousData(matrix, col);
				WindowPermuteSegmenter wps = new WindowPermuteSegmenter(data, mergeDistance);
				for(int i = 0; i < randomizations; i++) {
					List<Double> permutationData = wps.computeRandomScores(chr, window);
					dataDistribution.addAll(permutationData);
					minDistribution.add(Statistics.min(permutationData));
					maxDistribution.add(Statistics.max(permutationData));
				}
			} finally {
				if(br != null) {
					br.close();
				}
			}
			
			writeList(outdir + "/"+col+".permutations.txt", dataDistribution);
			writeList(outdir + "/"+col+".permutations.max.txt", maxDistribution);
			writeList(outdir + "/"+col+".permutations.min.txt", minDistribution);
		} else if ("8".equals(argMap.getTask())) {
			String regions = argMap.getMandatory("regions");
			String filterRegionFile = argMap.getMandatory("filterRegions");
			double maxPctInRegions = argMap.containsKey("maxPctInFilterRegions") ? argMap.getDouble("maxPctInFilterRegions") : 0.2;
			System.err.println("Filtering out regions with at most " + maxPctInRegions + " in filter regions ");
			BEDReader regionReader  = new BEDReader(regions);
			Set<LightweightGenomicAnnotation> passingRegions = new TreeSet<LightweightGenomicAnnotation>();
			Iterator<String> chrIt = regionReader.getChromosomeIterator();
			while(chrIt.hasNext()) {
				String chr = chrIt.next();
				BEDReader filterReader = loadChromosomeRegions(filterRegionFile, chr);
				System.err.println("Loaded filter regions for " + chr);
				IntervalTree<BED> chrAnnotationTree = regionReader.getChromosomeTree(chr);
				IntervalTree<BED> filterAnnotationTree = filterReader.getChromosomeTree(chr);
				if(filterAnnotationTree == null) {
					continue;
				}
				Iterator<BED> chrAnnotationIt = chrAnnotationTree.valueIterator();
				while(chrAnnotationIt.hasNext()) {
					BED annotation = chrAnnotationIt.next();
					Iterator<Node<BED>> overlapIt = filterAnnotationTree.overlappers(annotation.getStart(), annotation.getEnd());
					List<LightweightGenomicAnnotation> overlapList = new ArrayList<LightweightGenomicAnnotation> ();
					while(overlapIt.hasNext()) {
						overlapList.add(new BasicLightweightAnnotation(overlapIt.next().getValue()));
					}

					List<LightweightGenomicAnnotation> intersect = annotation.intersect(overlapList);
					int total = 0;
					for(LightweightGenomicAnnotation o : intersect) {
						total += o.length();
					}
					if(total/(double)annotation.length() < maxPctInRegions) {
						passingRegions.add(annotation);
					}
					
				}
			}

			
			BufferedWriter bw = argMap.getOutputWriter();
			write(passingRegions, bw);
			bw.close();
		} 
		else if("9".equals(argMap.getTask())) { //ADDED BY MORAN - make matrix of raw data (no floor or normalization)
			List<File> dataFiles = tau.getWorkFiles(argMap);
				
			if(dataFiles.isEmpty() ) {
				return;
			}
			String format = argMap.containsKey("format") ? argMap.get("format") : "GFF";
			File first = dataFiles.get(0);
			AnnotationReader<? extends GenomicAnnotation> reader = AnnotationReaderFactory.create(format); 
			LocationAwareMatrix matrix = new LocationAwareMatrix( );
			int rowsRead = reader.parse(first.getAbsolutePath(),  null, matrix);
		
			reader = null; //Help the GC by explicitly telling it now that the reader will no longer be referenced.
			List<String> colNames = new ArrayList<String>(dataFiles.size());
			for(File file : dataFiles) {
				colNames.add(CLUtil.replaceFileExtension(file.getName(),""));
			}
			matrix.setColumns(colNames);
			matrix.setDataDimensions(rowsRead, colNames.size());

			for(int i = 0; i < dataFiles.size(); i++) {
				tau.addRawData(dataFiles.get(i), matrix, format, colNames.get(i));
			}

			BufferedWriter bw = argMap.getOutputWriter();
			if(argMap.containsKey("gp")) {
				matrix.writeGCT(bw);
			} else {
				matrix.write(bw);
			}
			bw.close();
			
		}
		
		else if("10".equals(argMap.getTask())) { //ADDED BY MORAN ; extracts segmented regions (from task 2)
			//that overlap a list of other regions (such as lincRNAs). Outputs a list of segments and a mapping between the input regions and the 
			//extracted segments names.
			String set1In = argMap.getMandatory("genes");
			String set2In = argMap.getMandatory("regions");
			String outName = argMap.getMandatory("outPrefix");
			mapToAllOverlapingRegions(set1In,set2In,outName);
			
			
		}
		
		
		else {
			System.err.println("Unimplemented task " + argMap.getTask() + ", " + USAGE);
		}
	}


	private static void write(Set<LightweightGenomicAnnotation> passingRegions,
			BufferedWriter bw) throws IOException {
		for(LightweightGenomicAnnotation g : passingRegions) {
			bw.write((new BED(g)).toString());
			bw.newLine();
		}
		
	}


	private static BEDReader loadChromosomeRegions(String regions,
			final String chromosome) throws FileNotFoundException, IOException,
			ParseException {
		BEDReader reader = new BEDReader();
		BufferedReader fileReader = new BufferedReader(new FileReader(regions));
		reader.load(fileReader, AnnotationFactoryFactory.bedFactory, new GenomicAnnotationFilter<BED>() {
			
			public boolean accept(BED annotation){
				return chromosome.equals(annotation.getChromosome());
			}

			public boolean isEnough(BED annotation){
				return false;
			}
			
		});
		fileReader.close();
		return reader;
	}


	private static int[] getWindowSizes(ArgumentMap argMap) {
		int [] windowSizes = WindowPermuteSegmenter.DEFAULT_WINDOWS;
		if(argMap.containsKey("windows")){
			String [] windowSizesStr = argMap.getMandatory("windows").split(",");
			windowSizes = new int[windowSizesStr.length];
			for(int i = 0; i < windowSizesStr.length; i++) {
				windowSizes[i] = Integer.parseInt(windowSizesStr[i]);
			}
			
		}
		return windowSizes;
	}


	private static void writeList(String fileName, List<?> dataDistribution) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
		for(Object o : dataDistribution) {
			bw.write(o.toString());
			bw.newLine();
		}
		bw.close();
		
	}


	private  void addData(final File file, final MatrixWithHeaders matrix, String format,final String column, final double floor) throws ParseException, IOException {
		AnnotationReader<? extends GenomicAnnotation>  reader = AnnotationReaderFactory.create( format);
		reader.parse(file.getAbsolutePath(), null, new AnnotationHandler() {

			public void annotation(GenomicAnnotation a) {
				matrix.set(a.getLocationString(), column, a.getScore() < floor ? floor : a.getScore());
			}

			public void begin() {
				// Nothing to do
			}

			public void browserLine(String line) {
				// Nothing to do
			}

			public void eof() {
				System.err.println("Finished loading " + file);
			}

			public void track(String line) {
				// Not sure what to do if several track are provided. For now load them all.
			}
			
		});
		
	}
	
	private  void addRawData(final File file, final MatrixWithHeaders matrix, String format,final String column) throws ParseException, IOException {
		AnnotationReader<? extends GenomicAnnotation>  reader = AnnotationReaderFactory.create( format);
		reader.parse(file.getAbsolutePath(), null, new AnnotationHandler() {

			public void annotation(GenomicAnnotation a) {
				matrix.set(a.getLocationString(), column, a.getScore());
			}

			public void begin() {
				// Nothing to do
			}

			public void browserLine(String line) {
				// Nothing to do
			}

			public void eof() {
				System.err.println("Finished loading " + file);
			}

			public void track(String line) {
				// Not sure what to do if several track are provided. For now load them all.
			}
			
		});
		
	}

	private  List<File> getWorkFiles(ArgumentMap argMap) {
		List<File> dataFiles = new ArrayList<File>();
		if(argMap.hasInputFile()) {
			dataFiles.add(new File(argMap.getInput()));
		} else {
			String inputDir = argMap.getInputDir();
			List<String> extensions = CLUtil.listFromArray(argMap.getMandatory("ext").split(","));
			File inputDirFile = new File(inputDir);
			for(final String ext : extensions) {
				String [] files = inputDirFile.list(new FilenameFilter(){
					public boolean accept(File arg0, String arg1) {
						return arg1.endsWith(ext);
					}					
				});
				
				for(String file : files) {
					dataFiles.add(new File(inputDir + File.separator + file));
				}
			}
		}
		
		return dataFiles;
	}
	
	private static void mapToAllOverlapingRegions (String set1In, String set2In,String outname) throws IOException {
		BEDFileParser bed= new BEDFileParser();
		BEDFileParser set1= new BEDFileParser(set1In);
		BEDFileParser set2= new BEDFileParser(set2In);
		BufferedWriter bw = new BufferedWriter(new FileWriter(outname + "Map"));
		
		int multiNameCntr=0;
		List<RefSeqGene> set1genes=set1.GetGenes();
		String newName;
		
		for(RefSeqGene gene : set1genes)
		{
			IntervalTree<RefSeqGeneWithIsoforms> tree=set2.getOverlappers(gene);
			Iterator <RefSeqGeneWithIsoforms> it= tree.valueIterator();
			while (it.hasNext()){
				RefSeqGene overlapper =it.next();
				bed.addRefSeq(overlapper);
				bw.write(gene.getName()+"\t"+overlapper.getName()+"\n");
			}
		}
		
	  bw.close();
	  bed.writeFullBed(outname + "Segments");
		
	}
	
}
