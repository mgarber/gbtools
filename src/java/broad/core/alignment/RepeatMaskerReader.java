package broad.core.alignment;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import broad.core.annotation.AnnotationReader;
import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.GenomicAnnotationFilter;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;

public class RepeatMaskerReader  {
	private File source;
	private ArrayList<Repeat> repeatList;
	private HashMap<String, IntervalTree<Repeat>> chrRepeatTreeMap;
	private HashMap<String, List<Repeat>> repeatFamilyMap ;
	private TreeMap<String, RepeatStatistic> repeatStatistics = new TreeMap<String, RepeatStatistic>(); 
	private TreeMap<String, RepeatStatistic> repeatFamStatistics = new TreeMap<String, RepeatStatistic>(); 
	private HashMap<String, List<String>> familyMemberMap = new HashMap<String, List<String>>();
	
	private boolean useUCSCFormat;


	public RepeatMaskerReader() { 
		super();
		repeatList = new ArrayList<Repeat>();
		repeatFamilyMap  = new HashMap<String, List<Repeat>>();
		chrRepeatTreeMap     = new HashMap<String, IntervalTree<Repeat>>();
		this.useUCSCFormat = false;
	}
	
	public RepeatMaskerReader(boolean useUCSCFormat) {
		this();
		this.useUCSCFormat = useUCSCFormat;
	}
	
	public RepeatMaskerReader(String fileName) throws FileNotFoundException {
		this(false);
		init(new File(fileName), false);
	}
	
	public RepeatMaskerReader(File file, boolean cleanRepeats) throws FileNotFoundException {
		this();
		init(file, cleanRepeats);
	}
	
	public RepeatMaskerReader(String fileName, boolean cleanRepeats) throws FileNotFoundException {
		this();
		init(new File(fileName), cleanRepeats);
	}

	private void init(File source, boolean cleanRepeats) throws FileNotFoundException {
		this.source = source;
		repeatList = loadRepeats(source, cleanRepeats);
	}
	
	public ArrayList<Repeat> loadRepeats(File source,  GenomicAnnotationFilter<Repeat> filter, int shift) throws IOException {
		InputStream is = new FileInputStream(source);
		loadRepeats(is,filter, shift);
		is.close();
		return repeatList;
	}
	
	public ArrayList<Repeat> loadRepeats(InputStream is,  GenomicAnnotationFilter<Repeat> filter, int shift) throws FileNotFoundException {
		BufferedReader br = new BufferedReader(new InputStreamReader(is));
		String line;
		try {
			int lineNum = 0;
			while((line = br.readLine()) != null) {
				line = line.trim();
				if(line.startsWith("#") || lineNum++ < 3){
					continue;
				}	
				//System.out.println(line);
				Repeat repeat = null;
				try {
					if(!useUCSCFormat) {
						repeat = new Repeat(line.split("\\s+"));
					} else {
						repeat = new Repeat();
						repeat.initFromUCSCTableData(line.split("\\s+"));
					}
					repeat.setStart(repeat.getStart() + shift); 
					repeat.setEnd(repeat.getEnd() + shift);    
					
					if(filter.accept(repeat)) {
						add(repeat);						
					}
					
				} catch (NumberFormatException nfe) {
					System.err.println("Number format exception caught when parsing <"+line+"> "+nfe);
					throw new RuntimeException("Number Format Exception",nfe);
				}
				
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				br.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		return repeatList;		
	}
	
	public ArrayList<Repeat> loadRepeats(File source, boolean cleanUp, LightweightGenomicAnnotation inRegion, int shift, int totalQuerySize, String queryName) throws FileNotFoundException {
		System.out.println("Loading repeats from " + source);
		BufferedReader br = new BufferedReader(new FileReader(source));
		String line;
		try {
			int lineNum = 0;
			Repeat priorRepeat = null;
			while((line = br.readLine()) != null) {
				line = line.trim();
				if(line.startsWith("#") || lineNum++ < 3){
					continue;
				}	
				//System.out.println(line);
				Repeat repeat = null;
				try {
					if(!useUCSCFormat) {
						repeat = new Repeat(line.split("\\s+"));
					} else {
						repeat = new Repeat();
						repeat.initFromUCSCTableData(line.split("\\s+"));
					}
					repeat.setStart(repeat.getStart() + shift);
					repeat.setEnd(repeat.getEnd() + shift);
					
					if(queryName != null) {
						repeat.setQuery(queryName);
					}
					
					if(inRegion != null && !inRegion.overlaps(repeat)) {
						continue;
					}
					
					if(totalQuerySize > 0) {
						repeat.setQueryBasesLeft(totalQuerySize - repeat.getEnd());
					}
					if(cleanUp && priorRepeat != null && repeat.overlaps(priorRepeat)) {
						if(priorRepeat.contains(repeat)) {
							continue;
						} 
						mergeRepeats(priorRepeat, repeat);
					}
					priorRepeat = repeat;

					add(repeat);
				} catch (NumberFormatException nfe) {
					System.err.println("Number format exception caught when parsing <"+line+"> "+nfe);
					throw new RuntimeException("Number Format Exception",nfe);
				}
				
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				System.out.print("Closing "+source.getAbsolutePath());
				br.close();
				System.out.print(" ..... Closed\n");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		return repeatList;
	}
	public void loadOverlappingRepeats(File source, boolean cleanUp, AnnotationReader<? extends GenomicAnnotation> overlappingRegions, int shift) throws FileNotFoundException {
			System.err.println("Loading repeats from " + source);
			BufferedReader br = new BufferedReader(new FileReader(source));
			String line;
			try {
				int lineNum = 0;
				Repeat priorRepeat = null;
				while((line = br.readLine()) != null) {
					line = line.trim();
					if(line.startsWith("#") || lineNum++ < 3){
						continue;
					}	
					//System.out.println(line);
					Repeat repeat = null;
					try {
						if(!useUCSCFormat) {
							repeat = new Repeat(line.split("\\s+"));
						} else {
							repeat = new Repeat();
							repeat.initFromUCSCTableData(line.split("\\s+"));
						}
						repeat.setStart(repeat.getStart() + shift);
						repeat.setEnd(repeat.getEnd() + shift);
						
						List<? extends GenomicAnnotation> repeatOverlappers = overlappingRegions.getOverlappers(repeat);
						if(repeatOverlappers == null || repeatOverlappers.size() == 0) {
							continue;
						}
						
						if(cleanUp && priorRepeat != null && repeat.overlaps(priorRepeat)) {
							if(priorRepeat.contains(repeat)) {
								continue;
							} 
							mergeRepeats(priorRepeat, repeat);
						}
						priorRepeat = repeat;

						add(repeat);
					} catch (NumberFormatException nfe) {
						System.err.println("Number format exception caught when parsing <"+line+"> "+nfe);
						throw new RuntimeException("Number Format Exception",nfe);
					}
					
				}
			} catch (IOException e) {
				e.printStackTrace();
			} finally {
				try {
					System.out.print("Closing "+source.getAbsolutePath());
					br.close();
					System.out.print(" ..... Closed\n");
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		
		Iterator<String> familyIt = repeatFamilyMap.keySet().iterator();
		while(familyIt.hasNext()) {
			List<Repeat> repeats = repeatFamilyMap.get(familyIt.next());
			Collections.sort(repeats, new Comparator<Repeat>() {

				public int compare(Repeat arg0, Repeat arg1) {
					return arg0.getStart() == arg1.getStart() 
						? (int)(arg0.getEnd() - arg1.getEnd()) 
						: (int) (arg0.getStart() - arg1.getStart());
				}
				
			});
		}
		
		System.out.println("Repeats loaded so far " + repeatList.size());
	}
	/**
	 * @param fileName
	 * @throws FileNotFoundException 
	 */
	public ArrayList<Repeat> loadRepeats(File source, boolean cleanUp) throws FileNotFoundException {

		return loadRepeats(source, cleanUp, null);
	}
	
	public ArrayList<Repeat> loadRepeats(File source, boolean cleanUp, LightweightGenomicAnnotation inRegion) throws FileNotFoundException {
		return loadRepeats(source, cleanUp, inRegion, 0, 0, null);
	}
	
	public ArrayList<Repeat> loadRepeats(File source, boolean cleanUp, int shift, int newQuerySize, String queryName) throws FileNotFoundException {
		return loadRepeats(source, cleanUp, null, shift, newQuerySize, queryName);
	}
	
	public List<FlankingRepeatFamily> getFlankingRepeats(LightweightGenomicAnnotation region, int buffer, float maxDivergence) {
		List<FlankingRepeatFamily> flankingRepeats = new ArrayList<FlankingRepeatFamily>();
		Iterator<String> familyIt = repeatFamilyMap.keySet().iterator();
		String family = null;
		BasicGenomicAnnotation leftRegion = new BasicGenomicAnnotation("leftRegion");
		leftRegion.setStart(region.getStart() - buffer);
		leftRegion.setEnd(region.getStart() + buffer);
		
		BasicGenomicAnnotation rightRegion = new BasicGenomicAnnotation("rightRegion");
		rightRegion.setStart(region.getEnd() - buffer);
		rightRegion.setEnd(region.getEnd() + buffer);
		
		while(familyIt.hasNext()) {
			family = familyIt.next();
			List<Repeat> familyRepeats = repeatFamilyMap.get(family);
			Repeat r = null;
			FlankingRepeatFamily flankRepeatFam = new FlankingRepeatFamily(family);
			long atPosition = 0;
			Iterator<Repeat> repeatIt = familyRepeats.iterator();
			while(repeatIt.hasNext() && atPosition < region.getEnd() + buffer) {
				r = repeatIt.next();
				if(r.getPercentDivergence() <= maxDivergence && leftRegion.contains(r)) {
					flankRepeatFam.addLeftRepeat(r);
				}
				if(r.getPercentDivergence() <= maxDivergence && rightRegion.contains(r)) {
					flankRepeatFam.addRightRepeat(r);
				}
				atPosition = r.getStart();
			}
			if(flankRepeatFam.isFlanking()) {
				flankingRepeats.add(flankRepeatFam);
			}
		}
		
		return flankingRepeats;
	}
	
	private void mergeRepeats(Repeat priorRepeat, Repeat repeat) {
		if(priorRepeat.percentDivergence < repeat.percentDivergence) {
			repeat.reduceToDifference(priorRepeat);
		} else {
			priorRepeat.reduceToDifference(repeat);
		}
		
	}
	
	public void add(Repeat repeat) {
		List<Repeat> repeatFamily = repeatFamilyMap.get(repeat.getRepeatFamily());
		if(repeatFamily == null) {
			repeatFamily = new ArrayList<Repeat>();
			repeatFamilyMap.put(repeat.getRepeatFamily(), repeatFamily);
		}
		repeatFamily.add(repeat);
		repeatList.add(repeat);	
		
		IntervalTree<Repeat> chrRepeatTree = chrRepeatTreeMap.get(repeat.getChromosome());
		if(chrRepeatTree == null) {
			chrRepeatTree = new IntervalTree<Repeat>();
			chrRepeatTreeMap.put(repeat.getChromosome(), chrRepeatTree);
		}
		//System.out.println("Repeat " + repeat.getRepeatName()  + " chr" + repeat.getChromosome() + ":" + ":" + repeat.getStart() + "-" + repeat.getEnd());
		chrRepeatTree.put(repeat.getStart(), repeat.getEnd(), repeat);
	}

	public void writeHeader(BufferedWriter bw) throws IOException {
		bw.write(Repeat.rightJustify("SW", Repeat.MAX_SCORE_SIZE));
		bw.write(" perc ");
		bw.write("perc ");
		bw.write("perc ");
		bw.write(Repeat.leftJustify("Query", Repeat.MAX_SEQ_NAME_SIZE));
		bw.write(" ");
		bw.write(Repeat.leftJustify("position in query", 3*Repeat.MAX_GENOMIC_POS_SIZE + 2));
		bw.write("   ");
		bw.write(Repeat.leftJustify("matching", Repeat.MAX_SEQ_NAME_SIZE));
		bw.write(" ");
		bw.write(Repeat.leftJustify("repeat", Repeat.MAX_SEQ_NAME_SIZE));
		bw.write(" ");
		bw.write(Repeat.leftJustify("position in repeat", 3*Repeat.MAX_REPEAT_POS_SIZE + 2));
		bw.newLine();
		
		bw.write(Repeat.rightJustify("score", Repeat.MAX_SCORE_SIZE));
		bw.write(" div. ");
		bw.write("del. ");
		bw.write("ins. ");
		bw.write(Repeat.leftJustify("squence", Repeat.MAX_SEQ_NAME_SIZE));
		bw.write(" ");
		bw.write(Repeat.leftJustify("begin", Repeat.MAX_GENOMIC_POS_SIZE));
		bw.write(" ");
		bw.write(Repeat.leftJustify("end", Repeat.MAX_GENOMIC_POS_SIZE));
		bw.write(" ");
		bw.write(Repeat.leftJustify("(left)", Repeat.MAX_GENOMIC_POS_SIZE));
		bw.write("   ");
		bw.write(Repeat.leftJustify("repeat", Repeat.MAX_SEQ_NAME_SIZE));
		bw.write(" ");
		bw.write(Repeat.leftJustify("class/family", Repeat.MAX_SEQ_NAME_SIZE));
		bw.write(" ");
		bw.write(Repeat.leftJustify("begin", Repeat.MAX_REPEAT_POS_SIZE));
		bw.write(" ");
		bw.write(Repeat.leftJustify("end", Repeat.MAX_REPEAT_POS_SIZE));
		bw.write(" ");
		bw.write(Repeat.leftJustify("(left)", Repeat.MAX_REPEAT_POS_SIZE));
		bw.write(" ID  ");
		bw.newLine();
		bw.newLine();
	}
	
	public void writeRepeatList(List<Repeat> repeats, BufferedWriter bw) throws IOException {
		Iterator<Repeat> it = repeats.iterator();
		while(it.hasNext()) {
			bw.write(it.next().format());
			bw.newLine();
		}
		
	}
	
	public void writeRepeatList(BufferedWriter bw) throws IOException {
		writeRepeatList(repeatList, bw);
		bw.flush();
	}
	
	public void resetQuery(String newQueryName) {
		Iterator<Repeat> it = repeatList.iterator();
		while(it.hasNext()) {
			it.next().setQuery(newQueryName);
		}
	}
	
	public void computeDirectoryStatistics(String directory) throws IOException {
		File dir = new File(directory);
		
		if (!dir.isDirectory()) {
			throw new IOException(directory + " is not a directory");
		}
		
		String [] outFiles = dir.list(new FilenameFilter() {

			public boolean accept(File arg0, String arg1) {
				return arg1.endsWith(".out");
			}
			
		});
		
		for(int i = 0; i < outFiles.length; i++) {
			String repeatFile = dir.getAbsolutePath() + "/" + outFiles[i];
			System.out.println("Cumputing statistics for " + repeatFile);
			computeStatistics(repeatFile);
		}
	}
	
	public void computeStatistics(String fileName)  throws IOException {
		source = new File(fileName);
		
		BufferedReader br = new BufferedReader(new FileReader(source));
		String line = null;
		try {
			Repeat repeat = null;
			while((line = br.readLine()) != null) {
				line = line.trim();
				if(line.startsWith("#") || line.startsWith("SW") || line.startsWith("score") || line.trim().length() == 0){
					continue;
				}	
				
				//System.out.println(line);
				try {
					if(!useUCSCFormat) {
						repeat = new Repeat(line.split("\\s+"));
					} else {
						repeat = new Repeat();
						repeat.initFromUCSCTableData(line.split("\\s+"));
						//System.out.println("r: " + repeat.toString());
					}
				} catch (NumberFormatException nfe) {
					System.err.println("Number format exception caught when parsing <"+line+"> "+nfe);
					throw new RuntimeException("Number Format Exception",nfe);
				}
				RepeatStatistic repeatStat = repeatStatistics.get(repeat.getName());
				if(repeatStat == null) {
					repeatStat = new RepeatStatistic(repeat.getName());
					repeatStatistics.put(repeat.getName(), repeatStat);
				}
				repeatStat.add(repeat);
				//System.out.println("repeatStat " + repeatStat.toString());
				
				RepeatStatistic famRepeatStat = repeatFamStatistics.get(repeat.getRepeatFamily());
				if(famRepeatStat == null) {
					famRepeatStat = new RepeatStatistic(repeat.getRepeatFamily());
					repeatFamStatistics.put(repeat.getRepeatFamily(), famRepeatStat);
				}
				famRepeatStat.add(repeat);
				
				List<String> familyMemebers = familyMemberMap.get(repeat.getRepeatFamily());
				if(familyMemebers == null) {
					familyMemebers = new ArrayList<String>();
					familyMemberMap.put(repeat.getRepeatFamily(), familyMemebers);
				}
				
				if(!familyMemebers.contains(repeat.getRepeatName())) {
					familyMemebers.add(repeat.getRepeatName());
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				System.out.print("Closing "+fileName);
				br.close();
				System.out.print(" ..... Closed\n");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}		
		
	}
	
	public void computeStatistics() {
		Iterator<Repeat> it = repeatList.iterator();
		Repeat r = null;
		
		int lineNum = 0;
		while(it.hasNext()) {
			r = it.next();
			RepeatStatistic repeatStat = repeatStatistics.get(r.getName());
			if(repeatStat == null) {
				repeatStat = new RepeatStatistic(r.getName());
				repeatStatistics.put(r.getName(), repeatStat);
			}
			repeatStat.add(r);
			
			
			
			RepeatStatistic famRepeatStat = repeatFamStatistics.get(r.getRepeatFamily());
			if(famRepeatStat == null) {
				famRepeatStat = new RepeatStatistic(r.getRepeatFamily());
				repeatFamStatistics.put(r.getRepeatFamily(), famRepeatStat);
			}
			famRepeatStat.add(r);
			
			List<String> familyMemebers = familyMemberMap.get(r.getRepeatFamily());
			if(familyMemebers == null) {
				familyMemebers = new ArrayList<String>();
				familyMemberMap.put(r.getRepeatFamily(), familyMemebers);
			}
			
			if(!familyMemebers.contains(r.getRepeatName())) {
				familyMemebers.add(r.getRepeatName());
			}
			
			if(lineNum % 100000 == 0)
				System.out.println("Total Mememory after adding item: "+Runtime.getRuntime().totalMemory()+" ... ");	
			lineNum++;
		}
		
	}
	
	public Map<String,List<String>> getFamilyMemeberMap() { return familyMemberMap;}
	
	public RepeatStatistic getFamilyStatistic(String familyName) {return repeatFamStatistics.get(familyName);}
	public RepeatStatistic getRepeatStatistic(String repeatName) {return repeatStatistics.get(repeatName);}
	
	public List<Repeat> getRepeatList() { return repeatList; }
	
	public List<Repeat> getRepeatList(GenomicAnnotation target) {
		ArrayList<Repeat>  targetReps = new ArrayList<Repeat>();
		Iterator<Repeat> it = repeatList.iterator();
		int currStart = 0;
		while(it.hasNext() && currStart < target.getEnd()) {
			Repeat r = it.next();
			if(r.overlaps(target)) {
				targetReps.add(r);
			}
			currStart = r.getStart();
		}
		return targetReps;
	}
	
	public Map<String, RepeatStatistic> computeStatisticOverlapping(GenomicAnnotation target) {
		Repeat r = null;
		List<Repeat> repeats = getRepeatList(target);
		Iterator<Repeat> it = repeats.iterator();
		Map<String, RepeatStatistic> familyStatistics = new HashMap<String, RepeatStatistic>();
		
		while(it.hasNext()) {
			r = it.next();
			String repeatFamily = r.getGeneralFamily();
			
			RepeatStatistic famRepeatStat = familyStatistics.get(repeatFamily);
			if(famRepeatStat == null) {
				famRepeatStat = new RepeatStatistic(repeatFamily);
				familyStatistics.put(repeatFamily, famRepeatStat);
			}
			famRepeatStat.add(r);
		}
		
		return familyStatistics;
	}
	
	public Iterator<Repeat> overlappers(LightweightGenomicAnnotation region) {
		IntervalTree<Repeat> chrRepeatTree = chrRepeatTreeMap.get(region.getChromosome());
		
		return new IntervalTree.ValuesIterator<Repeat>(chrRepeatTree.overlappers(region.getStart(), region.getEnd()));
	}
	
	public void unchunkGenomeRepeats(File sourceDir, File outDir) throws Exception{
		HashMap<String, LinkedList<Repeat>> chrRepeatMap = new HashMap<String, LinkedList<Repeat>>();
		
		if(!sourceDir.isDirectory()) {
			throw new Exception("sourceDir<"+sourceDir+">"+" was not a directory");
		}
		if(!outDir.isDirectory()) {
			throw new Exception("outDir<"+outDir+">"+" was not a directory");
		}

		String [] files = sourceDir.list();
		ArrayList<Repeat> repeats = null;
		String chr = null;
		for(int i = 0; i < files.length; i++) {
			repeats = loadRepeats(new File(sourceDir + "/" + files[i]), true);
			if(repeats.size() == 0) {
				System.out.println("WARNING: file "+sourceDir + "/" + files[i]+" contained no repeats");
				continue;
			}
			Repeat repeat = null;
			repeat = repeats.get(0);
			String query = repeat.getQuery();
			chr = query.substring(0, query.indexOf(".") );
			
			Iterator<Repeat> it = repeats.iterator();
			int line = 0;
			while(it.hasNext()) {
				line++;
				repeat = it.next();
				String rawQuery = repeat.getQuery();
				repeat.setQuery("chr" + chr);
				String [] nameParts = rawQuery.split("/");
				int shift = 0;
				String baseChunk = nameParts[0].substring(nameParts[0].indexOf(".") + 1);
				if(nameParts.length == 3) {
					String [] chunkBaseCoords = baseChunk.split("-");
					String [] chunkCoords     = nameParts[2].split("-");
					shift = Integer.parseInt(chunkBaseCoords[0]) + Integer.parseInt(chunkCoords[0]) - 1;
				}else if(nameParts.length == 1) {
					String [] chunkBaseCoords = baseChunk.split("-");
					shift = Integer.parseInt(chunkBaseCoords[0].substring(2));
				} else {
					throw new Exception("query for repeat "+repeat.format()+"is not in the expected format line " + line+ " file " + files[i]);
				}

				repeat.setStart(repeat.getStart() + shift);
				repeat.setEnd(repeat.getEnd() + shift);
				repeat.setQueryBasesLeft(-1);
				LinkedList<Repeat> chrRepeat = chrRepeatMap.get(chr);
				if(chrRepeat == null)  {
					//System.out.println("chrRepeat was null for "+chr+" creating new");
					chrRepeat = new LinkedList<Repeat>();
					chrRepeatMap.put(chr, chrRepeat);
				}
				chrRepeat.add(repeat);
			}
		}
		Iterator<String> chrIt = chrRepeatMap.keySet().iterator();
		while(chrIt.hasNext()) {				
			String curChr = chrIt.next();
			//System.out.println("current chr " + curChr);
			List<Repeat> repList = chrRepeatMap.get(curChr); 
			Collections.sort(repList, new Comparator<Repeat>() {
				public int compare(Repeat arg0, Repeat arg1) {
					return (int) (arg0.getStart() == arg1.getStart() 
								? arg0.getEnd() - arg1.getEnd()
								: arg0.getStart() - arg1.getStart());
				}
				
			});
			
			File chrFile = new File(outDir+"/chr"+curChr+".fa.out");
			boolean newFile = !chrFile.exists();
			BufferedWriter bw = new BufferedWriter(new FileWriter(chrFile, true));
			if(newFile) {
				writeHeader(bw);
			} 
			writeRepeatList(repList, bw);
			bw.close();
		}			

	}
	
	public List<RepeatStatistic> getRepeatStatistics() {
		List<RepeatStatistic> repeatStats = new ArrayList<RepeatStatistic>(repeatStatistics.values());
		Collections.sort(repeatStats);
		return repeatStats;
	}
	
	public List<RepeatStatistic> getRepeatFamilyStatistics() {
		List<RepeatStatistic> repeatStats = new ArrayList<RepeatStatistic>(repeatFamStatistics.values());
		Collections.sort(repeatStats);
		return repeatStats;
	}

	public void clearRepeats() {
		repeatList.clear();
		repeatFamStatistics.clear();
		repeatStatistics.clear();
	}
	
	
	public static class RepeatStatistic implements Comparable<RepeatStatistic>{
		private String name;
		private int memberNumber;
		private long totalBases;
		private float averageDivergence;
		private float lengthWieightedAverageDivergence;
		
		public RepeatStatistic(String name) {
			this.name = name;
		}
		
		public void add(Repeat r) {
			int length = r.getRepeatEnd() - r.getRepeatStart();
			averageDivergence = (memberNumber * averageDivergence + r.getPercentDivergence())/(memberNumber+1);
			lengthWieightedAverageDivergence = (totalBases * lengthWieightedAverageDivergence + length * r.getPercentDivergence())/(totalBases + length);
			++memberNumber;
			totalBases += length;
			
		}
		
		public void add(RepeatStatistic rs) {
			if(!rs.getName().equals(name)) {
				throw new IllegalArgumentException("Traying to add a statistic of repeat <"+rs.getName()+"> to <"+name+">");
			}
			averageDivergence = ((memberNumber * averageDivergence) + (rs.getAverageDivergence() * rs.getMemberNumber())) / (rs.getMemberNumber() + memberNumber);
			lengthWieightedAverageDivergence = (totalBases * lengthWieightedAverageDivergence + rs.getTotalBases() * rs.getLengthWieightedAverageDivergence())/(totalBases + rs.getTotalBases());
			memberNumber += rs.getMemberNumber();
			totalBases   += rs.getTotalBases();
		}

		public float getAverageDivergence() {
			return averageDivergence;
		}

		public float getLengthWieightedAverageDivergence() {
			return lengthWieightedAverageDivergence;
		}

		public int getMemberNumber() {
			return memberNumber;
		}
		

		public String getName() {
			return name;
		}

		public long getTotalBases() {
			return totalBases;
		}
		
		public float getAverageLength() {
			return ((float) totalBases)/((float) memberNumber);
		}
		
		public boolean equals(RepeatStatistic stat) {
			return (name.equals(stat.name) && (getMemberNumber() == stat.getMemberNumber()) && (getTotalBases() == stat.getTotalBases()));
		}
		
		public String toString() {
			StringBuilder sb = new StringBuilder();
			sb.append(getName()).append("\t")
				.append(getTotalBases()).append("\t")
				.append(getMemberNumber()).append("\t")
				.append(getAverageDivergence()).append("\t")
				.append(getAverageLength());
			
			return sb.toString();
		}
		
		public int hashCode() {
			int hashcode = name.hashCode();
			hashcode = hashcode << 10000;
			hashcode += getTotalBases();
			
			return hashcode;
		}

		public int compareTo(RepeatStatistic arg0) {
			return (int) (arg0.getTotalBases() -  getTotalBases());
		}
		
	}
	
	public static class FlankingRepeatFamily {
		private String familyName;
		private List<Repeat> leftRepeats = new ArrayList<Repeat>();
		private List<Repeat> rightRepeats = new ArrayList<Repeat>();

		public FlankingRepeatFamily(String familyName){
			this.familyName = familyName;
		}
		
		public void addLeftRepeat(Repeat r) {
			leftRepeats.add(r);
		}
		
		public void addRightRepeat(Repeat r) {
			rightRepeats.add(r);
		}
		
		public boolean isFlanking() {
			return (leftRepeats.size() > 0) && (rightRepeats.size() > 0);
		}

		public String getFamilyName() {
			return familyName;
		}

		public List<Repeat> getLeftRepeats() {
			return leftRepeats;
		}

		public List<Repeat> getRightRepeats() {
			return rightRepeats;
		}
		
		public void write(BufferedWriter bw) throws IOException {
			bw.write("\t" + familyName );
			bw.newLine();
			bw.write("\t\tLEFT REPEATS\n");
			Repeat r = null;
			Iterator<Repeat> it = leftRepeats.iterator();
			while(it.hasNext()) {
				r = it.next();
				bw.write("\t\t");
				bw.write(r.getRepeatName());
				bw.write("\t");
				bw.write(String.valueOf(r.getPercentDivergence()));
				bw.write("\t");
				bw.write(String.valueOf(r.getStart()));
				bw.write("\t");
				bw.write(String.valueOf(r.getEnd()));
				bw.write("\t");
				bw.write(r.inReversedOrientation() ? "-" : "+");
				bw.newLine();
			}
			bw.write("\t\tRIGHT REPEATS\n");
			
			it = rightRepeats.iterator();
			while(it.hasNext()) {
				r = it.next();
				bw.write("\t\t");
				bw.write(r.getRepeatName());
				bw.write("\t");
				bw.write(String.valueOf(r.getPercentDivergence()));
				bw.write("\t");
				bw.write(String.valueOf(r.getStart()));
				bw.write("\t");
				bw.write(String.valueOf(r.getEnd()));
				bw.write("\t");
				bw.write(r.inReversedOrientation() ? "-" : "+");
				bw.newLine();
			}
		}
	}
	
	public static void main(String [] args) throws Exception {
		if(args.length != 2) {
			System.err.println("USAGE: RepeatMaskerReader inDir outDir");
			return;
		}
		
		RepeatMaskerReader rmr = new RepeatMaskerReader();
		
		rmr.unchunkGenomeRepeats(new File(args[0]), new File(args[1]));
	}



}
