package broad.core.alignment;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

public class ArachneQLTReader {
	private File source;
	private ArrayList<QLTAlignment> alignmentList;
	HashMap<String, Integer> chrIdMap = new HashMap<String, Integer>();
	public static final String HG17_ID_MAP="/wga/scr3/LOOKUP/build35/hg17.ids";
	public static final String HG16_ID_MAP="/wga/scr3/LOOKUP/human34/build34.ids";
	public static final String HG18_ID_MAP=" /wga/scr3/LOOKUP/human36/hg18.ids";
	
	public ArachneQLTReader(String fileName) throws FileNotFoundException {
		super();
		alignmentList = new ArrayList<QLTAlignment>();
		source = new File(fileName);
		if(!source.exists()) {
			throw new FileNotFoundException("File "+ fileName + " does not exist");
		}		
	}
	
	public List<QLTAlignment> getAlignments() { return alignmentList; } 
	
	public void load(String chrIdMapFile, String chrSymbol, int start, int end) throws IOException {
		if(chrIdMap.isEmpty()) {
			loadChrIdMap(chrIdMapFile);
		}
		
		Integer chrId = chrIdMap.get("chr"+chrSymbol);
		if(chrId == null) {
			throw new RuntimeException("chromosome " + chrSymbol + " has no scaffold  id");
		}
		System.out.println("chr " + chrSymbol + " has id " + chrId);
		
		load(chrId, start, end);
	}
	
	public void load(int subjectId, int subjectStartPos, int subjectEndPos) 
	throws FileNotFoundException {
		BufferedReader br = new BufferedReader(new FileReader(source));
		String line;
		try {
			while((line = br.readLine()) != null) {
				if(!line.startsWith("QUERY")){
					continue;
				}
				QLTAlignment summary = new QLTAlignment(line.split("\t"));
				if(String.valueOf(subjectId).equals(summary.getSubject()) &&
					summary.getSubjectStart() <= subjectEndPos &&
					summary.getSubjectEnd() >= subjectStartPos ) {
					alignmentList.add(summary);
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
				e.printStackTrace();
			}
		}	
	}

	public static class QLTAlignment extends AlignmentSummary {
		int queryLength;
		int subjectLength;
		
		public int getQueryLength() {
			return queryLength;
		}

		public int getSubjectLength() {
			return subjectLength;
		}

		public QLTAlignment(String[] strings) {
			int i = 1;
			super.setQuery(strings[i++]);
			super.setQueryStart(Integer.parseInt(strings[i++]));
			super.setQueryEnd(Integer.parseInt(strings[i++]));
			queryLength = Integer.parseInt(strings[i++]);
			super.setQueryOrientation(strings[i++]);
			super.setSubject(strings[i++]);
			super.setSubjectStart(Integer.parseInt(strings[i++]));
			super.setSubjectEnd(Integer.parseInt(strings[i++]));
			subjectLength = Integer.parseInt(strings[i++]);
			
		}
	}
	
	public static void main(String [] args) throws IOException {
		if(args.length != 5) {
			System.err.println("Ussage ArachneQLTReader <input QLT file> <chromosome symbol> <subject start> <subject end> <output file>");
			return;
		}
		String chr = args[1];
		int start = Integer.parseInt(args[2]);
		int end   = Integer.parseInt(args[3]);
		String out = args[4];
		
		ArachneQLTReader qltR = new ArachneQLTReader(args[0]);
		qltR.load(HG16_ID_MAP, chr, start, end);
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		Iterator<QLTAlignment> it = qltR.alignmentList.iterator();
		while(it.hasNext()) {
			bw.write(it.next().toString());
			bw.newLine();
		}
		bw.close();

	}
	private void loadChrIdMap(String mapFile) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(mapFile));
		String line = null;
		int i = 0;
		while((line = br.readLine()) != null) {
			line.trim();
			String [] info = line.split("\\.");
			chrIdMap.put(info[0], i++);
		}
	}
}
