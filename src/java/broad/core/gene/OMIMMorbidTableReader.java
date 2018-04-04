package broad.core.gene;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;

public class OMIMMorbidTableReader {
	private HashMap<Integer,MorbidEntry> morbidEntryList; 
	
	public OMIMMorbidTableReader(String fileName) throws FileNotFoundException {
		
		super();
		morbidEntryList = new HashMap<Integer,MorbidEntry>();
		File source = new File(fileName);
		BufferedReader br = new BufferedReader(new FileReader(source));
		String line;
		try {
			while((line = br.readLine()) != null) {
				if(line.startsWith("#")){
					continue;
				}
				MorbidEntry entry = new MorbidEntry(line.split("\\|"));
				morbidEntryList.put(entry.getOMIMId(),entry);
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
	
	public Iterator<MorbidEntry> getOMIMEntryIterator() {
		// TODO Auto-generated method stub
		return morbidEntryList.values().iterator();
	}
	
	public MorbidEntry getEntry(Integer id) {
		return morbidEntryList.get(id);
	}
	
	public static class MorbidEntry  {

		private Integer omimId;
		private String name;
		private String[] genes;
		private String chromosomeBand;
		

		public MorbidEntry(String[] strings) {
			name = strings[0];
			genes = strings[1].split(",");
			omimId = Integer.valueOf(strings[2]);
			chromosomeBand = strings[3];
		}

		public String getName() {
			return name;
		}

		public Integer getOMIMId() {
			return omimId;
		}

		public String getChromosomeBand() {
			return chromosomeBand;
		}

		public String[] getGenes() {
			return genes;
		}
		
		public String toString() {
			return "Morbid Entry <"+getName()+" id<"+getOMIMId();
		}
		
	}

}
