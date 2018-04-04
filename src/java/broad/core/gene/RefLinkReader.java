package broad.core.gene;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

public class RefLinkReader {

	private HashMap<Integer,String> omimToRefSeqTable;
	private HashMap<String,Integer> refSeqToOmimTable;
	private HashMap<String, String> refSeqNameTable;
	private HashMap<String,List<String>> symbolRefSeqTable;
	private HashMap<Integer, String> entrezToRefSeqTable;
	private HashMap<String, Integer> refSeqToEntrezTable;
	
	public RefLinkReader(String fileName) throws FileNotFoundException {
		super();
		omimToRefSeqTable = new HashMap<Integer,String>();
		refSeqToOmimTable = new HashMap<String,Integer>();
		refSeqNameTable   = new HashMap<String,String>();
		symbolRefSeqTable = new HashMap<String,List<String>>();
		entrezToRefSeqTable = new HashMap<Integer, String>();
		refSeqToEntrezTable = new HashMap<String, Integer>();
		
		
		File source = new File(fileName);
		BufferedReader br = new BufferedReader(new FileReader(source));
		String line;
		try {
			while((line = br.readLine()) != null) {
				if(line.startsWith("#")){
					continue;
				}
				String[] splitLine = line.split("\t");
				Integer omimId = Integer.parseInt(splitLine[7]);
				String  symbol = splitLine[0];
				Integer entrezId = Integer.parseInt(splitLine[6]);
				String mRNAId   =  splitLine[2];
				omimToRefSeqTable.put(omimId,mRNAId);
				refSeqToOmimTable.put(mRNAId,omimId);
				refSeqNameTable.put(mRNAId, symbol);
				List<String> refSeqListForSymbol = symbolRefSeqTable.get(splitLine[0]);
				if(refSeqListForSymbol == null) {
					refSeqListForSymbol = new ArrayList<String>();
					symbolRefSeqTable.put(symbol, refSeqListForSymbol);
				}
				refSeqListForSymbol.add(splitLine[2]);
				refSeqToEntrezTable.put(mRNAId, entrezId);
				entrezToRefSeqTable.put(entrezId, mRNAId);
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				System.out.print("Closing "+fileName);
				br.close();
				System.out.print(" ..... Closed\n");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}	
	}

	public Integer getOmimId(String refSeqId) {
		return refSeqToOmimTable.get(refSeqId);
	}
	
	public String getRefSeqId(Integer omimId) {
		return omimToRefSeqTable.get(omimId);
	}
	
	public String getRefSeqName(String refSeqId) {
		return refSeqNameTable.get(refSeqId);
	}
	
	public Iterator<Integer> getOMIMIdIterator() {
		return omimToRefSeqTable.keySet().iterator();
	}
	
	public Collection<String> getNMsForGeneSymbol(String symbol) {
		return symbolRefSeqTable.get(symbol);
	}

	public List<String> getRefSeqsForSymbol(String symbol) {
		List<String> refSeqs = symbolRefSeqTable.get(symbol);
		return refSeqs != null ? refSeqs : new ArrayList<String>();
	}
	
	public String getRefSeqForEntrezId(Integer entrezId) {return entrezToRefSeqTable.get(entrezId);}
	public Integer getEntrezIdForRefSeq(String refSeq) {return refSeqToEntrezTable.get(refSeq);}
}
