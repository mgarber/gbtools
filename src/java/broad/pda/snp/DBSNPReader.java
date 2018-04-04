package broad.pda.snp;

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
import java.util.TreeMap;

import broad.core.annotation.BasicGenomicAnnotation;

public class DBSNPReader  {
	public static File DB_SNP_DIR = new File("/seq/mgarber/dbSNP/");
	HashMap<String, DBSNP> snpMap;
	TreeMap<Integer, DBSNP> snpLocMap;
	
	public DBSNPReader(String chromosome) throws FileNotFoundException {
		snpMap = new HashMap<String, DBSNP>();
		snpLocMap = new TreeMap<Integer, DBSNP>();

		File chrDBSNPFile = new File(DB_SNP_DIR.getAbsolutePath() + "/chr" + chromosome+".snp");
		if(!chrDBSNPFile.exists()) {
			return;
		}
		BufferedReader br = new BufferedReader(new FileReader(chrDBSNPFile));
		

		String line;
		try {
			while((line = br.readLine()) != null) {
				if(line.startsWith("#")){
					continue;
				}
				//System.out.println(line.replace("\t", "-TAB-"));
				String [] lineSplit = line.split("\t");
				DBSNP snp = new DBSNP(lineSplit);
				snpMap.put(snp.getName(), snp);
				snpLocMap.put(snp.getStart(), snp);
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				System.out.print("Closing "+DB_SNP_DIR.getAbsolutePath() + "/chr" + chromosome+".snp");
				br.close();
				System.out.print(" ..... Closed\n");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}	
	}
	
	public Collection<DBSNP> getSNPList() { return snpLocMap.values();}
	
	public List<DBSNP> getSNPBetweenList(int start, int end) {
		Iterator<Integer> positionIt = snpLocMap.keySet().iterator();
		int position = 0;
		ArrayList<DBSNP> snps = new ArrayList<DBSNP>();
		
		while(position <= end && positionIt.hasNext()) {
			position = positionIt.next();
			if(start <= position && end >= position) {
				//System.out.println("SNP at pos " + position + ": " + snpLocMap.get(position));
				snps.add(snpLocMap.get(position));
			}
		}
		return snps;
	}
	
	public DBSNP getSNP(String dbSNPId) { return snpMap.get(dbSNPId);}
	
	public static class DBSNP extends BasicGenomicAnnotation {
		private String refNCBI;
		private String refUCSC;
		private String allele;
		private String moleculeType;
		private String sample;
		private String type;
		private float averageHeterozygocity;
		private float averageHeterozygocityStd;
		private String function;
		private String locationType;
		private String source;
		
		public DBSNP(String[] data) {
			int i = 1;
			setChromosome(data[i++].intern());
			setStart(data[i++]);
			setEnd(data[i++]);
			setName(data[i++]);
			setScore(Integer.parseInt(data[i++]));
			setReversedOrientation("-".equals(data[i++]));
			this.refNCBI = data[i++].intern();
			this.refUCSC = data[i++].intern();
			this.allele  = data[i++].intern();
			this.moleculeType = data[i++].intern();
			this.sample = data[i++].intern();
			this.type = data[i++].intern();
			this.averageHeterozygocity = Float.parseFloat(data[i++]);
			this.averageHeterozygocityStd = Float.parseFloat(data[i++]);
			this.function = data[i++].intern();
			this.locationType = data[i++].intern();
			this.source       = data[i++].intern();
		}
		
		public DBSNP(String id) {
			super(id);
		}

		public String getAlleles() {
			return allele;
		}

		public float getAverageHeterozygocity() {
			return averageHeterozygocity;
		}

		public float getAverageHeterozygocityStd() {
			return averageHeterozygocityStd;
		}

		public String getFunction() {
			return function;
		}

		public String getLocationType() {
			return locationType;
		}

		public String getMoleculeType() {
			return moleculeType;
		}

		public String getSource() {
			return source;
		}

		public String getType() {
			return type;
		}
		
		public String getSample() {
			return sample;
		}

		public String getRefNCBI() {
			return refNCBI;
		}

		public String getRefUCSC() {
			return refUCSC;
		}
		
	}



}
