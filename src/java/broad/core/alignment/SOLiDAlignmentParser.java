package broad.core.alignment;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import broad.core.sequence.AbstractFastaParser;
import broad.core.sequence.FastaHandler;
import broad.core.sequence.FastaParser;

public class SOLiDAlignmentParser implements FastaHandler{
	public static final Pattern SOLID_PAIRED_READ_NAME_PATTERN = Pattern.compile("([0-9]+_[0-9]+_[0-9]+)(_)(R|F)([0-9]+)");
	public static final Pattern SOLID_ALN_PATTERN              = Pattern.compile("(-?)([0-9]+)(\\.)([0-9]+)"); 
	public static final int SOLID_READ_SIZE = 25; //This is a hack to avoid parsing the read colors to figure out its size...

	//private HashMap<String, PairedEndAlignment> pairedReadAlignments;
	private int expectedReadDist;
	private int readDistVariance;
	private String parsingQueryHits;
	private HashMap<String, PairedEndAlignment> concordantAlignments;
	private HashMap<String, PairedEndAlignment> discordantAlignments;
	private int querySize;
	
	public SOLiDAlignmentParser(int expectedReadDist, int readDistVariance) {
		super();
		this.expectedReadDist = expectedReadDist;
		this.readDistVariance = readDistVariance;
		//pairedReadAlignments = new HashMap<String, PairedEndAlignment>();
		concordantAlignments = new HashMap<String, PairedEndAlignment>();
		discordantAlignments = new HashMap<String, PairedEndAlignment>();
	}
	
	public String getParsingQueryHits() {
		return parsingQueryHits;
	}

	public void setParsingQueryHits(String parsingQueryHits) {
		this.parsingQueryHits = parsingQueryHits;
	}

	public int getReadDistVariance() {
		return readDistVariance;
	}

	public void setReadDistVariance(int readDistVariance) {
		this.readDistVariance = readDistVariance;
	}

	public Map<String, PairedEndAlignment> getConcordantAlignments() {
		return concordantAlignments;
	}

	public Map<String, PairedEndAlignment> getDiscordantAlignments() {
		return discordantAlignments;
	}

	/*
	public Map<String, PairedEndAlignment> getPairedReadAlignments() {
		return pairedReadAlignments;
	}
	*/
	public void setExpectedReadDist(int expectedReadDist) {
		this.expectedReadDist = expectedReadDist;
	}

	public void load(File source, String query, int querySize) throws FileNotFoundException, IOException {
		parsingQueryHits = query;
		this.querySize = querySize;
		FastaParser fp = new FastaParser();
		fp.parse(new FileInputStream(source), this);
	}
	
	public void eof(AbstractFastaParser parser) throws IOException {
		//We do not care about the last bit of the file.
	}

	public void newBase(AbstractFastaParser parser) throws IOException {
		// Do nothing we do not care about read bases
		
	}

	public void newSequence(AbstractFastaParser parser) throws IOException {
		String sequence = parser.getCurrentSequenceId();
		String [] seqInfo = sequence.split(",");
		if (seqInfo.length > 1) {
			String readName = seqInfo[0];
			Matcher m = SOLID_PAIRED_READ_NAME_PATTERN.matcher(readName);
			if(m.find()) {
				
				String fr = m.group(3);
				String read = m.group(1) + m.group(2) + m.group(4);
				
				//System.out.println("Read " + read);
				
				PairedEndAlignment pra = concordantAlignments.get(read);
				if(pra == null) { discordantAlignments.get(read);} 

				if(pra == null) {
					pra = new PairedEndAlignment(read);
					concordantAlignments.put(read, pra);
				}
				
				for(int i = 1; i < seqInfo.length; i++) {
					String alignment = seqInfo[i];
					m = SOLID_ALN_PATTERN.matcher(alignment);
					if(m.find()) {
						boolean isReveserMatch = "-".equals(m.group(1));
						int refPosition = Integer.parseInt(m.group(2));
						int missmatches = Integer.parseInt(m.group(4));
						pra.addAlignment(parsingQueryHits, refPosition, missmatches, isReveserMatch, "F".equals(fr), querySize);
					}else {
						throw new IOException("Alignment " + alignment + " did not match SOLiD alignment pattern " + SOLID_ALN_PATTERN.toString());
					}
				}
				
				if(pra.isDiscordant(expectedReadDist, readDistVariance)) {
					concordantAlignments.remove(read);
					discordantAlignments.put(read, pra);
				} else {
					discordantAlignments.remove(read);
					concordantAlignments.put(read, pra);
				}
			}else {
				throw new IOException(readName + " did not match SOLiD read name pattern " + SOLID_PAIRED_READ_NAME_PATTERN.toString());
			}
		}
	}
	
	public static class PairedEndAlignment {
	    static final int MAX_ALNS = 5;
	    static final int MAX_MISSMATCH = 1;
		String name;
		List<AlignmentSummary> forwardAlignments;
		List<AlignmentSummary> reverseAlignments;
		
		public PairedEndAlignment(String name) {
			this.name = name;
			forwardAlignments = new ArrayList<AlignmentSummary>();
			reverseAlignments = new ArrayList<AlignmentSummary>();
		}
		public boolean isDiscordant(int expectedPairDistance, int pairDistanceVariance) {
			boolean isDiscordant = false;
			
			int dist = getDistanceCloseTo(expectedPairDistance);
			
			if(dist > 0 && Math.abs(dist - expectedPairDistance) > pairDistanceVariance) {
				isDiscordant = true;
			}
			
			return isDiscordant;
		}
		
		public int getDistanceCloseTo(int distanceToCompare) {
			int dist = 0;
			
			if(forwardAlignments.size() > 0 && reverseAlignments.size() > 0) {
				Iterator<AlignmentSummary> forwardHitIt = forwardAlignments.iterator();
				while(forwardHitIt.hasNext()) {
					Iterator<AlignmentSummary> reverseHitIt = reverseAlignments.iterator();
					AlignmentSummary fHit = forwardHitIt.next();
					while(reverseHitIt.hasNext()) {
						AlignmentSummary rHit = reverseHitIt.next();
						
						int hitDist = fHit.getA().getDistanceTo(rHit.getA());
						dist = dist == 0 ? hitDist : Math.min(Math.abs(distanceToCompare - hitDist), dist);
					}
				}
			}
			//System.out.println("Distance " + dist);
			return dist;
		}
		public void addAlignment(String query, int refPosition, int missmatches, boolean reverseMatch, boolean isForward, int querySize) {
			if(missmatches > MAX_MISSMATCH) {
				return;
			}
			AlignmentSummary aln = new AlignmentSummary();
			//System.out.println("Adding alignment query " + query + " refPos " + refPosition + " misses " + missmatches + " ForwardMatch? " + !reverseMatch + " is forward read " + isForward);
			aln.setQuery(query);
			aln.getA().setChromosome(query);
			aln.setMissmatches(missmatches);
			aln.setReversedOrientation(reverseMatch);
			aln.setQueryStart(reverseMatch ? querySize - refPosition : refPosition);
			aln.setQueryEnd(reverseMatch ? querySize - refPosition - SOLID_READ_SIZE : refPosition + SOLID_READ_SIZE);
			
			if(isForward) {
			    addAlignment(forwardAlignments,aln);
			} else {
			    addAlignment(reverseAlignments,aln);
			}
		}
		public List<AlignmentSummary> getForwardAlignments() {
			return forwardAlignments;
		}

	    private void addAlignment(List<AlignmentSummary> alns, AlignmentSummary aln) {
	    	if(alns.size() > 0) {
	    		if(alns.get(0).getMissmatches() == aln.getMissmatches()) {
	    			alns.add(aln);
	    		} else if(alns.get(0).getMissmatches() < aln.getMissmatches()) {
	    			alns.clear();
	    			alns.add(aln);
	    		}
	    	} else {
	    		alns.add(aln);
	    	}
	    }

		public List<AlignmentSummary> getReverseAlignments() {
			return reverseAlignments;
		}
		
		
		public String getName() {
			return name;
		}
		public AlignmentSummary getPivot() {
			AlignmentSummary pivot = null;
			if(isUnambiguous()) {
				AlignmentSummary fAln = forwardAlignments.get(0);
				AlignmentSummary rAln = reverseAlignments.get(0);
				
				if(fAln.getQuery().equals(rAln.getQuery())) {
					pivot = fAln.getA().getStart() < rAln.getA().getStart() ? fAln : rAln;
				} else {
					pivot = fAln.getQuery().compareTo(rAln.getQuery()) >= 0 ? rAln : fAln;
				}
			}else if(forwardAlignments.size() == 1) {
				pivot = forwardAlignments.get(0);
			} else if (reverseAlignments.size() == 1) {
				pivot = reverseAlignments.get(0);
			}
			
			return pivot;
		}
		
		public boolean isUnambiguous() {
			return forwardAlignments.size() == 1 && reverseAlignments.size() == 1;
		}
		
		public void write (BufferedWriter bw ) throws IOException {
			bw.write(getName());
			bw.write("\t");
			if(! isUnambiguous()) {
				bw.write(getPivot().getA().toString());
				bw.write("\t");
			}
			Iterator<AlignmentSummary> alnIt = getForwardAlignments().iterator();
			while(alnIt.hasNext()) {
				AlignmentSummary aln = alnIt.next();
				bw.write(aln.getA().toString());
				bw.write(";");
			}
			bw.write("\t");
			alnIt = getReverseAlignments().iterator();
			while(alnIt.hasNext()) {
				AlignmentSummary aln = alnIt.next();
				bw.write(aln.getA().toString());
				bw.write(";");
			}
		}

	}
	


}
