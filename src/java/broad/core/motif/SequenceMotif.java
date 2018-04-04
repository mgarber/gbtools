package broad.core.motif;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

public class SequenceMotif {
	static HashMap<Character, String> symbolToPatternMap;
	public static String USAGE = "Usage: SequenceMotif TASK=<task_num> <task_args>\n" +
	"\tTasks:\n" +
	"\t\t1. Get list of motif matches: IN=<fasta sequence> MOTIF=<Core Motif e.g. YRNRY> SUSCCESSIVE=<number of consecutive core motifs>\n";
	
	private Pattern motif;
	private Pattern [] coreMotifs;
	
	static {
		symbolToPatternMap = new HashMap<Character, String>();
		symbolToPatternMap.put('A',"A");
		symbolToPatternMap.put('C',"C");
		symbolToPatternMap.put('T',"T");
		symbolToPatternMap.put('G',"G");
		symbolToPatternMap.put('U',"U");
		
		symbolToPatternMap.put('R',"[AG]");
		symbolToPatternMap.put('Y',"[CTU]");
		
		symbolToPatternMap.put('M',"[AC]");
		symbolToPatternMap.put('K',"[GTU]");
		symbolToPatternMap.put('S',"[CG]");
		symbolToPatternMap.put('M',"[ATU]");
		
		symbolToPatternMap.put('B',"[^A]");
		symbolToPatternMap.put('D',"[^C]");
		symbolToPatternMap.put('H',"[^G]");
		symbolToPatternMap.put('V',"[^TU]");
		
		symbolToPatternMap.put('N',".");
	}


	public SequenceMotif(Pattern motif) {
		super();
		this.motif = motif;
	}
	 
	public SequenceMotif(String corePattern, int minsuccessive) throws SearchException {
		super();
		String [] patternList  = corePattern.split("-");
		coreMotifs = new Pattern[patternList.length];
		
		for(int j = 0; j < patternList.length; j++) {
			StringBuffer patternStr = new StringBuffer();
			char [] patternChars = patternList[j].toUpperCase().toCharArray();
			for(int i = 0; i < patternChars.length; i++) {
				String translatedPattern = symbolToPatternMap.get(patternChars[i]);
				if(translatedPattern == null) {
					throw new SearchException("The " +(i + 1) + "th character<" + patternChars[i] +"> of corePattern " + corePattern + " is not supported");
				}
				patternStr.append(translatedPattern);
			}
			coreMotifs[j] =  Pattern.compile(patternStr.toString());
		}
		
		if (coreMotifs.length == 1) {
			motif = Pattern.compile("(" + coreMotifs[0] + "){" + minsuccessive + ",}");
		} else {
			StringBuffer buf = new StringBuffer();
			for(int j = 0; j < coreMotifs.length; j++) {
				buf.append("(");
				buf.append(coreMotifs[j].pattern());
				buf.append("){" + minsuccessive + ",}");
				if(j < coreMotifs.length - 1) {
					buf.append("|");
				}
			}
			motif = Pattern.compile(buf.toString());
		}
	}
	
	protected void setMotif(Pattern motif) { this.motif = motif;}
	public Pattern getMotif() { return motif;}
	
	protected Pattern[] getCoreMotifs() { return coreMotifs;}
	
	public List<SequenceRegion> match(Sequence seq) {
		//System.out.println("Motif " + motif);
		ArrayList<SequenceRegion> matches = new ArrayList<SequenceRegion>();
		String seqBases = seq.getSequenceBases().toUpperCase();
		Matcher matcher = motif.matcher(seqBases);
		while(matcher.find()) {
			int start = matcher.start();
			int end   = matcher.end();
			matches.add(seq.getRegion(start, end));
		}
		
		return matches;
	}
	
	

	/**
	 * @param args
	 * @throws SearchException 
	 * @throws NumberFormatException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws NumberFormatException, SearchException, IOException {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE);
		if ("1".equals(argMap.getTask())) {
			SequenceMotif motif = new SequenceMotif(argMap.getMandatory("MOTIF"), argMap.getInteger("SUCCESSIVE"));
			FastaSequenceIO fsio = new FastaSequenceIO(argMap.getInput());
			Iterator<Sequence> seqIt = fsio.loadAll().iterator();
			Sequence seq = null;
			while(seqIt.hasNext()) {
				seq = seqIt.next();
				List<SequenceRegion> matches = motif.match(seq);
				//System.out.println("Sequence " + seq.getId());
				if(matches.size() == 0) {
					System.out.println("No matches");
				} else {
					Iterator<SequenceRegion> matchIt = matches.iterator();
					SequenceRegion seqReg = null;
					while(matchIt.hasNext()) {
						seqReg = matchIt.next();
						System.out.println((seqReg.getStart() + 1) + "\t" + (seqReg.getEnd())+ "\t" + seqReg.getSequenceBases());
					}
				}
				
			}
		} else {
			System.err.println(USAGE);
		}

	}

	public List<SequenceRegion> match(Sequence geneSeq, List<SequenceRegion> list, int length) {
		List<SequenceRegion> rtrn=new ArrayList<SequenceRegion>();
		
		for(SequenceRegion seq: list){
			Sequence subseq=geneSeq.getSubSequence(seq, length);
			rtrn.addAll(match(subseq));
		}
		
		return rtrn;
	}



}
