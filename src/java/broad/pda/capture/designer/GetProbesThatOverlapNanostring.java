package broad.pda.capture.designer;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixGenerator;

public class GetProbesThatOverlapNanostring {

	public GetProbesThatOverlapNanostring(String probeReport, String nanostringProbes, String save) throws IOException{
		Collection<String> probes=BEDFileParser.loadList(probeReport);
		Map<String, String> nanostring=parseNano(BEDFileParser.loadList(nanostringProbes));
		
		Collection<String> overlappers=new ArrayList<String>();
		FileWriter writer=new FileWriter(save);
		
		int i=0;
		for(String probe: probes){
			String[] tokens=probe.split("\t");
			String name=tokens[0];
			System.err.println(name+" "+i);
			String seq=tokens[8];
			boolean aligned=align(seq, nanostring, name);
			if(aligned){overlappers.add(probe); writer.write(probe+"\n"); System.err.println(name+" found in nanostring");}
			i++;
		}
		
		writer.close();
		
	}
	
	private boolean align(String seq, Map<String, String> nanostring, String name) {
		String nanoSeq=nanostring.get(name);
		if(nanoSeq==null){
			for(String n: nanostring.keySet()){
				nanoSeq=nanostring.get(n);
				if(align(seq, nanoSeq)){
					System.err.println(name+" Found in nanostring "+nanoSeq);
					return true;
				}
			}
		}
		else{
			return align(seq, nanoSeq);
		}
		
		return false;
	}

	private boolean align(String seq, String nanoSeq) {
		jaligner.Sequence s1=new jaligner.Sequence(seq);
		jaligner.Sequence s2=new jaligner.Sequence(nanoSeq);
		jaligner.Sequence s3=new jaligner.Sequence(Sequence.reverseSequence(nanoSeq));
		Matrix matrix=MatrixGenerator.generate(1.0f, -1.0f);
		jaligner.Alignment alignment=SmithWatermanGotoh.align(s1, s2, matrix, 1000000, 1000000);
		
		if(alignment.getNumberOfMatches()>60){return true;}
		alignment=SmithWatermanGotoh.align(s1, s3, matrix, 1000000, 1000000);
		if(alignment.getNumberOfMatches()>60){return true;}
		
		return false;
	}

	private Map<String, String> parseNano(Collection<String> list) {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		for(String line: list){
			String[] tokens=line.split("\t");
			rtrn.put(tokens[0], tokens[1]);
		}
		
		
		return rtrn;
	}

	public static void main(String[] args)throws IOException{
		if(args.length>2){
			String probes=args[0];
			String nano=args[1];
			String save=args[2];
			new GetProbesThatOverlapNanostring(probes, nano, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=probeReport \n args[1]=nanostringProbes \n args[2]=save";
	
}
