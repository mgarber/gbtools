package broad.pda.seq.segmentation;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.graph.ChromosomeWithBubbles2.MinReadFilter;

public class test {

	//use the rate of spliced reads globally to identify outliers on the negative tail of the distribution
	//lambda estimation should come from the number of spliced reads dividied by the number of spliced locations
	//also estimate locally, if passes both cutoffs then get rid of it
	
	//only filter if there are alternative isoforms in the graph
	//if only one path exists then keep it
	
	public test(AlignmentDataModel data, String save, String genomeDirectory, String chr)throws IOException{
		write(data, genomeDirectory, chr, save);
		
	}

	private void write(String save, Map<Alignments, double[]> map) throws IOException {
		FileWriter writer=new FileWriter(save);
				
		for(Alignments intron: map.keySet()){
			double[] vals=map.get(intron);
			writer.write(intron+"\t"+intron.getOrientation()+"\t"+vals[0]+"\t"+vals[1]+"\n");
		}
		
		writer.close();
	}
	
	private Sequence getChrSeq(String genomeDir, String chr) throws IOException{
		String seqFile=genomeDir+"/"+chr.replaceAll("chr", "")+"/"+chr+".fa";
		FastaSequenceIO fsio = new FastaSequenceIO(seqFile);
		List<Sequence> seqs = fsio.extractRecordsWithIDLike(chr, false);
		if(!seqs.isEmpty()) {
			return seqs.get(0);
		}
		else{return null;}
	}

	private void write(AlignmentDataModel data, String genomeDirectory, String chr, String save) throws IOException {
		//Map<Alignments, double[]> rtrn=new TreeMap<Alignments, double[]>();
		
		FileWriter writer=new FileWriter(save);
		//for(String chr: data.getChromosomeLengths().keySet()){
			System.err.println(chr);
			Sequence chrSeq=getChrSeq(genomeDirectory, chr);
			Alignments chrRegion=new Alignments(chr, 0, data.getChromosomeLengths().get(chr));
			List<ReadFilter> filters = new ArrayList<ReadFilter>();
			if(chrSeq != null && chrSeq.getLength() >0) {
				filters.add(new CannonicalSpliceFilter(chrSeq));
			}
			filters.add(new MinReadFilter(10));
			
			Map<Alignments, Integer> spliced=data.getSplicedReads(chrRegion, filters, 1, chrSeq);
			
			for(Alignments splice: spliced.keySet()){writer.write(splice+"\t"+splice.getOrientation()+"\n");}
		//}
		
		writer.close();
	}
	
	
	
	private double count(Map<Alignments, Integer> spliced) {
		double num=0;
		
		for(Alignments intron: spliced.keySet()){
			num+=spliced.get(intron);
		}
		
		return num;
	}

	
	
	private double poissonCDF(int k, double lambda){
		cern.jet.random.Poisson poiss=new cern.jet.random.Poisson(lambda, new cern.jet.random.engine.DRand());
		return poiss.cdf(k);
	}

	public static void main(String[] args)throws Exception{
		if(args.length>4){
			AlignmentDataModel data=new GenericAlignmentDataModel(args[0], args[1]);
			String save=args[2];
			String genomeDirectory=args[3];
			String chr=args[4];
			new test(data, save, genomeDirectory, chr);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=alignment file (SAM) \n args[1]=sizes \n args[2]=save \n args[3]=genomeDirectory \n args[4]=chr";
	
	
}
