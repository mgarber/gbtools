package broad.pda.tag;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.core.sequence.WindowSlider;
import broad.core.util.ParseGCTFile;

public class IdentifyHardToCloneSeqs {

	int windowSize=100;
	int overlap=99;
	double percentFilter=.9;
	
	
	
	public IdentifyHardToCloneSeqs(File proteinSequences, Map<String, String> refLink, Collection<String> genesOfInterest, File genes, String save) throws IOException {
		List<Sequence> proteinSequence = getSequences(proteinSequences);
		
		proteinSequence=filter(proteinSequence, genesOfInterest, refLink);
		
		//Find amino acid runs
		Map<String, int[]> runs=new TreeMap();
		for(Sequence seq: proteinSequence){
			int same6=0;
			int same7=0;
			WindowSlider slider=WindowSlider.getSlider(seq, 6, 5);
			while(slider.hasNext()){
				SequenceRegion region=slider.next();
				if(isSame(region)){same6++;}
			}
			slider=WindowSlider.getSlider(seq, 7, 6);
			while(slider.hasNext()){
				SequenceRegion region=slider.next();
				if(isSame(region)){same7++;}
			}
			//System.err.println(seq.getId()+"\t"+same6+"\t"+same7);
			int[] vals={same6, same7};
			runs.put(seq.getId(), vals);
		}
		
		Map<String, int[]> gc=new TreeMap();
		//Find high GC content
		List<Sequence> geneSequences=getSequences(genes);
		geneSequences=filter(geneSequences, genesOfInterest);
		
		for(Sequence seq: geneSequences){
			WindowSlider slider=WindowSlider.getSlider(seq, windowSize, overlap);
			int pct9=0;
			int pct8=0;
			int pct7=0;
			while(slider.hasNext()){
				SequenceRegion region=slider.next();
				double percentGC=region.percentGC();
				if(percentGC>.9){pct9++;}
				if(percentGC>.8){pct8++;}
				if(percentGC>.7){pct7++;}
			}
			int[] vals={pct7, pct8, pct9};
			gc.put(seq.getId().split("\\.")[0], vals);
		}
		
		write(save, runs, gc, refLink);
	}
	
	private void write(String save, Map<String, int[]> runs, Map<String, int[]> gc, Map<String, String> refLink) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("Gene\tRuns of 6\tRuns of 7\t70% GC\t80% GC\t90% GC\n");
		
		for(String prot: runs.keySet()){
			String gene=refLink.get(prot.split("\\.")[0]);
			int[] runVals=runs.get(prot);
			int[] gcVals=gc.get(gene);
			writer.write(gene+"\t"+runVals[0]+"\t"+runVals[1]+"\t"+gcVals[0]+"\t"+gcVals[1]+"\t"+gcVals[2]+"\n");
		}
		
		writer.close();
	}

	private boolean isSame(SequenceRegion region) {
		String seq=region.getSequenceBases();
		char[] seqs=seq.toCharArray();
		for(int i=1; i<seqs.length; i++){
			if(seqs[i]!=seqs[0]){return false;}
		}
		return true;
	}

	private List<Sequence> getSequences(File fa) throws IOException{
		FastaSequenceIO fsio = new FastaSequenceIO(fa);
		return fsio.loadAll();
	}
	
	/*public IdentifyHardToCloneSeqs(File geneSequence, Collection<String> genesOfInterest, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		FastaSequenceIO fsio = new FastaSequenceIO(geneSequence);
		List<Sequence> geneSequences = fsio.loadAll();
		
		geneSequences=filter(geneSequences, genesOfInterest);
		
		for(Sequence seq: geneSequences){
			WindowSlider slider=WindowSlider.getSlider(seq, windowSize, overlap);
			int pct9=0;
			int pct8=0;
			int pct7=0;
			while(slider.hasNext()){
				SequenceRegion region=slider.next();
				double percentGC=region.percentGC();
				if(percentGC>.9){pct9++;}
				if(percentGC>.8){pct8++;}
				if(percentGC>.7){pct7++;}
			}
			writer.write(seq.getId()+"\t"+pct9+"\t"+pct8+"\t"+pct7+"\n");
		}
		writer.close();
	}*/

	private List<Sequence> filter(List<Sequence> geneSequences,	Collection<String> genesOfInterest, Map<String, String> refLink) {
		List<Sequence> rtrn=new ArrayList<Sequence>();
		
		for(Sequence seq: geneSequences){
			String gene=refLink.get(seq.getId());
			if(gene==null){gene=refLink.get(seq.getId().split("\\.")[0]);}
			if(gene!=null){
				if(genesOfInterest.contains(gene) || genesOfInterest.contains(gene.toUpperCase()) || genesOfInterest.contains(gene.split("\\.")[0].toUpperCase())){
					System.err.println("Found: "+gene+" "+seq.getId());
					rtrn.add(seq);
				}
			}
		}
		
		return rtrn;
	}
	
	private List<Sequence> filter(List<Sequence> geneSequences,	Collection<String> genesOfInterest) {
		List<Sequence> rtrn=new ArrayList<Sequence>();
		
		for(Sequence seq: geneSequences){
			String gene=seq.getId();
			if(genesOfInterest.contains(gene) || genesOfInterest.contains(gene.toUpperCase()) || genesOfInterest.contains(gene.split("\\.")[0].toUpperCase())){
				rtrn.add(seq);
			}
		}
		
		return rtrn;
	}

	public static void main(String[] args)throws IOException{
		if(args.length>4){
			File proteinSeq=new File(args[0]);
			Map<String, String> refLink=ParseGCTFile.parseChipFile(new File(args[1]));
			Collection<String> genesOfIneterest=ParseGCTFile.parseGeneList(new File(args[2]));
			File genes=new File(args[3]);
			String save=args[4];
			new IdentifyHardToCloneSeqs(proteinSeq, refLink, genesOfIneterest, genes, save);
		}
		else{
			System.err.println(usage);
		}
	}
	
	static String usage=" args[0]=protein seq \n args[1]=refLink \n args[2]=genesOfInterest \n args[3]=gene fasta \n args[4]=save";
	
}
