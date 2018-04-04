package broad.pda.ribosome;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;

public class NovelSmallORFs {
	int minORF=10;
	int maxORF=100;

	/*public NovelSmallORFs(File testSetFile, File proteinFile, File nonCodingFile, String genomeDirectory, String save) throws IOException{
		Map<String, Collection<RefSeqGene>> testSetByChr=BEDFileParser.loadDataByChr(testSetFile);
		Map<String, Collection<RefSeqGene>> proteinsByChr=BEDFileParser.loadDataByChr(proteinFile);
		Map<String, Collection<RefSeqGene>> ncRNAByChr=BEDFileParser.loadDataByChr(nonCodingFile);
		
		//Step 1: generate Distribution for ncRNA
		
		EmpiricalDistribution negativeDist=generateDistributionForAllORFs(ncRNAByChr, genomeDirectory);
		
		//Step 2: generate Distribution for proteins
		
		EmpiricalDistribution positiveDist=generateDistributionForDefinedORFs(proteinsByChr, genomeDirectory);
		
		//Step 3: score test set
		assignScore(testSetByChr, negativeDist, positiveDist, genomeDirectory);
		
	}*/
	
	public NovelSmallORFs(File bedFile, String genomeDirectory, String save) throws IOException{
		Map<String, Collection<RefSeqGene>> genesByChr=BEDFileParser.loadDataByChr(bedFile);
		Map<RefSeqGene, Collection<RefSeqGene>> orfs=getAllORFs(genesByChr, genomeDirectory);
		
		Collection<RefSeqGene> candidateORFs=new TreeSet<RefSeqGene>();
		
		for(RefSeqGene transcript: orfs.keySet()){
			for(RefSeqGene orf: orfs.get(transcript)){
				int orfSize=orf.getSize();
				int aa=orfSize/3;
				if(aa>minORF && aa<maxORF){candidateORFs.add(orf);}
			}
		}
		
		//score candidate ORFs
		
		
		//write ORFs
		write(save, candidateORFs);
	}
	
	private void write(String save, Collection<RefSeqGene> candidateORFs) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(RefSeqGene orf: candidateORFs){
			writer.write(orf+"\n");
		}
		
		writer.close();
	}

	/*private EmpiricalDistribution generateDistributionForAllORFs(Map<String, Collection<RefSeqGene>> ncRNAByChr, String genomeDirectory) {
		Map<RefSeqGene, Collection<RefSeqGene>> orfs=getAllORFs(ncRNAByChr, genomeDirectory);
		
		File temp=new File("temp.bed");
		temp.deleteOnExit();
		Runtime run=Runtime.getRuntime();
		
		//TODO Make sure that if ORFs overlap we only take max
		for(RefSeqGene gene: orfs.keySet()){
			Collection<RefSeqGene> orf=orfs.get(gene);
			for(RefSeqGene orfRegion: orf){
				//write to temp
				orfRegion.writeBED(temp);
				String command=
			}
		}
		
		
		for(String chr: genesByChr.keySet()){
			String chrSeqFile=genomeDirectory+"/"+chr.replaceAll("chr", "")+"/"+chr+".fa";
			FastaSequenceIO fsio = new FastaSequenceIO(chrSeqFile);
			Sequence chrSeq=fsio.loadAll().get(0);
			
			
			Collection<RefSeqGene> genes=genesByChr.get(chr);
			
			for(RefSeqGene gene: genes){
				gene.setSequenceFromChromosome(chrSeq);
				Collection<RefSeqGene> rtrn=gene.findAllORFs();
				
				//TODO Score CSF for each ORF
			}
			
		}
		
		
	}*/

	private Map<RefSeqGene, Collection<RefSeqGene>> getAllORFs(Map<String, Collection<RefSeqGene>> genesByChr,	String genomeDirectory) throws IOException {
		Map<RefSeqGene, Collection<RefSeqGene>> orfs=new TreeMap<RefSeqGene, Collection<RefSeqGene>>();
		for(String chr: genesByChr.keySet()){
			Sequence chrSequence=makeChrSequence(chr, genomeDirectory);
			Collection<RefSeqGene> genes=genesByChr.get(chr);
			for(RefSeqGene gene: genes){
				gene.setSequenceFromChromosome(chrSequence);
				Collection<RefSeqGene> rtrn=gene.findAllORFs(false);
				orfs.put(gene, rtrn);
			}
		}
		return orfs;
	}

	private Sequence makeChrSequence(String chr, String genomeDirectory) throws IOException {
		String chrSeqFile=genomeDirectory+"/"+chr.replaceAll("chr", "")+"/"+chr+".fa";
		FastaSequenceIO fsio = new FastaSequenceIO(chrSeqFile);
		Sequence chrSeq=fsio.loadAll().get(0);
		return chrSeq;
	}

	public static void main(String[] args)throws IOException{
		if(args.length>2){
			File bedFile=new File(args[0]);
			String genomeDir=args[1];
			String save=args[2];
			new NovelSmallORFs(bedFile, genomeDir, save);
		}
		else{
			System.err.println(usage);
		}
	}
	
	static String usage=" args[0]=bedFile \n args[1]=genomeDirectory \n args[2]=save";
	
}
