package broad.pda.ribosome.misc;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.parser.CommandLineParser;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.GeneTools;
import broad.pda.gene.RefSeqGene;

public class FindORFs {

	private Map<RefSeqGene, Collection<RefSeqGene>> allORFs;
	private Map<RefSeqGene, Collection<RefSeqGene>> allNonORFs;
	private Map<String, IntervalTree<RefSeqGene>> orfTree;
	
	//This will find all ORFs and regions that are not ORFs
	public FindORFs(Map<String, Collection<RefSeqGene>> genesByChr, String seqDir, boolean trim3UTRs) throws IOException{
		this.allORFs=findAllORFs(genesByChr, seqDir, trim3UTRs);
		this.allNonORFs=findAllNonORFs(allORFs, genesByChr);
		buildORFTree();
	}
	
	private void buildORFTree() {
		this.orfTree=new TreeMap<String, IntervalTree<RefSeqGene>>();
		
		for(RefSeqGene gene: allORFs.keySet()){
			for(RefSeqGene orf: allORFs.get(gene)){
				IntervalTree<RefSeqGene> tree=new IntervalTree<RefSeqGene>();
				if(this.orfTree.containsKey(orf.getChr())){tree=this.orfTree.get(orf.getChr());}
				tree.put(orf.getCDSRegion().getStart(), orf.getCDSRegion().getEnd(), orf);
				this.orfTree.put(orf.getChr(), tree);
			}
		}
		
	}

	public FindORFs(Collection<RefSeqGene> genes, String genomeDir, boolean trim3UTRs) throws IOException{
		Map<String, Collection<RefSeqGene>> genesByChr=new TreeMap<String, Collection<RefSeqGene>>();
		
		for(RefSeqGene gene: genes){
			Collection<RefSeqGene> list=new TreeSet<RefSeqGene>();
			if(genesByChr.containsKey(gene.getChr())){
				list=genesByChr.get(gene.getChr());
			}
			list.add(gene);
			genesByChr.put(gene.getChr(), list);
		}
		
		this.allORFs=findAllORFs(genesByChr, genomeDir, trim3UTRs);
		this.allNonORFs=findAllNonORFs(allORFs, genesByChr);
		buildORFTree();
	}
	
	public Map<RefSeqGene, Collection<RefSeqGene>> getAllORFs(){return this.allORFs;}
	public Map<RefSeqGene, Collection<RefSeqGene>> getAllNonORFs(){return this.allNonORFs;}

	private Map<RefSeqGene, Collection<RefSeqGene>> findAllORFs(Collection<RefSeqGene> genes, Sequence chrSeq, boolean trim3UTRs) throws IOException {
		Map<RefSeqGene, Collection<RefSeqGene>> allORFs=new TreeMap<RefSeqGene, Collection<RefSeqGene>>();
		for(RefSeqGene gene: genes){
			Collection<RefSeqGene> orfs=findAllORFs(gene, chrSeq, trim3UTRs);
			allORFs.put(gene, orfs);
		}
		return allORFs;
	}
	
	private Map<RefSeqGene, Collection<RefSeqGene>> findAllORFs(Map<String, Collection<RefSeqGene>> genesByChr, String seqDir, boolean trim3UTRs) throws IOException {
		Map<RefSeqGene, Collection<RefSeqGene>> allORFs=new TreeMap<RefSeqGene, Collection<RefSeqGene>>();
		
		
		for(String chr: genesByChr.keySet()){
			//System.err.println(chr);
			try{
			String chrSeqFile=seqDir+"/"+chr+".fa";
			FastaSequenceIO fsio = new FastaSequenceIO(chrSeqFile);
			Sequence chrSeq=fsio.loadAll().get(0);
			for(RefSeqGene gene: genesByChr.get(chr)){
				Collection<RefSeqGene> orfs=findAllORFs(gene, chrSeq, trim3UTRs);
				allORFs.put(gene, orfs);
			}
			}catch(FileNotFoundException ex){
				//ex.printStackTrace();
				System.err.println("Skipping "+chr);
			}
		}
		return allORFs;
	}

	private Map<RefSeqGene, Collection<RefSeqGene>> mergeORFs(Map<RefSeqGene, Collection<RefSeqGene>> allORFs, Map<String, Collection<RefSeqGene>> genesByChr) throws IOException {
		Map<RefSeqGene, Collection<RefSeqGene>> rtrn=new TreeMap<RefSeqGene, Collection<RefSeqGene>>();
		Map<String, Collection<RefSeqGene>> map=new TreeMap<String, Collection<RefSeqGene>>();
		
		for(RefSeqGene gene: allORFs.keySet()){
			Collection<RefSeqGene> orfs=allORFs.get(gene);
			for(RefSeqGene orf: orfs){
				String chr=orf.getChr();
				Collection<RefSeqGene> list=new TreeSet<RefSeqGene>();
				if(map.containsKey(chr)){
					list=map.get(chr);
				}
				//System.err.println("Putting " + gene.toBED());
				list.add(orf.getCDS());
				map.put(chr, list);
			}
		}
		
		Map<String, IntervalTree<RefSeqGene>> trees=GeneTools.merge(map);
		
		for(String chr: genesByChr.keySet()){
			IntervalTree<RefSeqGene> tree=trees.get(chr);
			if(tree!=null){
				for(RefSeqGene gene: genesByChr.get(chr)){
					Iterator<Node<RefSeqGene>> iter=tree.overlappers(gene.getStart(), gene.getEnd());
					Collection<RefSeqGene> list=new TreeSet<RefSeqGene>();
					while(iter.hasNext()){
						RefSeqGene next=iter.next().getValue();
						if(next.getName().equalsIgnoreCase(gene.getName())){list.add(next);}
					}
					rtrn.put(gene, list);
				}
			}
		}
		
		return rtrn;
	}
	
	

	/**
	 * @param gene
	 * @param orfs
	 * @param chrSeq
	 * @return
	 * @throws IOException 
	 */
	private Map<RefSeqGene, Collection<RefSeqGene>> findAllNonORFs(Map<RefSeqGene, Collection<RefSeqGene>> allORFs, Map<String, Collection<RefSeqGene>> genesByChr) throws IOException {
		double sum=0;
		double counter=0;
		
		Map<RefSeqGene, Collection<RefSeqGene>> rtrn=new TreeMap<RefSeqGene, Collection<RefSeqGene>>();
		
		//Step 1: First collapse all ORFs
		Map<RefSeqGene, Collection<RefSeqGene>> collapsedORFs=mergeORFs(allORFs, genesByChr);
		
		
		//Step 2: For each gene
		for(RefSeqGene gene: collapsedORFs.keySet()){
			RefSeqGene previousORF=null;
			Collection<RefSeqGene> list=new TreeSet<RefSeqGene>();
			for(RefSeqGene orf: collapsedORFs.get(gene)){
				RefSeqGene nonORF;
				if(previousORF==null){
					nonORF=gene.trimAbsolute(gene.getStart(), orf.getStart());
					//System.err.println(gene.getChr()+":"+gene.getStart()+"-"+orf.getStart());
				}
				else{
					nonORF=gene.trimAbsolute(previousORF.getEnd(), orf.getStart());
					//System.err.println(gene.getChr()+":"+gene.getStart()+"-"+orf.getStart());
				}
				previousORF=orf;
				if(nonORF!=null){
					RefSeqGene testCDS=gene.copy();
					testCDS.setCDS(nonORF.getAlignment());
					list.add(testCDS);
				}
			}
			
			if(previousORF!=null){
				RefSeqGene nonORF=gene.trimAbsolute(previousORF.getEnd(), gene.getEnd());
				if(nonORF!=null){
					RefSeqGene testCDS=gene.copy();
					testCDS.setCDS(nonORF.getAlignment());
					list.add(testCDS);
					//list.add(nonORF);
				}
			}
			rtrn.put(gene, list);
		}
		
		//System.err.println(rtrn.size());
		
		return rtrn;
	}
	
	
	private Map<RefSeqGene, Collection<RefSeqGene>> findAllNonORFs(Map<RefSeqGene, Collection<RefSeqGene>> allORFs, Collection<RefSeqGene> genes) throws IOException {
		double sum=0;
		double counter=0;
		
		Map<RefSeqGene, Collection<RefSeqGene>> rtrn=new TreeMap<RefSeqGene, Collection<RefSeqGene>>();
		
		
		Map<String, Collection<RefSeqGene>> genesByChr=new TreeMap<String, Collection<RefSeqGene>>();
		genesByChr.put(genes.iterator().next().getChr(), genes);
		
		//Step 1: First collapse all ORFs
		Map<RefSeqGene, Collection<RefSeqGene>> collapsedORFs=mergeORFs(allORFs, genesByChr);
		
		
		//Step 2: For each gene
		for(RefSeqGene gene: collapsedORFs.keySet()){
			RefSeqGene previousORF=null;
			Collection<RefSeqGene> list=new TreeSet<RefSeqGene>();
			for(RefSeqGene orf: collapsedORFs.get(gene)){
				RefSeqGene nonORF;
				if(previousORF==null){
					nonORF=gene.trimAbsolute(gene.getStart(), orf.getStart());
					//System.err.println(gene.getChr()+":"+gene.getStart()+"-"+orf.getStart());
				}
				else{
					nonORF=gene.trimAbsolute(previousORF.getEnd(), orf.getStart());
					//System.err.println(gene.getChr()+":"+gene.getStart()+"-"+orf.getStart());
				}
				previousORF=orf;
				if(nonORF!=null){list.add(nonORF);}
			}
			
			if(previousORF!=null){
				RefSeqGene nonORF=gene.trimAbsolute(previousORF.getEnd(), gene.getEnd());
				if(nonORF!=null){list.add(nonORF);}
			}
			rtrn.put(gene, list);
		}
		
		//System.err.println(rtrn.size());
		
		return rtrn;
	}

	private void write(Map<RefSeqGene, Collection<RefSeqGene>> allORFs, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(RefSeqGene gene: allORFs.keySet()){
			Collection<RefSeqGene> orfs=allORFs.get(gene);
			for(RefSeqGene orf:orfs){
				if(orf!=null){
					orf.setName(gene.getName());
					writer.write(orf.toString()+"\n");
				}
				else{System.err.println(gene.getName()+" is null");}
			}
		}
		
		writer.close();
	}

	public static Collection<RefSeqGene> findAllORFs(RefSeqGene gene, Sequence chrSequence, boolean trim3UTRs) {
		
		gene.setSequenceFromChromosome(chrSequence);
		Collection<RefSeqGene> rtrn=gene.findAllORFs(trim3UTRs);
		return rtrn;
	}
	
	public static void main(String[] args)throws IOException{
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-g", "Gene bed file", true);
		p.addStringArg("-s", "Sequence directory", true);
		p.addStringArg("-o", "Output file prefix", true);
		p.addBooleanArg("-u", "Trim 3'UTRs to beginning of next ORF", false, Boolean.valueOf(false));
		p.parse(args);
		
		Map<String, Collection<RefSeqGene>> genesByChr = BEDFileParser.loadDataByChr(new File(p.getStringArg("-g")));
		String seqDir = p.getStringArg("-s");
		String save = p.getStringArg("-o");
		boolean trimUTR = p.getBooleanArg("-u").booleanValue();
		if(trimUTR) save += "_3UTRsTrimmed";
		FindORFs orfFinder=new FindORFs(genesByChr, seqDir, trimUTR);
		orfFinder.write(orfFinder.allORFs, save + ".bed");
		orfFinder.write(orfFinder.allNonORFs, save+".nonORF.bed");

	}
	
	public boolean overlapsORF(RefSeqGene window) {
		if(this.orfTree.get(window.getChr()).overlappers(window.getStart(), window.getEnd()).hasNext()){return true;}
		return false;
	}
	
}
