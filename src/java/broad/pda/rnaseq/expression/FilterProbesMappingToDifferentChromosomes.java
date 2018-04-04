package broad.pda.rnaseq.expression;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.util.ParseGCTFile;
import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;

public class FilterProbesMappingToDifferentChromosomes {

	//An attempt to cleanup the bad annotations on affy that are causing inconsistencies with genome aligned data
	public FilterProbesMappingToDifferentChromosomes(File chipFile, File mappingFile, String save)throws IOException{
		Map<String, String> chipMap=ParseGCTFile.parseChipFile(chipFile);
		Map<String, RefSeqGene> geneMap=BEDFileParser.loadDataByName(mappingFile);
		
		//find all genes with more than PID
		//map each of these PIDs
		Map<String, Collection<RefSeqGene>> geneToProbe=getGeneToProbeMap(chipMap, geneMap);
		
		//flag PIDs for the same gene that map to different chr
		writeFilteredMapping(save, geneToProbe);		
		
		//write(save, geneToProbe);
		
	}
	
	private void writeFilteredMapping(String save, Map<String, Collection<RefSeqGene>> geneToProbe) throws IOException{
		FileWriter writer=new FileWriter(save);	
		
		for(String geneName: geneToProbe.keySet()){
			Collection<RefSeqGene> genes=geneToProbe.get(geneName);
			//writer.write(geneName+"\t"+genes.size()+"\n");
			boolean write=true;
			if(genes.size()>1){
				String chr=genes.iterator().next().getChr();
				for(RefSeqGene gene: genes){
					if(!gene.getChr().equalsIgnoreCase(chr)){write=false;}
				}
			}
			if(write){write(writer, geneName, genes);}
			else{Collection<RefSeqGene> filtered=filter(genes); write(writer, geneName, filtered);}
		}
		
		writer.close();
	}
	
	
	
	private Collection<RefSeqGene> filter(Collection<RefSeqGene> genes) {
		Collection<RefSeqGene> rtrn=new TreeSet();
		
		// TODO Auto-generated method stub
		return rtrn;
	}

	private void write(String save, Map<String, Collection<RefSeqGene>> geneToProbe) throws IOException{
		FileWriter writer=new FileWriter(save);	
		for(String geneName: geneToProbe.keySet()){
			Collection<RefSeqGene> genes=geneToProbe.get(geneName);
			//writer.write(geneName+"\t"+genes.size()+"\n");
			boolean write=false;
			if(genes.size()>1){
				String chr=genes.iterator().next().getChr();
				for(RefSeqGene gene: genes){
					if(!gene.getChr().equalsIgnoreCase(chr)){write=true;}
				}
			}
			if(write){write(writer, geneName, genes);}
		}
		writer.close();
	}
	
	private void write(FileWriter writer, String name, Collection<RefSeqGene> genes) throws IOException{
		
		for(RefSeqGene gene: genes){
			writer.write(gene+"\n");
		}
		
	}

	private Map<String, Collection<RefSeqGene>> getGeneToProbeMap(Map<String, String> chipMap, Map<String, RefSeqGene> geneMap) {
		Map<String, Collection<RefSeqGene>> rtrn=new TreeMap();
		
		for(String probe: chipMap.keySet()){
			String gene=chipMap.get(probe);
			RefSeqGene geneAlign=geneMap.get(probe);
			if(geneAlign!=null){
			Collection<RefSeqGene> set=new TreeSet();
			if(rtrn.containsKey(gene)){set=(Collection<RefSeqGene>)rtrn.get(gene);}
			set.add(geneAlign);
			rtrn.put(gene, set);
			}
			//else{System.err.println(probe+" "+gene);}
		}
		
		return rtrn;
	}

	public static void main(String[] args)throws IOException{
		if(args.length>2){
			File chipFile=new File(args[0]);
			File mappingFile=new File(args[1]);
			String save=args[2];
			new FilterProbesMappingToDifferentChromosomes(chipFile, mappingFile, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=chip file \n args[1]=mapping file (BED File) \n args[2]=save";
	
}
