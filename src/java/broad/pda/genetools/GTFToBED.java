package broad.pda.genetools;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.pda.annotation.GTFFileParser;
import broad.pda.gene.RefSeqGene;
import broad.pda.gene.RefSeqGeneWithIsoforms;

public class GTFToBED {

	
	public GTFToBED(String file, String save) throws IOException {
		GTFFileParser gtf=new GTFFileParser(file);
		Map<String, IntervalTree<RefSeqGeneWithIsoforms>> trees=gtf.getIntervalTreeWithIsoforoms();
		
		FileWriter writer=new FileWriter(save);
		
		for(String chr: trees.keySet()){
			Iterator<Node<RefSeqGeneWithIsoforms>> iter=trees.get(chr).iterator();
			while(iter.hasNext()){
				RefSeqGeneWithIsoforms gene=iter.next().getValue();
				Collection<RefSeqGene> isoforms=gene.getAllIsoforms();
				for(RefSeqGene isoform: isoforms){
					writer.write(isoform+"\n");
				}
			}
		}
		
		writer.close();
		
	}

	public static void main(String[] args)throws IOException{
		if(args.length>1){
			String file=(args[0]);
			String save=args[1];
			new GTFToBED(file, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=GTF file \n args[1]=save";
	
}
