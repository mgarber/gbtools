package broad.pda.genetools;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.RefSeqGene;

public class FilterGenesByName {

	public FilterGenesByName(Collection<RefSeqGene> genes, Collection<String> names, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		for(RefSeqGene gene: genes){
			if(names.contains(gene.getName())){
				writer.write(gene+"\n");
			}
		}
		writer.close();
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>2){
			Collection<RefSeqGene> genes=BEDFileParser.loadData(new File(args[0]));
			Collection<String> names=BEDFileParser.loadList(args[1]);
			String save=args[2];
			new FilterGenesByName(genes, names, save);
		}
		else{
			System.err.println(usage);
		}
	}
	
	static String usage=" args[0]=full BED \n args[1]=names \n args[2]=save";
	
}
