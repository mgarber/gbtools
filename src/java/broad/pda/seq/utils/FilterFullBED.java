package broad.pda.seq.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;

import broad.pda.gene.RefSeqGene;

public class FilterFullBED {

	public FilterFullBED(File file, String save, int minNumExons)throws IOException{
		//Collection<RefSeqGene> genes=BEDFileParser.loadData(file);
		write(file, save, minNumExons);
	}
	
	private void write(File file, String save, int minNumExons)throws IOException{
		FileWriter writer=new FileWriter(save);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
            
			if(nextLine.startsWith("chr")){
        		String[] tokens=nextLine.split("\t");
        		int numExons=new Integer(tokens[9]);
				if(numExons>=minNumExons){writer.write(nextLine+"\n");}
			}
		
		}
    
    
		reader.close();
		writer.close();
	}
	
	private void write(String save, Collection<RefSeqGene> genes, int minNumExons)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(RefSeqGene gene: genes){
			int numExons=gene.getNumExons();
			if(numExons>=minNumExons){writer.write(gene.toBED()+"\n");}
		}
		
		writer.close();
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
		File bedFile=new File(args[0]);
		String save=args[1];
		int minNumExons=new Integer(args[2]);
		new FilterFullBED(bedFile, save, minNumExons);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=BED file \n args[1]=save file \n args[2]=min num exons";
	
}
