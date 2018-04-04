package broad.pda.rnaseq.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Map;
import java.util.TreeMap;

import broad.pda.gene.RefSeqGene;

public class MakeExpressionTable {

	public MakeExpressionTable(File[] files, String save) throws IOException{
		
		Map<RefSeqGene, double[]>[] maps=new Map[files.length];
		
		for(int i=0; i<files.length; i++){
			maps[i]=parse(files[i]);
		}
		
		write(save, maps, files);
	}

	private Map<RefSeqGene, double[]> parse(File file) throws IOException {
		Map<RefSeqGene, double[]> rtrn=new TreeMap<RefSeqGene, double[]>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		
    	String nextLine;
        int i=0;
        while ((nextLine = reader.readLine()) != null ) {
        	RefSeqGene track=new RefSeqGene(nextLine, false);
        	String[] tokens=nextLine.split("\t");
        	double[] array={new Double(tokens[12]), new Double(tokens[13])};
        	rtrn.put(track, array);
        }
            
            
       reader.close();
           
		
		
		return rtrn;
	}

	private void write(String save, Map<RefSeqGene, double[]>[] maps, File[] files) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("chr\tstart\tend\tname\tscore\tstrand\tCDS start\tCDS end\tRGB\tnum exons\texon sizes\texon starts");
		for(int i=0; i<files.length; i++){
			writer.write("\t"+files[i].getName()+"_pvalue");
		}
		for(int i=0; i<files.length; i++){
			writer.write("\t"+files[i].getName()+"_enrichment");
		}
		writer.write("\n");
		
		for(RefSeqGene gene: maps[0].keySet()){
			writer.write(gene.toBED());
			for(int i=0; i<maps.length; i++){
				double[] vals=maps[i].get(gene);
				writer.write("\t"+vals[0]);
			}
			for(int i=0; i<maps.length; i++){
				double[] vals=maps[i].get(gene);
				writer.write("\t"+vals[1]);
			}
			writer.write("\n");
		}
		
		writer.close();
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>1){
			File[] files=new File(args[0]).listFiles();
			String save=args[1];
			new MakeExpressionTable(files, save);
		}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=files \n args[1]=save";
}
