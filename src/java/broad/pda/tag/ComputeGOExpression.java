package broad.pda.tag;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.pda.annotation.BEDFileParser;

public class ComputeGOExpression {

	public ComputeGOExpression(File genelist, File refLink, File expression, String save) throws IOException {
		Collection<String> genes=parseGeneList(genelist);
		Map<String, Collection<String>> geneToRefSeq=getGenetoRefSeq(refLink);
		Map<String, double[]> refSeqExpression=parseExpression(expression);
		
		write(genes, geneToRefSeq, refSeqExpression, save);
		
	}

	private Map<String, double[]> parseExpression(File expression) throws IOException {
		Map<String, double[]> rtrn=new TreeMap<String, double[]>();
		
		Collection<String> lines=BEDFileParser.loadList(expression.getAbsolutePath());
		
		for(String line: lines){
			String[] tokens=line.split("\t");
			//System.err.println(line+" "+tokens.length);
			String name=tokens[3];
			double[] vals=new double[tokens.length-12];
			for(int i=12; i<tokens.length; i++){
				vals[i-12]=new Double(tokens[i]);
			}
			rtrn.put(name, vals);
		}
		
		return rtrn;
	}

	private Map<String, Collection<String>> getGenetoRefSeq(File refLink) throws IOException {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		Collection<String> lines=BEDFileParser.loadList(refLink.getAbsolutePath(), true);
		
		for(String line: lines){
			String[] tokens=line.split("\t");
			String name=tokens[0];
			String refSeq=tokens[2];
			Collection<String> set=new TreeSet<String>();
			if(rtrn.containsKey(name)){
				set=rtrn.get(name);
			}
			set.add(refSeq);
			rtrn.put(name, set);
		}
		
		return rtrn;
	}

	private Collection<String> parseGeneList(File genelist) throws IOException {
		Collection<String> rtrn=new TreeSet<String>();
		
		Collection<String> lines=BEDFileParser.loadList(genelist.getAbsolutePath());
		
		for(String line: lines){
			String[] tokens=line.split("\t");
			rtrn.add(tokens[2]);
		}
		
		
		return rtrn;
	}

	private void write(Collection<String> genes, Map<String, Collection<String>> geneToRefSeq, Map<String, double[]> refSeqExpression, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String gene: genes){
			Collection<String> refSeqs=geneToRefSeq.get(gene);
			if(refSeqs!=null){
			for(String refSeq: refSeqs){
				double[] vals=refSeqExpression.get(refSeq);
				if(vals!=null){
				writer.write(refSeq+"\t"+gene);
				for(int i=0; i<vals.length; i++){writer.write("\t"+vals[i]);}
				writer.write("\n");
				}
				else{System.err.println("Cant find "+refSeq+" "+gene);}
			}
			}
		}
		
		writer.close();
	}

	public static void main(String[] args) throws IOException{
		if(args.length>3){
			File genelist=new File(args[0]);
			File refLink=new File(args[1]);
			File expression=new File(args[2]);
			String save=args[3];
			new ComputeGOExpression(genelist, refLink, expression, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=Gene list (GO) \n args[1]=ref link (name to RefSeq) \n args[2]=RNA-Seq expression \n args[3]=save";
	
}
