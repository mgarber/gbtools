package broad.pda.miseq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;
import broad.pda.seq.fastq.FastqSequence;

public class AnalyzeBarcodes {

	public AnalyzeBarcodes(File[] fastqFiles, File barcodesUsed, String save) throws IOException{
		//Analyze barcode diversity
		
		Map<String, String> barcodeInfo=parseBarcodesUsed(barcodesUsed);
		
		write(save, barcodeInfo, fastqFiles);
		
		//Count ribosomal RNAs vs mRNA
		trimReads(fastqFiles);
	}
	
	private void trimReads(File[] fastqFiles) throws IOException {
		//go through each read and trim the barcodes
		for(int i=0; i<fastqFiles.length; i++){
			File file=fastqFiles[i];
			String save=fastqFiles[i].getAbsolutePath()+".trimmed.fq";
			trimReads(file, save);
		}
		
	}

	private void trimReads(File file, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	    	if(nextLine.startsWith("@")){
	    		FastqSequence seq=new FastqSequence(nextLine, reader.readLine(), reader.readLine(), reader.readLine());
	        	FastqSequence trimmed=seq.trimFirstNBPs(9);
	        	writer.write(trimmed+"\n");
	    	}
	    }
		writer.close();
	}

	private Map<String, String> parseBarcodesUsed(File barcodesUsed) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		Collection<String> lines=BEDFileParser.loadList(barcodesUsed.getAbsolutePath());
		
		for(String line: lines){
			String[] tokens=line.split("\t");
			rtrn.put(Sequence.reverseSequence(tokens[3]), line);
			//System.err.println(tokens[3]+" "+Sequence.reverseSequence(tokens[3]));
		}
		
		return rtrn;
	}

	private void write(String save, Map<String, String> barcodeInfo, File[] fastqFiles) throws IOException {
		Map<String, Integer>[] barcodeCount=countBarcodes(fastqFiles);
		
		FileWriter writer=new FileWriter(save);
		
		writer.write("Read\tPlate\tRow\tColumn\tBarcode");
		for(int i=0; i<fastqFiles.length; i++){
			writer.write("\t"+fastqFiles[i].getName());
		}
		writer.write("\n");
		
		
		Collection<String> barcodes=new TreeSet<String>();
		for(int i=0; i<barcodeCount.length; i++){
			barcodes.addAll(barcodeCount[i].keySet());
		}
				
		for(String barcode: barcodes){
			String info=barcodeInfo.get(barcode);
			writer.write(barcode);
			if(info!=null){writer.write("\t"+info);}
			else{writer.write("\tNA\tNA\tNA\tNA");}
			for(int i=0; i<barcodeCount.length; i++){
				int val=0;
				if(barcodeCount[i].containsKey(barcode)){val=barcodeCount[i].get(barcode);}
				int all=barcodeCount[i].get("all");
				double ratio=(double)val/(double)all;
				writer.write("\t"+ratio);
			}
			writer.write("\n");
		}
		
		writer.close();
	}

	private Map<String, Integer>[] countBarcodes(File[] fastqFiles) throws IOException{
		Map<String, Integer>[] rtrn=new Map[fastqFiles.length];
		for(int i=0; i<fastqFiles.length; i++){
			int all=0;
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fastqFiles[i])));
			String nextLine;
			Map<String, Integer> barcodeCounter=new TreeMap<String, Integer>();
		    while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
		    	if(nextLine.startsWith("@")){
		    		FastqSequence seq=new FastqSequence(nextLine, reader.readLine(), reader.readLine(), reader.readLine());
		        	String barcode=seq.getSequence().substring(0, 8);
		        	//System.err.println(barcode);
		        	int count=0;
		        	if(barcodeCounter.containsKey(barcode)){count=barcodeCounter.get(barcode);}
		        	count++;
		        	barcodeCounter.put(barcode, count);
		        	all++;
		    	}
		    }
		    barcodeCounter.put("all", all);
		    rtrn[i]=barcodeCounter;
		}
	   return rtrn;
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>2){
			File[] fastq=new File(args[0]).listFiles();
			File barcodes=new File(args[1]);
			String save=args[2];
			new AnalyzeBarcodes(fastq, barcodes, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=fastq file \n args[1]=barcodesUsed \n args[2]=save";
	
}
