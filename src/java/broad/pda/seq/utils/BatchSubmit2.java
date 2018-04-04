package broad.pda.seq.utils;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import broad.pda.annotation.BEDFileParser;

public class BatchSubmit2 {

	/*
	public static void main(String[] args)throws IOException{
		if(args.length>7){
		Runtime run=java.lang.Runtime.getRuntime();
		File file=new File(args[0]);
		String maskFiles=args[1];
		File saveFile=new File(args[2]);
		String sizes=args[3];
		String queue=args[4];
		String script=args[5];
		String memory=args[6];
		File junkFile=new File(args[7]);
		String genomeDir=args[8];
		
		
		Map<String, Integer> sizeMap=BEDFileParser.loadChrSizes(sizes);
		
		for(String chr: sizeMap.keySet()){
			String chrSeq=genomeDir+"/"+chr.replaceAll("chr", "")+"/"+chr+".fa";
			String save=saveFile.getParentFile().getAbsolutePath()+"/"+chr+"."+saveFile.getName();
			String junk=junkFile.getParentFile().getAbsolutePath()+"/"+chr+"."+junkFile.getName();
			String command="bsub -q "+queue+" -o "+ junk+" java -jar "+ memory+" "+script+" "+file.getAbsolutePath()+" "+maskFiles+" "+save+" "+sizes+" "+chr+" "+chrSeq;
			run.exec(command);
			System.err.println(command);
			System.gc();
		}
		
		}
		else{System.err.println(usage);}
	}
	
	
	static String usage=" args[0]=BAM File \n args[1]=mask file \n args[2]=save File \n args[3]=sizes \n args[4]=queue \n args[5]=script \n args[6]=memory \n args[7]=junkFile";
	*/
	
	
	static String usage=" args[0]=Alignment file \n args[1]=mask file \n args[2]=saveDir \n args[3]=sizes \n args[4]=window sizes (comma-seperated) \n args[5]=alpha \n args[6]=trim ends \n args[7]=queue \n args[8]=script \n args[9]=memory \n args[10]=junkFile";
	
	public static void main(String[] args)throws IOException{
		if(args.length>10){
		Runtime run=java.lang.Runtime.getRuntime();
		File BAMFile=new File(args[0]);
		String maskFiles=args[1];
		File saveFile=new File(args[2]);
		String sizes=args[3];
		String fixedWidths=(args[4]);
		double alpha=new Double(args[5]);
		boolean trimEnds=new Boolean(args[6]);
		
		String queue=args[7];
		String script=args[8];
		String memory=args[9];
		File junkFile=new File(args[10]);
		
		
		Map<String, Integer> sizeMap=BEDFileParser.loadChrSizes(sizes);
		
		for(String chr: sizeMap.keySet()){
			String save=saveFile.getParentFile().getAbsolutePath()+"/"+chr+"."+saveFile.getName();
			String junk=junkFile.getParentFile().getAbsolutePath()+"/"+chr+"."+junkFile.getName();
			String command="bsub -q "+queue+" -o "+ junk+" java -jar "+ memory+" "+script+" "+BAMFile.getAbsolutePath()+" "+maskFiles+" "+save+" "+sizes+" "+fixedWidths+" "+alpha+" "+trimEnds+" "+chr ;
			run.exec(command);
			System.err.println(command);
			System.gc();
		}
		
		}
		else{System.err.println(usage);}
	}
	
}
