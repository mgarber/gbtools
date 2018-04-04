package broad.pda.seq.utils;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import broad.pda.annotation.BEDFileParser;

public class BatchSubmit {

	public static void main(String[] args)throws IOException{
		if(args.length>12){
		Runtime run=java.lang.Runtime.getRuntime();
		File file=new File(args[0]);
		String maskFiles=args[1];
		File saveFile=new File(args[2]);
		String sizes=args[3];
		int windowSize=new Integer(args[4]);
		String queue=args[5];
		String script=args[6];
		String memory=args[7];
		File junkFile=new File(args[8]);
		String genomeDir=args[9];
		String paired=args[10];
		boolean upweight=new Boolean(args[11]);
		boolean hasPairs=new Boolean(args[12]);
		
		Map<String, Integer> sizeMap=BEDFileParser.loadChrSizes(sizes);
		
		for(String chr: sizeMap.keySet()){
			String chrSeq=genomeDir+"/"+chr.replaceAll("chr", "")+"/"+chr+".fa";
			String save=saveFile.getParentFile().getAbsolutePath()+"/"+chr+"."+saveFile.getName();
			String junk=junkFile.getParentFile().getAbsolutePath()+"/"+chr+"."+junkFile.getName();
			//TODO add upweight param
			String command="bsub -q "+queue+" -o "+ junk+" java -jar "+ memory+" "+script+" -alignment "+file.getAbsolutePath()+" -maskFileDir "+maskFiles+" -out "+save+" -sizeFile "+sizes+" -windows "+windowSize+" -chr "+chr +" -chrSequence "+chrSeq;
			if(hasPairs){command+=" -pairedEnd "+paired;}
			if(upweight){command+=" -upWeightSplices";}
			run.exec(command);
			System.err.println(command);
			System.gc();
		}
		
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=BAM File \n args[1]=mask file \n args[2]=save File \n args[3]=sizes \n args[4]=window size \n args[5]=queue \n args[6]=script \n args[7]=memory \n args[8]=junkFile \n args[9]=genome directory \n args[10]=paired ends \n args[11]=upweight splices \n args[12]=has Pairs";
	
}
