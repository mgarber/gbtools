package broad.pda.seq.utils;

import java.io.File;
import java.io.IOException;

public class TophatAlignBatch {
	private static final String quote="\"";


	public static void main(String[] args)throws IOException{
		if(args.length>3){
		Runtime run=java.lang.Runtime.getRuntime();
		File[] files=new File(args[0]).listFiles();
		String referenceLocation=args[1];
		String saveDir=args[2];
		String queue=args[3];
		//String gffFile=args[4];
		
		
		for(int i=0; i<files.length; i++){
			String save=saveDir+"/"+files[i].getName().split("//.")[0];
			new File(save).mkdir();
			String rusage="-R "+"rusage[mem=3]";
			String command="bsub -q "+queue+" -o "+save+"/junk.bsub "+rusage+" /seq/mguttman/scripts/TopHat/bin/tophat -o "+save +" "+referenceLocation+" "+files[i];
			if(args.length>4){
				command="bsub -q "+queue+" -o "+save+"/junk.bsub"+" /seq/mguttman/scripts/TopHat/bin/tophat -o "+save +" -j "+args[4]+" "+referenceLocation+" "+files[i];
			}
			run.exec(command);
			System.err.println(command);
		}
		
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=sequence files \n args[1]=reference location \n args[2]=save dir \n args[3]=queue \n args[4]=Junctions File (optional)";
	
}
