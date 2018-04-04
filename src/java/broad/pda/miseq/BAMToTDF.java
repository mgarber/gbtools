package broad.pda.miseq;

import java.io.File;

public class BAMToTDF {

	public BAMToTDF(File[] bamFiles, String saveDir){
		for(int i=0; i<bamFiles.length; i++){
			if(bamFiles[i].getName().endsWith("bam")){
				String command=" /xchip/igv/tool/igvtools count ";
			}
		}
		
	}
	
	public static void main(String[] args){
		if(args.length>1){
			File[] bamFiles=new File(args[0]).listFiles();
			String saveDir=args[1];
			new BAMToTDF(bamFiles, saveDir);
		}
		else{
			System.err.println(usage);
		}
	}
	
	static String usage=" args[0]=bamFiles \n args[1]=save directory";
	
}
