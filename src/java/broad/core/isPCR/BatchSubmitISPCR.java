package broad.core.isPCR;

import java.io.File;
import java.io.IOException;

public class BatchSubmitISPCR {

	
	public static void main(String[] args)throws IOException{
		if(args.length>3){
		File file=new File(args[0]);
		String genomeDirectory=args[1];
		String saveDir=args[2];
		String queue=args[3];
		
		Runtime run=java.lang.Runtime.getRuntime();
		File[] dirs=new File(genomeDirectory).listFiles();
		for(int i=0; i<dirs.length; i++){
			String chr="chr"+dirs[i].getName();
			String faFile=dirs[i].getAbsolutePath()+"/"+chr+".fa";
			String save=saveDir+"/"+chr+"."+file.getName()+".isPCR";
			String command="bsub -q "+queue+" /seq/mguttman/scripts/isPCR/isPcr "+faFile+" "+file.getAbsolutePath()+" "+save;
			run.exec(command);
			System.err.println(command);
		}
		}
		else{System.err.println(usage);}
		
	}
	
	static String usage=" args[0]=Primer file \n args[1]=genome directory \n args[2]=save directory \n args[3]=queue";
	
}
