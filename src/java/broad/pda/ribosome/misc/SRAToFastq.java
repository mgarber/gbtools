package broad.pda.ribosome.misc;

import java.io.File;
import java.io.IOException;

public class SRAToFastq {

	public SRAToFastq(String SRADir, String saveDir) throws IOException, InterruptedException{
		Runtime run=Runtime.getRuntime();
		File[] files=new File(SRADir).listFiles();
		for(int i=0; i<files.length; i++){
			File actual=files[i].listFiles()[0];
			
			//convert to fastq
			///seq/mgarber/tools/sratoolkit.2.1.7-centos_linux64/fastq-dump SRR315591.sra
			System.err.println(actual);
			String cmd="/seq/mgarber/tools/sratoolkit.2.1.7-centos_linux64/fastq-dump -O "+saveDir+" "+actual.getAbsolutePath();
			Process p=run.exec(cmd);
			p.waitFor();
		}
		
	}
	
	public static void main(String[] args)throws IOException, InterruptedException{
		if(args.length>1){
			String SRADir=args[0];
			String saveDir=args[1];
			new SRAToFastq(SRADir, saveDir);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=SRA files \n args[1]=saveDir";
	
}
