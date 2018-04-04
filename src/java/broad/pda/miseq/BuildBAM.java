package broad.pda.miseq;

import java.io.File;
import java.io.IOException;

public class BuildBAM {

	String jar="java -jar /seq/mgarber/tools/picard-tools-1.66/BuildBamIndex.jar";
	
	public BuildBAM(File[] files) throws IOException, InterruptedException{
		Runtime run=Runtime.getRuntime();
		for(int i=0; i<files.length; i++){
			if(files[i].getAbsolutePath().endsWith(".bam")){
				System.err.println(files[i].getAbsolutePath());
				String cmd=jar+" I="+files[i].getAbsolutePath()+" O="+files[i].getAbsolutePath()+".bai";
				Process p=run.exec(cmd);
				p.waitFor();
			}
		}	
	}
	
	public static void main(String[] args)throws IOException, InterruptedException{
		File[] files=new File(args[0]).listFiles();
		new BuildBAM(files);
	}
	
	
}
