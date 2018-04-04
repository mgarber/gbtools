package broad.pda.miseq;

import java.io.File;
import java.io.IOException;

public class ToTDF {

	public static void main(String[] args) throws IOException, InterruptedException{
		if(args.length>3){
			Runtime run=Runtime.getRuntime();
			File[] files=new File(args[0]).listFiles();
			String saveDir=args[1];
			String extension=args[2];
			String genome=args[3];
			
			for(int i=0; i<files.length; i++){
				File file=files[i];
				if(files[i].getName().endsWith(extension)){
					String output=saveDir+"/"+file.getName()+".tdf";
					String command="/xchip/igv/tools/igvtools count "+file.getAbsolutePath()+" "+output+" "+" "+genome;
					System.err.println(command);
					Process p=run.exec(command);
					p.waitFor();
					//[inputFile] [outputFile] [genome]
				}
					
			}
			
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=directory \n args[1]=save dir \n args[2]=extension \n args[3]=genome";
	
}
