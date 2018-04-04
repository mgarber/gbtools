package broad.pda.seq.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class MakeValidFastq {

	public MakeValidFastq(File file, String save) throws IOException{
		FileWriter writer =new FileWriter(save);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		boolean started=false;
		int counter=0;
        while ((nextLine = reader.readLine()) != null) {
        	if(counter %1000000 ==0){System.err.println(counter);}
        	if(nextLine.startsWith("@") || started){started=true; writer.write(nextLine+"\n");}
        	counter++;
        }
		
        reader.close();
		writer.close();
	}
	
	public static void main(String[] args)throws IOException{
		File file=new File(args[0]);
		String save=args[1];
		
		//for(int i=0; i<files.length; i++){
			//System.err.println(files[i]);
			//String save=saveDir+"/"+files[i].getName()+".valid.fq";
			new MakeValidFastq(file, save);
		//}
	}
	
}
