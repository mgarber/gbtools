package broad.pda.seq.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class CleanupFile {

	public CleanupFile(File file, String save)throws IOException{
		readWriteAndClean(file, save);
	}
	
	private void readWriteAndClean(File file, String save)throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		FileWriter writer=new FileWriter(save);
    	  String nextLine;
    	  int i=0;
            while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
               
					if(nextLine.startsWith("chr")){
                	
                	String[] tokens=nextLine.split("\t");
					String chr=(tokens[0]);
					int start=new Integer(tokens[1]);
					int end=new Integer(tokens[2]);
					
					if(start<end){writer.write(nextLine+"\n");}
                }
				i++;
				if(i%10000 ==0){System.err.println(i);}
            }
            
            
            reader.close();
            writer.close();
            
	}
	
	public static void main(String[] args)throws IOException{
		File file=new File(args[0]);
		String save=args[1];
		new CleanupFile(file, save);
	}
	
}
