package broad.pda.seq.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class AppendOrientation {

	public AppendOrientation(File file, String save)throws IOException{
		appendAndWrite(file, save);
		
	}
	
	private void appendAndWrite(File file, String save)throws IOException{
		FileWriter writer=new FileWriter(save);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
    	int i=0;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	writer.write(nextLine+"\t+"+"\n");
        	if(i%100000 ==0){System.err.println(i);}
        	i++;
        }
		writer.close();
		
	}
	
	public static void main(String[] args)throws IOException{
		File file=new File(args[0]);
		String save=args[1];
		new AppendOrientation(file, save);
	}
}
