package broad.pda.ribosome.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

import broad.pda.seq.alignment.sam.SAMUtils;

public class TrimAndSort {

	public TrimAndSort(File sam, String save) throws IOException{
		FileWriter writer=new FileWriter(save); 
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(sam)));
		 String nextLine;
	     while ((nextLine = reader.readLine()) != null) {
	    	 if(SAMUtils.isValid(nextLine)){
	    		 writer.write(nextLine+"\n");
	    	 }
	     }
	     reader.close();
	     writer.close();
	}
	
	public static void main(String[] args)throws IOException{
		new TrimAndSort(new File(args[0]), args[1]);
	}
	
}
