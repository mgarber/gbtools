package broad.pda.genetools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;

public class CountStrandMismatch {

	public CountStrandMismatch(File samFile) throws IOException{
		 BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(samFile)));
		 String nextLine;
		 double plus=0;
		 double minus=0;
		 double endPos=0;
		 double total=0;
		 while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String strand=tokens[1];
			if(strand.equalsIgnoreCase("0")){plus++;}
			else if(strand.equalsIgnoreCase("16")){minus++;}
			if(tokens.length>3 && (tokens[3].equalsIgnoreCase("115") || tokens[3].equalsIgnoreCase("114") || tokens[3].equalsIgnoreCase("113") || tokens[3].equalsIgnoreCase("116"))){endPos++;}
			total++;
		 }
		 reader.close();
		 System.err.println(plus+" "+minus+" "+total+" "+endPos);
	}
	
	public static void main(String[] args)throws IOException{
		File samFile=new File(args[0]);
		new CountStrandMismatch(samFile);
	}
	
}
