package broad.pda.seq.alignment;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

import broad.core.annotation.PSL;

public class PSLToSAM {
	
	public static void convertPSLToSAM(String file, String save) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		FileWriter writer=new FileWriter(save);
		
    	String nextLine;
        while ((nextLine = reader.readLine()) != null) {
        	PSL psl= PSL.create(nextLine);      
        	if(psl != null) {
        		writer.write(psl.toSAM()+"\n");
        	}
        }
        reader.close();
        writer.close();
	}
	
	public static void main(String[] args) throws IOException{
		convertPSLToSAM(args[0], args[1]);
	}

}
