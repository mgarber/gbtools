package broad.pda.seq.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

import broad.pda.datastructures.Alignments;

public class AlignedToSAM {

	public AlignedToSAM(File[] files, String save, boolean flipSign)throws IOException{
		
		for(int i=0; i<files.length; i++){
			readAndWrite(files[i], save, flipSign);
		}
		
	}
	
	private void readAndWrite(File file, String save, boolean flipSign) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		FileWriter writer=new FileWriter(save, true);
    	String nextLine;
       
       	int i=0;
	   while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	        	String[] tokens=nextLine.split("\t");
	        	Alignments align=new Alignments(tokens[0], tokens[1], tokens[2]);
	        	align.setOrientation(tokens[3]);
	        	if(flipSign){align.setOrientation(flip(align.getOrientation()));}
	        	writer.write(align.toSAM()+"\n");
	        	i++;
	        	if(i%1000000 ==0){System.err.println(i);}
	   }
	            
	            
	         reader.close();
	         writer.close();

	}

	private String flip(String s) {
		if(s.equalsIgnoreCase("+")){return "-";}
		if(s.equalsIgnoreCase("-")){return "+";}
		return "*";
	}

	public static void main(String[] args)throws IOException{
		if(args.length>2){
			File[] files=new File(args[0]).listFiles();
			String save=args[1];
			boolean flipSign=new Boolean(args[2]);
			new AlignedToSAM(files, save, flipSign);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=files \n args[1]=save \n args[2]=flip sign?";
	
}
