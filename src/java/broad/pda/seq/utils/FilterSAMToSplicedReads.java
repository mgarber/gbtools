package broad.pda.seq.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;


//This class will take a "Generic SAM" file and remove all reads that overlap the same position in the genome
public class FilterSAMToSplicedReads {

	public FilterSAMToSplicedReads(File SAMFile, String save)throws IOException{
	 	filterAndWrite(save, SAMFile);
	}
	
	

	private void filterAndWrite(String save, File SAMFile)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(SAMFile)));
		String nextLine;
        int i=0;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	String cigar=nextLine.split("\t")[5];
        	if(cigar.split("N").length>1){writer.write(nextLine+"\n");}
        	i++;
        	if(i% 100000 ==0){System.err.println(i);}
		}
		
		
		writer.close();
	}
	
	public static void main(String[] args)throws IOException{
		File file=new File(args[0]);
		String save=args[1];
		new FilterSAMToSplicedReads(file, save);
	}
	

	
}
