package broad.pda.seq.alignment;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

import broad.core.annotation.PSL;
import broad.pda.seq.alignment.sam.SAMRecord;

public class GetUniquePSLOnly {

	//assume the alignments are sorted
	//iterate through the alignments
	//retain the previous one
	//if the current is the same as the previous skip it
	//if different then write the previous and update the previous
	public static void getUniquePSLOnly(String string, String save)throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(string)));
		FileWriter writer=new FileWriter(save);
		
		PSL previousPSL=null;
		boolean flagged=false;
		
		String nextLine;
        while ((nextLine = reader.readLine()) != null) {
        	try{
        		PSL psl=new PSL(nextLine);              	
				if(previousPSL!=null && !psl.getName().equalsIgnoreCase(previousPSL.getName())){
					//write previous PSL
					if(!flagged){writer.write(previousPSL.toPSL()+"\n");}
					flagged=false;
				}
				else{flagged=true;}
				previousPSL=psl;
        	}catch(Exception ex){System.err.println("Skipping: "+nextLine); ex.printStackTrace();}
				
          }
        writer.close();
	}
	
	
	public static void main(String[] args)throws IOException{
		if(args.length>1){
			String pslFile=(args[0]);
			String save=args[1];
			uniqueToSAM(pslFile, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=PSL file \n args[1]=save";

	public static void uniqueToSAM(String string, String save) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(string)));
		FileWriter writer=new FileWriter(save);
		
		PSL previousPSL=null;
		boolean flagged=false;
		
		String nextLine;
        while ((nextLine = reader.readLine()) != null) {
        	try{
        		PSL psl=new PSL(nextLine);              	
				if(previousPSL!=null && !psl.getName().equalsIgnoreCase(previousPSL.getName())){
					//write previous PSL
					if(!flagged){
						String name=psl.getName().replaceAll("@", "");
						SAMRecord record=new SAMRecord(name, psl.toGene(), psl.getSequence());
						record.setNumMismatches(psl.getMismatches());
						writer.write(record.toString()+"\n");}
					flagged=false;
				}
				else{flagged=true;}
				previousPSL=psl;
        	}catch(NumberFormatException ex){System.err.println("Skipping: "+nextLine); ex.printStackTrace();}
				
          }
        writer.close();
		
	}
	
}
