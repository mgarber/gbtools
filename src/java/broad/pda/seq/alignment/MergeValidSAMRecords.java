package broad.pda.seq.alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class MergeValidSAMRecords {

	public MergeValidSAMRecords(File[] SAMFiles, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(int i=0; i<SAMFiles.length; i++){
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(SAMFiles[i])));
			
			int counter=0;
			int count=0;
			String nextLine;
			while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
				if(!nextLine.startsWith("@")){
					//try{
					String[] tokens=nextLine.split("\t");
					String match=tokens[2];
					if(!match.equalsIgnoreCase("*")){
						writer.write(nextLine+"\n");
						count++;
					}
					//}catch(Exception ex){}
				}
				counter++;
				if(counter%100000 ==0){System.err.println(counter+" Good: "+count);}
			}
						
			//parse SAM file
			//if valid SAM line write it
			//else skip it
		}
		
		writer.close();
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>1){
			File file=new File(args[0]);
			String save=args[1];
			File[] files={file};
			new MergeValidSAMRecords(files, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=sam files \n args[1]=save";
	
}
