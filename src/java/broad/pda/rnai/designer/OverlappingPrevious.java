package broad.pda.rnai.designer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.TreeSet;

public class OverlappingPrevious {

	
	public OverlappingPrevious(File[] files, String saveDir, File previousHairpinFile) throws IOException{
		Collection<String> previousHairpins=parseHairpinFile(previousHairpinFile);
		for(int i=0; i<files.length; i++){
			System.err.println(files[i]);
			readAndWriteNew(files[i], saveDir, previousHairpins);
		}
	}

	private void readAndWriteNew(File file, String saveDir,	Collection<String> previousHairpins) throws IOException {
		String save=saveDir+"/"+file.getName();
		FileWriter writer=new FileWriter(save);
	
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		int i=0;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			String[] tokens=nextLine.split("\t");
			String kmer=tokens[0];
			if(!previousHairpins.contains(kmer)){writer.write(nextLine+"\n");}
		}
		reader.close();
		
		
		writer.close();
	}

	private Collection<String> parseHairpinFile(File hairpinFile) throws IOException {
		Collection<String> rtrn=new TreeSet();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(hairpinFile)));
		String nextLine;
		int i=0;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			if(i>0){
				String[] tokens=nextLine.split("\t");
				String kmer=tokens[0];
				rtrn.add(kmer);
			}
			i++;
		}
		reader.close();
		
		return rtrn;
	}
	
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
		File[] files=new File(args[0]).listFiles();
		String saveDir=args[1];
		File previous=new File(args[2]);
		new OverlappingPrevious(files, saveDir, previous);
		}
		else{System.err.println(usage);}
	}

	static String usage=" args[0]=top files \n args[1]=save dir \n args[2]=previous hairpins";
	
}
