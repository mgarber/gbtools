package broad.pda.seq.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.TreeSet;

public class FilterSAMFile {

	public FilterSAMFile(File samFile, Collection<String> multimappers, String save)throws IOException{
		FileWriter writer=new FileWriter(save);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(samFile)));
		String nextLine;
		int i=0;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	String name=nextLine.split("\t")[0];
        	if(!multimappers.contains(name)){writer.write(nextLine+"\n");}
        	if(i%10000 ==0){System.err.println(i);}
        	i++;
        }
		writer.close();
	}
	
	private static Collection<String> parseMultimappers(String file) throws IOException{
		Collection<String> rtrn=new TreeSet();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	rtrn.add(nextLine);
        }
				
		return rtrn;
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
			File samFile=new File(args[0]);
			String save=args[1];
			Collection<String> multimappers=parseMultimappers(args[2]);
			new FilterSAMFile(samFile, multimappers, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage="v2 \n args[0]=samFile \n args[1]=save \n args[2]=multimappers";	
	
}
