package broad.pda.genetools;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import broad.pda.annotation.BEDFileParser;

public class CleanBEDFile {

	
	public static void main(String[] args)throws IOException{
		if(args.length>1){
			Collection<String> lines=BEDFileParser.loadList(args[0]);
			FileWriter writer=new FileWriter(args[1]);
			
			for(String line: lines){
				writer.write(line.replaceAll("\"", "")+"\n");
			}
			writer.close();
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=BED file \n args[1]=output";
	
}
