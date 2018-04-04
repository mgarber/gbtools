package broad.pda.arrays.tilling;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Map;
import java.util.TreeMap;

//Take an IGV file and generate all possible pairings

public class GenerateRatios {
	ArrayList<String> header;

	public GenerateRatios(File file, String saveDir)throws IOException{
		Map<String, Double>[] maps=parse(file);
		write(saveDir, maps);
	}
	
	
	private Map<String, Double>[] parse(File file)throws IOException{
		Map<String, Double>[] rtrn=null;
		this.header=new ArrayList();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		int counter=0;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			String[] tokens=nextLine.split("\t");
			String name=tokens[3];
			for(int i=4; i<tokens.length; i++){
				if(rtrn==null){
					rtrn=new Map[tokens.length-4];
					for(int j=0; j<rtrn.length; j++){rtrn[j]=new TreeMap();}
				}
				if(counter==0){header.add(tokens[i]);}
				else{
					Double val=new Double(tokens[i]);
					rtrn[i-4].put(name, val);
				}
			}
			counter++;
		}
		return rtrn;
	}
	
	private void write(String saveDir, Map<String, Double>[] maps)throws IOException{
		
		for(int i=0; i<maps.length; i++){
			for(int j=0; j<maps.length; j++){
				if(i!=j){
					System.err.println(header.get(i)+"_vs"+header.get(j)+".ratio");
					String save=saveDir+"/"+header.get(i)+"_vs_"+header.get(j)+".ratio";
					write(save, maps[i], maps[j]);					
				}
			}
		}
		
	}
	
	
	private void write(String save, Map<String, Double> map1, Map<String, Double> map2)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(String name: map1.keySet()){
			double val1=map1.get(name);
			double val2=map2.get(name);
			double ratio=val1-val2;
			writer.write(name+"\t"+ratio+"\n");
		}
		
		writer.close();
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>1){
		File file=new File(args[0]);
		String saveDir=args[1];
		new GenerateRatios(file, saveDir);
		}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=IGV File \n args[1]=save directory";
}
