package broad.projection.geo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Map;
import java.util.TreeMap;

public class FTPRetrieveFiles {

	public FTPRetrieveFiles(File file, String saveDir) throws IOException, InterruptedException{
		Map<String, String> map=parseGEO(file);
		retrieveViaFTP(map, saveDir);	
		//write(save, map);
	}
		
	private void write(String save, Map<String, String> map) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String nextLine: map.keySet()){
			writer.write(nextLine+"\n");
		}
		
		writer.close();
	}


	private void retrieveViaFTP(Map<String, String> map, String saveDir) throws IOException, InterruptedException {
		Runtime run=Runtime.getRuntime();
		
		for(String key: map.keySet()){
			String ftpAddress=map.get(key);
			String command="bsub -q priority -N wget -P "+saveDir+" "+ftpAddress;
			System.err.println(command);
			Process p=run.exec(command);
			//writeInputStream(p.getInputStream());
		}
		
	}

	
	private static void writeInputStream(InputStream stream, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		BufferedReader reader=	new BufferedReader(new InputStreamReader(stream));
		String rtrn="";
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			System.err.println(nextLine);
			writer.write(nextLine+"\n");
			rtrn+=nextLine;
		}
		writer.close();
	}

	private Map<String, String> parseGEO(File file) throws IOException {
		Map<String, String> rtrn=new TreeMap();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		
    	String nextLine;
    	while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
    		String[] tokens=nextLine.split("\t");
    		rtrn.put(tokens[0], tokens[3]);
    	}
    	
    	reader.close();
    	return rtrn;
	}


	public static void main(String[] args)throws IOException, InterruptedException{
		if(args.length>1){
			File file=new File(args[0]);
			String saveDir=args[1];
			new FTPRetrieveFiles(file, saveDir);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=GEO files \n args[1]=save directory";
	
}
