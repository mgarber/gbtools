package broad.pda.ribosome.misc;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import broad.pda.annotation.BEDFileParser;

public class MergeWindowsIntoTable {

	//Take a directory and merge all chromosomes by shared name
	public MergeWindowsIntoTable(File[] files, String saveDir) throws IOException{
		//Step 1: Sort files
		Map<String, Collection<File>> filesByName=sortFilesByName(files);
	
		//Step 2: Merge and write
		merge(filesByName, saveDir);
	}
	
	
	private void merge(Map<String, Collection<File>> filesByName, String saveDir) throws IOException {
		for(String name: filesByName.keySet()){
			FileWriter writer=new FileWriter(saveDir+"/"+name);
			System.err.println("writing "+name);
			Collection<File> files=filesByName.get(name);
			for(File file: files){
				write(writer, file);
			}
			writer.close();
		}
		
	}


	private void write(FileWriter writer, File file) throws IOException {
		Collection<String> lines=BEDFileParser.loadList(file.getAbsolutePath());
		for(String line: lines){
			writer.write(line+"\n");
		}
	}


	private Map<String, Collection<File>> sortFilesByName(File[] files) {
		Map<String, Collection<File>> rtrn=new TreeMap<String, Collection<File>>();
		
		for(int i=0; i<files.length; i++){
			String name=splitChr(files[i].getName());
			System.err.println(files[i].getName()+" "+name);
			Collection<File> list=new ArrayList<File>();
			if(rtrn.containsKey(name)){
				list=rtrn.get(name);
			}
			list.add(files[i]);
			rtrn.put(name, list);
		}
		
		return rtrn;
	}


	private String splitChr(String name) {
		String chr=name.split("\\.")[0];
		return name.replaceAll(chr+".", "");
	}


	public static void main(String[] args) throws IOException{
		if(args.length>1){
			File[] files=new File(args[0]).listFiles();
			String saveDir=args[1];
			new MergeWindowsIntoTable(files, saveDir);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=directory \n args[1]=saveDir";
	
}
