package broad.pda.ribosome;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import broad.pda.annotation.BEDFileParser;

public class MergeFeatureScores {

	public MergeFeatureScores(File[] files, String saveDir, String extension) throws IOException {
		Map<String, Collection<File>> filesByName=getFilesByName(files, extension);
		mergeAndWrite(filesByName, saveDir);
	}


	private void mergeAndWrite(Map<String, Collection<File>> filesByName, String saveDir) throws IOException {
		for(String baseName: filesByName.keySet()){
			Collection<File> files=filesByName.get(baseName);
			String save=saveDir+"/"+baseName+".scores";
			write(save, files);
		}
	}


	private void write(String save, Collection<File> files) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		boolean started=false;
		
		for(File file: files){
			Collection<String> lines=BEDFileParser.loadList(file.getAbsolutePath(), started);
			for(String line: lines){writer.write(line+"\n");}
			started=true;
		}
		
		writer.close();
	}


	private Map<String, Collection<File>> getFilesByName(File[] files, String extension) {
		Map<String, Collection<File>> rtrn=new TreeMap<String, Collection<File>>();
		
		for(int i=0; i<files.length; i++){
			File file=files[i];
			if(file.getName().endsWith(extension)){
				//grab base name
				String baseName=file.getName().split("\\.")[1];
				Collection<File> list=new ArrayList<File>();
				if(rtrn.containsKey(baseName)){
					list=rtrn.get(baseName);
				}
				list.add(file);
				rtrn.put(baseName, list);
			}
		}
		
		return rtrn;
	}


	public static void main(String[] args) throws IOException{
		if(args.length>2){
			File[] files=new File(args[0]).listFiles();
			String saveDir=args[1];
			String extension=args[2];
			new MergeFeatureScores(files, saveDir, extension);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=directory \n args[1]=save dir \n args[2]=extension to use";
	
}
