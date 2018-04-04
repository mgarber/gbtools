package broad.pda.ribosome.misc;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.util.PipelineUtils;
import broad.pda.annotation.BEDFileParser;

public class MergeSRAFiles {

	public MergeSRAFiles(File[] bamFiles, String infoFile, String saveDir) throws IOException, InterruptedException{
		Map<String, Collection<String>> info=parseInfoFile(infoFile);
		Map<String, File> files=assignFiles(bamFiles);
		
		mergeFilesBySRA(files, info, saveDir);	
	}
	
	private void mergeFilesBySRA(Map<String, File> files, Map<String, Collection<String>> info, String saveDir) throws IOException, InterruptedException {
		Runtime run=Runtime.getRuntime();
		String jobID=PipelineUtils.getJobID();
		
		for(String geo: info.keySet()){
			Collection<String> names=info.get(geo);
			String save=saveDir+"/"+geo+".bam";
			String output=saveDir+"/"+geo+".bsub";
			String cmd="java -jar /seq/mgarber/tools/picard-tools-1.66/MergeSamFiles.jar";
			for(String name: names){
				cmd+=" I="+files.get(name).getAbsolutePath();
			}
			cmd+=" O="+save;
			//run file
			System.err.println(cmd);
			//Process p=run.exec(cmd);
			//p.waitFor();
			
			String cmd2="java -jar /seq/mgarber/tools/picard-tools-1.66/BuildBamIndex.jar I="+save+" O="+save+".bai";
			//p=run.exec(cmd2);
			//p.waitFor();
			
			String[] commands={cmd, cmd2};
			PipelineUtils.bsubProcess(run, jobID, commands, output, "hour");
		}
		
	}

	private Map<String, File> assignFiles(File[] bamFiles) {
		Map<String, File> rtrn=new TreeMap<String, File>();
		
		for(int i=0; i<bamFiles.length; i++){
			File file=bamFiles[i];
			if(file.getName().endsWith(".bam")){
				String name=file.getName().split("\\.")[0];
				//System.err.println(name+" "+file);
				rtrn.put(name, file);
			}
		}
		
		return rtrn;
	}

	private Map<String, Collection<String>> parseInfoFile(String infoFile) throws IOException {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		Collection<String> lines=BEDFileParser.loadList(infoFile, true);
		
		for(String line: lines){
			String[] tokens=line.split("\t");
			String GEOID=tokens[0];
			String SRA=tokens[1];
			String name=tokens[2];
			Collection<String> list=new TreeSet<String>();
			if(rtrn.containsKey(GEOID)){list=rtrn.get(GEOID);}
			list.add(SRA);
			rtrn.put(GEOID, list);
		}
		
		return rtrn;
	}

	public static void main(String[] args)throws IOException, InterruptedException{
		if(args.length>2){
			File[] files=new File(args[0]).listFiles();
			String info=args[1];
			String saveDir=args[2];
			new MergeSRAFiles(files, info, saveDir);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam files \n args[1]=SRA info file \n args[2]=save dir";
	
}
