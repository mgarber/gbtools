 package broad.pda.seq.alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

public class BroadTopHatAlignmentPipeline {

	private static Collection<File> findReads(String barcode, String laneNum, Runtime run) throws IOException, InterruptedException{
		String command="/seq/software/bin/goSolexaRun.pl "+barcode+" gerald";
		Process p=run.exec(command);
		int completed=p.waitFor();
		int exitVal=p.exitValue();
		String path=getInputStream(p.getInputStream());
		writeError(p.getErrorStream());
					
		Collection<File> files=getValidFile(path, laneNum);
		
		return files;
	}
	
	private static Collection<File> getValidFile(String path, String laneNum) {
		Collection<File> rtrn=new ArrayList<File>();
		
		File[] files=new File(path).listFiles();
		
		for(int i=0; i<files.length; i++){
			if(files[i].getName().startsWith("s_"+laneNum) && files[i].getName().endsWith("sequence.txt")){rtrn.add(files[i]); System.err.println(files[i]);}
		}
		
		return rtrn;
	}

	private static Map<String, Collection<File>> findReads(File file, Runtime run) throws IOException, InterruptedException{
		Map<String, Collection<File>> rtrn=new TreeMap<String, Collection<File>>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	String[] tokens=nextLine.split("\t");
        	//Name, Flowcell, lane number
        	String name=tokens[0]+"_"+tokens[1]+"_s"+tokens[2];
        	String barcode=tokens[1];
        	Collection<File> path=findReads(barcode, tokens[2], run);
        	rtrn.put(name, path);
        }
		
		return rtrn;
	}
	
	private static int runTopHat(String save, String referenceLocation, String fastq, Runtime run) throws IOException, InterruptedException{
		String command="/seq/mguttman/scripts/TopHat/bin/tophat -o "+save +" "+referenceLocation+" "+fastq;
		Process p=run.exec(command);
		int completed=p.waitFor();
		int exitVal=p.exitValue();
		writeError(p.getErrorStream());
				
		return exitVal;
	}
	
	
	
	private static void writeError(InputStream errorStream) throws IOException {
		BufferedReader reader=	new BufferedReader(new InputStreamReader(errorStream));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			System.err.println(nextLine);
		}
		System.err.println();
	}
	
	private static String getInputStream(InputStream stream) throws IOException{
		BufferedReader reader=	new BufferedReader(new InputStreamReader(stream));
		String rtrn="";
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			System.err.println(nextLine);
			rtrn+=nextLine;
		}
		
		return rtrn;
	}
	
	private static int sortAlignmentFile(String alignmentFile, String save, Runtime run) throws IOException, InterruptedException{
		String command="/xchip/igv/tools/igvtools sort "+alignmentFile+" "+save;
		Process p=run.exec(command);
		int completed=p.waitFor();
		int exitVal=p.exitValue();
		writeError(p.getErrorStream());
		return exitVal;
	}
	
	private static int indexAlignmentFile(String alignmentFile, Runtime run) throws IOException, InterruptedException{
		String command="/xchip/igv/tools/igvtools index "+alignmentFile;
		Process p=run.exec(command);
		int completed=p.waitFor();
		int exitVal=p.exitValue();
		writeError(p.getErrorStream());
		return exitVal;
	}
	
	public static void main(String[] args)throws IOException, InterruptedException{
		if(args.length>1){
			File barcodeNames=new File(args[0]);
			String saveDir=args[1];
			
			Runtime run=java.lang.Runtime.getRuntime();
			Map<String, Collection<File>> fileMap=findReads(barcodeNames, run);
			for(String name: fileMap.keySet()){
				Collection<File> files=fileMap.get(name);
				for(File file: files){
				//Make a symbolic link
				String command="ln -s "+file.getAbsolutePath()+" "+saveDir;
				run.exec(command);
			}
			}
			
			
			
			//String alignmentFile=alignmentDir+"/accepted_hits.sam";
			
			//Map<String, Collection<File>> files=findReads(barcodeNames, run);
			//runTopHat(alignmentDir, bowtie, files, finalSave, run);
		}
		else{System.err.println(usage);}
	}
	
	private static void runTopHat(String alignmentDir, String bowtie, Map<String, Collection<File>> files, String finalSave, Runtime run) throws IOException, InterruptedException {
		for(String name: files.keySet()){
			Collection<File> names=files.get(name);
			for(File file: names){
				System.err.println(file);
				runTopHat(alignmentDir+"/"+file.getName(), bowtie, file.getAbsolutePath(), run);
				String alignmentFile=alignmentDir+"/"+file.getName()+"/accepted_hits.sam";
				String sortedName=finalSave+"/"+file.getName()+"."+name+".sorted.sam";
				sortAlignmentFile(alignmentFile, sortedName, run);
				indexAlignmentFile(sortedName, run);
			}
		}
	}

	static String usage=" args[0]=barcodes \n args[1]=save dir";
	
	
}
