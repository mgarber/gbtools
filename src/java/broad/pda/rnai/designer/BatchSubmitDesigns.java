package broad.pda.rnai.designer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class BatchSubmitDesigns {

	public static void main(String[] args)throws IOException{
		if(args.length>6){
		File file=new File(args[0]);
		File geneFastaFile=new File(args[1]);
		File geneCoordinateFile=new File(args[2]);
		String saveDir=args[3];
		int numSplit=new Integer(args[4]);
		
		File[] files=split(file, saveDir, numSplit);
		String script=args[5];
		String queue=args[6];
		String hairpinFile="";
		if(args.length>7){hairpinFile=(args[7]);}
		
		Runtime run=java.lang.Runtime.getRuntime();
		for(int i=0; i<files.length; i++){
			String junkFile=saveDir+"/"+files[i].getName()+".bsub";
			String command="bsub -q "+queue+" -o "+junkFile+" "+"java -jar -Xmx2000m "+script+" "+files[i].getAbsolutePath()+" "+geneFastaFile.getAbsolutePath()+" "+geneCoordinateFile.getAbsolutePath()+" "+saveDir+" "+hairpinFile;
			run.exec(command);
			System.err.println(command);
			System.gc();
		}
		}else{System.err.println(usage);}
	}
	
	private static File[] split(File file, String saveDir, int numSplit) throws IOException {
		new File(saveDir+"/temp/").mkdir();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		int i=0;
		String header="";
		FileWriter writer=new FileWriter(saveDir+"/temp/"+(i/numSplit)+".report");
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			if(i>0){
				if(i%numSplit ==0){writer.close(); writer=new FileWriter(saveDir+"/temp/"+(i/numSplit)+".report"); writer.write(header+"\n");}
				writer.write(nextLine+"\n");
			}
			else{header=nextLine; writer.write(header+"\n");}
			i++;
		}
		writer.close();
		return new File(saveDir+"/temp/").listFiles();
	}

	static String usage=" args[0]=File (RNAi Report) \n args[1]=gene fasta file \n args[2]=gene coordinate fie \n args[3]=save directory \n args[4]=num split \n args[5]=script \n args[6]=queue \n args[7]=already designed hairpins (optional)";
	
}
