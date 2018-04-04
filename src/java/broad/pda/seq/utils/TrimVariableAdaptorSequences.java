package broad.pda.seq.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

import broad.pda.seq.fastq.FastqSequence;

public class TrimVariableAdaptorSequences {

	String linkerSeq="ATCTCGTATGCCGTCTTCTGCTTG";
	
	//take fastq and trim adpator sequences and rewrite fastq
	
	public TrimVariableAdaptorSequences(File fastq, String save)throws IOException{
		FileWriter writer=new FileWriter(save);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fastq)));
		String nextLine;
    	int i=0;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	if(nextLine.startsWith("@")){
        		String firstLine=nextLine;
        		String secondLine=reader.readLine();
        		String thirdLine=reader.readLine();
        		String fourthLine=reader.readLine();
        		
        		FastqSequence seq=new FastqSequence(firstLine, secondLine, thirdLine, fourthLine);
        		FastqSequence trimmed=seq.trimLinker(linkerSeq);
        		if(i%100000 ==0){System.err.println(i+" "+seq.getSequence()+"  "+trimmed.getSequence());}
        		writer.write(trimmed+"\n");
        		i++;
        		//int num=10;
        		//String lastNBps=getLastBps(sequence, num);
    			//String polyN=polyN("T", num);
    			//if(lastNBps.equalsIgnoreCase(polyN)){System.err.println(sequence);}
        	}
        }
		writer.close();
	}
	
	
	
	public static void main(String[] args)throws IOException{
		File file=new File(args[0]);
		String save=args[1];
		
		new TrimVariableAdaptorSequences(file, save);
	}
	
}
