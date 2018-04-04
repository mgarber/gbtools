package broad.pda.seq.joshuaAnalyses;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class SampleReads {

	public SampleReads(String samFile, String save, int numReadsToSample, int numReadsInFile, int perm) throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(samFile)));
    	String nextLine;
        FileWriter[] writers=new FileWriter[perm];
		
        for(int i=0; i<writers.length; i++){
        	writers[i]=new FileWriter(save+"."+i+".sam");
        }
		double fraction=(double)numReadsToSample/(double)numReadsInFile;
		
		
		while ((nextLine = reader.readLine()) != null) {
			double rand=Math.random();
			double[] rands=new double[perm];
			for(int i=0; i<rands.length; i++){
				rands[i]=Math.random();
				if(rands[i]<fraction){writers[i].write(nextLine+"\n");}
			}
		}
		
		for(int i=0; i<writers.length; i++){writers[i].close();}
	}

	public static void main(String[] args)throws IOException{
		if(args.length>3){
			String sam=args[0];
			String save=args[1];
			int numToSample=new Integer(args[2]);
			int numInFile=new Integer(args[3]);
			int numPerm=new Integer(args[4]);
			new SampleReads(sam, save, numToSample, numInFile, numPerm);
		}
		else{System.err.println(usage);}
	}
	
	static String usage="args[0]=sam \n args[1]=save \n args[2]=numToSample \n args[3]=numInFile \n args[4]=numPerms";
	
}
