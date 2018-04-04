package broad.pda.arrays.tilling;

import java.io.File;
import java.io.IOException;
import java.util.Set;

import broad.pda.datastructures.ContinuousData;

public class RunAllChromosomes {

	public static void main(String[] args)throws IOException{
		if(args.length>0){
		Runtime run=java.lang.Runtime.getRuntime();
		File[] files=new File(args[0]).listFiles();
		String savedir=args[1];
		int numTails=new Integer(args[2]);
		int numPerms=new Integer(args[3]);
		double alpha=new Double(args[4]);
		
		String script=args[5];
		String queue=args[6];
		
		
		for(int i=0; i<files.length; i++){
		ContinuousData data=new ContinuousData(files[i]);
		Set<String> chromosomes=data.getChromosomes();
		for(String chr: chromosomes){
			String save=savedir+"/"+chr+"."+files[i].getName();
			String command="bsub -q "+queue+" -o  "+save+".bsub"+" java -jar -Xmx2000m ";
			command=command+script+" "+files[i].getAbsolutePath()+" "+save+" "+numTails+" "+numPerms+" "+alpha+" "+chr;
			System.err.println(command);
			run.exec(command);
		}
		}
		}
		else{
			System.out.println("Function: "+function+"\n"+"USAGE: "+usage);
		}
	}
	
	static String function="Runs window analysis on all chromosomes seperately on the LSF cluster";
	static String usage="args[0]=input directory \nargs[1]=directory to save to \nargs[2]=number of tails \nargs[3]=number of permutations \nargs[4]=alpha (.9,.95) \nargs[5]=script name \nargs[6]=queue";
	
}
