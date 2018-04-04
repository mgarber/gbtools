package broad.pda.seq.pairedend;

import java.io.File;
import java.io.IOException;
import java.util.Collection;

public class PairedEndDistributions {

	public PairedEndDistributions(File leftFile, File rightFile, String chr)throws IOException{
		Collection<PairedEndAlignment> pairs=PairedEndMapping.getAllPairs(leftFile, rightFile, chr);//make more efficient by caching chunks
		
	}
	
	public static void main(String[] args){
		
	}
}
