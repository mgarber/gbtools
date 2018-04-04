package broad.pda.rnai.designer;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

public class ConcatanateTopReports {

	
	public ConcatanateTopReports(File[] hairpinFiles, String save)throws Exception{
		Map<String, Collection<HairpinKmer>> hairpins=RNAiFileFormatUtils.parseHairpinDesignFiles(hairpinFiles);
		write(save, hairpins);
	}
	
	private void write(String save,	Map<String, Collection<HairpinKmer>> hairpins) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		
		
		for(String transcript: hairpins.keySet()){
			System.out.println(transcript+"\t"+hairpins.get(transcript).size());
			Collection<HairpinKmer> kmers=hairpins.get(transcript);
			for(HairpinKmer kmer: kmers){
				writer.write(transcript.split("_")[0]+"\t"+transcript+"\t"+kmer.getKmerStartPosition()+"\t"+kmer.getRS8Score()+"\t"+kmer.getKmerSequence()+"\n");
			}
		}
		
		writer.close();
	}

	public static void main(String[] args)throws Exception{
		if(args.length>1){
			File[] haiprinFiles=new File(args[0]).listFiles();
			String save=args[1];
			new ConcatanateTopReports(haiprinFiles, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=hairpin files \n args[1]=save";
	
}
