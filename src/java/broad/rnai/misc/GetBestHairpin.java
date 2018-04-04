package broad.rnai.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

//takes a validation file and orders the hairpins by gene name
//ranks them by percent knock down
//picks the best hairpin that has a percent knock down greater than some cutoff

public class GetBestHairpin {

	double kdCutoff=.5;
	double enStdCutoff=1;
	
	private Map<String, Collection<String>> hairpinsPerGene;
	private Map<String, Double> validationPerHairpin;
	private Map<String, String> hairpinFullInfo;
	Map<String, Double> bestHairpins;
	
	public GetBestHairpin(File validationFile, String save) throws IOException{
		this.hairpinsPerGene=new TreeMap();
		this.validationPerHairpin=new TreeMap();
		this.bestHairpins=new TreeMap();
		this.hairpinFullInfo=new TreeMap();
		
		parse(validationFile);
		this.bestHairpins=getBestHairpins();
		write(save);
	}
	
	private Map<String, Double> getBestHairpins() {
		Map<String, Double> rtrn=new TreeMap();
		
		for(String gene: hairpinsPerGene.keySet()){
			Collection<String> hairpins=hairpinsPerGene.get(gene);
			String hairpin=getBest(hairpins);
			double percentKD=this.validationPerHairpin.get(hairpin);
			if(percentKD<this.kdCutoff){rtrn.put(hairpin, this.validationPerHairpin.get(hairpin));}
		}
		
		return rtrn;
	}

	private String getBest(Collection<String> hairpins) {
		double min=Double.MAX_VALUE;
		String best="None";
		for(String hairpin: hairpins){
			if(this.validationPerHairpin.get(hairpin)<min){min=this.validationPerHairpin.get(hairpin); best=hairpin;}
		}
		return best;
	}

	private void write(String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write(this.hairpinFullInfo.get("header")+"\n");
		for(String hairpin: this.bestHairpins.keySet()){
			writer.write(this.hairpinFullInfo.get(hairpin)+"\n");
		}
		
		writer.close();
	}

	private void parse(File file) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		int i=0;
		while ((nextLine = reader.readLine()) != null) {
			if(i>0){
			String[] tokens=nextLine.split("\t");
			String hairpinName=tokens[0];
			String geneName=hairpinName.split("-")[0];
			double knockDown=1;
			try{knockDown=new Double(tokens[14]);}catch(Exception ex){}
			
			double variation=10;
			try{variation=new Double(tokens[15]);}catch(Exception ex){}
			if(variation<this.enStdCutoff){
				System.err.println(hairpinName+" "+knockDown);
				this.validationPerHairpin.put(hairpinName, knockDown);
				this.hairpinFullInfo.put(hairpinName, nextLine);
				Collection<String> hairpins=new TreeSet();
				if(this.hairpinsPerGene.containsKey(geneName)){
					hairpins=this.hairpinsPerGene.get(geneName);
				}
				hairpins.add(hairpinName);
				this.hairpinsPerGene.put(geneName, hairpins);
				}
			}
			else{this.hairpinFullInfo.put("header", nextLine);}
			i++;
		}
		reader.close();
	}

	public static void main(String[] args)throws IOException{
		if(args.length>1){
			File file=new File(args[0]);
			String save=args[1];
			new GetBestHairpin(file, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=validation file \n args[1]=save";
}
