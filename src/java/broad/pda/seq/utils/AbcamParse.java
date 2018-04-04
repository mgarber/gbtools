package broad.pda.seq.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

public class AbcamParse {

	public AbcamParse(File file, String save) throws IOException{
		Map<String, Collection<Antibody>> geneAntibodies=parse(file);
		Map<String,Antibody> best=findBestAntibody(geneAntibodies);
		write(save, best);
	}
	
	
	private Map<String, Collection<Antibody>> parse(File file) throws IOException {
		Map<String, Collection<Antibody>> rtrn=new TreeMap();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		int i=0;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	if(i>0){
        	String[] tokens=nextLine.split("\t");
        	if(tokens.length>7){
	        	String gene=tokens[0].trim().toUpperCase();
	        	Antibody ab=new Antibody(tokens[2], tokens[0], tokens[5].split(","), tokens[6].split(","), tokens[7]);
	        	Collection<Antibody> set=new TreeSet();
	        	if(rtrn.containsKey(gene)){set=rtrn.get(gene);}
	        	set.add(ab);
	        	rtrn.put(gene, set);
        	}
        	else{System.err.println("Skipped "+nextLine);}
        	}
        	i++;
        }
		
        reader.close();
        
		return rtrn;
	}


	private Map<String, Antibody> findBestAntibody(Map<String, Collection<Antibody>> geneAntibodies) {
		Map<String, Antibody> rtrn=new TreeMap();
		
		for(String gene: geneAntibodies.keySet()){
			Collection<Antibody> antibodies=geneAntibodies.get(gene);
			Antibody antibody=findBestAntibody(antibodies);
			rtrn.put(gene, antibody);
		}
		
		return rtrn;
	}


	private Antibody findBestAntibody(Collection<Antibody> antibodies) {
		return antibodies.iterator().next();
	}


	private void write(String save, Map<String, Antibody> best) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("Gene Name"+"\t"+"AbID"+"\t"+"Internal Target ID"+"\t"+"Species"+"\t"+"Tested Applications"+"\t"+"Priority\n");
		
		for(String gene: best.keySet()){
			writer.write(best.get(gene).toString()+"\n");
		}
		
		writer.close();
	}


	public static class Antibody implements Comparable{
		Collection<String> species;
		Collection<String> applications;
		String id;
		String gene;
		int priority;
		String internalID;
		
		public Antibody(String id, String gene, String[] species, String[] applications, String internalID) {
			this.species=toCollection(species);
			this.applications=toCollection(applications);
			this.id=id;
			this.gene=gene;
			this.priority=setPriority();
			this.internalID=internalID;
		}
		
		private Collection<String> toCollection(String[] species2) {
			Collection<String> rtrn=new ArrayList();
			for(int i=0; i<species2.length; i++){
				rtrn.add(species2[i].replaceAll("\"", "").trim().toUpperCase());
			}
			return rtrn;
		}

		public String toString(){
			return gene+"\t"+id+"\t"+internalID+"\t"+species+"\t"+applications+"\t"+this.priority;
		}

		public int setPriority(){
			//Use priorities set by ido
			//1) Mouse and ChIP
			//2) ChIP
			//3) Mouse
			//4) IF
			//5) WB
			
			if(this.species.contains("MS") && this.applications.contains("CHIP")){return 1;}
			if(this.applications.contains("CHIP")){return 2;}
			if(this.species.contains("MS") && this.applications.contains("IP")){return 3;}
			if(this.species.contains("MS") && this.applications.contains("ICC/IF")){return 4;}
			if(this.species.contains("MS") && this.applications.contains("IHC-P")){return 5;}
			if(this.species.contains("MS")){return 6;}
			if(this.applications.contains("IP")){return 7;}
			if(this.applications.contains("ICC/IF")){return 8;}
			if(this.applications.contains("IHC-P")){return 9;}
			else{return 10;}
			
		}
		
		public int compareTo(Object arg0) {
			Antibody ab=(Antibody) arg0;
			
			return new Integer(priority).compareTo(ab.priority); 
		}
		
	}

	public static void main(String[] args)throws IOException{
		if(args.length>1){
			File file=new File(args[0]);
			String save=args[1];
			new AbcamParse(file, save);
		}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=file \n args[1]=save";

}


