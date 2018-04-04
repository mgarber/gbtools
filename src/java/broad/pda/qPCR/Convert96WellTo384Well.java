package broad.pda.qPCR;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.pda.annotation.BEDFileParser;

public class Convert96WellTo384Well {

	public Convert96WellTo384Well(File file, String save) throws IOException{
		Map<String, String> plate96Format=parse96Format(file);
		//Map<String, String> plate384Format=convert96To384(plate96Format);
		write(save, plate96Format);
		
	}
	
	private void write(String save, Map<String, String> plate96Format) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("Location\tExperiment\n");
		
		for(String location: plate96Format.keySet()){
			Collection<String> plate384=convert96To384(location);
			for(String p: plate384){
				writer.write(p+"\t"+plate96Format.get(location)+"\n");
			}
			
		}
		
		writer.close();
		
	}

	private Map<String, String> convert96To384(Map<String, String> plate96Format) {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		for(String location: plate96Format.keySet()){
			String val=plate96Format.get(location);
			Collection<String> all=convert96To384(location);
			for(String n: all){rtrn.put(n, val);}
		}
		
		return rtrn;
	}

	private Collection<String> convert96To384(String location) {
		Collection<String> rtrn=new TreeSet<String>();
		
		String row=location.split("_")[0];
		String column=location.split("_")[1];
		
		int num=new Integer(column);
		int val1=num*2;
		int val2=val1-1;
		
		if(row.equalsIgnoreCase("A")){
			//then A and B
			rtrn.add("A"+val1);
			rtrn.add("A"+val2);
			rtrn.add("B"+val1);
			rtrn.add("B"+val2);
		}
		else if(row.equalsIgnoreCase("B")){
			rtrn.add("C"+val1);
			rtrn.add("C"+val2);
			rtrn.add("D"+val1);
			rtrn.add("D"+val2);
		}
		else if(row.equalsIgnoreCase("C")){
			//E, F
			rtrn.add("E"+val1);
			rtrn.add("E"+val2);
			rtrn.add("F"+val1);
			rtrn.add("F"+val2);
		}
		else if(row.equalsIgnoreCase("D")){
			//G, H
			rtrn.add("G"+val1);
			rtrn.add("G"+val2);
			rtrn.add("H"+val1);
			rtrn.add("H"+val2);
		}
		else if(row.equalsIgnoreCase("E")){
			//I,J
			rtrn.add("I"+val1);
			rtrn.add("I"+val2);
			rtrn.add("J"+val1);
			rtrn.add("J"+val2);
		}
		else if(row.equalsIgnoreCase("F")){
			//K, L
			rtrn.add("K"+val1);
			rtrn.add("K"+val2);
			rtrn.add("L"+val1);
			rtrn.add("L"+val2);
		}
		else if(row.equalsIgnoreCase("G")){
			//M, N
			rtrn.add("M"+val1);
			rtrn.add("M"+val2);
			rtrn.add("N"+val1);
			rtrn.add("N"+val2);
		}
		else if(row.equalsIgnoreCase("H")){
			//O,P
			rtrn.add("O"+val1);
			rtrn.add("O"+val2);
			rtrn.add("P"+val1);
			rtrn.add("P"+val2);
		}
	
		return rtrn;
	}

	private Map<String, String> parse96Format(File file) throws IOException{
		Map<String, String> plate96Format=new TreeMap<String, String>();
		Collection<String> lines=BEDFileParser.loadList(file.getAbsolutePath(), true);
		int rowNum=1;
		for(String line: lines){
			String[] tokens=line.split("\t");
			for(int columnNum=1; columnNum<tokens.length; columnNum++){
				String sample=tokens[columnNum];
				String position=convert(rowNum, columnNum);
				plate96Format.put(position, sample);
			}
			
			rowNum++;
		}
		return plate96Format;
	}

	private String convert(int rowNum, int columnNum) {
		String position="";
		if(rowNum==1){position+="A";}
		else if(rowNum==2){position+="B";}
		else if(rowNum==3){position+="C";}
		else if(rowNum==4){position+="D";}
		else if(rowNum==5){position+="E";}
		else if(rowNum==6){position+="F";}
		else if(rowNum==7){position+="G";}
		else if(rowNum==8){position+="H";}
		
		position+="_"+columnNum;
		return position;
	}
	
	public static void main(String[] args)throws IOException{
		File file=new File(args[0]);
		String save=args[1];
		new Convert96WellTo384Well(file, save);
	}
	
}
