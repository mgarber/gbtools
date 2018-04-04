package broad.pda.arrays.tilling;

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

import broad.core.math.Statistics;
import broad.pda.datastructures.Alignments;

public class MakeIGVFile {

	public MakeIGVFile(File[] files, String save, boolean medianNorm, boolean isLog)throws IOException{
		Map<Alignments, ArrayList<Double>> maps=parse(files);
		if(medianNorm){maps=medianNorm(maps, isLog);}
		write(save, maps, files);
	}
	
	
	private Map<Alignments, ArrayList<Double>> medianNorm(Map<Alignments, ArrayList<Double>> maps, boolean isLog) {
		Map<Alignments, ArrayList<Double>> rtrn=new TreeMap();
		
		//1) Compute median for each array
		Collection<ArrayList<Double>> values=maps.values();
		ArrayList<Double> median=getMedian(values);
		
		//2) divide through by the median (or subtract)
		for(Alignments align: maps.keySet()){
			ArrayList<Double> vals=maps.get(align);
			ArrayList<Double> norm=norm(vals, median, isLog);
			rtrn.put(align, norm);
		}
		return rtrn;
	}


	private ArrayList<Double> getMedian(Collection<ArrayList<Double>> values) {
		ArrayList<Double> rtrn=new ArrayList();
		int num=values.iterator().next().size();
		
		for(int i=0; i<num; i++){
			double median=median(values, i);
			rtrn.add(median);
		}
		
		return rtrn;
	}


	private double median(Collection<ArrayList<Double>> values, int i) {
		ArrayList<Double> temp=new ArrayList();
		for(ArrayList<Double> vals: values){
			temp.add(vals.get(i));
		}
		return Statistics.average(temp);
	}


	private ArrayList<Double> norm(ArrayList<Double> vals, ArrayList<Double> median, boolean isLog) {
		ArrayList<Double> rtrn=new ArrayList();

		int index=0;
		for(Double val: vals){
			double norm=val/median.get(index);
			if(isLog){norm=val-median.get(index);}
			else{norm=Math.log(norm)/Math.log(2);}
			rtrn.add(norm);
			index++;
		}
		
		return rtrn;
	}


	private Map<Alignments, ArrayList<Double>> parse(File[] files)throws IOException{
		Map<Alignments, ArrayList<Double>> rtrn=new TreeMap<Alignments, ArrayList<Double>>();
		
		for(int i=0; i<files.length; i++){
			System.err.println(files[i]);
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(files[i])));
			String nextLine;
			while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
				Alignments name=new Alignments(nextLine.split("\t")[0]);
				Double val=new Double(nextLine.split("\t")[1]);
				ArrayList<Double> list=new ArrayList<Double>();
				if(rtrn.containsKey(name)){list=rtrn.get(name);}
				list.add(val);
				rtrn.put(name, list);
			}
		}
		
		return rtrn;
	}
	
	private void write(String save, Map<Alignments, ArrayList<Double>> maps, File[] files)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		writer.write("Chromosome\tStart\tEnd\tName");
		for(int i=0; i<files.length; i++){
			writer.write("\t"+files[i].getName());
		}
		writer.write("\n");
		
		
		for(Alignments align: maps.keySet()){
			//Alignments align=new Alignments(name);
			ArrayList<Double> vals=maps.get(align);
			writer.write(align.getChr()+"\t"+align.getStart()+"\t"+align.getEnd()+"\t"+align.toUCSC());
			for(Double val: vals){writer.write("\t"+val);}
			writer.write("\n");
		}
		
		writer.close();
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>3){
		File[] files=new File(args[0]).listFiles();
		String save=args[1];
		boolean medianNorm=new Boolean(args[2]);
		boolean isLog=new Boolean(args[3]);
		new MakeIGVFile(files, save, medianNorm, isLog);
		}
		else{System.err.println(usage);}
	}

	static String usage=" args[0]=files \n args[1]=save file \n args[2]=median norm? \n args[3]=is log?";
}
