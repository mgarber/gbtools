package broad.core.conservationAnalysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.apache.commons.math.stat.descriptive.rank.Percentile;

import broad.core.siphy.EvolutionaryModel.OmegaFit;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;

public class FilterElementsByPValue {
	
	
	public FilterElementsByPValue(File BEDScores, String save, double alpha)throws IOException{
		parseAndWrite(BEDScores, save,alpha);
	}
	
	private void write(String save, Map<Alignments, List<OmegaFit>> omega, double alpha)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(Alignments align: omega.keySet()){
			List<OmegaFit> list=omega.get(align);
			for(OmegaFit fit: list){
				if(fit.getPVal()<alpha){writer.write(align+"\t"+fit.getRegion()+"\t"+fit.getOmega()+"\t"+fit.getLogOddsScore()+"\t"+fit.getPVal()+"\n");}
			}
		}
		
		writer.close();
	}
	
	
	private void parseAndWrite(File file, String save, double alpha)throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
        FileWriter writer=new FileWriter(save);
		
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	String[] tokens=nextLine.split("\t");
        	
        	
        	double p=new Double(tokens[6]);
        	if(p<alpha){writer.write(nextLine+"\n");}
        	
        }
		
        reader.close();
		writer.close();
	}
	

	public FilterElementsByPValue(File BEDFile, File modelFile, String alnDir, String alnFileFormat, String save, int windowSize, double percentile)throws Exception{
		Set<Alignments> set=BEDFileParser.loadAlignmentData(BEDFile);
		Map<Alignments, List<OmegaFit>> omega=this.computeOmegaPerExon(set, modelFile, alnDir, alnFileFormat, windowSize);
		Map<Alignments, List<OmegaFit>> percentileMap=getTopPercentileElements(omega, percentile);
		write(save, percentileMap);
	}
	
	private Map<Alignments, List<OmegaFit>> getTopPercentileElements(Map<Alignments, List<OmegaFit>> omega, double percentile){
		Map rtrn=new TreeMap();
		
		for(Alignments align: omega.keySet()){
			List<OmegaFit> list=omega.get(align);
			double crit=percentileCutoff(list, percentile);
			List<OmegaFit> filtered=filter(list, crit);
			rtrn.put(align, filtered);
		}
		
		return rtrn;
	}
	
	private List<OmegaFit> filter(List<OmegaFit> list, double crit){
		List<OmegaFit> rtrn=new ArrayList();
		
		for(OmegaFit fit: list){
			double score=fit.getLogOddsScore();
			if(score>crit){rtrn.add(fit);}
		}
		
		return rtrn;
	}
	
	private double percentileCutoff(List<OmegaFit> list, double percentile){
		double[] array=makeDoubleArray(list);
		Percentile per=new Percentile(percentile);
		double crit=per.evaluate(array, percentile);
		return crit;
	}
	
	private double[] makeDoubleArray(List<OmegaFit> list){
		double[] rtrn=new double[list.size()];
		
		int i=0;
		for(OmegaFit fit: list){
			rtrn[i++]=fit.getLogOddsScore();
		}
		
		return rtrn;
	}
	
	private void write(String save, Map<Alignments, List<OmegaFit>> omega)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(Alignments align: omega.keySet()){
			List<OmegaFit> fitList=omega.get(align);
			for(OmegaFit fit: fitList){
				writer.write(align+"\t"+fit.getRegion()+"\t"+fit.getOmega()+"\t"+fit.getLogOddsScore()+"\t"+fit.getPVal()+"\n");
			}
		}
		
		writer.close();
	}
	
	private Map<Alignments, List<OmegaFit>> computeOmegaPerExon(Set<Alignments> set, File modelFile, String alnDir, String alnFileFormat, int windowSize)throws Exception{
		Map rtrn=new TreeMap();
		for(Alignments align: set){
			System.err.println(align);
			String alnFile=alnDir+"/"+align.getChr()+".maf";
			List<OmegaFit> fit=EstimateOmegaPerExon.slideWindowComputeOmega(align, modelFile, alnFile, alnFileFormat, windowSize, windowSize-1);
			rtrn.put(align, fit);
		}
		return rtrn;
	}
	
	public static void main(String[] args)throws Exception{
		if(args.length>2){
			File bedScores=new File(args[0]);
			String save=args[1];
			double percentile=new Double(args[2]);
			//double alpha=new Double(args[3]);
			new FilterElementsByPValue(bedScores, save, percentile);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=BED Scores \n args[1]=save file \n args[2]=percentile";
	
	
}