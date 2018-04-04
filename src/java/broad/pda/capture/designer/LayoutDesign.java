package broad.pda.capture.designer;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;
import broad.pda.seq.alignment.Pair;

public class LayoutDesign {

	public LayoutDesign(String probeFile, String genesPerArrayFile, String primerFile, String save, boolean writeAntisense) throws IOException{
		Map<String, Collection<String>> probesByGene=parseReport(probeFile);
		Map<Integer, List<String>> genesPerArray=parseGenesPerArray(genesPerArrayFile);
		Map<String, Boolean> useAntisense=parseAntisenseInfo(genesPerArrayFile, writeAntisense);
		Map<Pair<String>, List<String>> primers=parsePrimers(primerFile);
		ArrayList<Pair<String>> mainPrimers=new ArrayList<Pair<String>>();
		
		for(Pair<String> pairs: primers.keySet()){mainPrimers.add(pairs);}
		
		FileWriter writer=new FileWriter(save);
		writer.write("Gene\tArray Number\tTiling Path\tCoverage Path\tGene Left Primer\tGene Right Primer\tCoverage Path Left Primer\tTiling Path Left Primer\tProbe Sequence\tFull Probe Sequence\n");
		
		FileWriter faWriter=new FileWriter(save+".fa");
		
		//assign primers
		for(Integer arrayNum: genesPerArray.keySet()){
			FileWriter faWriter2=new FileWriter(save+"."+arrayNum+".fa");
			int counter=0;
			List<String> list=genesPerArray.get(arrayNum);
			for(int index=0; index<list.size(); index++){
				String gene=list.get(index);
				Collection<String> probes=probesByGene.get(gene);
				for(String probe: probes){
					String[] tokens=probe.split("\t");
					String probeSeq=tokens[8];
					int tilingPath=new Integer(tokens[2]);
					int coveragePath=assignCoverageNumber(tilingPath);
					Pair<String> mainPrimer=mainPrimers.get(index);
					List<String> subprimers=primers.get(mainPrimer);
					String subprimer=subprimers.get(coveragePath);
					
					int tilingPrimerNum=getTilingPrimerNum(tilingPath);
					String subsub=Sequence.get3Prime(subprimers.get(tilingPrimerNum), 5);
					String subsubprimer=Sequence.get3Prime(mainPrimer.getValue1(), 10)+Sequence.get3Prime(subprimer, 5)+subsub;
					
					String fullProbeSeq=mainPrimer.getValue1()+Sequence.get3Prime(subprimer, 5)+subsub+probeSeq+Sequence.reverseSequence(mainPrimer.getValue2());
					
					
					writer.write(gene+"\t"+arrayNum+"\t"+tilingPath+"\t"+coveragePath+"\t"+mainPrimer.getValue1()+"\t"+mainPrimer.getValue2()+"\t"+subprimer+"\t"+subsubprimer+"\t"+probeSeq+"\t"+fullProbeSeq+"\n");
					faWriter.write(">"+gene+"_"+tilingPath+"_"+(counter++)+"_Sense"+"\n"+fullProbeSeq+"\n");
					boolean writeAntisenseLocal=useAntisense.get(gene);
					if(writeAntisenseLocal){
						faWriter.write(">"+gene+"_"+tilingPath+"_"+(counter++)+"_Antisense"+"\n"+Sequence.reverseSequence(fullProbeSeq)+"\n");
					}
					faWriter2.write(">"+gene+"_"+tilingPath+"_"+(counter++)+"_Sense"+"\n"+fullProbeSeq+"\n");
					if(writeAntisenseLocal){
						faWriter2.write(">"+gene+"_"+tilingPath+"_"+(counter++)+"_Antisense"+"\n"+Sequence.reverseSequence(fullProbeSeq)+"\n");
					}
				}
			}
			faWriter2.close();
		}
		
		writer.close();
		faWriter.close();
		
	}
	
	private Map<String, Boolean> parseAntisenseInfo(String genesPerArrayFile, boolean writeAntisense) throws IOException {
		Map<String, Boolean> rtrn=new TreeMap<String, Boolean>();
		
		Collection<String> list=BEDFileParser.loadList(genesPerArrayFile, true);
		
		for(String line: list){
			String[] tokens=line.split("\t");
			String gene=tokens[0];
			try{
			int num=new Integer(tokens[3]);
			boolean use=true;
			if(num==0){use=false;}
			rtrn.put(gene, use);
			}catch (Exception ex){ex.printStackTrace(); rtrn.put(gene, writeAntisense);}
		}
		
		return rtrn;
	}

	private int getTilingPrimerNum(int tilingPath) {
		if(tilingPath==0){return 5;}
		if(tilingPath==60){return 6;}
		if(tilingPath==30){return 7;}
		if(tilingPath==90){return 8;}
		if(tilingPath==15){return 9;}
		if(tilingPath==45){return 10;}
		if(tilingPath==75){return 11;}
		if(tilingPath==105){return 12;}
		
		return 13;
	}

	private Map<Pair<String>, List<String>> parsePrimers(String primerFile) throws IOException {
		Map<Pair<String>, List<String>> rtrn=new HashMap<Pair<String>, List<String>>();
		Collection<String> list=BEDFileParser.loadList(primerFile);
		
		for(String line: list){
			String[] tokens=line.split("\t");
			Collection<String> subprimer=new HashSet<String>();
			for(int i=2; i<tokens.length; i++){subprimer.add(tokens[i]);}
			List<String> set=new ArrayList<String>();
			set.addAll(subprimer);
			Pair<String> pair=new Pair<String>(tokens[0], tokens[1]);
			if(set.size()>14){
				rtrn.put(pair, set);
			}
		}
		
		
		return rtrn;
	}

	private Map<String, Collection<String>> parseReport(String probeFile) throws IOException {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		Collection<String> list=BEDFileParser.loadList(probeFile, true);
		
		for(String line: list){
			String[] tokens=line.split("\t");
			String gene=tokens[0];
			Collection<String> set=new ArrayList<String>();
			if(rtrn.containsKey(gene)){set=rtrn.get(gene);}
			set.add(line);
			rtrn.put(gene, set);
		}
		
		return rtrn;
	}

	private int assignCoverageNumber(int tilingPath) {
		if(tilingPath==0){return 0;}
		if(tilingPath==60){return 1;}
		if(tilingPath==30 || tilingPath==90){return 2;}
		if(tilingPath==15 || tilingPath==45 || tilingPath==75 || tilingPath==105){return 3;}
		
		System.err.println("Bad tiling path");
		return 4;
	}

	private Map<Integer, List<String>> parseGenesPerArray(String genesPerArrayFile) throws IOException {
		Map<Integer, List<String>> rtrn=new TreeMap<Integer, List<String>>();
		Collection<String> list=BEDFileParser.loadList(genesPerArrayFile, true);
		
		for(String line: list){
			String[] tokens=line.split("\t");
			String name=tokens[0];
			int num=new Integer(tokens[2]);
			List<String> set=new ArrayList<String>();
			if(rtrn.containsKey(num)){set=rtrn.get(num);}
			set.add(name);
			rtrn.put(num, set);
		}
		
		return rtrn;
	}

	public static void main(String[] args) throws IOException{
		if(args.length>4){
		String probeFile=args[0];
		String genesPerArray=args[1];
		String primerFile=args[2];
		String save=args[3];
		boolean writeAntisense=new Boolean(args[4]);
		new LayoutDesign(probeFile, genesPerArray, primerFile, save, writeAntisense);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=probe file \n args[1]=genes Per array \n args[2]=primer file \n args[3]=save \n args[4]=write antisense";
	
}
