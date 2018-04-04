package broad.pda.gene;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.MathException;

import analysisTools.functionalAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.math.EmpiricalDistribution;
import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.core.util.ParseGCTFile;
import broad.pda.annotation.BEDFileParser;

public class NeighborAnalysis {
	
	private BEDFileParser geneSet;
	private BEDFileParser refSet;
	private BEDFileParser unionSet;//the union of refSet and GeneSet across the genome
	private HashMap<RefSeqGene,RefSeqGeneWithIsoforms[] > neighborMap;
	private HashMap<RefSeqGene,double[] > neighborExpCorr;
	private HashMap<RefSeqGene,double[] > neighborExpCorrPval;
	private HashMap<RefSeqGene,double[] > neighborDistance;
	private HashMap<RefSeqGene,double[] > neighborJSDistance;
	private HashMap<RefSeqGene,double[] > isDivergentToNeighbor; //is the current gene transcribed in a divergent transcription pattern from its neighbor?
	private HashMap<RefSeqGene,double[] > isClosestToNeighbor; //is the current gene the closest to this neighbor? (several "lincs" can have the same neighbor but only one is the closest)
	private HashMap<RefSeqGene,double[] > isConvergentToNeighbor; //is the current gene transcribed in a convergent transcription pattern from its neighbor?
	private HashMap<RefSeqGene,double[] > isSameDirectToNeighbor; //is the current gene transcribed in a sameDirect transcription pattern from its neighbor?
	private ArrayList<Double> randomPvals;
	
	static String usage="Usage: NeighborAnalysis -task <task name> "+
	"\n\t GetNeighbors: get the neighbors from the ref set of the transcripts in set1 \n\t\t-set1 <set1 BED> \n\t\t-ref <set2 BED> \n\t\t-outPrefix <Output file prefix> \n\t\t -distThreshold <distance between transcripts to be called divergent>" +
	"\n\t NeighborsExpCorr: for each gene is set1, find its closest neighbor from the ref set and report thier correlation in expression,distance and if they are the closest to the neighbor. \n\t\t-set1 <set1 BED> \n\t\t-ref <set2 BED> \n\t\t-outPrefix <Output file prefix> \n\t\t-refSubSet <a sublist of the ref to be accounted as neighbors>  \n\t\t-gct1 <exp values of set1> \n\t\t-gct2 <exp values of ref>  \n\t\t-gct1ProbBed <mapping between set1 and gct1 probes> \n\t\t-toLogTransform <optional: log transform the data> \n\t\t-addfactor<optional: if log transform the data, add this value> \n\t\t -distThreshold <distance between transcripts to be called divergent> \n\t\t chrSizeF <optional, chrSize file for corr pval calculations>\n\t\t -numPerm <number of permutation for pval calculations>\n\t\t -addToUnion<A bed file to be added to the reference gene set, so all known annotation will be considered when asking  isClosestNeighbor>"+
	"\n\t DivergentExpCorr: for each gene is set1, find a divergent transcription neighbor from the ref set and report thier correlation in expression,distance and if they are the closest to the neighbor \n\t\t-set1 <set1 BED> \n\t\t-ref <set2 BED> \n\t\t-outPrefix <Output file prefix> \n\t\t-gct1 <exp values of set1> \n\t\t-gct2 <exp values of ref>  \n\t\t-gct1ProbBed <mapping between set1 and gct1 probes> \n\t\t-toLogTransform <optional: log transform the data> \n\t\t-addfactor<optional: if log transform the data, add this value> \n\t\t chrSizeF <optional, chrSize file for corr pval calculations>\n\t\t -distThreshold <distance between transcripts to be called divergent>"+
	"\n\t RandNeighborsExpCorr: calc correlation of random gene pairs. \n\t\t-set1 <set1 BED>  \n\t\t-outPrefix <Output file prefix> -gct1 <exp values of set1> \n\t\t -gct1ProbBed <mapping between set1 and gct1 probes> \n\t\t-toLogTransform <optional: log transform the data> \n\t\t-addfactor<optional: if log transform the data, add this value> \n\t\t  -numPerm <number of permutation for pval calculations>\n"+
	"\n\t GetDivergent: extract trascritps from set1 that are divergent to transcripts in the reference set   \n\t\t-set1 <set1 BED> \n\t\t-ref <set2 BED> \n\t\t-outPrefix <Output file prefix> \n\t\t -distThreshold <distance between transcripts to be called divergent>" +
	"\n\t GoAnalysis: given a set of genes, determine enrichment of their neighbors with specific annotations using random permutations of the neighbor set. \n\t\t -set1<set1 bed> \n\t\t -ref <ref set bed> -\n\t\t -goFile <Go annotFile> \n\t\t -goFormat <format of the GO file> \n\t\t goSizeThreshold <def= 5> \n\t\t -geneids <Optional, select only this subset of gene names> \n\t\t -batch <optional, geneIdsF is a list of geneIds files>\n\t\t-out <outPrefix>"+
	"\n\t GoTermStats: given a set of genes, report stats on the intersection of their neighbors with specific GO term . \n\t\t -set1<set1 bed> \n\t\t -ref <ref set bed> -\n\t\t -goFile <Go annotFile> \n\t\t -goFormat <format of the GO file> \n\t\t goSizeThreshold <def= 5> \n\t\t -GoTerm <GO:XXXX> -out <outPrefix>"+
	"\n\t pairedNeighbors : find sets of genes that have a linc neighbors on one side and a coding on the other side, mark Div/conv/uni . \n\t\t -set1 \n\t\t -set2 \n\t\t -outPrefix \n\t\t "+
	"\n";

	public static void main(String [] args) throws Exception  {
		ArgumentMap argmap = CLUtil.getParameters(args, usage, "orient");
		
		if("GetNeighbors".equalsIgnoreCase(argmap.getTask())) {
			String setname = argmap.getMandatory("set1");
			String refname = argmap.getMandatory("ref");
			String out = argmap.getMandatory("outPrefix");
			String subsetf=argmap.containsKey("refSubSet")? argmap.get("refSubSet"): null;
			double distThreshold=argmap.getDouble("distThreshold");
			
			NeighborAnalysis na= new NeighborAnalysis(setname,refname);
			//BEDFileParser set=new BEDFileParser(setname);
			//BEDFileParser ref=new BEDFileParser(refname);
			BEDFileParser neighborsBed = new BEDFileParser();
			BufferedWriter bw=new BufferedWriter(new FileWriter(out+"_neighborsMap"));
			
			if (subsetf==null)
				na.getNeighbors(neighborsBed,bw,false,distThreshold);
			else{
				Set<String> subset=ReadStringList(subsetf);
				na.getNeighbors(neighborsBed,bw,subset,false,distThreshold);//10 kb for divergent 
			}
			neighborsBed.writeFullBed(out+"_allNeighbors.bed");
		}
		
		else if("NeighborsExpCorr".equalsIgnoreCase(argmap.getTask())) {
			String setname = argmap.getMandatory("set1");
			String refname = argmap.getMandatory("ref");
			String out = argmap.getMandatory("outPrefix");
			String subsetf=argmap.containsKey("refSubSet")? argmap.get("refSubSet"): null;
			String set1gct = argmap.getMandatory("gct1");
			String set2gct = argmap.getMandatory("gct2");
			String probRegionMap = argmap.getMandatory("gct1ProbBed");
			boolean toLog= argmap.containsKey("toLogTransform")? true: false;
			double addFactor=argmap.containsKey("addFactor")? argmap.getDouble("addFactor"): 0;
			double distThreshold=argmap.containsKey("distThreshold")? argmap.getDouble("distThreshold"): 10000;
			boolean pearsonOnDensity = argmap.containsKey("pearsonOnDensity")? true:false;
			//Another bed file to be added to the union set (so we can say correctly if closest?)
			String addToUnion=argmap.containsKey("addToUnion")? argmap.get("addToUnion"): "";
			String chrSizeF= argmap.containsKey("chrSizeF")? argmap.get("chrSizeF"):"";
			int numPerm= argmap.containsKey("numPerm")? argmap.getInteger("numPerm"):100;
			NeighborAnalysis na= new NeighborAnalysis(setname,refname);
			if (!addToUnion.equalsIgnoreCase("")){
				BEDFileParser bed = new BEDFileParser(addToUnion);
				na.unionSet.addRefSeqSet(bed.GetGenes());
			}
			BEDFileParser neighborsBed = new BEDFileParser();
			if (subsetf==null)
				na.getNeighbors(neighborsBed,null,false,distThreshold);
			else{
				Set<String> subset=ReadStringList(subsetf);
				na.getNeighbors(neighborsBed,null,subset,false,distThreshold);
			}
			BEDFileParser probRegionBed =new BEDFileParser(probRegionMap);
			na.calcGeneNeighborExpCorr(set1gct,set2gct,setname,refname,probRegionBed,toLog,addFactor,chrSizeF,numPerm,pearsonOnDensity);
			na.printGeneNeighborDataTab(out+"_neighborsMap");
			na.printRandomCorr(out+"_randCorr");
		}
		else if("DivergentExpCorr".equalsIgnoreCase(argmap.getTask())) {
			String setname = argmap.getMandatory("set1");
			String refname = argmap.getMandatory("ref");
			String out = argmap.getMandatory("outPrefix");
			String set1gct = argmap.getMandatory("gct1");
			String set2gct = argmap.getMandatory("gct2");
			String probRegionMap = argmap.getMandatory("gct1ProbBed");
			boolean toLog= argmap.containsKey("toLogTransform")? true: false;
			double addFactor=argmap.containsKey("addFactor")? argmap.getDouble("addFactor"): 0;
			double distThreshold=argmap.getDouble("distThreshold");
			String chrSizeF= argmap.containsKey("chrSizeF")? argmap.get("chrSizeF"):"";
			int numPerm= argmap.containsKey("numPerm")? argmap.getInteger("numPerm"):100;
			boolean pearsonOnDensity = argmap.containsKey("pearsonOnDensity")? true:false;
			
			
			NeighborAnalysis na= new NeighborAnalysis(setname,refname);
			BEDFileParser neighborsBed = new BEDFileParser();
			na.getNeighbors(neighborsBed,null,false,distThreshold);
			
			BEDFileParser probRegionBed =new BEDFileParser(probRegionMap);
			na.calcGeneNeighborExpCorr(set1gct,set2gct,setname,refname,probRegionBed,toLog,addFactor,chrSizeF,numPerm,pearsonOnDensity);
			na.printDiveregentTranscriptionDataTab(out+"_DivergentTranscriptMap",distThreshold);
		}
		else if ("RandNeighborsExpCorr".equalsIgnoreCase(argmap.getTask())){
			
			String setname = argmap.getMandatory("set1");
			String out = argmap.getMandatory("outPrefix");
			String set1gct = argmap.getMandatory("gct1");
			boolean toLog= argmap.containsKey("toLogTransform")? true: false;
			double addFactor=argmap.containsKey("addFactor")? argmap.getDouble("addFactor"): 0;
			boolean pearsonOnDensity = argmap.containsKey("pearsonOnDensity")? true:false;
			int numPerm= argmap.containsKey("numPerm")? argmap.getInteger("numPerm"):10000;
			calcRandNeighborsExpCorr (set1gct,out,toLog,addFactor,pearsonOnDensity,numPerm);
		}
		else if("GetDivergent".equalsIgnoreCase(argmap.getTask())) {
			String setname = argmap.getMandatory("set1");
			String refname = argmap.getMandatory("ref");
			String out = argmap.getMandatory("outPrefix");
			double distThreshold=argmap.getDouble("distThreshold");
			
			NeighborAnalysis na= new NeighborAnalysis(setname,refname);
			BEDFileParser neighborsBed = new BEDFileParser();
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			na.getNeighbors(neighborsBed,bw,false,distThreshold);
			na.printDivergent(bw);
			bw.close();
			
		}
		else if ("GoAnalysis".equalsIgnoreCase(argmap.getTask())){
			String setname = argmap.getMandatory("set1");
			String refname = argmap.getMandatory("ref");
			String gofile =argmap.getMandatory("goFile");
			String goformat =argmap.getMandatory("goFormat");
			int numPerm= argmap.containsKey("numPerm")? argmap.getInteger("numPerm"):100;
			String geneidsF= argmap.containsKey("geneids")? argmap.get("geneids"):"";
			int goSizeThreshold= argmap.containsKey("goSizeThreshold")? argmap.getInteger("goSizeThreshold"):5;
			String outprefix= argmap.getOutput();
			Boolean  batchMode = argmap.containsKey("batch");
			if (batchMode)
				BatchGOAnalysis (setname,refname,gofile,goformat,outprefix,numPerm,goSizeThreshold,geneidsF);
			else
				GOAnalysis (setname,refname,gofile,goformat,outprefix,numPerm,goSizeThreshold,geneidsF);
		
		}
		else if ("GoTermStats".equalsIgnoreCase(argmap.getTask())){
			String setname = argmap.getMandatory("set1");
			String refname = argmap.getMandatory("ref");
			String gofile =argmap.getMandatory("goFile");
			String goformat =argmap.getMandatory("goFormat");
			int numPerm= argmap.containsKey("numPerm")? argmap.getInteger("numPerm"):100;
			int goSizeThreshold= argmap.containsKey("goSizeThreshold")? argmap.getInteger("goSizeThreshold"):5;
			String GoTerm =  argmap.getMandatory("GoTerm");
			String outprefix= argmap.getOutput();
			GOTermStats (setname,refname,gofile,goformat,outprefix,numPerm,goSizeThreshold,GoTerm );
		
		}
		else if ("testNull".equalsIgnoreCase(argmap.getTask())){
			String setname = argmap.getMandatory("set1");
			String refname = argmap.getMandatory("ref");
			String chrSizeF = argmap.getMandatory("chrSizeF");
			NeighborsNullModel na =new NeighborsNullModel (setname,refname,10,10);
			
			na.makeBinnedIntervalModel(true);
			//na.makeNaiveRandModel(chrSizeF);
			
			na.printRandNeighborsInfo();
			
		} else if ("NeighborsComparison".equalsIgnoreCase(argmap.getTask())) {
			String targets  = argmap.getMandatory("targets"); //List of annotation names here we do not really need any bed file
			String set1File = argmap.getMandatory("set1"); // Annotation list in BED format
			String set2File = argmap.getMandatory("set2"); // The other annotation list in BED format
			
			BEDFileParser set1 = new BEDFileParser(set1File);
			BEDFileParser set2 = new BEDFileParser(set2File);
			
			BufferedReader br = new BufferedReader(new FileReader(targets));
			ArrayList<String> targetList = new ArrayList<String>();
			String line = null;
			while  ( (line = br.readLine())!= null) {
				line.trim();
				targetList.add(line);
			}
			
			for (String t : targetList) {
				List<RefSeqGeneWithIsoforms> set1Annotations = set1.getAnnotations(t);
				for(RefSeqGeneWithIsoforms g : set1Annotations) {
					getLeftNeighbor(g, set1);
					getRightNeighbor(g, set1);
					
				}
			}
		} else {
			System.err.println(usage);
		}
		
	}
	
		
	


	



	private static void BatchGOAnalysis(String set, String refname, String gofile,
			String goformat, String outprefix, int numPerm,	int goSizeThreshold, String geneidsF) throws IOException, MathException {

		BufferedReader br = new BufferedReader(new FileReader(geneidsF));
		String line;
		while((line = br.readLine()) != null) {
			line = line.trim();
			String[] name=line.split("\t");
			GOAnalysis (set,refname,gofile,goformat,outprefix + name[1],numPerm,goSizeThreshold,name[0]);
		}		
		
	}





	//Ctor
	public NeighborAnalysis(){
		geneSet=new BEDFileParser();
		refSet=new BEDFileParser();
		unionSet=new BEDFileParser();
		neighborMap=new HashMap<RefSeqGene,RefSeqGeneWithIsoforms[] >();
		neighborDistance=new HashMap<RefSeqGene,double[]>();
		neighborJSDistance=new HashMap<RefSeqGene,double[]>();
		neighborExpCorr=new HashMap<RefSeqGene,double[]>();
		neighborExpCorrPval=new HashMap<RefSeqGene,double[]>();
		isClosestToNeighbor=new HashMap<RefSeqGene,double[]>();
		isDivergentToNeighbor=new HashMap<RefSeqGene,double[]>();
		isConvergentToNeighbor=new HashMap<RefSeqGene,double[]>();
		isSameDirectToNeighbor=new HashMap<RefSeqGene,double[]>();
		randomPvals =new ArrayList<Double> ();
	}

	public NeighborAnalysis(String setname,String refname) throws IOException{
		this();
		geneSet=new BEDFileParser(setname);
		refSet=new BEDFileParser(refname);
	}
	
	public NeighborAnalysis(BEDFileParser newGeneSet,BEDFileParser newRefSet) throws IOException{
		this();
		geneSet= newGeneSet;
		refSet=newRefSet;
	}
	
	private void updateNeighborMap(RefSeqGene g,RefSeqGeneWithIsoforms neighbor,
			RefSeqGeneWithIsoforms leftNeighbor, double rightDist,double leftDist, int isClosestRight, int isClosestLeft) {
		RefSeqGeneWithIsoforms[] arr=new RefSeqGeneWithIsoforms[2];
		arr[1]=leftNeighbor;
		arr[0]=neighbor;
		neighborMap.put(g,arr);
		double[] arr2=new double[2];
		arr2[1]=leftDist;
		arr2[0]=rightDist;
		neighborDistance.put(g,arr2);
		double[] arr3=new double[2];
		arr3[1]=isClosestLeft;
		arr3[0]=isClosestRight;
		isClosestToNeighbor.put(g,arr3);
	}

	
	private static Set<String> ReadStringList(String fname) throws IOException {
		HashSet<String> rtrn=new HashSet<String>();
		BufferedReader reader=new BufferedReader(new FileReader(fname));
		String nextLine;		
		while((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)){
			rtrn.add(nextLine);
		}
		return rtrn;
	}


	@SuppressWarnings("unchecked")
	private  void calcGeneNeighborExpCorr(String set1gct, String set2gct,String setname,String refname,
			BEDFileParser probRegionBed,boolean toLog,double addFactor,String chrSizeF, int numPerm, boolean pearsonOnDensity) throws IOException {
		
		
		Map<String, ArrayList<Double>> geneExp=ParseGCTFile.loadData(new File(set1gct));
		Map<String, ArrayList<Double>> refExp=ParseGCTFile.loadData(new File(set2gct));
		Map<String, ArrayList<Double>> geneExp2;
		Map<String, ArrayList<Double>> refExp2;
		
		int sampleNum=geneExp.get("header").size();
		refExp.remove("header");
		geneExp.remove("header");
		
		if (toLog){
			geneExp=ParseGCTFile.logTransform(geneExp,addFactor);
			refExp=ParseGCTFile.logTransform(refExp,addFactor);
			
			//For the JS
			if (addFactor<1){
				geneExp2=ParseGCTFile.loadData(new File(set1gct));
				refExp2=ParseGCTFile.loadData(new File(set2gct));
				refExp2.remove("header");
				geneExp2.remove("header");
				geneExp2=ParseGCTFile.logTransform(geneExp2,1.0);
				refExp2=ParseGCTFile.logTransform(refExp2,1.0);
			}
			else{
				geneExp2=geneExp;
				refExp2=refExp;
			}
		}
		else{
			geneExp2=geneExp;
			refExp2=refExp;
		}
		
		Map <RefSeqGene, String> geneProbMap=new HashMap <RefSeqGene, String>();
		//maybe we will remove
		for(RefSeqGene g:probRegionBed.GetGenes())
			geneProbMap.put(g,g.getName());
		
		if (this.neighborMap.size()==0)
			this.getNeighbors(new BEDFileParser(), null, false,-1);
		
		//Generate random-neighbors (control for distance between neighbors)
		NeighborsNullModel nullNeighbors=null;
		if (! chrSizeF.equals("")){
			BEDFileParser expressedRef= extractExpressedBed(refExp,this.refSet);
			nullNeighbors = new NeighborsNullModel(this.geneSet,expressedRef,0,numPerm);
			nullNeighbors.makeNaiveRandModel(chrSizeF);
		}
				
		
		for (RefSeqGene g:this.geneSet.GetGenes() ){
			if (!this.neighborMap.containsKey(g)){
				System.err.println("Not mapped " +g.getName());
				continue;
			}
			RefSeqGeneWithIsoforms left=this.neighborMap.get(g)[1];
			RefSeqGeneWithIsoforms right=this.neighborMap.get(g)[0];
			double[] corr =new double[2];
			double[] pval =new double[2];
			double[] JS =new double[2];
			if (pearsonOnDensity){
				updateCorrVal(g,left,geneExp2,refExp2,probRegionBed,geneProbMap,corr,pval,1,nullNeighbors);
				updateCorrVal(g,right,geneExp2,refExp2,probRegionBed,geneProbMap,corr,pval,0,nullNeighbors);
			}
			else{				
				updateCorrVal(g,left,geneExp,refExp,probRegionBed,geneProbMap,corr,pval,1,nullNeighbors);
				updateCorrVal(g,right,geneExp,refExp,probRegionBed,geneProbMap,corr,pval,0,nullNeighbors);
			}
			updateJSVal(g,left,geneExp2,refExp2,probRegionBed,geneProbMap,JS,1);
			updateJSVal(g,right,geneExp2,refExp2,probRegionBed,geneProbMap,JS,0);
			this.neighborExpCorr.put(g,corr);
			this.neighborExpCorrPval.put(g,pval);
			this.neighborJSDistance.put(g,JS);
		}
	
		
		
	}


	
	

	private BEDFileParser extractExpressedBed(
			Map<String, ArrayList<Double>> refExp, BEDFileParser refSet2) {

		BEDFileParser bed=new BEDFileParser();
		for (RefSeqGene g: refSet2.GetGenes())
		{
			if  (refExp.containsKey(g.getName()))
				bed.addRefSeq(g);
		}
		
		return bed;
	}


	private void updateCorrVal(RefSeqGene g, RefSeqGeneWithIsoforms neighbor,
			Map<String, ArrayList<Double>> geneExp,
			Map<String, ArrayList<Double>> refExp,
			BEDFileParser geneProbMapBed, Map<RefSeqGene, String> geneProbMap,
			double[] corr, double[] pval, int j, NeighborsNullModel nullNeighbors) {
		
		corr[j]=-2; pval[j]=-2;
		if (neighbor==null ){
			System.err.println(g.getName() +" does not have a neighbor\n");
			return;
		}
		if (!refExp.containsKey(neighbor.getName())){
			System.err.println(neighbor.getName() +" is not found in the refSet");
			return;
		}
		ArrayList<Double> neighborExp=refExp.get(neighbor.getName());
		String measureProbe="";
		IntervalTree<RefSeqGeneWithIsoforms> overlappers =geneProbMapBed.getOverlappers(g);
		if (!overlappers.isEmpty()){
			if (overlappers.size() >1)
				System.err.println(g.getName()+" overlaps " +overlappers.size()+ " different probes \n");
			Iterator <RefSeqGeneWithIsoforms> it=overlappers.valueIterator();
			while(it.hasNext()){
				RefSeqGene og= it.next();
				//String prob=geneProbMap.get(og);
				String prob=og.getName();
				prob=prob.replace("\"","");
				if(geneExp.containsKey(prob)){
					ArrayList<Double> gExp=geneExp.get(prob);
					double p=Statistics.pearsonDistance(gExp,neighborExp);
					if (new Double(p).isNaN())
						System.err.println ("Corr of "+g.getName()+ " "+neighbor.getName()+" is NAN");
					if (Math.abs(p)>Math.abs(corr[j]) || corr[j]==-2 ){
						corr[j]=p;
						measureProbe=prob;
					}
				}
				else{
					System.err.println(prob +" is not found in the gene exp mat");
				}
			}
			//Calculate pval for corr values
			if (corr[j]!= -2 & nullNeighbors!=null)
			{
				//pval[j]=calcCorrPval(corr[j],nullNeighbors,g,measureProbe,refExp,geneExp,j);
				double[] p =calcCorrPval(corr[j],nullNeighbors,g,measureProbe,refExp,geneExp,j);
				pval[j]=p[0];
				
			}
		}
		else
			System.err.println("The lincRNA "+g.getName() +" is not found in exp mat");
	}

	
	
	private double[] calcCorrPval(double corr, NeighborsNullModel nullNeighbors,RefSeqGene gene,
			String measureProbe,Map<String, ArrayList<Double>> refExp,
			Map<String, ArrayList<Double>> geneExp,int index) {

		double res=1.0;
		double total=0;
		double hits=0;
		double[] resArr=new double[2];
		LinkedList<Double> randP=new LinkedList<Double>();
		ArrayList<RefSeqGeneWithIsoforms[]> randPerm=nullNeighbors.getGeneRandNeighbors(gene);
		for (int i=0; i< randPerm.size(); i++){
			RefSeqGeneWithIsoforms n=randPerm.get(i)[index];
			ArrayList<Double> nExp=refExp.get(n.getName());
			ArrayList<Double> gExp=geneExp.get(measureProbe);
			double p=Statistics.pearsonDistance(gExp,nExp);
			
			if (! (new Double(p).isNaN())){
				total ++;
				if (p*corr > 0 & Math.abs(p) >= Math.abs(corr))
					hits++;
				this.randomPvals.add(p);//Maintain random corr distribution
				randP.add(p);
			}
		}
		res= hits/total;
		resArr[0]=res;
		
		//empirical distribution calc
		EmpiricalDistribution ED=new EmpiricalDistribution(randP);
		resArr[1]=1-(ED.getCummulativeProbability(corr));
		return resArr;
	}

	
	private static void calcRandNeighborsExpCorr(String set1gct, String out,
			boolean toLog, double addFactor, boolean pearsonOnDensity,
			int numPerm) {


		
		
	}





	private void updateJSVal(RefSeqGene g, RefSeqGeneWithIsoforms neighbor,
			Map<String, ArrayList<Double>> geneExp,
			Map<String, ArrayList<Double>> refExp,
			BEDFileParser geneProbMapBed, Map<RefSeqGene, String> geneProbMap,
			double[] JS, int j) {
		
		JS[j]=-2; 
		if (neighbor==null ){
			System.err.println(g.getName() +" does not have a neighbor\n");
			return;
		}
		if (!refExp.containsKey(neighbor.getName())){
			System.err.println(neighbor.getName() +" is not found in the refSet");
			return;
		}
		ArrayList<Double> neighborExp=refExp.get(neighbor.getName());
		double [] neighborExp_p= Statistics.normalizeToRelativeAbundance(Statistics.toDoubleArray(neighborExp));
		IntervalTree<RefSeqGeneWithIsoforms> overlappers =geneProbMapBed.getOverlappers(g);
		if (!overlappers.isEmpty()){
			if (overlappers.size() >1)
				System.err.println(g.getName()+" overlaps " +overlappers.size()+ " different probes \n");
			Iterator <RefSeqGeneWithIsoforms> it=overlappers.valueIterator();
			while(it.hasNext()){
				RefSeqGene og= it.next();
				String prob=og.getName();
				if(geneExp.containsKey(prob)){
					ArrayList<Double> gExp=geneExp.get(prob);
					double [] gExp_p= Statistics.normalizeToRelativeAbundance(Statistics.toDoubleArray(gExp));
					
					
					double p=Statistics.JSDist(gExp_p,neighborExp_p);
					if (new Double(p).isNaN())
						System.err.println ("JS of "+g.getName()+ " "+neighbor.getName()+" is NAN");
					if (Math.abs(p)<Math.abs(JS[j]) || JS[j]==-2 )
						JS[j]=p;
				}
				else{
					System.err.println(prob +" is not found in the gene exp mat");
				}
			}
		}
		else
			System.err.println("The lincRNA "+g.getName() +" is not found in exp mat");
	}


	public void getNeighbors(BEDFileParser neighborsBed,BufferedWriter bw,boolean toWrite,double distThreshold) throws IOException {
		Set<String> neighborsBedNames= new HashSet<String>(); 
		//if ((this.unionSet.GetGenes().isEmpty())){ //could be initialized in advance with other annotations
		this.unionSet.addRefSeqSet(this.refSet.GetGenes()); 
		this.unionSet.addRefSeqSet(this.geneSet.GetGenes()); 
		//}
		
		for(RefSeqGene g : this.geneSet.GetGenes()){
			RefSeqGeneWithIsoforms neighbor=getRightNeighbor(g, this.refSet);
			RefSeqGeneWithIsoforms leftNeighbor=getLeftNeighbor(g,this.refSet);
			
			
			if (neighbor!=null && !neighborsBedNames.contains(neighbor.getName())){
				neighborsBed.addRefSeq(neighbor);
				neighborsBedNames.add(neighbor.getName());
			}
			if (leftNeighbor!=null && !neighborsBedNames.contains(leftNeighbor.getName())){
				neighborsBed.addRefSeq(leftNeighbor);
				neighborsBedNames.add(leftNeighbor.getName());
			}
			int distRight=neighbor==null? 0 : (g.getStart()-neighbor.getEnd());
			int distLeft=leftNeighbor==null? 0: (leftNeighbor.getStart()-g.getEnd());
			String neighborName = neighbor==null? "" : neighbor.getName();
			String leftNeighborName = leftNeighbor==null? "" : leftNeighbor.getName();
			//Ask is closest based on the union set
			int isClosestRight=neighbor==null? 0 : this.isClosest(neighbor,g);
			int isClosestLeft=leftNeighbor==null? 0: this.isClosest(leftNeighbor,g);
						
			if (toWrite)
				bw.write (g.toBED()+"\t"+neighborName+"\t"+distRight+"\t"+leftNeighborName+"\t"+distLeft+"\n");	
			
			this.updateNeighborMap(g,neighbor,leftNeighbor,distRight,distLeft,isClosestRight,isClosestLeft);
			if (distThreshold>=0)
				this.updateDivergentPair(g,distThreshold);
		}
		
	}

	
	private  void getNeighbors(BEDFileParser neighborsBed,BufferedWriter bw,Set<String> subSet,boolean toWrite,double distThreshold) throws IOException {
		Set<String> neighborsBedNames=new HashSet<String>(); 
		boolean n1=false;
		boolean n2=false;
		int distLeft=0; int distRight = 0;
		int isClosestRight=0; int  isClosestLeft=0;
		for(RefSeqGene g : this.geneSet.GetGenes()){
			RefSeqGeneWithIsoforms neighbor=getRightNeighbor(g, this.refSet);
			RefSeqGeneWithIsoforms leftNeighbor=getLeftNeighbor(g,this.refSet);
			if (toWrite)
				bw.write (g.toBED());
			if (neighbor!=null && subSet.contains(neighbor.getName())){
				if (!neighborsBedNames.contains(neighbor.getName())){
					neighborsBed.addRefSeq(neighbor);
					neighborsBedNames.add(neighbor.getName());
				}
				distRight=neighbor==null? 0 : (g.getStart()-neighbor.getEnd());
				String neighborName = neighbor==null? "" : neighbor.getName();
				isClosestRight=neighbor==null? 0 : this.isClosest(neighbor,g);
				
				
				if (toWrite)
					bw.write ("\t"+neighborName+"\t"+distRight);
				n1=true;
			}
			else{
				if (toWrite)
					bw.write ("\t"+""+"\t"+0);
			}
			if (leftNeighbor!=null &&  subSet.contains(leftNeighbor.getName())){
				if (!neighborsBedNames.contains(leftNeighbor.getName())){
					neighborsBed.addRefSeq(leftNeighbor);
					neighborsBedNames.add(leftNeighbor.getName());
				}
				distLeft=leftNeighbor==null? 0: (leftNeighbor.getStart()-g.getEnd());
				String leftNeighborName = leftNeighbor==null? "" : leftNeighbor.getName();
				isClosestLeft=leftNeighbor==null? 0: this.isClosest(leftNeighbor,g);
				if (toWrite)
					bw.write ("\t"+leftNeighborName+"\t"+distLeft+"\n");
				n2=true;
			}
			else{
				if (toWrite)
					bw.write ("\t"+""+"\t"+0+"\n");
			}
			
			if (n1 || n2)
				this.updateNeighborMap(g,neighbor,leftNeighbor,distRight,distLeft,isClosestRight,isClosestLeft);

			if (distThreshold>=0)
				updateDivergentPair(g,distThreshold);
			
			
			
			n1=false;
			n2=false;
		}
		
	}

 private void updateDivergentPair(RefSeqGene g,double distThreshold) {
	 
	 	double[] arr=new double[2];
	    arr[0] =  isDivergentPair(g,0,distThreshold)? 1:0;
		arr[1] = isDivergentPair(g,1,distThreshold)? 1:0;
		this.isDivergentToNeighbor.put(g,arr);
		
		double[] arr2=new double[2];
	    arr2[0] =  isConvergentPair(g,0,distThreshold)? 1:0;
		arr2[1] = isConvergentPair(g,1,distThreshold)? 1:0;
		this.isConvergentToNeighbor.put(g,arr2);
		
		double[] arr3=new double[2];
	    arr3[0] =  isSameDirectPair(g,0,distThreshold)? 1:0;
		arr3[1] = isSameDirectPair(g,1,distThreshold)? 1:0;
		this.isSameDirectToNeighbor.put(g,arr3);
	}




//This was the original is closest that was based only on the gene set and not the union of the gene set and the reference set
	/*private int isClosest(RefSeqGeneWithIsoforms Neighbor, RefSeqGene gene) {
		
		RefSeqGeneWithIsoforms L= getLeftNeighbor(Neighbor,	this.geneSet);
		RefSeqGeneWithIsoforms R=getRightNeighbor(Neighbor,	this.geneSet);
		
		if(( L!= null && L.getName().equalsIgnoreCase(gene.getName()) ) || (R!= null && R.getName().equalsIgnoreCase(gene.getName())) )
			return 1;
		return 0;
	}*/

  private int isClosest(RefSeqGeneWithIsoforms Neighbor, RefSeqGene gene) {
		
		RefSeqGeneWithIsoforms L= getLeftNeighbor(Neighbor,	this.unionSet);
		RefSeqGeneWithIsoforms R=getRightNeighbor(Neighbor,	this.unionSet);
		
		if(( L!= null && L.getName().equalsIgnoreCase(gene.getName()) ) || (R!= null && R.getName().equalsIgnoreCase(gene.getName())) )
			return 1;
		return 0;
	}


	public static RefSeqGeneWithIsoforms getLeftNeighbor(RefSeqGene g,	BEDFileParser ref) {
		if (! ref.containChr(g.getChr()))
			return null;
		Iterator<Node<RefSeqGeneWithIsoforms>> iter=ref.getChrTree(g.getChr()).iterator(g.getStart(), g.getEnd());
        RefSeqGeneWithIsoforms rtrn=null;
		boolean overlap=true;
		while(iter.hasNext() && overlap){
			Node<RefSeqGeneWithIsoforms> node=(Node<RefSeqGeneWithIsoforms>)iter.next();
			rtrn=node.getValue();
			overlap=rtrn.overlapsGene(g);
		}
		return rtrn;
	}



	public static RefSeqGeneWithIsoforms getRightNeighbor(RefSeqGene g,BEDFileParser ref) {
		if (! ref.containChr(g.getChr()))
			return null;
		RefSeqGeneWithIsoforms rtrn=null;
		Iterator<Node<RefSeqGeneWithIsoforms>> iter=ref.getChrTree(g.getChr()).reverseIterator(g.getStart(), g.getEnd());
		boolean overlap=true;
		while(iter.hasNext() && overlap){
			Node<RefSeqGeneWithIsoforms> node=(Node<RefSeqGeneWithIsoforms>)iter.next();
			rtrn=node.getValue();
			overlap=rtrn.overlapsGene(g);
		}
		return rtrn;
	}


	public void updateDistanceToNeighbors() throws IOException {
		if (this.neighborMap.size()==0)
			this.getNeighbors(new BEDFileParser(), null, false,-1);
		for(RefSeqGene g:neighborMap.keySet()){
			double[]dist=new double[2];
			dist[0]= neighborMap.get(g)[0]==null? 0: (g.getStart() - neighborMap.get(g)[0].getEnd());
			dist[1]= neighborMap.get(g)[1]==null? 0: (neighborMap.get(g)[1].getStart()-g.getEnd());
			this.neighborDistance.put(g,dist);
		}
	}
	
  public RefSeqGeneWithIsoforms[]  getNeighbors (RefSeqGene g) throws IOException{
	  if (this.neighborMap.size()==0)
			this.getNeighbors(new BEDFileParser(), null, false,-1);
	  if (this.neighborMap.containsKey(g))
		  return neighborMap.get(g);
	  return null;
	
  }


public double[] getNeighborDistance(RefSeqGene g) throws IOException {
	if (this.neighborDistance.size()==0)
		this.updateDistanceToNeighbors();
	if(this.neighborDistance.containsKey(g))
		return this.neighborDistance.get(g);
	return null;
}

//for each gene report its neighbor distance to the distribution, use only the maximal value if flag 
//is set
private ArrayList<Integer> getNeighborDistances( boolean maxDist) {
	
	ArrayList<Integer> res= new ArrayList<Integer>();
	double maxd=0;
	for (RefSeqGene g:this.neighborDistance.keySet()){
		double[] d=this.neighborDistance.get(g);
		if (d==null)
			continue;
		for (int i=0; i<d.length;i++)
			if(new Double(d[i]).isNaN()) d[i]=0;
		if (maxDist){
			maxd=Math.max(d[0],d[1]);
			if (maxd>0)
				res.add((int)maxd);
		}
		else{
			for (int i=0; i<d.length;i++)
				if(d[i]>0) res.add((int)d[i]);
		}
	}
	return res;
}

private void printGeneNeighborDataTab(String fname) throws IOException {
	
	BufferedWriter bw=new BufferedWriter(new FileWriter(fname));
	//Header
	bw.write ("chr"+"\t"+"start"+"\t"+"end"+"\t"+"name"+"\t"+"score"+"\t"+"strand"+"\t"+"tStart"+"\t"+"tEnd"+"\t"+"rgb"+"\t"+"numBlocks"+"\t"+"blockSize"+"\t"+"blockStart"+"\t"+"rightNeighborName\tdistRight\tisClosestR\tcorrRight\tpvalRight\tJSR\tleftNeighborName\tdistLeft\tisClosestL\tcorrLeft\tpvalLeft\tJSL\tisDivR\tisDivL\tisConvR\tisConvL\tisSameDirectR\tisSameDirectL\n");	

	for(RefSeqGene g: this.neighborMap.keySet()){
		RefSeqGeneWithIsoforms left=this.neighborMap.get(g)[1];
		RefSeqGeneWithIsoforms right=this.neighborMap.get(g)[0];
		String leftNeighborName= left==null? "" :left.getName();
		String rightNeighborName= right==null? "" :right.getName();
		double distLeft=0;
		double distRight=0;
		double corrLeft=-2;
		double corrRight=-2;
		double pvalLeft=-1;
		double pvalRight=-1;
		double isClosestR=0;
		double isClosestL=0;
		double isDivR=0;
		double isDivL=0;
		double isConvR=0;
		double isConvL=0;
		double isSameDirectR=0;
		double isSameDirectL=0;
		double JSR=-1;
		double JSL=-1;
		if (this.neighborDistance.containsKey(g)){
			distLeft=this.neighborDistance.get(g)[1];
			distRight=this.neighborDistance.get(g)[0];
		}
		if (this.neighborExpCorr.containsKey(g)){
			corrLeft=this.neighborExpCorr.get(g)[1];
			corrRight=this.neighborExpCorr.get(g)[0];
		}
		if (this.neighborExpCorrPval.containsKey(g)){
			pvalLeft=this.neighborExpCorrPval.get(g)[1];
			pvalRight=this.neighborExpCorrPval.get(g)[0];
		}
		if (this.isClosestToNeighbor.containsKey(g)){
			isClosestR=this.isClosestToNeighbor.get(g)[0];
			isClosestL=this.isClosestToNeighbor.get(g)[1];
		}
		if (this.isDivergentToNeighbor.containsKey(g)){
			isDivR=this.isDivergentToNeighbor.get(g)[0];
			isDivL=this.isDivergentToNeighbor.get(g)[1];			
		}
		if (this.neighborJSDistance.containsKey(g)){
			JSR=this.neighborJSDistance.get(g)[0];
			JSL=this.neighborJSDistance.get(g)[1];
		}
		if (this.isConvergentToNeighbor.containsKey(g)){
			isConvR=this.isConvergentToNeighbor.get(g)[0];
			isConvL=this.isConvergentToNeighbor.get(g)[1];			
		}
		if (this.isSameDirectToNeighbor.containsKey(g)){
			isSameDirectR=this.isSameDirectToNeighbor.get(g)[0];
			isSameDirectL=this.isSameDirectToNeighbor.get(g)[1];			
		}
		bw.write (g.toBED()+"\t"+rightNeighborName+"\t"+distRight+"\t"+isClosestR+"\t"+corrRight+"\t"+pvalRight+"\t"+JSR+"\t"+leftNeighborName+"\t"+distLeft+"\t"+isClosestL+"\t"+corrLeft+"\t"+pvalLeft +"\t"+JSL+"\t"+isDivR+"\t"+isDivL+"\t"+isConvR+"\t"+isConvL+"\t"+isSameDirectR+"\t"+isSameDirectL+"\n");	
	}
	bw.close();
}

	
private void printDiveregentTranscriptionDataTab(String fname,double distThreshold) throws IOException {
	
	BufferedWriter bw=new BufferedWriter(new FileWriter(fname));
	bw.write ("chr"+"\t"+"start"+"\t"+"end"+"\t"+"name"+"\t"+"score"+"\t"+"strand"+"\t"+"tStart"+"\t"+"tEnd"+"\t"+"rgb"+"\t"+"numBlocks"+"\t"+"blockSize"+"\t"+"blockStart"+"\t"+"NeighborName\tdistRight\tisClosestR\tcorrRight\tpvalRight\tJS\n");	

	for(RefSeqGene g: this.neighborMap.keySet()){
		int ind=-1;
		if (isDivergentPair(g,0,distThreshold))
			ind=0;
		else if(isDivergentPair(g,1,distThreshold))
			ind=1;
		if (ind>-1){
			RefSeqGeneWithIsoforms nGene=this.neighborMap.get(g)[ind];
			double dist=this.neighborDistance.get(g)[ind];
			double isClosest = this.isClosestToNeighbor.get(g)[ind];
			double corr = this.neighborExpCorr.get(g)[ind];
			double pval = this.neighborExpCorrPval.get(g)[ind];
			double JS = this.neighborJSDistance.get(g)[ind];
			bw.write (g.toBED()+"\t"+nGene.getName()+"\t"+dist+"\t"+isClosest+"\t"+corr+"\t"+pval+"\t"+JS+"\n");	
		}
	}
	bw.close();
}

private void printRandomCorr(String fname) throws IOException {
	
	//Sample only 10,000 vals 
	int delta=(this.randomPvals.size()/10000);
	if(!this.randomPvals.isEmpty()){
		BufferedWriter bw=new BufferedWriter(new FileWriter(fname));
		for (int i=0; i<this.randomPvals.size(); i+=delta)
			bw.write(this.randomPvals.get(i)+"\n");
		bw.close();
	}
	
}

private boolean isDivergentPair(RefSeqGene g, int ind, double distThreshold) {
	boolean res=false;
	
	if((! this.neighborMap.containsKey(g)) || this.neighborMap.get(g)[ind] == null)
		return res;
	
	RefSeqGeneWithIsoforms nGene=this.neighborMap.get(g)[ind];
	double dist=this.neighborDistance.get(g)[ind];
	if ((g.getOrientation().equals("+") && nGene.getOrientation().equals("-")) && g.getStart() >= nGene.getStart() && dist<=distThreshold)
		res=true;
	if ((g.getOrientation().equals("-") && nGene.getOrientation().equals("+")) && g.getStart() <= nGene.getStart() && dist<=distThreshold)
		res=true;
	return res;
}

private boolean isConvergentPair(RefSeqGene g, int ind, double distThreshold) {
	boolean res=false;
	
	if((! this.neighborMap.containsKey(g)) || this.neighborMap.get(g)[ind] == null)
		return res;
	
	RefSeqGeneWithIsoforms nGene=this.neighborMap.get(g)[ind];
	double dist=this.neighborDistance.get(g)[ind];
	if ((g.getOrientation().equals("+") && nGene.getOrientation().equals("-")) && g.getEnd() <= nGene.getStart() && dist<=distThreshold)
		res=true;
	if ((g.getOrientation().equals("-") && nGene.getOrientation().equals("+")) && g.getStart() >= nGene.getEnd() && dist<=distThreshold)
		res=true;
	return res;
}

private boolean isSameDirectPair(RefSeqGene g, int ind, double distThreshold) {
	boolean res=false;
	
	if((! this.neighborMap.containsKey(g)) || this.neighborMap.get(g)[ind] == null)
		return res;
	
	RefSeqGeneWithIsoforms nGene=this.neighborMap.get(g)[ind];
	double dist=this.neighborDistance.get(g)[ind];
	if ((g.getOrientation().equals( nGene.getOrientation())) && dist<=distThreshold)
		res=true;
	
	return res;
}

private  void printDivergent(BufferedWriter bw) throws IOException {
	int cnt=0;
	int found=0;
	for (RefSeqGene g: this.isDivergentToNeighbor.keySet()){
		double a[]=this.isDivergentToNeighbor.get(g);
	   cnt++;
		if (a[0]==1 || a[1]==1){
			try {
				bw.write(g.toBED()+"\n");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				System.err.println(g.toBED());
			}
			
			found++;
		}
	}
	bw.flush();
	System.err.println(cnt);
	System.err.println(found);
	
}



	private static void GOAnalysis(String setname, String refname,
		String gofile, String goformat, String outprefix, int numPerm, int goSizeThreshold, String geneidsF ) throws IOException, MathException {
		
		BEDFileParser genes= new BEDFileParser(setname);
		BEDFileParser ref= new BEDFileParser(refname);
		
		//Extract just a subset of the genes
		if (!geneidsF.equals("")){
			genes=selectGenes(genes,geneidsF);
		}
		NeighborsNullModel na =new NeighborsNullModel (genes,ref,10,numPerm);
		na.makeBinnedIntervalModel(true);
		
		functionalAnnotation GO= new functionalAnnotation();
		GO.loadGoMap(gofile, goformat,0);
		
		GO.calcAllTermsEnrichment(na.getUniqueNeighbors(), makeRandomGeneLists(na));
		GO.writeTermsEnrichments(outprefix,goSizeThreshold,"",null);
		
		BufferedWriter bw=new BufferedWriter(new FileWriter(outprefix+"neighborList"));
		 for (RefSeqGene g:na.getUniqueNeighbors())
			 bw.write(g.getName()+"\n");
		 bw.close();
	}
	
	
	private static BEDFileParser selectGenes(BEDFileParser genes,String geneidsF) throws IOException {

		BEDFileParser res=new BEDFileParser();
		HashMap<String,RefSeqGene> map= genes.getNameGeneMap();
		BufferedReader reader=new BufferedReader(new FileReader(geneidsF));
		String name ;
		while ((name = reader.readLine()) != null ){
			if(map.containsKey(name))
				res.addRefSeq(map.get(name));
		}
		return res;
	}


	private static LinkedList<Collection<RefSeqGene>> makeRandomGeneLists(NeighborsNullModel na){
		
		LinkedList<Collection<RefSeqGene>> res=new LinkedList<Collection<RefSeqGene>>();
		
		for (int i=0; i<na.getPermNum(); i++ ){
			Set<RefSeqGene> s= new HashSet<RefSeqGene>();
			for (RefSeqGene g: na.getGeneRandNeighbors().keySet()){
				RefSeqGeneWithIsoforms[] a= na.getGeneRandNeighbors().get(g).get(i);
				if (a== null)
					continue;
				for (int j=0;j<a.length;j++){
					if(a[j] != null)
						s.addAll(a[j].getAllIsoforms());	
				}
			}
			res.add(s);
		}
		return res;
		
	}
	
	
	private static void GOTermStats(String setname, String refname,String gofile, String goformat,
			 String outprefix, int numPerm,	int goSizeThreshold, String goTerm) throws IOException {

		BEDFileParser genes= new BEDFileParser(setname);
		BEDFileParser ref= new BEDFileParser(refname);
		
		NeighborsNullModel na =new NeighborsNullModel (genes,ref,10,1);
		na.makeBinnedIntervalModel(true);
		
		functionalAnnotation GO= new functionalAnnotation();
		GO.loadGoMap(gofile, goformat,0);
		
		Set<RefSeqGene> goSet= GO.getTermGenes(goTerm,ref);
		Collection<RefSeqGene> lincNset = na.getUniqueNeighbors();
		Collection<RefSeqGene>  inter = GeneTools.intersectGeneSets(goSet, lincNset,false).GetGenes();
		
		
		BufferedWriter bw=new BufferedWriter(new FileWriter(outprefix+"neighborDist",true));
		bw.append(goTerm+"_DistFromNeighbors\t");
		ReportDistFromNeighbors(goSet,ref,bw);
		bw.append("\nlincsDistFromNeighbors\t");
		ReportDistFromNeighbors(lincNset,ref,bw);
		bw.append("\nlincsNeighborsGoIntersection_geneDistFromNeighbors\t");
		ReportDistFromNeighbors(inter,ref,bw);
		bw.append("\nlincsNeighborsGoIntersection_LincDistFromNeighbors\t");
		ReportDistFromNeighbors(inter,genes,ref,bw);//reports the distance between the lincRNA and the neighbors in inter
		bw.append("\n");
		bw.close();
		
	}


	private static void ReportDistFromNeighbors(Collection<RefSeqGene> geneset,BEDFileParser refset,	BufferedWriter bw) throws IOException {
	   
		BEDFileParser genebed=new BEDFileParser();
		genebed.addRefSeqSet(geneset);
		NeighborAnalysis na = new NeighborAnalysis(genebed,refset);
		
		na.getNeighbors(genebed, bw, false, 10000);
		ArrayList<Integer> closestDist=na.getNeighborDistances(true);
		for (int i=0; i<closestDist.size() ; i++){
			bw.append(closestDist.get(i)+",");
		}
	}


	private static void ReportDistFromNeighbors(Collection<RefSeqGene> inter,
			BEDFileParser genes, BEDFileParser ref, BufferedWriter bw) throws IOException {
		
		NeighborAnalysis na = new NeighborAnalysis(genes,ref);
		
		BEDFileParser interbed=new BEDFileParser();
		interbed.addRefSeqSet(inter);
		
		for (RefSeqGene linc: genes.GetGenes()){
			RefSeqGeneWithIsoforms[] narr= na.getNeighbors(linc);
			double[] dist=na.getNeighborDistance(linc);
			for (int i=0; i<narr.length;i++){
				if (narr[i]==null ) continue;
				if( (interbed.getExactIsoform(narr[i]) != null) & ! (new Double(dist[i]).isNaN()) & dist[i]>0 )
					bw.append((int)dist[i]+",");	
			}
		}
	}



	







	
	
	

	
	/*
	private static void getNeighbors(BEDFileParser set, BEDFileParser ref,BEDFileParser neighborsBed,BufferedWriter bw) throws IOException {
		Set<String> neighborsBedNames= new HashSet<String>(); 
		for(RefSeqGene g : set.GetGenes()){
			RefSeqGeneWithIsoforms neighbor=getRightNeighbor(g, ref);
			RefSeqGeneWithIsoforms leftNeighbor=getLeftNeighbor(g,ref);
			//System.err.println(g.toBED());
			if (neighbor!=null && !neighborsBedNames.contains(neighbor.getName())){
				neighborsBed.addRefSeq(neighbor);
				neighborsBedNames.add(neighbor.getName());
			}
			if (leftNeighbor!=null && !neighborsBedNames.contains(leftNeighbor.getName())){
				neighborsBed.addRefSeq(leftNeighbor);
				neighborsBedNames.add(leftNeighbor.getName());
			}
			int distRight=neighbor==null? 0 : (g.getStart()-neighbor.getEnd());
			int distLeft=leftNeighbor==null? 0: (leftNeighbor.getStart()-g.getEnd());
			String neighborName = neighbor==null? "" : neighbor.getName();
			String leftNeighborName = leftNeighbor==null? "" : leftNeighbor.getName();
			bw.write (g.toBED()+"\t"+neighborName+"\t"+distRight+"\t"+leftNeighborName+"\t"+distLeft+"\n");	
		}	
	}

	private static void getNeighbors(BEDFileParser set, BEDFileParser ref,BEDFileParser neighborsBed,BufferedWriter bw,Set<String> subSet) throws IOException {
		Set<String> neighborsBedNames=new HashSet<String>(); 
		for(RefSeqGene g : set.GetGenes()){
			RefSeqGeneWithIsoforms neighbor=getRightNeighbor(g, ref);
			RefSeqGeneWithIsoforms leftNeighbor=getLeftNeighbor(g,ref);
			bw.write (g.toBED());
			if (neighbor!=null && subSet.contains(neighbor.getName())){
				if (!neighborsBedNames.contains(neighbor.getName())){
					neighborsBed.addRefSeq(neighbor);
					neighborsBedNames.add(neighbor.getName());
				}
				int distRight=neighbor==null? 0 : (g.getStart()-neighbor.getEnd());
				String neighborName = neighbor==null? "" : neighbor.getName();
				bw.write ("\t"+neighborName+"\t"+distRight);
			}
			else
				bw.write ("\t"+""+"\t"+0);
			if (leftNeighbor!=null &&  subSet.contains(leftNeighbor.getName())){
				if (!neighborsBedNames.contains(leftNeighbor.getName())){
					neighborsBed.addRefSeq(leftNeighbor);
					neighborsBedNames.add(leftNeighbor.getName());
				}
				int distLeft=leftNeighbor==null? 0: (leftNeighbor.getStart()-g.getEnd());
				String leftNeighborName = leftNeighbor==null? "" : leftNeighbor.getName();
				bw.write ("\t"+leftNeighborName+"\t"+distLeft+"\n");
			}
			else
				bw.write ("\t"+""+"\t"+0+"\n");
		}
		
	}
*/
	
}
