package broad.pda.annotation;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeSet;

import broad.pda.datastructures.Alignments;
import broad.pda.gene.NeighborAnalysis;
import broad.pda.gene.RefSeqGene;
import broad.pda.gene.RefSeqGeneWithIsoforms;

public class Locus {

	private static final String SetIsoformsMap = null;
	private RefSeqGeneWithIsoforms refTranscript;
	//private IntervalTree<RefSeqGeneWithIsoforms> AllIsoforms;
	private BEDFileParser AllIsoformsBed; //This will store all the isoforms such that 2 isoforms that span the same genomic region will not override one another
	private HashMap<String,Collection<RefSeqGene>> SetIsoformMap ; 
	
	public Locus(RefSeqGeneWithIsoforms g){
		AllIsoformsBed=new BEDFileParser();
		SetIsoformMap=new HashMap<String,Collection<RefSeqGene>> ();
		refTranscript=g;
	}

	public void addIsoform(RefSeqGeneWithIsoforms iso,String set) {
		//AllIsoforms.put(iso.getStart(), iso.getEnd(), iso);
		Collection <RefSeqGene>allIso=iso.getAllIsoforms();
		AllIsoformsBed.addRefSeqSet(allIso);
		if (!SetIsoformMap.containsKey(set)){
			Collection<RefSeqGene> c=new LinkedList<RefSeqGene>();
			SetIsoformMap.put(set,c);
		}
		SetIsoformMap.get(set).addAll(allIso);
	}
	
	public RefSeqGeneWithIsoforms getReferenceTranscript(){
		return this.refTranscript;
	}
	public HashMap<String,Integer> getNumIsoformsPerSet(){
		
		HashMap<String,Integer> res=new HashMap<String,Integer>();
		for (String s: SetIsoformMap.keySet()){
			res.put(s,this.SetIsoformMap.get(s).size() );
		}
		
		return res;
	}
	
	public HashMap<String,Double> getMaxScorePerSet(){
		HashMap<String,Double> res=new HashMap<String,Double>();
		for (String s: SetIsoformMap.keySet()){
			Collection<RefSeqGene> isoforms=SetIsoformMap.get(s);
			Double maxScr=0.0;
			boolean flag=true;
			for (RefSeqGene iso: isoforms){
				if (flag) {maxScr=iso.getBedScore(); flag=false;}
				else {maxScr=Math.max(maxScr, iso.getBedScore());}
			}
			res.put(s,maxScr);
		}
		return res;
	}
	
	public HashMap<String,Integer> getIsComaptibleWithRefPerSet(){
		HashMap<String,Integer> res=new HashMap<String,Integer>();
		for (String s: SetIsoformMap.keySet()){
			Collection<RefSeqGene> isoforms=SetIsoformMap.get(s);
			Integer maxScr=0;
			for (RefSeqGene iso: isoforms){
				if (refTranscript.numOfCompatibleIntrons(iso)== refTranscript.getNumExons()-1)
					maxScr++;
			}
			res.put(s,maxScr);
		}
		return res;	
	}
	
	public HashMap<String,Integer> getIsPartiallyComaptibleWithRefPerSet(){
		HashMap<String,Integer> res=new HashMap<String,Integer>();
		for (String s: SetIsoformMap.keySet()){
			Collection<RefSeqGene> isoforms=SetIsoformMap.get(s);
			Integer maxScr=0;
			for (RefSeqGene iso: isoforms){
				if (refTranscript.numOfCompatibleIntrons(iso)>0)
					maxScr++;
			}
			res.put(s,maxScr);
		}
		return res;	
	}
	
	public int getIsComaptibleBetweenSubsets(ArrayList<String> set1Lst,ArrayList<String> set2Lst,boolean isFully){
		
		Iterator<String> set1It=set1Lst.iterator();
		int compatible=0;
		while (set1It.hasNext()){
			String set1name=set1It.next();
			if (! SetIsoformMap.containsKey(set1name))
				{ //System.err.println(set1name);
				continue;}
			Collection<RefSeqGene> isoforms1=SetIsoformMap.get(set1name);
			for (RefSeqGene iso1:isoforms1){
					Iterator<String> set2It=set2Lst.iterator();
					while (set2It.hasNext()){
						String set2name=set2It.next();
						if (! SetIsoformMap.containsKey(set2name))
						{ //System.err.println(set2name);
						continue;}
						Collection<RefSeqGene> isoforms2=SetIsoformMap.get(set2name);
						for (RefSeqGene iso2: isoforms2){
								if (isFully && isFullyComaptible(iso1,iso2,false))
									compatible++;
								if ((isFully==false) && (iso1.numOfCompatibleIntrons(iso2)>0))
									compatible++;
							}
						}
				}
			}
	return compatible;			
	}
	
	public int getIsComaptibleBetweenAny2Sets(ArrayList<String> setLst,	boolean isFully) {
		int compatible=0;
		for (int i=0; i< setLst.size()-1;i++){
			String set1name=setLst.get(i);
			if (! SetIsoformMap.containsKey(set1name))
				continue;
			Collection<RefSeqGene> isoforms1=SetIsoformMap.get(set1name);
			for (RefSeqGene iso1:isoforms1){
					for (int j=i+1; j<setLst.size(); j++ ){
						String set2name=setLst.get(j);
						if (! SetIsoformMap.containsKey(set2name))
							continue;
						Collection<RefSeqGene> isoforms2=SetIsoformMap.get(set2name);
						for (RefSeqGene iso2: isoforms2){
								if (isFully && isFullyComaptible(iso1,iso2,false))
									compatible++;
								if ((isFully==false) && (iso1.numOfCompatibleIntrons(iso2)>0))
									compatible++;
						}
					}
			}
				
		}
		
		return compatible;		
	}
	
	
	public int getIsComaptibleBetweenSubsets(ArrayList<String> set1Lst,
			ArrayList<String> set2Lst, BEDFileParser uniqIsoBed) {
		
		Iterator<RefSeqGeneWithIsoforms> unqRefsIt = uniqIsoBed.getOverlappers(this.refTranscript).valueIterator();
		int compatible=0;
		int i=0;
		ArrayList<Integer> compNum1=new ArrayList<Integer>();
		ArrayList<Integer> compNum2=new ArrayList<Integer>();
		//go over all unique reference isoforms
		while (unqRefsIt.hasNext()){
			Collection <RefSeqGene> refs = unqRefsIt.next().getAllIsoforms();
			for (RefSeqGene ref:refs){
				compNum1.add(i, numCompatibleInSet(set1Lst,ref));
				compNum2.add(i, numCompatibleInSet(set2Lst,ref));
				i++;
			}
		}
		
		for (int j=0; j<compNum1.size(); j++){
			int a=compNum1.get(j);
			int c= compNum2.get(j);
			if (a>0 & c >0){
				compatible=Math.max(compatible, a+c);
			}
		}
		return compatible;
		
	}
	
	public int getIsComaptibleBetweenAny2Sets(ArrayList<String> set3Lst,
			 BEDFileParser uniqIsoBed) {
		Iterator<RefSeqGeneWithIsoforms> unqRefsIt = uniqIsoBed.getOverlappers(this.refTranscript).valueIterator();
		int compatible=0;
		int i=0;
		ArrayList<Integer> compNum1=new ArrayList<Integer>();
		while (unqRefsIt.hasNext()){
			Collection <RefSeqGene> refs = unqRefsIt.next().getAllIsoforms();
			for (RefSeqGene ref:refs){
				compNum1.add(i, numCompatibleInSet(set3Lst,ref));
				i++;
			}
		}
		for (int j=0; j<compNum1.size(); j++)
			compatible=Math.max(compatible, compNum1.get(j));
		
		if (compatible >1)
			return compatible;
		else
			return 0;
	}
	
    private int numCompatibleInSet(ArrayList<String> set1Lst, RefSeqGene ref) {
    	
    	int res=0;
    	Iterator<String> set1It=set1Lst.iterator();
		while (set1It.hasNext()){
			String set1name=set1It.next();
			if (! SetIsoformMap.containsKey(set1name))
				continue;
			Collection<RefSeqGene> isoforms1=SetIsoformMap.get(set1name);
			for (RefSeqGene iso1:isoforms1){
				if (isFullyComaptible(ref,iso1,false))
					res++;	
			}
		}
		return res;
	}

	
	
		
	public RefSeqGene getExonUnion(){
		
		boolean first=true;
		RefSeqGene mergedElement=refTranscript;
		
		if (! this.SetIsoformMap.isEmpty()){
			Iterator<RefSeqGeneWithIsoforms> it=AllIsoformsBed.getChrTree(this.refTranscript.getChr()).valueIterator();
			
			RefSeqGene tmp;
			while(it.hasNext()){
				for (RefSeqGene iso: it.next().getAllIsoforms()){
					if (first) {mergedElement=iso; first=false;}
					else {mergedElement=  mergedElement.takeUnion(iso);}
				}
			}
		}
		return mergedElement;
	}
	
	//Output: First element specifies if genes are in the same orientation, and second element specifies the distance
	public double[] getDistanceTo3primeNeighbor(BEDFileParser geneSet)throws IOException{
		double[]res=new double[2];
		RefSeqGene consenzus=getExonUnion();
		BEDFileParser myBed=new BEDFileParser();
		myBed.addRefSeq(consenzus);
		NeighborAnalysis mNA=new NeighborAnalysis(myBed,geneSet);
		mNA.getNeighbors( new BEDFileParser(), null,false,-1);
		mNA.updateDistanceToNeighbors();
		int ix=1;
		if (consenzus.getOrientation().equalsIgnoreCase("-"))
			ix=0;
		RefSeqGene leftNeighbor= mNA.getNeighbors(consenzus)[ix];
		if (leftNeighbor!= null) 
			res[0]= leftNeighbor.getOrientation().equalsIgnoreCase( consenzus.getOrientation())? 1:0;
		else
			res[0]=0;
		res[1]=mNA.getNeighborDistance(consenzus)[ix];
		
		return res;		
		
	}
	
	//Output: First element specifies if genes are in the same orientation, and second element specifies the distance
	public double[] getDistanceTo5primeNeighbor(BEDFileParser geneSet)throws IOException{
		double[]res=new double[2];
		RefSeqGene consenzus=getExonUnion();
		BEDFileParser myBed=new BEDFileParser();
		myBed.addRefSeq(consenzus);
		NeighborAnalysis mNA=new NeighborAnalysis(myBed,geneSet);
		mNA.getNeighbors( new BEDFileParser(), null,false,-1);
		mNA.updateDistanceToNeighbors();
		int ix=0;
		if (consenzus.getOrientation().equalsIgnoreCase("-"))
			ix=1;
		RefSeqGene Neighbor= mNA.getNeighbors(consenzus)[ix];
		if (Neighbor!= null) 
			res[0]= Neighbor.getOrientation().equalsIgnoreCase( consenzus.getOrientation())? 1:0;
		else
			res[0]=0;
		res[1]=mNA.getNeighborDistance(consenzus)[ix];
		
		return res;		
		
	}
	
	public double[]distanceToOppositeStrand5primeNeighbor (BEDFileParser geneSet) throws IOException{
		
		double[]res=this.getDistanceTo5primeNeighbor(geneSet);
		if (res[0]!=0){ //Have the same orientation
			res[0]=0;
			res[1]=0;
		}
		return res;
	}

	public double getMaxScoreAcrossAllIso() {

		double maxScr=0;
		HashMap<String, Double> setMaxScoreMap=this.getMaxScorePerSet();
		Collection<Double> vals=setMaxScoreMap.values();
		for (Double d: vals){
			maxScr=Math.max(d,maxScr);
		}
		return maxScr;
	}

	//adds a suffix to the name of every isoform that specifies the loci position
	public void updateIsoformsNameWithLocusName() {
		
		String refName=";"+this.refTranscript.getChr()+":"+this.refTranscript.getStart()+"-"+this.refTranscript.getEnd();
		if (this.SetIsoformMap.isEmpty())
			return;
		
		/* This are supposed to be the same objects as in the sets themselves- with the exception of the transcript that initiates Super
		*/Iterator<RefSeqGeneWithIsoforms> it =this.AllIsoformsBed.getChrTree(this.refTranscript.getChr()).valueIterator();
		while(it.hasNext()){
			RefSeqGeneWithIsoforms g=it.next();
			g.addSuffixToName(refName);
			
		}
		
		
		for (String s:this.SetIsoformMap.keySet()){
			for (RefSeqGene g : this.SetIsoformMap.get(s))
				g.addSuffixToName(refName);
		}
		
	}

	public Collection<RefSeqGene> getSetIsoforms(String setName) {
		
		Collection<RefSeqGene> isoforms=new LinkedList<RefSeqGene> ();
		if ( SetIsoformMap.containsKey(setName))
			 isoforms=SetIsoformMap.get(setName);
		else
			System.err.println(setName);
	
		return isoforms;
	}

	public Collection<RefSeqGene> getAllIsoforms(String excludeSet) {
		
		LinkedList<RefSeqGene> lst= new LinkedList<RefSeqGene> ();
		if (! this.SetIsoformMap.isEmpty()){
			for (String s:this.SetIsoformMap.keySet()){
				if (!s.equalsIgnoreCase(excludeSet))
					lst.addAll(this.SetIsoformMap.get(s));
			}
		}
		return lst;
	}

	public Collection<RefSeqGene> getComaptibleBetweenAny2Sets(ArrayList<String> setLst, boolean isFully) {
		
		HashSet<RefSeqGene> res=new HashSet<RefSeqGene>();
		HashMap<String,RefSeqGene> map= new HashMap<String,RefSeqGene>();
		for (int i=0; i< setLst.size()-1;i++){
			String set1name=setLst.get(i);
			//System.err.println(set1name);
			if (! SetIsoformMap.containsKey(set1name))
				continue;
			Collection<RefSeqGene> isoforms1=SetIsoformMap.get(set1name);
			for (RefSeqGene iso1:isoforms1){
					for (int j=i+1; j<setLst.size(); j++ ){
						String set2name=setLst.get(j);
						if (! SetIsoformMap.containsKey(set2name))
							continue;
						Collection<RefSeqGene> isoforms2=SetIsoformMap.get(set2name);
						for (RefSeqGene iso2: isoforms2){
								if (isFully && isFullyComaptible(iso1,iso2,false)){
									map.put(iso1.getName(),iso1);
									//System.err.println("Found compatible " + iso1.getName()+ "From " +set1name +" and "+ set2name + "\n" + iso1.toBED());
								}
								if ((isFully==false) && (iso1.numOfCompatibleIntrons(iso2)>0)){
									map.put(iso1.getName(),iso1);
								}
						}
						
					}
			}
		}
		for(String s:map.keySet())
			res.add(map.get(s));
		return res;
	}

	public Collection<RefSeqGene> getUnqComaptibleBetweenAny2Sets(
			ArrayList<String> setLst, Collection<? extends RefSeqGene> unqIsoSet, boolean isFully) {
		
		LinkedList<RefSeqGene> candidates =new LinkedList <RefSeqGene>();
		for (RefSeqGene g:unqIsoSet){
			int scr=this.getNumberOfCompatibleSets(g,setLst,isFully);
			if (  scr>=2){
				g.setBedScore(scr);
				candidates.add(g);
			}
		}
				
		return candidates;
	}

	public Collection<RefSeqGene> getUnqComaptibleBetweenSubsets( ArrayList<String> set1Lst,ArrayList<String> set2Lst,
			Collection<RefSeqGene> unqIsoSet, boolean isFully) {
		
		LinkedList<RefSeqGene> candidates =new LinkedList <RefSeqGene>();
		for (RefSeqGene g:unqIsoSet){
			int scr1=this.getNumberOfCompatibleSets(g,set1Lst,isFully);
			int scr2=this.getNumberOfCompatibleSets(g,set2Lst,isFully);
			if ( scr1>0 && scr2>0){
				g.setBedScore(scr1+scr2);
				candidates.add(g);
			}
		}		
		return candidates;
	}
	
	public Collection<RefSeqGene> getComaptibleIsoforms(RefSeqGene g,boolean isFully) {
		
		Collection<RefSeqGene> lst=new LinkedList<RefSeqGene>();
		for (String setName:SetIsoformMap.keySet()){
			Collection<RefSeqGene> isoforms=SetIsoformMap.get(setName);
			for (RefSeqGene iso: isoforms){
				if (isFully && (g.numOfCompatibleIntrons(iso)==g.getNumExons()-1))
					lst.add(iso);
				if ((isFully==false) && (g.numOfCompatibleIntrons(iso)>0))
					lst.add(iso);
			}
		}
		return lst;
	}
	

	private int getNumberOfCompatibleSets(RefSeqGene g,ArrayList<String> setLst, boolean isFully) {
		
		int totalSet=0;
		for (String setName:setLst){
			if (! SetIsoformMap.containsKey(setName))
				continue;
			Collection<RefSeqGene> isoforms=SetIsoformMap.get(setName);
			boolean inSet=false;
			for (RefSeqGene iso: isoforms){
				if (isFully && (g.numOfCompatibleIntrons(iso)==g.getNumExons()-1))
					inSet=true;
				if ((isFully==false) && (g.numOfCompatibleIntrons(iso)>0))
					inSet=true;
			}
			if (inSet)
				totalSet++;
		}
		return totalSet;
		
		
	}
		
	
    public static Collection<RefSeqGene> SelectLongestIntronChainCandidate(Collection<RefSeqGene> betweenSubsetsIsoSet) {
			
    		int length=0;
    		Collection<RefSeqGene> lst= new LinkedList<RefSeqGene>();
			for (RefSeqGene g: betweenSubsetsIsoSet )
				length=Math.max(length,g.getNumExons());
			for (RefSeqGene g: betweenSubsetsIsoSet ){
				if (g.getNumExons()==length)
					lst.add(g);
			}
			
			return lst;
	}

	public static boolean isFullyComaptible(RefSeqGene iso1,RefSeqGene iso2,boolean stringent){
		int num1=iso1.numOfCompatibleIntrons(iso2);
		int num2=iso2.numOfCompatibleIntrons(iso1);
		boolean or =(num1 ==(iso1.getNumExons()-1) || num2 ==(iso2.getNumExons()-1) );
		boolean and=(num1 ==(iso1.getNumExons()-1) && num2 ==(iso2.getNumExons()-1) );
		if (stringent)
			return and;
		return or ; 
	}

	//RefSeq is a name of one set  stored in the locus that serves as 
	// the set of unique reference isoforms (like a cuff compare result)
	public RefSeqGene getLongestIntornChainIso(String refSet) {
		RefSeqGene res=null;
		/*if (! this.SetIsoformMap.containsKey(refSet)){
			System.err.println("No key , Ref set is: "+refSet + "keys are");
			Set<String> tmp=this.SetIsoformMap.keySet();
			for (String s:tmp)
				System.err.println(s);
		}*/
		if (refSet != null && this.SetIsoformMap.containsKey(refSet))
			res= longestIntornChainIso(this.SetIsoformMap.get(refSet));
		else{
			Collection <RefSeqGene> c=new LinkedList<RefSeqGene> ();
			for (String s:this.SetIsoformMap.keySet()){
				c.add(longestIntornChainIso(this.SetIsoformMap.get(s)));
			}
			res=longestIntornChainIso(c);
		}
		return res;
	}

	private RefSeqGene longestIntornChainIso(Collection<RefSeqGene> collection) {
		RefSeqGene res=null;
		HashMap<RefSeqGene,Double> map= new HashMap<RefSeqGene,Double> ();
		for (RefSeqGene g: collection)
			map.put(g,new Double(g.getNumExons()));
		Collection <RefSeqGene> longest =selectMaxScrIso (map,-1);
		HashMap<RefSeqGene,Double> map2= new HashMap<RefSeqGene,Double> ();
		for (RefSeqGene g: longest) //from the longest intron chain select the one that has the largest size
			map2.put(g,new Double(g.getTranscriptLength()));
		LinkedList <RefSeqGene> longest2 =(LinkedList<RefSeqGene>) selectMaxScrIso (map2,-1);
		if (longest2.size()>0)
			res=longest2.get(0);
		return res;
	}
		
	private static Collection <RefSeqGene> selectMaxScrIso(HashMap<RefSeqGene,Double> map,double initScr) {
		double scr=initScr;
		Collection <RefSeqGene>res=new LinkedList<RefSeqGene>(); 
		for (RefSeqGene g: map.keySet() )
			scr=Math.max(scr,map.get(g));
		for (RefSeqGene g: map.keySet() ){
			if (map.get(g)==scr)
				res.add(g);
		}
		return res;
	}

	public RefSeqGene getMerged() {
		
		List<RefSeqGene> genes= getAllMerged();
		return genes.get(0);
	}

	public List<RefSeqGene> getAllMerged() {
		
		BEDFileParser bed=new BEDFileParser();
		bed.addRefSeqSet(this.AllIsoformsBed.GetGenes());
		bed.merge();
		List<RefSeqGene> genes=bed.GetGenes();
		if (genes.size()>1){
			for (int i=0; i<genes.size();i++)
				System.err.println("merged in fragments: " + genes.get(i).toBED());
		}
		return genes;
	}


	public int numOfOverlapSets(ArrayList<String> setLst) {
		int res=0;
		for (int i=0; i<setLst.size() ;i++){
			boolean b=this.SetIsoformMap.containsKey(setLst.get(i));
			if (this.SetIsoformMap.containsKey(setLst.get(i)) &&  this.SetIsoformMap.get(setLst.get(i)).size()>0)
				res++;
		}
		return res;
	}

	public double getMaxPctCovered() {
		return getMaxCovered(0);
	}

	public double getMaxExonsCovered() {
		return getMaxCovered(1);
	}

	public double getMaxPctGenomeCovered() {
		return getMaxCovered(2);
	}
	
	private double getMaxCovered(int i){
		HashMap<RefSeqGene, ArrayList<Double>> overlapMap= pctOverlapWithRef();
		double maxOver=0;
		for (RefSeqGene g:overlapMap.keySet() )
			maxOver=Math.max(maxOver,overlapMap.get(g).get(i));
		return maxOver;
	}
	
	
	public HashMap<RefSeqGene, ArrayList<Double>> pctOverlapWithRef()
	{
		HashMap<RefSeqGene, ArrayList<Double>> res=new HashMap<RefSeqGene, ArrayList<Double>>();
	
		for (String s: SetIsoformMap.keySet()){
			Collection<RefSeqGene> isoforms=SetIsoformMap.get(s);
			ArrayList<Double> arr= new ArrayList<Double>();
			for (RefSeqGene iso: isoforms){
				arr.add(0,this.refTranscript.percentOverlapping(iso));
				arr.add(1,new Double(this.refTranscript.getMerged().numOverlappingExons(iso)));
				arr.add(2,this.refTranscript.getMerged().percentGenomeOverlapping(iso));
				res.put(iso,arr);
			}
			
		}
		return res;	
	
	}

	//TODO- decide how to implement this
	//Bias will be determined if only the first or last 50% of exons are covered 
	public int isTerminiCoverageBias() {
		int bias=0;
		RefSeqGene ref= this.getReferenceTranscript();
		int numEx=ref.getNumExons();
		Alignments[] aln = ref.getExons();
		for (String s: SetIsoformMap.keySet()){
			Collection<RefSeqGene> isoforms=SetIsoformMap.get(s);
			for (RefSeqGene iso: isoforms){
				int first=0;
				int last=0;
				for (int i=0;i<numEx; i++){
					
				}
				
			}
		}

		return bias;
	}

	public Collection<? extends RefSeqGene> getAllIsoforms() {
		//Don't exclude any set
		return getAllIsoforms("");
	}

	public Integer getNumberOfIsoforms() {
		
		return this.AllIsoformsBed.getNumberOfIsoforms();
	}

	public Collection<? extends RefSeqGene> selectRandIsoSubset(
			Integer maxIsoPerLoci) {
		
		
		Collection<RefSeqGene> rtrn = new TreeSet <RefSeqGene>();
		int i=0;
		for (RefSeqGene g: this.getAllIsoforms()){
			if (i< maxIsoPerLoci){
				rtrn.add(g); i++;
			}
			else
				break;
		}
			
		return rtrn;
	
	}

	
	
}
