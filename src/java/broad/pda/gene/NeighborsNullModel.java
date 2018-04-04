package broad.pda.gene;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Set;
import java.util.TreeSet;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.pda.annotation.BEDFileParser;

public class NeighborsNullModel {
	
	private BEDFileParser geneSet;
	private BEDFileParser refSet;
	private BEDFileParser unionSet; //the union of refSet and GeneSet across the genome
	private HashMap<RefSeqGene,ArrayList<RefSeqGeneWithIsoforms[]>> geneRandNeighbors;
	private HashMap<RefSeqGene,Integer> noDistanceControl;
	int permNum=1000;
	
	//Bin model
	private HashMap<RefSeqGene,Integer> neighborDistanceMap; //interval size between the 2 neighbors
	private HashMap<RefSeqGene,RefSeqGeneWithIsoforms[]> neighborMap;
	private HashMap<String,IntervalTree<NeighborsBin>>  binnedNeighborPairs;
	private IntervalTree<NeighborsBin>  binnedNeighborPairsNoChrCtrl;
	int binSize=50000;
	private int [] customBins = {0,20000,50000,100000,200000,500000,1000000,5000000,10000000,100000000,1000000000};
	
	
	//Dart model- Naive model
	private HashMap<String,Integer> chrSizes;
	private IntervalTree<chrSizeNode> chrProb;
	
	
	
	
	
	
	//CTORs
	public NeighborsNullModel(){
		geneSet=new BEDFileParser();
		refSet=new BEDFileParser();
		unionSet=new BEDFileParser();
		neighborDistanceMap=new HashMap<RefSeqGene,Integer>();
		neighborMap=new HashMap<RefSeqGene,RefSeqGeneWithIsoforms[]>();
		binnedNeighborPairs =new HashMap<String,IntervalTree<NeighborsBin>>();
		geneRandNeighbors = new HashMap<RefSeqGene,ArrayList<RefSeqGeneWithIsoforms[]>>();
		chrSizes=new HashMap<String,Integer>();
		chrProb=new IntervalTree<chrSizeNode> ();
		binnedNeighborPairsNoChrCtrl=new IntervalTree<NeighborsBin>();
	}
	
	public NeighborsNullModel(String setname,String refname,int binsize,int permnum) throws IOException{
		this();
		this.binSize=binsize;
		this.permNum=permnum;
		geneSet=new BEDFileParser(setname);
		refSet=new BEDFileParser(refname);
		unionSet.addRefSeqSet(this.geneSet.GetGenes());
		unionSet.addRefSeqSet(this.refSet.GetGenes());
		
	}
	
	public NeighborsNullModel(BEDFileParser setname,BEDFileParser refname,int binsize,int permnum) throws IOException{
		this();
		this.binSize=binsize;
		this.permNum=permnum;
		geneSet=setname;
		refSet=refname;
		unionSet.addRefSeqSet(this.geneSet.GetGenes());
		unionSet.addRefSeqSet(this.refSet.GetGenes());
		
	}
	

	//This will generate a null model where you pick neighbors that flank the same interval size as the neighbors
	//Across all genes in each permutation the distance distribution from neighbors is maintained
	//(suitable for GO permutation analysis)
	public void makeBinnedIntervalModel(boolean UseCustomBins){
		
		binRefNeighbors(UseCustomBins);
		setNeighborDistanceMap(); //updates the interval size between the neighbors of each gene, and the mapping between gene and its neighbors
		
		Integer[] errCtr= new Integer[1];
		errCtr[0]=0;
				
			for (int i=0; i<this.permNum; i++){
				for (RefSeqGene gene: this.geneSet.GetGenes() ){
				int dist=this.neighborDistanceMap.get(gene);
				if (dist<0)
					dist=1;
				//RefSeqGeneWithIsoforms[] randNeighbors= getRandNeighbors_fromBin(gene,dist);
				RefSeqGeneWithIsoforms[] randNeighbors= getRandNeighborsNoChrCtrl_fromBin(gene,dist,errCtr);
				if (! this.geneRandNeighbors.containsKey(gene))
					this.geneRandNeighbors.put(gene,new ArrayList <RefSeqGeneWithIsoforms[]>());
				this.geneRandNeighbors.get(gene).add(randNeighbors);
					
			}
		}
		System.err.println("Total number of non avaiable neighbors from correct bin in making binned interval model " + errCtr[0]);
	}
	
	//GetRandNeighbors: pick random neighboring genes that are flanking an interval in the same magnitude of the 
	//real neighbors. follow the rules:
	//1) if available- pick from the same chromosome. 2) else , pick a random chr
	//3) if the selected chromosome doesn't have intervals in this size range pick from the smallest size range within this chr
			
	private RefSeqGeneWithIsoforms[] getRandNeighbors_fromBin(RefSeqGene gene, int dist) {
		
		RefSeqGeneWithIsoforms[] res=null;
		IntervalTree<NeighborsBin> chrBins=null;
		
		if (this.binnedNeighborPairs.containsKey(gene.getChr()))
				chrBins=this.binnedNeighborPairs.get(gene.getChr());
		else{
			LinkedList<String> keys=new LinkedList<String> ( this.binnedNeighborPairs.keySet());
			int randVal=new Double(Math.random()* (keys.size()-1)).intValue();
			chrBins=this.binnedNeighborPairs.get(keys.get(randVal));
		}
		Iterator<Node<NeighborsBin>> binIt=chrBins.overlappers(dist, dist+1);
		NeighborsBin bin=null;
		if (binIt.hasNext())
			 bin=binIt.next().getValue();
		else
			bin=chrBins.findByIndex(0).getValue();
		
		if (bin.neighbors.size()==0){
			bin=chrBins.findByIndex(0).getValue();
			//System.err.println(gene.getName()+ " in distance " + dist +" dont have bin on  " + gene.getChr());
		}
		
		int index=new Double(Math.random()* bin.neighbors.size()-1).intValue();
		if (index>=0)
			res=bin.neighbors.get(index);
		else{
			//System.err.println(gene.getName()+ " in distance " + dist +" dont have bin on  "+ gene.getChr());
			
		}
		
		return res;
	}
	
	
    private RefSeqGeneWithIsoforms[] getRandNeighborsNoChrCtrl_fromBin(RefSeqGene gene, int dist, Integer[] errCtr) {
		
		RefSeqGeneWithIsoforms[] res=null;
		
		Iterator<Node<NeighborsBin>> binIt=this.binnedNeighborPairsNoChrCtrl.overlappers(dist, dist+1);
		NeighborsBin bin=null;
		if (binIt.hasNext())
			 bin=binIt.next().getValue();
		else
			{bin=this.binnedNeighborPairsNoChrCtrl.findByIndex(0).getValue(); errCtr[0]++;}
		
		if (bin.neighbors.size()==0){
			bin=this.binnedNeighborPairsNoChrCtrl.findByIndex(0).getValue(); errCtr[0]++;
		}
		
		int index=new Double(Math.random()* bin.neighbors.size()-1).intValue();
		if (index>=0)
			res=bin.neighbors.get(index);
		else{ 
			errCtr[0]++;}
		
		return res;
	}
	

	//Simple random gene model in which each gene randomly selects a gene from its chr
	public void makeNaiveRandModel(String chrSizeF) throws IOException{
		
		//if no chr size file is provided , return without making a random model
		if (chrSizeF.equals(""))
			return;
		//initialize chr sizes
		double genomeSize= initChrSize(chrSizeF);
		
		for (RefSeqGene gene: this.geneSet.GetGenes() ){
			ArrayList <RefSeqGeneWithIsoforms[]> randPerms= new ArrayList <RefSeqGeneWithIsoforms[]>();
		
			for (int i=0; i<this.permNum; i++){
				int index=new Double(Math.random()*(genomeSize-2)).intValue()+1;
				String currChr= gene.getChr();
				Iterator<Node<chrSizeNode>> t= this.chrProb.overlappers(index,index+1);
				if (  t.hasNext())
					currChr=t.next().getValue().getChr();
				IntervalTree<RefSeqGeneWithIsoforms> tree=null;
				if (! this.refSet.containChr(currChr))
					tree=refSet.getChrTree("chr1");
				else
					tree=refSet.getChrTree(currChr);
				index=new Double(Math.random()*(tree.size()-1)).intValue();
				Node<RefSeqGeneWithIsoforms> node=tree.findByIndex(index);
				RefSeqGeneWithIsoforms[] array={node.getValue(), node.getNext().getValue()};
				//System.err.println(array[0].getName()+"\t"+array[1].getName());
				randPerms.add(array);
			}
			this.geneRandNeighbors.put(gene, randPerms);
			//System.err.println();
		}
		
	}
			
	private double initChrSize(String chrSizeF) throws IOException {

	
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(chrSizeF)));
		String nextLine;
		int IsoformMissCntr=0;
		double sum=0;
		HashMap<String,Double> sizemap=  new HashMap<String,Double>();
		
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			String[] tokens=nextLine.split(" ");
			if (tokens[0].contains("rand") || tokens[0].contains("chrM")|| (! this.refSet.containChr(tokens[0])))
				continue;
			sizemap.put(tokens[0], Double.valueOf(tokens[1]));
			sum+=Double.valueOf(tokens[1]);
		}
		reader.close();
		HashMap<String,Integer> probmap= new HashMap<String,Integer>();
		for (String chr : sizemap.keySet()){
			int s=new Double((sizemap.get(chr)/sum)*100000).intValue(); ;
			probmap.put(chr, s);
		}
		int st=1;
		for (String chr : probmap.keySet()){
			int val=probmap.get(chr);
			int stop= (st +val)-1;
			chrSizeNode V=new chrSizeNode(chr,st,stop);
			this.chrProb.put (st,stop,V );
			//Iterator<Node<chrSizeNode>> testIt=this.chrProb.overlappers(st+1, st+5);
			//boolean b=testIt.hasNext();
			st=stop+1;
		}
		return (st-1);
	}

	//This will through darts on the genome to find random neighbors of each lincRNA 
	//and make null model that controls for the chances to fall in gene deserts
	public void makeDartThrowingModel(){
		
		for (int i=0; i<this.permNum; i++){
			//for every chr randomly place lincRNA but control that the  
			//TODO- complete this function if you want this to be functional
		}
	}
	
	private void binRefNeighbors(boolean UseCustomBins) {

		//Bin with no respect of chr 
		if (UseCustomBins)
			makeCustomBinsNoChrCtrl();
		//Bin by chr
		Iterator<String> chrIt=this.refSet.getChromosomeIterator();
		while(chrIt.hasNext()){
			String chr = chrIt.next();
			this.binnedNeighborPairs.put(chr,new IntervalTree<NeighborsBin>());
			//make bins
			if (UseCustomBins)
				makeCustomBins(chr);
			else
				makeEqualBins(chr);
			Iterator<RefSeqGeneWithIsoforms> treeIt=this.refSet.getChrTree(chr).valueIterator();
			if (! treeIt.hasNext())
				continue;
			RefSeqGeneWithIsoforms R=treeIt.next();
			while(treeIt.hasNext()){
				RefSeqGeneWithIsoforms L=treeIt.next();
				int dist= L.getStart()-R.getEnd();
				if (dist > 0 & ! L.overlaps(R)){
					addPairTobinnedNeighborPairs(chr,dist,R,L);
					addPairToBinnedNeighborPairsNoChrCtrl(dist,R,L);
				}
				R=L; //move to the next pair
			}
		}
		
		
	}

	
	
	private void makeCustomBins(String chr) {
		
		for (int i=0; i<this.customBins.length-1; i++){
			int start= this.customBins[i];
			int end=this.customBins[i+1];
			NeighborsBin bin= new NeighborsBin (start,end);
			this.binnedNeighborPairs.get(chr).put(start, end, bin);
		}
	}
	
	private void makeCustomBinsNoChrCtrl() {
		
		for (int i=0; i<this.customBins.length-1; i++){
			int start= this.customBins[i];
			int end=this.customBins[i+1];
			NeighborsBin bin= new NeighborsBin (start,end);
			this.binnedNeighborPairsNoChrCtrl.put(start, end, bin);
		}
	}


	private void makeEqualBins (String chr){
		for (int q=0; q<100; q++){
			int start= q*binSize +1;
			int end=(q+1) * binSize;
			NeighborsBin bin= new NeighborsBin (start,end);
			this.binnedNeighborPairs.get(chr).put(start, end, bin);
		}
	}
	
	private void addPairTobinnedNeighborPairs(String chr, int dist,
			RefSeqGeneWithIsoforms r, RefSeqGeneWithIsoforms l) {
		
		if (! this.binnedNeighborPairs.containsKey(chr))
			this.binnedNeighborPairs.put(chr,new IntervalTree<NeighborsBin>());
		Iterator<Node<NeighborsBin>> binIter=this.binnedNeighborPairs.get(chr).overlappers(dist,dist);
		if (! binIter.hasNext()){
			int q= dist / this.binSize;
			int start= q*binSize +1;
			int end=(q+1) * binSize;
			NeighborsBin bin= new NeighborsBin (start,end);
			this.binnedNeighborPairs.get(chr).put(start, end, bin);
		}
		NeighborsBin bin=binIter.next().getValue();
		bin.addPair(r,l);
		
	}

	private void addPairToBinnedNeighborPairsNoChrCtrl(int dist,
			RefSeqGeneWithIsoforms r, RefSeqGeneWithIsoforms l) {
		if (dist<0 )
			dist=0;
		Iterator<Node<NeighborsBin>> binIter=this.binnedNeighborPairsNoChrCtrl.overlappers(dist,dist);
		NeighborsBin bin=null;
		if (! binIter.hasNext())//If dist is too large, add to the last bin
			bin = this.binnedNeighborPairsNoChrCtrl.findByIndex(this.binnedNeighborPairsNoChrCtrl.size()-1).getValue();
		else
			bin=binIter.next().getValue();
		bin.addPair(r,l);
		
	}

	
	
	//updates the interval size between the neighbors of each gene, and the mapping between gene and its neighbors
	private void setNeighborDistanceMap() {
		for (RefSeqGene g: this.geneSet.GetGenes()){
			IntervalTree<RefSeqGeneWithIsoforms> tree=this.unionSet.getOverlappers(g);
			RefSeqGeneWithIsoforms L=NeighborAnalysis.getLeftNeighbor(g, this.refSet);
			RefSeqGeneWithIsoforms R=NeighborAnalysis.getRightNeighbor(g, this.refSet);
			RefSeqGeneWithIsoforms[] arr=new RefSeqGeneWithIsoforms[2];
			arr[0]=R; arr[1]=L;
			this.neighborMap.put(g, arr);
			if (L ==null)
				L=new RefSeqGeneWithIsoforms( g);
			if (R == null)
				R=new RefSeqGeneWithIsoforms( g);
			int dist = L.getStart()-R.getEnd();
			this.neighborDistanceMap.put(g,dist);
		}
		
	}

	public ArrayList<RefSeqGeneWithIsoforms[]> getGeneRandNeighbors(RefSeqGene gene) {
		if (this.geneRandNeighbors.containsKey(gene))
			return this.geneRandNeighbors.get(gene);
		return null;
	}
	
	class NeighborsBin{
		
		int start;
		int end;
		LinkedList<RefSeqGeneWithIsoforms[]> neighbors;
		
		public NeighborsBin (int st,int en){
			this.start=st;
			this.end=en;
			neighbors=new LinkedList<RefSeqGeneWithIsoforms[]>();
		}
		
		public boolean isProperSizePair(RefSeqGeneWithIsoforms right,RefSeqGeneWithIsoforms left ){
			int size=left.getStart()-right.getEnd();
			if ( size >= this.start & size <= this.end)
				return true;
			return false;
		}
		public void addPair (RefSeqGeneWithIsoforms right,RefSeqGeneWithIsoforms left ){
			RefSeqGeneWithIsoforms[] a=new RefSeqGeneWithIsoforms[2];
			a[0]=right;
			a[1]=left;
			neighbors.add(a);
		}
		
	}

	class chrSizeNode{
		String chr;
		int start;
		int stop;
		
		public chrSizeNode (String chr1,int st,int en){
			chr=chr1;
			start=st;
			stop=en;
		}
		
		public String getChr(){
			return chr;
		}
	}

	public void printRandNeighborsInfo() {

		for (RefSeqGene g:this.geneRandNeighbors.keySet()){
			ArrayList<RefSeqGeneWithIsoforms[]> a=this.geneRandNeighbors.get(g);
			for (int i=0; i< a.size(); i++)
			{	RefSeqGeneWithIsoforms[] n=a.get(i);
				int dist = n[1].getStart()- n[0].getEnd();
				System.err.println(n[0].getChr()+"\t"+n[0].getName()+"\t"+n[1].getName()+"\t" +dist);
			}
			
		}
	}

	public int getPermNum() {return this.permNum;}

	public HashMap<RefSeqGene, ArrayList<RefSeqGeneWithIsoforms[]>> getGeneRandNeighbors() {
		return this.geneRandNeighbors;
	}

	public Collection<RefSeqGene> getUniqueNeighbors() {
		Set<RefSeqGene> s= new TreeSet<RefSeqGene>();
		for (RefSeqGene g: this.neighborMap.keySet()){
			RefSeqGeneWithIsoforms[] a= this.neighborMap.get(g);
			if (a != null){
				for (int i=0; i<a.length;i++){
					if (a[i]!= null)
						s.addAll(a[i].getAllIsoforms());
				}
			}
		}
		return s;
	}
	
}
