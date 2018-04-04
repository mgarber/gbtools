package broad.pda.geneexpression.dge;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.annotation.BED;
import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.annotation.GTFFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.gene.RefSeqGeneWithIsoforms;
import broad.pda.seq.alignment.AlignmentUtils;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;

/**
 * 
 * @author skadri
 *
 */
public class APA {
	
	private static final int WINDOW = 5;
	private static final int NEARNESS_THRESHOLD = 25;
	private static final int MAGNITUDE_THRESHOLD = 5;
	private static final double QUANTILE = 0.8;
	static final String usage = "Usage: APA -task <task name> "+
			"\n\tcalculate: Computes the alternative 3' ends of annotated genes" +
			"\n\t\t-annotations <Specific genes to calculate APA [BED by default]> "+
			"\n\t\t-alignment <Alignment (mapped to genome) [bam or sam file]> "+
			"\n\t\t-out <Output Filename>";
	
	static Logger logger = Logger.getLogger(APA.class.getName());
	
	
	public APA(String[] args)throws IOException{
		
		Globals.setHeadless(true);
		
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"calculate");
		/*
		 * Read the annotation file
		 */
		String annotationsFile = argMap.getMandatory("annotations");
		/*
		 * Check the format of the annotation files and call the GTF or BED parser accordingly
		 */
		BEDFileParser annotations =  annotationsFile.endsWith(".gtf") || annotationsFile.endsWith(".GTF")? new GTFFileParser(annotationsFile) : new BEDFileParser(annotationsFile);
		/*
		 * Read the name of the alignment file
		 */
		String alignmentFile = argMap.getMandatory("alignment");
		/*
		 * Initialize the data model using the alignment file
		 * @param: <alignment flieName> <load_chromosome_stats> <minMappingQuality> <remove_PCR_duplicates> <weigh reads by NH flag>
		 */
		ContinuousDataAlignmentModel libDataModel = AlignmentUtils.loadAlignmentData(alignmentFile,false,5,false,true);
		
		/*
		 * Iterate through list of annotations. Iterator is over chromosome names.
		 */
		Iterator<String> chrom_iter = annotations.getChromosomeIterator();
		
		BufferedWriter bw = argMap.getOutputWriter();
		/*
		 * For each chromosome
		 */
		while(chrom_iter.hasNext()) {
			String chr = chrom_iter.next();
			/*
			 * If the alignment data has data from that chromosome
			 */
			if(libDataModel.hasDataForChromosome(chr)){
				logger.info("Processing " + chr);
				/*
				 * Obtain an iterator over the interval tree values built from the annotations (genes) on that particular chromosome.
				 * IntervalTree of RefSeqGeneWithIsoforms
				 */
				Iterator<RefSeqGeneWithIsoforms> annotation_iter = annotations.getChrTree(chr).valueIterator();
				/*
				 * Parse annotation tree
				 * While there is an annotated RefSeqGeneWithIsoforms in Interval tree to analyze
				 * Thus, for each gene
				 */
				
				while(annotation_iter.hasNext()){
					/*
					 * annotation = current node of RefSeqGeneWithIsoforms in interval tree
					 */
					RefSeqGeneWithIsoforms annotation_with_isoforms = annotation_iter.next();
					Iterator<RefSeqGene> isoform_iter = annotation_with_isoforms.getAllIsoforms().iterator();
					
					//System.out.println("Processing "+annotation_with_isoforms.getName());
					boolean firstIsoform=true;
					int farthestEnd = 0;
					int farthestStart =0;
					RefSeqGene farthestGene=null;
					//FOR EACH ANNOTATION (WITH ISOFORMS)
					while(isoform_iter.hasNext()){
						
						RefSeqGene annotation = isoform_iter.next();
						
						if(firstIsoform){
							
							if(annotation.getOrientation().equals("-")){
								farthestEnd = annotation.getStart();
								farthestStart = annotation.getEnd();
							}
							else{
								farthestEnd = annotation.getEnd();
								farthestStart = annotation.getStart();
							}
							farthestGene=annotation;
							firstIsoform=false;
						}
						// CALCULATE INTERVALS FOR THE ISOFORM ENDS
						//Get farthest end of isoform
						if(annotation.getOrientation().equals("-")){
							if(annotation.getStart()<farthestEnd){
								farthestEnd = annotation.getStart();
								farthestGene = annotation;
							}
							if(annotation.getEnd()>farthestStart){
								farthestStart = annotation.getEnd();
							}
						}
						else{
							if(annotation.getEnd()>farthestEnd){
								farthestEnd = annotation.getEnd();
								farthestGene = annotation;
							}
							if(annotation.getStart()<farthestStart){
								farthestStart = annotation.getStart();
							}
						}
					}
					//System.out.println("Farthest Gene is "+farthestGene.getName());
					
					PolyASite previousSite = new PolyASite();
					
					Vector<PolyASite> sites = findPolyASites(libDataModel,farthestGene,farthestStart, farthestEnd);
					//Starting from farthest 3' end of gene, move to farthest gene start
					//Find a ranked list of polyA sites with expression above a certain quantile
					//Negative strand
					//System.out.println(farthestStart+" "+farthestEnd);
				/*	if(farthestEnd<farthestStart){
						//GET DATA FOR ALIGNMENT BETWEEN FARTHEST START AND END
						int i=farthestEnd;
						while(i<farthestStart){
							PolyASite site = findNextHill(i,chr,libDataModel,farthestStart);
							
							//If the previous hill is sufficiently far from current hill
							if(!passNearnessTest(previousSite,site)){
								//Add the site to the vector of sites
								sites.add(site);
							}
							else{
								//Dont add the site
								//For now, we will consider the start of the first hill in the closely places hills as the site
							}
							previousSite = site;
							i = site.getPosition();
						}
					}
					// positive strand
					//WRONG
					else{
						double[] data = libDataModel.getDataForAlignment(new Alignments(farthestGene.getChr(),farthestStart,farthestEnd));
						
						int i=farthestStart;
						while(i<farthestEnd){
							PolyASite site = findNextHill(i,chr,libDataModel,farthestEnd);
							
							//If the previous hill is sufficiently far from current hill
							if(!passNearnessTest(previousSite,site)){
								//Add the site to the vector of sites
								sites.add(site);	
							}
							else{
								//Dont add the site
								//For now, we will consider the start of the first hill in the closely places hills as the site
							}
							previousSite = site;
							i = site.getPosition();
						}
					}*/
					
					writeToBedFile(bw,sites,farthestGene);
				}
			}
		}
		bw.close();
	}

	public void writeToBedFile(BufferedWriter bw,Vector<PolyASite> sites,RefSeqGene g)  throws IOException{
		
		int rank=1;
		//CONVERT POLYA SITES TO A BED ENTRY
		
		for(PolyASite s:sites){
			//if(s.getMagnitude()>MAGNITUDE_THRESHOLD){
				BED entry = new BED((g.getName()+"_"+rank),g.getChr(),s.getPosition(),s.getPosition());
				entry.setScore(s.getMagnitude());
				rank++;
				bw.write(entry.toString());
				bw.newLine();
				
		//	}
			
			//System.out.println(entry.toString());
		}
	}
	
	private Vector<PolyASite> findPolyASites(ContinuousDataAlignmentModel model,RefSeqGene gene, int start, int end) throws IOException{
		
		Vector<PolyASite> sites = new Vector<PolyASite>();
		PolyASite previousSite = new PolyASite();
		
		List<Double> gene_list = model.getDataForGene(gene);
		//System.out.println("Data for gene imported");
		Collections.sort(gene_list);
		//System.out.println("gene length = "+gene_list.size());
		//System.out.println("Calculating cutoff");
		
		//Quantile calculated based on gene transcript not genomic region
		double cutoff = Math.max(MAGNITUDE_THRESHOLD, Statistics.quantile(gene_list, QUANTILE));
		//System.out.println("Cutoff calculated");
		
		//negative strand
		//go from end to start
		if(end<start){
			double[] data_arr = model.getDataForAlignment(new Alignments(gene.getChr(),end,start));
			//System.out.println(data_arr.length+" "+gene.getGenomicLength()+" "+gene.getTranscriptLength()+" "+(start-end));
			//List<Double> data_list=array2List(data_arr);
			data_arr = subtract(data_arr, Math.max(model.getLambda(gene.getChr()),cutoff));
			int i=end;
			int index=0;
			while(i<start){ 
				//System.out.println(end+" index= "+index+" "+start);
				/**
				 * index: current position in the array
				 * data_arr: genomic region to parse
				 * end: the starting point of the gene: here, end because negative strand
				 * "downstream": because on negative strand going from 3' end to 5'end of gene
				 */
				PolyASite site = findNextHill (index,data_arr,end,"downstream");
				//If the previous hill is sufficiently far from current hill
				if(!passNearnessTest(previousSite,site)&& site.getMagnitude()>MAGNITUDE_THRESHOLD){
					//Add the site to the vector of sites
					sites.add(site);	
				}
				else{
					//Dont add the site
					//For now, we will consider the start of the first hill in the closely places hills as the site
				}
				previousSite = site;
				i = site.getPosition();
				index = i-end+1;
				//System.out.println("Number of sites:"+sites.size());
			}
		}
		//positive strand. GO from end to start but higher to lower position
		else{
			double[] data_arr = model.getDataForAlignment(new Alignments(gene.getChr(),start,end));
			//List<Double> data_list=array2List(data_arr);
			//System.out.println("Data imported");
			double lambda = model.getLambda(gene.getChr());
			//System.out.println(lambda);
			double max = Math.max(lambda,cutoff);
			//System.out.println(max);
			data_arr = subtract(data_arr, max);
			//System.out.println("Subtraction done");
			int i=end;
			int index=data_arr.length;
			//System.out.println(data_arr.length+" i= "+i+" index="+index);
			index = index -1;
			//System.out.println(data_arr.length+" i= "+i+" index="+index);
			while(i>start){
				//System.out.println(end+" i= "+index+" "+start);
				PolyASite site = findNextHill (index,data_arr,start,"upstream");
				//If the previous hill is sufficiently far from current hill
				if(!passNearnessTest(previousSite,site)&& site.getMagnitude()>MAGNITUDE_THRESHOLD){
					//Add the site to the vector of sites
					sites.add(site);
				}
				else{
					//Dont add the site
					//For now, we will consider the start of the first hill in the closely places hills as the site
				}
				previousSite = site;
				i = site.getPosition();
				index = i-start-1;
				//System.out.println("Number of sites:"+sites.size());
			}
		}
		//System.out.println("Exit");
		return sites;
	}
	

	/**
	 * This function will find the next hill, that is, next polyA site given a data array, direction to search in the array and a starting index
	 * @param start_arr
	 * @param data
	 * @param region_start
	 * @param direction
	 * @return
	 * @throws IOException
	 */
	private PolyASite findNextHill (int start_arr,double[] data,int region_start,String direction) throws IOException{
		
		//IF gene is on negative strand, and we are going from gene end to start which is downstream in genomic orientation
		if("downstream".equalsIgnoreCase(direction)){
			double prev = Double.NEGATIVE_INFINITY;
			double mag = 0.0;
			for(int i=start_arr;i<data.length;i++){
				if(data[i]<prev && data[i]>MAGNITUDE_THRESHOLD){
					//region_start is the gene_end
					return (new PolyASite(region_start+i-1,prev));
				}
				prev = data[i];
			}
			return (new PolyASite(region_start+data.length,prev));
		}
		//IF gene is on positive strand, and we are going from gene end to start which is upstream in genomic orientation
		else if("upstream".equalsIgnoreCase(direction)){
			double prev = Double.NEGATIVE_INFINITY;
			double mag = 0.0;
			
			//Here i is the index of data towards end of array
			for(int i=start_arr;i>=0;i--){
				if(data[i]<prev && data[i]>MAGNITUDE_THRESHOLD){
					//region_start is the gene start
					return (new PolyASite(region_start+i+1,prev));
				}
				prev = data[i];
			}
			//System.out.println("Nothing happened");
			return (new PolyASite(region_start,prev));
		}
		else{
			logger.error("Wrong direction supplied");
			return(null);
		}
			
	}
	
	/**
	 * This function will find the next hill, that is, next polyA site given a chormosome start position and the alignment model
	 * @param start Start of the window
	 * @param chr Chromosome
	 * @param model
	 * @param end
	 * @return
	 * @throws IOException
	 */
	private PolyASite findNextHill (int start,String chr,ContinuousDataAlignmentModel model, int end) throws IOException{
		
		double prev = 0.0;
		double mag = 0.0;
		
		while(start<end-WINDOW){
			mag = model.scoreGene(new RefSeqGene(new Alignments(chr,start,(start+WINDOW))))[2];
			//mag = model.getCount(new Alignments(chr,start,(start+WINDOW)));
			if(prev>mag){
				return (new PolyASite(start+WINDOW,mag)); 
			}
			prev = mag;
			start = start+WINDOW;
		}
		return (new PolyASite(end,mag));
	}
	
	/**
	 * This function converts an array of double to List of type Double
	 * @param double[]
	 * @return
	 */
	private List<Double> array2List(double[] array){
		List<Double> lt=new ArrayList<Double>();

		for(int i=0;i<array.length;i++){
			lt.add(array[i]);
		}
	
		return lt;
	}
	
	/**
	 * This functions returns an array after subtracting "factor" from each of the array elements.
	 * @copied from ContinuousDataAlignmentModel
	 * @param array
	 * @param factor
	 * @return
	 */
	private double[] subtract(double[] array, double factor){
		double[] rtrn=new double[array.length];
	
		for(int i=0; i<array.length; i++){
			//System.out.println(i);
			rtrn[i]=array[i]-factor;
		}
		//System.out.println("one");
		return rtrn;
	}

	/**
	 * This function returns true if the two polyA sites are closer than the NEARNESS_THRESHOLD
	 * @param site1
	 * @param site2
	 * @return
	 */
	private boolean passNearnessTest(PolyASite site1,PolyASite site2){
		
		if((Math.abs(site1.getPosition()-site2.getPosition()))<NEARNESS_THRESHOLD){
			return true;
		}
		else{
			return false;
		}
	}
	
	private class PolyASite{
		
		int position;
		double magnitude;
		int rank;
		
		PolyASite(int p,double m){
			position = p;
			magnitude = m;
		}
		
		PolyASite(){
			position = 0;
			magnitude = -1;
		}
		
		int getPosition(){
			return position;
		}
		
		double getMagnitude(){
			return magnitude;
		}
		 
		void setMagnitude(double m){
			magnitude =m;
		}
		
		void setPosition(int p){
			position=p;
		}

	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException{
		

		APA dummy = new APA(args);
	}

}
