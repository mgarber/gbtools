package broad.pda.rnai;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.pda.annotation.BEDFileParser;
import broad.pda.chromosome.Chromosome;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.rnai.designer.RNAiFileFormatUtils;

public class ProduceRNAiReport {
	
	int numNs=0;
	boolean repeatMask=true;
	int incrementNumber;
	private String species="10090";
	private String ncbiBuild="37";
	

	public ProduceRNAiReport(Collection<Alignments> lincs, Collection<RefSeqGene> transcripts, String save, String genomeDirectory, File okRepeatFile, boolean repeatMask, int incrementNumber)throws Exception{
		this.incrementNumber=incrementNumber;
		this.repeatMask=repeatMask;
		transcripts=setSequence(transcripts, genomeDirectory, okRepeatFile);
		
		//get names based on the lincRNA that overlaps it
		Map<Alignments, Collection<RefSeqGene>> overlaps=getOverlaps(lincs, transcripts);
		
		//output info and sequence for each version
		writeForRNAi(save, overlaps);
		System.err.println("Made "+total+" unique linc reports");
	}
	
	//Augment the current report
	public ProduceRNAiReport(Collection<Alignments> lincs, Collection<RefSeqGene> transcripts, String save, String genomeDirectory, File okRepeatFile, boolean repeatMask, Map<String, RNAiGeneAnnotation> rnaiReport, String format, boolean includePrevious)throws Exception{
		this.incrementNumber=getLastUsedNumber(rnaiReport);
		this.repeatMask=repeatMask;
		transcripts=setSequence(transcripts, genomeDirectory, okRepeatFile);
		
		//get names based on the lincRNA that overlaps it
		Map<Alignments, Collection<RefSeqGene>> overlaps=getOverlaps(lincs, transcripts);
		
		Collection<RNAiGeneAnnotation> lincRNAAnnotations=makeLincRNAAnnotations(overlaps, rnaiReport);
		if(includePrevious){
			lincRNAAnnotations.addAll(rnaiReport.values());
		}
		
		//output info and sequence for each version
		write(save, lincRNAAnnotations, format);
		System.err.println("Made "+total+" unique linc reports");
	}
	
	

	private int getLastUsedNumber(Map<String, RNAiGeneAnnotation> rnaiReport) {
		int rtrn=0;
		
		for(RNAiGeneAnnotation annotation: rnaiReport.values()){
			rtrn=Math.max(rtrn, Integer.parseInt(annotation.getGeneName().replaceAll("linc", "")));
		}
		
		return rtrn;
	}

	private Collection<RNAiGeneAnnotation> makeLincRNAAnnotations(Map<Alignments, Collection<RefSeqGene>> overlaps, Map<String, RNAiGeneAnnotation> rnaiReport) {
		Collection<RNAiGeneAnnotation> rtrn=new HashSet();
		
		Map<String, IntervalTree<Collection<RNAiGeneAnnotation>>> previousReport=makeTree(rnaiReport);
		
		//There are 2 forms of names
		//If an alignment overlaps a previous RNAiGeneAnnotation then its the same
		//Need to check if sequence matches, if so, then ignore, if not then update number and add
		
		for(Alignments align: overlaps.keySet()){
			String plusName="";
			String minusName="";
			int plusVersion=0;
			int minusVersion=0;
			//first check if novel or exists
			Iterator<Node<Collection<RNAiGeneAnnotation>>> overlappers=previousReport.get(align.getChr()).overlappers(align.getStart(), align.getEnd());
			boolean hasPrevious=overlappers.hasNext();
			while(overlappers.hasNext()){
				Collection<RNAiGeneAnnotation> annotations=overlappers.next().getValue();
				//check if its same strand
				for(RNAiGeneAnnotation annotation: annotations){
					String strand=annotation.geneSourceStrand;
					if(strand.equalsIgnoreCase("+")){plusName=annotation.getGeneName(); plusVersion=version(plusVersion, annotation.getTranscriptName());}
					if(strand.equalsIgnoreCase("-")){minusName=annotation.getGeneName(); minusVersion=version(minusVersion, annotation.getTranscriptName());}
				}
			}
			
			if(hasPrevious){
				Collection<RefSeqGene> genes=overlaps.get(align);
				//check if gene version is already contained
				for(RefSeqGene gene: genes){
					if(rnaiReport.containsKey(gene.getSequence().toUpperCase())){rtrn.add(rnaiReport.get(gene.getSequence().toUpperCase()));}
					else{
						//is different and therefore need to classify
						//make new report and increment version number
						String geneName=plusName;
						int versionNumber=plusVersion;
						if(gene.getOrientation().equalsIgnoreCase("-")){
							if(minusName.isEmpty()){minusName=increment();}
							geneName=minusName; versionNumber=(++minusVersion);
						}
						else{
							versionNumber=++plusVersion;
							if(plusName.isEmpty()){plusName=increment();}
							geneName=plusName;
						}
						RNAiGeneAnnotation annotation=makeRNAiGeneAnnotation(gene, align, geneName, versionNumber);
						rtrn.add(annotation);
					}
				}
			}
			else{
				//doesnt have a name so make one
				Collection<RefSeqGene> genes=overlaps.get(align);
				for(RefSeqGene gene: genes){
					String name="";
					int version=0;
					if(gene.getOrientation().equalsIgnoreCase("+")){
						if(plusName.isEmpty()){plusName=increment();}
						name=plusName;
						version=(++plusVersion);
					}
					else{
						if(minusName.isEmpty()){minusName=increment();}
						name=minusName;
						version=(++minusVersion);
					}
					//is different and therefore need to classify
					//make new report and increment version number
					RNAiGeneAnnotation annotation=makeRNAiGeneAnnotation(gene, align, name, version);
					rtrn.add(annotation);
					
				}
			}
		}
		
		
		return rtrn;
	}

	private Map<String, IntervalTree<Collection<RNAiGeneAnnotation>>> makeTree(Map<String, RNAiGeneAnnotation> rnaiReport) {
		Map<String, IntervalTree<Collection<RNAiGeneAnnotation>>> rtrn=new TreeMap();
		
		for(RNAiGeneAnnotation annotation: rnaiReport.values()){
			String chr=annotation.getRegion().getChr();
			IntervalTree<Collection<RNAiGeneAnnotation>> tree=new IntervalTree();
			if(rtrn.containsKey(chr)){tree=rtrn.get(chr);}
			//TODO add to tree
			Alignments region=annotation.getRegion();
			Node<Collection<RNAiGeneAnnotation>> node=tree.find(region.getStart(), region.getEnd());
			if(node!=null){
				Collection<RNAiGeneAnnotation> temp=node.getValue();
				temp.add(annotation);
				tree.put(region.getStart(), region.getEnd(), temp);
			}
			else{
				Collection<RNAiGeneAnnotation> temp=new ArrayList();
				temp.add(annotation);
				tree.put(region.getStart(), region.getEnd(), temp);
			}
			rtrn.put(chr, tree);
		}
		
		return rtrn;
	}

	private int version(int plusVersion, String transcriptName) {
		int num=new Integer(transcriptName.split("_")[1]);
		return Math.max(plusVersion, num);
	}

	//TODO Make a new transcript version annotation
	private RNAiGeneAnnotation makeRNAiGeneAnnotation(RefSeqGene gene, Alignments linc, String name, int version) {
		RNAiGeneAnnotation annotation=new RNAiGeneAnnotation(gene, linc, name, version, gene.getSequence(), species, ncbiBuild);
		// TODO Auto-generated method stub
		return annotation;
	}

	private String increment() {
		this.incrementNumber++;
		return "linc"+this.incrementNumber;
	}

	private void write(String save,	Collection<RNAiGeneAnnotation> lincRNAAnnotations, String format) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write(RNAiGeneAnnotation.nanostringHeader+"\n");
		
		for(RNAiGeneAnnotation linc: lincRNAAnnotations){
			/*if(format.equalsIgnoreCase("nanostring")){writer.write(linc.toNanostring());}
			else if(format.equalsIgnoreCase("RNAi")){writer.write(linc.toRNAi());}
			else{writer.write(linc.toDB());}*/
			writer.write(linc.toNanostring()+"\n");
		}
		
		writer.close();
	}

	

	private Map<Alignments, Collection<RefSeqGene>> getOverlaps(Collection<Alignments> lincs, Collection<RefSeqGene> transcripts){
		Map<Alignments, Collection<RefSeqGene>> rtrn=new TreeMap();
		
		Map<String, IntervalTree<RefSeqGene>> trees=makeIntervalTrees(transcripts);
		
		for(Alignments linc: lincs){
			Iterator<Node<RefSeqGene>> iter=trees.get(linc.getChr()).overlappers(linc.getStart(), linc.getEnd());
			rtrn.put(linc, toCollection(iter));
		}
		
		return rtrn;
	}
	
	private Map<String, IntervalTree<RefSeqGene>> makeIntervalTrees(Collection<RefSeqGene> alignments){
		Map<String, IntervalTree<RefSeqGene>> rtrn=new TreeMap();
		
		for(RefSeqGene align: alignments){
			IntervalTree<RefSeqGene> tree=new IntervalTree();
			String chr=align.getChr();
			if(rtrn.containsKey(chr)){tree=rtrn.get(chr);}
			tree.put(align.getAlignment().getStart(), align.getAlignment().getEnd(), align);
			rtrn.put(align.getChr(), tree);
		}
		
		return rtrn;
	}
	
	private Collection<RefSeqGene> toCollection(Iterator<Node<RefSeqGene>> iter){
		Collection rtrn=new TreeSet();
		while(iter.hasNext()){
			rtrn.add(iter.next().getValue());
		}
		return rtrn;
	}
	
	private Collection<RefSeqGene> setSequence(Collection<RefSeqGene> alignments, String genomeDirectory, File okRepeatFile)throws Exception{
		Map<String, IntervalTree<Alignments>> okRepeats=BEDFileParser.loadAlignmentDataToTree(okRepeatFile);
		Collection<RefSeqGene> rtrn=new TreeSet();
		Map<String, Set> alignmentsByChr=ExtractSequence.splitByChr(alignments);
		
		for(String chr: alignmentsByChr.keySet()){
			System.err.println(chr);
			Set<RefSeqGene> set=alignmentsByChr.get(chr);
			String sequenceFile=genomeDirectory+"/"+chr.replaceAll("chr", "").trim()+"/"+chr+".agp";
			Chromosome chrom = new Chromosome(sequenceFile);
			chrom.loadSequence();
			for(RefSeqGene align: set){
				String seq=ExtractSequence.getSequenceForGene(align, chrom, repeatMask, okRepeats);
				align.setSequence(seq);
				rtrn.add(align);
			}
		}
		
		return rtrn;
	}
	
	private void writeForRNAi(String save, Map<Alignments, Collection<RefSeqGene>> geneVersions)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		writer.write(getRNAiHeader()+"\n");
		
		int i=1;
		for(Alignments linc: geneVersions.keySet()){
			//System.err.println("index "+i+" "+linc.toUCSC());
			Collection<RefSeqGene> transcripts=geneVersions.get(linc);
			getRNAiString(linc, transcripts, i, writer);
			i++;
		}
		
		writer.close();
	}
	
	
	private String getRNAiHeader(){
		String rtrn="GENE.SOURCEID\tGENE.SOURCEVERSION\tGENE.SYMBOL\tGENE.ENTREZGENEID\tGENE.TAXONID\tGENE.CHROMOSOME\tGENE.MAPLOCATION\tGENE.NCBI_BUILDID\tGENE.SOURCECONTIG\tGENE.SOURCESTART\tGENE.SOURCEEND\tGENE.SOURCESTRAND\tGENE.CREATEDATE\tGENE.UPDATEDATE\tTRANS.SOURCEID\tTRANS.SOURCEVERSION\tTRANS.SEQ\tTRANS.LENGTH\tTRANS.SEQCREATEDATE\tTRANS.REFSEQ\tTRANS.GENBANKID\tTotal.length\tUnique.length";
		return rtrn;
	}
	
	
	private void getRNAiString(Alignments linc, Collection<RefSeqGene> transcripts, int index, FileWriter writer)throws IOException{
		String rtrn="";
		
		index=(index+this.incrementNumber);
		
		String sourceID="linc"+index;
		String sourceVersion="1";
		String symbol="linc"+index;
		String entrez="NA";
		String taxonID="10090"; //Hardcoded for mouse
		String build="37";
		
		boolean used=false;
		
		if(linc.getName()!=null && !linc.getName().isEmpty()){sourceID=linc.getName(); symbol=linc.getName();}
		
		String general=(sourceID+"\t"+sourceVersion+"\t"+symbol+"\t"+entrez+"\t"+taxonID);
		
		int subindex=1;
		for(RefSeqGene transcript: transcripts){
			String add="\t"+transcript.getChr()+"\tNA\t"+build+"\t"+linc.getChr()+"\t"+linc.getStart()+"\t"+linc.getEnd()+"\t"+transcript.getOrientation()+"\t"+getDateTime()+"\tNA\t"+sourceID+"_"+subindex+"\t"+"1"+"\t"+transcript.getSequence()+"\t"+transcript.getSequence().toCharArray().length+"\t"+getDateTime()+"\tNA\tNA";
			writer.write(general+add);
			int[] counts=getCounts(transcript.getSequence());
			writer.write("\t"+counts[0]+"\t"+counts[1]+"\n");
			subindex++;
			used=true;
		}
		if(used){this.total++;}
	}
	
	private int[] getCounts(String seq){
		int total=0;
		int upper=0;
		
		char[] chars=seq.toCharArray();
		
		for(int i=0; i<chars.length; i++){
			if(chars[i]=='A' | chars[i]=='C' || chars[i]=='G' || chars[i]=='T'){upper++;}
			total++;
		}
		
		int[] array={total, upper};
		return array;
	}
	
	 private String getDateTime() {
	        DateFormat dateFormat = new SimpleDateFormat("MM/dd/yyyy");
	        Date date = new Date();
	        return dateFormat.format(date);
	    }
	
	 
	 static String usage=" args[0]=linc List (K4-K36 bed) \n args[1]=transcripts (full BED) \n args[2]=save file \n args[3]=genomeDirectory \n args[4]=ok repeats \n args[5]=repeatMask \n args[6]=RNAi Report";
		int total=0;
		
		public static void main(String[] args)throws Exception{
			if(args.length>7){
				Collection<Alignments> lincFile=BEDFileParser.loadAlignmentData((new File(args[0])));
				Collection<RefSeqGene> transcriptFile=BEDFileParser.loadData(new File(args[1]));
				String save=args[2];
				String genomeDirectory=args[3];
				File okRepeats=new File(args[4]);
				boolean repeatMask=new Boolean(args[5]);
				Map<String, RNAiGeneAnnotation> report=RNAiFileFormatUtils.parseRNAiReportFileBySequence(new File(args[6]));
				boolean includePrevious=new Boolean(args[7]);
				new ProduceRNAiReport(lincFile, transcriptFile, save, genomeDirectory, okRepeats, repeatMask, report, "RNAi", includePrevious);
			}
			else{System.err.println(usage);}
		}
	
}
