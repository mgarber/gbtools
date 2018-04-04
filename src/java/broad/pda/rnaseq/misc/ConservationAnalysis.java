package broad.pda.rnaseq.misc;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.math.Statistics;
import broad.core.siphy.EvolutionaryModel.PiFit;
import broad.core.siphy.StationaryDistributionIO;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.assembly.ChunkResolver;
import broad.pda.assembly.ChunkResolver.Chunk;
import broad.pda.chromosome.Chromosome;
import broad.pda.chromosome.GenericOrganism;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.RefSeqGene;
import broad.pda.gene.RefSeqGeneWithIsoforms;

public class ConservationAnalysis {
	
	static final String USAGE = "\n\tsplice: Get PI scores around splice junctions \n\t-genes <Gene set> \n\t-assemblyDir <Top level directory of assembly>\n\t-pidir <PI Call directory> \n\t-radius <Radious around the splice site to compute extract pis> \n\t-extract <[LOD (log odds ratio score)],BL (branch length)> \n\t-piChunkSize <Pi call chunk size, default is 1000000>"+
	"\n\texons: Get PI average stats for exons \n\t-genes <Gene set> \n\t-assemblyDir <Top level directory of assembly>\n\t-pidir <PI Call directory>  \n\t-piChunkSize <Pi call chunk size, default is 1000000>"+
	"\n\tCrossWithSiPhyOutput: extract from a Siphy output file only rows that overlap in the input bed \n\t -bed <BED> \n\t -siphy <siphy file> \n\t -best3 <optional, for each region in the bed file extract only the 3 exons that have the lowest omega score and overlap the ref seq> \n\t -exactExons <flag to include only exons that exactly appear in the bed> -out <siphy output file>"+
	"\n\textractConservedTranscripts: extract from bed only the transcripts that overlap a conserved region in the Siphy output file \n\t -bed <BED> \n\t -siphy <siphy file> \n\t -omega <upper threshold of omega value> \n\t\t -bl <minimal branch length> \n\\t\t -out <siphy output file>"+
	"\n";
	
	public static void main(String [] args) throws Exception  {
		ArgumentMap argmap = CLUtil.getParameters(args, USAGE, "splice");

		if("splice".equalsIgnoreCase(argmap.getTask())) {
			String geneFile = argmap.getMandatory("genes");
			String piDir    = argmap.getMandatory("pidir");
			int    radius   = argmap.getInteger("radius");
			String assemblyDir = argmap.getMandatory("assemblyDir");
			int    chunkSize= argmap.containsKey("piChunkSize") ? argmap.getInteger("piChunkSize") : 1000000;
			boolean useLOD  = argmap.containsKey("use") && !"BL".equalsIgnoreCase(argmap.get("use"));
			Map<String, Collection<RefSeqGene>> regionMap =BEDFileParser.loadDataByChr(new File(geneFile));
			
			GenericOrganism go = new GenericOrganism(new File(assemblyDir));
			ChunkResolver cr = new ChunkResolver(new File(piDir), chunkSize, ".pi" );
			StationaryDistributionIO currentChunkData = null;
			Chunk currentChunk = null;
			Collection<PiFit []> data = new ArrayList<PiFit []>();
			for(String chr : regionMap.keySet()) {
				Chromosome chrObj = go.getChromosome(chr.replace("chr", ""));
				Collection<RefSeqGene> chrGenes = regionMap.get(chr);
				for(RefSeqGene g : chrGenes) {
					Alignments [] gExons = g.getExons();
					for(int i = 0; i < gExons.length; i++) {
						Alignments exon = gExons[i];
						if(currentChunk == null || !currentChunk.contains(exon)) {
							System.err.println("Loading PI data for exon " + exon.toUCSC() + "  because " + (currentChunk == null ? "null": currentChunk.toUCSC() + "  contains exon " + exon.toUCSC() + " " + currentChunk.contains(exon) ));
							currentChunk = cr.getChunk(chr, exon.getStart(), chrObj.getSize());
							currentChunkData = new StationaryDistributionIO();
							currentChunkData.load(currentChunk.getLocation().getAbsolutePath());
						}
						if(exon.length() > radius) {
							if(i>0) {
								LightweightGenomicAnnotation prime5 = new BasicLightweightAnnotation(chr, exon.getStart() - radius, exon.getStart());
								PiFit [] prime5Data = getPis(currentChunkData, prime5);
								if(prime5Data != null) {
									data.add(prime5Data);
								} else {
									System.err.println("\tNo 5' data pi for " + exon.toUCSC());
								}
							}
							if(i< gExons.length - 1) {
								LightweightGenomicAnnotation prime3 = new BasicLightweightAnnotation(chr, exon.getEnd(), exon.getEnd() + radius);
								PiFit [] prime3Data = getPis(currentChunkData, prime3);
								if(prime3Data != null) {
									data.add(prime3Data);
								}else {
									System.err.println("\tNo 3' data pi for " + exon.toUCSC());
								}
							}
						}	
					}
				}	
			}
			
			BufferedWriter bw = argmap.getOutputWriter();
			writePiData(useLOD, data, bw);
			bw.close();
		} else if("exons".equalsIgnoreCase(argmap.getTask())) {
			String geneFile = argmap.getMandatory("genes");
			String piDir    = argmap.getMandatory("pidir");
			String assemblyDir = argmap.getMandatory("assemblyDir");
			int    chunkSize= argmap.containsKey("piChunkSize") ? argmap.getInteger("piChunkSize") : 1000000;
			Map<String, Collection<RefSeqGene>> regionMap =BEDFileParser.loadDataByChr(new File(geneFile));
			
			GenericOrganism go = new GenericOrganism(new File(assemblyDir));
			ChunkResolver cr = new ChunkResolver(new File(piDir), chunkSize, ".pi" );
			StationaryDistributionIO currentChunkData = null;
			Chunk currentChunk = null;
			Map<LightweightGenomicAnnotation, double []> data = new LinkedHashMap<LightweightGenomicAnnotation, double []>();
			for(String chr : regionMap.keySet()) {
				Chromosome chrObj = go.getChromosome(chr.replace("chr", ""));
				Collection<RefSeqGene> chrGenes = regionMap.get(chr);
				for(RefSeqGene g : chrGenes) {
					Alignments [] gExons = g.getExons();
					for(int i = 0; i < gExons.length; i++) {
						Alignments exon = gExons[i];
						if(currentChunk == null || !currentChunk.contains(exon)) {
							System.err.println("Loading PI data for exon " + exon.toUCSC() + "  because " + (currentChunk == null ? "null": currentChunk.toUCSC() + "  contains exon " + exon.toUCSC() + " " + currentChunk.contains(exon) ));
							currentChunk = cr.getChunk(chr, exon.getStart(), chrObj.getSize());
							currentChunkData = new StationaryDistributionIO();
							currentChunkData.load(currentChunk.getLocation().getAbsolutePath());
						}
						PiFit [] exonFits = getPis(currentChunkData, exon);
						if(exonFits != null) {
							data.put(exon, summaryze(exonFits));
						}
					}
				}	
			}
			
			BufferedWriter bw = argmap.getOutputWriter();
			writeSummaryPiData( data, bw);
			bw.close();
		}
		else if("CrossWithSiPhyOutput".equalsIgnoreCase(argmap.getTask())){
			String bedFile = argmap.getMandatory("bed");
			String siphyFile    = argmap.getMandatory("siphy");
			boolean  best3 = argmap.containsKey("best3")? true : false ;
			BufferedWriter bw = argmap.getOutputWriter();
			boolean exactExons = argmap.containsKey("exactExons")? true: false;
			if (exactExons)
				CrossWithSiPhyOutputByExactExons(bedFile,siphyFile,best3,bw);
			else
				CrossWithSiPhyOutput(bedFile,siphyFile,best3,bw);
			bw.close();
		}
		else if("extractConservedTranscripts".equalsIgnoreCase(argmap.getTask())){
			String bedFile = argmap.getMandatory("bed");
			String siphyFile    = argmap.getMandatory("siphy");
			double  omega = argmap.getDouble("omega") ;
			double  bl= argmap.getDouble("bl") ;
			BufferedWriter bw = argmap.getOutputWriter();
			extractConservedTranscripts(bedFile,siphyFile,omega,bl,bw);
			bw.close();
		}
		else{System.err.println(USAGE);}
	}

	

	



	private static double[] summaryze(PiFit[] fits) {
		List<Double> lods = new ArrayList<Double>(fits.length);
		List<Double> bls  = new ArrayList<Double>(fits.length);
		
		for(PiFit fit : fits) {
			lods.add(fit == null ? 0d : fit.getLogLikelihoodRatio());
			bls.add(fit == null ? 0d : fit.getTreeLength());
		}
		
		double [] summary = {Statistics.mean(lods), Statistics.mean(bls)};
		return summary;
	}

	private static void writeSummaryPiData(Map<LightweightGenomicAnnotation, double []> data, BufferedWriter bw) throws IOException {
		for(LightweightGenomicAnnotation region : data.keySet()) {
			double [] regionData = data.get(region);
			bw.write(region.getName()+"\t" + region.toUCSC());
			for(double datapoint : regionData){
				bw.write("\t" + datapoint);
			}
			bw.newLine();
		}
		
	}

	private static void writePiData(boolean useLOD, Collection<PiFit[]> data,BufferedWriter bw) throws IOException {
		for(PiFit [] pis : data) {
			for(int i = 0; i < pis.length; i++) {
				if(pis[i] == null) {
					bw.write("0");
				} else {
					bw.write(String.valueOf(useLOD ? pis[i].getLogLikelihoodRatio() : pis[i].getTreeLength()));
				}
				if(i < pis.length - 1) {
					bw.write("\t");
				}
			}	
			bw.newLine();
		}
	}

	private static PiFit[] getPis(StationaryDistributionIO currentChunkData, LightweightGenomicAnnotation region) {
		Iterator<PiFit> piIt = currentChunkData.getOverlappers(region);
		PiFit [] data = new PiFit[region.length()];
		
		int i = 0;
		//System.err.println("PIs for region " + region.toUCSC());
		while(piIt.hasNext() && i < region.length()) {
			PiFit fit = piIt.next();
			while(fit.getPosition() > region.getStart() + i) {
				data[i++] = null;
			}
			//System.err.println("\t" + fit.getPosition() + " - " + fit.getLogLikelihoodRatio() + " - " + fit.getTreeLength());
			data[i++] = fit;
		}
		return i<region.length() ?  null : data;
	}

	
	
	private static void CrossWithSiPhyOutput(String bedFile, String siphyFile,boolean best3, BufferedWriter bw) throws IOException {
		
		BEDFileParser bed =new BEDFileParser(bedFile);
		BEDFileParser siphy=BEDFileParser.readSiphyOutToBed (siphyFile);
		double totalGene=0.0;
		double totalAdded=0.0;
		Set<RefSeqGeneWithIsoforms> outSet=new HashSet<RefSeqGeneWithIsoforms>(); //Set such that every element will be included only once
		for (RefSeqGene g :bed.GetGenes()){
			IntervalTree<RefSeqGeneWithIsoforms> overlaps=siphy.getOverlappers(g);
			//Add the gene Name as a last column to the siphy ; June 2012
			for (RefSeqGeneWithIsoforms currIso: overlaps.toCollection())
				currIso.setName(g.getName());
			outSet.addAll(overlaps.toCollection());
			totalGene++;
			if (! overlaps.isEmpty())
				totalAdded++;
		}
		for(RefSeqGeneWithIsoforms g: outSet){
			for(RefSeqGene iso:g.getAllIsoforms()){
				bw.write(iso.getChr()+"\t"+iso.getStart()+"\t"+iso.getEnd()+"\t"+iso.getExtraFields(0)+"\t"+iso.getExtraFields(1)+"\t"+iso.getExtraFields(2)+"\t"+iso.getExtraFields(3)+"\t"+iso.getName()+"\n");
			}
		}
		System.err.println("Total genes is: "+totalGene +" total added: "+totalAdded+ "precentage "+ totalAdded/totalGene);
	}
	
   private static void CrossWithSiPhyOutputByExactExons(String bedFile, String siphyFile,boolean best3, BufferedWriter bw) throws IOException {
		
		BEDFileParser bed =new BEDFileParser(bedFile);
		BEDFileParser siphy=BEDFileParser.readSiphyOutToBed (siphyFile);
		double totalGene=0.0;
		double totalAdded=0.0;
		Set<RefSeqGeneWithIsoforms> outSet=new HashSet<RefSeqGeneWithIsoforms>(); //Set such that every element will be included only once
		for (RefSeqGene g :bed.GetGenes()){
			IntervalTree<RefSeqGeneWithIsoforms> overlaps=siphy.getOverlappers(g);
			Set<Alignments> exons=g.getExonSet();
			Set<RefSeqGeneWithIsoforms> exactOverlaps= new  TreeSet<RefSeqGeneWithIsoforms>();
			Iterator <RefSeqGeneWithIsoforms> it = overlaps.valueIterator();
			while(it.hasNext()){
				RefSeqGeneWithIsoforms currEx=it.next();
				for (Alignments ex:exons){
					if (ex.getStart()==currEx.getStart() & ex.getEnd()==currEx.getEnd() )
					{ exactOverlaps.add(currEx);	break;}		
				}
			}
			
			for (RefSeqGeneWithIsoforms currIso: exactOverlaps)
				currIso.setName(g.getName());
			outSet.addAll(exactOverlaps);
			totalGene++;
			if (! exactOverlaps.isEmpty())
				totalAdded++;
		}
		for(RefSeqGeneWithIsoforms g: outSet){
			for(RefSeqGene iso:g.getAllIsoforms()){
				bw.write(iso.getChr()+"\t"+iso.getStart()+"\t"+iso.getEnd()+"\t"+iso.getExtraFields(0)+"\t"+iso.getExtraFields(1)+"\t"+iso.getExtraFields(2)+"\t"+iso.getExtraFields(3)+"\t"+iso.getName()+"\n");
			}
		}
		System.err.println("Total genes is: "+totalGene +" total added: "+totalAdded+ "precentage "+ totalAdded/totalGene);
	}
	
	
	
	private static void extractConservedTranscripts(String bedFile,String siphyFile, double omega,double bl, BufferedWriter bw) throws IOException {
		
		BEDFileParser bed =new BEDFileParser(bedFile);
		BEDFileParser siphy=BEDFileParser.readSiphyOutToBed (siphyFile);
		BEDFileParser outSet=new BEDFileParser(); //Set such that every element will be included only once
		for (RefSeqGene g :bed.GetGenes()){
			g.setBedScore(10);
			Iterator<RefSeqGeneWithIsoforms> it=siphy.getOverlappers(g).valueIterator();
			boolean flag=false;
			while(it.hasNext()){
				RefSeqGene overlapper=it.next();
				double curr_scr=Double.valueOf(overlapper.getExtraFields(0));
				if (g.overlaps(overlapper) && curr_scr <= omega && Double.valueOf(overlapper.getExtraFields(3)) >= bl){
					flag=true;
					if(g.getBedScore()> curr_scr)
						g.setBedScore(curr_scr);
				}
			}
			if (flag)
				outSet.addRefSeq(g);
		}
		outSet.writeFullBed(bw);
		
	}
	
}
