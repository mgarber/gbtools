package broad.pda.gene;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.broad.igv.Globals;

import broad.core.datastructures.IntervalTree;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.annotation.GTFFileParser;

public class IsoformComparer {

	static String usage="Usage: IsoformComparer -task <task name> "+
						"\n-compare";
	
	private double MIN_OVERLAP = 0.4;
	
	public IsoformComparer(BEDFileParser annotations,ArgumentMap argmap) throws IOException{
		
		BufferedWriter bw =argmap.getOutputWriter();
		Map<String, IntervalTree<RefSeqGeneWithIsoforms>> mergedGenes = annotations.getMergedAnnotationMap(MIN_OVERLAP);
		Iterator<String> chrIt=mergedGenes.keySet().iterator();

		while(chrIt.hasNext()){
			String chr=chrIt.next();
			IntervalTree<RefSeqGeneWithIsoforms> tree=mergedGenes.get(chr);
			Iterator<RefSeqGeneWithIsoforms> geneIt=tree.valueIterator();
			while(geneIt.hasNext()) {
				RefSeqGeneWithIsoforms gene = geneIt.next();
				List<Integer> starts = new ArrayList<Integer>();
				List<Integer> ends = new ArrayList<Integer>();
				boolean sameStart = true;
				boolean sameEnd = true;
				
				gene.setName("gene");
				gene.cleanIsoforms();
				int isoformsAdded = 0;
				Iterator<RefSeqGeneWithIsoforms> overlappers = annotations.getOverlappers(gene).valueIterator();
				while(overlappers.hasNext()) {
					RefSeqGeneWithIsoforms overlapper = overlappers.next();
					if(BEDFileParser.isOverlapCompatible(gene, overlapper, MIN_OVERLAP)) {
						//System.err.println("Adding overlapper " + overlapper.toBED());
						Collection<RefSeqGene> isoforms = overlapper.getAllIsoforms();
						if(isoforms.size()>1){
							
							for(RefSeqGene iso: isoforms) {
								gene.setName(gene.getName()+"_"+iso.getName());
								boolean couldAdd = gene.addContainedIsoform(iso);
								if(!couldAdd) {
									//System.err.println("WARN: Could not add isoform " + overlapper.toBED() + " to " + gene.toBED());
								} else {
									isoformsAdded++;
									//System.err.println("Added isoform " + overlapper.getName() + " to " + gene.getName());
								}
								//Check if the starts and ends are the same
								starts.add(iso.getOrientedStart());
								ends.add(iso.getOrientedEnd());
							}
						}
					}
				}
				
				if(isoformsAdded ==0) {
					System.err.println("ERROR: Gene " + gene.getName() + " " + gene.toBED() + "  had no overlapping isoforms");
				}
				else{
					int s = starts.get(0);
					int e = ends.get(0);
					for(int i=1;i<starts.size();i++){
						if(starts.get(i)!=s){
							sameStart=false;
						}
						if(ends.get(i)!=e)
							sameEnd = false;
					}
					bw.write(gene.getName()+"\t"+isoformsAdded+"\t"+sameStart+"\t"+sameEnd+"\n");

				}
			}
		}		
		
	}
	
	public static void main(String [] args) throws IOException {
		
		/*
		 * Gives a log4j error. Check later.
		 */
		Globals.setHeadless(true);
		ArgumentMap argmap = CLUtil.getParameters(args, usage, "compare");
		
		String annotationFile = argmap.getMandatory("annotations");
		BEDFileParser annotationParser =  annotationFile.endsWith(".gtf") || annotationFile.endsWith(".GTF")? new GTFFileParser(annotationFile) : new BEDFileParser(annotationFile);
		
		IsoformComparer IC = new IsoformComparer(annotationParser,argmap);
	}
}
