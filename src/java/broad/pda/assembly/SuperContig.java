package broad.pda.assembly;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;


public class SuperContig extends Sequence{
	String name;
	int size;
	List<AgpEntry> gaps;
	List<AgpEntry> contigs;
	
	public SuperContig(String name) {
		super(name);
		this.name = name;
		gaps = new ArrayList<AgpEntry>();
		contigs = new ArrayList<AgpEntry>();
	}
	
	public void addAGPEntry(AgpEntry entry) {
		//System.out.println("Adding entry: " + entry.getContainingSequenceId() + "\t" + entry.getName());
		if(entry.getType() == AgpEntry.GAP_TYPE) {
			addGap(entry);
		} else {
			addContig(entry);
		}
	}
	
	
	public void addContig(AgpEntry contig) {
		contigs.add(contig);
		size += contig.getLength();
	}
	
	public void addGap(AgpEntry gap) {
		gaps.add(gap);
		size += gap.getLength();
	}


	public List<AgpEntry> getContigs() {
		return contigs;
	}

	public List<AgpEntry> getGaps() {
		return gaps;
	}

	public String getName() {
		return name;
	}

	public int getSize() {
		return size;
	}

	public void loadContigsFromContigFile(File sequenceFile) throws Exception {
		FastaSequenceIO fsio = new FastaSequenceIO(sequenceFile);
		fsio.extractRecordsIntoSequenceList(contigs);
	}

	/**
	 * Assumes that all contigs' sequence has been loaded, methods assumes contig and gap lists
	 * are ordered.
	 * @see loadContigsFromContigFile
	 * @param output file to append super contig.
	 */
	public String getSequenceBases() {
		String sequenceBases = super.getSequenceBases();
		if(sequenceBases == null || sequenceBases.length() == 0) {
		    //System.out.println("Processing supercontig " + getName() + " with " + gaps.size() + " gaps and " + contigs.size() + " contigs");
			AgpEntry firstGap = gaps.size() > 0 ? gaps.get(0) : null;
			AgpEntry firstContig = contigs.get(0);
		
			List<AgpEntry> first = (firstGap != null) && firstGap.getStart() < firstContig.getStart() ? gaps : contigs;
			List<AgpEntry> second = (firstGap != null) && firstGap.getStart() < firstContig.getStart() ? contigs : gaps;

			int i = 0;
			for(i = 0; i < Math.min(first.size(), second.size()); i++) {
				if(first.get(i).getSequenceBases() == null) {
					throw new RuntimeException("contig " + first.get(i) + " has no sequence");
				} else if (second.get(i).getSequenceBases() == null) {
					throw new RuntimeException("contig " + second.get(i) + " has no sequence");
				}
				appendToSequence(first.get(i).getSequenceBases());
				appendToSequence(second.get(i).getSequenceBases());
			}
		
			if(first.size() > second.size()) {
				appendToSequence(first.get(i).getSequenceBases());
			}
			sequenceBases = super.getSequenceBases();
		}
		return sequenceBases;
	}

	public void unloadSequence() {
		Iterator<AgpEntry> contigIt = contigs.iterator();
		while(contigIt.hasNext()) {
			contigIt.next().unloadSequence();
		}
		super.unloadSequence();
	}
	
	public void loadContigsFromContigFileByName(File sequenceFile, final String agpContigPrefix) throws IOException {
		FastaSequenceIO fsio = new FastaSequenceIO(sequenceFile);
		List<String> regExList = new ArrayList<String>(contigs.size());
		Iterator<AgpEntry> contigIt = contigs.iterator();
		while(contigIt.hasNext()){
			AgpEntry contig = contigIt.next();
			regExList.add("_" + contig.getName().replace(agpContigPrefix, "") + " ");
			//System.out.println("Adding " + "_" + contig.getName().replace(agpContigPrefix, ""));
		}		
		addSequenceToContigs(fsio, regExList, new Comparator< Sequence>() {

			public int compare(Sequence arg0,Sequence arg1){
				
				return arg0.getId().contains("_" + arg1.getId().replace(agpContigPrefix, "") + " ") ? 0 : 1;
			}
			
		});
	}

	public void loadContigsFromContigFileBySize(File sequenceFile) throws Exception {
		FastaSequenceIO fsio = new FastaSequenceIO(sequenceFile);
		List<String> regExList = new ArrayList<String>(contigs.size());
		Iterator<AgpEntry> contigIt = contigs.iterator();
		while(contigIt.hasNext()){
			AgpEntry contig = contigIt.next();
			regExList.add(contig.getLength() + " nt");
		}
		addSequenceToContigs(fsio, regExList, new Comparator<Sequence> (){

			public int compare(Sequence arg0, Sequence arg1) {
				return arg0.getLength() - arg1.getLength();
			}
			
		});
	}

	private void addSequenceToContigs(FastaSequenceIO fsio, List<String> regExList, Comparator<Sequence> comparator) throws IOException {
		Iterator<AgpEntry> contigIt;
		List<Sequence> contigSeqs = fsio.extractRecordsWithIDsMatching(regExList, false);

		contigIt = contigs.iterator();
		while(contigIt.hasNext() ) {
			Iterator<Sequence> contigSeqIt = contigSeqs.iterator();
			AgpEntry contig = contigIt.next();
			boolean found = false;
			System.out.println("Trying to find sequence for " + contig.getId() + " is in reversed orientation " + contig.inReversedOrientation());
			while(contigSeqIt.hasNext() && !found) {
				Sequence contigSeq = contigSeqIt.next();

				if(comparator.compare(contigSeq, contig) ==  0) {
					if(contig.inReversedOrientation()) {
						System.out.println("contig " + contig.getName() + " is inverted ");
						contigSeq.reverse();
					}
					contig.setSequenceBases(contigSeq.getSequenceBases());
					contigSeq.unloadSequence();
					System.out.println("\tFOUND!Adding sequence to " + contig.getName() + " size " + contig.getLength() + " name " + contig.getName());
					found = true;
				}
			}
			if(!found) {
				throw new RuntimeException("\tCould not find " + contig.getId() + " length " + contig.getLength() );
			}
		}
	}

	public void addInitialGap() {
	    //System.out.println("Processing supercontig " + getName() + " with " + gaps.size() + " gaps and " + contigs.size() + " contigs");
		AgpEntry firstGap = gaps.size() > 0 ? gaps.get(0) : null;
		AgpEntry firstContig = contigs.get(0);
	
		List<AgpEntry> first = (firstGap != null) && firstGap.getStart() < firstContig.getStart() ? gaps : contigs;
		
		String chr = first.get(0).getChromosome() ;
		if(first.get(0).getStart() > 1) {
			System.out.println("chr" + chr + " agp does not start at 1, adding initial gap");
			Gap initialGap = new Gap(chr);
			initialGap.setStart(1);
			initialGap.setEnd(first.get(0).getStart() - 1);
			initialGap.setChromosome(chr);
			
			List<AgpEntry> tmpGaps = new ArrayList<AgpEntry>(gaps.size() + 1);
			tmpGaps.add(initialGap);
			Iterator<AgpEntry> gapIt = gaps.iterator();
			while(gapIt.hasNext()) {
				tmpGaps.add(gapIt.next());
			}
			gaps = tmpGaps;
		}
	
	}


}
