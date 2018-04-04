package broad.pda.assembly;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import broad.core.alignment.RepeatMaskerReader;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

public class Assembly {
	public static String USAGE = "Usage: Assembly TASK=<task_num> <task_args>\n" +
	"\tTasks:\n" +
	"\t\t1. Convert contig assembly to supercontig: OUT=<super contig assembly file name> IN=agp file SEQFILE=<contig multifasta>\n" +
	"\t\t2. Install broad assembly in current directory IN=<multifasta file with all chromosomes> AGP=<One file agp> [-unchunked <if the chromosomes are unchunked in the multifasta file>]\n" +
	"\t\t3. Install only Un from Broad assembly in current dir -in <assembly>\n" +
	"\t\t4. Install contig based assembly -in <agp file> -seqfile <contig based fasta file> \n"+
	"\t\t\t[-usesize <in case the agp does not have unique contig names> OR -agpContigNamePrefix <prefix of contig names in the agp file]\n";
	
	private final static int NUM_CONTIG_TO_GROUP = 5000;
	
	String organismName;
	String version;
	File source;
	int size;
	File sequenceFile;
	ArrayList<SuperContig> supers;
	RepeatMaskerReader repeatReader;
	

	public Assembly(String agpFile) throws IOException {
		super();
		size = 0;
		File source = new File(agpFile);
		BufferedReader br = new BufferedReader(new FileReader(source));
		String agpWithoutExt = agpFile.substring(0,agpFile.lastIndexOf("."));
		sequenceFile = new File(agpWithoutExt + ".fa");
		String line;
		supers = new ArrayList<SuperContig>();
		repeatReader = new RepeatMaskerReader();
		
		try {
			String currentSuper = "";
			while((line = br.readLine()) != null) {
				if(line.startsWith("#")){
					continue;
				}
				String[] lineSplit = line.split("\t");
				AgpEntry entry = AgpEntryFactory.createEntry(lineSplit);
				//System.out.println("Entry chr " + entry.getChromosome());
				
				if(!currentSuper.equals(entry.getChromosome())) {
					if(supers.size() > 0) {
						size += supers.get(supers.size() - 1).getSize();
					}
					SuperContig superContig = new SuperContig(entry.getChromosome());					
					supers.add(superContig);
					currentSuper = superContig.getName();
					//System.out.println("Containing seq " + entry.getChromosome() + " entry " + entry.getName());
				}
				
				supers.get(supers.size() - 1).addAGPEntry(entry);
			}
			Iterator<SuperContig> it = supers.iterator();
			while(it.hasNext()) {
				it.next().addInitialGap();
			}
		}  finally {
			try {
				System.out.print("Closing "+agpFile);
				br.close();
				System.out.print(" ..... Closed\n");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

	public void setSequenceFile(String pathToSequenceFile) {
		this.sequenceFile = new File(pathToSequenceFile);
	}
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws Exception {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE);
		if ("1".equals(argMap.getTask())) {
			String output = argMap.getOutput();
			File outFile = new File(output);
			if (outFile.exists()) {
				outFile.delete();
			}
			Assembly assembly = new Assembly(argMap.getInput());
			assembly.setSequenceFile(argMap.getMandatory("SEQFILE"));
			assembly.createSuperContigFasta(output);
		}else if ("2".equals(argMap.getTask())) {
			String fastaFile = argMap.getInput();
			String agpFile   = argMap.getMandatory("AGP");
			boolean unchunked = argMap.containsKey("unchunked");
			String curChr = null;
			BufferedReader agpBR = new BufferedReader(new FileReader(agpFile));
			String line = null;
			String lineChr = null;
			BufferedWriter agpBW = null;
			while((line = agpBR.readLine()) != null) {
				lineChr = line.substring(3, line.indexOf("\t") );
				if(!lineChr.equals(curChr)) {
					if(agpBW != null) {
						System.out.println("Closing agp file for chr " + curChr);
						agpBW.close();						
					}
					File chrDir = new File(lineChr);


					if(!lineChr.startsWith("Un")){
						if(!chrDir.exists()) {
							chrDir.mkdir();
						}
						agpBW = new BufferedWriter(new FileWriter(chrDir + "/chr" + chrDir + ".agp"));
						FastaSequenceIO fsio = new FastaSequenceIO(fastaFile);
						Sequence wholeChr = new Sequence("chr" + lineChr);
						if(unchunked) {
							System.out.println("unchunked! looking for chr" + lineChr);
							ArrayList<String> onelist = new ArrayList<String>(1);
							onelist.add("chr" + lineChr);
							fsio.write(fsio.extractRecords(onelist), chrDir + "/chr" + chrDir + ".fa");
						} else {
							fsio.mergeRecords(wholeChr, "^"+lineChr + "\\.");
							fsio.write(wholeChr, chrDir + "/chr" + chrDir + ".fa");
						}
						
						wholeChr.unloadSequence();
						agpBW.write(line);
						agpBW.newLine();
					}
					curChr = lineChr;
				}

			}
			agpBW.close();
			agpBR.close();
			
			installUn(fastaFile);
		}else if ("3".equals(argMap.getTask())) {
			String fastaFile = argMap.getInput();
			installUn(fastaFile);

		}else if ("4".equals(argMap.getTask())) {
			String agp = argMap.getInput();
			String seqFile = argMap.getMandatory("seqfile");
			Assembly assembly = new Assembly(agp);
			assembly.setSequenceFile(seqFile);
			if(argMap.containsKey("usesize")) {
				assembly.createChromosomeBasedFastaByContigSize();
			} else if(argMap.isPresent("agpContigNamePrefix")){
				assembly.createChromosomeBasedFasta(argMap.getMandatory("agpContigNamePrefix"));
			}else {
				System.err.println("Must specify either usesize or agpContigNamePrefix. "+ USAGE);
			}
			

		}else {
			System.err.println(USAGE);
		}
	}

	private static void installUn(String fastaFile) throws IOException {
		FastaSequenceIO fsio = new  FastaSequenceIO(fastaFile);
		List<Sequence> unList = fsio.extractRecordsWithIDLike("^chrUn", false);
		
		fsio.write(unList, "Un/chrUn.fa");
		/*
		Iterator<Sequence> unSequenceChunkIt = unList.iterator();	
		if(!unSequenceChunkIt.hasNext()) {
			return;
		}
		
		Sequence currentUn = unSequenceChunkIt.next();
		resetChrName(currentUn);
		while(unSequenceChunkIt.hasNext()) {
			Sequence chunk = unSequenceChunkIt.next();
			if(chunk.getId().startsWith(currentUn.getId())) {
				currentUn.appendToSequence(chunk.getSequenceBases());
				chunk.unloadSequence();
			} else {
				fsio.write(currentUn, currentUn.getId() + "/chr" + currentUn.getId() + ".fa");
				resetChrName(chunk);
				currentUn = chunk;
			}
		}
		fsio.write(currentUn, currentUn.getId() + "/chr" + currentUn.getId() + ".fa");
		*/
	}
	
	private static void resetChrName(Sequence unChunk) {
		String chunkId = unChunk.getId();
		unChunk.setId(chunkId.substring(0, chunkId.indexOf(".")));
	}

	private void createSuperContigFasta(String output) throws Exception {
		Iterator<SuperContig> superIt = supers.iterator();

		SuperContig superContig = null;
		FastaSequenceIO fsio = new FastaSequenceIO(sequenceFile);
		while(superIt.hasNext()) {
			ArrayList<AgpEntry> contigs = new ArrayList<AgpEntry>();
			ArrayList<SuperContig> tmpSupers = new ArrayList<SuperContig>(NUM_CONTIG_TO_GROUP);
			int i = 0;
			while(i++ < NUM_CONTIG_TO_GROUP && superIt.hasNext()) {
				superContig = superIt.next();
				contigs.addAll(superContig.getContigs());
				tmpSupers.add(superContig);
			}
			fsio.extractRecordsIntoSequenceList(contigs);
			for(int j = 0; j < tmpSupers.size(); j++) {
				fsio.append(tmpSupers.get(j), output);
				tmpSupers.get(j).unloadSequence();				
			}
		}
	}
	
	private void createChromosomeBasedFastaByContigSize() throws Exception {
		Iterator<SuperContig> superIt = supers.iterator();

		FastaSequenceIO fsio = new FastaSequenceIO(sequenceFile);
		while(superIt.hasNext()) {
			SuperContig sc  = superIt.next();
			sc.loadContigsFromContigFileBySize(sequenceFile);
			System.out.println("Generating chr" + sc.getName() + " sequence bases length " + sc.getSequenceBases().length());
			System.out.print("..... done memory " + Runtime.getRuntime().totalMemory());
			File superDir = new File(sc.getName());
			if(!superDir.exists()) {
				superDir.mkdir();
			}
			BufferedWriter bw = new BufferedWriter(new FileWriter(superDir.getAbsolutePath() + "/chr"+sc.getName() + ".fa"));
			fsio.write(sc, bw);
			bw.close();
			sc.unloadSequence();
		}
	}
	
	private void createChromosomeBasedFasta(String agpContigPrefix) throws IOException {
		Iterator<SuperContig> superIt = supers.iterator();

		FastaSequenceIO fsio = new FastaSequenceIO(sequenceFile);
		while(superIt.hasNext()) {
			SuperContig sc  = superIt.next();
			sc.loadContigsFromContigFileByName(sequenceFile, agpContigPrefix);
			System.out.println("Generating chr" + sc.getName() + " sequence bases length " + sc.getSequenceBases().length());
			System.out.print("..... done memory " + Runtime.getRuntime().totalMemory());
			File superDir = new File(sc.getName());
			if(!superDir.exists()) {
				superDir.mkdir();
			}
			BufferedWriter bw = new BufferedWriter(new FileWriter(superDir.getAbsolutePath() + "/chr"+sc.getName() + ".fa"));
			fsio.write(sc, bw);
			bw.close();
			sc.unloadSequence();
		}
		
	}

}
