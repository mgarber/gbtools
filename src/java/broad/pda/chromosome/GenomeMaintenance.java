package broad.pda.chromosome;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.Stack;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import broad.core.alignment.RepeatMaskerReader;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.assembly.AgpEntry;



public class GenomeMaintenance {
	public static String USAGE = "Usage: GenomeMaintenance TASK=<task_num> <task_args>\n" +
	"\tTasks:\n" +
	"\t1. Generate size file SEQDIR=<sequence directory >\n" +
	"\t2. Generate agp files SEQDIR=<sequence directory for the assembly to work on\n" +
	"\t3. Generate artificial AGP file for a given sequence IN=<Fasta file with sequence to analyze>\n" +
	"\t4. Chunk -org organism -chunk <chunk size> -overlap <chunk overlap> -outdir <output directory> \n" +
	"\t\t-seqdir <sequence directory if not the standard version>]\n" +
	"\t5. Unchunk RepeatMasker chunked output -indir <chunks input directory> -chunk <chunk size> -outdir <output directory> -seqdir <org sequence directory>\n" +
	"\t6. Introduce Rearrangements -org <Organism> \n\t\t\t" +
		"-dups <# of duplications to introduce> -translocations <# translocations to introduce> -deletions <# deletions> -inversions <# inversions>\n\t\t\t"+
		"outdir <output directory for the rearranged genome> -orgDir <Organism sequence dir if not the default>\n" +
	"\t7. Generate random insert ends -org <organism> -insertSize <average Insert Size> -readLength <average read length> \n" + 
		"\t\t\t-readLengthVar <read length variance> -coverage <desired coverage eg. 2> [-orgDir <org sequence directory if not the default>]\n" +
	"\t8. Combined a multifasta file into a single fasta using Ns as separators (At this point the whole file need to be loaded into memory so it requires lots of it) -seqId <Sequence name for the combined data> -separatorSize <Number of Ns separating chunks> -in <Multifasta file (Ussually the multifasta Un chromosome) standard in if none specified> -out <Output file or standard out if non is psecified> ";
	
	public static void main(String [] args) throws Exception {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE);
		
		if("1".equals(argMap.getTask())) {
			Mammal m = new GenericOrganism(new File(argMap.get("SEQDIR")));
			m.writeSizes(m.getSequenceDir() + "/sizes");
		} else if ("2".equals(argMap.getTask())) {
			File  sequenceDirectory = new File(argMap.getMandatory("SEQDIR")); 
			String [] sequenceDirList = sequenceDirectory.list();
			FastaSequenceIO fsio = new FastaSequenceIO();
			for(int i = 0; i < sequenceDirList.length; i++) {
				File f =  new File(sequenceDirectory.getAbsolutePath() + "/" + sequenceDirList[i]);
				if(f.isDirectory()) {
					String [] chrList = f.list();
					for(int j = 0; j < chrList.length;  j++) {
						if(chrList[j].endsWith(".fa")){
							System.out.println("Doing " + f.getAbsolutePath() + "/" + chrList[j]);
							List<AgpEntry> agp = getAGP(new File(f.getAbsolutePath() + "/" + chrList[j]));
							Iterator<AgpEntry> it = agp.iterator();
							String fileName = chrList[j].substring(0,chrList[j].lastIndexOf("."));
							String agpFile = fileName + ".agp";
							BufferedWriter agpBW = new BufferedWriter(new FileWriter(f.getAbsolutePath() + "/" + agpFile));
							while(it.hasNext()) {
								agpBW.write(it.next().toString());
								agpBW.newLine();
							}
							agpBW.close();
						}
					}
				}
			}
		} else if("3".equals(argMap.getTask())) { 
			String input = argMap.getInput();
			FastaSequenceIO fsio = new FastaSequenceIO();
			String output = CLUtil.replaceFileExtension(input, "agp");
			
			List<AgpEntry> agp = getAGP(new File(input));
			Iterator<AgpEntry> it = agp.iterator();
			BufferedWriter agpBW = new BufferedWriter(new FileWriter(output));
			while(it.hasNext()) {
				agpBW.write(it.next().toString());
				agpBW.newLine();
			}
			agpBW.close();			
	    }else if("4".equals(argMap.getTask())) {
			String org = argMap.getMandatory("org");
			String sequenceDirectory = argMap.getMandatory("seqdir");
			String outdir = argMap.getOutputDir();
			int chunkSize = argMap.getInteger("chunk");
			int chunkOverlap = argMap.getInteger("overlap");
			
			Mammal m = new GenericOrganism(new File(sequenceDirectory));
			
			FastaSequenceIO fsio = new FastaSequenceIO();
			Iterator<Chromosome> chrIt = m.getAllChromosomes().iterator();
			while(chrIt.hasNext()) {
				Chromosome c = chrIt.next();
				c.loadSequence();
				System.out.println("Processing chromosome " + c.getSymbol());
				File chrDir = new File(outdir + "/" + c.getSymbol());
				if(!chrDir.exists()) {
					chrDir.mkdir();
				}
				Iterator<SequenceRegion> chrChunkIt = c.chunk(chunkSize, chunkOverlap).iterator();
				while(chrChunkIt.hasNext()) {
					SequenceRegion chunk = chrChunkIt.next();
					fsio.write(chunk, chrDir.getAbsolutePath() + "/chr" + c.getSymbol() + "_" + chunk.getRegionStart() +"_" + chunk.getRegionEnd() + ".fa");
				}
				c.unloadSequence();
			}
		} else if("5".equals(argMap.getTask())) { 
			int chunkSize = argMap.getInteger("chunksize");
			final String chunkExt = "fa.out";
			final Pattern extPat = Pattern.compile("\\." + chunkExt.replace("\\.", "\\.") + "$");
			File inDir = new File(argMap.getInputDir());
			String outDir = argMap.getOutputDir();
			Mammal m = new GenericOrganism(new File(argMap.get("seqdir")));

			System.out.println("Ext matching pattern " + extPat);
			
			Iterator<Chromosome> chrIt = m.getAllChromosomes().iterator(); 
			//Iterator<Chromosome> chrIt = org.getUn().iterator();
			while(chrIt.hasNext()) {
				Chromosome c = chrIt.next();
				System.out.println("processing chr" + c.getSymbol());
				File chrDir = new File(inDir.getAbsolutePath() + "/" + c.getSymbol());
				File [] inDirList = chrDir.listFiles(new FilenameFilter() {
	
					public boolean accept(File dir, String fileName) {
						Matcher m = extPat.matcher(fileName);
						return m.find();
					}
					
				});
				
				Arrays.sort(inDirList, new Comparator<File>(){
	
					public int compare(File arg0, File arg1) {
						String fName0 = arg0.getName().replace("."+chunkExt, "");
						String fName1 = arg1.getName().replace("."+chunkExt, "");
						
						String [] arg0Info = fName0.split("_");
						int start0 = Integer.parseInt(arg0Info[arg0Info.length - 2]);
						String [] arg1Info = fName1.split("_");
						int start1 = Integer.parseInt(arg1Info[arg1Info.length - 2]);
						return start0 - start1;
					}
					
				});
				
				// since we assume a naming convention we know how to get the orginal sequence name
				if(inDirList.length == 0 ) {
					System.err.println("No files in " + chrDir.getAbsolutePath() + " matched the given extension " + chunkExt);
					return;
				}
				String seqName = inDirList[0].getName().split("_")[0];			
				
				RepeatMaskerReader rmr = new RepeatMaskerReader();
				long totalSize = c.getSize();
	
				for(int i = 0; i < inDirList.length; i++) {
					File file = inDirList[i];
					String [] nameInfo = file.getName().split("_");
					int chunkStart = Integer.parseInt(nameInfo[nameInfo.length - 2].replace("."+chunkExt, "").split("-")[0]);
					//Lets check the chunks to make sure we are not missing any.
					if(chunkStart != (chunkSize * i )) {
						throw new RuntimeException("Expected chunk starting at " + (chunkSize * i) + " but got " + inDirList[i].getName() + " which starts at " + chunkStart);
					}
					rmr.loadRepeats(file, false, chunkStart, (int) totalSize, seqName);
				}
				
				File outChrDir = new File(outDir + "/" + c.getSymbol());
				if(!outChrDir.exists()) {
					outChrDir.mkdir();
				}
				BufferedWriter bw = new BufferedWriter(new FileWriter(outChrDir.getAbsolutePath() + "/chr" + c.getSymbol() + "." + chunkExt));
				rmr.writeHeader(bw);
				rmr.writeRepeatList(bw);
				bw.close();
			}
		}else if("6".equals(argMap.getTask())) { 
			String org = argMap.getMandatory("org");
			String outdir = argMap.getOutputDir();
			Mammal m = new GenericOrganism(new File(argMap.get("orgDir")));
			int numOfDups = argMap.isPresent("dups") ? argMap.getInteger("dups") : 0;
			int numOfTrans = argMap.isPresent("translocations") ? argMap.getInteger("translocations") : 0;
			int numOfInv   = argMap.isPresent("inversions") ? argMap.getInteger("inversions") : 0;
			int numOfDel   = argMap.isPresent("deletions") ? argMap.getInteger("deletions") : 0;
			
			Random randomizer = new Random(Math.round(Math.random()*1000000) );
			
			Iterator<Chromosome> chrIt = m.getAllChromosomes().iterator();
			while(chrIt.hasNext()) {
				chrIt.next().loadSequence();
			}
			for(int i = 0; i < numOfDups; i++) {
				SequenceRegion regionToDuplicate = m.drawRandomRegion(randomizer.nextInt(100) * 1000);
				SequenceRegion insertionPoint = m.drawRandomRegion(1);
				m.insert(regionToDuplicate, insertionPoint);
				System.out.println("Duplication " + regionToDuplicate + " at " + insertionPoint);
			}
			
			for(int i = 0; i < numOfTrans; i++) {
				SequenceRegion region1 = m.drawRandomRegion(randomizer.nextInt(100) * 1000);
				SequenceRegion region2 = m.drawRandomRegion(randomizer.nextInt(100) * 1000);
				
				m.translocate(region1, region2);
				System.out.println("translocation " + region1 + " and " + region2);
			}
			
			for(int i = 0; i < numOfInv; i++) {
				SequenceRegion region = m.drawRandomRegion(randomizer.nextInt(100) * 1000);
				m.invert(region);
				System.out.println("Inversion " + region);
			}
			
			for(int i = 0; i < numOfDel; i++) {
				SequenceRegion region = m.drawRandomRegion(randomizer.nextInt(100) * 1000);
				m.delete(region);
				System.out.println("Deletion " + region);
			}
			
			
			chrIt = m.getAllChromosomes().iterator();
			FastaSequenceIO fsio = new FastaSequenceIO();
			while(chrIt.hasNext()) {
				Chromosome c = chrIt.next();
				File outChrDir = new File(outdir + "/" + c.getSymbol());
				if(!outChrDir.exists()) {
					outChrDir.mkdir();
				}
				fsio.write(c.getSequence(), outChrDir.getAbsolutePath() + "/chr"+c.getSymbol() + ".fa");
			}
	    }else if("7".equals(argMap.getTask())) { 
			String org        = argMap.getMandatory("org");
			String outdir     = argMap.getOutputDir();
			Mammal m = new GenericOrganism(new File(argMap.get("orgDir")));
			int readLength    = argMap.getInteger("readLength");
			int readLengthVar = argMap.getInteger("readLengthVar");
			int insertSize    = argMap.getInteger("insertSize");
			float coverage      = argMap.getFloat("coverage");
			
			Random randomizer = new Random(Math.round(Math.random()*1000000) );
			int totalInserts = Math.round(m.getGenomeSize() /(float)insertSize * coverage);
			
			//System.out.println("Genome size, " + m.getGenomeSize() + " Total inserts to generate " + totalInserts);
			HashMap<String, List<SequenceRegion>> chrInsertEnds = new HashMap<String, List<SequenceRegion>>();

			for (int i = 0; i <= totalInserts; i++ ) {
				SequenceRegion r = m.drawRandomRegion(insertSize);
				List<SequenceRegion> insertEnds = chrInsertEnds.get(r.getChromosome());
				if(insertEnds == null) {
					insertEnds = new ArrayList<SequenceRegion>();
					chrInsertEnds.put(r.getChromosome(), insertEnds);
				}
				int leftReadSize = randomizer.nextInt(2*readLengthVar + 1) + (readLength - readLengthVar);
				SequenceRegion leftRead = new SequenceRegion(r.getContainingSequenceId(), r);
				leftRead.setEnd(r.getStart() + leftReadSize );
				leftRead.setName(i + "_left");
				
				int rightReadSize = randomizer.nextInt(2*readLengthVar + 1) + (readLength - readLengthVar);
				SequenceRegion rightRead = new SequenceRegion(r.getContainingSequenceId(), r);
				rightRead.setName(i + "_right");
				rightRead.setStart(r.getEnd() - rightReadSize );
				
				insertEnds.add(leftRead);
				insertEnds.add(rightRead);
				
			}
			
			System.out.println("Generated random inserts got " + chrInsertEnds.values().iterator().next().size() + " for one of them");
			
			Iterator<Chromosome> chrIt = m.getAllNonRandomChromosomes().iterator();

			FastaSequenceIO fsio = new FastaSequenceIO();
			while(chrIt.hasNext()) {
				Chromosome c = chrIt.next();
				System.out.print("Filling end reads for chr" + c.getSymbol() + " ... ");
				c.loadSequence();
				List<SequenceRegion> ends = chrInsertEnds.get(c.getSymbol());
				Iterator<SequenceRegion> endsIt = ends.iterator();
				
				while(endsIt.hasNext()) {
					c.getRegion(endsIt.next());
				}
				c.unloadSequence();
				System.out.println("done now writing the out");
				fsio.write(ends, outdir + "/chr" + c.getSymbol() +"_reads.fa");
			}
			
	    } else if("8".equals(argMap.getTask())) { 
	    	int separatorSize = argMap.getInteger("separatorSize");
	    	String seqId      = argMap.getMandatory("seqId");
			FastaSequenceIO fsio = new FastaSequenceIO();
			InputStream is = argMap.getInputStream();
			List<Sequence> chunks = fsio.loadAll(is);
			is.close();
			
			StringBuilder separatorBuilder = new StringBuilder(separatorSize);
			for(int i = 0; i < separatorSize; i++) {
				separatorBuilder.append("N");
			}
			String separator = separatorBuilder.toString();
			
			
			Iterator<Sequence> chunkIt = chunks.iterator();
			Sequence oneSeq = new Sequence(seqId);
			while(chunkIt.hasNext() ) {
				Sequence chunk = chunkIt.next();
				oneSeq.appendToSequence(chunk.getSequenceBases());
				if(chunkIt.hasNext()) {
					oneSeq.appendToSequence(separator);
				}
				chunk.unloadSequence();
			}
			BufferedWriter bw = argMap.getOutputWriter();
			fsio.write(oneSeq, bw);
			bw.close();			
	    } else {
			System.err.println(USAGE);
		}
	}
	
	public static List<AgpEntry> getAGP(File file) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line = br.readLine();
		String seqName = line.substring(1);
		if(seqName.startsWith("chr")) {
			seqName = seqName.substring(3);
		}
		Stack<AgpEntry> agp = new Stack<AgpEntry>();

		boolean inGap = true;

		int currEntrySize = 0;
		int position = 1;
		char [] lineChars = null;
		int num = 1;
		int contig = 1;
		while((line = br.readLine()) != null) {
			lineChars = line.toCharArray();
			for(int i = 0; i < lineChars.length; i++) {
				char c = Character.toUpperCase(lineChars[i]);
				if(c == 'N' && (!inGap || agp.isEmpty() )) {
					if(!agp.isEmpty()) {
						AgpEntry lastContig = agp.peek();
						lastContig.setEnd(position - 1);
					}
					AgpEntry newGap = new AgpEntry(seqName);
					newGap.setStart(position);
					newGap.setType(AgpEntry.GAP_TYPE);
					newGap.setNumber(num++);
					newGap.setName("fragment");
					agp.push(newGap);
					currEntrySize = 1;
					inGap = true;
				} else if(c == 'N' && inGap) {
					currEntrySize++;
				} else if (c != 'N' && !inGap) {
					currEntrySize++;
				} else if (c != 'N' && inGap  ) {
					if(!agp.isEmpty()) {
						AgpEntry lastGap = agp.peek();
						lastGap.setEnd(position - 1);
					}
					AgpEntry newContig = new AgpEntry(seqName);
					newContig.setStart(position);
					newContig.setType(AgpEntry.CONTIG_TYPE);
					newContig.setNumber(num++);
					newContig.setName("contig_" + contig++);
					agp.push(newContig);					
					currEntrySize = 1;
					inGap = false;
				}
				position++;
			}
			if(!agp.isEmpty()) {
				AgpEntry last = agp.peek();
				last.setEnd(position - 1);
			}
		}
		br.close();
		return agp;
	}

}
