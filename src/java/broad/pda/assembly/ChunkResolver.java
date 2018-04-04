package broad.pda.assembly;

import java.io.File;
import java.io.FilenameFilter;
import java.util.Iterator;

import broad.core.annotation.BasicGenomicAnnotation;

public class ChunkResolver {
	private File dir;
	private int chunkSize;
	private String suffix;

	public ChunkResolver(File topLevelDirectory, int chunkSize, String suffix) {
		super();
		this.dir = topLevelDirectory;
		this.chunkSize = chunkSize;
		this.suffix = suffix;
	}
	
	public Chunk getChunk(String chr, int position, int maxSize) {
		if (chr.length() < 3) {
			chr = "chr" + chr;
		}
		
		int chunkStart = (position / chunkSize) * chunkSize;
		
		Chunk chunk = new Chunk(chr, chunkStart, chunkSize, dir.getAbsolutePath(), suffix, maxSize);
		
		return chunk;
	}
	
	public Iterator<Chunk> iterator(String chr, int l) {
		return new ChunkIterator(chr, l, this);	
	}
	
	public static class Chunk extends BasicGenomicAnnotation{
		File location;
		public Chunk(String chr, int start, int chunkSize, String path, final String suffix, int maxSize) {
			setStart(start);
			setEnd(Math.min(start + chunkSize - 1, maxSize));
			setChromosome(chr);
			String name = "chr"+chr.replace("chr", "") + "_" + getStart() + "-" + getEnd();
			setName(name);
			location = new File(path + "/"/*chr"*/ + chr.replace("chr", "") + "/" + name + suffix);
			if(!location.exists()) {
				final String startName = "chr"+chr.replace("chr", "") + "_" + getStart() + "-";
				System.err.println("WARN: Next chunk " + location.getName() + " NOT FOUND trying finding chunk by start " + startName);
				File dir = new File(path+"/"/*chr"*/ + chr.replace("chr", "") + "/");
				String [] fileList = dir.list(new FilenameFilter() {
					public boolean accept(File dir, String name) {
						
						return name.startsWith(startName) && name.endsWith(suffix);
					}
				});
				//System.err.println("List for files starting with startName has " + fileList.length + " files");
				if(fileList.length == 1) {
					location = new File(path + "/" +  chr.replace("chr", "") + "/" +fileList[0]);
					setName(fileList[0]);
					String [] fileInfo = fileList[0].split("-");
					setEnd(Integer.valueOf(fileInfo[1].replace(suffix, "")));
				}
			}
		}
		
		public File getLocation() {
			return location;
		}
		
		public int getOffset(int genomicLocation) {
			int offset = 0;
			if(genomicLocation >= getStart() && genomicLocation <= getEnd()){
				offset =  genomicLocation - getStart();
			} else {
				throw new IllegalArgumentException("Location " + genomicLocation + " does not fall with chunk " + toString());
			}
			
			return offset;
		}
		
		
	}
	
	public static class ChunkIterator implements Iterator<Chunk> {
		ChunkResolver cr;
		int maxSize;
		int nextStart;
		String chr;
		
		ChunkIterator(String chr, int maxSize, ChunkResolver cr) {
			this.chr = chr;
			this.maxSize = maxSize;
			this.cr = cr;
		}
		
		public boolean hasNext() {
			return nextStart < maxSize;
		}

		public Chunk next() {
			Chunk next = cr.getChunk(chr, nextStart, maxSize);
			nextStart = next.getEnd() + 1;
			return next;
		}

		public void remove() {
			//Do nothink, can't really remove a chunk;			
		}
		
	}
}
