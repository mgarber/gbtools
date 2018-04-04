package broad.pda.feature.genome;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.ChromosomeImpl;
import org.broad.igv.feature.genome.ChromosomeCoordinate;
import org.broad.igv.feature.genome.Genome;

import broad.core.annotation.ShortBEDReader;
import broad.core.error.ParseException;
import broad.pda.annotation.BEDFileParser;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.slide.Slide;

/**
 * TODO:  integrate better with org.broad.igv.feature.genome;
 * @author engreitz
 *
 */
public class SimpleGenome implements Genome {
    private static Logger log = Logger.getLogger(Genome.class);

    private String id;
    private String displayName;
    private List<String> chromosomeNames;
    private LinkedHashMap<String, Chromosome> chromosomeMap;
    private long length = -1;
    private Map<String, Long> cumulativeOffsets;   
	private Map<String, Integer> chromosomeSizes;
	
	private Map<String, ShortBEDReader> maskedRegions;   // TODO: store this as some other class .. like GenericAlignmentDataModel
														 // TODO: this should be a single ShortBEDReader
	private Map<String, Integer> maskFileData;
	
	
	public SimpleGenome(Map<String, Integer> sizes, File[] maskFiles) throws IOException, ParseException {
		initialize(sizes, maskFiles);
	}
	
	public SimpleGenome(Map<String, Integer> sizes) throws IOException, ParseException {
		this(sizes, null);
	}
	
	public SimpleGenome(String sizeFile, File[] maskFiles) throws IOException, ParseException {
		Map<String, Integer> sizes = BEDFileParser.loadChrSizes(sizeFile);
		initialize(sizes, maskFiles);
	}
	
	public SimpleGenome(String sizeFile) throws IOException, ParseException {
		this(sizeFile, null);
	}
	
	
	public void initialize(Map<String, Integer> sizes, File[] maskFiles) throws IOException, ParseException {
		this.chromosomeNames = new ArrayList<String>();
		this.chromosomeMap = new LinkedHashMap<String, Chromosome>();
		this.maskedRegions = new LinkedHashMap<String, ShortBEDReader>();
		this.cumulativeOffsets = new HashMap<String, Long>();   
		
		for (String name : sizes.keySet()) {
			Chromosome chr = new ChromosomeImpl(name, sizes.get(name));
			chromosomeMap.put(name, chr);
			chromosomeNames.add(name);
			
			if (maskFiles != null) {
				maskedRegions.put(name, Slide.loadMaskedRegions(maskFiles, name));
			}
		}
		
		maskFileData = ContinuousDataAlignmentModel.parseMaskFiles(maskFiles);
	}
	
	
	public Map<String, Integer> getSizes() { return chromosomeSizes; }
	public Map<String, Integer> getMaskFileData() { return maskFileData; }
	public Map<String, ShortBEDReader> getMaskedRegions() { return maskedRegions; }
	
	public long getNumMaskedBases() { 
		long total = 0;
		for (Integer value : maskFileData.values()) total += value;
		return total;
	}
	
	public long getNumUnmaskedBases() { return getLength() - getNumMaskedBases(); }
	
	@Override
    public String getId() { return id; };
    
    @Override
    public String getHomeChromosome() { return chromosomeNames.get(0); };
    
    @Override
    public Chromosome getChromosome(String chrName) {
    	return chromosomeMap.get(chrName);
    }

    @Override
    public List<String> getChromosomeNames() {
        return chromosomeNames;
    }

    @Override
    public Collection<Chromosome> getChromosomes() {
        return chromosomeMap.values();
    }

    @Override
    public String getChromosomeAlias(String str) {
    	return str;   // this is 
    }

    @Override
    public long getLength() {
        if (length < 0) {
            length = 0;
            for (Chromosome chr : chromosomeMap.values()) {
                length += chr.getLength();
            }
        }
        return length;
    }
    

    @Override
    public long getCumulativeOffset(String chr) {
        Long cumOffset = cumulativeOffsets.get(chr);
        if (cumOffset == null) {
            long offset = 0;
            for (String c : getChromosomeNames()) {
                if (chr.equals(c)) {
                    break;
                }
                offset += getChromosome(c).getLength();
            }
            cumOffset = new Long(offset);
            cumulativeOffsets.put(chr, cumOffset);
        }
        return cumOffset.longValue();
    }


    /**
     * Covert the chromosome coordinate in BP to genome coordinates in KBP
     *
     * @param chr
     * @param locationBP
     * @return
     */
    @Override
    public int getGenomeCoordinate(String chr, int locationBP) {
        return (int) ((getCumulativeOffset(chr) + locationBP) / 1000);
    }

    /**
     * Convert the genome coordinates in KBP to a chromosome coordinate
     */
    @Override
    public ChromosomeCoordinate getChromosomeCoordinate(int genomeKBP) {

        long cumOffset = 0;
        for (String c : chromosomeNames) {
            int chrLen = getChromosome(c).getLength();
            if ((cumOffset + chrLen) / 1000 > genomeKBP) {
                int bp = (int) (genomeKBP * 1000 - cumOffset);
                return new ChromosomeCoordinate(c, bp);
            }
            cumOffset += chrLen;
        }

        String c = chromosomeNames.get(chromosomeNames.size() - 1);
        int bp = (int) (genomeKBP - cumOffset) * 1000;
        return new ChromosomeCoordinate(c, bp);
    }

    @Override
    public String getNextChrName(String chr) {
        List<String> chrList = getChromosomeNames();
        for (int i = 0; i < chrList.size() - 1; i++) {
            if (chrList.get(i).equals(chr)) {
                return chrList.get(i + 1);
            }
        }
        return null;
    }

    @Override
    public String getPrevChrName(String chr) {
        List<String> chrList = getChromosomeNames();
        for (int i = chrList.size() - 1; i > 0; i--) {
            if (chrList.get(i).equals(chr)) {
                return chrList.get(i - 1);
            }
        }
        return null;
    }

    /**
     * Return the nucleotide sequence on the + strand for the genomic interval.  This method can return null
     * if sequence is not available.
     *
     * @param chr
     * @param start  start position in "zero-based" coordinates
     * @param end  end position
     * @return  sequence, or null if not available
     */
    @Override
    public byte[] getSequence(String chr, int start, int end) {
    	return null; 
    }

    @Override
    public String getDisplayName() {
    	return displayName;
    }

    @Override
    public byte getReference(String chr, int pos) {
    	return 0;
    }
    
    @Override
    public void addChrAliases(Map<String, String> aliases) {}
    
    @Override
    public void setChromosomeMap(LinkedHashMap<String,Chromosome> map, boolean x) {}
    
    @Override
    public void loadUserDefinedAliases() {}

    
    public static void main(String[] args) {
    	// For testing purposes
    }
    
    
    /**
     * TODO:
     * getRandomRegion(size) // account for masked regions
     * getRandomRegion(size, pctMaskedAllowed)
     * getRandomRegion(size, chromosome)
     * getRandomRegion(size, region)
     */
}
