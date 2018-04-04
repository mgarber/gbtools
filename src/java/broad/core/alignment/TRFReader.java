package broad.core.alignment;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import broad.core.datastructures.IntervalTree;
import broad.core.error.ParseException;

public class TRFReader implements TRFHandler {
	private Map<String, IntervalTree<TRFHit>> sequenceTRFMap;
	private TRFHitFilter filter;
	
	public TRFReader() {
		sequenceTRFMap = new HashMap<String, IntervalTree<TRFHit>>();
	}
	
	public TRFReader(TRFHitFilter filter) {
		this();
		this.filter = filter;
	}


	public void load(InputStream is) throws IOException, ParseException {
		BufferedReader br = new BufferedReader(new InputStreamReader(is));
		TRFParser parser = new TRFParser();
		parser.parse(br, this);
	}
	
	public void sequence(String sequenceName) {
		IntervalTree<TRFHit> sequenceTree = new IntervalTree<TRFHit>();
		sequenceTRFMap.put(sequenceName, sequenceTree);
	}

	public void tandemRepeat( TRFHit tandemRepeat) throws ParseException {
		IntervalTree<TRFHit> seqTree = sequenceTRFMap.get(tandemRepeat.getContainingSequenceId());
		if(seqTree == null) {
			throw new ParseException("Sequence for this Tandem Repeat has not been seen for sequence " + tandemRepeat.getContainingSequenceId() + ". sequence() must be called before adding TRFHits, check file");
		}
		if(filter == null || filter.accept(tandemRepeat) ){
			seqTree.put(tandemRepeat.getStart(), tandemRepeat.getEnd(), tandemRepeat);
		}
	}
	
	public List<TRFHit> getSequenceTandemRepeats(String sequenceName) {
		IntervalTree<TRFHit> sequenceTRFs = sequenceTRFMap.get(sequenceName);
		if(sequenceTRFs == null) {
			return new ArrayList<TRFHit>();
		}
		
		List<TRFHit> sequenceTRFList = new ArrayList<TRFHit>(sequenceTRFs.size());
		Iterator<TRFHit> TRFHitIt = sequenceTRFs.valueIterator();
		while(TRFHitIt.hasNext()) {
			sequenceTRFList.add(TRFHitIt.next());
		}
		return sequenceTRFList;
	}
	
	

}
