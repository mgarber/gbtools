package broad.core.alignment;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import broad.core.annotation.LightweightGenomicAnnotation;


public class PHTabularReader extends AlignmentList<PHOneLineAlignment>{

	
	public PHTabularReader(String alignmentFile) throws IOException {
		super();
		loadAll(alignmentFile);
	}
	
	public PHTabularReader() {
		super();
	}
	
	public void loadAll(String alignmentFile) throws IOException {
		loadAll(alignmentFile, 0);
	}
	
	public void loadAll(String alignmentFile, int shiftQuery) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(alignmentFile));
		String line;

		while((line = br.readLine()) != null) {
			if(line.startsWith("#") || line.trim().length() == 0  || line.contains("****")){
				continue;
			}
			//System.out.println(line);
			PHOneLineAlignment summary = new PHOneLineAlignment(line.split("\t+"));
			LightweightGenomicAnnotation query = summary.getA();
			query.setStart(query.getStart() + shiftQuery);
			query.setEnd(query.getEnd() + shiftQuery);
			//System.out.println("Loaded " + summary);
			addAlignment(summary);
		} 
		br.close();
		
	}

}

