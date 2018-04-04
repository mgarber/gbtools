package broad.core.alignment;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;


public class AlignmentTableSummaryReader extends AlignmentSummaryReader {

	
	public AlignmentTableSummaryReader(String fileName) throws FileNotFoundException{
		super(fileName);
		BufferedReader br = new BufferedReader(new FileReader(fileName));
		String line;
		try {
			while((line = br.readLine()) != null) {
				if(line.startsWith("#") || line.startsWith("Query\tSubject")){
					continue;
				}
				AlignmentTableSummary summary = new AlignmentTableSummary(line.split("\t"));
				addAlignmentSummary(summary);
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				System.out.print("Closing "+fileName);
				br.close();
				System.out.print(" ..... Closed\n");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}


	public static class AlignmentTableSummary extends AlignmentSummary {
		protected AlignmentTableSummary(String [] rawData) {
			setQuery(rawData[0]);
			setSubject(rawData[1]);
			setQueryStart(Integer.parseInt(rawData[2]));
			setQueryEnd(Integer.parseInt(rawData[3]));
			setSubjectStart(Integer.parseInt(rawData[4]));
			setSubjectEnd(Integer.parseInt(rawData[5]));
			setPid(rawData[6].length() == 0 ? 0 : Float.parseFloat(rawData[6]));
			//setScore(rawData[7].length() == 0 ? 0 : Integer.parseInt(rawData[7]));
			//setQueryOrientation(rawData[8]);
			setSubjectOrientation("+");
		}
	}

}
