package broad.core.alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;


public class AlignmentSummaryReader extends AlignmentList<AlignmentSummary>  {
	private File source;


	
	public AlignmentSummaryReader(String fileName) throws FileNotFoundException{
		super();
		source = new File(fileName);
		BufferedReader br = new BufferedReader(new FileReader(source));

		String line;
		try {
			while((line = br.readLine()) != null) {
				if(line.startsWith("#")){
					continue;
				}
				AlignmentSummary summary = new AlignmentSummary(line.split("\t"));
				addAlignment(summary);
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

}
