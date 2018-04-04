package broad.core.alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public class MegablastReader extends AlignmentList<MegablastAlignment>{
	private File source;
	
	public MegablastReader(String fileName) throws FileNotFoundException {
		super();
		source = new File(fileName);
		BufferedReader br = new BufferedReader(new FileReader(source));
		String line;
		try {
			while((line = br.readLine()) != null) {
				if(line.startsWith("#")){
					continue;
				}
				MegablastAlignment summary = new MegablastAlignment(line.split("\t"));
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
