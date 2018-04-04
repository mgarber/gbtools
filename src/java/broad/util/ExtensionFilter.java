package broad.util;

import java.io.File;
import java.io.FilenameFilter;

public class ExtensionFilter implements FilenameFilter {
	// use:  File myFiles[] = new File(dirname).listFiles(new ExtensionFilter(".fq"));
	// This will list only .fq files.
	
	String ext;
	public ExtensionFilter(String ext) {
		this.ext = ext;
	}
	public boolean accept(File dir, String name) {
		return (name.endsWith(this.ext));
	}
}
