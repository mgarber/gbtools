package broad.pda.chromosome;

import java.io.File;


public class GenericOrganism extends Mammal {

	public GenericOrganism(File sequenceDirectory) throws Exception {
		super(sequenceDirectory);
		// TODO Auto-generated constructor stub
	}

	@Override
	public String getLetterCode() {
		return "Gr";
	}

	@Override
	public String getName() {
		return "Generic";
	}

	@Override
	public File getPreferredSequenceDir() {
		// TODO Auto-generated method stub
		return null;
	}

}
