package broad.pda.assembly;

public class AgpEntryFactory {

	public static AgpEntry createEntry(String [] rawInfo) {
		AgpEntry entry = null;
		String parentName = rawInfo[0];
		
		entry = new AgpEntry(parentName);
		entry.setName(rawInfo[5]);
		entry.setInReverseOrientation(rawInfo.length == 9 && "-".equals(rawInfo[8]));	
		entry.setStart(Integer.parseInt(rawInfo[1]));
		entry.setEnd(Integer.parseInt(rawInfo[2]));
		entry.setNumber(Integer.parseInt(rawInfo[3]));
		if(rawInfo[0].startsWith("Un") || rawInfo[0].startsWith("un")) {
			entry.setChromosome("Un");
		} else {
			entry.setChromosome(rawInfo[0].length() < 4 ? rawInfo[0] : rawInfo[0].substring(3)); //try to handle both chrNN and NN notations for chromosomes
			//System.out.println("Entry's chromosome " + entry.getChromosome());
			if(entry.getChromosome().startsWith("0")) {
				entry.setChromosome(entry.getChromosome().substring(1));
			}
		}
		if ("F".equals(rawInfo[4])) {			
			entry.setName(rawInfo[5]);
			entry.setType(AgpEntry.CLONE_TYPE);
		} else if ("centromere".equalsIgnoreCase(rawInfo[6])) {
			entry.setName("centromere");
			entry.setType(AgpEntry.CENTROMERE_TYPE);	
		} else if ("short_arm".equalsIgnoreCase(rawInfo[6])) {
			entry.setName("short_arm");
			entry.setType(AgpEntry.SHORT_ARM_TYPE);	
		} else if ("N".equals(rawInfo[4])) {
			entry = new Gap(entry);
		} else if ("W".equalsIgnoreCase(rawInfo[4])) {
			entry.setType(AgpEntry.CONTIG_TYPE);
		} else {
			entry.setName("other");
			entry.setType(AgpEntry.OTHER_TYPE);
			//System.out.print("Can't determine what this AGP entry is type<" + rawInfo[4] +"> name <  rawInfo<"+rawInfo[5]);
		}
		
		//System.out.println("created entry " + entry);
		return entry;
	}


}
