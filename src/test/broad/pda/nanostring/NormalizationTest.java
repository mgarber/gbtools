package broad.pda.nanostring;

import java.io.FileInputStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import junit.framework.TestCase;

public class NormalizationTest extends TestCase {
	public static final String testData = "plate_13_a1_h8_96.txt";
	
	public void testSpikedInNormalization() throws Exception {
		URL trfDataUrl = getClass().getResource(testData);
		FileInputStream fis = new FileInputStream(trfDataUrl.getFile());
		NanostringReader reader = new NanostringReader();
		try {
			reader.load(fis);
		} finally {
			fis.close();
		}
		
		
		assertEquals("Unexpected number of experiments loaded ", 96, reader.getExperiments().size());
		Experiment first = reader.getExperiments().iterator().next();
		assertEquals("Unexpected number of observations",128, first.observations().size());
		
		reader.normalizeBySpikeIns();
		reader.writeSpikedInNormalized(System.out, 0, false,0,0);
	}
	
	public void testExperimentalControlNormalization() throws Exception {
		URL trfDataUrl = getClass().getResource(testData);
		FileInputStream fis = new FileInputStream(trfDataUrl.getFile());
		NanostringReader reader = new NanostringReader();
		try {
			reader.load(fis);
		} finally {
			fis.close();
		}
		
		List<String> controlTranscripts = new ArrayList<String>();
		controlTranscripts.add("Ik");
		controlTranscripts.add("Tbca");
		controlTranscripts.add("Ndufa7");
		controlTranscripts.add("Ndufs5");
		controlTranscripts.add("Ywhaz");
		controlTranscripts.add("Mea1");
		controlTranscripts.add("Shfm1");
		controlTranscripts.add("Tomm7");
		
		reader.setControlTranscripts(controlTranscripts);
		
		assertEquals("Unexpected number of experiments loaded ", 96, reader.getExperiments().size());
		Experiment first = reader.getExperiments().iterator().next();
		assertEquals("Unexpected number of observations",120, first.observations().size());
		
		reader.normalizeBySpikeIns();
		reader.normalizeByNormalizationTranscripts();
		
		reader.writeSpikedInNormalized(System.out, 0, false,0,0);
	}
	
	public void testZScore() throws Exception {
		URL trfDataUrl = getClass().getResource(testData);
		FileInputStream fis = new FileInputStream(trfDataUrl.getFile());
		NanostringReader reader = new NanostringReader();
		try {
			reader.load(fis);
		} finally {
			fis.close();
		}
		
		List<String> controlTranscripts = new ArrayList<String>();
		controlTranscripts.add("Ik");
		controlTranscripts.add("Tbca");
		controlTranscripts.add("Ndufa7");
		controlTranscripts.add("Ndufs5");
		controlTranscripts.add("Ywhaz");
		controlTranscripts.add("Mea1");
		controlTranscripts.add("Shfm1");
		controlTranscripts.add("Tomm7");
		
		reader.setControlTranscripts(controlTranscripts);
		
		assertEquals("Unexpected number of experiments loaded ", 96, reader.getExperiments().size());
		Experiment first = reader.getExperiments().iterator().next();
		assertEquals("Unexpected number of observations",120, first.observations().size());
		
		reader.normalizeBySpikeIns();
		reader.normalizeByNormalizationTranscripts();
		
		int controls [] = { 5,11,17,23,29,35,41,47,48,54,60,66,72,78,84,90};
		reader.setControls(controls);
		reader.computeZScores();
		
		reader.writeSpikedInNormalized(System.out, 2, false,0,0);
	}
	
	public void testZScoreConfidence() throws Exception {
		URL trfDataUrl = getClass().getResource(testData);
		FileInputStream fis = new FileInputStream(trfDataUrl.getFile());
		NanostringReader reader = new NanostringReader();
		try {
			reader.load(fis);
		} finally {
			fis.close();
		}
		
		List<String> controlTranscripts = new ArrayList<String>();
		controlTranscripts.add("Ik");
		controlTranscripts.add("Tbca");
		controlTranscripts.add("Ndufa7");
		controlTranscripts.add("Ndufs5");
		controlTranscripts.add("Ywhaz");
		controlTranscripts.add("Mea1");
		controlTranscripts.add("Shfm1");
		controlTranscripts.add("Tomm7");
		
		reader.setControlTranscripts(controlTranscripts);
		
		assertEquals("Unexpected number of experiments loaded ", 96, reader.getExperiments().size());
		Experiment first = reader.getExperiments().iterator().next();
		assertEquals("Unexpected number of observations",120, first.observations().size());
		
		reader.normalizeBySpikeIns();
		reader.normalizeByNormalizationTranscripts();
		
		int controls [] = { 5,11,17,23,29,35,41,47,48,54,60,66,72,78,84,90};
		reader.setControls(controls);
		reader.computeZScores();
		reader.computeConfidenceScores(2, 1);
		
		reader.writeSpikedInNormalized(System.out, 1,false,0.8,0.8);
	}
	
	public void testZScoreConfidenceByControlTranscripts() throws Exception {
		URL trfDataUrl = getClass().getResource(testData);
		FileInputStream fis = new FileInputStream(trfDataUrl.getFile());
		NanostringReader reader = new NanostringReader();
		try {
			reader.load(fis);
		} finally {
			fis.close();
		}
		
		List<String> controlTranscripts = new ArrayList<String>();
		controlTranscripts.add("Ik");
		controlTranscripts.add("Tbca");
		controlTranscripts.add("Ndufa7");
		controlTranscripts.add("Ndufs5");
		controlTranscripts.add("Ywhaz");
		controlTranscripts.add("Mea1");
		controlTranscripts.add("Shfm1");
		controlTranscripts.add("Tomm7");
		
		reader.setControlTranscripts(controlTranscripts);
		
		String [] fdrTranscripts = new String[controlTranscripts.size() + 1];
		for(int i = 0; i < controlTranscripts.size(); i++) {
			fdrTranscripts[i] = controlTranscripts.get(i);
		}
		fdrTranscripts[controlTranscripts.size()] = "Gapdh";
  		
		assertEquals("Unexpected number of experiments loaded ", 96, reader.getExperiments().size());
		Experiment first = reader.getExperiments().iterator().next();
		assertEquals("Unexpected number of observations",120, first.observations().size());
		
		reader.normalizeBySpikeIns();
		reader.normalizeByNormalizationTranscripts();
		
		int controls [] = { 5,11,17,23,29,35,41,47,48,54,60,66,72,78,84,90};
		reader.setControls(controls);
		reader.computeZScores();
		reader.computeConfidenceScoresUsingTranscripts(fdrTranscripts, 1);
		reader.writeSpikedInNormalized(System.out, 1,false,0.8,0.8);
	}
}
