package broad.pda.genomicplots;

import java.awt.BasicStroke;
import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.jibble.epsgraphics.EpsGraphics2D;

import broad.core.alignment.AlignmentSummary;
import broad.core.alignment.AlignmentSummaryReader;
import broad.core.alignment.AlignmentTableSummaryReader;
import broad.core.alignment.Repeat;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.chromosome.Chromosome;

public class ManyPairwiseAlignmentPlot {
	private static final long serialVersionUID = 31344562L;
	static final int RESOLUTION = 10000;	
	public static int WIDTH = 5000;
	public static int HIGHT = 5000;
	public static int RIGHT_MARGIN = 50;
	public static int VERTICAL_SPACING = 100;
	public static int INTER_VERTICAL_SPACING = 40;
	
	public static String USAGE = "Usage: ManyPairwiseAlignmentPlot TASK=<task_num> <task_args>\n" +
	"\tTasks:\n" +
	"\t\t1. Draw Alignments: OUT=output file SCAFFOLD=scaffold AGP file SEQ1=aligned sequence <SEQ2= .... SEQN=> ALN1=alignment of seq1 and scaffold <ALN2= ... ALNN=>\n";

	private Chromosome scaffold;
	private List<AlignmentTableSummaryReader> scaffoldAligns;
	private List<Chromosome> alignedSeqs;
	
	
	public ManyPairwiseAlignmentPlot(Chromosome scaffold, 
			List<AlignmentTableSummaryReader> scaffoldAligns,
			List<Chromosome> alignedSeqs) {
		super();
		this.scaffold = scaffold;
		this.scaffoldAligns = scaffoldAligns;
		this.alignedSeqs = alignedSeqs;
	}
	
	public void plot(String outputFile) throws FileNotFoundException, IOException {
		Color alignedRegionColor = Color.GREEN;
		Color repeatColor        = Color.LIGHT_GRAY;
		Color gapColor           = Color.RED;

		// FIRST plot reference
		int yStart = 150;
		int length = (int)Math.min(WIDTH, scaffold.getSize());
		float scale = ((float)length) / scaffold.getSize();
		
		double powerOfTen = Math.floor(Math.log10(scaffold.getSize()));		
		float tickSpacing = (int) Math.pow(10, powerOfTen)/4;
		int tickNumber  = (int) Math.floor(scaffold.getSize()/tickSpacing);
		
		EpsGraphics2D g2d = new EpsGraphics2D("Multiple Pairwise Plot", new FileOutputStream(outputFile), 0, 0, (int) 1.5*WIDTH,(int) 1.2*HIGHT);
		java.awt.Font psFont = new java.awt.Font("arial", 0, 18);
		g2d.setFont(psFont);
		g2d.setStroke(new BasicStroke(8.0f));

		g2d.drawLine(RIGHT_MARGIN, yStart , length + RIGHT_MARGIN, yStart );		
		for(int i = 0 ; i <= tickNumber; i++) {
			int pos = (int) (i * tickSpacing);
			int plotPos = (int) (pos * scale);
			g2d.drawLine(RIGHT_MARGIN + plotPos, yStart - 20, RIGHT_MARGIN + plotPos, yStart);
			g2d.drawString(String.valueOf(pos), RIGHT_MARGIN + plotPos - 50, yStart - 40);
		}
		
		g2d.setColor(gapColor);
		Iterator<? extends GenomicAnnotation> gapIt = scaffold.Gaps();
		while(gapIt.hasNext()) {
			LightweightGenomicAnnotation annot = gapIt.next();
			g2d.drawLine(RIGHT_MARGIN + (int) (annot.getStart() * scale), yStart, RIGHT_MARGIN + (int) (annot.getEnd() * scale), yStart);
		}
		 
		Iterator<Repeat> repeatIt = scaffold.getRepeatReader().getRepeatList().iterator();
		System.out.println("repeats in scaffold " + scaffold.getRepeatReader().getRepeatList());
		g2d.setColor(repeatColor);
		while(repeatIt.hasNext()) {
			LightweightGenomicAnnotation annot = repeatIt.next();
			g2d.drawLine(RIGHT_MARGIN + (int) (annot.getStart() * scale), yStart, RIGHT_MARGIN + (int) (annot.getEnd() * scale), yStart);
			//System.out.println("Plotted reference repeat" + annot);
		}
		g2d.setColor(Color.BLACK);			

		// Plot aligned sequences
		for(int j = 0; j < scaffoldAligns.size(); j++) {
			System.out.println("Processing assembly " + j);
			int yBase = yStart + (j+2)*VERTICAL_SPACING;
			
			Chromosome seq = alignedSeqs.get(j);
			gapIt = seq.Gaps();
			repeatIt = seq.getRepeatReader().getRepeatList().iterator();

			AlignmentSummaryReader asr = scaffoldAligns.get(j);
			List<AlignmentSummary> alns = asr.getAlignmentList();
			Iterator<AlignmentSummary> alnIt = alns.iterator();
			long start = alns.get(0).getSubjectStart();
			long end   = alns.get(alns.size() - 1).getSubjectEnd();
			g2d.setColor(Color.BLACK);
			g2d.drawLine(RIGHT_MARGIN + (int) (start * scale), 
					yBase,
					RIGHT_MARGIN + (int) (end * scale), 
					yBase);
			
			/*
			g2d.setColor(gapColor);
			while(gapIt.hasNext()) {
				GenomicAnnotation annot = gapIt.next();
				g2d.drawLine(RIGHT_MARGIN + (int) (annot.getStart() * scale), yBase, RIGHT_MARGIN + (int) (annot.getEnd() * scale), yBase);
				System.out.println("Plotted " + annot);
			}
			 
			g2d.setColor(repeatColor);
			while(repeatIt.hasNext()) {
				GenomicAnnotation annot = repeatIt.next();
				g2d.drawLine(RIGHT_MARGIN + (int) (annot.getStart() * scale), yBase, RIGHT_MARGIN + (int) (annot.getEnd() * scale), yBase);
				System.out.println("Plotted " + annot);
			}
			*/
			while(alnIt.hasNext()) {								
				AlignmentSummary al = alnIt.next();
				long algnStart = al.getQueryStart();
				long algnEnd   = al.getQueryEnd();
				System.out.println("Processing " + al);
				// Now annotate the alignments...
				g2d.setColor(alignedRegionColor);
				g2d.drawLine(RIGHT_MARGIN + (int) (algnStart * scale), 
						yBase, RIGHT_MARGIN + (int) (algnEnd * scale), 
						yBase);
			}
		}
		g2d.flush();
		g2d.close();
	}

	
	public static void main(String [] args) throws Exception {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE);
		
		String scaffoldAGPFile = argMap.get("SCAFFOLD");
		Chromosome scaffold = new Chromosome(scaffoldAGPFile);
		scaffold.loadSequence();
		String scaffoldRepeatFile = scaffoldAGPFile.substring(0,scaffoldAGPFile.lastIndexOf(".")) + ".fa.out";
		scaffold.loadRepeatInfo(new File(scaffoldRepeatFile), false);
		
		int i = 1;
		ArrayList<Chromosome> alignedSeqs = new ArrayList<Chromosome>();
		List<AlignmentTableSummaryReader> alignments = new ArrayList<AlignmentTableSummaryReader>();
		
		while(true) { 
			String alignedSeqAGP = argMap.get("SEQ"+i);
			String alignedSeqAlnFile = argMap.get("ALN"+i);
			if(alignedSeqAlnFile == null) {
				break;
			}
			AlignmentTableSummaryReader asr = new AlignmentTableSummaryReader(alignedSeqAlnFile);
			Chromosome seq = new Chromosome(alignedSeqAGP);
			seq.loadSequence();
			String seqRepeatInfoFile = alignedSeqAGP.substring(0,alignedSeqAGP.lastIndexOf(".")) + ".fa.out";
			seq.loadRepeatInfo(new File(seqRepeatInfoFile), false);
			alignments.add(asr);
			alignedSeqs.add(seq);
			i++;
		}
		
		ManyPairwiseAlignmentPlot mpap = new ManyPairwiseAlignmentPlot(scaffold, alignments, alignedSeqs);
		
		mpap.plot(argMap.get("OUT"));
	}

}
