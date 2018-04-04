package broad.core.gene;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import broad.core.annotation.AnnotationFactoryFactory;
import broad.core.annotation.AnnotationFactoryFactory.GFFFactory;
import broad.core.annotation.AnnotationHandler;
import broad.core.annotation.AnnotationReader;
import broad.core.annotation.GFF;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.GenomicAnnotationFilter;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.error.ParseException;



public class GFFGeneReader extends AnnotationReader<GeneAnnotation> {
	public GFFGeneReader() {
		super();
	}
	
	public GFFGeneReader(BufferedReader br) throws ParseException, IOException {
		super();
		load(br);
	}
	
	public GFFGeneReader(String input) throws ParseException, IOException {
		super();
		BufferedReader br = new BufferedReader(new FileReader(input));
		try {
			load(br);
		} finally {
			try {
				//System.out.print("Closing "+input);
				br.close();
				//System.out.print(" ..... Closed\n");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	public void load(BufferedReader br) throws IOException, ParseException {
		String line;
		String currentGene = null;
		List<GFF> currentGeneExons = new ArrayList<GFF>();
		GFFFactory exonFactory = AnnotationFactoryFactory.gffFactory;
		try {
			int tracks = 0;
			while((line = br.readLine()) != null) {
				line = line.trim();
				if(line.startsWith("#") || line.length() == 0) {
					continue;
				}
				
				if( line.toLowerCase().startsWith("browser")){
					addBrowserLine(line);
				} else if ( line.toLowerCase().startsWith("track")){ //Assuming only one track per annotation set, TODO: check!
					startTrack( line, tracks);
					tracks++;
				} else {
					String [] lineSplit = line.split("\t");
					GFF exon = exonFactory.create(lineSplit);
					List<String> transcriptAttributes = exon.getAttributes().get("transcript_id");
					String exonTranscript = transcriptAttributes == null || transcriptAttributes.size() == 0 ? exon.getName() : transcriptAttributes.get(0);
					if(!exonTranscript.equals(currentGene)) {
						if(currentGene!=null) {
							GeneAnnotation gene = new GeneAnnotation(currentGene);
							for(GFF e : currentGeneExons) {
								gene.addExon(e, e.getFrame());
							}
							gene.setStart(currentGeneExons.get(0).getStart());
							LightweightGenomicAnnotation lastExon = currentGeneExons.get(currentGeneExons.size() - 1);
							gene.setEnd(lastExon.getEnd());
							gene.setOrientation(lastExon.getOrientation());
							gene.setChromosome(lastExon.getChromosome());
							super.addAnnotationToLastTrack(gene);
						}
						currentGene = exonTranscript;
						currentGeneExons.clear();
						currentGeneExons.add(exon);
					} else {
						currentGeneExons.add(exon);
					}

				
		
				}
				
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				//System.out.print("Closing "+input);
				br.close();
				//System.out.print(" ..... Closed\n");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	public GeneAnnotation createAnnotation(GenomicAnnotation a) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int parse(String file,
			GenomicAnnotationFilter<GeneAnnotation> filter,
			AnnotationHandler handler) throws ParseException, IOException {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public int parse(BufferedReader br,
			GenomicAnnotationFilter<GeneAnnotation> filter,
			AnnotationHandler handler) throws ParseException, IOException {
		// TODO Auto-generated method stub
		return 0;
	}

}
