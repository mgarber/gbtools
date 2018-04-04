package broad.pda.seq.alignment;

import java.io.File;


//Writing my own alignment pipeline using Bowtie and BLAT
public class BowtieBLATPipeline {

	public BowtieBLATPipeline(File fastqFile){
		//Bowtie to genome
		//Bowtie to all junctions
		bowtieToJunctions();
		
		//BLAT remainder
		//Tophat remainder
		//Bowtie remaining reads to all new BLAT/Tophat junctions
	}

	private void bowtieToJunctions() {
		String command="bowtie /seq/rinnscratch/mguttman/ReferenceGenomes/MM9ESTs/junctions unaligned.fq ES.BowtieJunctions.sam --sam -k 3 -m 1 --un unaligned.junctions.fq --max multimappers.junctions.fq";
		
	}
	
	private void gsnap(){
		String command="/seq/annotation/bio_tools/GMAP/gmap-gsnap-pmap-2010-02-03/gmap-2010-02-03/bin/gsnap -d mouse -D /seq/annotation/bio_tools/GMAP/gmap-gsnap-pmap-2010-02-03/gmap-2010-02-03/share/mouse/ -T 1 -N 1 /broad/shptmp/mguttman/ES/BLAT/42RH8_s8.0.fq/42RH8_s8.0.fq.fa > ES.0.gsnap";
	}
	
	String blatCommand="blat /seq/genome/mouse/mouse_Mm9/19/chr19.fa seq.mdust temp2/chr19.maskLower.blat -minScore=50 -minIdentity=93 -mask=lower -repeats=lower";

	String bowtie="bowtie ES.BLAT.junctions /seq/rinnscratch/mguttman/RNASeq/ES/Illumina/fastq/42RH8_s8.fq BowtieV2/ES.BLATJunctions.bowtie -m 1 -a --best --strata --un BowtieV2/ES.BLATJunctions.unaligned.fq --max BowtieV2/ES.BLATJunctions.multimappers.fq";
	String nameSort="sort BowtieV2/ES.All.sam -o BowtieV2/ES.All.nameSort.sam -T /broad/shptmp/";
}
