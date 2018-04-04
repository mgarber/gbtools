package broad.pda.ribosome;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.broad.tribble.FeatureSource;
import org.broad.tribble.bed.BEDCodec;
import org.broad.tribble.bed.BEDFeature;
import org.broad.tribble.source.BasicFeatureSource;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import broad.pda.datastructures.Alignments;

/**
   * @author Jim Robinson
   * @date 12/7/11
   */
public class VCFExample {

      public VCFExample(String vcfFile, String bed, String save) throws IOException {
    	  boolean requiresIndex = true;
          FeatureSource<VariantContext> source = BasicFeatureSource.getFeatureSource(vcfFile, new VCFCodec(), requiresIndex);

         Iterator<VariantContext> iter = source.query("17", 35647006, 35647114);// chr17:35,647,006-35,647,114
	     while (iter.hasNext()) {
	    	 VariantContext variant = iter.next();
	    	 System.err.println(new Alignments("17", 35647006, 35647114).fullyContained(new Alignments(variant.getChr(), variant.getStart(), variant.getEnd())));
	        Map<String, Genotype> genotypes= variant.getGenotypes();
	        for(String sample: genotypes.keySet()){
	        	Genotype genotype=genotypes.get(sample);
	        	List<Allele> alleles=genotype.getAlleles();
	        	for(Allele allele: alleles){
	        		if(allele.isNonReference()){
	        			System.err.println(sample+" "+allele+" "+variant.getChr()+" "+variant.getStart()+" "+variant.getEnd());	
	        		}
	        	}
	        
	        }
	        
	      /*  for(Allele allele: alleles){
	        	System.err.println(allele);
	        	 System.err.println(allele.isNonReference()+" "+allele.getBaseString());
	         }
	         */
	     }
	     
	     
	     
         source.close();
	}


	public void testBED() throws IOException {

          String bedFile = "path to your indexed bed";
          FeatureSource<BEDFeature> source = BasicFeatureSource.getFeatureSource(bedFile, new BEDCodec(), true);

          Iterator<BEDFeature> iter = source.query("chr1", 100000, 200000);

          while (iter.hasNext()) {
              BEDFeature f = iter.next();
              //assertTrue(f.getEnd() >= 0 && f.getStart() <= 200);
          }

          source.close();
      }


      public void testVCF() throws IOException {

          String vcfFile = "path to your indexed VCF";
          boolean requiresIndex = true;
          FeatureSource<VariantContext> source = BasicFeatureSource.getFeatureSource(vcfFile, new VCFCodec(), true);

          Iterator<VariantContext> iter = source.query("chr1", 100000, 200000);

          while (iter.hasNext()) {
              VariantContext variant = iter.next();
              variant.getStart();
              variant.getAlleles();
              // etc etc
          }
          source.close();
      }
      
      public static void main(String[] args) throws IOException{
    	  if(args.length>2){
    		  String vcfFile=args[0];
    		  String bed=args[1];
    		  String save=args[2];
    		  new VCFExample(vcfFile, bed, save);
    	  }
    	  else{System.err.println(usage);}
      }

      static String usage=" args[0]=vcfFile \n args[1]=BED file \n args[2]=save";
      
}


