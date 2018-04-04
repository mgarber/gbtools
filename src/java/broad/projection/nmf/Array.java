package broad.projection.nmf;


import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Set;

import broad.core.math.EmpiricalDistribution;
import broad.projection.math.Matrix;
import broad.projection.math.NMFCostFunction;
import broad.projection.utils.ClusterHierarchical;
import broad.projection.utils.ColorBlend;
import broad.projection.utils.Gene;
import broad.projection.utils.Inc;
import broad.projection.utils.Infile;
import broad.projection.utils.Outfile;
import broad.projection.utils.Phenotype;
import broad.projection.utils.Sample;
import broad.projection.utils.Text;
import broad.projection.utils.TextTree;
import broad.projection.utils.Texts;
import broad.projection.utils.Tree;
import broad.projection.utils.Url;
import broad.projection.utils.svm;
import broad.projection.utils.svm_model;
import broad.projection.utils.svm_node;
import broad.projection.utils.svm_parameter;
import broad.projection.utils.svm_problem;

/*--------------------------------------------------------------------
 TODO
 ---------------------------------------------------------------------
 matching two non-matching genesets still returns 1, bug...
 remove all plotting stuff once the classes have been made

 get easy way to get sample color during radial plot. (probably generic color blend would help)
 color setter for phenotypes, list...
 make the pre/post example using abbreviations
 make a function that loads a samples stuff, first by phenotype, then by sample

 Get cls stuff working and get all the functions from microarray
 Fix maphomolog later to be able to handle multiple genes of same name
   in gct
 Make a phenotype selector that sorts by phenotype, hc?
 make cls reader be able to handle older versions and 0,1,2, etc...
 function to copy an array
 ---------------------------------------------------------------------
 IDEAS
 ---------------------------------------------------------------------
 ---------------------------------------------------------------------
 BUGS
 ---------------------------------------------------------------------
 clusterHierarchicalGene ghc is screwing up
--------------------------------------------------------------------*/

/**
 * Microarray data analysis
 *
 * @author Daniel Scanfeld
 * @version 1.0
 */

public class Array extends Inc {
/*--------------------------------------------------------------------
 Class Variables
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 Attributes
--------------------------------------------------------------------*/

    /**
     * Unique experiment id
     */
    public String id;

    /**
     * Experiment name
     */
    public Text name;

    /**
     * Short description of the experiment
     */
    public Text shortDescription;

    /**
     * Long description of the experiment
     */
    public Text longDescription;

/*--------------------------------------------------------------------
 Samples
--------------------------------------------------------------------*/

    /**
     * Number of Samples
     */
    public int numSample;

    /**
     * Sample Names
     */
    public Texts sampleName;

    /**
     * Sample objects
     */
    public Sample[] sample;

    /**
     * Sample clustering
     */
    public ClusterHierarchical shc;

/*--------------------------------------------------------------------
 Genes
--------------------------------------------------------------------*/

    /**
     * Number of Genes
     */
    public int numGene;

    /**
     * Gene Names
     */
    public Texts geneName;
    
    public Texts descriptions;

    /**
     * Gene Tree
     */
    public TextTree geneTree;

    /**
     * Gene objects
     */
    public Gene[] gene;

    /**
     * Missing values
     */
    public boolean[] missing;

    /**
     * Gene clustering
     */
    public ClusterHierarchical ghc;

    /**
     * Gene espression levels
     */
    public Matrix x;

/*--------------------------------------------------------------------
 Phenotypes
--------------------------------------------------------------------*/

    /**
     * Number of Phenotypes
     */
    public int numPhenotype;

    /**
     * Phenotype Names
     */
    public Texts phenotypeName;

    /**
     * Phenotype objects
     */
    public Phenotype[] phenotype;

/*--------------------------------------------------------------------
 Support Vector Machine
--------------------------------------------------------------------*/

    /**
     * Support Vector Machine
     */
    public svm sv;

    /**
     * SVM Parameters
     */
    public svm_parameter svp;

    /**
     * SVM models
     */
    public svm_model[] svmm;

    /**
     * SVM model indices
     */
    public int[] svmi;

    /**
     * SVM model max
     */
    public double svm_max;

/*--------------------------------------------------------------------
 Output
--------------------------------------------------------------------*/

    /**
     * The current PDF
     */
    //public PDF curPDF;

    /**
     * NMF Error string
     */
    public double nmferr;

    /**
     * NMF Trim Error string
     */
    public double nmferrtr;

    /**
     * Color Blend
     */
    public ColorBlend cb;

/*--------------------------------------------------------------------
 Output Flags
--------------------------------------------------------------------*/

    /**
     * Show samples below the score cutoff
     */
    public boolean showBelowCutoff;

    /**
     * Shade samples to gray
     */
    public boolean shadeToGray;

    /**
     * Show the heatmap on the plot
     */
    public boolean showHeatMap;

/*--------------------------------------------------------------------
 CONSTRUCTOR FUNCTIONS
--------------------------------------------------------------------*/

    /**
     * Null constructor
     */

    public Array() {
        init();  // Initializes the object
    }

    /**
     * Initializes with a GCT and CLS filename
     *
     * @param gctName GCT file name
     * @param clsName CLS fiel name
     */

    public Array(String gctName, String clsName) {
        init(); // Initializes the object
        if(clsName != null && clsName.length() > 0) {
        	load(gctName, clsName);  // Loads the GCT and CLS files
        } else {
        	load(gctName);
        }
    }

/*--------------------------------------------------------------------
 Array - GCT constructor
--------------------------------------------------------------------*/

    public Array(String gctName) {
        init();
        load(gctName);  // Load the GCT file
    }

/*--------------------------------------------------------------------
 Array - Create an array with a's genes, new samples, and phenotypes
--------------------------------------------------------------------*/

    public Array(Array a, Matrix ex, Texts ns, Texts np, Texts ph) {
        init();
        geneName = a.geneName.copy();
        descriptions=a.descriptions.copy();
        geneTree = a.geneTree;  // asdf - change to copy
        gene = a.gene;          // asdf - change to copy
        sampleName = ns;
        numGene = gene.length;
        numSample = ns.size;
        phenotypeName = np;
        numPhenotype = np.size;
        x = ex;

        sample = new Sample[numSample];
        missing = new boolean[numGene];
        phenotype = new Phenotype[numPhenotype];

        // LOOP THROUGH THE PHENOTYPES
        for (int i = 0; i < numPhenotype; i++) {
            phenotype[i] = new Phenotype(phenotypeName.get(i), this);  // Make the ith phenotype
        }

        // LOOP THROUGH THE SAMPLES
        for (int i = 0; i < numSample; i++) {
            sample[i] = new Sample(sampleName.get(i), i, this);  // Make the ith sample
        }

        // LOOP THROUGH THE GENES
        for (int i = 0; i < numGene; i++) {
            missing[i] = false;
        }

        // LOOP THROUGH THE SAMPLE
        for (int i = 0; i < numSample; i++) {
            sample[i].addPhenotype(ph.get(i));  // Set the phenotype
        }
    }

/*--------------------------------------------------------------------
 Array - Create an array with a's genes, new samples, and phenotypes
--------------------------------------------------------------------*/

    public Array(Matrix a, Texts gnme, Texts descriptions, Texts snme) {
        init();

        x = a.copy();
        geneName = gnme.copy();
        this.descriptions=descriptions.copy();
        sampleName = snme.copy();

        numGene = gnme.size;
        numSample = snme.size;

        phenotypeName = new Texts();
        phenotypeName.rpush("none");

        numPhenotype = phenotypeName.size;

        gene = new Gene[numGene];
        sample = new Sample[numSample];
        missing = new boolean[numGene];
        phenotype = new Phenotype[numPhenotype];

        // LOOP THROUGH THE PHENOTYPES
        for (int i = 0; i < numPhenotype; i++) {
            phenotype[i] = new Phenotype(phenotypeName.get(i), this);  // Make the ith phenotype
        }

        // LOOP THROUGH THE SAMPLES
        for (int i = 0; i < numSample; i++) {
            sample[i] = new Sample(sampleName.get(i), i, this);  // Make the ith sample
        }

        // LOOP THROUGH THE GENES
        for (int i = 0; i < numGene; i++) {
            gene[i] = new Gene(geneName.get(i), geneName.get(i), this);
            missing[i] = false;
        }

        // LOOP THROUGH THE SAMPLE
        for (int i = 0; i < numSample; i++) {
            sample[i].addPhenotype(phenotypeName.get(0));
        }
    }
    
    
   
    
    /*--------------------------------------------------------------------
    Array - Create an array with a's genes, new samples, and phenotypes
   --------------------------------------------------------------------*/

       public Array(Matrix a, Texts gnme, Texts snme) {
          this(a, gnme, gnme, snme);
       }
    

/*--------------------------------------------------------------------
 init - Initialize the array
--------------------------------------------------------------------*/

    public void init() {
        showBelowCutoff = false;
        shadeToGray = false;
        showHeatMap = true;
    }

/*--------------------------------------------------------------------
 copy - Copy the array
--------------------------------------------------------------------*/

    public Array copy() {
        Array a = new Array();
        a.id = id;
        a.name = name;
        a.shortDescription = shortDescription;
        a.longDescription = longDescription;
        a.sampleName = sampleName.copy();
        a.geneName = geneName.copy();
        a.descriptions=descriptions.copy();
        a.phenotypeName = phenotypeName.copy();

        a.numGene = numGene;
        a.numSample = numSample;
        a.numPhenotype = numPhenotype;

        a.gene = new Gene[numGene];
        a.sample = new Sample[numSample];
        a.phenotype = new Phenotype[numPhenotype];

        // LOOP THROUGH THE GENES
        for (int i = 0; i < a.numGene; i++) {
            a.gene[i] = new Gene(gene[i].name, gene[i].description, a);
        }

        // LOOP THROUGH THE SAMPLES
        for (int i = 0; i < a.numSample; i++) {
            a.sample[i] = new Sample(a.sampleName.get(i), i, a);  // Make the ith sample
            a.sample[i].phenotype = sample[i].phenotype.copy();
        }

        // LOOP THROUGH THE PHENOTYPES
        for (int i = 0; i < a.numPhenotype; i++) {
            a.phenotype[i] = new Phenotype(a.phenotypeName.get(i), a);  // Make the ith phenotype

            if (phenotype[i].c != null) {
                a.phenotype[i].c = copy(phenotype[i].c);
            }

            if (phenotype[i].shape != null) {
                a.phenotype[i].shape = "" + phenotype[i].shape;
            }

            a.phenotype[i].size = phenotype[i].size;
        }

        a.x = x.copy();

        if(x.H!=null && x.W!=null){
        	a.x.H=this.x.H.copy();
        	a.x.W=this.x.W.copy();
        }
        return a;
    }

/*--------------------------------------------------------------------
 addDescription - Add descriptions
--------------------------------------------------------------------*/

    public void addDescription(TextTree m) {
        for (int i = 0; i < gene.length; i++) {
            Text z = m.get(geneName.get(i));
            if (z != null) {
                gene[i].description = z;
            }
        }
    }

/*--------------------------------------------------------------------
 load - Load the array from a gct and cls file
--------------------------------------------------------------------*/

    public void load(String gctName, String clsName) {
        loadGCT(gctName);  // Load the GCT file
        loadCLS(clsName);  // Load the CLS file
    }

/*--------------------------------------------------------------------
 load - Load the array from a gct and cls file
--------------------------------------------------------------------*/

    public void load(String gctName) {
        loadGCT(gctName);  // Load the GCT file
        loadCLS();           //No CLS file
    }

/*--------------------------------------------------------------------
 write - Write the array to a gct and cls file
--------------------------------------------------------------------*/

    public void write(String gctName, String clsName) {
        writeGCT(gctName);  // Write the GCT file
        writeCLS(clsName);  // Write the CLS file
    }

    /*--------------------------------------------------------------------
    write - Write the array to a gct and cls file
   --------------------------------------------------------------------*/

       public void write(String gctName, String clsName, boolean newWriter) {
           writeGCT(gctName, newWriter);  // Write the GCT file
           writeCLS(clsName, newWriter);  // Write the CLS file
       }
    
/*--------------------------------------------------------------------
 write - Write the array to a gct file
--------------------------------------------------------------------*/

    public void write(String gctName) {
        writeGCT(gctName);  // Write the GCT file
    }

    /*--------------------------------------------------------------------
    write - Write the array to a gct file
   --------------------------------------------------------------------*/

       public void write(String gctName, boolean newWriter) {
           writeGCT(gctName, newWriter);  // Write the GCT file
       }    
    
/*--------------------------------------------------------------------
 loadGCT - Load the array from a gct file
--------------------------------------------------------------------*/

    public void loadGCT(String fname) {
        this.id=fname;
    	Infile a = new Infile(fname);  // Open the GCT file
        a.skip(1);                     // Skip the first line
        a.readLine();                    // Read the second line
        Text[] b = a.splitTab();         // Split on tabs

        numGene = num(b, 0);    // Number of genes
        numSample = num(b, 1);  // Number of samples

        a.readLine();      // Read the next line
        b = a.splitTab();  // Split on tabs

        x = new Matrix(numGene, numSample);  // Make the expression matrix
        geneName = new Texts(numGene);       // Make a gene name array
        descriptions=new Texts(numGene);
        gene = new Gene[numGene];             // Make a gene array
        missing = new boolean[numGene];       // Make a missing value array
        sampleName = new Texts(numSample);   // Make a sample name array
        sample = new Sample[numSample];       // Make a sample array

        // LOOP THROUGH THE SAMPLES
        for (int i = 0; i < numSample; i++) {
            sample[i] = new Sample(txt(b, i + 2), i, this);  // Make the ith sample
            sampleName.set(i, sample[i].name);                  // ith sample name
        }

        // LOOP THROUGH THE GENES
        for (int i = 0; i < numGene; i++) {
            a.readLine();      // Read the next line
            b = a.splitTab();  // Split on tabs

            missing[i] = false;
            geneName.set(i, txt(b, 0));                                 // ith gene name
            descriptions.set(i, txt(b,1));
            gene[i] = new Gene(geneName.get(i), txt(b, 1), this);  // Make the ith gene

            Text zero = new Text(0);

            // LOOP THROUGH THE SAMPLES
            for (int j = 0; j < numSample; j++) {
                float val = (float) (dbl(b, j + 2));  // Read the expression for gene i in sample j

                if (val == 0 && txt(b, j + 2).find(zero) < 0) {
                    missing[i] = true;
                }

                x.x[i][j] = val;
            }
        }

        a.close();  // Close the GCT file

        phenotypeName = new Texts();
        phenotypeName.rpush("none");
    }

/*--------------------------------------------------------------------
 abbrNumber - Number the abbreviations
--------------------------------------------------------------------*/

    public void abbrNumber() {
        for (int i = 0; i < sample.length; i++) {
            sample[i].abbr = "" + (i + 1);
        }
    }

/*--------------------------------------------------------------------
 abbr - Set some abbreviations
--------------------------------------------------------------------*/

    public void abbr(String[] z) {
        for (int i = 0; i < z.length; i += 2) {
            int j = sampleName.find(z[i]);

            if (j >= 0) {
                sample[j].abbr = z[i + 1];
            }
        }
    }


/*--------------------------------------------------------------------
 writeGCT - Write the array to a gct file
--------------------------------------------------------------------*/

    public void writeGCT(String fname) {
        Outfile a = new Outfile(fname);           // Open the GCT file
        a.writeLine("#1.2");                      // Version number
        a.writeLine(numGene + "\t" + numSample);  // Version number
        a.write("Name\tDescription\t");           // Name and description

        // WRITE THE SAMPLES
        for (int i = 0; i < numSample; i++) {
            if (i > 0) {
                a.write("\t");
            }

            a.write(sample[i].makeName());
        }

        a.write("\n");

        // WRITE THE ROWS OF THE MATRIX
        for (int i = 0; i < numGene; i++) {
            a.write(geneName.get(i) + "\t");
            // a.write( gene[ i ].name + " " );
            a.write(gene[i].description);  // Write the gene name

            // WRITE EACH SAMPLE
            for (int j = 0; j < numSample; j++) {
                a.write("\t" + (double) Math.round(x.x[i][j] * 1000.0) / 1000.0);
            }

            a.write("\n");  // Go to the next line
        }

        a.close();  // Close the GCT file
    }
    
    
    /*--------------------------------------------------------------------
    writeGCT - Write the array to a gct file
   --------------------------------------------------------------------*/

       public void writeGCT(String fname, boolean newWriter) {
          try{
           System.err.println("Array's count: "+numGene+" "+numSample);
    	   FileWriter a = new FileWriter(fname);           // Open the GCT file
           a.write("#1.2\n");                      // Version number
           a.write(numGene + "\t" + numSample+"\n");  // Version number
           a.write("Name\tDescription\t");           // Name and description

           // WRITE THE SAMPLES
           for (int i = 0; i < numSample; i++) {
               if (i > 0) {
                   a.write("\t");
               }

               a.write(sample[i].makeName());
           }

           a.write("\n");

           // WRITE THE ROWS OF THE MATRIX
           for (int i = 0; i < numGene; i++) {
               a.write(geneName.get(i) + "\t");
               // a.write( gene[ i ].name + " " );
               a.write(gene[i].description.toString());  // Write the gene name

               // WRITE EACH SAMPLE
               for (int j = 0; j < numSample; j++) {
            	  // System.err.println(x.x[i][j]);
                   a.write("\t" + x.x[i][j]);
               }

               a.write("\n");  // Go to the next line
           }

           a.close();  // Close the GCT file
          }catch(IOException ex){}
       }

/*--------------------------------------------------------------------
 loadCLS - Initialize an empty cls file
--------------------------------------------------------------------*/

   public void loadCLS() {
        numPhenotype = numSample;                 // Number of phenotypes
        phenotypeName = new Texts(numPhenotype);  // Make a phenotype name array
        phenotype = new Phenotype[numPhenotype];   // Make a phenotype array

        // LOOP THROUGH THE PHENOTYPES
        for (int i = 0; i < numPhenotype; i++) {
        	phenotypeName.set(i, sample[i].name);                         // ith phenotype name
            phenotype[i] = new Phenotype(phenotypeName.get(i), this);  // Make the ith phenotype
        }

        // LOOP THROUGH THE SAMPLE
        for (int i = 0; i < numSample; i++) {
        	sample[i].addPhenotype(sample[i].name);  // Set the phenotype
        }


        setPhenotypeSize();  // Get the size of each phenotype
        colorPhenotype();
    }

/*--------------------------------------------------------------------
 loadCLS - Load the phenotype information from a cls file
--------------------------------------------------------------------*/

    public void loadCLS(String fname) {
        Infile a = new Infile(fname);  // Open the CLS file
        a.readLine();                    // Read the first line
        Text[] b = a.splitSpace();       // Split on spaces

        numPhenotype = num(b, 1);                 // Number of phenotypes
        phenotypeName = new Texts(numPhenotype);  // Make a phenotype name array
        phenotype = new Phenotype[numPhenotype];   // Make a phenotype array

        a.readLine();        // Read the second line
        b = a.splitSpace();  // Split on spaces

        // LOOP THROUGH THE PHENOTYPES
        for (int i = 0; i < numPhenotype; i++) {
        	phenotypeName.set(i, txt(b, i + 1));                         // ith phenotype name
            phenotype[i] = new Phenotype(phenotypeName.get(i), this);  // Make the ith phenotype
        }

        a.readLine();        // Read the third line
        b = a.splitSpace();  // Split on spaces

        // LOOP THROUGH THE SAMPLE
        for (int i = 0; i < numSample; i++) {
        	//System.err.println(b[i]);
        	sample[i].addPhenotype(b[i]);  // Set the phenotype
        }

        a.close();  // Close the CLS file

        setPhenotypeSize();  // Get the size of each phenotype
        colorPhenotype();
    }

/*--------------------------------------------------------------------
 setPhenotypeSize - Find the size of the phenotypes
--------------------------------------------------------------------*/

    public void setPhenotypeSize() {
        for (int i = 0; i < phenotype.length; i++) {
            int s = 0;

            for (int j = 0; j < numSample; j++) {
                if (sample[j].isPhenotype(phenotypeName.get(i))) {
                    s++;
                }
            }

            phenotype[i].size = s;
        }
    }

/*--------------------------------------------------------------------
 writeCLS - Write the array to a cls file
--------------------------------------------------------------------*/

    public void writeCLS(String fname) {
        writeCLS(fname, 0);
    }

/*--------------------------------------------------------------------
 writeCLS - Write the array to a cls file
--------------------------------------------------------------------*/

    public void writeCLS(String fname, int j) {
        Texts p = new Texts();  // Create a temporary com.scanfeld.core.Texts object

        // LOOP THROUGH THE SAMPLES
        for (int i = 0; i < numSample; i++) {
            p.rpush(sample[i].phenotype.get(j));  // Push the phenotype name
        }

        p = p.removeDuplicate();  // Remove duplicate names

        Outfile a = new Outfile(fname);                // Open the CLS file
        a.writeLine(numSample + " " + p.size + " 1");  // Number of samples and phenotypes
        a.writeLine("# " + p.join(" "));            // Phenotypes

        // LOOP THROUGH THE SAMPLES
        for (int i = 0; i < numSample; i++) {
            // IF NOT THE FIRST SAMPLE
            if (i > 0) {
                a.write(" ");  // Include a space
            }

            a.write(sample[i].phenotype.get(j));  // Write the phenotype
        }

        a.write("\n");  // Add a new line

        a.close();  // Close the CLS file
    }

    /*--------------------------------------------------------------------
    writeCLS - Write the array to a cls file
   --------------------------------------------------------------------*/

       public void writeCLS(String fname, boolean newWriter) {
    	   try{
    	   int j=0;
           Texts p = new Texts();  // Create a temporary com.scanfeld.core.Texts object

           // LOOP THROUGH THE SAMPLES
           for (int i = 0; i < numSample; i++) {
               p.rpush(sample[i].phenotype.get(j));  // Push the phenotype name
           }

           p = p.removeDuplicate();  // Remove duplicate names

           FileWriter a = new FileWriter(fname);                // Open the CLS file
           a.write(numSample + " " + p.size + " 1\n");  // Number of samples and phenotypes
           a.write("# " + p.join(" ")+"\n");            // Phenotypes

           // LOOP THROUGH THE SAMPLES
           for (int i = 0; i < numSample; i++) {
               // IF NOT THE FIRST SAMPLE
               if (i > 0) {
                   a.write(" ");  // Include a space
               }

               a.write(sample[i].phenotype.get(j).toString());  // Write the phenotype
           }

           a.write("\n");  // Add a new line

           a.close();  // Close the CLS file
    	   }catch(IOException ex){}
       }

    
/*--------------------------------------------------------------------
 writeINF - Write information about an array
--------------------------------------------------------------------*/

    public void writeINF(String fname) {
        Outfile a = new Outfile(fname);  // Open the INF file

        a.writeLine("ARRAY INFORMATION\n");                     // Section Header
        a.writeLine("Number of Genes:\t" + numGene);            // Number of genes
        a.writeLine("Number of Samples:\t" + numSample);        // Number of samples
        a.writeLine("Number of Phenotypes:\t" + numPhenotype);  // Number of phenotypes

        a.writeLine("\nSTATISTICS\n");                                            // Section Header
        a.writeLine("Index\tName\tClass\tMean\tVariance\tSD\tMin\tMax\tRange");  // Column Statistics

        Matrix mean = x.colMean();
        Matrix var = x.colVar();
        Matrix sd = x.colSD();
        Matrix min = x.colMin();
        Matrix max = x.colMax();
        Matrix range = x.colRange();
        Matrix cor = x.spearman();

        // LOOP THROUGH THE SAMPLES
        for (int i = 0; i < numSample; i++) {
            a.write((i + 1) + "\t" + sampleName.get(i));  // Sample index and name
            a.write("\t" + sample[i].phenotype.join(", "));    // Output the phenotype information
            a.write("\t" + mean.x[0][i]);                  // Mean expression
            a.write("\t" + var.x[0][i]);                   // Column variance
            a.write("\t" + sd.x[0][i]);                    // Column standard deviation
            a.write("\t" + min.x[0][i]);                   // Column minimum
            a.write("\t" + max.x[0][i]);                   // Column maximum
            a.write("\t" + range.x[0][i]);                 // Column range
            a.writeLine();                                       // End the line
        }

        a.writeLine("\nSPEARMAN CORRELATION\n");  // Section Header

        // LOOP THROUGH THE SAMPLES
        for (int i = 0; i < numSample; i++) {
            a.write("\t" + sampleName.get(i));
        }
        a.writeLine();

        // LOOP THROUGH THE SAMPLES
        for (int i = 0; i < numSample; i++) {
            a.write(sampleName.get(i));

            // LOOP THROUGH THE SAMPLES
            for (int j = 0; j < numSample; j++) {
                a.write("\t" + cor.x[i][j]);  // Correlation
            }

            a.writeLine();  // End the line
        }

        a.close();  // Close the INF file
    }

/*--------------------------------------------------------------------
 markerGene - Find the marker genes between phenotype a and b
--------------------------------------------------------------------*/

    public int[] markerGene(String a, int num) {
        int[] z = new int[num * 2];

        int[] ai = selectPhenotypeIndex(a);

        Matrix y = new Matrix(ai.length, 1);
        Matrix w = x.ntranspose();
        w = w.selectRow(ai);
        w.unit();
        y = y.mergeRight(w);

        for (int i = 0; i < ai.length; i++) {
            y.x[i][0] = 1; // + 0.00001 * randDouble ();
        }

        int sze = y.nc - 1;

        double[] dis = new double[sze];
        for (int i = 1; i < sze + 1; i++) {
            dis[i - 1] = y.euclidean(0, i);
        }

        int ind[] = sort(dis);

        for (int i = 0; i < num; i++) {
            z[i] = ind[i];
        }

        for (int i = 0; i < num; i++) {
            z[num + i] = ind[sze - num + i];
        }

        y = y.selectCol(seq(0, 10));
        return z;
    }

/*--------------------------------------------------------------------
 markerGene - Find the marker genes between phenotype a and b
--------------------------------------------------------------------*/

    public int[] markerGene(String a, String b, int num) {
        int[] z = new int[num * 2];

        int[] ai = selectPhenotypeIndex(a);
        int[] bi;

        if (b.equals("rest")) {
            bi = new int[numSample - ai.length];

            int k = 0;
            for (int i = 0; i < numSample; i++) {
                boolean found = false;

                for (int j = 0; j < ai.length; j++) {
                    if (ai[j] == i) {
                        found = true;
                        break;
                    }
                }

                if (!found) {
                    bi[k] = i;
                    k++;
                }
            }
        } else {
            bi = selectPhenotypeIndex(b);
        }

        int[] ab = c(ai, bi);

        Matrix y = new Matrix(ab.length, 1);
        Matrix w = x.ntranspose();
        w = w.selectRow(ab);
        w.unit();
        y = y.mergeRight(w);

        for (int i = 0; i < ai.length; i++) {
            y.x[i][0] = 1;
        }

        for (int i = ai.length; i < ai.length + bi.length; i++) {
            y.x[i][0] = 0;
        }

        int sze = y.nc - 1;

        double[] dis = new double[sze];
        for (int i = 1; i < sze + 1; i++) {
            dis[i - 1] = y.pearson(0, i);
        }

        int ind[] = rev(sort(dis));
        dis = rev(dis);

        for (int i = 0; i < num; i++) {
            z[i] = ind[i];
        }

        for (int i = 0; i < num; i++) {
            z[num + i] = ind[sze - num + i];
        }

        return z;
    }

/*--------------------------------------------------------------------
 writeTopGene - Write the top genes per sample
--------------------------------------------------------------------*/

    public void writeTopGene(String fname, int top) {
        Outfile a = new Outfile(fname);  // Open the INF file

        for (int i = 0; i < numSample; i++) {
            Texts set = new Texts();
            int ord[] = new int[numGene];

            for (int j = 0; j < numGene; j++) {
                for (int k = 0; k < numGene; k++) {
                    if (x.x[j][i] > x.x[k][i]) {
                        ord[j]++;
                    }
                }
            }

            for (int j = numGene - 1; j > numGene - 100; j--) {
                for (int k = 0; k < numGene; k++) {
                    if (ord[k] == j) {
                        set.rpush(geneName.get(k));
                    }
                }
            }

            set = set.removeDuplicate();
            a.write(sampleName.get(i) + "\tBLACK\t" + set.join("\t") + "\n");
        }

        a.close();
    }

/*--------------------------------------------------------------------
 threshold - Threshold the genes
--------------------------------------------------------------------*/

    public void threshold(double min, double max) {
        x.threshold(min, max);
    }

/*--------------------------------------------------------------------
 variation - Variation filter
--------------------------------------------------------------------*/

    public void variation(double fold, double delta) {
        int ind[] = x.variation(fold, delta);  // Find the indices of the satisfying genes;
        selectGene(ind);                       // Select the genes
    }

/*--------------------------------------------------------------------
 rank - Rank the genes
--------------------------------------------------------------------*/

    public void rank() {
        x = x.rank();  // Create a rank matrix
    }

/*--------------------------------------------------------------------
 rank - Rank the genes
--------------------------------------------------------------------*/

    public void rank(double max) {
        x = x.rank(max);  // Create a rank matrix
    }

    
    /*--------------------------------------------------------------------
    max - get the max value in matrix
   --------------------------------------------------------------------*/

       public double HMax() {
    	   double max=-Double.MAX_VALUE;
           for(int i=0; i<x.H.nr; i++){
        	   for(int j=0; j<x.H.nc; j++){
        		   max=Math.max(max, x.H.x[i][j]);
        	   }
           }
           return max;
       }
    
       public double[] HMaxByMetagene() {
    	   double[] max=new double[x.H.nr];
           for(int i=0; i<x.H.nr; i++){
        	   max[i]=-Double.MAX_VALUE;
        	   for(int j=0; j<x.H.nc; j++){
        		   max[i]=Math.max(max[i], x.H.x[i][j]);
        	   }
           }
           return max;
       }     
       
  public Array assignPValues(EmpiricalDistribution dist){
	  Array rtrn=this.copy();
	  
	  for(int i=0; i<rtrn.x.H.nr; i++){
		  for(int j=0; j<rtrn.x.H.nc; j++){
			  double p=1-dist.getCummulativeProbability(x.H.x[i][j]);
			  rtrn.x.H.x[i][j]=p;
		  }
	  }
	  
	  return rtrn;
  }     
  
    
  
  
  public Array assignPValues(EmpiricalDistribution[] distsByMetagene){
	  Array rtrn=this.copy();
	  
	  for(int i=0; i<rtrn.x.H.nr; i++){
		  EmpiricalDistribution dist=distsByMetagene[i];
		  for(int j=0; j<rtrn.x.H.nc; j++){
			  double p=1-dist.getCummulativeProbability(x.H.x[i][j]);
			  rtrn.x.H.x[i][j]=p;
		  }
	  }
	  
	  return rtrn;
  }     
  
  
  public Array assignPValues(ArrayList[] distsByMetagene){
	  Array rtrn=this.copy();
	  
	  for(int i=0; i<rtrn.x.H.nr; i++){
		  ArrayList dist=distsByMetagene[i];
		  for(int j=0; j<rtrn.x.H.nc; j++){
			  double p=1-getCummulativeProbability(dist, x.H.x[i][j]);
			  rtrn.x.H.x[i][j]=p;
		  }
	  }
	  
	  return rtrn;
  }     
  
  private double getCummulativeProbability(ArrayList<Double> vals, double val){
	  double counterTotal=0;
	  double countLessThan=0;
	  
	  for(Double observed: vals){
		  if(observed<val){countLessThan++;}
		  counterTotal++;
	  }
	  
	  return countLessThan/counterTotal;
  }
  
/*--------------------------------------------------------------------
 rank - Rank the genes
--------------------------------------------------------------------*/

    public void rankGene() {
        x = x.ntranspose().rank().ntranspose();  // Create a rank matrix
    }

/*--------------------------------------------------------------------
 rank - Rank the genes
--------------------------------------------------------------------*/

    public void rankGene(double max) {
        x = x.ntranspose().rank(max).ntranspose();  // Create a rank matrix
    }

/*--------------------------------------------------------------------
 minPhenotype - Get the number of samples in the smallest phenotype
--------------------------------------------------------------------*/

    public int minPhenotype() {
        int min = selectPhenotypeIndex(0).length;

        for (int i = 1; i < phenotype.length; i++) {
            int tsze = selectPhenotypeIndex(i).length;

            if (tsze < min) {
                min = tsze;
            }
        }

        return min;
    }

/*--------------------------------------------------------------------
 getPhenotype - Get the sample indices for the phenotype
--------------------------------------------------------------------*/

    public int[] getPhenotype(String s) {
        return getPhenotype(new Text(s));
    }

/*--------------------------------------------------------------------
 getPhenotype - Get the sample indices for the phenotype
--------------------------------------------------------------------*/

    public int[] getPhenotype(Text s) {
        Texts num = new Texts();

        for (int i = 0; i < numSample; i++) {
            // IF THE SAMPLE HAS THE PHENOTYPE
            if (sample[i].isPhenotype(s)) {
                num.rpush(i);
            }
        }

        return num.num();
    }

/*--------------------------------------------------------------------
 selectGeneObject - Select the genes
--------------------------------------------------------------------*/

    public void selectGeneObject(int[] si) {
        Gene[] ns = new Gene[si.length];  // Create a new gene array

        // LOOP THROUGH THE GENE ARRAY
        for (int i = 0; i < si.length; i++) {
            ns[i] = gene[si[i]];  // Set the ith gene
        }

        gene = ns;               // Set to the new gene array
        x = x.selectRow(si);  // Set to the new columns
    }

/*--------------------------------------------------------------------
 selectGene - Select a set of genes
--------------------------------------------------------------------*/

    public void selectGene(Texts z) {
        int[] si = geneName.select(z);  // Get the gene indices
        selectGene(si);                 // Select the genes
    }

/*--------------------------------------------------------------------
 selectGene - Select a set of genes
--------------------------------------------------------------------*/

    public void selectGene(Text[] z) {
        int[] si = geneName.select(z);  // Get the gene indices
        selectGene(si);                 // Select the genes
    }

/*--------------------------------------------------------------------
 selectGene - Select a set of genes
--------------------------------------------------------------------*/

    public void selectGene(String[] z) {
        int[] si = geneName.select(z);  // Get the gene indices
        selectGene(si);                 // Select the genes
    }

/*--------------------------------------------------------------------
 selectGene - Select a set of genes
--------------------------------------------------------------------*/

    public void selectGene(boolean[] si) {
        Texts ind = new Texts();
        for (int i = 0; i < si.length; i++) {
            if (si[i]) {
                ind.rpush(i);
            }
        }

        selectGene(ind.num());
    }

/*--------------------------------------------------------------------
 selectGene - Select a set of genes
--------------------------------------------------------------------*/

    public void selectGene(int[] si) {
        geneName = geneName.select(si);  // Select the gene names
        selectGeneObject(si);            // Select the gene objects
        numGene = geneName.size;            // Update the number of genes
    }

/*--------------------------------------------------------------------
 geneTree - Return a gene tree
--------------------------------------------------------------------*/

    public TextTree geneTree() {
        if (geneTree != null) {
            return geneTree;
        } else {
            geneTree = new TextTree();

            for (int i = 0; i < geneName.size; i++) {
                geneTree.set(geneName.get(i), gene[i].description);
            }

            return geneTree;
        }
    }

/*--------------------------------------------------------------------
 merge - Merge with another dataset
--------------------------------------------------------------------*/

    public void merge(Array b) {
        merge(b,true);
    }

/*--------------------------------------------------------------------
 merge - Merge with another dataset
--------------------------------------------------------------------*/

    public void merge(Array b, boolean removeFlag) {
        Texts ga = geneName;
        Texts gb = b.geneName;

        if (removeFlag) {
            ga = ga.match(gb);
            gb = gb.match(ga);

            selectGene(ga);
            b.selectGene(gb);
            x = x.mergeRight(b.x);

        } else {
            Texts genes = new Texts(ga);
            genes.rpush(gb);
            genes.sort();
            genes = genes.removeDuplicate();

            Tree t = new Tree();
            for (int i = 0; i < geneName.size; i++) {
                t.set(geneName.get(i), 1);
            }

            Tree ft = new Tree();
            for (int i = 0; i < genes.size; i++) {
                ft.set(genes.get(i), i + 1);
            }

            Matrix y = new Matrix(genes.size, numSample + b.numSample);

            for (int i = 0; i < numGene; i++) {
                int gind = ft.get(geneName.get(i)) - 1;
                if (gind > 0) {
                    for (int j = 0; j < numSample; j++) {
                        y.x[gind][j] = x.x[i][j];
                    }
                }
            }

            for (int i = 0; i < b.numGene; i++) {
                int gind = ft.get(b.geneName.get(i)) - 1;
                if (gind > 0) {
                    for (int j = 0; j < b.numSample; j++) {
                        y.x[gind][numSample + j] = b.x.x[i][j];
                    }
                }
            }

            geneName = genes;
            x = y;
        }

        sampleName.rpush(b.sampleName);

        Sample[] newsample = new Sample[numSample + b.numSample];

        for (int i = 0; i < numSample; i++) {
            newsample[i] = sample[i];
            sample[i].array = this;
        }

        for (int i = 0; i < b.numSample; i++) {
            newsample[numSample + i] = b.sample[i];
            newsample[numSample + i].array = this;
        }

        sample = newsample;
        numSample += b.numSample;

        phenotypeName.rpush(b.phenotypeName);

        Phenotype[] p = new Phenotype[phenotype.length + b.phenotype.length];

        for (int i = 0; i < phenotype.length; i++) {
            p[i] = phenotype[i];
        }

        for (int i = 0; i < b.phenotype.length; i++) {
            p[phenotype.length + i] = b.phenotype[i];
            p[phenotype.length + i].array = this;
        }

        numPhenotype += b.numPhenotype;
        phenotype = p;
    }

/*--------------------------------------------------------------------
 clusterHierarchical - A hierarchical clustering of the samples
--------------------------------------------------------------------*/

    public void clusterHierarchical() {
        clusterHierarchical("pearson", "complete");
    }

/*--------------------------------------------------------------------
 clusterHierarchical - A hierarchical clustering of the samples
--------------------------------------------------------------------*/

    public void clusterHierarchical(String metric, String type) {
        shc = x.clusterHierarchical(metric, type);
        selectSample(shc.order);
    }

/*--------------------------------------------------------------------
 clusterHierarchicalGene - A hierarchical clustering of the genes
--------------------------------------------------------------------*/

    public void clusterHierarchicalGene() {
        clusterHierarchicalGene("pearson", "complete");
    }

/*--------------------------------------------------------------------
 clusterHierarchicalGene - A hierarchical clustering of the genes
--------------------------------------------------------------------*/

    public void clusterHierarchicalGene(String metric, String type) {
        Matrix xc = x.copy();
        xc.transpose();
        ghc = xc.clusterHierarchical(metric, type);
        selectGene(ghc.order);
    }

/*--------------------------------------------------------------------
 phenotypeShape - Set a phenotype shape
--------------------------------------------------------------------*/

    public void phenotypeShape(String n, String s) {
        for (int i = 0; i < numPhenotype; i++) {
            if (phenotypeName.get(i).eq(n)) {
                phenotype[i].shape = s;
                break;
            }
        }
    }

/*--------------------------------------------------------------------
 colorPhenotypeRGB - Color a phenotype
--------------------------------------------------------------------*/

    public void colorPhenotypeRGB(String n, double r, double g, double b) {
        for (int i = 0; i < numPhenotype; i++) {
            if (phenotypeName.get(i).eq(n)) {
                phenotype[i].c = rgb( r, g, b );
                break;
            }
        }
    }

/*--------------------------------------------------------------------
 colorPhenotypeHSV - Color a phenotype
--------------------------------------------------------------------*/

    public void colorPhenotypeHSV(String n, double h, double s, double v) {
        for (int i = 0; i < numPhenotype; i++) {
            if (phenotypeName.get(i).eq(n)) {
                phenotype[i].c = hsv( h, s, v );
                break;
            }
        }
    }

/*--------------------------------------------------------------------
 colorPhenotype - Color a phenotype
--------------------------------------------------------------------*/

    public void colorPhenotype(String n, String h) {
        for (int i = 0; i < numPhenotype; i++) {
            if (phenotypeName.get(i).eq(n)) {
                phenotype[i].c = rgb( h );
                break;
            }
        }
    }

/*--------------------------------------------------------------------
 colorPhenotype - Color the phenotypes
--------------------------------------------------------------------*/

    public void colorPhenotype() {
        ColorBlend pb = new ColorBlend(1.0, 0.0, 0.0);
        pb.add(1.0, 0.5, 0.0, 1.0 / 6.0);
        pb.add(1.0, 1.0, 0.0, 2.0 / 6.0);
        pb.add(0.0, 1.0, 0.0, 3.0 / 6.0);
        pb.add(0.0, 0.0, 1.0, 4.0 / 6.0);
        pb.add(1.0, 0.0, 1.0, 5.0 / 6.0);

        for (int i = 0; i < numPhenotype; i++) {
            phenotype[i].c = pb.get((double) (i) / (double) (numPhenotype));
        }
    }

/*--------------------------------------------------------------------
 colorPhenotype - Color the phenotypes
--------------------------------------------------------------------*/

    public void colorPhenotype(double col[][]) {
        for (int i = 0; i < numPhenotype; i++) {
            phenotype[i].c = col[i];
        }
    }

/*--------------------------------------------------------------------
 getPhenotype - Get the phenotype for each sample
--------------------------------------------------------------------*/

    public int[] getPhenotype() {
        int[] p = new int[numSample];

        for (int i = 0; i < numSample; i++) {
            p[i] = getPhenotype(i);
        }

        return p;
    }

/*--------------------------------------------------------------------
 getPhenotype - Get the phenotype for a sample
--------------------------------------------------------------------*/

    public int getPhenotype(int k) {
       for (int i = 0; i < numPhenotype; i++) {
            if (sample[k].isPhenotype(phenotypeName.get(i))) {
                return i;
            }
        }

        return -1;
    }

/*--------------------------------------------------------------------
 getPhenotypeShape - Get the shape for the phenotype
--------------------------------------------------------------------*/

    public String getPhenotypeShape(int k)

    {
        return phenotype[k].shape;
    }

/*--------------------------------------------------------------------
 getSampleShape- Get the shape for the sample
--------------------------------------------------------------------*/

    public String getSampleShape(int k) {
        String phenotypeShape = getPhenotypeShape(getPhenotype(k));
        String sampleShape = sample[k].shape;

        if (sampleShape != null) {
            return sampleShape;
        }

        if (phenotypeShape != null) {
            return phenotypeShape;
        }

        return "";
    }

/*--------------------------------------------------------------------
 radial - Create a radial plot
--------------------------------------------------------------------*/

   /* public void radial(String title, String fname) {
        double[] xdis = x.euclidean();
        double maxdis = max(xdis);

        double[] cbc = new double[numPhenotype];
        cbc = add(cbc, 0.0);

        colorPhenotype();
        radial(title, fname, maxdis, cbc);
    }*/

/*--------------------------------------------------------------------
 radial - Create a radial plot
--------------------------------------------------------------------*/

    /*public void radial(String title, String fname, double maxdis, double[] cbc) {
        radial(title, fname, maxdis, cbc, 0, "euclidean", null);
    }*/

/*--------------------------------------------------------------------
 radial - Create a radial plot
--------------------------------------------------------------------*/

    /*public void radial(String title, String fname, double maxdis, double cbc[], int type, String zOrder, Array b) {
        if (type == 0 || type == 1) {
            curPDF = new PDF(fname);
        }

        PDF p = curPDF;

        if (showHeatMap) {
            p.title(title);
            p.push(p.ys, 0, p.xs, p.ys);
            heatmap(p, seq(0, numGene - 1), true, "redblue");
            // heatmap ( p, seq ( 0, numSample - 1 ), true, "redblue" );
            p.pop();
            p.push(0, 0, p.ys, p.ys);
        } else {
            p.squareTitle(title);
        }
        // p.push ( 10, 10, -10, -10 );

        int k = numGene;
        double kc[][] = new double[k][2];
        double sc[][] = new double[numSample][2];

        double[] ang = new double[k];
        double[] spc = new double[k];

        boolean spreadFlag = false;
        int[] order = seq(0, k - 1);

        Matrix y = x.ntranspose();
        if (spreadFlag) {
            ClusterHierarchical ch = y.clusterHierarchical("pearson", "complete");
            order = ch.order;
            x = x.selectRow(order);
        }

        for (int i = 0; i < k; i++) {
            int j = i + 1;

            if (j == k) {
                j = 0;
            }

            if (spreadFlag) {
                spc[i] = (y.pearson(i, j) + 1) / 2;
            } else {
                spc[i] = 1.0 / (double) (k);
            }
        }

        double spcsum = sum(spc);
        ang[0] = 0;
        for (int i = 1; i < k; i++) {
            ang[i] = ang[i - 1] + spc[i] / spcsum;
        }

        for (int i = 0; i < k; i++) {
            double r = ang[i];
            kc[i][0] = Math.cos(r * 2.0 * Math.PI);
            kc[i][1] = Math.sin(r * 2.0 * Math.PI);
        }

        double d = 35.0 / (double) (p.xs);

        p.beginText();
        p.stroke(0.0, 0.0, 0.0);
        p.fill(0.0, 0.0, 0.0);
        p.setFontSize(14);

        for (int i = 0; i < k; i++) {
            p.centerText("" + (order[i] + 1), n2p(kc[i][0] * (1 - d)), n2p(kc[i][1] * (1 - d)));
        }

        p.push(d, d, -d, -d);

        ColorBlend c = new ColorBlend(0.95, 0.95, 0.95);
        c.add(1.0, 1.0, 1.0, 1.0);

        for (int i = 0; i < 10; i++) {
            double r = (double) (i) / 10.0;

            if (i == 0) {
                p.stroke(0.5, 0.5, 0.5);
            } else {
                p.stroke(0.75, 0.75, 0.75);
            }

            p.fill(c.get(r));
            p.circle(0.5, 0.5, 0.5 - r * 0.5);
        }

        p.stroke(0.5, 0.5, 0.5);

        for (int i = 0; i < k; i++) {
            p.moveTo(0.5, 0.5);
            p.lineTo(n2p(kc[i][0]), n2p(kc[i][1]));
        }

        p.stroke();

        // draw  circle for the actual mean

        double mn = calcRadial(kc, sc, maxdis, 1.0);  // ASDF - Write code to calculate the proper e.

       

        double pwid = 7.5 / (double) (p.xs);

        double br[];
        if (b == null) {
            b = this;
        }

        if (zOrder.equals("brier")) {
            br = b.brier();
        } else {
            br = new Matrix(sc).ntranspose().euclidean();
        }

        int[] ind = sort(copy(br));
        int[] maxind = b.colMaxInd();

        ColorBlend[] pc = new ColorBlend[b.numPhenotype];

        for (int i = 0; i < b.numPhenotype; i++) {
            double tcbc = cbc[i];

            pc[i] = new ColorBlend(b.phenotype[i].c);

           

            if (!shadeToGray) {
                pc[i].add(b.phenotype[i].c, (1 - tcbc));
            }

            pc[i].add(0.75, 0.75, 0.75, (1 - tcbc) + 0.0001);
            pc[i].add(0.75, 0.75, 0.75, 1.0);
        }

        // String abbr[] = new String[] { "1", "2", "3" };

        for (int i = 0; i < numSample; i++) {
            int j = ind[i];
            // int ph = maxind[ j ]; // ASDF getPhenotype is off here, so switched to maxind
            int ph = getPhenotype(j); // ASDF getPhenotype is off here, so switched to maxind

            // p.stroke ( 0.1, 0.1, 0.1 );

            double[] col = pc[ph].get(1 - br[j]);
            // col = phenotype[ getPhenotype ( j ) ].c;

            double tcbc = cbc[ph];

            if (showBelowCutoff) {
                tcbc = 0.0;
            }

            // if ( br[ j ] >= tcbc && b.getPhenotype ( j ) >= 0 )
            // {

            double sh = 1 - pwid * 3;

            if (br[j] >= tcbc && b.getPhenotype(j) >= 0) {
                if (sample[j].abbr != null && sample[j].abbr.length() > 0) {
                    p.stroke(mul(col, 0.5));
                    p.fill(col);

                    String shape = b.getSampleShape(j);

                    if (shape.equals("circle")) {
                        p.circle(n2p(sc[j][0] * sh), n2p(sc[j][1] * sh), pwid);
                    } else {
                        p.rect(n2p(sc[j][0] * sh) - pwid, n2p(sc[j][1] * sh) - pwid, n2p(sc[j][0] * sh) + pwid, n2p(sc[j][1] * sh) + pwid);
                    }

                    if (sample[j].abbr != null && sample[j].abbr.length() > 0) {
                        p.beginText();
                        p.stroke(0.0, 0.0, 0.0);
                        p.fill(0.0, 0.0, 0.0);
                        p.setFontSize(12);
                        p.centerText(sample[j].abbr, n2p(sc[j][0] * sh), n2p(sc[j][1] * sh));
                    }
                }
            }

           
        }

        for (int i = 0; i < 1; i++) {
            double r = (double) (i) / 10.0;

            if (i == 0) {
                p.stroke(0.5, 0.5, 0.5);
            } else {
                p.stroke(0.75, 0.75, 0.75);
            }

            p.circleLine(0.5, 0.5, 0.5 - r * 0.5);
        }

        p.stroke(0.5, 0.5, 0.5);

        if (type == 0 || type == 3) {
            p.close();
        } else {
            p.next();
        }
    }*/

/*--------------------------------------------------------------------
 writeDoubleHC - Write a double hc plot
--------------------------------------------------------------------*/

    /**
     * Creates a double hierarchical cluster plot
     * <a href="http://www.google.com">google</a>
     *
     * @param title
     * @param fname
     * @param b
     */

    /*
    public void writeDoubleHC(String title, String fname, Array b) {
        // clusterHierarchical ();
        clusterHierarchicalGene();
        b.clusterHierarchicalGene();

        PDF p = new PDF(fname + ".pdf");
        p.showBackground = true;
        p.squareTitle(title);
        p.push(10, 10, -10, -10);

        double hcWidth = (50.0 / p.xs);
        double spaceWidth = (5.0 / p.xs);
        double hcHeight = (50.0 / p.ys);
        double spaceHeight = (5.0 / p.ys);

        double squareWidth = (1.0 - hcWidth - spaceWidth) / (double) (b.numGene + numSample);
        double squareHeight = (1.0 - hcHeight - spaceHeight) / (double) (numGene + b.numSample);

        double awid = squareWidth * numSample;
        double bhei = squareHeight * numSample;

        Texts bs = b.sampleName;
        Texts bg = b.geneName;
        b.geneName = bs;
        b.sampleName = bg;
        b.numGene = bs.size;
        b.numSample = bg.size;
        b.x = b.x.ntranspose();
        b.shc = b.ghc;
        b.ghc = null;

        p.push(0, 0, hcWidth + awid, 1 - hcHeight - bhei - spaceHeight);
        heatmap(p, seq(0, numGene - 1), true, "redblue", 50, true);
        p.pop();

        p.push(hcWidth + awid + spaceWidth, 1 - hcHeight - bhei, 1, 1);
        b.heatmap(p, seq(0, b.numGene - 1), true, "redblue", 50, false);
        p.pop();

        ColorBlend cb = new ColorBlend("whiteblack");

        p.push(hcWidth + awid + spaceWidth, 0, 1, 1 - hcHeight - bhei - spaceHeight);
        p.lineWidth(1.0);
        p.stroke(0, 0, 0);
        p.fill(cb.get(0));
        p.rect(0, 0, 1.0, 1.0);

        int ycut1 = 0;
        int ycut2 = 100;  // 498
        int zcut1 = 0;
        int zcut2 = 100;  // 498

        ClusterHierarchical y = ghc;
        int yNumSet = y.setd.length - 1;
        int[] yord = order2rank(y.order);

        ClusterHierarchical z = b.shc;
        int zNumSet = z.setd.length - 1;
        int[] zord = order2rank(z.order);

        double[] y1 = new double[yNumSet];
        double[] y2 = new double[yNumSet];

        for (int i = 0; i < y.nc; i++) {
            y1[i] = (double) (yord[i]) / (double) (y.nc);
            y2[i] = (double) (yord[i] + 1) / (double) (y.nc);
        }

        for (int i = y.nc; i < yNumSet; i++) {
            y1[i] = y1[y.seta[i]];
            y2[i] = y2[y.setb[i]];
        }

        double[] z1 = new double[zNumSet];
        double[] z2 = new double[zNumSet];

        for (int i = 0; i < z.nc; i++) {
            z1[i] = (double) (zord[i]) / (double) (z.nc);
            z2[i] = (double) (zord[i] + 1) / (double) (z.nc);
        }

        for (int i = z.nc; i < zNumSet; i++) {
            z1[i] = z1[z.seta[i]];
            z2[i] = z2[z.setb[i]];
        }

        Texts t = new Texts();
        HCEnrich(t, y, y.setd.length - 2, y.setd.length - 2 - ycut1, y.setd.length - 2 - ycut2, z, z.setd.length - 2, z.setd.length - 2 - zcut1, z.setd.length - 2 - zcut2, numGene, 1.0, 0.1);

        int[] bi = new int[t.size];
        int[] bj = new int[t.size];
        double[] be = new double[t.size];

        for (int i = 0; i < t.size; i++) {
            Text[] x = t.get(i).splitSpace();
            bi[i] = x[0].num();
            bj[i] = x[1].num();
            be[i] = x[2].dbl();
        }

        int ind[] = sort(be);
        ind = rev(ind);
        be = rev(be);
        bi = sel(bi, ind);
        bj = sel(bj, ind);

        double bemin = 0; // -log ( max ( be ) );
        double bemax = -log(min(be));

        for (int k = 0; k < t.size; k++) {
            int i = bi[k];
            int j = bj[k];
            double e = (-log(be[k]) - bemin) / (bemax - bemin);

            if (e >= 0.0) {
                p.lineWidth(0.0);
                p.stroke(cb.get(e));
                p.fill(cb.get(e));
                p.rect(z1[j], 1.0 - y1[i], z2[j], 1.0 - y2[i]);
            }
        }

        p.lineWidth(1.0);
        p.stroke(0, 0, 0);
        p.rect();

        

        p.pop();

        p.close();
    }*/

/*--------------------------------------------------------------------
 HCEnrich - Show HC enrichment
--------------------------------------------------------------------*/

    public void HCEnrich(Texts p, ClusterHierarchical y, int i, int ycut1, int ycut2, ClusterHierarchical z, int j, int zcut1, int zcut2, int size, double e, double min) {
        double d = 1;

        if (i <= ycut1 && j <= zcut1) {
            Texts ym = y.sort(i);
            ym.lpop(2);
            Texts zm = z.sort(j);
            zm.lpop(2);
            Texts mm = ym.match(zm);

            int yct = ym.size;
            int zct = zm.size;
            int mct = mm.size;

            if (yct <= zct) {
                d = hyperGeometric(mct, yct, zct, size);
            } else {
                d = hyperGeometric(mct, zct, yct, size);
            }
        }

        if (d < e && d < min) {
            p.rpush(i + " " + j + " " + d);
            e = d;
        }

        if (z.seta[j] >= zcut2) {
            HCEnrich(p, y, i, ycut1, ycut2, z, z.seta[j], zcut1, zcut2, size, e, min);
        }

        if (z.setb[j] >= zcut2) {
            HCEnrich(p, y, i, ycut1, ycut2, z, z.setb[j], zcut1, zcut2, size, e, min);
        }

        if (j == z.setd.length - 2) {
            if (y.seta[i] >= ycut2) {
                HCEnrich(p, y, y.seta[i], ycut1, ycut2, z, j, zcut1, zcut2, size, e, min);
            }

            if (y.setb[i] >= ycut2) {
                HCEnrich(p, y, y.setb[i], ycut1, ycut2, z, j, zcut1, zcut2, size, e, min);
            }
        }
    }

/*--------------------------------------------------------------------
 HCEnrich2 - Show HC enrichment
--------------------------------------------------------------------*/

    public void HCEnrich2(Texts p, ClusterHierarchical y, int i, int ycut1, int ycut2, ClusterHierarchical z, int j, int zcut1, int zcut2, int size, double e) {
        double d = 1;

        if (i <= ycut1 && j <= zcut1) {
            Texts ym = y.sort(i);
            ym.lpop(2);
            Texts zm = z.sort(j);
            zm.lpop(2);
            Texts mm = ym.match(zm);

            int yct = ym.size;
            int zct = zm.size;
            int mct = mm.size;

            if (yct <= zct) {
                d = hyperGeometric(mct, yct, zct, size);
            } else {
                d = hyperGeometric(mct, zct, yct, size);
            }
        }

        if (d < e) {
            p.rpush(i + " " + j + " " + e);
            e = d;
        }

        if (z.seta[j] >= zcut2) {
            HCEnrich2(p, y, i, ycut1, ycut2, z, z.seta[j], zcut1, zcut2, size, e);
        }

        if (z.setb[j] >= zcut2) {
            HCEnrich2(p, y, i, ycut1, ycut2, z, z.setb[j], zcut1, zcut2, size, e);
        }

        if (j == z.setd.length - 2) {
            if (y.seta[i] >= ycut2) {
                HCEnrich2(p, y, y.seta[i], ycut1, ycut2, z, j, zcut1, zcut2, size, e);
            }

            if (y.setb[i] >= ycut2) {
                HCEnrich2(p, y, y.setb[i], ycut1, ycut2, z, j, zcut1, zcut2, size, e);
            }
        }
    }

/*--------------------------------------------------------------------
 writeWH - Write out the heatmaps for W, H, and this matrix
--------------------------------------------------------------------*/

    /*public void writeWH(String title, String fname, Array w, Array h, boolean sampleFlag, boolean geneFlag) {
        writeWH(title, fname, w, h, null, sampleFlag, geneFlag);
    }*/

/*--------------------------------------------------------------------
 writeWH - Write out the heatmaps for W, H, and this matrix
--------------------------------------------------------------------*/

    /*public void writeWH(String title, String fname, Array w, Array h, Array svm, boolean sampleFlag, boolean geneFlag) {
        Array a = copy();

        w.selectGene(a.geneName);
        a.selectGene(w.geneName);

        if (geneFlag) {
            a.clusterHierarchicalGene();
        }

        w.selectGene(a.geneName);

        if (sampleFlag) {
            a.sortByPhenotypeHierarchical();
        }

        // a.clusterHierarchical ();
        h.selectSample(a.sampleName);

        if (svm != null) {
            svm.selectSample(a.sampleName);
        }

        a.ghc = null;
        a.shc = null;

        PDF p = new PDF(fname + ".wh.pdf");
        p.title(title);
        h.write(fname + ".wh.H.gct");
        w.write(fname + ".wh.W.gct");
        a.write(fname + ".wh.gct", fname + ".wh.cls");
        a.writeWHPDF(p, w, h, svm, true, "redblue");
        p.close();

        Matrix a2 = a.x.copy();
        a.x = w.x.nmatmul(h.x);

        p = new PDF(fname + ".wh.m.pdf");
        p.title(title + " WxH");
        a.write(fname + ".wh.A.gct", fname + ".wh.A.cls");
        a.writeWHPDF(p, w, h, svm, true, "redblue");
        p.close();

        a.x.sub(a2);
        a.x.abs();
        w.selectSample(new int[]{0});
        h.selectGene(new int[]{0});

        w.x = a.x.rowMean();
        w.x.transpose();
        h.x = a.x.colMean();

        p = new PDF(fname + ".wh.d.pdf");
        p.title(title + " WxH Diff");
        a.write(fname + ".wh.d.gct", fname + ".wh.d.cls");
        a.writeWHPDF(p, w, h, svm, false, "whiteblack");
        p.close();
    }*/

/*--------------------------------------------------------------------
 writeWHPDF - Write out the heatmaps for W, H, and this matrix
--------------------------------------------------------------------*/

    /*public void writeWHPDF(PDF p, Array w, Array h) {
        writeWHPDF(p, w, h, null, true, "redblue");
    }*/

/*--------------------------------------------------------------------
 writeWHPDF - Write out the heatmaps for W, H, and this matrix
--------------------------------------------------------------------*/

    /*public void writeWHPDF(PDF p, Array w, Array h, Array svm, boolean norm, String color) {
        p.push(10, 10, -10, -10);
        double totalSample = w.numSample + numSample;
        double wwid = (double) (w.numSample) / totalSample;

        double hhei = 10 * h.numGene;
        double hhei2 = 0;
        double hhei3 = 5;

        if (svm != null) {
            hhei2 += 10 * svm.numGene + hhei3;
            hhei += 10 * svm.numGene + hhei3;
        }

        p.push(wwid, 0, 1, 1);  // Leave space for the W matrix
        p.push(10, -30, 0, 0);   // The top rectangle

        p.setFontSize(p.textFontSize);

        double lp = 0;
        double rp = 0;

        // LOOP THROUGH THE PHENOTYPES
        for (int i = 0; i < numPhenotype; i++) {
            // FIND THE COORDINATES OF THE PHENOTYPE BAR
            lp = rp;
            rp = lp + (double) (phenotype[i].size) / (double) (numSample);

            p.push(lp, 0, rp, 1);  // Focus on the phenotype width
            p.push(0, 25, 0, 30);  // The top 5 pixels of the header
            p.stroke(0, 0, 0);
            p.fill(0, 0, 0);
            p.centerText(phenotypeName.get(i).str(), 0.5, 0);  // Output the phenotype name
            p.pop();

            p.push(0, 5, 0, 12);  // The color strip for the phenotype
            p.stroke(0, 0, 0);
            p.fill(phenotype[i].c);
            p.rect(0, 0, 1.0, 1.0);
            p.pop();
            p.pop();
        }

        p.pop();
        p.pop();

        p.push(0, 0, 0, -30);  // The rest of the plot

        p.push(0, 0, wwid, 1.0 - hhei);
        p.push(0, 0, 0, -10);
        w.heatmap(p, seq(0, w.numGene - 1), norm, color);
        p.pop();
        p.pop();

        if (svm != null) {
            p.push(wwid, 1.0 - hhei2 + hhei3, 1.0, 1.0);
            p.push(10, 0, 0, 0);
            svm.heatmap(p, seq(0, svm.numGene - 1), false, "whiteblack");
            p.pop();
            p.pop();
        }

        p.push(wwid, 1.0 - hhei, 1.0, 1.0 - hhei2);
        p.push(10, 0, 0, 0);
        h.heatmap(p, seq(0, h.numGene - 1), norm, color);
        p.pop();
        p.pop();

        p.push(wwid, 0, 1.0, 1.0 - hhei);
        p.push(10, 0, 0, -10);
        heatmap(p, seq(0, numGene - 1), norm, color);
        p.pop();
        p.pop();

        p.pop();
    }*/

/*--------------------------------------------------------------------
 maxFactor - Find the maximum factor for each phenotype
--------------------------------------------------------------------*/

    public int[] maxFactor() {
        Matrix th = x.copy();
        int numFactor = th.nr;
        int[] cm = new int[numFactor];
        int[] ct = new int[numPhenotype];
        int mu = numFactor / numPhenotype;

        for (int i = 0; i < numFactor; i++) {
            int[] cmi = d2i(th.colMaxInd().getRow(0));
            int[] cmc = new int[numFactor * numPhenotype];

            for (int j = 0; j < cmi.length; j++) {
                cmc[getPhenotype(j) * numFactor + cmi[j]]++;
            }

            int[] ind = rev(sort(cmc));
            int cmin = min(ct);

            for (int j = 0; j < ind.length; j++) {
                int tc = ind[j] % numFactor;
                int tr = (ind[j] - tc) / numFactor;

                if (ct[tr] == cmin) {
                    cm[tr * mu + ct[tr]] = tc;
                    ct[tr]++;
                    th.zeroRow(tc);
                    break;
                }
            }
        }

        return cm;
    }

/*--------------------------------------------------------------------
 maxFactorH - Find the maximum factor for each phenotype
--------------------------------------------------------------------*/

    public int[] maxFactorH() {
        Matrix th = x.H.copy(); //original H matrix
        int numFactor = th.nr; //number of metagenes
        int[] cm = new int[numFactor]; // array with metagene columns
        int[] ct = new int[numPhenotype]; //array with phenotype columns
        int mu = numFactor / numPhenotype; //metagenes/phenotypes

        for (int i = 0; i < numFactor; i++) {
            int[] cmi = d2i(th.colMaxInd().getRow(0));
            int[] cmc = new int[numFactor * numPhenotype];

            for (int j = 0; j < cmi.length; j++) {
            	//System.err.println(getPhenotype(j));
                cmc[getPhenotype(j) * numFactor + cmi[j]]++;
            }

            int[] ind = rev(sort(cmc));
            int cmin = min(ct);

            for (int j = 0; j < ind.length; j++) {
                int tc = ind[j] % numFactor;
                int tr = (ind[j] - tc) / numFactor;

                if (ct[tr] == cmin) {
                    cm[tr * mu + ct[tr]] = tc;
                    ct[tr]++;
                    th.zeroRow(tc);
                    break;
                }
            }
        }

        return cm;
    }

/*--------------------------------------------------------------------
 maxClass - Create classes for the factors
--------------------------------------------------------------------*/

    public void maxClass() {
        numPhenotype = numGene;
        phenotypeName = new Texts();
        phenotype = new Phenotype[numPhenotype];

        for (int i = 0; i < numPhenotype; i++) {
            phenotypeName.rpush((i + 1));
            phenotype[i] = new Phenotype(phenotypeName.get(i), this);
        }

        Matrix d = x.colMaxInd();
        for (int i = 0; i < numSample; i++) {
            sample[i].addPhenotype(phenotypeName.get((int) (d.x[0][i])));
        }
    }

/*--------------------------------------------------------------------
 maxClass - Create classes for the factors
--------------------------------------------------------------------*/

    public void maxClass(int s) {
        numPhenotype = s;
        phenotypeName = new Texts();
        phenotype = new Phenotype[numPhenotype];
        int mu = s / numGene;

        for (int i = 0; i < numPhenotype; i++) {
            phenotypeName.rpush((i + 1));
            phenotype[i] = new Phenotype(phenotypeName.get(i), this);
        }

        Matrix d = x.colMaxInd();
        for (int i = 0; i < numSample; i++) {
            sample[i].addPhenotype(phenotypeName.get((int) (d.x[0][i]) * mu));
        }
    }

/*--------------------------------------------------------------------
 calcRadial - Calculate the radial coordinates
--------------------------------------------------------------------*/

    public double calcRadial(double[][] kc, double[][] sc, double maxdis, double e) {
        double[] xdis = x.euclidean();

        int k = numGene;
        Matrix cs = x.colSum();

        Matrix xc = x.copy();
        xc.div(cs);

        if (e != 1.0) {
            xc.apow(e);
        }

        cs = xc.colSum();
        xc.div(cs);

        Matrix yc = xc.copy();

        xc.mul(new Matrix(kc).selectCol(new int[]{0}));
        yc.mul(new Matrix(kc).selectCol(new int[]{1}));

        for (int i = 0; i < numSample; i++) {
            double d = min(1, xdis[i] / maxdis);

            if (e != 1.0) {
                d = 1 - pow((1 - d), e);
            }

            // d = 1.0;  // if stretching, can't have distance adjustment
            sc[i][0] = xc.colSum(i) * d;
            sc[i][1] = yc.colSum(i) * d;
        }

        double[] a = angle(kc);

        for (int i = 0; i < numSample; i++) {
            double sa = angle(sc[i][0], sc[i][1]);

            int j = 0;
            for (j = 0; j < numGene; j++) {
                if (sa < a[j]) {
                    break;
                }
            }

            int fa = j - 1;
            int fb = j;

            if (fb >= numGene) {
                fb = 0;
            }

            double dx = Math.cos(sa);
            double dy = Math.sin(sa);

            double[] mp = intersect(kc[fa][0], kc[fa][1], kc[fb][0], kc[fb][1], dx, dy, 0.0, 0.0);
            double mx = mp[0];
            double my = mp[1];

            double md = 1 / pow(mx * mx + my * my, 0.5);
            sc[i][0] *= md;
            sc[i][1] *= md;
        }

        Matrix d = new Matrix(sc).ntranspose();
        double[] sdis = d.euclidean();
        double mn = mean(sdis);

        return mn;
    }

/*--------------------------------------------------------------------
 heatmapDesc - Create a heatmap
--------------------------------------------------------------------*/

   /* public void heatmapDesc(String title, String fname, String tname, String dname) {
        boolean norm = true;
        String color = "redblue";
        Array a = copy();

        int headerHeight = 15;

        Text[][] tab = splitTabLine(tname);
        Texts setSizet = new Texts(tab[2]);
        setSizet.lpop();
        Texts leSizet = new Texts(tab[3]);
        leSizet.lpop();

        double[] setSize = setSizet.dbl();
        double[] leSize = leSizet.dbl();

        Texts des = new Texts();
        Texts nme = new Texts();

        for (int i = 0; i < numGene; i++) {
            des.rpush(a.gene[i].description);
            nme.rpush(a.geneName.get(i));
        }

        int ind[] = des.sort();
        nme = nme.select(ind);
        a.selectGene(nme);

        Text[] desord = splitTabLine(dname)[0];
        Text[] desnme = splitTabLine(dname)[0];
        int[][] dessmp = new int[desord.length][];
        int[] desa = new int[desord.length];
        int[] desb = new int[desord.length];
        int sze = 0;

        PDF p = new PDF(fname + ".pdf", false);
        p.title();
        // p.push ( 10, 10, -10, -10 );

        p.push(0, 0, 0, -headerHeight);
        upperFirst(desnme);

        double dw = p.strWidth(desnme);
        int dwi = (int) (dw);

        int numShown = 0;
        for (int i = 0; i < desord.length; i++) {
            dessmp[i] = des.finds(desord[i]);
            numShown += dessmp[i].length;
        }

        int space = 5;
        int totalHeight = p.ys - space * (desord.length - 1);
        double geneHeight = ((double) (totalHeight) / (double) (numShown));

        int setSizeWidth = 20;
        int setSizeSpace = 10;

        p.push(0, 0, -dwi - (setSizeWidth * 2 + setSizeSpace * 2), 0);

        for (int i = 0; i < desord.length; i++) {
            sze += dessmp[i].length;

            desa[i] = (int) (p.ys - (space * (i)) - (double) (sze) * geneHeight);
            desb[i] = desa[i] + (int) ((double) (dessmp[i].length) * geneHeight);
        }

        for (int i = 0; i < desord.length; i++) {
            p.push(0, desa[i], 0, desb[i]);
            a.heatmap(p, dessmp[i], norm, color);
            p.pop();
        }

        p.pop();

        setSize = log10(setSize);
        leSize = log10(leSize);

        double setSizeMax = max(setSize);
        double leSizeMax = max(leSize);

        // PLOT THE LEADING EDGE SET SIZE
        p.lineWidth(0.0);

        p.push(-dwi - (setSizeWidth * 2 + setSizeSpace), 0, -dwi - (setSizeWidth + setSizeSpace), 0);
        for (int i = 0; i < desord.length; i++) {
            p.push(0, desa[i], 0, desb[i]);

            p.stroke(0.75, 0.75, 0.75);
            p.fill(0.75, 0.75, 0.75);
            p.rect(0, 0, 1.0, 1.0);

            p.stroke(0.25, 0.25, 0.25);
            p.fill(0.25, 0.25, 0.25);

            for (int j = 0; j < dessmp[i].length; j++) {
                double y1 = 1.0 - (double) (j) / (double) (dessmp[i].length);
                double y2 = 1.0 - (double) (j + 1) / (double) (dessmp[i].length);
                double wid = (double) (leSize[dessmp[i][j]]) / (double) (leSizeMax);
                p.rect(1.0 - wid, y1, 1.0, y2);
            }
            p.pop();
        }
        p.pop();

        // PLOT THE GENE SET SIZE
        p.lineWidth(0.0);

        p.push(-dwi - (setSizeWidth + setSizeSpace), 0, -dwi - setSizeSpace, 0);
        for (int i = 0; i < desord.length; i++) {
            p.push(0, desa[i], 0, desb[i]);

            p.stroke(0.75, 0.75, 0.75);
            p.fill(0.75, 0.75, 0.75);
            p.rect(0, 0, 1.0, 1.0);

            p.stroke(0.5, 0.5, 0.5);
            p.fill(0.5, 0.5, 0.5);

            for (int j = 0; j < dessmp[i].length; j++) {
                double y1 = 1.0 - (double) (j) / (double) (dessmp[i].length);
                double y2 = 1.0 - (double) (j + 1) / (double) (dessmp[i].length);
                double wid = (double) (setSize[dessmp[i][j]]) / (double) (setSizeMax);
                p.rect(0.0, y1, wid, y2);
            }
            p.pop();
        }
        p.pop();

        p.lineWidth(1.0);
        p.push(-dwi, 0, 0, 0);
        p.stroke(0.25, 0.25, 0.25);
        p.fill(0.25, 0.25, 0.25);
        p.leftText(desnme, mean(desa, desb));
        p.pop();

        p.pop();

        p.push(0, -headerHeight + 2, 0, 0);

        p.push(-dwi - (setSizeWidth * 2 + setSizeSpace), 0, -dwi - setSizeSpace, 0);

        p.push(0, 0, 0, 3);
        p.stroke(0.25, 0.25, 0.25);
        p.fill(0.25, 0.25, 0.25);
        p.lineWidth(1.0);
        p.moveTo(0.0, 1.0);
        p.lineTo(0.0, 0.0);
        p.lineTo(1.0, 0.0);
        p.lineTo(1.0, 1.0);
        p.line(0.5, 0.0, 0.5, 1.0);
        p.stroke();
        p.pop();

        p.push(0, 3, 0, 0);
        p.setFontSize(8);
        p.centerText((int) (round(pow(10, leSizeMax), 0)) + "", 0.0, 0.5);
        p.centerText("0", 0.5, 0.5);
        p.centerText((int) (round(pow(10, setSizeMax), 0)) + "", 1.0, 0.5);
        p.pop();

        p.pop();

        p.pop();

        p.pop();
        p.close();

        a.writeGCT(fname + ".gct");
    }*/

/*--------------------------------------------------------------------
 Heatmaps
--------------------------------------------------------------------*/

    /**
     * Create a heatmap (redblue and normalized)
     *
     * @param title Title to be displayed
     * @param fname File to save the pdf
     */

    /*public void heatmap(String title, String fname) {
        heatmap(title, fname, true, "redblue");  // Call
    }*/

/*--------------------------------------------------------------------
 heatmap - Create a heatmap
--------------------------------------------------------------------*/

    /*public void heatmap(String title, String fname, boolean norm) {
        heatmap(title, fname, norm, "redblue");
    }*/

/*--------------------------------------------------------------------
 heatmap - Create a heatmap
--------------------------------------------------------------------*/

    /*public void heatmap(String title, String fname, boolean norm, String color) {
        heatmap(title, fname, norm, color, "right");
    }*/

/*--------------------------------------------------------------------
 heatmap - Create a heatmap
--------------------------------------------------------------------*/

    /*public void heatmap(String title, String fname, boolean norm, String color, String dir) {
        PDF p = new PDF(fname);
        p.title(title);
        p.push(10, 10, -10, -10);

        if (dir == "top") {
            p.push(0, 0, 0, -30);
            heatmap(p, seq(0, numGene - 1), norm, color);
            p.pop();
            p.push(0, -20, 0, 0);
            p.rect(cb);
            p.pop();
        } else if (dir == "bottom") {
            p.push(0, 30, 0, 0);
            heatmap(p, seq(0, numGene - 1), norm, color);
            p.pop();
            p.push(0, 0, 0, 20);
            p.rect(cb);
            p.pop();
        } else if (dir == "left") {
            p.push(30, 0, 0, 0);
            heatmap(p, seq(0, numGene - 1), norm, color);
            p.pop();
            p.push(0, 0, 20, 0);
            p.rect(cb);
            p.pop();
        } else {
            p.push(0, 0, -30, 0);
            heatmap(p, seq(0, numGene - 1), norm, color);
            p.pop();
            p.push(-20, 0, 0, 0);
            p.rect(cb);
            p.pop();
        }

        p.pop();
        p.close();
    }*/

/*--------------------------------------------------------------------
 dendrogram - Show the dendrogram
--------------------------------------------------------------------*/

    /*public void dendrogram(PDF p, ClusterHierarchical hc, boolean top) {
        int numSet = hc.setd.length - 1;
        int[] ord = order2rank(hc.order);

        double[] sx = new double[numSet];
        double[] sy = new double[numSet];

        if (top) {
            for (int i = 0; i < hc.nc; i++) {
                sx[i] = ((double) (ord[i]) + 0.5) / (double) (hc.nc);
                sy[i] = 0.0;
            }

            double maxd = max(hc.setd);
            for (int i = hc.nc; i < numSet; i++) {
                sx[i] = (sx[hc.seta[i]] + sx[hc.setb[i]]) / 2;
                sy[i] = 0.1 + (hc.setd[i] / maxd) * 0.9;
                p.moveTo(sx[hc.seta[i]], sy[hc.seta[i]]);
                p.lineTo(sx[hc.seta[i]], sy[i]);
                p.lineTo(sx[hc.setb[i]], sy[i]);
                p.lineTo(sx[hc.setb[i]], sy[hc.setb[i]]);
            }
        } else {
            for (int i = 0; i < hc.nc; i++) {
                sx[i] = 1 - ((double) (ord[i]) + 0.5) / (double) (hc.nc);
                sy[i] = 1.0;
            }

            double maxd = max(hc.setd);
            for (int i = hc.nc; i < numSet; i++) {
                sx[i] = (sx[hc.seta[i]] + sx[hc.setb[i]]) / 2;
                sy[i] = 1 - (0.1 + (hc.setd[i] / maxd) * 0.9);
                p.moveTo(sy[hc.seta[i]], sx[hc.seta[i]]);
                p.lineTo(sy[i], sx[hc.seta[i]]);
                p.lineTo(sy[i], sx[hc.setb[i]]);
                p.lineTo(sy[hc.setb[i]], sx[hc.setb[i]]);
            }
        }

        p.stroke();
    }*/

/*--------------------------------------------------------------------
 heatmap - Create a subset of a heatmap
--------------------------------------------------------------------*/

   /* public void heatmap(PDF p, int[] z, boolean norm) {
        heatmap(p, z, norm, "redblue");
    }*/

/*--------------------------------------------------------------------
 heatmap - Create a subset of a heatmap
--------------------------------------------------------------------*/

    /*public void heatmap(PDF p, int[] z, boolean norm, String color) {
        heatmap(p, z, norm, color, 50, true);
    }*/

/*--------------------------------------------------------------------
 heatmap - Create a subset of a heatmap
--------------------------------------------------------------------*/

   /* public void heatmap(PDF p, int[] z, boolean norm, String color, int hcwid, boolean geneFlag) {
        int hcw = hcwid;
        int hch = hcwid;

        if (shc == null) {
            hch = 0;
        }

        if (ghc == null) {
            hcw = 0;
        }

        if (shc != null) {
            p.lineWidth(1.0);
            p.stroke(0.5, 0.5, 0.5);
            p.push(hcw, -hch, p.xs, p.ys);
            dendrogram(p, shc, true);
            p.pop();
        }

        if (ghc != null) {
            p.lineWidth(1.0);
            p.stroke(0.5, 0.5, 0.5);
            p.push(0, 0, hcw, -hch);
            dendrogram(p, ghc, false);
            p.pop();
        }

        p.push(hcw, 0, p.xs, p.ys - hch);

        cb = new ColorBlend(color);

        double min = x.min();
        double mean = x.mean();
        double max = x.max();

        double rowHeight = p.ys / (double) (z.length);
        double colWidth = p.xs / (double) (numSample);

        p.lineWidth(0.0);
        p.stroke(0.25, 0.25, 0.25);

        if (geneFlag) {
            for (int i = 0; i < z.length; i++) {
                if (norm) {
                    min = x.rowMin(i);
                    mean = x.rowMean(i);
                    max = x.rowMax(i);
                }

                double mp = (mean - min) / (max - min);
                double nm = 0.5 / mp;
                double pm = 0.5 / (1 - mp);

                for (int j = 0; j < numSample; j++) {
                    double x1 = (double) (j) / (double) (numSample);
                    double y1 = (double) (z.length - i - 1) / (double) (z.length);
                    double x2 = (double) (j + 1) / (double) (numSample);
                    double y2 = (double) (z.length - i) / (double) (z.length);

                    double col = (x.x[z[i]][j] - min) / (max - min);

                    if (col < mp) {
                        col = col * nm;
                    } else {
                        col = 1 - ((1 - col) * pm);
                    }

                    p.stroke(cb.get(col));
                    p.fill(cb.get(col));
                    p.rect(x1, y1, x2, y2);
                }
            }
        } else {
            for (int j = 0; j < numSample; j++) {
                if (norm) {
                    min = x.colMin(j);
                    mean = x.colMean(j);
                    max = x.colMax(j);
                }

                double mp = (mean - min) / (max - min);
                double nm = 0.5 / mp;
                double pm = 0.5 / (1 - mp);

                for (int i = 0; i < z.length; i++) {
                    double x1 = (double) (j) / (double) (numSample);
                    double y1 = (double) (z.length - i - 1) / (double) (z.length);
                    double x2 = (double) (j + 1) / (double) (numSample);
                    double y2 = (double) (z.length - i) / (double) (z.length);

                    double col = (x.x[z[i]][j] - min) / (max - min);

                    if (col < mp) {
                        col = col * nm;
                    } else {
                        col = 1 - ((1 - col) * pm);
                    }

                    p.stroke(cb.get(col));
                    p.fill(cb.get(col));
                    p.rect(x1, y1, x2, y2);
                }
            }
        }

        if (rowHeight >= 5) {
            for (int i = 1; i < z.length; i++) {
                double x1 = 0;
                double x2 = 1;
                double y1 = (double) (z.length - i) / (double) (z.length);

                p.lineWidth(1.0);
                p.stroke(0.25, 0.25, 0.25);

                p.moveTo(x1, y1);
                p.lineTo(x2, y1);
                p.stroke();
            }
        }

        for (int j = 0; j <= numSample; j++) {
            double x1 = (double) (j) / (double) (numSample);
            double y1 = 0;
            double y2 = 1;

            boolean flag = false;

            if (j > 0 && j < numSample) {
                if (getPhenotype(j) != getPhenotype(j - 1)) {
                    flag = true;
                }
            }

            if (flag || colWidth >= 5) {
                if (flag) {
                    p.lineWidth(2.0);
                    p.stroke(0.0, 0.0, 0.0);
                } else {
                    p.lineWidth(1.0);
                    p.stroke(0.25, 0.25, 0.25);
                }

                p.moveTo(x1, y1);
                p.lineTo(x1, y2);
                p.stroke();
            }
        }

        p.lineWidth(1.0);
        p.stroke(0.0, 0.0, 0.0);

        p.moveTo(0, 0);
        p.lineTo(0, 1);
        p.lineTo(1, 1);
        p.lineTo(1, 0);
        p.lineTo(0, 0);

        p.stroke();

        p.pop();
    }*/

/*--------------------------------------------------------------------
 sortGene - Sort the genes
--------------------------------------------------------------------*/

    public void sortGene() {
        int[] si = geneName.sort();  // Sort the sample names
        selectGeneObject(si);      // Select the gene objects
    }

/*--------------------------------------------------------------------
 sortSample - Sort the samples
--------------------------------------------------------------------*/

    public void sortSample() {
        int[] si = sampleName.sort();  // Sort the sample names
        selectSampleObject(si);      // Select the samples
    }

/*--------------------------------------------------------------------
 sortPhenotype - Sort the phenotypes
--------------------------------------------------------------------*/

    public void sortPhenotype() {
        int ind[] = phenotypeName.sort();  // Sort the phenotype names
        selectPhenotypeObject(ind);
    }

/*--------------------------------------------------------------------
 sortByPhenotype - Sort by the phenotypes
--------------------------------------------------------------------*/

    public void sortByPhenotype() {
        int ind[] = phenotypeName.sort();  // Sort the phenotype names
        selectPhenotypeObject(ind);

        Text lst = new Text();
        for (int i = 0; i < phenotypeName.size; i++) {
            int s[] = selectPhenotypeIndexObject(new int[]{i});

            for (int j = 0; j < s.length; j++) {
                lst.rpush(s[j] + " ");
            }
        }

        ind = new Texts(lst.trimWhite().splitSpace()).num();
        selectSample(ind);
    }

/*--------------------------------------------------------------------
 sortByPhenotypeHierarchical - Sort by the phenotypes
--------------------------------------------------------------------*/

    public void sortByPhenotypeHierarchical() {
        // int indp[] = phenotypeName.sort ();  // Sort the phenotype names
        // selectPhenotypeObject ( indp );

        Text lst = new Text();
        for (int i = 0; i < phenotypeName.size; i++) {
            int s[] = selectPhenotypeIndexObject(new int[]{i});

            Matrix y = x.selectCol(s);
            ClusterHierarchical c = y.clusterHierarchical("pearson", "complete");

            for (int j = 0; j < s.length; j++) {
                lst.rpush(s[c.order[j]] + " ");
            }
        }

        int[] ind = new Texts(lst.trimWhite().splitSpace()).num();
        selectSample(ind);
    }

/*--------------------------------------------------------------------
 selectSampleObject - Select the samples
--------------------------------------------------------------------*/

    public void selectSampleObject(int[] si) {
        Sample[] ns = new Sample[si.length];  // Create a new sample array

        // LOOP THROUGH THE SAMPLE ARRAY
        for (int i = 0; i < si.length; i++) {
            ns[i] = sample[si[i]];  // Set the ith sample
        }

        sample = ns;            // Set to the new sample array
        x = x.selectCol(si);  // Set to the new columns
    }

/*--------------------------------------------------------------------
 selectSample - Select a set of samples
--------------------------------------------------------------------*/

    public void selectSample(Texts z) {
        int[] si = sampleName.select(z);  // Get the sample indices
        selectSample(si);                 // Select the samples
    }

/*--------------------------------------------------------------------
 selectSample - Select a set of samples
--------------------------------------------------------------------*/

    public void selectSample(Text[] z) {
        int[] si = sampleName.select(z);  // Get the sample indices
        selectSample(si);                 // Select the samples
    }

/*--------------------------------------------------------------------
 selectSample - Select a set of samples
--------------------------------------------------------------------*/

    public void selectSample(String[] z) {
        int[] si = sampleName.select(z);  // Get the sample indices
        selectSample(si);                 // Select the samples
    }

/*--------------------------------------------------------------------
 selectSample - Select a set of samples
--------------------------------------------------------------------*/

    public void selectSample(boolean[] si) {
        Texts ind = new Texts();
        for (int i = 0; i < si.length; i++) {
            if (si[i]) {
                ind.rpush(i);
            }
        }

        selectSample(ind.num());
    }

/*--------------------------------------------------------------------
 selectSample - Select a set of samples
--------------------------------------------------------------------*/

    public void selectSample(int[] si) {
        sampleName = sampleName.select(si);   // Select the sample names
        selectSampleObject(si);               // Select the sample objects
        numSample = sampleName.size;             // Update the number of samples
    }

/*--------------------------------------------------------------------
 selectPhenotypeObject - Select the phenotypes
--------------------------------------------------------------------*/

    public void selectPhenotypeObject(int[] si) {
        Phenotype[] ns = new Phenotype[si.length];  // Create a new phenotype array

        // LOOP THROUGH THE PHENOTYPE ARRAY
        for (int i = 0; i < si.length; i++) {
            ns[i] = phenotype[si[i]];  // Set the ith phenotype
        }

        phenotype = ns;  // Set to the new phenotype array

        Texts num = new Texts();

        // LOOP THROUGH THE SAMPLES
        for (int i = 0; i < numSample; i++) {
            // LOOP THORUGH THE PHENOTYPES
            for (int j = 0; j < numPhenotype; j++) {
                // IF THE SAMPLE HAS THE PHENOTYPE
                if (sample[i].isPhenotype(phenotypeName.get(j))) {
                    num.rpush(" " + i);
                    break;
                }
            }
        }

        selectSample(num.num());  // Select the samples
    }

/*--------------------------------------------------------------------
 selectPhenotype - Select a set of phenotypes
--------------------------------------------------------------------*/

    public void selectPhenotype(Texts z) {
        int[] si = phenotypeName.select(z);  // Get the phenotype indices
        selectPhenotype(si);                 // Select the phenotypes
    }

/*--------------------------------------------------------------------
 selectPhenotype - Select a set of phenotypes
--------------------------------------------------------------------*/

    public void selectPhenotype(Text[] z) {
        int[] si = phenotypeName.select(z);  // Get the phenotype indices
        selectPhenotype(si);                 // Select the phenotypes
    }

/*--------------------------------------------------------------------
 selectPhenotype - Select a set of phenotypes
--------------------------------------------------------------------*/

    public void selectPhenotype(String[] z) {
        int[] si = phenotypeName.select(z);  // Get the phenotype indices
        selectPhenotype(si);                 // Select the phenotypes
    }

/*--------------------------------------------------------------------
 selectPhenotype - Select a set of phenotypes
--------------------------------------------------------------------*/

    public void selectPhenotype(String z) {
        int[] si = phenotypeName.select(new String[]{z});  // Get the phenotype indices
        selectPhenotype(si);                                  // Select the phenotypes
    }

/*--------------------------------------------------------------------
 selectPhenotype - Select a phenotype
--------------------------------------------------------------------*/

    public void selectPhenotype(int i) {
        int[] si = new int[]{i};
        phenotypeName = phenotypeName.select(si);  // Select the phenotype names
        numPhenotype = phenotypeName.size;            // Update the number of phenotypes
        selectPhenotypeObject(si);                 // Select the phenotype objects
    }

/*--------------------------------------------------------------------
 selectPhenotype - Select a set of phenotypes
--------------------------------------------------------------------*/

    public void selectPhenotype(int[] si) {
        phenotypeName = phenotypeName.select(si);  // Select the phenotype names
        numPhenotype = phenotypeName.size;            // Update the number of phenotypes
        selectPhenotypeObject(si);                 // Select the phenotype objects
    }

/*--------------------------------------------------------------------
 selectPhenotypeIndexObject - Select the phenotypes
--------------------------------------------------------------------*/

    public int[] selectPhenotypeIndexObject(int[] si) {
        Texts num = new Texts();

        // LOOP THROUGH THE SAMPLES
        for (int i = 0; i < numSample; i++) {
            // LOOP THORUGH THE PHENOTYPES
            for (int j = 0; j < si.length; j++) {
                // IF THE SAMPLE HAS THE PHENOTYPE
                if (sample[i].isPhenotype(phenotypeName.get(si[j]))) {
                    num.rpush(" " + i);
                    break;
                }
            }
        }

        return num.num();  // Return the indices
    }

/*--------------------------------------------------------------------
 selectPhenotypeIndex - Select a set of phenotypes
--------------------------------------------------------------------*/

    public int[] selectPhenotypeIndex(Texts z) {
        int[] si = phenotypeName.select(z);  // Get the phenotype indices
        return selectPhenotypeIndex(si);     // Select the phenotypes
    }

/*--------------------------------------------------------------------
 selectPhenotypeIndex - Select a set of phenotypes
--------------------------------------------------------------------*/

    public int[] selectPhenotypeIndex(Text[] z) {
        int[] si = phenotypeName.select(z);  // Get the phenotype indices
        return selectPhenotypeIndex(si);     // Select the phenotypes
    }

/*--------------------------------------------------------------------
 selectPhenotypeIndex - Select a set of phenotypes
--------------------------------------------------------------------*/

    public int[] selectPhenotypeIndex(int i) {
        int[] si = phenotypeName.select(new String[]{phenotypeName.get(i).str()});  // Get the phenotype indices
        return selectPhenotypeIndex(si);     // Select the phenotypes
    }

/*--------------------------------------------------------------------
 selectPhenotypeIndex - Select a set of phenotypes
--------------------------------------------------------------------*/

    public int[] selectPhenotypeIndex(String z) {
        int[] si = phenotypeName.select(new String[]{z});  // Get the phenotype indices
        return selectPhenotypeIndex(si);     // Select the phenotypes
    }

/*--------------------------------------------------------------------
 selectPhenotypeIndex - Select a set of phenotypes
--------------------------------------------------------------------*/

    public int[] selectPhenotypeIndex(String[] z) {
        int[] si = phenotypeName.select(z);  // Get the phenotype indices
        return selectPhenotypeIndex(si);     // Select the phenotypes
    }

/*--------------------------------------------------------------------
 selectPhenotypeIndex - Select a set of phenotypes
--------------------------------------------------------------------*/

    public int[] selectPhenotypeIndex(int[] si) {
        return selectPhenotypeIndexObject(si);   // Select the phenotype objects
    }

/*--------------------------------------------------------------------
 mapHomolog - com.scanfeld.pr.Map genes to their homologs
--------------------------------------------------------------------*/

    public void mapHomolog(String fname) {

        Infile a = new Infile();         // Create a new com.scanfeld.io.Infile
        a.load(fname);                 // Load the homolog chart
        Text[][] z = a.splitTabLine();   // Split it up
        Texts g1 = new Texts(z[0]);  // Create the first gene list
        Texts g2 = new Texts(z[1]);  // Create the second gene list
        map(g1, g2);                   // com.scanfeld.pr.Map the homologs
    }

/*--------------------------------------------------------------------
 des - com.scanfeld.pr.Map descriptions onto the genes
--------------------------------------------------------------------*/

    public void des(String fname) {
        TextTree d = new TextTree(fname);

        for (int i = 0; i < geneName.size; i++) {
            Text td = d.get(geneName.get(i));

            if (td != null) {
                gene[i].description = td;
            }
        }
    }

/*--------------------------------------------------------------------
 map - com.scanfeld.pr.Map genes to a list of genes
--------------------------------------------------------------------*/

    public void map(String fname) {
        Infile a = new Infile();         // Create a new com.scanfeld.io.Infile
        a.load(fname);                 // Load the homolog chart
        Text[][] z = a.splitTabLine();   // Split it up
        Texts g1 = new Texts(z[0]);  // Create the first gene list
        Texts g2 = new Texts(z[1]);  // Create the second gene list
        map(g1, g2);                   // com.scanfeld.pr.Map the homologs
    }

/*--------------------------------------------------------------------
 map - com.scanfeld.pr.Map genes to a list of genes
--------------------------------------------------------------------*/

    public void map(Texts y, Texts z) {
        int ord[] = y.sort();  // Order the y gene list
        z = z.select(ord);   // Then order the z gene list

        selectGene(y);       // Select the desired genes

        int j = 0;              // Gene name counter

        // LOOP THROUGH THE GENE NAMES
        for (int i = 0; i < y.size && j < geneName.size; i++) {
            // WHILE GENE j IS LESS THAN GENE i
            while (i < y.size && y.get(i).lt(geneName.get(j))) {
                i++;  // Increment the lookup counter
            }

            // IF A MATCH IS FOUND
            if (j < geneName.size && i < y.size && y.get(i).eq(geneName.get(j))) {
                geneName.set(j, z.get(i));  // Set the new gene name
            }

            // IF NO MATCH WAS FOUND
            else {
                i--;  // Decrement the lookup counter
            }

            j++;  // Increment the gene name counter
        }
    }

/*--------------------------------------------------------------------
 removeDescription - Remove genes with this description
--------------------------------------------------------------------*/

    public void removeDescription(String z) {
        boolean b[] = new boolean[geneName.size];  // Boolean flag array

        // LOOP THROUGH THE ARRAY
        for (int i = 0; i < geneName.size; i++) {
            if (gene[i].description.find(z) >= 0) {
                b[i] = false;
            } else {
                b[i] = true;
            }
        }

        selectGene(b);
    }

/*--------------------------------------------------------------------
 randomMix - Create a random mix of two phenotypes
--------------------------------------------------------------------*/

    public Array randomMix(String a, String b, int groupSize, int groupDiv) {
        int numGroup = groupDiv - 1;

        Matrix ex = new Matrix(numGene, groupSize * numGroup);
        Texts ns = new Texts();
        Texts np = new Texts();
        Texts ph = new Texts();

        for (int i = 0; i < numGroup; i++) {
            int n = (int) Math.round(100 * (double) (i + 1) / (double) (groupDiv));
            np.rpush("m" + n);
        }

        int[] ai = selectPhenotypeIndex(a);
        int[] bi = selectPhenotypeIndex(b);

        int alen = ai.length;
        int blen = bi.length;

        int[] aa = randInt(0, alen - 1, groupSize);
        int[] bb = randInt(0, blen - 1, groupSize);

        for (int i = 0; i < numGroup; i++) {
            for (int j = 0; j < groupSize; j++) {
                double n = (double) (i + 1) / (double) (groupDiv);
                ns.rpush("S" + (j + 1) + "." + (int) Math.round(100 * n) + "." + ai[aa[j]] + "." + bi[bb[j]]);

                for (int l = 0; l < numGene; l++) {
                    ex.x[l][i * groupSize + j] = (1 - n) * x.x[l][ai[aa[j]]];
                    ex.x[l][i * groupSize + j] += n * x.x[l][bi[bb[j]]];
                }

                ph.rpush(np.get(i));
            }
        }

        Array r = new Array(this, ex, ns, np, ph);

        return r;
    }

/*--------------------------------------------------------------------
 normalizeRow - Normalize each row
--------------------------------------------------------------------*/

    public void normalizeRow() {
        x.normalizeRow();
    }

/*--------------------------------------------------------------------
 normalizeRowMean - Normalize each row
--------------------------------------------------------------------*/

    public void normalizeRowMean() {
        x.normalizeRowMean();
    }

/*--------------------------------------------------------------------
 randomMax - Create a random max of two phenotypes
--------------------------------------------------------------------*/

    public Array randomMax(String a, String b, int groupSize) {
        Matrix ex = new Matrix(numGene, groupSize);
        Texts ns = new Texts();
        Texts np = new Texts();
        Texts ph = new Texts();

        np.rpush("max");

        int[] ai = selectPhenotypeIndex(a);
        int[] bi = selectPhenotypeIndex(b);

        int alen = ai.length;
        int blen = bi.length;

        int[] aa = randInt(0, alen - 1, groupSize);
        int[] bb = randInt(0, blen - 1, groupSize);

        for (int i = 0; i < groupSize; i++) {
            ns.rpush("S" + (i + 1) + "." + ai[aa[i]] + "." + bi[bb[i]]);

            for (int l = 0; l < numGene; l++) {
                ex.x[l][i] = max(x.x[l][ai[aa[i]]], x.x[l][bi[bb[i]]]);
            }

            ph.rpush(np.get(0));
        }

        Array r = new Array(this, ex, ns, np, ph);

        return r;
    }

/*--------------------------------------------------------------------
 sizePhenotypeDiverse - Resize each phenotype to be equal size
--------------------------------------------------------------------*/

    public void sizePhenotypeDiverse() {
        int min = 0;

        Texts ni = new Texts();
        Matrix[] m = new Matrix[numPhenotype];
        int[][] ind = new int[numPhenotype][];

        for (int i = 0; i < numPhenotype; i++) {
            ind[i] = selectPhenotypeIndex(i);
            m[i] = x.copy();
            m[i] = m[i].selectCol(ind[i]);

            if (ind[i].length < min || min == 0) {
                min = ind[i].length;
            }
        }

        for (int i = 0; i < numPhenotype; i++) {
            Matrix p = m[i].pearson();
            p.add(1);

            int[] ch = new int[min];
            for (int j = 0; j < min; j++) {
                double[] s;

                if (j == 0) {
                    s = p.rowSum().getRow(0);
                } else {
                    Matrix p2 = p.selectCol(ch);
                    p2 = p2.selectCol(seq(0, j - 1));
                    s = p2.rowSum().getRow(0);
                    // out ( s );
                }

                int[] si = sort(s);
                // out ( j + " " + si[ 0 ] );
                ni.rpush(ind[i][si[0]]);
                ch[j] = si[0];

                for (int k = 0; k < p.nr; k++) {
                    p.x[si[0]][k] = 10;
                }
            }
        }

        // out ( ni.joinTab () );
        selectSample(ni.num());
        // out ( sampleName.joinTab () );
    }

/*--------------------------------------------------------------------
 sizePhenotype - Resize each phenotype to be equal
--------------------------------------------------------------------*/

    public void sizePhenotype() {
        int min = 0;

        for (int i = 0; i < numPhenotype; i++) {
            int[] ind = selectPhenotypeIndex(i);
            if (ind.length < min || min == 0) {
                min = ind.length;
            }
        }

        randomSample(min);
    }

/*--------------------------------------------------------------------
 randomGene - Get k random genes
--------------------------------------------------------------------*/

    public void randomGene(int k) {
        int[] ind = randIntArray(numGene, k);
        selectGene(ind);
    }

/*--------------------------------------------------------------------
 random - Add noise to the samples
--------------------------------------------------------------------*/

    public void random(double z) {
        x.addRandom(z);
    }

/*--------------------------------------------------------------------
 randomGene - Get k random genes
--------------------------------------------------------------------*/

    public void randomGene(double k) {
        if (k < 1) {
            k = numGene * k;
        }

        if (k > 0) {
            randomGene((int) (k));
        }
    }

/*--------------------------------------------------------------------
 randomSample - Get k random samples from each phenotype
--------------------------------------------------------------------*/

    public void randomSample(int k) {
        Texts s = new Texts();

        for (int i = 0; i < numPhenotype; i++) {
            int[] ind = selectPhenotypeIndex(i);
            int[] rnd = randInt(ind.length);

            for (int j = 0; j < k && j < ind.length; j++) {
                s.rpush(ind[rnd[j]]);
            }
        }

        selectSample(s.num());
    }

/*--------------------------------------------------------------------
 randomSampleReplace - Get k random samples with replacement from
                       each phenotype
--------------------------------------------------------------------*/

    public void randomSampleReplace(int k) {
        Texts s = new Texts();

        for (int i = 0; i < numPhenotype; i++) {
            int[] ind = selectPhenotypeIndex(i);
            int[] rnd = randInt(0, (ind.length) - 1, k);

            for (int j = 0; j < k; j++) {
                s.rpush(ind[rnd[j]]);
            }
        }

        selectSample(s.num());
    }

/*--------------------------------------------------------------------
 randomSampleReplace - Get k random samples with replacement from
                       each phenotype
--------------------------------------------------------------------*/

    public void randomSampleReplaceFull(int k) {
        Texts s = new Texts();

        for (int i = 0; i < numPhenotype; i++) {
            int[] ind = selectPhenotypeIndex(i);
            int[] rnd = randIntArray(ind.length, k);

            for (int j = 0; j < k; j++) {
                s.rpush(ind[rnd[j]]);
            }
        }

        selectSample(s.num());
    }
    
    
    /*--------------------------------------------------------------------
    randomizeArray - For each gene randomize the values by sample
   --------------------------------------------------------------------*/

       public Array randomizeArray() {
           Array rtrn=this.copy();
           
           for(int i=0; i<rtrn.x.nr; i++){
        	   double[] array=rtrn.x.x[i];
        	   double[] randomized=randomize(array);
        	   rtrn.x.x[i]=randomized;
           	}
           return rtrn;
       }   
       
       private double[] randomize(double[] array){
    	  ArrayList temp=array2List(array);
    	   Collections.shuffle(temp);
    	   double[] rtrn=list2Array(temp);
    	   return rtrn;
       }

       private ArrayList array2List(double[] array){
    	   ArrayList rtrn=new ArrayList();
    	   for(int i=0; i<array.length; i++){rtrn.add(array[i]);}
    	   return rtrn;
       }
       
       private double[] list2Array(ArrayList<Double> list){
    	   double[] rtrn=new double[list.size()];
    	   for(int i=0; i<list.size(); i++){rtrn[i]=list.get(i);}
    	   return rtrn;
       }
       
/*--------------------------------------------------------------------
 removeDuplicate - Remove duplicate genes
--------------------------------------------------------------------*/

    public void removeDuplicate() {
        Matrix m = x.foldChange();  // Find the indices of the satisfying genes;
        double v[] = m.col(0);
        int ord[] = sort(v);
        ord = rev(ord);
        selectGene(ord);

        Tree z = new Tree();                      // Temporary lookup tree
        boolean b[] = new boolean[geneName.size];  // Boolean flag array

        // LOOP THROUGH THE ARRAY
        for (int i = 0; i < geneName.size; i++) {
            // IF THE TEXT HAS ALREADY BEEN ENTERED
            if (z.get(geneName.x[geneName.k + i]) > 0) {
                b[i] = false;  // Set the flag to false
            }

            // IF THE TEXT HAS NOT BEEN ENTERED
            else {
                z.set(geneName.x[geneName.k + i], (i + 1));  // Add the ith element
                b[i] = true;                    // Set the flag to true
            }
        }

        selectGene(b);
        sortGene();
    }

/*--------------------------------------------------------------------
 removeNull - Remove null genes
--------------------------------------------------------------------*/

    public void removeNull() {
        Texts ind = new Texts();

        for (int i = 0; i < numGene; i++) {
            if (!missing[i] && geneName.get(i).ne("NULL") && geneName.get(i).ne("null")) {
                ind.rpush(i);
            }
        }

        selectGene(ind.num());  // Select the genes
        missing = null;             // Empty the missing array
    }

/*--------------------------------------------------------------------
 trimSampleName - Remove extra letters from sample names
--------------------------------------------------------------------*/

    public void trimSampleName() {
        for (int i = 0; i < numSample; i++) {
            sampleName.set(i, sampleName.get(i).alphaNumSign().removeExtraWhite());
            sample[i].name = sample[i].name.alphaNumSign().removeExtraWhite();
        }
    }

/*--------------------------------------------------------------------
 renameSample - Rename all samples to their phenotypes with an index
--------------------------------------------------------------------*/

    public void renameSample() {
        Tree t = new Tree();

        // LOOP THROUGH THE SAMPLE ARRAY
        for (int i = 0; i < numSample; i++) {
            Text p = sample[i].phenotype.get(0);
            int j = t.get(p);

            if (j > 0) {
                j++;
                t.set(p, j);
            } else {
                j = 1;
                t.set(p, j);
            }

            Text newName = c(p, s2t("."), i2t(j), s2t("."), sampleName.get(i));
            sample[i].name = newName;
            sampleName.set(i, newName);
        }
    }

/*--------------------------------------------------------------------
 renameSample - Rename all samples
--------------------------------------------------------------------*/

    public void renameSample(String[] z) {
        // LOOP THROUGH THE SAMPLE ARRAY
        for (int i = 0; i < z.length; i++) {
            Text newName = s2t(z[i]);
            sample[i].name = newName;
            sampleName.set(i, newName);
        }
    }

/*--------------------------------------------------------------------
 renamePhenotype - Rename the phenotypes
--------------------------------------------------------------------*/

    public void renamePhenotype(String s) {
        Text z = new Text(s);
        Text[] y = z.splitSpace();

        for (int i = 0; i < y.length; i++) {
            Text[] x = y[i].split("=");

            if (x.length >= 2) {
                renamePhenotype(x[0].str(), x[1].str());
            }
        }
    }

/*--------------------------------------------------------------------
 prefixPhenotype - Add a prefix to all phenotypes
--------------------------------------------------------------------*/

    public void prefixPhenotype(String a) {
        for (int i = 0; i < numPhenotype; i++) {
            renamePhenotype(phenotypeName.get(i).str(), a + phenotypeName.get(i).str());
        }
    }

/*--------------------------------------------------------------------
 renamePhenotype - Rename a phenotype
--------------------------------------------------------------------*/

    public void renamePhenotype(String a, String b) {
        Text ta = s2t(a);
        Text tb = s2t(b);

        for (int i = 0; i < numPhenotype; i++) {
            if (phenotypeName.get(i).eq(ta)) {
                phenotypeName.set(i, tb);
                phenotype[i].name = tb;
            }
        }

        for (int i = 0; i < numSample; i++) {
            if (sample[i].isPhenotype(ta)) {
                for (int j = 0; j < sample[i].phenotype.size; j++) {
                    if (sample[i].phenotype.get(j).eq(ta)) {
                        sample[i].phenotype.set(j, tb);
                    }
                }
            }
        }
    }

/*--------------------------------------------------------------------
 removeBadSample - Remove any incorrect samples
--------------------------------------------------------------------*/

    public void removeBadSample(double bc) {
        double[] br = brier();
        int[] ind = colMaxInd();
        int[] ph = getPhenotype();

        boolean[] ns = new boolean[numSample];

        for (int i = 0; i < numSample; i++) {
            if (br[i] < bc || ph[i] != ind[i]) {
                ns[i] = false;
            } else {
                ns[i] = true;
            }
        }

        selectSample(ns);
    }

/*--------------------------------------------------------------------
 countBadSample - Count the incorrect samples
--------------------------------------------------------------------*/

    public int countBadSample(double bc) {
        double[] br = brier();
        int[] ind = colMaxInd();
        int[] ph = getPhenotype();

        int ct = 0;
        for (int i = 0; i < numSample; i++) {
            if (br[i] < bc || ph[i] != ind[i]) {
                ct++;
            }
        }

        return ct;
    }

/*--------------------------------------------------------------------
 buildModel - Build an nmf model
--------------------------------------------------------------------*/

    public void buildModel(String name, String fname, boolean evalFlag) {
        buildModel(name, fname, evalFlag, 0, 0, 3, 0.5, 3);
    }

/*--------------------------------------------------------------------
 buildModel - Build an nmf model
--------------------------------------------------------------------*/

    public void buildModel(String name, String fname, boolean evalFlag, int svm_type, int kernel_type, double c, double gamma, int degree) {
        //write(fname + ".gct", fname + ".cls");
        write(fname + ".gct", fname + ".cls", true);

        nmf(numPhenotype);
        nmferr = x.err[x.erri];
       // writeNMF(name + " All Samples", fname + ".mod");
        writeNMF(fname + ".mod", true);

        Array n = new Array(fname + ".gct", fname + ".cls");
        Array p = n.project(fname + ".mod.W.gct");
        p.x.hypersphere(p.x.max());
        //p.write(fname + ".mod.proj.gct", fname + ".mod.proj.cls");
        p.write(fname + ".mod.proj.gct", fname + ".mod.proj.cls", true);
        Array s = p.evaluate(p, svm_type, kernel_type, c, gamma, degree);
        //s.write(fname + ".mod.eval.gct");
        s.write(fname + ".mod.eval.gct", true);

        /*
        out(x.H);
        //x.H = p.x;
        if ( evalFlag )
        {
          evaluateH ( svm_type, kernel_type, c, gamma, degree );
        }

        else
        {
          x.H.nonNegUnit ();
        }

        sizeForNMF ( minPhenotype () );
        */

        //remove bad sample with hardcoded score less than 0.3
        s.removeBadSample(0.3);
        
        
        Texts lst = new Texts();
        double[] br = s.brier();

        /*
        for ( int i = 0; i < s.numSample; i++ )
        {
          out ( sampleName.get ( i ) + " " + br[ i ] + " " + s.x.x[ 0 ][ i ] + " " + s.x.x[ 1 ][ i ] + " " + s.x.x[ 2 ][ i ] );
        }
        */

        int k = s.minPhenotype();
        for (int i = 0; i < s.phenotype.length; i++) {
            int[] ind = s.selectPhenotypeIndex(i);
            int[] pi = sort(sel(br, ind));

            /*
            // PICK THE BEST
            for ( int j = ind.length - 1; j >= ind.length - k; j-- )
            {
              lst.rpush ( ind[ pi[ j ] ] );
            }
            */

            // PICK THE WORST
            for (int j = 0; j < k; j++) {
                lst.rpush(ind[pi[j]]);
            }
        }

        int[] sel = lst.num();
        s.selectSample(sel);
        selectSample(s.sampleName);

        nmf(numPhenotype);
        nmferrtr = x.err[x.erri];

        //writeNMF(name + " Trimmed", fname + ".mod.tr");
        writeNMF(fname + ".mod.tr", true);
       // write(fname + ".tr.gct", fname + ".tr.cls");
        write(fname + ".tr.gct", fname + ".tr.cls", true);

        n = new Array(fname + ".gct", fname + ".cls");
        p = n.project(fname + ".mod.tr.W.gct");
        p.x.hypersphere(p.x.max());
       // p.write(fname + ".mod.tr.proj.gct", fname + ".mod.tr.proj.cls");
        p.write(fname + ".mod.tr.proj.gct", fname + ".mod.tr.proj.cls", true);
        s = p.evaluate(p, svm_type, kernel_type, c, gamma, degree);
       // s.write(fname + ".mod.tr.eval.gct");
        s.write(fname + ".mod.tr.eval.gct", true);
    }

/*--------------------------------------------------------------------
 buildModel - Build an nmf model
--------------------------------------------------------------------*/

    public void buildModel(String fname, boolean evalFlag, int svm_type, int kernel_type, double c, double gamma, int degree, double badCutoff, String sampleSelect) {
        write(fname + ".gct", fname + ".cls");

        nmf(numPhenotype);
        nmferr = x.err[x.erri];
        writeNMF("", fname + ".mod");

        Array n = new Array(fname + ".gct", fname + ".cls");
        Array p = n.project(fname + ".mod.W.gct");
        p.x.hypersphere(p.x.max());
        p.write(fname + ".mod.proj.gct", fname + ".mod.proj.cls");
        Array s = p.evaluate(p, svm_type, kernel_type, c, gamma, degree);
        s.write(fname + ".mod.eval.gct");

        /*
        out(x.H);
        //x.H = p.x;
        if ( evalFlag )
        {
          evaluateH ( svm_type, kernel_type, c, gamma, degree );
        }

        else
        {
          x.H.nonNegUnit ();
        }

        sizeForNMF ( minPhenotype () );
        */

        if (badCutoff > 0.0) {
            s.removeBadSample(badCutoff);
        }

        Texts lst = new Texts();
        double[] br = s.brier();

        /*
        for ( int i = 0; i < s.numSample; i++ )
        {
          out ( sampleName.get ( i ) + " " + br[ i ] + " " + s.x.x[ 0 ][ i ] + " " + s.x.x[ 1 ][ i ] + " " + s.x.x[ 2 ][ i ] );
        }
        */

        if (sampleSelect.equals("best") || sampleSelect.equals("worst")) {
            int k = s.minPhenotype();
            for (int i = 0; i < s.phenotype.length; i++) {
                int[] ind = s.selectPhenotypeIndex(i);
                int[] pi = sort(sel(br, ind));

                if (sampleSelect.equals("best")) {
                    // PICK THE BEST
                    for (int j = ind.length - 1; j >= ind.length - k; j--) {
                        lst.rpush(ind[pi[j]]);
                    }
                } else {
                    // PICK THE WORST
                    for (int j = 0; j < k; j++) {
                        lst.rpush(ind[pi[j]]);
                    }
                }
            }

            int[] sel = lst.num();
            s.selectSample(sel);

            selectSample(s.sampleName);
        }

        nmf(numPhenotype);
        nmferrtr = x.err[x.erri];

        writeNMF(name + " Trimmed", fname + ".mod.tr");
        write(fname + ".tr.gct", fname + ".tr.cls");

        n = new Array(fname + ".gct", fname + ".cls");
        p = n.project(fname + ".mod.tr.W.gct");
        p.x.hypersphere(p.x.max());
        p.write(fname + ".mod.tr.proj.gct", fname + ".mod.tr.proj.cls");
        s = p.evaluate(p, svm_type, kernel_type, c, gamma, degree);
        s.write(fname + ".mod.tr.eval.gct");
    }

    /*--------------------------------------------------------------------
    globalShift - Reduce all negative values to positive linear addition of min
   --------------------------------------------------------------------*/

       public void globalShift() {
           //get min of all samples
    	   double min=min(x.x);
    	   //System.err.println(min);
    	   
    	   //add min to all values
    	   globalShift(-min);
           
       }
       
   /*--------------------------------------------------------------------
   globalShift - Reduce all negative values to positive linear addition of min
   --------------------------------------------------------------------*/

     public void globalShift(double addFactor) {
    	for(int i=0; i<x.x.length; i++){
    		 for(int j=0; j<x.x[i].length; j++){
    			 x.x[i][j]+=addFactor;
    		 }
    	 }
     }

     /*--------------------------------------------------------------------
     min - get minimum value in matrix
     --------------------------------------------------------------------*/

       public double min(double[][] y) {
    	   double min=Double.MAX_VALUE;
    	   for(int i=0; i<y.length; i++){
    		   for(int j=0; j<y[i].length; j++){
    			   min=Math.min(min, y[i][j]);
    		   }
    	   }
    	   return min;
       }
     
/*--------------------------------------------------------------------
 nmf - Run the nmf algorithm
--------------------------------------------------------------------*/

    public double nmf(int k) {
        x.NMFdiv(k);
        //int[] ind = maxFactorH();
        //x.H = x.H.selectRow(ind);
        //x.W = x.W.selectCol(ind);
        return x.getNMFError();
    }
    
    public double nmf(int k, NMFCostFunction costFunction) {
        x.NMFdiv(k, costFunction);
        //int[] ind = maxFactorH();
        //x.H = x.H.selectRow(ind);
        //x.W = x.W.selectCol(ind);
        return x.getNMFError();
    }
    
    public double[] nmf(int k, int numIterations, NMFCostFunction costFunction) {
        return nmf(k, numIterations, costFunction, 0.0);
    }
    
    public double[] nmf(int k, int numIterations, NMFCostFunction costFunction, double precision) {
        double minCost=Double.MAX_VALUE;
    	
        Matrix minH=null;
        Matrix minW=null;
        
        double[] costs=new double[numIterations];
        
        for(int i=0; i<numIterations; i++){
        	x.NMFdiv(k, costFunction, precision);
    		double cost = x.getNMFError();
//    		if(i == 0) {
//    			double newPrecision = cost/(double)10000;
//    			x.setPrecision(newPrecision);
//    			System.err.println("\tAdjusted precision " + newPrecision);
//    		}
    		if(cost<minCost){minCost=cost; minH=x.H.copy(); minW=x.W.copy();}
    		costs[i]=cost;
        }
        
        x.H=minH;
        x.W=minW;
        
        return costs;
    }
    
    public double[] nmf(Matrix Hi, int numIterations, NMFCostFunction costFunction, double precision, int divFactor, String save) {
        double minCost=Double.MAX_VALUE;
    	
        Matrix minH=null;
        Matrix minW=null;
        
        double[] costs=new double[numIterations];
        
        for(int i=0; i<numIterations; i++){
        	Matrix HRandom=addRandom(Hi, divFactor);
        	//try{HRandom.write(save+".Iteration"+i+".gct");}catch(IOException ex){ex.printStackTrace();}
        	x.NMFdiv(HRandom, costFunction, precision);
    		double cost = x.getNMFError();
    		if(cost<minCost){minCost=cost; minH=x.H.copy(); minW=x.W.copy();}
    		costs[i]=cost;
        }
        
        x.H=minH;
        x.W=minW;
        
        return costs;
    }
    
    //Add random variable to each row
    private Matrix addRandom(Matrix H, int divFactor){
    	double[][] vals=H.x;
    	
    	for(int i=0; i<vals.length; i++){
    		for(int j=0; j<vals[i].length; j++){
    			if(divFactor>0){vals[i][j]+=(Math.random()/divFactor);}
    		}
    	}
    	
    	return new Matrix(vals);
    }
    
    
    public double[] nmf(int k, int numIterations) {        
        return nmf(k,numIterations, null);
    }
   
    public double[] nmf(int k, int numIterations, double precision) {        
        return nmf(k,numIterations, null, precision);
    }
    
    
    public double nmf(Matrix Hi, NMFCostFunction costFunction) {
    	x.NMFdiv(Hi, costFunction);
        double cost = x.getNMFError();
       // int[] ind = maxFactorH();
       // x.H = x.H.selectRow(ind);
       // x.W = x.W.selectCol(ind);
        return cost;
    }
    
    public double nmf(Matrix Hi) {
    	return nmf(Hi, Matrix.DEFAULT_COST_FUNCTION);
    }
    

/*--------------------------------------------------------------------
 sizeForNMF - Resize for the NMF set
--------------------------------------------------------------------*/

    public void sizeForNMF(int k) {
        Texts lst = new Texts();
        double[] br = x.H.brier();

        /*
        out ( br );

        if ( fflag )
        {
          for ( int i = 0; i < numSample; i++ )
          {
            // UPDATE FOR MULTIPLE ROWS PER PHENOTYPE
            br[ i ] = x.H.x[ getPhenotype ( i ) ][ i ];
          }
        }
        */

        // out ( br );
        // out ( x.H );

        for (int i = 0; i < phenotype.length; i++) {
            int[] ind = selectPhenotypeIndex(i);
            int[] pi = sort(sel(br, ind));

            for (int j = ind.length - 1; j >= ind.length - k; j--) {
                lst.rpush(ind[pi[j]]);
            }
        }

        int[] sel = lst.num();
        selectSample(sel);
        // out ( sampleName.joinLine () );
    }

/*--------------------------------------------------------------------
 writeNMF - Run the nmf algorithm
--------------------------------------------------------------------*/

    public void writeNMF(String name, String fname) {
        Outfile a = new Outfile(fname + ".W.gct");
        a.writeLine("#1.2");                      // Version number
        a.writeLine(x.W.nr + "\t" + x.W.nc);      // Size information
        a.write("Name\tDescription");             // Name and description

        for (int i = 0; i < x.W.nc; i++) {
            a.write("\t" + "F" + (i + 1));
        }

        a.write("\n");

        for (int i = 0; i < x.W.nr; i++) {
            a.write(geneName.get(i) + "\t");

            for (int j = 0; j < x.W.nc; j++) {
                a.write("\t" + x.W.x[i][j]);
            }

            a.write("\n");
        }

        a.close();

        a = new Outfile(fname + ".H.gct");
        a.writeLine("#1.2");                      // Version number
        a.writeLine(x.H.nr + "\t" + x.H.nc);      // Size information
        a.write("Name\tDescription");             // Name and description

        for (int i = 0; i < x.H.nc; i++) {
            a.write("\t" + sampleName.get(i));
        }

        a.write("\n");

        for (int i = 0; i < x.H.nr; i++) {
            a.write("F" + (i + 1) + "\t");

            for (int j = 0; j < x.H.nc; j++) {
                a.write("\t" + x.H.x[i][j]);
            }

            a.write("\n");
        }

        a.close();

        /*  ASDF - uncomment when you have hierarchical clustering working
        Array W = new Array ( fname + ".W.gct" );
        W.clusterHierarchical ( "pearson", "complete" );
        W.clusterHierarchicalGene ( "pearson", "complete" );
        W.heatmap ( name + " W Heatmap", fname + ".W.pdf", false );

        Array H = new Array ( fname + ".H.gct" );
        H.clusterHierarchical ( "pearson", "complete" );
        H.clusterHierarchicalGene ( "pearson", "complete" );
        H.heatmap ( name + " H Heatmap", fname + ".H.pdf", false );
        */

        // H.radial ( name + " H Radial", fname + ".H.RAD.pdf" );
    }
    
    /*--------------------------------------------------------------------
    getNMF - Get NMF decomposition matrix
   --------------------------------------------------------------------*/
    public Matrix getNMF(){
    	return x;
    }
    
    
    
    /*--------------------------------------------------------------------
    writeNMF - Run the nmf algorithm
   --------------------------------------------------------------------*/

       public void writeNMF(String fname) {
           Outfile a = new Outfile(fname + ".W.gct");
           a.writeLine("#1.2");                      // Version number
           a.writeLine(x.W.nr + "\t" + x.W.nc);      // Size information
           a.write("Name\tDescription");             // Name and description

           for (int i = 0; i < x.W.nc; i++) {
               a.write("\t" + "F" + (i + 1));
           }

           a.write("\n");

           for (int i = 0; i < x.W.nr; i++) {
               a.write(geneName.get(i) + "\t");

               for (int j = 0; j < x.W.nc; j++) {
                   a.write("\t" + x.W.x[i][j]);
               }

               a.write("\n");
           }

           a.close();

           a = new Outfile(fname + ".H.gct");
           a.writeLine("#1.2");                      // Version number
           a.writeLine(x.H.nr + "\t" + x.H.nc);      // Size information
           a.write("Name\tDescription");             // Name and description

           for (int i = 0; i < x.H.nc; i++) {
               a.write("\t" + sampleName.get(i));
           }

           a.write("\n");

           for (int i = 0; i < x.H.nr; i++) {
               a.write("F" + (i + 1) + "\t");

               for (int j = 0; j < x.H.nc; j++) {
                   a.write("\t" + x.H.x[i][j]);
               }

               a.write("\n");
           }

           a.close();

           /*  ASDF - uncomment when you have hierarchical clustering working
           Array W = new Array ( fname + ".W.gct" );
           W.clusterHierarchical ( "pearson", "complete" );
           W.clusterHierarchicalGene ( "pearson", "complete" );
           W.heatmap ( name + " W Heatmap", fname + ".W.pdf", false );

           Array H = new Array ( fname + ".H.gct" );
           H.clusterHierarchical ( "pearson", "complete" );
           H.clusterHierarchicalGene ( "pearson", "complete" );
           H.heatmap ( name + " H Heatmap", fname + ".H.pdf", false );
           */

           // H.radial ( name + " H Radial", fname + ".H.RAD.pdf" );
       }
       
       
       
       
       public void  writeNMF (String fname, boolean newWriter) {
           try{
    	   //System.err.println("test");
    	   FileWriter writer = new FileWriter(fname + ".W.gct");
           writer.write("#1.2\n");                      // Version number
           writer.write(x.W.nr + "\t" + x.W.nc+"\n");      // Size information
           writer.write("Name\tDescription");             // Name and description

           for (int i = 0; i < x.W.nc; i++) {
               writer.write("\t" + "F" + (i + 1));
           }

           writer.write("\n");

           for (int i = 0; i < x.W.nr; i++) {
               writer.write(geneName.get(i) + "\t"+descriptions.get(i));

               for (int j = 0; j < x.W.nc; j++) {
                   writer.write("\t" + x.W.x[i][j]);
               }

               writer.write("\n");
           }

           writer.close();
           
           if(x.C!=null){
           writer = new FileWriter(fname + ".C.gct");
           writer.write("#1.2\n");                      // Version number
           writer.write(x.C.nr + "\t" + x.C.nc+"\n");      // Size information
           writer.write("Name\tDescription");             // Name and description

           for (int i = 0; i < x.C.nc; i++) {
               writer.write("\t" + "F" + (i + 1));
           }

           writer.write("\n");

           for (int i = 0; i < x.C.nr; i++) {
               writer.write(geneName.get(i) + "\t"+descriptions.get(i));

               for (int j = 0; j < x.C.nc; j++) {
                   writer.write("\t" + x.C.x[i][j]);
               }

               writer.write("\n");
           }

           writer.close();
           }

           FileWriter a = new FileWriter(fname + ".H.gct");
           a.write("#1.2\n");                      // Version number
           a.write(x.H.nr + "\t" + x.H.nc+"\n");      // Size information
           a.write("Name\tDescription");             // Name and description

           for (int i = 0; i < x.H.nc; i++) {
               a.write("\t" + sampleName.get(i));
           }

           a.write("\n");

           for (int i = 0; i < x.H.nr; i++) {
               a.write("F" + (i + 1) + "\t");

               for (int j = 0; j < x.H.nc; j++) {
                   a.write("\t" + x.H.x[i][j]);
               }

               a.write("\n");
           }

           a.close();

           /*  ASDF - uncomment when you have hierarchical clustering working
           Array W = new Array ( fname + ".W.gct" );
           W.clusterHierarchical ( "pearson", "complete" );
           W.clusterHierarchicalGene ( "pearson", "complete" );
           W.heatmap ( name + " W Heatmap", fname + ".W.pdf", false );

           Array H = new Array ( fname + ".H.gct" );
           H.clusterHierarchical ( "pearson", "complete" );
           H.clusterHierarchicalGene ( "pearson", "complete" );
           H.heatmap ( name + " H Heatmap", fname + ".H.pdf", false );
           */

           // H.radial ( name + " H Radial", fname + ".H.RAD.pdf" );
           }catch(IOException ex){}
       }

       
       public void  writeNMF (String fname, Set<Integer> badMetagenes, boolean newWriter) {
           try{
    	   //System.err.println("test");
    	   FileWriter writer = new FileWriter(fname + ".W.gct");
           writer.write("#1.2\n");                      // Version number
           writer.write(x.W.nr + "\t" + x.W.nc+"\n");      // Size information
           writer.write("Name\tDescription");             // Name and description

           for (int i = 0; i < x.W.nc; i++) {
               writer.write("\t" + "F" + (i + 1));
           }

           writer.write("\n");

           for (int i = 0; i < x.W.nr; i++) {
               writer.write(geneName.get(i) + "\t");

               for (int j = 0; j < x.W.nc; j++) {
                   writer.write("\t" + x.W.x[i][j]);
               }

               writer.write("\n");
           }

           writer.close();

           FileWriter a = new FileWriter(fname + ".H.gct");
           a.write("#1.2\n");                      // Version number
           a.write(x.H.nr + "\t" + x.H.nc+"\n");      // Size information
           a.write("Name\tDescription");             // Name and description

           for (int i = 0; i < x.H.nc; i++) {
               a.write("\t" + sampleName.get(i));
           }

           a.write("\n");

           for (int i = 0; i < x.H.nr; i++) {
               a.write("F" + (i + 1) + "\t");

               for (int j = 0; j < x.H.nc; j++) {
                   a.write("\t" + x.H.x[i][j]);
               }

               a.write("\n");
           }

           a.close();
           
           
           a = new FileWriter(fname + ".H.filteredBAD.gct");
           a.write("#1.2\n");                      // Version number
           a.write((x.H.nr-badMetagenes.size()) + "\t" + x.H.nc+"\n");      // Size information
           a.write("Name\tDescription");             // Name and description

           for (int i = 0; i < x.H.nc; i++) {
               a.write("\t" + sampleName.get(i));
           }

           a.write("\n");

           for (int i = 0; i < x.H.nr; i++) {
        	   if(!badMetagenes.contains(i)){
               a.write("F" + (i + 1) + "\t");

               for (int j = 0; j < x.H.nc; j++) {
                   a.write("\t" + x.H.x[i][j]);
               }

               a.write("\n");
        	   }
           }

           a.close();

           /*  ASDF - uncomment when you have hierarchical clustering working
           Array W = new Array ( fname + ".W.gct" );
           W.clusterHierarchical ( "pearson", "complete" );
           W.clusterHierarchicalGene ( "pearson", "complete" );
           W.heatmap ( name + " W Heatmap", fname + ".W.pdf", false );

           Array H = new Array ( fname + ".H.gct" );
           H.clusterHierarchical ( "pearson", "complete" );
           H.clusterHierarchicalGene ( "pearson", "complete" );
           H.heatmap ( name + " H Heatmap", fname + ".H.pdf", false );
           */

           // H.radial ( name + " H Radial", fname + ".H.RAD.pdf" );
           }catch(IOException ex){}
       }
       
       public ArrayList<String> getHeader(){
    	   ArrayList header=new ArrayList();
    	   for (int i = 0; i < x.nc; i++) {
    		   header.add(id+"_"+sampleName.get(i).toString());
           }
    	   return header;
       }
       
       public ArrayList getGenes(){
    	   ArrayList genes=new ArrayList();
    	   for (int i = 0; i < x.nr; i++) {
              genes.add(geneName.get(i));
    	   	}
    	   return genes;
       }
       
     /*  public ArrayList<double[]>  filterH (Set<Integer> badMetagenes) {
    	   ArrayList list=new ArrayList();
           //need to remove bad rows
    	   
    	   for (int i = 0; i < x.H.nr; i++) {
        	   if(!badMetagenes.contains(i)){
        		   list.add(x.H.x[i]);
        		}
           }

        return list;   
       }*/
       
       public Array  filterH (Set<Integer> badMetagenes) {
    	   ArrayList list=new ArrayList();
           //need to remove bad rows
    	   
    	   
    	   x.H.removeRows(badMetagenes);
    	   x.W.removeColumns(badMetagenes);
        	   //if(badMetagenes.contains(i)){
        		 //  System.err.println("removing row "+i+" of total "+x.H.nr);
        		  // x.H.removeRow(i);
        		   //x.W.removeColumn(i);
        	   //}
           

        return this;   
       }
       
       public Matrix  filterH (Set<Integer> badMetagenes, boolean test) {
    	   //need to remove bad rows
    	   
    	 for(Integer rowIndex: badMetagenes){
    		 x.H.removeRow(rowIndex);
    	 }

        return x.H;   
       }
       
       public ArrayList<double[]>  filterW (Set<Integer> badMetagenes) {
    	   ArrayList list=new ArrayList();
         
    	   //need to remove bad columns
    	   for(int j=0; j<x.W.nr; j++){
    		   double[] array=new double[x.W.nc-badMetagenes.size()];
    		   int count=0;
    		   for (int i = 0; i < x.W.nc; i++) {
    			   
    		   		if(!badMetagenes.contains(i)){
        		   array[count++]=x.W.x[j][i];
    			   
    		   	}
    		   		list.add(array);
           	}
    	   }

        return list;   
       }
       
/*--------------------------------------------------------------------
 brier - Calculate the brier score
--------------------------------------------------------------------*/

    public double[] brier() {
        return x.brier();
    }

/*--------------------------------------------------------------------
 brierNoNorm - Calculate the brier score
--------------------------------------------------------------------*/

    public double[] brierNoNorm() {
        int[] p = getPhenotype();
        return x.brierNoNorm(p);
    }

/*--------------------------------------------------------------------
 colMaxInd - Calculate the indices of the maxiumum row in each column
--------------------------------------------------------------------*/

    public int[] colMaxInd() {
        return d2i(x.colMaxInd().x[0]);
    }

/*--------------------------------------------------------------------
 colMaxIndH - Calculate the indices of the maxiumum row in each column
--------------------------------------------------------------------*/

    public int[] colMaxIndH() {
        return d2i(x.H.colMaxInd().x[0]);
    }

/*--------------------------------------------------------------------
 evaluateH - Evaluate the H matrix
--------------------------------------------------------------------*/

    public void evaluateH() {
        evaluateH(0, 0, 3, 0.5, 3);
    }

/*--------------------------------------------------------------------
 evaluateH - Evaluate the H matrix
--------------------------------------------------------------------*/

    public void evaluateH(int svm_type, int kernel_type, double c, double gamma, int degree) {
        if (svmm == null) {
            model(svm_type, kernel_type, c, gamma, degree);
        }

        Matrix m = x.H.copy();
        m.scale(m.max() / svm_max);  // Think this fixes the bug
        svm sv = new svm();

        /*
        // ORIGINAL STATE
        pa.svm_type = 0;  // C-SVC
        pa.kernel_type = 0;  // linear
        pa.cache_size = 40;
        pa.C = 2 ^ 3;
        pa.gamma = 2 ^ -5;
        //pa.degree = 3;
        pa.eps = 1e-5;
        pa.probability = 1;
        */

        /*
        // ORIGINAL STATE
        pa.svm_type = 0;  // C-SVC
        pa.kernel_type = 0;  // linear
        pa.cache_size = 40;
        pa.C = 3;
        //pa.gamma = 1 / 2;
        //pa.degree = 3;
        pa.eps = 1e-5;
        pa.probability = 1;
        */

        double[] pe = new double[2];
        for (int i = 0; i < numPhenotype; i++) {
            svm_node[][] tv = new svm_node[numSample][numPhenotype];

            for (int j = 0; j < numSample; j++) {
                for (int k = 0; k < numPhenotype; k++) {
                    tv[j][k] = new svm_node();
                    tv[j][k].index = k + 1;
                    tv[j][k].value = x.H.x[k][j];
                }
            }

            svm_problem pr = new svm_problem();

            pr.l = numSample;
            pr.y = new double[numSample];
            pr.x = tv;

            svm_parameter pa = new svm_parameter();

            pa.svm_type = svm_type;
            pa.kernel_type = kernel_type;
            pa.cache_size = 40;
            pa.C = c;
            pa.gamma = gamma;
            pa.degree = degree;
            pa.eps = 1e-5;
            pa.probability = 1;

            Text s = phenotypeName.get(i);
            int ind = 1;

            if (sample[0].isPhenotype(s)) {
                ind = 0;
            }

            for (int j = 0; j < numSample; j++) {
                if (sample[j].isPhenotype(s)) {
                    pr.y[j] = 1;
                } else {
                    pr.y[j] = -1;
                }
            }

            svm_model mo = svm.svm_train(pr, pa);

            for (int j = 0; j < numSample; j++) {
                double pred = svm.svm_predict_probability(mo, tv[j], pe);
                m.x[i][j] = pe[ind];
            }
        }

        m.div(m.colSum());  // asdf - do we want to add up to 1?
        x.H = m;
    }

/*--------------------------------------------------------------------
 model - Model
--------------------------------------------------------------------*/

    public void model() {
        model(0, 0, 3, 0.5, 3);
    }

/*--------------------------------------------------------------------
 model - Model
--------------------------------------------------------------------*/

    public void model(int svm_type, int kernel_type, double c, double gamma, int degree) {
        Matrix m = x.copy();
        svm_max = m.max();
        m.scale();

        svm_node[][] tv = new svm_node[numSample][numGene];

        for (int i = 0; i < numSample; i++) {
            for (int j = 0; j < numGene; j++) {
                tv[i][j] = new svm_node();
                tv[i][j].index = j + 1;
                tv[i][j].value = m.x[j][i];
            }
        }

        sv = new svm();

        svm_problem pr = new svm_problem();
        pr.l = numSample;
        pr.y = new double[numSample];
        pr.x = tv;

        svp = new svm_parameter();

        svp.svm_type = svm_type;
        svp.kernel_type = kernel_type;
        svp.cache_size = 40;
        svp.C = c;
        svp.gamma = gamma;
        svp.degree = degree;

        /*
        svp.svm_type = 0;  // C-SVC
        svp.kernel_type = 0;  // linear
        svp.cache_size = 40;
        svp.C = 3;
        //svp.gamma = 1 / 2;
        //svp.degree = 3;
        */

        svp.eps = 1e-5;
        svp.probability = 1;

        svmm = new svm_model[numPhenotype];
        svmi = new int[numPhenotype];

        for (int i = 0; i < numPhenotype; i++) {
            Text s = phenotypeName.get(i);
            svmi[i] = 1;

            if (sample[0].isPhenotype(s)) {
                svmi[i] = 0;
            }

            for (int j = 0; j < numSample; j++) {
                if (sample[j].isPhenotype(s)) {
                    pr.y[j] = 1;
                } else {
                    pr.y[j] = -1;
                }
            }

            svmm[i] = svm.svm_train(pr, svp);
        }
    }

/*--------------------------------------------------------------------
 evaluate - Evaluate the projection
--------------------------------------------------------------------*/

    public Array evaluate(Array a) {
        return evaluate(a, 0, 0, 3, 0.5, 3);
    }

/*--------------------------------------------------------------------
 evaluate - Evaluate the projection
--------------------------------------------------------------------*/

    public Array evaluate(Array a, int svm_type, int kernel_type, double c, double gamma, int degree) {
        if (svmm == null) {
            model(svm_type, kernel_type, c, gamma, degree);
        }

        Matrix t = a.x.copy();
        t.scale(t.max() / svm_max);

        Array m = a.copy();
        m.x = new Matrix(numPhenotype, a.numSample);
        svm_node[][] rv = new svm_node[a.numSample][a.numGene];

        for (int i = 0; i < a.numSample; i++) {
            for (int j = 0; j < a.numGene; j++) {
                rv[i][j] = new svm_node();
                rv[i][j].index = j + 1;
                rv[i][j].value = t.x[j][i];
            }
        }

        double[] pe = new double[2];
        for (int i = 0; i < numPhenotype; i++) {
            for (int j = 0; j < a.numSample; j++) {
                double pred = svm.svm_predict_probability(svmm[i], rv[j], pe);
                m.x.x[i][j] = pe[svmi[i]];
            }
        }

        m.x.div(m.x.colSum());  // asdf - do we want to add up to 1?
        m.selectGene(seq(0, numPhenotype - 1));

        return m;
    }

/*--------------------------------------------------------------------
 project - Project a datset by W
--------------------------------------------------------------------*/

    public Array project(String oname, Array W, String fname) {
        return project(oname, W, fname, 0, 0, 3, 0.5, 3);
    }

/*--------------------------------------------------------------------
 project - Project a datset by W
--------------------------------------------------------------------*/

    public Array project(String oname, Array W, String fname, int svm_type, int kernel_type, double c, double gamma, int degree) {
        return project(oname, W, fname, -1, svm_type, kernel_type, c, gamma, degree);
    }

/*--------------------------------------------------------------------
 project - Project a datset by W
--------------------------------------------------------------------*/

    public Array project(String oname, Array W, String fname, int num) {
        return project(oname, W, fname, num, 0, 0, 3, 0.5, 3);
    }

/*--------------------------------------------------------------------
 project - Project a datset by W
--------------------------------------------------------------------*/

    public Array project(String oname, Array W, String fname, int num, int svm_type, int kernel_type, double c, double gamma, int degree) {
        Texts[] pr = new Texts[numGene];
        Texts[] pre = new Texts[numPhenotype];

        for (int i = 0; i < numGene; i++) {
            pr[i] = new Texts();
        }

        for (int i = 0; i < numPhenotype; i++) {
            pre[i] = new Texts();
        }

        Texts nme = new Texts();
        Array a, p, s;

        // remove this afterward
        // com.scanfeld.core.Tree gt = new com.scanfeld.core.Tree ();
        // com.scanfeld.core.Texts geneList = new com.scanfeld.core.Texts ();

        if (num <= 1) {
            a = new Array(fname + ".gct");
            nme.rpush(a.sampleName);

            // for ( int j = 0; j < a.numGene; j++ )
            // {
            //  if ( gt.get ( a.geneName.get ( j ) ) < 1 )
            //  {
            //    gt.set ( a.geneName.get ( j ), 1 );
            //    geneList.rpush ( a.geneName.get ( j ) );
            //  }
            // }

            p = a.project(W.copy());

            for (int j = 0; j < a.numSample; j++) {
                for (int k = 0; k < numGene; k++) {
                    pr[k].rpush(p.x.x[k][j]);
                }
            }

            s = evaluate(p, svm_type, kernel_type, c, gamma, degree);

            for (int j = 0; j < a.numSample; j++) {
                for (int k = 0; k < numGene; k++) {
                    pre[k].rpush(s.x.x[k][j]);
                }
            }
        } else {
            for (int i = 0; i < num; i++) {
                a = new Array(fname + i + ".gct");
                nme.rpush(a.sampleName);

                // for ( int j = 0; j < a.numGene; j++ )
                // {
                //  if ( gt.get ( a.geneName.get ( j ) ) < 1 )
                //  {
                //    gt.set ( a.geneName.get ( j ), 1 );
                //    geneList.rpush ( a.geneName.get ( j ) );
                //  }
                // }

                p = a.project(W.copy());

                for (int j = 0; j < a.numSample; j++) {
                    for (int k = 0; k < numGene; k++) {
                        pr[k].rpush(p.x.x[k][j]);
                    }
                }

                s = evaluate(p, svm_type, kernel_type, c, gamma, degree);

                for (int j = 0; j < a.numSample; j++) {
                    for (int k = 0; k < numPhenotype; k++) {
                        pre[k].rpush(s.x.x[k][j]);
                    }
                }
            }
        }

        // geneList.sort ();
        // geneList = geneList.removeDuplicate ();

        // a = new Array ( "c:/com.scanfeld/pr/ma/gs/jd.pre.sort.gct", "c:/com.scanfeld/pr/ma/gs/jd.pre.sort.cls" );
        // a.selectGene ( geneList );
        // a.write ( "c:/com.scanfeld/pr/ma/gs/jd.pre.sort.yeast.gct", "c:/com.scanfeld/pr/ma/gs/jd.pre.sort.yeast.cls" );

        Matrix px = new Matrix(numGene, pr[0].size);
        Matrix sx = new Matrix(numPhenotype, pre[0].size);

        Texts cnme = new Texts();
        for (int i = 0; i < numGene; i++) {
            cnme.rpush("F" + (i + 1));

            for (int j = 0; j < pr[0].size; j++) {
                px.x[i][j] = pr[i].get(j).dbl();
            }
        }

        Array pxa = new Array(px, cnme, nme);
        pxa.writeGCT(oname + ".proj.gct");

        cnme = new Texts();
        for (int i = 0; i < numPhenotype; i++) {
            cnme.rpush("F" + (i + 1));

            for (int j = 0; j < pre[0].size; j++) {
                sx.x[i][j] = pre[i].get(j).dbl();
            }
        }

        Array sxa = new Array(sx, cnme, nme);
        sxa.writeGCT(oname + ".eval.gct");

        return pxa;
    }

/*--------------------------------------------------------------------
 mpEvalTab - Create a mp eval tab file
--------------------------------------------------------------------*/

    public void mpEvalTab(String fname, String aname, double[] brierCuts, String[] brierNames) {
        mpEvalTab(fname, aname, -1, brierCuts, brierNames);
    }

/*--------------------------------------------------------------------
 mpEvalTab - Create an mp eval tab file
--------------------------------------------------------------------*/

    public void mpEvalTab(String fname, String aname, int act, double[] brierCuts, String[] brierNames) {
        TextTree map = new TextTree();
        Texts cat = new Texts();

        int asze = act;

        if (asze <= 0) {
            asze = 1;
        }

        for (int i = 0; i < asze; i++) {
            String afn = aname;

            if (act > 0) {
                afn += i;
            }

            afn += ".map";

            Infile mapf = new Infile(afn);
            mapf.in();

            while (!mapf.eof) {
                Text[] z = mapf.splitTab();
                cat.rpush(z[1]);
                map.add(z[1], "\t" + z[0].alphaNumSign().removeExtraWhite());
                mapf.in();
            }

            mapf.close();

            cat.sort();
            cat = cat.removeDuplicate();
        }

        Texts[] sl = new Texts[numSample];
        int[] cct = new int[cat.size];
        int[][] pct = new int[cat.size][numGene];
        double[][] hsc = new double[cat.size][numGene];
        Tree st = new Tree();

        for (int i = 0; i < numSample; i++) {
            sl[i] = new Texts();
            st.set(sampleName.get(i).alphaNumSign().removeExtraWhite(), i + 1);
        }

        for (int i = 0; i < cat.size; i++) {
            Text[] z = map.get(cat.get(i)).trimWhite().splitTab();

            for (int j = 0; j < z.length; j++) {
                sl[st.get(z[j]) - 1].rpush(i);
            }
        }

        for (int i = 0; i < sl.length; i++) {
            sl[i] = sl[i].removeDuplicate();
        }

        double[] br = brier();
        int mi[] = colMaxInd();
        int bind[] = rev(sort(br));
        br = rev(br);
        mi = sel(mi, bind);

        for (int i = 0; i < brierCuts.length; i++) {
            String bn = brierNames[i];
            double bc = brierCuts[i];

            Outfile o = new Outfile(fname + "." + bc + ".set.score.txt");

            for (int j = 0; j < cat.size; j++) {
                cct[j] = 0;

                for (int k = 0; k < numGene; k++) {
                    pct[j][k] = 0;
                }
            }

            int numGeneAbove = 0;
            int[] numGeneAbovePhenotype = new int[numGene];

            for (int j = 0; j < numSample; j++) {
                if (br[j] < bc) {
                    break;
                }

                numGeneAbove++;
                numGeneAbovePhenotype[mi[j]]++;

                //  Text sname = sampleName.get ( bind[ j ] ).removeExtraWhite ();

                for (int k = 0; k < sl[bind[j]].size; k++) {
                    int cind = sl[bind[j]].get(k).num();
                    cct[cind]++;
                    pct[cind][mi[j]]++;
                }
            }

            for (int j = 0; j < cat.size; j++) {
                for (int k = 0; k < numGene; k++) {
                    hsc[j][k] = hyperGeometric(pct[j][k], cct[j], numGeneAbovePhenotype[k], numGeneAbove);
                }
            }

            String z = "Name\tTotal";

            for (int k = 0; k < numGene; k++) {
                z += "\t" + (k + 1);
            }

            for (int k = 0; k < numGene; k++) {
                z += "\t" + (k + 1) + " Score";
            }

            o.writeLine(z);

            for (int j = 0; j < cat.size; j++) {
                z = cat.get(j) + "\t" + cct[j];

                for (int k = 0; k < numGene; k++) {
                    z += "\t" + pct[j][k];
                }

                for (int k = 0; k < numGene; k++) {
                    z += "\t" + hsc[j][k];
                }

                o.writeLine(z);
            }

            o.close();

            o = new Outfile(fname + "." + bc + ".matrix.txt");

            for (int j = 0; j < numSample; j++) {
                if (br[j] < bc) {
                    break;
                }

                int[] d = b2b(i2b(sl[bind[j]].num(), cat.size));
                o.writeLine(joinTab(d));
            }

            o.close();

            o = new Outfile(fname + "." + bc + ".sample.txt");

            Texts ns = sampleName.select(bind);
            for (int j = 0; j < numSample; j++) {
                if (br[j] < bc) {
                    break;
                }

                o.writeLine(ns.get(j));
            }

            o.close();

            o = new Outfile(fname + "." + bc + ".class.txt");

            for (int j = 0; j < numSample; j++) {
                if (br[j] < bc) {
                    break;
                }

                o.writeLine("" + (mi[j] + 1));
            }

            o.close();

            o = new Outfile(fname + "." + bc + ".sample.score.txt");

            String s = "Name";

            for (int k = 0; k < numGene; k++) {
                s += "\t" + (k + 1);
            }

            s += "\tBrier\tClass\tDescription";

            o.writeLine(s);

            for (int j = 0; j < numSample; j++) {
                if (br[j] < bc) {
                    break;
                }

                s = "" + ns.get(j);

                for (int k = 0; k < numGene; k++) {
                    s += "\t" + round(x.x[k][bind[j]], -2);
                }

                s += "\t" + round(br[j], -2) + "\t" + (mi[j] + 1) + "\t";

                for (int k = 0; k < sl[bind[j]].size; k++) {
                    int cind = sl[bind[j]].get(k).num();

                    if (k > 0) {
                        s += ", ";
                    }

                    s += cat.get(cind);
                }

                o.writeLine(s);
            }

            o.close();

            o = new Outfile(fname + "." + bc + ".set.txt");
            o.write(cat.joinLine());
            o.close();
        }
    }

/*--------------------------------------------------------------------
 plotProjection - Plot the projection
--------------------------------------------------------------------*/

    /*public void plotProjection(String name, String oname, String aname, boolean evalFlag, Array c) {
        double[] brierCuts = new double[]{0.0, 0.5, 0.75, 0.9, 10, 20, 50};
        String[] brierNames = new String[]{">= 0.0", ">= 0.5", ">= 0.75", ">= 0.9", "Top 10", "Top 20", "Top 50"};

        // brierCuts = new double[] { 0.75 };
        // brierNames = new String[] { ">= 0.75" };

        plotProjection(name, oname, aname, evalFlag, c, brierCuts, brierNames);
    }*/

/*--------------------------------------------------------------------
 plotProjection - Plot the projection
--------------------------------------------------------------------*/

   /* public void plotProjection(String name, String oname, String aname, boolean evalFlag, Array c, double[] brierCuts, String[] brierNames) {
        double maxdis = 1.0;

        if (evalFlag) {
            oname += ".eval";
        } else {
            oname += ".proj";
        }

        Array b = this.copy();
        b.trimSampleName();

        if (!evalFlag) {
            b.x.nonNegUnit();
            maxdis = max(b.x.euclidean());
        }

        Array a = b.copy();
        if (c == null) {
            c = a.copy();
        }

        c = c.copy();
        c.trimSampleName();

        if (c.numGene > numPhenotype) {
            int sze = c.numGene / numPhenotype;

            for (int i = 0; i < c.numSample; i++) {
                for (int j = 0; j < numPhenotype; j++) {
                    double max = 0;

                    for (int k = 0; k < sze; k++) {
                        if (c.x.x[j * sze + k][i] > max) {
                            max = c.x.x[j * sze + k][i];
                        }
                    }

                    for (int k = 0; k < sze; k++) {
                        c.x.x[j * sze + k][i] = max;
                    }
                }
            }

            int[] gind = new int[numPhenotype];

            for (int i = 0; i < numPhenotype; i++) {
                gind[i] = i * sze;
            }

            c.selectGene(gind);
        }

        for (int bi = 0; bi < brierCuts.length; bi++) {
            double brierCutOff = brierCuts[bi];

            double[] brs = c.brier();
            int[] ind = sort(copy(brs));

            double[] bc = new double[c.numPhenotype];
            for (int i = 0; i < c.numPhenotype; i++) {
                bc[i] = brierCutOff;
            }

            if (brierCutOff > 1) {
                for (int i = 0; i < c.numPhenotype; i++) {
                    double pct = 0;

                    for (int m = 0; m < c.numSample; m++) {
                        int j = ind[c.numSample - m - 1];
                        int ph = c.getPhenotype(j);

                        if (ph == i) {
                            pct++;

                            if (abs(pct - brierCutOff) < 0.001) {
                                bc[i] = brs[j];
                                break;
                            }
                        }
                    }

                    if (bc[i] > 1) {
                        bc[i] = 0;
                    }
                }
            }

            TextTree map = new TextTree();
            Texts cat = new Texts();
            Infile mapf = new Infile(aname);
            mapf.in();

            while (!mapf.eof) {
                Text[] z = mapf.splitTab();
                cat.rpush(z[1]);
                map.add(z[1], "\t" + z[0].alphaNumSign().removeExtraWhite());
                mapf.in();
            }

            cat.sort();
            cat = cat.removeDuplicate();

            Tree[] maps = new Tree[cat.size];
            for (int i = 0; i < cat.size; i++) {
                maps[i] = new Tree();
                Text[] z = map.get(cat.get(i)).trimWhite().splitTab();

                for (int j = 0; j < z.length; j++) {
                    maps[i].set(z[j], 1);
                }
            }

            mapf.close();

            cat.size = 0;  // ASDF - remove this to see the categories
            if (cat.size == 0) {
                a.radial(name + " (" + brierNames[bi] + ")", oname + ".RAD." + brierCutOff + ".pdf", 1.0, bc, 0, "brier", c);
            } else {
                a.radial(name + " (" + brierNames[bi] + ")", oname + ".RAD." + brierCutOff + ".pdf", 1.0, bc, 1, "brier", c);
            }

            PDF curPDF = a.curPDF;

            for (int i = 0; i < cat.size; i++) {
                a = b.copy();
                a.curPDF = curPDF;
                Texts z = new Texts(map.get(cat.get(i)).trimWhite().splitTab());

                a.selectSample(z);
                Array c2 = c.copy();
                c2.selectSample(z);

                if (i < cat.size - 1) {
                    a.radial(name + " " + cat.get(i) + " (" + brierNames[bi] + ")", oname + ".RAD." + brierCutOff + "." + (i + 1) + ".pdf", 1.0, bc, 2, "brier", c2);
                } else {
                    a.radial(name + " " + cat.get(i) + " (" + brierNames[bi] + ")", oname + ".RAD." + brierCutOff + "." + (i + 1) + ".pdf", 1.0, bc, 3, "brier", c2);
                }
            }

            // FIX BRIER SCORE STUFF AFTER HERE
            // why two brier score runs?

            a = b.copy();
            Outfile proj = new Outfile(oname + "." + brierCutOff + ".tab");

      

            double[] br = a.brier();
            int mi[] = a.colMaxInd();
            int bind[] = rev(sort(br));
            br = rev(br);

            Tree[] grps = new Tree[a.numGene];
            int gsze[] = new int[a.numGene];

            for (int i = 0; i < a.numGene; i++) {
                grps[i] = new Tree();
            }

            for (int j = 0; j < a.numSample; j++) {
                int i = bind[j];

                if ((brierCutOff < 1 && br[j] >= brierCutOff) || (brierCutOff > 1 && gsze[mi[i]] < brierCutOff)) {
                    gsze[mi[i]]++;

                    for (int k = 0; k < cat.size; k++) {
                        if (maps[k].get(a.sampleName.get(i)) > 0) {
                            // Increment the category's count in that group
                            grps[mi[i]].inc(cat.get(k));
                        }
                    }
                }
            }

            int[] catsze = new int[cat.size];

            for (int i = 0; i < a.numGene; i++) {
                // proj.write ( "Class " + ( i + 1 ) + " Counts\t\t\t\t\t\t" );

                for (int j = 0; j < cat.size; j++) {
                    // proj.write ( "\t" + grps[ i ].get ( cat.get ( j ) ) );
                    catsze[j] += grps[i].get(cat.get(j));
                }

                // proj.writeLine ();
            }

            Matrix catp = new Matrix(a.numGene, cat.size);

            for (int i = 0; i < a.numGene; i++) {
                // proj.write ( "Class " + ( i + 1 ) + " %\t\t\t\t\t\t" );

                for (int j = 0; j < cat.size; j++) {
                    // BRIER SCORE ON CATEGORY SIZE
                    // if ( catsze[ i ] > 0 )
                    // {
                    //   catp.x[ i ][ j ] = ( double ) ( grps[ i ].get ( cat.get ( j ) ) ) / ( double ) ( catsze[ i ] );
                    // }

                    // BRIER SCORE ON CATEGORY PERCENTAGE
                    if (gsze[i] > 0) {
                        catp.x[i][j] = (double) (grps[i].get(cat.get(j))) / gsze[i];
                    } else {
                        catp.x[i][j] = 0;
                    }

                    // proj.write ( "\t" + round ( catp.x[ i ][ j ], -3 ) );
                }

                // proj.writeLine ();
            }

            Matrix catm = catp.colMaxInd();
            double[] catb = catp.brier();

          

            proj.write("Name\tDescription\tBrier Score\tClass");
            for (int j = 0; j < a.numGene; j++) {
                proj.write("\t" + (j + 1));
            }

            for (int j = 0; j < cat.size; j++) {
                proj.write("\t" + cat.get(j));
            }

            proj.writeLine();

            gsze = new int[a.numGene];

            for (int pi = 0; pi < a.numGene; pi++) {
                for (int j = 0; j < a.numSample; j++) {
                    int i = bind[j];

                    if (mi[i] == pi) {
                        if ((brierCutOff < 1 && br[j] >= brierCutOff) || (brierCutOff > 1 && gsze[mi[i]] < brierCutOff)) {
                            gsze[mi[i]]++;

                            proj.write(a.sampleName.get(i) + "\t");
                            proj.write("\t" + br[j] + "\t" + (mi[i] + 1));

                            for (int k = 0; k < a.numGene; k++) {
                                proj.write("\t" + a.x.x[k][i]);
                            }

                            for (int k = 0; k < cat.size; k++) {
                                if (maps[k].get(a.sampleName.get(i)) > 0) {
                                    grps[mi[i]].inc(cat.get(k));
                                    proj.write("\t1");
                                } else {
                                    proj.write("\t0");
                                }
                            }

                            proj.writeLine();
                        }
                    }
                }
            }

            proj.close();
        }
    }*/

/*--------------------------------------------------------------------
 project - Project a datset by W
--------------------------------------------------------------------*/

    public Array project(String fname) {
        Array w = new Array(fname);
        return project(w);
    }

/*--------------------------------------------------------------------
 project - Project an array
--------------------------------------------------------------------*/

    public Array project(Array W) {
        Texts Z = geneName.match(W.geneName);
        Z.sort();

        W.selectGene(Z);
        selectGene(Z);

        /*
        com.scanfeld.io.Outfile g1 = new com.scanfeld.io.Outfile ( "c:/temp/g1.txt" );
        for ( int i = 0; i < W.geneName.size; i++ )
        {
          g1.writeLine ( W.geneName.get ( i ) );
        }
        g1.close ();

        com.scanfeld.io.Outfile g2 = new com.scanfeld.io.Outfile ( "c:/temp/g2.txt" );
        for ( int i = 0; i < geneName.size; i++ )
        {
          g2.writeLine ( geneName.get ( i ) );
        }
        g2.close ();

        com.scanfeld.io.Outfile g3 = new com.scanfeld.io.Outfile ( "c:/temp/g3.txt" );
        for ( int i = 0; i < Z.size; i++ )
        {
          g3.writeLine ( Z.get ( i ) );
        }
        g3.close ();
        */

        Matrix Wi = W.x.ninverse();
        Matrix n = Wi.nmatmul(x);
        // out( freemem() + " " + totalmem() + " " + maxmem() );

        Array a = new Array();
        a.x = n;
        a.sampleName = new Texts(sampleName);
        a.geneName = new Texts();
        a.gene = new Gene[Wi.nr];

        for (int i = 0; i < Wi.nr; i++) {
            a.geneName.rpush("F" + (i + 1));
            a.gene[i] = new Gene(s2t("F" + (i + 1)), s2t(""), a);
        }

        a.numSample = a.sampleName.size;
        a.numGene = a.geneName.size;

        a.numPhenotype = numPhenotype;
        a.phenotypeName = phenotypeName.copy();
        a.phenotype = new Phenotype[a.numPhenotype];

        a.numSample = numSample;
        a.sampleName = sampleName.copy();
        a.sample = new Sample[a.numSample];

        // LOOP THROUGH THE SAMPLES
        for (int i = 0; i < a.numSample; i++) {
            a.sample[i] = new Sample(a.sampleName.get(i), i, a);  // Make the ith sample
            a.sample[i].phenotype = sample[i].phenotype.copy();
        }

        // LOOP THROUGH THE PHENOTYPES
        for (int i = 0; i < a.numPhenotype; i++) {
            a.phenotype[i] = new Phenotype(a.phenotypeName.get(i), a);  // Make the ith phenotype
        }

        return a;
    }
    
    public Array projectFixedW(Array W, double precision, NMFCostFunction cost) {
        Texts Z = geneName.match(W.geneName);
        Z.sort();

        W.selectGene(Z);
        selectGene(Z);

        int k=W.x.nc; 
        Matrix WMatrix=W.x;
        
        System.err.println(Z.size);
        System.err.println(WMatrix.nr+" "+WMatrix.nc);
        System.err.println(x.nr+" "+x.nc);
        
        x.NMFdiv(k, cost, precision, WMatrix);
        return this;
    }

/*--------------------------------------------------------------------
 homolog - com.scanfeld.pr.Map to homolog file
--------------------------------------------------------------------*/

    public void homolog(String fname) {
        /*
        com.scanfeld.io.Infile a = new com.scanfeld.io.Infile ();
        a.load ( fname );

        com.scanfeld.core.Text[] b = a.splitLine ();
        com.scanfeld.core.Text[] y = new com.scanfeld.core.Text[ b.length ];
        com.scanfeld.core.Text[] z = new com.scanfeld.core.Text[ b.length ];

        for ( int i = 0; i < b.length; i++ )
        {
          com.scanfeld.core.Text[] c = b[ i ].splitTab ();
          y[ i ] = c[ 0 ];
          z[ i ] = c[ 1 ];
        }

        Array na = new Array();
        na.numGene = b.length;
        na.numSample = numSample;

        na.x = new float[ na.numGene ][ na.numSample ];

        int m = 0;
        int n = 0;

        while( m < b.length || n <for ( int i = 0; i < b.length; i++ )
        {
          if( z[ i ].eq( gene[ i ].name ) )
          {

          }
        }

        out ( y[ 0 ] + " " + z[ 0 ] );
        */
    }

/*--------------------------------------------------------------------
 loadCLS - Load the phenotype information from a cls file
--------------------------------------------------------------------*/

    public void homolog(int dir, String curOrg, String newOrg, int a, int b, String fname) {
        Url u;                                  // Hold the webpages
        Text t;                                 // Temporary com.scanfeld.core.Text object
        Text[] s;                               // Temporary com.scanfeld.core.Text objects
        String newGene;                         // Temporary gene for the new organism
        Outfile f = new Outfile(fname);      // Open the output file

        // WRITE THE FILE HEADER
        f.writeLine("index\t" + curOrg + "\tdesc\t" + newOrg + "\tdesc\tlen\tSW-score\tmargin\tbits\tidentity\toverlap\tbest");

        // LOOP THROUGH THE GENES
        for (int i = a; i < b; i++) {
            String test = "best_best";

            if (dir < 0) {
                test = "rev_best";
            } else if (dir > 0) {
                test = "best";
            }

            u = new Url("http://www.genome.jp/ssdb-bin/ssdb_" + test + "?org_gene=" + curOrg + ":" + gene[i].name);

            t = u.findBetween(newOrg + ":", "<INPUT");

            // IF NO MATCH WAS FOUND
            if (t.length() > 0) {
                t = t.remove(new String[]{"\"", "</A>"});
                t = t.remove(new String[]{"-", "log>"});
                t = t.replace("</a>", "");
                t = t.replace("(", "");
                t = t.replace(")", "");
                t = t.replace("&lt;", "");
                t = t.replace("&gt;", "");
                t = t.removeExtraWhite();
                s = t.splitSpace();

                t.set("");
                t.rpush(txt(s, 0));
                t.rpush("\t");
                t.rpush(txt(s, 1));

                for (int j = 2; j < s.length; j++) {
                    if (j < s.length - 7) {
                        t.rpush(" ");
                        t.rpush(s[j]);
                    } else {
                        t.rpush("\t");
                        t.rpush(s[j]);
                    }
                }

                f.writeLine((i + 1) + "\t" + gene[i].name + "\t" + gene[i].description + "\t" + t);
            }
        }

        f.close();  // Close the output file
    }

/*--------------------------------------------------------------------
 square - Set the phenotype shapes to square
--------------------------------------------------------------------*/

    public void square() {
        for (int i = 0; i < numPhenotype; i++) {
            phenotype[i].shape = "square";
        }
    }

/*--------------------------------------------------------------------
 circle - Set the phenotype shapes to circle
--------------------------------------------------------------------*/

    public void circle() {
        for (int i = 0; i < numPhenotype; i++) {
            phenotype[i].shape = "circle";
        }
    }

 
 
}
