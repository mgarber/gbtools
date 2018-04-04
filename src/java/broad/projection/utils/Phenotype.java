package broad.projection.utils;

import broad.projection.nmf.Array;


/*--------------------------------------------------------------------
 com.scanfeld.bio.Phenotype
 ---------------------------------------------------------------------
 Models a com.scanfeld.bio.Phenotype in a microarray experiment
 ---------------------------------------------------------------------
 Author: Daniel Scanfeld
   Date: 08.23.06
 ---------------------------------------------------------------------

--------------------------------------------------------------------*/

public class Phenotype extends Inc
{
/*--------------------------------------------------------------------
 Class Variables
--------------------------------------------------------------------*/

  public Array array;  // com.scanfeld.bio.Array the com.scanfeld.bio.Phenotype belongs to
  public int index;    // Index in the array
  public int size;     // Number of samples
  public Text name;    // Name of the phenotype
  public double[] c;
  public String shape; // Shape for plots

/*--------------------------------------------------------------------
 CONSTRUCTOR FUNCTIONS
 ---------------------------------------------------------------------
 Initialize com.scanfeld.bio.Phenotype objects
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 com.scanfeld.bio.Phenotype - Null constructor
--------------------------------------------------------------------*/

  public Phenotype ()
  {
    c = new double[] { 1.0, 0.0, 0.0 };
    shape = "circle";
  }

/*--------------------------------------------------------------------
 com.scanfeld.bio.Phenotype - com.scanfeld.core.Text constructor
--------------------------------------------------------------------*/

  public Phenotype ( Text name, Array array )
  {
    this.name = name;
    this.array = array;
    shape = "circle";
  }
}