package broad.projection.utils;

import broad.projection.nmf.Array;



/*--------------------------------------------------------------------
 com.scanfeld.bio.Gene
 ---------------------------------------------------------------------
 Models a gene in a microarray experiment
 ---------------------------------------------------------------------
 Author: Daniel Scanfeld
   Date: 08.23.06
 ---------------------------------------------------------------------

--------------------------------------------------------------------*/

public class Gene extends Inc
{
/*--------------------------------------------------------------------
 Class Variables
--------------------------------------------------------------------*/

  public Array array;       // com.scanfeld.bio.Array the gene belongs to
  public int index;         // Index in the array
  public Text name;         // Name
  public Text description;  // Description

/*--------------------------------------------------------------------
 CONSTRUCTOR FUNCTIONS
 ---------------------------------------------------------------------
 Initialize com.scanfeld.bio.Gene objects
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 com.scanfeld.bio.Gene - Null constructor
--------------------------------------------------------------------*/

  public Gene ()
  {
  }

/*--------------------------------------------------------------------
 com.scanfeld.bio.Gene - com.scanfeld.core.Text constructor
--------------------------------------------------------------------*/

  public Gene ( Text name, Text description, Array array )
  {
    this.name = name;
    this.description = description;
    this.array = array;
  }
}