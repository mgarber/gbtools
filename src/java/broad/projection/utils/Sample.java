package broad.projection.utils;

import broad.projection.nmf.Array;



/*--------------------------------------------------------------------
 TODO
 ---------------------------------------------------------------------
 update index in the array when samples get sorted
 ---------------------------------------------------------------------
 IDEAS
 ---------------------------------------------------------------------
 ---------------------------------------------------------------------
 BUGS
 ---------------------------------------------------------------------
--------------------------------------------------------------------*/

/**
 * A sample in a microarray experiment
 *
 * @author Daniel Scanfeld
 * @version 1.0
 */

public class Sample extends Inc
{
/*--------------------------------------------------------------------
 Class Variables
--------------------------------------------------------------------*/

  /** Array the sample belongs to */
  public Array array;

  /** Index in the array */
  public int index;

  /** Sample name */
  public Text name;

  /** Sample abbreviation */
  public String abbr;

  /** Shape for plots */
  public String shape;

  /** Color for plots */
  public double[] color;

  /** Set of phenotypes */
  public Texts phenotype;

/*--------------------------------------------------------------------
 CONSTRUCTOR FUNCTIONS
--------------------------------------------------------------------*/

  /** Null constructor */

  public Sample ()
  {
    init ( new Text ( "" ), 0, null, "", null, "" );  // Initialize the sample
  }

  /**
   * Full constructor
   *
   * @param name  Name of the sample
   * @param index Index in the Array
   * @param array The Array
   */

  public Sample ( Text name, int index, Array array )
  {
    init ( name, index, array, "", null, "" );
  }

  /**
   * Initializes the object
   *
   * @param name  Name of the sample
   * @param index Index in the Array
   * @param array The Array
   */

  public void init ( Text name, int index, Array array, String abbr, double[] color, String shape )
  {
    this.name = name;    // Set the name
    this.index = index;  // Set the index
    this.array = array;  // Set the Array
    this.abbr = abbr;    // Set the abbreviation
    this.color = color;  // Set the color
    this.shape = shape;  // Set the shape

    parseName ( name );
    phenotype = new Texts ();
  }

  /**
   * Parses a name for other variables
   *
   * @param name Name string from the gct file
   */

  public void parseName ( Text name )
  {
    Text[] sampleArr = name.split ( "'" );  // Split the name string
    this.name = sampleArr[ 0 ];             // Set the name

    // IF AN ABBREVIATION WAS INCLUDED
    if ( sampleArr.length > 1 )
    {
      abbr = sampleArr[ 1 ].str ();  // Set the abbreviation
    }

    // IF A SHAPE WAS INCLUDED
    if ( sampleArr.length > 2 )
    {
      shape = sampleArr[ 2 ].str ();  // Set the shape
    }

    // IF A COLOR WAS INCLUDED
    if ( sampleArr.length > 4 )
    {
      // IF SETTING AN HSV COLOR
      if ( sampleArr[ 3 ].str ().equals ( "h" ) )
      {
        double h = sampleArr[ 4 ].dbl ();  // Set the hue
        double s = 1.0;                    // Initial saturation
        double v = 1.0;                    // Initial value

        // IF THE SATURATION WAS INCLUDED
        if ( sampleArr.length > 5 )
        {
          s = sampleArr[ 5 ].dbl ();  // Set the saturation
        }

        // IF THE VALUE WAS INCLUDED
        if ( sampleArr.length > 6 )
        {
          v = sampleArr[ 6 ].dbl ();  // Set the value
        }

        color = hsv ( h, s, v );  // Set the color
      }

      // IF SETTING AN RGB COLOR
      else
      {
        double r = sampleArr[ 4 ].dbl ();  // Set the red
        double g = r;                      // Initial green
        double b = r;                      // Initial blue

        // IF THE GREEN WAS INCLUDED
        if ( sampleArr.length > 5 )
        {
          g = sampleArr[ 5 ].dbl ();  // Set the green
        }

        // IF THE BLUE WAS INCLUDED
        if ( sampleArr.length > 6 )
        {
          b = sampleArr[ 5 ].dbl ();  // Set the blue
        }

        color = rgb ( r, g, b );  // Set the color
      }
    }
  }

  /**
   * Makes a name for the GCT file
   *
   * @return string filled with the sample variables
   */

  public String makeName ()
  {
    String z = "";  // Name string
    z += "'";       // Add a delimiter

    // IF AN ABBREVIATION EXISTS
    if ( abbr != null )
    {
      z += abbr; // Add the abbreviation
    }

    z += "'";  // Add a delimiter

    // IF A SHAPE EXISTS
    if ( shape != null )
    {
      z += shape;  // Add the shape
    }

    z += "'";  // Add a delimiter

    // IF A COLOR EXISTS
    if ( color != null )
    {
      z += "r'";              // Specify RGB
      z += color[ 0 ] + "'";  // Add red
      z += color[ 1 ] + "'";  // Add green
      z += color[ 2 ];        // Add blue
    }

    // IF NO PARAMETERS WERE ADDED
    if ( z.length () == 3 )
    {
      z = "";  // Remove the suffix
    }

    z = name + z;  // Add the name

    return z;  // Return the name string
  }

  /**
   * Checks if this sample has this phenotype
   *
   * @param p Phenotype to search for
   *
   * @return true if the phenotype matches, false if it does not match
   */

  public boolean isPhenotype ( Text p )
  {
    // IF THE PHENOTYPE IS IN THE LIST
    if ( phenotype.find ( p ) >= 0 )
    {
      return true;
    }

    // IF THE PHENOTYPE IS NOT IN THE LIST
    else
    {
      return false;
    }
  }

  /**
   * Adds a phenotype to the list
   *
   * @param p Phenotype name
   */

  public void addPhenotype ( Text p )
  {
    // IF IT IS NOT THE PHENOTYPE ALREADY
    if ( !isPhenotype ( p ) )
    {
      phenotype.rpush ( p );  // Add the phenotype
    }
  }
}