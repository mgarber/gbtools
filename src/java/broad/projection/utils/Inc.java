package broad.projection.utils;/*--------------------------------------------------------------------
 com.scanfeld.core.Inc
 ---------------------------------------------------------------------
 Helper functions for all java programs
 ---------------------------------------------------------------------
 Author: Daniel Scanfeld
   Date: 08.23.06
 ---------------------------------------------------------------------
 Elaborate Date Output
 Function to add zeros to the left
 Search through code and trim file loads when necessary
 Get com.scanfeld.core.RC4 implemented
 Get com.scanfeld.math.Matrix implemented as well
 pdf and ps output
 log load and save, etc...
 make a log file
 random matrix function
--------------------------------------------------------------------*/

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Calendar;
import java.util.Enumeration;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import broad.projection.math.Matrix;


public class Inc
{
/*--------------------------------------------------------------------
 Class Variables
--------------------------------------------------------------------*/

  public static RC4 rc;              // Random number generator
  public static long seconds;        //
  public static long delay;          //
  public static boolean statusFlag;  //
  public static boolean dbgFlag;     // Debug flag
  public static int bufSize;         // Size of the buffer
  public static int indentAmount;    // Number of spaces to indent
  public static int dashLength;      // Dash length

  public static Texts unzipped;      // Unzipped Files
  public static Texts created;       // Created Files

/*--------------------------------------------------------------------
 RANDOM FUNCTIONS
 ---------------------------------------------------------------------
 Uses the com.scanfeld.core.RC4 random number class
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 initInc - Initialize the com.scanfeld.core.Inc variables
--------------------------------------------------------------------*/

  public static void initInc ()
  {
    rc = new RC4 ();   // Initialize the random number generator
    dbgFlag = true;    // Debug by default
    delay = 1;         // Set the status delay
    bufSize = 1024;    // Buffer size
    indentAmount = 1;  // Number of spaces to indent
    dashLength = 40;   // Default dash length
    startStatus ();     // Start the status checker

    unzipped = new Texts();
    created = new Texts();
  }

/*--------------------------------------------------------------------
 rand - Return a random byte
--------------------------------------------------------------------*/

  public char rand ()
  {
    return rc.rand ();
  }

/*--------------------------------------------------------------------
 randString - Return a random string
--------------------------------------------------------------------*/

  public String randString ( int size )
  {
    return rc.randString ( size );
  }

/*--------------------------------------------------------------------
 randNormal - Create a random independent normally distributed number
   with mean m and standard deviation s
--------------------------------------------------------------------*/

  public double randNormal ( double m, double s )
  {
    return rc.randNormal ( m, s );
  }

/*--------------------------------------------------------------------
 randNormal - Create an array of normal doubles of size n with mean m
   and standard deviation s
--------------------------------------------------------------------*/

  public double[] randNormal ( double m, double s, int n )
  {
    return rc.randNormal ( m, s, n );
  }

/*--------------------------------------------------------------------
 randBoolean- Create a random boolean
--------------------------------------------------------------------*/

  public boolean randBoolean ()
  {
    return rc.randDouble () < 0.5;
  }

/*--------------------------------------------------------------------
 randBoolean- Create a random boolean array
--------------------------------------------------------------------*/

  public boolean[] randBoolean ( int n )
  {
    boolean[] z = new boolean[n];

    for ( int i = 0; i < n; i++ )
    {
      z[ i ] = rc.randDouble () < 0.5;
    }

    return z;
  }

/*--------------------------------------------------------------------
 randDouble- Create a random double between 0 and 1.  doubleSize is
   the nuber of bytes used to encode the double.
--------------------------------------------------------------------*/

  public double randDouble ()
  {
    //return rc.randDouble ();
	  return Math.random();
  }

/*--------------------------------------------------------------------
 randDouble - Create a random double array of size n between 0 and 1
--------------------------------------------------------------------*/

  public double[] randDouble ( int n )
  {
    return rc.randDouble ( n );
  }

/*--------------------------------------------------------------------
 randInt - Create a random integer between min and max
--------------------------------------------------------------------*/

  public int randInt ( int min, int max )
  {
    return rc.randInt ( min, max );
  }

/*--------------------------------------------------------------------
 randInt - Create a random array of integers 0 to k - 1
--------------------------------------------------------------------*/

  public int[] randInt ( int k )
  {
    return rc.randInt ( k );
  }

/*--------------------------------------------------------------------
 randIntArray - Create a random array of n integers 0 to k - 1
--------------------------------------------------------------------*/

  public int[] randIntArray ( int k, int n )
  {
    return rc.randIntArray ( k, n );
  }

/*--------------------------------------------------------------------
 randIntArray - Create a random array of n integers 0 to k - 1
--------------------------------------------------------------------*/

  public int[] randIntArray ( int k, int n, boolean f )
  {
    if ( f )
    {
      return rc.randInt ( 0, k - 1, n );
    }

    else
    {
      return rc.randIntArray ( k, n );
    }
  }

/*--------------------------------------------------------------------
 randInt - Create a random array of n integers min to max
--------------------------------------------------------------------*/

  public int[] randInt ( int min, int max, int n )
  {
    return rc.randInt ( min, max, n );
  }

/*--------------------------------------------------------------------
 SEQUENCE FUNCTIONS
 ---------------------------------------------------------------------
 Manipulates sequences
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 seq - Create a sequence from a to b
--------------------------------------------------------------------*/

  public int[] seq ( int a, int b )
  {
    int[] z;

    if ( b < a )
    {
      z = new int[a - b + 1];

      // LOOP THROUGH THE SEQUENCE
      for ( int i = 0; i < a - b + 1; i++ )
      {
        z[ i ] = a - i;
      }
    }

    else
    {
      z = new int[b - a + 1];

      // LOOP THROUGH THE SEQUENCE
      for ( int i = 0; i < b - a + 1; i++ )
      {
        z[ i ] = i + a;
      }
    }

    return z;
  }

/*--------------------------------------------------------------------
 seq - Create a sequence from a to b, incrementing by k
--------------------------------------------------------------------*/

  public int[] seq ( int a, int b, int k )
  {
    int[] z;

    if ( b < a )
    {
      z = new int[( a - b ) / k + 1];

      // LOOP THROUGH THE SEQUENCE
      for ( int i = 0; i < a - b + 1; i += k )
      {
        z[ i / k ] = a - i;
      }
    }

    else
    {
      z = new int[( b - a ) / k + 1];

      // LOOP THROUGH THE SEQUENCE
      for ( int i = 0; i < b - a + 1; i += k )
      {
        z[ i / k ] = i + a;
      }
    }

    return z;
  }

/*--------------------------------------------------------------------
 angle - Calculate the angle of a set of points about the origin
--------------------------------------------------------------------*/

  public double[] angle ( double[][] z )
  {
    double[] a = new double[z.length];

    for ( int i = 0; i < z.length; i++ )
    {
      a[ i ] = angle ( z[ i ][ 0 ], z[ i ][ 1 ] );
    }

    return a;
  }

/*--------------------------------------------------------------------
 intersect - Find the intersection of two lines
--------------------------------------------------------------------*/

  public double[] intersect ( double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4 )
  {
    double[] a = new double[2];

    double ua = ( ( x4 - x3 ) * ( y1 - y3 ) - ( y4 - y3 ) * ( x1 - x3 ) ) / ( ( y4 - y3 ) * ( x2 - x1 ) - ( x4 - x3 ) * ( y2 - y1 ) );
    a[ 0 ] = x1 + ua * ( x2 - x1 );
    a[ 1 ] = y1 + ua * ( y2 - y1 );

    return a;
  }

/*--------------------------------------------------------------------
 angle - Calculate the angle of a point about the origin
--------------------------------------------------------------------*/

  public double angle ( double x, double y )
  {
    double a = 0;

    if ( x == 0 )
    {
      if ( y < 0 )
      {
        a = 3.0 * Math.PI / 2.0;
      }

      else
      {
        a = Math.PI / 2.0;
      }
    }

    else
    {
      a = Math.atan ( abs ( y ) / abs ( x ) );

      if ( x < 0 )
      {
        if ( y < 0 )
        {
          a = Math.PI + a;
        }

        else
        {
          a = Math.PI - a;
        }
      }

      else
      {
        if ( y < 0 )
        {
          a = 2.0 * Math.PI - a;
        }

        else
        {
          a = a;
        }
      }
    }

    if ( a < 0 )
    {
      a += Math.PI * 2;
    }

    if ( a > Math.PI * 2 )
    {
      a -= Math.PI * 2;
    }

    return a;
  }

/*--------------------------------------------------------------------
 seq - Create a sequence from a to b
--------------------------------------------------------------------*/

  public char[] seq ( char a, char b )
  {
    char[] z;

    if ( b < a )
    {
      z = new char[a - b + 1];

      // LOOP THROUGH THE SEQUENCE
      for ( int i = 0; i < a - b + 1; i++ )
      {
        z[ i ] = ( char ) ( a - i );
      }
    }

    else
    {
      z = new char[b - a + 1];

      // LOOP THROUGH THE SEQUENCE
      for ( int i = 0; i < b - a + 1; i++ )
      {
        z[ i ] = ( char ) ( a + i );
      }
    }

    return z;
  }

/*--------------------------------------------------------------------
 entropy - Find the entropy of a set of doubles
--------------------------------------------------------------------*/

  public double entropy ( double[] x )
  {
    double h = 0;

    // LOOP THROUGH THE SEQUENCE
    for ( int i = 0; i < x.length; i++ )
    {
      if ( x[ i ] > 0 )
      {
        h += x[ i ] * log2 ( x[ i ] );
      }
    }

    return -h;
  }

/*--------------------------------------------------------------------
 entropy - Find the entropy of a set of ints
--------------------------------------------------------------------*/

  public double entropy ( int[] y )
  {
    double[] x = new double[y.length];
    double sum = sum ( y );

    // LOOP THROUGH THE SEQUENCE
    for ( int i = 0; i < x.length; i++ )
    {
      x[ i ] = ( double ) ( y[ i ] ) / sum;
    }

    return entropy ( x );
  }

/*--------------------------------------------------------------------
 entropy01 - Find the entropy of a set of doubles
--------------------------------------------------------------------*/

  public double entropy01 ( double[] x )
  {
    double h = 0;

    // LOOP THROUGH THE SEQUENCE
    for ( int i = 0; i < x.length; i++ )
    {
      if ( x[ i ] > 0 )
      {
        h += x[ i ] * log2 ( x[ i ] );
      }
    }

    double rnd = 1.0 / ( double ) ( x.length );
    h /= ( double ) ( x.length ) * rnd * log2 ( rnd );

    return h;
  }

/*--------------------------------------------------------------------
 entropy01 - Find the entropy of a set of ints from 0 to 1
--------------------------------------------------------------------*/

  public double entropy01 ( int[] y )
  {
    double[] x = new double[y.length];
    double sum = sum ( y );

    // LOOP THROUGH THE SEQUENCE
    for ( int i = 0; i < x.length; i++ )
    {
      x[ i ] = ( double ) ( y[ i ] ) / sum;
    }

    return entropy01 ( x );
  }

/*--------------------------------------------------------------------
 rev - Reverse an int array
--------------------------------------------------------------------*/

  public int[] rev ( int[] x )
  {
    int[] z = new int[x.length];

    // LOOP THROUGH THE SEQUENCE
    for ( int i = 0; i < x.length; i++ )
    {
      z[ x.length - i - 1 ] = x[ i ];
    }

    return z;
  }

/*--------------------------------------------------------------------
 rev - Reverse a double array
--------------------------------------------------------------------*/

  public double[] rev ( double[] x )
  {
    double[] z = new double[x.length];

    // LOOP THROUGH THE SEQUENCE
    for ( int i = 0; i < x.length; i++ )
    {
      z[ x.length - i - 1 ] = x[ i ];
    }

    return z;
  }

/*--------------------------------------------------------------------
 rev - Reverse an int array
--------------------------------------------------------------------*/

  public char[] rev ( char[] x )
  {
    char[] z = new char[x.length];

    // LOOP THROUGH THE SEQUENCE
    for ( int i = 0; i < x.length; i++ )
    {
      z[ x.length - i - 1 ] = x[ i ];
    }

    return z;
  }

/*--------------------------------------------------------------------
 sel - Select a set of integers
--------------------------------------------------------------------*/

  public int[] sel ( int[] x, int[] y )
  {
    int[] z = new int[y.length];

    // LOOP THROUGH THE SEQUENCE
    for ( int i = 0; i < y.length; i++ )
    {
      z[ i ] = x[ y[ i ] ];
    }

    return z;
  }

/*--------------------------------------------------------------------
 sel - Select a set of integers
--------------------------------------------------------------------*/

  public int[][] sel ( int[][] x, int[] y )
  {
    int[][] z = new int[x.length][y.length];

    // LOOP THROUGH THE SEQUENCE
    for ( int i = 0; i < y.length; i++ )
    {
      z[ i ] = x[ y[ i ] ];
    }

    return z;
  }

/*--------------------------------------------------------------------
 sel - Select a set of doubles
--------------------------------------------------------------------*/

  public double[] sel ( double[] x, int[] y )
  {
    double[] z = new double[y.length];

    // LOOP THROUGH THE SEQUENCE
    for ( int i = 0; i < y.length; i++ )
    {
      z[ i ] = x[ y[ i ] ];
    }

    return z;
  }

/*--------------------------------------------------------------------
 sel - Select a set of doubles
--------------------------------------------------------------------*/

  public double[][] sel ( double[][] x, int[] y )
  {
    double[][] z = new double[x.length][y.length];

    // LOOP THROUGH THE SEQUENCE
    for ( int i = 0; i < y.length; i++ )
    {
      z[ i ] = x[ y[ i ] ];
    }

    return z;
  }

/*--------------------------------------------------------------------
 sel - Select a set of booleans
--------------------------------------------------------------------*/

  public boolean[] sel ( boolean[] x, int[] y )
  {
    boolean[] z = new boolean[y.length];

    // LOOP THROUGH THE SEQUENCE
    for ( int i = 0; i < y.length; i++ )
    {
      z[ i ] = x[ y[ i ] ];
    }

    return z;
  }

/*--------------------------------------------------------------------
 sel - Select a set of com.scanfeld.core.Texts
--------------------------------------------------------------------*/

  public Texts[] sel ( Texts[] x, int[] y )
  {
    Texts[] z = new Texts[y.length];

    // LOOP THROUGH THE SEQUENCE
    for ( int i = 0; i < y.length; i++ )
    {
      z[ i ] = x[ y[ i ] ];
    }

    return z;
  }

/*--------------------------------------------------------------------
 sel - Select a set of characters
--------------------------------------------------------------------*/

  public char[] sel ( char[] x, int[] y )
  {
    char[] z = new char[y.length];

    // LOOP THROUGH THE SEQUENCE
    for ( int i = 0; i < y.length; i++ )
    {
      z[ i ] = x[ y[ i ] ];
    }

    return z;
  }

/*--------------------------------------------------------------------
 sel - Select a set of strings
--------------------------------------------------------------------*/

  public String[] sel ( String[] x, int[] y )
  {
    String[] z = new String[y.length];

    // LOOP THROUGH THE SEQUENCE
    for ( int i = 0; i < y.length; i++ )
    {
      z[ i ] = x[ y[ i ] ];
    }

    return z;
  }

/*--------------------------------------------------------------------
 copy - Copy an int array into another one
--------------------------------------------------------------------*/

  public void copy ( int[] a, int ai, int[] b, int bi )
  {
    copy ( a, ai, b, bi, b.length - bi );
  }

/*--------------------------------------------------------------------
 copy - Copy an int array into another one
--------------------------------------------------------------------*/

  public void copy ( int[] a, int ai, int[] b, int bi, int len )
  {
    int al = a.length;
    int bl = b.length;

    for ( int i = 0; i < len; i++ )
    {
      if ( ai + i < al && bi + i < bl )
      {
        a[ ai + i ] = b[ bi + i ];
      }
    }
  }

/*--------------------------------------------------------------------
 booleanFalse - Create an array of boolean falses
--------------------------------------------------------------------*/

  public boolean[] booleanFalse ( int a )
  {
    boolean[] b = new boolean[a];

    for ( int i = 0; i < a; i++ )
    {
      b[ i ] = false;
    }

    return b;
  }

/*--------------------------------------------------------------------
 booleanTrue - Create an array of boolean trues
--------------------------------------------------------------------*/

  public boolean[] booleanTrue ( int a )
  {
    boolean[] b = new boolean[a];

    for ( int i = 0; i < a; i++ )
    {
      b[ i ] = true;
    }

    return b;
  }

/*--------------------------------------------------------------------
 OUTPUT FUNCTIONS
 ---------------------------------------------------------------------
 Handle output to the screen
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 out - Output a string
--------------------------------------------------------------------*/

  public void out ( String z )
  {
    System.out.println ( z );  // Write to the system output function
  }

/*--------------------------------------------------------------------
 outDash - Output a string with dashes
--------------------------------------------------------------------*/

  public void outDash ( String z )
  {
    outDash ();
    out ( z );
    outDash ();
  }

/*--------------------------------------------------------------------
 outDash - Output a dash
--------------------------------------------------------------------*/

  public void outDash ()
  {
    out ( dash () );
  }

/*--------------------------------------------------------------------
 out - Output a new line
--------------------------------------------------------------------*/

  public void out ()
  {
    out ( "" );  // Write a new line
  }

/*--------------------------------------------------------------------
 out - Output a character
--------------------------------------------------------------------*/

  public void out ( char z )
  {
    out ( "" + z );
  }

/*--------------------------------------------------------------------
 out - Output a character array
--------------------------------------------------------------------*/

  public void out ( char[] z )
  {
    out ( c2s ( z ) );
  }

/*--------------------------------------------------------------------
 out - Output a boolean
--------------------------------------------------------------------*/

  public void out ( boolean z )
  {
    // If true
    if ( z )
    {
      out ( "true" );  // Output true
    }

    // If false
    else
    {
      out ( "false" );  // Output false
    }
  }

/*--------------------------------------------------------------------
 out - Output an integer
--------------------------------------------------------------------*/

  public void out ( int i )
  {
    out ( "" + i );
  }

/*--------------------------------------------------------------------
 out - Output a double array
--------------------------------------------------------------------*/

  public void out ( double[] z )
  {
    out ( "Double Array (" + z.length + ")" );

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < z.length; i++ )
    {
      out ( i + ": " + z[ i ] );  // Output the ith element
    }

    out ( "" );  // Add a blank line at the end
  }

/*--------------------------------------------------------------------
 showOne - Output an int array
--------------------------------------------------------------------*/

  public String showOne ( int[] z )
  {
    String a = new String ( "( " );

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < z.length; i++ )
    {
      if ( i > 0 )
      {
        a += ", ";
      }

      a += z[ i ];  // Output the ith element
    }

    a += " )";

    return a;
  }

/*--------------------------------------------------------------------
 showOne - Output a double array
--------------------------------------------------------------------*/

  public String showOne ( double[] z )
  {
    String a = new String ( "( " );

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < z.length; i++ )
    {
      if ( i > 0 )
      {
        a += ", ";
      }

      a += z[ i ];  // Output the ith element
    }

    a += " )";

    return a;
  }

/*--------------------------------------------------------------------
 showOne - Output a String array
--------------------------------------------------------------------*/

  public String showOne ( String[] z )
  {
    String a = new String ( "( " );

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < z.length; i++ )
    {
      if ( i > 0 )
      {
        a += ", ";
      }

      a += z[ i ];  // Output the ith element
    }

    a += " )";

    return a;
  }

/*--------------------------------------------------------------------
 out - Output a boolean array
--------------------------------------------------------------------*/

  public void out ( boolean[] z )
  {
    out ( "Boolean Array (" + z.length + ")" );

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < z.length; i++ )
    {
      out ( i + ": " + z[ i ] );  // Output the ith element
    }

    out ( "" );  // Add a blank line at the end
  }

/*--------------------------------------------------------------------
 out - Output an integer array
--------------------------------------------------------------------*/

  public void out ( int[] z )
  {
    out ( "Integer Array (" + z.length + ")" );

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < z.length; i++ )
    {
      out ( i + ": " + z[ i ] );  // Output the ith element
    }

    out ( "" );  // Add a blank line at the end
  }

/*--------------------------------------------------------------------
 out - Output a Text array
--------------------------------------------------------------------*/

  public void out ( String[] z )
  {
    out ( "String Array (" + z.length + ")" );

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < z.length; i++ )
    {
      out ( i + ": " + z[ i ] );  // Output the ith element
    }

    out ( "" );  // Add a blank line at the end
  }

/*--------------------------------------------------------------------
 out - Output a Text array
--------------------------------------------------------------------*/

  public void out ( Text[] z )
  {
    out ( "Text Array (" + z.length + ")" );

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < z.length; i++ )
    {
      out ( i + ": " + z[ i ] );  // Output the ith element
    }

    out ( "" );  // Add a blank line at the end
  }

/*--------------------------------------------------------------------
 out - Output a double
--------------------------------------------------------------------*/

  public void out ( double d )
  {
    out ( "" + d );
  }

/*--------------------------------------------------------------------
 out - Output text
--------------------------------------------------------------------*/

  public void out ( Text t )
  {
    out ( "" + t );
  }

/*--------------------------------------------------------------------
 out - Output text
--------------------------------------------------------------------*/

  public void out ( Texts t )
  {
    out ( "" + t );
  }

/*--------------------------------------------------------------------
 out - Output two texts as columns
--------------------------------------------------------------------*/

  public void out ( Texts a, Texts b )
  {
    int len = a.maxLength () + 1;

    for ( int i = 0; i < a.size && i < b.size; i++ )
    {
      out ( spaceRight ( a.get ( i ) + ":", len ) + " " + b.get ( i ) );
    }
  }

/*--------------------------------------------------------------------
 rep - Repeat a string n times
--------------------------------------------------------------------*/

  public String rep ( String s, int n )
  {
    String z = "";

    for ( int i = 0; i < n; i++ )
    {
      z += s;
    }

    return z;
  }

/*--------------------------------------------------------------------
 rep - Repeat a character n times
--------------------------------------------------------------------*/

  public String rep ( char s, int n )
  {
    String z = "";

    for ( int i = 0; i < n; i++ )
    {
      z += s;
    }

    return z;
  }

/*--------------------------------------------------------------------
 fillRight - Output with fills to the right
--------------------------------------------------------------------*/

  public String fillRight ( String s, char c, int n )
  {
    return s + rep ( c, n - s.length () );
  }

/*--------------------------------------------------------------------
 fillRight - Output with fills to the right
--------------------------------------------------------------------*/

  public Text fillRight ( Text s, char c, int n )
  {
    Text z = s.copy ();
    z.rpush ( rep ( c, n - s.length () ) );
    return z;
  }

/*--------------------------------------------------------------------
 fillLeft - Output with fills to the left
--------------------------------------------------------------------*/

  public String fillLeft ( String s, char c, int n )
  {
    return rep ( c, n - s.length () ) + s;
  }

/*--------------------------------------------------------------------
 fillLeft - Output with fills to the left
--------------------------------------------------------------------*/

  public Text fillLeft ( Text s, char c, int n )
  {
    Text z = s.copy ();
    z.lpush ( rep ( c, n - s.length () ) );
    return z;
  }

/*--------------------------------------------------------------------
 spaceRight - Output with spaces to the right
--------------------------------------------------------------------*/

  public String spaceRight ( String s, int n )
  {
    return fillRight ( s, ' ', n );
  }

/*--------------------------------------------------------------------
 spaceRight - Output with spaces to the right
--------------------------------------------------------------------*/

  public Text spaceRight ( Text s, int n )
  {
    return fillRight ( s, ' ', n );
  }

/*--------------------------------------------------------------------
 spaceLeft - Output with spaces to the left
--------------------------------------------------------------------*/

  public String spaceLeft ( String s, int n )
  {
    return fillLeft ( s, ' ', n );
  }

/*--------------------------------------------------------------------
 spaceLeft - Output with spaces to the left
--------------------------------------------------------------------*/

  public Text spaceLeft ( Text s, int n )
  {
    return fillLeft ( s, ' ', n );
  }

/*--------------------------------------------------------------------
 zeroFill - Fill zeros to the left of a number
--------------------------------------------------------------------*/

  public String zeroFill ( String s, int n )
  {
    return fillLeft ( s, '0', n );
  }

/*--------------------------------------------------------------------
 zeroFill - Fill zeros to the left of a number
--------------------------------------------------------------------*/

  public Text zeroFill ( Text s, int n )
  {
    return fillLeft ( s, '0', n );
  }

/*--------------------------------------------------------------------
 zeroFill - Fill zeros to the left of a number
--------------------------------------------------------------------*/

  public Text zeroFill ( int s, int n )
  {
    return fillLeft ( i2t ( s ), '0', n );
  }

/*--------------------------------------------------------------------
 out - Output matrix
--------------------------------------------------------------------*/

  public void out ( Matrix m )
  {
    m.out ();
  }

/*--------------------------------------------------------------------
 err - Output an error string
--------------------------------------------------------------------*/

  public void err ( String z )
  {
    out ( z );          // Output the string
    System.exit ( 0 );  // Exit the program
  }

/*--------------------------------------------------------------------
 exit - Exit the program
--------------------------------------------------------------------*/

  public void exit ()
  {
    System.exit ( 0 );  // Exit the program
  }

/*--------------------------------------------------------------------
 dbg - Output a debug string
--------------------------------------------------------------------*/

  public void dbg ( String z )
  {
    // IF DEBUGGING
    if ( dbgFlag )
    {
      out ( z );  // Output the string
    }
  }

/*--------------------------------------------------------------------
 sts - Output a status message
--------------------------------------------------------------------*/

  public void sts ( String z )
  {
    if ( seconds > 0 )
    {
      Calendar d = Calendar.getInstance ();  // A calendar instance
      long c = d.getTimeInMillis ();         // Current time

      // IF DELAY SECONDS HAVE TRANSPIRED
      if ( ( c - ( delay * 1000 ) ) > seconds )
      {
        out ( z );          // Output the message
        seconds = c;        // Save the current time
        statusFlag = true;  // Set the status flag
      }
    }
  }

/*--------------------------------------------------------------------
 startStatus - Start the status checker
--------------------------------------------------------------------*/

  public static void startStatus ()
  {
    Calendar d = Calendar.getInstance ();  // A calendar instance
    seconds = d.getTimeInMillis ();        // Save the current time
    statusFlag = false;                    // Unset the status flag
  }

/*--------------------------------------------------------------------
 endStatus - End the status checker
--------------------------------------------------------------------*/

  public static void endStatus ()
  {
    // IF A STATUS MESSAGE HAS BEEN SENT
    if ( statusFlag )
    {
      System.out.println ( "" );          // Output a new line
      statusFlag = false;  // Unset the status flag
    }
  }

/*--------------------------------------------------------------------
 dash - Create a dash
--------------------------------------------------------------------*/

  public String dash ( int len )
  {
    String z = "";

    for ( int i = 0; i < len; i++ )
    {
      z += "-";
    }

    return z;
  }

/*--------------------------------------------------------------------
 dash - Create a dash
--------------------------------------------------------------------*/

  public String dash ()
  {
    return dash ( dashLength );
  }

/*--------------------------------------------------------------------
 date - Return the date
--------------------------------------------------------------------*/

  public String date ()
  {
    Calendar d = Calendar.getInstance ();  // A calendar instance
    String z = zeroFill ( d.get ( Calendar.MONTH ), 2 ) + ".";
    z += zeroFill ( d.get ( Calendar.DAY_OF_MONTH ), 2 ) + ".";
    z += zeroFill ( d.get ( Calendar.YEAR ), 4 ) + " ";
    z += zeroFill ( d.get ( Calendar.HOUR_OF_DAY ), 2 ) + ":";
    z += zeroFill ( d.get ( Calendar.MINUTE ), 2 ) + ":";
    z += zeroFill ( d.get ( Calendar.SECOND ), 2 );
    return z;
  }

/*--------------------------------------------------------------------
 FILE FUNCTIONS
 ---------------------------------------------------------------------
 Ensure variables satisfy certain conditions
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 exists - Check if a file exists
--------------------------------------------------------------------*/

  public boolean exists ( String z )
  {
    File f = new File ( z );  // Create the file
    return f.exists ();
  }

/*--------------------------------------------------------------------
 delete - Delete a file
--------------------------------------------------------------------*/

  public void delete ( String z )
  {
    File f = new File ( z );  // Create the file

    // IF THE FILE EXISTS
    if ( f.exists () )
    {
      f.delete ();  // Delete the file
    }
  }

/*--------------------------------------------------------------------
 dir - Retrieve a directory listing
--------------------------------------------------------------------*/

  public String[] dir ( String fname )
  {
    String[] fn = null;  // Initialize a string array

    // TRY TO READ THE DIRECTORY
    try
    {
      File f = new File ( fname );  // Create a file object for the directory
      fn = f.list ();               // List the files
      f = null;                     // Set the file to null
    }

    // IF AN ERROR OCCURRED
    catch ( Exception e )
    {
      err ( "Error opening: \"" + fname + "\" for input." );  // Output an error message
    }

    return fn;  // Return the string array
  }

/*--------------------------------------------------------------------
 url - Load a url into a file
--------------------------------------------------------------------*/

  public void url ( String iname, String oname )
  {
    Url u = new Url ( iname, oname );
  }

/*--------------------------------------------------------------------
 upperFirst - Capitalize the first letter
--------------------------------------------------------------------*/

  public void upperFirst ( Text[] z )
  {
    for ( int i = 0; i < z.length; i++ )
    {
      z[ i ] = z[ i ].upperFirst ();
    }
  }

/*--------------------------------------------------------------------
 mkdirs - Make directories if they don't exist
--------------------------------------------------------------------*/

  public void mkdirs ( String fname )
  {
      File f;
      Text z = new Text ( fname );  // Create a text version of the filename
      z = z.replace ( "\\", "/" );  // Ensure only forward slashse
      int[] d = z.finds ( "/" );    // Find all forward slashes

      // If a directory was found
      if ( d.length > 0 )
      {
        String dir = z.sub ( 0, d[ d.length - 1 ] ).str ();  // Create the directory string
        f = new File ( dir );                                // Create the directory file

        // IF THE DIRECTORY DOES NOT EXIST
        if ( !f.exists () )
        {
          f.mkdirs ();  // Make all required directories to make it exist
        }
      }
  }


/*--------------------------------------------------------------------
 unzip - Unzip a zip file
--------------------------------------------------------------------*/

  public void unzip ( String iname, String oname )
  {
    try
    {
      ZipFile zf = new ZipFile(iname);
      Enumeration e = zf.entries();
      while (e.hasMoreElements()) {
	ZipEntry ze = (ZipEntry) e.nextElement();
        mkdirs( oname + ze.getName() );

        unzipped.rpush( oname + ze.getName() );
        if(!ze.isDirectory()) {
	 FileOutputStream fout = new FileOutputStream(oname + ze.getName());
	 InputStream in = zf.getInputStream(ze);
	 for (int c = in.read(); c != -1; c = in.read()) {
	   fout.write(c);
	 }
	 in.close();
	 fout.close();
        }
      }
    }

    // IF ANY ERROR HAS OCCURRED
    catch ( Exception e )
    {
      err ( e.getMessage () );  // Output the error message
    }
  }

/*--------------------------------------------------------------------
 deleteUnzipped - Delete the unzipped files
--------------------------------------------------------------------*/

  public void deleteUnzipped ( )
  {
    try
    {
      while( unzipped.size>0 )
      {
       String fname = t2s( unzipped.rpop() );
       File f = new File( fname );

       if(f.isDirectory())
       {
        String[] files = f.list();

        if (files.length == 0)
	{
	 f.delete();
	}
       }

       else
       {
        String pname = f.getParent();
        f.delete();

        if( pname!=null && pname.length() > 0 )
	{
	 f = new File( pname );

	 if(f.isDirectory())
	 {
	  String[] files = f.list();

	  if (files.length == 0)
	  {
	   f.delete();
	  }
	 }
	}
       }
      }
    }

    // IF ANY ERROR HAS OCCURRED
    catch ( Exception e )
    {
      err ( e.getMessage () );  // Output the error message
    }
  }

/*--------------------------------------------------------------------
 deleteCreated - Delete the unzipped files
--------------------------------------------------------------------*/

  public void deleteCreated ( )
  {
    try
    {
      while( created.size>0 )
      {
       String fname = t2s( created.rpop() );
       File f = new File( fname );

       if(f.isDirectory())
       {
        String[] files = f.list();

        if (files.length == 0)
	{
	 f.delete();
	}
       }

       else
       {
        String pname = f.getParent();
        f.delete();

        if( pname!=null && pname.length() > 0 )
	{
	 f = new File( pname );

	 if(f.isDirectory())
	 {
	  String[] files = f.list();

	  if (files.length == 0)
	  {
	   f.delete();
	  }
	 }
	}
       }
      }
    }

    // IF ANY ERROR HAS OCCURRED
    catch ( Exception e )
    {
      err ( e.getMessage () );  // Output the error message
    }
  }

/*--------------------------------------------------------------------
 zipDirectory - Zip up a directory
--------------------------------------------------------------------*/

  /** Zip the contents of the directory, and save it in the zipfile */
  public void zipDirectory(String dir, String zipfile, Texts t) 
  {
   try
   {
    if( dir.length()==0 )
    {
     Text z = new Text ( System.getProperty("user.dir") );
     z = z.replace ( "\\", "/" );  // Ensure only forward slashse
     dir = t2s(z) + "/";
    }

    // Create a stream to compress data and write it to the zipfile
    ZipOutputStream out = new ZipOutputStream(new FileOutputStream(zipfile));
    zipDirectoryRec( out, "", dir, t );
    out.close();
   }

    // IF ANY ERROR HAS OCCURRED
    catch ( Exception e )
    {
      err ( e.getMessage () );  // Output the error message
    }
  }

/*--------------------------------------------------------------------
 zipDirectoryRec - Zip up a directory
--------------------------------------------------------------------*/

  public void zipDirectoryRec(ZipOutputStream out, String curDir, String dir, Texts t)
  {
   try
   {
    File d = new File(dir);
    if (!d.isDirectory())
      throw new IllegalArgumentException("Compress: not a directory:  " + dir);
    String[] entries = d.list();
    byte[] buffer = new byte[4096];  // Create a buffer for copying 
    int bytes_read;

    // Loop through all entries in the directory
    for(int i = 0; i < entries.length; i++) {
      if( curDir.equals("") )
      {
       if( t.find(entries[i])>=0 )
       {
        continue;
       }
      }

      File f = new File(d, entries[i]);

      if (f.isDirectory())
      {
       zipDirectoryRec( out, curDir + entries[i] + "/", dir + entries[i] + "/", t );
       continue;
      }

      FileInputStream in = new FileInputStream(f); // Stream to read file
      ZipEntry entry = new ZipEntry(curDir + entries[i]);  // Make a ZipEntry
      out.putNextEntry(entry);                     // Store entry in zipfile
      while((bytes_read = in.read(buffer)) != -1)  // Copy bytes to zipfile
        out.write(buffer, 0, bytes_read);
      in.close();                                  // Close input stream
    }
   }

    // IF ANY ERROR HAS OCCURRED
    catch ( Exception e )
    {
      err ( e.getMessage () );  // Output the error message
    }
  }

/*--------------------------------------------------------------------
 ungzip - Un gzip a file
--------------------------------------------------------------------*/

  public void ungzip ( String iname, String oname )
  {
    try
    {
      dbg ( "Unzip \"" + iname + "\" into \"" + oname + "\"" );  // Status message

      // Open the compressed file
      GZIPInputStream in = new GZIPInputStream ( new FileInputStream ( iname ) );

      // Open the output file
      OutputStream out = new FileOutputStream ( oname );

      // Transfer bytes from the compressed file to the output file
      byte[] buf = new byte[1024];
      int len;
      while ( ( len = in.read ( buf ) ) > 0 )
      {
        out.write ( buf, 0, len );
      }

      // Close the file and stream
      in.close ();
      out.close ();
    }

    // IF ANY ERROR HAS OCCURRED
    catch ( Exception e )
    {
      err ( e.getMessage () );  // Output the error message
    }
  }

/*--------------------------------------------------------------------
 gse - Load and parse a gse file
--------------------------------------------------------------------*/

  public void gse ( String id, String dir, String pa, String pb, String sa, String sb )
  {
    String wfile = "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/by_series/GSE" + id + "/GSE" + id + "_family.soft.gz";
    String gfile = dir + "GSE" + id + ".soft.gz";
    String ufile = dir + "GSE" + id + ".soft";

    // GET THIS PART TO WORK WITH FTP FILES
    if ( !exists ( gfile ) )
    {
      url ( wfile, gfile );
    }

    if ( !exists ( ufile ) )
    {
      ungzip ( gfile, ufile );
    }

    Infile a = new Infile ( ufile );
    a.readLine ();

    Outfile ad = new Outfile ();
    Outfile pd = new Outfile ();
    Outfile sd = new Outfile ();

    int ac = 0;
    int pc = 0;
    String platformid = "";
    String sampleid = "";

    Texts sample = new Texts ();

    while ( !a.eof )
    {
      char c = a.get ( 0 );

      if ( c == '^' )
      {
        if ( a.find ( "=" ) >= 0 )
        {
          Text[] x = a.right ( -1 ).split ( " = " );

          if ( x[ 0 ].eq ( "SERIES" ) )
          {
            sd = new Outfile ( dir + "GSE" + id + ".des" );
            sd.out ( dash ( 70 ) + "\n" + x[ 1 ] + "\n" + dash ( 70 ) + "\n" );
            sd.out ( "id\t" + x[ 1 ] );
          }

          else if ( x[ 0 ].eq ( "PLATFORM" ) )
          {
            platformid = x[ 1 ].str ();

            pc++;
            sd.out ( "platform\t" + pc );
            sd.out ( "id\t" + x[ 1 ] );

            if ( pd != null )
            {
              pd.close ();
            }

            pd = new Outfile ( dir + x[ 1 ] + ".des" );
            pd.out ( "id\t" + x[ 1 ] );
          }

          else if ( x[ 0 ].eq ( "SAMPLE" ) )
          {
            sampleid = x[ 1 ].str ();
            sample.rpush ( x[ 1 ] );

            if ( pd != null )
            {
              pd.close ();
              pd = null;
            }

            ac++;
            sd.out ( "sample\t" + ac );
            sd.out ( "id\t" + x[ 1 ] );

            if ( ad != null )
            {
              ad.close ();
            }

            ad = new Outfile ( dir + x[ 1 ] + ".des" );
            ad.out ( "id\t" + x[ 1 ] );
          }
        }
      }

      if ( c == '!' )
      {
        if ( a.find ( "=" ) >= 0 )
        {
          Text[] x = a.right ( -1 ).split ( " = " );

          Text[] y = x[ 0 ].split ( "_" );

          for ( int i = 2; i < y.length; i++ )
          {
            y[ 1 ].rpush ( " " + y[ i ] );
          }

          if ( y[ 0 ].eq ( "Database" ) )
          {
          }

          else if ( y[ 0 ].eq ( "Series" ) )
          {
            sd.out ( y[ 1 ] + "\t" + x[ 1 ] );
          }

          else if ( y[ 0 ].eq ( "Platform" ) )
          {
            pd.out ( y[ 1 ] + "\t" + x[ 1 ] );
          }

          else if ( y[ 0 ].eq ( "Sample" ) )
          {
            ad.out ( y[ 1 ] + "\t" + x[ 1 ] );
          }

          else
          {
            // out ( y[ 0 ] + "\t" + y[ 1 ] + "\t" + x[ 1 ] );
          }
        }

        else if ( a.eq ( "!platform_table_begin" ) )
        {
          Texts paa = new Texts ();
          Texts pba = new Texts ();

          a.readLine ();
          Texts pl = new Texts ( a.splitTab () ).lower ();
          int pai = pl.find ( new Text ( pa ).lower () );
          int pbi = pl.find ( new Text ( pb ).lower () );

          a.readLine ();
          while ( !a.eof && a.ne ( "!platform_table_end" ) )
          {
            Text[] x = a.splitTab ();
            Text na = x[ pai ].splitSpace ()[ 0 ];
            Text nb = x[ pbi ].splitSpace ()[ 0 ];

            if ( na != null && na.size > 0 && nb != null && nb.size > 0 )
            {
              paa.rpush ( na );
              pba.rpush ( nb );
            }

            // add in code to find platform for map

            a.readLine ();
          }

          int ind[] = paa.sort ();
          pba = pba.select ( ind );

          Outfile pm = new Outfile ( dir + platformid + ".map" );

          for ( int i = 0; i < paa.size; i++ )
          {
            pm.out ( paa.get ( i ) + "\t" + pba.get ( i ) );
          }

          pm.close ();
        }

        else if ( a.eq ( "!sample_table_begin" ) )
        {
          Texts paa = new Texts ();
          Texts pba = new Texts ();

          a.readLine ();
          Texts pl = new Texts ( a.splitTab () ).lower ();
          int pai = pl.find ( new Text ( sa ).lower () );
          int pbi = pl.find ( new Text ( sb ).lower () );

          a.readLine ();
          while ( !a.eof && a.ne ( "!sample_table_end" ) )
          {
            Text[] x = a.splitTab ();
            paa.rpush ( x[ pai ] );
            pba.rpush ( x[ pbi ] );
            a.readLine ();
          }

          int ind[] = paa.sort ();
          pba = pba.select ( ind );

          Outfile pm = new Outfile ( dir + sampleid + ".gct" );

          pm.out ( "#1.2" );
          pm.out ( paa.size + "\t1" );
          pm.out ( "Name\tDescription\t" + sampleid );

          for ( int i = 0; i < paa.size; i++ )
          {
            pm.out ( paa.get ( i ) + "\t\t" + pba.get ( i ) );
          }

          pm.close ();

          pm = new Outfile ( dir + sampleid + ".cls" );
          pm.out ( "1 1 1" );
          pm.out ( "# NONE" );
          pm.out ( "NONE" );
        }
      }

      if ( c == '#' )
      {
        if ( a.find ( "=" ) >= 0 )
        {
          Text[] x = a.right ( -1 ).split ( " = " );

          if ( pd != null )
          {
            pd.out ( x[ 0 ] + "\t" + x[ 1 ] );
          }

          else
          {
            ad.out ( x[ 0 ] + "\t" + x[ 1 ] );
          }
        }
      }

      a.readLine ();
    }

    if ( ad != null )
    {
      ad.close ();
    }

    if ( pd != null )
    {
      pd.close ();
      pd = null;
    }

    a.close ();

    sample = sample.removeDuplicate ();

    a = new Infile ( dir + sample.get ( 0 ) + ".des" );
    a.readLine ();

    sd.out ( "\n" + dash ( 70 ) + "\n" + sample.get ( 0 ) + "\n" + dash ( 70 ) + "\n" );
    while ( !a.eof )
    {
      sd.out ( a );
      a.readLine ();
    }

    a.close ();

    if ( sd != null )
    {
      sd.close ();
    }
  }

/*--------------------------------------------------------------------
 gse2gct - Convert a set of description files to a gct
--------------------------------------------------------------------*/

  public void gse2gct ( String id, String dir )
  {
    gse2gct ( id, dir, null );
  }

/*--------------------------------------------------------------------
 gse2gct - Convert a set of description files to a gct
--------------------------------------------------------------------*/

  public void gse2gct ( String id, String dir, Texts sample )
  {
    if ( sample == null )
    {
      sample = new Texts ();

      Infile a = new Infile ( dir + "GSE" + id + ".des" );

      a.readLine ();
      while ( !a.eof )
      {
        Text[] x = a.splitTab ();
        if ( x.length >= 2 && x[ 1 ].find ( "GSM" ) == 0 )
        {
          sample.rpush ( x[ 1 ] );
        }
        a.readLine ();
      }
      a.close ();

      sample = sample.removeDuplicate ();
    }

    //Array sar = new Array ( dir + sample.get ( 0 ) + ".gct", dir + sample.get ( 0 ) + ".cls" );

    TextTree sdes = new TextTree ( dir + sample.get ( 0 ) + ".des" );
    Text gpl = sdes.get ( "platform id" );
    out ( gpl );

    // delete ( dir + sample.get ( 0 ) + ".gct" );
    // delete ( dir + sample.get ( 0 ) + ".cls" );
    // delete ( dir + sample.get ( 0 ) + ".des" );

    for ( int i = 1; i < sample.size; i++ )
    {
      //Array sar2 = new Array ( dir + sample.get ( i ) + ".gct", dir + sample.get ( i ) + ".cls" );
      //sar.merge ( sar2 );

      sdes = new TextTree ( dir + sample.get ( i ) + ".des" );
      gpl = sdes.get ( "platform id" );
      out ( gpl );

      // delete ( dir + sample.get ( i ) + ".gct" );
      // delete ( dir + sample.get ( i ) + ".cls" );
      // delete ( dir + sample.get ( i ) + ".des" );
    }

    // sar.write ( dir + "GSE" + id + ".gct", dir + "GSE" + id + ".cls" );
  }

/*--------------------------------------------------------------------
 copy - Copy part of a file into another
--------------------------------------------------------------------*/

  public void copy ( String iname, String oname, int a )
  {
    copy ( iname, oname, 0, a - 1 );
  }

/*--------------------------------------------------------------------
 copy - Copy part of a file into another
--------------------------------------------------------------------*/

  public void copy ( String iname, String oname, int a, int b )
  {
    Infile d = new Infile ( iname );
    Outfile e = new Outfile ( oname );

    d.in ();
    int i = 0;
    while ( !d.eof && i < a )
    {
      i++;
      d.in ();
    }

    while ( !d.eof && i <= b )
    {
      e.writeLine ( d );
      i++;
      d.in ();
    }

    e.close ();
    d.close ();
  }

/*--------------------------------------------------------------------
 count - Count the number of instances
--------------------------------------------------------------------*/

  public int[] count ( int[] z )
  {
    Tree t = new Tree ();
    Texts n = new Texts ();

    for ( int i = 0; i < z.length; i++ )
    {
      int j = z[ i ];

      t.inc ( j );

      if ( t.get ( j ) == 1 )
      {
        n.rpush ( j );
      }
    }

    n = n.removeDuplicate ();

    int[] y = new int[n.size];

    for ( int i = 0; i < n.size; i++ )
    {
      y[ i ] = t.get ( n.get ( i ) );
    }

    return y;
  }

/*--------------------------------------------------------------------
 aliasToMap - Map an alias file to a map file
--------------------------------------------------------------------*/

  public void aliasToMap ( String aname, String mname )
  {
    Infile a = new Infile ();           // Create a new Infile
    a.load ( aname );                   // Load the homolog chart
    Text[][] z = a.splitTabLine ();     // Split it up
    int zlen = z[ 0 ].length;           // Number of genes
    Outfile b = new Outfile ( mname );  // Create a new Outfile

    // LOOP THROUGH THE GENES
    for ( int i = 1; i < zlen; i++ )
    {
      Text[] x = z[ 0 ][ i ].splitSpace ();  // Split the text on the spaces
      Text[] y = z[ 1 ][ i ].splitSpace ();  // Split the text on the spaces

      // LOOP THROUGH THE IDS
      for ( int j = 0; j < x.length; j++ )
      {
        b.writeLine ( x[ j ] + "\t" + x[ 0 ] );  // Write the matching gene
      }

      // IF THERE ARE SOME GENES
      if ( z[ 1 ][ i ].size > 0 )
      {
        // LOOP THROUGH THE IDS
        for ( int j = 0; j < y.length; j++ )
        {
          b.writeLine ( y[ j ] + "\t" + x[ 0 ] );  // Write the matching gene
        }
      }
    }

    b.close ();  // Close the outfile
  }

/*--------------------------------------------------------------------
 ssdbToMap - Convert an ssdb file to a map
--------------------------------------------------------------------*/

  public void ssdbToMap ( String sname, String mname, String oa, String ob )
  {
    Infile a = new Infile ();
    a.load ( sname );
    Text[] s = a.split ( "<\n>" );

    Texts ga = new Texts ();
    Texts gb = new Texts ();

    for ( int i = 1; i < s.length; i++ )
    {
      Text genea = s[ i ].findBetween ( "\n> " + oa + ":", "(" );
      Text geneb = s[ i ].findBetween ( "\n> " + ob + ":", "(" );

      if ( genea != null && genea.size > 0 && geneb != null && geneb.size > 0 )
      {
        ga.rpush ( genea.splitTab ()[ 0 ] );
        gb.rpush ( geneb.splitTab ()[ 0 ] );
      }
    }

    int ind[] = ga.sort ();
    gb = gb.select ( ind );

    Outfile b = new Outfile ( mname );

    for ( int i = 0; i < ga.size; i++ )
    {
      b.writeLine ( ga.get ( i ) + "\t" + gb.get ( i ) );
    }

    b.close ();
  }

/*--------------------------------------------------------------------
 mapBest - Convert two maps to a best - best
--------------------------------------------------------------------*/

  public void mapBest ( String aname, String bname, String cname )
  {
    TextTree a = new TextTree ( aname );
    TextTree b = new TextTree ( bname );

    Outfile c = new Outfile ( cname );
    Text[][] m = splitTabLine ( aname );

    for ( int i = 0; i < m[ 0 ].length; i++ )
    {
      Text bg = a.get ( m[ 0 ][ i ] );

      if ( bg != null && bg.size > 0 )
      {
        Text ag = b.get ( bg );

        if ( ag != null && ag.size > 0 )
        {
          if ( m[ 0 ][ i ].eq ( ag ) )
          {
            c.writeLine ( ag + "\t" + bg );
          }
        }
      }
    }

    c.close ();
  }

/*--------------------------------------------------------------------
 mapRemoveDuplicate - Remove the duplicates from a map
--------------------------------------------------------------------*/

  public void mapRemoveDuplicate ( String aname, String bname )
  {
    Text[][] z = splitTabLine ( aname );   // Load the map file
    Texts a = new Texts ( z[ 0 ] );
    Texts b = new Texts ( z[ 1 ] );
    boolean[] u = a.unique ();
    a = a.select ( u );
    b = b.select ( u );
    joinTabLine ( bname, a, b );
  }

/*--------------------------------------------------------------------
 joinMap - Join two maps
--------------------------------------------------------------------*/

  public void joinMap ( String aname, String bname, String cname )
  {
    TextTree t = new TextTree ();
    Texts g = new Texts ();

    Infile a = new Infile ( aname );
    a.in ();

    while ( !a.eof )
    {
      Text[] x = a.splitTab ();
      t.set ( x[ 0 ], x[ 1 ] );
      g.rpush ( x[ 0 ] );
      a.in ();
    }

    a.close ();

    a = new Infile ( bname );
    a.in ();

    while ( !a.eof )
    {
      Text[] x = a.splitTab ();
      t.set ( x[ 0 ], x[ 1 ] );
      g.rpush ( x[ 0 ] );
      a.in ();
    }

    a.close ();

    g = g.removeDuplicate ();

    Outfile b = new Outfile ( cname );

    for ( int i = 0; i < g.size; i++ )
    {
      b.writeLine ( g.get ( i ) + "\t" + t.get ( g.get ( i ) ) );
    }

    b.close ();
  }

/*--------------------------------------------------------------------
 mapMap - Map a map by two mapping files
--------------------------------------------------------------------*/

  public void mapMap ( String mname, String nname, String aname, String bname )
  {
    Text[][] z = splitTabLine ( mname );   // Load the map file to update
    int zlen = z[ 0 ].length;              // Number of genes

    TextTree ta = new TextTree ();
    TextTree tb = new TextTree ();

    if ( aname != null )
    {
      Text[][] la = splitTabLine ( aname );  // Load the map for gene 1
      ta.add ( la[ 0 ], la[ 1 ] );           // Set up the map for gene 1
    }

    if ( bname != null )
    {
      Text[][] lb = splitTabLine ( bname );  // Load the map for gene 2
      tb.add ( lb[ 0 ], lb[ 1 ] );           // Set up the map gene 2
    }

    Outfile b = new Outfile ( nname );     // Create a new Outfile

    // LOOP THROUGH THE GENES
    for ( int i = 0; i < zlen; i++ )
    {
      Text g1 = z[ 0 ][ i ];  // Find gene i.1
      Text g2 = z[ 1 ][ i ];  // Find gene i.2

      if ( aname != null )
      {
        g1 = ta.get ( g1 );     // Get the gene list for 1
      }

      if ( bname != null )
      {
        g2 = tb.get ( g2 );     // Get the gene list for 2
      }

      // IF BOTH LISTS EXIST
      if ( g1 != null && g2 != null )
      {
        Text[] gn1 = g1.removeExtraWhite ().splitSpace ();  // Create the list for 1
        Text[] gn2 = g2.removeExtraWhite ().splitSpace ();  // Create the list for 2

        // LOOP THROUGH THE FIRST LIST OF GENES
        for ( int j = 0; j < gn1.length; j++ )
        {
          // LOOP THROUGH THE SECOND LIST OF GENES
          for ( int k = 0; k < gn2.length; k++ )
          {
            b.writeLine ( gn1[ j ] + "\t" + gn2[ k ] );  // Write the gene pair into the file
          }
        }
      }
    }

    b.close ();  // Close the outfile
  }

/*--------------------------------------------------------------------
 splitTabLine - Load a file, split into columns
--------------------------------------------------------------------*/

  public Text[][] splitTabLine ( String fname )
  {
    Infile m = new Infile ();        // Create a new Infile
    m.load ( fname );                // Load the homolog chart
    Text[][] z = m.splitTabLine ();  // Split it up
    return z;                        // Return the arrays
  }

/*--------------------------------------------------------------------
 joinTabLine - Save a file in tabs
--------------------------------------------------------------------*/

  public void joinTabLine ( String fname, Texts a, Texts b )
  {
    Outfile m = new Outfile ( fname );  // Create a new Outfile

    // LOOP THROUGH THE ROWS
    for ( int i = 0; i < a.size && i < b.size; i++ )
    {
      m.writeLine ( a.get ( i ) + "\t" + b.get ( i ) );  // Write the row
    }

    m.close ();  // Close the Outfile
  }

/*--------------------------------------------------------------------
 splitSpaceLine - Load a file, split into columns
--------------------------------------------------------------------*/

  public Text[][] splitSpaceLine ( String fname )
  {
    Infile m = new Infile ();        // Create a new Infile
    m.load ( fname );                // Load the homolog chart
    Text[][] z = m.splitSpaceLine ();  // Split it up
    return z;                        // Return the arrays
  }

/*--------------------------------------------------------------------
 splitLineSpace -
--------------------------------------------------------------------*/

  public Text[][] splitLineSpace ( String fname )
  {
    Infile m = new Infile ();        // Create a new Infile
    m.load ( fname );                // Load the homolog chart
    Text[][] z = m.splitLineSpace ();  // Split it up
    return z;                        // Return the arrays
  }

/*--------------------------------------------------------------------
 splitLineTab -
--------------------------------------------------------------------*/

  public Text[][] splitLineTab ( String fname )
  {
    Infile m = new Infile ();        // Create a new Infile
    m.load ( fname );                // Load the homolog chart
    Text[][] z = m.splitLineTab ();  // Split it up
    return z;                        // Return the arrays
  }

/*--------------------------------------------------------------------
 mapDescription - Create a map descripion file
--------------------------------------------------------------------*/

  public void mapDescription ( String iname, String oname, String m1name, String m2name, String d1name, String d2name )
  {
    TextTree d1 = new TextTree ( d1name );
    TextTree d2 = new TextTree ( d2name );
    TextTree map1 = new TextTree ( m1name );

    TextTree map2 = null;

    if ( m2name.length () > 0 )
    {
      map2 = new TextTree ( m2name );
    }

    Infile d = new Infile ( iname );
    Outfile e = new Outfile ( oname );
    d.readLine ();

    while ( !d.eof )
    {
      Text z = new Text ( "" );
      Text mapid = map1.get ( d );

      if ( mapid != null && mapid.size > 0 && map2 != null )
      {
        mapid = map2.get ( mapid );
      }

      if ( mapid != null && mapid.size > 0 )
      {
        z.rpush ( d );
        z.rpush ( "\t" );

        Text des = d1.get ( d );
        if ( des != null )
        {
          z.rpush ( des );
        }

        z.rpush ( "\n" );

        des = d2.get ( mapid );
        if ( des != null )
        {
          z.rpush ( mapid );
          z.rpush ( "\t" );
          z.rpush ( des );
        }

        z.rpush ( "\n" );
        e.out ( z );
      }

      d.readLine ();
    }

    e.close ();
    d.close ();
  }

/*--------------------------------------------------------------------
 convertPred - Convert a pred file to other output types
--------------------------------------------------------------------*/

  public void convertPredTab ( String pname, String tname, String fname )
  {
    Infile d = new Infile ( pname );
    d.readLine ();
    d.readLine ();

    Texts p = new Texts ();
    while ( !d.eof )
    {
      if ( d.eq ( "" ) )
      {
        break;
      }

      p.rpush ( new Text ( d ) );
      d.readLine ();
    }
    d.close ();

    TextTree t = new TextTree ();
    d = new Infile ( tname );
    d.readLine ();

    int found = 0;

    while ( !d.eof )
    {
      Text z = new Text ( d );

      Text[] x = z.split ( ". " );

      Text y = x[ 1 ];
      for ( int j = 2; j < x.length; j++ )
      {
        y.rpush ( x[ j ] );
      }

      y = y.alphaNum ().replace ( " ", "" );

      Text y2 = new Text ( "" );

      d.readLine ();
      while ( d.ne ( "" ) )
      {
        if ( y2.size > 0 )
        {
          y2.rpush ( " :: " );
        }

        y2.rpush ( d.trimWhite () );
        d.readLine ();
      }

      t.set ( y, y2 );
      d.readLine ();
    }

    d.close ();

    Texts samp = new Texts ();
    Texts conf = new Texts ();
    Texts phen = new Texts ();

    Outfile e = new Outfile ( fname + ".pred.txt" );
    for ( int i = 0; i < p.size; i++ )
    {
      Text z = new Text ();

      Text[] x = p.get ( i ).splitTab ();
      Text y = x[ 0 ].alphaNum ().replace ( " ", "" );
      Text attr = t.get ( y );

      z.rpush ( x[ 0 ] );
      samp.rpush ( x[ 0 ] );

      if ( attr != null )
      {
        z.rpush ( "\t" );
        z.rpush ( attr );
      }

      else
      {
        z.rpush ( "\t" );
      }

      z.rpush ( "\t" );
      z.rpush ( x[ 2 ] );
      phen.rpush ( x[ 2 ] );

      z.rpush ( "\t" );
      z.rpush ( x[ 4 ] );

      z.rpush ( "\t" );
      z.rpush ( x[ 5 ] );
      conf.rpush ( x[ 5 ] );

      z.rpush ( "\t" );
      z.rpush ( x[ 9 ] );

      z.rpush ( "\t" );
      z.rpush ( x[ 10 ] );

      z.rpush ( "\t" );
      z.rpush ( x[ 11 ] );

      e.writeLine ( z );
    }

    e.close ();

    Texts phenName = phen.removeDuplicate ();
    phenName.sort ();

    int numPhen = phenName.size;
    Tree[] an = new Tree[numPhen];
    Texts alist = new Texts ();

    int numToCheck = 100;
    double confCutoff = 0.0;

    e = new Outfile ( fname + ".combo.txt" );

    for ( int i = 0; i < numPhen; i++ )
    {
      int ct = 0;
      int ct2 = 0;

      int[] pi = phen.select ( phenName.get ( i ) );
      Texts csamp = samp.select ( pi );
      Texts cconf = conf.select ( pi );
      Texts cphen = phen.select ( pi );

      int[] ind = cconf.sortd ();
      cphen = cphen.select ( ind );
      csamp = csamp.select ( ind );

      an[ i ] = new Tree ();
      for ( int j = 0; j < cphen.size; j++ )
      {
        if ( j >= numToCheck && numToCheck > 0 )
        {
          break;
        }

        if ( cconf.get ( j ).dbl () < confCutoff )
        {
          break;
        }

        Text y = csamp.get ( j ).alphaNum ().replace ( " ", "" );
        Text attr = t.get ( y );

        if ( attr != null )
        {
          Text x[] = attr.split ( " :: " );

          for ( int k = 0; k < x.length; k++ )
          {
            if ( x[ k ].size > 0 )
            {
              an[ i ].inc ( x[ k ] );
              alist.rpush ( x[ k ] );
            }
          }

          ct++;
        }

        ct2++;
        //e.out ( csamp.get ( j ) + "\t" + attr + "\t" + cconf.get ( j ) + "\t" + cphen.get ( j ) );
      }

      out ( phenName.get ( i ) + " " + ct + " " + ct2 + " " + csamp.size );
    }

    alist = alist.removeDuplicate ();
    alist.sort ();

    Text z = new Text ( "Sample" );

    for ( int i = 0; i < numPhen; i++ )
    {
      z.rpush ( "\t" );
      z.rpush ( phenName.get ( i ) );
    }

    for ( int i = 0; i < numPhen; i++ )
    {
      z.rpush ( "\t" );
      z.rpush ( phenName.get ( i ) );
    }

    z.rpush ( "\tSum" );
    e.out ( z );

    for ( int i = 0; i < alist.size; i++ )
    {
      int[] score = new int[numPhen];
      int sum = 0;

      for ( int j = 0; j < numPhen; j++ )
      {
        score[ j ] = an[ j ].get ( alist.get ( i ) );
        sum += score[ j ];
      }

      z = alist.get ( i );

      for ( int j = 0; j < numPhen; j++ )
      {
        z.rpush ( "\t" );
        z.rpush ( ( double ) ( score[ j ] ) / ( double ) ( sum ) );
      }

      for ( int j = 0; j < numPhen; j++ )
      {
        z.rpush ( "\t" );
        z.rpush ( score[ j ] );
      }

      z.rpush ( "\t" );
      z.rpush ( sum );

      e.out ( z );
    }

    e.close ();
  }

/*--------------------------------------------------------------------
 copy - Copy a file into another
--------------------------------------------------------------------*/

  public void copy ( String iname, String oname )
  {
    // TRY TO COPY THE FILE
    try
    {
      File f = new File ( oname );  // Create a file object

      // IF THE FILE DOES NOT ALREADY EXIST
      if ( !f.exists () )
      {
        Outfile a = new Outfile ( oname );                // Create the outfile
        f = new File ( iname );                           // Create a file object
        FileInputStream fin = new FileInputStream ( f );  // Create a file input stream
        byte[] buffer = new byte[bufSize];              // Create a byte buffer for the input
        int numRead;                                      // Track the number of bytes read

        // LOOP UNTIL THE END OF THE FILE REACHED
        while ( ( numRead = fin.read ( buffer ) ) != -1 )
        {
          a.fout.write ( buffer, 0, numRead );  // Output to the file
        }

        fin.close ();                                            // Close the input stream
        a.close ();                                              // Close the outfile
        dbg ( "COPY \"" + iname + " \" TO \"" + oname + "\"" );  // Status message
      }
    }

    // IF ANY ERROR HAS OCCURRED
    catch ( Exception e )
    {
      err ( e.getMessage () );  // Output the error message
    }
  }

/*--------------------------------------------------------------------
 CHECK FUNCTIONS
 ---------------------------------------------------------------------
 Ensure variables satisfy certain conditions
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 lchk - Less than check
--------------------------------------------------------------------*/

  public int lchk ( int z, int y )
  {
    // IF NOT LESS THAN
    if ( z >= y )
    {
      z = y - 1;  // Update the value
    }

    return z;  // Return the updated value
  }

/*--------------------------------------------------------------------
 lechk - Less than or equal to check
--------------------------------------------------------------------*/

  public int lechk ( int z, int y )
  {
    // IF GREATER THAN
    if ( z > y )
    {
      z = y;  // Update the value
    }

    return z;  // Return the updated value
  }

/*--------------------------------------------------------------------
 gchk - Greater than check
--------------------------------------------------------------------*/

  public int gchk ( int z, int y )
  {
    // IF NOT GREATER THAN
    if ( z <= y )
    {
      z = y + 1;  // Update the value
    }

    return z;  // Return the updated value
  }

/*--------------------------------------------------------------------
 nchk - Null check
--------------------------------------------------------------------*/

  public Text nchk ( Text a, Text b )
  {
    if ( a == null )
    {
      return b;
    }

    else
    {
      return a;
    }
  }

/*--------------------------------------------------------------------
 gechk - Greater than or equal to check
--------------------------------------------------------------------*/

  public int gechk ( int z, int y )
  {
    // IF LESS THAN
    if ( z < y )
    {
      z = y;  // Update the value
    }

    return z;  // Return the updated value
  }

/*--------------------------------------------------------------------
 eq - is a equal to this character array
--------------------------------------------------------------------*/

  public boolean eq ( char[] a, char[] b )
  {
    return ( cmp ( a, b ) == 0 );
  }

/*--------------------------------------------------------------------
 ne - is a not equal to this character array
--------------------------------------------------------------------*/

  public boolean ne ( char[] a, char[] b )
  {
    return ( cmp ( a, b ) != 0 );
  }

/*--------------------------------------------------------------------
 lt - is a less than this character array
--------------------------------------------------------------------*/

  public boolean lt ( char[] a, char[] b )
  {
    return ( cmp ( a, b ) < 0 );
  }

/*--------------------------------------------------------------------
 le - is a less than or equal to this character array
--------------------------------------------------------------------*/

  public boolean le ( char[] a, char[] b )
  {
    return ( cmp ( a, b ) <= 0 );
  }

/*--------------------------------------------------------------------
 gt - is a greater than this character array
--------------------------------------------------------------------*/

  public boolean gt ( char[] a, char[] b )
  {
    return ( cmp ( a, b ) > 0 );
  }

/*--------------------------------------------------------------------
 ge - is a greater than or equal to this character array
--------------------------------------------------------------------*/

  public boolean ge ( char[] a, char[] b )
  {
    return ( cmp ( a, b ) >= 0 );
  }

/*--------------------------------------------------------------------
 le - count the doubles less than or equal to b
--------------------------------------------------------------------*/

  public double le ( double[] a, double b )
  {
    double z = 0;

    for ( int i = 0; i < a.length; i++ )
    {
      if ( a[ i ] <= b )
      {
        z++;
      }
    }

    return z;
  }

/*--------------------------------------------------------------------
 lt - count the doubles less than to b
--------------------------------------------------------------------*/

  public double lt ( double[] a, double b )
  {
    double z = 0;

    for ( int i = 0; i < a.length; i++ )
    {
      if ( a[ i ] < b )
      {
        z++;
      }
    }

    return z;
  }

/*--------------------------------------------------------------------
 ge - count the doubles greater than or equal to b
--------------------------------------------------------------------*/

  public double ge ( double[] a, double b )
  {
    double z = 0;

    for ( int i = 0; i < a.length; i++ )
    {
      if ( a[ i ] >= b )
      {
        z++;
      }
    }

    return z;
  }

/*--------------------------------------------------------------------
 gt - count the doubles greater than or equal to b
--------------------------------------------------------------------*/

  public double gt ( double[] a, double b )
  {
    double z = 0;

    for ( int i = 0; i < a.length; i++ )
    {
      if ( a[ i ] > b )
      {
        z++;
      }
    }

    return z;
  }

/*--------------------------------------------------------------------
 gtInd - return the indices greater than b
--------------------------------------------------------------------*/

  public int[] gtInd ( double[] a, double b )
  {
    Texts z = new Texts ();

    for ( int i = 0; i < a.length; i++ )
    {
      if ( a[ i ] > b )
      {
        z.rpush ( i );
      }
    }

    return z.num ();
  }

/*--------------------------------------------------------------------
 geInd - return the indices greater than or equal to b
--------------------------------------------------------------------*/

  public int[] geInd ( double[] a, double b )
  {
    Texts z = new Texts ();

    for ( int i = 0; i < a.length; i++ )
    {
      if ( a[ i ] >= b )
      {
        z.rpush ( i );
      }
    }

    return z.num ();
  }

/*--------------------------------------------------------------------
 ltInd - return the indices less than b
--------------------------------------------------------------------*/

  public int[] ltInd ( double[] a, double b )
  {
    Texts z = new Texts ();

    for ( int i = 0; i < a.length; i++ )
    {
      if ( a[ i ] < b )
      {
        z.rpush ( i );
      }
    }

    return z.num ();
  }

/*--------------------------------------------------------------------
 leInd - return the indices less than or equal to b
--------------------------------------------------------------------*/

  public int[] leInd ( double[] a, double b )
  {
    Texts z = new Texts ();

    for ( int i = 0; i < a.length; i++ )
    {
      if ( a[ i ] <= b )
      {
        z.rpush ( i );
      }
    }

    return z.num ();
  }

/*--------------------------------------------------------------------
 cmp - Compare two character arrays
--------------------------------------------------------------------*/

  public int cmp ( char[] a, char[] b )
  {
    int ans = 0;

    for ( int i = 0; i < a.length && i < b.length; i++ )
    {
      if ( a[ i ] != b[ i ] )
      {
        if ( a[ i ] < b[ i ] )
        {
          ans = -1;
        }

        else
        {
          ans = 1;
        }

        break;
      }
    }

    if ( ans == 0 )
    {
      if ( a.length < b.length )
      {
        ans = -1;
      }

      else if ( b.length < a.length )
      {
        ans = 1;
      }
    }

    return ans;
  }

/*--------------------------------------------------------------------
 eq - is a equal to this character array
--------------------------------------------------------------------*/

  public boolean eq ( char[] a, int c, char[] b, int d )
  {
    return ( cmp ( a, c, b, d ) == 0 );
  }

/*--------------------------------------------------------------------
 ne - is a not equal to this character array
--------------------------------------------------------------------*/

  public boolean ne ( char[] a, int c, char[] b, int d )
  {
    return ( cmp ( a, c, b, d ) != 0 );
  }

/*--------------------------------------------------------------------
 lt - is a less than this character array
--------------------------------------------------------------------*/

  public boolean lt ( char[] a, int c, char[] b, int d )
  {
    return ( cmp ( a, c, b, d ) < 0 );
  }

/*--------------------------------------------------------------------
 le - is a less than or equal to this character array
--------------------------------------------------------------------*/

  public boolean le ( char[] a, int c, char[] b, int d )
  {
    return ( cmp ( a, c, b, d ) <= 0 );
  }

/*--------------------------------------------------------------------
 gt - is a greater than this character array
--------------------------------------------------------------------*/

  public boolean gt ( char[] a, int c, char[] b, int d )
  {
    return ( cmp ( a, c, b, d ) > 0 );
  }

/*--------------------------------------------------------------------
 ge - is a greater than or equal to this character array
--------------------------------------------------------------------*/

  public boolean ge ( char[] a, int c, char[] b, int d )
  {
    return ( cmp ( a, c, b, d ) >= 0 );
  }

/*--------------------------------------------------------------------
 cmp - Compare two text objects
--------------------------------------------------------------------*/

  public int cmp ( char[] a, int c, char[] b, int d )
  {
    int ans = 0;

    for ( int i = 0; i < a.length - c && i < b.length - d; i++ )
    {
      if ( a[ i + c ] != b[ i + d ] )
      {
        if ( a[ i + c ] < b[ i + d ] )
        {
          ans = -1;
        }

        else
        {
          ans = 1;
        }

        break;
      }
    }

    if ( ans == 0 )
    {
      if ( a.length - c < b.length - d )
      {
        ans = -1;
      }

      else if ( b.length - d < a.length - c )
      {
        ans = 1;
      }
    }

    return ans;
  }

/*--------------------------------------------------------------------
 IS FUNCTIONS
 ---------------------------------------------------------------------
 Handle queries about variables
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 isEmpty - Is the text empty?
--------------------------------------------------------------------*/

  public boolean isEmpty ( Text t )
  {
    // IF IT IS NULL OR ZERO LENGTH
    return t == null || t.size == 0;
  }

/*--------------------------------------------------------------------
 notEmpty - Is the text not empty?
--------------------------------------------------------------------*/

  public boolean notEmpty ( Text t )
  {
    // IF IT IS NULL OR ZERO LENGTH
    return !isEmpty ( t );
  }

/*--------------------------------------------------------------------
 isWhite - Is the character whitespace?
--------------------------------------------------------------------*/

  public boolean isWhite ( char c )
  {
    // IF IT IS A SPACE, TAB, OR NEW LINE
    return c == ' ' || c == '\t' || c == '\r' || c == '\n';
  }

/*--------------------------------------------------------------------
 isLower - Is the character lowercase?
--------------------------------------------------------------------*/

  public boolean isLower ( char c )
  {
    // IF A LOWERCASE CHARACTER
    return c >= 'a' && c <= 'z';
  }

/*--------------------------------------------------------------------
 isUpper - Is the character uppercase?
--------------------------------------------------------------------*/

  public boolean isUpper ( char c )
  {
    // IF AN UPPERCASE CHARACTER
    return c >= 'A' && c <= 'Z';
  }

/*--------------------------------------------------------------------
 isAlpha- Is the character a letter?
--------------------------------------------------------------------*/

  public boolean isAlpha ( char c )
  {
    // IF A LETTER
    return isLower ( c ) || isUpper ( c );
  }

/*--------------------------------------------------------------------
 isNum - Is the character a number?
--------------------------------------------------------------------*/

  public boolean isNum ( char c )
  {
    // IF A NUMBER
    return c >= '0' && c <= '9';
  }

/*--------------------------------------------------------------------
 isAlphaNum - Is the character alphanumeric?
--------------------------------------------------------------------*/

  public boolean isAlphaNum ( char c )
  {
    // IF A LETTER OR A NUMBER
    return isAlpha ( c ) || isNum ( c );
  }

/*--------------------------------------------------------------------
 isAlphaNumSign - Is the character alphanumeric or a sign?
--------------------------------------------------------------------*/

  public boolean isAlphaNumSign ( char c )
  {
    // IF A LETTER OR A NUMBER
    return isAlpha ( c ) || isNum ( c ) || c == '+' || c == '-';
  }

/*--------------------------------------------------------------------
 TEXT FUNCTIONS
 ---------------------------------------------------------------------
 Handle text manipulation
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 joinSpace - Join a text array a with spaces
--------------------------------------------------------------------*/

  public Text joinSpace ( Text a[] )
  {
    return join ( a, new Text ( " " ) );
  }

/*--------------------------------------------------------------------
 joinTab - Join a text array a with tabs
--------------------------------------------------------------------*/

  public Text joinTab ( Text a[] )
  {
    return join ( a, new Text ( "\t" ) );
  }

/*--------------------------------------------------------------------
 joinTab - Join an int array a with tabs
--------------------------------------------------------------------*/

  public Text joinTab ( int a[] )
  {
    return join ( a, new Text ( "\t" ) );
  }

/*--------------------------------------------------------------------
 join - Join a text array a with string b
--------------------------------------------------------------------*/

  public Text join ( Text a[], String b )
  {
    return join ( a, new Text ( b ) );
  }

/*--------------------------------------------------------------------
 join - Join a text array a with text b
--------------------------------------------------------------------*/

  public Text join ( Text a[], Text b )
  {
    int asize = a.length;  // Size of a
    Text c = new Text ();  // Text to join into

    // IF A CONTAINS ANY TEXT
    if ( asize > 0 )
    {
      c.rpush ( a[ 0 ] );  // Push the first element

      // LOOP THROUGH THE ELEMENTS
      for ( int i = 1; i < asize; i = i + 1 )
      {
        c.rpush ( b );       // Push b
        c.rpush ( a[ i ] );  // Push the ith element
      }
    }

    return c;  // Return the text
  }

/*--------------------------------------------------------------------
 join - Join an int array a with text b
--------------------------------------------------------------------*/

  public Text join ( int a[], Text b )
  {
    int asize = a.length;  // Size of a
    Text c = new Text ();  // Text to join into

    // IF A CONTAINS ANY TEXT
    if ( asize > 0 )
    {
      c.rpush ( a[ 0 ] );  // Push the first element

      // LOOP THROUGH THE ELEMENTS
      for ( int i = 1; i < asize; i = i + 1 )
      {
        c.rpush ( b );       // Push b
        c.rpush ( a[ i ] );  // Push the ith element
      }
    }

    return c;  // Return the text
  }

/*--------------------------------------------------------------------
 c - Concatenate
--------------------------------------------------------------------*/

  public Text c ( Text a, Text b )
  {
    Text z = new Text ( a );  // Set to a
    z.rpush ( b );            // Push b
    return z;                 // Return the text
  }

/*--------------------------------------------------------------------
 c - Concatenate
--------------------------------------------------------------------*/

  public Text c ( Text a, Text b, Text c )
  {
    Text z = new Text ( a );  // Set to a
    z.rpush ( b );            // Push b
    z.rpush ( c );            // Push c
    return z;                 // Return the text
  }

/*--------------------------------------------------------------------
 c - Concatenate
--------------------------------------------------------------------*/

  public Text c ( Text a, Text b, Text c, Text d )
  {
    Text z = new Text ( a );  // Set to a
    z.rpush ( b );            // Push b
    z.rpush ( c );            // Push c
    z.rpush ( d );            // Push d
    return z;                 // Return the text
  }

/*--------------------------------------------------------------------
 c - Concatenate
--------------------------------------------------------------------*/

  public Text c ( Text a, Text b, Text c, Text d, Text e )
  {
    Text z = new Text ( a );  // Set to a
    z.rpush ( b );            // Push b
    z.rpush ( c );            // Push c
    z.rpush ( d );            // Push d
    z.rpush ( e );            // Push e
    return z;                 // Return the text
  }

/*--------------------------------------------------------------------
 n - Concatenate
--------------------------------------------------------------------*/

  public Text n ( Text a, Text b )
  {
    Text z = new Text ( a );  // Set to a
    z.rpush ( "\n" );         // Push a new line
    z.rpush ( b );            // Push b
    return z;                 // Return the text
  }

/*--------------------------------------------------------------------
 n - Concatenate
--------------------------------------------------------------------*/

  public Text n ( Text a, Text b, Text c )
  {
    Text z = new Text ( a );  // Set to a
    z.rpush ( "\n" );         // Push a new line
    z.rpush ( b );            // Push b
    z.rpush ( "\n" );         // Push a new line
    z.rpush ( c );            // Push c
    return z;                 // Return the text
  }

/*--------------------------------------------------------------------
 n - Concatenate
--------------------------------------------------------------------*/

  public Text n ( Text a, Text b, Text c, Text d )
  {
    Text z = new Text ( a );  // Set to a
    z.rpush ( "\n" );         // Push a new line
    z.rpush ( b );            // Push b
    z.rpush ( "\n" );         // Push a new line
    z.rpush ( c );            // Push c
    z.rpush ( "\n" );         // Push a new line
    z.rpush ( d );            // Push d
    return z;                 // Return the text
  }

/*--------------------------------------------------------------------
 n - Concatenate
--------------------------------------------------------------------*/

  public Text n ( Text a, Text b, Text c, Text d, Text e )
  {
    Text z = new Text ( a );  // Set to a
    z.rpush ( "\n" );         // Push a new line
    z.rpush ( b );            // Push b
    z.rpush ( "\n" );         // Push a new line
    z.rpush ( c );            // Push c
    z.rpush ( "\n" );         // Push a new line
    z.rpush ( d );            // Push d
    z.rpush ( "\n" );         // Push a new line
    z.rpush ( e );            // Push e
    return z;                 // Return the text
  }

/*--------------------------------------------------------------------
 STRING FUNCTIONS
 ---------------------------------------------------------------------
 Handle string manipulation
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 chrat - Return the specified character in a string
--------------------------------------------------------------------*/

  public char chrat ( String z, int i )
  {
    return z.charAt ( i );
  }

/*--------------------------------------------------------------------
 strlen - Find the length of a string
--------------------------------------------------------------------*/

  public int strlen ( String z )
  {
    return z.length ();
  }

/*--------------------------------------------------------------------
 copy - Copy a matrix
--------------------------------------------------------------------*/

  public double[][] copy ( double[][] z )
  {
    double[][] y = new double[z.length][];

    for ( int i = 0; i < z.length; i++ )
    {
      y[ i ] = new double[z[ i ].length];

      for ( int j = 0; j < z[ i ].length; j++ )
      {
        y[ i ][ j ] = z[ i ][ j ];
      }
    }

    return y;
  }

/*--------------------------------------------------------------------
 copy - Copy a matrix
--------------------------------------------------------------------*/

  public double[] copy ( double[] z )
  {
    double[] y = new double[z.length];

    for ( int i = 0; i < z.length; i++ )
    {
      y[ i ] = z[ i ];
    }

    return y;
  }

/*--------------------------------------------------------------------
 c - Concatenate
--------------------------------------------------------------------*/

  public int[] c ( int[] a, int[] b )
  {
    int[] z = new int[a.length + b.length];

    for ( int i = 0; i < a.length; i++ )
    {
      z[ i ] = a[ i ];
    }

    for ( int i = 0; i < b.length; i++ )
    {
      z[ a.length + i ] = b[ i ];
    }

    return z;
  }

/*--------------------------------------------------------------------
 c - Concatenate
--------------------------------------------------------------------*/

  public String c ( String a, String b )
  {
    return a + b;
  }

/*--------------------------------------------------------------------
 c - Concatenate
--------------------------------------------------------------------*/

  public String c ( String a, String b, String c )
  {
    return a + b + c;
  }

/*--------------------------------------------------------------------
 c - Concatenate
--------------------------------------------------------------------*/

  public String c ( String a, String b, String c, String d )
  {
    return a + b + c + d;
  }

/*--------------------------------------------------------------------
 c - Concatenate
--------------------------------------------------------------------*/

  public String c ( String a, String b, String c, String d, String e )
  {
    return a + b + c + d + e;
  }

/*--------------------------------------------------------------------
 n - Concatenate
--------------------------------------------------------------------*/

  public String n ( String a, String b )
  {
    return a + "\n" + b;
  }

/*--------------------------------------------------------------------
 n - Concatenate
--------------------------------------------------------------------*/

  public String n ( String a, String b, String c )
  {
    return a + "\n" + b + "\n" + c;
  }

/*--------------------------------------------------------------------
 n - Concatenate
--------------------------------------------------------------------*/

  public String n ( String a, String b, String c, String d )
  {
    return a + "\n" + b + "\n" + c + "\n" + d;
  }

/*--------------------------------------------------------------------
 n - Concatenate
--------------------------------------------------------------------*/

  public String n ( String a, String b, String c, String d, String e )
  {
    return a + "\n" + b + "\n" + c + "\n" + d + "\n" + e;
  }

/*--------------------------------------------------------------------
 CONVERSION FUNCTIONS
 ---------------------------------------------------------------------
 Handle string manipulation
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 h2s - Convert a character to a string
--------------------------------------------------------------------*/

  public String h2s ( char z )
  {
    return new String ( z + "" );
  }

/*--------------------------------------------------------------------
 b2c - Convert a byte array into a character array
--------------------------------------------------------------------*/

  public char[] b2c ( byte[] z, int length )
  {
    char[] y = new char[length];  // Create a new character array

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < length; i++ )
    {
      y[ i ] = ( char ) ( z[ i ] );  // Convert the character
    }

    return y;  // Return the array
  }

/*--------------------------------------------------------------------
 c2s - Convert a character array to a string
--------------------------------------------------------------------*/

  public String c2s ( char[] z )
  {
    return new String ( z );
  }

/*--------------------------------------------------------------------
 s2c - Convert a string to a character array
--------------------------------------------------------------------*/

  public char[] s2c ( String z )
  {
    char[] y = new char[z.length ()];   // Create the array
    z.getChars ( 0, z.length (), y, 0 );  // Fill the array
    return y;                             // Return the array
  }

/*--------------------------------------------------------------------
 d2s - Convert a double to a string
--------------------------------------------------------------------*/

  public String d2s ( double z )
  {
    return new String ( z + "" );
  }

/*--------------------------------------------------------------------
 d2i - Convert a double to an integer
--------------------------------------------------------------------*/

  public int d2i ( double z )
  {
    return ( int ) ( z );
  }

/*--------------------------------------------------------------------
 d2i - Convert a double array to an integer array
--------------------------------------------------------------------*/

  public int[] d2i ( double[] z )
  {
    int[] y = new int[z.length];

    for ( int i = 0; i < y.length; i++ )
    {
      y[ i ] = ( int ) ( z[ i ] );
    }

    return y;
  }

/*--------------------------------------------------------------------
 i2s - Convert an integer to a string
--------------------------------------------------------------------*/

  public String i2s ( int z )
  {
    return new String ( z + "" );
  }

/*--------------------------------------------------------------------
 i2t - Convert an integer to a Text
--------------------------------------------------------------------*/

  public Text i2t ( int z )
  {
    return new Text ( z );
  }

/*--------------------------------------------------------------------
 i2tr - Convert an integer array to a Tree
--------------------------------------------------------------------*/

  public Tree i2tr ( int z[] )
  {
    Tree t = new Tree ();

    for ( int i = 0; i < z.length; i++ )
    {
      t.set ( z[ i ] + "", 1 );
    }

    return t;
  }

/*--------------------------------------------------------------------
 s2i - Convert a string to an integer
--------------------------------------------------------------------*/

  public int s2i ( String z )
  {
    Text t = new Text ( z );  // Create a Text version
    return t.num ();          // Return the integer
  }

/*--------------------------------------------------------------------
 h2c - Convert a character to a character array
--------------------------------------------------------------------*/

  public char[] h2c ( char z )
  {
    char[] y = new char[1];  // Create a character array
    y[ 0 ] = z;                // Store the character
    return y;                  // Return the character array
  }

/*--------------------------------------------------------------------
 h2i - Convert a character to an integer
--------------------------------------------------------------------*/

  public int h2i ( char z )
  {
    return ( int ) ( z - '0' );
  }

/*--------------------------------------------------------------------
 h2l - Convert a character to lowercase
--------------------------------------------------------------------*/

  public char h2l ( char z )
  {
    // IF THE CHARACTER IS UPPERCASE
    if ( isUpper ( z ) )
    {
      return ( char ) ( ( char ) ( z - 'A' ) + 'a' );
    }

    // IF THE CHARACTER IS NOT UPPERCASE
    else
    {
      return z;
    }
  }

/*--------------------------------------------------------------------
 h2u - Convert a character to uppercase
--------------------------------------------------------------------*/

  public char h2u ( char z )
  {
    // IF THE CHARACTER IS LOWERCASE
    if ( isLower ( z ) )
    {
      return ( char ) ( ( char ) ( z - 'a' ) + 'A' );
    }

    // IF THE CHARACTER IS NOT LOWERCASE
    else
    {
      return z;
    }
  }

/*--------------------------------------------------------------------
 s2t - Convert a string to a Text
--------------------------------------------------------------------*/

  public Text s2t ( String z )
  {
    return new Text ( z );  // Return the new text
  }

/*--------------------------------------------------------------------
 s2t - Convert an array of strings to an array of Texts
--------------------------------------------------------------------*/

  public Text[] s2t ( String[] z )
  {
    Text[] y = new Text[z.length];  // Create the new array

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < z.length; i++ )
    {
      y[ i ] = s2t ( z[ i ] );  // Convert the string to Text
    }

    return y;  // Return the new Text array
  }

/*--------------------------------------------------------------------
 t2s - Convert a Text to a string
--------------------------------------------------------------------*/

  public String t2s ( Text z )
  {
    if ( z == null )
    {
      return "";
    }

    else
    {
      return z.str ();  // Return the new string
    }
  }

/*--------------------------------------------------------------------
 s2t - Convert an array of Texts to an array of strings
--------------------------------------------------------------------*/

  public String[] t2s ( Text[] z )
  {
    String[] y = new String[z.length];  // Create the new array

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < z.length; i++ )
    {
      y[ i ] = t2s ( z[ i ] );  // Convert the Text to string
    }

    return y;  // Return the new string array
  }

/*--------------------------------------------------------------------
 b2i - Convert a boolean array to an integer index array
--------------------------------------------------------------------*/

  public int[] b2i ( boolean[] b )
  {
    int asize = 0;         // Create a size variable for a
    int bsize = b.length;  // Create a size bariable for b

    // FIND THE NUMBER OF INDICES TO SELECT
    for ( int i = 0; i < bsize; i++ )
    {
      // IF THE INDEX IS SELECTED
      if ( b[ i ] )
      {
        asize++;  // Increment the size
      }
    }

    int[] a = new int[asize];  // New select array
    int j = 0;                   // Counter variable

    // LOOP THROUGH THE ITEMS
    for ( int i = 0; i < bsize; i++ )
    {
      // IF THE INDEX IS SELECTED
      if ( b[ i ] )
      {
        a[ j ] = i;  // Set it in the select array
        j++;         // Increment the counter
      }
    }

    return a;  // Return the index array
  }

/*--------------------------------------------------------------------
 rgb - create an rgb color
--------------------------------------------------------------------*/

  public double[] rgb ( double r, double g, double b )
  {
    return new double[] { r, g, b };
  }

/*--------------------------------------------------------------------
 hsv - convert hsv to rgb
--------------------------------------------------------------------*/

  public double[] hsv ( double h, double s, double v )
  {
    return hsv2rgb ( h, s, v );
  }

/*--------------------------------------------------------------------
 rgb - rgb
--------------------------------------------------------------------*/

  public double[] rgb ( String h )
  {
   double[] c = new double[3];

   if(h.length()>6)
   {
    h = h.substring( h.length()-6, h.length() );
   }

   if(h.length()==6)
   {
    c[0] = (double) (hex2num( h.substring(0,2) ) ) / 255.0;
    c[1] = (double) (hex2num( h.substring(2,4) ) ) / 255.0;
    c[2] = (double) (hex2num( h.substring(4,6) ) ) / 255.0;
   }

   return c;
  }

/*--------------------------------------------------------------------
 hex2int - convert hex to integer
--------------------------------------------------------------------*/

  public int hex2num ( String h )
  {
   int n = 0;
   String a = "0123456789ABCDEF";

   for(int i=0; i<h.length(); i++)
   {
    n *= 16;
    n += a.indexOf(h.charAt(i));
   }

   return n;
  }

/*--------------------------------------------------------------------
 hsv2rgb - convert hsv to rgb
--------------------------------------------------------------------*/

  public double[] hsv2rgb ( double h, double s, double v )
  {
    // If grey
    if ( s == 0 )
    {
      return new double[] { v, v, v };
    }

    if ( h < 1.0 / 3.0 )
    {
      h *= ( 1.0 / 2.0 );
    }

    else if ( h < 1.0 / 2.0 )
    {
      h -= ( 1.0 / 6.0 );
    }

    else if ( h < 2.0 / 3.0 )
    {
      h = h * 2.0 - ( 2.0 / 3.0 );
    }

    h *= 6.0;
    if ( h >= 6.0 )
    {
      h = 0;
    }

    int i = ( int ) ( Math.floor ( h ) );

    double f = h - i;
    double p = v * ( 1 - s );
    double q = v * ( 1 - s * f );
    double t = v * ( 1 - s * ( 1 - f ) );

    if ( i == 0 )
    {
      return new double[] { v, t, p };
    }

    else if ( i == 1 )
    {
      return new double[] { q, v, p };
    }

    else if ( i == 2 )
    {
      return new double[] { p, v, t };
    }

    else if ( i == 3 )
    {
      return new double[] { p, q, v };
    }

    else if ( i == 4 )
    {
      return new double[] { t, p, v };
    }

    else
    {
      return new double[] { v, p, q };
    }
  }

/*--------------------------------------------------------------------
 b2b - Convert a boolean array to an integer boolean array
--------------------------------------------------------------------*/

  public int[] b2b ( boolean[] b )
  {
    int[] c = new int[b.length];

    for ( int i = 0; i < c.length; i++ )
    {
      if ( b[ i ] )
      {
        c[ i ] = 1;
      }

      else
      {
        c[ i ] = 0;
      }
    }

    return c;
  }

/*--------------------------------------------------------------------
 i2b - Convert a boolean array to an integer index array
--------------------------------------------------------------------*/

  public boolean[] i2b ( int[] b )
  {
    return i2b ( b, -1 );
  }

/*--------------------------------------------------------------------
 i2b - Convert a boolean array to an integer index array
--------------------------------------------------------------------*/

  public boolean[] i2b ( int[] b, int size )
  {
    int max = max ( b );

    if ( size <= max )
    {
      size = max + 1;
    }

    boolean[] c = new boolean[size];

    for ( int i = 0; i < c.length; i++ )
    {
      c[ i ] = false;
    }

    for ( int i = 0; i < b.length; i++ )
    {
      c[ b[ i ] ] = true;
    }

    return c;
  }

/*--------------------------------------------------------------------
 order2rank - Change from order to rank
--------------------------------------------------------------------*/

  public int[] order2rank ( int[] z )
  {
    int[] y = new int[z.length];

    for ( int i = 0; i < z.length; i++ )
    {
      y[ z[ i ] ] = i;
    }

    return y;  // Return the index array
  }

/*--------------------------------------------------------------------
 next - Find the next item in an array after a
--------------------------------------------------------------------*/

  public int next ( int[] z, int a )
  {
    int b = -1;  // Holds the index

    // LOOP THROUGH THE ELEMENTS OF THE ARRAY
    for ( int i = 0; i < z.length; i++ )
    {
      // IF THE CURRENT INDEX IS NOT LESS THAN a
      if ( z[ i ] >= a )
      {
        b = z[ i ];  // Store the index
        break;       // Break from the loop
      }
    }

    return b;
  }

/*--------------------------------------------------------------------
 num - Find the ith item in a Text array and convert to an integer
--------------------------------------------------------------------*/

  public int num ( Text[] z, int i )
  {
    // IF THE ARRAY IS NOT NULL AND THE PROPER LENGTH
    if ( z != null && z.length > i )
    {
      return z[ i ].num ();  // Convert the item to an integer
    }

    // OTHERWISE RETURN 0
    else
    {
      return 0;
    }
  }

/*--------------------------------------------------------------------
 dbl - Find the ith item in a Text array and convert to a double
--------------------------------------------------------------------*/

  public double dbl ( Text[] z, int i )
  {
    // IF THE ARRAY IS NOT NULL AND THE PROPER LENGTH
    if ( z != null && z.length > i )
    {
      return z[ i ].dbl ();  // Convert the item to a double
    }

    // OTHERWISE RETURN 0
    else
    {
      return 0;
    }
  }

/*--------------------------------------------------------------------
 txt - Find the ith item in a Text array
--------------------------------------------------------------------*/

  public Text txt ( Text[] z, int i )
  {
    // IF THE ARRAY IS NOT NULL AND THE PROPER LENGTH
    if ( z != null && z.length > i )
    {
      return z[ i ];  // Return the item
    }

    // OTHERWISE RETURN 0
    else
    {
      return new Text ();  // Return a blank text
    }
  }

/*--------------------------------------------------------------------
 str - Find the ith item in a Text array and convert it to a string
--------------------------------------------------------------------*/

  public String str ( Text[] z, int i )
  {
    // IF THE ARRAY IS NOT NULL AND THE PROPER LENGTH
    if ( z != null && z.length > i )
    {
      return z[ i ].str ();  // Convert the item to a string
    }

    // OTHERWISE RETURN 0
    else
    {
      return "";  // Return a blank string
    }
  }

/*--------------------------------------------------------------------
 sub - Select a substring of a character array
--------------------------------------------------------------------*/

  public char[] sub ( char[] z, int c, int d )
  {
    d = lechk ( d, z.length );  // Ensure d <= size
    d = gechk ( d, 0 );         // Ensure d >= 0
    c = lechk ( c, d );         // Ensure c <= d
    c = gechk ( c, 0 );         // Ensure c >= 0
    int newsize = d - c;        // Size of the new array

    char[] a = new char[newsize];  // Create the new array

    // LOOP THROUGH THE ARRAY
    for ( int i = c; i < d; i++ )
    {
      a[ i - c ] = z[ i ];  // Set the ith character
    }

    return a;  // Return the new array
  }

/*--------------------------------------------------------------------
 MATH FUNCTIONS
 ---------------------------------------------------------------------
 Handle the standard mathematical operators
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 sum - Return the sum
--------------------------------------------------------------------*/

  public int sum ( int[] a )
  {
    int sum = 0;

    for ( int i = 0; i < a.length; i++ )
    {
      sum += a[ i ];
    }

    return sum;
  }

/*--------------------------------------------------------------------
 sum - Return the sum
--------------------------------------------------------------------*/

  public double sum ( double[] a )
  {
    double sum = 0;

    for ( int i = 0; i < a.length; i++ )
    {
      sum += a[ i ];
    }

    return sum;
  }

/*--------------------------------------------------------------------
 mean - Return the mean
--------------------------------------------------------------------*/

  public double mean ( double[] a )
  {
    double sum = 0;

    for ( int i = 0; i < a.length; i++ )
    {
      sum += a[ i ];
    }

    return sum / a.length;
  }

/*--------------------------------------------------------------------
 mean - Return the mean of two integer arrays
--------------------------------------------------------------------*/

  public int[] mean ( int[] a, int[] b )
  {
    int[] c = new int[a.length];

    for ( int i = 0; i < a.length; i++ )
    {
      c[ i ] = ( a[ i ] + b[ i ] ) / 2;
    }

    return c;
  }

/*--------------------------------------------------------------------
 median - Return the median
--------------------------------------------------------------------*/

  public double median ( double[] a )
  {
    double[] b = copy ( a );
    sort ( b );
    int i = floor ( a.length / 2 );
    return a[ i ];
  }

/*--------------------------------------------------------------------
 var - Return the variance
--------------------------------------------------------------------*/

  public double var ( double[] a )
  {
    double sum = 0;
    double m = mean ( a );

    for ( int i = 0; i < a.length; i++ )
    {
      double d = a[ i ] - m;
      sum += d * d;
    }

    return ( sum / ( a.length - 1 ) );
  }

/*--------------------------------------------------------------------
 sd - Return the standard deviation
--------------------------------------------------------------------*/

  public double sd ( double[] a )
  {
    return pow ( var ( a ), 0.5 );
  }

/*--------------------------------------------------------------------
 max - Return the maximum value
--------------------------------------------------------------------*/

  public int max ( int a, int b )
  {
    // IF b IS LARGER THAN a
    if ( b > a )
    {
      return b;
    }

    // IF b IS NOT LARGER THAN a
    else
    {
      return a;
    }
  }

/*--------------------------------------------------------------------
 max - Return the maximum value
--------------------------------------------------------------------*/

  public double max ( double a, double b )
  {
    // IF b IS LARGER THAN a
    if ( b > a )
    {
      return b;
    }

    // IF b IS NOT LARGER THAN a
    else
    {
      return a;
    }
  }

/*--------------------------------------------------------------------
 max - Return the maximum value
--------------------------------------------------------------------*/

  public int max ( int[] a )
  {
    int max = 0;

    if ( a.length > 0 )
    {
      max = a[ 0 ];

      for ( int i = 1; i < a.length; i++ )
      {
        if ( a[ i ] > max )
        {
          max = a[ i ];
        }
      }
    }

    return max;
  }

/*--------------------------------------------------------------------
 max - Return the maximum value
--------------------------------------------------------------------*/

  public double max ( double[] a )
  {
    double max = 0;

    if ( a.length > 0 )
    {
      max = a[ 0 ];

      for ( int i = 1; i < a.length; i++ )
      {
        if ( a[ i ] > max )
        {
          max = a[ i ];
        }
      }
    }

    return max;
  }

/*--------------------------------------------------------------------
 maxInd - Return the index of the maximum value
--------------------------------------------------------------------*/

  public int maxInd ( int[] a )
  {
    int max = 0;

    for ( int i = 1; i < a.length; i++ )
    {
      if ( a[ i ] > a[ max ] )
      {
        max = i;
      }
    }

    return max;
  }

/*--------------------------------------------------------------------
 maxInd - Return the index of the maximum value
--------------------------------------------------------------------*/

  public int maxInd ( double[] a )
  {
    int max = 0;

    for ( int i = 1; i < a.length; i++ )
    {
      if ( a[ i ] > a[ max ] )
      {
        max = i;
      }
    }

    return max;
  }

/*--------------------------------------------------------------------
 min - Return the minimum value
--------------------------------------------------------------------*/

  public int min ( int a, int b )
  {
    // IF b IS LESS THAN a
    if ( b < a )
    {
      return b;
    }

    // IF b IS NOT LESS THAN a
    else
    {
      return a;
    }
  }

/*--------------------------------------------------------------------
 min - Return the minimum value
--------------------------------------------------------------------*/

  public double min ( double a, double b )
  {
    // IF b IS LESS THAN a
    if ( b < a )
    {
      return b;
    }

    // IF b IS NOT LESS THAN a
    else
    {
      return a;
    }
  }

/*--------------------------------------------------------------------
 min - Return the minimum value
--------------------------------------------------------------------*/

  public int min ( int[] a )
  {
    int min = a[ 0 ];

    for ( int i = 1; i < a.length; i++ )
    {
      if ( a[ i ] < min )
      {
        min = a[ i ];
      }
    }

    return min;
  }

/*--------------------------------------------------------------------
 min - Return the minimum value
--------------------------------------------------------------------*/

  public double min ( double[] a )
  {
    double min = a[ 0 ];

    for ( int i = 1; i < a.length; i++ )
    {
      if ( a[ i ] < min )
      {
        min = a[ i ];
      }
    }

    return min;
  }

/*--------------------------------------------------------------------
 minInd - Return the index of the minimum value
--------------------------------------------------------------------*/

  public int minInd ( int[] a )
  {
    int min = 0;

    for ( int i = 1; i < a.length; i++ )
    {
      if ( a[ i ] < a[ min ] )
      {
        min = i;
      }
    }

    return min;
  }

/*--------------------------------------------------------------------
 minInd - Return the index of minimum value
--------------------------------------------------------------------*/

  public int minInd ( double[] a )
  {
    int min = 0;

    for ( int i = 1; i < a.length; i++ )
    {
      if ( a[ i ] < a[ min ] )
      {
        min = i;
      }
    }

    return min;
  }

/*--------------------------------------------------------------------
 log - Calculate the log of a
--------------------------------------------------------------------*/

  public double log ( double a )
  {
    return Math.log ( a );
  }

/*--------------------------------------------------------------------
 log - Calculate the log of a
--------------------------------------------------------------------*/

  public double[] log ( double[] a )
  {
    double[] b = new double[a.length];

    for ( int i = 0; i < a.length; i++ )
    {
      b[ i ] = log ( a[ i ] );
    }

    return b;
  }

/*--------------------------------------------------------------------
 exp - Calculate the reverse of the log of a
--------------------------------------------------------------------*/

  public double exp ( double a )
  {
    return Math.exp ( a );
  }

/*--------------------------------------------------------------------
 exp - Calculate the reverse of the log of a
--------------------------------------------------------------------*/

  public double[] exp ( double[] a )
  {
    double[] b = new double[a.length];

    for ( int i = 0; i < a.length; i++ )
    {
      b[ i ] = exp ( a[ i ] );
    }

    return b;
  }

/*--------------------------------------------------------------------
 log2 - Calculate the log base 2 of a
--------------------------------------------------------------------*/

  public double log2 ( double a )
  {
    return Math.log ( a ) / Math.log ( 2 );
  }

/*--------------------------------------------------------------------
 log2 - Calculate the log base 2 of a
--------------------------------------------------------------------*/

  public double[] log2 ( double[] a )
  {
    double[] b = new double[a.length];

    for ( int i = 0; i < a.length; i++ )
    {
      b[ i ] = log2 ( a[ i ] );
    }

    return b;
  }

/*--------------------------------------------------------------------
 log10 - Calculate the log base 10 of a
--------------------------------------------------------------------*/

  public double log10 ( double a )
  {
    return Math.log ( a ) / Math.log ( 10 );
  }

/*--------------------------------------------------------------------
 log10 - Calculate the log base 10 of a
--------------------------------------------------------------------*/

  public double[] log10 ( double[] a )
  {
    double[] b = new double[a.length];

    for ( int i = 0; i < a.length; i++ )
    {
      b[ i ] = log10 ( a[ i ] );
    }

    return b;
  }

/*--------------------------------------------------------------------
 pow - Calculate a to the bth power pow
--------------------------------------------------------------------*/

  public double pow ( double a, double b )
  {
    return Math.pow ( a, b );
  }

/*--------------------------------------------------------------------
 apow - Calculate a to the bth power pow
--------------------------------------------------------------------*/

  public double apow ( double a, double b )
  {
    if ( a < 0 )
    {
      return -Math.pow ( -a, b );
    }

    else
    {
      return Math.pow ( a, b );
    }
  }

/*--------------------------------------------------------------------
 abs - Calculate the absolute value
--------------------------------------------------------------------*/

  public double abs ( double a )
  {
    return Math.abs ( a );
  }

/*--------------------------------------------------------------------
 logChoose - Log version of choose
--------------------------------------------------------------------*/

  public double logChoose ( int a, int b )
  {
    double n = 0;
    double d = 0;

    if ( b == 0 )
    {
      return 0;
    }

    for ( int i = a - b + 1; i <= a; i++ )
    {
      n += log ( i );
    }

    for ( int i = 1; i <= b; i++ )
    {
      d += log ( i );
    }

    return n - d;
  }

/*--------------------------------------------------------------------
 hyperGeometric - Probability of choosing c of category C
                  in n items from total set N
--------------------------------------------------------------------*/

  public double hyperGeometric ( int c, int C, int n, int N )
  {
    double p = 0;
    int max = min ( C, n );

    for ( int i = c; i <= max; i++ )
    {
      p += exp ( logChoose ( C, i ) + logChoose ( N - C, n - i ) - logChoose ( N, n ) );
    }

    // out ( c + " of " + C + ", " + n + " of " + N + " = " + p );

    return p;
  }

/*--------------------------------------------------------------------
 add - Add an integer to an integer array
--------------------------------------------------------------------*/

  public int[] add ( int[] z, int a )
  {
    int[] y = new int[z.length];  // Create a new array

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < z.length; i++ )
    {
      y[ i ] = z[ i ] + a;  // Add the integer
    }

    return y;  // Return the array
  }

/*--------------------------------------------------------------------
 add - Add two integer arrays
--------------------------------------------------------------------*/

  public int[] add ( int[] z, int[] a )
  {
    int[] y = new int[z.length];  // Create a new array

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < z.length; i++ )
    {
      y[ i ] = z[ i ] + a[ i ];  // Add the integer
    }

    return y;  // Return the array
  }

/*--------------------------------------------------------------------
 add - Add a double to an integer array
--------------------------------------------------------------------*/

  public double[] add ( int[] z, double a )
  {
    double[] y = new double[z.length];  // Create a new array

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < z.length; i++ )
    {
      y[ i ] = z[ i ] + a;  // Add the integer
    }

    return y;  // Return the array
  }

/*--------------------------------------------------------------------
 add - Add a double to a double array
--------------------------------------------------------------------*/

  public double[] add ( double[] z, double a )
  {
    double[] y = new double[z.length];  // Create a new array

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < z.length; i++ )
    {
      y[ i ] = z[ i ] + a;  // Add the integer
    }

    return y;  // Return the array
  }

/*--------------------------------------------------------------------
 mul - Multiply an integer array by an integer
--------------------------------------------------------------------*/

  public int[] mul ( int[] z, int d )
  {
    int[] y = new int[z.length];  // Create a new array

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < z.length; i++ )
    {
      y[ i ] = z[ i ] * d;  // Multiply
    }

    return y;  // Return the array
  }

/*--------------------------------------------------------------------
 mul - Multiply an integer array by a double
--------------------------------------------------------------------*/

  public double[] mul ( int[] z, double d )
  {
    double[] y = new double[z.length];  // Create a new array

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < z.length; i++ )
    {
      y[ i ] = ( double ) ( z[ i ] ) * d;  // Divide
    }

    return y;  // Return the array
  }

/*--------------------------------------------------------------------
 mul - Multiply a double array by a double
--------------------------------------------------------------------*/

  public double[] mul ( double[] z, double d )
  {
    double[] y = new double[z.length];  // Create a new array

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < z.length; i++ )
    {
      y[ i ] = ( double ) ( z[ i ] ) * d;  // Divide
    }

    return y;  // Return the array
  }

/*--------------------------------------------------------------------
 div - Divide an integer array by a double
--------------------------------------------------------------------*/

  public double[] div ( int[] z, double d )
  {
    double[] y = new double[z.length];  // Create a new array

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < z.length; i++ )
    {
      y[ i ] = ( double ) ( z[ i ] ) / d;  // Divide
    }

    return y;  // Return the array
  }

/*--------------------------------------------------------------------
 div - Divide a double array by a double
--------------------------------------------------------------------*/

  public double[] div ( double[] z, double d )
  {
    double[] y = new double[z.length];  // Create a new array

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < z.length; i++ )
    {
      y[ i ] = ( double ) ( z[ i ] ) / d;  // Divide
    }

    return y;  // Return the array
  }

/*--------------------------------------------------------------------
 ceil - Ceiling
--------------------------------------------------------------------*/

  public int ceil ( double z )
  {
    return ( int ) ( Math.ceil ( z ) );
  }

/*--------------------------------------------------------------------
 ceil - Ceiling
--------------------------------------------------------------------*/

  public int[] ceil ( double[] z )
  {
    int[] y = new int[z.length];  // Create a new array

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < z.length; i++ )
    {
      y[ i ] = ( int ) ( Math.ceil ( z[ i ] ) );
    }

    return y;  // Return the array
  }

/*--------------------------------------------------------------------
 floor - Floor
--------------------------------------------------------------------*/

  public int floor ( double z )
  {
    return ( int ) ( Math.floor ( z ) );
  }

/*--------------------------------------------------------------------
 floor - Floor
--------------------------------------------------------------------*/

  public int[] floor ( double[] z )
  {
    int[] y = new int[z.length];  // Create a new array

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < z.length; i++ )
    {
      y[ i ] = ( int ) ( Math.floor ( z[ i ] ) );
    }

    return y;  // Return the array
  }

/*--------------------------------------------------------------------
 round - Round a number
--------------------------------------------------------------------*/

  public double round ( double z )
  {
    return Math.round ( z );
  }

/*--------------------------------------------------------------------
 round - Round a number to the specified number of digits
--------------------------------------------------------------------*/

  public double round ( double z, double d )
  {
    double n = pow ( 10, -d );
    return Math.round ( z * n ) / n;
  }

/*--------------------------------------------------------------------
 SORT FUNCTIONS
 ---------------------------------------------------------------------
 Sort arrays of various datatypes
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 find - Find the index of the item closest to the left of a
--------------------------------------------------------------------*/

  public int find ( double[] z, double a )
  {
    int b = -1;

    for ( int i = 0; i < z.length && z[ i ] < a; i++ )
    {
      b++;
    }

    if ( z[ z.length - 1 ] < a )
    {
      b = z.length;
    }

    return b;
  }

/*--------------------------------------------------------------------
 findd - Find the index of the item closest to the left of a
--------------------------------------------------------------------*/

  public int findd ( double[] z, double a )
  {
    int b = -1;

    for ( int i = 0; i < z.length && z[ i ] >= a; i++ )
    {
      b++;
    }

    if ( z[ z.length - 1 ] >= a )
    {
      b = z.length;
    }

    return b;
  }

/*--------------------------------------------------------------------
 sort - Sort the indices of an array of integers
--------------------------------------------------------------------*/

  public int[] sort ( int[] z )
  {
    int size = z.length;        // Size of the array
    int[] a = new int[size];  // First index array
    int[] b = new int[size];  // Second index array

    // INITIALIZE THE ARRAYS
    for ( int i = 0; i < size; i = i + 1 )
    {
      a[ i ] = i;  // Set the ith item to i
      b[ i ] = i;  // Set the ith item to i
    }

    int[] y = new int[size];                 // Temporary integer array
    mergeSort ( a, b, z, y, 0, size - 1 );  // Mergse sort the array
    return a;                                  // Return the indices
  }

/*--------------------------------------------------------------------
 mergeSort - Merge sort the indices of an array of integers
--------------------------------------------------------------------*/

  public void mergeSort ( int[] a, int[] b, int[] p, int[] q, int min, int max )
  {
    // FIND THE SIZE OF THE CURRENT SUBARRAY
    int size = max - min + 1;

    // IF THERE ARE MORE THAN TWO ELEMENTS
    if ( size > 2 )
    {
      // SPLIT THE ARRAY
      int s = min + ( size / 2 );               // Find the split point
      mergeSort ( a, b, p, q, min, s - 1 );  // Merge sort the lower half
      mergeSort ( a, b, p, q, s, max );      // Merge sort the upper half

      // INITIALIZE THE LOOPING INDICES
      int i = min;  // Index for the first split
      int j = s;    // Index for the second split
      int k = min;  // Index for the new ordering

      // WHILE THERE ARE COMPARISONS TO MAKE
      while ( ( i < s ) || ( j <= max ) )
      {
        // IF THE FIRST AND SECOND ARRAYS EACH HAVE ITEMS LEFT
        if ( ( i < s ) && ( j <= max ) )
        {
          // IF THE FIRST IS LESS THAN OR EQUAL TO THE SECOND
          if ( ( q[ i ] ) <= ( q[ j ] ) )
          {
            p[ k ] = q[ i ];  // Store the first value
            a[ k ] = b[ i ];  // Store the first index
            i = i + 1;        // Increment the i index
          }

          // IF THE FIRST IS GREATER THAN THE SECOND
          else
          {
            p[ k ] = q[ j ];  // Store the second value
            a[ k ] = b[ j ];  // Store the second index
            j = j + 1;        // Increment the j index
          }
        }

        // IF I IS LESS THAN THE SPLIT
        else if ( i < s )
        {
          p[ k ] = q[ i ];  // Store the first value
          a[ k ] = b[ i ];  // Store the first index
          i = i + 1;        // Increment the i index
        }

        // IF I IS AT THE SPLIT
        else
        {
          p[ k ] = q[ j ];  // Store the second value
          a[ k ] = b[ j ];  // Store the second index
          j = j + 1;        // Increment the j index
        }

        k = k + 1;  // Increment the k index
      }

      // LOOP FROM MIN TO MAX
      for ( int h = 0; h <= max; h = h + 1 )
      {
        q[ h ] = p[ h ];  // Store the value in q
        b[ h ] = a[ h ];  // Store the index in b
      }
    }

    // IF THERE ARE TWO ITEMS IN THE ARRAY
    else if ( size == 2 )
    {
      // IF THE FIRST ITEM IS LESS THAN OR EQUAL TO THE SECOND
      if ( ( p[ min ] ) <= ( p[ min + 1 ] ) )
      {
        q[ min ] = p[ min ];          // Store the first value in q
        q[ min + 1 ] = p[ min + 1 ];  // Store the second value in q
        b[ min ] = a[ min ];          // Store the first index in b
        b[ min + 1 ] = a[ min + 1 ];  // Store the second index in b
      }

      // IF THE FIRST ITEM IS GREATER THAN THE SECOND
      else
      {
        // SWITCH THE VALUES IN Q
        q[ min ] = p[ min + 1 ];
        q[ min + 1 ] = p[ min ];

        // STORE THE VALUES BACK IN P
        p[ min ] = q[ min ];
        p[ min + 1 ] = q[ min + 1 ];

        // SWITCH THE INDICES IN B
        b[ min ] = a[ min + 1 ];
        b[ min + 1 ] = a[ min ];

        // STORE THE INDICES BACK IN A
        a[ min ] = b[ min ];
        a[ min + 1 ] = b[ min + 1 ];
      }
    }

    // IF THERE IS ONLY ONE ITEM IN THE ARRAY
    else if ( size == 1 )
    {
      q[ min ] = p[ min ];  // Store the value in q
      b[ min ] = a[ min ];  // Store the index in b
    }
  }

/*--------------------------------------------------------------------
 sort - Sort the indices of an array of doubles
--------------------------------------------------------------------*/

  public int[] sort ( double[] z )
  {
    int size = z.length;        // Size of the array
    int[] a = new int[size];  // First index array
    int[] b = new int[size];  // Second index array

    // INITIALIZE THE ARRAYS
    for ( int i = 0; i < size; i = i + 1 )
    {
      a[ i ] = i;  // Set the ith item to i
      b[ i ] = i;  // Set the ith item to i
    }

    double[] y = new double[size];           // Temporary double array
    mergeSort ( a, b, z, y, 0, size - 1 );  // Mergse sort the array
    return a;                                  // Return the indices
  }

/*--------------------------------------------------------------------
 sortd - Sort the indices of an array of doubles descending
--------------------------------------------------------------------*/

  public int[] sortd ( double[] z )
  {
    int[] ind = sort ( z );
    ind = rev ( ind );

    double t = 0;
    for ( int i = 0; i < z.length / 2; i++ )
    {
      t = z[ z.length - i - 1 ];
      z[ z.length - i - 1 ] = z[ i ];
      z[ i ] = t;
    }

    return ind;
  }

/*--------------------------------------------------------------------
 sum - Sum an array of booleans
--------------------------------------------------------------------*/

  public int sum ( boolean[] z )
  {
    int sum = 0;

    for ( int i = 0; i < z.length; i++ )
    {
      if ( z[ i ] )
      {
        sum++;
      }
    }

    return sum;
  }

/*--------------------------------------------------------------------
 mergeSort - Merge sort the indices of an array of doubles
--------------------------------------------------------------------*/

  public void mergeSort ( int[] a, int[] b, double[] p, double[] q, int min, int max )
  {
    // FIND THE SIZE OF THE CURRENT SUBARRAY
    int size = max - min + 1;

    // IF THERE ARE MORE THAN TWO ELEMENTS
    if ( size > 2 )
    {
      // SPLIT THE ARRAY
      int s = min + ( size / 2 );               // Find the split point
      mergeSort ( a, b, p, q, min, s - 1 );  // Merge sort the lower half
      mergeSort ( a, b, p, q, s, max );      // Merge sort the upper half

      // INITIALIZE THE LOOPING INDICES
      int i = min;  // Index for the first split
      int j = s;    // Index for the second split
      int k = min;  // Index for the new ordering

      // WHILE THERE ARE COMPARISONS TO MAKE
      while ( ( i < s ) || ( j <= max ) )
      {
        // IF THE FIRST AND SECOND ARRAYS EACH HAVE ITEMS LEFT
        if ( ( i < s ) && ( j <= max ) )
        {
          // IF THE FIRST IS LESS THAN OR EQUAL TO THE SECOND
          if ( ( q[ i ] ) <= ( q[ j ] ) )
          {
            p[ k ] = q[ i ];  // Store the first value
            a[ k ] = b[ i ];  // Store the first index
            i = i + 1;        // Increment the i index
          }

          // IF THE FIRST IS GREATER THAN THE SECOND
          else
          {
            p[ k ] = q[ j ];  // Store the second value
            a[ k ] = b[ j ];  // Store the second index
            j = j + 1;        // Increment the j index
          }
        }

        // IF I IS LESS THAN THE SPLIT
        else if ( i < s )
        {
          p[ k ] = q[ i ];  // Store the first value
          a[ k ] = b[ i ];  // Store the first index
          i = i + 1;        // Increment the i index
        }

        // IF I IS AT THE SPLIT
        else
        {
          p[ k ] = q[ j ];  // Store the second value
          a[ k ] = b[ j ];  // Store the second index
          j = j + 1;        // Increment the j index
        }

        k = k + 1;  // Increment the k index
      }

      // LOOP FROM MIN TO MAX
      for ( int h = 0; h <= max; h = h + 1 )
      {
        q[ h ] = p[ h ];  // Store the value in q
        b[ h ] = a[ h ];  // Store the index in b
      }
    }

    // IF THERE ARE TWO ITEMS IN THE ARRAY
    else if ( size == 2 )
    {
      // IF THE FIRST ITEM IS LESS THAN OR EQUAL TO THE SECOND
      if ( ( p[ min ] ) <= ( p[ min + 1 ] ) )
      {
        q[ min ] = p[ min ];          // Store the first value in q
        q[ min + 1 ] = p[ min + 1 ];  // Store the second value in q
        b[ min ] = a[ min ];          // Store the first index in b
        b[ min + 1 ] = a[ min + 1 ];  // Store the second index in b
      }

      // IF THE FIRST ITEM IS GREATER THAN THE SECOND
      else
      {
        // SWITCH THE VALUES IN Q
        q[ min ] = p[ min + 1 ];
        q[ min + 1 ] = p[ min ];

        // STORE THE VALUES BACK IN P
        p[ min ] = q[ min ];
        p[ min + 1 ] = q[ min + 1 ];

        // SWITCH THE INDICES IN B
        b[ min ] = a[ min + 1 ];
        b[ min + 1 ] = a[ min ];

        // STORE THE INDICES BACK IN A
        a[ min ] = b[ min ];
        a[ min + 1 ] = b[ min + 1 ];
      }
    }

    // IF THERE IS ONLY ONE ITEM IN THE ARRAY
    else if ( size == 1 )
    {
      q[ min ] = p[ min ];  // Store the value in q
      b[ min ] = a[ min ];  // Store the index in b
    }
  }

/*--------------------------------------------------------------------
 sort - Sort the indices
--------------------------------------------------------------------*/

  public int[] sort ( Text[] z )
  {
    int size = z.length;        // Size of the array
    int[] a = new int[size];  // First index array
    int[] b = new int[size];  // Second index array

    // INITIALIZE THE ARRAYS
    for ( int i = 0; i < size; i = i + 1 )
    {
      a[ i ] = i;  // Set the ith item to i
      b[ i ] = i;  // Set the ith item to i
    }

    Text[] y = new Text[size];            // Temporary Text array
    mergeSort ( a, b, z, y, 0, size - 1 );  // Mergse sort the array
    return a;                               // Return the indices
  }

/*--------------------------------------------------------------------
 mergeSort - Merge sort the indices of an array of Texts
--------------------------------------------------------------------*/

  public void mergeSort ( int[] a, int[] b, Text[] p, Text[] q, int min, int max )
  {
    // FIND THE SIZE OF THE CURRENT SUBARRAY
    int size = max - min + 1;

    // IF THERE ARE MORE THAN TWO ELEMENTS
    if ( size > 2 )
    {
      // SPLIT THE ARRAY
      int s = min + ( size / 2 );            // Find the split point
      mergeSort ( a, b, p, q, min, s - 1 );  // Merge sort the lower half
      mergeSort ( a, b, p, q, s, max );      // Merge sort the upper half

      // INITIALIZE THE LOOPING INDICES
      int i = min;  // Index for the first split
      int j = s;    // Index for the second split
      int k = min;  // Index for the new ordering

      // WHILE THERE ARE COMPARISONS TO MAKE
      while ( ( i < s ) || ( j <= max ) )
      {
        // IF THE FIRST AND SECOND ARRAYS EACH HAVE ITEMS LEFT
        if ( ( i < s ) && ( j <= max ) )
        {
          // IF THE FIRST IS LESS THAN OR EQUAL TO THE SECOND
          if ( q[ i ].le ( q[ j ] ) )
          {
            p[ k ] = q[ i ];  // Store the first value
            a[ k ] = b[ i ];  // Store the first index
            i = i + 1;        // Increment the i index
          }

          // IF THE FIRST IS GREATER THAN THE SECOND
          else
          {
            p[ k ] = q[ j ];  // Store the second value
            a[ k ] = b[ j ];  // Store the second index
            j = j + 1;        // Increment the j index
          }
        }

        // IF I IS LESS THAN THE SPLIT
        else if ( i < s )
        {
          p[ k ] = q[ i ];  // Store the first value
          a[ k ] = b[ i ];  // Store the first index
          i = i + 1;        // Increment the i index
        }

        // IF I IS AT THE SPLIT
        else
        {
          p[ k ] = q[ j ];  // Store the second value
          a[ k ] = b[ j ];  // Store the second index
          j = j + 1;        // Increment the j index
        }

        k = k + 1;  // Increment the k index
      }

      // LOOP FROM MIN TO MAX
      for ( int h = 0; h <= max; h = h + 1 )
      {
        q[ h ] = p[ h ];  // Store the value in q
        b[ h ] = a[ h ];  // Store the index in b
      }
    }

    // IF THERE ARE TWO ITEMS IN THE ARRAY
    else if ( size == 2 )
    {
      // IF THE FIRST ITEM IS LESS THAN OR EQUAL TO THE SECOND
      if ( p[ min ].le ( p[ min + 1 ] ) )
      {
        q[ min ] = p[ min ];          // Store the first value in q
        q[ min + 1 ] = p[ min + 1 ];  // Store the second value in q
        b[ min ] = a[ min ];          // Store the first index in b
        b[ min + 1 ] = a[ min + 1 ];  // Store the second index in b
      }

      // IF THE FIRST ITEM IS GREATER THAN THE SECOND
      else
      {
        // SWITCH THE VALUES IN Q
        q[ min ] = p[ min + 1 ];
        q[ min + 1 ] = p[ min ];

        // STORE THE VALUES BACK IN P
        p[ min ] = q[ min ];
        p[ min + 1 ] = q[ min + 1 ];

        // SWITCH THE INDICES IN B
        b[ min ] = a[ min + 1 ];
        b[ min + 1 ] = a[ min ];

        // STORE THE INDICES BACK IN A
        a[ min ] = b[ min ];
        a[ min + 1 ] = b[ min + 1 ];
      }
    }

    // IF THERE IS ONLY ONE ITEM IN THE ARRAY
    else if ( size == 1 )
    {
      q[ min ] = p[ min ];  // Store the value in q
      b[ min ] = a[ min ];  // Store the index in b
    }
  }

/*--------------------------------------------------------------------
 SYSTEM FUNCTIONS
 ---------------------------------------------------------------------
 Calls various system functions
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 gc - Garbage collect
--------------------------------------------------------------------*/

  public void gc ()
  {
    Runtime.getRuntime ().gc ();  // Garbage collection
  }

/*--------------------------------------------------------------------
 freemem - Return the free memory in megabytes
--------------------------------------------------------------------*/

  public double freemem ()
  {
    return Runtime.getRuntime ().freeMemory () / ( pow ( 2, 20 ) );
  }

/*--------------------------------------------------------------------
 totalmem - Return the total memory in megabytes
--------------------------------------------------------------------*/

  public double totalmem ()
  {
    return Runtime.getRuntime ().totalMemory () / ( pow ( 2, 20 ) );
  }

/*--------------------------------------------------------------------
 maxmem - Return the maximum memory in megabytes
--------------------------------------------------------------------*/

  public double maxmem ()
  {
    return Runtime.getRuntime ().maxMemory () / ( pow ( 2, 20 ) );
  }

/*--------------------------------------------------------------------
 n2p - Negative to positive points
--------------------------------------------------------------------*/

  public double n2p ( double z )
  {
    return 0.5 + z * 0.5;
  }

  /**
   * Blend two colors
   */

  /**
   * Blend two colors
   *
   * @param a first color
   * @param b second color
   * @param i blend amount
   *
   * @return blended color
   */

  public double[] blend ( double[] a, double[] b, double i )
  {
    double[] c = new double[3];

    if ( i < 0 )
    {
      i = 0;
    }

    if ( i > 1 )
    {
      i = 1;
    }

    for ( int j = 0; j < 3; j++ )
    {
      c[ j ] = i * ( b[ j ] - a[ j ] ) + a[ j ];
    }

    return c;
  }
}
