package broad.projection.utils;/*--------------------------------------------------------------------
 com.scanfeld.core.RC4
 ---------------------------------------------------------------------
 Output random numbers following the com.scanfeld.core.RC4 algorithm.
 ---------------------------------------------------------------------
 Author: Daniel Scanfeld
   Date: 08.23.06
 ---------------------------------------------------------------------

--------------------------------------------------------------------*/

public class RC4 extends Inc
{
/*--------------------------------------------------------------------
 Class Variables
--------------------------------------------------------------------*/

  public char[] RC4S;            // Main character array
  public int RC4i;               // i index
  public int RC4j;               // j index
  public double saveNormal;      // Save a normal
  public boolean flagNormal;     // Flag for normals
  public String seed;            // Seed for the algorithm
  public static int doubleSize;  // Size of randoms
  public int ct;                 // Number of randoms

/*--------------------------------------------------------------------
 CONSTRUCTOR FUNCTIONS
 ---------------------------------------------------------------------
 Constructors for the possible com.scanfeld.core.RC4 objects
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 com.scanfeld.core.RC4 - Null constructor
--------------------------------------------------------------------*/

  public RC4 ()
  {
    init ();
  }

/*--------------------------------------------------------------------
 com.scanfeld.core.RC4 - Initialize with a string
--------------------------------------------------------------------*/

  public RC4 ( String s )
  {
    init ( s );
  }

/*--------------------------------------------------------------------
 init - Initialization function
--------------------------------------------------------------------*/

  public void init ()
  {
    init ( date () );          // Initialize to the date
    init ( randString ( 10 ) );  // Then initialize to a random 10 character word
  }

/*--------------------------------------------------------------------
 init - Initialization function
--------------------------------------------------------------------*/

  public void init ( String s )
  {
    RC4S = new char[256];  // Set up the character array
    doubleSize = 3;          // Number of bytes in random doubles
    set ( s );               // Initialize to the seed s
  }

/*--------------------------------------------------------------------
 set - Set the seed for the com.scanfeld.core.RC4 object with an index
--------------------------------------------------------------------*/

  public void set ( String s, int i )
  {
    set ( s );
    
    while ( ct < i )
    {
      rand ();
    }
  }

/*--------------------------------------------------------------------
 seedRandom - Seed the java random number generator
--------------------------------------------------------------------*/

  public void seedRandom ( )
  {
    // Math.random("dan");
  }

/*--------------------------------------------------------------------
 set - Set the seed for the com.scanfeld.core.RC4 object
--------------------------------------------------------------------*/

  public void set ( String s )
  {
    seed = s;                    // SET THE SEED TO S
    flagNormal = false;          // NOT CURRENTLY MAKING RANDOM NUMBERS FROM A NORMAL DISTRIBUTION
    saveNormal = 0.0;            // THE LAST NORMAL NUMBER CREATED
    char[] k = new char[256];    // CREATE A TEMPORARY CHARACTER ARRAY
    char temp;                   // CREATE A TEMPORARY CHARACTER VARIABLE

    // DECLARE LOOP VARIABLES
    int i;
    int j = 0;

    // LOOP THROUGH THE RC4S ARRAY
    for ( i = 0; i < 256; i = i + 1 )
    {
      RC4S[ i ] = ( char ) ( i );  // INITIALIZE THE ITH ITEM TO I
      k[ i ] = seed.charAt ( j );  // FILL THE K ARRAY WITH REPEATS OF THE SEED
      j = j + 1;                   // INCREMENT J

      // IF AT THE END OF THE SEED
      if ( j >= seed.length () )
      {
        j = 0;  // GO TO THE FIRST CHARACTER OF THE SEED
      }
    }

    j = 0;

    // LOOP THROUGH THE RC4S ARRAY
    for ( i = 0; i < 256; i = i + 1 )
    {
      j = ( j + ( RC4S[ i ] ) + ( k[ i ] ) ) & 0xff;  // Update j

      // SWITCH TWO ITEMS
      temp = RC4S[ i ];
      RC4S[ i ] = RC4S[ j ];
      RC4S[ j ] = temp;
    }

    RC4i = 0;
    RC4j = 0;

    ct = 0;  // Initialize the count to 0;

    out ( "Random: " + seed );
  }

/*--------------------------------------------------------------------
 rand - Return a random byte
--------------------------------------------------------------------*/

  public char rand ()
  {
    ct++;

    char temp;                                  // Temporary character
    RC4i = ( RC4i + 1 ) & 0xff;                 // Update RC4i
    RC4j = ( RC4j + ( RC4S[ RC4i ] ) ) & 0xff;  // Update RC4j

    // SWITCH TWO ITEMS
    temp = RC4S[ RC4i ];
    RC4S[ RC4i ] = RC4S[ RC4j ];
    RC4S[ RC4j ] = temp;

    // RETURN A RANDOM BYTE
    return RC4S[ ( ( RC4S[ RC4i ] ) + ( RC4S[ RC4j ] ) ) & 0xff ];
  }

/*--------------------------------------------------------------------
 randString - Return a random string
--------------------------------------------------------------------*/

  public String randString ( int size )
  {
    // CREATE A STRING VARIABLE
    String s = "";

    // CREATE A CHARACTER LOOKUP STRING
    String c = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

    // FOR EACH NEW CHARACTER
    for ( int i = 0; i < size; i = i + 1 )
    {
      // ADD A CHARACTER TO THE RANDOM STRING
      s = s + c.charAt ( ( int ) ( randDouble () * ( double ) ( c.length () ) ) );
    }

    // RETURN THE RANDOM STRING
    return s;
  }

/*--------------------------------------------------------------------
 randNormal - Create a random independent normally distributed number
   with mean m and standard deviation s
--------------------------------------------------------------------*/

  public double randNormal ( double m, double s )
  {
    // IF THERE IS A SAVED NORMAL
    if ( flagNormal )
    {
      flagNormal = false;
      return ( saveNormal * s ) + m;
    }

    // OTHERWISE, CREATE TWO RANDOM NORMALS
    else
    {
      // MAKE SOME TEMPORARY VARIABLES
      double x = 0.0;
      double y = 0.0;
      double w = 1.0;

      // WHILE W IS GREATER THAN OR EQUAL TO ONE
      while ( w >= 1.0 )
      {
        // INCREMENT THE VARIABLES
        x = ( 2.0 * randDouble () ) - 1.0;
        y = ( 2.0 * randDouble () ) - 1.0;
        w = ( x * x ) + ( y * y );
      }

      w = pow ( -2.0 * log ( w ) / w, 0.5 );  // Update w
      saveNormal = x * w;                     // Save the next normal
      flagNormal = true;                      // Set the flag that a normal is available
      return ( y * w * s ) + m;               // Return the random normal
    }
  }

/*--------------------------------------------------------------------
 randNormal - Create an array of normal doubles of size n with mean m
   and standard deviation s
--------------------------------------------------------------------*/

  public double[] randNormal ( double m, double s, int n )
  {
    // CREATE THE RANDOM ARRAY
    double[] d = new double[n];

    // LOOP THROUGH THE N NUMBERS
    for ( int i = 0; i < n; i = i + 1 )
    {
      d[ i ] = randNormal ( m, s );  // Find the ith random double
    }

    // RETURN THE RANDOM ARRAY
    return d;
  }

/*--------------------------------------------------------------------
 randDouble- Create a random double between 0 and 1.  doubleSize is
   the nuber of bytes used to encode the double.
--------------------------------------------------------------------*/

  public double randDouble ()
  {
    double randDouble = 0.0;  // Hold the random double

    // LOOP THROUGH THE ALLOTED DOUBLESIZE
    for ( int i = 0; i < doubleSize; i = i + 1 )
    {
      randDouble = randDouble + ( rand () << ( i * 8 ) );  // Add the next byte of double
    }

    return randDouble / ( pow ( 256, doubleSize ) );  // Return the random double
  }

/*--------------------------------------------------------------------
 randDouble - Create a random double array of size n between 0 and 1
--------------------------------------------------------------------*/

  public double[] randDouble ( int n )
  {
    double[] d = new double[n];  // The random array

    // LOOP THROUGH THE n NUMBERS
    for ( int i = 0; i < n; i = i + 1 )
    {
      d[ i ] = randDouble ();  // Find the ith random double
    }

    return d;  // Return the random array
  }

/*--------------------------------------------------------------------
 randInt - Create a random integer between min and max
--------------------------------------------------------------------*/

  public int randInt ( int min, int max )
  {
    double size = max - min + 1;                         // The range of the integers
    int n = ( ( int ) ( randDouble () * size ) ) + min;  // Find the random integer
    n = lechk ( n, max );                                 // Ensure n is less than max
    return n;                                            // Return the random integer
  }

/*--------------------------------------------------------------------
 randInt - Create a random array of integers 0 to k - 1
--------------------------------------------------------------------*/

  public int[] randInt ( int k )
  {
    int[] d = new int[k];          // Initialize a random array
    boolean[] b = new boolean[k];  // Initialize a flag array

    // LOOP THROUGH ALL THE FLAGS
    for ( int i = 0; i < k; i = i + 1 )
    {
      b[ i ] = false;  // Set the flag to false
    }

    int size = k;  // Initialize size to k
    int a = 0;     // Temporary int holder

    // LOOP THROUGH RANDOM INTEGERS
    for ( int i = 0; i < k; i = i + 1 )
    {
      a = randInt ( 0, size - 1 );  // Create a random integer

      // LOOP THROUGH RANDOM INTEGERS
      for ( int j = 0; j <= a; j = j + 1 )
      {
        // IF THE J FLAG IS SET
        if ( b[ j ] )
        {
          a++;
        }
      }

      b[ a ] = true;  // Set the flag for this integer
      d[ i ] = a;     // Insert the random integer into the array
      size--;         // Decrement the size
    }

    return d;  // Return the random integer array
  }

/*--------------------------------------------------------------------
 randIntArray - Create a random array of n integers 0 to k - 1
--------------------------------------------------------------------*/

  public int[] randIntArray ( int k, int n )
  {
    int[] z = new int[n];  // Initialize a random array
    int l = 0;

    while ( l < n )
    {
      int[] y = randInt ( k );
      copy ( z, l, y, 0 );
      l += k;
    }

    return z;  // Return the random integer array
  }

/*--------------------------------------------------------------------
 randInt - Create a random array of n integers min to max
--------------------------------------------------------------------*/

  public int[] randInt ( int min, int max, int n )
  {
    int[] d = new int[n];  // Initialize a random array

    // LOOP THROUGH N INTEGERS
    for ( int i = 0; i < n; i = i + 1 )
    {
      d[ i ] = randInt ( min, max );  // Set the ith integer
    }

    return d;  // Return the random integer array
  }
}
