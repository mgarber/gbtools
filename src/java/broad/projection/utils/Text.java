package broad.projection.utils;

/*--------------------------------------------------------------------
 TODO
 ---------------------------------------------------------------------
 output exp doubles
 output formatted doubles
 left, right
 subleft, subright
 findrev, rev
 switch to ml, mr, set them from current inc settings
 trim, trimWhite
 alpha, numeric, alphaNumeric
 removeExtraWhite, removeWhites
 get dbl to handle exponents and stuff 1.34e34
 make a double split function that returns com.scanfeld.core.Text[][] or com.scanfeld.core.Texts[]
 add ability to remove things from com.scanfeld.core.Text
 ---------------------------------------------------------------------
 IDEAS
 ---------------------------------------------------------------------
 ---------------------------------------------------------------------
 BUGS
 ---------------------------------------------------------------------
--------------------------------------------------------------------*/

/**
 * Stores and manipulates a large amount of text
 *
 * @author Daniel Scanfeld
 * @version 1.0
 */

public class Text extends Inc implements Compare
{
/*--------------------------------------------------------------------
 Class Variables
--------------------------------------------------------------------*/

  /** Character storage */
  public char[] x;

  /** Index of the first character */
  public int k;

  /** Size of the text */
  public int size;

  /** Size of the text array */
  public int space;

  /** Amount to multiply space on size increase */
  public double m;

/*--------------------------------------------------------------------
 CONSTRUCTOR FUNCTIONS
--------------------------------------------------------------------*/

  /** Null constructor */

  public Text ()
  {
    init ();     // Initialize the Text object
    set ( "" );  // Set to the empty string
  }

  /**
   * String constructor
   *
   * @param z String to set to
   */

  public Text ( String z )
  {
    init ();     // Initialize the Text object
    set ( z );   // Set to the supplied string
  }

  /**
   * Text constructor
   *
   * @param z Text to set to
   */

  public Text ( Text z )
  {
    init ();     // Initialize the Text object
    set ( z );   // Set to the supplied Text object
  }

  /**
   * Integer constructor
   *
   * @param z int to set to
   */

  public Text ( int z )
  {
    init ();     // Initialize the Text object
    set ( z );   // Set to the supplied integer
  }

  /**
   * Double constructor
   *
   * @param z double to set to
   */

  public Text ( double z )
  {
    init ();     // Initialize the Text object
    set ( z );   // Set to the supplied double
  }

  /**
   * Character array constructor
   *
   * @param z character array to set to
   */

  public Text ( char[] z )
  {
    init ();     // Initialize the Text object
    set ( z );   // Set to the supplied double
  }

/*--------------------------------------------------------------------
 TEXT COMPARISON FUNCTIONS
--------------------------------------------------------------------*/

  /**
   * Compares two text objects
   *
   * @param a Text to compare to
   *
   * @return -1 if less than a, 0 if equal to a, 1 if greater than a
   */

  public int cmp ( Text a )
  {
    int ans = 0;

    // LOOP THROUGH THE CHARACTERS
    for ( int i = 0; i < a.size && i < size; i++ )
    {
      // IF A CHARACTER DOES NOT MATCH
      if ( x[ k + i ] != a.x[ a.k + i ] )
      {
        // IF THE CHARACTER IS LESS THAN IN a
        if ( x[ k + i ] < a.x[ a.k + i ] )
        {
          ans = -1;  // Return less than
        }

        // IF THE CHARACTER IS GREATER THAN IN a
        else
        {
          ans = 1;  // Return greater than
        }

        break;  // Break out of the loop
      }
    }

    // IF NO DIFFERENCE WAS FOUND
    if ( ans == 0 )
    {
      // IF SMALLER THAN a
      if ( size < a.size )
      {
        ans = -1;  // Return less than
      }

      // IF LARGER THAN a
      else if ( a.size < size )
      {
        ans = 1;  // Return greater than
      }
    }

    return ans;  // Return the answer
  }

  /**
   * Is a equal to this Text object?
   *
   * @param a Text to compare to
   *
   * @return true if equal to a, or false if not equal to a
   */

  public boolean eq ( Text a )
  {
    return ( cmp ( a ) == 0 );
  }

  /**
   * Is a not equal to this Text object?
   *
   * @param a Text to compare to
   *
   * @return true if not equal to a, or false if equal to a
   */

  public boolean ne ( Text a )
  {
    return ( cmp ( a ) != 0 );
  }

  /**
   * Is a less than this Text object?
   *
   * @param a Text to compare to
   *
   * @return true if less than a, or false if not less than a
   */

  public boolean lt ( Text a )
  {
    return ( cmp ( a ) < 0 );
  }

  /**
   * Is a less than or equal to this Text object?
   *
   * @param a Text to compare to
   *
   * @return true if less than or equal to a, or false if not less than or equal to a
   */

  public boolean le ( Text a )
  {
    return ( cmp ( a ) <= 0 );
  }

  /**
   * Is a greater than this Text object?
   *
   * @param a Text to compare to
   *
   * @return true if greater than a, or false if not greater than a
   */

  public boolean gt ( Text a )
  {
    return ( cmp ( a ) > 0 );
  }

  /**
   * Is a greater than or equal to this Text object?
   *
   * @param a Text to compare to
   *
   * @return true if greater than or equal to a, or false if not greater than or equal to a
   */

  public boolean ge ( Text a )
  {
    return ( cmp ( a ) >= 0 );
  }

/*--------------------------------------------------------------------
 STRING COMPARISON FUNCTIONS
--------------------------------------------------------------------*/

  /**
   * Compares two string objects
   *
   * @param a String to compare to
   *
   * @return -1 if less than a, 0 if equal to a, 1 if greater than a
   */

  public int cmp ( String a )
  {
    return cmp ( s2t ( a ) );  // Compare using a Text object
  }

  /**
   * Is a equal to this Text object?
   *
   * @param a String to compare to
   *
   * @return true if equal to a, or false if not equal to a
   */

  public boolean eq ( String a )
  {
    return ( cmp ( a ) == 0 );
  }

  /**
   * Is a not equal to this Text object?
   *
   * @param a String to compare to
   *
   * @return true if not equal to a, or false if equal to a
   */

  public boolean ne ( String a )
  {
    return ( cmp ( a ) != 0 );
  }

  /**
   * Is a less than this Text object?
   *
   * @param a String to compare to
   *
   * @return true if less than a, or false if not less than a
   */

  public boolean lt ( String a )
  {
    return ( cmp ( a ) < 0 );
  }

  /**
   * Is a less than or equal to this Text object?
   *
   * @param a String to compare to
   *
   * @return true if less than or equal to a, or false if not less than or equal to a
   */

  public boolean le ( String a )
  {
    return ( cmp ( a ) <= 0 );
  }

  /**
   * Is a greater than this Text object?
   *
   * @param a String to compare to
   *
   * @return true if greater than a, or false if not greater than a
   */

  public boolean gt ( String a )
  {
    return ( cmp ( a ) > 0 );
  }

  /**
   * Is a greater than or equal to this Text object?
   *
   * @param a String to compare to
   *
   * @return true if greater than or equal to a, or false if not greater than or equal to a
   */

  public boolean ge ( String a )
  {
    return ( cmp ( a ) >= 0 );
  }

/*--------------------------------------------------------------------
 HELPER FUNCTIONS
--------------------------------------------------------------------*/

  /** Initialize the object */

  public void init ()
  {
    m = 2.0;                // Multiply amount
    k = 10;                 // Initialize to center of the array
    space = 20;             // Initialize to 20 characters
    x = new char[space];    // The text array
    clear ();               // Clear the text object
  }

  /** Clear the Text object */

  public void clear ()
  {
    size = 0;  // Start with an empty array

    // FILL THE ARRAY WITH SPACES
    for ( int i = 0; i < space; i = i + 1 )
    {
      x[ i ] = ' ';  // Set the ith character to a space
    }
  }

  /**
   * Copy the Text object
   *
   * @return a copy of the Text object
   */

  public Text copy ()
  {
    return new Text ( this );  // Return a copy
  }

  /**
   * Return the size of the Text object
   *
   * @return the size of the Text object
   */

  public int length ()
  {
    return size;  // Return the size
  }

/*--------------------------------------------------------------------
 isEmpty - Check if the array is empty
--------------------------------------------------------------------*/

  /**
   * Is it empty?
   *
   * @return true if the size is 0, false if size > 0
   */

  public boolean isEmpty ()
  {
    // IF THE ARRAY IS EMPTY
    if ( size == 0 )
    {
      return true;
    }

    // IF THE ARRAY IS NOT EMPTY
    else
    {
      return false;
    }
  }

  /**
   * Resize the text
   *
   * @param newspace new space
   * @param newk     new offset
   */

  public void resize ( int newspace, int newk )
  {
    // IF THE ARRAY IS TOO SMALL
    if ( newspace * ( 1.0 + m * 0.5 ) > space )
    {
      int oldspace = space;                          // Remember the size of the original array
      space = ( int ) ( newspace * ( 1.0 + m ) );    // Update to the new size
      char[] y = new char[space];                    // Create the new array
      int nk = ( space / 2 ) - ( size + newk / 2 );  // Find the new k
      nk = nk + newk;                                // Add the offset

      // COPY THE PREVIOUS ARRAY TO THE NEW ONE
      for ( int i = 0; i < size; i = i + 1 )
      {
        y[ nk + i ] = x[ k + i ];  // Set the ith character
      }

      // FILL THE REST OF THE ARRAY WITH SPACES
      for ( int i = 0; i < nk; i = i + 1 )
      {
        y[ i ] = ' ';  // Set the ith character to a space
      }

      // FILL THE REST OF THE ARRAY WITH SPACES
      for ( int i = nk + size; i < space; i = i + 1 )
      {
        y[ i ] = ' ';  // Set the ith character to a space
      }

      x = y;   // Save the new array
      k = nk;  // Save the new start
    }

    // IF THE TEXT IS GETTING TOO CLOSE TO THE LEFT
    if ( k < ( newspace - size ) )
    {
      int nk = ( space / 2 ) - ( size / 2 );  // Find the new k

      // MOVE THE ARRAY TO THE RIGHT
      for ( int i = size - 1; i >= 0; i = i - 1 )
      {
        x[ nk + i ] = x[ k + i ];  // Copy the ith character
      }

      // FILL THE REST OF THE ARRAY WITH SPACES
      for ( int i = nk + size; i < k + size; i = i + 1 )
      {
        x[ i ] = ' ';  // Set the ith character to a space
      }

      k = nk;  // Save the new start
    }

    // IF THE TEXT IS GETTING TOO CLOSE TO THE RIGHT
    if ( k > ( space - newspace ) )
    {
      int nk = ( space / 2 ) - ( size / 2 );  // Find the new k

      // MOVE THE ARRAY TO THE LEFT
      for ( int i = 0; i < size; i = i + 1 )
      {
        x[ nk + i ] = x[ k + i ];  // Copy the ith character
      }

      // FILL THE REST OF THE ARRAY WITH SPACES
      for ( int i = k; i < nk; i = i + 1 )
      {
        x[ i ] = ' ';  // Set the ith character to a space
      }

      k = nk;  // Save the new start
    }
  }

/*--------------------------------------------------------------------
 CONVERSION FUNCTIONS
--------------------------------------------------------------------*/

  /**
   * Convert to an integer
   *
   * @return an integer version of the object
   */

  public int num ()
  {
    return ( int ) ( Math.floor ( dbl () ) );
  }

  /**
   * Convert to a double
   *
   * @return a double version of the object
   */

  public double dbl ()
  {
    int s = 1;     // Track the sign
    double d = 0;  // Temporary number storage
    char c = ' ';  // Temporary character variable
    int sze = 0;   // Track the number of digits
    int dec = -1;  // Track the decimal
    int exp = 0;   // Extra exponent

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < size; i = i + 1 )
    {
      c = x[ k + i ];  // Get the ith character

      // IF AT THE FIRST CHARACTER AND IT IS A MINUS SIGN
      if ( ( i == 0 ) && ( c == '-' ) )
      {
        s = -1;  // Set the sign to negative
      }

      // OTHERWISE
      else
      {
        // IF THE CHARACTER IS A DIGIT
        if ( isNum ( c ) )
        {
          d = d * 10;         // Move the number to the right
          d = d + h2i ( c );  // Add the digit
          sze = sze + 1;      // Update the size
        }

        // IF THE CHARACTER IS THE DECIMAL POINT
        else if ( c == '.' )
        {
          dec = sze;          // Save the position
        }

        // IF THE CHARACTER IS THE DECIMAL POINT
        else if ( c == 'e' || c == 'E' )
        {
          exp = sub ( i + 1, size ).num ();  // Find the exponent number
          break;  // Break out of the loop
        }
      }
    }

    d = d * s;  // multiply by the sign

    // IF A DECIMAL WAS FOUND
    if ( dec >= 0 )
    {
      d = d / ( pow ( 10, sze - dec ) );  // Put in the decimal
    }

    d = d * pow ( 10, exp );  // Incorporate the exponent if included

    return d;  // Return the number
  }

  /**
   * Convert to lowercase
   *
   * @return a lowercase version of the object
   */

  public Text lower ()
  {
    // COPY THE TEXT
    Text a = copy ();

    // LOOP THROUGH THE TEXT
    for ( int i = 0; i < size; i = i + 1 )
    {
      a.x[ a.k + i ] = h2l ( x[ k + i ] );  // Get the ith character
    }

    return a;  // Return the lowercase Text
  }

  /**
   * Convert to uppercase
   *
   * @return an uppercase version of the object
   */

  public Text upper ()
  {
    // COPY THE TEXT
    Text a = copy ();

    // LOOP THROUGH THE TEXT
    for ( int i = 0; i < size; i = i + 1 )
    {
      a.x[ a.k + i ] = h2u ( x[ k + i ] );  // Get the ith character
    }

    return a;  // Return the uppercase Text
  }

  /**
   * Convert the first letter to uppercase
   *
   * @return a new Text with an uppercase first letter
   */

  public Text upperFirst ()
  {
    // COPY THE TEXT
    Text a = copy ();

    if ( a.size > 0 )
    {
      a.x[ a.k ] = h2u ( x[ k ] );  // Get the ith character
    }

    return a;  // Return the new Text
  }

  /**
   * Convert to only aphabet characters and whitespace
   *
   * @return a alphabet and whitespace only Text
   */

  public Text alpha ()
  {
    Text a = new Text ( "" );

    // LOOP THROUGH THE TEXT
    for ( int i = 0; i < size; i = i + 1 )
    {
      // IF THE CHARACTER IS FROM THE ALPHABET OR WHITESPACE
      if ( isAlpha ( x[ k + i ] ) || isWhite ( x[ k + i ] ) )
      {
        a.rpush ( x[ k + i ] );  // Add the character
      }

      // IF NON-ALPHABETIC AND NOT WHITESPACE
      else
      {
        a.rpush ( " " );  // Push a space
      }
    }

    return a.removeExtraWhite ();  // Return the new Text
  }

  /**
   * Convert to only numeric characters and whitespace
   *
   * @return a numeric and whitespace only Text
   */

  public Text numeric ()
  {
    Text a = new Text ( "" );

    // LOOP THROUGH THE TEXT
    for ( int i = 0; i < size; i = i + 1 )
    {
      if ( isNum ( x[ k + i ] ) || isWhite ( x[ k + i ] ) )
      {
        a.rpush ( x[ k + i ] );
      }

      else
      {
        a.rpush ( " " );
      }
    }

    return a.removeExtraWhite ();  // Return the new Text
  }

/*--------------------------------------------------------------------
 alphaNum - Convert to alpha-numeric with sign
--------------------------------------------------------------------*/

  public Text alphaNumSign ()
  {
    Text a = new Text ( "" );

    // LOOP THROUGH THE TEXT
    for ( int i = 0; i < size; i++ )
    {
      if ( isAlphaNumSign ( x[ k + i ] ) || isWhite ( x[ k + i ] ) )
      {
        a.rpush ( x[ k + i ] );
      }

      else
      {
        a.rpush ( " " );
      }
    }

    return a.removeExtraWhite ();  // Return the new Text
  }

/*--------------------------------------------------------------------
 alphaNum - Convert to alpha-numeric
--------------------------------------------------------------------*/

  public Text alphaNum ()
  {
    Text a = new Text ( "" );

    // LOOP THROUGH THE TEXT
    for ( int i = 0; i < size; i++ )
    {
      if ( isAlphaNum ( x[ k + i ] ) || isWhite ( x[ k + i ] ) )
      {
        a.rpush ( x[ k + i ] );
      }

      else
      {
        a.rpush ( " " );
      }
    }

    return a.removeExtraWhite ();  // Return the new Text
  }

/*--------------------------------------------------------------------
 alphaNum - Convert to alpha-numeric
--------------------------------------------------------------------*/

  public Text alphaNum ( String z )
  {
    return alphaNum ( s2t ( z ) );
  }

/*--------------------------------------------------------------------
 alphaNum - Convert to alpha-numeric
--------------------------------------------------------------------*/

  public Text alphaNum ( Text z )
  {
    Text a = new Text ( "" );

    // LOOP THROUGH THE TEXT
    for ( int i = 0; i < size; i = i + 1 )
    {
      if ( isAlphaNum ( x[ k + i ] ) || isWhite ( x[ k + i ] ) || z.find ( x[ k + i ] ) >= 0 )
      {
        a.rpush ( x[ k + i ] );
      }

      else
      {
        a.rpush ( " " );
      }
    }

    return a.removeExtraWhite ();  // Return the new Text
  }

/*--------------------------------------------------------------------
 EDIT FUNCTIONS
 ---------------------------------------------------------------------
 These functions edit the text
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 indent - Indent
--------------------------------------------------------------------*/

  public Text indent ()
  {
    return indent ( indentAmount );
  }

/*--------------------------------------------------------------------
 indent - Indent by s
--------------------------------------------------------------------*/

  public Text indent ( int s )
  {
    String y = "";

    for ( int i = 0; i < s; i++ )
    {
      y += " ";
    }

    Text z = join ( split ( "\n" ), "\n" + y );
    z.lpush ( y );

    y = "\n" + y;
    if ( z.sub ( z.size - y.length (), z.size ).eq ( y ) )
    {
      z = z.sub ( 0, z.size - y.length () );
    }

    return z;
  }

/*--------------------------------------------------------------------
 GET / SET FUNCTIONS
 ---------------------------------------------------------------------
 These functions manipulate the letters in the string
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 get - Get a character
--------------------------------------------------------------------*/

  public char get ( int i )
  {
    // IF i IS WITHIN RANGE
    if ( ( i >= 0 ) && ( i < size ) )
    {
      return x[ k + i ];  // Return the ith character
    }

    // IF i IS OUT OF RANGE
    else
    {
      return 0;  // Return the null character
    }
  }

/*--------------------------------------------------------------------
 toString - Get a string
--------------------------------------------------------------------*/

  public String toString ()
  {
    return str ();  // Return the text as a string
  }

/*--------------------------------------------------------------------
 get - Get a string
--------------------------------------------------------------------*/

  public String get ()
  {
    return str ();
  }

/*--------------------------------------------------------------------
 chr - Get a character array
--------------------------------------------------------------------*/

  public char[] chr ()
  {
    char[] y = new char[size];  // Create a new character array

    // LOOP THROUGH THE ARRAY
    for ( int i = k; i < k + size; i = i + 1 )
    {
      y[ i - k ] = x[ i ];  // Set the ith character
    }

    return y;  // Return the new array
  }

/*--------------------------------------------------------------------
 bytes - Get a byte array
--------------------------------------------------------------------*/

  public byte[] bytes ()
  {
    byte[] y = new byte[size];  // Create a new character array

    // LOOP THROUGH THE ARRAY
    for ( int i = k; i < k + size; i = i + 1 )
    {
      y[ i - k ] = ( byte ) ( x[ i ] );  // Set the ith character
    }

    return y;  // Return the new array
  }

/*--------------------------------------------------------------------
 str - Get a string
--------------------------------------------------------------------*/

  public String str ()
  {
    return c2s ( chr () );  // Return a string version of the character array
  }

/*--------------------------------------------------------------------
 set - Set a character
--------------------------------------------------------------------*/

  public void set ( int i, char z )
  {
    resize ( i, 0 );   // Ensure the array is a desired size
    x[ k + i ] = z;  // Set the ith character to z

    // IF i IS OUTSIDE THE CURRENT SIZE
    if ( i >= size )
    {
      size = i + 1;  // Set the size to fit i
    }
  }

/*--------------------------------------------------------------------
 set - Set to a character array
--------------------------------------------------------------------*/

  public void set ( char[] z )
  {
    clear ();     // Clear the Text object
    rpush ( z );  // Set to the supplied character array
  }

/*--------------------------------------------------------------------
 set - Set to a string
--------------------------------------------------------------------*/

  public void set ( String z )
  {
    clear ();     // Clear the Text object
    rpush ( z );  // Set to the supplied string
  }

/*--------------------------------------------------------------------
 set - Set to an integer
--------------------------------------------------------------------*/

  public void set ( int z )
  {
    clear ();     // Clear the Text object
    rpush ( z );  // Set to the supplied integer
  }

/*--------------------------------------------------------------------
 set - Set to a double
--------------------------------------------------------------------*/

  public void set ( double z )
  {
    clear ();     // Clear the Text object
    rpush ( z );  // Set to the supplied double
  }

/*--------------------------------------------------------------------
 set - Set to a Text object
--------------------------------------------------------------------*/

  public void set ( Text z )
  {
    clear ();     // Clear the Text object
    rpush ( z );  // Set to the supplied Text
  }

/*--------------------------------------------------------------------
 PUSH / POP FUNCTIONS
 ---------------------------------------------------------------------
 Add and remove letters from the object
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 lpop - Pop a character from the left side
--------------------------------------------------------------------*/

  public char lpop ()
  {
    // IF THE ARRAY IS NOT EMPTY
    if ( size > 0 )
    {
      char c = x[ k ];  // Save the first character
      x[ k ] = ' ';     // Set the character to a space
      k = k + 1;        // Increment the start index
      size = size - 1;  // Update the size of the array
      return c;         // Return the first character
    }

    // IF THE ARRAY IS EMPTY
    else
    {
      return 0;  // Return the null character
    }
  }

/*--------------------------------------------------------------------
 rpop - Pop a character from the right side
--------------------------------------------------------------------*/

  public char rpop ()
  {
    // IF THE ARRAY IS NOT EMPTY
    if ( size > 0 )
    {
      char c = x[ k + size - 1 ];  // Save the last character
      x[ k + size - 1 ] = ' ';     // Set the charatcter to a space
      size = size - 1;             // Update the size of the array
      return c;                    // Return the last character
    }

    // IF THE ARRAY IS EMPTY
    else
    {
      return 0;  // Return the null character
    }
  }

/*--------------------------------------------------------------------
 lpush - Add a character to the left side
--------------------------------------------------------------------*/

  public void lpush ( char z )
  {
    resize ( size + 1, 1 );   // Ensure the array is the desired size
    x[ k - 1 ] = z;         // Set the character
    k = k - 1;              // Update k
    size = size + 1;        // Store the new size
  }

/*--------------------------------------------------------------------
 lpush - Add a character array to the left side
--------------------------------------------------------------------*/

  public void lpush ( char[] z )
  {
    int length = z.length;        // Length of the character array
    int newsize = size + length;  // New size of the array
    resize ( newsize, length );     // Ensure the array is the desired size

    // ADD THE CHARACTER ARRAY TO THE LEFT SIDE OF THE ARRAY
    for ( int i = 0; i < length; i = i + 1 )
    {
      x[ k - length + i ] = z[ i ];  // Set the ith character
    }

    k = k - length;        // Update k
    size = size + length;  // Store the new size
  }

/*--------------------------------------------------------------------
 lpush - Add a double to the left side
--------------------------------------------------------------------*/

  public void lpush ( double z )
  {
    lpush ( d2s ( z ) );  // Push the double as a string
  }

/*--------------------------------------------------------------------
 lpush - Add an integer to the left side
--------------------------------------------------------------------*/

  public void lpush ( int z )
  {
    lpush ( i2s ( z ) );  // Push the integer as a string
  }

/*--------------------------------------------------------------------
 lpush - Add a string to the left side
--------------------------------------------------------------------*/

  public void lpush ( String z )
  {
    int length = strlen ( z );    // Length of the string
    int newsize = size + length;  // New size of the array
    resize ( newsize, length );     // Ensure the array is the desired size

    // ADD THE STRING TO THE LEFT SIDE OF THE ARRAY
    for ( int i = 0; i < length; i = i + 1 )
    {
      x[ k - length + i ] = chrat ( z, i );  // Set the ith character
    }

    k = k - length;        // Update k
    size = size + length;  // Store the new size
  }

/*--------------------------------------------------------------------
 lpush - Add text to the left side
--------------------------------------------------------------------*/

  public void lpush ( Text z )
  {
    int length = z.length ();     // Length of the text
    int newsize = size + length;  // New size of the array
    resize ( newsize, length );     // Ensure the array is the desired size

    // ADD TEXT TO THE LEFT SIDE OF THE ARRAY
    for ( int i = 0; i < length; i = i + 1 )
    {
      x[ k - length + i ] = z.get ( i );  // Set the ith character
    }

    k = k - length;        // Update k
    size = size + length;  // Store the new size
  }

/*--------------------------------------------------------------------
 nlpush - Copy and then add a character to the left side
--------------------------------------------------------------------*/

  public Text nlpush ( char z )
  {
    Text b = this.copy ();  // Copy the text
    b.lpush ( z );          // Push z onto the left side
    return b;               // Return the new Text object
  }

/*--------------------------------------------------------------------
 nlpush - Copy and then add a character array to the left side
--------------------------------------------------------------------*/

  public Text nlpush ( char[] z )
  {
    Text b = this.copy ();  // Copy the text
    b.lpush ( z );          // Push z onto the left side
    return b;               // Return the new Text object
  }

/*--------------------------------------------------------------------
 nlpush - Copy and then add a double to the left side
--------------------------------------------------------------------*/

  public Text nlpush ( double z )
  {
    Text b = this.copy ();  // Copy the text
    b.lpush ( z );          // Push z onto the left side
    return b;               // Return the new Text object
  }

/*--------------------------------------------------------------------
 nlpush - Copy and then add an integer to the left side
--------------------------------------------------------------------*/

  public Text nlpush ( int z )
  {
    Text b = this.copy ();  // Copy the text
    b.lpush ( z );          // Push z onto the left side
    return b;               // Return the new Text object
  }

/*--------------------------------------------------------------------
 nlpush - Copy and then add a string to the left side
--------------------------------------------------------------------*/

  public Text nlpush ( String z )
  {
    Text b = this.copy ();  // Copy the text
    b.lpush ( z );          // Push z onto the left side
    return b;               // Return the new Text object
  }

/*--------------------------------------------------------------------
 nlpush - Copy and then add text to the left side
--------------------------------------------------------------------*/

  public Text nlpush ( Text z )
  {
    Text b = this.copy ();  // Copy the text
    b.lpush ( z );          // Push z onto the left side
    return b;               // Return the new Text object
  }

/*--------------------------------------------------------------------
 rpush - Add a character to the right side
--------------------------------------------------------------------*/

  public void rpush ( char z )
  {
    resize ( size + 1, -1 );  // Ensure the array is the desired size
    x[ k + size ] = z;  // Set the character
    size = size + 1;        // Store the new size
  }

/*--------------------------------------------------------------------
 rpush - Add a character array to the right side
--------------------------------------------------------------------*/

  public void rpush ( char[] z )
  {
    int length = z.length;        // Length of the character array
    int newsize = size + length;  // New size of the array
    resize ( newsize, -length );    // Ensure the array is the desired size

    // ADD THE CHARACTER ARRAY TO THE RIGHT SIDE OF THE ARRAY
    for ( int i = 0; i < length; i = i + 1 )
    {
      x[ k + i + size ] = z[ i ];  // Set the ith character
    }

    size = size + length;  // Store the new size
  }

/*--------------------------------------------------------------------
 rpush - Add a double to the right side
--------------------------------------------------------------------*/

  public void rpush ( double z )
  {
    rpush ( d2s ( z ) );  // Push the double as a string
  }

/*--------------------------------------------------------------------
 rpush - Add an integer to the right side
--------------------------------------------------------------------*/

  public void rpush ( int z )
  {
    rpush ( i2s ( z ) );  // Push the integer as a string
  }

/*--------------------------------------------------------------------
 rpush - Add a string to the right side
--------------------------------------------------------------------*/

  public void rpush ( String z )
  {
    int length = strlen ( z );    // Length of the string
    int newsize = size + length;  // New size of the array
    resize ( newsize, -length );    // Ensure the array is the desired size

    // ADD THE STRING TO THE RIGHT SIDE OF THE ARRAY
    for ( int i = 0; i < length; i = i + 1 )
    {
      x[ k + i + size ] = chrat ( z, i );  // Set the ith character
    }

    size = size + length;  // Store the new size
  }

/*--------------------------------------------------------------------
 rpush - Add text to the right side
--------------------------------------------------------------------*/

  public void rpush ( Text z )
  {
    int length = z.length ();     // Length of the text
    int newsize = size + length;  // New size of the array
    resize ( newsize, -length );    // Ensure the array is the desired size

    // ADD TEXT TO THE RIGHT SIDE OF THE ARRAY
    for ( int i = 0; i < length; i = i + 1 )
    {
      x[ k + i + size ] = z.get ( i );  // Set the ith character
    }

    size = size + length;  // Store the new size
  }

/*--------------------------------------------------------------------
 nrpush - Copy and then add a character to the right side
--------------------------------------------------------------------*/

  public Text nrpush ( char z )
  {
    Text b = this.copy ();  // Copy the text
    b.rpush ( z );          // Push z onto the right side
    return b;               // Return the new Text object
  }

/*--------------------------------------------------------------------
 nrpush - Copy and then add a character array to the right side
--------------------------------------------------------------------*/

  public Text nrpush ( char[] z )
  {
    Text b = this.copy ();  // Copy the text
    b.rpush ( z );          // Push z onto the right side
    return b;               // Return the new Text object
  }

/*--------------------------------------------------------------------
 nrpush - Copy and then add a double to the right side
--------------------------------------------------------------------*/

  public Text nrpush ( double z )
  {
    Text b = this.copy ();  // Copy the text
    b.rpush ( z );          // Push z onto the right side
    return b;               // Return the new Text object
  }

/*--------------------------------------------------------------------
 nrpush - Copy and then add an integer to the right side
--------------------------------------------------------------------*/

  public Text nrpush ( int z )
  {
    Text b = this.copy ();  // Copy the text
    b.rpush ( z );          // Push z onto the right side
    return b;               // Return the new Text object
  }

/*--------------------------------------------------------------------
 nrpush - Copy and then add a string to the right side
--------------------------------------------------------------------*/

  public Text nrpush ( String z )
  {
    Text b = this.copy ();  // Copy the text
    b.rpush ( z );          // Push z onto the right side
    return b;               // Return the new Text object
  }

/*--------------------------------------------------------------------
 nrpush - Copy and then add text to the right side
--------------------------------------------------------------------*/

  public Text nrpush ( Text z )
  {
    Text b = this.copy ();  // Copy the text
    b.rpush ( z );          // Push z onto the right side
    return b;               // Return the new Text object
  }

/*--------------------------------------------------------------------
 SUBSTRING FUNCTIONS
 ---------------------------------------------------------------------
 Create substrings of the text
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 left - Search for a string, and grab text to the left
--------------------------------------------------------------------*/

  public Text left ( String z )
  {
    return left ( s2t ( z ) );
  }

/*--------------------------------------------------------------------
 left - Search for a string, and grab text to the left
--------------------------------------------------------------------*/

  public Text left ( Text z )
  {
    Text y = new Text ();  // Return text
    int i = find ( z );    // Find the first instance of the text

    // If an instance was found
    if ( i > 0 )
    {
      y = left ( i );  // Return the substring to the left of the Text
    }

    return y;  // Return the text
  }

/*--------------------------------------------------------------------
 onLeft - Check if a text is on the left
--------------------------------------------------------------------*/

  public boolean onLeft ( Text a )
  {
    return left ( a.size ).eq ( a );
  }

/*--------------------------------------------------------------------
 onRight - Check if a text is on the left
-------------------------------------------------------------------*/

  public boolean onRight ( Text a )
  {
    return right ( a.size ).eq ( a );
  }

/*--------------------------------------------------------------------
 left - Create a left substring
--------------------------------------------------------------------*/

  public Text left ( int a )
  {
    Text n;

    if ( a < 0 )
    {
      n = sub ( size + a, size );
    }

    else
    {
      n = sub ( 0, a );
    }

    return n;
  }

/*--------------------------------------------------------------------
 right - Create a right substring
--------------------------------------------------------------------*/

  public Text right ( int a )
  {
    Text n;

    if ( a < 0 )
    {
      n = sub ( -a, size );
    }

    else
    {
      n = sub ( size - a, size );
    }

    return n;
  }

/*--------------------------------------------------------------------
 sub - Select a substring of the array
--------------------------------------------------------------------*/

  public Text sub ( int c, int d )
  {
    d = lechk ( d, size );  // Ensure d <= size
    d = gechk ( d, 0 );     // Ensure d >= 0
    c = lechk ( c, d );     // Ensure c <= d
    c = gechk ( c, 0 );     // Ensure c >= 0
    int newsize = d - c;    // Size of the new array

    Text a = new Text ();
    a.resize ( newsize, 0 );

    // LOOP THROUGH THE ARRAY
    for ( int i = newsize - 1; i >= 0; i = i - 1 )
    {
      // SET THE ITH CHARACTER
      a.x[ a.k + i ] = x[ k + i + c ];
    }

    a.size = newsize;  // Size of the new array
    return a;          // Return the new array
  }

/*--------------------------------------------------------------------
 replace - Replace character a with b
--------------------------------------------------------------------*/

  public Text replace ( char a, char b )
  {
    return replace ( new Text ( a ), new Text ( b ) );
  }

/*--------------------------------------------------------------------
 replace - Replace character array a with b
--------------------------------------------------------------------*/

  public Text replace ( char[] a, char[] b )
  {
    return replace ( new Text ( a ), new Text ( b ) );
  }

/*--------------------------------------------------------------------
 replace - Replace double a with b
--------------------------------------------------------------------*/

  public Text replace ( double a, double b )
  {
    return replace ( new Text ( a ), new Text ( b ) );
  }

/*--------------------------------------------------------------------
 replace - Replace integer a with b
--------------------------------------------------------------------*/

  public Text replace ( int a, int b )
  {
    return replace ( new Text ( a ), new Text ( b ) );
  }

/*--------------------------------------------------------------------
 replace - Replace string a with b
--------------------------------------------------------------------*/

  public Text replace ( String a, String b )
  {
    return replace ( new Text ( a ), new Text ( b ) );
  }

/*--------------------------------------------------------------------
 replace - Replace a with b
--------------------------------------------------------------------*/

  public Text replace ( Text a, Text b )
  {
    int num = counts ( a );  // Number of instances of a
    int asize = a.size;     // Size of a
    int bsize = b.size;     // Size of b
    int offset = 0;         // Character offset

    // IF THERE ARE INSTANCES OF A
    if ( num > 0 )
    {
      int[] y = finds ( a );  // Find the instances of a
      Text c = new Text ();  // Create a new Text objext

      // LOOP THROUGH THE INSTANCES
      for ( int i = 0; i < num; i = i + 1 )
      {
        c.rpush ( sub ( offset, y[ i ] ) );  // Push the string to the left
        c.rpush ( b );                          // Push the replacement
        offset = y[ i ] + asize;                // Update the offset
      }

      c.rpush ( sub ( offset, size ) );  // Push the rest of the string
      return c;                             // Return the new Text object
    }

    // IF NO INSTANCES WERE FOUND
    else
    {
      return this.copy ();  // Return a copy of the Text object
    }
  }

/*--------------------------------------------------------------------
 split - Split with string a
--------------------------------------------------------------------*/

  public Text[] split ( String a )
  {
    return split ( new Text ( a ) );
  }

/*--------------------------------------------------------------------
 split - Split with a
--------------------------------------------------------------------*/

  public Text[] split ( Text a )
  {
    int num = counts ( a );        // Number of instances of a
    int asize = a.size;            // Size of a
    int offset = 0;                // Character offset
    Text[] c = new Text[num + 1];  // Result array

    // IF THERE ARE INSTANCES OF A
    if ( num > 0 )
    {
      int[] y = finds ( a );  // Find the instances

      // LOOP THROUGH THE INSTANCES
      for ( int i = 0; i < num; i = i + 1 )
      {
        c[ i ] = sub ( offset, y[ i ] );  // Save the text
        offset = ( y[ i ] ) + asize;      // Update the offset
      }

      c[ num ] = sub ( offset, size );  // Save the final string
    }

    // IF NO INSTANCES WERE FOUND
    else
    {
      c[ 0 ] = this.copy ();  // Copy the array
    }

    return c;  // Return the final array
  }

/*--------------------------------------------------------------------
 splitLine - Split with a new line
--------------------------------------------------------------------*/

  public Text[] splitLine ()
  {
    return split ( new Text ( "\n" ) );  // Split with a new line
  }

/*--------------------------------------------------------------------
 splitSpace - Split with a space
--------------------------------------------------------------------*/

  public Text[] splitSpace ()
  {
    return split ( new Text ( " " ) );  // Split with a space
  }

/*--------------------------------------------------------------------
 splitTab - Split with a tab
--------------------------------------------------------------------*/

  public Text[] splitTab ()
  {
    return split ( new Text ( "\t" ) );  // Split with a tab
  }

/*--------------------------------------------------------------------
 splitLineTab - Split with a tab on the lines
--------------------------------------------------------------------*/

  public Text[][] splitLineTab ()
  {
    return new Texts ( trimWhite ().splitLine () ).splitTab ();
  }

/*--------------------------------------------------------------------
 splitTabLine - Split with a tab on the lines
--------------------------------------------------------------------*/

  public Text[][] splitTabLine ()
  {
    int s = left ( "\n" ).splitTab ().length;  // Get the number of items per line
    Text b = replace ( "\n", "\t" );           // Switch new lines to tabs
    Text[] z = b.splitTab ();                  // Split on the tabs
    int zlen = z.length / s;                   // z length
    Text[][] y = new Text[s][zlen];            // Double text array
    int j = -1;                                // Line counter
    int k = 0;                                 // Item counter
    zlen = min ( zlen * s, z.length );          // Find the real z length

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < zlen; i++ )
    {
      k = i % s;  // Find the correct item

      // IF AT THE BEGINNING OF A LINE
      if ( k == 0 )
      {
        j++;  // Increment the line
      }

      y[ k ][ j ] = z[ i ];  // Set the item
    }

    return y;  // Return the Text arrays
  }

/*--------------------------------------------------------------------
 splitLineSpace - Split with a Space on the lines
--------------------------------------------------------------------*/

  public Text[][] splitLineSpace ()
  {
    return new Texts ( trimWhite ().splitLine () ).splitSpace ();
  }

/*--------------------------------------------------------------------
 splitSpaceLine - Split with a Space on the lines
--------------------------------------------------------------------*/

  public Text[][] splitSpaceLine ()
  {
    int s = left ( "\n" ).splitSpace ().length;  // Get the number of items per line
    Text b = replace ( "\n", " " );              // Switch new lines to tabs
    Text[] z = b.splitSpace ();                  // Split on the tabs
    int zlen = z.length / s;                     // z length
    Text[][] y = new Text[s][zlen];              // Double text array
    int j = -1;                                  // Line counter
    int k = 0;                                   // Item counter
    zlen = min ( zlen * s, z.length );           // Find the real z length

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < zlen; i++ )
    {
      k = i % s;  // Find the correct item

      // IF AT THE BEGINNING OF A LINE
      if ( k == 0 )
      {
        j++;  // Increment the line
      }

      y[ k ][ j ] = z[ i ];  // Set the item
    }

    return y;  // Return the Text arrays
  }

/*--------------------------------------------------------------------
 SEARCH FUNCTIONS
 ---------------------------------------------------------------------
 Search the text
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 contains - Does it contain z
--------------------------------------------------------------------*/

  public boolean contains ( String z )
  {
    return count ( z ) > 0;
  }

/*--------------------------------------------------------------------
 contains - Does it contain z
--------------------------------------------------------------------*/

  public boolean contains ( Text z )
  {
    return count ( z ) > 0;
  }

/*--------------------------------------------------------------------
 counts - Number of instances of z
--------------------------------------------------------------------*/

  public int counts ( String z )
  {
    return counts ( new Text ( z ), 0, size );
  }

/*--------------------------------------------------------------------
 counts - Number of instances of z starting with a
--------------------------------------------------------------------*/

  public int counts ( String z, int a )
  {
    return counts ( new Text ( z ), a, size );
  }

/*--------------------------------------------------------------------
 counts - Number of instances of z between a and b
--------------------------------------------------------------------*/

  public int counts ( String z, int a, int b )
  {
    return counts ( new Text ( z ), a, b );
  }

/*--------------------------------------------------------------------
 counts - Number of instances of z
--------------------------------------------------------------------*/

  public int counts ( Text z )
  {
    return counts ( z, 0, size );
  }

/*--------------------------------------------------------------------
 counts - Number of instances of z starting with a
--------------------------------------------------------------------*/

  public int counts ( Text z, int a )
  {
    return counts ( z, a, size );
  }

/*--------------------------------------------------------------------
 counts - Number of instances of z between a and b
--------------------------------------------------------------------*/

  public int counts ( Text z, int a, int b )
  {
    a = gechk ( a, 0 );     // Ensure a >= 0
    b = lechk ( b, size );  // Ensure b <= size
    int ct = 0;            // Counter
    int zsize = z.size;    // Size of z
    boolean flag = false;  // Instance flag

    // LOOP THROUGH THE ARRAY
    for ( int i = a; i < b - zsize + 1; i = i + 1 )
    {
      flag = true;  // Assume an instance

      // LOOP THROUGH Z
      for ( int j = 0; j < zsize; j = j + 1 )
      {
        // IF A CHARACTER DOES NOT MATCH
        if ( ( x[ k + i + j ] ) != ( z.x[ z.k + j ] ) )
        {
          flag = false;  // Set the flag to false
          break;         // Break out of the loop
        }
      }

      // IF A MATCH WAS FOUND
      if ( flag )
      {
        ct = ct + 1;  // Increment the counter
      }
    }

    return ct;  // Return the counter
  }

/*--------------------------------------------------------------------
 count - Number of instances of z
--------------------------------------------------------------------*/

  public int count ( String z )
  {
    return count ( new Text ( z ), 1, 0, size );
  }

/*--------------------------------------------------------------------
 count - Number of instances of z
--------------------------------------------------------------------*/

  public int count ( String z, int n )
  {
    return count ( new Text ( z ), n, 0, size );
  }

/*--------------------------------------------------------------------
 count - Number of instances of z starting with a
--------------------------------------------------------------------*/

  public int count ( String z, int n, int a )
  {
    return count ( new Text ( z ), n, a, size );
  }

/*--------------------------------------------------------------------
 count - Number of instances of z between a and b
--------------------------------------------------------------------*/

  public int count ( String z, int n, int a, int b )
  {
    return count ( new Text ( z ), n, a, b );
  }

/*--------------------------------------------------------------------
 count - Number of instances of z
--------------------------------------------------------------------*/

  public int count ( Text z )
  {
    return count ( z, 1, 0, size );
  }

/*--------------------------------------------------------------------
 count - Number of instances of z
--------------------------------------------------------------------*/

  public int count ( Text z, int n )
  {
    return count ( z, n, 0, size );
  }

/*--------------------------------------------------------------------
 count - Number of instances of z starting with a
--------------------------------------------------------------------*/

  public int count ( Text z, int n, int a )
  {
    return count ( z, n, a, size );
  }

/*--------------------------------------------------------------------
 count - Number of instances of z between a and b
--------------------------------------------------------------------*/

  public int count ( Text z, int n, int a, int b )
  {
    a = gechk ( a, 0 );     // Ensure a >= 0
    b = lechk ( b, size );  // Ensure b <= size
    int ct = 0;            // CountFirster
    int zsize = z.size;    // Size of z
    boolean flag = false;  // Instance flag

    // LOOP THROUGH THE ARRAY
    for ( int i = a; i < b - zsize + 1; i = i + 1 )
    {
      flag = true;  // Assume an instance

      // LOOP THROUGH Z
      for ( int j = 0; j < zsize; j = j + 1 )
      {
        // IF A CHARACTER DOES NOT MATCH
        if ( ( x[ k + i + j ] ) != ( z.x[ z.k + j ] ) )
        {
          flag = false;  // Set the flag to false
          break;         // Break out of the loop
        }
      }

      // IF A MATCH WAS FOUND
      if ( flag )
      {
        ct = ct + 1;  // Increment the counter

        // IF DESIRED NUMBER REACHED
        if ( ct == n )
        {
          break;  // Break out of the loop
        }
      }
    }

    return ct;  // Return the counter
  }

/*--------------------------------------------------------------------
 finds - Find the instances of z
--------------------------------------------------------------------*/

  public int[] finds ( String z )
  {
    return finds ( new Text ( z ), 0, size );
  }

/*--------------------------------------------------------------------
 finds - Find the instances of z starting at a
--------------------------------------------------------------------*/

  public int[] finds ( String z, int a )
  {
    return finds ( new Text ( z ), a, size );
  }

/*--------------------------------------------------------------------
 finds - Find the instances of z between a and b
--------------------------------------------------------------------*/

  public int[] finds ( String z, int a, int b )
  {
    return finds ( new Text ( z ), a, b );
  }

/*--------------------------------------------------------------------
 finds - Find the instances of z
--------------------------------------------------------------------*/

  public int[] finds ( Text z )
  {
    return finds ( z, 0, size );
  }

/*--------------------------------------------------------------------
 finds - Find the instances of z starting at a
--------------------------------------------------------------------*/

  public int[] finds ( Text z, int a )
  {
    return finds ( z, a, size );
  }

/*--------------------------------------------------------------------
 finds - Find the instances of z between a and b
--------------------------------------------------------------------*/

  public int[] finds ( Text z, int a, int b )
  {
    a = gechk ( a, 0 );                       // Ensure a >= 0
    b = lechk ( b, size );                    // Ensure b <= size
    int ct = 0;                              // Counter
    int[] y = new int[counts ( z, a, b )];  // Int array for the indices
    int zsize = z.size;                      // Size of z
    boolean flag = false;                    // Instance flag

    // LOOP THROUGH THE ARRAY
    for ( int i = a; i < b - zsize + 1; i = i + 1 )
    {
      flag = true;  // Assume an instance

      // LOOP THROUGH Z
      for ( int j = 0; j < zsize; j = j + 1 )
      {
        // IF A CHARACTER DOES NOT MATCH
        if ( ( x[ k + i + j ] ) != ( z.x[ z.k + j ] ) )
        {
          flag = false;  // Set the flag to false
          break;         // Break out of the loop
        }
      }

      // IF A MATCH WAS FOUND
      if ( flag )
      {
        y[ ct ] = i;  // Store the index
        ct = ct + 1;  // Increment the counter
      }
    }

    return y;  // Return the array
  }

/*--------------------------------------------------------------------
 find - Find the first n instance of z
--------------------------------------------------------------------*/

  public int find ( char z )
  {
    int[] y = find ( new Text ( z + "" ), 1, 0, size );  // Find the first instance

    // IF A RESULT WAS FOUND
    if ( y.length > 0 )
    {
      return y[ 0 ];  // Return the instance
    }

    else
    {
      return -1;  // Return a negative one
    }
  }

/*--------------------------------------------------------------------
 find - Find the first n instance of z
--------------------------------------------------------------------*/

  public int find ( String z )
  {
    int[] y = find ( new Text ( z ), 1, 0, size );  // Find the first instance

    // IF A RESULT WAS FOUND
    if ( y.length > 0 )
    {
      return y[ 0 ];  // Return the instance
    }

    else
    {
      return -1;  // Return a negative one
    }
  }

/*--------------------------------------------------------------------
 find - Find the first n instances of z
--------------------------------------------------------------------*/

  public int[] find ( String z, int n )
  {
    return find ( new Text ( z ), n, 0, size );
  }

/*--------------------------------------------------------------------
 find - Find the first n instances of z starting at a
--------------------------------------------------------------------*/

  public int[] find ( String z, int n, int a )
  {
    return find ( new Text ( z ), n, a, size );
  }

/*--------------------------------------------------------------------
 find - Find the first n instances of z between a and b
--------------------------------------------------------------------*/

  public int[] find ( String z, int n, int a, int b )
  {
    return find ( new Text ( z ), n, a, b );
  }

/*--------------------------------------------------------------------
 find - Find the first n instance of z
--------------------------------------------------------------------*/

  public int find ( Text z )
  {
    int[] y = find ( z, 1, 0, size );  // Find the first instance

    // IF A RESULT WAS FOUND
    if ( y.length > 0 )
    {
      return y[ 0 ];  // Return the instance
    }

    else
    {
      return -1;  // Return a negative one
    }
  }

/*--------------------------------------------------------------------
 find - Find the first n instances of z
--------------------------------------------------------------------*/

  public int[] find ( Text z, int n )
  {
    return find ( z, n, 0, size );
  }

/*--------------------------------------------------------------------
 find - Find the first n instances of z starting at a
--------------------------------------------------------------------*/

  public int[] find ( Text z, int n, int a )
  {
    return find ( z, n, a, size );
  }

/*--------------------------------------------------------------------
 find - Find the first n instances of z between a and b
--------------------------------------------------------------------*/

  public int[] find ( Text z, int n, int a, int b )
  {
    a = gechk ( a, 0 );                       // Ensure a >= 0
    b = lechk ( b, size );                    // Ensure b <= size
    int ct = 0;                              // Counter
    int[] y = new int[min ( n, count ( z, n, a, b ) )];  // Int array for the indices
    int zsize = z.size;                      // Size of z
    boolean flag = false;                    // Instance flag

    // LOOP THROUGH THE ARRAY
    for ( int i = a; i < b - zsize + 1; i = i + 1 )
    {
      flag = true;  // Assume an instance

      // LOOP THROUGH Z
      for ( int j = 0; j < zsize; j = j + 1 )
      {
        // IF A CHARACTER DOES NOT MATCH
        if ( ( x[ k + i + j ] ) != ( z.x[ z.k + j ] ) )
        {
          flag = false;  // Set the flag to false
          break;         // Break out of the loop
        }
      }

      // IF A MATCH WAS FOUND
      if ( flag )
      {
        y[ ct ] = i;  // Store the index
        ct = ct + 1;  // Increment the counter

        // IF ALREADY FOUND THE REQUIRED NUMBER
        if ( ct == y.length )
        {
          break; // Break out of the loop
        }
      }
    }

    return y;  // Return the array
  }

/*--------------------------------------------------------------------
 findPosts - Find the instances of z
--------------------------------------------------------------------*/

  public int[] findPosts ( String z )
  {
    return findPosts ( new Text ( z ), 0, size );
  }

/*--------------------------------------------------------------------
 findPosts - Find the instances of z starting at a
--------------------------------------------------------------------*/

  public int[] findPosts ( String z, int a )
  {
    return findPosts ( new Text ( z ), a, size );
  }

/*--------------------------------------------------------------------
 findPosts - Find the instances of z between a and b
--------------------------------------------------------------------*/

  public int[] findPosts ( String z, int a, int b )
  {
    return findPosts ( new Text ( z ), a, b );
  }

/*--------------------------------------------------------------------
 findPosts - Find the instances of z
--------------------------------------------------------------------*/

  public int[] findPosts ( Text z )
  {
    return findPosts ( z, 0, size );
  }

/*--------------------------------------------------------------------
 findPosts - Find the instances of z starting at a
--------------------------------------------------------------------*/

  public int[] findPosts ( Text z, int a )
  {
    return findPosts ( z, a, size );
  }

/*--------------------------------------------------------------------
 findPosts - Find the instances of z between a and b
--------------------------------------------------------------------*/

  public int[] findPosts ( Text z, int a, int b )
  {
    int[] y = finds ( z, a, b ); // Find the instances of z
    y = add ( y, z.length () );  // Add the length of z
    return y;                    // Return the results
  }

/*--------------------------------------------------------------------
 findPost - Find the first n instance of z
--------------------------------------------------------------------*/

  public int findPost ( String z )
  {
    int[] y = findPost ( new Text ( z ), 1, 0, size );  // Find the first instance

    // IF A RESULT WAS FOUND
    if ( y.length > 0 )
    {
      return y[ 0 ];  // Return the instance
    }

    else
    {
      return -1;  // Return a negative one
    }
  }

/*--------------------------------------------------------------------
 findPost - Find the first n instances of z
--------------------------------------------------------------------*/

  public int[] findPost ( String z, int n )
  {
    return findPost ( new Text ( z ), n, 0, size );
  }

/*--------------------------------------------------------------------
 findPost - Find the first n instances of z starting at a
--------------------------------------------------------------------*/

  public int[] findPost ( String z, int n, int a )
  {
    return findPost ( new Text ( z ), n, a, size );
  }

/*--------------------------------------------------------------------
 findPost - Find the first n instances of z between a and b
--------------------------------------------------------------------*/

  public int[] findPost ( String z, int n, int a, int b )
  {
    return findPost ( new Text ( z ), n, a, b );
  }

/*--------------------------------------------------------------------
 findPost - Find the first n instance of z
--------------------------------------------------------------------*/

  public int findPost ( Text z )
  {
    int[] y = findPost ( z, 1, 0, size );  // Find the first instance

    // IF A RESULT WAS FOUND
    if ( y.length > 0 )
    {
      return y[ 0 ];  // Return the instance
    }

    else
    {
      return -1;  // Return a negative one
    }
  }

/*--------------------------------------------------------------------
 findPost - Find the first n instances of z
--------------------------------------------------------------------*/

  public int[] findPost ( Text z, int n )
  {
    return findPost ( z, n, 0, size );
  }

/*--------------------------------------------------------------------
 findPost - Find the first n instances of z starting at a
--------------------------------------------------------------------*/

  public int[] findPost ( Text z, int n, int a )
  {
    return findPost ( z, n, a, size );
  }

/*--------------------------------------------------------------------
 findPost - Find the first n instances of z between a and b
--------------------------------------------------------------------*/

  public int[] findPost ( Text z, int n, int a, int b )
  {
    int[] y = find ( z, n, a, b );  // Find the instances of z
    y = add ( y, z.length () );          // Add the length of z
    return y;                            // Return the results
  }

/*--------------------------------------------------------------------
 findBetween - Find a string between two strings
--------------------------------------------------------------------*/

  public Text findBetween ( String y, String z )
  {
    return findBetween ( y, z, 0 );  // Call the finds function with offset 0
  }

/*--------------------------------------------------------------------
 findBetween - Find a string between y and z
--------------------------------------------------------------------*/

  public Text findBetween ( String y, String z, int j )
  {
    Text a = findBetween ( new Text ( y ), new Text ( z ), j );  // Call the finds function with offset j

    // IF A RESULT WAS RETURNED
    if ( a != null )
    {
      return a;  // Return a as a string
    }

    // IF NO RESULT WAS FOUND
    else
    {
      return null;
    }
  }

/*--------------------------------------------------------------------
 findBetween - Find a string between two strings
--------------------------------------------------------------------*/

  public Text findBetween ( Text y, Text z )
  {
    return findBetween ( y, z, 0 );  // Call the finds function with offset 0
  }

/*--------------------------------------------------------------------
 findBetween - Find a Text between y and z
--------------------------------------------------------------------*/

  public Text findBetween ( Text y, Text z, int j )
  {
    Text t = new Text ();                // A Text object to hold the result
    int a[] = findPost ( y, 1, j );  // Find the first string

    // IF THE STRING WAS FOUND
    if ( a.length > 0 && a[ 0 ] >= 0 )
    {
      int b[] = find ( z, 1, a[ 0 ] );  // Find the second string

      // IF THE STRING WAS FOUND
      if ( b.length > 0 && b[ 0 ] >= 0 )
      {
        t = sub ( a[ 0 ], b[ 0 ] );  // Get the string
      }
    }

    // RETURN THE TEXT
    return t;
  }

/*--------------------------------------------------------------------
 findBetweens - Find a string between two strings
--------------------------------------------------------------------*/

  public Text[] findBetweens ( String y, String z )
  {
    return findBetweens ( y, z, 0 );  // Call the finds function with offset 0
  }

/*--------------------------------------------------------------------
 findBetweens - Find a string between y and z
--------------------------------------------------------------------*/

  public Text[] findBetweens ( String y, String z, int j )
  {
    return findBetweens ( new Text ( y ), new Text ( z ), j );
  }

/*--------------------------------------------------------------------
 findBetweens - Find a string between two strings
--------------------------------------------------------------------*/

  public Text[] findBetweens ( Text y, Text z )
  {
    return findBetweens ( y, z, 0 );  // Call the finds function with offset 0
  }

/*--------------------------------------------------------------------
 findBetweens - Find a Text between y and z
--------------------------------------------------------------------*/

  public Text[] findBetweens ( Text y, Text z, int j )
  {
    int js = j;                           // Save j
    int a[] = findPost ( y, 1, j );       // Find the first string
    int n = 0;                            // result counter

    // IF THE STRING WAS FOUND
    while ( a.length > 0 && a[ 0 ] >= 0 )
    {
      int b[] = find ( z, 1, a[ 0 ] );  // Find the second string

      // IF THE STRING WAS FOUND
      if ( b.length > 0 && b[ 0 ] >= 0 )
      {
        n++;                            // Increment the result counter
        j = b[ 0 ] + b.length;      // Increment the start index
        a = findPost ( y, 1, j );       // Find the next string
      }

      // IF THE STRING WAS NOT FOUND
      else
      {
        break;  // Break out of the loop
      }
    }

    j = js;                         // Reset index
    Text t[] = new Text[n];       // Create a Text array for the results
    int i = 0;                      // Result counter
    a = findPost ( y, 1, j );       // Find the first string

    // IF THE STRING WAS FOUND
    while ( a.length > 0 && a[ 0 ] >= 0 )
    {
      int b[] = find ( z, 1, a[ 0 ] );  // Find the second string

      // IF THE STRING WAS FOUND
      if ( b.length > 0 && b[ 0 ] >= 0 )
      {
        t[ i ] = sub ( a[ 0 ], b[ 0 ] );  // Get the string
        i++;                                 // Increment the counter
        j = b[ 0 ] + b.length;           // Increment the start index
        a = findPost ( y, 1, j );            // Find the next string
      }

      // IF THE STRING WAS NOT FOUND
      else
      {
        break;  // Break out of the loop
      }
    }

    return t;  // Return the text array
  }

/*--------------------------------------------------------------------
 tag - Get the XML tag
--------------------------------------------------------------------*/

  public Text tag ( String z )
  {
    return tag ( z, 0 );  // Return the tag
  }

/*--------------------------------------------------------------------
 tag - Get the XML tag starting at i
--------------------------------------------------------------------*/

  public Text tag ( String z, int i )
  {
    return findBetween ( "<" + z + ">", "</" + z + ">", i );  // Return the tag
  }

/*--------------------------------------------------------------------
 tag - Get the XML tag
--------------------------------------------------------------------*/

  public Text[] tags ( String z )
  {
    return tags ( z, 0 );  // Return the tag
  }

/*--------------------------------------------------------------------
 tag - Get the XML tag starting at i
--------------------------------------------------------------------*/

  public Text[] tags ( String z, int i )
  {
    return findBetweens ( "<" + z + ">", "</" + z + ">", i );  // Return the tag
  }

/*--------------------------------------------------------------------
 remove - Remove the string z
--------------------------------------------------------------------*/

  public Text remove ( String z )
  {
    return replace ( z, "" );  // Return the new removed text
  }

/*--------------------------------------------------------------------
 remove - Remove the Text z
--------------------------------------------------------------------*/

  public Text remove ( Text z )
  {
    return replace ( z, new Text ( "" ) );  // Return the new removed text
  }

/*--------------------------------------------------------------------
 remove - Remove any text with the start and end text
--------------------------------------------------------------------*/

  public Text remove ( String[] z )
  {
    return remove ( s2t ( z ) );  // Return the new removed text
  }

/*--------------------------------------------------------------------
 remove - Remove any text with the start and end text
--------------------------------------------------------------------*/

  public Text remove ( Text[] z )
  {
    Text t, a, b, u;        // Temporary Text objects
    int ai[], bi[];         // Index arrays
    int an, bn, pr;         // Index holders
    int asize;              // Text size
    t = new Text ( this );  // Create a new Text object

    // LOOP THROUGH THE PAIRS OF TEXTS
    for ( int i = 0; i < z.length; i += 2 )
    {
      u = new Text ( t );  // Copy the Text
      t = new Text ();     // Create a new Text object

      a = z[ i ];           // Get the first Text a
      b = z[ i + 1 ];       // Get the second Text b
      asize = a.length ();  // Size of a

      ai = u.finds ( a );      // Find the instances of a
      bi = u.findPosts ( b );  // Find the instances of b

      an = next ( ai, 0 );  // Find the next item
      pr = 0;               // Previous index

      // LOOP UNTIL NO OTHER ITEMS FOUND
      while ( an >= 0 )
      {
        // IF THE NEXT ITEM WAS FOUND
        if ( an >= 0 )
        {
          bn = next ( bi, an + asize );  // Find the next b

          // IF THE NEXT ITEM WAS FOUND
          if ( bn >= 0 )
          {
            t.rpush ( u.sub ( pr, an ) );  // Push the non-removed text
            pr = bn;                          // Update the previous index
            an = next ( ai, bn );             // Find the next a index
          }

          // IF THE NEXT ITEM WAS NOT FOUND
          else
          {
            break;  // Break from the loop
          }
        }

        // IF THE NEXT ITEM WAS NOT FOUND
        else
        {
          break;  // Break from the loop
        }
      }

      t.rpush ( u.sub ( pr, u.length () ) );   // Push the non-removed text
    }

    return t;  // Return the Text object
  }

/*--------------------------------------------------------------------
 removeTag - Remove any html tags
--------------------------------------------------------------------*/

  public Text removeTag ()
  {
    Text t = new Text ();       // Create a new Text object
    boolean whiteFlag = false;  // Whitespace flag

    // LOOP THROUGH THE CHARACTER
    for ( int i = 0; i < size; i++ )
    {
      // IF THE CURRENT CHARACTER IS A TAG OPENER
      if ( x[ k + i ] == '<' )
      {
        // LOOP UNTIL THE TAG CLOSER
        while ( i < size && x[ k + i ] != '>' )
        {
          i++;  // Increment i
        }

        i++;  // Increment past the tag closer
      }

      // IF NOT AT THE END
      if ( i < size )
      {
        t.rpush ( x[ k + i ] );  // Push the character
      }
    }

    return t.trimWhite ();  // Return the Text object
  }

/*--------------------------------------------------------------------
 removeExtraWhite - Remove any multiple whitespace
--------------------------------------------------------------------*/

  public Text removeExtraWhite ()
  {
    Text t = new Text ();  // Create a new Text object
    boolean whiteFlag = false;   // Whitespace flag

    // LOOP THROUGH THE CHARACTER
    for ( int i = 0; i < size; i++ )
    {
      // IF THE CURRENT CHARACTER IS WHITESPACE
      if ( isWhite ( x[ k + i ] ) )
      {
        // IF THE LAST CHARACTER WAS NOT WHITESPACE
        if ( !whiteFlag )
        {
          whiteFlag = true;  // Set the whitespace flag
          t.rpush ( " " );    // Push a space
        }
      }

      // IF THE CURRENT CHARACTER IS NOT WHITESPACE
      else
      {
        t.rpush ( x[ k + i ] );  // Push the character
        whiteFlag = false;       // Unset the whitespace flag
      }
    }

    return t.trimWhite ();  // Return the Text object
  }

/*--------------------------------------------------------------------
 trimWhite - Remove any outside whitespace
--------------------------------------------------------------------*/

  public Text trimWhite ()
  {
    Text t = new Text ();  // Create a new Text object
    int i, j;              // Indices

    // LOOP THROUGH THE CHARACTERS
    for ( i = 0; i < size; i++ )
    {
      // IF THE CURRENT CHARACTER IS NOT WHITESPACE
      if ( !isWhite ( x[ k + i ] ) )
      {
        break;  // Break out of the loop
      }
    }

    // LOOP THROUGH THE CHARACTERS
    for ( j = size - 1; j >= 0; j-- )
    {
      // IF THE CURRENT CHARACTER IS NOT WHITESPACE
      if ( !isWhite ( x[ k + j ] ) )
      {
        j++;    // Need to choose the next character
        break;  // Break out of the loop
      }
    }

    return sub ( i, j );  // Return the Text object
  }

/*
// FIND THE URLS IN A TEXT

public Text[][] findUrlTitles()
{
  // FIND ALL HREFS
  Text[] tu = findBetweens( "href=\"", "</a>" );

  // CREATE THE URL ARRAY
  Text[][] u = new Text[ tu.length ][ 2 ];

  // LOOP THROUGH THE URLS
  for ( int i = 0; i < tu.length; i++ )
  {
    // FIND THE CLOSING QUOTE
    int j = tu[ i ].finds( "\"" );

    // FIND THE CLOSING BRACKET
    int k = tu[ i ].finds( ">" );

    // IF THE CLOSING QUOTE WAS FOUND
    if ( j >= 0 )
    {
      // FIND THE URL
      u[ i ][ 0 ] = tu[ i ].sub( 0, j );
    }

    // IF THE CLOSING QUOTE WAS NOT FOUND
    else
    {
      // SET IT TO THE EMPTY STRING
      u[ i ][ 0 ] = new Text( "" );
    }

    // IF THE CLOSING BRACKET WAS FOUND
    if ( k >= 0 )
    {
      // FIND THE URL
      u[ i ][ 1 ] = tu[ i ].sub( k + 1, tu[ i ].length() );
    }

    // IF THE CLOSING BRACKET WAS NOT FOUND
    else
    {
      // SET IT TO THE EMPTY STRING
      u[ i ][ 1 ] = new Text( "" );
    }
  }

  // RETURN THE HREFS
  return u;
}

// FIND THE URLS IN A TEXT

public Text[] findUrls()
{
  // FIND ALL HREFS
  Text[] u = findBetweens( "href=\"", "\"" );

  // RETURN THE HREFS
  return u;
}

// FIX A LINK

public Text fixlink()
{
  // CREATE A COPY OF THE TEXT
  Text t = new Text( this );

  // RETURN THE TEXT
  return t;
}

// UNFIX A LINK

public Text unfixlink()
{
  // CREATE A COPY OF THE TEXT
  Text t = new Text( this );

  // FIX THE ESCAPED CHARACTERS
  t = t.replace( "%3F", "?" );
  t = t.replace( "%3D", "=" );
  t = t.replace( "%26", "&" );

  // RETURN THE TEXT
  return t;
}

// REMOVE A TAG

public Text removeTag( String z )
{
  return replace( "<" + z + ">", "" ).replace( "</" + z + ">", "" );
}

// GET THE XML TAG

public Text[] tags( String z )
{
  // RETURN THE TAG
  return tags( z, 0 );
}

// GET THE XML TAG

public Text[] tags( String z, int i )
{
  // RETURN THE TAG
  return findBetweens( "<" + z + ">", "</" + z + ">", i );
}

// STRIP A TEXT OF NON-ALPHANUMERIC CHARACTERS

public Text strip()
{
  // CREATE A TEXT OBJECT
  Text n = new Text( this );

  // REMOVE SOME HTML TAGS
  n = n.replace( "&hellip;", "" );
  n = n.replace( "&nbsp;", "" );

  // CREATE A TEMPORARY CHARACTER
  char c = ' ';

  // CREATE A WHITESPACE FLAG
  boolean spaceFlag = false;

  // CREATE A TAG FLAG
  boolean tagFlag = false;

  // CREATE AN OUTPUT STRING
  char[] o = new char[ length() ];

  // CREATE AN INDEX
  int j = 0;

  // LOOP THROUGH THE CHARACTERS
  for ( int i = 0; i < n.length(); i++ )
  {
    // RETRIEVE THE CHARACTER
    c = n.chr( i );

    // IF IN A TAG
    if ( tagFlag )
    {
      // IF AT THE FINAL CHARACTER IN A TAG
      if ( c == '>' )
      {
        // UNSET THE TAG FLAG
        tagFlag = false;

        // IF PREVIOUS CHARACTER WAS NOT WHITESPACE
        if ( ! spaceFlag )
        {
          // ADD THE CHARACTER
          o[ j ] = ' ';

          // INCREMENT THE INDEX
          j++;

          // SET THE WHITESPACE FLAG
          spaceFlag = true;
        }
      }
    }

    // IF AT THE START OF A TAG
    else if ( c == '<' )
    {
      // SET THE TAG FLAG
      tagFlag = true;
    }

    // IF THE LETTER IS ALPHANUMERIC OR PUNCTUATION
    else if ( isAlphaNum( c ) || ( c == '.' ) || ( c == ',' ) || ( c == '-' ) )
    {
      // ADD THE CHARACTER
      o[ j ] = c;

      // INCREMENT THE INDEX
      j++;

      // UNSET THE WHITESPACE FLAG
      spaceFlag = false;
    }

    // IF THE CHARACTER IS NOT ALPHANUMERIC
    else
    {
      // IF PREVIOUS CHARACTER WAS NOT WHITESPACE
      if ( ! spaceFlag )
      {
        // ADD THE CHARACTER
        o[ j ] = ' ';

        // INCREMENT THE INDEX
        j++;

        // SET THE WHITESPACE FLAG
        spaceFlag = true;
      }
    }
  }

  // CREATE A TEXT OBJECT
  Text t = new Text( new String( o, 0, j ) );

  // RETURN THE STRIPPED TEXT
  return t.trim();
}

// TRIM A TEXT OF WHITESPACE

public Text trim()
{
  // CREATE A TEMPORARY CHARACTER
  char c = ' ';

  // CREATE AN INDEX FOR THE START OF THE STRING
  int i = 0;

  // CREATE AN INDEX FOR THE END OF THE STRING
  int j = 0;

  // LOOP THROUGH THE CHARACTERS
  for ( i = 0; i < length(); i++ )
  {
    // RETRIEVE THE CHARACTER
    c = chr( i );

    // IF A LETTER, IS NOT WHITESPACE
    if ( ! isWhite( c ) )
    {
      // BREAK OUT OF THE LOOP
      break;
    }
  }

  // LOOP THROUGH THE CHARACTERS
  for ( j = length() - 1; j >= 0; j-- )
  {
    // RETRIEVE THE CHARACTER
    c = chr( j );

    // IF A LETTER, IS NOT WHITESPACE
    if ( ! isWhite( c ) )
    {
      // BREAK OUT OF THE LOOP
      break;
    }
  }

  // RETURN THE TRIMMED TEXT
  return sub( i, j + 1 );
}

// CONVERT TO DIGITS

public Text digit()
{
  // CREATE A TEXT OBJECT
  Text t = new Text( "" );

  // CREATE A TEMPORARY CHARACTER
  char c = ' ';

  // LOOP THROUGH THE CHARACTERS
  for ( int i = 0; i < size; i++ )
  {
    // RETRIEVE THE CHARACTER
    c = get( i );

    // IF AN UPPERCASE LETTER
    if ( isNum( c ) || c == '-' || c == '.' )
    {
      // ADD THE CHARACTER
      t.rpush( c );
    }
  }

  // RETURN THE DIGIT TEXT
  return t;
}
*/

  public Texts split ()
  {
    Text z = trimWhite ();
    Text[] y = z.splitLine ();  // Split into lines
    Texts x = new Texts ();     // Holds tags sequentially

    // LOOP THROUGH THE LINES
    for ( int i = 0; i < y.length; i++ )
    {
      // PUSH THE TAGS AND A NEW LINE
      x.rpush ( y[ i ].splitSpace () );
    }

    return x;
  }
}
