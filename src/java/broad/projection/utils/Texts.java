package broad.projection.utils;

import java.util.ArrayList;



/*--------------------------------------------------------------------
 com.scanfeld.core.Texts
 ---------------------------------------------------------------------
 Models a list of Texts
 ---------------------------------------------------------------------
 Author: Daniel Scanfeld
   Date: 08.31.06
 ---------------------------------------------------------------------
 Optimize the memory stuff mr, ml, ...
--------------------------------------------------------------------*/

public class Texts extends Inc implements Compare
{
/*--------------------------------------------------------------------
 Class Variables
--------------------------------------------------------------------*/

  public Text[] x;   // Text array
  public Tree t;     // Lookup tree for the list
  public int k;      // The index of the first item
  public int size;   // The size of the text
  public int space;  // The size of the text array
  public double ml;  // Left side multiply amount
  public double mr;  // Right side multiply amount

/*--------------------------------------------------------------------
 CONSTRUCTOR FUNCTIONS
 ---------------------------------------------------------------------
 Constructors for the possible text objects
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 Texts - Null constructor
--------------------------------------------------------------------*/

  public Texts ()
  {
    init ();
  }
/*--------------------------------------------------------------------
 Texts - size constructor
--------------------------------------------------------------------*/

  public Texts ( int sze )
  {
    init ();
    size ( sze, 0 );
  }

/*--------------------------------------------------------------------
 Texts - Texts constructor
--------------------------------------------------------------------*/

  public Texts ( Texts z )
  {
    init ();

    size ( z.size, -z.size / 2 );  // Size the array for the new samples

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < z.size; i++ )
    {
      rpush ( z.get ( i ) );  // Push the ith item
    }
  }

/*--------------------------------------------------------------------
 Texts - Text array constructor
--------------------------------------------------------------------*/

  public Texts ( Text[] z )
  {
    init ();

    size ( z.length, -z.length / 2 );  // Size the array for the new samples

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < z.length; i++ )
    {
      rpush ( z[ i ] );  // Push the ith item
    }
  }

/*--------------------------------------------------------------------
 Texts - Text array constructor
--------------------------------------------------------------------*/

  public Texts ( Text[] z, int a, int b )
  {
    init ();

    if ( a < 0 )
    {
      a = 0;
    }

    if ( b >= z.length )
    {
      b = z.length - 1;
    }

    int len = b - a;
    size ( len, -len / 2 );  // Size the array for the new samples

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < len; i++ )
    {
      rpush ( z[ a + i ] );  // Push the ith item
    }
  }

/*--------------------------------------------------------------------
 Texts - String array constructor
--------------------------------------------------------------------*/

  public Texts ( String[] z )
  {
    init ();

    size ( z.length, -z.length / 2 );  // Size the array for the new samples

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < z.length; i++ )
    {
      rpush ( z[ i ] );  // Push the ith item
    }
  }

  
  public Texts ( ArrayList<String> z )
  {
    init ();

    size ( z.size(), -z.size() / 2 );  // Size the array for the new samples

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < z.size(); i++ )
    {
      rpush ( z.get(i)  );  // Push the ith item
    }
  }
  
/*--------------------------------------------------------------------
 Texts - Int array constructor
--------------------------------------------------------------------*/

  public Texts ( int[] z )
  {
    init ();

    size ( z.length, -z.length / 2 );  // Size the array for the new samples

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < z.length; i++ )
    {
      rpush ( z[ i ] );  // Push the ith item
    }
  }

/*--------------------------------------------------------------------
 HELPER FUNCTIONS
 ---------------------------------------------------------------------
 These are the main helper functions for manipulating the object
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 contains - does it contain a
--------------------------------------------------------------------*/

  public boolean[] contains ( String a )
  {
    return contains ( new Text ( a ) );
  }

/*--------------------------------------------------------------------
 contains - does it contain a
--------------------------------------------------------------------*/

  public boolean[] contains ( Text a )
  {
    boolean b[] = new boolean[size];  // Boolean flag array

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < size; i++ )
    {
      // IF THE TEXT IS EQUAL
      if ( get ( i ).contains ( a ) )
      {
        b[ i ] = true;  // Set the flag to true
      }

      // IF THE TEXT HAS NOT BEEN ENTERED
      else
      {
        b[ i ] = false;  // Set the flag to false
      }
    }

    return b;
  }

/*--------------------------------------------------------------------
 eq - is a equal to this Texts object
--------------------------------------------------------------------*/

  public boolean eq ( Texts a )
  {
    return ( cmp ( a ) == 0 );
  }

/*--------------------------------------------------------------------
 eq - is a equal to this Text object
--------------------------------------------------------------------*/

  public boolean[] eq ( Text a )
  {
    boolean b[] = new boolean[size];  // Boolean flag array

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < size; i++ )
    {
      // IF THE TEXT IS EQUAL
      if ( get ( i ).eq ( a ) )
      {
        b[ i ] = true;  // Set the flag to true
      }

      // IF THE TEXT HAS NOT BEEN ENTERED
      else
      {
        b[ i ] = false;  // Set the flag to false
      }
    }

    return b;
  }

/*--------------------------------------------------------------------
 ne - is a not equal to this Texts object
--------------------------------------------------------------------*/

  public boolean ne ( Texts a )
  {
    return ( cmp ( a ) != 0 );
  }

/*--------------------------------------------------------------------
 ne - is a not equal to this Text object
--------------------------------------------------------------------*/

  public boolean[] ne ( Text a )
  {
    boolean b[] = new boolean[size];  // Boolean flag array

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < size; i++ )
    {
      // IF THE TEXT IS EQUAL
      if ( get ( i ).ne ( a ) )
      {
        b[ i ] = true;  // Set the flag to true
      }

      // IF THE TEXT HAS NOT BEEN ENTERED
      else
      {
        b[ i ] = false;  // Set the flag to false
      }
    }

    return b;
  }

/*--------------------------------------------------------------------
 lt - is a less than this Texts object
--------------------------------------------------------------------*/

  public boolean lt ( Texts a )
  {
    return ( cmp ( a ) < 0 );
  }

/*--------------------------------------------------------------------
 le - is a less than or equal to this Texts object
--------------------------------------------------------------------*/

  public boolean le ( Texts a )
  {
    return ( cmp ( a ) <= 0 );
  }

/*--------------------------------------------------------------------
 gt - is a greater than this Texts object
--------------------------------------------------------------------*/

  public boolean gt ( Texts a )
  {
    return ( cmp ( a ) > 0 );
  }

/*--------------------------------------------------------------------
 ge - is a greater than or equal to this Texts object
--------------------------------------------------------------------*/

  public boolean ge ( Texts a )
  {
    return ( cmp ( a ) >= 0 );
  }

/*--------------------------------------------------------------------
 cmp - com.scanfeld.cluster.Compare two texts objects
--------------------------------------------------------------------*/

  public int cmp ( Texts a )
  {
    int ans = 0;

    for ( int i = 0; i < a.size && i < size; i++ )
    {
      if ( x[ k + i ].ne ( a.x[ a.k + i ] ) )
      {
        if ( x[ k + i ].lt ( a.x[ a.k + i ] ) )
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
      if ( size < a.size )
      {
        ans = -1;
      }

      else if ( a.size < size )
      {
        ans = 1;
      }
    }

    return ans;
  }

/*--------------------------------------------------------------------
 init - Initialization function
--------------------------------------------------------------------*/

  public void init ()
  {
    t = new Tree ();        // Create a lookup tree
    ml = 2.0;               // Multiply amount
    mr = 2.0;               // Multiply amount
    k = 10;                 // Initialize to center of the array
    space = 20;             // Initialize to 20 characters
    x = new Text[space];  // The Text array
    clear ();               // Clear the Text object
  }

/*--------------------------------------------------------------------
 clear - Clear the text object
--------------------------------------------------------------------*/

  public void clear ()
  {
    size = 0;  // Start with an empty array

    // FILL THE ARRAY WITH SPACES
    for ( int i = 0; i < space; i++ )
    {
      x[ i ] = null;  // Set the ith character to a space
    }
  }

/*--------------------------------------------------------------------
 copy - Copy the texts object
--------------------------------------------------------------------*/

  public Texts copy ()
  {
    return new Texts ( this );
  }

/*--------------------------------------------------------------------
 alphaNum - Return just alphanumeric characters
--------------------------------------------------------------------*/

  public Texts alphaNum ()
  {
    Texts n = new Texts ();

    for ( int i = 0; i < size; i++ )
    {
      n.rpush ( get ( i ).alphaNum () );
    }

    return n;
  }

/*--------------------------------------------------------------------
 alphaNumSign - Return just alphanumeric characters and sign
--------------------------------------------------------------------*/

  public Texts alphaNumSign ()
  {
    Texts n = new Texts ();

    for ( int i = 0; i < size; i++ )
    {
      n.rpush ( get ( i ).alphaNumSign () );
    }

    return n;
  }

/*--------------------------------------------------------------------
 trimWhite - Trim the whitespace
--------------------------------------------------------------------*/

  public Texts trimWhite ()
  {
    Texts n = new Texts ();

    for ( int i = 0; i < size; i++ )
    {
      n.rpush ( get ( i ).trimWhite () );
    }

    return n;
  }

/*--------------------------------------------------------------------
 removeTag - Remove any enclosing tag
--------------------------------------------------------------------*/

  public Texts removeTag ()
  {
    Texts n = new Texts ();

    for ( int i = 0; i < size; i++ )
    {
      n.rpush ( get ( i ).removeTag () );
    }

    return n;
  }

/*--------------------------------------------------------------------
 sub - Select a substring of the array
--------------------------------------------------------------------*/

  public Texts sub ( int c, int d )
  {
    d = lechk ( d, size );  // Ensure d <= size
    d = gechk ( d, 0 );     // Ensure d >= 0
    c = lechk ( c, d );     // Ensure c <= d
    c = gechk ( c, 0 );     // Ensure c >= 0
    int newsize = d - c;    // Size of the new array

    Texts a = new Texts ();

    // LOOP THROUGH THE ARRAY
    for ( int i = c; i < d; i++ )
    {
      a.rpush ( get ( i ) );
    }

    return a;          // Return the new array
  }

/*--------------------------------------------------------------------
 left - Return to the left of a string
--------------------------------------------------------------------*/

  public Texts left ( String z )
  {
    return left ( new Text ( z ) );
  }

/*--------------------------------------------------------------------
 left - Return to the left of a text
--------------------------------------------------------------------*/

  public Texts left ( Text z )
  {
    Texts n = new Texts ();

    for ( int i = 0; i < size; i++ )
    {
      n.rpush ( get ( i ).left ( z ) );
    }

    return n;
  }

/*--------------------------------------------------------------------
 upper - Convert to uppercase
--------------------------------------------------------------------*/

  public Texts upper ()
  {
    Texts n = new Texts ();

    for ( int i = 0; i < size; i++ )
    {
      n.rpush ( get ( i ).upper () );
    }

    return n;
  }

/*--------------------------------------------------------------------
 lower - Convert to lowercase
--------------------------------------------------------------------*/

  public Texts lower ()
  {
    Texts n = new Texts ();

    for ( int i = 0; i < size; i++ )
    {
      n.rpush ( get ( i ).lower () );
    }

    return n;
  }

/*--------------------------------------------------------------------
 removeExtraWhite - Remove white space
--------------------------------------------------------------------*/

  public Texts removeWhite ()
  {
    Texts n = new Texts ();

    for ( int i = 0; i < size; i++ )
    {
      n.rpush ( get ( i ).removeExtraWhite () );
    }

    return n;
  }

/*--------------------------------------------------------------------
 split - Split texts by z
--------------------------------------------------------------------*/

  public Text[][] split ( String z )
  {
    return split ( new Text ( z ) );
  }

/*--------------------------------------------------------------------
 splitSpace - Split texts by a space
--------------------------------------------------------------------*/

  public Text[][] splitSpace ()
  {
    return split ( new Text ( " " ) );
  }

/*--------------------------------------------------------------------
 splitTab - Split texts by a tab
--------------------------------------------------------------------*/

  public Text[][] splitTab ()
  {
    return split ( new Text ( "\t" ) );
  }

/*--------------------------------------------------------------------
 split - Split texts by a space
--------------------------------------------------------------------*/

  public Text[][] splitLine ()
  {
    return split ( new Text ( "\n" ) );
  }

/*--------------------------------------------------------------------
 split - Split texts by z
--------------------------------------------------------------------*/

  public Text[][] split ( Text z )
  {
    Text[][] n = new Text[size][];

    for ( int i = 0; i < size; i++ )
    {
      n[ i ] = get ( i ).split ( z );
    }

    return n;
  }

/*--------------------------------------------------------------------
 reverse - Reverse the texts object
--------------------------------------------------------------------*/

  public void reverse ()
  {
    // FIND THE MIDPOINT
    int mid = size / 2;

    if ( mid * 2 == size )
    {
      mid--;
    }

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i <= mid; i++ )
    {
      Text temp = x[ k + i ];
      x[ k + i ] = x[ k + size - 1 - i ];
      x[ k + size - 1 - i ] = temp;
    }
  }

/*--------------------------------------------------------------------
 length - Return the size of the Texts object
--------------------------------------------------------------------*/

  public int length ()
  {
    return size;
  }

/*--------------------------------------------------------------------
 lengths - Return the sizes of the Text objects
--------------------------------------------------------------------*/

  public int[] lengths ()
  {
    int[] z = new int[size];

    for ( int i = 0; i < size; i++ )
    {
      z[ i ] = get ( i ).size;
    }

    return z;
  }

/*--------------------------------------------------------------------
 maxLength - Return the maximum length in the Texts
--------------------------------------------------------------------*/

  public int maxLength ()
  {
    return max ( lengths () );
  }

/*--------------------------------------------------------------------
 minLength - Return the minimum length in the Texts
--------------------------------------------------------------------*/

  public int minLength ()
  {
    return min ( lengths () );
  }

/*--------------------------------------------------------------------
 isEmpty - Check if the array is empty
--------------------------------------------------------------------*/

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

/*--------------------------------------------------------------------
 sim - Calculate the overlap between a and this texts object
--------------------------------------------------------------------*/

  public double sim ( Texts a )
  {
    Tree z = a.tree ();  // Lookup tree for a
    int j = 0;          // Match count

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < size; i++ )
    {
      if ( z.get ( x[ k + i ] ) > 0 )
      {
        j++;  // Increment the match count
      }
    }

    return ( ( double ) ( j ) / ( double ) ( size ) );  // Return the match percentage
  }

/*--------------------------------------------------------------------
 size - Resize the array
--------------------------------------------------------------------*/

  public void size ( int newspace )
  {
    size ( newspace, 0 );
  }

/*--------------------------------------------------------------------
 size - Resize the array
--------------------------------------------------------------------*/

  public void size ( int newspace, int offset )
  {
    // IF THE ARRAY IS TOO SMALL
    if ( newspace * ( 1.0 + ml * 0.5 + mr * 0.5 ) > space )
    {
      int oldspace = space;                              // Remember the size of the original array
      space = ( int ) ( newspace * ( 1.0 + ml + mr ) );  // Update to the new size
      Text[] y = new Text[space];                      // Create the new array
      int nk = ( space / 2 ) - ( size + offset / 2 );    // Find the new k
      nk = nk + offset;                                  // Add the offset

      // COPY THE PREVIOUS ARRAY TO THE NEW ONE
      for ( int i = 0; i < size; i = i + 1 )
      {
        y[ nk + i ] = x[ k + i ];  // Set the ith character
      }

      // FILL THE REST OF THE ARRAY WITH SPACES
      for ( int i = 0; i < nk; i = i + 1 )
      {
        y[ i ] = null;  // Set the ith character to a space
      }

      // FILL THE REST OF THE ARRAY WITH SPACES
      for ( int i = nk + size; i < space; i++ )
      {
        y[ i ] = null;  // Set the ith character to a space
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
        x[ i ] = null;  // Set the ith character to a space
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
        x[ i ] = null;  // Set the ith character to a space
      }

      k = nk;  // Save the new start
    }
  }

/*--------------------------------------------------------------------
 PUSH / POP FUNCTIONS
 ---------------------------------------------------------------------
 Add and remove letters from the object
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 lpop - Pop a Text from the left side
--------------------------------------------------------------------*/

  public Text lpop ( int n )
  {
    for ( int i = 0; i < n - 1; i++ )
    {
      lpop ();
    }

    return lpop ();
  }

/*--------------------------------------------------------------------
 lpop - Pop a Text from the left side
--------------------------------------------------------------------*/

  public Text lpop ()
  {
    // IF THE ARRAY IS NOT EMPTY
    if ( size > 0 )
    {
      Text c = x[ k ];  // Save the first Text
      x[ k ] = null;    // Set the Text to a space
      k = k + 1;        // Increment the start index
      size = size - 1;  // Update the size of the array
      return c;         // Return the first Text
    }

    // IF THE ARRAY IS EMPTY
    else
    {
      return null;  // Return the null Text
    }
  }

/*--------------------------------------------------------------------
 rpop - Pop a Text from the right side
--------------------------------------------------------------------*/

  public Text rpop ( int n )
  {
    for ( int i = 0; i < n - 1; i++ )
    {
      rpop ();
    }

    return rpop ();
  }

/*--------------------------------------------------------------------
 rpop - Pop a Text from the right side
--------------------------------------------------------------------*/

  public Text rpop ()
  {
    // IF THE ARRAY IS NOT EMPTY
    if ( size > 0 )
    {
      Text c = x[ k + size - 1 ];  // Save the last Text
      x[ k + size - 1 ] = null;    // Set the charatcter to a space
      size = size - 1;             // Update the size of the array
      return c;                    // Return the last Text
    }

    // IF THE ARRAY IS EMPTY
    else
    {
      return null;  // Return the null Text
    }
  }

/*--------------------------------------------------------------------
 tree - Create a lookup tree
--------------------------------------------------------------------*/

  public Tree tree ()
  {
    Tree z = new Tree ();  // Lookup tree

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < size; i++ )
    {
      z.set ( x[ k + i ], ( i + 1 ) );
    }

    return z; // Return the tree
  }

/*--------------------------------------------------------------------
 lpush - Add a Text to the left side
--------------------------------------------------------------------*/

  public void lpush ( Text z )
  {
    size ( size + 1, 1 );   // Ensure the array is the desired size
    x[ k - 1 ] = z;         // Set the Text
    k = k - 1;              // Update k
    size = size + 1;        // Store the new size
  }

/*--------------------------------------------------------------------
 lpush - Add a Text array to the left side
--------------------------------------------------------------------*/

  public void lpush ( Text[] z )
  {
    int length = z.length;        // Length of the Text array
    int newsize = size + length;  // New size of the array
    size ( newsize, length );     // Ensure the array is the desired size

    // ADD THE TEXT ARRAY TO THE LEFT SIDE OF THE ARRAY
    for ( int i = 0; i < length; i = i + 1 )
    {
      x[ k - length + i ] = z[ i ];  // Set the ith Text
    }

    k = k - length;        // Update k
    size = size + length;  // Store the new size
  }

/*--------------------------------------------------------------------
 lpush - Add a Texts to the left side
--------------------------------------------------------------------*/

  public void lpush ( Texts z )
  {
    int length = z.size;          // Length of the Text array
    int newsize = size + length;  // New size of the array
    size ( newsize, length );     // Ensure the array is the desired size

    // ADD THE TEXT ARRAY TO THE LEFT SIDE OF THE ARRAY
    for ( int i = 0; i < length; i = i + 1 )
    {
      x[ k - length + i ] = z.get ( i );  // Set the ith Text
    }

    k = k - length;        // Update k
    size = size + length;  // Store the new size
  }

/*--------------------------------------------------------------------
 lpush - Add a character to the left side
--------------------------------------------------------------------*/

  public void lpush ( char z )
  {
    lpush ( new Text ( z ) );  // Push the character as a Text
  }

/*--------------------------------------------------------------------
 lpush - Add a character array to the left side
--------------------------------------------------------------------*/

  public void lpush ( char[] z )
  {
    lpush ( new Text ( z ) );  // Push the character array as a Text
  }

/*--------------------------------------------------------------------
 lpush - Add a double to the left side
--------------------------------------------------------------------*/

  public void lpush ( double z )
  {
    lpush ( new Text ( z ) );  // Push the double as a Text
  }

/*--------------------------------------------------------------------
 lpush - Add an integer to the left side
--------------------------------------------------------------------*/

  public void lpush ( int z )
  {
    lpush ( new Text ( z ) );  // Push the integer as a Text
  }

/*--------------------------------------------------------------------
 lpush - Add a string to the left side
--------------------------------------------------------------------*/

  public void lpush ( String z )
  {
    lpush ( new Text ( z ) );  // Push the string as a Text
  }

/*--------------------------------------------------------------------
 rpush - Add a Text to the right side
--------------------------------------------------------------------*/

  public void rpush ( Text z )
  {
    size ( size + 1, -1 );  // Ensure the array is the desired size
    x[ k + size ] = z;      // Set the Text
    size++;                 // Store the new size
  }

/*--------------------------------------------------------------------
 rpush - Add a Text array to the right side
--------------------------------------------------------------------*/

  public void rpush ( Text[] z )
  {
    int length = z.length;        // Length of the Text array
    int newsize = size + length;  // New size of the array
    size ( newsize, length );     // Ensure the array is the desired size

    // ADD THE TEXT ARRAY TO THE RIGHT SIDE OF THE ARRAY
    for ( int i = 0; i < length; i = i + 1 )
    {
      x[ k + size + i ] = z[ i ];  // Set the ith Text
    }

    size = size + length;  // Store the new size
  }

/*--------------------------------------------------------------------
 rpush - Add an int array to the right side
--------------------------------------------------------------------*/

  public void rpush ( int[] z )
  {
    int length = z.length;        // Length of the Text array
    int newsize = size + length;  // New size of the array
    size ( newsize, length );     // Ensure the array is the desired size

    // ADD THE TEXT ARRAY TO THE RIGHT SIDE OF THE ARRAY
    for ( int i = 0; i < length; i = i + 1 )
    {
      x[ k + size + i ] = i2t ( z[ i ] );  // Set the ith Text
    }

    size = size + length;  // Store the new size
  }

/*--------------------------------------------------------------------
 rpush - Add a Texts to the right side
--------------------------------------------------------------------*/

  public void rpush ( Texts z )
  {
    size ( size + z.size, -z.size );   // Ensure the array is the desired size

    // LOOP THROUGH THE TEXTS OBJECT
    for ( int i = 0; i < z.size; i++ )
    {
      x[ k + size + i ] = z.x[ z.k + i ];  // Set the item
    }

    size += z.size;                    // Store the new size
  }

/*--------------------------------------------------------------------
 rpush - Add a character to the right side
--------------------------------------------------------------------*/

  public void rpush ( char z )
  {
    rpush ( new Text ( z ) );  // Push the character as a Text
  }

/*--------------------------------------------------------------------
 rpush - Add a character array to the right side
--------------------------------------------------------------------*/

  public void rpush ( char[] z )
  {
    rpush ( new Text ( z ) );  // Push the character array as a Text
  }

/*--------------------------------------------------------------------
 rpush - Add a double to the right side
--------------------------------------------------------------------*/

  public void rpush ( double z )
  {
    rpush ( new Text ( z ) );  // Push the double as a Text
  }

/*--------------------------------------------------------------------
 rpush - Add an integer to the right side
--------------------------------------------------------------------*/

  public void rpush ( int z )
  {
    rpush ( new Text ( z ) );  // Push the integer as a Text
  }

/*--------------------------------------------------------------------
 rpush - Add a string to the right side
--------------------------------------------------------------------*/

  public void rpush ( String z )
  {
    rpush ( new Text ( z ) );  // Push the string as a Text
  }

/*--------------------------------------------------------------------
 GET / SET FUNCTIONS
 ---------------------------------------------------------------------
 These functions manipulate the Text objects in the array
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 get - Get a Text
--------------------------------------------------------------------*/

  public Text get ( int i )
  {
    // IF i IS WITHIN RANGE
    if ( ( i >= 0 ) && ( i < size ) )
    {
      return x[ k + i ];  // Return the ith character
    }

    // IF i IS OUT OF RANGE
    else
    {
      return null;  // Return the null Text
    }
  }

/*--------------------------------------------------------------------
 find - Find a string
--------------------------------------------------------------------*/

  public int find ( String z )
  {
    return find ( new Text ( z ) );
  }

/*--------------------------------------------------------------------
 find - Find an integer
--------------------------------------------------------------------*/

  public int find ( int z )
  {
    return find ( i2t ( z ) );
  }

/*--------------------------------------------------------------------
 find - Find a text
--------------------------------------------------------------------*/

  public int find ( Text z )
  {
    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < size; i++ )
    {
      // IF THE TEXT IS FOUND
      if ( x[ k + i ].eq ( z ) )
      {
        return i;  // Return the index
      }
    }

    return -1;  // If not found, return -1
  }

/*--------------------------------------------------------------------
 finds - Find a set of strings
--------------------------------------------------------------------*/

  public int[] finds ( String z )
  {
    return finds ( new Text ( z ) );
  }

/*--------------------------------------------------------------------
 finds - Find a set of integers
--------------------------------------------------------------------*/

  public int[] finds ( int z )
  {
    return finds ( i2t ( z ) );
  }

/*--------------------------------------------------------------------
 finds - Find a set of texts
--------------------------------------------------------------------*/

  public int[] finds ( Text z )
  {
    Texts ind = new Texts ();

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < size; i++ )
    {
      // IF THE TEXT IS FOUND
      if ( x[ k + i ].eq ( z ) )
      {
        ind.rpush ( i );
      }
    }

    return ind.num ();  // If not found, return -1
  }

/*--------------------------------------------------------------------
 toString - Get a string
--------------------------------------------------------------------*/

  public String toString ()
  {
    return txt ().str ();  // Return the Texts as a string
  }

/*--------------------------------------------------------------------
 out - Show the Text array
--------------------------------------------------------------------*/

  public void out ()
  {
    out ( "Text com.scanfeld.bio.Array (" + size + ")" );

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < size; i++ )
    {
      out ( i + ": " + x[ k + i ] );  // Output the ith element
    }

    out ( "" );  // Add a blank line at the end
  }

/*--------------------------------------------------------------------
 outLine - Output the text array to a file line by line
--------------------------------------------------------------------*/

  public void outLine ( String f )
  {
    Outfile a = new Outfile ( f );
    a.write ( joinLine () );
    a.close ();
  }

/*--------------------------------------------------------------------
 get - Get a string
--------------------------------------------------------------------*/

  public String get ()
  {
    return txt ().str ();
  }

/*--------------------------------------------------------------------
 chr - Get a character array
--------------------------------------------------------------------*/

  public char[] chr ()
  {
    int fullSize = 0;  // Store the full size of the set of Texts

    // LOOP THROUGH THE ARRAY
    for ( int i = k; i < k + size; i = i + 1 )
    {
      fullSize += x[ i ].size;  // Increment the total size counter
    }

    char[] y = new char[fullSize];  // Create a new character array
    int yind = 0;                     // Create a y index

    // LOOP THROUGH THE ARRAY
    for ( int i = k; i < k + size; i = i + 1 )
    {
      // LOOP THROUGH THE TEXT
      for ( int j = x[ i ].k; j < x[ i ].k + x[ i ].size; j++ )
      {
        y[ yind ] = x[ i ].x[ j ];  // Set the y index character
        yind++;                     // Increment the y index
      }
    }

    return y;  // Return the new array
  }

/*--------------------------------------------------------------------
 bytes - Get a byte array
--------------------------------------------------------------------*/

  public byte[] bytes ()
  {
    int fullSize = 0;  // Store the full size of the set of Texts

    // LOOP THROUGH THE ARRAY
    for ( int i = k; i < k + size; i = i + 1 )
    {
      fullSize += x[ i ].size;  // Increment the total size counter
    }

    byte[] y = new byte[fullSize];  // Create a new character array
    int yind = 0;                     // Create a y index

    // LOOP THROUGH THE ARRAY
    for ( int i = k; i < k + size; i++ )
    {
      // LOOP THROUGH THE TEXT
      for ( int j = x[ i ].k; j < x[ i ].k + x[ i ].size; j++ )
      {
        y[ yind ] = ( byte ) x[ i ].x[ j ];  // Set the yind character
        yind++;                              // Increment the y index
      }
    }

    return y;  // Return the new array
  }

/*--------------------------------------------------------------------
 join - Join into one Text
--------------------------------------------------------------------*/

  public Text join ( String y )
  {
    return join ( new Text ( y ) );
  }

/*--------------------------------------------------------------------
 join - Join into one Text
--------------------------------------------------------------------*/

  public Text join ( Text y )
  {
    Text z = new Text ();  // Text to hold the result

    // IF THERE IS AT LEAST ONE ITEM
    if ( size > 0 )
    {
      z.rpush ( x[ k ] );  // Push the first item

      // LOOP THROUGH THE ELEMENTS
      for ( int i = 1; i < size; i++ )
      {
        z.rpush ( y );           // Push the joiner
        z.rpush ( x[ k + i ] );  // Push the ith element
      }
    }

    return z;  // Return the Text object
  }

/*--------------------------------------------------------------------
 joinSpace - Join into one Text
--------------------------------------------------------------------*/

  public Text joinSpace ()
  {
    return join ( " " );
  }

/*--------------------------------------------------------------------
 joinTab - Join into one Text
--------------------------------------------------------------------*/

  public Text joinTab ()
  {
    return join ( "\t" );
  }

/*--------------------------------------------------------------------
 joinLine - Join into one Text
--------------------------------------------------------------------*/

  public Text joinLine ()
  {
    return join ( "\n" );
  }

/*--------------------------------------------------------------------
 str - Get a string array
--------------------------------------------------------------------*/

  public String[] str ()
  {
    String[] z = new String[size];  // Create a string array

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < size; i++ )
    {
      z[ i ] = x[ k + i ].str ();  // Set the ith element
    }

    return z;  // Return the string array
  }

/*--------------------------------------------------------------------
 num - Get an integer array
--------------------------------------------------------------------*/

  public int[] num ()
  {
    int[] z = new int[size];  // Create an integer array

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < size; i++ )
    {
      z[ i ] = x[ k + i ].num ();  // Set the ith element
    }

    return z;  // Return the integer array
  }

/*--------------------------------------------------------------------
 dbl - Get a double array
--------------------------------------------------------------------*/

  public double[] dbl ()
  {
    double[] z = new double[size];  // Create an double array

    // LOOP THROUGH THE ELEMENTS
    for ( int i = 0; i < size; i++ )
    {
      z[ i ] = x[ k + i ].dbl ();  // Set the ith element
    }

    return z;  // Return the double array
  }

/*--------------------------------------------------------------------
 txt - Get a Text
--------------------------------------------------------------------*/

  public Text txt ()
  {
    Text y = new Text ( "" );  // New Text object
    y.set ( chr () );          // Set to the character array
    return y;                  // Return the Text object
  }

/*--------------------------------------------------------------------
 set - Set a Text
--------------------------------------------------------------------*/

  public void set ( int i, Text z )
  {
    size ( i, 0 );   // Ensure the array is a desired size
    x[ k + i ] = z;  // Set the ith character to z

    // IF i IS OUTSIDE THE CURRENT SIZE
    if ( i >= size )
    {
      size = i + 1;  // Set the size to fit i
    }
  }

/*--------------------------------------------------------------------
 select - Select a subset
--------------------------------------------------------------------*/

  public Texts select ( int[] set )
  {
    Texts y = new Texts ();  // New Texts object

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < set.length; i++ )
    {
      if ( i >= 0 && set[ i ] < size )
      {
        y.rpush ( x[ set[ i ] + k ] );  // Push the ith element
      }
    }

    return y;  // Return the Texts object
  }

/*--------------------------------------------------------------------
 select - Select a subset
--------------------------------------------------------------------*/

  public Texts select ( boolean[] set )
  {
    Texts y = new Texts ();  // New Texts object

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < set.length; i++ )
    {
      if ( i >= 0 && i < size && set[ i ] )
      {
        y.rpush ( x[ k + i ] );  // Push the ith element
      }
    }

    return y;  // Return the Texts object
  }

/*--------------------------------------------------------------------
 removeDuplicate - Remove the duplicates
--------------------------------------------------------------------*/

  public Texts removeDuplicate ()
  {
    Tree z = new Tree ();               // Temporary lookup tree
    boolean b[] = new boolean[size];  // Boolean flag array

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < size; i++ )
    {
      // IF THE TEXT HAS ALREADY BEEN ENTERED
      if ( z.get ( x[ k + i ] ) > 0 )
      {
        b[ i ] = false;  // Set the flag to false
      }

      // IF THE TEXT HAS NOT BEEN ENTERED
      else
      {
        z.set ( x[ k + i ], ( i + 1 ) );  // Add the ith element
        b[ i ] = true;                    // Set the flag to true
      }
    }

    return select ( b );  // Return the Texts object
  }

/*--------------------------------------------------------------------
 removeDuplicateInd - Remove the duplicates return indices
--------------------------------------------------------------------*/

  public int[] removeDuplicateInd ()
  {
    Tree z = new Tree ();               // Temporary lookup tree
    boolean b[] = new boolean[size];  // Boolean flag array

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < size; i++ )
    {
      // IF THE TEXT HAS ALREADY BEEN ENTERED
      if ( z.get ( x[ k + i ] ) > 0 )
      {
        b[ i ] = false;  // Set the flag to false
      }

      // IF THE TEXT HAS NOT BEEN ENTERED
      else
      {
        z.set ( x[ k + i ], ( i + 1 ) );  // Add the ith element
        b[ i ] = true;                    // Set the flag to true
      }
    }

    return b2i( b );  // Return the indices
  }

/*--------------------------------------------------------------------
 unique - Get the unique
--------------------------------------------------------------------*/

  public boolean[] unique ()
  {
    Tree z = new Tree ();               // Temporary lookup tree
    boolean b[] = new boolean[size];  // Boolean flag array

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < size; i++ )
    {
      // IF THE TEXT HAS ALREADY BEEN ENTERED
      if ( z.get ( x[ k + i ] ) > 0 )
      {
        b[ i ] = false;  // Set the flag to false
      }

      // IF THE TEXT HAS NOT BEEN ENTERED
      else
      {
        z.set ( x[ k + i ], ( i + 1 ) );  // Add the ith element
        b[ i ] = true;                    // Set the flag to true
      }
    }

    return b;  // Return the unique indices
  }

/*--------------------------------------------------------------------
 select - Select the indices of a set of Strings
--------------------------------------------------------------------*/

  public int[] select ( String[] z )
  {
    Texts y = new Texts ( z );  // Create a Texts array of the strings
    return select ( y );        // Select the strings
  }

/*--------------------------------------------------------------------
 select - Select the indices of a set of Strings
--------------------------------------------------------------------*/

  public int[] select ( Text[] z )
  {
    Texts y = new Texts ( z );  // Create a Texts array of the strings
    return select ( y );        // Select the strings
  }

/*--------------------------------------------------------------------
 select - Select the indices of a string
--------------------------------------------------------------------*/

  public int[] select ( String y )
  {
    String[] z = new String[1];
    z[ 0 ] = y;
    return select ( z );
  }

/*--------------------------------------------------------------------
 select - Select the indices of a text
--------------------------------------------------------------------*/

  public int[] select ( Text y )
  {
    Text[] z = new Text[1];
    z[ 0 ] = y;
    return select ( z );
  }

/*--------------------------------------------------------------------
 select - Select the indices of a set of Texts
--------------------------------------------------------------------*/

  public int[] select ( Texts y )
  {
    Tree z = new Tree ();                // Temporary lookup tree
    Text[] f = new Text[y.length ()];    // Find Text
    int[][] g = new int[y.length ()][];  // Lookups

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < y.size; i++ )
    {
      int j = z.get ( y.x[ y.k + i ] );

      if ( j > 0 )
      {
        if ( f[ j - 1 ] == null )
        {
          f[ j - 1 ] = new Text ();
          f[ j - 1 ].rpush ( i );
        }

        else
        {
          f[ j - 1 ].rpush ( "\t" );
          f[ j - 1 ].rpush ( i );
        }
      }

      else
      {
        z.set ( y.x[ y.k + i ], ( i + 1 ) );  // Add the ith element
      }
    }

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < y.size; i++ )
    {
      if ( f[ i ] != null )
      {
        Texts tt = new Texts ( f[ i ].splitTab () );
        g[ i ] = tt.num ();
      }
    }

    return select ( z, g );  // Return the indices of the match
  }

/*--------------------------------------------------------------------
 select - Select the indices of a set of Texts
--------------------------------------------------------------------*/

  public int[] select ( Tree z, int[][] y )
  {
    Text ft = new Text ();  // Found text
    Text fi = new Text ();  // Found indices
    int f = 0;              // Find integer

    // LOOP THROUGH THE ARRAY
    for ( int i = 0; i < size; i++ )
    {
      f = z.get ( x[ k + i ] );  // Look for the string

      // IF THE TEXT MATCHES
      if ( f > 0 )
      {
        if ( ft.size > 0 )
        {
          ft.rpush ( "\t" );
        }

        ft.rpush ( i );  // Push the index

        if ( fi.size > 0 )
        {
          fi.rpush ( "\t" );
        }

        fi.rpush ( f );  // Push the order

        if ( y[ f - 1 ] != null )
        {
          for ( int j = 0; j < y[ f - 1 ].length; j++ )
          {
            // out ( "match: " + x[ k + i ] );
            ft.rpush ( "\t" );
            ft.rpush ( i );
            fi.rpush ( "\t" );
            fi.rpush ( y[ f - 1 ][ j ] );
          }
        }
      }
    }

    Texts ftt = new Texts ( ft.splitTab () );
    Texts fit = new Texts ( fi.splitTab () );
    int[] fitord = sort ( fit.num () );
    ftt = ftt.select ( fitord );

    return ftt.num ();  // Return the indices
  }

/*--------------------------------------------------------------------
 match - Match to a Texts
--------------------------------------------------------------------*/

  public Texts match ( Texts z )
  {
    int[] ind = select ( z );
    return select ( ind );
  }

/*--------------------------------------------------------------------
 map - Map using a TextTree
--------------------------------------------------------------------*/

  public Texts map ( TextTree t )
  {
    Texts z = new Texts ();

    for ( int i = 0; i < size; i++ )
    {
      Text n = t.get ( get ( i ) );

      if ( n != null )
      {
        z.rpush ( n );
      }
    }

    return z;
  }

/*--------------------------------------------------------------------
 map - Map using a TextTree
--------------------------------------------------------------------*/

  public Texts map ( String f )
  {
    TextTree t = new TextTree ( f );

    Texts z = new Texts ();

    for ( int i = 0; i < size; i++ )
    {
      Text n = t.get ( get ( i ) );

      if ( n != null )
      {
        z.rpush ( n );
      }
    }

    return z;
  }

/*--------------------------------------------------------------------
 mapInd - Map using a TextTree and return indices
--------------------------------------------------------------------*/

  public int[] mapInd ( String f )
  {
    TextTree t = new TextTree ( f );

    Texts z = new Texts ();

    for ( int i = 0; i < size; i++ )
    {
      Text n = t.get ( get ( i ) );

      if ( n != null )
      {
        z.rpush ( i );
      }
    }

    return z.num ();
  }

/*--------------------------------------------------------------------
 mapMult - Map using multiple hits
--------------------------------------------------------------------*/

  public Texts mapMult ( String f )
  {
    TextTree t = new TextTree ();
    t.mapMult ( f );

    Texts z = new Texts ();

    for ( int i = 0; i < size; i++ )
    {
      Text[] n = t.getMult ( get ( i ) );

      if ( n.length > 0 )
      {
        z.rpush ( n );
      }
    }

    return z;
  }

/*--------------------------------------------------------------------
 mapMultInd - Map using multiple hits and return indices
--------------------------------------------------------------------*/

  public int[] mapMultInd ( String f )
  {
    TextTree t = new TextTree ();
    t.mapMult ( f );

    Texts z = new Texts ();

    for ( int i = 0; i < size; i++ )
    {
      Text[] n = t.getMult ( get ( i ) );

      for ( int j = 0; j < n.length; j++ )
      {
        z.rpush ( i );
      }
    }

    return z.num ();
  }

/*--------------------------------------------------------------------
 sortd - Reverse sort the indices of an array of Texts
--------------------------------------------------------------------*/

  public int[] sortd ()
  {
    int[] n = sort ();
    reverse ();
    return rev ( n );
  }

/*--------------------------------------------------------------------
 sort - Sort the indices of an array of Texts
--------------------------------------------------------------------*/

  public int[] sort ()
  {
    int[] a = new int[size];  // First index array
    int[] b = new int[size];  // Second index array

    // INITIALIZE THE ARRAYS
    for ( int i = 0; i < size; i = i + 1 )
    {
      a[ i ] = i;  // Set the ith item to i
      b[ i ] = i;  // Set the ith item to i
    }

    Text[] y = new Text[size];          // Temporary Text array
    mergeSort ( a, b, y, 0, size - 1 );   // Merge sort the array
    return a;                             // Return the indices
  }

/*--------------------------------------------------------------------
 mergeSort - Merge sort the indices of an array of Texts
--------------------------------------------------------------------*/

  public void mergeSort ( int[] a, int[] b, Text[] q, int min, int max )
  {
    // FIND THE SIZE OF THE CURRENT SUBARRAY
    int sze = max - min + 1;

    // IF THERE ARE MORE THAN TWO ELEMENTS
    if ( sze > 2 )
    {
      // SPLIT THE ARRAY
      int s = min + ( sze / 2 );          // Find the split point
      mergeSort ( a, b, q, min, s - 1 );  // Merge sort the lower half
      mergeSort ( a, b, q, s, max );      // Merge sort the upper half

      // INITIALIZE THE LOOPING INDICES
      int i = min;  // Index for the first split
      int j = s;    // Index for the second split
      int l = min;  // Index for the new ordering

      // WHILE THERE ARE COMPARISONS TO MAKE
      while ( ( i < s ) || ( j <= max ) )
      {
        // IF THE FIRST AND SECOND ARRAYS EACH HAVE ITEMS LEFT
        if ( ( i < s ) && ( j <= max ) )
        {
          // IF THE FIRST IS LESS THAN OR EQUAL TO THE SECOND
          if ( q[ i ].le ( q[ j ] ) )
          {
            x[ k + l ] = q[ i ];  // Store the first value
            a[ l ] = b[ i ];  // Store the first index
            i = i + 1;        // Increment the i index
          }

          // IF THE FIRST IS GREATER THAN THE SECOND
          else
          {
            x[ k + l ] = q[ j ];  // Store the second value
            a[ l ] = b[ j ];  // Store the second index
            j = j + 1;        // Increment the j index
          }
        }

        // IF I IS LESS THAN THE SPLIT
        else if ( i < s )
        {
          x[ k + l ] = q[ i ];  // Store the first value
          a[ l ] = b[ i ];  // Store the first index
          i = i + 1;        // Increment the i index
        }

        // IF I IS AT THE SPLIT
        else
        {
          x[ k + l ] = q[ j ];  // Store the second value
          a[ l ] = b[ j ];  // Store the second index
          j = j + 1;        // Increment the j index
        }

        l = l + 1;  // Increment the k index
      }

      // LOOP FROM MIN TO MAX
      for ( int h = 0; h <= max; h = h + 1 )
      {
        q[ h ] = x[ h + k ];  // Store the value in q
        b[ h ] = a[ h ];  // Store the index in b
      }
    }

    // IF THERE ARE TWO ITEMS IN THE ARRAY
    else if ( sze == 2 )
    {
      // IF THE FIRST ITEM IS LESS THAN OR EQUAL TO THE SECOND
      if ( x[ k + min ].le ( x[ k + min + 1 ] ) )
      {
        q[ min ] = x[ k + min ];          // Store the first value in q
        q[ min + 1 ] = x[ k + min + 1 ];  // Store the second value in q
        b[ min ] = a[ min ];              // Store the first index in b
        b[ min + 1 ] = a[ min + 1 ];      // Store the second index in b
      }

      // IF THE FIRST ITEM IS GREATER THAN THE SECOND
      else
      {
        // SWITCH THE VALUES IN Q
        q[ min ] = x[ k + min + 1 ];
        q[ min + 1 ] = x[ k + min ];

        // STORE THE VALUES BACK IN P
        x[ k + min ] = q[ min ];
        x[ k + min + 1 ] = q[ min + 1 ];

        // SWITCH THE INDICES IN B
        b[ min ] = a[ min + 1 ];
        b[ min + 1 ] = a[ min ];

        // STORE THE INDICES BACK IN A
        a[ min ] = b[ min ];
        a[ min + 1 ] = b[ min + 1 ];
      }
    }

    // IF THERE IS ONLY ONE ITEM IN THE ARRAY
    else if ( sze == 1 )
    {
      q[ min ] = x[ k + min ];  // Store the value in q
      b[ min ] = a[ min ];      // Store the index in b
    }
  }
}