package broad.projection.utils;



/*--------------------------------------------------------------------
 TextTree
 ---------------------------------------------------------------------
 Models a tree for text lookups
 ---------------------------------------------------------------------
 Author: Daniel Scanfeld
   Date: 08.31.06
 ---------------------------------------------------------------------

--------------------------------------------------------------------*/

public class TextTree extends Inc
{
/*--------------------------------------------------------------------
 Class Variables
--------------------------------------------------------------------*/

  public char[] x;       // Remainder of the string
  public Text v;         // Value
  public char[] y;       // Lookups
  public TextTree[] t;   // TextTree
  public TextTree next;  // Next TextTree
  public TextTree prev;  // Previous TextTree
  public int s;          // Size of array

/*--------------------------------------------------------------------
 CONSTRUCTOR FUNCTIONS
 ---------------------------------------------------------------------
 Constructors for the possible text objects
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 TextTree - Null constructor
--------------------------------------------------------------------*/

  public TextTree ()
  {
    s = 0;  // Initial size of 0
  }

/*--------------------------------------------------------------------
 TextTree - File constructor
--------------------------------------------------------------------*/

  public TextTree ( String fname )
  {
    s = 0;  // Initial size of 0
    map ( fname );
  }

/*--------------------------------------------------------------------
 getMult - GetMult the value for z
--------------------------------------------------------------------*/

  public Text[] getMult ( char z )
  {
    return getMult ( h2c ( z ) );  // Run getMult with a character array
  }

/*--------------------------------------------------------------------
 getMult - GetMult the value for z
--------------------------------------------------------------------*/

  public Text[] getMult ( Text z )
  {
    return getMult ( z.chr () );  // Run getMult with a character array
  }

/*--------------------------------------------------------------------
 getMult - GetMult the value for z
--------------------------------------------------------------------*/

  public Text[] getMult ( String z )
  {
    return getMult ( s2c ( z ) );  // Run getMult with a character array
  }

/*--------------------------------------------------------------------
 getMult - Get the value for z
--------------------------------------------------------------------*/

  public Text[] getMult ( char[] z )
  {
    Text t = get ( z, 0 );

    if ( t == null )
    {
      return new Text[0];
    }

    else
    {
      t = t.right ( -1 );
      return t.split ( "|" );
    }
  }

/*--------------------------------------------------------------------
 get - Get the value for z
--------------------------------------------------------------------*/

  public Text get ( char z )
  {
    return get ( h2c ( z ), 0 );  // Run get with a character array
  }

/*--------------------------------------------------------------------
 get - Get the value for z
--------------------------------------------------------------------*/

  public Text get ( Text z )
  {
    return get ( z.chr (), 0 );  // Run get with a character array
  }

/*--------------------------------------------------------------------
 get - Get the value for z
--------------------------------------------------------------------*/

  public Text get ( String z )
  {
    return get ( s2c ( z ), 0 );  // Run get with a character array
  }

/*--------------------------------------------------------------------
 get - Get the value for z
--------------------------------------------------------------------*/

  public Text get ( char[] z, int ind )
  {
    // IF NO LETTERS LEFT
    if ( ind == z.length )
    {
      // IF X IS NULL OR EMPTY
      if ( x == null || x.length == 0 )
      {
        return v;  // Return the value
      }
    }

    // ELSE IF THE TREE IS NULL
    else if ( t == null )
    {
      // IF THE STRINGS ARE EQUAL
      if ( x != null && eq ( z, ind, x, 0 ) )
      {
        return v;  // Return the value
      }

      // IF THE STRINGS ARE NOT EQUAL
      else
      {
        return null;  // Return null
      }
    }

    // IF THE TREE EXISTS
    else
    {
      int i = find ( z[ ind ] );  // Find the branch

      // IF A BRANCH WAS FOUND
      if ( i < s )
      {
        return t[ i ].get ( z, ind + 1 );  // Return the value
      }
    }

    return null;
  }

/*--------------------------------------------------------------------
 find - Find the character
--------------------------------------------------------------------*/

  public int find ( char a )
  {
    int min = 0;                // Start min at 0
    int max = s;                // Start max at size
    int i = ( max - min ) / 2;  // Initial midpoint

    // WHILE WITHIN THE STRING
    while ( i >= min && i <= max )
    {
      // IF THE LETTER IS FOUND
      if ( a == y[ i ] )
      {
        break;  // Break out of the loop
      }

      // IF THE LETTER IS LESS
      else if ( a < y[ i ] )
      {
        max = i;                      // Update max
        i = i - ( i - min + 1 ) / 2;  // Move i left

        // IF AT THE BOTTOM
        if ( i < min || max == min )
        {
          i = s + 1 + ( min - 1 );  // Set to out to the left
          break;                    // Break out of the loop
        }
      }

      // IF THE LETTER IS MORE
      else
      {
        min = i;                      // Update min
        i = i + ( max - i + 1 ) / 2;  // Move i right

        // IF AT THE TOP
        if ( i >= max || max == min )
        {
          i = s + max;  // Set to out to the right
          break;        // Break out of the loop
        }
      }
    }

    return i;  // Return the index
  }

/*--------------------------------------------------------------------
 set - Set the tree to a Text and a value
--------------------------------------------------------------------*/

  public void set ( String z, String v )
  {
    set ( s2t ( z ).chr (), 0, s2t ( v ) );  // Run set with a character array
  }

/*--------------------------------------------------------------------
 set - Set the tree to a Text and a value
--------------------------------------------------------------------*/

  public void set ( Text z, String v )
  {
    set ( z.chr (), 0, s2t ( v ) );  // Run set with a character array
  }

/*--------------------------------------------------------------------
 set - Set the tree to a Text and a value
--------------------------------------------------------------------*/

  public void set ( Text z, Text v )
  {
    set ( z.chr (), 0, v );  // Run set with a character array
  }

/*--------------------------------------------------------------------
 set - Set the tree to a string and a value
--------------------------------------------------------------------*/

  public void set ( String z, Text v )
  {
    set ( s2c ( z ), 0, v );  // Run set with a character array
  }

/*--------------------------------------------------------------------
 set - Set the tree to a character and a value
--------------------------------------------------------------------*/

  public void set ( char z, Text v )
  {
    set ( h2c ( z ), 0, v );  // Run set with a character array
  }

/*--------------------------------------------------------------------
 set - Set the tree to a character array, an index, and a value
--------------------------------------------------------------------*/

  public void set ( char[] z, int ind, Text v )
  {
    // IF THE TREE IS EMPTY
    if ( t == null )
    {
      char[] r = sub ( z, ind, z.length );  // The remainder of the string

      // IF x IS NULL
      if ( x == null )
      {
        x = r;       // Set the lookup
        this.v = v;  // Set the value
      }

      // IF x IS EMPTY
      else if ( x.length == 0 )
      {
        // IF THE STRING IS EMPTY
        if ( z.length == ind )
        {
          x = r;       // Set the lookup
          this.v = v;  // Set the value
        }

        // IF THE STRING IS NOT EMPTY
        else
        {
          TextTree a = new TextTree ();     // Create the a tree
          a.set ( z, ind + 1, v );  // Set the a tree

          y = new char[1];      // Initialize y
          t = new TextTree[1];  // Initialize t
          s = 1;                // Initialize size to 1
          y[ 0 ] = z[ ind ];    // Set a's first character
          t[ 0 ] = a;           // Set a's tree
        }
      }

      // IF THE STRINGS ARE EQUAL
      else if ( eq ( x, r ) )
      {
        this.v = v;  // Set the value
      }

      // IF THE STRINGS ARE NOT EQUAL
      else
      {
        // IF THE STRING IS EMPTY
        if ( z.length == ind )
        {
          TextTree a = new TextTree ();    // Create the a tree
          a.set ( x, 1, this.v );  // Set the a tree

          y = new char[1];      // Initialize y
          t = new TextTree[1];  // Initialize t
          s = 1;                // Initialize size to 1
          y[ 0 ] = x[ 0 ];      // Set a's first character
          t[ 0 ] = a;           // Set a's tree
          this.v = v;           // Set the value
          x = new char[0];      // Zero out x   NEW
        }

        // IF THE STRING IS NOT EMPTY
        else
        {
          char a1 = r[ 0 ];  // First letter of a
          char b1 = x[ 0 ];  // First letter of b

          // IF A IS NOT EQUAL TO B
          if ( a1 != b1 )
          {
            TextTree a = new TextTree ();  // Create the a tree
            a.set ( r, 1, v );     // Set the a tree

            TextTree b = new TextTree ();    // Create the b tree
            b.set ( x, 1, this.v );  // Set the b tree

            y = new char[2];      // Initialize y
            t = new TextTree[2];  // Initialize t
            s = 2;                // Initialize size to 2

            // IF A IS LESS THAN B
            if ( a1 < b1 )
            {
              y[ 0 ] = a1;  // Set a's first character
              y[ 1 ] = b1;  // Set b's first character
              t[ 0 ] = a;   // Set a's tree
              t[ 1 ] = b;   // Set b's tree
            }

            // IF B IS LESS THAN A
            else
            {
              y[ 0 ] = b1;  // Set b's first character
              y[ 1 ] = a1;  // Set a's first character
              t[ 0 ] = b;   // Set b's tree
              t[ 1 ] = a;   // Set a's tree
            }
          }

          // IF THE CHARACTERS ARE THE SAME
          else
          {
            TextTree a = new TextTree ();  // Create the a tree
            a.set ( r, 1, v );             // Set the a tree
            a.set ( x, 1, this.v );        // Add b
            y = new char[1];               // Initialize y
            t = new TextTree[1];           // Initialize t
            s = 1;                         // Initialize size
            y[ 0 ] = a1;                   // Set the character
            t[ 0 ] = a;                    // Set the tree
          }

          this.x = new char[0];  // Ensure a blank x
          this.v = null;         // Null out v
        }
      }
    }

    // IF THE TREE EXISTS
    else
    {
      // IF THE STRING IS EMPTY
      if ( z.length == ind )
      {
        x = new char[0];  // Set x to an empty character array
        this.v = v;       // Set the value
      }

      // IF THE STRING IS NOT EMPTY
      else
      {
        int i = find ( z[ ind ] );  // Find the branch

        // IF IN THE TREE
        if ( i < s )
        {
          t[ i ].set ( z, ind + 1, v );  // Set the branch
        }

        // IF NOT IN THE TREE
        else
        {
          TextTree a = new TextTree ();  // Create the a tree
          a.set ( z, ind + 1, v );       // Set the a tree
          i = i - s;                     // Fix i

          // IF THE TREE IS FULL
          if ( t.length == s )
          {
            char[] ny = new char[t.length * 2];          // Double the size
            TextTree[] nt = new TextTree[t.length * 2];  // Double the size

            // LOOP THROUGH THE ELEMENTS
            for ( int j = 0; j < i; j++ )
            {
              ny[ j ] = y[ j ];  // Store the jth lookup item
              nt[ j ] = t[ j ];  // Store the jth tree item
            }

            // ADD THE VARIABLE
            ny[ i ] = z[ ind ];  // Store the ith lookup item
            nt[ i ] = a;         // Store the ith tree item

            // LOOP THROUGH THE ELEMENTS
            for ( int j = i; j < t.length; j++ )
            {
              ny[ j + 1 ] = y[ j ];  // Store the jth lookup item
              nt[ j + 1 ] = t[ j ];  // Store the jth tree item
            }

            y = ny;  // Update keys
            t = nt;  // Update trees

            s++;  // Increase the size
          }

          // IF THE TREE IS NOT FULL
          else
          {
            // LOOP THROUGH THE ELEMENTS
            for ( int j = s - 1; j >= i; j-- )
            {
              y[ j + 1 ] = y[ j ];  // Store the jth lookup item
              t[ j + 1 ] = t[ j ];  // Store the jth tree item
            }

            // ADD THE VARIABLE
            y[ i ] = z[ ind ];  // Store the ith lookup item
            t[ i ] = a;         // Store the ith tree item

            s++;  // Increase the size
          }
        }
      }
    }
  }

/*--------------------------------------------------------------------
 del - Delete the tree to a string and a value
--------------------------------------------------------------------*/

  public void del ( String z )
  {
    del ( s2c ( z ), 0 );  // Run del with a character array
  }

/*--------------------------------------------------------------------
 del - Delete the tree to a character and a value
--------------------------------------------------------------------*/

  public void del ( char z )
  {
    del ( h2c ( z ), 0 );  // Run del with a character array
  }

/*--------------------------------------------------------------------
 del - Delete the tree to a character array, an index, and a value
--------------------------------------------------------------------*/

  public boolean del ( char[] z, int ind )
  {
    // IF THE TREE IS EMPTY
    if ( t == null )
    {
      // if strings match, then empty value and string
    }

    // IF THE TREE EXISTS
    else
    {
      // IF THE STRING IS EMPTY
      if ( z.length == ind )
      {
        // zero out the value, if necessary move a branch to here
        // actually, not checking for 0, checking for 1
      }

      // IF THE STRING IS NOT EMPTY
      else
      {
        int i = find ( z[ ind ] );  // Find the branch

        // IF IN THE TREE
        if ( i < s )
        {
          // IF THE SUBSTRINGS MATCH UP
          if ( eq ( z, ind + 1, t[ i ].x, 0 ) )
          {
            t[ i ].x = null;  // Unset the string
            t[ i ].v = null;  // Null out the value

            // IF IT DOES NOT HAVE A TREE OF ITS OWN
            if ( t[ i ].t == null )
            {
              // REMOVE THE BRANCH
              for ( int j = i; j < s - 1; j++ )
              {
                y[ j ] = y[ j + 1 ];  // Save the jth item
                t[ j ] = t[ j + 1 ];  // Save the jth item
              }

              s--;  // Decrement size
            }

            // IF SIZE IS 1
            if ( s == 1 )
            {
              // IF THE LAST BRANCH HAS NO TREE
              if ( t[ 0 ].t == null && x == null )
              {
                String tz = y[ 0 ] + c2s ( x );  // Find the new string
                this.x = s2c ( tz );             // Set the string
                this.v = t[ 0 ].v;               // Set the value
                y = null;                        // Remove the lookups
                t = null;                        // Remove the branch array
                s = 0;                           // Set size to zero
                return true;                     // Possible remove this tree
              }
            }
          }

          // IF THE SUBSTRINGS DO NOT MATCH UP
          else
          {
            t[ i ].del ( z, ind + 1 );  // Delete the branch?
          }
        }
      }
    }

    return false;  // No need to delete anything more
  }

/*--------------------------------------------------------------------
 add - Add two Text arrays
--------------------------------------------------------------------*/

  public void add ( Text[] a, Text[] b )
  {
    for ( int i = 0; i < a.length; i++ )
    {
      add ( a[ i ], " " + b[ i ] );
    }
  }

/*--------------------------------------------------------------------
 add - Add the text to a spot in the tree
--------------------------------------------------------------------*/

  public void add ( String z, Text v )
  {
    add ( s2t ( z ), v );
  }

/*--------------------------------------------------------------------
 add - Add the text to a spot in the tree
--------------------------------------------------------------------*/

  public void add ( Text z, String v )
  {
    add ( z, s2t ( v ) );
  }

/*--------------------------------------------------------------------
 add - Add the text to a spot in the tree
--------------------------------------------------------------------*/

  public void add ( String z, String v )
  {
    add ( s2t ( z ), s2t ( v ) );
  }

/*--------------------------------------------------------------------
 add - Add the text to a spot in the tree
--------------------------------------------------------------------*/

  public void add ( Text z, Text v )
  {
    Text a = get ( z );  // Get the text

    // IF THE TEXT DOES NOT EXIST
    if ( a == null )
    {
      a = new Text ( v );  // Make a new text
      set ( z, v );        // Set the text
    }

    // IF THE TEXT DOES EXIST
    else
    {
      a.rpush ( v );  // Push the new text
    }
  }

/*--------------------------------------------------------------------
 map - Load a map
--------------------------------------------------------------------*/

  public void map ( String fname )
  {
    Infile d = new Infile ( fname );
    d.readLine ();
    while ( !d.eof )
    {
      Text[] y = d.splitTab ();
      set ( y[ 0 ], y[ 1 ] );
      d.readLine ();
    }

    d.close ();
  }

/*--------------------------------------------------------------------
 mapMult - Load a multiple map
--------------------------------------------------------------------*/

  public void mapMult ( String fname )
  {
    Infile d = new Infile ( fname );
    d.readLine ();
    while ( !d.eof )
    {
      Text[] y = d.splitTab ();
      add ( y[ 0 ], "|" + y[ 1 ] );
      d.readLine ();
    }

    d.close ();
  }
}