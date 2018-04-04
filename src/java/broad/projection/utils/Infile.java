package broad.projection.utils;

/*--------------------------------------------------------------------
 com.scanfeld.io.Infile
 ---------------------------------------------------------------------
 Models an input file
 ---------------------------------------------------------------------
 Author: Daniel Scanfeld
   Date: 08.23.06
 ---------------------------------------------------------------------

--------------------------------------------------------------------*/

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;



public class Infile extends Text
{
/*--------------------------------------------------------------------
 Class Variables
--------------------------------------------------------------------*/

  public File f;               // Main file object
  public FileInputStream fin;  // Stream to read in characters
  public BufferedReader br;    // Buffer for the characters
  public boolean eof;          // End of file flag
  public int bufSize;          // Size of the buffer

/*--------------------------------------------------------------------
 CONSTRUCTOR FUNCTIONS
 ---------------------------------------------------------------------
 Constructors for the possible infile objects
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 Infile - Create an infile
 --------------------------------------------------------------------*/

  public Infile ()
  {
    bufSize = 1024;  // Buffer size
  }

/*--------------------------------------------------------------------
 Infile - Load a file fname
--------------------------------------------------------------------*/

  public Infile ( String fname )
  {
    bufSize = 1024;      // Buffer size
    open ( fname );  // Open fname
  }

/*--------------------------------------------------------------------
 FILE MANIPULATION FUNCTIONS
 ---------------------------------------------------------------------
 Provide file read functionality
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 load - Load the entire contents of a file
--------------------------------------------------------------------*/

  public void load ( String fname )
  {
    // TRY TO LOAD THE FILE
    try
    {
      f = new File ( fname );               // Create  a file object
      fin = new FileInputStream ( f );      // Create a file input stream
      byte[] buffer = new byte[bufSize];  // Create a byte buffer for the input
      int numRead;                          // Track the number of bytes read

      // LOOP UNTIL THE END OF THE FILE REACHED
      while ( ( numRead = fin.read ( buffer ) ) != -1 )
      {
        rpush ( b2c ( buffer, numRead ) );  // Push the buffer
      }

      fin.close ();  // Close the input stream
    }

    // IF ANY ERROR HAS OCCURRED
    catch ( Exception e )
    {
      err ( e.getMessage () );  // Output the error message
    }

    Text a = replace ( "\r", "" );

    for ( int i = a.size - 1; i >= 0; i-- )
    {
      if ( a.get ( i ) != '\n' )
      {
        if ( i != ( a.size - 1 ) )
        {
          a = a.sub ( 0, i + 1 );
        }
        
        break;
      }
    }

    set ( a );

    dbg ( "LOAD \"" + fname + "\"" );  // Status message
  }

/*--------------------------------------------------------------------
 open - Open a file for reading
--------------------------------------------------------------------*/

  public void open ( String fname )
  {
    // TRY TO OPEN THE FILE
    try
    {
      f = new File ( fname );                                     // Create  a file object
      fin = new FileInputStream ( f );                            // Create a file input stream
      br = new BufferedReader ( new InputStreamReader ( fin ) );  // Create a buffered reader
      eof = false;                                                // Unset the end of file flag
    }

    // IF ANY ERROR HAS OCCURRED
    catch ( Exception e )
    {
      err ( "Error opening: \"" + fname + "\" for input." );  // Output the error message
    }

    dbg ( "LOAD \"" + fname + "\"" );  // Status message
  }

/*--------------------------------------------------------------------
 close - Close a file for reading
--------------------------------------------------------------------*/

  public void close ()
  {
    // IF THE FILE OBJECT IS NOT NULL
    if ( f != null )
    {
      // TRY TO CLOSE THE FILE
      try
      {
        // IF THE BUFFERED READER IS NOT NULL
        if ( br != null )
        {
          br.close ();  // Close the buffered reader
          br = null;    // Set the reader to null
        }

        // IF THE FILE INPUT STREAM IS NOT NULL
        if ( fin != null )
        {
          fin.close ();  // Close the file input stream
          fin = null;    // Set the stream to null
        }

        f = null;    // Set the file to null
        eof = true;  // Set the end of file flag
      }

      // IF AN ERROR OCCURRED
      catch ( Exception e )
      {
        err ( "Error closing input file." );  // Output an error message
      }
    }
  }

/*--------------------------------------------------------------------
 readLine - Read in a line
--------------------------------------------------------------------*/

  public void readLine ()
  {
    String z = null;  // Create a temporary string

    // IF THE FILE AND BUFFER EXISTS AND NOT AT THE END
    if ( f != null && br != null && !eof )
    {
      // TRY TO READ A LINE
      try
      {
        z = br.readLine ();  // Save the string

        // IF NOT AT THE END OF THE FILE
        if ( z != null )
        {
          set ( z );  // Set to the line
        }

        // IF AT THE END OF THE FILE
        else
        {
          set ( "" );     // Set to the empty string
          eof = true;  // Set the end of file string
        }
      }

      // IF AN ERROR OCCURRED
      catch ( Exception e )
      {
        err ( "Error inputting from file." );  // Output an error message
      }
    }
  }

/*--------------------------------------------------------------------
 in - Read in a line
--------------------------------------------------------------------*/

  public void in ()
  {
    readLine ();
  }

/*--------------------------------------------------------------------
 in - Read in a line
--------------------------------------------------------------------*/

  public void in ( int n )
  {
    skip ( n );
  }

/*--------------------------------------------------------------------
 skip - Skip a number of lines
--------------------------------------------------------------------*/

  public void skip ( int n )
  {
    // LOOP THROUGH THE NUMBER OF LINES
    for ( int i = 0; i < n; i++ )
    {
      readLine ();  // Skip the line
    }
  }
}