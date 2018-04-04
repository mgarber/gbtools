package broad.projection.utils;

/*--------------------------------------------------------------------
 com.scanfeld.io.Outfile
 ---------------------------------------------------------------------
 Models an outfile file
 ---------------------------------------------------------------------
 Author: Daniel Scanfeld
   Date: 08.23.06
 ---------------------------------------------------------------------

--------------------------------------------------------------------*/

import java.io.File;
import java.io.FileOutputStream;


public class Outfile extends Text
{
/*--------------------------------------------------------------------
 Class Variables
--------------------------------------------------------------------*/

  public File f;                 // Main file object
  public FileOutputStream fout;  // Stream to write characters

/*--------------------------------------------------------------------
 CONSTRUCTOR FUNCTIONS
 ---------------------------------------------------------------------
 Constructors for the possible infile objects
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 com.scanfeld.io.Outfile - Create an outfile
 --------------------------------------------------------------------*/

  public Outfile ()
  {
  }

/*--------------------------------------------------------------------
 com.scanfeld.io.Outfile - Load a file fname
--------------------------------------------------------------------*/

  public Outfile ( String fname )
  {
    open ( fname );  // Open fname
  }

/*--------------------------------------------------------------------
 FILE MANIPULATION FUNCTIONS
 ---------------------------------------------------------------------
 Provide file write functionality
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 save - Write the text to a file
--------------------------------------------------------------------*/

  public void save ( String fname )
  {
    // TRY TO SAVE THE FILE
    try
    {
      File f = new File ( fname );                         // Create  a file object
      FileOutputStream fout = new FileOutputStream ( f );  // Create a file output stream
      fout.write ( bytes () );                             // Write the bytes
      fout.close ();                                       // Close the output stream
    }

    // IF ANY ERROR HAS OCCURRED
    catch ( Exception e )
    {
      err ( e.getMessage () );  // Output the error message
    }

    dbg ( "SAVE \"" + fname + "\"" );  // Status message
    created.rpush( fname );
  }

/*--------------------------------------------------------------------
 open - Open a file for reading
--------------------------------------------------------------------*/

  public void open ( String fname )
  {
    // TRY TO OPEN THE FILE
    try
    {
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
          created.rpush( dir );
        }
      }

      f = new File ( fname );             // Create  a file object
      fout = new FileOutputStream ( f );  // Create a file output stream
    }

    // IF ANY ERROR HAS OCCURRED
    catch ( Exception e )
    {
      err ( "Error opening: \"" + fname + "\" for output." );  // Output the error message
    }

    dbg ( "SAVE \"" + fname + "\"" );  // Status message
    created.rpush(fname);
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
        // IF THE FILE OUTPUT STREAM IS NOT NULL
        if ( fout != null )
        {
          fout.close ();  // Close the file input stream
          fout = null;    // Set the stream to null
        }

        f = null;    // Set the file to null
      }

      // IF AN ERROR OCCURRED
      catch ( Exception e )
      {
        err ( "Error closing output file." );  // Output an error message
      }
    }
  }

/*--------------------------------------------------------------------
 writeLine - Write a line
--------------------------------------------------------------------*/

  public void writeLine ()
  {
    write ( "\n" );  // Write a new line character
  }

/*--------------------------------------------------------------------
 writeLine - Write a line
--------------------------------------------------------------------*/

  public void writeLine ( String z )
  {
    write ( z + "\n" );  // Write the string plus a new line character
  }

/*--------------------------------------------------------------------
 writeLine - Write a line
--------------------------------------------------------------------*/

  public void writeLine ( Text z )
  {
    write ( z );     // Write the com.scanfeld.core.Text
    write ( "\n" );  // Write the new line
  }

/*--------------------------------------------------------------------
 writeLine - Write a line
--------------------------------------------------------------------*/

  public void writeLine ( byte[] z )
  {
    write ( z );     // Write the byte array
    write ( "\n" );  // Write the new line
  }

/*--------------------------------------------------------------------
 write - Write a string
--------------------------------------------------------------------*/

  public void write ( String z )
  {
    write ( z.getBytes () );  // Write the bytes of the string
  }

/*--------------------------------------------------------------------
 write - Write a com.scanfeld.core.Text
--------------------------------------------------------------------*/

  public void write ( Text z )
  {
    write ( z.bytes () );  // Write the bytes of the com.scanfeld.core.Text
  }

/*--------------------------------------------------------------------
 write - Write a com.scanfeld.core.Texts
--------------------------------------------------------------------*/

  public void write ( Texts z )
  {
    write ( z.bytes () );  // Write the bytes of the com.scanfeld.core.Texts
  }

/*--------------------------------------------------------------------
 write - Write a byte array
--------------------------------------------------------------------*/

  public void write ( byte b[] )
  {
    // IF THE FILE AND STREAM EXISTS
    if ( f != null && fout != null )
    {
      // TRY TO WRITE THE BYTES
      try
      {
        fout.write ( b );  // Write the bytes
      }

      // IF AN ERROR OCCURRED
      catch ( Exception e )
      {
        err ( "Error outputting to file." );  // Output an error message
      }
    }
  }

/*--------------------------------------------------------------------
 out - Write a line
--------------------------------------------------------------------*/

  public void out ()
  {
    writeLine ( "" );
  }

/*--------------------------------------------------------------------
 out - Write a line
--------------------------------------------------------------------*/

  public void out ( byte[] z )
  {
    writeLine ( z );
  }

/*--------------------------------------------------------------------
 out - Write a line
--------------------------------------------------------------------*/

  public void out ( String z )
  {
    writeLine ( z );
  }

/*--------------------------------------------------------------------
 out - Write a line
--------------------------------------------------------------------*/

  public void out ( Text z )
  {
    writeLine ( z );
  }

/*--------------------------------------------------------------------
 dbg - Output a debug string
--------------------------------------------------------------------*/

  public void dbg ( String z )
  {
    // IF DEBUGGING
    if ( dbgFlag )
    {
      System.out.println ( z );  // Write to the system output function
    }
  }

/*--------------------------------------------------------------------
 err - Output an error string
--------------------------------------------------------------------*/

  public void err ( String z )
  {
    System.out.println ( z );  // Write to the system output function
    System.exit ( 0 );         // Exit the program
  }
}
