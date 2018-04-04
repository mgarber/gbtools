package broad.projection.utils;/*--------------------------------------------------------------------
 com.scanfeld.io.web.Url
 ---------------------------------------------------------------------
 Models a url object
 ---------------------------------------------------------------------
 Author: Daniel Scanfeld
   Date: 08.23.06
 ---------------------------------------------------------------------
--------------------------------------------------------------------*/

import java.io.InputStream;
import java.net.URL;
import java.net.URLConnection;



public class Url extends Text
{
/*--------------------------------------------------------------------
 Class Variables
--------------------------------------------------------------------*/

  public int bufSize;  // Size of the buffer

/*--------------------------------------------------------------------
 CONSTRUCTOR FUNCTIONS
 ---------------------------------------------------------------------
 Constructors for the possible text objects
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 com.scanfeld.io.web.Url - Main constructor
--------------------------------------------------------------------*/

  public Url ()
  {
    init ();         // Initialize the com.scanfeld.core.Text object
    bufSize = 1024;  // Buffer size
  }

/*--------------------------------------------------------------------
 com.scanfeld.io.web.Url - Load a website
--------------------------------------------------------------------*/

  public Url ( String url )
  {
    init ();         // Initialize the com.scanfeld.core.Text object
    bufSize = 1024;  // Buffer size
    load ( url );    // Load the url
  }

/*--------------------------------------------------------------------
 com.scanfeld.io.web.Url - Load a website
--------------------------------------------------------------------*/

  public Url ( String url, String fname )
  {
    init ();         // Initialize the com.scanfeld.core.Text object
    bufSize = 1024;  // Buffer size
    load ( url, fname );    // Load the url
  }

/*--------------------------------------------------------------------
 load - Load the full website
--------------------------------------------------------------------*/

  public void load ( String url )
  {
    // TRY TO LOAD THE WEBSITE
    try
    {
      URL u = new URL ( url );                                        // Create a url object
      URLConnection connection = u.openConnection ();                 // Create a url connection
      connection.setRequestProperty ( "User-agent", "Mozilla/4.0" );  // Allow google searching
      InputStream urlData = connection.getInputStream ();             // Create an input stream to read the data
      byte[] buffer = new byte[bufSize];                            // Create a byte buffer for the input
      int numRead;                                                    // Track the number of bytes read

      // LOOP UNTIL THE END OF THE FILE REACHED
      while ( ( numRead = urlData.read ( buffer ) ) != -1 )
      {
        rpush ( b2c ( buffer, numRead ) );  // Push the buffer
      }

      urlData.close ();  // Close the input stream
    }

    // IF ANY ERROR HAS OCCURRED
    catch ( Exception e )
    {
      err ( e.getMessage () );  // Output the error message
    }

    dbg ( "LOAD \"" + url + "\"" );  // Status message
  }

/*--------------------------------------------------------------------
 load - Load the full website into a file
--------------------------------------------------------------------*/

  public void load ( String url, String fname )
  {
    Outfile a = new Outfile ( fname );  // Create the outfile

    // TRY TO LOAD THE WEBSITE
    try
    {
      URL u = new URL ( url );                                        // Create a url object
      URLConnection connection = u.openConnection ();                 // Create a url connection
      connection.setRequestProperty ( "User-agent", "Mozilla/4.0" );  // Allow google searching
      InputStream urlData = connection.getInputStream ();             // Create an input stream to read the data
      byte[] buffer = new byte[bufSize];                            // Create a byte buffer for the input
      int numRead;                                                    // Track the number of bytes read

      // LOOP UNTIL THE END OF THE FILE REACHED
      while ( ( numRead = urlData.read ( buffer ) ) != -1 )
      {
        a.fout.write ( buffer, 0, numRead );  // Output to the file
      }

      urlData.close ();                       // Close the input stream
    }

    // IF ANY ERROR HAS OCCURRED
    catch ( Exception e )
    {
      err ( e.getMessage () );  // Output the error message
    }

    a.close ();                                            // Close the outfile
    dbg ( "SAVE \"" + url + " \" TO \"" + fname + "\"" );  // Status message
  }
}
