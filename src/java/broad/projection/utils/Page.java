package broad.projection.utils;



/*--------------------------------------------------------------------
 com.scanfeld.io.web.Page
 ---------------------------------------------------------------------
 Models a webpage
 ---------------------------------------------------------------------
 Author: Daniel Scanfeld
   Date: 08.23.06
 ---------------------------------------------------------------------

--------------------------------------------------------------------*/

public class Page extends Text
{
/*--------------------------------------------------------------------
 Class Variables
--------------------------------------------------------------------*/

  public String title;  // Title of the webpage
  public Text br;       // Pagebreak

/*--------------------------------------------------------------------
 CONSTRUCTOR FUNCTIONS
 ---------------------------------------------------------------------
 Constructors for the possible page objects
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 com.scanfeld.io.web.Page - Null constructor
--------------------------------------------------------------------*/

  public Page ()
  {
    br = s2t ( "<br>" );
    set ( htm ( "com.scanfeld.pr.Test", "CCCCFF", cen ( "dan" ) ) );
  }

/*--------------------------------------------------------------------
 save - Save the web page to a directory and a file name
--------------------------------------------------------------------*/

  public void save ( String dir, String name )
  {
    copy ( "inc/s.gif", dir + "s.gif" );  // Copy the spacer gif to this directory
    copy ( "inc/s.css", dir + "s.css" );  // Copy the style sheet to this directory
    copy ( "inc/i.js", dir + "i.js" );    // Copy the javascript include file to this directory

    Outfile a = new Outfile ( dir + name );  // Create an outfile for the page
    a.write ( this );                        // Write the page
    a.close ();                              // Close the file
  }

/*--------------------------------------------------------------------
 cen - Center
--------------------------------------------------------------------*/

  public Text cen ( String z )
  {
    return cen ( s2t ( z ) );  // Call the com.scanfeld.core.Text version
  }

/*--------------------------------------------------------------------
 cen - Center
--------------------------------------------------------------------*/

  public Text cen ( Text y )
  {
    Text z = y.copy ();
    z.lpush ( "<center>" );  // Push the opening tag
    z.rpush ( "</center>" ); // Push the closing tag
    return z;                // Return the new com.scanfeld.core.Text
  }

/*--------------------------------------------------------------------
 nb - No line break
--------------------------------------------------------------------*/

  public Text nb ( String z )
  {
    return nb ( s2t ( z ) );  // Call the com.scanfeld.core.Text version
  }

/*--------------------------------------------------------------------
 nb - No line break
--------------------------------------------------------------------*/

  public Text nb ( Text y )
  {
    Text z = y.copy ();
    z.lpush ( "<nobr>" );  // Push the opening tag
    z.rpush ( "</nobr>" ); // Push the closing tag
    return z;              // Return the new com.scanfeld.core.Text
  }

/*--------------------------------------------------------------------
 tab - Default table function
--------------------------------------------------------------------*/

  public Text tab ( String z )
  {
    return tab ( s2t ( z ) );  // Call the com.scanfeld.core.Text version
  }

/*--------------------------------------------------------------------
 tab - Default table function
--------------------------------------------------------------------*/

  public Text tab ( Text y )
  {
    Text z = y.copy ();
    z.lpush ( "<table cellspacing=\"0\" cellpadding=\"0\" border=\"0\">" );  // Push the opening tag
    z.rpush ( "</table>" );                                                  // Push the closing tag
    return z;                                                                // Return the new com.scanfeld.core.Text
  }

/*--------------------------------------------------------------------
 tabw - Table with width
--------------------------------------------------------------------*/

  public Text tabw ( int w, String z )
  {
    return tabw ( i2t ( w ), s2t ( z ) );  // Call the com.scanfeld.core.Text version
  }

/*--------------------------------------------------------------------
 tabw - Table with width
--------------------------------------------------------------------*/

  public Text tabw ( Text w, Text y )
  {
    Text z = y.copy ();
    z.lpush ( "<table cellspacing=\"0\" cellpadding=\"0\" border=\"0\" width=\"" + w + "\">" );  // Push the opening tag
    z.rpush ( "</table>" );                                                                      // Push the closing tag
    return z;                                                                                    // Return the new com.scanfeld.core.Text
  }

/*--------------------------------------------------------------------
 tabp - Table with padding
--------------------------------------------------------------------*/

  public Text tabp ( int p, String z )
  {
    return tabp ( i2t ( p ), s2t ( z ) );  // Call the com.scanfeld.core.Text version
  }

/*--------------------------------------------------------------------
 tabp - Table with padding
--------------------------------------------------------------------*/

  public Text tabp ( Text p, Text y )
  {
    Text z = y.copy ();
    z.lpush ( "<table cellspacing=\"0\" cellpadding=\"" + p + "\" border=\"0\">" );  // Push the opening tag
    z.rpush ( "</table>" );                                                          // Push the closing tag
    return z;                                                                        // Return the new com.scanfeld.core.Text
  }

/*--------------------------------------------------------------------
 tabwp - Table with width and padding
--------------------------------------------------------------------*/

  public Text tabwp ( int w, int p, String z )
  {
    return tabwp ( i2t ( w ), i2t ( p ), s2t ( z ) );  // Call the com.scanfeld.core.Text version
  }

/*--------------------------------------------------------------------
 tabwp - Table with width and padding
--------------------------------------------------------------------*/

  public Text tabwp ( Text w, Text p, Text y )
  {
    Text z = y.copy ();
    z.lpush ( "<table cellspacing=\"0\" cellpadding=\"" + p + "\" border=\"0\" width=\"" + w + "\">" );  // Push the opening tag
    z.rpush ( "</table>" );                                                                              // Push the closing tag
    return z;                                                                                            // Return the new com.scanfeld.core.Text
  }

/*--------------------------------------------------------------------
 tr - Table row
--------------------------------------------------------------------*/

  public Text tr ( String z )
  {
    return tr ( s2t ( z ) );  // Call the com.scanfeld.core.Text version
  }

/*--------------------------------------------------------------------
 tr - Table row
--------------------------------------------------------------------*/

  public Text tr ( Text y )
  {
    Text z = y.copy ();
    z.lpush ( "<tr>" );  // Push the opening tag
    z.rpush ( "</tr>" ); // Push the closing tag
    return z;            // Return the new com.scanfeld.core.Text
  }

/*--------------------------------------------------------------------
 td - Table cell
--------------------------------------------------------------------*/

  public Text td ( String z )
  {
    return td ( s2t ( z ) );  // Call the com.scanfeld.core.Text version
  }

/*--------------------------------------------------------------------
 td - Table cell
--------------------------------------------------------------------*/

  public Text td ( Text y )
  {
    Text z = y.copy ();
    z.lpush ( "<td>" );  // Push the opening tag
    z.rpush ( "</td>" ); // Push the closing tag
    return z;            // Return the new com.scanfeld.core.Text
  }

/*--------------------------------------------------------------------
 tdt - Table cell aligned to the top
--------------------------------------------------------------------*/

  public Text tdt ( String z )
  {
    return tdt ( s2t ( z ) );  // Call the com.scanfeld.core.Text version
  }

/*--------------------------------------------------------------------
 tdt - Table cell aligned to the top
--------------------------------------------------------------------*/

  public Text tdt ( Text y )
  {
    Text z = y.copy ();
    z.lpush ( "<td valign=\"top\">" );  // Push the opening tag
    z.rpush ( "</td>" );                // Push the closing tag
    return z;                           // Return the new com.scanfeld.core.Text
  }

/*--------------------------------------------------------------------
 tdl - Table cell with colspan
--------------------------------------------------------------------*/

  public Text tdl ( int n, String z )
  {
    return tdl ( i2t ( n ), s2t ( z ) );  // Call the com.scanfeld.core.Text version
  }

/*--------------------------------------------------------------------
 tdl - Table cell with colspan
--------------------------------------------------------------------*/

  public Text tdl ( Text n, Text y )
  {
    Text z = y.copy ();
    z.lpush ( "<td colspan=" + n + ">" );  // Push the opening tag
    z.rpush ( "</td>" );                   // Push the closing tag
    return z;                              // Return the new com.scanfeld.core.Text
  }

/*--------------------------------------------------------------------
 tds - Table cell with a spacer image
--------------------------------------------------------------------*/

  public Text tds ( int c, int w, int h )
  {
    return tds ( i2t ( c ), i2t ( w ), i2t ( h ) );  // Call the com.scanfeld.core.Text version
  }

/*--------------------------------------------------------------------
 tds - Table cell with a spacer image
--------------------------------------------------------------------*/

  public Text tds ( Text c, Text w, Text h )
  {
    return new Text ( "<td bgcolor=\"" + c + "\" width=\"" + w + "\">" + sp ( w, h ) + "</td>" );  // Build the tag
  }

/*--------------------------------------------------------------------
 tdc - Table cell with color
--------------------------------------------------------------------*/

  public Text tdc ( String c, String z )
  {
    return tdc ( s2t ( c ), s2t ( z ) );  // Call the com.scanfeld.core.Text version
  }

/*--------------------------------------------------------------------
 tdc - Table cell with color
--------------------------------------------------------------------*/

  public Text tdc ( Text c, Text y )
  {
    Text z = y.copy ();
    z.lpush ( "<td bgcolor=\"" + c + "\">" );  // Push the opening tag
    z.rpush ( "</td>" );                       // Push the closing tag
    return z;                                  // Return the new com.scanfeld.core.Text
  }

/*--------------------------------------------------------------------
 sp - Spacer image
--------------------------------------------------------------------*/

  public Text sp ( int w, int h )
  {
    return sp ( i2t ( w ), i2t ( h ) );  // Call the com.scanfeld.core.Text version
  }

/*--------------------------------------------------------------------
 sp - Spacer image
--------------------------------------------------------------------*/

  public Text sp ( Text w, Text h )
  {
    return new Text ( "<img src=\"s.gif\" width=\"" + w + "\" height=\"" + h + "\">" );  // Build the tag
  }

/*--------------------------------------------------------------------
 brd - Bordered table
--------------------------------------------------------------------*/

  public Text brd ( Text w, Text bgc, Text bw, Text brc, Text pad, Text y )
  {
    Text rw = s2t ( "100%" );
    Text iw = s2t ( "100%" );

    if ( w.ne ( "100%" ) )
    {
      iw.set ( w.num () - pad.num () * 2 );
      rw.set ( iw.num () - pad.num () * 2 );
    }

    Text rowa = tr ( c ( tds ( brc, bw, bw ), tds ( brc, iw, bw ), tds ( brc, bw, bw ) ) );
    Text rowb = tr ( c ( tds ( brc, bw, bw ), tdc ( bgc, tabwp ( iw, pad, tdc ( bgc, y ) ) ), tds ( brc, bw, bw ) ) );

    return tabw ( w, c ( rowa, rowb, rowa ) );
  }

/*--------------------------------------------------------------------
 brdt - Bordered table top
--------------------------------------------------------------------*/

  public Text brdt ( Text bgc, Text bw, Text brc, Text pad, Text y )
  {
    Text rowa = tr ( c ( tds ( brc, bw, bw ), tds ( brc, i2t ( 1 ), bw ), tds ( brc, bw, bw ) ) );
    Text rowb = tr ( c ( tds ( brc, bw, bw ), tdc ( bgc, tabp ( pad, tdc ( bgc, y ) ) ), tds ( brc, bw, bw ) ) );

    return tab ( c ( rowa, rowb ) );
  }

/*--------------------------------------------------------------------
 pan - Panel
--------------------------------------------------------------------*/

  public Text pan ( Text w, Text tbc, Text bgc, Text bw, Text brc, Text pad, Text x, Text y )
  {
    return tab ( tr ( td ( c ( brdt ( tbc, bw, brc, pad, x ), brd ( w, bgc, bw, brc, pad, y ) ) ) ) );
  }

/*--------------------------------------------------------------------
 pan - Panel
--------------------------------------------------------------------*/

  public Text pand ( Text w, Text tbc, Text bgc, Text pdc, Text bw, Text brc, Text ibc, Text pad, Text pad2, Text x, Text y )
  {
    Text a = brd ( i2t ( 10 ), tbc, bw, ibc, pad2, x );
    Text b = brd ( i2t ( w.num () - pad.num () * 2 - bw.num () * 2 ), bgc, bw, ibc, pad2, y );

    return pan ( w, pdc, pdc, bw, brc, pad, a, b );
  }

/*--------------------------------------------------------------------
 pane - Panel
--------------------------------------------------------------------*/

  public Text pane ( Text w, Text tbc, Text bgc, Text pdc, Text x, Text y )
  {
    return pand ( w, tbc, bgc, pdc, i2t ( 1 ), s2t ( "999999" ), s2t ( "CCCCCC" ), i2t ( 2 ), i2t ( 2 ), nb ( txt ( "pt", x ) ), y );
  }

/*--------------------------------------------------------------------
 txt - com.scanfeld.core.Text
--------------------------------------------------------------------*/

  public Text txt ( Text n, Text y )
  {
    Text z = y.copy ();
    z.lpush ( "<font class=\"" + n + "\">" );  // Push the opening tag
    z.rpush ( "</font>" );                     // Push the closing tag
    return z;                                  // Return the new com.scanfeld.core.Text
  }

/*--------------------------------------------------------------------
 txt - com.scanfeld.core.Text
--------------------------------------------------------------------*/

  public Text txt ( String n, Text z )
  {
    return txt ( s2t ( n ), z );
  }

/*--------------------------------------------------------------------
 tit - Title
--------------------------------------------------------------------*/

  public Text tit ( String z )
  {
    return tit ( s2t ( z ) );  // Call the com.scanfeld.core.Text version
  }

/*--------------------------------------------------------------------
 tit - Title
--------------------------------------------------------------------*/

  public Text tit ( Text y )
  {
    Text z = y.copy ();
    z.lpush ( "<title>" );  // Push the opening tag
    z.rpush ( "</title>" ); // Push the closing tag
    return z;               // Return the new com.scanfeld.core.Text
  }

/*--------------------------------------------------------------------
 htm - HTML com.scanfeld.io.web.Page
--------------------------------------------------------------------*/

  public Text htm ( String title, String bgc, String z )
  {
    return htm ( s2t ( title ), s2t ( bgc ), s2t ( z ) );  // Call the com.scanfeld.core.Text version
  }

/*--------------------------------------------------------------------
 htm - HTML com.scanfeld.io.web.Page
--------------------------------------------------------------------*/

  public Text htm ( String title, String bgc, Text z )
  {
    return htm ( s2t ( title ), s2t ( bgc ), z );  // Call the com.scanfeld.core.Text version
  }

/*--------------------------------------------------------------------
 htm - HTML com.scanfeld.io.web.Page
--------------------------------------------------------------------*/

  public Text htm ( Text title, Text bgc, Text y )
  {
    Text z = y.copy ().trimWhite ().indent ();
    z = n ( hea ( n ( tit ( title ), css ( "s.css" ), scr ( "i.js" ) ) ), bdy ( bgc, z ) );
    z.lpush ( "<html>\n" );   // Push the opening tag
    z.rpush ( "\n</html>" );  // Push the closing tag
    return z;                 // Return the new com.scanfeld.core.Text
  }

/*--------------------------------------------------------------------
 bdy - Body
--------------------------------------------------------------------*/

  public Text bdy ( String c, String z )
  {
    return bdy ( s2t ( c ), s2t ( z ) );  // Call the com.scanfeld.core.Text version
  }

/*--------------------------------------------------------------------
 bdy - Body
--------------------------------------------------------------------*/

  public Text bdy ( String c, Text z )
  {
    return bdy ( s2t ( c ), z );  // Call the com.scanfeld.core.Text version
  }

/*--------------------------------------------------------------------
 bdy - Body
--------------------------------------------------------------------*/

  public Text bdy ( Text c, Text y )
  {
    Text z = new Text ( "" );
    z.rpush ( "<body bgcolor=\"" );
    z.rpush ( c );
    z.rpush ( "\" marginwidth=\"0\" marginheight=\"0\" leftmargin=\"0\" topmargin=\"0\" onload=\"init();\">\n" );   // Push the opening tag
    z.rpush ( c ( sp ( 1, 10 ), br ) );
    z.rpush ( "\n" );
    z.rpush ( y.trimWhite () );
    z.rpush ( "\n</body>" );  // Push the closing tag
    return z;                 // Return the new com.scanfeld.core.Text
  }

/*--------------------------------------------------------------------
 hea - Head
--------------------------------------------------------------------*/

  public Text hea ( String z )
  {
    return hea ( s2t ( z ) );  // Call the com.scanfeld.core.Text version
  }

/*--------------------------------------------------------------------
 hea - Head
--------------------------------------------------------------------*/

  public Text hea ( Text y )
  {
    Text z = y.copy ().trimWhite ().indent ();
    z.lpush ( "<head>\n" );   // Push the opening tag
    z.rpush ( "\n</head>" );  // Push the closing tag
    return z;                 // Return the new com.scanfeld.core.Text
  }

/*--------------------------------------------------------------------
 scr - Javascript
--------------------------------------------------------------------*/

  public Text scr ( String z )
  {
    return scr ( s2t ( z ) );  // Call the com.scanfeld.core.Text version
  }

/*--------------------------------------------------------------------
 scr - Javascript
--------------------------------------------------------------------*/

  public Text scr ( Text y )
  {
    Text z = y.copy ();
    z.lpush ( "<script language=\"javascript\" src=\"" );  // Push the opening tag
    z.rpush ( "\">" );                                     // Push the closing tag
    return z;                                              // Return the new com.scanfeld.core.Text
  }

/*--------------------------------------------------------------------
 css - Style sheet
--------------------------------------------------------------------*/

  public Text css ( String z )
  {
    return css ( s2t ( z ) );  // Call the com.scanfeld.core.Text version
  }

/*--------------------------------------------------------------------
 css - Style sheet
--------------------------------------------------------------------*/

  public Text css ( Text y )
  {
    Text z = y.copy ();
    z.lpush ( "<link rel=\"stylesheet\" type=\"text/css\" href=\"" );  // Push the opening tag
    z.rpush ( "\">" );                                                 // Push the closing tag
    return z;                                                          // Return the new com.scanfeld.core.Text
  }
}