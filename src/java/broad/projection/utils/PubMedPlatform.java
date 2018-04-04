package broad.projection.utils;



public class PubMedPlatform extends Inc implements Web
{
/*--------------------------------------------------------------------
 Class Variables
--------------------------------------------------------------------*/

  public String title;
  public String summary;
  public String longDescription;
  public String id;
  public String dir;
  public String pid;
  public String gid;
  public String did;
  public TextTree t;

/*--------------------------------------------------------------------
 CONSTRUCTOR FUNCTIONS
 ---------------------------------------------------------------------
 Initialize com.scanfeld.io.web.PubMed objects
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 com.scanfeld.io.web.PubMed - Main constructor, requires a pubmed id
--------------------------------------------------------------------*/

  public PubMedPlatform ( String id, String dir )
  {
    this.id = id;    // Set the id
    this.dir = dir;  // Set the directory
    this.pid = "";
    this.gid = "";
    this.did = "";
    webLoad ();      // Load the com.scanfeld.io.web.PubMed from the web
  }

/*--------------------------------------------------------------------
 com.scanfeld.io.web.PubMed - Main constructor, requires a pubmed id
--------------------------------------------------------------------*/

  public PubMedPlatform ( String id, String dir, String pid, String gid, String did )
  {
    this.id = id;    // Set the id
    this.dir = dir;  // Set the directory
    this.pid = pid;
    this.gid = gid;
    this.did = did;
    webLoad ();      // Load the com.scanfeld.io.web.PubMed from the web
  }

/*--------------------------------------------------------------------
 com.scanfeld.io.web.PubMedSample - Main constructor, requires a pubmed id
--------------------------------------------------------------------*/

  public PubMedPlatform ( Text id, Text dir )
  {
    this.id = id.str ();    // Set the id
    this.dir = dir.str ();  // Set the directory
    this.pid = "";
    this.gid = "";
    this.did = "";
    webLoad ();             // Load the com.scanfeld.io.web.PubMed from the web
  }

/*--------------------------------------------------------------------
 webLoad - Load the geneset from the web
--------------------------------------------------------------------*/

  public void webLoad ()
  {
    // CHECK IF THE PLATFORM HAS BEEN LOADED
    if ( !exists ( dir + "GPL" + id + ".txt" ) )
    {
      // LOAD THE PLATFORM WEBSITE
      Url u = new Url ( "http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?view=data&acc=GPL" + id, dir + "GPL" + id + ".txt" );
    }

    // CHECK IF THE PLATFORM MAP HAS BEEN CREATED
    if ( !exists ( dir + "GPL" + id + ".map" ) )
    {
      Infile a = new Infile ( dir + "GPL" + id + ".txt" );    // com.scanfeld.io.Infile for the website
      Outfile b = new Outfile ( dir + "GPL" + id + ".map" );  // com.scanfeld.io.Outfile for the map

      Texts g = new Texts ();  // com.scanfeld.bio.Gene list
      Texts v = new Texts ();  // Values list

      a.readLine ();  // Read in the first line

      // LOOP TO THE FIRST GENE
      while ( !a.eof && a.sub ( 0, 3 ).ne ( "<st" ) )
      {
        a.readLine ();  // Read in the next line
      }

      Texts tags = new Texts ( a.tags ( "strong" ) );  // Find all values between strong tags
      tags = tags.trimWhite ().lower ();               // Remove whitespace

      int id = 0;
      int tid = tags.find ( new Text ( pid ).trimWhite ().lower () );  // Find the location of the probe id

      if ( tid >= 0 )
      {
        id = tid;
      }

      int mid = 0;
      int tmid = tags.find ( new Text ( gid ).trimWhite ().lower () );

      if ( tmid >= 0 )
      {
        mid = tmid;
      }

      Text t1, t2;
      a.readLine ();

      while ( !a.eof )
      {
        if ( a.size == 0 )
        {
          break;
        }

        Texts x = new Texts ( a.splitTab () ).removeTag ();

        if ( x.size > id && x.size > mid )
        {
          t1 = x.get ( id );
          t2 = x.get ( mid );

          if ( t1 != null && t1.size > 0 && t2 != null && t2.size > 0 )
          {
            g.rpush ( x.get ( id ) );
            v.rpush ( x.get ( mid ) );
          }
        }

        a.readLine ();
      }

      int[] ind = g.sort ();
      v = v.select ( ind );

      for ( int i = 0; i < g.size; i++ )
      {
        b.writeLine ( g.get ( i ) + "\t" + v.get ( i ) );
      }

      b.close ();
      a.close ();
    }

    if ( !exists ( dir + "GPL" + id + ".des" ) )
    {
      Infile a = new Infile ( dir + "GPL" + id + ".txt" );
      Outfile b = new Outfile ( dir + "GPL" + id + ".des" );
      a.readLine ();

      TextTree d = new TextTree ();
      Texts g = new Texts ();

      a.readLine ();
      while ( !a.eof && a.sub ( 0, 3 ).ne ( "<st" ) )
      {
        a.readLine ();
      }

      Texts tags = new Texts ( a.tags ( "strong" ) );
      tags = tags.trimWhite ().lower ();

      int id = 0;
      int tid = tags.find ( new Text ( gid ).trimWhite ().lower () );

      if ( tid >= 0 )
      {
        id = tid;
      }

      int mid = 0;
      int tmid = tags.find ( new Text ( did ).trimWhite ().lower () );

      if ( tmid >= 0 )
      {
        mid = tmid;
      }

      a.readLine ();
      while ( !a.eof )
      {
        if ( a.size == 0 )
        {
          break;
        }

        Text t1, t2;
        Texts x = new Texts ( a.splitTab () ).removeTag ();

        if ( x.size > id && x.size > mid )
        {
          t1 = x.get ( id );
          t2 = x.get ( mid );

          if ( t1 != null && t1.size > 0 && t2 != null && t2.size > 0 )
          {
            d.set ( t1, t2 );
            g.rpush ( t1 );
          }
        }

        a.readLine ();
      }

      g.sort ();
      g = g.removeDuplicate ();

      for ( int i = 0; i < g.size; i++ )
      {
        b.writeLine ( g.get ( i ) + "\t" + d.get ( g.get ( i ) ) );
      }

      b.close ();
      a.close ();
    }

    t = new TextTree ();
    t.map ( dir + "GPL" + id + ".map" );
  }

/*--------------------------------------------------------------------
 get - Get a gene
--------------------------------------------------------------------*/

  public Text get ( String name )
  {
    return t.get ( name );
  }

/*--------------------------------------------------------------------
 get - Get a gene
--------------------------------------------------------------------*/

  public Text get ( Text name )
  {
    return t.get ( name );
  }
}
