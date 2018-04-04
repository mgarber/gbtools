package broad.projection.utils;

public class PubMedSeries extends Inc implements Web
{
/*--------------------------------------------------------------------
 Class Variables
--------------------------------------------------------------------*/

  public String title;
  public String summary;
  public String longDescription;
  public String id;
  public String pubMedId;
  public String platformId;
  public String platformName;
  public String dir;
  public String pid;
  public String gid;
  public String did;
  public PubMedPlatform p;

/*--------------------------------------------------------------------
 CONSTRUCTOR FUNCTIONS
 ---------------------------------------------------------------------
 Initialize com.scanfeld.io.web.PubMed objects
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 com.scanfeld.io.web.PubMed - Main constructor, requires a pubmed id
--------------------------------------------------------------------*/

  public PubMedSeries ( String id )
  {
    this.id = id;    // Set the id
    this.dir = "c:/com.scanfeld/pr/da/GSE" + this.id + "/";
    this.pid = "ID";
    this.gid = "com.scanfeld.bio.Gene Symbol";
    this.did = "com.scanfeld.bio.Gene Title";
    webLoad ();    // Load the com.scanfeld.io.web.PubMed from the web
  }

/*--------------------------------------------------------------------
 com.scanfeld.io.web.PubMed - Main constructor, requires a pubmed id
--------------------------------------------------------------------*/

  public PubMedSeries ( String id, String dir )
  {
    this.id = id;    // Set the id
    this.dir = dir;  // Set the directory
    this.pid = "ID";
    this.gid = "com.scanfeld.bio.Gene Symbol";
    this.did = "com.scanfeld.bio.Gene Title";
    webLoad ();    // Load the com.scanfeld.io.web.PubMed from the web
  }

/*--------------------------------------------------------------------
 com.scanfeld.io.web.PubMed - Main constructor, requires a pubmed id
--------------------------------------------------------------------*/

  public PubMedSeries ( String id, String pid, String gid, String did )
  {
    this.id = id;    // Set the id
    this.dir = "c:/com.scanfeld/pr/da/GSE" + this.id + "/";
    this.pid = pid;
    this.gid = gid;
    this.did = did;
    webLoad ();    // Load the com.scanfeld.io.web.PubMed from the web
  }

/*--------------------------------------------------------------------
 com.scanfeld.io.web.PubMed - Main constructor, requires a pubmed id
--------------------------------------------------------------------*/

  public PubMedSeries ( String id, String pid, String gid, String did, String dir )
  {
    this.id = id;    // Set the id
    this.dir = dir;  // Set the directory
    this.pid = pid;
    this.gid = gid;
    this.did = did;
    webLoad ();    // Load the com.scanfeld.io.web.PubMed from the web
  }

/*--------------------------------------------------------------------
 com.scanfeld.io.web.PubMed - Main constructor, requires a pubmed id
--------------------------------------------------------------------*/

  public PubMedSeries ( Text id, Text dir )
  {
    this.id = id.str ();    // Set the id
    this.dir = dir.str ();  // Set the directory
    this.pid = "";
    this.gid = "";
    this.did = "";
    webLoad ();             // Load the com.scanfeld.io.web.PubMed from the web
  }

public void webLoad() {
	// TODO Auto-generated method stub
	
}

/*--------------------------------------------------------------------
 webLoad - Load the geneset from the web
--------------------------------------------------------------------*/

  /*public void webLoad ()
  {
    Url u = new Url ( "http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE" + id );

    title = u.findBetween ( "Title", "</tr>" ).findBetween ( "<td>", "</td>" ).str ();
    summary = u.findBetween ( "Summary", "</tr>" ).findBetween ( "<td>", "</td>" ).replace ( "<br>", "\n" ).str ();
    pubMedId = u.findBetween ( "com.scanfeld.io.web.PubMed ID", "</tr>" ).findBetween ( "term=", "\">" ).str ();

    Text t = u.findBetween ( "acc=GPL", "</tr>" );
    platformId = t.findBetween ( ">GPL", "</a>" ).str ();
    platformName = t.findBetween ( "<td valign=\"top\">", "</td>" ).str ();

    Array na = new Array ();
    p = new PubMedPlatform ( platformId, dir, pid, gid, did );

    Text[] si = u.findBetweens ( "acc=GSM", "</tr>" );

    for ( int i = 0; i < si.length; i++ )
    {
      si[ i ] = si[ i ].findBetween ( ">GSM", "</a>" );

      if ( si[ i ].length () > 0 )
      {
        PubMedSample s = new PubMedSample ( si[ i ].str (), si[ i ].str (), "none", dir, p, na );
      }
    }

    Outfile b = new Outfile ( dir + "GSE" + id + "GSM.txt" );

    for ( int i = 0; i < si.length; i++ )
    {
      b.writeLine ( "GSM" + si[ i ] );
    }

    b.close ();

    b = new Outfile ( dir + "GSE" + id + ".txt" );
    b.writeLine ( title );
    b.writeLine ( pubMedId );
    b.writeLine ( platformName );
    b.writeLine ( summary );
    b.close ();

    na.map ( dir + "GPL" + platformId + ".map" );
    na.removeDuplicate ();

    TextTree des = new TextTree ();
    des.map ( dir + "GPL" + platformId + ".des" );
    na.addDescription ( des );

    na.write ( dir + "GSE" + id + ".gct", dir + "GSE" + id + ".cls" );
  }*/
}