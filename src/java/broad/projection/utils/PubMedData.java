package broad.projection.utils;


public class PubMedData extends Inc implements Web
{
/*--------------------------------------------------------------------
 Class Variables
--------------------------------------------------------------------*/

  public String title;
  public String summary;
  public String longDescription;
  public String id;
  public String dir;
  public PubMedSample[] s;

/*--------------------------------------------------------------------
 CONSTRUCTOR FUNCTIONS
 ---------------------------------------------------------------------
 Initialize com.scanfeld.io.web.PubMed objects
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 com.scanfeld.io.web.PubMed - Main constructor, requires a pubmed id
--------------------------------------------------------------------*/

  public PubMedData ( String id, String dir )
  {
    this.id = id;    // Set the id
    this.dir = dir;  // Set the directory
    webLoad ();    // Load the com.scanfeld.io.web.PubMed from the web
  }

/*--------------------------------------------------------------------
 com.scanfeld.io.web.PubMed - Main constructor, requires a pubmed id
--------------------------------------------------------------------*/

  public PubMedData ( Text id, Text dir )
  {
    this.id = id.str ();    // Set the id
    this.dir = dir.str ();  // Set the directory
    webLoad ();             // Load the com.scanfeld.io.web.PubMed from the web
  }

/*--------------------------------------------------------------------
 webLoad - Load the geneset from the web
--------------------------------------------------------------------*/

  public void webLoad ()
  {
    // http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gds&id=4394&report=docsum&mode=text
    // http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=GSM18462&id=10247&db=GeoDb_blob02&token=

    // PLATFORM
    // http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=GPL90&id=23603&db=GeoDb_blob03&token=

    // LOAD THE PUBMED WEBSITE
    Url u = new Url ( "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gds&id=" + id + "&report=docsum&mode=text" );

    out ( u );
    Text[] si = u.findBetweens ( "GSM", ":" );
    s = new PubMedSample[si.length];

    for ( int i = 0; i < 1; i++ ) // s.length; i++ )
    {
      //s[ i ] = new com.scanfeld.io.web.PubMedSample ( si[ i ].str (), si[ i ].str (), "none", dir, null );
    }

    // FIND THE TITLE
    // title = u.tag ( "ArticleTitle" ).get ();

    // FIND THE TITLE
    // summary = u.tag ( "AbstractText" ).get ();
  }
}
