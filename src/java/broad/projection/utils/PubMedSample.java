package broad.projection.utils;

import broad.projection.nmf.Array;


public class PubMedSample extends Inc implements Web
{
/*--------------------------------------------------------------------
 Class Variables
--------------------------------------------------------------------*/

  public String name;
  public String phen;
  public String title;
  public String summary;
  public String longDescription;
  public String id;
  public String dir;
  public PubMedPlatform p;
  public Array dat;

/*--------------------------------------------------------------------
 CONSTRUCTOR FUNCTIONS
 ---------------------------------------------------------------------
 Initialize com.scanfeld.io.web.PubMed objects
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 com.scanfeld.io.web.PubMed - Main constructor, requires a pubmed id
--------------------------------------------------------------------*/

  public PubMedSample ( String id, String name, String phen, String dir, PubMedPlatform p, Array dat )
  {
    this.dat = dat;
    this.name = name;
    this.phen = phen;
    this.p = p;
    this.id = id;    // Set the id
    this.dir = dir;  // Set the directory
    webLoad ();      // Load the com.scanfeld.io.web.PubMed from the web
  }

/*--------------------------------------------------------------------
 com.scanfeld.io.web.PubMedSample - Main constructor, requires a pubmed id
--------------------------------------------------------------------*/

  public PubMedSample ( Text id, Text name, Text phen, Text dir, PubMedPlatform p, Array dat )
  {
    this.dat = dat;
    this.name = name.str ();
    this.phen = phen.str ();
    this.p = p;
    this.id = id.str ();    // Set the id
    this.dir = dir.str ();  // Set the directory
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
    if ( !exists ( dir + "GSM" + id + ".txt" ) )
    {
      // LOAD THE PUBMED WEBSITE // why accession number and id?
      Url u = new Url ( "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=GSM" + id + "&db=GeoDb_blob02&token=", dir + "GSM" + id + ".txt" );
    }

    if ( !exists ( dir + "GSM" + id + ".exp" ) || !exists ( dir + "GSM" + id + "gene.exp" ) || !exists ( dir + "GSM" + id + ".gct" ) || !exists ( dir + "GSM" + id + ".cls" ) )
    {
      Infile a = new Infile ( dir + "GSM" + id + ".txt" );
      a.readLine ();

      while ( !a.eof )
      {
        if ( a.find ( "<strong>" ) >= 0 )
        {
          break;
        }

        a.readLine ();
      }

      Texts g = new Texts ();
      Texts v = new Texts ();

      a.readLine ();
      while ( !a.eof )
      {
        Text[] x = a.splitTab ();
        if ( x.length > 3 && x[ 0 ] != null && x[ 0 ].size > 0 && x[ 1 ] != null && x[ 1 ].size > 0 )
        {
          g.rpush ( x[ 0 ] );
          v.rpush ( x[ 1 ] );
        }

        a.readLine ();
      }
      a.close ();

      int[] ind = g.sort ();
      v = v.select ( ind );

      Outfile b = new Outfile ( dir + "GSM" + id + ".exp" );

      for ( int i = 0; i < g.size; i++ )
      {
        b.writeLine ( g.get ( i ) + "\t" + v.get ( i ) );
      }

      b.close ();

      if ( dat.numGene == 0 )
      {
        dat.loadEXP ( g, v, name, phen );
        dat.write ( dir + "GSM" + id + ".gct", dir + "GSM" + id + ".cls" );
      }

      else
      {
        Array z = new Array ();
        z.loadEXP ( g, v, name, phen );
        z.write ( dir + "GSM" + id + ".gct", dir + "GSM" + id + ".cls" );
        dat.merge ( z );
      }

      TextTree t = new TextTree ();
      Texts ng = new Texts ();

      Text t1;
      for ( int i = 0; i < g.size; i++ )
      {
        t1 = p.get ( g.get ( i ) );

        if ( t1 != null )
        {
          ng.rpush ( t1 );
          Text val = t.get ( g.get ( i ) );

          if ( val == null || val.dbl () < v.get ( i ).dbl () )
          {
            t.set ( t1, v.get ( i ) );
          }
        }
      }

      ng.sort ();
      ng = ng.removeDuplicate ();

      b = new Outfile ( dir + "GSM" + id + "gene.exp" );

      v = new Texts ();
      for ( int i = 0; i < ng.size; i++ )
      {
        b.writeLine ( ng.get ( i ) + "\t" + t.get ( ng.get ( i ) ) );
        v.rpush ( t.get ( ng.get ( i ) ) );
      }

      b.close ();
    }
  }*/
}
