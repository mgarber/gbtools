package broad.projection.utils;

import java.util.Collection;
import java.util.Set;
import java.util.TreeSet;


public class PubMed extends Inc implements Web
{
/*--------------------------------------------------------------------
 Class Variables
--------------------------------------------------------------------*/

  public String title;
  public String summary;
  public String longDescription;
  public String id;
  public String db;
  public Text webEnv;
  public String referenceSeriesNum;
  public String experimentSummary;
  public Url url;

/*--------------------------------------------------------------------
 CONSTRUCTOR FUNCTIONS
 ---------------------------------------------------------------------
 Initialize com.scanfeld.io.web.PubMed objects
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 com.scanfeld.io.web.PubMed - Main constructor, requires a pubmed id
--------------------------------------------------------------------*/

  public PubMed ( String db, String id)
  {
    this.webEnv=webEnv;
	  this.db = db;  // Set the database
    this.id = id;  // Set the id
    webLoad ();    // Load the com.scanfeld.io.web.PubMed from the web
  }

/*--------------------------------------------------------------------
 com.scanfeld.io.web.PubMed - Main constructor, requires a pubmed id
--------------------------------------------------------------------*/

  public PubMed ( String db, Text id )
  {
	  this.webEnv=webEnv;
    this.db = db;         // Set the database
    this.id = id.str ();  // Set the id
    webLoad ();           // Load the com.scanfeld.io.web.PubMed from the web
  }

/*--------------------------------------------------------------------
 webLoad - Load the geneset from the web
--------------------------------------------------------------------*/

  private String getReferenceNumber(String str){
	  String[] lines=str.split("\n");
	  String refNum;
	  boolean has=false;
	    for(int i=0; i<lines.length; i++){
	    	if(lines[i].startsWith("Reference")){
	    		try{
	    		refNum=lines[i].split(":")[1].trim();
	    		return refNum;
	    		}catch(Exception ex){System.err.println(lines[i]);}
	    	}
	    }
	    for(int i=0; i<lines.length; i++){
	    	if(lines[i].startsWith("1:")){
	    		refNum=lines[i].split(":")[1].replaceAll("record", "").trim();
	        	return refNum;
	        }
	    }
	    
	    return null;
  }
  
  private String getExperimentSummary(String str){
	  String[] lines=str.split("\n");
	  String rtrn;
	  boolean has=false;
	    for(int i=0; i<lines.length; i++){
	    	if(lines[i].startsWith("Summary")){
	    		rtrn=lines[i].split(":")[1].trim();
	    		return rtrn;
	    	}
	    }
	 return null;
  }
  
  public void webLoad ()
  {
    // http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gds&id=4394&report=docsum&mode=text
    // http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=GSM18462&id=10247&db=GeoDb_blob02&token=

    // LOAD THE PUBMED WEBSITE
    Url u = new Url ( "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=" + db + "&id=" + id  + "&report=docsum&mode=text"); //Text object
   Url t=new Url("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=" + db + "&id=" + id  + "&report=brief&mode=text");
   //System.err.println(t);
   this.url=t;
    String str=u.toString();
    this.referenceSeriesNum=this.getReferenceNumber(str);
    this.experimentSummary=t.toString().replaceAll("\n", " ").trim();
   
    System.err.println(referenceSeriesNum+" "+this.experimentSummary);
   
   
   // System.err.println(u);
    
    // FIND THE TITLE
    title = u.tag ( "ArticleTitle" ).get ();

    // FIND THE TITLE
    summary = u.tag ( "AbstractText" ).get ();

  }
  
  public String toString(){return this.url.toString();}
  
  public String getFtpCommand(){
	  return "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/series/"+this.referenceSeriesNum+"/"+this.referenceSeriesNum+"_RAW.tar";
  }
  
  private Collection<String> getSet(String speciesString, String delim){
	  Set<String> rtrn=new TreeSet();
	  
	  String[] tokens=speciesString.split(delim);
	  for(int i=0; i<tokens.length; i++){rtrn.add(tokens[i]);}
	  
	  return rtrn;
  }
  
  public String getSpecies(){
	  String speciesString=experimentSummary.split("\\[")[1].split("\\]")[0];
	  return speciesString;
  }
  
}
