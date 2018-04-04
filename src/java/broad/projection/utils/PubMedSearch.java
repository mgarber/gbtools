package broad.projection.utils;

import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;


public class PubMedSearch extends Url
{
/*--------------------------------------------------------------------
 Class Variables
--------------------------------------------------------------------*/

  public String db;        // pubmed, gds
  public String search;
  public int max;
  public int size;
  public int fullSize;
  public Text[] id;
  public PubMed[] p;
  
/*--------------------------------------------------------------------
 CONSTRUCTOR FUNCTIONS
 ---------------------------------------------------------------------
 Initialize com.scanfeld.io.web.PubMed objects
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 com.scanfeld.io.web.PubMedSearch - Main constructor, requires a search string and a db
--------------------------------------------------------------------*/

  public PubMedSearch ( String db, String search, int max )
  {
    this.db = db;        // Set the database
    this.max = max;      // Set the maximum number of results
    loadWeb ( search );  // Load the com.scanfeld.io.web.PubMedSearch from the web
  }

/*--------------------------------------------------------------------
 com.scanfeld.io.web.PubMedSearch - Main constructor, requires a search string and a db
--------------------------------------------------------------------*/

  public PubMedSearch ( Text db, Text search )
  {
    this.db = db.get ();        // Set the database
    loadWeb ( search.get () );  // Load the com.scanfeld.io.web.PubMedSearch from the web
  }

/*--------------------------------------------------------------------
 loadWeb - Load a com.scanfeld.io.web.PubMedSearch from the web
--------------------------------------------------------------------*/

  public void loadWeb ( String search )
  {
    this.search = search;  // Save the search string

    // CREATE A NEW URL
    Url u = new Url ();
    
    // LOAD THE PUBMEDSEARCH WEBSITE
    u.load ( "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=" + db + "&usehistory=no&retmode=xml&retmax=" + max + "&term=" + search );

    // GET THE ID LIST
    Text idList = u.tag ( "IdList" );

    // GET THE IDS
    id = idList.tags ( "Id" );

    // THE NUMBER OF IDS
    size = id.length;

    // GET THE FULL SIZE
    fullSize = u.tag ( "Count" ).num ();
    Text webEnv=u.tag("WebEnv");
    //System.err.println(webEnv);
    
    //System.err.println(size+" "+fullSize);

    p = new PubMed[size];  // Create a com.scanfeld.io.web.PubMed array

    for ( int i = 0; i < size; i++ )
    {
    	//System.err.println(id[i]);
      p[ i ] = new PubMed ( db, id[ i ] );
     
    }
  }
  
  public void launchFTPGet(String speciesFilter){
	  String ftp="";
	  for(int i=0; i<p.length; i++){
		 String speciesString= p[i].experimentSummary.split("\\[")[1].split("\\]")[0];
		 Set<String> species=getSet(speciesString, ";");
		 if(species.contains(speciesFilter)){
			 ftp+=(" "+p[i].getFtpCommand());
		 }
	  }
	  System.err.println(ftp);
  }
  
  private Set getSet(String speciesString, String delim){
	  Set rtrn=new TreeSet();
	  
	  String[] tokens=speciesString.split(delim);
	  for(int i=0; i<tokens.length; i++){rtrn.add(tokens[i]);}
	  
	  return rtrn;
  }
  
  public void writeSummary(String save)throws IOException{
	  FileWriter writer=new FileWriter(save);
	  
	  for(int i=0; i<this.p.length; i++){
		 String speciesString= p[i].experimentSummary.split("\\[")[1].split("\\]")[0];
		 if(p[i].referenceSeriesNum!=null){ writer.write(p[i].referenceSeriesNum+"\t"+speciesString+"\t"+p[i].experimentSummary+"\t"+p[i].getFtpCommand()+"\n");}
		 else{System.err.println(p[i]);}
	  }
	  
	  writer.close();
  }
  
  
  public static void main(String[] args)throws IOException{
	  if(args.length>4){
		  String db=args[0];
		  String[] search=args[1].split(",");
		  int max=new Integer(args[2]);
		  String save=args[3];
		  
		  Map<String, PubMed> set=new HashMap<String, PubMed>();
		  for(int i=0; i<search.length; i++){
			  PubMedSearch p=new PubMedSearch(db, search[i], max);
			  set.putAll(p.getFTPs());
		  }
		  
		  writeSearch(set, save);
		  
	  }
	  else{System.err.println(usage);}
  }
 
  private static void writeSearch(Map<String, PubMed> set, String save) throws IOException {
	FileWriter writer=new FileWriter(save);
	
	for(String ftp: set.keySet()){
		writer.write(ftp+"\t"+set.get(ftp).experimentSummary+"\t"+set.get(ftp).getSpecies()+"\n");
	}
	
	writer.close();
  }

private Map<String, PubMed> getFTPs() {
	 Map<String, PubMed> rtrn=new HashMap(); 
	  for(int i=0; i<p.length; i++){
		  String ftp=p[i].getFtpCommand();
		  rtrn.put(ftp, p[i]);
	  }
	return rtrn;
  }

static String usage=" args[0]=database (gds) \n args[1]=search phrases (comma seperated) \n args[2]=maxNumResults \n args[3]=save";
  
}
