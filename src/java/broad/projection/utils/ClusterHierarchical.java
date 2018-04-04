package broad.projection.utils;

import broad.projection.math.Matrix;


/*--------------------------------------------------------------------
 ClusterHierarchical
 ---------------------------------------------------------------------
 Models a hierarchical clustering
 ---------------------------------------------------------------------
 Author: Daniel Scanfeld
   Date: 08.23.06
 ---------------------------------------------------------------------
 test to ensure the correct results
 could speed up the algorithm a bit
 make it not a matrix, but a more sensible data type
 function to display as a pdf
 implement all the distance metrics from gene pattern
   uncentered correlation
   pearson correlation
   uncentered correlation, absolute value
   pearson correlation, absolute value
   spearman's rank correlation
   kendall's tau
   euclidean distance
   city-block distance

   pairwise complete linkage
   pairwise single linkage
   pairwise centroid linkage
   pairwise average linkage
--------------------------------------------------------------------*/

public class ClusterHierarchical extends Inc
{
/*--------------------------------------------------------------------
 Class Variables
--------------------------------------------------------------------*/

  public Matrix m;     // Cluster com.scanfeld.math.Matrix
  public int order[];  // Ordering
  public Conga[] c;    // An array of Congas
  public int cc;       // Current Conga
  public int logCut;   // Log cutoff

  public int nc;
  public String metric;
  public String type;
  public int[] seta;      // Left Branch
  public int[] setb;      // Right Branch
  public double[] setd;   // Distance
  public boolean[] setf;  // Set Existence Flag

/*--------------------------------------------------------------------
 CONSTRUCTOR FUNCTIONS
 ---------------------------------------------------------------------
 Initialize ClusterHierarchical objects
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 ClusterHierarchical - Null constructor
--------------------------------------------------------------------*/

  public ClusterHierarchical ()
  {
  }

/*--------------------------------------------------------------------
 ClusterHierarchical - com.scanfeld.math.Matrix constructor
--------------------------------------------------------------------*/

  //TODO Implement a clustering method that starts with groups rather than individual samples
  public ClusterHierarchical ( Matrix m, String metric, String type )
  {
    cluster ( m, metric, type );
  }

/*--------------------------------------------------------------------
 dis - Measure the distance between two sets
--------------------------------------------------------------------*/

  public double dis ( int a, int b )
  {
    double d = 0;

    // COMPLETE LINKAGE - MAXIMUM DISTANCE
    if ( type.equals ( "complete" ) )
    {
      if ( a < nc )
      {
        if ( b < nc )
        {
          d = -m.pearson ( a, b ) + 1;
        }

        else
        {
          d = max ( dis ( a, seta[ b ] ), dis ( a, setb[ b ] ) );
        }
      }

      else
      {
        d = max ( dis ( seta[ a ], b ), dis ( setb[ a ], b ) );
      }
    }

    // SINGLE LINKAGE - MINIMUM DISTANCE
    else
    {
      if ( a < nc )
      {
        if ( b < nc )
        {
          d = -m.pearson ( a, b ) + 1;
        }

        else
        {
          d = min ( dis ( a, seta[ b ] ), dis ( a, setb[ b ] ) );
        }
      }

      else
      {
        d = min ( dis ( seta[ a ], b ), dis ( setb[ a ], b ) );
      }
    }

    return d;
  }

/*--------------------------------------------------------------------
 cluster - Cluster the matrix
--------------------------------------------------------------------*/

  public void cluster ( Matrix m, String metric, String type )
  {
    if ( metric.equals ( "spearman" ) )
    {
      m = m.rank ();
    }

    out ( "Hierarchical Clustering: " + metric + " " + type );

    this.m = m;
    this.metric = metric;
    this.type = type;
    nc = m.nc;

    seta = new int[nc * 2];
    setb = new int[nc * 2];
    setd = new double[nc * 2];

    
    //Initialize the sets so that each sample is 1 set
    for ( int i = 0; i < nc; i++ ){
      seta[ i ] = i;
      setb[ i ] = -1;
    }

    int cs = nc;
    for ( int i = nc; i < nc * 2; i++ ){
      seta[ i ] = -1;
      setb[ i ] = -1;
    }

    double[][] d = new double[nc * 2][];
    int[] md = new int[nc * 2];

    for ( int i = 0; i < nc; i++ ){
      sts ( "Calculating distance " + i + " of " + nc );
      d[ i ] = new double[nc * 2];
      md[ i ] = -1;

      
      //Loop through each column
      for ( int j = 0; j < nc; j++ ){
        if ( j < i ){
          d[ i ][ j ] = d[ j ][ i ];
        }
        //Compute distance between samples i and j
        else if ( j > i ){
          d[ i ][ j ] = dis ( i, j );
          // out ( i + "\t" + j + "\t" + d[ i ][ j ] );
        }

        else{continue;}

        if ( md[ i ] < 0 || d[ i ][ j ] < d[ i ][ md[ i ] ] ){md[ i ] = j;}
      }
    }

    for ( int i = 0; i < nc - 1; i++ )  {
      sts ( "Calculating merge " + i + " of " + ( nc - 1 ) );
      int mx1 = -1;
      int mx2 = -1;

      for ( int j = 0; j < md.length; j++ )
      {
        if ( d[ j ] != null && ( mx1 < 0 || d[ j ][ md[ j ] ] < d[ mx1 ][ mx2 ] ) )
        {
          mx1 = j;
          mx2 = md[ j ];
        }
      }

      seta[ cs ] = mx1;
      setb[ cs ] = mx2;
      setd[ cs ] = d[ mx1 ][ mx2 ];
      md[ cs ] = -1;

      for ( int j = 0; j < md.length; j++ ){
        if ( mx1 != j && mx2 != j && d[ j ] != null )
        {
          d[ mx1 ][ j ] = max ( d[ mx1 ][ j ], d[ mx2 ][ j ] );

          if ( md[ cs ] < 0 || d[ mx1 ][ j ] < d[ mx1 ][ md[ cs ] ] )
          {
            md[ cs ] = j;
          }
        }
      }

      d[ cs ] = d[ mx1 ];
      d[ mx1 ] = null;
      d[ mx2 ] = null;

      for ( int j = 0; j < md.length; j++ )
      {
        if ( d[ j ] != null && j != cs )
        {
          d[ j ][ cs ] = d[ cs ][ j ];

          if ( md[ j ] == mx1 || md[ j ] == mx2 )
          {
            md[ j ] = -1;
            for ( int k = 0; k < md.length; k++ )
            {
              if ( j != k && d[ k ] != null )
              {
                if ( md[ j ] < 0 || d[ j ][ k ] < d[ j ][ md[ j ] ] )
                {
                  md[ j ] = k;
                }
              }
            }
          }

          else if ( d[ j ][ cs ] < d[ j ][ md[ j ] ] )
          {
            md[ j ] = cs;
          }
        }
      }

      cs++;
    }

    Texts ord = sort ( ( nc * 2 ) - 2 );
    ord.lpop ();
    ord.lpop ();
    order = ord.num ();
  }

/*--------------------------------------------------------------------
 cluster3 - Cluster the matrix
--------------------------------------------------------------------*/

  public void cluster3 ( Matrix m, String metric, String type )
  {
    if ( metric.equals ( "spearman" ) )
    {
      m = m.rank ();
    }

    out ( "Hierarchical Clustering: " + metric + " " + type );
    this.m = m;
    this.metric = metric;
    this.type = type;
    nc = m.nc;

    logCut = ceil ( log2 ( nc ) );
    c = new Conga[logCut * 3];
    cc = 0;

    seta = new int[nc * 2];
    setb = new int[nc * 2];
    setd = new double[nc * 2];
    setf = new boolean[nc * 2];

    for ( int i = 0; i < nc; i++ )
    {
      seta[ i ] = i;
      setb[ i ] = -1;
      setf[ i ] = true;
    }

    for ( int i = nc; i < nc * 2; i++ )
    {
      seta[ i ] = -1;
      setb[ i ] = -1;
      setf[ i ] = false;
    }

    Conga g = next ();
    g.setAll ();
    g.draw ();

    int numSet = b2i ( setf ).length;

    while ( numSet > 1 )
    {
      sts ( numSet + " sets to be clustered." );

      double d = -1;
      int ci = 0;
      int cj = 0;

      // out ( "Find the smallest distance" );

      // FIND THE SMALLEST DISTANCE
      // com.scanfeld.core.Texts pts = new com.scanfeld.core.Texts ();
      for ( int i = 0; i < c.length; i++ )
      {
        if ( c[ i ] == null )
        {
          continue;
        }

        int[] pt = c[ i ].pt.num ();
        double[] dis = c[ i ].dis.dbl ();

        // out ( dis );
        // pts.rpush ( c[ i ].pt );

        int mn = minInd ( dis );
        if ( d < 0 || dis[ mn ] < d )
        {
          d = dis[ mn ];
          ci = i;
          cj = mn;
        }
      }

      // RECORD THE TWO POINTS
      int ra = c[ ci ].pt.get ( cj ).num ();
      int rb = c[ ci ].pt.get ( cj + 1 ).num ();

      /*
      // TEST THE CLUSTERING
      
      double min = 10000;
      int mina = -1;
      int minb = -1;

      for ( int i = 0; i < setf.length; i++ )
      {
        if ( setf[ i ] )
        {
          for ( int j = i + 1; j < setf.length; j++ )
          {
            if ( setf[ j ] )
            {
              double ijd = dis ( i, j );
              if ( ijd < min )
              {
                min = ijd;
                mina = i;
                minb = j;
              }
            }
          }
        }
      }

      if ( mina + minb != ra + rb )
      {
        out ( "oops" );
      }
      */

      /*
      pts.sort ();
      pts = pts.removeDuplicate ();
      int[] ptsi = pts.num ();
      sort ( ptsi );
      */

      // FIND THE NEXT SPOT TO STORE THE SET
      int rc = 0;
      for ( rc = 0; rc < seta.length; rc++ )
      {
        if ( seta[ rc ] < 0 )
        {
          seta[ rc ] = ra;
          setb[ rc ] = rb;
          setd[ rc ] = d;
          break;
        }
      }

      setf[ rc ] = true;

      // out ( "Remove the first point" );
      // REMOVE THE FIRST POINT
      remove ( ra );

      for ( int i = 0; i < c.length; i++ )
      {
        if ( c[ i ] == null )
        {
          continue;
        }
      }

      // out ( "Remove the second point" );
      // REMOVE THE SECOND POINT
      remove ( rb );

      // out ( "Create the next conga line" );
      // CREATE THE NEXT CONGA LINE
      g = next ();
      g.setOne ( rc );
      g.draw ();

      numSet--;
    }

    Texts ord = sort ( ( nc * 2 ) - 2 );
    ord.lpop ();
    ord.lpop ();
    order = ord.num ();
  }

/*--------------------------------------------------------------------
 sort - Sort a set
--------------------------------------------------------------------*/

  public Texts sort ( int k )
  {
    /*
    In hierarchical cluster displays, a decision is needed at each merge to specify which subtree
    should go on the left and which on the right. Since, for n observations there are n-1 merges,
    there are 2^{(n-1)} possible orderings for the leaves in a cluster tree, or dendrogram.
    The algorithm used in hclust is to order the subtree so that the tighter cluster is on the left
    (the last, i.e., most recent, merge of the left subtree is at a lower value than the last merge
    of the right subtree). Single observations are the tightest clusters possible, and merges
    involving two observations place them in order by their observation sequence number.
     */
    Texts a = new Texts ();
    boolean single = false;

    if ( setb[ k ] < 0 )
    {
      a.rpush ( k );
      a.rpush ( k );
      a.rpush ( k );
    }

    else
    {
      a = sort ( seta[ k ] );
      int p1 = a.lpop ().num ();
      int p2 = a.lpop ().num ();
      Texts b = sort ( setb[ k ] );
      int p3 = b.lpop ().num ();
      int p4 = b.lpop ().num ();

      /*
      double[] d = new double[4];
      d[ 0 ] = dis ( p1, p3 );
      d[ 1 ] = dis ( p1, p4 );
      d[ 2 ] = dis ( p2, p3 );
      d[ 3 ] = dis ( p2, p4 );

      // if ( ( p1 == 0 && p1 == 13 ) || ( p2 == 0 && p3 == 13 ) )
      // {
      // out ( p1 + " " + p2 + " " + p3 + " " + p4 );
      // }

      int mx = minInd ( d );
      // out(mx+ " " + showOne( d ));

      if ( mx == 0 )
      {
        a.reverse ();
      }

      else if ( mx == 1 )
      {
        a.reverse ();
        b.reverse ();
      }

      else if ( mx == 3 )
      {
        b.reverse ();
      }
      */

      if( setb[k] < seta[k] )
      {
        Texts c = a;
        a = b;
        b = c;
      }

      Text astart = b.get ( b.size-1 );
      Text bend = a.get( 0 );

      a.rpush ( b );

      // Compare to reset of the set
      // a.lpush ( setb[ k ] );
      // a.lpush ( seta[ k ] );

      // Compare to just the edges
      a.lpush ( bend );
      a.lpush ( astart );
    }

    // out ( k + ": " + a.joinSpace () );
    return a;
  }

/*--------------------------------------------------------------------
 remove - Remove point k
--------------------------------------------------------------------*/

  public void remove ( int k )
  {
    Texts newSet = new Texts ();

    setf[ k ] = false;
    for ( int i = 0; i < c.length; i++ )
    {
      if ( c[ i ] == null )
      {
        continue;
      }

      int ind = c[ i ].remove ( k );

      if ( ind >= 0 )
      {
        if ( ind > 0 )
        {
          newSet.rpush ( c[ i ].pt.get ( ind - 1 ) );
        }

        delete ( i );
      }
    }

    if ( newSet.size > 0 )
    {
      Conga g = next ();

      for ( int i = 0; i < newSet.size; i++ )
      {
        g.s[ newSet.get ( i ).num () ] = true;
      }

      g.draw ();
    }

    merge ();

    Tree all = new Tree ();
    for ( int i = 0; i < c.length; i++ )
    {
      if ( c[ i ] != null )
      {
        for ( int j = 0; j < c[ i ].s.length; j++ )
        {
          if ( c[ i ].s[ j ] )
          {
            all.set ( j, i + 1 );
          }
        }
      }
    }
  }

/*--------------------------------------------------------------------
 next - Get a new Conga line
--------------------------------------------------------------------*/

  public Conga next ()
  {
    for ( int i = 0; i < c.length; i++ )
    {
      if ( c[ i ] == null )
      {
        c[ i ] = new Conga ( this );
        cc++;
        return c[ i ];
      }
    }

    return null;
  }

/*--------------------------------------------------------------------
 delete - Delete the Conga line
--------------------------------------------------------------------*/

  public void delete ( int i )
  {
    c[ i ] = null;
    cc--;
  }

/*--------------------------------------------------------------------
 merge - See if we're over the limit, if so merge
--------------------------------------------------------------------*/

  public void merge ()
  {
    while ( cc > logCut )
    {
      // out ( "merge" );
      Texts setInd = new Texts ();
      Texts setSize = new Texts ();

      for ( int i = 0; i < c.length; i++ )
      {
        if ( c[ i ] != null )
        {
          setInd.rpush ( i );
          setSize.rpush ( c[ i ].setCount );
        }
      }

      int[] intSetSize = setSize.num ();
      int ind[] = sort ( intSetSize );
      int[] intSetInd = sel ( setInd.num (), ind );

      int min = 0;
      for ( int i = 1; i < intSetSize.length - 1; i++ )
      {
        if ( ( double ) ( intSetSize[ i + 1 ] ) / ( double ) ( intSetSize[ i ] ) < ( double ) ( intSetSize[ min + 1 ] ) / ( double ) ( intSetSize[ min ] ) )
        {
          min = i;
        }
      }
      min = 0;

      int ca = intSetInd[ min ];
      int cb = intSetInd[ min + 1 ];

      for ( int i = 0; i < c[ ca ].s.length; i++ )
      {
        c[ ca ].s[ i ] = c[ ca ].s[ i ] || c[ cb ].s[ i ];
      }

      delete ( cb );
      // out ( "draw " + b2i ( c[ ca ].s ).length );
      c[ ca ].draw ();
    }
  }

/*--------------------------------------------------------------------
 cluster2 - Original method to cluster the matrix
--------------------------------------------------------------------*/

  public void cluster2 ( Matrix m, String metric, String type )
  {
    Matrix r = new Matrix ( 3, m.nc - 1 );

    Matrix d = new Matrix ( m.nc, m.nc );
    d.add ( 1 );

    if ( metric.equals ( "pearson" ) )
    {
      d.sub ( m.pearson () );
    }

    else
    {
      d.sub ( m.spearman () );
    }

    Matrix d2 = d.copy ();

    double max = d.max () + 1;

    Matrix y = new Matrix ( m.nc - 1, 4 );
    int[] ind = seq ( 0, m.nc - 1 );

    int p = 0;
    int[] z = new int[2];

    out ( "start cluster" );
    for ( int i = 0; i < m.nc - 1; i++ )
    {
      double min = d.nonZeroMin ( z );
      int a = z[ 0 ];
      int b = z[ 1 ];

      y.x[ p ][ 0 ] = ind[ a ];
      y.x[ p ][ 1 ] = ind[ b ];
      y.x[ p ][ 2 ] = min;

      int newid = m.nc + p;
      int aid = ind[ a ];
      int bid = ind[ b ];

      for ( int j = 0; j < m.nc; j++ )
      {
        if ( ind[ j ] == aid || ind[ j ] == bid )
        {
          ind[ j ] = newid;
        }
      }

      for ( int j = 0; j < m.nc; j++ )
      {
        if ( ind[ j ] == newid )
        {
          // COMPLETE LINKAGE - MAXIMUM DISTANCE
          if ( type.equals ( "complete" ) )
          {
            for ( int k = 0; k < m.nc; k++ )
            {
              if ( ind[ k ] != newid && d2.x[ j ][ k ] > 0 && d2.x[ j ][ k ] > d.x[ b ][ k ] )
              {
                d.x[ b ][ k ] = d2.x[ j ][ k ];
              }
            }
          }

          // SINGLE LINKAGE - MINIMUM DISTANCE
          else
          {
            for ( int k = 0; k < m.nc; k++ )
            {
              if ( ind[ k ] != newid && d2.x[ j ][ k ] > 0 && d2.x[ j ][ k ] < d.x[ b ][ k ] )
              {
                d.x[ b ][ k ] = d2.x[ j ][ k ];
              }
            }
          }
        }
      }

      for ( int j = 0; j < m.nc; j++ )
      {
        if ( ind[ j ] == newid && j != b )
        {
          d.zeroCol ( j );
          d.zeroRow ( j );
        }
      }

      p++;
    }

    out ( "end cluster" );
    Texts[] ord = new Texts[m.nc * 2];
    Texts top = new Texts ();

    for ( int i = 0; i < m.nc; i++ )
    {
      ord[ i ] = new Texts ();
      ord[ i ].rpush ( i );
    }

    p = 0;
    for ( int i = 0; i < m.nc - 1; i++ )
    {
      int a = ( int ) ( y.x[ i ][ 0 ] );
      int b = ( int ) ( y.x[ i ][ 1 ] );
      int newid = m.nc + p;

      ord[ newid ] = new Texts ();
      ord[ newid ].rpush ( ord[ a ] );

      int a1 = ord[ a ].get ( 0 ).num ();
      int a2 = ord[ a ].get ( ord[ a ].size - 1 ).num ();
      int b1 = ord[ b ].get ( 0 ).num ();
      int b2 = ord[ b ].get ( ord[ b ].size - 1 ).num ();

      double dis1 = d2.x[ a1 ][ b1 ];
      double dis2 = d2.x[ a1 ][ b2 ];
      double dis3 = d2.x[ a2 ][ b1 ];
      double dis4 = d2.x[ a2 ][ b2 ];

      if ( dis1 < dis2 && dis1 < dis3 && dis1 < dis4 )
      {
        ord[ b ].reverse ();
        ord[ newid ].lpush ( ord[ b ] );
      }

      else if ( dis2 < dis3 && dis2 < dis4 )
      {
        ord[ newid ].lpush ( ord[ b ] );
      }

      else if ( dis3 < dis4 )
      {
        ord[ newid ].rpush ( ord[ b ] );
      }

      else
      {
        ord[ b ].reverse ();
        ord[ newid ].rpush ( ord[ b ] );
      }

      ord[ a ] = new Texts ();
      ord[ b ] = new Texts ();

      top = ord[ newid ];
      p++;
    }

    this.m = y;
    order = top.num ();
  }

	public int[] getOrder() {
		return this.order;
	}
}

