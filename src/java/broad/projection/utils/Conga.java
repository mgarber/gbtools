package broad.projection.utils;


/*--------------------------------------------------------------------
 com.scanfeld.cluster.Conga
 ---------------------------------------------------------------------
 Creates the conga data structure from:
 Fast Hierarchical Clustering and Other Applications of Dynamic Closest Pairs
 David Eppstein
 ---------------------------------------------------------------------
 Author: Daniel Scanfeld
   Date: 08.31.06
 ---------------------------------------------------------------------

--------------------------------------------------------------------*/

public class Conga extends Inc
{
/*--------------------------------------------------------------------
Class Variables
--------------------------------------------------------------------*/

  public boolean[] s;
  public Texts pt;
  public Texts dis;

  public Tree setTree;
  public int setCount;
  public int pointCount;
  public Tree path;

  public ClusterHierarchical c;

/*--------------------------------------------------------------------
 CONSTRUCTOR FUNCTIONS
 ---------------------------------------------------------------------
 Constructors for the possible text objects
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 com.scanfeld.cluster.Conga - Null constructor
--------------------------------------------------------------------*/

  public Conga ()
  {
  }

/*--------------------------------------------------------------------
 com.scanfeld.cluster.Conga - With a matrix, metric, and type
--------------------------------------------------------------------*/

  public Conga ( ClusterHierarchical c )
  {
    this.c = c;
    s = booleanFalse ( c.nc * 2 );
  }

/*--------------------------------------------------------------------
 setAll - Set to all points
--------------------------------------------------------------------*/

  public void setAll ()
  {
    for ( int i = c.nc; i < c.nc * 2; i++ )
    {
      s[ i ] = false;
    }

    for ( int i = 0; i < c.nc; i++ )
    {
      s[ i ] = true;
    }
  }

/*--------------------------------------------------------------------
 setOne - Set to one point
--------------------------------------------------------------------*/

  public void setOne ( int j )
  {
    for ( int i = 0; i < c.nc * 2; i++ )
    {
      s[ i ] = false;
    }

    s[ j ] = true;
  }

/*--------------------------------------------------------------------
 draw - Draw a line using the set
--------------------------------------------------------------------*/

  public void draw ()
  {
    pt = new Texts ();
    dis = new Texts ();

    int[] setArr = b2i ( s );
    setCount = setArr.length;
    pointCount = c.nc;

    setTree = i2tr ( setArr );
    path = new Tree ();

    int k = setArr[ 0 ];
    int n = -1;

    path.set ( k, 1 );
    pt.rpush ( k );
    setCount--;
    pointCount--;

    if ( setCount == 0 && setTree.get ( k ) == 1 )
    {
      if ( pointCount > 0 )
      {
        n = closest ( k );
      }
    }

    while ( setCount > 0 && k >= 0 )
    {
      if ( setTree.get ( k ) == 1 )
      {
        n = closest ( k );
      }

      else
      {
        n = closestIn ( k );
      }

      k = n;

      if ( setCount == 0 && setTree.get ( k ) == 1 )
      {
        if ( pointCount > 0 )
        {
          n = closest ( k );
        }
      }
    }

    setCount = setArr.length;
  }

/*--------------------------------------------------------------------
 closest - Get the closest not in the path
--------------------------------------------------------------------*/

  public int closest ( int a )
  {
    int k = -1;
    double d = -1;
    double td = 0;

    for ( int i = 0; i < c.setf.length; i++ )
    {
      if ( !( !c.setf[ i ] || a == i || path.get ( i ) == 1 ) )
      {
        td = c.dis ( a, i );

        if ( td < d || d < 0 )
        {
          d = td;
          k = i;
        }
      }
    }

    if ( k >= 0 )
    {
      pt.rpush ( k );
      dis.rpush ( d );
      path.set ( k, 1 );

      if ( setTree.get ( k ) == 1 )
      {
        setCount--;
      }
      pointCount--;
    }

    return k;
  }

/*--------------------------------------------------------------------
 closestIn - Get the closest in the set, not in the path
--------------------------------------------------------------------*/

  public int closestIn ( int a )
  {
    int k = -1;
    double d = -1;
    double td = 0;

    for ( int i = 0; i < c.setf.length; i++ )
    {
      if ( !( !c.setf[ i ] || a == i || setTree.get ( i ) != 1 || path.get ( i ) == 1 ) )
      {
        td = c.dis ( a, i );

        if ( td < d || d < 0 )
        {
          d = td;
          k = i;
        }
      }
    }

    if ( k >= 0 )
    {
      pt.rpush ( k );
      dis.rpush ( d );
      path.set ( k, 1 );

      if ( setTree.get ( k ) == 1 )
      {
        setCount--;
      }
      pointCount--;
    }

    return k;
  }

/*--------------------------------------------------------------------
 merge - Merge two sets
--------------------------------------------------------------------*/

  public void merge ()
  {
  }

/*--------------------------------------------------------------------
 remove - Remove an item
--------------------------------------------------------------------*/

  public int remove ( int k )
  {
    int ind = pt.find ( k );

    if ( ind >= 0 )
    {
      int ind2 = ind + 1;
      while ( ind2 < pt.size && setTree.get ( pt.get ( ind2 ) ) != 1 )
      {
        ind2++;
      }

      if ( ind > 1 )
      {
        Conga g = add ( 0, ind - 1, ind - 1 );
      }

      if ( ind2 < pt.size )
      {
        Conga g = add ( ind2, pt.size - 1, -1 );

        if ( ind2 == pt.size - 1 )
        {
          g.draw ();
        }
      }

      if ( ind > 0 )
      {
        if ( setTree.get ( pt.get ( ind - 1 ) ) != 1 )
        {
          ind = 0;
        }
      }
    }

    return ind;
  }

/*--------------------------------------------------------------------
 add - add a new com.scanfeld.cluster.Conga
--------------------------------------------------------------------*/

  public Conga add ( int a, int b, int x )
  {
    int[] psel = seq ( a, b );
    int[] dsel = seq ( a, b - 1 );

    Conga g = c.next ();
    g.pt = pt.select ( psel );
    g.dis = dis.select ( dsel );

    g.path = new Tree ();
    g.setTree = new Tree ();
    g.setCount = 0;

    for ( int i = 0; i < g.pt.size; i++ )
    {
      int j = g.pt.get ( i ).num ();
      g.path.set ( j, 1 );

      if ( setTree.get ( j ) == 1 && i != x )
      {
        g.s[ j ] = true;
        g.setCount++;
        g.setTree.set ( j, 1 );
      }
    }

    return g;
  }
}
