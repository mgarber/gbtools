package broad.projection.utils;

/*--------------------------------------------------------------------
 ColorBlend
 ---------------------------------------------------------------------
 Blends a set of colors
 ---------------------------------------------------------------------
 Author: Daniel Scanfeld
   Date: 08.23.06
 ---------------------------------------------------------------------

--------------------------------------------------------------------*/

public class ColorBlend extends Inc
{
/*--------------------------------------------------------------------
 Class Variables
--------------------------------------------------------------------*/

  public double[] r;  // Reds
  public double[] g;  // Greens
  public double[] b;  // Blues
  public double[] p;  // Percentage
  public int numCol;  // Number of colors

/*--------------------------------------------------------------------
 CONSTRUCTOR FUNCTIONS
 ---------------------------------------------------------------------
 Initialize ColorBlend objects
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 ColorBlend - Null constructor
--------------------------------------------------------------------*/

  public ColorBlend ()
  {
    numCol = 0;
    r = new double[100];
    g = new double[100];
    b = new double[100];
    p = new double[100];
  }

/*--------------------------------------------------------------------
 ColorBlend - Null constructor
--------------------------------------------------------------------*/

  public ColorBlend ( double rv, double gv, double bv )
  {
    numCol = 0;
    r = new double[100];
    g = new double[100];
    b = new double[100];
    p = new double[100];
    add ( rv, gv, bv, 0.0 );
  }

/*--------------------------------------------------------------------
 ColorBlend - Null constructor
--------------------------------------------------------------------*/

  public ColorBlend ( double v[] )
  {
    numCol = 0;
    r = new double[100];
    g = new double[100];
    b = new double[100];
    p = new double[100];
    add ( v[ 0 ], v[ 1 ], v[ 2 ], 0.0 );
  }

/*--------------------------------------------------------------------
 ColorBlend - String constructor
--------------------------------------------------------------------*/

  public ColorBlend ( String color )
  {
    ColorBlend cb;

    if ( color.equals ( "redwhite" ) )
    {
      cb = new ColorBlend ( 1.0, 1.0, 1.0 );
      cb.add ( 214 / 255.0, 12 / 255.0, 0 / 255.0, 1.0 );
    }

    else if ( color.equals ( "blackwhite" ) )
    {
      cb = new ColorBlend ( 1.0, 1.0, 1.0 );
      cb.add ( 0, 0, 0, 1.0 );
    }

    else if ( color.equals ( "whiteblack" ) )
    {
      cb = new ColorBlend ( 0.0, 0.0, 0.0 );
      cb.add ( 1, 1, 1, 1.0 );
    }

    else if ( color.equals ( "rainbow" ) )
    {
      double hue = 0;
      double twothirds = 50 * ( 2.0 / 3.0 );
      cb = new ColorBlend ();

      for ( int i = 0; i < 50; i++ )
      {
        if ( i >= twothirds )
        {
          hue = ( 0.5 - 0.5 * ( double ) ( i - twothirds ) / ( 50.0 - twothirds ) );
        }

        else
        {
          hue = ( 1.0 - 0.5 * ( double ) ( i ) / twothirds );
        }

        boolean pos = ( hue >= 0.5 );
        hue = 1.0 - pow ( ( 1.0 - Math.abs ( ( hue - 0.5 ) * 2 ) ), 2 );

        if ( !pos )
        {
          hue = -hue;
        }

        hue = hue / 2.0 + 0.5;

        hue = hue * ( 2.0 / 3.0 );
        cb.addHSV ( hue, 1.0, 1.0, ( double ) ( i ) / 49.0 );
      }
    }

    else
    {
      cb = new ColorBlend ( 69 / 255.0, 0 / 255.0, 173 / 255.0 );
      cb.add ( 39 / 255.0, 0 / 255.0, 209 / 255.0, 1.0 / 11.0 );
      cb.add ( 107 / 255.0, 88 / 255.0, 239 / 255.0, 2.0 / 11.0 );
      cb.add ( 136 / 255.0, 136 / 255.0, 255 / 255.0, 3.0 / 11.0 );
      cb.add ( 199 / 255.0, 193 / 255.0, 255 / 255.0, 4.0 / 11.0 );
      cb.add ( 213 / 255.0, 213 / 255.0, 255 / 255.0, 5.0 / 11.0 );
      cb.add ( 255 / 255.0, 192 / 255.0, 229 / 255.0, 6.0 / 11.0 );
      cb.add ( 255 / 255.0, 137 / 255.0, 137 / 255.0, 7.0 / 11.0 );
      cb.add ( 255 / 255.0, 112 / 255.0, 128 / 255.0, 8.0 / 11.0 );
      cb.add ( 255 / 255.0, 90 / 255.0, 90 / 255.0, 9.0 / 11.0 );
      cb.add ( 239 / 255.0, 64 / 255.0, 64 / 255.0, 10.0 / 11.0 );
      cb.add ( 214 / 255.0, 12 / 255.0, 0 / 255.0, 1.0 );
    }

    r = cb.r;
    g = cb.g;
    b = cb.b;
    p = cb.p;
    numCol = cb.numCol;
  }

/*--------------------------------------------------------------------
 add - Add a color
--------------------------------------------------------------------*/

  public void add ( double rv, double gv, double bv, double pv )
  {
    r[ numCol ] = rv;
    g[ numCol ] = gv;
    b[ numCol ] = bv;
    p[ numCol ] = pv;
    numCol++;
  }

/*--------------------------------------------------------------------
 add - Add an HSV color
--------------------------------------------------------------------*/

  public void addHSV ( double hv, double sv, double vv, double pv )
  {
    double[] h = hsv2rgb ( hv, sv, vv );
    r[ numCol ] = h[ 0 ];
    g[ numCol ] = h[ 1 ];
    b[ numCol ] = h[ 2 ];
    p[ numCol ] = pv;
    numCol++;
  }

/*--------------------------------------------------------------------
 add - Add a color
--------------------------------------------------------------------*/

  public double[] add ( double v[], double pv )
  {
    r[ numCol ] = v[ 0 ];
    g[ numCol ] = v[ 1 ];
    b[ numCol ] = v[ 2 ];
    p[ numCol ] = pv;
    numCol++;

    return null;
  }

/*--------------------------------------------------------------------
 get - Get a color
--------------------------------------------------------------------*/

  public double[] get ( double z )
  {
    if ( z < 0 )
    {
      z = 0;
    }

    if ( z > 1 )
    {
      z = 1;
    }

    int i = 0;
    for ( i = 0; i < numCol - 1; i++ )
    {
      if ( p[ i ] <= z && p[ i + 1 ] >= z )
      {
        break;
      }
    }

    double[] c = new double[3];
    z = ( z - p[ i ] ) / ( p[ i + 1 ] - p[ i ] );
    c[ 0 ] = r[ i ] + z * ( r[ i + 1 ] - r[ i ] );
    c[ 1 ] = g[ i ] + z * ( g[ i + 1 ] - g[ i ] );
    c[ 2 ] = b[ i ] + z * ( b[ i + 1 ] - b[ i ] );

    return c;
  }
}