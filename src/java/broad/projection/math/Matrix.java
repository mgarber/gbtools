package broad.projection.math;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import broad.projection.utils.ClusterHierarchical;
import broad.projection.utils.Inc;
import broad.projection.utils.Text;
import broad.projection.utils.Texts;


/*--------------------------------------------------------------------
 com.scanfeld.math.Matrix
 ---------------------------------------------------------------------
 com.scanfeld.math.Matrix operations
 ---------------------------------------------------------------------
 Author: Daniel Scanfeld
   Date: 08.31.06
 ---------------------------------------------------------------------
 mod all matrix operations
 make all add, etc... work with varying matrix sizes by repeats
 make all matrix operations work on their own cells
 finish the matrix stuff that stores in this, did add and sub
 hierarchical clustering
 http://www.elet.polimi.it/upload/matteucc/Clustering/tutorial_html/hierarchical.html
 merge matrices functions
 change all col functions to non-prefixed
--------------------------------------------------------------------*/

public class Matrix extends Inc {
	private static final double e = Math.pow(10, -14);  // Save the minimal double
	private static final double SMALL_NUMBER = Math.pow(10, -6);
	public static final NMFCostFunction DEFAULT_COST_FUNCTION = new PoissonNMFCost();
/*--------------------------------------------------------------------
Class Variables
--------------------------------------------------------------------*/

    public double[][] x;
    private double  precision = SMALL_NUMBER;
    public int nr;
    public int nc;
    public double[] err;
    public int erri;
    public List<Double> errors;
    public Matrix W;
    public Matrix Wi;
    public Matrix H;
    public Matrix C;

/*--------------------------------------------------------------------
 CONSTRUCTOR FUNCTIONS
 ---------------------------------------------------------------------
 Constructors for the possible text objects
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 Matrix - Null constructor
--------------------------------------------------------------------*/

    public Matrix() {
        init(3, 3);  // Initialize with a 3 by 3 matrix
    }

/*--------------------------------------------------------------------
 Matrix - Width and height constructor
--------------------------------------------------------------------*/

    public Matrix(int nr, int nc) {
        init(nr, nc);
    }

/*--------------------------------------------------------------------
 Matrix - Matrix constructor
--------------------------------------------------------------------*/

    public Matrix(Matrix a) {
        init(a.nr, a.nc);
        set(a);
    }

/*--------------------------------------------------------------------
 Matrix - Matrix constructor
--------------------------------------------------------------------*/

    public Matrix(double[][] a) {
        init(a.length, a[0].length);
        set(a);
    }

/*--------------------------------------------------------------------
 Matrix - Matrix constructor
--------------------------------------------------------------------*/

    public Matrix(double[] a) {
        init(1, a.length);
        set(a);
    }

/*--------------------------------------------------------------------
 hypersphere - Project to the hypersphere
--------------------------------------------------------------------*/

    public void hypersphere() {
        hypersphere(1.0);
    }

/*--------------------------------------------------------------------
 hypersphere - Project to the hypersphere
--------------------------------------------------------------------*/

    public void hypersphere(double max) {
        double[] e = euclidean();

        for (int j = 0; j < nr; j++) {
            for (int k = 0; k < nc; k++) {
                x[j][k] /= e[k];
            }
        }

        mul(max);
    }

    /*--------------------------------------------------------------------
    clusterHierarchical - Hierarchically cluster the columns
   --------------------------------------------------------------------*/

       public ClusterHierarchical clusterHierarchical(String metric, String type) {
           return new ClusterHierarchical(this, metric, type);
       }

    
 /*--------------------------------------------------------------------
 init - Initialize the matrix
--------------------------------------------------------------------*/

    public void init(int nr, int nc) {
        this.nr = 0;             // Start with zero rows
        this.nc = 0;             // Start with zero columns
        this.resize(nr, nc);  // Resize the matrix
    }

/*--------------------------------------------------------------------
 inverse - Find the inverse
--------------------------------------------------------------------*/

    public void inverse() {
        set(ninverse());
    }

/*--------------------------------------------------------------------
 ninverse - Find the inverse
--------------------------------------------------------------------*/

    public Matrix ninverse() {
        Matrix b;

        if (nr > nc && 1==0) {  // Makes it rank deficient
            b = new Matrix(nr, nc);
            //weka.core.matrix.Matrix m = new weka.core.matrix.Matrix(ntranspose().x);
           broad.projection.weka.matrix.Matrix m = new broad.projection.weka.matrix.Matrix(x);
            m.inverse();
            
            b.x = m.inverse().getArray();
            //b.transpose();
        } else {
            b = new Matrix(nc, nr);
            broad.projection.weka.matrix.Matrix m = new broad.projection.weka.matrix.Matrix(x);
            m.inverse();
            b.x = m.inverse().getArray();
        }

        b.nr = b.x.length;
        b.nc = b.x[0].length;
        return b;                                 // Return the new matrix
    }

/*--------------------------------------------------------------------
 merge - Merge two matrices
--------------------------------------------------------------------*/

    public Matrix mergeRight(Matrix a) {
        Matrix b = new Matrix(nr, nc + a.nc);  // Create the new matrix
        b.set(this, 0, 0);                     // Add this matrix
        b.set(a, 0, nc);                       // Add a to the right
        return b;                                 // Return the new matrix
    }

/*--------------------------------------------------------------------
 set - Set to a matrix
--------------------------------------------------------------------*/

    public void set(Matrix a, int r, int c) {
        // LOOP THROUGH THE ROWS
        for (int i = 0; i < a.nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < a.nc; j++) {
                x[r + i][c + j] = a.x[i][j];
            }
        }
    }

/*--------------------------------------------------------------------
 set - Set to a matrix
--------------------------------------------------------------------*/

    public void set(double[][] a, int r, int c) {
        // LOOP THROUGH THE ROWS
        for (int i = 0; i < a.length; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < a[0].length; j++) {
                x[r + i][c + j] = a[i][j];
            }
        }
    }

    
    
/*--------------------------------------------------------------------
 standardCorrelation - Correlation between two columns
--------------------------------------------------------------------*/

    public double standardCorrelation(int a, int b) {
        double c = 0;
        double d = 0;
        double e = 0;

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            c += x[i][a] * x[i][b];
            d += x[i][a] * x[i][a];
            e += x[i][b] * x[i][b];
        }

        return c / pow(d * e, 0.5);
    }

/*--------------------------------------------------------------------
 standardCorrelation - Correlation between all columns
--------------------------------------------------------------------*/

    public Matrix standardCorrelation() {
        Matrix c = new Matrix(nc, nc);
        Matrix d = new Matrix(nc, nc);
        Matrix e = new Matrix(nc, nc);
        Matrix mean = colMean();

        // LOOP THROUGH THE COLUMNS
        for (int a = 0; a < nc; a++) {
            // LOOP THROUGH THE COLUMNS
            for (int b = 0; b < nc; b++) {
                // LOOP THROUGH THE ROWS
                for (int i = 0; i < nr; i++) {
                    c.x[a][b] += x[i][a] * x[i][b];
                    d.x[a][b] += x[i][a] * x[i][a];
                    e.x[a][b] += x[i][b] * x[i][b];
                }
            }
        }

        d.mul(e);
        d.pow(0.5);
        c.div(d);

        return c;
    }

/*--------------------------------------------------------------------
 pearson - Correlation between two columns
--------------------------------------------------------------------*/

    public Matrix pearson(Matrix n) {
        Matrix m = new Matrix(nc, 1);

        for (int i = 0; i < nc; i++) {
            m.x[i][0] = pearson(i, n);
        }

        return m;
    }

/*--------------------------------------------------------------------
 pearson - Correlation between two columns
--------------------------------------------------------------------*/

    public double pearson(int a, Matrix n) {
        double c = 0;
        double d = 0;
        double e = 0;

        double ma = colMean(a);
        double mb = n.colMean(0);

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            c += (x[i][a] - ma) * (n.x[i][0] - mb);
            d += (x[i][a] - ma) * (x[i][a] - ma);
            e += (n.x[i][0] - mb) * (n.x[i][0] - mb);
        }

        return c / pow(d * e, 0.5);
    }

/*--------------------------------------------------------------------
 pearson - Correlation between two columns
--------------------------------------------------------------------*/

    public double pearson(int a, int b) {
        double c = 0;
        double d = 0;
        double e = 0;

        double ma = colMean(a);
        double mb = colMean(b);

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            c += (x[i][a] - ma) * (x[i][b] - mb);
            d += (x[i][a] - ma) * (x[i][a] - ma);
            e += (x[i][b] - mb) * (x[i][b] - mb);
        }

        return c / pow(d * e, 0.5);
    }

/*--------------------------------------------------------------------
 pearson - Correlation between all columns
--------------------------------------------------------------------*/

    public Matrix pearson() {
        double ma, mb;
        Matrix c = new Matrix(nc, nc);

        for (int i = 0; i < nc; i++) {
            for (int j = 0; j < nc; j++) {
                c.x[i][j] = pearson(i, j);
            }
        }

        /*
        // eventually fix this up

        Matrix d = new Matrix ( nc, nc );
        Matrix e = new Matrix ( nc, nc );
        Matrix mean = colMean ();

        // LOOP THROUGH THE COLUMNS
        for ( int a = 0; a < nc; a++ )
        {
          // LOOP THROUGH THE COLUMNS
          for ( int b = a + 1; b < nc; b++ )
          {
            // LOOP THROUGH THE ROWS
            for ( int i = 0; i < nr; i++ )
            {
              ma = mean.x[ 0 ][ a ];
              mb = mean.x[ 0 ][ b ];

              c.x[ a ][ b ] += ( x[ i ][ a ] - ma ) * ( x[ i ][ b ] - mb );
              d.x[ a ][ b ] += ( x[ i ][ a ] - ma ) * ( x[ i ][ a ] - ma );
              e.x[ a ][ b ] += ( x[ i ][ b ] - mb ) * ( x[ i ][ b ] - mb );
            }

            c.x[ b ][ a ] = c.x[ a ][ b ];
            d.x[ b ][ a ] = d.x[ a ][ b ];
            e.x[ b ][ a ] = e.x[ a ][ b ];
          }
        }

        d.mul ( e );
        d.pow ( 0.5 );
        c.div ( d );
        */

        return c;
    }
    
    
    /*--------------------------------------------------------------------
    KLDivergence - Kulbeck-Leibler divergence between two columns
   --------------------------------------------------------------------*/

       public double KLDivergence(int a, int b) {
           double[] v=this.getCol(a);
           double[] w=this.getCol(b);
    	   
           double rtrn=0;
           for(int i=0; i<v.length; i++){
        	   rtrn+=(v[i]*Math.log(v[i]/w[i]))-v[i]+w[i];
           }
           
    	  return rtrn;
       }
    
    /*--------------------------------------------------------------------
    KLDivergence - Kulbeck-Leibler divergence between columns
   --------------------------------------------------------------------*/

       public Matrix KLDivergence() {
           double ma, mb;
           Matrix c = new Matrix(nc, nc);

           for (int i = 0; i < nc; i++) {
               for (int j = 0; j < nc; j++) {
                   c.x[i][j] = Math.min(KLDivergence(i, j), KLDivergence(j,i));
               }
           }

           /*
           // eventually fix this up

           Matrix d = new Matrix ( nc, nc );
           Matrix e = new Matrix ( nc, nc );
           Matrix mean = colMean ();

           // LOOP THROUGH THE COLUMNS
           for ( int a = 0; a < nc; a++ )
           {
             // LOOP THROUGH THE COLUMNS
             for ( int b = a + 1; b < nc; b++ )
             {
               // LOOP THROUGH THE ROWS
               for ( int i = 0; i < nr; i++ )
               {
                 ma = mean.x[ 0 ][ a ];
                 mb = mean.x[ 0 ][ b ];

                 c.x[ a ][ b ] += ( x[ i ][ a ] - ma ) * ( x[ i ][ b ] - mb );
                 d.x[ a ][ b ] += ( x[ i ][ a ] - ma ) * ( x[ i ][ a ] - ma );
                 e.x[ a ][ b ] += ( x[ i ][ b ] - mb ) * ( x[ i ][ b ] - mb );
               }

               c.x[ b ][ a ] = c.x[ a ][ b ];
               d.x[ b ][ a ] = d.x[ a ][ b ];
               e.x[ b ][ a ] = e.x[ a ][ b ];
             }
           }

           d.mul ( e );
           d.pow ( 0.5 );
           c.div ( d );
           */

           return c;
       }

/*--------------------------------------------------------------------
 spearman - Correlation between all columns
--------------------------------------------------------------------*/

    public Matrix spearman() {
        return rank().pearson();
    }


/*--------------------------------------------------------------------
 select - Select part of the array
--------------------------------------------------------------------*/

    public Matrix select(int[] rs, int[] cs) {
        Matrix m = new Matrix(rs.length, cs.length);

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < m.nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < m.nc; j++) {
                m.x[i][j] = x[rs[i]][cs[j]];  // Copy the i,j item
            }
        }

        return m;
    }

/*--------------------------------------------------------------------
 selectRow - Select part of the array
--------------------------------------------------------------------*/

    public Matrix selectRow(int[] rs) {
        Matrix m = new Matrix(rs.length, nc);

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < m.nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < m.nc; j++) {
                m.x[i][j] = x[rs[i]][j];  // Copy the i,j item
            }
        }

        return m;
    }

/*--------------------------------------------------------------------
 sortColInd - Return a column's indices
--------------------------------------------------------------------*/

    public int[] sortColInd(int c) {
        double[] nc = new double[nr];

        for (int i = 0; i < nr; i++) {
            nc[i] = x[i][c];
        }

        return sort(nc);
    }

/*--------------------------------------------------------------------
 selectCol - Select part of the array
--------------------------------------------------------------------*/

    public Matrix selectCol(int c) {
        return selectCol(new int[]{c});
    }

/*--------------------------------------------------------------------
 selectCol - Select part of the array
--------------------------------------------------------------------*/

    public Matrix selectCol(int[] cs) {
        Matrix m = new Matrix(nr, cs.length);

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < m.nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < m.nc; j++) {
                m.x[i][j] = x[i][cs[j]];  // Copy the i,j item
            }
        }

        return m;
    }

/*--------------------------------------------------------------------
 set - Set to a matrix
--------------------------------------------------------------------*/

    public void set(Matrix m) {
        if (nr != m.nr || nc != m.nc) {
            nr = m.nr;               // Copy row count
            nc = m.nc;               // Copy column count
            x = new double[nr][nc];  // Create the space for the matrix
        }

        // LOOP THROUGH THE ROWS
        {
            for (int i = 0; i < nr; i++) {
                // LOOP THROUGH THE COLUMNS
                for (int j = 0; j < nc; j++) {
                    x[i][j] = m.x[i][j];  // Copy the i,j item
                }
            }
        }
    }

/*--------------------------------------------------------------------
 set - Set to a matrix
--------------------------------------------------------------------*/

    public void set(double[][] a) {
        if (nr != a.length || nc != a[0].length) {
            nr = a.length;               // Copy row count
            nc = a[0].length;               // Copy column count
            x = new double[nr][nc];  // Create the space for the matrix
        }

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < a.length; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < a[0].length; j++) {
                x[i][j] = a[i][j];
            }
        }
    }

/*--------------------------------------------------------------------
 set - Set to a matrix
--------------------------------------------------------------------*/

    public void set(double[] a) {
        if (nr != 1 || nc != a.length) {
            nr = 1;                   // Copy row count
            nc = a.length;           // Copy column count
            x = new double[nr][nc];  // Create the space for the matrix
        }

        // LOOP THROUGH THE COLUMNS
        for (int j = 0; j < a.length; j++) {
            x[0][j] = a[j];
        }
    }

/*--------------------------------------------------------------------
 dup - Copy a matrix
--------------------------------------------------------------------*/

    public Matrix dup() {
        return new Matrix(this);
    }

/*--------------------------------------------------------------------
 copy - Copy a matrix
--------------------------------------------------------------------*/

    public Matrix copy() {
        return new Matrix(this);
    }

/*--------------------------------------------------------------------
 resize - Resize the matrix
--------------------------------------------------------------------*/

    public void resize(int nr, int nc) {
        // IF THE NEW SIZE IS DIFFERENT THAN THE OLD SIZE
        if ((nr != (this.nr)) || (nc != (this.nc))) {
            this.nr = nr;                // Set the new nr
            this.nc = nc;                // Set the new nc
            x = new double[nr][nc];  // Create a new matrix
        }
    }

/*--------------------------------------------------------------------
 out - Output the matrix
--------------------------------------------------------------------*/

    public void out() {
        Text t = new Text();  // Create an ouptut line

        // Output the size of the matrix
        out("Matrix ( " + nr + " x " + nc + " )");

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            t.set("");  // Empty the string

            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                // IF NOT THE FIRST ELEMENT
                if (j != 0) {
                    t.rpush("\t");  // Output a tab
                }

                t.rpush(x[i][j]);  // Output the element
            }

            out(t);  // Output the line
        }

        out("");  // Add a blank line at the end
    }

/*--------------------------------------------------------------------
 write - Output the matrix
--------------------------------------------------------------------*/

    public void write(String fname) throws IOException{
        FileWriter writer = new FileWriter(fname);
        Text t = new Text();

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            t.set("");  // Empty the string

            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                // IF NOT THE FIRST ELEMENT
                if (j != 0) {
                    t.rpush("\t");  // Output a tab
                }

                t.rpush(x[i][j]);  // Output the element
            }

            writer.write(t+"\n");  // Output the line
        }
        writer.close();
    }
    
  
    
    public void write(String fname, ArrayList genes) throws IOException{
        FileWriter writer = new FileWriter(fname);
        Text t = new Text();

        writer.write("PID\tName");
        for(Object gene: genes){
        	writer.write("\t"+gene.toString());
        }
        writer.write("\n");
        
        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
        	writer.write(genes.get(i)+"\t"+genes.get(i)+"\t");
            t.set("");  // Empty the string

            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                // IF NOT THE FIRST ELEMENT
                if (j != 0) {
                    t.rpush("\t");  // Output a tab
                }

                t.rpush(x[i][j]);  // Output the element
            }

            writer.write("F"+i+"\t"+t.toString()+"\n");  // Output the line
        }
        writer.close();
    }

/*--------------------------------------------------------------------
pow - Take each element to the nth power
--------------------------------------------------------------------*/

    public void pow(double n) {
        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                x[i][j] = pow(x[i][j], n);  // Take the element to the nth power
            }
        }
    }

/*--------------------------------------------------------------------
apow - Take each element to the nth power
--------------------------------------------------------------------*/

    public void apow(double n) {
        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                x[i][j] = apow(x[i][j], n);  // Take the element to the nth power
            }
        }
    }

/*--------------------------------------------------------------------
rank - Create a rank matrix
--------------------------------------------------------------------*/

    public Matrix rank() {
        Matrix a = new Matrix(nr, nc);
        double z[] = new double[nr];

        // LOOP THROUGH THE COLUMNS
        for (int j = 0; j < nc; j++) {
            // LOOP THROUGH THE ROWS
            for (int i = 0; i < nr; i++) {
                z[i] = x[i][j];  // Get the jth item
            }

            int[] o = sort(z);

            // LOOP THROUGH THE ROWS
            for (int i = 0; i < nr; i++) {
                a.x[o[i]][j] = i + 1;
            }
        }

        return a;
    }

/*--------------------------------------------------------------------
rank - Create a rank matrix
--------------------------------------------------------------------*/

    public Matrix rank(double max) {
        Matrix a = rank();
        a.sub(1);
        a.mul(max / (nr - 1));
        return a;
    }

/*--------------------------------------------------------------------
add - Add a double
--------------------------------------------------------------------*/

    public void add(double a) {
        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                x[i][j] = x[i][j] + a;  // Add a to each element
            }
        }
    }

/*--------------------------------------------------------------------
add - Add a matrix to a double
--------------------------------------------------------------------*/

    public void add(Matrix b, double a) {
        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                x[i][j] = b.x[i][j] + a;  // Add a to each element
            }
        }
    }

/*--------------------------------------------------------------------
nlog - Take the log
--------------------------------------------------------------------*/

    public Matrix nlog() {
        Matrix b = new Matrix(nr, nc);  // New matrix

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                b.x[i][j] = log(x[i][j]);  // Log
            }
        }

        return b;  // Return the new solution
    }

    
   
/*--------------------------------------------------------------------
nadd - Add a double
--------------------------------------------------------------------*/

    public Matrix nadd(double a) {
        Matrix b = new Matrix(nr, nc);  // New matrix

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                b.x[i][j] = x[i][j] + a;  // Add a to each element
            }
        }

        return b;  // Return the new solution
    }

/*--------------------------------------------------------------------
add - Add a matrix
--------------------------------------------------------------------*/

    public void add(Matrix a) {
        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                x[i][j] = x[i][j] + a.x[i][j];  // Add the element
            }
        }
    }

/*--------------------------------------------------------------------
add - Add two matrices
--------------------------------------------------------------------*/

    public void add(Matrix a, Matrix b) {
        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                x[i][j] = a.x[i][j] + b.x[i][j];  // Add the element
            }
        }
    }

/*--------------------------------------------------------------------
nadd - Add a matrix
--------------------------------------------------------------------*/

    public Matrix nadd(Matrix a) {
        Matrix b = new Matrix(nr, nc);  // new Matrix

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                b.x[i][j] = x[i][j] + a.x[i][j];  // Add the element
            }
        }

        return b;  // Return the new solution
    }

/*--------------------------------------------------------------------
sub - Subtract a double
--------------------------------------------------------------------*/

    public void sub(double a) {
        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                x[i][j] = x[i][j] - a;  // Subtract a from each element
            }
        }
    }

/*--------------------------------------------------------------------
sub - Subtract a double from a matrix
--------------------------------------------------------------------*/

    public void sub(Matrix b, double a) {
        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                x[i][j] = b.x[i][j] - a;  // Subtract a from each element
            }
        }
    }

/*--------------------------------------------------------------------
nsub - Subtract a double
--------------------------------------------------------------------*/

    public Matrix nsub(double a) {
        Matrix b = new Matrix(nr, nc);  // New matrix

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                b.x[i][j] = x[i][j] - a;  // Subtract a from each element
            }
        }

        return b;  // Return the new solution
    }

/*--------------------------------------------------------------------
sub - Subtract a matrix
--------------------------------------------------------------------*/

    public void sub(Matrix a) {
        // IF THE MATRICES ARE OF DIFFERENT SIZES
        if (nr != a.nr || nc != a.nc) {
            // LOOP THROUGH THE ROWS
            for (int i = 0; i < nr; i++) {
                // LOOP THROUGH THE COLUMNS
                for (int j = 0; j < nc; j++) {
                    x[i][j] = x[i][j] - a.x[i % a.nr][j % a.nc];  // Subtract the element
                }
            }
        }

        // IF THE MATRICES ARE OF THE SAME SIZE
        else {
            // LOOP THROUGH THE ROWS
            for (int i = 0; i < nr; i++) {
                // LOOP THROUGH THE COLUMNS
                for (int j = 0; j < nc; j++) {
                    x[i][j] = x[i][j] - a.x[i][j];  // Subtract the element
                }
            }
        }
    }

/*--------------------------------------------------------------------
sub - Subtract a matrix from a matrix
--------------------------------------------------------------------*/

    public void sub(Matrix a, Matrix b) {
        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                x[i][j] = a.x[i][j] - b.x[i][j];  // Subtract the element
            }
        }
    }

/*--------------------------------------------------------------------
nsub - Subtract a matrix
--------------------------------------------------------------------*/

    public Matrix nsub(Matrix a) {
        Matrix b = new Matrix(nr, nc);  // new Matrix

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                b.x[i][j] = x[i][j] - a.x[i][j];  // Subtract the element
            }
        }

        return b;  // Return the new solution
    }

/*--------------------------------------------------------------------
mul - Multiply a double
--------------------------------------------------------------------*/

    public void mul(double a) {
        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                x[i][j] = x[i][j] * a;  // Multiply a to each element
            }
        }
    }

/*--------------------------------------------------------------------
nmul - Multiply a double
--------------------------------------------------------------------*/

    public Matrix nmul(double a) {
        Matrix b = new Matrix(nr, nc);  // New matrix

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                b.x[i][j] = x[i][j] * a;  // Multiply a to each element
            }
        }

        return b;  // Return the new solution
    }

/*--------------------------------------------------------------------
mul - Multiply a matrix
--------------------------------------------------------------------*/

    public void mul(Matrix a) {
        // IF THE MATRICES ARE OF DIFFERENT SIZES
        if (nr != a.nr || nc != a.nc) {
            // LOOP THROUGH THE ROWS
            for (int i = 0; i < nr; i++) {
                // LOOP THROUGH THE COLUMNS
                for (int j = 0; j < nc; j++) {
                    x[i][j] = x[i][j] * a.x[i % a.nr][j % a.nc];  // Multiply a
                }
            }
        }

        // IF THE MATRICES ARE OF THE SAME SIZE
        else {
            // LOOP THROUGH THE ROWS
            for (int i = 0; i < nr; i++) {
                // LOOP THROUGH THE COLUMNS
                for (int j = 0; j < nc; j++) {
                    x[i][j] = x[i][j] * a.x[i][j];  // Multiply a
                }
            }
        }
    }

/*--------------------------------------------------------------------
nmul - Multiply a matrix
--------------------------------------------------------------------*/

    public Matrix nmul(Matrix a) {
        Matrix b = new Matrix(nr, nc);  // new Matrix

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                b.x[i][j] = x[i][j] * a.x[i][j];  // Multiply
            }
        }

        return b;  // Return the new solution
    }

/*--------------------------------------------------------------------
max - Take the maximum
--------------------------------------------------------------------*/

    public void max(Matrix a) {
        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                x[i][j] = max(x[i][j], a.x[i][j]);  // Take the max
            }
        }
    }

/*--------------------------------------------------------------------
matmul - Matrix multiply a matrix
--------------------------------------------------------------------*/

    public void matmul(Matrix a) {
        set(nmatmul(a));  // Set to the new matrix
    }

/*--------------------------------------------------------------------
normalizeRow - Normalize the rows
--------------------------------------------------------------------*/

    public void normalizeRow() {
        Matrix max = rowMax();
        Matrix min = rowMin();
        max.sub(min);
        sub(min);
        div(max);
    }

/*--------------------------------------------------------------------
normalizeRowMean01 - Normalize the rows by the mean and between 0 and 1
--------------------------------------------------------------------*/

    public void normalizeRowMean01() {
        for (int i = 0; i < nr; i++) {
            double mean = mean(x[i]);
            double min = min(x[i]) - mean;
            double max = max(x[i]) - mean;

            for (int j = 0; j < nc; j++) {
                x[i][j] = x[i][j] - mean;

                if (x[i][j] > 0) {
                    x[i][j] /= max;
                } else {
                    x[i][j] /= -min;
                }

                x[i][j] = (x[i][j] * 0.5) + 0.5;
            }
        }
    }

/*--------------------------------------------------------------------
normalizeRowMean - Normalize the rows with mean and stdev
--------------------------------------------------------------------*/

    public void normalizeRowMean() {
        Matrix me = rowMean();
        Matrix sd = rowSD();
        sub(me.ntranspose());
        div(sd.ntranspose());
    }

/*--------------------------------------------------------------------
normalizeColMean - Normalize the cols with mean and stdev
--------------------------------------------------------------------*/

    public void normalizeColMean() {
        Matrix me = colMean();
        Matrix sd = colSD();
        sub(me);
        div(sd);
    }

/*--------------------------------------------------------------------
nonNegUnit - normalize to the unit matrix, remove negatives
--------------------------------------------------------------------*/

    public void nonNegUnit() {
        double amin = min();
        double amax = max();

        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nc; j++) {
                if (x[i][j] < 0) {
                    x[i][j] = 0;
                }

                x[i][j] = x[i][j] / amax;
            }
        }
    }

/*--------------------------------------------------------------------
unit - normalize to the unit matrix
--------------------------------------------------------------------*/

    public void unit() {
        double amin = min();
        double amax = max();

        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nc; j++) {
                x[i][j] = (x[i][j] - amin) / (amax - amin);
            }
        }
    }

/*--------------------------------------------------------------------
scale - scale to the unit matrix
--------------------------------------------------------------------*/

    public void scale() {
        double amax = max();

        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nc; j++) {
                x[i][j] = x[i][j] / amax;
            }
        }
    }

/*--------------------------------------------------------------------
scale - scale to the unit matrix * s
--------------------------------------------------------------------*/

    public void scale(double s) {
        double amax = max();

        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nc; j++) {
                x[i][j] = s * x[i][j] / amax;
            }
        }
    }

/*--------------------------------------------------------------------
nmatmul - Matrix multiply a matrix
--------------------------------------------------------------------*/

    public Matrix nmatmul(Matrix a) {
        // CREATE A NEW MATRIX TO STORE THE SOLUTION
        Matrix b = new Matrix(nr, a.nc);  // new Matrix

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < a.nc; j++) {
                b.x[i][j] = 0;  // Initialize to zero

                // LOOP THROUGH THE COLUMNS
                for (int k = 0; k < nc; k = k + 1) {
                    // MULTIPLY THE TWO ELEMENTS TOGETHER AND ADD THEM TO THE CURRENT ELEMENT IN B
                    b.x[i][j] = b.x[i][j] + (x[i][k] * a.x[k][j]);
                }
            }
        }

        return b;  // Return the solution
    }

/*--------------------------------------------------------------------
abs - Absolute value
--------------------------------------------------------------------*/

    public void abs() {
        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                x[i][j] = abs(x[i][j]);  // Absolute value
            }
        }
    }

/*--------------------------------------------------------------------
div - Divide a double
--------------------------------------------------------------------*/

    public void div(double a) {
        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                x[i][j] = x[i][j] / a;  // Divide a from each element
            }
        }
    }

/*--------------------------------------------------------------------
ndiv - Divide a double
--------------------------------------------------------------------*/

    public Matrix ndiv(double a) {
        Matrix b = new Matrix(nr, nc);  // New matrix

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                b.x[i][j] = x[i][j] / a;  // Divide a from each element
            }
        }

        return b;  // Return the new solution
    }

/*--------------------------------------------------------------------
div - Divide a matrix
--------------------------------------------------------------------*/

    public void div(Matrix a) {
        // IF THE MATRICES ARE OF DIFFERENT SIZES
        if (nr != a.nr || nc != a.nc) {
            // LOOP THROUGH THE ROWS
            for (int i = 0; i < nr; i++) {
                // LOOP THROUGH THE COLUMNS
                for (int j = 0; j < nc; j++) {
                    x[i][j] = x[i][j] / a.x[i % a.nr][j % a.nc];  // Divide a from the element
                }
            }
        }

        // IF THE MATRICES ARE OF THE SAME SIZE
        else {
            // LOOP THROUGH THE ROWS
            for (int i = 0; i < nr; i++) {
                // LOOP THROUGH THE COLUMNS
                for (int j = 0; j < nc; j++) {
                    x[i][j] = x[i][j] / a.x[i][j];  // Divide a from the element
                }
            }
        }
    }

/*--------------------------------------------------------------------
ndiv - Divide a matrix
--------------------------------------------------------------------*/

    public Matrix ndiv(Matrix a) {
        Matrix b = new Matrix(nr, nc);  // new Matrix

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                b.x[i][j] = x[i][j] / a.x[i][j];  // Divide a from  the element
            }
        }

        return b;  // Return the new solution
    }

/*--------------------------------------------------------------------
transpose - Transpose the elements in a matrix
--------------------------------------------------------------------*/

    public void t() {
        transpose();
    }

/*--------------------------------------------------------------------
ntranspose - Transpose the elements in a matrix
--------------------------------------------------------------------*/

    public Matrix nt() {
        return ntranspose();
    }

/*--------------------------------------------------------------------
transpose - Transpose the elements in a matrix
--------------------------------------------------------------------*/

    public void transpose() {
        Matrix b = new Matrix(nc, nr);  // new Matrix

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                // TRANSPOSE THE ELEMENT INTO B
                b.x[j][i] = x[i][j];
            }
        }

        set(b);  // Set to the new matrix
    }

/*--------------------------------------------------------------------
ntranspose - Transpose the elements in a matrix
--------------------------------------------------------------------*/

    public Matrix ntranspose() {
        Matrix b = new Matrix(nc, nr);  // new Matrix

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                // TRANSPOSE THE ELEMENT INTO B
                b.x[j][i] = x[i][j];
            }
        }

        return b;  // Return the transposed matrix
    }

/*--------------------------------------------------------------------
 threshold - Threshold the values
--------------------------------------------------------------------*/

    public void threshold(double min, double max) {
        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                // IF THE ELEMENT IS LESS THAN THE MINIMUM
                if (x[i][j] < min) {
                    x[i][j] = min;  // Store the minimum
                }

                // IF THE ELEMENT IS GREATER THAN THE MAXIMUM
                if (x[i][j] > max) {
                    x[i][j] = max;  // Store the maximum
                }
            }
        }
    }

/*--------------------------------------------------------------------
 foldChange - calculate the fold change
--------------------------------------------------------------------*/

    public Matrix foldChange() {
        Matrix min = rowMin();
        Matrix max = rowMax();
        max.div(min);
        max.abs();
        return max;
    }

/*--------------------------------------------------------------------
 variation - Variation filter
--------------------------------------------------------------------*/

    public int[] variation(double fold, double delta) {
        Texts ind = new Texts();
        Matrix min = rowMin();
        Matrix max = rowMax();

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // IF ABOVE THE FOLD CHANGE AND DELTA THRESHOLDS
            if ((max.x[i][0] - min.x[i][0] >= delta) && (abs(max.x[i][0] / min.x[i][0]) >= fold)) {
                ind.rpush(i);  // Push the index
            }
        }

        return ind.num();  // Return the indices
    }

/*--------------------------------------------------------------------
min - Find the minimum location in the matrix
--------------------------------------------------------------------*/

    public double min() {
        double y = x[0][0];  // Initial value

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                // IF THE ELEMENT IS LESS THAN THE CURRENT MAX
                if (x[i][j] < y) {
                    y = x[i][j];  // Store the new maximum
                }
            }
        }

        return y;  // Return the minimum
    }

/*--------------------------------------------------------------------
min - Find the minimum location in the matrix
--------------------------------------------------------------------*/

    public double minItem(int[] z) {
        z[0] = 0;
        z[1] = 0;
        double y = x[0][0];  // Initial value

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                // IF THE ELEMENT IS LESS THAN THE CURRENT MIN
                if (x[i][j] < y) {
                    y = x[i][j];  // Store the new maximum
                    z[0] = i;
                    z[1] = j;
                }
            }
        }

        return y;  // Return the minimum
    }

/*--------------------------------------------------------------------
min - Find the minimum location in the matrix
--------------------------------------------------------------------*/

    public double nonZeroMin() {
        double y = 0;  // Initial value

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                // IF THE ELEMENT IS LESS THAN THE CURRENT MAX
                if (((y == 0) || (x[i][j] < y)) && x[i][j] != 0) {
                    y = x[i][j];  // Store the new maximum
                }
            }
        }

        return y;  // Return the minimum
    }

/*--------------------------------------------------------------------
min - Find the minimum location in the matrix
--------------------------------------------------------------------*/

    public double nonZeroMin(int[] z) {
        z[0] = -1;
        z[1] = -1;
        double y = 0;  // Initial value

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                // IF THE ELEMENT IS LESS THAN THE CURRENT MIN
                if (((y == 0) || (x[i][j] < y)) && x[i][j] != 0) {
                    y = x[i][j];  // Store the new minimum
                    z[0] = i;
                    z[1] = j;
                }
            }
        }

        return y;  // Return the minimum
    }

/*--------------------------------------------------------------------
mean - Find the mean of the matrix
--------------------------------------------------------------------*/

    public double mean() {
        double m = 0;

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            m += rowMean(i);
        }

        return m / nr;
    }

/*--------------------------------------------------------------------
max - Find the maximum value in the matrix
--------------------------------------------------------------------*/

    public double max() {
        double y = x[0][0];  // Initial value

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                // IF THE ELEMENT IS GREATER THAN THE CURRENT MAX
                if (x[i][j] > y) {
                    y = x[i][j];  // Store the new maximum
                }
            }
        }

        return y;  // Return the maximum
    }

/*--------------------------------------------------------------------
max - Find the maximum location in the matrix
--------------------------------------------------------------------*/

    public double maxItem(int[] z) {
        z[0] = 0;
        z[1] = 0;
        double y = x[0][0];  // Initial value

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                // IF THE ELEMENT IS GREATER THAN THE CURRENT MAX
                if (x[i][j] > y) {
                    y = x[i][j];  // Store the new maximum
                    z[0] = i;
                    z[1] = j;
                }
            }
        }

        return y;  // Return the maximum
    }

/*--------------------------------------------------------------------
 sum - Sum the matrix
--------------------------------------------------------------------*/

    public double sum() {
        double z = 0;  // Sum value

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                // SUM THE ELEMENTS IN THE COLUMNS
                z += x[i][j];
            }
        }

        return z;  // Return the sum
    }

/*--------------------------------------------------------------------
rowSum - Sum a row of a matrix
--------------------------------------------------------------------*/

    public double rowSum(int i) {
        double z = 0;  // Sum value

        // LOOP THROUGH THE COLUMNS
        for (int j = 0; j < nc; j++) {
            // SUM THE ELEMENTS IN THE COLUMNS
            z += x[i][j];
        }

        return z;  // Return the sum
    }

/*--------------------------------------------------------------------
rowSum - Sum the rows of a matrix
--------------------------------------------------------------------*/

    public Matrix rowSum() {
        Matrix a = new Matrix(1, nr);  // Column sum matrix

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // INITIALIZE THE ELEMENT TO 0
            a.x[0][i] = 0;

            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                // SUM THE ELEMENTS IN THE COLUMNS
                a.x[0][i] = a.x[0][i] + x[i][j];
            }
        }

        return a;  // Return the matrix
    }

/*--------------------------------------------------------------------
rowMean - Average the row of a matrix
--------------------------------------------------------------------*/

    public double rowMean(int i) {
        return rowSum(i) / nc;
    }

/*--------------------------------------------------------------------
rowMean - Average the rows of a matrix
--------------------------------------------------------------------*/

    public Matrix rowMean() {
        Matrix a = rowSum();  // Sum the rows
        a.div(nc);          // Divide by the number of columns
        return a;              // Return the matrix
    }

/*--------------------------------------------------------------------
rowVar - Calculate the variance of a row of a matrix
--------------------------------------------------------------------*/

    public double rowVar(int i) {
        double z = 0;
        double m = rowMean(i);  // Find the mean of each row

        // LOOP THROUGH THE COLS
        for (int j = 0; j < nc; j++) {
            z += pow((x[i][j] - m), 2);
        }

        return z / (nc - 1);  // Return the matrix
    }

/*--------------------------------------------------------------------
rowVar - Calculate the variance of the rows of a matrix
--------------------------------------------------------------------*/

    public Matrix rowVar() {
        Matrix a = copy();         // Copy this matrix

        Matrix b = rowMean();      // Find the mean of each row
        a.sub(b.ntranspose());  // Subtract the mean from each element
        a.pow(2);                // Square each element

        Matrix c = a.rowSum();     // Sum the rows
        c.div(nc - 1);           // Divide by the number of columns

        return c;                   // Return the matrix
    }

/*--------------------------------------------------------------------
rowSD - Calculate the standard deviation of a row of a matrix
--------------------------------------------------------------------*/

    public double rowSD(int i) {
        return pow(rowVar(i), 0.5);
    }

/*--------------------------------------------------------------------
rowSD - Calculate the standard deviation of the rows of a matrix
--------------------------------------------------------------------*/

    public Matrix rowSD() {
        Matrix a = rowVar();  // Calculate the variance
        a.pow(0.5);         // Take the square root of the variance
        return a;              // Return the matrix
    }

/*--------------------------------------------------------------------
colSum - Sum the column of a matrix
--------------------------------------------------------------------*/

    public double colSum(int j) {
        double z = 0;  // Sum value

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // SUM THE ELEMENTS IN THE COLUMNS
            z += x[i][j];
        }

        return z;  // Return the sum
    }

/*--------------------------------------------------------------------
colSum - Sum the columns of a matrix
--------------------------------------------------------------------*/

    public Matrix colSum() {
        Matrix a = new Matrix(1, nc);  // Column sum matrix

        // LOOP THROUGH THE COLUMNS
        for (int j = 0; j < nc; j++) {
            // INITIALIZE THE ELEMENT TO 0
            a.x[0][j] = colSum(j);
        }

        return a;  // Return the matrix
    }

/*--------------------------------------------------------------------
colMean - Average the column of a matrix
--------------------------------------------------------------------*/

    public double colMean(int j) {
        return colSum(j) / nr;
    }

/*--------------------------------------------------------------------
colMean - Average the columns of a matrix
--------------------------------------------------------------------*/

    public Matrix colMean() {
        Matrix a = colSum();  // Sum the columns
        a.div(nr);          // Divide by the number of rows
        return a;              // Return the matrix
    }

/*--------------------------------------------------------------------
colVar - Calculate the variance of a column of a matrix
--------------------------------------------------------------------*/

    public double colVar(int j) {
        double z = 0;
        double m = colMean(j);  // Find the mean of each column

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            z += pow((x[i][j] - m), 2);
        }

        return z / (nr - 1);  // Return the matrix
    }

/*--------------------------------------------------------------------
colVar - Calculate the variance of the columns of a matrix
--------------------------------------------------------------------*/

    public Matrix colVar() {
        Matrix a = copy();      // Copy this matrix

        Matrix b = colMean();   // Find the mean of each column
        a.sub(b);             // Subtract the mean from each element
        a.pow(2);             // Square each element

        Matrix c = a.colSum();  // Sum the columns
        c.div(nr - 1);        // Divide by the number of rows

        return c;                // Return the matrix
    }

/*--------------------------------------------------------------------
colSD - Calculate the standard deviation of a column of a matrix
--------------------------------------------------------------------*/

    public double colSD(int j) {
        return pow(colVar(j), 0.5);
    }

/*--------------------------------------------------------------------
colSD - Calculate the standard deviation of the columns of a matrix
--------------------------------------------------------------------*/

    public Matrix colSD() {
        Matrix a = colVar();  // Calculate the variance
        a.pow(0.5);         // Take the square root of the variance
        return a;              // Return the matrix
    }

/*--------------------------------------------------------------------
colMin - Find the minimum value for a column
--------------------------------------------------------------------*/

    public double colMin(int j) {
        double z = x[0][j];  // Sum value

        // LOOP THROUGH THE ROWS
        for (int i = 1; i < nr; i++) {
            // IF THE ELEMENT IS LESS THAN THE CURRENT MIN
            if (x[i][j] < z) {
                z = x[i][j];  // Store the new minimum
            }
        }

        return z;  // Return the minimum
    }

/*--------------------------------------------------------------------
 getCol - Get a column
--------------------------------------------------------------------*/

    public double[] getCol(int j) {
        double[] z = new double[nr];

        for (int i = 0; i < z.length; i++) {
            z[i] = x[i][j];
        }

        return z;
    }
    
    
    /*--------------------------------------------------------------------
    replaceRow - Replace a row
   --------------------------------------------------------------------*/

       public void replaceRow(int i, double[] newRow) {
    	   x[i]=newRow;
       }

/*--------------------------------------------------------------------
 getRow - Get a row
--------------------------------------------------------------------*/

    public double[] getRow(int i) {
        double[] z = new double[nc];

        for (int j = 0; j < z.length; j++) {
            z[j] = x[i][j];
        }

        return z;
    }

/*--------------------------------------------------------------------
colMin - Find the minimum value for each column
--------------------------------------------------------------------*/

    public Matrix colMin() {
        Matrix a = new Matrix(1, nc);  // Column min matrix

        // LOOP THROUGH THE COLUMNS
        for (int j = 0; j < nc; j++) {
            // INITIALIZE THE ELEMENT TO 0
            a.x[0][j] = colMin(j);
        }

        return a;  // Return the matrix
    }

/*--------------------------------------------------------------------
colMax - Find the maximum value for a column
--------------------------------------------------------------------*/

    public double colMax(int j) {
        double z = x[0][j];  // Sum value

        // LOOP THROUGH THE ROWS
        for (int i = 1; i < nr; i++) {
            // IF THE ELEMENT IS LESS THAN THE CURRENT MAX
            if (x[i][j] > z) {
                z = x[i][j];  // Store the new maximum
            }
        }

        return z;  // Return the maximum
    }

/*--------------------------------------------------------------------
colMaxInd - Find the maximum value for a column
--------------------------------------------------------------------*/

    public int colMaxInd(int j) {
        int z = 0;

        // LOOP THROUGH THE ROWS
        for (int i = 1; i < nr; i++) {
            // IF THE ELEMENT IS GREATER THAN THE CURRENT MAX
            if (x[i][j] > x[z][j]) {
                z = i;  // Store the new maximum
            }
        }

        return z;  // Return the maximum
    }

/*--------------------------------------------------------------------
 brier - Calculate the brier score
--------------------------------------------------------------------*/

    public double[] brier() {
        double[] z = new double[nc];

        Matrix m = copy();
        m.div(m.colSum());
        int[] ind = d2i(m.colMaxInd().x[0]);

        double ng = (double) (nr);

        double avg = pow((1 - (1 / ng)), 2) + (ng - 1) * pow((1 / ng), 2);

        for (int i = 0; i < nc; i++) {
            z[i] = pow((1 - m.x[ind[i]][i]), 2);

            for (int j = 0; j < nr; j++) {
                if (j != ind[i]) {
                    z[i] += pow(m.x[j][i], 2);
                }
            }

            z[i] = 1 - (z[i] / avg);
        }

        return z;
    }

    public double[] brierNoNorm() {
        double[] z = new double[nc];

        Matrix m = copy();
        m.div(m.colSum());
        int[] ind = d2i(m.colMaxInd().x[0]);

        double ng = (double) (nr);

        for (int i = 0; i < nc; i++) {
            z[i] = pow((1 - m.x[ind[i]][i]), 2);

            for (int j = 0; j < nr; j++) {
                if (j != ind[i]) {
                    z[i] += pow(m.x[j][i], 2);
                }
            }
        }

        return z;
    }

    public double[] brierNoNorm(int[] ind) {
        double[] z = new double[nc];

        Matrix m = copy();
        m.div(m.colSum());
        int[] ind2 = d2i(m.colMaxInd().x[0]);

        double ng = (double) (nr);

        for (int i = 0; i < nc; i++) {
            z[i] = pow((1 - m.x[ind[i]][i]), 2);

            for (int j = 0; j < nr; j++) {
                if (j != ind[i]) {
                    z[i] += pow(m.x[j][i], 2);
                }
            }
        }

        return z;
    }

/*--------------------------------------------------------------------
euclidean - Calculate the Euclidean distance betwen columns
--------------------------------------------------------------------*/

    public double euclidean(int a, int b) {
        double z = 0;  // Distance holder

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            z += (pow(x[i][a] - x[i][b], 2));  // ADD EACH SQUARED DIFFERENCE
        }

        return pow(z, 0.5);  // Return the distance
    }

/*--------------------------------------------------------------------
euclidean - Calculate the Euclidean distances from zero
--------------------------------------------------------------------*/

    public double[] euclidean() {
        double[] d = new double[nc];

        // LOOP THROUGH THE COLUMNS
        for (int j = 0; j < nc; j++) {
            // LOOP THROUGH THE ROWS
            for (int i = 0; i < nr; i++) {
                d[j] += (pow(x[i][j], 2));  // ADD EACH SQUARED DIFFERENCE
            }

            d[j] = pow(d[j], 0.5);
        }

        return d;
    }

/*--------------------------------------------------------------------
maximum - Calculate the Maximum distance betwen columns
--------------------------------------------------------------------*/

    public double maximum(int a, int b) {
        double z = 0;    // Distance holder
        double max = 0;  // Maximum

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            z = (pow(x[i][a] - x[i][b], 2));  // SQUARED DISTANCE

            // IF LARGER THAN THE CURRENT MAXIMUM
            if (z > max) {
                max = z;  // Set z equal to max
            }
        }

        return pow(max, 0.5);  // Return the distance
    }

/*--------------------------------------------------------------------
 manhattan - Calculate the Manhattan distance betwen columns
--------------------------------------------------------------------*/

    public double manhattan(int a, int b) {
        double z = 0;     // Distance holder

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            z += abs(x[i][a] - x[i][b]);  // Add the term
        }

        return z;  // Return the distance
    }

/*--------------------------------------------------------------------
 canberra - Calculate the Canberra distance betwen columns
--------------------------------------------------------------------*/

    public double canberra(int a, int b) {
        double z = 0;    // Distance holder
        double num = 0;  // Numerator
        double den = 0;  // Denominator

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            num = abs(x[i][a] - x[i][b]);  // Calculate the denominator
            den = abs(x[i][a] + x[i][b]);  // Calculate the numerator

            // IF THE NUMERATOR AND DENOMINATOR ARE > 0
            if (num > 0 && den > 0) {
                z += num / den;  // Add to the sum
            }
        }

        return z;  // Return the distance
    }

/*--------------------------------------------------------------------
 minkowski - Calculate the Minkowski distance betwen columns
--------------------------------------------------------------------*/

    public double minkowski(int a, int b, double p) {
        p = max(p, 1);  // Ensure it is >= 1
        double z = 0;     // Distance holder

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            z += pow(abs(x[i][a] - x[i][b]), p);  // Add the term
        }

        return pow(z, (1 / p));  // Return the distance
    }

/*--------------------------------------------------------------------
zeroCol - Zero out a column
--------------------------------------------------------------------*/

    public void zeroCol(int j) {
        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // INITIALIZE THE ELEMENT TO 0
            x[i][j] = 0.0;
        }
    }
    
/*--------------------------------------------------------------------
 *  Probabilistic objective distance (Lee-Seung MNF objective function)   
 -------------------------------------------------------------------*/
    public double probabilistic(Matrix y) {
        Matrix Mn = nadd(e).ndiv(y.nadd(e)).nlog();
        Mn = nmul(Mn);
        Mn.sub(this);
        Mn.add(y);
        Mn.div(nr * nc);

        return Mn.sum();
    }
    
/*--------------------------------------------------------------------
 * Euclidean distance between this and another matrix
 -------------------------------------------------------------------*/
    public double euclidean(Matrix y) {
    	double dist = 0;
    	// LOOP THROUGH EACH ROW IN Z
    	for (int i = 0; i < nr; i = i + 1) {
    		// LOOP THROUGH EACH COLUMN IN Z
    		for (int j = 0; j < nc; j = j + 1) {
    			dist = dist + abs(pow((y.x[i][j]) - (x[i][j]), 2));
    		}
    	}
    	return dist;
    }
/*--------------------------------------------------------------------
zeroRow - Zero out a row
--------------------------------------------------------------------*/

    public void zeroRow(int i) {
        // LOOP THROUGH THE ROWS
        for (int j = 0; j < nc; j++) {
            // INITIALIZE THE ELEMENT TO 0
            x[i][j] = 0.0;
        }
    }

/*--------------------------------------------------------------------
col - Return a column
--------------------------------------------------------------------*/

    public double[] col(int j) {
        double[] z = new double[nr];

        // LOOP THROUGH THE COLUMNS
        for (int i = 0; i < nr; i++) {
            // SAVE THE ITEM
            z[i] = x[i][j];
        }

        return z;  // Return the row
    }

/*--------------------------------------------------------------------
colMaxInd - Find the index of the maximum value for each column
--------------------------------------------------------------------*/

    public Matrix colMaxInd() {
        Matrix a = new Matrix(1, nc);  // Column max matrix

        // LOOP THROUGH THE COLUMNS
        for (int j = 0; j < nc; j++) {
            // INITIALIZE THE ELEMENT TO 0
            a.x[0][j] = colMaxInd(j);
        }

        return a;  // Return the matrix
    }

/*--------------------------------------------------------------------
colMax - Find the maximum value for each column
--------------------------------------------------------------------*/

    public Matrix colMax() {
        Matrix a = new Matrix(1, nc);  // Column max matrix

        // LOOP THROUGH THE COLUMNS
        for (int j = 0; j < nc; j++) {
            // INITIALIZE THE ELEMENT TO 0
            a.x[0][j] = colMax(j);
        }

        return a;  // Return the matrix
    }

/*--------------------------------------------------------------------
colRange - Find the range for a column
--------------------------------------------------------------------*/

    public double colRange(int j) {
        return colMax(j) - colMin(j);
    }

/*--------------------------------------------------------------------
colRange - Find the range for each column
--------------------------------------------------------------------*/

    public Matrix colRange() {
        Matrix a = colMax();  // Find the column maximums
        a.sub(colMin());     // Subtract out the column minimums
        return a;              // Return the matrix
    }

/*--------------------------------------------------------------------
rowMin - Find the minimum value for a row
--------------------------------------------------------------------*/

    public double rowMin(int i) {
        double z = x[i][0];  // Initial value

        // LOOP THROUGH THE COLUMNS
        for (int j = 1; j < nc; j++) {
            // IF THE ELEMENT IS LESS THAN THE CURRENT MIN
            if (x[i][j] < z) {
                z = x[i][j];  // Store the new minimum
            }
        }

        return z;  // Return the maximum
    }

/*--------------------------------------------------------------------
rowMin - Find the minimum value for each row
--------------------------------------------------------------------*/

    public Matrix rowMin() {
        Matrix a = new Matrix(nr, 1);  // Row min matrix

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // INITIALIZE THE ELEMENT TO 0
            a.x[i][0] = rowMin(i);
        }

        return a;  // Return the matrix
    }

/*--------------------------------------------------------------------
rowMax - Find the maximum value for a row
--------------------------------------------------------------------*/

    public double rowMax(int i) {
        double z = x[i][0];  // Initial value

        // LOOP THROUGH THE COLUMNS
        for (int j = 1; j < nc; j++) {
            // IF THE ELEMENT IS GREATER THAN THE CURRENT MAX
            if (x[i][j] > z) {
                z = x[i][j];  // Store the new maximum
            }
        }

        return z;  // Return the maximum
    }

/*--------------------------------------------------------------------
row - Return a row
--------------------------------------------------------------------*/

    public double[] row(int i) {
        double[] z = new double[nc];

        // LOOP THROUGH THE ROWS
        for (int j = 0; j < nc; j++) {
            // SAVE THE ITEM
            z[j] = x[i][j];
        }

        return z;  // Return the row
    }

/*--------------------------------------------------------------------
rowMax - Find the maximum value for each row
--------------------------------------------------------------------*/

    public Matrix rowMax() {
        Matrix a = new Matrix(nr, 1);  // Row max matrix

        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // INITIALIZE THE ELEMENT TO 0
            a.x[i][0] = rowMax(i);
        }

        return a;  // Return the matrix
    }

/*--------------------------------------------------------------------
rowRange - Find the range for a row
--------------------------------------------------------------------*/

    public double rowRange(int i) {
        return rowMax(i) - rowMin(i);
    }

/*--------------------------------------------------------------------
rowRange - Find the range for each row
--------------------------------------------------------------------*/

    public Matrix rowRange() {
        Matrix a = rowMax();  // Find the row maximums
        a.sub(rowMin());   // Subtract out the row minimums
        return a;              // Return the matrix
    }

/*--------------------------------------------------------------------
random - Randomize the matrix
--------------------------------------------------------------------*/

    public void random() {
    	// LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                // SET THE ELEMENT TO A RANDOM VALUE
                //x[i][j] = randDouble();
               // x[i][j]=rand.nextDouble();
               //if(j==0){ System.err.println(x[i][j]);}
                x[i][j] = Math.random();
            }
        }
    }

/*--------------------------------------------------------------------
addRandom - Add random perturbations
--------------------------------------------------------------------*/

    public void addRandom() {
        addRandom(0.0000001);
    }

/*--------------------------------------------------------------------
addRandom - Add random perturbations
--------------------------------------------------------------------*/

    public void addRandom(double d) {
        // LOOP THROUGH THE ROWS
        for (int i = 0; i < nr; i++) {
            // LOOP THROUGH THE COLUMNS
            for (int j = 0; j < nc; j++) {
                // SET THE ELEMENT TO A RANDOM VALUE
                x[i][j] += randDouble() * d;
            }
        }
    }

    /*--------------------------------------------------------------------
    Remove column - remove a column from matrix
    --------------------------------------------------------------------*/

        public void removeColumn(int columnIndex) {
        	double[][] newMatrix=new double[nr][nc-1];
	        // LOOP THROUGH THE ROWS
	        for (int i = 0; i < nr; i++) {
	        	int counter=0;
	        	for(int j=0; j<nc; j++){
	        		if(j!=columnIndex){newMatrix[i][counter++]=x[i][j];}
	            }
	        }           
	        this.x=newMatrix;
	        this.nc=nc-1;
      }
        
        public void removeColumns(Set<Integer> columns) {
        	double[][] newMatrix=new double[nr][nc-columns.size()];
	        // LOOP THROUGH THE ROWS
	        for (int i = 0; i < nr; i++) {
	        	int counter=0;
	        	for(int j=0; j<nc; j++){
	        		if(!columns.contains(j)){newMatrix[i][counter++]=x[i][j];}
	            }
	        }           
	        this.x=newMatrix;
	        this.nc=nc-columns.size();
      }
	
     /*--------------------------------------------------------------------
	 removeRow - Remove a row to the current matrix
	 --------------------------------------------------------------------*/
	
	 public void removeRow(int rowIndex) {
		 double[][] newMatrix=new double[nr-1][nc];
	     
		 
		 
		 // LOOP THROUGH THE ROWS
	     int counter=0;
	     for (int i = 0; i < nr; i++) {
	    	 if(i!=rowIndex){newMatrix[counter++]=x[i];}
	     }
	     this.x=newMatrix;
	     this.nr=nr-1;
	 }
	 
	 public void removeRows(Set<Integer> rows) {
		 double[][] newMatrix=new double[nr-rows.size()][nc];
	     
		  // LOOP THROUGH THE ROWS
	     int counter=0;
	     for (int i = 0; i < nr; i++) {
	    	 if(!rows.contains(i)){newMatrix[counter++]=x[i];}
	     }
	     this.x=newMatrix;
	     this.nr=nr-rows.size();
	 }
    
    /*--------------------------------------------------------------------
    addRow - Append a row to the current matrix
    --------------------------------------------------------------------*/

        public void addRow(double[] newRow) {
            double[][] newMatrix=new double[nr+1][nc];
        	// LOOP THROUGH THE ROWS
            for (int i = 0; i < nr; i++) {
                newMatrix[i]=x[i];
            }
            newMatrix[nr]=newRow;
            this.x=newMatrix;
            this.nr=nr+1;
      }
       
        public void addRow(ArrayList<Double> newRow) {
            double[] row=new double[newRow.size()];
            for(int i=0; i<newRow.size(); i++){
            	row[i]=newRow.get(i);
            }
            addRow(row);
      }
        
        public void addRows(Matrix matrix) {
           for(int i=0; i<matrix.nr; i++){
        	  addRow(matrix.getRow(i));
           }
      }
        
        public void addColumns(Matrix matrix){
        	for(int i=0; i<matrix.nc; i++){addColumn(matrix.getCol(i));}
        }
        
        /*--------------------------------------------------------------------
	    addColumn - Append a column to the current matrix
	    --------------------------------------------------------------------*/
	
	    public void addColumn(double[] newColumn) {
	    	double[][] newMatrix=new double[nr][nc+1];
	        // LOOP THROUGH THE ROWS
	        for (int i = 0; i < nr; i++) {
	        	for(int j=0; j<nc; j++){
	        		newMatrix[i][j]=x[i][j];
	            }
	          newMatrix[i][nc]=newColumn[i];
	        }           
	        this.x=newMatrix;
	        this.nc=nc+1;
	     }
    
/*--------------------------------------------------------------------
NMFdiv - Run the NMF algorithm on the matrix
--------------------------------------------------------------------*/
    public void NMFdiv(int k, NMFCostFunction cost) {
       
        NMFdiv(k, cost, 0.0);
    }
    
    public void NMFdiv(int k, NMFCostFunction cost, double precision) {
        H = new Matrix(k, nc);          // New matrix for H
        H.random();                       // Randomize the matrix
        NMFdiv(H, cost == null ? DEFAULT_COST_FUNCTION : cost, precision);
    }
    
    
    
    public void NMFdiv(int k) {
    	NMFdiv(k, null);

    }
    
    /*--------------------------------------------------------------------
    NMFdiv - Run the NMF algorithm on the matrix with specific initial H
    --------------------------------------------------------------------*/
    	public void NMFdiv(Matrix Hi) {
    		NMFdiv(Hi, DEFAULT_COST_FUNCTION);
    	}
    	
    	public void NMFdiv(Matrix Hi, NMFCostFunction cost) {
    		NMFdiv(Hi, cost, 0.0);
    	}
    	
    	
    	 public void NMFdiv(int k, NMFCostFunction cost, double precision, Matrix WFixed) {
    		 H = new Matrix(k, nc);          // New matrix for H
    	     H.random();                       // Randomize the matrix
    		 int num = 2000;  // Number of iterations
            
             //W = new Matrix(nr, k);          // New matrix for W
             //W.random();                       // Randomize the matrix
             W=WFixed.copy();
             
             Matrix I = new Matrix(nr, k);   // New matrix for I
             //H = new Matrix(k, nc);          // New matrix for H
             //H.random();                       // Randomize the matrix
             
             Matrix J = new Matrix(k, nc);   // New matrix for J
             Matrix K = new Matrix(nr, nc);  // New matrix for K
             Matrix V = new Matrix(nr, nc);  // New matrix for V
             Matrix Y = new Matrix(k, 1);    // New matrix for Y
             Matrix t = new Matrix(k, nr);   // New matrix for t
             Matrix u = new Matrix(nc, k);   // New matrix for u

             List<Double> errors = new ArrayList<Double>();

             //out("\nNMF Algorithm");

             erri = num - 1;

             // ITERATE THE ALGORITHM num TIMES
             int a = 0;
             
             double deltaError=Double.MAX_VALUE;
             double error=Double.MAX_VALUE;
             while( a < num && deltaError>precision/*&& (a < 2 || abs(errors.get(a - 1) - errors.get(a-2)) > precision)*/) {
                  V = W.nmatmul(H);                    // Multiply W and H
                 t = W.ntranspose();                    // Transpose W

                 // CREATE A NEW ESTIMATE FOR H
                 K = ndiv(V);       // Divide this matrix by V
                 J = t.nmatmul(K);  // Multiply t by K
                 J = H.nmul(J);     // Multiply H by J
                 H = J.nadd(e);     // Add e to J

                 // CALCULATE THE SUM FOR EACH COLUMN
                 Y = W.colSum();  // Sum the columns of W

                 // LOOP THROUGH EACH ROW IN H
                 for (int i = 0; i < k; i = i + 1) {
                     // LOOP THROUGH EACH COLUMN IN H
                     for (int j = 0; j < nc; j = j + 1) {
                         // NORMALIZE THE ROW
                         H.x[i][j] = (H.x[i][j]) / (Y.x[0][i]);
                     }
                 }

                 V = W.nmatmul(H);  // Multiply W and H
                 u = H.ntranspose();  // Transpose u

                 // CREATE A NEW ESTIMATE FOR W
                 K = ndiv(V);       // Divide this matrix by V
                 I = K.nmatmul(u);  // Multiply K by u
                 I = W.nmul(I);     // Multuple W by I
                 //W = I.nadd(e);     // Add e to I //NOT UODATING W

                 // CALCULATE THE SUM FOR EACH ROW
                 Y = H.rowSum();  // Sum the rows of H

                 // LOOP THROUGH EACH COLUMN IN W
                 for (int j = 0; j < k; j = j + 1) {
                     // LOOP THROUGH EACH ROW IN W
                     for (int i = 0; i < nr; i = i + 1) {
                         // NORMALIZE THE ROW
                         //W.x[i][j] = (W.x[i][j]) / (Y.x[0][j]);//NOT UPDATING W
                     }
                 }

                 // CALCULATE THE ERROR OF THE CURRENT ESTIMATE.
                 errors.add(cost.cost(this, V));
                 double newError=cost.cost(this, V);
                 deltaError=error-newError;
                 error=newError;
                 System.err.println("Iteration "+a+" Error "+error+" Delat Error: "+deltaError);
                 //errors.add(probabilistic(V));
                 a++;  
                 //System.err.println("\tNMF iteration " +a + " cost " + errors.get(errors.size() - 1));
                 //Tune "smallNumber" to the magnitude of the error.
                 //if(smallNumber == SMALL_NUMBER) {
                 	//smallNumber = errors.get(errors.size() - 1)/(double)1000;
                 	//System.err.println("Adjusted smallNumber " + smallNumber);
                 //}
              }
              err = new double[errors.size()];
              for(int i = 0;i < errors.size(); i++){err[i] = errors.get(i);}
              
              //System.err.println("Finished NMF iterations " + a + " cost " + getNMFError() + " second but last cost " + err[err.length - 2] + " ecuclidean cost " + euclidean(V) + " probabilistic cost " + probabilistic(V));
         }
    	
       /* public void NMFdiv(Matrix Hi, NMFCostFunction cost, double precision) {
        	int k=Hi.nr;
            int num = 2000;  // Number of iterations
           
            W = new Matrix(nr, k);          // New matrix for W
            W.random();                       // Randomize the matrix

            Matrix I = new Matrix(nr, k);   // New matrix for I
            //H = new Matrix(k, nc);          // New matrix for H
            //H.random();                       // Randomize the matrix
            H=Hi; //specify to initial H
            
            Matrix J = new Matrix(k, nc);   // New matrix for J
            Matrix K = new Matrix(nr, nc);  // New matrix for K
            Matrix V = new Matrix(nr, nc);  // New matrix for V
            Matrix Y = new Matrix(k, 1);    // New matrix for Y
            Matrix t = new Matrix(k, nr);   // New matrix for t
            Matrix u = new Matrix(nc, k);   // New matrix for u

            List<Double> errors = new ArrayList<Double>();

            //out("\nNMF Algorithm");

            erri = num - 1;

            // ITERATE THE ALGORITHM num TIMES
            int a = 0;
            
            double deltaError=Double.MAX_VALUE;
            double error=Double.MAX_VALUE;
            while( a < num && deltaError>precision) {
                 V = W.nmatmul(H);                    // Multiply W and H
                t = W.ntranspose();                    // Transpose W

                // CREATE A NEW ESTIMATE FOR H
                K = ndiv(V);       // Divide this matrix by V
                J = t.nmatmul(K);  // Multiply t by K
                J = H.nmul(J);     // Multiply H by J
              //  H = J.nadd(e);     // Add e to J
                H=J;
                
                // CALCULATE THE SUM FOR EACH COLUMN
                Y = W.colSum();  // Sum the columns of W

                // LOOP THROUGH EACH ROW IN H
                for (int i = 0; i < k; i = i + 1) {
                    // LOOP THROUGH EACH COLUMN IN H
                    for (int j = 0; j < nc; j = j + 1) {
                        // NORMALIZE THE ROW
                        H.x[i][j] = (H.x[i][j]) / (Y.x[0][i]);
                    }
                }

                V = W.nmatmul(H);  // Multiply W and H
                u = H.ntranspose();  // Transpose u

                // CREATE A NEW ESTIMATE FOR W
                K = ndiv(V);       // Divide this matrix by V
                I = K.nmatmul(u);  // Multiply K by u
                I = W.nmul(I);     // Multuple W by I
               // W = I.nadd(e);     // Add e to I
                W=I;
                
                // CALCULATE THE SUM FOR EACH ROW
                Y = H.rowSum();  // Sum the rows of H

                // LOOP THROUGH EACH COLUMN IN W
                for (int j = 0; j < k; j = j + 1) {
                    // LOOP THROUGH EACH ROW IN W
                    for (int i = 0; i < nr; i = i + 1) {
                        // NORMALIZE THE ROW
                        W.x[i][j] = (W.x[i][j]) / (Y.x[0][j]);
                    }
                }

                // CALCULATE THE ERROR OF THE CURRENT ESTIMATE.
                errors.add(cost.cost(this, V));
                double newError=cost.cost(this, V);
                deltaError=error-newError;
                error=newError;
                System.err.println("Iteration "+a+" Error "+error+" Delat Error: "+deltaError);
                //errors.add(probabilistic(V));
                a++;  
                //System.err.println("\tNMF iteration " +a + " cost " + errors.get(errors.size() - 1));
                //Tune "smallNumber" to the magnitude of the error.
                //if(smallNumber == SMALL_NUMBER) {
                	//smallNumber = errors.get(errors.size() - 1)/(double)1000;
                	//System.err.println("Adjusted smallNumber " + smallNumber);
                //}
             }
             err = new double[errors.size()];
             for(int i = 0;i < errors.size(); i++){err[i] = errors.get(i);}
             
             System.err.println("Finished NMF iterations " + a + " cost " + getNMFError() + " second but last cost " + err[err.length - 2] + " ecuclidean cost " + euclidean(V) + " probabilistic cost " + probabilistic(V));
        }*/
        
        
       public void NMFdiv(Matrix Hi, NMFCostFunction cost, double precision) {
        	//System.err.println(cost.getName());
    	   int k=Hi.nr;
        	W = new Matrix(nr, k);          // New matrix for W
            W.random();                       // Randomize the matrix
            H=Hi; //specify to initial H
            int num = 2000;  // Number of iterations
            List<Double> errors = new ArrayList<Double>();

            
            int a = 0;
            double deltaError=Double.MAX_VALUE;
            double error=Double.MAX_VALUE;
            while( a < num && deltaError>precision) {
            	Matrix[] HW=cost.update(this, W, H);
            	W=HW[0];
            	H=HW[1];
            	Matrix VNew=W.nmatmul(H);
            	
               	// CALCULATE THE ERROR OF THE CURRENT ESTIMATE.
                errors.add(cost.cost(this, VNew));
                double newError=cost.cost(this, VNew);
                deltaError=error-newError;
                error=newError;
                System.err.println("Iteration "+a+" Error "+error+" Delta Error: "+deltaError);
            	a++;
            }
            this.errors=errors;
        }
                
        public String toString(){
      	   String rtrn="";
      	   for(int i=0; i<this.nr; i++){
      		   for(int j=0; j<this.nc; j++){
      			   rtrn+=this.x[i][j]+" ";
      		   }
      		   rtrn+="\n";
      	   }
      	   return rtrn;
         }

		public double getNMFError() {
			return this.errors.get(errors.size()-1);
			//return err != null ? err[err.length - 1] : 0;
		}

		public void setPrecision(double d) {
			this.precision = d;
			
		}
    
}
