/* hb_filters.h
................................................................................
                          *******  HydroBase2  *******
               ...................................................
                          
                    author:  Ruth  Curry
                             Woods Hole Oceanographic Institution
                             2000
................................................................................
*/

#ifndef HB_FILTERS
#define HB_FILTERS 1


#ifndef TOOBIG
#define TOOBIG 0.9e+35     /* value larger than any oceanographic property */
#endif

#ifndef ABS
#define    ABS(x)       (((x) < 0) ? -(x) : (x))
#endif

/*  prototypes for functions in filters.c */

void filter(char *, double *, int, int);  
   /* interface to 1-d filters */
   
double filter2d(double *, double *, double , double, int, int, int);
    /*  Smooths value in middle of array xin (at index xrad, yrad)
      which represents a 2D array with dimensions 
                ncols = xrad * 2 + 1 
		nrows = yrad * 2 + 1 
       using weighted arithmetic average of surrounding gridnodes.
       Array element 0 corresponds to row0, col0;
       elements are stored in row order (column varying fastest). */
  
void gauss(double *, int, int);   
   /*  Uses a weighted arithmetic average to smooth an array of values.
       The filter width is the total number of pts which contribute to
       any one smoothed value and should be an odd number to incorporate an
       equal number of points above and below each point */


double *laplacian(double *, int, int, double, double, double, int);

  /* performs a 5-pt laplacian filter on the input array (dimension nx * ny
     stored in row-order:  sq = row * ncols + col). Filter is:
     
         out[ij] = in[ij] + weight * sumof(in[neighbor] - in[ij]) / 4  
     
           in :  input array
        nx, ny:  x, y dimensions
        weight:  weight factor 
     empty_val:  marks empty nodes
      mask_val:   marks masked nodes
        report:  BOOLEAN flag to print info to stderr
     
     returns a pointer to the filtered array
     
  */
#endif /* ifndef HB_FILTERS */
