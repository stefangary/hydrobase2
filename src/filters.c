/*  filters.c
................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
                             1993
................................................................................
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "hb_filters.h"
#include "hb_memory.h"

void filter(char *type, double *x, int nx, int width)
  /* type:    filter identifier
     x:       array to filter 
     nx:      size of the array 
     width:   number of pts in filter 
  */
{
   void gauss(double *, int, int);
   
   if (width <= 1)
      return;
    
   switch (type[0]) {
      case 'g':
      case 'G':
            gauss(x,nx,width);
            break;
      default:
            fprintf(stderr,"\nUnknown filter type <%c> requested in call to filter()\n", type);
            break;
   } /*end switch */
   return;
}  /* end filter() */
/*****************************************************************************/
void gauss(double *x, int nx, int width)

   /*  Uses a weighted arithmetic average to smooth an array of values.
       The filter width is the total number of pts which contribute to
       any one smoothed value and should be an odd number to incorporate an
       equal number of points above and below each point */
{
   double *tmp, *wght;
   double pi = 3.1415926;
   double alpha, sum;
   int i, j, k, halfwidth, window, end;
      
   tmp = (double *) calloc((size_t)nx, sizeof(double));
   
   halfwidth = width >> 1;      /* divide by 2 */
   width = (halfwidth << 1) + 1;   /* ensure that width is an odd number */
   wght = (double *) calloc((size_t) width, sizeof(double));
   
   alpha = (pi / (double) halfwidth);
   alpha = alpha * alpha / 4.5;
   sum = 0.0;
   
   /* first , calculate weights for bottom half of filter window... */
   
   i = 0;
   while (i < halfwidth) {
       wght[i] = halfwidth - i;
       wght[i] = exp(-alpha * wght[i] * wght[i]);
       sum += wght[i];
       ++i;
   }
   wght[halfwidth] = 1.0;    /* weight at central element */
   sum = sum + sum + 1.0;    /* sum of weights for entire filter */
   
   /* normalize the weights */
   
   for (i = 0; i <= halfwidth; ++i) 
      wght[i] /= sum;
   
   /* now set the weights in the top half of the filter to mirror the bottom */
   
   j = halfwidth - 1;
   for (i = halfwidth+1; i < width; ++i ) {
      wght[i] = wght[j--];
   }
   
   /* copy the data array into temporary space */
   
   for (i = 0; i < nx; ++i) {
      tmp[i] = x[i];
   }
   
   /* apply the filter -- this leaves the values at the top and bottom
      of the array (distance of halfwidth) unchanged */
   
   end = nx - halfwidth;
   for (i = halfwidth; i < end; ++i) {
      x[i] = 0.0;
      window = i + halfwidth;
      k = 0;
      for (j = i-halfwidth; j <= window; ++j) 
         x[i] += tmp[j] * wght[k++];
   }
   free((void *)tmp);
   free((void *)wght);
   return;
   
}  /* end gauss() */
/******************************************************************/
double filter2d(double *xin, double *wghts, double empty_val, double mask_val, int xrad, int yrad, int ncols)
   /*  Smooths value in middle of array xin (at index xrad, yrad)
      which represents a 2D array with dimensions 
                ncols = xrad * 2 + 1 
		nrows = yrad * 2 + 1 
       using weighted values in surrounding gridnodes. Array element 0 corresponds to row0, col0;
       elements are stored in row order (column varying fastest). 

      Corresponding weights stored in similar grid (wghts) are applied for smoothing.
      Searches in each of 8 directions from central gridnode for gridnodes to incorporate. 
      Missing values are denoted by empty_val (large negative), masked values (below seafloor) 
      are large positive values (mask_val).
      If mask_val is encountered, the search is suspended in that direction to prevent
      mixing watermasses across topographic barriers.  
      Returned values:
          xoutptr:   weighted mean
   */

{
   int row, col, stop, sq, nobs, xradius, yradius;
   double  wsum, xwsum, toosmall, toobig;  

   xwsum = 0.0;
   wsum = 0.0;
   nobs = 0;
   
   toobig = mask_val / 10.0;    /* for testing */
   toosmall = empty_val / 10.0;
   
/* Start with middle gridnode */  
   row = yrad;
   col = xrad;
   sq = row * ncols + col;
   if (xin[sq] > toosmall)  {  /* either a real value or a mask */
         if (xin[sq] < toobig) {
            xwsum += xin[sq] * wghts[sq];
	   wsum += wghts[sq];
	   ++nobs;
         }
    }  
     
/*search west */
   xradius = 0;
   stop = 0;
   while (++xradius <= xrad && !stop) {
      --col;
      sq = row * ncols + col;
      if (xin[sq] > toosmall)  {  /* either a real value or a mask */
         if (xin[sq] >= toobig) 
	   stop = 1;
	 else {
	   xwsum += xin[sq] * wghts[sq];
	   wsum += wghts[sq];
	   ++nobs;
	 }
      }
   } /* end while */
   
/*search east */
   row = yrad;
   col = xrad;
   xradius = 0;
   stop = 0;
   while (++xradius <= xrad && !stop) {
      ++col;
      sq = row * ncols + col;
      if (xin[sq] > toosmall)  {  /* either a real value or a mask */
         if (xin[sq] >= toobig) 
	   stop = 1;
	 else {
	   xwsum += xin[sq] * wghts[sq];
	   wsum += wghts[sq];
	   ++nobs;
	 }
      }
   } /* end while */
  
   
/*search north */
   row = yrad;
   col = xrad;
   yradius = 0;
   stop = 0;
   while (++yradius <= yrad && !stop) {
      ++row;
      sq = row * ncols + col;
      if (xin[sq] > toosmall)  {  /* either a real value or a mask */
         if (xin[sq] >= toobig) 
	   stop = 1;
	 else {
	   xwsum += xin[sq] * wghts[sq];
	   wsum += wghts[sq];
	   ++nobs;
	 }
      }
   } /* end while */
   
/*search south */
   row = yrad;
   col = xrad;
   yradius = 0;
   stop = 0;
   while (++yradius <= yrad && !stop) {
      --row;
      sq = row * ncols + col;
      if (xin[sq] > toosmall)  {  /* either a real value or a mask */
         if (xin[sq] >= toobig) 
	   stop = 1;
	 else {
	   xwsum += xin[sq] * wghts[sq];
	   wsum += wghts[sq];
	   ++nobs;
	 }
      }
   } /* end while */
   
/*search northwest */
   row = yrad;
   col = xrad;
   yradius = 0;
   xradius = 0;
   stop = 0;
   while ((++xradius < xrad) && (++yradius < yrad) && !stop) {
      ++row;
      --col;
      sq = row * ncols + col;
      if (xin[sq] > toosmall)  {  /* either a real value or a mask */
         if (xin[sq] >= toobig) 
	   stop = 1;
	 else {
	   xwsum += xin[sq] * wghts[sq];
	   wsum += wghts[sq];
	   ++nobs;
	 }
      }
   } /* end while */
   
/*search northeast */
   row = yrad;
   col = xrad;
   yradius = 0;
   xradius = 0;
   stop = 0;
   while ((++xradius < xrad) && (++yradius < yrad) && !stop) {
      ++row;
      ++col;
      sq = row * ncols + col;
      if (xin[sq] > toosmall)  {  /* either a real value or a mask */
         if (xin[sq] >= toobig) 
	   stop = 1;
	 else {
	   xwsum += xin[sq] * wghts[sq];
	   wsum += wghts[sq];
	   ++nobs;
	 }
      }
   } /* end while */

/*search southeast */
   row = yrad;
   col = xrad;
   yradius = 0;
   xradius = 0;
   stop = 0;
   while ((++xradius < xrad) && (++yradius < yrad) && !stop) {
      --row;
      ++col;
      sq = row * ncols + col;
      if (xin[sq] > toosmall)  {  /* either a real value or a mask */
         if (xin[sq] >= toobig) 
	   stop = 1;
	 else {
	   xwsum += xin[sq] * wghts[sq];
	   wsum += wghts[sq];
	   ++nobs;
	 }
      }
   } /* end while */
   
/*search southwest */
   row = yrad;
   col = xrad;
   yradius = 0;
   xradius = 0;
   stop = 0;
   while ((++xradius < xrad) && (++yradius < yrad) && !stop) {
      --row;
      --col;
      sq = row * ncols + col;
      if (xin[sq] > toosmall)  {  /* either a real value or a mask */
         if (xin[sq] >= toobig) 
	   stop = 1;
	 else {
	   xwsum += xin[sq] * wghts[sq];
	   wsum += wghts[sq];
	   ++nobs;
	 }
      }
   } /* end while */
   
   if (nobs == 0) 
      return (empty_val);   
   
   
    return (xwsum / wsum);    
  

} /* end smooth2d() */
 
/******************************************************************/
double *laplacian(double *in, int nx, int ny, double weight, double empty_val, double mask_val, int report)

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
{
   int row, col, sq, nsq, index;
   int nadj[3], okay;
   double coef1, zadj, z1m, zij, z1p;
   double adj_sum[3], z_del;
   double *out;
   double toobig;

/* define a value which can be used to test for empty/masked nodes */

   toobig = ABS(empty_val) > ABS(mask_val) ? ABS(mask_val) : ABS(empty_val);
   toobig -= 1.0;
   
   nsq = nx * ny;
   out = (double *) get_memory((void *)NULL, (size_t)nsq, sizeof(double));
   coef1 = weight * 0.25;
   nadj[0] = nadj[1] = nadj[2] = 0;
   adj_sum[0] = adj_sum[1] = adj_sum[2] = 0.0;
   
   for (col =0; col < nx; ++col) {
      for (row = 0; row < ny; ++row) {
      
        sq = row * nx + col;
         zij = in[sq];
	 out[sq] = in[sq];
      
         /* Is point missing? */
         if (ABS(zij) >= toobig) 
            continue;
	 
         /* Does the filter stay inside the x boundary? 
          * If not use a forward second difference. 
          */
         z_del = 0.0;
         if (col == 0) {         /* first column */
           index = sq;
           z1m = in[index];
	   zij = in[++index];
	   z1p = in[++index];
         }
         else if (col == nx-1) {  /* last column */
           index = sq;
           z1p = in[index];
           zij = in[--index];
           z1m = in[--index];
         }
         else {    /* interior column */
           index = sq;
           z1m = in[index -1];
           z1p = in[index +1];
         }
      
        /* If any x filter points are missing, set them to point being filtered */
     
        if (ABS(zij) >= toobig)
          zij = in[sq];
        if (ABS(z1p) >= toobig )
          z1p = in[sq];
        if (ABS(z1m) >= toobig )
          z1m = in[sq];
     
	  zadj = z1m + z1p - 2.0*zij;
	  z_del = zadj;
	  ++nadj[0];
	  adj_sum[0] += zadj;
	
         /* Does the filter stay inside the y boundary? 
          * If not use a forward second difference. 
          */
	
         zij = in[sq];
	 if (row == 0) {          /* first row */
	   index = sq;
	   z1m = in[index];
	   zij = in[index + nx];
	   z1p = in[index + nx + nx];
	 }
	 else if (row == ny-1) {
	   index = sq;           /* last row */
	   z1p = in[index];
	   zij = in[index - nx];
	   z1m = in[index - 2*nx];
	 }
	 else {
	   index = sq;
	   z1m = in[index - nx];
	   z1p = in[index + nx];
	 
	 }
        /* If any y filter points are missing, set them to point being filtered */
     
        if (ABS(zij) >= toobig )
          zij = in[sq];
        if (ABS(z1p) >= toobig )
          z1p = in[sq];
        if (ABS(z1m) >= toobig )
          z1m = in[sq];
     
	  zadj = z1m + z1p - 2.0*zij;
	  z_del += zadj;
	  ++nadj[1];
	  adj_sum[1] += zadj;
	
	  ++nadj[2];
	  adj_sum[2] += z_del;
	  
	  out[sq] = in[sq] + coef1 * z_del;
      
      }  /* end for row */
   } /* end for col */

  if (report && nadj[0] >0. && nadj[1] > 0.0) {
    fprintf(stderr,"\n%d points smoothed with average delta = %.5lg\n", nadj[2], (double)(adj_sum[2]*coef1/nadj[2]));
    fprintf(stderr,"avg x-delta: %.5lg   avg y-delta: %.5lg\n",  (double)(adj_sum[0]*coef1/nadj[0]), (double)(adj_sum[1]*coef1/nadj[1]));
    
  }
  return(out);
}
