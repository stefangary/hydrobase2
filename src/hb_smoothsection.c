/* hb_smoothsection.c

................................................................................
                          *******  HydroBase2  *******
               ...................................................
                          
                    author:  Ruth  Curry
                             Woods Hole Oceanographic Institution
                             2000
................................................................................
--------------------------------------------------------------------

 * Smooths a 2-dimensional gridded surface using a distance weighted 
 * gaussian filter with variable x:y weighting ratio.
 *
 *
 *-------------------------------------------------------------------- 
*/

#include <stdio.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "hb_memory.h"
#include "hb_filters.h"
#include "hb_fit.h"

/* globally defined variables */

/* Define how input empty/masked nodes are designated. */

int xradius, yradius;
double xwght, ywght;
double empty_val, mask_val;     /* for output empty/masked nodes */

/*  prototypes for locally declared functions */

void print_usage(char *);
int read_data(int, char **, struct GRID_INFO *, double *, int);
void get_weights( double *, struct GRID_INFO *, int, int);

main (int argc, char **argv)
{
  int bflag, iflag, nfiles;
  int i, j, n;
  int nsq, sq, ix, iy;
  int nwsq, wsq, wrow, wcol;

  BOOLEAN error;
  BOOLEAN set_empty, latfirst;
  BOOLEAN use_mask;
  
  float  x, y, wx, wy;

  double  zout, flagged;
  double xoffset, yoffset;
  double X1, Y1;
  double *z, *zwork, *weights;
  
  char *st, **sptr;
  char *mask;

  FILE *outfile;
  FILE *maskfile;
  
  struct GRID_INFO h, hwork;
  

/* Arguments? */
  
  if (argc < 1) {
    print_usage(argv[0]);
    exit(0);
  }
  
/* Set some default values. */


  empty_val = -TOOBIG;    /* output values for empty/masked nodes */       
  mask_val = TOOBIG;
  flagged = TOOBIG - 10.0;
  use_mask  = FALSE;
  error = 0;
  latfirst = FALSE;
  bflag = iflag = 0;
  outfile = stdout;
  h.node_offset = 0;  /* default gridnode registration */
  nfiles = 0;
  xwght = ywght = 1.0;


  /* parse command line arguments */

  
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      
        case 'B':  /* get grid bounds */
	  bflag = 1;
          st = &argv[i][2];
          if (*st == '/')
             ++st;
          error = (sscanf(st,"%lf", &h.x_min) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%lf", &h.x_max) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%lf", &h.y_min) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%lf", &h.y_max) != 1);
                        	     
	  if (h.x_min > h.x_max) {
	    fprintf(stderr,"\nWest bound must be numerically <= east bound");
	    exit(1);
	  }
	  
	  if (h.y_min > h.y_max) {
	    fprintf(stderr,"\nNorth bound must be numerically <= south bound");
	    exit(1);
	  }
          break;
	  
        case 'G':
	  h.node_offset = 0;
	  break;

        case 'I':
	  iflag = 1;
          error = (sscanf(&argv[i][2], "%lf", &h.x_inc) == 1) ? 0 : 1;
	  h.y_inc = h.x_inc;
          st = strchr(&argv[i][2],'/');
          if (st != NULL) {
             sscanf(++st,"%lf", &h.y_inc);
          }
          break;

        case 'M':
          maskfile = fopen(&argv[i][2],"r");
          if (maskfile == NULL) {
                fprintf(stderr,"\nError opening %s for reading.\n",&argv[i][2]);
                exit(1);
          }
          use_mask = TRUE;
          break;
	  
        case 'O':
          outfile = fopen(&argv[i][2],"w");
          if (outfile == NULL) {
                fprintf(stderr,"\nError opening %s for output.\n",&argv[i][2]);
                exit(1);
          }
          break;
	  
        case 'P':
	  h.node_offset = 1;
	  break;

        case 'S':
          error = sscanf(&argv[i][2], "%d", &xradius) != 1;
	  st = strchr(&argv[i][2], '/');
	  yradius = xradius;
	  if (st != NULL) 
	      sscanf(++st, "%d", &yradius);
          break;

        case 'W':
          error = sscanf(&argv[i][2], "%d", &xwght) != 1;
	  st = strchr(&argv[i][2], '/');
	  ywght = xwght;
	  if (st != NULL) 
	      sscanf(++st, "%d", &ywght);
          break;

	case 'h':  
	   print_usage(argv[0]);
	   exit(0);
	   
	case ':':
	   latfirst = 1;
	   break;

        default:
          error = TRUE;
          break;
        } /* end switch */
	  
       if (error ) {
         fprintf(stderr,"\nError parsing command line args.\n");
         fprintf(stderr,"     in particular: '%s'\n", argv[i]);
         exit(1);
       }

     }
     else 
        ++nfiles;
	
    }  /* end for  */
    
    
   if (!bflag || !iflag ) {
      fprintf(stderr,"\nYou must specify -B<bounds> and -I<gridspacing>");
      fprintf(stderr,"\nUse -h for complete usage info. \n");
      exit(1);
   }
    
  /* Set X1, Y1 to  lower/left gridnode.   */

  xoffset = 0.5*h.x_inc;
  yoffset = 0.5*h.y_inc;
  h.nx = (int) NINT((h.x_max - h.x_min) / h.x_inc);
  h.ny = (int) NINT((h.y_max - h.y_min) / h.y_inc);

/* for gridnode registration... */
 
  if (! h.node_offset) {   /*xmin/xmax coincides with lower/left gridpt */
    ++h.nx;
    ++h.ny;
    xoffset = 0.0;
    yoffset = 0.0;
  }

  X1 = h.x_min + xoffset;           /* west gridnode */
  Y1 = h.y_min + yoffset;           /* south gridnode */

  fprintf (stderr, "\nGrid dimensions are nx = %d, ny = %d", h.nx, h.ny);
  if (h.node_offset)
     fprintf (stderr, "\n using %s ", "pixel");
  else
     fprintf (stderr, "\n using %s ", "gridnode");
  fprintf (stderr, "registration.");


/* allocate space, initialize to empty ... */

  nsq = h.nx * h.ny;
  z = (double *) get_memory ((void *)NULL, (size_t)(nsq), sizeof(double));
  
  for (i = 0; i < nsq; ++i)
     z[i] = empty_val;
  
/* get masked points  ... */

  if (use_mask) {
     mask = (char *) get_memory((void *)NULL, (size_t)nsq, sizeof(char *));
     fprintf(stderr,"\nReading mask file....");
     n = get_mask(maskfile, &h, mask, latfirst);
     fprintf(stderr,"   %d gridnodes will be masked.", n);
     for (i = 0; i < nsq; i++)  {
       if (mask[i] )
            z[i] = mask_val;
     }
     free((void *)mask);
  }
  /* Read in xyz files */  
  n = read_data(nfiles, argv, &h, z, latfirst);  
  fprintf (stderr, "\nFinished input:  inserted %d non-empty points", n);
    
 /* set up work grid and weights*/
 
   hwork.nx = xradius * 2 + 1;
   hwork.ny = yradius * 2 + 1; 
   nwsq = hwork.nx * hwork.ny;
   hwork.x_inc = h.x_inc;
   hwork.y_inc = h.y_inc;
   hwork.node_offset = h.node_offset;
   
   zwork = (double *) get_memory((void *) NULL, (size_t)nwsq, sizeof(double));
   weights = (double *) get_memory((void *) NULL, (size_t)nwsq, sizeof(double));
   
   get_weights(weights, &hwork, xradius, yradius);
  
  /*visit each square in big grid, smooth nodes with actual values and output x,y,z triplets */
  
  for (i = 0; i < h.nx; ++i)  {
     for (j = 0; j < h.ny; ++j ) {
	sq = j * h.nx + i;
	/* reset these at each gridnode */
        error = ij2xy(&h, i, j, &x, &y);
	hwork.x_min =  x - xradius * hwork.x_inc;
	hwork.x_max =  x + xradius * hwork.x_inc;
	hwork.y_min =  y - yradius * hwork.y_inc;
	hwork.y_max =  y + yradius * hwork.y_inc;
	
	for (wrow = 0; wrow < hwork.ny; ++wrow) {
	    for (wcol = 0; wcol < hwork.nx; ++wcol) {
	       wsq = wrow * hwork.nx + wcol;
	       error = ij2xy(&hwork, wcol, wrow, &wx, &wy);
	          error = xy2ij(&h, wx, wy, &ix, &iy);
		  if (error < 0)
		     zwork[wsq] = mask_val;
		  else {
		     zwork[wsq] = z[iy * h.nx + ix];
		  }
	    }  /* for wcol */
	} /* for wrow */
	
	zout = filter2d(zwork, weights, empty_val, mask_val, xradius, yradius, hwork.nx);
	if (ABS(zout) < flagged)
	    fprintf(outfile,"%f %f %lf\n", x, y, zout);
     
     }  /* end for j */
  }  /* end for i*/

  free ((void *)z);
  fprintf(stderr,"\nEnd of %s\n", argv[0]);
  exit(0);
}  /* end main */

/****************************************************************/
void print_usage(char *program)
{
  
    fprintf (stderr, "%s - Smooths a 2D vertical section using a distance weighted gaussian filter\n\n", program);
    fprintf(stderr, "USAGE: %s [xyzfile(s)] -B<xmin/xmax/ymin/ymax> -I<dx>[/<dy>]", program);
    fprintf(stderr, " [-O<output_file>]  [-G] ");
    fprintf(stderr, "  [-M<mask_file> ] [-O<output_file>] [-P] [-S<xrad/yrad> ] [-W<weight_ratio>");
    fprintf(stderr, "  [-:] [-h]\n\n");
    fprintf(stderr, "-B   sets the grid bounds in user units: xmin/xmax/ymin/ymax.\n");
    fprintf(stderr, "-I   sets the grid spacing for the x/y directions.\n");
    fprintf(stderr, "\n\tOPTIONS:\n");

    fprintf(stderr, "-G  force gridnode registration (default ).\n");
    fprintf(stderr, "-M name of file containing masking info.\n");
    fprintf(stderr, "      Each line of file may specify individual points to mask\n");
    fprintf(stderr, "      or multiple polygons separated by a '>' character. \n");
    fprintf(stderr, "      Specify 'I' or 'O' immediately after the '>' character\n");
    fprintf(stderr, "      to mask INSIDE or OUTSIDE of polygon.\n");
    fprintf(stderr, "      No '>' or a 'C' immediately after the '>'  \n");
    fprintf(stderr, "      signifies cell mode masking.  In this mode\n");
    fprintf(stderr, "      an {x,y} pair OR {x,y,mask} triplet is given.\n");
    fprintf(stderr, "      The <mask> value can be 1 to mask the point\n");
    fprintf(stderr, "      or 0 (zero) to unmask the point.\n");
    fprintf(stderr, "-O  name of output file  (default is stdout).\n");
    fprintf(stderr, "-P  force pixel registration (default is gridnode registration).\n");
    fprintf(stderr, "-S xradius/yradius of smoothing ellipse\n");
    fprintf(stderr, "-W x:y weight ratio.  Ex:  -W3/1 for weighting x direction 3 times heavier than y\n");
    fprintf(stderr, "      Default: [%1.1lf/%1.1lf]\n", xwght, ywght);
    fprintf(stderr, "-: input data are ordered y x z  \n");
    fprintf(stderr, "       [default order is x y z]   \n");
    fprintf(stderr, "-h help....prints this message.  \n");
    return;

} /* end print_usage() */
/****************************************************************/
int read_data(int nfiles, char **argv, struct GRID_INFO *hptr, double *z, int latfirst)

/*  On input, z is dimensioned hptr->nx * ny and must be initialized 
   everywhere to 0.  Data are read into the appropriate grid square.
   Multiple entries for the same square are averaged. */
{  
  int n, nsq, curfile, error;
  int n_fields, ix, iy;
  int row, col, sq;
  int *count;
  
  double in[3], X1, Y1;
  double flagged;
  
  char *st, line[BUFSIZ];
  FILE *infile;
  
  ix = latfirst; 
  iy = 1 - ix;
  flagged = mask_val - 10.0;
  
  infile = stdin;
  curfile = 1;

  X1 = hptr->x_min;
  Y1 = hptr->y_min;
    
  if(hptr->node_offset) {    /* lower/left gridnode for pixel grid */
    X1 += hptr->x_inc * 0.5;
    Y1 += hptr->y_inc * 0.5;
  } 
   
  nsq = hptr->nx * hptr->ny;
  count = (int *) get_memory((void *)NULL, (size_t)nsq, sizeof(int));     
  n = 0;  /* counts number of points read in */

 /* loop for each input file */
 
   do {
   
  
     if (! nfiles)
        fprintf(stderr,"\nExpecting input from stdin ...");
      
     else {
        infile = fopen(argv[curfile],"r");
        if (infile == NULL) {
         fprintf(stderr, "\nUnable to open %s for reading.", argv[curfile]);
         goto NEXTFILE;
        }
        fprintf(stderr,"\nOpened %s ", argv[curfile]);
     }
   

     while (fscanf(infile,"%[^\n]", line) != EOF) { 
        error = getc(infile);   /* move past newline */
	
	if ((st = strchr(line,'N')) != NULL) continue; /* check for NaN */
	
        n_fields = sscanf(line,"%lf %lf %lf", &in[0], &in[1], &in[2]);

        if (n_fields != 3) {
           fprintf(stderr, "Mismatch between actual (%d) and expected (3) fields near line %d \n", n_fields, n);
           continue;
        }

               /* weed out empty/mask flags */
        if (ABS(in[2]) >= flagged) continue; 
	  
        row = (int) NINT((in[iy] - Y1) / hptr->y_inc);  
	col = (int) NINT((in[ix] - X1) / hptr->x_inc); 
	
	if (row >= 0 && row < hptr->ny && col >= 0 && col < hptr->nx) {
	   sq = row * hptr->nx + col;
	   if (count[sq] == 0)
	     z[sq] = 0;
	   z[sq] += in[2];
	   ++count[sq];
           ++n;
        }  
	
    }  /* End while  */
      
NEXTFILE:
    if (nfiles) fclose (infile);

  } while (++curfile < nfiles);      /* End input phase */


/* If multiple entries for a node, find the average.  If no
   entries, mark it empty. */
     
  for (sq = 0; sq < nsq; ++sq) { 
     if (count[sq] > 1) 
        z[sq] /= count[sq];
  }
  
  free((void *)count);
  return(n);
  
}  /* end read_data() */

/**********************************************************/
void get_weights(double *weight, struct GRID_INFO *hptr, int xcntr, int ycntr)
  /* sets weight based on distance (grid_node) from center node
        weight[n] = e^[- (pi/halfwidth)^2 /4.5 * dist^2 ].  If xwght/ywght ratio is specified,
	multiplies weight in x or y direction by ratio*/
{

int ix, iy, sq;
double ratio, dist, halfwidth, alpha;

   halfwidth = xcntr;
   if (ycntr > xcntr)
      halfwidth = ycntr;
   for (ix = 0; ix < hptr->nx; ++ix) {
      for (iy = 0; iy < hptr->ny; ++iy) {
          sq = iy * hptr->nx + ix;
	  dist = sqrt((iy-ycntr) * (iy - ycntr) + (ix - xcntr) * (ix - xcntr));
	  alpha = 3.14 / halfwidth;
	  alpha  *= alpha;
	  alpha /= 4.5;
	  weight[sq] = exp( -alpha * dist * dist);
      }
   }
   ratio = xwght / ywght;
   if (xwght > ywght) {
      for (ix = 0; ix < hptr->nx; ++ix) {
         sq = ycntr * hptr->nx + ix;
         weight[sq] *= ratio;
      }
      return;
   }
   if (ywght > xwght) {
      for (iy = 0; iy < hptr->nx; ++iy) {
         sq = iy * hptr->nx + xcntr;
         weight[sq] *= ratio;
      }
      return;
   
   }

   return;
} /* end get_weights() */
