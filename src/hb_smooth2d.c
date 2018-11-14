/* hb_smooth2d.c

................................................................................
                          *******  HydroBase2  *******
               ...................................................
                          
                    author:  Ruth  Curry
                             Woods Hole Oceanographic Institution
                             2000
................................................................................
--------------------------------------------------------------------

 * Smooths a 2-dimensional gridded surface by iteratively applying 
 *   a Laplacian 5-pt filter.
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

double HB_d_mask = TOOBIG;      /* a very large  value */
double HB_d_empty = -(TOOBIG);  /* a very negative value */
double flagged;               /* test value for input masked/empty nodes */


double weight;     
int n_iter;

/*  prototypes for locally declared functions */

void print_usage(char *);
int read_data(int, char **, struct GRID_INFO *, double *, int);

main (int argc, char **argv)
{
  int bflag, iflag, nfiles;
  int i, n;
  int ix, iy;
  int nsq, row, col;

  BOOLEAN error;
  BOOLEAN set_empty, latfirst;
  BOOLEAN use_mask;
  BOOLEAN mask_NaN, empty_NaN;
  
  double empty_val, mask_val;     /* for output empty/masked nodes */
  double xout, yout;
  double xoffset, yoffset;
  double X1, Y1;
  double *Z, *zout;
  
  char *st, **sptr;
  char *mask;

  FILE *outfile;
  FILE *maskfile;
  FILE *new_maskfile;
  
  struct GRID_INFO h;
  

/* Arguments? */
  
  if (argc < 1) {
    print_usage(argv[0]);
    exit(0);
  }
  
/* Set some default values. */


  flagged = HB_d_mask - 10.0;
  empty_val = -TOOBIG;    /* output values for empty/masked nodes */       
  mask_val = TOOBIG;
  use_mask  = FALSE;
  set_empty = FALSE;
  mask_NaN = empty_NaN = FALSE;
  error = 0;
  latfirst = FALSE;
  bflag = iflag = 0;
  outfile = stdout;
  new_maskfile = (FILE *)NULL;
  h.node_offset = 1;  /* default pixel registration */
  
  weight = 0.25;
  n_iter = 3;
  nfiles = 0;


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
	  
        case 'E':
	  st = &argv[i][2];  
          set_empty = TRUE;
	  
	  if (*st == '\0')     /* no options specified */
	     break;
	  
	  if (*st != '/') {     /* set option for empty values */
	     if ((*st == 'N') || (*st == 'n')) {
	       empty_NaN = TRUE;
	     }
	     else {
	       fprintf(stderr,"\nUnknown qualifier -E%s", st);
	       fprintf(stderr,"\nuse -Enan  or -E/nan or -Enan/nan");
	       fprintf(stderr,"\nto use NaN instead of ");
	       fprintf(stderr,"\n [%.3e] for empty values ", empty_val);
	       fprintf(stderr,"\n [%.3e] for mask values ", mask_val);
	       exit(1);
	     }
	     
	  }
	  
	  
	  if ((st = strchr(&argv[i][2],'/')) != NULL) {
	    ++st;
	     if ((*st == 'N') || (*st == 'n')) {
	       mask_NaN = TRUE;
	     }
	    else {
	       fprintf(stderr,"\nUnknown mask qualifier %s", argv[i]);
	       fprintf(stderr,"\nuse -Enan  or -E/nan or -Enan/nan");
	       fprintf(stderr,"\nto use NaN instead of ");
	       fprintf(stderr,"\n [%.3e] for empty values ", empty_val);
	       fprintf(stderr,"\n [%.3e] for mask values ", mask_val);
	       exit(1);
	    }
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
	  
        case 'N':
          new_maskfile = fopen(&argv[i][2],"w");
          if (new_maskfile == NULL) {
                fprintf(stderr,"\nError opening %s for output.\n",&argv[i][2]);
                exit(1);
          }
          break;
                    
        case 'O':
          outfile = fopen(&argv[i][2],"w");
          if (outfile == NULL) {
                fprintf(stderr,"\nError opening %s for output.\n",&argv[i][2]);
                exit(1);
          }
          break;

        case 'S':
          error = sscanf(&argv[i][2],"%d", &n_iter) != 1;
	  
          break;

        case 'W':
          error = sscanf(&argv[i][2],"%lf", &weight) != 1;
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
    
   error = 0;
  
  if (n_iter < 1) {
    fprintf (stderr, "SYNTAX ERROR -S option.  Must specify interations >= 1\n");
    error++;
    }
    
  if (weight <= 0.0) {
    fprintf (stderr, "SYNTAX ERROR -W option.  Must specify weight > 0\n");
    error++;
    }
    
  if (error) 
        exit(1);
    
    
  /* Set X1, Y1 to  lower/left gridnode.  
     default is for pixel grids: 
     Leave the struct GRID_INFO as is. */

  xoffset = 0.5*h.x_inc;
  yoffset = 0.5*h.y_inc;
  h.nx = (int) NINT((h.x_max - h.x_min) / h.x_inc);
  h.ny = (int) NINT((h.y_max - h.y_min) / h.y_inc);

/* for gridnode registration... */
 
  if (! h.node_offset) {   /*xmin/xmax coincides with lower/left gridpt */         ++h.nx;
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


/* allocate space and get input xyz data ... */

  nsq = h.nx * h.ny;
  Z = (double *) get_memory ((void *)NULL, (size_t)(nsq), sizeof(double));
  
  n = read_data(nfiles, argv, &h, Z, latfirst);  
  fprintf (stderr, "\nFinished input:  inserted %d non-empty points", n);
    
  
/* get masked points  ... */

  if (use_mask) {
     mask = (char *) get_memory((void *)NULL, (size_t)(nsq), sizeof(char *));
     fprintf(stderr,"\nReading mask file....");
     n = get_mask(maskfile, &h, mask, latfirst);
     fprintf(stderr,"   %d gridnodes will be masked.", n);
     for (i = 0; i < nsq; i++)  {
       if (mask[i] )
            Z[i] = HB_d_mask;
     }
     free((void *)mask);
  }
  
/* Iteratively smooth the Z array */

  for (i = 0; i < n_iter; ++i) {
    fprintf(stderr,"\nSmoothing pass #%d ", i);
    zout = laplacian(Z, h.nx, h.ny, weight, HB_d_empty, HB_d_mask, TRUE);
    free((void *)Z);
    Z = zout;
  }
  
/* Output the new grid, setting empty and masked gridnodes to 
   appropriate values ... */


  for (i = 0;  i < nsq; i++)  {
         
    row = i / h.nx;
    col = i - row * h.nx;
    xout =  col * h.x_inc + X1;
    yout =  row * h.y_inc + Y1;
    
    if (ABS(Z[i]) >= flagged )  {  /* Is it somehow flagged? */
       if (set_empty) {
         if (Z[i] >= flagged) {
	    if (mask_NaN) 
	       fprintf(outfile,"%8.3lf %8.3lf NaN \n", xout, yout);
	    else
	       fprintf(outfile,"%8.3lf %8.3lf %.1e \n", xout, yout, mask_val);
	 }
	 else {
	    if (empty_NaN) 
	       fprintf(outfile,"%8.3lf %8.3lf NaN \n", xout, yout);
	    else
	       fprintf(outfile,"%8.3lf %8.3lf %.1e \n", xout, yout, empty_val);
	 }
       }
       
       if (new_maskfile)
         fprintf(new_maskfile, "%8.3lf %8.3lf \n", xout, yout);
    }
    else {
       fprintf(outfile,"%8.3lf %8.3lf %.9g \n", xout, yout, Z[i]);
    
    }
  }  /* end for */

  free ((void *)Z);
  fprintf(stderr,"\nEnd of %s\n", argv[0]);
  exit(0);
}  /* end main */

/****************************************************************/
void print_usage(char *program)
{
  
    fprintf (stderr, "%s - An iterative 5-pt Laplacian filter for 2-D gridded data.\n\n", program);
    fprintf(stderr, "USAGE: %s [xyzfile(s)] -B<xmin/xmax/ymin/ymax> -I<dx>[/<dy>]", program);
    fprintf(stderr, " [-O<output_file>] [-E<empty>/<mask_val>] [-G] ");
    fprintf(stderr, "  [-M<mask_file> ] [-N<new_mask_file> ] [-O<output_file>]  [-S<#_smoothing_iterations> ] [-W<weight>]");
    fprintf(stderr, "  [-:] [-h]\n\n");
    fprintf(stderr, "-B   sets the grid bounds in user units: xmin/xmax/ymin/ymax.\n");
    fprintf(stderr, "-I   sets the grid spacing for the x/y directions.\n");
    fprintf(stderr, "\n\tOPTIONS:\n");
    fprintf(stderr, "-E  Explicity write out nodes that are empty or masked.\n");
    fprintf(stderr, "       Without -E, only non-empty, non-masked gridnodes are written out.\n");
    fprintf(stderr, "       Optionally use NaN for empty nodes ");
    fprintf(stderr, "instead of [%.3e] \n", HB_d_empty);
    fprintf(stderr, "       and/or for masked nodes instead of [%.3e]\n", HB_d_mask);
    fprintf(stderr, "       Specify -Enan  for empty nodes\n");
    fprintf(stderr, "               -E/nan  for masked nodes\n");
    fprintf(stderr, "               -Enan/nan for both empty and masked.\n");

    fprintf(stderr, "-G  force gridnode registration (default is pixel registration).\n");
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
    fprintf(stderr, "-N name of new file for x,y pairs of empty or masked gridnodes after smoothing.\n");
    fprintf(stderr, "-O  name of output file  (default is stdout).\n");
    fprintf(stderr, "-S #of iterations for smoothing\n");
    fprintf(stderr, "      Default is [-S%1d]\n", n_iter);
    fprintf(stderr, "-W weighting factor. Default:[%4.2lf]\n", weight);
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
  
  char *st, line[BUFSIZ];
  FILE *infile;
  
  ix = latfirst; 
  iy = 1 - ix;
  
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
     else if (count[sq] == 0)
        z[sq] = HB_d_empty;
  }
  
  free((void *)count);
  return(n);
  
}  /* end read_data() */
