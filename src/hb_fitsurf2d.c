/*  hb_fitsurf2d.c

................................................................................
                          *******  HydroBase 2 *******
                             Ruth Curry
                             Woods Hole Oceanographic Institution
                             Dec 2000
...................................................

*  Interpolates for missing gridnodes in a 2-dimensional surface using an
*  iterative Laplacian/spline fitting algorithm.  Reads xyz files.  If a 
*  masking file is specified, those gridnodes are excluded from the fit.                     
   
 */

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include "hb_fit.h"
#include "hb_memory.h"

/* global variables */
 
int max_iterations;		
int verbose;
int xrad, yrad;	           /* Search radius for empty nodes  */
double converge_limit;	   /* Convergence limit */
double	tension = 0.0;	    /* Tension parameter on the surface  */
double HB_d_mask;      
double HB_d_empty;
double flagged;        /* test for input xyz data values */ 


/* prototypes for locally defined functions */	

void print_usage(char *);
struct POINT *read_data(int, char **, struct GRID_INFO *, int *, int);

main(int argc, char **argv)
{
  int i, j, npoints;
  int nsqout, nsqwork;
  int n_set, n_mask, n_empty;
  int nfiles = 0;
  int row, col, error;
  int startrow, endrow, startcol, endcol;
  int bflag, iflag;
  int xadd, yadd; 
  BOOLEAN set_empty, latfirst;
  BOOLEAN use_mask, pixel_grid;
  BOOLEAN mask_NaN, empty_NaN;
  double empty_val, mask_val;
  double xoffset, yoffset;
  double xradius, yradius;
  double lat, lon;
  double toobig;          /* zgrid() representation of empty/masked nodes */
  double bigger;
  double HB_zgrid_mask;  /* zgrid() expects this on input to mark masked nodes */
  int *knxt, *imnew;
  double *Z, *zpij, *new_Z;
  char *st;
  char *mask;
  struct GRID_INFO hwork, hout;
  FILE	*outfile, *new_maskfile, *maskfile;
  struct POINT *data; 

/* Set these default values. */

  toobig = 1.0E34;        /* flag to test for (<) zgrid() ouput flag */	
  bigger = 1.0E36;        /* larger than empty flag/smaller than mask flag */
  flagged = 1.0E34;       /* flag to test input: ABS(input) < flag if valid */
  HB_zgrid_mask = (double) DBL_MAX;
  HB_d_mask = TOOBIG;     /* a very large value */
  HB_d_empty = -(TOOBIG); /* a very negative value */
  empty_val = HB_d_empty; /* value to write out */            
  mask_val = HB_d_mask;   /* value to write out */
  mask_NaN = empty_NaN = FALSE;  /* write out NaN instead of above */
  use_mask  = FALSE;
  set_empty = FALSE;
  
  verbose = FALSE;
  error = 0;
  latfirst = FALSE;
  bflag = iflag = 0;
  max_iterations = 1;
  n_mask = 0;
  tension = 1.0;
  converge_limit = 0.0000001;   /* something very small */
  xrad = yrad = 1;     /* Search radii in gridnodes */
  outfile = stdout;
  new_maskfile = (FILE *)NULL;
  pixel_grid = TRUE;          /* default pixel registration */
 
  
  if (argc < 1) {
    print_usage(argv[0]);
    exit(0);
  }


  /* parse command line arguments */

  
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      
        case 'B':  /* get grid bounds */
	  bflag = 1;
          st = &argv[i][2];
          if (*st == '/')
             ++st;
          error = (sscanf(st,"%lf", &hout.x_min) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%lf", &hout.x_max) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%lf", &hout.y_min) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%lf", &hout.y_max) != 1);
                        	     
	  if (hout.x_min > hout.x_max) {
	    fprintf(stderr,"\nWest bound must be numerically <= east bound");
	    exit(1);
	  }
	  
	  if (hout.y_min > hout.y_max) {
	    fprintf(stderr,"\nNorth bound must be numerically <= south bound");
	    exit(1);
	  }
          break;
	  
        case 'C':
          error = sscanf(&argv[i][2],"%lf", &converge_limit) != 1;
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
	  pixel_grid = FALSE;
	  break;

        case 'I':
	  iflag = 1;
          error = (sscanf(&argv[i][2], "%lf", &hout.x_inc) == 1) ? 0 : 1;
	  hout.y_inc = hout.x_inc;
          st = strchr(&argv[i][2],'/');
          if (st != NULL) {
             sscanf(++st,"%lf", &hout.y_inc);
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

        case 'Q':
          error = sscanf(&argv[i][2],"%d", &max_iterations) != 1;
          break;

        case 'S':
          error = sscanf(&argv[i][2],"%d", &xrad) != 1;
          yrad = xrad;
          st = &argv[i][2];
          while (*(st++) != '\0') {
            if (*st == '/') {
               ++st;
               error = (sscanf(st,"%d", &yrad) == 1) ? 0 : 1;
               break;
            }
          }
          break;
	  
        case 'T':
	  tension = atof(&argv[i][2]);
          break;
	  
        case 'V':
	  verbose = TRUE;
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
    
/*  Check syntax of options */    
                              
error = 0;

if (tension < 0.0 || tension > 1.0) {
    fprintf(stderr,"ERROR:  specify tension parameter between 0->1\n");
    ++error;
}


if (max_iterations < 1) {
    fprintf (stderr, "SYNTAX ERROR -Q option.  Must specify interations >= 1\n");
    error++;
    }
if (converge_limit <= 0.0) {
    fprintf (stderr, "SYNTAX ERROR -C option.  Must specify convergence limit > 0\n");
    error++;
    }
if (xrad < 1 || yrad < 1) {
    fprintf (stderr, "SYNTAX ERROR -S option: search radius (gridnodes) must be >= 1.\n");
    error++;
    }
    
if (error) 
        exit(1);
 
        
/* set xmin/ymin to coordinates of lower/left gridnode 
   and xmax/ymax to upper/left gridnode so we don't have
   to keep checking for grid registration.   */

xoffset = 0.5*hout.x_inc;    /* default is for pixel registration */
yoffset = 0.5*hout.y_inc;
hout.nx = (int) NINT((hout.x_max - hout.x_min) / hout.x_inc);
hout.ny = (int) NINT((hout.y_max - hout.y_min) / hout.y_inc);

if (!pixel_grid) {  /* for forced gridnode registration */
  xoffset = 0.0;
  yoffset = 0.0;
  ++hout.nx;
  ++hout.ny;
}

hout.x_min +=  xoffset;                       /* west gridnode */
hout.x_max = hout.x_min + (hout.nx -1) * hout.x_inc;  /* east  */
hout.y_min += yoffset;                                /* south */
hout.y_max = hout.y_min + (hout.ny -1) * hout.y_inc;  /* north */

hout.node_offset = 0;   /* now normalized to gridnode registration */

nsqout = hout.nx * hout.ny;

/* Set up work grid :   xmin/ymin is lower/left gridnode
    and xmax/ymax is upper/right gridnode determined by the value of 
    node_offset.  expand work area by <xrad/yrad> gridnode units. 
*/

   xradius = xrad * hout.x_inc;    /* radius in x-units */
   yradius = yrad * hout.y_inc;    /* radius in x-units */
        
   hwork.x_min = hout.x_min - xradius;
   hwork.x_max = hout.x_max + xradius;
   hwork.y_min = hout.y_min - yradius;
   hwork.y_max = hout.y_max + yradius;
   hwork.x_inc = hout.x_inc;
   hwork.y_inc = hout.y_inc;
   
   xadd = xrad;   
   yadd = yrad;  
   hwork.nx = hout.nx + 2 * xadd;
   hwork.ny = hout.ny + 2 * yadd;
   nsqwork =  hwork.nx * hwork.ny; 
   hwork.node_offset = 0;      /* normalized to gridnode registration */

   fprintf (stderr, "\nOutput xmin/ymin gridnode: %8.3lf/%8.3lf ", hout.x_min, hout.y_min);
   
   fprintf (stderr, " grid dimensions are nx = %d, ny = %d\n", hout.nx, hout.ny);

  if (pixel_grid)
     fprintf (stderr, "\n using pixel registration");
  else
     fprintf (stderr, "\n using gridnode registration");

   fprintf (stderr, "\nNumber of gridnodes in search radii: %d/%d", xrad, yrad);
	
/*-------------------------------------------*/
   
  data = read_data(nfiles, argv, &hwork, &npoints, latfirst);
  

  /* allocate space for work arrays */
  
  knxt = (int *) get_memory ((void *)NULL, npoints+1, sizeof (int));
  zpij = (double *) get_memory ((void *)NULL, npoints+1, sizeof (double));
  imnew = (int *) get_memory ((void *)NULL, (int) hwork.ny, sizeof (int));
  Z = (double *) get_memory ((void *)NULL, (size_t)(nsqwork), sizeof(double));

  if (use_mask) {
     mask = (char *) get_memory((void *)NULL, (size_t)(nsqwork), sizeof(char *));
     fprintf(stderr,"\nReading mask file....");
     n_mask = get_mask(maskfile, &hwork, mask, latfirst);
     for (i = 0; i < nsqwork; i++)  {
       if (mask[i] )
          Z[i] = HB_zgrid_mask;
     }
     free((void *)mask);
     fclose(maskfile);
  }

  fprintf(stderr,"\ngridding....");
  new_Z = zgrid (converge_limit, max_iterations, npoints, xrad, yrad, tension, &hwork, knxt, imnew, zpij, Z, data, &n_mask, &n_empty, &n_set);
  
  if (n_set < 0)  {
    fprintf(stderr,"\nzgrid() returned error.\n");
    exit (-1);
  }

  free((void *)Z);
  Z = new_Z;

/* Finished with surface fitting.
   Output the new grid, setting empty and masked gridnodes to 
   appropriate values ... */

  n_set = 0;
  n_mask = 0;
  n_empty = 0;

  startrow = yadd;
  endrow = hwork.ny - yadd - 1;
  startcol = xadd;
  endcol = hwork.nx - xadd - 1;
      
  for (row = startrow; row <= endrow; ++row) { 
    lat = hwork.y_min + row * hwork.y_inc;  
    for (col = startcol; col <= endcol; ++col) {
    
      i =  row * hwork.nx + col;           /* index to work grid */
      lon =  hwork.x_min + col * hwork.x_inc;
    
      if (ABS(Z[i]) >= toobig )  {  /* Is it somehow flagged? */
        if (Z[i] > bigger) {
	   ++n_mask;
           if (set_empty) {
	       if (mask_NaN) 
	          fprintf(outfile,"%8.3lf %8.3lf NaN \n", lon, lat);
	       else
	          fprintf(outfile,"%8.3lf %8.3lf %.1e \n", lon, lat, mask_val);
	   }
	}
	else {
	   ++n_empty;
           if (set_empty) {
	       if (empty_NaN) 
	          fprintf(outfile,"%8.3lf %8.3lf NaN \n", lon, lat);
	       else
	          fprintf(outfile,"%8.3lf %8.3lf %.1e \n", lon, lat, empty_val);
	   }
         }
       
         if (new_maskfile)
         fprintf(new_maskfile, "%8.3lf %8.3lf \n", lon, lat);
      }
      else {
       fprintf(outfile,"%8.3lf %8.3lf %.9g \n", lon, lat, Z[i]);
       ++n_set;
      }
    }  /* end for col */ 
  } /* end for row */
  
  free ((void *)Z);
  free ((void *)data);
  free ((void *)knxt);
  free ((void *)zpij);
  free ((void *)imnew);
  
  
  if (verbose) {
     fprintf(stderr,"\nNumber of work grid nodes: <%d>", nsqout);
     fprintf(stderr,"\nNumber filled: <%d>   empty: <%d>   masked: <%d>\n", n_set, n_empty, n_mask );
  
  }
  fprintf(stderr,"\nEnd of %s.\n", argv[0]);
exit(0);
	
}  /* end main */

/************************************************************************/
void print_usage(char *program)
{

  fprintf(stderr,"\nUses a Laplacian/spline algorithm to fit a surface to xyz data.");
  fprintf(stderr,"\n\nUSAGE:  %s xyzfile(s) -B<w/e/s/n> -I<x_inc>[/<y_inc>] [-C<converge_limit>] [-E[nan[/nan]]] [-G] [-M<maskfile>] [-N<new_maskfile>] [-Ooutfile] [-Q<max_iterations>] [-S<xradius>[/<yradius>]] [-T<tension>] [-V] [-h] [-:]\n\n", program);
  fprintf(stderr,"-B  sets the grid bounds w/e/s/n.\n");
  fprintf(stderr,"-I  sets the grid spacing for the x/y directions.\n");
  fprintf(stderr, "\n\tOPTIONS:\n");
  fprintf(stderr,"-C  set the convergence criteria. Default is %lf\n", converge_limit);
  fprintf(stderr,"-E  Explicity write out nodes that are empty or masked.\n");
  fprintf(stderr,"       Without -E, only non-empty, non-masked gridnodes are written out.\n");
  fprintf(stderr,"       Optionally use NaN instead of [%.3e] for empty nodes\n", HB_d_empty);
  fprintf(stderr,"       and/or NaN instead of [%.3e] for masked nodes\n", HB_d_mask);
  fprintf(stderr,"       Specify -Enan  for empty nodes\n");
  fprintf(stderr,"               -E/nan  for masked nodes\n");
  fprintf(stderr,"               -Enan/nan for both empty and masked.\n");

  fprintf(stderr,"-G force gridnode registration (default is pixel registration).\n");
  fprintf(stderr,"-M name of input file containing masking info.\n");
  fprintf(stderr,"      Each line of file may specify individual points to mask\n");
  fprintf(stderr,"      or multiple polygons separated by a '>' character. \n");
  fprintf(stderr,"      Specify 'I' or 'O' immediately after the '>' character\n");
  fprintf(stderr,"      to mask INSIDE or OUTSIDE of polygon.\n");
  fprintf(stderr,"      No '>' or a 'C' immediately after the '>'  \n");
  fprintf(stderr,"      signifies cell mode masking for individual gridnodes.\n");
  fprintf(stderr,"      In this mode an [x,y] pair OR [x,y,mask] triplet\n");
  fprintf(stderr,"      is given. The <mask> value can be 1 to mask the point\n");
  fprintf(stderr,"      or 0 (zero) to unmask the point.\n");
  fprintf(stderr,"-N name of new file for lon/lat pairs of masked gridnodes after the fitting.\n");
  fprintf(stderr,"-O name of output file  (default is stdout).\n");
  fprintf(stderr,"-Q sets the max iterations to achieve convergence.  [%d]\n", max_iterations);
  fprintf(stderr,"-S search radius in integer grid increments\n");
  fprintf(stderr,"      for x and y directions.  If no data are \n");
  fprintf(stderr,"      within range of a node it is set to empty.\n");
  fprintf(stderr,"      If -E is set, the grid points are written to outfile\n");
  fprintf(stderr,"      with z set to -9e34 or NaN. If -E is not set\n");
  fprintf(stderr,"      the grid point is omitted from outfile. If -N is set, the x,y\n");
  fprintf(stderr,"      coordinates are listed in the new_maskfile.\n");
  fprintf(stderr,"      Default is [-S%1d/%1d]\n", xrad, yrad);
  fprintf(stderr,"-T tension parameter -- range [0..1]\n");
  fprintf(stderr,"      A value of 1 produces a pure laplacian solution,\n");
  fprintf(stderr,"      while a 0 value gives a harmonic spline\n");
  fprintf(stderr,"      solution with a smoother field but the possibility of \n");
  fprintf(stderr,"      spurious peaks or valleys.  Default: [%.2lf]\n", tension);
  fprintf(stderr,"-V verbose.  \n");
  fprintf(stderr,"-: all input files (xyz and mask) are ordered lat/lon  \n");
  fprintf(stderr,"      [default order is lon/lat]   \n");
  fprintf(stderr,"-h help...... prints this message. \n");
  return;
  
} /* end print_usage() */
/************************************************************************/
struct POINT * read_data(int nfiles, char **argv, struct GRID_INFO *hptr, int *npts_addr, int yfirst)
{
  int	i, j, nalloc, n_fields, ix, iy;
  int   curfile, error;
  double in[3];
  char	line[BUFSIZ], *st;
  FILE *infile;
  struct POINT *p, *data;
	
	/* Read in xyz data and store it in a structure */

  ix = yfirst; 
  iy = 1 - ix;
  nalloc = 10000;	
  data = (struct POINT *) get_memory((void *)NULL,(size_t)nalloc, sizeof(struct POINT));
  *npts_addr = 0;
  curfile = 1;
  infile = stdin;

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

	if (n_fields <= 2) 
	       continue;
	
	if (ABS(in[2]) > flagged) 
	      continue;   /* signifies an empty or masked node */
		
	i = (int)NINT((in[ix]- hptr->x_min) / hptr->x_inc);
	if (i < 0 || i >= hptr->nx) continue;
	
	j = (int)NINT((in[iy]- hptr->y_min) / hptr->y_inc);
	if (j < 0 || j >= hptr->ny) continue;
        p = &data[*npts_addr];
	p->XP = in[ix];
	p->YP = in[iy];
	p->ZP = in[2];
	if (++(*npts_addr) == nalloc) {
           nalloc += 1000;
           data = (struct POINT *) get_memory ((void *)data, (size_t)nalloc, sizeof(struct POINT));
	}
	
    }  /* End while  */
      
NEXTFILE:
    if (nfiles) fclose (infile);

  } while (++curfile < nfiles);      /* End input phase */
  
  if (*npts_addr == 0) {
	fprintf (stderr, " No datapoints inside region.... exiting.\n");
	exit(1);
  }
	
  data = (struct POINT *) get_memory ((void *)data, (size_t)*npts_addr, sizeof(struct POINT));

  if (verbose) {
     fprintf(stderr, "\n Read in %d xyz points.", *npts_addr);
  }
  return(data);

} /*end read_data() */

