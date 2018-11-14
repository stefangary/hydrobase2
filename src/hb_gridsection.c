/* hb_gridsection.c
 ................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth G Curry
                             Woods Hole Oceanographic Instition
                             2000  ANSI compliant
...............................................................................
--------------------------------------------------------------------
.
 * Interpolates xyz values onto a grid by linear interpolation in the vertical (y)
 * direction first, followed by interpolation in the horizontal (x) direction.   
 *
 *-------------------------------------------------------------------- */

#include <stdio.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "hydrobase.h"
#include "hb_fit.h"
#include "hb_memory.h"


/* globally defined variables */

int  xradius, yradius; 

/* Define internal representation of masked / flagged values    
   The value for empty/masked nodes written out is defined 
   in hb_fit.h as TOOBIG */
   
double empty_val = -TOOBIG;         /* large negative value */
double mask_val  = TOOBIG;          /* large positive value */
double flagged = -TOOBIG +1;        /* check for input missing values  */
double masked = TOOBIG -1;        /* check for input masked values  */

/*  prototypes for locally declared functions */

void print_usage(char *);

main (int argc, char **argv)
{
  int bflag, iflag;
  int i, j,  nread, n_fields;
  int ix, iy, n_set, n_mask, n_empty;
  int nfiles = 0, curfile = 1;
  int nsq, row, col, sq;
  int rad, found;
  int *zcount;
  

  BOOLEAN error;
  BOOLEAN  yfirst;
  BOOLEAN  skip;
  BOOLEAN use_mask;

  float xflt, yflt;
  double in[3], x_left, x_right, y_top, y_bottom; 
  double xoffset, yoffset;
  double x1, x2, y1, y2, z1, z2;
  double *zout, *zwork, *xwork, *ywork;    


  char line[BUFSIZ];
  char *st, **sptr;
  char *mask;
  
  FILE *outfile;
  FILE *infile;
  FILE *maskfile;
  FILE *new_maskfile;
  
  struct GRID_INFO header;
  
  /* initialize these values */

  use_mask  = FALSE;
  error = 0;
  yfirst = 0;
  bflag = iflag = 0;
  n_mask = 0;
  xradius  = 1;
  yradius = 1;
  outfile = stdout;
  infile = stdin;
  new_maskfile = (FILE *)NULL;
  header.node_offset = 0;  /* default gridnode registration */
  mask = (char *)NULL;  
  
  
  if (argc < 1) {
    print_usage(argv[0]);
    exit(1);
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
          error = (sscanf(st,"%lf", &header.x_min) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%lf", &header.x_max) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%lf", &header.y_min) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%lf", &header.y_max) != 1);
                        	     
	  if (header.x_min > header.x_max) {
	    fprintf(stderr,"\nxmin bound must be numerically < xmax bound");
	    exit(1);
	  }
	  
	  if (header.y_min > header.y_max) {
	    fprintf(stderr,"\nymin bound must be numerically < ymax bound");
	    exit(1);
	  }
          break;
	  
        case 'I':
	  iflag = 1;
          error = (sscanf(&argv[i][2], "%lf", &header.x_inc) == 1) ? 0 : 1;
	  header.y_inc = header.x_inc;
          st = strchr(&argv[i][2],'/');
          if (st != NULL) {
             sscanf(++st,"%lf", &header.y_inc);
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

        case 'P':
	  header.node_offset = 1;
	  break;
	  

        case 'S':
          error = sscanf(&argv[i][2],"%d", &xradius) != 1;
          st = strchr(&argv[i][2],'/');
	  yradius = xradius;
          if (st != NULL) {
             sscanf(++st,"%d", &yradius);
          }
          break;
	  
	case 'h':  
	   print_usage(argv[0]);
	   exit(0);
	   
	case ':':
	   yfirst = 1;
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
   
  if (xradius < 0 || yradius < 0) {
    fprintf (stderr, "SYNTAX ERROR -S option.  Search radius must be positive.\n");
    error++;
    }
  if (error) 
        exit(1);
 
  ix = yfirst; 
  iy = 1 - ix;
        
   fprintf (stderr, "search xradius = %d   yradius = %d\n", xradius, yradius);
   xoffset = 0.5*header.x_inc;
   yoffset = 0.5*header.y_inc;
     
  /*  for default gridnode registration */
  
   header.nx = 1 + (int) NINT((header.x_max - header.x_min) / header.x_inc);
   header.ny = 1 + (int) NINT((header.y_max - header.y_min) / header.y_inc);

      /* compute data bounds */
   
   x_left = header.x_min - xoffset;   
   x_right = header.x_min + (header.nx-1) * header.x_inc + xoffset;
   y_bottom = header.y_min - yoffset;
   y_top = header.y_min + (header.ny-1) * header.y_inc + yoffset;  
  
 /* for pixel registration */
 
   if ( header.node_offset == 1) {        
     --header.nx;
     --header.ny;
     x_left = header.x_min;    /* data boundaries */ 
     x_right = header.x_min  + header.nx * header.x_inc;
     y_bottom = header.y_min;
     y_top = header.y_min + header.ny * header.y_inc;  
  }


  fprintf (stderr, "Grid dimensions are nx = %d, ny = %d\n", header.nx, header.ny);
  if (header.node_offset)
     fprintf (stderr, "\n using pixel registration. ");
  else
     fprintf (stderr, "\n using gridnode registration. ");

/* allocate space for z arrays */

   nsq = header.nx * header.ny;
   
   zout = (double *) get_memory((void *)NULL, nsq, sizeof(double));
   zcount = (int *) get_memory((void *)NULL, nsq, sizeof(int));
   
  nread = 0;  /* counts number of points read in */

 /* loop for each input file :  attach each point read in to its
    nearest grid point.  Sum values and find average in case multiple points fall on a grid node*/
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
           fprintf(stderr, "Mismatch between actual (%d) and expected (3) fields near line %d\n", n_fields, nread);
           exit(1);
        }
	
        skip = FALSE;
	xflt = in[ix];
	yflt = in[iy];
        if (xflt < x_left || xflt > x_right) skip = TRUE;
        if (yflt < y_bottom || yflt > y_top) skip = TRUE;
        if (ABS(in[2]) > masked) skip = TRUE;   /* weed out empty/mask flags */
    
        if (!skip) {
	   error = xy2ij(&header, xflt, yflt, &col, &row);
	   sq = row * header.nx + col;
	   zout[sq] += in[2];
	   ++zcount[sq];
           ++nread;
         }  /* End if skip. */
	
    }  /* End while  */
      
NEXTFILE:
    if (nfiles) fclose (infile);

  } while (++curfile < nfiles);      /* End input phase */

/*  Find average values for each gridnode that has data attached to it */
   for (sq = 0; sq < nsq; ++sq) {
      if (zcount[sq] > 0)
          zout[sq] /= zcount[sq];
      else
          zout[sq] = empty_val;
   }  
 
  if (use_mask) {
     mask = (char *) get_memory((void *)NULL, (size_t)(nsq), sizeof(char *));
     fprintf(stderr,"\nReading mask file....");
     n_mask = get_mask(maskfile, &header, mask, yfirst);
     for (i = 0; i < nsq; i++)  {
       if (mask[i] )
          zout[i] = mask_val;
     }
     free((void *)mask);
  }
  
/* At this point, every element of zout is either a value, empty or masked */

  ywork = (double *) calloc(header.ny, sizeof(double));
  zwork = (double *) calloc(header.ny, sizeof(double));
  n_set = 0;

/* explicitly generate y-vector corresponding to grid*/  
  for (row = 0; row < header.ny; ++ row) {
     ij2xy(&header, 0, row, &xflt, &yflt);
     ywork[row] = (double) yflt;
  }
/* interpolate vertically column by column... */
  
  for (col = 0; col < header.nx; ++col) {
      for (row = 0; row < header.ny; ++row) {
         sq = row * header.nx + col;
	 zwork[row] = zout[sq];
      }
      
      for (row= 0; row < header.ny; ++row) {
          if (zwork[row] <= flagged ) {
	     /* search backward and forward in zwork array to find neighbors within yradius;
	     Suspend search if a masked neighbor is encountered.*/
	     found = 0;
	     rad = 1;
	     i = row -1;
	     while (i > 0 && rad <= yradius && !found) {

	         if (zwork[i] >= masked) {
		     rad = yradius +1;   /* suspend search */
		 }
		  else if (zwork[i] <= flagged) {
		      --i;
		      ++rad;
		  }
		  else  {
		       found = 1;
		  }
	     } /* end while */
	     
	     if (found) {
	         z1 = zwork[i];
	         y1 = ywork[i];
	     
	        found = 0;
	        rad = 1;
	        i = row + 1;
	        while (i < header.ny && rad <= yradius && !found) {
	            if (zwork[i] >= masked) {
		        rad = yradius +1;   /* suspend search */
		     }
		     else if (zwork[i] <= flagged) {
		         ++i;
			 ++rad;
		     }
		     else  {
		       found = 1;
	               z2 = zwork[i];
	               y2 = ywork[i];
		       zwork[row] = z1 + (z2 - z1) * (ywork[row] - y1) / (y2 - y1);
		       ++n_set;
		    }
	        }  /* end while */
	     }  /* end if found */
	  } /*end if zwork[row] */
      }/* end for row */
      
      /* load zwork back into zout */
      
       for (row = 0; row < header.ny; ++row) {
         sq = row * header.nx + col;
	 zout[sq] = zwork[row];
      }
     
  }  /* end for col */

  free(zwork);
  
/* interpolate horizontally row by row... */
  xwork = (double *) calloc(header.nx, sizeof(double));
  zwork = (double *) calloc(header.nx, sizeof(double));
  
/* explicitly generate x-vector corresponding to grid*/  
  for (col = 0; col < header.nx; ++ col) {
     ij2xy(&header, col, 0, &xflt, &yflt);
     xwork[col] = (double) xflt;
  }
  
   for (row = 0; row < header.ny; ++row) {
       for (col = 0; col < header.nx; ++col) {
         sq = row * header.nx + col;
	 zwork[col] = zout[sq];
       } 
       
       for (col = 0; col < header.nx; ++col) {
         if (zwork[col] <= flagged) {  /* search backward and forward in array to find neighbors within xradius */
	    found = 0;
	    rad = 1;
	    i = col-1;
	    while (i > 0 && rad <= xradius && !found) {
	         if (zwork[i] >= masked)  {
		    rad = xradius + 1;   /* suspend search */
		 }
		 else if (zwork[i] <= flagged) {
		    --i;
		    ++rad;
		 }
		 else {
		    found = 1;
		 }
	    } /* end while */
	    
	    if  (found) {
	       z1 = zwork[i];
	       x1 = xwork[i];
	       
	      found = 0;
	      rad = 1;
	      i = col + 1;
	      while (i < header.nx && rad <= xradius && !found ) {
	         if (zwork[i] >= masked)  {
		     rad = xradius + 1;  /* suspend search */
		  }
		 else if (zwork[i] <= flagged) {
		     ++i;
		     ++rad;
		 }
		 else  {
		    found = 1;
		    z2 = zwork[i];
		    x2 = xwork[i];
		    zwork[col] = z1 + (z2 - z1) * (xwork[col] - x1) / (x2 - x1);
		    ++n_set;
		 }
	      }  /* end while */ 
	    }  /* end if found */
	 }  /* end if zwork[col] */
       }/* end for col */
       
             /* load zwork back into zout */
      
       for (col = 0; col < header.nx; ++col) {
         sq = row * header.nx + col;
	 zout[sq] = zwork[col];
      }

   } /*end for row */ 
  
  free(zwork);

   /* output xyz triplets */
   
  n_empty = 0; 
   for (row = 0; row < header.ny; ++row) {
      for (col = 0; col < header.nx; ++col) {
         sq = row * header.nx + col;
         if ( zout[sq] > flagged && zout[sq] < masked) {
	    fprintf(outfile, "%8.3lf %8.3lf %.9g \n", xwork[col], ywork[row], zout[sq]);
	 }
	 
	 if (new_maskfile && zout[sq] >= masked)
	    fprintf(new_maskfile, "%8.3lf %8.3lf \n", xwork[col], ywork[row]);
      
         if (zout[sq] <= flagged )
	    ++n_empty;
      } /* end for col */
   } /* end for row */
  
  fprintf(stderr,"\nNumber of gridnodes: <%d>", nsq);
  fprintf(stderr,"\nNumber read in: <%d>  filled: <%d>   empty: <%d>   masked: <%d>\n", nread, n_set, n_empty, n_mask );
 
  fprintf(stderr,"\nEnd of %s.\n", argv[0]);
  
  exit(0);
  
}


/******************************************************************/
void print_usage(char *program)
{
  
    fprintf (stderr, "%s - Interpolates vertical section profiles onto a regular grid.\n\n", program);
    fprintf(stderr, "USAGE: %s [xyzfile(s)] -I<dx>[/<dy>]", program);
    fprintf(stderr, " -B<xmin/xmax/ymin/ymax> [-O<output_file>] [-P]");
    fprintf(stderr, "  [-M<mask_file> ] [-N<new_mask_file> ] [-S<xradius>/<yradius> ");
    fprintf(stderr, "  [-:] [-h]\n\n");
    fprintf(stderr, "-B   sets the grid bounds in user units: xmin/xmax/ymin/ymax.\n");
    fprintf(stderr, "-I   sets the grid spacing for the x/y directions.\n");
    fprintf(stderr, "\n\tOPTIONS:\n");
    fprintf(stderr, "-P  force pixel registration (default is gridnode registration).\n");
    fprintf(stderr, "-M name of file containing masking info.\n");
    fprintf(stderr, "      Each line of file may specify individual points to mask\n");
    fprintf(stderr, "      or multiple polygons separated by a '>' character. \n");
    fprintf(stderr, "      Specify 'I' or 'O' immediately after the '>' character\n");
    fprintf(stderr, "      to mask INSIDE or OUTSIDE of polygon.\n");
    fprintf(stderr, "      No '>' or a 'C' immediately after the '>'  \n");
    fprintf(stderr, "      signifies cell mode masking.  In this mode}\n");
    fprintf(stderr, "      an {x,y} pair OR {x,y,mask} triplet is given.\n");
    fprintf(stderr, "      The <mask> value can be 1 to mask the point\n");
    fprintf(stderr, "      or 0 (zero) to unmask the point.\n");
    fprintf(stderr, "-N name of new file for x/y pairs of empty or masked gridnodes after the interpolation.\n");
    fprintf(stderr, "-O  name of output file  (default is stdout).\n");
    fprintf(stderr, "-S set a search radius in integer grid increments\n");
    fprintf(stderr, "      for x and y directions.  If no data are \n");
    fprintf(stderr, "      within range of a node it is set to empty.\n");
    fprintf(stderr, "      Default is [-S%1d/%1d]\n", xradius, yradius);
    fprintf(stderr, "-: input data are ordered y x z  \n");
    fprintf(stderr, "       [default order is x y z]   \n");
    fprintf(stderr, "-h help....prints this message.  \n");
    return;
}  /* end print_usage() */

