/* hb_read_topo.c
................................................................................
                         *******  HydroBase2 *******
................................................................................
Based on the Hydrobase 2 framework written by Ruth Curry.
Framework modified by Stefan Gary
................................................................................
................................................................................
.  Reads in a file of global topography and spits output in xyz format
.  to stdout for error checking.
.
................................................................................
................................................................................

*/
#include <stdio.h>
#include <stdlib.h>   /*sfg added to deal with exit and malloc warnings*/
#include <ctype.h>
#include <string.h>
#include "hydrobase.h"
#include "hb_paths.h"

#define    PRINT_MSG 1


     /* boundaries for grid */

double   xmin, xmax, ymin, ymax, delta_x, delta_y;
float  *topo_lat, *topo_lon;
int 	lon0to360, depth_is_neg;   
int     ncols, nrows;


  /* prototypes for functions declared locally */
  
void print_usage(char *);
short **get_topo(FILE *);

main (int argc, char **argv)
{
  int     i, j;
  int     i_flag, b_flag, t_flag;
  int     error = 0;
  short **ztopo, z;
  char   *st;
  char    fname[200];
  FILE   *topofile;
  
/*  set these default values */

    i_flag = b_flag = t_flag = 0;
    lon0to360 = 1;
    depth_is_neg = 1;
    xmin = 0;
    xmax = 360;
    ymin = -90.0;
    ymax = 90.0;
    delta_x = delta_y = 0.1;

/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'B':                    /* get grid bounds */
                        b_flag = 1;
                        st = &argv[i][2];
                           if (*st == '/')
                               ++st;
                        error = (sscanf(st,"%lf", &xmin) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%lf", &xmax) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%lf", &ymin) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%lf", &ymax) != 1);

                        if (xmin > xmax)  
                          fprintf(stderr,"\nW bound cannot exceed E bound.\n");
                                               
                        if (ymin > ymax)  
                          fprintf(stderr,"\nS bound cannot exceed N bound.\n");
                          
                        if (xmin < 0)
                           lon0to360 = 0;
                           
                        break;

               case 'I' :
                        i_flag = 1;
                        error = (sscanf(&argv[i][2],"%lf", &delta_x) == 1) ? 0 : 1;
                        delta_y = delta_x;
                        st = strchr(&argv[i][2],'/');
                        if (st != NULL) {
                          sscanf(++st,"%lf", &delta_y);
                        }
                        break;

               case 'N' :
                        depth_is_neg = 0;
                        break;

               case 'T':
                        topofile = fopen(&argv[i][2],"r");
                        if (topofile == NULL) {
                          fprintf(stderr,"\nError opening %s for input\n", &argv[i][2]);
                          exit(1);
                        }
                        t_flag = 1;
                        break;

               case 'h':
                        print_usage(argv[0]);
                        exit(0);
			break;

               default  :
                        fprintf(stderr,"\nError parsing command line");
                        fprintf(stderr,"\n in particular: %s\n", argv[i]);
                        exit(1);
            }  /* end switch */
       }  /* end if */
       else  {
	 /* There should be nothing w/o prepended "-" on command line.*/
	 fprintf(stderr,"\nError parsing command line");
	 fprintf(stderr,"\n in particular: %s\n", argv[i]);
	 exit(1);
       }
   }  /* end for */

   if (! t_flag ) {
   
      topofile = fopen(BATHPATH,"r");
      if (topofile == NULL) {
          fprintf(stderr,"\nUnable to open default topography file\n", BATHPATH);
          fprintf(stderr,"You can specify the name of a binary topography file for input with [-T]\n"); 
          exit(1);
      }
    
   } 
   if (! i_flag) {
      fprintf(stderr,"\nIncrement for topofile: [%.2lf/%.2lf]", delta_x, delta_y);
   } 
   
/* read in topography file ... */

   ztopo = get_topo(topofile);
   fclose(topofile);

   for ( j = 0; j < nrows; j++ ) {
     for ( i = 0; i < ncols; i++ ) {
       fprintf(stderr,"\rAt node location [%5d,%5d] of [%5d,%5d].",j,i,nrows,ncols);
       fprintf(stdout,"\n%6f %6f %6f",topo_lon[i],topo_lat[j],(float)ztopo[j][i]);
     }
   }

   fprintf(stderr,"\n\nEnd of %s.\n\n", argv[0]);
   exit(0);

} /* end main() */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\n%s Reads topography file and posts it to stderr.", program);
   fprintf(stderr,"\nThe default topography file is %s", BATHPATH);
   fprintf(stderr,"\nFor an alternate topography file, use -T to specify its name,");
   fprintf(stderr,"\n-B to specify its bounds, and -I for its grid increment");
   fprintf(stderr,"\nUse -N if the alternate topography has positive seafloor depths.");
   
   fprintf(stderr,"\n\nUsage:  %s", program);
   fprintf(stderr," [-T<topofile>] [-B<west/east/south/north>] [-I<xincr/yincr>] [-N] [-h]");
   fprintf(stderr,"\n");
   
   fprintf(stderr,"\tOPTIONS:");
   fprintf(stderr,"\n[-B] : specifies bounds of topography grid. Default [0/360/-90/90]");
   fprintf(stderr,"\n[-I] : specify x-increment [/y-increment] for topography grid");
   fprintf(stderr,"\n         Default [%.2lf/%.2lf]", delta_x, delta_y);
   fprintf(stderr,"\n[-T] : name of topography file. Default is [%s]", BATHPATH);
   fprintf(stderr,"\n[-N] : specify this flag if topography file has positive seafloor depths.");
   fprintf(stderr,"\n[-h] : help -- prints this message");
  fprintf(stderr,"\n\n");  
   return;
}

/****************************************************************************/

short **get_topo(FILE *fptr)
   /* allocates memory and reads in topography values.  Returns a pointer to
      the start of the memory block. An error causes an error message to be
      printed and an exit.   */
{

   int row, n;
   short **z;
   
   /* Set some globally defined variables */
   
   /* Compute the number of rows and columns. */
   nrows = (int)NINT((ymax - ymin) / delta_y) + 1;
   ncols = (int)NINT((xmax - xmin) / delta_x) + 1;
   
   /* Allocate space for coordinate arrays
    * based on number of rows and columns. */
   topo_lat = (float *) malloc(nrows * sizeof(float));
   topo_lon = (float *) malloc(ncols * sizeof(float));
   
   /* Fill in the values for the coordiante
    * arrays for topography grid. */
   for (n = 0; n < nrows; ++n) {
     topo_lat[n] = (float) ((double) n * delta_y + ymin);
   }
   
   for (n = 0; n < ncols; ++n) {
     topo_lon[n] = (float) ((double) n * delta_x + xmin);
   }
   
   /* Allocate space for the topomap based
    * on the number of rows and columns. */
   z = (short **) malloc(nrows * sizeof(short *));
   if (z == NULL) {
      fprintf(stderr,"\nError allocating memory.\n");
      exit(1);
   }
   
   for (row = 0; row < nrows; ++row) {
     z[row] = (short *) malloc(ncols * sizeof(short));
     if (z[row] == NULL) {
         fprintf(stderr,"\nError allocating memory.\n");
         exit(1);
      }
   }
   
   fprintf(stderr,"\nReading in topography values ");
   
   for (row = 0; row < nrows; ++row) {
     n = fread((short *)z[row], sizeof(short), ncols, fptr);
     if (n != ncols) {
         fprintf(stderr,"\nError reading the topofile at row %d\n", row);
         exit(1);
     }
     if ((row % 10) == 0)   /* signal progress */
        fprintf(stderr,".");
   }
   
   fprintf(stderr,"\nFinished reading topofile.\n");
   
   return(z);

} /* end get_topo() */
/****************************************************************************/
