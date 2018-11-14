/*  hb_findblanks.c  
................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
                             original: 1993
			     updated to ANSI-C standards Feb 2000
................................................................................

   Reads a file of gridded xyz values, and identifies the gridpoints 
  which contain no data and writes lon/lat out to the stdout device.
................................................................................
*/

#include <stdio.h>
#include <stdlib.h>   /*sfg added to deal with exit and malloc warnings*/
#include <math.h>
#include <string.h>



     /* boundaries for grid */

float   xmin, xmax, ymin, ymax, delta_x, delta_y;     
int     ncols, nrows;

int main(int argc, char **argv)
{
   FILE  *outfile;
   short bopt,  iopt;
   char *st, *name[2];
   int   error, n, i, row, col;
   float lat, lon;
   float **x, x1; 
   void print_usage(char *);
   int readprop(char *, float **);

/* check for command line arguments */

   if (argc < 2) {
     print_usage(argv[0]);
     exit(1);
   }

/* set these default values... */

   bopt  = iopt  = 0;
   error = 0;
   outfile = stdout;
   n = 0;


/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'B':                    /* get grid bounds */
                        bopt = 1;
                        st = &argv[i][2];
                           if (*st == '/')
                               ++st;
                        error = (sscanf(st,"%f", &xmin) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%f", &xmax) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%f", &ymin) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%f", &ymax) != 1);
                        break;

               case 'I':
                        iopt = 1;
                        error = (sscanf(&argv[i][2],"%f", &delta_x) == 1) ? 0 : 1;
                        delta_y = delta_x;
                        st = strchr(&argv[i][2],'/');
                        if (st != NULL) {
                          sscanf(++st,"%f", &delta_y);
                        }
                        break;

      	       case 'h':  
	           print_usage(argv[0]);
	           exit(0);


               default:
                        error = 1;

          }    /* end switch */

          if (error ) {
             fprintf(stderr,"\nError parsing command line args.\n");
             fprintf(stderr,"     in particular: '%s'\n", argv[i]);
             exit(1);
          }

       }  /* end if */
       else {
          if (n < 1) 
             name[n++] = argv[i];
          
          else {
             fprintf(stderr,"\nToo many input files specified!\n");
             fprintf(stderr,"Ignoring %s\n", argv[i]);
          }
       }

   }  /* end for */

   if (!bopt || !iopt || (n < 1) ) {
       fprintf(stderr,"\nYou must specify input file, bounds, and gridspacing!\n");
       exit(1);
   }
   /* compute dimensions of matrix formed by grid  */

   nrows = (int) (ceil((double)((ymax - ymin) / delta_y)) + .0001);
   ncols = (int) (ceil((double)((xmax - xmin) / delta_x)) + .0001);
      

/*   allocate space for gridded values and initialize... */

   x = (float **) malloc(nrows * sizeof(float *));
   for (i = 0; i < nrows; ++i ) {
      x[i] = (float *) malloc(ncols * sizeof(float));
   }
   for (row = 0; row < nrows; ++row) {
      for (col = 0; col < ncols; ++col) {
          x[row][col] = -99999.0;
      }
   }

/*   read in values from  file ... */
 
   n = readprop(name[0], x);

/* now  output lon/lat pairs where there is no data ... */


   for (row = 0; row < nrows; ++row) {
       for (col = 0; col < ncols; ++col) {
        if (x[row][col] < -9998.) {
           lon = xmin +  (col + .5) * delta_x;
           lat = ymin + (row + .5) * delta_y;
           fprintf(outfile,"%8.3f %8.3f \n", lon, lat);
        }
      }
   }
   
   fflush(outfile);
   fprintf(stderr, "hb_findblanks done.\n");
   exit(0);
}  /* end main */

/****************************************************************************/

void print_usage(char *program)
{
 fprintf(stderr,"***********************************************************"); 
 fprintf(stderr,"\nReads a file of gridded xyz values, identifies gridpoints"); 
 fprintf(stderr,"\nwhich contain no data, and writes lon/lat of those points to stdout. ");
 fprintf(stderr,"\n***********************************************************"); 
  fprintf(stderr,"\nUsage:  %s xyz_file(s) -Bwest/east/south/north -Ideltax/deltay  > output_file", program);
   fprintf(stderr,"\n   -B : specifies grid bounds. Ex: -B-80/0/0/65");
   fprintf(stderr,"\n   -I : specifies grid increments. Ex: -I1/.5");
  fprintf(stderr,"\n\n");  
  return;
}
/***********************************************************************/
int readprop(char *filename, float **x)

{
   FILE *infile;
   float  y;
   float lon, lat;
   int row, col;


   if ((infile = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "\n Unable to open %s for input\n\n", filename);
        exit(1);
   }
   fprintf(stderr," Opened %s ...\n", filename);

   while (fscanf(infile,"%f%f%f", &lon, &lat, &y) != EOF) {
     row = (int) (.0001 + (lat - ymin) / delta_y);
     col = (int) (.0001 + (lon - xmin) / delta_x);


     /* check that it is within the bounds */

     if ((row >= 0) && (row < nrows) && (col >= 0) && (col < ncols)) 
            x[row][col] = y;
   }

   close (infile);

   return (0);

}  /* end readprop() */




