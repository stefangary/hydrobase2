/* hb_stats1d.c
................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                    author:  Ruth Curry
                             Woods Hole Oceanographic Institution
                             June 2000
................................................................................
..................................................................  
.  reads a  file of x,y pairs and determines the mean 
.  and std deviation of values y for specified bins of x. 
..................................................................  
*/ 
#include <stdio.h> 
#include <stdlib.h>   /*sfg added to deal with exit and malloc warnings*/
#include <math.h>
#include <string.h>

#define    ABS(x)       (((x) < 0) ? -(x) : (x))

/*  prototypes for locally defined functions */
void print_usage(char *);
           
int main (int argc, char **argv)
{ 
   short ropt; 
   char *s,*fname; 
   int   curfile = 1, nfiles = 0; 
   FILE *infile, *outfile; 
   int  line, i, ix, *n, error, nxbins; 
   float xmin, xmax, xdelta; 
   float y, x;
   double *sum, *sumsq; 
   float mean, stddev, var;

   if (argc < 1) { 
      print_usage(argv[0]); 
      exit(1); 
   }

   infile = stdin; 
   outfile = stdout; 
   ropt = 0; 
   line = 0;

/* parse the command line arguments */

   for (i = 1; i < argc; i++) { 
      if (argv[i][0] == '-') {
         s = &argv[i][1]; 
         switch (*s) { 
            case 'R':  
               ++s; 
               ropt = 1; 
               if (*s == '/')
                  ++s; 
               error += (sscanf(s,"%f", &xmin) == 1); 
               s = strchr(s,'/'); /* check for another delimiter*/ 
               error = (s == NULL); 
               if (s != NULL) { 
                  ++s; /* move past delimiter */ 
                  error += (sscanf(s,"%f", &xmax) != 1); 
               } 
               s = strchr(s,'/'); /* check for another delimiter*/ 
               if (s != NULL) { 
                  ++s; /* move past delimiter */ 
                  error += (sscanf(s,"%f", &xdelta) != 1); 
               } 
               if (error) 
                  break; 
               if (xmax < xmin) { 
                  fprintf(stderr,"\n xmin exceeds xmax\n"); 
                  exit(1); 
               } 
               nxbins = (xmax - xmin)/ xdelta + 1; 
               break; 
               
            case 'O':  
               ++s; 
               if ((outfile = fopen(s, "w")) == NULL) { 
                  error = 1;
                  fprintf(stderr, "\nError opening %s for output\n", s); 
               } 
               break; 
               
            case 'h': 
               print_usage(argv[0]);
               exit(0);
               
            default:  
               error = 1;
               
         } /* end switch */
             

         if (error ) { 
           fprintf(stderr,"\nError parsing command line args.\n");
           fprintf(stderr," in particular:  '%s'\n", argv[i]); 
           exit(1); 
         }
      }
      else  {
        ++nfiles;
      }
      
   } /* end for */

   if (! ropt ) { 
      fprintf(stderr,"\nYou must specify -R arguments to define xmin/xmax"); 
      fprintf(stderr,"\nThe default delta-x is 1.0, but a different value");
      fprintf(stderr,"\ncan be specified.  ex:  -R1900/1999/4\n"); 
      exit(1); 
   }

  /* allocate memory and initialize arrays...  */

   sum = (double *) malloc(nxbins * sizeof(double)); 
   sumsq = (double *) malloc(nxbins * sizeof(double)); 
   n = (int *) malloc(nxbins * sizeof(int));

   for (i = 0; i < nxbins; ++i) { 
      sum[i] = 0.0; 
      sumsq[i] = 0.0; 
      n[i] = 0; 
   }

 /* loop for each input file */
 
   do {
   
   
     if (nfiles > 0) {
       infile = fopen(argv[curfile],"r");
       if (infile == NULL) {
          fprintf(stderr,"\nUnable to open %s for input.\n", argv[curfile]);
          goto NEXTFILE;
       }
     }

     while ((i = fscanf(infile,"%f%f", &x, &y)) != EOF) { 
        ++line; 
        if (i != 2) {
           fprintf(stderr,"\n Error reading infile at line #%d\n", line); 
           exit(1); 
        }
        ix = (int) ((x - xmin) / xdelta); /* ix now contains index to arrays */
        if (ix >=0 && ix < nxbins) {   /* check that this point is in range */
           sum[ix] += y; 
           sumsq[ix] += (y * y); 
           ++n[ix]; 
        }
     }

/* compute means and std devs and output results...*/

     for (i = 0; i < nxbins; ++i) { 
        if (n[i] > 0) { 
           mean = (float)sum[i]; 
           stddev = 0.0;
           if (n[i] > 1) { 
              mean = (float)(sum[i] / (double) n[i]); 
              var = (sumsq[i] - sum[i] * sum[i] / (double) n[i]) / (double) (n[i]-1); 
              stddev = (float) sqrt(ABS(var)); 
           } 
           fprintf(outfile,"%f %10.4f %10.4f %5d\n", (i * xdelta + xmin), mean, stddev, n[i]);
        }
     }
     
NEXTFILE:
     if (nfiles)
         fclose(infile);
      
   } while (++curfile < nfiles);

   fclose(outfile); 
   fprintf(stderr,"\n End of %s.\n", argv[0]);
   
   exit(0);


} /* end main */


/****************************************************************************/

void print_usage(char *program) 
{ 
   fprintf(stderr,"\nReads x-y pairs and computes mean and standard deviation", program);
   fprintf(stderr,"\nof y for the x binned according to min/max/increment ");
   fprintf(stderr,"\nspecified with the -R argument. Input and output files can be");
   fprintf(stderr,"\nspecified -- or the stdin and stdout devices will be used");
   fprintf(stderr,"\n\nUsage:  %s input_file_list [-O<output_file>]  -R<min/max/incr> ", program);
   fprintf(stderr,"\n [-O] : specifies optional output file");
   fprintf(stderr,"\n        default is stdout ");
   fprintf(stderr,"\n  -R  : specify minimum/maximum/[delta] to create x bins. "); 
   fprintf(stderr,"\n        ex: -R1920/1995/5    default delta is 1 "); 
   fprintf(stderr,"\n\n");
   return;
}
   

