/* hb_smooth1d.c
................................................................................
                          *******  HydroBase 2 *******
                             Ruth Curry
                             Woods Hole Oceanographic Institution
                             May 2001
..................................................................  
.  Reads a file of x,y pairs and applies the specified filter
.  for the specified window.  It will interpolate for y within x
.  domain, but will not extrapolate.  x,y pairs are output.  
..................................................................  
*/ 
#include<stdio.h> 
#include <math.h>
#include "hb_filters.h"

/* prototypes for locally defined functions */	
void print_usage(char *);
           
int main(int argc, char **argv)
{
   short iopt, yopt, wopt; 
   char *s,*fname; 
   char stuff[100];
   FILE *infile, *outfile; 
   int yr, line, i, j, *n, nn, nnn, error; 
   int minyr, maxyr,nyrs;
   int window, half, end; 
   float y;
   double *x, xx; 
   void print_usage();

   if (argc < 2) { 
      print_usage(argv[0]); 
      exit(1); 
   }

   infile = NULL; 
   outfile = stdout; 
   iopt = yopt = wopt = 0; 
   line = 0;

/* parse the command line arguments */

   for (i = 1; i < argc; i++) { 
      error = (argv[i][0] == '-') ?  0 :  1; 
      if (!error) { 
         s = &argv[i][1]; 
         switch (*s) { 
            case 'Y':  
               ++s; 
               yopt = 1; 
               if (*s == '/')
                  ++s; 
               error += (sscanf(s,"%d", &minyr) != 1); 
               s = strchr(s,'/'); /* check for another delimiter*/ 
               error = (s == NULL); 
               if (s != NULL) { 
                  ++s; /* move past delimiter */ 
                  error += (sscanf(s,"%d", &maxyr) != 1); 
               } 
               if (error) 
                  break; 
               if (maxyr < minyr) { 
                  fprintf(stderr,"\n minyr exceeds maxyr\n"); 
                  exit(1); 
               } 
               nyrs = maxyr - minyr + 1; 
               break; 
               
            case 'I':  
               iopt = 1; 
               ++s; 
               if ((infile = fopen(s, "r")) == NULL) { 
                  error = 1; 
                  fprintf(stderr, "\nError opening %s\n", s); 
               }
               break; 
            case 'O':  
               ++s; 
               if ((outfile = fopen(s, "w")) == NULL) { 
                  error = 1;
                  fprintf(stderr, "\nError opening %s\n", s); 
               } 
               break; 
            case 'W':  
               ++s; 
               wopt = 1; 
               if (*s == '/')
                  ++s; 
               error += (sscanf(s,"%d", &window) != 1); 
               if ((window % 2) == 0)
                  ++window;
               break; 
            default:  
               error = 1;
         } /* end switch */
             
      } /* end if */


      if (error ) { 
        fprintf(stderr,"\nError parsing command line args.\n");
        fprintf(stderr," in particular:  '%s'\n", argv[i]); 
        exit(1); 
     }


   } /* end for */

   if (!  (iopt && yopt && wopt) ) { 
      print_usage(argv[0]); 
      fprintf(stderr,"\nYou must specify -Y, -W and -I options!!!\n"); 
      exit(1); 
   }

  /* allocate memory and initialize arrays...  */

   x = (double *) calloc(nyrs, sizeof(double)); 
   n = (int *) calloc(nyrs, sizeof(int));


/* loop to read input file */

   while ((i = fscanf(infile,"%d", &yr)) != EOF) { 
      ++line; 
      if (i != 1) {
         fprintf(stderr,"\n Error reading infile at line #%d\n", line); 
         exit(1); 
      } 
      yr -= minyr; /* yr now contains index to arrays */ 
      i = fscanf(infile,"%f%[^\n]", &y, stuff);
      if (yr >=0 && yr < nyrs) {   /* check that this point is in range */
         x[yr] += y; 
         ++n[yr]; 
      }
   }

/* first make sure each year is an average...*/

   for (i = 0; i < nyrs; ++i) { 
      if (n[i] > 0) { 
         x[yr] /= n[yr];
         n[yr] = 1;
      }
   }
   
/* now apply the filter...*/

   half = window >> 1;   /* divide by 2 */
   end = nyrs - half;
   
   for (i = 0; i < half;  ++i) {
      if (n[i])
       fprintf(outfile,"%5d %10.4lf\n", (i+minyr), x[i]);
   }
   
   for (i = half; i < end; ++i) {
      xx = 0;
      nn = 0;
      for (j = i-half; j <= i+half; ++j) {
        if (j == i){
           if (!nn) {
              nn = n[i];  /* no obs in 1st half -- don't average*/
              xx = x[i];
              break;
           }
           nnn = nn + n[i];  /* save nobs in first half */
        }
        if (n[j]) {
          xx += x[j];
          ++nn;
        }
        if (j == (i+half)) {
           if ((nn - nnn) == 0) {  /* no obs in 2nd half -- don't average*/
              nn = n[i];
              xx = x[i];
           }
        }
      }
      if (nn > 0)
       fprintf(outfile,"%5d %10.4lf\n", (i+minyr), (xx / (double) nn));
   }
   
    for (i = end; i < nyrs;  ++i) {
      if (n[i])
       fprintf(outfile,"%5d %10.4lf\n", (i+minyr), x[i]);
   }
  

   fclose(infile); 
   fclose(outfile); 
   fprintf(stderr,"\n End of %s.\n", argv[0]);
   
   exit(0);


} /* end main */


/****************************************************************************/

void print_usage(program) 
char *program; 
{ 
   fprintf(stderr,"\nUsage:  %s -I<input_filename> -O<output_file_root> -Y<minyear/maxyear> -W<window>", program);
   fprintf(stderr,"\n\n  -I  :  specifies input file_name"); 
   fprintf(stderr,"\n [-O] : specifies optional output file");
   fprintf(stderr,"\n        default is stdout ");
   fprintf(stderr,"\n  -Y  : specifies min/max years "); 
   fprintf(stderr,"\n  -W  : specifies window in # of pts "); 
   fprintf(stderr,"\n\n");
   return;
}
   

