/* hb_msq10bounds.c
   
   Accepts a 10-deg MSq and (optional) increment as arguments and
   returns bounds in the proper format: w/e/s/n  .  Specifying
   an increment will enlarge the bounds by half the increment
   -- useful to generate node-grids as opposed to pixel-grids.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void print_usage(char *);

main(int argc, char **argv)
{

   int msq10, i, hem, rem, error, no_offset;
   float incr;
   float minlon, maxlon, minlat, maxlat;
   float lat, lon;
   char *s;
   
   
   if (argc < 2) { 
      print_usage(argv[0]); 
      exit(1); 
   }
  
  incr = 0.0;
  no_offset = 1;
  error = 0;
  
/* parse the  command line arguments */
   
    for (i = 1; i < argc; i++) { 
      if (argv[i][0] == '-') { 
         s = &argv[i][1]; 
         switch (*s) { 
            case 'I': 
               ++s;
               incr = (float) atof(s);
               incr *= .5;
               no_offset = 0; 
               break;
            case 'h': 
               print_usage(argv[0]); 
               exit(0);
	       
            default:  
               error = 1;
         } /* end switch */
             
      } /* end if */
      else {

         if (! (msq10 = atoi(argv[i]))) {
          fprintf(stderr,"\nUnable to parse MSQ10: %s \n", argv[i]);
          exit(1);
         }
      }

      if (error ) { 
        fprintf(stderr,"\nError parsing command line args.\n");
        fprintf(stderr," in particular:  '%s'\n", argv[i]); 
        exit(1); 
      }
   } /* end for */
 
   hem = msq10 / 1000; 
   rem = msq10 - (hem * 1000);
   lat = (float)(rem / 100 * 10);
   lon = (float) (rem % 100 * 10);

   if (no_offset) {   
      if (hem == 7) {
        fprintf(stdout,"%.3f/%.3f/%.3f/%.3f \n", -(lon +10), -lon, lat, lat+10);
        exit(0);
      }
      if (hem == 5) {
        fprintf(stdout,"%.3f/%.3f/%.3f/%.3f \n", -(lon +10), -lon, -(lat+10), -lat);
        exit(0);
      }
      if (hem == 3) {
        fprintf(stdout,"%.3f/%.3f/%.3f/%.3f \n", lon, lon+10, -(lat+10), -lat);
        exit(0);
      }
      if (hem == 1) {
        fprintf(stdout,"%.3f/%.3f/%.3f/%.3f \n", lon, lon+10, lat, lat+10);
        exit(0);
      }
   }
   
   /* to offset bounds so first pixel is centered on minlat/minlon...*/
   
      if (hem == 7) {
        minlon = -lon - 10 - incr;
        maxlon = -lon + incr;
        minlat = lat - incr;
        maxlat = lat + 10 + incr;
        
        fprintf(stdout,"%.3f/%.3f/%.3f/%.3f \n", minlon,maxlon, minlat, maxlat);
        exit(0);
      }
      if (hem == 5) {
        minlon = -lon - 10 - incr;
        maxlon = -lon + incr;
        minlat = -lat -10 - incr;
        maxlat = -lat + incr;
         fprintf(stdout,"%.3f/%.3f/%.3f/%.3f \n", minlon,maxlon, minlat, maxlat);
        exit(0);
      }
      if (hem == 3) {
        minlon = lon - incr;
        maxlon = lon + 10 + incr;
        minlat = -lat -10 - incr;
        maxlat = -lat + incr;
         fprintf(stdout,"%.3f/%.3f/%.3f/%.3f \n", minlon,maxlon, minlat, maxlat);
        exit(0);
      }
      if (hem == 1) {
        minlon = (lon -incr);
        maxlon = lon + 10 + incr;
        minlat = lat - incr;
        maxlat = lat + 10 + incr;
       
        fprintf(stdout,"%.3f/%.3f/%.3f/%.3f \n", minlon,maxlon, minlat, maxlat);
        exit(0);
      }
    
   fprintf(stderr, "\nCould not recognize the hemisphere in MSQ #%d\n", msq10);
   exit(2);  

} /* end main */
/****************************************************************************/

void print_usage(char *program) 
{ 
   fprintf(stderr,"\n Accepts a 10-deg MSq and (optional) increment as arguments and");
   fprintf(stderr,"\nreturns bounds in the proper format: w/e/s/n  .  Specifying");
   fprintf(stderr,"\nan increment will enlarge the bounds by half the increment");
   fprintf(stderr,"\n-- useful to generate node-grids as opposed to pixel-grids.");

   fprintf(stderr,"\n\nUsage:  %s msq_10 [-I<incr> ", program);
   fprintf(stderr,"\n\n msq_10 is the 4-digit WMO (Marsden Square)."); 
   fprintf(stderr,"\n \n  OPTIONS:"); 
   fprintf(stderr,"\n[-I]  : specifies grid increment."); 
   fprintf(stderr,"\n        Bounds will be offset by half the increment."); 
   fprintf(stderr,"\n[-h]  : help -- prints this message."); 
  fprintf(stderr,"\n\n");
   return;
}
   
   
