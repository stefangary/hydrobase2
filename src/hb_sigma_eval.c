/*  hb_sigma_eval.c

................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             1995
                             Updated Mar 2000 to ANSI
................................................................................
____________________________________________________________________________
  USAGE:  

 hb_sigma_eval filename(_roots)   [-D<dirname>] [-E<file_extent>]

  list of filenames MUST be first argument!

 -D : specifies directory for input files (default is current dir)  
      ex: -D/home/pete/data

 -E : specifies input_file extent (default is no extent) ex: -E.sum

____________________________________________________________________________
sigma_eval determines the min, max,and mean potential density
for a list of standard levels.                                                    ____________________________________________________________________________
*/


#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include "hydrobase.h"

#define   MISSINGVAL    -99999. 


/* input file pathnames */

#define    EXTENT   ""
#define    DIR      ""


   /* global variables to store station ... */

struct HYDRO_DATA station;
struct HYDRO_HDR hdr;

#define NLEV0 31
double lev0[NLEV0] = 
{0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
120, 140, 160, 180, 200,
220, 240, 260, 280, 300,
320, 340, 360, 380, 400,
420, 440, 460, 480, 500};

#define NLEV1 51
double lev1[NLEV1] = 
{500, 520, 540, 560, 580,
600, 620, 640, 660, 680,
700, 720, 740, 760, 780,
800, 820, 840, 860, 880,
900, 920, 940, 960, 980,
1000, 1020, 1040, 1060, 1080,
1100, 1120, 1140, 1160, 1180,
1200, 1220, 1240, 1260, 1280,
1300, 1320, 1340, 1360, 1380,
1400, 1420, 1440, 1460, 1480,
1500};

#define NLEV2 21
double lev2[NLEV2] = 
{1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950,
 2000, 2050, 2100, 2150, 2200, 2250, 2300, 2350, 2400, 2450,
 2500};
 
#define NLEV3 21
double lev3[NLEV3] = 
{2500, 2550, 2600, 2650, 2700, 2750, 2800, 2850, 2900, 2950,
 3000, 3050, 3100, 3150, 3200, 3250, 3300, 3350, 3400, 3450,
 3500};

#define NLEV4 46
double lev4[NLEV4] = 
{3500, 3550, 3600, 3650, 3700, 3750, 3800, 3850, 3900, 3950,
 4000, 4050, 4100, 4150, 4200, 4250, 4300, 4350, 4400, 4450,
 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400,
 5500, 5600, 5700, 5800, 5900, 6000, 6100, 6200, 6300, 6400,
 6500, 6600, 6700, 6800, 6900, 7000};
 
  /* global vars for stats at each standard level */
double min0[NLEV0], max0[NLEV0], sum0[NLEV0];
double min1[NLEV1], max1[NLEV1], sum1[NLEV1];
double min2[NLEV2], max2[NLEV2], sum2[NLEV2];
double min3[NLEV3], max3[NLEV3], sum3[NLEV3];
double min4[NLEV4], max4[NLEV4], sum4[NLEV4];
 int count0[NLEV0];
 int count1[NLEV1];
 int count2[NLEV2];
 int count3[NLEV3];
 int count4[NLEV4];

  /* prototypes for locally defined functions */
  
void    print_usage(char *);
void    get_hydro_data(int);


main (int argc, char **argv)
{
   int     index, nprops = 0;
   int     nstdlevs;
   int     i;
   int     curfile = 1, nfiles = 0; 
   int     print_msg = 1;       /* print message in file open routines*/
   int     error, prompt;
   char    *extent, *dir, *st;
   int     infile;
   double  last, mean;


/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
 
/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    error = 0;
    
/* initialize these ... */

   for (i = 0; i < MAXPROP; ++i) {
      station.observ[i] = (double *) NULL;
   }
   hdr.prop_id = (int *) NULL;

   for (i = 0; i < NLEV0; ++i) {
     min0[i] = 999;
     max0[i] = 0;
     sum0[i] = 0;
     count0[i] = 0;
   }
   
   for (i = 0; i < NLEV1; ++i) {
     min1[i] = 999;
     max1[i] = 0;
     sum1[i] = 0;
     count1[i] = 0;
   }
   
   for (i = 0; i < NLEV2; ++i) {
     min2[i] = 999;
     max2[i] = 0;
     sum2[i] = 0;
     count2[i] = 0;
   }

   for (i = 0; i < NLEV3; ++i) {
     min3[i] = 999;
     max3[i] = 0;
     sum3[i] = 0;
     count3[i] = 0;
   }

   for (i = 0; i < NLEV4; ++i) {
     min4[i] = 999;
     max4[i] = 0;
     sum4[i] = 0;
     count4[i] = 0;
   }

/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'D':                   /* get input dir */
                        dir = &argv[i][2];
                        break;
               case 'E':                    /* get file extent */
                        extent = &argv[i][2];
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

       else  {
           ++nfiles;
       }
   }  /* end for */



/* loop for each input file */

   do {

   if ( !nfiles) {
      infile = STDIN;
      fprintf(stderr,"\n Expecting data from stdin....  ");
   }
   else {
       infile = open_hydro_file(dir, argv[curfile], extent, print_msg);
       if (infile < 0)
          goto NEXTFILE;
   }

     
            /* read each file completely */

 
         get_hydro_data(infile);

NEXTFILE:

         close(infile);

   } while (curfile++ < nfiles );
   
/* Print out accumulated statistics for each standard level */
   
   fprintf(stdout, "\n%13s %6s %7s %7s %7s  %7s %7s\n", "StdLev (pref)", "N", "Min", "Max", "Range", "Mean", "dMean" ); 
   
   last = 0; 
   for (i = 0; i < NLEV0; ++i) {
      if (count0[i] > 0) {
         mean = sum0[i]/(double)count0[i];
         fprintf(stdout, "%7.1lf (0)     %6d %7.4lf %7.4lf %7.4lf %7.4lf %7.4lf\n", 
             lev0[i], count0[i], min0[i], max0[i], max0[i]-min0[i], mean, mean-last);
         last = mean;
       }
   }

   fprintf(stdout, " \n");  
   last = 0; 
   for (i = 0; i < NLEV1; ++i) {
      if (count1[i] > 0) {
         mean = sum1[i]/(double)count1[i];
         fprintf(stdout, "%7.1lf (1000)  %6d %7.4lf %7.4lf %7.4lf %7.4lf %7.4lf\n", 
             lev1[i], count1[i], min1[i], max1[i], max1[i]-min1[i], mean, mean-last);
         last = mean;
       }
   }
   
   fprintf(stdout, " \n");  
   last = 0; 
   for (i = 0; i < NLEV2; ++i) {
      if (count2[i] > 0) {
         mean = sum2[i]/(double)count2[i];
         fprintf(stdout, "%7.1lf (2000)  %6d %7.4lf %7.4lf %7.4lf %7.4lf %7.4lf\n", 
             lev2[i], count2[i], min2[i], max2[i], max2[i]-min2[i], mean, mean-last);
         last = mean;
       }
   }
   
   fprintf(stdout, " \n");  
   last = 0; 
   for (i = 0; i < NLEV3; ++i) {
      if (count3[i] > 0) {
         mean = sum3[i]/(double)count3[i];
         fprintf(stdout, "%7.1lf (3000)  %6d %7.4lf %7.4lf %7.4lf %7.4lf %7.4lf\n", 
             lev3[i], count3[i], min3[i], max3[i], max3[i]-min3[i], mean, mean-last);
         last = mean;
       }
   }
   
   fprintf(stdout, " \n");  
   last = 0; 
   for (i = 0; i < NLEV4; ++i) {
      if (count4[i] > 0) {
         mean = sum4[i]/(double)count4[i];
         fprintf(stdout, "%7.1lf (4000)  %6d %7.4lf %7.4lf %7.4lf %7.4lf %7.4lf\n", 
             lev4[i], count4[i], min4[i], max4[i], max4[i]-min4[i], mean, mean-last);
         last = mean;
       }
   }
   
   fprintf(stderr,"\n\nEnd of hb_sigma_eval.\n");
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{   

   fprintf(stderr,"\n%s determines the min, max, mean, and delta values of potential density at predefined pressure levels.\n", program);

   fprintf(stderr,"\nUsage:  %s filename_root(s) [-D<dirname>] [-E<file_extent>] [-h] ", program);

   fprintf(stderr,"\n\n  List of filenames MUST be first argument!");
   fprintf(stderr,"\n\n   OPTIONS: ");
   fprintf(stderr,"\n-D : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n        ex: -D../data/ ");
   fprintf(stderr,"\n-E : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n        ex: -E.dat ");
   fprintf(stderr,"\n-h :  help...... prints this message. ");
   fprintf(stderr,"\n\n");  
   return;
}

/*****************************************************************************/
void get_hydro_data(int file)
   /*  Reads each station in a HydroBase file and computes sigma values
       at each standard level.  This module requires that the HydroBase
       file contains a minimum of pr, te, sa observations.  */
{
   int error, i, j, npts, datagap;
   short too_deep;
   double  z;
   double  *x;


/* read each station in file ... */

   while ((error = get_station(file, &hdr, &station)) == 0) {

  /* ensure that pr, de, te, and sa are available ... */

       if (!(available(PR, &hdr) &&  available(TE, &hdr) && available(SA, &hdr))) {
         fprintf(stderr,"\n>>>>> WARNING!!!  ");
         fprintf(stderr,"Station does not include PR, TE, and SA.");
         break;
       }


 /* compute potential density properties at each level in station ... */

      npts = hdr.nobs;
      free_and_alloc(&station.observ[(int)S0], hdr.nobs);
      compute_sigma(0., hdr.nobs, station.observ[(int)S0], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);

      free_and_alloc(&station.observ[(int)S1], hdr.nobs);
      compute_sigma(1000., hdr.nobs, station.observ[(int)S1], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);

      free_and_alloc(&station.observ[(int)S2], hdr.nobs);
      compute_sigma(2000., hdr.nobs, station.observ[(int)S2], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);

      free_and_alloc(&station.observ[(int)S3], hdr.nobs);
      compute_sigma(3000., hdr.nobs, station.observ[(int)S3], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);

      free_and_alloc(&station.observ[(int)S4], hdr.nobs);
      compute_sigma(4000., hdr.nobs, station.observ[(int)S4], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);

/*  linearly interpolate density onto standard levels */

      x = station.observ[(int)PR];
      too_deep = 0;
    
            /* SIG0 */    
      i = 0;
      while (!too_deep && i < NLEV0) {
         z =  hb_linterp(lev0[i], x, station.observ[(int)S0], npts);
        
         if (z > -998. && npts > 1) {                  /* check for vertical datagaps */
            j = 0;
            while (x[++j] < lev0[i]) 
               ;

            if (((int)lev0[i] - (int)x[j-1] < 10) 
              || ( (int)x[j] - (int)lev0[i] < 10)) 
                  datagap = 0;

            else
                  datagap = (x[j] - x[j-1]) > 220;
                  
            if (!datagap) {
              if (z < min0[i])
                 min0[i] = z;
              if (z > max0[i])
                 max0[i] = z;
              sum0[i] += z;
              ++count0[i];   
            }
         }
       
         if (++i < NLEV0) {
            too_deep = lev0[i] > station.observ[(int)PR][npts-1];
         } 
      } /* end while */
   
 /* SIG1 */    
  
      i = 0;
      while (!too_deep && i < NLEV1) {
         z =  hb_linterp(lev1[i], x, station.observ[(int)S1], npts);
         
         if (z > -998.  && npts > 1) {                  /* check for vertical datagaps */
            j = 0;
            while (x[++j] < lev1[i]) 
               ;

            if (((int)lev1[i] - (int)x[j-1] < 10) 
              || ( (int)x[j] - (int)lev1[i] < 10)) 
                  datagap = 0;

            else if (lev1[i] < 1005)
                  datagap = (x[j] - x[j-1]) > 220;
            else       
                   datagap = (x[j] - x[j-1]) > 620;
                   
           if (!datagap) {
              if (z < min1[i])
                 min1[i] = z;
              if (z > max1[i])
                 max1[i] = z;
              sum1[i] += z;
              ++count1[i];   
            }
         }
       
         if (++i < NLEV1) {
            too_deep = lev1[i] > station.observ[(int)PR][npts-1];
         } 
      } /* end while */
     
    
  /* SIG2 */    
  
      i = 0;
      while (!too_deep && i < NLEV2) {
         z =  hb_linterp(lev2[i], x, station.observ[(int)S2], npts);
       
         if (z > -998.  && npts > 1) {                  /* check for vertical datagaps */
            j = 0;
            while (x[++j] < lev2[i]) 
               ;

            if ( ((int)x[j] - (int)lev2[i] < 10) 
              || ((int)lev2[i] - (int)x[j-1] < 10)) 
                  datagap = 0;

            else       
                   datagap = (x[j] - x[j-1]) > 620;
                   
           if (!datagap) {
              if (z < min2[i])
                 min2[i] = z;
              if (z > max2[i])
                 max2[i] = z;
              sum2[i] += z;
              ++count2[i];   
            }
         }
       
         if (++i < NLEV2) {
            too_deep = lev2[i] > station.observ[(int)PR][npts-1];
         } 
      } /* end while */


    
  /* SIG3 */    
  
      i = 0;
      while (!too_deep && i < NLEV3) {
         z =  hb_linterp(lev3[i], x, station.observ[(int)S3], npts);
       
         if (z > -998. && npts > 1) {                  /* check for vertical datagaps */
            j = 0;
            while (x[++j] < lev3[i]) 
               ;

            if ( ((int)x[j] - (int)lev3[i] < 10) 
              || ((int)lev3[i] - (int)x[j-1] < 10)) 
                  datagap = 0;
            else       
                   datagap = (x[j] - x[j-1]) > 620;
                   
           if (!datagap) {
              if (z < min3[i])
                 min3[i] = z;
              if (z > max3[i])
                 max3[i] = z;
              sum3[i] += z;
              ++count3[i];   
            }
         }
       
         if (++i < NLEV3) {
            too_deep = lev3[i] > station.observ[(int)PR][npts-1];
         } 
      } /* end while */


    
  /* SIG4 */    
  
      i = 0;
      while (!too_deep && i < NLEV4) {
         z =  hb_linterp(lev4[i], x, station.observ[(int)S4], npts);
       
         if (z > -998.  && npts > 1) {                  /* check for vertical datagaps */
            j = 0;
            while (x[++j] < lev4[i]) 
               ;

            if ( ((int)x[j] - (int)lev4[i] < 10) 
              || ((int)lev4[i] - (int)x[j-1] < 10)) 
                  datagap = 0;
            else       
                   datagap = (x[j] - x[j-1]) > 620;
                   
           if (!datagap) {
              if (z < min4[i])
                 min4[i] = z;
              if (z > max4[i])
                 max4[i] = z;
              sum4[i] += z;
              ++count4[i];   
            }
         }
       
         if (++i < NLEV4) {
            too_deep = lev4[i] > station.observ[(int)PR][npts-1];
         } 
      } /* end while */
    
      if ( station.observ[(int)PR][npts-1] > lev4[NLEV4-1]){
       fprintf(stdout, "%7.3f %8.3f %8.1lf %8.3lf %8.3lf %8.3lf\n", hdr.lat, hdr.lon,
               station.observ[(int)PR][npts-1], station.observ[(int)TE][npts-1],
               station.observ[(int)SA][npts-1], station.observ[(int)S4][npts-1]);
      }

   }  /* end while get_station */

   if (error > 0)
      report_status(error, stderr);

   return;

} /* end get_hydro_data() */

