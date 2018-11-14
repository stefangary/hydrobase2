/*  mmp_ctdconvert.c
................................................................................
 Converts ascii files of profiles from moored data (MMP or VGR) that are
 of the form:
 
  date pr sa te gn vdo vcr
  
  into HydroBase station format.
.................................................................................
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "hydrobase.h"
#include "hb_memory.h"

#define   NPOUT     7

#define   DIR    ""
#define   EXTENT ""

  /* globally defined variables */
  
struct HYDRO_HDR  hdr;
struct HYDRO_DATA data;
double *p, *d, *t, *s, *g, *vdo, *vcr;

  /* prototypes for locally define functions */
  
void print_usage(char *);
FILE *openfile(char *, char *, char *);
int readdata(FILE *);

int main ( int argc, char **argv)
{
   int error, nobs;
   int outfile;
   int  i, j, curfile, nfiles;
   char *dir, *extent, *st;
   double dlat;
   int  npts;
   FILE *infile;
   
/* check for command line arguments */

   if (argc < 2) {
     print_usage(argv[0]);
     exit(1);
   }
 
/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    error = 0;
    outfile = STDOUT;
    dlat = -99.0;
    curfile = 1;
    nfiles = 0;
    
    /* set header info */
    hdr.cruise = 999;
    strncpy(hdr.ship, "MP", 3);
    strncpy(hdr.country, "32", 3);
    hdr.instrument = 'm';
    hdr.origin = '9';
    hdr.station = 1;
    hdr.pdr = 0;
    hdr.qual[0] = '0';
    hdr.qual[1] = '0';
    hdr.qual[2] = '0';
    hdr.qual[3] = '1';
    hdr.prop_id = (int *) calloc((size_t)NPOUT, sizeof(int));
    hdr.prop_id[0] = (int)PR;
    hdr.prop_id[1] = (int)DE;
    hdr.prop_id[2] = (int)TE;
    hdr.prop_id[3] = (int)SA;
    hdr.prop_id[4] = (int)GN;
    hdr.prop_id[5] = (int)VE;
    hdr.prop_id[6] = (int)VN;
 
/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'C':                    /* get cruise number */
                        st = &argv[i][2];
                        if (*st == '/' )
                          ++st;
                        error = (sscanf(st,"%d", &hdr.cruise) != 1);
                        break;
               case 'D':                   /* get input dir */
                        dir = &argv[i][2];
                        break;
               case 'E':                    /* get file extent */
                        extent = &argv[i][2];
                        break;
               case 'I':                    /* get ship code  */
                        st = &argv[i][2];
                        if (*st == '/' )
                          ++st;
                        error = (sscanf(st,"%2s", &hdr.ship) != 1);
                        break;
               case 'L':                    /* get lat/lon */
                        error = (sscanf(&argv[i][2],"%f/%f", &hdr.lat, &hdr.lon) != 2);
			dlat = (double) hdr.lat;
			hdr.ms10 = ms10(hdr.lat, hdr.lon, &hdr.ms1);
                        break;

                case 'O':                    /* get output file  */
                        outfile = create_hydro_file(&argv[i][2], OVERWRITE);
                        if (outfile < 1) {
                           fprintf(stderr,"\nError opening output file: %s\n", 
                                  &argv[i][2]);
                           exit(1);
                        }
                        break;

               case 'S':                    /* get starting station number */
                        st = &argv[i][2];
                        error = (sscanf(st,"%d", &hdr.station) != 1);
                        break;
               case 'W':                    /* get water_depth */
                        st = &argv[i][2];
                        error = (sscanf(st,"%d", &hdr.pdr) != 1);
                        break;

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

   if (! (nfiles)) {
       fprintf(stderr,"\nYou must specify input file(s)  \n");
       exit(1);
   }
    
     /*  allocate space */
     nobs = 5000;  
     p = (double *) get_memory((void *)NULL, (size_t)nobs, sizeof(double));
     d = (double *) get_memory((void *)NULL, (size_t)nobs, sizeof(double));
     t = (double *) get_memory((void *)NULL, (size_t)nobs, sizeof(double));
     s = (double *) get_memory((void *)NULL, (size_t)nobs, sizeof(double));
     g = (double *) get_memory((void *)NULL, (size_t)nobs, sizeof(double));
     vdo = (double *) get_memory((void *)NULL, (size_t)nobs, sizeof(double));
     vcr = (double *) get_memory((void *)NULL, (size_t)nobs, sizeof(double));

     data.observ[(int)PR] = p;   /* set these pointers */
     data.observ[(int)DE] = d;
     data.observ[(int)TE] = t;
     data.observ[(int)SA] = s;
     data.observ[(int)GN] = g;
     data.observ[(int)VE] = vdo;
     data.observ[(int)VN] = vcr;
     
     --hdr.station;  

/* loop for each input file */

   do {
   
     infile = openfile(dir, argv[curfile], extent);
     if (infile == NULL) 
           goto NEXTFILE;
      
     npts = readdata(infile);
     
     if (npts  > 0 )
         for (i = 0; i < npts; ++i) {
             d[i] = hb_depth(p[i], dlat);
	     vcr[i] *= 0.01;
	     vdo[i] *= 0.01;  /* convert units to m/s */
         }
     
     /* adjust some of the header info ...   */
         ++hdr.station;
         hdr.nobs = data.nobs = npts;
         data.nprops = hdr.nprops = NPOUT; 
     
         write_hydro_station(outfile, &hdr, &data);
         
NEXTFILE:
         fclose(infile);

   } while (curfile++ < nfiles );

   fprintf(stderr,"End of %s\n", argv[0]);
   exit(0);
        
} /* end main() */ 
  
/*****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nUsage:  %s filelist  -I<mp_id> -L<lat>/<lon> [-C<cruise_id>]   [-D<dirname>] [-E<file_extent>] [-O<outfile>]  [-S<starting_station>]  [-W<water_depth>]", program);
   fprintf(stderr,"\n\n  List of filenames must be first argument.");
   fprintf(stderr,"\n    -I  :  lD (2-char) of this moored profiler  (e.g. -IW1)");
   fprintf(stderr,"\n    -L  :  lat/lon for this moored profiler ");
   fprintf(stderr,"\n OPTIONAL Arguments: ");
   fprintf(stderr,"\n   [-C]  : specifies (deployment) cruise id.  default = 999 ");
   fprintf(stderr,"\n   [-D] : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n   [-E] : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n   [-O] : specifies output file ");
   fprintf(stderr,"\n   [-S] : specifies starting station # [1] which is incremented for each file");
   fprintf(stderr,"\n   [-W] : specifies water depth of mooring");

   fprintf(stderr,"\n\n");  
   return;
}
   /*****************************************************************************/
FILE *openfile(char *dir, char *root, char *extent)

{
   char st[80];
   int i;
   FILE *infile;
   
   strcpy(st, dir);
   if ((i=strlen(dir)) != 0) {
      if (dir[i-1] != '/')
         strncat(st,"/",1);
   }
   strcat(st, root);
   strcat(st, extent);
   infile = fopen(st,"r");
   if (infile == NULL)
         fprintf(stderr,"\n Unable to open %s \n", st);
   else
         fprintf(stderr,"\n   opened %s ...\n", st);
   
   return(infile);
   
}  /* end openfile() */
/*****************************************************************************/
int readdata(FILE *fptr)

/* Reads  pressure-sorted property values for an entire profile and returns the 
   number of observations. */
 
{
   int i, n, error;
   char line[1000];
   char datestring[12];
   char month[4];
   
   n = 0;
   i = 0;
   
   while ( fscanf(fptr, "%[^\n]", line) != EOF) {
      ++i;
      getc(fptr);
      
      if (strstr(line, "NaN") == NULL)  {
           error = sscanf(line, "%s %lf %lf %lf %lf %lf %lf", datestring, &p[n],&s[n], &t[n],&g[n],&vdo[n], &vcr[n]) != 7;
     
         if (error) {
            fprintf(stderr,"\nError parsing data line #%d\n %s\n", i, line);
            exit(1);
         }
         ++n;
      } /* end if */
   } /* end while */
   
   if (sscanf(datestring, "%2d-%3s-%4d", &hdr.day, month, &hdr.year) != 3 ) {
            fprintf(stderr,"\nError parsing date %s\n", datestring);
	    exit(1);
   }
   
   switch (month[0]) {
      case 'A':
          switch (month[1]) {
	     case 'u':
	        hdr.month = 8;
		break;
	     case 'p':
	        hdr.month = 4;
		break;
	     default:
	        fprintf(stderr, "Unrecognized month %s\n", month);
		exit(1);
	  }
          break;
      case 'D':
         switch (month[1]){
	     case 'e':
	        hdr.month = 12;
		break;
	     default:
	        fprintf(stderr, "Unrecognized month %s\n", month);
		exit(1);
	  }
          break;
      case 'F':
         switch (month[1]) {
	     case 'e':
	        hdr.month = 2;
		break;
	     default:
	        fprintf(stderr, "Unrecognized month %s\n", month);
		exit(1);
	  }
          break;
       case 'J':
          switch (month[1]) {
	     case 'a':
	        hdr.month = 1;
		break;
	     case 'u':
	        hdr.month = 7;
		if (month[2] == 'n')
		    hdr.month = 6;
		break;
	     default:
	        fprintf(stderr, "Unrecognized month %s\n", month);
		exit(1);
	  }
          break;
      case 'M':
          switch (month[2]) {
	     case 'r':
	        hdr.month = 3;
		break;
	     case 'y':
	        hdr.month = 5;
		break;
	     default:
	        fprintf(stderr, "Unrecognized month %s\n", month);
		exit(1);
	  }
          break;
      case 'N':
         switch (month[1]) {
	     case 'o':
	        hdr.month = 11;
		break;
	     default:
	        fprintf(stderr, "Unrecognized month %s\n", month);
		exit(1);
	  }
          break;
      case 'O':
         switch (month[1]) {
	     case 'c':
	        hdr.month = 10;
		break;
	     default:
	        fprintf(stderr, "Unrecognized month %s\n", month);
		exit(1);
	  }
          break;
      case 'S':
         switch (month[1]) {
	     case 'e':
	        hdr.month = 9;
		break;
	     default:
	        fprintf(stderr, "Unrecognized month %s\n", month);
		exit(1);
	  }
          break;
      default:
	  fprintf(stderr, "Unrecognized month %s\n", month);
          exit(1);
    } /* end switch */
    
   return(n);
   
}  /* end readdata() */
