/*  awi_convert.c
................................................................................
.   Reads Alfred Wegener Institute 's CTDA data files
.   extracts :    header info
.                 p,t,s,o2 
.   and outputs p,d,t,s and o2 at the specified prs_int to the output
.   file.
.   The names of ctd files and header info are contained in a single summary file which must be supplied as an argument.
.
.   USAGE: awi_convert -Fsummary_file -Ooutfile -Sshipcode -Ncountrycode -Ccruise_code -P<prs_int>
           [-D<dir>] 
................................................................................
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "hydrobase.h"
#include "hb_memory.h"

#define   NPOUT     5
#define   DELTAP    1     /* predetermined pressure interval  */

#define   DIR    ""
#define   EXTENT ""



  /* globally defined variables */
  
struct HYDRO_HDR  hdr;
struct HYDRO_DATA data;
double *p, *d, *t, *s, *o;
int ox_available;

  /* prototypes for locally define functions */
  
void print_usage(char *);
FILE *openfile(char *, char *, char *);
int readdata(FILE *, int *);
int parse_header(FILE *);
int parse_date(char *);


int main ( int argc, char **argv)
{
   int error, nobs, nprops, npout;
   int outfile;
   int  i, j;
   int nread; 
   int sflag, oflag, nflag, cflag; 
   int n_ox_available; 
   char *dir, *extent, *st, buffer[5000];
   char fname[1000];
   char datestr[12];
   double dlat, prsint;
   int  npts;
   
   FILE *ctdfile, *sumfile;
   
/* check for command line arguments */

   if (argc < 2) {
     print_usage(argv[0]);
     exit(1);
   }
   
 
/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    error = 0;
    cflag = nflag = sflag = oflag  = 0;
    npout = NPOUT;
    prsint = 10.0;
    sumfile = NULL;
     
/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'C':                    /* get cruise number */
                        cflag = 1;
                        st = &argv[i][2];
                        if (*st == '/' )
                          ++st;
                        error = (sscanf(st,"%d", &hdr.cruise) != 1);
                        break;
               case 'F':                    /* get summary file*/
	                sumfile = fopen(&argv[i][2], "r");
			if (sumfile == NULL) {
			    fprintf(stderr,"Unable to open %s for input\n", &argv[i][2]);
			    exit(1);
			}
			    fprintf(stderr,"Opened %s \n", &argv[i][2]);
                        break;

               case 'N':                    /* get country code */
                        nflag = 1;
                        st = &argv[i][2];
                        if (*st == '/' )
                          ++st;
                        error = (sscanf(st,"%2s", hdr.country) != 1);
                        break;
                case 'O':                    /* get output file  */
                        oflag = 1;
                        outfile = create_hydro_file(&argv[i][2], OVERWRITE);
                        if (outfile < 1) {
                           fprintf(stderr,"\nError opening output file: %s\n", 
                                  &argv[i][2]);
                           exit(1);
                        }
                        break;

               case 'P':                    /* get pressure interval */
                        error = (sscanf(&argv[i][2],"%lf", &prsint) != 1);
                        break;

               case 'S':                    /* get ship code */
                        sflag = 1;
                        st = &argv[i][2];
                        if (*st == '/' )
                          ++st;
                        error = (sscanf(st,"%2s", hdr.ship) != 1);
                        break;

               case 'h':                    /* help */
	          print_usage(argv[0]);
		  exit(1);
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
             fprintf(stderr,"\nError parsing command line args.\n");
             fprintf(stderr,"     in particular: '%s'\n", argv[i]);
             exit(1);
       }

   }  /* end for */

   if (! (oflag )) {
       fprintf(stderr,"\nYou must specify an output file \n");
       exit(1);
   }
   if (sumfile == NULL) {
       fprintf(stderr,"\nYou must specify an input file containing summary of ctd files/positions \n");
       exit(1);
   }
   if (! sflag) {
       fprintf(stderr,"\nYou must specify a 2-char ship code \n");
       exit(1);
   }
   if (! nflag) {
       fprintf(stderr,"\nYou must specify a 2-char country code \n");
       exit(1);
   }
   if (! cflag) {
       fprintf(stderr,"\nYou must specify a cruise number \n");
       exit(1);
   }
   
   fprintf(stderr,"\nOutput pressure series increment: %4.1lf db\n", prsint);
    
/* loop for each input data file */

  while (fscanf(sumfile, "%[^\n]", buffer) != EOF) {
     getc(sumfile);   /* move past LF */
     
     nread = sscanf(buffer, "%s%d%*d%f%f%d%*d%s", fname, &hdr.station, &hdr.lat, &hdr.lon, &hdr.pdr, datestr); 
     ctdfile = openfile(dir, fname, extent);
     if (ctdfile == NULL) goto NEXTFILE;
      
       
     
     /*  allocate space */
     
     p = (double *) get_memory((void *)NULL, (size_t)7000, sizeof(double));
     d = (double *) get_memory((void *)NULL, (size_t)7000, sizeof(double));
     t = (double *) get_memory((void *)NULL, (size_t)7000, sizeof(double));
     s = (double *) get_memory((void *)NULL, (size_t)7000, sizeof(double));
     o = (double *) get_memory((void *)NULL, (size_t)7000, sizeof(double));
    
    /* read each property */    
        
     n_ox_available = readdata(ctdfile, &nobs);
     npout = NPOUT;
     if (n_ox_available == 0)
        npout = 4;

     /* decimate the arrays according to specified prsint */
     
     npts = NINT(prsint / DELTAP);  /*  get # of pts per interval*/
     if (npts == 0) {
        fprintf(stderr,"Bad number of points per interval: %d \n", npts);
        exit(1);
     }
    
     j = 0;
     for (i = 0; i < nobs; i += npts) {
       p[j] = p[i];
       d[j] = hb_depth(p[i], (double)hdr.lat);
       t[j] = t[i];
       s[j] = s[i];
       if (ox_available)
          o[j] = o[i];
       ++j;
     }
     
     /* add bottom observation  ... */
     
     if ((i - npts) != (--nobs)) {
        p[j] = p[nobs];
        d[j] = hb_depth(p[nobs], dlat);
        t[j] = t[nobs];
        s[j] = s[nobs];
        if (ox_available)
            o[j] = o[nobs];
        ++j;
     }
     
     /* adjust some of the header info ...   */
     
     hdr.ms10 = ms10(hdr.lat, hdr.lon, &hdr.ms1);
     hdr.prop_id = (int *) calloc((size_t)npout, sizeof(int));
     hdr.prop_id[0] = (int)PR;
     hdr.prop_id[1] = (int)DE;
     hdr.prop_id[2] = (int)TE;
     hdr.prop_id[3] = (int)SA;
     if (npout == 5)
        hdr.prop_id[4] = (int)OX;
     hdr.nprops = npout;
     
     hdr.origin = '9';
     hdr.instrument = 'c';
     hdr.qual[0] = 0;
     hdr.qual[1] = 0;
     hdr.qual[2] = 0;
     hdr.qual[3] = 1;
     hdr.nobs = data.nobs = j;
     data.nprops = hdr.nprops; 
     
     error = parse_date(datestr);
     if (error)
        fprintf(stderr, "\n Error parsing date string: %s\n", datestr);
     
     data.observ[(int)PR] = p;   /* set these pointers */
     data.observ[(int)DE] = d;
     data.observ[(int)TE] = t;
     data.observ[(int)SA] = s;
     if (npout > 4)
        data.observ[(int)OX] = o;
        
     if (hdr.nobs > 0 )
        write_hydro_station(outfile, &hdr, &data);
     
     free(hdr.prop_id);   
     free(p);
     free(d);
     free(t);
     free(s);
     free(o);

NEXTFILE:
         fclose(ctdfile);

 
   } /* end while */

   fclose(sumfile);
   fprintf(stderr,"End of %s\n", argv[0]);
   exit(0);
        
} /* end main() */ 
  
/*****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nUsage:  %s -F<sum_file_name> -C<cruise> -N<nation> -O<outfile>  -P<prs_int> -S<shipcode> ", program);
   fprintf(stderr,"\n    -C  : specifies cruise number ");
   fprintf(stderr,"\n    -F : specifies file containing list of input files and station positions");
   fprintf(stderr,"\n    -N  : 2-char country code ");
   fprintf(stderr,"\n    -O  : specifies output file ");
   fprintf(stderr,"\n    -P  : specifies pressure interval of output series ");
   fprintf(stderr,"\n    -S  : 2-char ship code ");

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
int parse_date(char *str)
    /* Returns 0 for successful parsing, -1 for an error */
{
  
   sscanf(str, "%d", &hdr.day);
   sscanf(&str[7], "%d", &hdr.year);
  
  
  switch (str[3]) {
     case 'A':
     case 'a':
         switch (str[4]) {
	    case 'p':
	    case 'P':
	      hdr.month = 4;
	      break;
	    case 'u':
	    case 'U':
	      hdr.month = 8;
	      break;
	    default:
	      return(-1);
	 } /* end switch */  
	    
         break;
	 
     case 'D':
     case 'd':
         hdr.month = 12;
	 break;
	 
     case 'F':
     case 'f':
         hdr.month = 2;
	 break;
     case 'J':
     case 'j':
         switch (str[5]) {
	    case 'l':
	    case 'L':
	      hdr.month = 7;
	      break;
	    case 'n':
	    case 'N':
	      hdr.month = 6;
	      if (str[4] == 'a' || str[4] == 'A')
	         hdr.month = 1;
	      break;
	    default:
	      return(-1);
	 } /* end switch */  
	    
         break;
     case 'M':
     case 'm':
         switch (str[5]) {
	    case 'r':
	    case 'R':
	      hdr.month = 3;
	      break;
	    case 'y':
	    case 'Y':
	      hdr.month = 5;
	      break;
	    default:
	      return(-1);
	 } /* end switch */  
	    
         break;
     case 'N':
     case 'n':
         hdr.month = 11;
	 break;
	 
     case 'O':
     case 'o':
         hdr.month = 10;
	 break;
	 
     case 'S':
     case 's':
         hdr.month = 9;
	 break;
	 
     default:
        return(-1);
  }
  return(0);

} /*end parse_date() */
/*****************************************************************************/
int readdata(FILE *fptr, int *nobs_addr)

/* Reads ctd p,t,s,ox for an entire station and returns the 
   number of oxygen observations. */
 
{
   int i, n, error, nobs, nox_avail;
   char line[120];
   
   n = 0;
   i = 0;
   nox_avail = 0;

   /* skip first two header records */
   
    fscanf(fptr, "%[^\n]", line);
      error = getc(fptr);
    fscanf(fptr, "%[^\n]", line);
      error = getc(fptr);
     
   while ( fscanf(fptr, "%[^\n]", line) != EOF) {
      error = getc(fptr);
     
             error = sscanf(line, "%d%lf%lf%*lf%lf%*lf%*lf%lf", &i, &p[n], &t[n], &s[n], &o[n]) != 5;
             if (o[n] > 0)
                 ++nox_avail;
     
      if (error) {
         fprintf(stderr,"\nError parsing data line #%d.\n", i);
         exit(1);
      }
      
      ++n;
   } /* end while */
   
   *nobs_addr = n;
   return(nox_avail);
   
}  /* end readdata() */
