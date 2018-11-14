/*  primer_ctdconvert.c
................................................................................
.   Reads WHOI CTD Group's CTDA data files
.   extracts :    header info
.                 p,t,s,ox 
.   applies a nominal filter to each (5 pts for p,t,s;)
.   and outputs p,d,t,s and ox at the specified prs_int to the output
.   file.
.
.   USAGE: ctda_convert infile_list -Ooutfile -Sshipcode -Ncountrycode -Ccruise_code -P<prs_int>
          [-F<filtwidth>[/<ox_filtwidth>]] [-D<dir>] [-E<extent>]
................................................................................
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "hydrobase.h"
#include "hb_memory.h"
#include "hb_filters.h"

#define   FILTWIDTH  5     /* # of pts in gaussian filter */
#define   OFILTWIDTH  11     /* # of pts in gaussian filter */
#define   NPOUT     5
#define   DELTAP    2     /* predetermined pressure interval  */

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


int main ( int argc, char **argv)
{
   int error, nobs, nprops, npout;
   int outfile;
   int  i, j, curfile, nfiles;
   int filtwidth, ofiltwidth; 
   int sflag, oflag, nflag, cflag; 
   int n_ox_available; 
   char *dir, *extent, *st;
   double dlat, prsint;
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
    cflag = nflag = sflag = oflag  = 0;
    npout = NPOUT;
    prsint = 10.0;
    filtwidth = FILTWIDTH;
    ofiltwidth = OFILTWIDTH;
    curfile = 1;
    nfiles = 0;
 
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
               case 'D':                   /* get input dir */
                        dir = &argv[i][2];
                        break;
               case 'E':                    /* get file extent */
                        extent = &argv[i][2];
                        break;
               case 'F':                    /* get filter width */
                        error = (sscanf(&argv[i][2],"%d", &filtwidth) != 1);
                        st = &argv[i][2];
                        while (*(st++) != '\0') {
                         if (*st == '/') {
                           ++st;
                           error = (sscanf(st,"%d", &ofiltwidth) == 1) ? 0 : 1;
                           break;
                         }
                        }
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

   if (! (oflag && nfiles)) {
       fprintf(stderr,"\nYou must specify input file(s) and an output file \n");
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
    
/* loop for each input file */

   do {
   
     infile = openfile(dir, argv[curfile], extent);
     if (infile == NULL) goto NEXTFILE;
      
      /* parse the header info */
      
     nobs = parse_header(infile);
       
     
     /*  allocate space */
     
     p = (double *) get_memory((void *)NULL, (size_t)nobs, sizeof(double));
     d = (double *) get_memory((void *)NULL, (size_t)nobs, sizeof(double));
     t = (double *) get_memory((void *)NULL, (size_t)nobs, sizeof(double));
     s = (double *) get_memory((void *)NULL, (size_t)nobs, sizeof(double));
     o = (double *) get_memory((void *)NULL, (size_t)nobs, sizeof(double));
    
    /* read and filter each property */    
        
     n_ox_available = readdata(infile, &nobs);
     npout = NPOUT;
     if (n_ox_available == 0)
        npout = 4;
     filter("g",t, nobs, filtwidth);
     filter("g",s, nobs, filtwidth);
     if (n_ox_available == nobs)
       filter("g",o, nobs, ofiltwidth);
        

     /* decimate the arrays according to specified prsint */
     
     npts = NINT(prsint / DELTAP);  /*  get # of pts per interval*/
     if (npts == 0) {
        fprintf(stderr,"Bad number of points per interval: %d \n", npts);
        exit(1);
     }
    
     j = 0;
     dlat = (double) hdr.lat;
     for (i = 0; i < nobs; i += npts) {
       p[j] = p[i];
       d[j] = hb_depth(p[i], dlat);
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
     
     hdr.origin = '6';
     hdr.instrument = 'c';
     hdr.qual[0] = 0;
     hdr.qual[1] = 0;
     hdr.qual[2] = 0;
     hdr.qual[3] = 1;
     hdr.nobs = data.nobs = j;
     data.nprops = hdr.nprops; 
     
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
         fclose(infile);

 
   } while (curfile++ < nfiles );

   fprintf(stderr,"End of %s\n", argv[0]);
   exit(0);
        
} /* end main() */ 
  
/*****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nUsage:  %s filelist -C<cruise> -N<nation> -O<outfile>  -P<prs_int> -S<shipcode> [-F<filtwidth>[/<ox_filtwidth>]] [-D<dirname>] [-E<file_extent>]", program);
   fprintf(stderr,"\n\n  List of filenames must be first argument.");
   fprintf(stderr,"\n    -C  : specifies cruise number ");
   fprintf(stderr,"\n    -N  : 2-char country code ");
   fprintf(stderr,"\n    -O  : specifies output file ");
   fprintf(stderr,"\n    -P  : specifies pressure interval of output series ");
   fprintf(stderr,"\n    -S  : 2-char ship code ");
   fprintf(stderr,"\n   [-F] : specifies filterwidth (npts) or t,s and ox (default is 5/11) ");
   fprintf(stderr,"\n   [-D] : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n   [-E] : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");

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
int parse_header(FILE *fptr)
   /* reads the header info from an already opened ctd data file and returns
      the number of observation records.  The info is used to fill in
      appropriate values in the struct HYDRO_HDR hdr, a global variable.  */
{
   int n,  nobs, i, error, cast;
   char *line, *st;

   line = (char *) calloc(256, sizeof(char));
   
/* line 1 */

   if (fscanf(fptr,"%[^\n]", line) != 1) {
        fprintf(stderr,"\n Error attempting to read 1st header line. \n");
        exit(1);   
   } 
   error = getc(fptr);
   st = strstr( line, "STAT" );
   if (st == NULL) {
        fprintf(stderr,"\n error in sscanf attempt to find station #\n");
        exit(1);
   }
   if (sscanf(st+5,"%d", &hdr.station) != 1) {
        fprintf(stderr,"\n error in sscanf attempt to read station #\n");
        exit(1);
   } 
   
/* line 2 */

   if (fscanf(fptr,"%[^\n]", line) != 1) {
        fprintf(stderr,"\n Error attempting to read 2nd header line. \n");
        exit(1);   
   }    
   error = getc(fptr);
      
   st = &line[5];
   if (sscanf(st,"%d-%d-%d", &hdr.year, &hdr.month, &hdr.day) != 3) {
        fprintf(stderr,"\n error in sscanf attempt to read Date \n");
        exit(1);
   }
   if (hdr.year < 100)
     hdr.year += 1900;
         
/* line 3 */
   error = fscanf(fptr,"%[^\n]", line);
   error = getc(fptr);
        
/* line 4 */
   error = fscanf(fptr,"%[^\n]", line);
   error = getc(fptr);
   if (sscanf(&line[9],"%f", &hdr.lat) != 1) {
        fprintf(stderr,"\n error in sscanf attempt to read decimal lat \n");
        exit(1);
   }
   if (sscanf(&line[23],"%f", &hdr.lon) != 1) {
        fprintf(stderr,"\n error in sscanf attempt to read decimal lon \n");
        exit(1);
   }
   
/* line 5.. until property codes */
   while ((st = strstr(line,"OBS=")) == NULL) {
      error = fscanf(fptr,"%[^\n]", line);
      error = getc(fptr);
   }
   
   if (sscanf(&line[5],"%d", &nobs) != 1) {
        fprintf(stderr,"\n error in sscanf attempt to read OBS= \n");
        exit(1);
   }
   
   
/* Get line with property codes */
   error = fscanf(fptr,"%[^\n]", line);
   error = getc(fptr);
       
   
   if (strstr(line, "PR") == NULL) {
        fprintf(stderr,"\n cannot find PR in header \n");
        exit(1);
   }   
   if (strstr(line, "TE") == NULL) {
        fprintf(stderr,"\n cannot find TE in header \n");
        exit(1);
   }   
   if (strstr(line, "SA") == NULL) {
        fprintf(stderr,"\n cannot find SA in header \n");
        exit(1);
   } 
   ox_available = 1;  
   if (strstr(line, "OX") == NULL) {
        ox_available = 0;
    }   
  
/* ctd data file is now positioned at start of data records. */

   return(nobs);   
}  /* end parse_header() */
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
   nobs = *nobs_addr;
   
   while ( fscanf(fptr, "%[^\n]", line) != EOF) {
      ++i;
      error = getc(fptr);
     
      if (ox_available)  {
             error = sscanf(line, "%lf %lf %lf %lf", &p[n],&t[n],&s[n],&o[n]) != 4;
             if (o[n] > 0)
                 ++nox_avail;
      }	 
      else {
             error = sscanf(line, "%lf %lf %lf", &p[n],&t[n],&s[n]) != 3;
     }
     
      if (error) {
         fprintf(stderr,"\nError parsing data line #%d.\n", i);
         exit(1);
      }
      
      ++n;
   } /* end while */
   
   *nobs_addr = n;
   return(nox_avail);
   
}  /* end readdata() */
