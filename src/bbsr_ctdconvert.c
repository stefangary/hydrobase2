/*  bbsr_ctdconvert.c
................................................................................
                          *******  HydroBase *******
               ...................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
                             1993
................................................................................
................................................................................
.   Reads Bermuda Bio Station ctd format file(s)
.   extracts :    header info
.                 p,t,s,ox,  
.   applies a nominal filter to each (5 pts for p,t,s; 11 pts for ox)
.   and outputs p,d,t,s and o at the specified prs_int to the output
.   file.
.
.   USAGE: ctdconvert infile_list -Ooutfile -P<prs_int> [-D<dir>] [-E<extent>]
................................................................................
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "hydrobase.h"


#define   LF        0x0a   /* ascii code for linefeed char */
#define   MISSING   -9.0   /* missing value flag */

#define   DIR    ""
#define   EXTENT ""

  /* define possible modes for output files ... */
#define   OVERWRITE 0
#define   NO_CLOB   1
#define   APPEND    2


  /* globally defined variables */
  
struct HYDRO_HDR  hdr;
struct HYDRO_DATA data;
double *p, *d, *t, *s, *o;
char pr_avail, te_avail, sa_avail, ox_avail;
int deltap;          /* pressure interval of ctd cast */
int temp90;

   /* prototypes for locally define functions */

void print_usage(char *);
int readdata(FILE *);
void parse_header(FILE *);
FILE *openfile(char *, char *, char *);

int main (int argc, char **argv)
{
   int error, nobs, nprops;
   int outfile;
   int  i, j, curfile = 1, nfiles = 0; 
   short oflag, pflag; 
   char *dir, *extent;
   double dlat;
   int prsint, npts;
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
    pflag = oflag  = 0;
    deltap = 2.0;
 
 
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
                        pflag = 1;
                        error = (sscanf(&argv[i][2],"%d", &prsint) != 1);
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

   if (! (oflag && nfiles && pflag)) {
       fprintf(stderr,"\nYou must specify input file(s), output file and pressure interval.\n");
       exit(1);
   }

/* loop for each input file */

   do {
   
     infile = openfile(dir, argv[curfile], extent);
     if (infile == NULL) goto NEXTFILE;
      
      /* parse the header info */
      
     parse_header(infile);
     
     nprops = 1;   /* assume pr is available. */
     if (te_avail)   
         ++nprops;
     if (sa_avail)   
         ++nprops;
     if (ox_avail)   
         ++nprops;
     hdr.nprops = nprops + 1;    /* add depth to the output vars */
         
       
     
     /* adjust some of the header info ...   */
     
     dlat = (double) hdr.lat;
     strncpy(hdr.country,"32",3);
     strncpy(hdr.ship,"0G",3);
     for (i = 0; i < NQUAL; ++i)
       hdr.qual[i] = '0';
     hdr.origin = '5';
     hdr.instrument = 'c';
     hdr.pdr = 3015;
     hdr.ms10 = ms10(hdr.lat, hdr.lon, &hdr.ms1);
     hdr.prop_id = (int *) malloc(hdr.nprops * sizeof(int));
     hdr.prop_id[0] = (int)PR;
     hdr.prop_id[1] = (int)DE;
     hdr.prop_id[2] = (int)TE;
     hdr.prop_id[3] = (int)SA;
     if (ox_avail) 
        hdr.prop_id[4] = (int)OX;
     
     /*  allocate space */
     
     nobs = 3000;
     
     p = (double *) calloc(nobs, sizeof(double));
     d = (double *) calloc(nobs, sizeof(double));
     t = (double *) calloc(nobs, sizeof(double));
     s = (double *) calloc(nobs, sizeof(double));
     if (ox_avail) 
        o = (double *) calloc(nobs, sizeof(double));
        
    /* read and filter each property */    
        
     nobs = readdata(infile);

     /* decimate the arrays according to specified prsint */
     
     npts = prsint / deltap;  /*  get # of pts per interval*/
     if (prsint % deltap)
        fprintf(stderr,"WARNING: prsint (%d) requested is not an even multiple of the pressure sorted ctd file: output will be a %d db series\n", prsint, (npts*deltap));
     
     if (npts == 0) {
        fprintf(stderr,"Bad number of points per interval: %d \n", npts);
        exit(1);
     }
    
     j = 0;
     for (i = 0; i < nobs; i += npts) {
       p[j] = p[i];
       d[j] = hb_depth(p[i], dlat);
       t[j] = t[i];
       s[j] = s[i];
       if (ox_avail) {
          if (o[i] >= 0)
             o[j] = ox_kg2l(o[i], p[i], t[i], s[i]);
          else
             o[j] = MISSING;   
       }
       ++j;
     }
     
     /* add bottom observation  ... */
     
     if ((i - npts) != (--nobs)) {
        p[j] = p[nobs];
        d[j] = hb_depth(p[nobs], dlat);
        t[j] = t[nobs];
        s[j] = s[nobs];
        if (ox_avail) {
          if (o[nobs] >= 0)
             o[j] = ox_kg2l(o[nobs], p[nobs], t[nobs], s[nobs]);
          else
             o[j] = MISSING; 
        }  
        ++j;
     }
     
     hdr.nobs = data.nobs = j;
     data.nprops = hdr.nprops; 
     
     data.observ[(int)PR] = p;   /* set these pointers */
     data.observ[(int)DE] = d;
     data.observ[(int)TE] = t;
     data.observ[(int)SA] = s;
     if (ox_avail)
        data.observ[(int)OX] = o;
        
     if (hdr.nobs > 0 )
        write_hydro_station(outfile, &hdr, &data);
     
     free(hdr.prop_id);   
     free(p);
     free(d);
     free(t);
     free(s);
     if (ox_avail)
       free(o);

NEXTFILE:
         fclose(infile);

 
   } while (curfile++ < nfiles );

   exit(0);
        
} /* end main() */ 
  
/*****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nUsage:  %s filelist -Ooutfile -Pprs_int [-D<dirname>] [-E<file_extent>]", program);
   fprintf(stderr,"\n\n  List of filenames MUST be first argument!");
   fprintf(stderr,"\n    -D  : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n    -E  : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n    -O  : specifies output filename");
   fprintf(stderr,"\n    -P  : specifies pressure interval for output file");  
   fprintf(stderr,"\n            ex: -P10 ");

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
void parse_header(FILE *fptr)

   /* reads the header info from an already opened ctd data file 
     The info is used to fill in appropriate values in
      the struct HYDRO_HDR hdr, a global variable.  */
{
   int error, i;
   char line[120], *st;
   
   if (fscanf(fptr,"%[^\n]", line) != 1) {
        fprintf(stderr,"\n Error attempting to read 1st header line. \n");
        exit(1);   
   } 
   error = getc(fptr);
   
   /* search for header info ... */ 
   
   while ((st = strstr(line, "temperature")) == NULL) {
      fscanf(fptr,"%[^\n]", line);
      error = getc(fptr);
   }
   
   temp90 = 1;
   if (strstr(line,"90") == NULL) {
       if (strstr(line,"68") == NULL)
           fprintf(stderr,"\n Cannot determine whether temperature is ITS 90 or IPTS 68 \n");
       temp90 = 0;
   }

   while ((st = strstr(line, "cruise")) == NULL) {
      fscanf(fptr,"%[^\n]", line);
      error = getc(fptr);
   }
   st = strchr(st, '=');
   ++st;
   if (sscanf(st,"%d", &hdr.cruise) != 1) {
        fprintf(stderr,"\n error in sscanf attempt to read cruise \n");
        exit(1);
   }
   
   /* search for  cast number ... */ 
   fscanf(fptr,"%[^\n]", line);
   error = getc(fptr);
   while ((st = strstr(line, "cast")) == NULL) {
       fscanf(fptr,"%[^\n]", line);
       error = getc(fptr);
   }
   st = strchr(st, '=');
   ++st;
   if (sscanf(st,"%d", &hdr.station) != 1) {
        fprintf(stderr,"\n error in sscanf attempt to read cast # \n");
        exit(1);
   }
   
   /* search for start info ... */ 
   fscanf(fptr,"%[^\n]", line);
   error = getc(fptr);
   while ((st = strstr(line, "start")) == NULL) {
       fscanf(fptr,"%[^\n]", line);
       error = getc(fptr);
   }
   /* get date ... */
  
   if ((st = strstr(line, "ddmmyy")) != NULL ) {
      while (*st != '=')
         ++st;
    
    ++st; 
      if (sscanf(st,"%2d%2d%2d", &hdr.day, &hdr.month, &hdr.year) != 3) {
        fprintf(stderr,"\n error in sscanf attempt to read date \n");
        exit(1);
      } 
      hdr.year += 1900;
      if (hdr.year < 1950)
         hdr.year += 100;
   }
   else {
      if ((st = strstr(line, "yyyymmdd")) == NULL ) {
        fprintf(stderr,"\n error in sscanf attempt to read date # \n");
        exit(1);
      }	
       while (*st != '=')
         ++st;

      ++st; 
      if (sscanf(st,"%4d%2d%2d", &hdr.year, &hdr.month, &hdr.day) != 3) {
        fprintf(stderr,"\n error in sscanf attempt to read date \n");
        exit(1);
      } 
   }
      
   /* get position ... */
   st = strstr(st, " lat "); 
       while (*st != '=')
         ++st;
   ++st;
   if (sscanf(st,"%f", &hdr.lat) != 1) {
        fprintf(stderr,"\n error in sscanf attempt to read lat \n");
        exit(1);
   }
   st = strstr(st, "long "); 
       while (*st != '=')
         ++st;
   ++st;
   if (sscanf(st,"%f", &hdr.lon) != 1) {
        fprintf(stderr,"\n error in sscanf attempt to read lon \n");
        exit(1);
   }
   if (hdr.lon > 0)
      hdr.lon = - hdr.lon;
      
      
   /* search for variable descriptions...*/  
      
   fscanf(fptr,"%[^\n]", line);
   error = getc(fptr);
   
   while (strstr(line, "/variable") == NULL) {
        fscanf(fptr,"%[^\n]", line);
        error = getc(fptr);
   }
   
   fscanf(fptr,"%[^\n]", line);
   error = getc(fptr);
   
   pr_avail = (strstr(line, "pres") == NULL) && (strstr(line,"P,") == NULL)? 0 : 1;
   te_avail = (strstr(line, "temp") == NULL)  && (strstr(line,"t,") == NULL) ? 0 : 1;
   sa_avail = (strstr(line, "sal") == NULL)  && (strstr(line,"s,") == NULL)  ? 0 : 1;
   ox_avail = (strstr(line, " ox") == NULL)  && (strstr(line,"O2,") == NULL)   ? 0 : 1;
 
    if (!pr_avail)
       fprintf(stderr,"\nWARNING:  no PR in this file.");
    if (!te_avail)
       fprintf(stderr,"\nWARNING:  no TE in this file.");
    if (!sa_avail)
       fprintf(stderr,"\nWARNING:  no SA in this file.");
       
   while (strstr(line, "/data") == NULL) {
        fscanf(fptr,"%[^\n]", line);
        error = getc(fptr);
   }
        fscanf(fptr,"%[^\n]", line);   /* file now positioned at start of data */
        error = getc(fptr);

   return;
}  /* end parse_header() */
/*****************************************************************************/
int readdata(FILE *fptr)
{
   int nobs, i, error, propcount, nread;
   char line[80];
   double buffer[4];
   
   propcount = 1;
   if (te_avail) 
      ++propcount;
   if (sa_avail) 
      ++propcount;
   if (ox_avail) 
      ++propcount;
   
   nobs = 0;
   while ((nread = fscanf(fptr, "%[^\n]", line)) != EOF) {
      if (nread != 1){
         fprintf(stderr,"\nError reading properties at data line: \n%s\n", line);
         exit(1);
      }
      error = getc(fptr);
      
      nread = sscanf(line,"%lf%lf%lf%lf",&buffer[0],&buffer[1],&buffer[2],&buffer[3]);
      
      i = 0;
      p[nobs] = buffer[i];
      if (te_avail) {
         t[nobs] = buffer[++i];
      } 

      if (sa_avail){ 
         s[nobs] = buffer[++i];
      }  
      if (ox_avail) {
         o[nobs] = buffer[++i];
      }
      
      if (nread < propcount){
         fprintf(stderr,"\nOnly able to parse %d properties at data line: \n%s\n", nread, line);
         exit(1);
      }
      if ((p[nobs] > -8) && (t[nobs] > -8) && (s[nobs] > -8) ) {
         if (temp90) {
             t[nobs] *= 1.00024;
	 }
         ++nobs;
      }
   }
   return (nobs);
}  /* end readdata() */
