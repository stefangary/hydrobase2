/*  a20_ctdconvert.c
................................................................................
                          *******  HydroBase2 *******
               ...................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
                             updated to ANSI C July 2001
................................................................................
................................................................................
.   Reads ctd   file(s) (from Scripps ODF CTD group)
.   Header info stored in separate file
.
.   USAGE: a20_ctdconvert infile_list  -Ooutfile -P<prs_int> -Sshipcode -Nnationcode -Ccruise_no [-D<dir>] [-E<extent>] [-Hhdr.file]
................................................................................
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "hydrobase.h"


#define   DELTAP 2

#define   DIR    ""
#define   EXTENT ""
#define   HEADERPATH  "/harpo/Shared/a20stations.txt"

#define NCTDPROPS 7
#define MAXHDRS 100000



  /* globally defined variables */
  
struct HYDRO_HDR  hdr;
struct HYDRO_DATA data;
double *p, *d, *t, *s, *o;
int ngood_ox;
int deltap;          /* pressure interval of ctd cast */
int nscans, nctdprops;
int cruiseflag;

struct INFO   {
    int year, month, day;
    int depth;
    float lat, lon;
};
struct INFO info[MAXHDRS];

  /* prototypes for locally define functions */
  
void print_usage(char *);
FILE *openfile(char *, char *, char *);
int readdata(FILE *);
void parse_header(FILE *);
   

int main (int argc, char **argv)
{
   int error, nobs, nprops, index;
   int outfile;
   int sta, cast;
   int  i, j, curfile = 1, nfiles = 0; 
   short oflag, pflag, sflag, nflag, hflag; 
   char *dir, *extent, *st;
   double dlat;
   int prsint, npts;
   FILE *infile, *hdrfile;
   
/* check for command line arguments */

   if (argc < 2) {
     print_usage(argv[0]);
     exit(1);
   }
   
 
/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    error = 0;
    pflag = oflag =sflag = nflag =   hflag = 0;
    deltap = 2.0;
    cruiseflag = 0;
    nscans = 6000;  /* dimension for ctd profile arrays */
 
 
/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'C':                    /* get cruise number */
                        cruiseflag = 1;
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

               case 'H':                    /* get header file  */
                        hflag = 1;
                        hdrfile = fopen(&argv[i][2],"r");
                        if (hdrfile == NULL) {
                           fprintf(stderr,"\nError opening header file: %s\n", 
                                  &argv[i][2]);
                           exit(1);
                        }
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
               case 'N':                    /* get country code */
                        nflag = 1;
                        st = &argv[i][2];
                        if (*st == '/' )
                          ++st;
                        error = (sscanf(st,"%2s", hdr.country) != 1);
                        break;

               case 'P':                    /* get pressure interval */
                        pflag = 1;
                        error = (sscanf(&argv[i][2],"%d", &prsint) != 1);
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

   if (! (oflag && nfiles && pflag && nflag && sflag)) {
       fprintf(stderr,"\nYou must specify input file(s), output file, nation_code, ship_code, and pressure interval.\n");
       exit(1);
   }
   
     /* read in the header file */
 
     if ( !hflag ) {
    
           hdrfile = fopen(HEADERPATH,"r");
           if (hdrfile == NULL) {
                fprintf(stderr,"\nError opening default header file: %s\n", 
                                  HEADERPATH);
                fprintf(stderr,"Specify an alternate file with -H option.\n");
               exit(1);
           }
     }     
     parse_header(hdrfile);
     
     fclose(hdrfile);

/* loop for each input file */

   do {
   
     infile = openfile(dir, argv[curfile], extent);
     if (infile == NULL) goto NEXTFILE;
      
     if (strlen(argv[curfile]) == 5) {
        st = argv[curfile];
     }
     else if ((st = strstr(argv[curfile],".ctd")) != NULL) {
        st -=  5;
     }
     else {
        st = argv[curfile];
        st += strlen(argv[curfile]);
        st -= 5;
     }
     
     sscanf(st, "%3d%2d", &sta, &cast);
     index = sta * 100 + cast;
     
     nprops = 4;   /* assume pr,te,sa,ox is available. */
     hdr.nprops = nprops + 1;    /* add depth to the output vars */
         
              
     
     /* adjust some of the header info ...   */
     
     hdr.station = sta;
     hdr.year = info[index].year;
     hdr.month = info[index].month;
     hdr.day = info[index].day;
     hdr.lat = info[index].lat;
     hdr.lon = info[index].lon;
     hdr.pdr = info[index].depth;
     
     dlat = (double) hdr.lat;
     hdr.ms10 = ms10(hdr.lat, hdr.lon, &hdr.ms1);
     hdr.qual[0] = hdr.qual[1] = '0';
     hdr.qual[2] = hdr.qual[3] = '1';
     hdr.origin = '8';
     hdr.instrument = 'c';
     hdr.prop_id = (int *) malloc(hdr.nprops * sizeof(int));
     hdr.prop_id[0] = (int)PR;
     hdr.prop_id[1] = (int)DE;
     hdr.prop_id[2] = (int)TE;
     hdr.prop_id[3] = (int)SA;
     hdr.prop_id[4] = (int)O2;
     
     /*  allocate space -- here nscans is set to a value at start of main() */
          
     p = (double *) calloc(nscans, sizeof(double));
     d = (double *) calloc(nscans, sizeof(double));
     t = (double *) calloc(nscans, sizeof(double));
     s = (double *) calloc(nscans, sizeof(double));
     o = (double *) calloc(nscans, sizeof(double));
        
    /* read each property */    
        
     nobs = readdata(infile);

     /* decimate the arrays according to specified prsint */
    
     deltap = DELTAP;
     npts = (int) (prsint / deltap);  /*  get # of pts per interval*/
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
       if (o[i] >= 0)
             o[j] = o[i];
       else
             o[j] = HB_MISSING;   
       ++j;
     }
     
     /* add bottom observation  ... */
     
     if ((i - npts) != (--nobs)) {
        p[j] = p[nobs];
        d[j] = hb_depth(p[nobs], dlat);
        t[j] = t[nobs];
        s[j] = s[nobs];
        if (o[nobs] >= 0)
             o[j] = o[nobs];
        else
             o[j] = HB_MISSING;   
        ++j;
     }
     
     hdr.nobs = data.nobs = j;
     data.nprops = hdr.nprops; 
     
     data.observ[(int)PR] = p;   /* set these pointers */
     data.observ[(int)DE] = d;
     data.observ[(int)TE] = t;
     data.observ[(int)SA] = s;
     data.observ[(int)O2] = o;
        
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

   exit(0);
        
} /* end main() */ 
  
/*****************************************************************************/

void print_usage(program)
char *program;
{   
   fprintf(stderr,"\nReads ctd files from Scripps ODF (A20 cruise) and converts to HydroBase format.\n");

   fprintf(stderr,"\nUsage:  %s filelist -Ooutfile -Pprs_int -N<nation_code> -S<ship_code>  [-D<dirname>] [-E<file_extent>]  [-H<header_file>]  [-X]", program);
   fprintf(stderr,"\n\n  List of filenames MUST be first argument!");
   fprintf(stderr,"\n    -O   : specifies output_file ");  
   fprintf(stderr,"\n    -P   : specifies pressure interval ");  
   fprintf(stderr,"\n    -N   : specifies 2char country_code ");  
   fprintf(stderr,"\n    -S   : specifies 2 char shipcode ");  
   fprintf(stderr,"\n   [-D]  : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n   [-E]  : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n   [-H]  : specifiy file of header info.  Default:  %s", HEADERPATH);  

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
   /* reads the header info from an already opened file and fills in an array of records.   */
{
   int error, sta, cast, index, count;
   char line[100], blah[10];

   count = 0;
   while ((error = fscanf(fptr, "%[^\n]", line)) != EOF) {
   
        ++count;
        if (error != 1) {
          fprintf(stderr,"Error reading line #%d from header file.\n", count);
         exit(1);
        }
       
        fgetc(fptr);                                     /* move past LF */
        
       if ( sscanf(line, "%d%d", &sta, &cast) != 2) {
         fprintf(stderr,"Error parsing station and cast from header file:  \n%s \n(line printed above)\n", line);
         exit(1);
       }
       
       index = sta * 100 + cast;


/*     fprintf(stderr,"sta=%d, cast=%d\n", sta, cast);   */
      
       if (sscanf(&line[0], "%d%d%4d-%2d-%2d%9s%f%f%d",  &sta, &cast, &info[index].year, &info[index].month, &info[index].day, &blah, &info[index].lat, &info[index].lon,&info[index].depth) != 9) {
       
        fprintf(stderr,"Error parsing info from header file:  \n%s \n(line printed above)\n", line);
        fprintf(stderr,"index= %d, year=%d, month=%d, blah=%s\n", index, info[index].year, info[index].month, blah);
        
         exit(1);
       }
       
    }   /* end while */
    return;
}  /* end parse_header() */
/*****************************************************************************/
int readdata(FILE *fptr)
{
   int   j, error;
   char line[200];
   
   ngood_ox = 0;
   
   j = 0;     /* line counter */
   while ( (error = fscanf(fptr, "%[^\n]", line)) != EOF) {
        if (error != 1) {
          fprintf(stderr,"Error reading line #%d from data file.\n", j);
         exit(1);
        }
       
        fgetc(fptr);                                     /* move past LF */
        
        
         if ( sscanf(line, "%lf%lf%lf%lf", &p[j], &t[j], &s[j], &o[j]) != 4) {
         fprintf(stderr,"Error parsing pr, te, sa, ox from data file: \n%s \n(line printed above)\n", line);
         exit(1);
     }
       
      if (o[j] > 0) 
         ++ngood_ox;
         
       ++j;
   } /* end while */
   
   return(j);
}  /* end readdata() */
