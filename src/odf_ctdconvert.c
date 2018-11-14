/*  odf_ctdconvert.c
................................................................................
                          *******  HydroBase2 *******
               ...................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
                             updated to ANSI C July 2001
................................................................................
................................................................................
.   Reads odf  format file(s) (from Bedford Inst. of Oceanography)
.   extracts :    header info
.                 p,t,s,ox,  
.   applies a nominal filter to each (5 pts for p,t,s; 11 pts for ox)
.   and outputs p,d,t,s and o at the specified prs_int to the output
.   file.
.
.   USAGE: odf_ctdconvert infile_list -Ooutfile -P<prs_int> [-D<dir>] [-E<extent>]
................................................................................
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "hydrobase.h"
#include "hb_filters.h"


#define   FILTWIDTH  5     /* # of pts in gaussian filter */
#define   FILTOX    11     /* # of pts in gaussian filter */

#define   DIR    ""
#define   EXTENT ""


#define P_POS 1
#define T_POS 0
#define S_POS 3
#define O_POS -9
#define NCTDPROPS 8



  /* globally defined variables */
  
struct HYDRO_HDR  hdr;
struct HYDRO_DATA data;
double *p, *d, *t, *s, *o;
char pr_avail, te_avail, sa_avail, ox_avail;
int p_pos, t_pos, s_pos, o_pos;
int ngood_ox;
int deltap;          /* pressure interval of ctd cast */
int nscans, nctdprops;
int cruiseflag;

  /* prototypes for locally define functions */
  
void print_usage(char *);
FILE *openfile(char *, char *, char *);
int readdata(FILE *);
void parse_header(FILE *);
   

int main (int argc, char **argv)
{
   int error, nobs, nprops;
   int outfile;
   int  i, j, curfile = 1, nfiles = 0; 
   short oflag, pflag, sflag, nflag, xflag; 
   char *dir, *extent, *st;
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
    pflag = oflag =sflag = nflag =  xflag = 0;
    deltap = 2.0;
    cruiseflag = 0;
 
 
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

               case 'X':                    /* preset data positions */
                        xflag = 1;
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

/* loop for each input file */

   do {
   
     infile = openfile(dir, argv[curfile], extent);
     if (infile == NULL) goto NEXTFILE;
      
      /* parse the header info */
      
     parse_header(infile);
     
     nprops = 3;   /* assume pr,te,sa is available. */
     if (ox_avail)   
         ++nprops;
     hdr.nprops = nprops + 1;    /* add depth to the output vars */
         
     if (xflag) {
       p_pos = P_POS;
       t_pos = T_POS;
       s_pos = S_POS;
       o_pos = O_POS;
       nctdprops = NCTDPROPS;
     }
     
     if ( o_pos < 0 )
       ox_avail = 0;
             
     
     /* adjust some of the header info ...   */
     
     dlat = (double) hdr.lat;
     hdr.ms10 = ms10(hdr.lat, hdr.lon, &hdr.ms1);
     for (i = 0; i < NQUAL; ++i)
       hdr.qual[i] = '0';
     hdr.origin = '8';
     hdr.instrument = 'c';
     hdr.prop_id = (int *) malloc(hdr.nprops * sizeof(int));
     hdr.prop_id[0] = (int)PR;
     hdr.prop_id[1] = (int)DE;
     hdr.prop_id[2] = (int)TE;
     hdr.prop_id[3] = (int)SA;
     if (ox_avail) 
        hdr.prop_id[4] = (int)OX;
     
     /*  allocate space */
          
     p = (double *) calloc(nscans, sizeof(double));
     d = (double *) calloc(nscans, sizeof(double));
     t = (double *) calloc(nscans, sizeof(double));
     s = (double *) calloc(nscans, sizeof(double));
     if (ox_avail) 
        o = (double *) calloc(nscans, sizeof(double));
        
    /* read each property */    
        
     nobs = readdata(infile);

/*******  COMMENTED OUT  filter 
     if (te_avail)
        filter('g',t, nobs, FILTWIDTH);
     if (sa_avail)
        filter('g',s, nobs, FILTWIDTH);
     if (ox_avail ) {
        if (ngood_ox == nobs)
          filter('g',o, nobs, FILTOX);
     }   
************/

     /* decimate the arrays according to specified prsint */
    
     deltap = p[2] - p[1];
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
       if (ox_avail) {
          if (o[i] >= 0)
             o[j] = o[i];
          else
             o[j] = HB_MISSING;   
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
             o[j] = o[nobs];
          else
             o[j] = HB_MISSING;   
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

void print_usage(program)
char *program;
{   
   fprintf(stderr,"\nReads odf files from Bedford Inst. of Oceanography and converts to HydroBase format.\n");

   fprintf(stderr,"\nUsage:  %s filelist -Ooutfile -Pprs_int [-D<dirname>] [-E<file_extent>]", program);
   fprintf(stderr,"\n\n  List of filenames MUST be first argument!");
   fprintf(stderr,"\n    -D  : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n    -E  : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n    -O  : specifies output_file ");  
   fprintf(stderr,"\n    -P  : specifies pressure interval ");  
   fprintf(stderr,"\n    -N  : specifies 2char country_code ");  
   fprintf(stderr,"\n    -S  : specifies 2 char shipcode ");  
   fprintf(stderr,"\n    -X  : use preset values for data arrangement");  

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
   /* reads the header info from an already opened ctd data file and returns
      the number of obs.  The info is used to fill in appropriate values in
      the struct HYDRO_HDR hdr, a global variable.  */
{
   int error, i, found, rflag, index;
   char line[256], *st, mon[4], desc[20];

   found = 0;
   rflag = 0;   
   while (!found ){
      error = fscanf(fptr,"%[^\n]", line);
      getc(fptr);
      
      if (error == EOF) {
         if (rflag) {
           fprintf(stderr,"\n Unable to find START_DATE. \n");
           exit(1);   
         }
         rflag = 1;
         rewind(fptr);
      }
      else {
        if (error != 1) {  /* couldn't read in any chars */
          fprintf(stderr,"\n Error attempting to read line:\n%s\n", line);
          exit(1);
        }   
        found = ((st = strstr(line,"START_DATE_TIME")) != NULL); 
      } 
   }  /* end while */ 
      
   st = strstr(line, "=");
   ++st;
   ++st;
   if (sscanf(st,"%d-%3s-%4d", &hdr.day, mon, &hdr.year) != 3) {
      fprintf(stderr,"\nError parsing START_DATE in this line:\n%s\n", line);
      exit(1);   
   }
   switch (*mon) {
      case 'A':
         hdr.month = 4;
         if (mon[1] == 'U') {
           hdr.month = 8;
           break;
        }
        break;
       case 'J':
        if (mon[1] == 'A') {
           hdr.month = 1;
           break;
        }
        if (mon[1] == 'U') {
           hdr.month = 7;
           if (mon[2] == 'N')
              hdr.month = 6;
           break;
        }
        break;
      case 'F':
          hdr.month = 2;
      case 'M':
          hdr.month = 3;
        if (mon[2] == 'Y') {
           hdr.month = 5;
           break;
        }
      case 'S':
          hdr.month = 9;
           break;
      case 'O':
          hdr.month = 10;
           break;
      case 'N':
          hdr.month = 11;
           break;
      case 'D':
          hdr.month = 12;
           break;
      default:
          break;
          
   } /* end switch */
   
      
   found = 0;
   rflag = 0;   
   while (!found ){
      error = fscanf(fptr,"%[^\n]", line);
      getc(fptr);
      
      if (error == EOF) {
         if (rflag) {
           fprintf(stderr,"\n Unable to find INITIAL_LATITUDE. \n");
           exit(1);   
         }
         rflag = 1;
         rewind(fptr);
      }
      else {
        if (error != 1) {  /* couldn't read in any chars */
          fprintf(stderr,"\n Error attempting to read line:\n%s\n", line);
          exit(1);
        }   
        found = ((st = strstr(line,"INITIAL_LATITUDE")) != NULL); 
      } 
    } /* end while */  
   st = strstr( line, "=");
   ++st;
   if (sscanf(st,"%f", &hdr.lat) != 1) {
      fprintf(stderr,"\nError parsing INITIAL_LATITUDE in this line:\n%s\n", line);
      exit(1);   
   }
   
   
   found = 0;
   rflag = 0;   
   while (!found ){
      error = fscanf(fptr,"%[^\n]", line);
      getc(fptr);
      
      if (error == EOF) {
         if (rflag) {
           fprintf(stderr,"\n Unable to find INITIAL_LONGITUDE. \n");
           exit(1);   
         }
         rflag = 1;
         rewind(fptr);
      }
      else {
        if (error != 1) {  /* couldn't read in any chars */
          fprintf(stderr,"\n Error attempting to read line:\n%s\n", line);
          exit(1);
        }   
        found = ((st = strstr(line,"INITIAL_LONGITUDE")) != NULL); 
      } 
    }  
      
   st = strstr(line, "=");
   ++st;
   if (sscanf(st,"%f", &hdr.lon) != 1) {
      fprintf(stderr,"\nError parsing INITIAL_LONGITUDE in this line:\n%s\n", line);
      exit(1);   
   }
 
   found = 0;
   rflag = 0;   
   while (!found ){
      error = fscanf(fptr,"%[^\n]", line);
      getc(fptr);
      
      if (error == EOF) {
         if (rflag) {
           fprintf(stderr,"\n Unable to find SOUNDING. \n");
           exit(1);   
         }
         rflag = 1;
         rewind(fptr);
      }
      else {
        if (error != 1) {  /* couldn't read in any chars */
          fprintf(stderr,"\n Error attempting to read line:\n%s\n", line);
          exit(1);
        }   
        found = ((st = strstr(line, "SOUNDING")) != NULL); 
      } 
    }  
      
   st = strstr(line, "=");
   ++st;
   if (sscanf(st,"%d", &hdr.pdr) != 1) {
      fprintf(stderr,"\nError parsing SOUNDING in this line:\n%s\n", line);
      exit(1);   
   }
 
/* now find first line labelled PROCESS= */
   found = 0;
   rflag = 0;   
   while (!found ){
      error = fscanf(fptr,"%[^\n]", line);
      getc(fptr);
      
      if (error == EOF) {
         if (rflag) {
           fprintf(stderr,"\n Unable to find PROCESS=. \n");
           exit(1);   
         }
         rflag = 1;
         rewind(fptr);
      }
      else {
        if (error != 1) {  /* couldn't read in any chars */
          fprintf(stderr,"\n Error attempting to read line:\n%s\n", line);
          exit(1);
        }   
        found = ((st = strstr(line, "PROCESS=")) != NULL); 
      } 
    } 
    
   if (!cruiseflag) {    
      found = 0;
      rflag = 0;   
      while (!found ){
         error = fscanf(fptr,"%[^\n]", line);
         getc(fptr);
      
         if (error == EOF) {
            if (rflag) {
              fprintf(stderr,"\n Unable to find Cruise. \n");
              exit(1);   
            }
            rflag = 1;
            rewind(fptr);
         }
         else {
           if (error != 1) {  /* couldn't read in any chars */
             fprintf(stderr,"\n Error attempting to read line:\n%s\n", line);
             exit(1);
           }   
           found = ((st = strstr(line, "Cruise")) != NULL); 
         } 
       } 
      st += 7;
      if (sscanf(st,"%d", &hdr.cruise) != 1) {
         fprintf(stderr,"\nError parsing Cruise in this line:\n%s\n", line);
         exit(1);   
      }
   }
      
   found = 0;
   rflag = 0;   
   while (!found ){
      error = fscanf(fptr,"%[^\n]", line);
      getc(fptr);
      
      if (error == EOF) {
         if (rflag) {
           fprintf(stderr,"\n Unable to find Station. \n");
           exit(1);   
         }
         rflag = 1;
         rewind(fptr);
      }
      else {
        if (error != 1) {  /* couldn't read in any chars */
          fprintf(stderr,"\n Error attempting to read line:\n%s\n", line);
          exit(1);
        }   
        found = ((st = strstr(line, "Station")) != NULL); 
      } 
    } 
   st += 8;
   if (sscanf(st,"%d", &hdr.station) != 1) {
      fprintf(stderr,"\nError parsing Station in this line:\n%s\n", line);
      exit(1);   
   }
   
   found = 0;
   rflag = 0;   
   while (!found ){
      error = fscanf(fptr,"%[^\n]", line);
      getc(fptr);
      
      if (error == EOF) {
         if (rflag) {
           fprintf(stderr,"\n Unable to find nquan. \n");
           exit(1);   
         }
         rflag = 1;
         rewind(fptr);
      }
      else {
        if (error != 1) {  /* couldn't read in any chars */
          fprintf(stderr,"\n Error attempting to read line:\n%s\n", line);
          exit(1);
        }   
        found = ((st = strstr(line, "nquan")) != NULL); 
      } 
    } 
   
   while (*st != '=' )
       ++st;
   ++st;
   if (sscanf(st,"%d", &nctdprops) != 1) {
      fprintf(stderr,"\nError parsing nquan in this line:\n%s\n", line);
      exit(1);   
   }
   nctdprops;     /* there is one less than this property says  */
    
   found = 0;
   rflag = 0;   
   while (!found ){
      error = fscanf(fptr,"%[^\n]", line);
      getc(fptr);
      
      if (error == EOF) {
         if (rflag) {
           fprintf(stderr,"\n Unable to find nvalues. \n");
           exit(1);   
         }
         rflag = 1;
         rewind(fptr);
      }
      else {
        if (error != 1) {  /* couldn't read in any chars */
          fprintf(stderr,"\n Error attempting to read line:\n%s\n", line);
          exit(1);
        }   
        found = ((st = strstr(line,"nvalues")) != NULL); 
      } 
    } 
   
   while (*st != '=')
       ++st;
   ++st;
   if (sscanf(st,"%d", &nscans) != 1) {
      fprintf(stderr,"\nError parsing nvalues in this line:\n%s\n", line);
      exit(1);   
   }
  
  /* now find the name label so we can get the order of properties in each scan */
  
   found = 0;
   rflag = 0;   
   while (!found ){
      error = fscanf(fptr,"%[^\n]", line);
      getc(fptr);
      
      if (error == EOF) {
         if (rflag) {
           fprintf(stderr,"\n Unable to find PROCESS= # name. \n");
           exit(1);   
         }
         rflag = 1;
         rewind(fptr);
      }
      else {
        if (error != 1) {  /* couldn't read in any chars */
          fprintf(stderr,"\n Error attempting to read line:\n%s\n", line);
          exit(1);
        }   
        found = ((st = strstr(line, "# name")) != NULL); 
      } 
    }  /*end while */
   
   pr_avail = 0;
   te_avail = 0;
   sa_avail = 0;
   ox_avail = 0;
   
   do {
      st += 6;
      if ((i = sscanf(st,"%d = %s", &index, desc)) != 2) {
        fprintf(stderr,"\nError parsing property descriptor:\n%s\n", line);
        exit(1);
      }
      
      switch (desc[0]) {
        case 'p':
            if (strncmp(desc,"pr:",3) == 0) {
               pr_avail = 1;
               p_pos = index;
            }
            break;
        case 't':
            if (strncmp(desc,"t068:",5) == 0) {
               te_avail = 1;
               t_pos = index;
            }
            break;
        case 's':
            if (strncmp(desc,"sal00",5) == 0) {
               sa_avail = 1;
               s_pos = index;
            }
            break;
        case 'o':
            if (strncmp(desc,"oxML/L",6) == 0) {
               ox_avail = 1;
               o_pos = index;
            }
            break;
        default:
            /*get next line */
            break;
      } /* end switch */
   
      error = fscanf(fptr,"%[^\n]", line);
      getc(fptr);
      
   } while ((st = strstr(line, "# name")) != NULL);
   
    if (!pr_avail)
       fprintf(stderr,"\nWARNING:  no pr in this file.");
    if (!te_avail)
       fprintf(stderr,"\nWARNING:  no t068 in this file.");
    if (!sa_avail)
       fprintf(stderr,"\nWARNING:  no sal00 in this file.");
      
   found = 0;
   rflag = 0;   
   while (!found ){
      error = fscanf(fptr,"%[^\n]", line);
      getc(fptr);
      
      if (error == EOF) {
         if (rflag) {
           fprintf(stderr,"\n Unable to find -- DATA. \n");
           exit(1);   
         }
         rflag = 1;
         rewind(fptr);
      }
      else {
        if (error != 1) {  /* couldn't read in any chars */
          fprintf(stderr,"\n Error attempting to read line:\n%s\n", line);
          exit(1);
        }   
        found = ((st = strstr(line, "-- DATA")) != NULL); 
      } 
    } 
    /*  file should now be positioned at start of data */

   return;
}  /* end parse_header() */
/*****************************************************************************/
int readdata(FILE *fptr)
{
   int nobs, i, j, nread;
   char line[512];
   double buffer[100];
   
   nobs = 0;
   ngood_ox = 0;
   
   j = 0;     /* line counter */
   while (1) {
      for (i = 0; i < nctdprops; ++i) {
         nread = fscanf(fptr, "%lf", &buffer[i]);
         if (nread == EOF){                        /* terminates loop */
            if (i > 0) {
              fprintf(stderr,"\nEOF but datascan was incomplete: \n expecting %d values in scan, read %d values, %d scans were previously read.\n", nctdprops, i-1, nobs);
              exit(1);
            }
            return (nobs);
         }
         
         if (nread != 1){
            fprintf(stderr,"\nError reading %dth prop at data line#%d\n", i, j);
            exit(1);
         }
      } /* end for */
      
      p[nobs] = buffer[p_pos];
      if (te_avail) {
         t[nobs] = buffer[t_pos];
      } 

      if (sa_avail){ 
         s[nobs] = buffer[s_pos];
      }  
      if (ox_avail) {
         o[nobs] = buffer[o_pos];
      }
      
      
      
      if ((p[nobs] > -8) && (t[nobs] > -8) && (s[nobs] > -8) ) {

         if (nobs > 0) {
            if (p[nobs] > p[nobs-1]) {
              if (o[nobs] > 0)  /* check this before ++ nobs */
                 ++ngood_ox;
              ++nobs;
            }
         }
         else {                 /* first scan */
              if (o[nobs] > 0)
                 ++ngood_ox;
          ++nobs;              /* check this before ++ nobs */
         }
      }
      ++j;     
      
   } /* end while */
}  /* end readdata() */
