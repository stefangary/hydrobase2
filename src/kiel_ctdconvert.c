/*  bio_ctdconvert.c
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
.   USAGE: bio_ctdconvert infile_list -Ooutfile -P<prs_int> [-D<dir>] [-E<extent>] -Ncountry_code -Sship_code -Ccruise_code
................................................................................
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "hydrobase.h"

#define   DIR    ""
#define   EXTENT ""


  /* globally defined variables */
  
struct HYDRO_HDR  hdr;
struct HYDRO_DATA data;
double *p, *d, *t, *s, *o;
char  ox_avail;
int p_pos, t_pos, s_pos, o_pos;
int ngood_ox;
int deltap;          /* pressure interval of ctd cast */
int nscans, nctdprops;
int cruiseflag, lflag;

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
   short oflag,  sflag, nflag, xflag; 
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
    oflag = lflag = sflag = nflag = xflag = 0;
    cruiseflag = 0;
    prsint = 10.0;  /* output pressure interval */
    nscans = 6000;  /* for dimensioning arrays */
 
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

               case 'L':                    /* lat/lon format is degHmin */
                        lflag = 1;
                        break;
               case 'S':                    /* get ship code */
                        sflag = 1;
                        st = &argv[i][2];
                        if (*st == '/' )
                          ++st;
                        error = (sscanf(st,"%2s", hdr.ship) != 1);
                        break;

               case 'X':                    /* specify data columns: p, t, s, ox */
			xflag = 1;
                        error = (sscanf(&argv[i][2],"%d/%d/%d/%d", &p_pos, &t_pos, &s_pos, &o_pos) < 4);
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

       else  {
           ++nfiles;
       }
   }  /* end for */

   if (! (oflag && nfiles  && nflag && sflag && cruiseflag)) {
       fprintf(stderr,"\nYou must specify input file(s), output file, nation_code, ship_code, and cruise id.\n");
       exit(1);
   }

   if (!xflag) {
       fprintf(stderr,"\nSpecify the column of input variables pr, te, sa, ox -X0/1/3/5 \n");
       fprintf(stderr,"  Use -1 if no oxygen \n");
       exit(1);
   }
   if (lflag) 
       fprintf(stderr,"\nInput positions are formatted like 45N66.3 \n");
   else
       fprintf(stderr,"\nInput positions are formatted like 45 66.3 N \n");
   
              

/* loop for each input file */

   do {
   
     infile = openfile(dir, argv[curfile], extent);
     if (infile == NULL) goto NEXTFILE;
      
      /* parse the header info */
      
     parse_header(infile);
     ox_avail = 1;
     if ( o_pos < 0 )
       ox_avail = 0;
     
     nprops = 3;   /* assume pr,te,sa is available. */
     if (ox_avail)   
         ++nprops;
     hdr.nprops = nprops + 1;    /* add depth to the output vars */
         
     
     /* adjust some of the header info ...   */
     
     dlat = (double) hdr.lat;
     hdr.ms10 = ms10(hdr.lat, hdr.lon, &hdr.ms1);
     for (i = 0; i < NQUAL; ++i)
       hdr.qual[i] = '0';
     hdr.origin = '9';
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
     
     if (hdr.pdr <= 0 )
        hdr.pdr = d[j-1] + 10.0;
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
   fprintf(stderr,"\nReads odf files from Kiel IfM and converts to HydroBase format.\n");

   fprintf(stderr,"\nUsage:  %s filelist -Ooutfile -Pprs_int [-D<dirname>] [-E<file_extent>]", program);
   fprintf(stderr,"\n    -C  : specifies cruise id ");  
   fprintf(stderr,"\n    -N  : specifies 2-char country_code ");  
   fprintf(stderr,"\n    -L  : lat/lon is in format like 45N36.0 (default is 45 36.0 N) ");  
   fprintf(stderr,"\n    -O  : specifies output_file ");  
   fprintf(stderr,"\n    -P  : specifies pressure interval ");  
   fprintf(stderr,"\n    -S  : specifies 2 char shipcode ");  
   fprintf(stderr,"\n    -X  : specifies column position for pr/te/sa/ox arrangement [-X0/1/3/-1]");  
   fprintf(stderr,"\n    -D  : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n    -E  : specifies input_file extent (default is no extent)");  
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
void parse_header(FILE *fptr)
   /* reads the header info from an already opened ctd data file and returns
      the number of obs.  The info is used to fill in appropriate values in
      the struct HYDRO_HDR hdr, a global variable.  */
{
   int error, i, found, rflag, index;
   float deg, min;
   char hem;
   char line[256], *st;

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
      
   st = strstr(line, "=");
   ++st;
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
           fprintf(stderr,"\n Unable to find Date. \n");
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
        found = ((st = strstr(line,"Date")) != NULL); 
      } 
   }  /* end while */ 
      
   st = strstr(line, "=");
   ++st;
   if (sscanf(st,"%d/%d/%d", &hdr.year, &hdr.month, &hdr.day) != 3) {
      fprintf(stderr,"\nError parsing Date in this line:\n%s\n", line);
      exit(1);   
   }
   found = 0;
   rflag = 0;   
   while (!found ){
      error = fscanf(fptr,"%[^\n]", line);
      getc(fptr);
      
      if (error == EOF) {
         if (rflag) {
           fprintf(stderr,"\n Unable to find Lat. \n");
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
        found = ((st = strstr(line,"Lat")) != NULL); 
      } 
    } /* end while */  
   st = strstr( line, "=");
   ++st;
   
   
   if (lflag) {
      if (sscanf(st,"%f%c%f", &deg, &hem, &min ) != 3) {
         fprintf(stderr,"\nError parsing Lat in this line:\n%s\n", line);
         exit(1);   
      }
   
   }
   else {
      if (sscanf(st,"%f %f %c", &deg, &min, &hem) != 3) {
         fprintf(stderr,"\nError parsing Lat in this line:\n%s\n", line);
         exit(1);   
      }
   }
   
   hdr.lat = deg + min / 60.0;
   if ((hem == 'S') || (hem == 's'))
       hdr.lat = -hdr.lat;
   
   found = 0;
   rflag = 0;   
   while (!found ){
      error = fscanf(fptr,"%[^\n]", line);
      getc(fptr);
      
      if (error == EOF) {
         if (rflag) {
           fprintf(stderr,"\n Unable to find Lon. \n");
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
        found = ((st = strstr(line,"Lon")) != NULL); 
      } 
    }  
      
   st = strstr(line, "=");
   ++st;
   
   if (lflag) {
      if (sscanf(st,"%f%c%f", &deg, &hem, &min ) != 3) {
         fprintf(stderr,"\nError parsing Lon in this line:\n%s\n", line);
         exit(1);   
      }
   
   }
   else {
      if (sscanf(st,"%f %f %c", &deg, &min, &hem) != 3) {
         fprintf(stderr,"\nError parsing Lon in this line:\n%s\n", line);
         exit(1);   
      }
   }
   
   hdr.lon = deg + min/60.0;
   if ((hem == 'W') || (hem == 'w'))
      hdr.lon = -hdr.lon;
 
   found = 0;
   rflag = 0;   
   while (!found ){
      error = fscanf(fptr,"%[^\n]", line);
      getc(fptr);
      
      if (error == EOF) {
         if (rflag) {
           fprintf(stderr,"\n Unable to find WaterDepth. \n");
	   hdr.pdr = 0;
           goto BYPASS;  
         }
         rflag = 1;
         rewind(fptr);
      }
      else {
        if (error != 1) {  /* couldn't read in any chars */
          fprintf(stderr,"\n Error attempting to read line:\n%s\n", line);
          exit(1);
        }   
        found = ((st = strstr(line, "WaterDepth")) != NULL); 
      } 
    }  
      
   st = strstr(line, "=");
   ++st;
   if (sscanf(st,"%d", &hdr.pdr) != 1) {
      fprintf(stderr,"\nError parsing WaterDepth in this line:\n%s\n", line);
      exit(1);   
   }

BYPASS: 
   found = 0;
   rflag = 0;   
   while (!found ){
      error = fscanf(fptr,"%[^\n]", line);
      getc(fptr);
      
      if (error == EOF) {
         if (rflag) {
           fprintf(stderr,"\n Unable to find Columns. \n");
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
        found = ((st = strstr(line, "Columns")) != NULL); 
      } 
    } 
    st = strstr(line, "=");
   ++st;

   nctdprops = 1;   
   while ((st = strstr(st, ":")) != NULL ) {
      ++nctdprops;
      ++st;
   }
   fprintf(stderr,"Station %d has %d properties\n", hdr.station, nctdprops);
      
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
      t[nobs] = buffer[t_pos];
      s[nobs] = buffer[s_pos];
      if (ox_avail) {
         o[nobs] = buffer[o_pos];
      }
      
      if ((p[nobs] > -8) && (t[nobs] > -8) && (s[nobs] > -8) ) {

         if (nobs > 0) {
            if (p[nobs] > p[nobs-1]) {
	       if (ox_avail) {
                 if (o[nobs] > 0)  /* check this before ++ nobs */
                    ++ngood_ox;
	       }
              ++nobs;
            }
         }
         else {                 /* first scan */
	   if (ox_avail) {
              if (o[nobs] > 0)
                 ++ngood_ox;
	   }
          ++nobs;              /* check this before ++ nobs */
         }
      }
      ++j;     
      
   } /* end while */
}  /* end readdata() */
