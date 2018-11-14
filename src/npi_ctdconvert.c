/*  npi_ctdconvert.c
................................................................................
                          *******  HydroBase *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             1993
................................................................................
................................................................................
.   Reads Norwegian Polar Institute  ctd  files from ship and helicopter(SBE9)
.   extracts :    header info
.                 p,t,s (o) 
.   and outputs properties at the specified prs_int to the output
.   file.
.
.   USAGE: npi_ctdconvert infile_list -Ooutfile -P<prs_int> -N<country_code> -C<cruise> -S<ship> -X<ox_code> [-D<dir>] [-E<extent>]
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
int nprops_per_scan;
int p_pos, t_pos, s_pos, o_pos;  /* position of parameters in each scan */
int deltap;          /* pressure interval of ctd cast */
int temp90;
int ox_code;          /* 0=no oxygen; 1=ml/l; 2=micromoles/kg */
int helicopter;      /* flag ctd format for helicopter stations*/
int is_downcast;
int binavg;

   /* prototypes for locally define functions */

void print_usage(char *);
int readdata(FILE *);
int avgdata(FILE *);
void parse_header(FILE *);
FILE *openfile(char *, char *, char *);
int get_scan(FILE *, double *);
void reverse_order(double *, int);

int main (int argc, char **argv)
{
   int error, nobs, nprops;
   int outfile;
   int  i, j, curfile = 1, nfiles = 0; 
   short oflag, pflag, nflag, sflag, cflag; 
   char *dir, *extent;
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
    cflag = nflag = sflag = 0;
    binavg = 0;
    ox_code = 2;  /* default is micromoles/kg*/
    deltap = 1.0;
    helicopter = 0;
 
 
/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'B':                    /* get cruise  */
                        binavg = 1;
                        break;
               case 'C':                    /* get cruise  */
                        cflag = 1;
                        error = sscanf(&argv[i][2], "%d", &hdr.cruise) != 1;
                        break;
               case 'D':                   /* get input dir */
                        dir = &argv[i][2];
                        break;
               case 'E':                    /* get file extent */
                        extent = &argv[i][2];
                        break;

               case 'H':                    /* helicopter ctd format  */
                        helicopter = 1;
                        break;
               case 'N':                    /* get country code  */
                        nflag = 1;
                        strncpy(hdr.country, &argv[i][2],2);
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
               case 'S':                    /* get ship  */
                        sflag = 1;
                        strncpy(hdr.ship, &argv[i][2],2);
                        break;

               case 'X':                    /*  oxygen flag  */
	                error = (sscanf(&argv[i][2], "%d", &ox_code) != 1);
                        break;
	       case 'h':
	                print_usage(argv[0]);	
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

   if (! (oflag && nfiles && pflag && cflag && sflag && nflag)) {
       fprintf(stderr,"\nYou must specify input file(s), output file, pressure interval, country_code, ship, and cruise\n");
       exit(1);
   }
   
   if (helicopter)
       fprintf(stderr,"\nFiles are helicopter format...");

   switch (ox_code) {
       case 0:
            fprintf(stderr,"\nNo oxygen data acquired");
            break;
       case 1:
            fprintf(stderr,"\nOxygen in ml/liter");
            break;
       case 2:
            fprintf(stderr,"\nOxygen in micromoles/kg");
            break;
       default:
             ;
       
   } /* end switch */   

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
     
     for (i = 0; i < NQUAL; ++i)
       hdr.qual[i] = '0';
     hdr.origin = '5';
     hdr.instrument = 'c';
     hdr.ms10 = ms10(hdr.lat, hdr.lon, &hdr.ms1);
     hdr.prop_id = (int *) malloc(hdr.nprops * sizeof(int));
     hdr.prop_id[0] = (int)PR;
     hdr.prop_id[1] = (int)DE;
     hdr.prop_id[2] = (int)TE;
     hdr.prop_id[3] = (int)SA;
     switch (ox_code) {
        case 1: 
           hdr.prop_id[4] = (int)OX;
	   break;
        case 2: 
           hdr.prop_id[4] = (int)O2;
	   break;
	default:
	   ;   /* no oxygen */
     } /* end switch */
     
     /*  allocate space */
     
        
    /* read each property */  
      
     nobs = 7000;
     p = (double *) calloc(nobs, sizeof(double));
     d = (double *) calloc(nobs, sizeof(double));
     t = (double *) calloc(nobs, sizeof(double));
     s = (double *) calloc(nobs, sizeof(double));
     if (ox_avail) 
        o = (double *) calloc(nobs, sizeof(double));
	
     if (binavg)
        nobs = avgdata(infile);	
     else   
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
       d[j] = hb_depth(p[i], (double) hdr.lat);
       t[j] = t[i];
       s[j] = s[i];
       if (ox_avail) {
          if (o[i] >= 0)
             o[j] = o[i];
          else
             o[j] = MISSING;   
       }
       ++j;
     }
     
     /* add bottom observation  ... */
     
     if ((i - npts) != (--nobs)) {
        p[j] = p[nobs];
        d[j] = hb_depth(p[nobs], (double)hdr.lat);
        t[j] = t[nobs];
        s[j] = s[nobs];
        if (ox_avail) {
          if (o[nobs] >= 0)
             o[j] = o[nobs];
          else
             o[j] = MISSING; 
        }  
        ++j;
     }
     
     hdr.pdr = d[j-1]+ 10.0;
     
     hdr.nobs = data.nobs = j;
     data.nprops = hdr.nprops; 
     
     data.observ[(int)PR] = p;   /* set these pointers */
     data.observ[(int)DE] = d;
     data.observ[(int)TE] = t;
     data.observ[(int)SA] = s;
     if (ox_avail)
        data.observ[(int)O2] = o;
	if (ox_code = 1)
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
   fprintf(stderr,"\nUsage:  %s filelist  -C<cruise_id> -N<country_code> -S<ship_code> -Ooutfile -Pprs_int [-B] [-D<dirname>] [-E<file_extent>] [-X<ox_code>]", program);
   fprintf(stderr,"\n\n  List of filenames MUST be first argument!");
   fprintf(stderr,"\n    -B  : Bin average the properties to form 1 db pressure series.  ");
   fprintf(stderr,"\n          The data file should be edited beforehand to limit the data to either the downcast or the upcast.  Otherwise, both will be included in the averaging.");
   fprintf(stderr,"\n    -C  : cruise number (5-char max) ");
   fprintf(stderr,"\n    -D  : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n    -E  : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n    -H  : file is helicopter ctd format  ");
   fprintf(stderr,"\n    -N  : country code (2-char)  ");
   fprintf(stderr,"\n    -O  : output filename");
   fprintf(stderr,"\n    -P  : specifies pressure interval for output file");  
   fprintf(stderr,"\n            ex: -P10 ");
   fprintf(stderr,"\n    -S  : ship code (2-char) ");
   fprintf(stderr,"\n    -X  : oxygen code: [default = 2]   ");
   fprintf(stderr,"\n            0 = no oxygen ");
   fprintf(stderr,"\n            1 =  oxygen  (ml/l)");
   fprintf(stderr,"\n            2 =  oxygen  (micromoles/kg)");

   fprintf(stderr,"\n\n");  
   return;
}
   
/*****************************************************************************/
FILE *openfile(char *dir, char *root, char *extent)
{
   char st[1000];
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
   int error, i, pos, eof;
   float deg, min;
   char line[2000], *st;
   char  hem;
   char month[4];
   
   if (fscanf(fptr,"%[^\n]", line) != 1) {
        fprintf(stderr,"\n Error attempting to read 1st header line. \n");
        exit(1);   
   } 
   error = getc(fptr);
   
/* Latitude */ 
   
   while ((st = strstr(line, "Lat")) == NULL && !(eof = feof(fptr)) ) {
      fscanf(fptr,"%[^\n]", line);
      error = getc(fptr);
   }
   
   if (eof) {
        fprintf(stderr,"\n Unable to find keyword 'Lat'. \n");
        exit(1);   
   }
   
   if (helicopter) {
      st = strchr(line, (int)':');
      if (st != NULL) {
         ++st;
         error = (sscanf(st, "%f %c", &hdr.lat, &hem) != 2);
      }
   }
   else {
      st = strchr(line, (int)'=');
      if (st == NULL) 
        st = strchr(line, (int)':');
      ++st;
      error = (sscanf(st, "%f %f %c", &deg, &min, &hem) != 3);
      hdr.lat =  deg + min / 60.0;  
   }

   if (error) {
        fprintf(stderr,"\n Error parsing Latitude from: %s \n", line);
        exit(1);   
   }
   
   if (hem == 'S' || hem == 's' )
       hdr.lat = -hdr.lat;
       
       
   rewind(fptr); 

   
 /* Longitude */      
   fscanf(fptr,"%[^\n]", line);
   error = getc(fptr);
   
   while ((st = strstr(line, "Lon")) == NULL && !(eof = feof(fptr)) ) {
      fscanf(fptr,"%[^\n]", line);
      error = getc(fptr);
   }
   
   if (eof) {
        fprintf(stderr,"\n Unable to find keyword 'Lon'. \n");
        exit(1);   
   }
   
   if (helicopter) {
      st = strchr(line, (int)':');
      if (st != NULL) {
         ++st;
         error = (sscanf(st, "%f %c", &hdr.lon, &hem) != 2);
      }
   }
   else {
      st = strchr(line, (int)'=');
      if (st == NULL) 
        st = strchr(line, (int)':');
      ++st;
      error = (sscanf(st, "%f %f %c", &deg, &min, &hem) != 3);
      hdr.lon =  deg + min / 60.0;  
   }

   if (error) {
        fprintf(stderr,"\n Error parsing Longitude from: %s \n", line);
        exit(1);   
   }
   
   if (hem == 'W' || hem == 'w')
      hdr.lon = -hdr.lon;
      
   rewind(fptr); 

   
 /* Date */      
 
   fscanf(fptr,"%[^\n]", line);
   error = getc(fptr);
   
   while ((st = strstr(line, "start_time")) == NULL && !(eof = feof(fptr)) ) {
      fscanf(fptr,"%[^\n]", line);
      error = getc(fptr);
   }
   
   if (eof) {
        fprintf(stderr,"\n Unable to find keyword 'start_time'. \n");
        exit(1);   
   }
   
   st = strchr(line, (int)'=');
   ++st;
   error = (sscanf(st, "%s %d %d", &month, &hdr.day, &hdr.year) != 3);
   if (error) {
        fprintf(stderr,"\n Error parsing start_time. \n");
        exit(1);   
   }

   switch (*month) {
      case 'A':
      case 'a':
          switch (*(month+1)) {
	      case 'P':
	      case 'p':
	           hdr.month = 4;
	           break;
	      case 'U':
	      case 'u':
	           hdr.month = 8;
	           break;
	      default:
	           fprintf(stderr,"\nUnrecognized month: %s \n", month);
		   exit(1);
	  }
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
          switch (*(month+1)) {
	      case 'A':
	      case 'a':
	           hdr.month = 1;
	           break;
	      case 'U':
	      case 'u':
	           hdr.month = 7;
	           if (*(month+2) == 'n')
	               hdr.month = 6;
	           break;
	      default:
	           fprintf(stderr,"\nUnrecognized month: %s \n", month);
		   exit(1);
	  }
	  break;
      case 'M':
      case 'm':
          switch (*(month+2)) {
	      case 'R':
	      case 'r':
	           hdr.month = 3;
	           break;
	      case 'Y':
	      case 'y':
	           hdr.month = 5;
	           break;
	      default:
	           fprintf(stderr,"\nUnrecognized month: %s \n", month);
		   exit(1);
	  }
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
	  fprintf(stderr,"\nUnrecognized month: %s \n", month);
          exit(1);
   } /* end switch */
   
   rewind(fptr); 
   
   
   
 /* Station # */      
   fscanf(fptr,"%[^\n]", line);
   error = getc(fptr);
   
   while ((st = strstr(line, "Station")) == NULL && !(eof = feof(fptr)) ) {
      fscanf(fptr,"%[^\n]", line);
      error = getc(fptr);
   }
   
   if (eof) {
        fprintf(stderr,"\n Unable to find keyword 'Station'. \n");
        exit(1);   
   }
   
   st = strchr(line, (int)':');
   ++st;
   
   if (helicopter) {
      while ( !isdigit((int)*st) )
         ++st;
   }
   
   error = (sscanf(st, "%d", &hdr.station) != 1);
   if (error) {
        fprintf(stderr,"\n Error parsing Station number. \n");
        exit(1);   
   }
         
   /* search for variable descriptions...*/ 
   
   while ((st = strstr(line, "nquan")) == NULL  & !(eof = feof(fptr)) ) {
        fscanf(fptr,"%[^\n]", line);
        error = getc(fptr);
   }
   if (eof) {
        fprintf(stderr,"\n Unable to find keyword 'nquan'. \n");
        exit(1);   
   }
   st = strchr(line, (int) '=');
   ++st;
   error = (sscanf(st, "%d", &nprops_per_scan) != 1);
   if (error) {
        fprintf(stderr,"\n Error parsing number of props per scan: 'nquan'. \n");
        exit(1);   
   }
   
   pr_avail = te_avail = sa_avail = ox_avail = 0;
   
   fscanf(fptr,"%[^\n]", line);
   error = getc(fptr);
   
   while ((st = strstr(line, "name")) == NULL) {
        fscanf(fptr,"%[^\n]", line);
        error = getc(fptr);
   }
   
   do {
      error = sscanf(st+4,"%d", &pos) != 1;
      if (error) {
        fprintf(stderr,"\n Error parsing variable position from: %s \n", line);
        exit(1);   
      } 
   
      if (strstr(line, "Pres") != NULL && !pr_avail) {
         p_pos = pos;
	 pr_avail = 1;
      }
      if (strstr(line, "Temp") != NULL && !te_avail) {
         t_pos = pos;
	 te_avail = 1;
  
         temp90 = 1;
         if (strstr(line,"90") == NULL) {
             if (strstr(line,"68") == NULL)
                fprintf(stderr,"\n Cannot determine whether temperature is ITS 90 or ITS 68 \n");
             else
	        temp90 = 0;
         }
       
      }
      
      if (strstr(line, "Sal") != NULL && !sa_avail) {
         s_pos = pos;
	 sa_avail = 1;
      }
      
      if (strstr(line, "Oxy") != NULL && !ox_avail) {
         o_pos = pos;
	 ox_avail = 1;
         ox_code = 2;
         if (strstr(line,"kg") == NULL) {
	    ox_code = 1;
	 }   
      }
      
      fscanf(fptr,"%[^\n]", line);
      error = getc(fptr);
   } while ((st = strstr(line, "name")) != NULL);
   
   if (!pr_avail)
       fprintf(stderr,"\nWARNING:  no PR in this file.");
   if (!te_avail)
       fprintf(stderr,"\nWARNING:  no TE in this file.");
   if (!sa_avail)
       fprintf(stderr,"\nWARNING:  no SA in this file.");
       
   while (strstr(line, "END") == NULL) {
        fscanf(fptr,"%[^\n]", line);
        error = getc(fptr);
   }
   
 /* File is now positioned at start of data. */
 
   return;
}  /* end parse_header() */
/*****************************************************************************/
int avgdata(FILE *fptr)
{
   int nobs, i, error,  nread, index;
   int *count, *ocount;
   char  *st, line[1000];
   double *buffer;
   double maxp;

   buffer = (double *) calloc(nprops_per_scan, sizeof(double));
   count = (int *) calloc(7000, sizeof(int));
   if (ox_avail)
      ocount = (int *) calloc(7000, sizeof(int));

   /* Cycle through all scans, summing properties within pressure bins */
   maxp = 0;
   while ((nread = get_scan(fptr, buffer)) != EOF) {
   
      if (nread != 1) {
         fprintf(stderr,"\nError reading data scan #%d\n", nobs);
	 exit(1);
      }
      
      index = NINT(buffer[p_pos]);
      if (index >= 0) {
      
          if (index > maxp)
	     maxp = (double) index;
	     
          if (buffer[t_pos] > -2.1 && buffer[s_pos] >= 0) {
	  
                 if (te_avail ) {
                    t[index] += buffer[t_pos];
		    ++count[index];
                 } 

                 if (sa_avail){ 
                    s[index] += buffer[s_pos];
                 }  
                 if (ox_avail && (buffer[o_pos] >= 0)) {
                    o[index] += buffer[o_pos];
		    ++ocount[index];
                 }
         }
     }
   }  /* end while */
   
   
   for (index = 0; index < maxp; ++index) { 
   
      p[index] = (double) index;
   
      if (count[index] > 0)  {
          t[index] /= count[index];
 	 s[index] /= count[index];
	 
	 if (ox_avail && ocount[index] > 0)
	   o[index] /= ocount[index];
	   
      } /* end if count[index] */
   }  /* end for */
   
   
   nobs = 0;
   for (i = 0; i < maxp; ++i ) {
       if (count[i] > 0) {
          p[nobs] = p[i];
	  t[nobs] = t[i];
          if (temp90) 
             t[nobs] *= 1.00024;
	  s[nobs] = s[i];
	  if (ox_avail) {
	     o[nobs] = o[i];
	  } 
	  
	  ++nobs;
       
       }
   }

   free(buffer);
   free(count);
   if (ox_avail)
      free(ocount);
   
   	 
   return (nobs);
}  /* end avgdata() */


/*****************************************************************************/
int readdata(FILE *fptr)
{
   int nobs, i, error,  nread;
   char  *st, line[1000];
   double *buffer;
   double last_pr;

   buffer = (double *) calloc(nprops_per_scan, sizeof(double));


   /* Read first 2 lines of data to determine 
      if pressure is increasing or decreasing  */
      
 /* Data Line 1 */        
   if ( (error = get_scan(fptr, buffer)) != 1){ 
      if (error  == EOF) {
         fprintf(stderr, "\nNo data in this file.\n");
         return(0);
      }
      
      fprintf(stderr, "\nError reading data line 1.\n");
      exit(1);
   }
   
   last_pr = buffer[p_pos];
     

 /* Data Line 2 */   

    is_downcast = 0;
    
    if (get_scan(fptr, buffer) == 1) { 
      if ((buffer[p_pos] - last_pr) > 0) {
         is_downcast = 1;
	 last_pr = -1;
      }
      else
         last_pr = 99999;
   }
    
   rewind(fptr);
   do { 
        fscanf(fptr,"%[^\n]", line);
        error = getc(fptr);
   }  while (strstr(line, "END") == NULL); 
   
 /***** File is now re-positioned at start of data. */
 
   nobs = 0;
   while ((nread = get_scan(fptr, buffer)) != EOF) {
   
      if (nread != 1) {
         fprintf(stderr,"\nError reading data scan #%d\n", nobs);
	 exit(1);
      }
      
      if  ( (is_downcast && (buffer[p_pos] > last_pr))
         || (! is_downcast && (buffer[p_pos] < last_pr)) ) {
           
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
      
         if ((p[nobs] > 0) && (t[nobs] > -8) && (s[nobs] > 0) ) {
            if (temp90) {
               t[nobs] *= 1.00024;
	    }
	    
	    last_pr = p[nobs];
            ++nobs;
         }
      }
   }  /* end while */
   
   free(buffer);
   
   
   if (!is_downcast) { 
      reverse_order(p, nobs);
      reverse_order(s, nobs);
      reverse_order(t, nobs);
      if (ox_avail)
         reverse_order(o, nobs);
   }
   	 
   return (nobs);
}  /* end readdata() */

/*****************************************************************************/
int get_scan(FILE *fptr, double *buffer)

{
   int status, eoscan, i;
   char line[1000], *st;
   
   eoscan = nprops_per_scan - 1;
         
   status = fscanf(fptr, "%[^\n]", line);
   getc(fptr);
   if (status != 1) 
          return (status);
   
   st = &line[0];
   for (i = 0; i < nprops_per_scan; ++i) {
         if (sscanf(st,"%lf",&buffer[i]) != 1) {
            fprintf(stderr,"\nError parsing property %2d at data line: \n%s\n", i, line);
            exit(1);
	 }
	 if (i < eoscan) {
	    while (*st == ' ')  /* first read white space between columns */
	       ++st;
	    while (*st != ' ')  /* then read column */
	       ++st;
	 }
   } /* end for */
   
   return(status);
} /* end get_scan() */

/*****************************************************************************/
void reverse_order(double *x, int n)
/* Reverse the order of array elements in x */
{
   int i;
   double *xtmp;
   
   xtmp = (double *) calloc(n, sizeof(double));
   
   
   for (i = 0; i < n; ++i) {
      xtmp[i] = x[i];
   }
   for (i = 0; i < n; ++i) {
      x[i] = xtmp[n-i-1];
   }
   
   free(xtmp);
   return;

}  /* end reverse_order */

