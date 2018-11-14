/* seabird_btl_convert.c
................................................................................
.   Reads bottle files of observed level hydrographic data files
.   extracts :    header info
.                 p,t,s,ox,n2,n3,si,p4
.   
................................................................................
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "hydrobase.h"

#define   BLANK     0x20   /* ascii code for blank */
#define   MISSING   -9.0   /* missing value flag */
#define   BUFSIZE   512   /* buffer for a line of data */
#define   DIR    ""
#define   EXTENT ".nut"
#define   NPROPS    6   /* pr,de,te,sa,ox,si*/
#define   MAXOBS    50



  /* globally defined variables */
  
struct HYDRO_HDR  hdr;
struct HYDRO_DATA data;
double *p, *d, *t, *o, *s, *si;
double *p4, *n2, *n3;   
int primary_t;      /* flag to use primary temperature sensor */     

void print_usage(char *);
int read_data(FILE *);
void read_hdr(FILE *);


main ( int argc, char **argv)
{
   int status, error;
   int outfile, nobs;
   int staread, staout;
   int  i, j, curfile = 1, nfiles = 0;
   short staOK;
   short sflag, cflag, nflag, oflag, pflag; 
   char *dir, *extent, *st, *prefix; 
   char  fname[100];
   double dlat;
   int  npts;
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
    staout = staread = 0;
    outfile = STDOUT;
    cflag = sflag = nflag = oflag = pflag = 0;
    primary_t =1;
  
    for (i = 0; i < MAXPROP; ++i) {
       data.observ[i] = NULL;
    }

/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
              case 'C':                    /* get cruise  */
                          cflag = 1;
                      error = (sscanf(&argv[i][2],"%d", &hdr.cruise) != 1);
                        break;
              case 'D':                   /* get input dir */
                        dir = &argv[i][2];
                        break;
              case 'E':                    /* get file extent */
                        extent = &argv[i][2];
                        break;
              case 'N':                    /* get nation code  */
                       nflag = 1;
                       error = (sscanf(&argv[i][2],"%s", &hdr.country) != 1);
                        break;
              case 'O':                    /* get output file  */
                        outfile = create_hydro_file(&argv[i][2], OVERWRITE);
                        if (outfile < 1) {
                           fprintf(stderr,"\nError opening output file: %s\n", 
                                  &argv[i][2]);
                           exit(1);
                        }
                        break;
              case 'P':                    /* get file_prefix for station_filenames */
                        pflag = 1;
                        prefix = &argv[i][2];
                        break;
              case 'S':                    /* get ship  */
                       sflag = 1;
                       error = (sscanf(&argv[i][2],"%s", &hdr.ship) != 1);
                        break;

              case 'T':
                       error = (sscanf(&argv[i][2], "%d", &primary_t) != 1);
                       if (primary_t > 1)
                          primary_t = 0;        /* turn off primary_t switch if secondary sensor is requested */
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

   if (! nfiles) {
       fprintf(stderr,"\nMust specify input station files. \n");
       exit(1);
   }
    if (! (sflag && nflag && cflag )) {
       fprintf(stderr,"\nYou must specify ship, country and cruise info.\n");
       exit(1);
   }
  
   /* set up some assumed structures */

   hdr.nprops = NPROPS;   
   hdr.prop_id = (int *) calloc((size_t) hdr.nprops, sizeof(int));
   
   hdr.prop_id[0] = (int)PR;
   hdr.prop_id[1] = (int)DE;
   hdr.prop_id[2] = (int)TE;
   hdr.prop_id[3] = (int)SA;
   hdr.prop_id[4] = (int)OX;
   hdr.prop_id[5] = (int)SI;
/*  hdr.prop_id[6] = (int)P4;  */
/* hdr.prop_id[7] = (int)N2; */
/*  hdr.prop_id[8] = (int)N3; */
   
   p = (double *) calloc ((size_t) MAXOBS, sizeof(double));
   d = (double *) calloc ((size_t) MAXOBS, sizeof(double));
   t = (double *) calloc ((size_t) MAXOBS, sizeof(double));
   s = (double *) calloc ((size_t) MAXOBS, sizeof(double));
   o = (double *) calloc ((size_t) MAXOBS, sizeof(double));
   si = (double *) calloc ((size_t) MAXOBS, sizeof(double));
/*   p4 = (double *) calloc ((size_t) MAXOBS, sizeof(double));   */
/*   n2 = (double *) calloc ((size_t) MAXOBS, sizeof(double));   */
/*   n3 = (double *) calloc ((size_t) MAXOBS, sizeof(double));   */
      
  /* loop for each file ...*/   
  
  do {
  
     strcpy(fname, dir); 
     strncat(fname,"/", 1);
     if (pflag)
        strcat(fname, prefix);
     strcat(fname, argv[curfile]); 
     strcat(fname, extent);
  
     infile = fopen(fname, "r");
     if (infile == NULL) {
        fprintf(stderr,"\nUnable to open %s for input\n\n.", fname);
	goto NEXTFILE;
     }
     fprintf(stderr,"Opened %s\n", fname);

     strcpy(fname, dir); 
     if (dir != "") 
        strncat(fname,"/", 1);
     if (pflag)
        strcat(fname, prefix);
     strcat(fname, argv[curfile]); 
     strcat(fname, ".hdr");
  
     hdrfile = fopen(fname, "r");
     if (hdrfile == NULL) {
        fprintf(stderr,"\nUnable to open header file: %s \n\n", fname);
        exit(1);
     }
     fprintf(stderr,"Opened %s\n", fname);

     read_hdr(hdrfile); 
     fclose(hdrfile);
         
          /* adjust some of the header info ...   */
     
     nobs = read_data(infile);
     ++staread;
    
     
     /* adjust some of the header info ...   */
 
     if (sscanf(argv[curfile], "%d", &hdr.station) != 1) {
        fprintf(stderr,"\nError parsing station # from this string <%s> \n", argv[curfile]);
        exit(1);
     }   

     dlat = (double) hdr.lat;
     hdr.nobs = data.nobs = nobs;
     data.nprops = hdr.nprops; 
     
     hdr.ms10 = ms10(hdr.lat, hdr.lon, &hdr.ms1);
     data.nobs = hdr.nobs;
     data.nprops = hdr.nprops;
            
          /* the data were loaded starting at the end of the arrays
             so set the pointers to the shallowest observation */  
            
      i = MAXOBS - hdr.nobs;
        
     data.observ[(int)PR] = &p[i];
     data.observ[(int)DE] = &d[i];
     data.observ[(int)TE] = &t[i];
     data.observ[(int)SA] = &s[i];
     data.observ[(int)OX] = &o[i];
     data.observ[(int)SI] = &si[i];
  
     for (j = i; j< MAXOBS; ++j) {
        d[j] = hb_depth(p[j], dlat) ;
     }

     if (hdr.nobs > 0 )  {
        write_hydro_station(outfile, &hdr, &data);
          ++staout;
    }  /* end if */
 
NEXTFILE:
      fclose(infile);       
        
   } while (curfile++ < nfiles);       

   fprintf(stderr,"\nEnd of conversion.\n");
   fprintf(stderr,"\n  %d stations read in\n", staread);
   fprintf(stderr,"\n  %d stations converted\n", staout);
   exit(0);
        
} /* end main() */ 
  
/*****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\n%s converts Seabird bottle files (merged with salt, ox, nuts) to HydroBase\n", program);
   fprintf(stderr,"\nUsage:  %s filelist -Ccruise -Sship -Nnation_code -Ooutfile -Ttemp_sensor_# [-D<dirname>] [-E<file_extent>] [-Pfilename_prefix]", program);

   fprintf(stderr,"\n\n  List of station ids MUST be first argument!");
   fprintf(stderr,"\n    -C  : specify cruise  ");
   fprintf(stderr,"\n   [-D] : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n   [-E] : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n    -N  : specify country code  ");
   fprintf(stderr,"\n   [-O] : name of outfile  ");
   fprintf(stderr,"\n    -P  : prefix for input files;  ex: -P161-7");
   fprintf(stderr,"\n    -S  : specify ship  ");
   fprintf(stderr,"\n    -T  : specify temperature sensor 1 or 2;  ex: -T1 (default is 1)");

   fprintf(stderr,"\n\n");  
   return;
}
   

/*****************************************************************************/
void read_hdr(FILE *fptr)
{
  char  line[100], *st;
  int  n, error;
  float deg, min;
  char  latfound, lonfound, depfound, datefound, hem, month[15];

  latfound = lonfound = depfound = datefound = 0;

  while (! (latfound && lonfound && depfound && datefound)   && ! feof(fptr)) {
      if ( fscanf(fptr, "%[^\n]", line) != 1) {
          fprintf(stderr, "\nError reading .hdr file \n");
          exit(1);
      }

    fgetc(fptr);   /* move past LF */
     if (! latfound) {
       if ((st = strstr(line, "Latitude")) != NULL) {
          latfound = 1;
          st += 11;           /* move past descriptor */
          if (sscanf(st, "%f %f %c", &deg, &min, &hem) != 3) {
             fprintf(stderr, "\nError parsing latitude \n");
             exit(1);
          }
          hdr.lat = deg + min/60;
         if (hem == 'S' || hem == 's') 
            hdr.lat = -hdr.lat;      
       }
    }

    if (! lonfound) {
        if ((st = strstr(line, "Longitude")) != NULL) {
          lonfound = 1;
          st += 12;           /* move past descriptor */
          if (sscanf(st, "%f %f %c", &deg, &min, &hem) != 3) {
             fprintf(stderr, "\nError parsing longitude \n");
             exit(1);
          }
          hdr.lon = deg + min/60;
         if (hem == 'W' || hem == 'w') 
            hdr.lon = -hdr.lon;      
       }
    }

    if (! depfound) {
       if ((st = strstr(line, "Depth")) != NULL) {
          depfound = 1;
          st += 6;           /* move past descriptor */
          if (sscanf(st, "%d", &hdr.pdr) != 1) {
             fprintf(stderr, "\nError parsing pdr depth \n");
             exit(1);
          }
       }
    }

    if (! datefound) {
       if ((st = strstr(line, "UpLoad Time")) != NULL) {
          datefound = 1;
          st += 14;           /* move past descriptor */
          if (sscanf(st, "%s %d %d", &month[0], &hdr.day, &hdr.year) != 3) {
             fprintf(stderr, "\nError parsing date in hdr \n");
             exit(1);
          }
         error = 0;
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
                       error = 1;
                       break;                
               }
              break;
            case 'D':
               switch (month[1]) {
                   case 'e':
                       hdr.month = 12;
                       break;
                   default:
                       error = 1;
                       break;                
               }
              break;
            case 'F':
               switch (month[1]) {
                   case 'e':
                       hdr.month = 2;
                       break;
                   default:
                       error = 1;
                       break;                
               }
               break;
           case 'J':
               switch (month[1]) {
                  case 'a':
                      hdr.month = 1;
                      break;
                  case 'u':
                      if (month[2] == 'l')
                            hdr.month = 7;
                       if (month[2] == 'n')
                            hdr.month = 6;
                      break;
                   default:
                       error = 1;
                       break;                
               }
              break;

            case 'M':
               switch (month[1]) {
                   case 'a':
                       if (month[2] == 'r')
                          hdr.month = 3;
                       
                       if (month[2] == 'y')
                          hdr.month = 5;
                      break;

                   default:
                       error = 1;
                       break; 
              }
              break;
            case 'N':
                   hdr.month = 11;
                   break;
 
            case 'O':
                   hdr.month = 10;
                   break;

            case 'S':
                   hdr.month = 9;
                   break;

            default:
                    error = 1;
         } /* end switch */


          if (error) {
             fprintf(stderr,"\nUnable to parse month: %s (%d)\n", month, hdr.month);
             exit(1);
         }
          
      }  /* end if strstr  */
   } /* end if ! datefound */

  } /* end while */

  return;

}  /* end read_hdr() */

/*****************************************************************************/
int read_data(FILE *fptr)

/* Reads bottle data and quality codes for an entire station and returns the number
    of observations which have a minimum of pr, te, sa measurements.
    It is assumed that
   the stations are in ascending order in the data file.  Within each station
   the observations are stacked from deepest to shallowest and so are loaded
   into the arrays backward starting at the end of each array. */
 
{
   int  n, i, nobs;
   char flag[20], staOK, sflag, oflag;
   char buffer[BUFSIZE];


  /* Read 4 header records */

   for (i = 0; i < 4; ++i)  {
         fscanf(fptr, "%[^\n]", buffer);
         fgetc(fptr);   /* move past LF */
   }

   i = MAXOBS;
   --i;   /* set this pointer to end of arrays*/ 
   nobs = 0;
    
  /* loop through data scans in the file */
   
  while (  fscanf(fptr, "%[^\n]", buffer) != EOF) {

         fgetc(fptr);   /* move past LF */

      if (primary_t)
          n = sscanf(buffer, "%*d %lf %lf %*lf %*lf %*lf %*lf %*lf %*lf %lf %lf %*lf %*lf %lf %*lf %s",  &p[i], &t[i], &s[i], &o[i], &si[i], &flag[0]);
     else 
          n = sscanf(buffer, "%*d %lf %*lf %lf %*lf %*lf %*lf %*lf %*lf %lf %lf %*lf %*lf %lf %*lf %s",  &p[i], &t[i], &s[i], &o[i], &si[i], &flag[0]);

  
     if (n != 6) {
       fprintf(stderr,"\n error in sscanf attempt to parse data\n%s\n", buffer);
       exit(1);
     } 
     
     /* convert ITS-90 temperature to ITS-68 */
         t[i] = 1.00024 * t[i];

   /* check flags and decrement index */
    
     sflag = flag[8];
     oflag = flag[9];     
     staOK = (p[i] > -8.) && (s[i] > -8.) && (t[i] > -8.);
     staOK = staOK && (sflag == '2');
     
     
     if (staOK) {
        if (oflag != '2')
           o[i] = -9.;
        if (o[i] > 20)
           o[i] = ox_kg2l(o[i], p[i], t[i], s[i]);
        --i;
        ++nobs;
     }

    }  /* end while */  
   
   return(nobs);
}  /* end readdata() */
