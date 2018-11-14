/*  seabird_dcc_convert.c
................................................................................
                          *******  HydroBase2 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Institution
                             May 2000
................................................................................
................................................................................
.   Reads  .dcc and .hdr files 
.   to produce HydroBase station files. 
.   USAGE: seabird_dcc_convert station_list -Pprefix -Ooutfile -Ncountry_code
           -Sship_code -Ccruise_id [-D<dir>] [-E<extent>]
................................................................................
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "hydrobase.h"

#define   DELTAP        2.0   /* pressure interval of original file */

#define   MISSING   -9.0   /* missing value flag */

#define   DIR    ""
#define   EXTENT ".dcc"


   /* prototypes for locally define functions */

void print_usage(char *);
int read_data(FILE *);
void read_hdr(FILE *);


  /* globally defined variables */
  
struct HYDRO_HDR  hdr;
struct HYDRO_DATA data;
double p[5000], d[5000], t[5000], s[5000], o[5000];
double t2[5000], s2[5000];
int primary_t;      /* flag to use primary temperature sensor */     

main (int argc, char **argv)
{
   int error, nobs, nprops;
   int outfile;
   int  i, j, n, sta, curfile = 1, nfiles = 0; 
   short sflag, cflag, nflag, oflag, pflag; 
   char *dir, *extent, fname[100], *prefix;
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
                        error = sscanf(&argv[i][2], "%d", &hdr.cruise) != 1;
                        break;
               case 'D':                   /* get input dir */
                        dir = &argv[i][2];
                        break;
               case 'E':                    /* get file extent */
                        extent = &argv[i][2];
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

               case 'P':                    /* get prefix for station_filenames */
                        pflag = 1;
                        prefix = &argv[i][2];
                        break;

               case 'S':                    /* get ship  */
                        sflag = 1;
                        strncpy(hdr.ship, &argv[i][2],2);
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

   if (! (sflag && nflag && cflag && oflag &&  nfiles )) {
       fprintf(stderr,"\nYou must specify station id(s),  output file ship, country and cruise.\n");
       exit(1);
   }



  hdr.nprops =5;
  hdr.prop_id = (int *) malloc(hdr.nprops * sizeof(int));
  hdr.prop_id[0] = (int)PR;
  hdr.prop_id[1] = (int)DE;
  hdr.prop_id[2] = (int)TE;
  hdr.prop_id[3] = (int)SA;
  hdr.prop_id[4] = (int)OX;

     /* loop for each station */
  
  do {
  
     strcpy(fname, dir); 
     if (dir != "")
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
  
     nobs = read_data(infile);
     fclose(infile);  
    
     
     /* adjust some of the header info ...   */
 
     if (sscanf(argv[curfile], "%d", &hdr.station) != 1) {
        fprintf(stderr,"\nError parsing station # from this string <%s> \n", argv[curfile]);
        exit(1);
     }   

     dlat = (double) hdr.lat;
     hdr.ms10 = ms10(hdr.lat, hdr.lon, &hdr.ms1);
                
     hdr.nobs = data.nobs = nobs;
     data.nprops = hdr.nprops; 
         
     
     for (i=0; i<nobs; ++i) {
       d[i] = hb_depth(p[i], dlat) ;
     }
     
     data.observ[(int)PR] = p;   /* set these pointers */
     data.observ[(int)DE] = d;
     data.observ[(int)TE] = t;
     data.observ[(int)SA] = s;
     if (! primary_t) {
         data.observ[(int)TE] = t2;
         data.observ[(int)SA] = s2;
     }
     data.observ[(int)OX] = o;
     if (hdr.nobs > 0 )
        write_hydro_station(outfile, &hdr, &data);
     
NEXTFILE:
     ; 
   } while (++curfile < nfiles);

   fprintf(stderr,"\nEnd of seabird_dcc_convert\n");  
   exit(0);
        
} /* end main() */ 
  
/*****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\n%s converts SeaBird .dcc (post-cal) files to HydroBase\n", program);
   fprintf(stderr,"\nUsage:  %s filelist -P<prefix> -C<cruise> -S<ship> -N<nation_code> -T<sensor_channel> -Ooutfile  [-D<dirname>] [-E<file_extent>]", program);
   fprintf(stderr,"\n\n  List of station ids MUST be first argument!");
   fprintf(stderr,"\n    -O  : name of outfile  ");
   fprintf(stderr,"\n    -P  : prefix for input files;  ex: -P161-7");
   fprintf(stderr,"\n    -C  : specify cruise  ");
   fprintf(stderr,"\n    -N  : specify country code  ");
   fprintf(stderr,"\n    -S  : specify ship  ");
    fprintf(stderr,"\n    -T  : specify temp/salt sensor (1 or 2) ex: -T2  (default is 1)  ");
  fprintf(stderr,"\n    -D  : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n    -E  : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");

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
{
   int nobs, n;
   double blah;
   char line[100];

   nobs = 0;   
   while ((n = fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf%[^\n]",  &p[nobs],  &t[nobs], &t2[nobs], &blah, &blah, &blah,  &blah, &blah,  &o[nobs], &s[nobs], &s2[nobs], &line)) != EOF) {
      if (n < 11 ){
         fprintf(stderr,"\nError reading data line line #%d.\n", nobs);
         exit(1);
      }

     /* convert ITS-90 temperature to ITS-68 */
         t[nobs] = 1.00024 * t[nobs];
         t2[nobs] = 1.00024 * t2[nobs];


      ++nobs;

      fgetc(fptr);   /* move past LF */
      
   } /* end while */

   return (nobs);
}  /* end readdata() */
