/*  igor_ctdconvert.c
................................................................................
.   Reads Igor Yashayaev's data files
.   extracts :    header info
.                 p,theta,s,ox 
.   and outputs p,d,t,s and ox at to the output
.   file.
.
.   USAGE: igor_ctdconvert infile_list -Ooutfile  
          [-D<dir>] [-E<extent>] -N<nation_code> -S<ship_code> -C<cruise_no>
................................................................................
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "hydrobase.h"
#include "hb_memory.h"


#define   DIR    ""
#define   EXTENT ""
#define   NPOUT 4


  /* globally defined variables */
  
struct HYDRO_HDR  hdr;
struct HYDRO_DATA data;
double *p, *d, *t,  *s;
char buffer[300];  /* store next line in file */

  /* prototypes for locally defined functions */
  
void print_usage(char *);
FILE *openfile(char *, char *, char *);
int readdata(FILE *, int *);


int main ( int argc, char **argv)
{
   int error, nobs, nprops, npout, more;
   int outfile;
   int  i, j, curfile, nfiles, nsta;
   int sflag, oflag, nflag, cflag; 
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
    cflag = nflag = sflag = oflag  = 0;
    npout = NPOUT;
    curfile = 1;
    nfiles = 0;
    buffer[0] = '\0';
 
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
   
     
     /*  allocate space */
     
     nobs = 6000;       
     p = (double *) get_memory((void *)NULL, (size_t)nobs, sizeof(double));
     d = (double *) get_memory((void *)NULL, (size_t)nobs, sizeof(double));
     t = (double *) get_memory((void *)NULL, (size_t)nobs, sizeof(double));
     s = (double *) get_memory((void *)NULL, (size_t)nobs, sizeof(double));

     hdr.prop_id = (int *) calloc((size_t)npout, sizeof(int));
     hdr.prop_id[0] = (int)PR;
     hdr.prop_id[1] = (int)DE;
     hdr.prop_id[2] = (int)TE;
     hdr.prop_id[3] = (int)SA;
     hdr.nprops = npout;
     nsta = 0;
    
   
    
/* loop for each input file */

   do {
   
     infile = openfile(dir, argv[curfile], extent);
     if (infile == NULL) goto NEXTFILE;
      
     do {
    
    /* read each station */    
        
       more =  readdata(infile, &nobs);
       if (nobs > 0) {       
         dlat = (double) hdr.lat;
         for (i = 0; i < nobs; ++i) {
          d[i] = hb_depth(p[i], dlat);
          t[i] = hb_theta( s[i], t[i], 0.0, p[i]);
         }
     
     
         /* adjust some of the header info ...   */
     
         hdr.ms10 = ms10(hdr.lat, hdr.lon, &hdr.ms1);
     
         hdr.nobs = data.nobs = nobs;
         data.nprops = hdr.nprops; 
     
         data.observ[(int)PR] = p;   /* set these pointers */
         data.observ[(int)DE] = d;
         data.observ[(int)TE] = t;
         data.observ[(int)SA] = s;
        
         write_hydro_station(outfile, &hdr, &data);
	 ++nsta;
	 
      } /* end if nobs > 0 */
    } while (more);

NEXTFILE:
         fclose(infile);

 
   } while (curfile++ < nfiles );

   fprintf(stderr,"%d stations written out.\n", nsta);
   fprintf(stderr,"End of %s\n", argv[0]);
   exit(0);
        
} /* end main() */ 
  
/*****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nUsage:  %s filelist -C<cruise> -N<nation> -O<outfile>  -P<prs_int> -S<shipcode> [-D<dirname>] [-E<file_extent>]", program);
   fprintf(stderr,"\n\n  List of filenames must be first argument.");
   fprintf(stderr,"\n    -C  : specifies cruise number ");
   fprintf(stderr,"\n    -N  : 2-char country code ");
   fprintf(stderr,"\n    -O  : specifies output file ");
   fprintf(stderr,"\n    -P  : specifies pressure interval of output series ");
   fprintf(stderr,"\n    -S  : 2-char ship code ");
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
/*****************************************************************************/
int readdata(FILE *fptr, int *nobs_addr)

/* Reads ctd p,t,s for an entire station. Returns 1 if there is more data, 
   or 0 for EOF. */
 
{
   int i, n, error, samesta;
   
   n = 0;

   if (*buffer == '\0') {
     if ((i = fscanf(fptr, "%[^\n]", buffer)) == EOF) {
         *nobs_addr = 0;
         return (0);
     }
     if (i != 1) {
         fprintf(stderr,"\nError reading line into buffer:\n%s\n", buffer);
         exit(1);
     }
     error = getc(fptr);
   }

  /* parse out the header info */
  
   error = sscanf(buffer, "%f,%f,%*f,%*f,%*f,%*f,%*f,%d,%d,%d", &hdr.lat, &hdr.lon, &hdr.day, &hdr.month, &hdr.year) != 5;
   if (error) {
      fprintf(stderr,"\nError parsing header info:\n%s\n", buffer);
      exit(1);
   }
   hdr.lon = -hdr.lon;
   samesta = 1;
      
   while (samesta) {  
     
      error = sscanf(buffer, "%*f,%*f,%*f,%lf,%lf,%lf", &p[n],&s[n],&t[n]) != 3;
      if (error) {
         fprintf(stderr,"\nError parsing data:\n%s\n", buffer);
         exit(1);
      }
      p[n] = -p[n];
      ++n;
       
      /* read next line into buffer */
      
      if ((i = fscanf(fptr, "%[^\n]", buffer)) == EOF) {
         *buffer = '\0';
         *nobs_addr = n;
         return (0);
      }
      
     if (i != 1) {
         fprintf(stderr,"\nError reading line into buffer:\n%s\n", buffer);
         exit(1);
     }
     error = getc(fptr);
     
     error = sscanf(buffer, "%*f,%*f,%*f,%lf", &p[n]) != 3;
     p[n] = -p[n];
     samesta = p[n] > p[n-1];
     
        
  } /* end while */
   
   *nobs_addr = n;
   return (1);
   
}  /* end readdata() */
