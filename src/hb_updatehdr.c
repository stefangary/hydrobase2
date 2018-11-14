/* hb_updatehdr.c
................................................................................
                         *******  HydroBase2 *******
................................................................................
                          
                    author:  Ruth  Curry
                             Woods Hole Oceanographic Inst
                             July 2000 
................................................................................
................................................................................
.  reads stations and applies changes to header 
.
................................................................................
................................................................................

*/
#include <stdio.h>
#include <stdlib.h>   /*sfg added to deal with exit and malloc warnings*/
#include <ctype.h>
#include <string.h>
#include "hydrobase.h"

#define    DIR     ""
#define    EXTENT   ""
#define    PRINT_MSG 1

void print_usage(char *);

main (int argc, char **argv)
{
   int     curfile = 1, nfiles = 0; 
   int     i, infile, m;
   int     m_flag, y_flag, msq_flag;
   int     n_flag, s_flag, c_flag;
   int     o_flag, q_flag, i_flag;
   char    *ship, *country;
   char    orig, instr;
   char    *qual_mask;
   int     status;
   int     error = 0;
   int     cru, year, month;
   struct HYDRO_HDR h;
   struct HYDRO_DATA data;
   char   *dir, *extent;

/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
   
/* allocate char arrays */

       ship = (char *) calloc(3, (size_t) sizeof(char));
       country = (char *) calloc(3, (size_t) sizeof(char));

/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    cru = 0;
    m_flag = y_flag = n_flag = s_flag = msq_flag = 0;
    c_flag = o_flag = q_flag = i_flag = 0;

/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'C' :
                        c_flag = 1;
                        error = sscanf(&argv[i][2],"%d", &cru) != 1;
                        break;
               case 'D':
                        dir = &argv[i][2];
                        break;
               case 'E':
                        extent = &argv[i][2];
                        break;
               case 'I' :
                        i_flag = 1;
                        instr = argv[i][2];
                        break;
               case 'M' :
                        m_flag = 1;
                        error = sscanf(&argv[i][2],"%d", &month) != 1;
                        break;
               case 'N' :
                        n_flag = 1;
                        error = (sscanf(&argv[i][2],"%2s", country) != 1);
                        break;
               case 'O' :
                        o_flag = 1;
                        orig = argv[i][2];
                        break;
               case 'Q' :
                        q_flag = 1;
                        qual_mask = &argv[i][2];
                        if (strlen(qual_mask) != 4) {
                          fprintf(stderr,"Quality mask must have 4 characters.");
                          fprintf(stderr,"Ex: -Q0011");
                          fprintf(stderr," a '0' means no quality control");
                          fprintf(stderr," a '1' means the property has been checked");
                          fprintf(stderr,"the first byte or character (from left) is for cfs's");
                          fprintf(stderr,"the second is for nutients");
                          fprintf(stderr,"the third is for oxygen");
                          fprintf(stderr,"the fourth is for temperature and salinity");
                          fprintf(stderr,"In the above example, the mask will set the quality bytes for T,S and Ox to '1' for all stations.  The nutrient and cfc bytes will remain whatever they are already.");
                          exit(1);
                        }
                        break;
               case 'S' :
                        s_flag = 1;
                        error = (sscanf(&argv[i][2],"%2s", ship) != 1);
                        break;
               case 'W' :
                        msq_flag = 1;
                        break;
               case 'Y' :
                        y_flag = 1;
                        error = sscanf(&argv[i][2],"%d", &year) != 1;
                        break;
               case 'h':
                        print_usage(argv[0]);
                        exit(0);
               default  :
                        fprintf(stderr,"\nError parsing command line");
                        fprintf(stderr,"\n in particular: %s\n", argv[i]);
                        exit(1);
            }  /* end switch */
       }  /* end if */
       else  {
           ++nfiles;
       }
   }  /* end for */

   if (!nfiles) {
       print_usage(argv[0]);
       fprintf(stderr,"\n\nYou must specify input file(s)!\n");
       exit(1);
   }

/* initialize these... */

   for (i = 0; i < MAXPROP; ++i)
      data.observ[i] = (double *) NULL;
   h.prop_id = (int *) NULL;

 /* loop for each input file */

   do {

     infile = open_hydro_file(dir, argv[curfile], extent, PRINT_MSG);
     if (infile  < 0) 
       goto NEXTFILE;
     

     /* loop for each station */

     while ((status = get_station(infile, &h, &data)) == 0) { 

          
       if (msq_flag) {
          h.ms10 = ms10(h.lat, h.lon, &h.ms1);
       }
       if (c_flag) {
          h.cruise = cru;
        }
       if (s_flag ) {
          strncpy(h.ship, ship, 3);
        }
       if (n_flag ) {
           strncpy(h.country, country, 3);
        }

       if (y_flag) {
          h.year = year;
        }
        
       if (m_flag) {
          h.month = month;
        }
	
       if (i_flag) {
          h.instrument = instr;
       }
       
       if (o_flag) {
          h.origin = orig;
       }
       
       if (q_flag) {
          for (i = 0; i < NQUAL; ++i) {
            h.qual[i] = (h.qual[i] - '0') | (qual_mask[i] - '0');  /*bitwise or */
            h.qual[i] += '0';
          }
       }

        error = write_hydro_station(STDOUT, &h, &data);
        if (error) {
          fprintf(stderr,"\nError in write_hydro_station(): error code = %d.\n", error);
          exit(1);
        }

     }  /* end while */

     report_status(status, stderr);
     close(infile);

NEXTFILE:
     ;
   } while (curfile++ < nfiles);

   fprintf(stderr,"\n\nEnd of %s.\n\n", argv[0]);
   exit(0);

} /* end main() */




/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nUsage:  %s filename_root(s)", program);

   fprintf(stderr," [-Ddirname] [-Eextent] [-Mmonth] [-Ccruise] [-Sship] [-Nnation] [-Yyear] [-W][-Iinstrument_code] [-Oorigin_code] [-Qquality_mask]  ");
   fprintf(stderr,"\n");
   fprintf(stderr,"\n    -D  : specifies dirname (default is ./) ");
   fprintf(stderr,"\n        ex: -D../data/ ");
   fprintf(stderr,"\n    -E  : specifies file extent (default is no extent)");  
   fprintf(stderr,"\n        ex: -E.dat ");
   fprintf(stderr,"\n    -W  : update WMO square designation");
   fprintf(stderr,"\n    -M  : update month");
   fprintf(stderr,"\n    -C  : update cruise #");
   fprintf(stderr,"\n    -S  : update ship code");
   fprintf(stderr,"\n    -N  : update country code");
   fprintf(stderr,"\n    -Y  : update year");
   fprintf(stderr,"\n    -I  : update instrument type: (b,c,f,s,m, or u");
   fprintf(stderr,"\n    -O  : update origination code: ('0' .. '9')");
   fprintf(stderr,"\n    -Q  : update quality code with a 4-byte mask");
   fprintf(stderr,"\n          ex: -Q1000  :  a bitwise OR of this mask is applied to each station's quality bytes.");
   fprintf(stderr,"\n    -h  : help -- prints this message");
  fprintf(stderr,"\n\n");  
   return;
}

/****************************************************************************/
