/*  ices_convert.c
................................................................................
                          *******  HydroBase *******
               ...................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
                             1993
................................................................................
................................................................................
.   Reads ICES "Punch Card" format file(s)
.  
.                  
.
.   USAGE: ices_convert infile_list -Ooutfile [-D<dir>] [-E<extent>]
................................................................................
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "hydrobase.h"
#include "hb_filters.h"


#define   LF        0x0a   /* ascii code for linefeed char */
#define   MISSING   -9.0   /* missing value flag */
#define   MAXOBS    6000

#define   DIR    ""
#define   EXTENT ""

  /* define possible modes for output files ... */
#define   OVERWRITE 0
#define   NO_CLOB   1
#define   APPEND    2


  /* globally defined variables */
  
struct HYDRO_HDR  hdr;
struct HYDRO_DATA data;
double *p, *d, *t, *s, *o, *ph, *ni, *si, *na;
int n_p, n_o, n_ph, n_ni, n_si, n_na;
char *buffer;
int suppress;


   /* prototypes for locally defined functions */

void print_usage(char *);
FILE *openfile(char *, char *, char *);
void parse_header(char *);
int parse_hydro_rec(char *, int);
int parse_chem_rec(char *, char *,int);
double convert_kg2l(double, double, double, double);
double convert_l2kg(double, double, double, double);
int pressure_sort(int);

int main (int argc, char **argv)
{
   int error, maxobs, nsta, nout, nerr;
   int eof, nextsta, nfilesout;
   int outfile, errfile, merror, pr_int;
   int  i, j, curfile = 1, nfiles = 0; 
   short oflag; 
   char *dir, *extent;
   char *hline, *cline;
   char *rootout, fname[100];
   FILE *infile;
   double dlast;
   
/* check for command line arguments */

   if (argc < 2) {
     print_usage(argv[0]);
     exit(1);
   }
   
 
/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    error = 0;
    oflag  = 0;
    cline = hline = buffer = NULL;
    maxobs = 60000;
    pr_int = 10;  /* default pressure series interval for output */
    nfilesout;
    suppress = 0;   /* write out WARNING only once */
 
 
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
			rootout = &argv[i][2];
                        if (outfile < 1) {
                           fprintf(stderr,"\nError opening output file: %s\n", 
                                  &argv[i][2]);
                           exit(1);
                        }
                        strcpy(fname, rootout);
                        strncat(fname, ".err",  4);
                        errfile = create_hydro_file(fname, OVERWRITE);
                        break;
               case 'P':                    /* get pressure interval for output */
                        error = sscanf(&argv[i][2],"%d", &pr_int) != 1;
                        break;

               default:
                        error = 1;

          }    /* end switch */

          if (error) {
             fprintf(stderr,"\nError parsing command line args.\n");
             fprintf(stderr,"     in particular: '%s'\n", argv[i]);
             exit(1);
          }

       }  /* end if */

       else  {
           ++nfiles;
       }
   }  /* end for */

   if (! (oflag && nfiles )) {
       fprintf(stderr,"\nYou must specify input file(s)and an output file.\n");
       exit(1);
   }
   
   nsta = 0;
   nout = 0;
   nerr = 0;

/* loop for each input file */

   do {
     infile = openfile(dir, argv[curfile], extent);
     if (infile == NULL) goto NEXTFILE;
     
      
     /*  loop to read all stations in file */
     
     suppress = 0;
     eof = 0;
     buffer = (char *) calloc(100, sizeof(char));
     error = fscanf(infile, "%[^\n]", buffer);
       
     if (error == EOF)
          eof = 1;
     else {
        if (error != 1) {
	     fprintf(stderr,"\nERROR reading input file at station #%d\n", nsta+1);
	     exit(1);
	}
	
	error = getc(infile);  /* move past LF */
	
        if (!(buffer[79] == 'J') ) {
	     fprintf(stderr,"\nERROR reading header record for station #%d\n%s\n", nsta+1, buffer);
	     exit(1);
	}
     }
       
     do {
     
        /*  allocate space and initialize bottle props*/
     
        hdr.instrument = 'c';
     
        p = (double *) calloc(maxobs, sizeof(double));
        d = (double *) calloc(maxobs, sizeof(double));
        t = (double *) calloc(maxobs, sizeof(double));
        s = (double *) calloc(maxobs, sizeof(double));
        o = (double *) calloc(maxobs, sizeof(double));
        ph = (double *) calloc(maxobs, sizeof(double));
        si = (double *) calloc(maxobs, sizeof(double));
        ni = (double *) calloc(maxobs, sizeof(double));
        na = (double *) calloc(maxobs, sizeof(double));

        for (i = 0; i < maxobs; ++i) {
	   o[i] = MISSING;
	   ph[i] = MISSING;
	   si[i] = MISSING;
	   na[i] = MISSING;
	   ni[i] = MISSING;
	}
	
        n_p = n_o = n_ph = n_si = n_ni = n_na = 0;
	        
	  
        while (buffer[79] != 'J') {
	
	       fprintf(stderr,"Expecting Master Record but code does not match (skip for now):\n%s\n", buffer);
               error = fscanf(infile, "%[^\n]", buffer);
	       error = getc(infile);  /* move past LF */
	       if (error == EOF)
	           eof = 1;
	}  /* end while */
        
	parse_header(buffer);
	free(buffer);
	buffer = NULL;
	hdr.nobs = 0;
        nextsta = 0;
        ++nsta;
	    
       /* loop to read all records for a single station */
	  
     
       while (! (eof || nextsta)) {
       
         buffer = (char *) calloc(100, sizeof(char));
         error = fscanf(infile, "%[^\n]", buffer);
         if (error != 1) {
	     if (error == EOF)
	        eof = 1;
	     else {
	       fprintf(stderr,"\nERROR reading input file at station #%d\nContents of buffer:\n%s\n", nsta, buffer);
	       exit(1);
	     }
	 }
	 else {
            error = getc(infile);   /* move past LF char*/
	 
	    switch (buffer[79]) {
	    
	       case 0:
	         eof = 1;
		 break;
	 
	       case '3': 
	          if (hline != NULL)
		     free(hline);
	          hline = buffer;
	          hdr.nobs = parse_hydro_rec(hline, hdr.nobs);
	          break;
	    
	       case '6':
	          cline = buffer;
	          hdr.nobs = parse_chem_rec(cline, hline, hdr.nobs);
                  hdr.instrument = 'b';
	          free((char *)cline);
	          cline = (char *)NULL;
	          if (hline != NULL) {
	            free((char *)hline);
		    hline = (char *)NULL;
	          }
	          break;
	    
	       case 'J':
	       case '2':
	          nextsta = 1;
	          break;
	    
	       case 'Z':
	          fprintf(stderr,"\nSkipping:  %s",buffer);
	          free((char *)buffer);
	          buffer = (char *)NULL;	    
	          break;
	       
	       
		  
	       default:

                  	    
	           fprintf(stderr,"\nERROR parsing line from input file:\n%s\n", buffer);
	           fprintf(stderr,"\nUnknown record type\n");
	           exit(1);
	     
	 
	    } /* end switch */
	 } /* end else */
    
       } /* end while */
       
       
       /* output the station */
       
       if (hdr.nobs > 0 ) {
          for (i = 0; i < NQUAL; ++i)
            hdr.qual[i] = '0';
          hdr.origin = '2';
       
          hdr.ms10 = ms10(hdr.lat, hdr.lon, &hdr.ms1);

        
          hdr.nprops = 4;
          if (n_o > 0) 
             hdr.nprops++;
          if (n_ph > 0) 
             hdr.nprops++;
          if (n_si > 0) 
             hdr.nprops++;
          if (n_na > 0) 
             hdr.nprops++;
          if (n_ni > 0) 
             hdr.nprops++;
	  
          hdr.prop_id = (int *) malloc(hdr.nprops * sizeof(int));
          hdr.prop_id[0] = (int)PR;
          hdr.prop_id[1] = (int)DE;
          hdr.prop_id[2] = (int)TE;
          hdr.prop_id[3] = (int)SA;
	  
          i = 4;
          if (n_o > 0) {
             hdr.prop_id[i++] = (int)OX;
	  }
          if (n_ph > 0) { 
             hdr.prop_id[i++] = (int)P4;
	  }
          if (n_si > 0) { 
             hdr.prop_id[i++] = (int)SI;
	  }
          if (n_na > 0) { 
             hdr.prop_id[i++] = (int)N2;
	  }
          if (n_ni > 0) { 
             hdr.prop_id[i++] = (int)N3;
	  }
	  
/* interpolate CTD data onto 10-db pressure series */

          if (hdr.instrument == 'c' ) {
	    hdr.nobs = pressure_sort(pr_int );
	  }
	  
	  data.observ[(int)PR] = p;
	  data.observ[(int)DE] = d;
	  data.observ[(int)TE] = t;
	  data.observ[(int)SA] = s;
	  
          if (n_o > 0) {
	     data.observ[(int)OX] = o;
	  }
          if (n_ph > 0) { 
	     data.observ[(int)P4] = ph;
	  }
          if (n_si > 0) { 
	     data.observ[(int)SI] = si;
	  }
          if (n_na > 0) { 
	     data.observ[(int)N2] = na;
	  }
          if (n_ni > 0) { 
	     data.observ[(int)N3] = ni;
	  }
          data.nprops = hdr.nprops;
          data.nobs = hdr.nobs;	
          
          
          /* check for monotonically increasing depth as indicator of decimal place problem */
          
         dlast = 0.0;      
          for (i = 0; i < hdr.nobs; ++i) {
             if (data.observ[(int)DE][i] < (dlast -2) && data.observ[(int)DE][i] < 500 ) {
                data.observ[(int)DE][i]   *= 10.0;
                data.observ[(int)PR][i]   *= 10.0;
             }
             dlast = data.observ[(int)DE][i];   
          }
        
        /* now check again */
        
           dlast = 0.0;
           merror = 0;
           for (i = 0; i < hdr.nobs; ++i) {
             if (data.observ[(int)DE][i] < (dlast-2)) 
                  ++merror;
           }
         
         
         if ( ! merror) {
             write_hydro_station(outfile, &hdr, &data);
	     ++nout;
          }  
          else {
               write_hydro_station(errfile, &hdr, &data);
               ++nerr;
          }
	  
              free(hdr.prop_id); 
         if ((nout % 1000) == 0) {
             close(outfile);
	     ++nfilesout;
	  
	     strcpy(fname, rootout);
	     sprintf(fname, "%s.%04d", rootout, nfilesout);
             outfile = create_hydro_file(fname, OVERWRITE);
             if (outfile < 1) {
                fprintf(stderr,"\nError opening output file: %s\n", fname);
                exit(1);
             }
                 fprintf(stderr,"\nOpening output file: %s\n", fname);
         }
       }

       /* clean up space... */ 
           
       free(p);
       free(d);
       free(t);
       free(s);
       free(o);
       free(ph);
       free(si);
       free(na);
       free(ni);
       
      
     } while (!eof);
     

NEXTFILE:
         fclose(infile);

 
   } while (curfile++ < nfiles );
   
   fprintf(stderr,"\n%d stations read in", nsta);
   fprintf(stderr,"\n%d stations written out in %d output files", nout, nfilesout);
   fprintf(stderr,"\nREMINDER:  ICES format does not include cruise id.");
   fprintf(stderr,"\nNumber of non-monotonic profiles: %d", nerr);
   if (nerr > 0)
      fprintf(stderr,"\nCHECK ERROR FILE");
   
   fprintf(stderr,"\nEnd of %s\n", argv[0]);
   
   exit(0);
        
} /* end main() */ 
  
/*****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nUsage:  %s filelist -Ooutfile [-D<dirname>] [-E<file_extent>]", program);
   fprintf(stderr,"\n\n  List of filenames MUST be first argument!");
   fprintf(stderr,"\n    -O  : specifies output filename "); 
   fprintf(stderr,"\n OPTIONS:");   
   fprintf(stderr,"\n   [-D]  : specifies directory for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n   [-E]  : specifies input_file extent (default is no extent)");  
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
void parse_header(char *line)

   /* reads the header info containe in buffered line and fills in global 
      variables.  An error causes an exit. */
{
   int deci, min, deg, error;
   char xx;
   
   hdr.country[0] = line[0];
   hdr.country[1] = line[1];
   hdr.country[2] = '\0';
   
   hdr.ship[0] = line[2];
   hdr.ship[1] = line[3];
   hdr.ship[2] = '\0';

   xx = line[8];
   line[8] = '\0';  
   error = sscanf(&line[4],"%d", &hdr.station) != 1;
   if (error) {
       hdr.station = 9999;
   }
   if (hdr.station < 0)
      hdr.station = 9999;
   while (hdr.station > 9999)
      hdr.station -= 10000;
      
   line[8] = xx;
   
   error = sscanf(&line[8],"%2d%2d",  &deg, &min) != 2;
   if (error) {
      fprintf(stderr,"\nError in parse_header() [latdeg, min]\n");
      exit(1);
   }
      
   xx = line[66];
   line[66] = '\0';  
   error = sscanf(&line[64],"%2d", &deci) != 1;
   if (error) {
      deci = 0;
   }   
   line[66] = xx;  
   
   hdr.lat = (float)( deg + ((min + (float)deci/100.) / 60.0));
   if (line[17] == '2' || line[17] == '3')
      hdr.lat = -hdr.lat;
      
   error = sscanf(&line[12],"%3d%2d",&deg, &min) != 2;
   if (error) {
      fprintf(stderr,"\nError in parse_header()[londeg,min]\n");
      exit(1);
   }   
   xx = line[68];
   line[68] = '\0';  
   error = sscanf(&line[66],"%2d", &deci) != 1;
   if (error) {
      deci = 0;
   }   
   line[68] = xx;  
   hdr.lon = (float)( deg + ((min + (float)deci/100.) / 60.0));
   if (line[17] == '1' || line[17] == '3')
      hdr.lon = -hdr.lon;
      
   
   
   error = sscanf(&line[18],"%3d%2d%2d", &hdr.year, &hdr.month, &hdr.day) != 3;
   if (error) {
      fprintf(stderr,"\nError in parse_header() [year,month,day]\n");
      exit(1);
   }   
   hdr.year += 1000;
   if (hdr.year < 1800)   /* for year 2000 and beyond */
        hdr.year += 1000;

   error = sscanf(&line[27],"%4d", &hdr.pdr) != 1;
   if (error && !suppress) {
      fprintf(stderr,"\nWARNING: no PDR depth\n");
      hdr.pdr = 0.0;
      suppress = 1;
   }   
   
   /*  ICES does not store cruise id; use year/month as default.  Can 
        use hb_updatehdr to replace with appropriate cruise #  */
   
   hdr.cruise = (hdr.year % 100) * 100 + hdr.month;
   if (hdr.lat > 90. || hdr.lat < -90)
      hdr.lat = 90.0;
   return;
   
}  /* end parse_header() */
/*****************************************************************************/

/*****************************************************************************/
int parse_hydro_rec(char *line, int nobs)
{
   int error, i, deci;
   char xx;
   double mult;



   /* get pressure/depth */
      
   xx = line[31];
   line[31] = '\0';
        
   switch (line[40]) {
      case 'p':
      case '+':
         error = (sscanf(&line[27],"%lf", &p[nobs]) != 1);
         line[31] = xx; 
         xx = line[43];
         line[43] = '\0';
         error = (sscanf(&line[41],"%d", &deci) != 1);
	 if (error)
	     deci = 0;
         line[43] = xx;
	 
	 mult = .01;
	   
 	 p[nobs] += (double) deci *mult;
	 d[nobs] = hb_depth(p[nobs], (double) hdr.lat);
         break;
	 
      default:
      
         error = (sscanf(&line[27],"%lf", &d[nobs]) != 1);
         line[31] = xx; 
         xx = line[43];
         line[43] = '\0';
         error = (sscanf(&line[41],"%d", &deci) != 1);
	 if (error)
	     deci = 0;
	 mult = .01;
         line[43] = xx;
	 d[nobs] += (double) deci * mult;
	 p[nobs] = hb_p80(d[nobs], (double) hdr.lat);
         break;
      
          
   }  /* end switch */
   
   /* get temperature */
   
   xx = line[35];
   line[35] = '\0';
        
   if (line[31] == '}')
        error = (sscanf(&line[32],"%lf", &t[nobs]) != 1);
   else 
        error = (sscanf(&line[31],"%lf", &t[nobs]) != 1);
   
   line[35] = xx;
   if (error) 
       t[nobs] = MISSING;
   else { 
      t[nobs] /= 100.0;
      
      xx = line[46];
      line[46] = '\0';
      if  (sscanf(&line[44],"%d", &deci) == 1) {
          mult = 0.0001;
          t[nobs] += (double) deci * mult;
      }
      line[46] = xx;
      if (line[31] == '}')
          t[nobs] = -t[nobs];
   }
   
   if(t[nobs] < -2.0)
       t[nobs] = MISSING;
   
   /* get salinity */
   
   xx = line[40];
   line[40] = '\0';
   error = (sscanf(&line[35],"%lf", &s[nobs]) != 1);
   line[40] = xx;
   if (error) 
       s[nobs] = MISSING;
   
   else { 
      s[nobs] /= 1000.0;
	    
      xx = line[49];
      line[49] = '\0';
      if ((sscanf(&line[47],"%1d", &deci) == 1)) {
          mult = 0.0001;
	  if (deci >= 10)
	    mult = 0.00001;
          s[nobs] += (double) deci * mult;
      }
      line[49] = xx;
   }
   
   if (s[nobs] < 1.0)
       s[nobs] = MISSING;
     

   /* if pr, te, or sa is missing, skip this record */
      
   if ((s[nobs] < -3.) || (t[nobs] < -3.) || (p[nobs] < -3))
       return (nobs);
       
   ++n_p;    /* increment counter */
   
   /* check for oxygen */   
   
   if (isdigit(line[57])) {
      xx = line[60];
      line[60] = '\0';
      error = (sscanf(&line[57],"%lf", &o[nobs]) != 1);
      line[60] = xx;
      if (error)
         o[nobs] = MISSING;
      else{
         o[nobs] /= 100.0;
         ++n_o;
      }
	 
      if (line[77] == 'K') 
         o[nobs] = convert_kg2l(o[nobs], p[nobs], t[nobs], s[nobs]);
      
   } /* end if isdigit()*/
   
   return(++nobs);
       
}  /* end parse_hydro_rec() */

/*****************************************************************************/
int parse_chem_rec(char *line, char *ctd_data, int nobs)
{
   int error, got_pts, deci;
   char xx;
   double mult;
   
   got_pts = 0;

   if (ctd_data) {
      got_pts = ((p[nobs] > -3.) && (t[nobs] > -3.) && (s[nobs] > -3.));
   }
   
   if (!got_pts) {
   

   /* get pressure/depth */
      
      xx = line[31];
      line[31] = '\0';
      switch (line[40]) {
         case 'p':
         case '+':
            error = (sscanf(&line[27],"%lf", &p[nobs]) != 1);
            line[31] = xx;
	    d[nobs] = hb_depth(p[nobs], (double) hdr.lat);
            break;
	 
         case 'd':  /* fall through */
         case ' ':      
            error = (sscanf(&line[27],"%lf", &d[nobs]) != 1);
            line[31] = xx;
	    p[nobs] = hb_p80(d[nobs], (double) hdr.lat);
            break;
      
         default:
          fprintf(stderr,"Unknown code at offset 40 in data record (expecting 'p','d', or ' ')\n%s\n", line);
          
      }  /* end switch */
   
      /* get temperature */
   
      xx = line[35];
      line[35] = '\0';
      
      if (line[31] == '}')
        error = (sscanf(&line[32],"%lf", &t[nobs]) != 1);
      else 
        error = (sscanf(&line[31],"%lf", &t[nobs]) != 1);
   
      line[35] = xx;
      if (error ) 
          t[nobs] = MISSING;
      else {  
          t[nobs] /= 100.0;
         if (line[31] == '}')
             t[nobs] = -t[nobs];
      }
      
      if (t[nobs] < -2.0)
         t[nobs] = MISSING;
   
   /* get salinity */
   
      xx = line[39];
      line[39] = '\0';
      error = (sscanf(&line[35],"%lf", &s[nobs]) != 1);
      line[39] = xx;
      if (error) 
          s[nobs] = MISSING;
      else  
         s[nobs] /=  100.0; 
	      
      if (s[nobs] < 1.0)
         s[nobs] = MISSING;

      /* if pr, te, or sa is missing, skip this record */
      
      if ((s[nobs] < -3.) || (t[nobs] < -3.) || (p[nobs] < -3))
          return (nobs);
       
      ++n_p;    /* increment counter */
   
   } /* end if !got_pts */
   
   if (o[nobs] < -3.) {
   
      if (isdigit(line[57])) {
         xx = line[60];
         line[60] = '\0';
         error = (sscanf(&line[57],"%lf", &o[nobs]) != 1);
         line[60] = xx;
         if (error)
            o[nobs] = MISSING;
         else {
            ++n_o;
	    o[nobs] /= 100.0;
	 }
	 
         if (line[77] == 'K') 
            o[nobs] = convert_kg2l(o[nobs], p[nobs], t[nobs], s[nobs]);
      }
   }

    /* phosphate */
       
   if (isdigit(line[42])) {
         xx = line[45];
         line[45] = '\0';
         error = (sscanf(&line[42],"%lf", &ph[nobs]) != 1);
         line[45] = xx;
         if (error)
            ph[nobs] = MISSING;
         else {
           ph[nobs] /= 100.0;;
           ++n_ph;
	 
            if (line[77] == 'K') 
               ph[nobs] = convert_kg2l(ph[nobs], p[nobs], t[nobs], s[nobs]);
	    if (line[78] == 'P')
	       ph[nobs] *=  10.0; 
	  } 
   }
    
    /* silicate */
       
   if (isdigit(line[48])) {
         xx = line[51];
         line[51] = '\0';
         error = (sscanf(&line[48],"%lf", &si[nobs]) != 1);
         line[51] = xx;
         if (error)
            si[nobs] = MISSING;
         else {
            si[nobs] /= 10.0; 
	    ++n_si;
	 
            if (line[77] == 'K') 
               si[nobs] = convert_kg2l(si[nobs], p[nobs], t[nobs], s[nobs]);        
	    if (line[78] == 'P')
	       si[nobs] *= 10.0;
	 }
   }
   
    /* nitrate */
       
   if (isdigit(line[51])) {
         xx = line[54];
         line[54] = '\0';
         error = (sscanf(&line[51],"%lf", &na[nobs]) != 1);
         line[54] = xx;
         if (error)
            na[nobs] = MISSING;
         else {
            na[nobs] /= 10.0; 
            ++n_na;
	 
            if (line[77] == 'K') 
               na[nobs] = convert_kg2l(na[nobs], p[nobs], t[nobs], s[nobs]);        
	    if (line[78] == 'P')
	       na[nobs] *= 10.0;
	 }
   }
   
    /* nitrite */
       
   if (isdigit(line[54])) {
         xx = line[57];
         line[57] = '\0';
         error = (sscanf(&line[54],"%lf", &ni[nobs]) != 1);
         line[57] = xx;
         if (error)
            ni[nobs] = MISSING;
         else {
            ni[nobs] /= 10.0; 
            ++n_ni;
	 
            if (line[77] == 'K') 
               ni[nobs] = convert_kg2l(ni[nobs], p[nobs], t[nobs], s[nobs]);        
	    if (line[78] == 'P')
	       ni[nobs] *= 10.0;
	 }
   }
   
   
   
   return (nobs);
}  /* end parse_chem_rec() */
/*****************************************************************************/
double convert_kg2l(double x, double pin, double tin, double sin)
 /* convert units/kg to units/liter by dividing by density */
 
{

   double pt, sig, sva;
   
   if (x < 0 )
      return (MISSING);
      
    pt = hb_theta(sin, tin, pin, 0.0);
    sva = hb_svan(sin, pt, 0.0, &sig);
    return (x * (1 + sig/1000.));
    

} /* end convert_kg2l()*/ 
/*****************************************************************************/

double convert_l2kg(double x, double pin, double tin, double sin)
 /* convert units/kg to units/liter by dividing by density */
 
{

   double pt, sig, sva;
   
   if (x < 0 )
      return (MISSING);
      
    pt = hb_theta(sin, tin, pin, 0.0);
    sva = hb_svan(sin, pt, 0.0, &sig);
    return (x / (1 + sig/1000.));
    

} /* end convert_kg2l()*/ 

/*****************************************************************************/
int pressure_sort(int prsint)
 /* interpolates profiles onto pressure series */
 
{
   int nout, i, j, filter_width, filter_it;
   double *pout, *yin, *yout;
   double delta_p, maxp, lastp;
   
   
   if (hdr.nobs <= 1)
      return (hdr.nobs);
 
 /* find average delta-p */
      
   delta_p = 0;  
   for (i = 1; i < hdr.nobs; ++i) 
      delta_p += p[i] - p[i-1];
   
    delta_p = delta_p / (hdr.nobs-1);
    
    filter_it = 0;
    
    if (delta_p < 2) {
       filter_it = 1;
       filter_width = NINT(prsint / delta_p) / 2 * 2 + 1;
    }
      
    pout = (double *) calloc(MAXOBS, sizeof(double));
    pout[0] = p[0];
    pout[1] = prsint;
    while (pout[1] <= pout[0] )
        pout[1] += prsint;
	
    nout = 1;
    maxp = p[hdr.nobs-1];
    
    while (pout[nout] < maxp) {
       lastp = pout[nout];
       pout[++nout] = lastp + prsint;
    }
    
    pout[nout] = maxp;
    ++nout;
    
    if (filter_it) {
       filter("g", t, hdr.nobs, filter_width);
       filter("g", s, hdr.nobs, filter_width);
       if (n_o == hdr.nobs )
          filter("g", o, hdr.nobs, filter_width);
       
    }
    
    /* interpolate arrays onto pressure series */
    
    for (j = 0; j < hdr.nprops; ++j) {
       if (hdr.prop_id[j] != (int)PR) {
       
          switch (hdr.prop_id[j]) {
	     case (int) TE:
	        yin = t;
		break;
	     case (int) DE:
	        yin = d;
		break;
	     case (int) SA:
	        yin = s;
		break;
	  
	     case (int) OX:
	     case (int) O2:
	        yin = o;
		break;

	     default:
	     
	        fprintf(stderr,"\nUnknown property passed to pressure_sort()\n");
		exit(1);
	  } /* end switch */
	  
          yout = (double *) calloc(nout, sizeof(double));
	  yout[0] = yin[0];
	  --nout;  /* for now to avoid unnecessary arithmetic */
          for (i = 1; i < nout; ++i) {
             yout[i] = hb_linterp(pout[i], p, yin, hdr.nobs);
	     if (yout[i] < -8.9)
	        yout[i] = HB_MISSING;
          }
	  yout[nout] = yin[hdr.nobs-1];
	  ++nout;  /* adjust back */
	 
	  free(yin);
     
          switch (hdr.prop_id[j]) {
	     case (int) TE:
		t = yout;
		break;
	     case (int) DE:
		d = yout;
		break;
	     case (int) SA:
		s = yout;
		break;
	  
	     case (int) OX:
	     case (int) O2:
		o = yout;
		break;

	  } /* end switch */
       } /* end if hdr.prop_id */
    } /* end for j */
    
    free(p);
    p = pout;
    return (nout);
} /* end pressure_sort() */
