/*  sac_convert.c
................................................................................
.   Reads SAC gridded standard level hydrographic data files
.   and creates HydroBase *.cdf gridded files.
.   
.   USAGE: sac_convert -O<outdir> -D<indir> 
...............................................................................
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "hydrobase.h"
#include "hydro_cdf.h"


#define NO_PREFILL 0
#define NSD  45      /* No of standard depths in SAC */
#define NPROPS 5    /* No of properties not including de */


float stdd[NSD] = { 0, 10, 20, 30, 40, 50, 75, 100, 125, 150, 
                  175, 200, 250, 300, 350, 400, 500, 600, 700, 800, 
		  900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 2000, 2250,
		  2500,2750, 3000, 3250, 3500, 3750, 4000, 4250, 4500, 4750,
		  5000, 5250, 5500, 5750, 6000 };
		  

float HB_f_mask = (float)HBMASK;       /* float values are used since the output */
float HB_f_empty = (float)HBEMPTY;     /* cdf values are of this type */

struct CDF_HDR cdf;

 /*  prototypes for locally defined functions */

void print_usage(char *);
FILE *openfile(char *, char *, char *);

int main (int argc, char **argv)
{
   int error, dummy, nread, any_data;
   int i, j, k, ip, ii, jj, nsta;
   int nsq, indx, msq10;
   int row, col, zlev;
   int nobs;
   int cdf_file;
   int include_counts;
   int prop_indx[NPROPS];
   FILE *infile;
   char *indir, *outdir;
   char filename[200];
   char line[200];
   float latmin, lonmin;
   float lat, lon, fdum;
   float *xout;
   float *bottom;
   float ***p;
   float ***t;
   float ***s;
   float ***o;
   float ***g;
   float ****data;
   double dlev;
   
 
/*  set these default values */

    error = 0;
    indir = "";
    outdir = "";
   
    
/*  parse the command line arguments */

   if (argc > 1) {
      for (i = 1; i < argc; i++) {
         if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'D':                   /* get input dir */
                        indir = &argv[i][2];
                        break;
               case 'O':                    /* get output dir  */
                        outdir = &argv[i][2];
                        break;
     	       case 'h':  
	           print_usage(argv[0]);
	           exit(0);

             default:
                        error = 1;

          }    /* end switch */

          if (error ) {
             fprintf(stderr,"\nError parsing command line args.\n");
             fprintf(stderr,"     in particular: '%s'\n", argv[i]);
             exit(1);
          }

        }  /* end if */
      }  /* end for */
   } /* end if argc > 1*/
   
/* allocate memory for global property arrays */

   fprintf(stderr,"Allocating memory...");

   p = (float ***) calloc((size_t)180, sizeof(float **));
   t = (float ***) calloc((size_t)180, sizeof(float **));
   s = (float ***) calloc((size_t)180, sizeof(float **));
   o = (float ***) calloc((size_t)180, sizeof(float **));
   g = (float ***) calloc((size_t)180, sizeof(float **));
   
   for (i = 0; i < 180; ++i)  {
      p[i] = (float **) calloc((size_t)360, sizeof(float *));
      t[i] = (float **) calloc((size_t)360, sizeof(float *));
      s[i] = (float **) calloc((size_t)360, sizeof(float *));
      o[i] = (float **) calloc((size_t)360, sizeof(float *));
      g[i] = (float **) calloc((size_t)360, sizeof(float *));
      
      for (j = 0; j < 360; ++j) {
        p[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
        t[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
        s[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
        o[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
        g[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
	
        if (g[i][j] == NULL) {
	 fprintf(stderr,"\nError allocating memory.\n");
	 exit(1);
	}
	
        for (k = 0; k < NSD; ++k) {
	  p[i][j][k] = HB_f_mask;
	  t[i][j][k] = HB_f_mask;
	  s[i][j][k] = HB_f_mask;
	  o[i][j][k] = HB_f_mask;
	  g[i][j][k] = HB_f_mask;
	}     
      }
   }
  
/* hardwire the properties that will be stored */

   prop_indx[0] = (int)PR;
   prop_indx[1] = (int)TE;
   prop_indx[2] = (int)SA;
   prop_indx[3] = (int)OX;
   prop_indx[4] = (int)GN;

   data = (float ****)calloc((size_t)NPROPS, sizeof(float ***));
   data[0] = p;   
   data[1] = t;
   data[2] = s;
   data[3] = o;
   data[4] = g;
  
/* initialize standard depths */

  NSTDLEVS = NSD + 1;
  for (zlev = 0; zlev < NSD; ++zlev)  /* loop for SAC std levs */
     std_depth[zlev] = (double) stdd[zlev];
     
  std_depth[zlev] = -99.0;   /* add the variable bottom depth */  
  std_depth_initialized = 1;

/*  open file and read gridded profiles  */
   
   infile = openfile(indir,"sac_wtso.asc", "");
   if (infile == NULL) 
      exit(1);

   nsq = 0;
   while (fscanf(infile,"%[^\n]", line) != EOF) {
     
     error = getc(infile);  /* move past LF */
      
     nread = sscanf(line,"%d", &nobs);
     nread = sscanf(&line[3],"%f%f%f%f%d", &lon, &lat, &fdum, &fdum, &dummy);
      if (nread != 5) {
	  fprintf(stderr,"\nError parsing header:\n%s\n", line);
	  exit(1);
      }
      
      if ((++nsq % 161) == 0) 
         fprintf(stderr,".");  /* show progress */
         
      if (nobs > 0) {
         row = NINT(lat - -90.0);
	 if (lon < 0)
	   lon += 360;
	 col = NINT(lon);
	 
	 for (i = 0; i < nobs; ++i) {
	    fscanf(infile,"%[^\n]", line);
	    error = getc(infile);   /* move past LF */
	    
	    nread = sscanf(line,"%lf", &dlev);
	    zlev = d2stdlev(dlev);
	    if (zlev < 0) {
	       fprintf(stderr,"\nError parsing observation level at lon:%.4f  lat:%.4f\n%s\n", lon, lat, line);
	       fprintf(stderr," d2stdlev(%.1lf) returned %d\n", dlev, zlev);
	       exit(1);
	    }
	    nread = sscanf(line,"%*lf %f %f %f %f", &g[row][col][zlev],
	     &t[row][col][zlev], &s[row][col][zlev], &o[row][col][zlev]);
	     
	    if (nread != 4) {
	       fprintf(stderr,"\nError parsing data values at lon:%.4f  lat:%.4f \n%s\n", lon, lat, line);
	       exit(1);
	    }
	    
	    if (g[row][col][zlev] < -9.) {
	        g[row][col][zlev] = HB_f_empty;
	        t[row][col][zlev] = HB_f_empty;
	        s[row][col][zlev] = HB_f_empty;
	        o[row][col][zlev] = HB_f_empty;
	    }
	    else {
	       p[row][col][zlev] = hb_p80(dlev, (double)lat);
	    }
	 } /* end for i */
      
      
      }
   
   } /* end while */
   
   fprintf(stderr,"\nFinished input.\n");
   
   fclose(infile);
   
   
/******************* End of Input Phase ************/

/* Construct a cdf header */

  cdf.nx = 10;
  cdf.ny = 10;
  cdf.nz = NSD+1;     /* includes a variable bottom depth */
  cdf.nprops = NPROPS;
  cdf.counts_included = 0;
  cdf.node_offset = 1;
  cdf.fill_value = HB_f_empty;
  cdf.tmin = (int *) calloc(1, sizeof(int)); 
  cdf.tmin[0] = 0;  
  cdf.tmax = (int *) calloc(1, sizeof(int));   
  cdf.tmax[0] = 9999;
  cdf.xincr = 1.0;
  cdf.yincr = 1.0;
  strncpy(cdf.x_units, "degrees", 8);
  strncpy(cdf.y_units, "degrees", 8);
  strncpy(cdf.z_units, "meters", 7);
  strncpy(cdf.title,"SAC :  1-degree gridded data", 10);
  strcpy(cdf.command, *argv);
  cdf.prop_id = (char **) malloc(cdf.nprops * sizeof(char *));
  cdf.prop_units = (char **) malloc(cdf.nprops * sizeof(char *));
  for (i = 0; i < cdf.nprops; ++i) {
      cdf.prop_id[i] = (char *) malloc(3);
      cdf.prop_units[i] = (char *) malloc(50);
      strncpy(cdf.prop_id[i], get_prop_mne(prop_indx[i]),3);
      strcpy(cdf.prop_units[i], get_prop_units(prop_indx[i]));
   }
  
  

/* Since everything is standard depths, mark the bottom depth empty */

  nsq = cdf.nx * cdf.ny;
  bottom = (float *) calloc((size_t)nsq, sizeof(float));
  
  for (i = 0; i < nsq; ++i)
     bottom[i] = HB_f_empty;
  
/* allocate memory to output a latitude band at a time */

  xout = (float *) calloc((size_t)(cdf.nz*cdf.nx), sizeof(float));

/*  loop for each output cdf file */

   latmin = -80.5;
   for (jj = 1; jj < 17; ++jj) {
   
      lonmin = -0.5;

      for (ii = 0; ii < 36; ++ii) {
      
        msq10 = ms10(latmin+1.0, lonmin+1.0, &dummy);
	if (*outdir != '\0') 
	   sprintf(filename, "%s/%4d.sac.cdf", outdir, msq10);
	else
	   sprintf(filename, "%4d.sac.cdf", msq10);
	   
        cdf_file = cdf_init(filename);
	fprintf(stderr,"\nOpened %s for output", filename); 
	   
        cdf.xmin = lonmin;
        cdf.xmax = lonmin + 10.0;
        cdf.ymin = latmin;
        cdf.ymax = latmin + 10.0;
	   
        error = cdf_define(cdf_file, &cdf, NO_PREFILL, cdf.counts_included);
        error = write_std_depths_cdf(cdf_file, &cdf);
        error = write_time_bins_cdf(cdf_file, &cdf);
	error = write_bottom_depth_cdf(cdf_file, 0, 0, 0, cdf.ny, cdf.nx, 1, bottom);
	
	
	/* file is now ready for writing data 
	
	  ip = index to property
	  [col][row] are indices to input (global) arrays   
	  i, j are indices to 10-deg files
	  indx is for the xout array.
	  first point in cdf_file is at {lonmin,latmax} 
	  so start at top row and work down.        */
	
	/* do each property  */
	
	any_data = 0;
        for (ip = 0; ip < cdf.nprops; ++ip) {
	   row = jj * 10 + 10; 
	             
	   for (j = 0; j < cdf.ny; ++j) {    
	      --row;
	      col = ii * 10;
	      indx = 0;
	      
	      for (i = 0; i < cdf.nx; ++i) {
	      
	         for (zlev = 0; zlev < NSD; ++zlev) {  /* standard levels */
		    xout[indx] = data[ip][row][col][zlev];
		    if (ip == 2 && !any_data){
		      if (xout[indx] < 40 && xout[indx] > 0)
		        any_data = 1;
		    }
		    ++indx;
		 } 
		 
		 xout[indx++] = HB_f_empty;   /* mark bottom level empty */
	         ++col;
	      }  /* end for i */
	      
	      write_prop_cdf(cdf_file, xout, cdf.prop_id[ip], j, 0, 0, 0, 1, cdf.nx, 1, cdf.nz);
	      
	   }  /* end for j */
	}  /* end for ip */
	
        cdf_close(cdf_file);
	
	lonmin += 10.0;
        if (lonmin >= 179.0)
           lonmin -= 360;
	
      } /* end for ii */
      
      latmin += 10.0;
   
   } /* end for jj */
    
  fprintf(stderr,"\n\nThat's it.... we're done.\n");
  exit(0);
        
} /* end main() */ 
  
/*****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nConverts SAC 1-deg analyzed fields");
   fprintf(stderr,"\nto HydroBase *.cdf files.  Output files are named by");
   fprintf(stderr,"\nto 4-digit WMOsq with extent .sac.cdf");
   fprintf(stderr,"\nUsage:  %s [-Ooutdir] [-D<dirname>] [-E<file_extent>] [-h]", program);
   fprintf(stderr,"\n-D : directory for input file sac_wtso.asc ");
   fprintf(stderr,"\n     (default is current directory) ex: -D../data/ ");
   fprintf(stderr,"\n-O : directory for output files.");
   fprintf(stderr,"\n-h : help...... prints this message. ");

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
         fprintf(stderr,"\n   opened %s ... ", st);
   
   return(infile);
   
}  /* end openfile() */

/*****************************************************************************/
