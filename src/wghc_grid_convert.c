/*  wghc_grid_convert.c
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
#define NPROPS 13    /* No of properties not including de */
#define NFILES 6    /* # of input files */


float stdd[NSD] = { 0, 10, 20, 30, 40, 50, 75, 100, 125, 150, 
                  175, 200, 250, 300, 350, 400, 500, 600, 700, 800, 
		  900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 2000, 2250,
		  2500,2750, 3000, 3250, 3500, 3750, 4000, 4250, 4500, 4750,
		  5000, 5250, 5500, 5750, 6000 };
		  

float HB_f_mask = (float)HBMASK;       /* float values are used since the output */
float HB_f_empty = (float)HBEMPTY;     /* cdf values are of this type */
char fname[NFILES][15];          /* hard-wired names of input files */
struct CDF_HDR cdf;
float grid_inc;

 /*  prototypes for locally defined functions */

void print_usage(char *);
FILE *openfile(char *, char *, char *);

int main (int argc, char **argv)
{
   int error, dummy, nread, any_data;
   int i, j, k, ip, ii, jj, nsta, ifile;
   int nsq, indx, msq10;
   int row, col, zlev;
   int nobs;
   int cdf_file;
   int include_counts;
   int prop_indx[NPROPS];
   int ***tcount;
   int ***ocount;
   int ***silcount;
   int ***phoscount;
   int ***nitcount;
   int ****count;
   FILE *infile;
   char *indir, *outdir;
   char filename[200];
   char line[2000];
   short *nout;
   float latmin, lonmin;
   float lat, lon, fdum, bdepth;
   float radius;
   float *xout;
   float *bottom;
   float ***p;
   float ***t, ***t_std;
   float ***s, ***s_std;
   float ***o, ***o_std;
   float ***sil, ***sil_std;
   float ***phos, ***phos_std;
   float ***nit, ***nit_std;
   float ****data;
   double dlev;
   
 
/*  set these default values */

    error = 0;
    grid_inc = 0.5;
    indir = "";
    outdir = "";
    strcpy(fname[0], "WGHC0000_0595");
    strcpy(fname[1], "WGHC0600_1195");
    strcpy(fname[2], "WGHC1200_1795");
    strcpy(fname[3], "WGHC1800_2395");
    strcpy(fname[4], "WGHC2400_2995");
    strcpy(fname[5], "WGHC3000_3595");

    
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

   p = (float ***) calloc((size_t)360, sizeof(float **));
   t = (float ***) calloc((size_t)360, sizeof(float **));
   s = (float ***) calloc((size_t)360, sizeof(float **));
   o = (float ***) calloc((size_t)360, sizeof(float **));
   sil = (float ***) calloc((size_t)360, sizeof(float **));
   phos = (float ***) calloc((size_t)360, sizeof(float **));
   nit = (float ***) calloc((size_t)360, sizeof(float **));
    t_std = (float ***) calloc((size_t)360, sizeof(float **));
   s_std = (float ***) calloc((size_t)360, sizeof(float **));
   o_std = (float ***) calloc((size_t)360, sizeof(float **));
   sil_std = (float ***) calloc((size_t)360, sizeof(float **));
   phos_std = (float ***) calloc((size_t)360, sizeof(float **));
   nit_std = (float ***) calloc((size_t)360, sizeof(float **));

   tcount = (int ***) calloc((size_t)360, sizeof(int **));
   ocount = (int ***) calloc((size_t)360, sizeof(int **));
   silcount = (int ***) calloc((size_t)360, sizeof(int **));
   phoscount = (int ***) calloc((size_t)360, sizeof(int **));
   nitcount = (int ***) calloc((size_t)360, sizeof(int **));
  
   for (i = 0; i < 360; ++i)  {
      p[i] = (float **) calloc((size_t)720, sizeof(float *));
      t[i] = (float **) calloc((size_t)720, sizeof(float *));
      s[i] = (float **) calloc((size_t)720, sizeof(float *));
      o[i] = (float **) calloc((size_t)720, sizeof(float *));
      sil[i] = (float **) calloc((size_t)720, sizeof(float *));
      phos[i] = (float **) calloc((size_t)720, sizeof(float *));
      nit[i] = (float **) calloc((size_t)720, sizeof(float *));
      t_std[i] = (float **) calloc((size_t)720, sizeof(float *));
      s_std[i] = (float **) calloc((size_t)720, sizeof(float *));
      o_std[i] = (float **) calloc((size_t)720, sizeof(float *));
      sil_std[i] = (float **) calloc((size_t)720, sizeof(float *));
      phos_std[i] = (float **) calloc((size_t)720, sizeof(float *));
      nit_std[i] = (float **) calloc((size_t)720, sizeof(float *));
      tcount[i] = (int **) calloc((size_t)720, sizeof(int *));
      ocount[i] = (int **) calloc((size_t)720, sizeof(int *));
      silcount[i] = (int **) calloc((size_t)720, sizeof(int *));
      phoscount[i] = (int **) calloc((size_t)720, sizeof(int *));
      nitcount[i] = (int **) calloc((size_t)720, sizeof(int *));
      
      for (j = 0; j < 720; ++j) {
        p[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
        t[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
        s[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
        o[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
        sil[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
        phos[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
        nit[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
        t_std[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
        s_std[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
        o_std[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
        sil_std[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
        phos_std[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
        nit_std[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
        tcount[i][j] = (int *) calloc((size_t)NSD, sizeof(int ));
        ocount[i][j] = (int *) calloc((size_t)NSD, sizeof(int ));
        silcount[i][j] = (int *) calloc((size_t)NSD, sizeof(int ));
        phoscount[i][j] = (int *) calloc((size_t)NSD, sizeof(int ));
        nitcount[i][j] = (int *) calloc((size_t)NSD, sizeof(int ));
	
        if (nit_std[i][j] == NULL) {
	 fprintf(stderr,"\nError allocating memory.\n");
	 exit(1);
	}
	
        for (k = 0; k < NSD; ++k) {
	  p[i][j][k] = HB_f_mask;
	  t[i][j][k] = HB_f_mask;
	  s[i][j][k] = HB_f_mask;
	  o[i][j][k] = HB_f_mask;
	  sil[i][j][k] = HB_f_mask;
	  phos[i][j][k] = HB_f_mask;
	  nit[i][j][k] = HB_f_mask;
	  t_std[i][j][k] = HB_f_mask;
	  s_std[i][j][k] = HB_f_mask;
	  o_std[i][j][k] = HB_f_mask;
	  sil_std[i][j][k] = HB_f_mask;
	  phos_std[i][j][k] = HB_f_mask;
	  nit_std[i][j][k] = HB_f_mask;
	}     
      }
   }
  
/* hardwire the properties that will be stored */

   prop_indx[0] = (int)PR;
   prop_indx[1] = (int)TE;
   prop_indx[2] = (int)SA;
   prop_indx[3] = (int)OX;
   prop_indx[4] = (int)SI;
   prop_indx[5] = (int)P4;
   prop_indx[6] = (int)N3;
   prop_indx[7] = (int)F1;   /* use tracer variables to store std dev props for now */
   prop_indx[8] = (int)F2;
   prop_indx[9] = (int)F3;
   prop_indx[10] = (int)HE;
   prop_indx[11] = (int)TU;
   prop_indx[12] = (int)N2;

   data = (float ****)calloc((size_t)NPROPS, sizeof(float ***));
   data[0] = p;   
   data[1] = t;
   data[2] = s;
   data[3] = o;
   data[4] = sil;
   data[5] = phos;
   data[6] = nit;
   data[7] = t_std;
   data[8] = s_std;
   data[9] = o_std;
   data[10] = sil_std;
   data[11] = phos_std;
   data[12] = nit_std;
   
   count = (int ****)calloc((size_t)NPROPS, sizeof(int ***));
   count[0] = tcount;   
   count[1] = tcount;
   count[2] = tcount;
   count[3] = ocount;
   count[4] = silcount;
   count[5] = phoscount;
   count[6] = nitcount;
   count[7] = tcount;   /* std dev variables */
   count[8] = tcount;
   count[9] = ocount;
   count[10] = silcount;
   count[11] = phoscount;
   count[12] = nitcount;
  
  
/* initialize standard depths */

  NSTDLEVS = NSD + 1;
  for (zlev = 0; zlev < NSD; ++zlev)  /* loop for WGHC std levs */
     std_depth[zlev] = (double) stdd[zlev];
     
  std_depth[zlev] = -99.0;   /* add the variable bottom depth */  
  std_depth_initialized = 1;

/*  open file and read gridded profiles  */
for (ifile = 0; ifile < NFILES; ++ifile) {   
   infile = openfile(indir,fname[ifile], "");
   if (infile == NULL) 
      exit(1);

   nsq = 0;
   while (fscanf(infile,"%[^\n]", line) != EOF) {
     
     error = getc(infile);  /* move past LF */
      
     nread = sscanf(line,"%d %f", &nobs, &radius);
     
     if ((++nsq % 101) == 0) 
           fprintf(stderr,".");      /*show progress */
         
     for (i = 0; i < nobs; ++i) {
	 fscanf(infile,"%[^\n]", line);
	 error = getc(infile);   /* move past LF */
	    
	 if (radius > 0) { 
	    nread = sscanf(line,"%f%f%f%lf", &lon, &lat, &bdepth, &dlev);
	    if (nread != 4) {
	       fprintf(stderr,"\nError parsing observation level at lon:%.4f  lat:%.4f\n%s\n", lon, lat, line);
	       exit(1);
	    }
		
	    if (i == 0) {
                row = NINT((lat - -90.0) / grid_inc);
	        col = NINT(lon / grid_inc);
	    }
	    
	    if (dlev >= 0  && row < 360)  {   /* skip North Pole and any non-standard depths */
	       zlev = d2stdlev(dlev);
	       if (zlev < 0)  {
	            fprintf(stderr,"ERROR:  d2stdlev(%.1lf) returned %d\n", dlev, zlev);
		    exit(1);
	        }
	   
	       nread = sscanf(line,"%*f%*f%*f%*lf%f%f%f%f%f%f%f%f%*f%*f%*f%*f%*f%d%d%d%d%d%f%f%f%f%f%f", &p[row][col][zlev], &t[row][col][zlev], &fdum, &s[row][col][zlev], &o[row][col][zlev], &sil[row][col][zlev], &nit[row][col][zlev], &phos[row][col][zlev],  &tcount[row][col][zlev],  &ocount[row][col][zlev], &silcount[row][col][zlev],  &nitcount[row][col][zlev], &phoscount[row][col][zlev],  &t_std[row][col][zlev], &s_std[row][col][zlev], &o_std[row][col][zlev], &sil_std[row][col][zlev], &nit_std[row][col][zlev], &phos_std[row][col][zlev]  );
	     
	       if (nread != 19) {
	          fprintf(stderr,"\nError parsing data values at lon:%.4f  lat:%.4f \n%s\n", lon, lat, line);
	          exit(1);
	       }
	    
	       if (p[row][col][zlev] < -8.) 
	           p[row][col][zlev] = HB_f_empty;
		   
	       if (t[row][col][zlev] < -8.) {
	           t[row][col][zlev] = HB_f_empty;
	           t_std[row][col][zlev] = HB_f_empty;
		}   
	       if (s[row][col][zlev] < -8.)  {
	           s[row][col][zlev] = HB_f_empty;
	           s_std[row][col][zlev] = HB_f_empty;
		}   
	       if (o[row][col][zlev] < -8.)  {
	           o[row][col][zlev] = HB_f_empty;
	           o_std[row][col][zlev] = HB_f_empty;
		}   
	       if (sil[row][col][zlev] < -8.)  {
	           sil[row][col][zlev] = HB_f_empty;
	           sil_std[row][col][zlev] = HB_f_empty;
		}   
	       if (nit[row][col][zlev] < -8.)  {
	           nit[row][col][zlev] = HB_f_empty;
	           nit_std[row][col][zlev] = HB_f_empty;
		}   
	       if (phos[row][col][zlev] < -8.)  {
	           phos[row][col][zlev] = HB_f_empty;
	           phos_std[row][col][zlev] = HB_f_empty;
		}
		
		if (tcount [row][col][zlev] <= 1) {
	           t_std[row][col][zlev] = HB_f_empty;
	           s_std[row][col][zlev] = HB_f_empty;
		}
		if (ocount [row][col][zlev] <= 1) {
	           o_std[row][col][zlev] = HB_f_empty;
		}   
		if (silcount [row][col][zlev] <= 1) {
	           sil_std[row][col][zlev] = HB_f_empty;
		}   
		if (nitcount [row][col][zlev] <= 1) {
	           nit_std[row][col][zlev] = HB_f_empty;
		}   
		if (phoscount [row][col][zlev] <= 1) {
	           phos_std[row][col][zlev] = HB_f_empty;
		}   
		
	       
	    } /* end dlev */
         } /* if radius */
      } /* end for i */
   } /* end while */

   fclose(infile);
} /*end for ifile */   
   fprintf(stderr,"\nFinished input.\n");
   
/******************* End of Input Phase ************/

/* Construct a cdf header */

  cdf.nx = 20;
  cdf.ny = 20;
  cdf.nz = NSD+1;     /* includes a variable bottom depth */
  cdf.nprops = NPROPS;
  cdf.counts_included = 1;
  cdf.node_offset = 1;
  cdf.fill_value = HB_f_empty;
  cdf.mask_value = HB_f_mask;
  cdf.tmin = (int *) calloc(1, sizeof(int)); 
  cdf.tmin[0] = 0;  
  cdf.tmax = (int *) calloc(1, sizeof(int));   
  cdf.tmax[0] = 9999;
  cdf.xincr = 0.5;
  cdf.yincr = 0.5;
  strncpy(cdf.x_units, "degrees", 8);
  strncpy(cdf.y_units, "degrees", 8);
  strncpy(cdf.z_units, "meters", 7);
  strncpy(cdf.title,"WGHC :  half-degree gridded data", 10);
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
  nout = (short *) calloc((size_t)(cdf.nz*cdf.nx), sizeof(short));

/*  loop for each output cdf file */

   latmin = -80.25;
   for (jj = 1; jj < 18; ++jj) {
   
      lonmin = -0.25;

      for (ii = 0; ii < 36; ++ii) {
      
        msq10 = ms10(latmin+1.0, lonmin+1.0, &dummy);
	if (*outdir != '\0') 
	   sprintf(filename, "%s/%4d.wghc.cdf", outdir, msq10);
	else
	   sprintf(filename, "%4d.wghc.cdf", msq10);
	   
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
	   row = jj * cdf.ny + cdf.ny; 
	             
	   for (j = 0; j < cdf.ny; ++j) {    
	      --row;
	      col = ii * cdf.nx;
	      indx = 0;
	      
	      for (i = 0; i < cdf.nx; ++i) {
	      
	         for (zlev = 0; zlev < NSD; ++zlev) {  /* standard levels */
		    xout[indx] = data[ip][row][col][zlev];
		    nout[indx] = count[ip][row][col][zlev];
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
	      write_prop_count_cdf(cdf_file, nout, cdf.prop_id[ip], j, 0, 0, 0, 1, cdf.nx, 1, cdf.nz);
	      
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
   fprintf(stderr,"\nConverts WGHC half-deg analyzed fields");
   fprintf(stderr,"\nto HydroBase *.cdf files.  Output files are named by");
   fprintf(stderr,"\nto 4-digit WMOsq with extent .wghc.cdf");
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
