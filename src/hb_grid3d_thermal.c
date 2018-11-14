/*  hb_grid3d_thermal.c
................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             June 2005
................................................................................
..........................................................................
. Computes an average potential temperature profile for each property specified at each lat/lon
. gridnode from randomly spaced stations and outputs a netCDF file of potential temperature and pressure
. gridded as a function of lat, lon, and depth.  The depth
. dimension consists of a series of standard depths, but these values also
. can be optionally specified by the user with -Z.  The averaging is done 
. on isothermal surfaces and interpolated back onto the depth levels.  
..........................................................................
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hydro_cdf.h"
#include "hydrobase.h"


#define    UI   unsigned int
#define    PREFILL      0        /* turn off prefilling of cdf variables */
#define    ADD_PROP_COUNTS  1    /* turn on option to include # of obs info */
#define    PRINT_MSG  1          /* turn on message printing to stderr */
#define    MAX_TIME_BINS  1      /* anachronistic time dimension stores just one matrix  */
#define    MIX_LAYER_DEF 0.1    /* default temperature range for mixed layer */
#define    MAX_CROSSINGS  5     /* max number of crossings of each theta level */

/* input_file pathnames */

#define    EXTENT   ""
#define    DIR      ""

/***************************************************************/
/*  define the standard theta levels on which data will be averaged...*/

#define  MAXTLEVS 800    /* max size of temperature series -- arbitrarily large number */

double tlevs[MAXTLEVS];   /* stores the theta series  */
int  ntlevs;                /* number of elements in full theta series */


/***************************************************************/

/* globally referenced variables for input station data */

struct HYDRO_DATA sta;
struct HYDRO_HDR hdr;
double *pr, *de, *te, *th, *sa;

double mix_layer_def;      /* option to specify mixed layer definition */

/************ structures to represent grid node ****************/

struct surfrec {
          double  sst;
          double  depth;
          double  tavg;
              UI  n;
  struct surfrec *next;
};

struct deepestrec {
          double depth, pressure;
	  double theta;
          struct deepestrec *next;
};

struct gridnode {
         double **d;
	 int **gradient;
         UI  **nobs;
 struct deepestrec *deepest;
 struct surfrec  *mix_layer;
};
/***************************************************************/
/*  prototypes for functions defined within this module... */

struct gridnode **alloc_grid(int, int, int);

int get_theta_series(int, FILE *,  double *);
int get_time_bins(int, int, struct CDF_HDR *);
int do_std_depth(struct gridnode *, int,  int, int, float **, short **, float *);
int mixed_layer_vals(int, int, struct surfrec *, double *, double *, double *, int, int);

int monotonic(double **, UI **, int);

int sort_by_depth(double *, UI *, double *, double *, int);
void quicksort(double *, int *, int);
int find_pivot(double *, int, double *);
int partition(double *, int *, int, double);
void swap_d(double *, double *);
void swap_i(int *, int *);

void free_grid(int, int,  struct gridnode **);
void insert_data(double *,int, double **,UI **, int **, int);
void insert_surf_vals(struct gridnode *);
void insert_deepest_vals(struct gridnode *);
void compute_avg(double **,UI **, int);
void sum_levels(double *, double *, int, double *, UI *);
void print_usage(char *);

struct surfrec *get_surfrec(struct surfrec *, double);
struct surfrec *create_surfrec();
struct surfrec *get_sst_rec(int, struct surfrec *);
struct surfrec *define_avg_mixed_layer(struct surfrec *);
struct deepestrec *create_deepestrec();

double interpolate(double, double *, double *, int);


/***************************************************************/

int main (int argc, char **argv)
{
   short   bopt, iopt, popt, copt;
   int     curfile = 1, nfiles = 0; 
   int     print_mess, xdateline; 
   int     nlevs, nprops, prop_indx[MAXPROP];
   int     i, j, error, test_outcrop, prop_avail;
   int     nstdlevs, tmin, tmax, n_filled;
   int     include_counts;
   FILE   *tlev_file;
   FILE   *z_file;
   int     row, col, nrows, ncols, tbin, ntbins;
   int     out_of_bounds=0;
   char   *extent, *dir, *st;
   char   *cdf_filename;
   int     infile, cdf_file;
   float  **data, *bottomdepth;  
   short  **count; 
   struct CDF_HDR  h;
   struct gridnode **grid[MAX_TIME_BINS];


/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
 
/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    bopt = iopt  = copt = 0;
    z_file = tlev_file = NULL;
    include_counts = ADD_PROP_COUNTS;
    error = 0;
    print_mess = PRINT_MSG;
    xdateline = 0;
    mix_layer_def = MIX_LAYER_DEF;
    tmin = 0;
    tmax = 9999;
    
/*  initialize default theta series... */
    ntlevs = 0;
    tlevs[ntlevs] = -2.0;
    while (tlevs[ntlevs++] < 40.0){
      tlevs[ntlevs] = tlevs[ntlevs-1] + 0.1;
    }
    
   
/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'B':                    /* get grid bounds */
                        bopt = 1;
                        st = &argv[i][2];
                           if (*st == '/')
                               ++st;
                        error = (sscanf(st,"%f", &h.xmin) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%f", &h.xmax) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%f", &h.ymin) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%f", &h.ymax) != 1);
                        
                        if (h.xmin > 0 && h.xmax < 0)
                           h.xmax += 360.;
                        if (h.xmax > 180)
                           xdateline = 1;
                           
                        break;

               case 'C':
                        cdf_file = cdf_init(&argv[i][2]);
			cdf_filename = &argv[i][2];
                        copt = 1;
                        break;

               case 'D':                   /* get input directory */
                        dir = &argv[i][2];
                        break;
               case 'E':                    /* get file extent */
                        extent = &argv[i][2];
                        break;

               case 'I':
                        iopt = 1;
                        error = (sscanf(&argv[i][2],"%f", &h.xincr) == 1) ? 0 : 1;
                        h.yincr = h.xincr;
                        st = &argv[i][2];
                        while (*(st++) != '\0') {
                            if (*st == '/') {
                              ++st;
                              error = (sscanf(st,"%f", &h.yincr) == 1) ? 0 : 1;
                              break;
                            }
                        }
                        
                        break;
                        
               case 'M':
                        error = (sscanf(&argv[i][2],"%lf", &mix_layer_def) == 1) ? 0 : 1;
                        break;
                        

               case 'T':
                        
                        if (argv[i][2] == '\0') {
                           fprintf(stdout,"\nStandard theta series: \n");
                           for (j = 0; j < ntlevs; ++j) {
                              fprintf(stdout,"  %.1lf", tlevs[j]);
                           }
                           fprintf(stdout,"\n");
                           exit(0);
                        }
			tlev_file = fopen(&argv[i][2],"r");
                        if (tlev_file == NULL) {
                           fprintf(stderr,"\nError opening theta-series file: %s\n",&argv[i][2]);
                           exit(1);
                        }
                        break;
               case 'Z':
                        if (argv[i][2] == '\0') {
                           nstdlevs = std_depth_init(z_file);
                           fprintf(stdout,"\nStandard depth levels: \n");
                           for (j = 0; j < nstdlevs; ++j) {
                              fprintf(stdout,"  %.1lf", std_depth[j]);
                           }
                           fprintf(stdout,"\n");
                           exit(0);
                        }
                        z_file = fopen(&argv[i][2],"r");
                        if (z_file == NULL) {
                           fprintf(stderr,"\nError opening stddep file: %s\n",&argv[i][2]);
                           exit(1);
                        }
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
             print_usage(argv[0]);
             exit(1);
          }

       }  /* end if */

       else  {
           ++nfiles;
       }
   }  /* end for */

   if (!bopt || !iopt || !nfiles ||  !copt) {
       fprintf(stderr,"\nYou must specify input file(s), bounds, cdf_output_file and gridspacing!\n");
       print_usage(argv[0]);
       exit(1);
   }


/* initialize global array of std_depths */

     nstdlevs = std_depth_init(z_file);


/* get optionally specified values for theta series ... */

   if (tlev_file == NULL) {
        fprintf(stderr,"\n Using standard theta levels ranging %.2ld to %.2ld", tlevs[0], tlevs[ntlevs-1]);
   }
   else {
     ntlevs = get_theta_series(tlev_file, tlevs);
      if (ntlevs >= MAXTLEVS) {
          fprintf(stderr,"\n # of theta levels exceeds space allocated");
          fprintf(stderr,"\n Recompile program with larger MAXTLEVS.\n");
          exit(1);
      }
     fprintf(stderr,"\n total number of theta levels: %d\n", ntlevs);
   }


/* compute dimensions of grid and allocate space for computation ...*/

   nrows = NINT( ceil((double)((h.ymax - h.ymin) / h.yincr)));
   ncols = NINT(ceil((double)((h.xmax - h.xmin) / h.xincr)));
   h.node_offset = 1;       /* 0 for node grid, 1 for pixel grid */

/*  Time bins are an anachronism (part of the header record) so just do it ... */

   ntbins = get_time_bins(tmin, tmax, &h);

/*  set these cdf header values */
   h.nx = ncols;
   h.ny = nrows;
   h.nz = nstdlevs;
   h.nt = ntbins;

   fprintf(stderr,"\n allocating memory for grid ... " );

   for (tbin = 0; tbin < ntbins; ++tbin) {
     grid[tbin] = alloc_grid(ncols, nrows, ntlevs);
   }
   fprintf(stderr,"\n");

/* initialize these... */

   for (i = 0; i < MAXPROP; ++i)
      sta.observ[i] = (double *) NULL;
   hdr.prop_id = (int *) NULL;

/* loop for each input_file */

   do {
     if ((infile = open_hydro_file(dir, argv[curfile], extent, print_mess)) 
          < 0) {
       goto NEXTFILE;
     }

     
     /* loop for each station */

    while ((error = get_station(infile, &hdr, &sta)) == 0) {

       /* is station within bounds ?   */
       
       if (xdateline && hdr.lon < 0)
          hdr.lon += 360;

       if (get_indices(&h, hdr.lat, hdr.lon, &row, &col) < 0) {
           fprintf(stderr,"Out of bounds: %.3f %.3f\n", hdr.lat, hdr.lon);
           ++out_of_bounds;
       }
       else {

         pr = sta.observ[(int)PR];   /* set frequently used pointers */
         de = sta.observ[(int)DE];
         te = sta.observ[(int)TE];
         th = sta.observ[(int)TH];
         sa = sta.observ[(int)SA];

       /* compute pr, de, or th if not available ... */
       	 
         if (pr == NULL && de == NULL) {
             fprintf(stderr, "\nNo pressure or depth available at this station.\n");
	      continue;
	 }
	 if ((th == NULL ) && ( te == NULL || sa == NULL)){
           fprintf(stderr, "\nUnable to compute theta at this station.\n");
           continue;
	 }
	 
	   
	 if (pr == NULL) {
           free_and_alloc(&sta.observ[(int)PR], hdr.nobs);
	   for (i = 0; i < hdr.nobs; ++i) {
	       pr[i] = hb_p80(de[i], (double)hdr.lat);
	   }    
	 }
	 
         if (de == NULL) {
           free_and_alloc(&sta.observ[(int)DE], hdr.nobs);
	   for (i = 0; i < hdr.nobs; ++i) {
	       de[i] = hb_depth(pr[i], (double)hdr.lat);
	   }    
	 }
	 
	 if (th == NULL) {
           free_and_alloc(&sta.observ[(int)TH], hdr.nobs);
           compute_theta(hdr.nobs, sta.observ[(int)TH], pr, te, sa);
         }
	  
         tbin = 0;
         test_outcrop =  0;
         insert_data(de, hdr.nobs, grid[tbin][row][col].d, grid[tbin][row][col].nobs, grid[tbin][row][col].gradient, test_outcrop);
         insert_surf_vals(&grid[tbin][row][col]);
         insert_deepest_vals(&grid[tbin][row][col]);
	 
       } /* end else */
       

       /* free up space... */

       for (i = 0; i < MAXPROP; ++i) {
          if (sta.observ[i] != NULL) {
             free((void *) sta.observ[i]);
             sta.observ[i] = (double *) NULL;
          }
       }

     }  /*end while !eof */ 
     if (error > 0) {
         report_status(error, stderr);
         exit(1);
   }

NEXTFILE:
     close(infile);

   } while (curfile++ < nfiles );             


/******** end of first input phase *********/

   if (out_of_bounds) {
      fprintf(stderr,"\n  %d out-of-bounds stations were read and ignored.\n", out_of_bounds);
   }

   fprintf(stderr,"\n  constructing header ...\n");


/* construct the cdf header  */

   h.nprops = nprops = 2;   /*just pr and th */
   prop_indx[0] = (int) PR;
   prop_indx[1] = (int) TH;
   
   h.fill_value = HBEMPTY;
   h.mask_value = HBMASK;
   strncpy(h.x_units, "degrees", 8);
   strncpy(h.y_units, "degrees", 8);
   strncpy(h.z_units, "meters", 7);
   strncpy(h.title,"HydroBase", 10);
   strcpy(h.command, "hb_grid3d_thermal");

   h.prop_id = (char **) malloc(nprops * sizeof(char *));
   h.prop_units = (char **) malloc(nprops * sizeof(char *));
   for (i = 0; i < nprops; ++i) {
      h.prop_id[i] = (char *) malloc(3);
      h.prop_units[i] = (char *) malloc(50);
      strncpy(h.prop_id[i], get_prop_mne(prop_indx[i]),3);
      strcpy(h.prop_units[i], get_prop_units(prop_indx[i]));
   }
   if ( h.prop_units[nprops-1] == NULL) {
      fprintf(stderr, "\nUnable to malloc memory for the cdf header.\n");
      exit(1);
   }

   error = cdf_define(cdf_file, &h, PREFILL, 1 );
   error = write_std_depths_cdf(cdf_file, &h);
   error = write_time_bins_cdf(cdf_file, &h);


   data = (float **) malloc ((UI) nprops  * sizeof(float *));
   count = (short **) malloc ((UI) (nprops ) * sizeof(short *));
   bottomdepth = (float *) malloc ((UI) ncols * sizeof(float));

   for (i = 0; i < nprops; ++i) {
      data[i] = (float *) malloc((UI) (nstdlevs * ncols) * sizeof(float));
      count[i] = (short *) malloc((UI) (nstdlevs * ncols) * sizeof(short));
   }

   if (data[nprops-1] == NULL) {
      fprintf(stderr, "\nUnable to allocate memory for data & count arrays.\n");
      exit(1);
   }

   fprintf(stderr,"  computing averages ...\n");

/* for each gridnode, compute means at all sigma levels ... */

   for (row = 0; row < nrows; ++row) {
      for (col = 0; col < ncols; ++col) {
         for (i = 0; i < ntlevs; ++i) {
            compute_avg(grid[tbin][row][col].d[i], grid[tbin][row][col].nobs[i], MAX_CROSSINGS);
	 }
      }
   }


/* interpolate the theta series back onto the standard depth levels and 
   output the property data to the netcdf file... */

   fprintf(stderr,"  writing data ...\n");

   for (tbin = 0; tbin < ntbins; ++tbin) {
      for (row = 0; row < nrows; ++row) {
         for (col = 0; col < ncols; ++col) {
            get_lat_lon(&h, row, col, &latitude, &longitude);
            nlevs = do_std_depth(&grid[tbin][row][col], ntlevs, 
                  ncols, col, data, count, &bottomdepth[col]);
         }

         for (i = 0; i < nprops; ++i) {
            write_prop_cdf(cdf_file, data[i], h.prop_id[i], row, 0, tbin, 0,
                           1, ncols, 1, nlevs);
            write_prop_count_cdf(cdf_file, count[i], h.prop_id[i],
                           row, 0, tbin, 0, 1, ncols, 1, nlevs);
         }
         write_bottom_depth_cdf(cdf_file, row, 0, tbin, 1, ncols, 1,
                                bottomdepth);
      }
   }
   cdf_close(cdf_file);
   free((void *) data);
   free((void *) count);
   free((void *) bottomdepth);
   for (tbin = 0; tbin < ntbins; ++tbin) 
      free_grid(ncols, nrows, grid[tbin]);
   
/*  Now that memory has been freed up, search for vertical data
    gaps in isopycnally averaged grid and determine which result from
    pycnostads.  
    
    Strategy:
    Re-read all the station files again, this time averaging data
    on stddepth surfaces.  Then re-open the cdf file, visit each
    grid node and see if there are vertical data gaps which may
    have resulted from weak vertical gradients.  Replace empty levels
    with isobarically averaged datapoints if they exist. Rewrite the
    updated property values to the cdf file. */
    
    
    
   fprintf(stderr,"\nChecking for pycnostads....");
   fprintf(stderr,"\n    allocating memory for grid ... " );

   for (tbin = 0; tbin < ntbins; ++tbin) 
     grid[tbin] = alloc_grid(ncols, nrows, nstdlevs, nprops);
   

/* initialize these... */

   for (i = 0; i < MAXPROP; ++i)
      sta.observ[i] = (double *) NULL;
   hdr.prop_id = (int *) NULL;

/* loop for each input_file */

   fprintf(stderr,"\n    summing ...");
   curfile = 1;
   
   do {
     if ((infile = open_hydro_file(dir, argv[curfile], extent, print_mess)) 
          < 0) {
       goto NEXTFILE2;
     }

      
     /* loop for each station */

     while ((error = get_station(infile, &hdr, &sta)) == 0) {

       /* is station within bounds ?   */
       
       if (xdateline && hdr.lon < 0)
          hdr.lon += 360;

       if (get_indices(&h, hdr.lat, hdr.lon, &row, &col) == 0) {

         pr = sta.observ[(int)PR];   /* set frequently used pointers */
         de = sta.observ[(int)DE];
         te = sta.observ[(int)TE];
         th = sta.observ[(int)TH];
         sa = sta.observ[(int)SA];
	 
       /* compute pr, de, or th if not available ... */
       	 
         if (pr == NULL && de == NULL) {
             fprintf(stderr, "\nNo pressure or depth available at this station.\n");
	      continue;
	 }
	 if ((th == NULL ) && ( te == NULL || sa == NULL)){
           fprintf(stderr, "\nUnable to compute theta at this station.\n");
           continue;
	 }
	 
	   
	 if (pr == NULL) {
           free_and_alloc(&sta.observ[(int)PR], hdr.nobs);
	   for (i = 0; i < hdr.nobs; ++i) {
	       pr[i] = hb_p80(de[i], (double)hdr.lat);
	   }    
	 }
	 
         if (de == NULL) {
           free_and_alloc(&sta.observ[(int)DE], hdr.nobs);
	   for (i = 0; i < hdr.nobs; ++i) {
	       de[i] = hb_depth(pr[i], (double)hdr.lat);
	   }    
	 }
	 
	 if (th == NULL) {
           free_and_alloc(&sta.observ[(int)TH], hdr.nobs);
           compute_theta(hdr.nobs, sta.observ[(int)TH], pr, te, sa);
         }
	 
	 

         sum_levels(sta.observ[int(PR)], de, hdr.nobs, grid[tbin][row][col].prop[0], grid[tbin][row][col].count[0]);
         sum_levels(sta.observ[int(TH)], de, hdr.nobs, grid[tbin][row][col].prop[1], grid[tbin][row][col].count[1]);
	 
       } /* end if */

       /* free up space... */

       for (i = 0; i < MAXPROP; ++i) {
          if (sta.observ[i] != NULL) {
             free((void *) sta.observ[i]);
             sta.observ[i] = (double *) NULL;
          }
       }

     }  /*end while !eof */ 
     if (error > 0) {
         report_status(error, stderr);
         exit(1);
   }

NEXTFILE2:
     close(infile);

   } while (curfile++ < nfiles );             

/******** end of second input phase *********/

   cdf_file = cdf_update("", cdf_filename, "", print_mess);
   
   data = (float **) calloc ((size_t) nprops, sizeof(float *));
   count = (short **) calloc ((size_t) nprops, sizeof(short *));
   bottomdepth = (float *) malloc (sizeof(float));

   for (i = 0; i < nprops; ++i) {
      data[i] = (float *) calloc((size_t) nstdlevs, sizeof(float));
      count[i] = (short *) calloc((size_t) nstdlevs, sizeof(short));
   }
   
   nlevs = nstdlevs - 1;

   fprintf(stderr,"\n  filling thermostads ...");
   n_filled = 0;

/* check each thermally averaged gridnode.
   If depth level is flagged as missing, check
   the isobarically averaged data to differentiate a vertical datagap
   from a thermostad ... */

   
      for (row = 0; row < nrows; ++row) {
         for (col = 0; col < ncols; ++col) {
	    read_cdf_bottom_depth(cdf_file, bottomdepth, row, col, tbin);
            for (i = 0; i < nprops; ++i) {
               read_cdf_prop(cdf_file, h.prop_id[i], data[i], row, col, tbin, 0, h.nz);
	       read_cdf_prop_count(cdf_file, h.prop_id[i], count[i], row, col, tbin, 0, h.nz);
	       j = 0;
	       while (j < nlevs && std_depth[j] < *bottomdepth) {
		  if (is_flagged(data[i][j], (float)HBEMPTY)) {
		     if (grid[tbin][row][col].count[i][j] > 0) {
		       count[i][j] = (short) grid[tbin][row][col].count[i][j];
		       data[i][j] = (float) (grid[tbin][row][col].prop[i][j] / (double) count[i][j]);
		       if (i == 0) ++n_filled;
		     }
		  }
		  ++j;
	       }
               write_prop_cdf(cdf_file, data[i], h.prop_id[i], row, col, tbin, 0,
                           1, 1, 1, nlevs);
               write_prop_count_cdf(cdf_file, count[i], h.prop_id[i],
                           row, col, tbin, 0, 1, 1, 1, nlevs);
	    } /* end for i */
	 } /* end for col */
      } /* end for row */
 
   fprintf(stderr,"\n  %d levels were identified and filled", n_filled);
   fprintf(stderr,"\n End of %s\n", argv[0]);
  
   exit(0);
} /* end main */


/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nComputes an average profile for pressure and potential temperature");
   fprintf(stderr,"\nat each lat/lon gridnode from randomly spaced stations.");
   fprintf(stderr,"\nOutputs a netCDF file of hydro properties gridded as a ");
   fprintf(stderr,"\nfunction of lat, lon, and depth.  The depth dimension ");
   fprintf(stderr,"\nconsists of a series of standard depths, optionally ");
   fprintf(stderr,"\nspecified with -Z. Averaging is done on isotherms (delta=0.1degC)");
   fprintf(stderr,"\nand interpolated back onto the depth levels.  The user");
   fprintf(stderr,"\nmay specify different theta surfaces with -T<filename>.");
   fprintf(stderr,"\n");
   
   fprintf(stderr,"\nUsage:  %s filename_root(s)  -B/west/east/south/north -C<cdf_output_file> -I<gridspacing>  [-Z<std_depth_file>] [-T<theta_series_file>] [-M<mixed_layer_def>]  [-D<dirname>] [-E<file_extent>] \n", program);
   fprintf(stderr,"\n -B   specifies grid bounds");
   fprintf(stderr,"\n -C   name of netCDF output file.");

   fprintf(stderr,"\n -I   specifies grid increment in degrees;  ex: -I0.5");
   fprintf(stderr,"\n          OR specify separate x,y increments with a slash");
   fprintf(stderr,"\n          to delimit xincr/yincr;  ex: -I2.0/0.5\n");
   fprintf(stderr,"\n\n OPTIONS:");
   fprintf(stderr,"\n -D  specifies directory for input data files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/\n ");
   fprintf(stderr,"\n -E  specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat \n");
   fprintf(stderr,"\n -M  specify the definition of the mixed layer as ");
   fprintf(stderr,"\n          theta(sea surface) + this value.  default: %f ", MIX_LAYER_DEF);
   fprintf(stderr,"\n          ex:  -M.02");
   fprintf(stderr,"\n -T  filename containing optionally specified theta levels for averaging on.");
   fprintf(stderr,"\n -Z  file containing list of standard depths.");
   fprintf(stderr,"\n          Values MUST be monotonically INCREASING.\n");
   fprintf(stderr,"\n -h  help ... prints this message.");
   fprintf(stderr,"\n\n");  
   return;
} /* end print_usage() */
/****************************************************************************/

/*****************************************************************************/
struct gridnode **alloc_grid(int nx, int ny, int nz)
{
   int i, j,  n;
   struct gridnode **g;

   g = (struct gridnode **) malloc(ny * sizeof(struct gridnode *));
   if (g == NULL) {
      fprintf(stderr,"\nUnable to allocate memory for grid.\n");
      exit(1);
   }
   for (i = 0; i < ny; ++i) {
         g[i] = (struct gridnode *) malloc(nx * sizeof(struct gridnode));
         if (g[i] == NULL) {
           fprintf(stderr,"\nUnable to allocate memory for grid[%d]\n", i);
           exit(1);
         }
         for (j = 0; j < nx; ++j) {

            }
            g[i][j].d = (double **) malloc(nz * sizeof(double *));
            if (g[i][j].d == NULL) {
               fprintf(stderr,"\nUnable to allocate memory for g[%d][%d].d\n", i, j);
               exit(1);
            }

            g[i][j].nobs = (UI **) malloc(nz * sizeof(UI *));
            if (g[i][j].nobs == NULL) {
               fprintf(stderr,"\nUnable to allocate memory for g[%d][%d].nobs\n", i, j);
               exit(1);
            }

            g[i][j].gradient = (int **) malloc(nz * sizeof(int *));
            if (g[i][j].gradient == NULL) {
               fprintf(stderr,"\nUnable to allocate memory for g[%d][%d].gradient\n", i, j);
               exit(1);
            }
	    
	    for (n = 0; n < nz; ++n ) {
	    
	       g[i][j].d[n] = (double *) malloc(MAX_CROSSINGS * sizeof(double));
	       if (g[i][j].d[n] == NULL) {
                  fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].d\n", i, j);
                  exit(1);
	      }
	    

              g[i][j].nobs[n] = (UI *) malloc(MAX_CROSSINGS * sizeof(UI));
                if (g[i][j].nobs[n] == NULL) {
                   fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].nobs\n", i, j);
                   exit(1);
              }

              g[i][j].gradient[n] = (int *) malloc(MAX_CROSSINGS * sizeof(int));
                if (g[i][j].gradient[n] == NULL) {
                   fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].gradient\n", i, j);
                   exit(1);
              }
	    
	    }
	    
            g[i][j].mix_layer = (struct surfrec *) NULL;
            g[i][j].deepest = (struct deepestrec *)NULL;

            for (n = 0; n < nz; ++n) {       /* initialize  arrays */
	       for (ix = 0; i < MAX_CROSSINGS; ++ix) {
                   g[i][j].d[n][ix] = 0.0;
                   g[i][j].nobs[n][ix] = 0;
                   g[i][j].gradient[n][ix] = 0;
	       }   /* end for ix */
            } /* end for n */

         } /* end for j*/
   }  /* end for i */

   return g;
} /* end alloc_grid */
/*****************************************************************************/
void free_grid(int nx, int ny, struct gridnode **gridptr)
  /* frees all the memory for ny * nx struct gridnodes */
{
   int row, col, k;
   struct gridnode *g;
   
   for (row = 0; row < ny; ++row) {
      for (col = 0; col < nx; ++col) {
         g = &gridptr[row][col];
	 for (k = 0; k < ntlevs; ++k) {
	    free((void *)g->d[k]);
	    free((void *)g->nobs[k]);
	    free((void *)g->gradient[k]);
	 }
	 free((void *)g->nobs);
	 free((void *)g->d);
	 free((void *)g->gradient);
	 if (g->mix_layer != NULL)
	   free((void *)g->mix_layer);
	 if (g->deepest != NULL)
	   free((void *)g->deepest);
      } /* end for col */
      
      free((void *) gridptr[row]);
   } /* end for row */
   
   free((void *) gridptr);
   
   return;
}  /* end free_grid() */
/*****************************************************************************/
int get_theta_series(FILE *fptr, double *thetalist)
  
    /*  the file will be read for (min, max, incr) triplets
       from which a sigma series will be generated and inserted at siglist.
       The number of sigma values is returned. 
       
       arguments:
       ref_id:   defines the sigma variable: (int) enum property 
         fptr:   pointer to file containing list of sigma values OR nil 
      thetalist:   starting addr of array to insert values 
   */
{
   double  min, max, incr;
   double  next;
   int i, count;


     count = 0;
     while( fscanf(fptr,"%lf %lf %lf", &min, &max, &incr) == 3) {
         if (++count >= MAXTLEVS) {
           return (count);
         }
         *thetalist = min;
         while ( ( next = *thetalist + incr) <= max) {
           if (++count >= MAXTLEVS) { /* avoid a SEG FAULT */
              return (count);
           }
           *(++thetalist) = next; 
         }
         ++thetalist;
     }

     fclose(fptr);
     return (count);

} /* end get_theta_series */
/*****************************************************************************/
int get_time_bins(int minyr, int maxyr, struct CDF_HDR *hptr)
    /* Sets time bin min/max values in the CDF header */
{

   hptr->tmin = (int *) malloc( sizeof(int));
   hptr->tmax = (int *) malloc( sizeof(int));

   hptr->tmin[0] = minyr;        
   hptr->tmax[0] = maxyr;

   return (1);

 } /* end get_time_bins */

/*****************************************************************************/
void insert_data(double *y, int ny, double **ysum, UI **count,int **tgradient, int test_outcrop)

/* Interpolate
   to find the depth (y) at each level in the theta series,  add each
   yval to its appropriate sum, and increment the counter. Handles multiple crossings of an isotherm.
   Checks whether theta is monotonic, then for each crossing finds the depth level and whether or not theta is increasing or decreasing in the immediate depth range.  Adds information to the array with the closest depth level if multiple crossings occur in profile.  Checks for vertical
   data gaps using the globally defined de array and does not interpolate over
   gaps which exceed 200 m in the thermocline (upper 1000 m) or 1000 m 
   elsewhere.  The theta series is specified in the global array tlevs.
   
           y:    array of  observed depth values
          ny:    dimension of y 
        ysum:    sum of yvalues on each theta surface
       count:    counts # of observations on each surface 
   tgradient:   -1 if theta is decreasing, 1 if increasing, 0 if max/min point
test_outcrop:    0 or 1: if set, tests for outcropping surfaces
*/
{
   int  i, j, datagap, is_monotonic, n;
   double *xtmp,  *dtmp, z;
   double reflev;
   
   xtmp = (double *) malloc((UI) (ny * sizeof(double)));
   dtmp = (double *) malloc((UI) (ny * sizeof(double)));
   if (dtmp == NULL) {
       fprintf(stderr,"\nUnable to allocate memory in insert_data()\n");
       exit(1);
   }

/* Ensure continuous (no missing values) x,  y, and depth arrays */
   n = 0;
   for (i = 0; i < ny; ++i) {
      if ( y[i] > -8.9)  {
         xtmp[n] = th[i];
         dtmp[n++] = y[i];
      }
   }
   
   if (n <= 1) {               /* not enough values */
     free((void *)xtmp);
     free((void *)dtmp);
     return;
   }


   for (i = 0; i < ntlevs; ++i) {


        /* see if temperature level exists at this station ... */

      if ((z = interpolate(tlevs[i], xtmp, dtmp, n)) > -998.) {

            /* now check for datagaps around this depth...*/
            j = 0;
            while ((++j < n) &&( dtmp[j] < z)  ) 
                ;
	    if (j == n)
	        --j;
		
            if ((dtmp[j-1] == z) || (dtmp[j] == z) )
                datagap = 0;
            else if (z < 1001)
                datagap = (dtmp[j] - dtmp[j-1]) > GAP_SHALLOW;
            else
                datagap = (dtmp[j] - dtmp[j-1]) > GAP_DEEP;
                
         
  
            if (!datagap) {
                  ysum[i] += z;
                  ++count[i];
            }
      }
   }

   free((void *)xtmp);
   free((void *)ytmp);
   free((void *)dtmp);
   return;
} /* end insert_data() */

/*****************************************************************************/
void insert_deepest_vals(struct gridnode *g_ptr)
   /* Maintains a record containing the deepest observed level. 
       
         g_ptr:   ptr to current gridnode  
   */
{
   int i, j, reflev;
   int depth_greater, depth_diff;
   double key;
   struct deepestrec  *dr_ptr, *r1, *r2;
 
   dr_ptr = g_ptr->deepest; 
   j = hdr.nobs - 1;
 
 /*  empty list case ...*/ 
     
   if (dr_ptr == NULL) {
      dr_ptr = create_deepestrec();
      
      dr_ptr->depth = sta.observ[(int)DE][j];
      dr_ptr->pressure = sta.observ[(int)PR][j];
      dr_ptr->theta = sta.observ[(int)TH][j];
      dr_ptr->next = NULL;
      g_ptr->deepest = dr_ptr;
      return;
   }
   
   
/* Compare to existing deepest record.... */   
   
   dr_ptr = g_ptr->deepest;
   depth_greater = sta.observ[(int)DE][j] > dr_ptr->depth ? 1 : 0;

   
   /* if it is  deeper than the deepest observation so far, 
       replace the deepest depth/pressure/theta
       into the linked list */
   
   if (depth_greater) {
     r1 = g_ptr->deepest;
     r1->depth = sta.observ[(int)DE][j];
     r1->pressure = sta.observ[(int)PR][j];
     r1->theta = sta.observ[(int)TH][j];
   }
   
   return;
   
} /* end insert_deepest_vals () */
/*****************************************************************************/
struct deepestrec * create_deepestrec()
{
   struct deepestrec *r1;
   
   r1 = (struct deepestrec *) malloc(sizeof(struct deepestrec));
   if (r1 == NULL) {
       fprintf(stderr,"\nOut of memory in create_deepestrec() \n");
       exit(1);
   }
  
   return(r1);
   
} /* end create_deepestrec() */
/*****************************************************************************/
/*****************************************************************************/
void insert_surf_vals(struct gridnode *g)
   /* defines the depth of a mixed layer at the sea surface, determines 
      the depth-averaged value of theta within it, and inserts the information
      into a linked list of records sorted by sst. The station data are 
      accessed through the global variable: struct HYDRO_DATA sta 
              g:   ptr to gridnode  
   */
{
   int i, j, n, weight;
   double x, v, vprev, dprev;
   double sst, depth;
   struct surfrec *m;

   if (de[0] > 10.0)
       return;

   /* round off surface temperature to nearest tenth */

   sst =  (double) (NINT(sta.observ[(int)TH][0] * 10.)) / 10.0;

  /* bottom of mixed layer is defined where theta = SST + mix_layer_def */  
     
   depth = interpolate((sta.observ[(int)TH][0]+ mix_layer_def), sta.observ[(int)TH], de, hdr.nobs);

   if (depth < -999.)  /* all observations are within defined mixed layer range of temperature */
       depth = de[hdr.nobs-1];

   m = get_surfrec(g->mix_layer, sst);

   if (g->mix_layer == NULL ) {   /* empty list */           
      g->mix_layer = m;
   }

   /* add depth info to this record */

   m->depth += depth;
   ++m->n;

  /* compute average theta for mixed layer 
     The average value is computed by summing the observed values weighted by
     the depth between observations.*/

         j = 0;
         x = sta.observ[(int)TH][0];
         n = 1;               /* weight the observation at sea surface */
         if (x < -8.9 ) {     /* unless it is a missing value */
            x = 0.0;
            n = 0;
         }
         dprev = 0.0;
         vprev = sta.observ[(int)TH][0];
         while ( (j < hdr.nobs) && (de[j] <= depth)) {
            if ( (v = sta.observ[(int)TH][j]) > -8.9) {
               if (vprev < -8.9) 
                   vprev = v;
               weight = (int) (de[j] - dprev);
               x += (v + vprev) * .5 * weight;
               n += weight;
               dprev = de[j];
               vprev = v;
            }
            ++j;
         } 
         if (n > 0) {
            m->tavg += (x / (float) n);
         }
      }
   }
   return;

}  /* end insert_surf_vals() */

/*****************************************************************************/
struct surfrec *get_surfrec(struct surfrec *rptr, double key)

   /* Recursively searches a linked list of surfrecs to:
        a) find an existing record for the specified surface temperature (sst);
     or b) create a new record and insert it into the ordered list.

      Returns a pointer to the appropriate record.
      
  arguments: 
       rptr:    pointer to start of list 
       key:    temperature to key on 
   */
{
    double key;
    struct surfrec *r1ptr;

    if (rptr == NULL) {         /* empty list */
       r1ptr = create_surfrec();
       r1ptr->sst = key;
       return(r1ptr);
    }
    if (NINT(key * 10) == NINT(rptr->sst * 10)) {  /* current rec */
        return (rptr);
    }

    if (key < (rptr->sst - .00001)) {   /* comes before the current rec */

       r1ptr = create_surfrec();

         /* copy all fields from rptr into r1ptr */
       r1ptr->sst = rptr->sst;
       r1ptr->depth = rptr->depth;
       r1ptr->tavg = rptr->tavg;
       r1ptr->n = rptr->n;
       r1ptr->next = rptr->next;

        /* zero the fields of rptr and link it to r1ptr */
       rptr->sst = key;
       rptr->depth = 0;
       rptr->tavg = 0;
       rptr->n = 0;
       rptr->next = r1ptr;
       
       return(rptr);
    }

    r1ptr = get_surfrec(rptr->next, key);  /* search rest of list */
    if (rptr->next == NULL)
          rptr->next = r1ptr;
    return (r1ptr);

}   /* end get_surfrec() */

/*****************************************************************************/
struct surfrec *create_surfrec()

   /* Allocates memory and initializes the fields of a struct surfrec.
      Returns a pointer to the record */
{
   struct surfrec *r;
   int i;

   r = (struct surfrec *) malloc(sizeof(struct surfrec));
   if (r == NULL) {
      fprintf(stderr,"\nUnable to allocate memory in create_surfrec()\n");
      exit(1);
   }
   r->depth = 0;
   r->sst = 0;
   r->tavg = 0;
   r->n = 0;
   r->next = NULL;
   r->n = (UI *) malloc(nprops * sizeof(UI));


   return (r);
}   /* end create_surfrec() */
/*****************************************************************************/
void compute_avg(double *sum, UI *nobs, int nlevs)

   /*   sum:    array containing summed values 
       nobs:    array containing count of values in each sum 
      nlevs:    # of elements in array 
   */
{
   int j;

   for (j = 0; j < nlevs; ++j) {
      sum[j] = (nobs[j] > 0) ?  sum[j] / nobs[j] : (double) HBEMPTY;
   }
   return;

}  /* end compute_avg() */
/*****************************************************************************/
/*****************************************************************************/
struct surfrec *define_avg_mixed_layer(struct surfrec *listptr)

   /* Traverses a linked list of surfrecs sorted by sst and finds 
      the average sst , then searches the linked list
      for that sst or 2  records bracketing the average sst and
      returns a pointer to a surfrec containing the property values
      associated with the average sst. 
      
      listptr:   address of first record in linked list
       nprops:   number of properties being averaged 
   */ 
{
   struct surfrec *r1;
   double x;
   UI n;
   int key;
   
   if (listptr == NULL)                         /* empty list */
       return ((struct surfrec *) NULL);

   r1 = listptr;
   x = 0.0;
   n = 0;
   while (r1 != NULL) {
      x += r1->sst ; 
      ++n;   
      r1 = r1->next;    
   }
   key = NINT((x / (double)n) * 10.);     /* average sst * 10  */
   
   /* Now search list for the key value ... */

   r1 = get_sst_rec(key, listptr);

   return (r1);

}  /* end define_avg_mixed_layer() */

/***********************************************************************/
struct surfrec *get_sst_rec(int key, struct surfrec *rptr)
    /* Recursively searches a linked list sorted by sst for the record
       corresponding to key or for two records bracketing that value.
       Returns a record containing averaged property values for that density. 
       
          key:  sst  being searched for *10  
         rptr:  ptr to element of linked list 
 
    */
{
   struct surfrec *r1, *r2, *new;
   double x[2], y[2];
   UI n;
   int  key1, key2, i;

   if (rptr == NULL) {           /* YIKES! */
       fprintf(stderr,"\nError in get_sst_rec(): ");
       fprintf(stderr,"\nEnd of linked list encountered before finding avg sst!");
       fprintf(stderr,"\nError in program logic....exiting.\n");
       exit(1);
   }

   r1 = rptr;
   r2 = rptr->next;

   if (r2 == NULL) {             /* end of list, so r1 must be the sst */
      r1->tavg /= (double)r1->n;
      r1->depth /= (double)r1->n;
      return (r1);
   }

   key2 = NINT(r2->sst * 10);

   if (key > key2)                /* recursive part */
       return( get_sst_rec(key, r2));
   

   if (key == key2) {             /* exact match! */
       r2->tavg /= (double)r2->n;
       r2->depth /= (double)r2->n;
       return(r2);
   }
     
   key1 = NINT(r1->sst * 10);

   if (key < key1) {              /* YIKES! */
       fprintf(stderr,"\nError in get_sst_rec(): ");
       fprintf(stderr,"\nAvg density is less than first sst in list!");
       fprintf(stderr,"\nThis is a bug in the program logic....exiting.\n");
       exit(1);
   }
   if (key == key1)  {            /* exact match! */
       r1->tavg /= (double)r1->n;
       r1->depth /= (double)r1->n;
       return(r1);
   }

/* if we get this far, then the key must lie between r1 & r2... so 
      choose r1 or r2-- whichever is closest.   */


   if ((key - key1) < (key2 - key) ) {
       r1->tavg /= (double)r1->n;
       r1->depth /= (double)r1->n;
       return(r1);
   }
   else {
       r2->tavg /= (double)r2->n;
       r2->depth /= (double)r2->n;
       return(r2);
   }


}  /* end get_sst_rec() */

/***********************************************************************/
int do_std_depth(struct gridnode *g, int npts, int ncols, int icol, float **dataptr, short **countptr, float *bottom)

/*  For this gridnode, interpolates potential temperature onto the std depth levels 
    and computes pressure at each level.  The arrays at dataptr and countptr
    are assumed to have dimensions data[2][nstdlevs*ncols].  This
    function deals with one gridpoint at a time, therefore icol must be
    specified to determine the starting position within the array for
    each property; i.e. starting address = data[iprop][nstdlevs * icol].  
    (The arrays have been setup to include an entire row of
    data to improve the efficiency of writing out to the netcdf file).
    THE CALLING PROGRAM IS RESPONSIBLE FOR ALLOCATING THE SPACE AT
     *dataptr AND *countptr: 
   
    The # of observations at each standard level and for each property 
    are approximated and returned starting at the address pointed to by 
    countptr.

    The function returns the # of std levels (including the bottom) . 
    std_depth[MAXSTDLEVS] is globally defined and already initialized.

            g:   address of this gridnode  
         npts:   size of arrays for this gridnode  
        ncols:   # of cols in data and count arrays  
         icol:   column # associated with this gridpt  
      dataptr:   ptr to start of data arrays   
     countptr:   ptr to start of count arrays  
       bottom:   returned depth of deepest observation 

*/
{
   int i, ii, j, jj, k, datagap, room, npropsout;
   int *nsort;
   double *t2,*d2,*n2, z, r;
   double *tsort, *dsort;
   float  *d;
   short  *n;
   int size, start;
   int deepestlev;
   struct surfrec *m, *m2;
   struct deepestrec *btmptr, *b1;
   extern double std_depth[];
   
   
   
   /*  Strategy:  
         construct depth, theta, n arrays for all levels found (including
                  multiple crossings of tlevs)
         sort them by depth
	 merge mixed layer, theta-series, bottom level 
	 interpolate theta onto standard depth levels
	 return them ready to write to cdf file
   */
		  
		  
   npropsout = 2;   /* pressure and theta are output properties */		 

   tsort = (double *) malloc(MAX_CROSSINGS*npts * sizeof(double));
   dsort = (double *) malloc(MAX_CROSSINGS*npts * sizeof(double));
   nsort = (double *) malloc(MAX_CROSSINGS*npts * sizeof(double));
   size = sort_by_depth(g->d, g->nobs, dsort, tsort, nsort);

   if (size <= 1) {
      for (i = 0; i < npropsout; ++i) {
           d = &dataptr[i][icol * NSTDLEVS];
           n = &countptr[i][icol * NSTDLEVS];
           for (j = 0; j < NSTDLEVS; ++j) {
               *(d++) = (float) HBEMPTY;
               *(n++) = 0;
           }
      }
      *bottom = (float) HBEMPTY;
      free((void *)tsort);
      free((void *)dsort);
      free((void *)nsort);
      return(NSTDLEVS);
   }

/* construct temporary arrays, t2,d2,n2 to store depth, theta and 
   # of obs .  Enlarge the size of the arrays 
   to accommodate the mixed layer and deepest observations.  */

   room = 50;
   t2 = (double *) malloc((size+room) * sizeof(double));
   d2 = (double *) malloc((size+room) * sizeof(double));
   n2 = (double *) malloc((size+room) * sizeof(double));
   if (n2 == NULL) {
     fprintf(stderr,"\nUnable to allocate memory in do_std_depth()\n");
     exit(1);
   }

   m = define_avg_mixed_layer(g->mix_layer);

/* find index of deepest level in sorted arrays for which there is data */

   deepestlev = 0;
   for (k = 0; k < size; ++k) {
      if (nsort[k] > 0)
           deepestlev = k;
   }
   
/* Define the starting point in the sorted series which is deeper than
    the bottom of the average mixed layer. Check the temperature gradient for correct sign */
    
   start = 0;
   if (m != NULL) {
   
      if (m->depth > dsort[deepestlev])
          m->depth = dsort[deepestlev];


      while ((m->depth > dsort[start]) && (start < size))
         ++start;
   }



/* Eliminate missing values from each y-property. Interpolate the property
   onto std depth surfaces.  Check for vertical datagaps and flag a standard
   depth as missing if the gap is too large.  Also interpolate to approximate 
   the # of observations at each std level.  The interpolation routine returns 
   HBEMPTY when appropriate, so that each level is assigned a value
   whether or not any data is present.  The value of count at a level with 
   no data is assigned zero. */



   for (j = 0; j < nprops; ++j) {
         npts = 0;
         if (m != NULL) {
            if (m->count[j] > 0) {   
              npts = mixed_layer_vals(j, prop_id[j], m, x, y, w, i_te, i_sa);
              hb_svan(m->prop[i_sa], m->prop[i_te], 0.0, &sig[0]);
              for (i = 1; i < npts; ++i)
                  sig[i] = sig[0]; 
              sig[npts-1] = sig_bml;
            }  
         }

         if ((npts - start ) > room)  {
            fprintf(stderr, "\nIncrease the value of 'room' to at least %d", npts-start);
            fprintf(stderr, "\n in do_std_depth()\n");
            exit(1);
         }

         for (k = start; k < size; ++k) {
            if (g->count[j][k] > 0) {
               y[npts] = g->prop[j][k];
               w[npts] = (double) g->count[j][k];
               x[npts] = g->d[k];
               sig[npts] = sigval[k];
               ++npts;
            }  
         } 
         
         /* add deepest observation to arrays.  If there was no
         measurement at the deepest level, check the linked list
         for a measurement */


         btmptr = g->deepest;
        *bottom = (float) btmptr->depth;
	
	/* check that sigma of deepest observation is denser than deepest sig value */
	
	reflev = 0;
	if (g->d[deepestlev] >= zmin[4])
	    reflev = 4;
	else if (g->d[deepestlev] >= zmin[3])
	    reflev = 3;
	else if (g->d[deepestlev] >= zmin[2])
	    reflev = 2;
	else if (g->d[deepestlev] >= zmin[1])
	    reflev = 1;
	
	switch (reflev) {
	    case 0:
	       is_denser = sigval[deepestlev] <= btmptr->sig_0 ? 1 : 0;
	       sig[npts] = btmptr->sig_0;
	       pref = 0.0;
	       break;
	    case 1:
	       is_denser = sigval[deepestlev] <= btmptr->sig_1 ? 1 : 0;
	       pref = 1000.0;
	       sig[npts] = btmptr->sig_1;
	       break;
	    case 2:
	       is_denser = sigval[deepestlev] <= btmptr->sig_2 ? 1 : 0;
	       pref = 2000.0;
	       sig[npts] = btmptr->sig_2;
	       break;
	    case 3:
	       is_denser = sigval[deepestlev] <= btmptr->sig_3 ? 1 : 0;
	       pref = 3000.0;
	       sig[npts] = btmptr->sig_3;
	       break;
	    case 4:
	       is_denser = sigval[deepestlev] <= btmptr->sig_4 ? 1 : 0;
	       pref = 4000.0;
	       sig[npts] = btmptr->sig_4;
	       break;
	} /*end switch */
   
        if (is_denser) {    /* add a property value and increment npts */
           y[npts] = btmptr->prop[j];
           
           if (prop_id[j] == (int) PR)       /* explicitly insert pressure if it is the y prop */
                y[npts] = btmptr->pressure;
                
                
           b1 = btmptr;
           while ((b1->next != NULL) && (y[npts] < -8.)) {
              b1 = b1->next;
              y[npts] = b1->prop[j];
           }
           if (y[npts] > -8.) {  /* found an actual obs */
             w[npts] = 1.0;
             x[npts] = btmptr->depth;
             btmptr->prop[j]  = y[npts];  
                 ++npts;
           }  
         }
	 
         d = &dataptr[j][icol * NSTDLEVS];
         n = &countptr[j][icol * NSTDLEVS];
         
        for (i = 0; i < NSTDLEVS-1; ++i) {   /* interpolate all but bottom */
           if (npts <= 1) {
               z = (double) HBEMPTY;
               k = 0;
           }
           else {
             z =  interpolate(std_depth[i], x, y, npts);
             k = 0;
           }

           if (z > -998.) {                  /* check for vertical datagaps */
              jj = 0;
              while (x[++jj] < std_depth[i]) 
                   ;

              if (((int)x[jj-1] == (int)std_depth[i]) 
                  || ((int)x[jj] == (int)std_depth[i])) 
                  datagap = 0;

              else if (std_depth[i] < 1001) 
                  datagap = (x[jj] - x[jj-1]) > GAP_SHALLOW;

              else
                  datagap = (x[jj] - x[jj-1]) > GAP_DEEP;

              if (datagap) {
                  z = (double) HBEMPTY;
                  if (j == 0 && report_gaps)
                    fprintf(stderr," datagap:  %.1lf  %.1lf\n", x[jj-1], x[jj]);
              }
              else {
                   r = interpolate(std_depth[i], x, w, npts);
                   if (r < 0)
                       k = 0;
                   else {
                       k = (r - (int)r) > 0 ? (int) r + 1 : (int) r;
                   }
              }
           }

           *(d++) = (float) z;
           *(n++) = (short) k;
         }

/* add the deepest observation */

       *(d++) = (float) btmptr->prop[j];
       *(n++) = (btmptr->prop[j] > -8.) ? 1 : 0;
   }

/* clean up space no longer needed ... */

   free((void *)sigval);
   free((void *)sig);
   free ((void *)x);
   free ((void *)y);
   free ((void *)w);
   m = g->mix_layer;
   while (m != NULL) {
      m2 = m->next;
      free((void *) m->prop);
      free((void *) m->count);
      free((void *)m);
      m = m2;
   }
   g->mix_layer = NULL;

   return(NSTDLEVS);

} /* end do_std_depth() */
/****************************************************************************/

int mixed_layer_vals(int iprop, int prop_id, struct surfrec *m, double *x, double *y, double *w, int i_te, int i_sa)

   /* Creates arrays of depth, property, and counts at
      100 m intervals for the requested property.  Returns
      the number of points in each array. 
      
         iprop:          indicates the ith property of nprops 
       prop_id:       indicates the ith property of MAXPROPS 
             m:      info defining the mixed layer 
             x:       depth array 
             y:       property array 
             w:       count array 
    i_te, i_sa:       index to temperature, salt (or -1) 
 
      */
      
{
   int i, npts;
   double  *t, *s, pref;
   float deltap;

/* Define top and bottom only unless the mixed layer depth
   exceeds 200 m ... */

   npts = 2;       
   if (m->depth > 199.)
        npts = 1 + (int) m->depth / 100;

/*  Assign depth and count values ... */

   for (i = 0; i < npts-1; ++i) {
         w[i] = m->count[iprop];
         x[i] = i * 100.0;           
   }         
   w[npts-1] = m->count[iprop];
   x[npts-1] = m->depth;


   return (npts);

}  /* end mixed_layer_vals() */
/****************************************************************************/
int sort_by_depth(double *d, UI  *nobs, double *thetas, int npts)

/* sorts the arrays into increasing order by depth, 
   ensuring there are no depth inversions.  Returns the 
   # of depth levels containing data. 
   
           d    starting addr of depth  
        nobs    starting addr of nobs array  
      thetas    array of theta values 
        npts    # of cols (points) in each property array 

   */
     
{
   int i, j, k, mono;
   int size;
   double  *diff;

   diff = (double *) calloc(npts,sizeof(double));
   if (diff == NULL) {
     fprintf(stderr, "\nUnable to allocate memory in monotonic()\n");
     exit(1);
   }

/* first, eliminate any depth levels with no observations */
   k = 0;
   diff[0] = 0.0;
   mono = 1;
   for (i = 0; i < npts; ++i) {
      if (nobs[i] > 0) {
         d[k] = d[i];
         nobs[k] = nobs[i];
         sigmas[k] = siglevs[i];
         for (j = 0; j < nprops; ++j) {
           xx[j][k] = xx[j][i];
           count[j][k] = count[j][i];
         }
	 if (k > 0) {
	    diff[k] = d[k] - d[k-1];
	    mono = mono && (diff[k] >= 0);
	 }
         ++k;
      }
   }
   
   
   size = k;
   if ( mono || (size < 2)) {
       free((void *) diff);
       return(size);
   }

   
   /* Remove levels where depth inversions occur */

   while (!mono) {
     k = 0;
     for (i = 0; i < size; ++i) {
       if (diff[i] > -1 ) {
          d[k] = d[i];
          nobs[k] = nobs[i];
          sigmas[k] = siglevs[i];
          for (j = 0; j < nprops; ++j) {
             xx[j][k] = xx[j][i];
             count[j][k] = count[j][i];
          }
	  ++k;
	}
	else if (i > 2) {
	     if ((d[i] > d[i-2]) && k > 0) {
	        --k;
                 d[k] = d[i];
                nobs[k] = nobs[i];
                sigmas[k] = siglevs[i];
                for (j = 0; j < nprops; ++j) {
                   xx[j][k] = xx[j][i];
                   count[j][k] = count[j][i];
                }
	        ++k;
	     }
	 
	     else {
	     
	       if (k >= 1) {
	          --k;
	          if ( i < (size-1)) {
	            diff[i+1] = d[i+1] - d[k-1];
	          }
	       }
	     }  
	 }
        
     } /* end for */
     
     size = k;
     mono = 1;
     diff[0] = 0.0;
     if (size > 2) {
        for (i = 1; i < size; ++i) {
	 diff[i] = d[i] - d[i-1];
	 mono = mono && (diff[i] >= -1);
        }
      }
   } /* end while !mono */
   
   free((void *) diff);

   return (size);
   

} /* end sort_by_depth() */ 

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/* Performs a linear interpolation to find the position of xval in array, x,
   and returns the corresponding value in the array, y.  If xval does not
   appear in array x, the value HBEMPTY is returned.  This routine assumes 
   that the x array is monotonic and continuous (no missing values); and it
   assumes the y array is continuous.    */

double interpolate(double xval, double *x, double *y, int nypts)
{
   int    k;
   double v1, v2;

   for (k = 0; k < nypts-1; ++k) {

      v1 = xval - x[k];
      v2 = xval - x[k+1];

      if (v1 == 0)             /* x[k] == xval */
          return (y[k]);
      if (v2 == 0)             /* x[k+1] == xval */
          return (y[k+1]);
      if (v1 < 0. && v2 < 0.)  /* xval not between x1 and x2 */  
          continue;
      if (v1 > 0. && v2 > 0.) 
          continue;

      return ( y[k] + (y[k+1] - y[k]) * v1 / (x[k+1] - x[k]) );
   }

   return ((double)HBEMPTY);

}   /* end interpolate() */

/****************************************************************************/

void sum_levels(double *y, double *d, int ny, double *ysum, UI *count)
  /* Interpolates the y and d arrays with ny observation levels onto the stddepth levels
  (globally defined).  Checks for vertical data gaps and does not interpolate over
   gaps which exceed 200 m in the thermocline (upper 1000 m) or 1000 m 
   elsewhere. The interpolated values are added to the appropriate
  levels of the array defined by ysum and the count array is incremented. */ 
{
   int  i, j, nz, datagap, n;
   double *ytmp, *dtmp, yint, flag;

   ytmp = (double *) calloc((size_t)ny, sizeof(double));
   dtmp = (double *) calloc((size_t)ny, sizeof(double));
   if (dtmp == NULL) {
       fprintf(stderr,"\nUnable to allocate memory in sum_levels()\n");
       exit(1);
   }

   nz = NSTDLEVS - 1;   /* NSTDLEVS is globally defined and includes a bottom depth*/
   flag = HBEMPTY + 1.0;  
   
/* Ensure continuous (no missing values) for  y and d arrays */
   n = 0;
   for (i = 0; i < ny; ++i) {
      if ( y[i] > -8.9)  {
         ytmp[n] = y[i];
         dtmp[n++] = d[i];
      }
   }
   
   if (n <= 1) {               /* not enough values */
     free((void *)ytmp);
     free((void *)dtmp);
     return;
   }


   for (i = 0; i < nz; ++i) { 
      if ((yint = interpolate(std_depth[i], dtmp, ytmp, n)) > flag) {

            /* now check for datagaps around this depth...*/
	    
            j = 0;
            while (dtmp[++j] < std_depth[i]) 
                ;
		
            if ((dtmp[j-1] == std_depth[i]) || (dtmp[j] == std_depth[i]) )
                datagap = 0;
            else if ( std_depth[i] < 1001)
                datagap = (dtmp[j] - dtmp[j-1]) > GAP_SHALLOW;
            else
                datagap = (dtmp[j] - dtmp[j-1]) > GAP_DEEP;
                           
  
            if (!datagap) {
                  ysum[i] += yint;
                  ++count[i];
            }
       } /* end if */
   } /* end for i */
   
   free((void *)ytmp);
   free((void *)dtmp);
   return;

} /* end sum_levels() */
