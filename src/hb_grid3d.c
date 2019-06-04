/*  hb_grid3d.c
...............................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             original: 1993
                             Updated to ANSI standards Feb 2000
...............................................................................
..........................................................................
. Computes an average profile for each property specified at each lat/lon
. gridnode from randomly spaced stations and outputs a netCDF file of hydro 
. properties gridded as a function of lat, lon, and depth.  The depth
. dimension consists of a series of standard depths, but these values also
. can be optionally specified by the user with -Z.  The averaging is done 
. on density surfaces and interpolated back onto the depth levels.  The user
. supplies the sigma values with the -S option.  Optimum results will be
. obtained by customizing the sigma series for a particular part of the
. ocean (see hb_sigma_eval).  The goal is to specify sigma values which 
. span the entire water column with appropriate sampling as a function 
. of depth.  
..........................................................................
*/

/***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hydro_cdf.h"
#include "hydrobase.h"

/***************************************************************/

#define   RADperDEG 0.017453292  /* pi/180 */
#define   EarthRadius  6371.0    /* km */

#define    UI   unsigned int
#define    PREFILL      0        /* turn off prefilling of cdf variables */
#define    ADD_PROP_COUNTS  1    /* turn on option to include # of obs info */
#define    ADD_PROP_ERRORS  1    /* turn on option to include std error */
#define    PRINT_MSG  1          /* turn on message printing to stderr */
#define    MAX_TIME_BINS  1      /* time dimension stores just one matrix  */
#define    MIX_LAYER_DEF 0.02    /* default sigma range for mixed layer */

#define    EXTENT   ""           /* Input_file extension and pathnames */
#define    DIR      ""

/***************************************************************/

/*  Define the standard sigma levels on which data will be averaged;
    sigma series is divided by reference levels. */

#define  MAXSIGLEVS 800    /* max size of sigma series -- arbitrary number */
#define  NREFLEVS  5       /* # of ref levs defining sigma values
                            * corresponding to sigma 0,1,2,3, and 4. */

/***************************************************************/

double siglevs[MAXSIGLEVS];   /* stores the sigma series  */
int  nsiglevs;                /* number of elements in full sigma series */

int ref_prop_id[NREFLEVS] = {(int)S0,(int)S1,(int)S2,(int)S3,(int)S4};
int nsig[NREFLEVS];       /* number elements of sigseries for each ref lev */
int zmin[NREFLEVS] = {-1,500,1500,2500,3500}; /* depth ranges for ref levs */
int zmax[NREFLEVS] = {500,1500,2500,3500,99999};
double *sigptr[NREFLEVS];   /* ptr to station arrays of ref sigma values */

/***************************************************************/

/* Globally referenced variables for input station data */

struct HYDRO_DATA sta;
struct HYDRO_HDR hdr;
double *pr, *de, *te, *sa;
double s_pref;
double pe_pref;
double ht_pref;

int report_gaps;           /* switch turns on reporting of vertical datagaps */
double mix_layer_def;      /* option to specify mixed layer definition */
float lengthscale;         /* L for weight function = e^[-(d/L)^2] */

/* Global variables used in computing potential vorticity. */
double dlat;
float latitude, longitude;
int window, w_incr;     /* used to specify pr window for vertical
                         * gradient properties*/

/************ structures to represent grid node ****************/
struct surfrec {
  double  density, depth;
  double *prop;
  double *wghtsum;
  UI  n;
  struct surfrec *next;
};

struct deepestrec {
  double  depth, pressure;
  double sig_0, sig_1, sig_2, sig_3, sig_4;
  double *prop;
  struct deepestrec *next;
};

struct gridnode {
  double **prop;        /* Property values at node - prop and depth dims. */
  double **weightsum;   /* Weightsum, one for each property value. */
  double  *d;           /* Depth of node (depth dim only). */
  double  *dweightsum;  /* Weightsum specific to depth (depth dim only). */
  UI **count;           /* Count, prop and depth dims */
  UI  *nobs;            /* observation count (depth dim only) */
  struct deepestrec *deepest; /* Prop dim only */
  struct surfrec  *mix_layer; /* Prop dim only */
};
/***************************************************************/

/*  prototypes for functions defined within this module... */

struct gridnode **alloc_grid(int, int, int, int);

int get_sig_series(int, FILE *,  double *);
int get_time_bins(int, int, struct CDF_HDR *);
int parse_p_option(char *, int *);
int do_std_depth(struct gridnode *, int, int *, int, int, int, float **, short **, float *);
int mixed_layer_vals(int, int, struct surfrec *, double *, double *, double *, int, int);


/*  Remove these functions:
int sort_by_depth(double *, UI *, double *, double **, UI **, int, int);
int find_pivot(double *, int, double *);
int partition(double *, int *, int, double);
*/

void free_grid(int, int, int, struct gridnode **);
void parse_s_option(char *, FILE **);
void insert_data(double *,int, double *, double *, double, UI *, int);
void insert_surf_vals(struct gridnode *,int *,int);
void insert_deepest_vals(struct gridnode *,int *,int);
void compute_avg(double *,double *, int);
void sum_levels(double *, double *, int, double *, UI *);
void print_usage(char *);
void define_sigma_transitions(struct gridnode *, int);
void quicksort(double *, int *, int);
void swap_d(double *, double *);
void swap_i(int *, int *);
void delete_list(struct deepestrec *);
void delete_surfrec_list(struct surfrec *);

struct surfrec *get_surfrec(struct surfrec *, double, int);
struct surfrec *create_surfrec(int);
struct surfrec *get_density_rec(int, struct surfrec *, int);
struct surfrec *define_avg_mixed_layer(struct surfrec *, int);
struct deepestrec *create_deepestrec(int);

double interpolate(double, double *, double *, int);
double get_weight(float, float, float, float, float); 


/***************************************************************/
int main (int argc, char **argv) {

  /* Flags marking presence of options on command
   * line (and run mode). */
  short   bopt, iopt, popt, copt, no_distance_weight;

  /* File counters. */
  int     curfile = 1, nfiles = 0;

  int     print_mess, xdateline, ratio_done;

  /* nlevs = # depth levels
   * nprops = number of requested properties.
   * prop_indx = property ID indices. */
  int     nlevs, nprops = 0, prop_indx[MAXPROP];

  int     i, j, error, test_outcrop, prop_avail;
  int     nstdlevs, tmin, tmax, n_filled;
  int     include_counts;
  int     include_errors;
  FILE   *sigfile[NREFLEVS]; /* List of files containing density surfaces. */
  FILE   *z_file;            /* File containing list of depths to interp. */
  int     row, col, nrows, ncols, tbin, ntbins;
  int     in_bounds, out_of_bounds=0;
  char   *extent, *dir, *st;
  char   *cdf_filename;
  int     infile, cdf_file;
  float  **data, *bottomdepth;  
  float  latc, lonc;
  short  **count; 
  struct CDF_HDR  h;
  struct gridnode **grid[MAX_TIME_BINS];
  double weight;


  /* Are there command line arguments? */
  if (argc < 2) {
    print_usage(argv[0]);
    exit(1);
  }
 
  /* Set default values. */
  dir = DIR;
  extent = EXTENT;
  bopt = iopt = popt = copt = 0;
  z_file = NULL;
  for (j = 0; j < NREFLEVS; ++j) {
    sigfile[j] = NULL;
  }
  include_counts = ADD_PROP_COUNTS;
  include_errors = ADD_PROP_ERRORS;
  error = 0;
  print_mess = PRINT_MSG;
  report_gaps = 0;
  xdateline = 0;
  window = 100;
  w_incr = 10;
  mix_layer_def = MIX_LAYER_DEF;
  tmin = 0;
  tmax = 9999;
  no_distance_weight = 0;
  lengthscale = 100.0;   /* L in km for distance weighting = e^-(d/L)^2 */
    
   
  /* Parse the command line arguments */
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
	       /* Get grid bounds and copy values
		* into the cdf output file header, h. */
               case 'B':                    
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

               case 'L':
		 error = (sscanf(&argv[i][2],"%f", &lengthscale) == 1) ? 0 : 1;
		 break;
                        
               case 'M':
		 error = (sscanf(&argv[i][2],"%lf", &mix_layer_def) == 1) ? 0 : 1;
		 break;
                        
               case 'P':
		 popt = 1;
		 nprops = parse_p_option(&argv[i][2], prop_indx);
		 break;

               case 'S':
		 parse_s_option(&argv[i][2], sigfile);
		 break;

               case 'T':
		 st = &argv[i][2];
		 if (*st == '/') 
		   ++st;
		 error = (sscanf(st,"%d/%d", &tmin, &tmax) == 2) ? 0 : 1;
		 break;

               case 'W':
		 error = (sscanf(&argv[i][2],"%d", &window) == 1) ? 0 : 1;
		 st = &argv[i][2];
		 while (*(st++) != '\0') {
		   if (*st == '/') {
		     ++st;
		     error = (sscanf(st,"%d", &w_incr) == 1) ? 0 : 1;
		     break;
		   }
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

          } /* End command line flag/option parsing switch. */

          if (error) {
             fprintf(stderr,"\nError parsing command line args.\n");
             fprintf(stderr,"     in particular: '%s'\n", argv[i]);
             print_usage(argv[0]);
             exit(1);
          }

       }  /* End command line flag/option check. */

       else  {
           ++nfiles;
       }
   }  /* End for loop over all command line entries. */

  if (!bopt || !iopt || !nfiles || !popt || !copt) {
    fprintf(stderr,"\nYou must specify input file(s), bounds, props, output file, and gridspacing!\n");
    print_usage(argv[0]);
    exit(1);
   }

  if (lengthscale == 0)
    no_distance_weight = 1;

  /* Initialize global array of std_depths */
  nstdlevs = std_depth_init(z_file);

  /* Define values for sigma series ... */
  nsiglevs = 0;
  for (i = 0; i < NREFLEVS; ++i ) {
    nsig[i] = get_sig_series(ref_prop_id[i], sigfile[i], &siglevs[nsiglevs]);
    nsiglevs += nsig[i];
    if (nsiglevs >= MAXSIGLEVS) {
      fprintf(stderr,"\n # of sigma levels exceeds space allocated");
      fprintf(stderr,"\n Recompile program with larger MAXSIGLEVS.\n");
      exit(1);
    }
  }
  fprintf(stderr,"\n Total number of sigma levels: %d\n", nsiglevs);

  /* Compute dimensions of grid and allocate space for computation. */
   nrows = NINT(ceil((double)((h.ymax - h.ymin) / h.yincr)));
   ncols = NINT(ceil((double)((h.xmax - h.xmin) / h.xincr)));
   h.node_offset = 1;       /* 0 for node grid, 1 for pixel grid */
  
   ntbins = get_time_bins(tmin, tmax, &h);

   /* Set cdf header values */
   h.nx = ncols;
   h.ny = nrows;
   h.nz = nstdlevs;
   h.nt = ntbins;

   fprintf(stderr,"\n allocating memory for grid ... " );
   
   for (tbin = 0; tbin < ntbins; ++tbin) {
     grid[tbin] = alloc_grid(ncols, nrows, nsiglevs, nprops);
   }
   fprintf(stderr,"\n");

   if (no_distance_weight)
      fprintf(stderr,"Distance weighting will not be applied\n");
   else 
      fprintf(stderr,"Lengthscale for distance weighting [= e^-(d/L)^2]  is %.1lf km\n", lengthscale);
   
   /* Initialize structures to hold current station header + data. */
   for (i = 0; i < MAXPROP; ++i)
      sta.observ[i] = (double *) NULL;
   hdr.prop_id = (int *) NULL;

   /* Loop for each input_file. */
   do {
     if ((infile = open_hydro_file(dir, argv[curfile], extent, print_mess)) < 0) {
       goto NEXTFILE;
     }

     /* Loop for each station. */

     while ((error = get_station(infile, &hdr, &sta)) == 0) {

       /* Is the station within bounds ? */
       if (xdateline && hdr.lon < 0)
          hdr.lon += 360;
	  
       in_bounds = 1;
       if (get_indices(&h, hdr.lat, hdr.lon, &row, &col) < 0) {
	 /* If it is on the border, accept it */
	 if (row == -1)  {
	   row = 0;
	   if (col == -1)
	     col = 0;
	   else if (col == h.nx)
	     --col;
	 }
	 else if (row == h.ny) {
	   --row;
	   if (col == -1)
	     col = 0;
	   else if (col == h.nx)
	     --col;
	 }
	 else if (col == -1) {
	   col = 0;
	   if (row == -1)
	     row = 0;
	   else if (row == h.ny)
	     --row ;
	 }
	 else if (col == h.nx) {
	   --col;
	   if (row == -1)
	     row = 0;
	   else if (row == h.ny)
	     --row ;
	 }
	 else {
	   fprintf(stderr,"Out of bounds: %.3f %.3f\n", hdr.lat, hdr.lon);
	   ++out_of_bounds;
	   in_bounds = 0;
	 }
       }
       
       if (in_bounds) {
         ratio_done = 0;
         pr = sta.observ[(int)PR];   /* set frequently used pointers */
         de = sta.observ[(int)DE];
         te = sta.observ[(int)TE];
         sa = sta.observ[(int)SA];
         if (pr == NULL || te == NULL || sa == NULL || de == NULL) {
           fprintf(stderr, "\nRequired parameters: pr, te, sa, de are not available at this station.\n");
           write_hydro_hdr(STDERR, &hdr);
           continue;
         }

         /* Allocate space for and compute depth
	  * referenced sigmas for this station. */
         free_and_alloc(&sta.observ[(int)S0], hdr.nobs);
         compute_sigma(0., hdr.nobs, sta.observ[(int)S0], pr, te, sa);
         free_and_alloc(&sta.observ[(int)S1], hdr.nobs);
         compute_sigma(1000., hdr.nobs, sta.observ[(int)S1], pr, te, sa);
         free_and_alloc(&sta.observ[(int)S2], hdr.nobs);
         compute_sigma(2000., hdr.nobs, sta.observ[(int)S2], pr, te, sa);
         free_and_alloc(&sta.observ[(int)S3], hdr.nobs);
         compute_sigma(3000., hdr.nobs, sta.observ[(int)S3], pr, te, sa);
         free_and_alloc(&sta.observ[(int)S4], hdr.nobs);
         compute_sigma(4000., hdr.nobs, sta.observ[(int)S4], pr, te, sa);

         for (i = 0; i < NREFLEVS; ++i) {
            sigptr[i] = sta.observ[ref_prop_id[i]];
         }

         /* Compute weight as a function of distance
	  * from center of grid node. */
	 if (no_distance_weight ) 
	   weight = 1.0;
	 else {
	   get_lat_lon(&h, row, col, &latc, &lonc);
	   weight = get_weight(hdr.lat, hdr.lon, latc, lonc, lengthscale);
	 }

         /* Compute other appropriate properties for this station ... */
         for (i = 0; i < nprops; ++i) {
	   test_outcrop = 0;
	   prop_avail = 1;

	   /* !**! Special cases for individual properties.
	    * Additional computations are only performed if
	    * we have requested a derived property - i.e.
	    * temp and salt and pressure and depth are already
	    * in the profile and will not be computed here. */
	   switch ((enum property) prop_indx[i]) {

               case PR:
		 /* If the depth is less than 51 meters,
		  * then we assume that the surface has
		  * outcropped by setting test_outcrop to 1. */
                  test_outcrop = (de[0] <= 51.) ? 1 : 0;
                  break;

               case OX:  
                  prop_avail = available((enum property) prop_indx[i], &hdr);
                  if (!prop_avail) {
                    if (available(O2, &hdr)) {
                      free_and_alloc(&sta.observ[(int)OX], hdr.nobs);
                      for (j=0; j < hdr.nobs; ++j) {
                        sta.observ[OX][j] = ox_kg2l(sta.observ[(int)O2][j], pr[j], te[j], sa[j]);
                      }
                      prop_avail = 1;
                    }
                  }
                  break;
               
               case O2:  
                  prop_avail = available((enum property) prop_indx[i], &hdr);
                  if (!prop_avail) {
                    if (available(OX, &hdr)) {
                      free_and_alloc(&sta.observ[(int)O2], hdr.nobs);
                      for (j=0; j < hdr.nobs; ++j) {
                        sta.observ[OX][j] = ox_l2kg(sta.observ[(int)OX][j], pr[j], te[j], sa[j]);
                      }
                      prop_avail = 1;
                    }
                  }
                  break;

               case S_:
                  free_and_alloc(&sta.observ[(int)S_], hdr.nobs);
                  compute_sigma(s_pref, hdr.nobs, sta.observ[(int)S_],pr,te,sa);
                  break;

               case SD:
                  free_and_alloc(&sta.observ[(int)SD], hdr.nobs);
                  compute_dref_sigma(de, hdr.nobs, sta.observ[(int)SD],pr,te,sa);
                  break;

               case TH:
                  free_and_alloc(&sta.observ[(int)TH], hdr.nobs);
                  compute_theta(hdr.nobs, sta.observ[(int)TH], pr, te, sa);
                  break;

       	       case TP:
		  free_and_alloc(&sta.observ[(int)TP], hdr.nobs);
		  compute_dtdp(sta.observ[(int)TP], pr, te, sa, hdr.nobs, window, w_incr);
		  break;

               case TZ:
		  free_and_alloc(&sta.observ[(int)TZ], hdr.nobs);
		  compute_dtdz(sta.observ[(int)TZ], pr, te, sa, hdr.nobs, window, w_incr, (double)hdr.lat);
		  break;

               case HC:
                  free_and_alloc(&sta.observ[(int)HC], hdr.nobs);
                  compute_heat_capacity(hdr.nobs, sta.observ[(int)HC], pr, te, sa);
                  break;

               case HT:
                  free_and_alloc(&sta.observ[(int)HT], hdr.nobs);
                  compute_height(hdr.nobs, pr, te, sa, ht_pref, sta.observ[(int)HT]);
                  break;

               case PE:
                  free_and_alloc(&sta.observ[(int)PE], hdr.nobs);
                  compute_energy(hdr.nobs, pr, te, sa, pe_pref, sta.observ[(int)PE]);
                  break;

               case DR:
	       case AL:
	       case BE:
	          if (!ratio_done) {
                     free_and_alloc(&sta.observ[(int)DR], hdr.nobs);
                     free_and_alloc(&sta.observ[(int)AL], hdr.nobs);
                     free_and_alloc(&sta.observ[(int)BE], hdr.nobs);
                     compute_ratio(hdr.nobs, sta.observ[(int)DR], pr, te, sa, sta.observ[(int)AL], sta.observ[(int)BE]);
		     ratio_done = 1;
		  }
                  break;
               case SV:
                  free_and_alloc(&sta.observ[(int)SV], hdr.nobs);
                  compute_sp_vol(hdr.nobs, sta.observ[(int)SV], pr, te, sa);
                  break;

               case VA:
                  free_and_alloc(&sta.observ[(int)VA], hdr.nobs);
                  compute_svan(hdr.nobs, sta.observ[(int)VA], pr, te, sa);
                  break;

               case VS:
                   free_and_alloc(&sta.observ[(int)VS], hdr.nobs);
                  compute_sound_vel( sta.observ[(int)VS], pr, te, sa, hdr.nobs);
                break;


               case BF:   
                  free_and_alloc(&sta.observ[(int)BF], hdr.nobs);
                    buoy_freq(sta.observ[(int)BF],pr,te,sa,hdr.nobs,window, w_incr);
                  for (j = 0; j < hdr.nobs; ++j) {
                    if (sta.observ[(int)BF][j] > -9998.)
                      sta.observ[(int)BF][j] *= 1.0e5;
                  }
                  break;

               case PV:
                  free_and_alloc(&sta.observ[(int)PV], hdr.nobs);
                  buoy_freq(sta.observ[(int)PV],pr,te,sa,hdr.nobs,window, w_incr);
                  for (j = 0; j < hdr.nobs; ++j) 
                    sta.observ[(int)PV][j] *= sta.observ[(int)PV][j];
                  po_vort(sta.observ[(int)PV],sta.observ[(int)PV], hdr.nobs, (double)hdr.lat);
                  break;
		  
               case RR:
                  free_and_alloc(&sta.observ[(int)RR], hdr.nobs);
		  dlat = (double) hdr.lat;
		  compute_approx_rossby_radius(sta.observ[i], hdr.nobs, hdr.pdr, sta.observ[(int)DE], sta.observ[(int)PR], sta.observ[(int)TE], sta.observ[(int)SA], dlat, window, w_incr);
                  break;


		  
               default:
		   /* Default checks if the property is available
		    * and if it is not, overrides the prop_avail=1
		    * setting at the top of this code block. This
		    * is done especially for other measured props
		    * such as tracers and nutrients because many of
		    * the derived quantities (i.e. pv or bf) can 
		    * always be computed. */
                   prop_avail = available((enum property) prop_indx[i], &hdr);
                   break;
		   
	   } /* End property selection switch */
	   
	   /* Sort station into appropriate time bin
	    * and insert the data into averages. This
	    * section applies to any property in the
	    * station, even depth. */
	   for (tbin = 0; tbin < ntbins; ++tbin) {
	     if (hdr.year >= h.tmin[tbin] && hdr.year <= h.tmax[tbin]) {
	       if (prop_avail) {
		 insert_data(sta.observ[prop_indx[i]], hdr.nobs, grid[tbin][row][col].prop[i], grid[tbin][row][col].weightsum[i], weight, grid[tbin][row][col].count[i], test_outcrop);
	       }
	     } /* End time bin sort check. */
	   } /* End for loop over all time bins. */ 
         } /* End for loop over all possible requested properties. */

	 /* Insert depth values and deal with the sea surface/bottom ... */
         test_outcrop = (de[0] <= 51.) ? 1 : 0;
         for (tbin = 0; tbin < ntbins; ++tbin) {
	   if (hdr.year >= h.tmin[tbin] && hdr.year <= h.tmax[tbin]) {

	     /* Water column data */
	     insert_data( de, hdr.nobs, grid[tbin][row][col].d, grid[tbin][row][col].dweightsum, weight,
			  grid[tbin][row][col].nobs, test_outcrop);

	     /* Surface data */
	     insert_surf_vals(&grid[tbin][row][col], prop_indx, nprops);

	     /* Bottom data */
	     insert_deepest_vals(&grid[tbin][row][col], prop_indx, nprops);

	   } /* End if check on time bin sort. */
         } /* End for loop over all time bins. */
       } /* End if check on station location in_bounds */
       
       /* Free up space by removing space allocated for
	* properties which are not being requested. */
       for (i = 0; i < MAXPROP; ++i) {
          if (sta.observ[i] != NULL) {
             free((void *) sta.observ[i]);
             sta.observ[i] = (double *) NULL;
          } /* End if check empty space in observations. */
       } /* End of loop over all possible properties.*/

     }  /* End while !eof - loop over all stations in file. */ 
     if (error > 0) {
         report_status(error, stderr);
         exit(1);
     }

   NEXTFILE:
     close(infile);

   } while (curfile++ < nfiles ); /* End of loop over all input files. */


   /******** End of first input phase *********/

   if (out_of_bounds) {
     fprintf(stderr,"\n  %d out-of-bounds stations were read and ignored.\n", out_of_bounds);
   }

   /* Construct the cdf header  */
   fprintf(stderr,"\n  Constructing header ...\n");
   h.nprops = nprops;
   h.fill_value = HBEMPTY;
   h.mask_value = HBMASK;
   strncpy(h.x_units, "degrees", 8);
   strncpy(h.y_units, "degrees", 8);
   strncpy(h.z_units, "meters", 7);
   strncpy(h.title,"HydroBase", 10);
   strcpy(h.command, *argv);
   for (i = 1; i < argc; ++i) {
     if (strlen(h.command) < 900) {
       strncat(h.command, " ", 1);
       strcat(h.command, argv[i]);
     }
   }
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

   /* Create arrays to write to cdf. */
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
   
   /* For each gridnode, compute means at all sigma levels. */
   fprintf(stderr,"  Computing averages ...\n");
   for (tbin = 0; tbin < ntbins; ++tbin) {
     for (row = 0; row < nrows; ++row) {
       for (col = 0; col < ncols; ++col) {
	 for (i = 0; i < nprops; ++i) {
	   compute_avg(grid[tbin][row][col].prop[i], 
		       grid[tbin][row][col].weightsum[i], nsiglevs);
	 } /* End of loop over all requested properties. */
	 compute_avg(grid[tbin][row][col].d, grid[tbin][row][col].dweightsum, 
		     nsiglevs);
	 define_sigma_transitions(&grid[tbin][row][col], nprops);
       } /* End of loop over all columns. */
     } /* End of loop over all rows. */
   } /* End of loop over all time bins. */


   /* Interpolate the sigma series back onto
    * the standard depth levels and output the
    * property data to the netcdf file... */
   fprintf(stderr,"  writing data ...\n");
   for (tbin = 0; tbin < ntbins; ++tbin) {
     for (row = 0; row < nrows; ++row) {
       for (col = 0; col < ncols; ++col) {
	 get_lat_lon(&h, row, col, &latitude, &longitude);
	 nlevs = do_std_depth(&grid[tbin][row][col], nsiglevs, prop_indx, 
			      nprops, ncols, col, data, count, &bottomdepth[col]);
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
   fprintf(stderr,"writing std depths a second time.\n");
   error = write_std_depths_cdf(cdf_file, &h);
   cdf_close(cdf_file);
   free((void *) data);
   free((void *) count);
   free((void *) bottomdepth);
   for (tbin = 0; tbin < ntbins; ++tbin) 
     free_grid(ncols, nrows, nprops, grid[tbin]);
   
   /* Now that memory has been freed up, search
    * for vertical data gaps in isopycnally averaged
    * grid and determine which result from pycnostads.  
    *
    * Strategy:
    * Re-read all the station files again, this time
    * averaging data on stddepth surfaces.  Then
    * re-open the cdf file, visit each grid node and
    * see if there are vertical data gaps which may
    * have resulted from weak vertical gradients.
    * Replace empty levels with isobarically averaged
    * datapoints if they exist. Rewrite the
    * updated property values to the cdf file. */
   fprintf(stderr,"\nChecking for pycnostads....");
   fprintf(stderr,"\n    allocating memory for grid ... " );
   for (tbin = 0; tbin < ntbins; ++tbin) 
     grid[tbin] = alloc_grid(ncols, nrows, nstdlevs, nprops);

   /* Initialize a data/header storage structure. */
   for (i = 0; i < MAXPROP; ++i)
     sta.observ[i] = (double *) NULL;
   hdr.prop_id = (int *) NULL;

   /* Loop for each input_file */
   fprintf(stderr,"\n    summing ...");
   curfile = 1;
   
   do {
     if ((infile = open_hydro_file(dir, argv[curfile], extent, print_mess)) 
	 < 0) {
       goto NEXTFILE2;
     }
      
     /* Loop for each station. */
     while ((error = get_station(infile, &hdr, &sta)) == 0) {

       /* Is station within bounds ?   */
       if (xdateline && hdr.lon < 0)
	 hdr.lon += 360;
       
       if (get_indices(&h, hdr.lat, hdr.lon, &row, &col) == 0) {

         pr = sta.observ[(int)PR];   /* set frequently used pointers */
         de = sta.observ[(int)DE];
         te = sta.observ[(int)TE];
         sa = sta.observ[(int)SA];
	 
	 ratio_done = 0;
         if (pr == NULL || te == NULL || sa == NULL || de == NULL) 
           continue;
         
         /* compute other appropriate properties for this station ... */
         for (i = 0; i < nprops; ++i) {
	 
	   /* !**! Special cases for individual properties... */
	   switch ((enum property) prop_indx[i]) {
               case OX:  
		 prop_avail = available((enum property) prop_indx[i], &hdr);
		 if (!prop_avail) {
		   if (available(O2, &hdr)) {
		     free_and_alloc(&sta.observ[(int)OX], hdr.nobs);
		     for (j=0; j < hdr.nobs; ++j) {
		       sta.observ[OX][j] = ox_kg2l(sta.observ[(int)O2][j], pr[j], te[j], sa[j]);
		     }
		     prop_avail = 1;
		   }
		 }
		 break;
               
               case O2:  
		 prop_avail = available((enum property) prop_indx[i], &hdr);
		 if (!prop_avail) {
		   if (available(OX, &hdr)) {
		     free_and_alloc(&sta.observ[(int)O2], hdr.nobs);
		     for (j=0; j < hdr.nobs; ++j) {
		       sta.observ[OX][j] = ox_l2kg(sta.observ[(int)OX][j], pr[j], te[j], sa[j]);
		     }
		     prop_avail = 1;
		   }
		 }
		 break;

               case S_:
		 free_and_alloc(&sta.observ[(int)S_], hdr.nobs);
		 compute_sigma(s_pref, hdr.nobs, sta.observ[(int)S_],pr,te,sa);
		 break;

               case S0:
		 free_and_alloc(&sta.observ[(int)S0], hdr.nobs);
		 compute_sigma(0.0, hdr.nobs, sta.observ[(int)S0],pr,te,sa);
		 break;

               case S1:
		 free_and_alloc(&sta.observ[(int)S1], hdr.nobs);
		 compute_sigma(1000.0, hdr.nobs, sta.observ[(int)S1],pr,te,sa);
		 break;

               case S2:
		 free_and_alloc(&sta.observ[(int)S2], hdr.nobs);
		 compute_sigma(2000.0, hdr.nobs, sta.observ[(int)S2],pr,te,sa);
		 break;

               case S3:
		 free_and_alloc(&sta.observ[(int)S3], hdr.nobs);
		 compute_sigma(3000.0, hdr.nobs, sta.observ[(int)S3],pr,te,sa);
		 break;

               case S4:
		 free_and_alloc(&sta.observ[(int)S4], hdr.nobs);
		 compute_sigma(4000.0, hdr.nobs, sta.observ[(int)S4],pr,te,sa);
		 break;

               case SD:
		 free_and_alloc(&sta.observ[(int)SD], hdr.nobs);
		 compute_dref_sigma(de, hdr.nobs, sta.observ[(int)SD],pr,te,sa);
		 break;

               case TH:
		 free_and_alloc(&sta.observ[(int)TH], hdr.nobs);
		 compute_theta(hdr.nobs, sta.observ[(int)TH], pr, te, sa);
		 break;

               case HC:
		 free_and_alloc(&sta.observ[(int)HC], hdr.nobs);
		 compute_heat_capacity(hdr.nobs, sta.observ[(int)HC], pr, te, sa);
		 break;

               case HT:
		 free_and_alloc(&sta.observ[(int)HT], hdr.nobs);
		 compute_height(hdr.nobs, pr, te, sa, ht_pref, sta.observ[(int)HT]);
		 break;

               case PE:
		 free_and_alloc(&sta.observ[(int)PE], hdr.nobs);
		 compute_energy(hdr.nobs, pr, te, sa, pe_pref, sta.observ[(int)PE]);
		 break;
		 
               case DR:
		 if (!ratio_done) {
		   free_and_alloc(&sta.observ[(int)DR], hdr.nobs);
		   free_and_alloc(&sta.observ[(int)AL], hdr.nobs);
		   free_and_alloc(&sta.observ[(int)BE], hdr.nobs);
		   compute_ratio(hdr.nobs, sta.observ[(int)DR], pr, te, sa, sta.observ[(int)AL], sta.observ[(int)BE]);
		   ratio_done = 1;
		 }
		 break;

	       case SV:
		 free_and_alloc(&sta.observ[(int)SV], hdr.nobs);
		 compute_sp_vol(hdr.nobs, sta.observ[(int)SV], pr, te, sa);
		 break;

               case VA:
		 free_and_alloc(&sta.observ[(int)VA], hdr.nobs);
		 compute_svan(hdr.nobs, sta.observ[(int)VA], pr, te, sa);
		 break;

               case VS:
		 free_and_alloc(&sta.observ[(int)VS], hdr.nobs);
		 compute_sound_vel( sta.observ[(int)VS], pr, te, sa, hdr.nobs);
		 break;


               case BF:   
                  free_and_alloc(&sta.observ[(int)BF], hdr.nobs);
		  buoy_freq(sta.observ[(int)BF],pr,te,sa,hdr.nobs,window, w_incr);
                  for (j = 0; j < hdr.nobs; ++j) {
                    if (sta.observ[(int)BF][j] > -9998.)
                      sta.observ[(int)BF][j] *= 1.0e5;
                  }
                  break;

               case PV:
		 free_and_alloc(&sta.observ[(int)PV], hdr.nobs);
		 buoy_freq(sta.observ[(int)PV],pr,te,sa,hdr.nobs,window, w_incr);
		 for (j = 0; j < hdr.nobs; ++j) 
		   sta.observ[(int)PV][j] *= sta.observ[(int)PV][j];
		 po_vort(sta.observ[(int)PV],sta.observ[(int)PV], hdr.nobs, (double)hdr.lat);
		 break;
		 
	       default:
		 prop_avail = available((enum property) prop_indx[i], &hdr);
		 break;

	   } /* End property selection switch. */

	   for (tbin = 0; tbin < ntbins; ++tbin) {
	     if (hdr.year >= h.tmin[tbin] && hdr.year <= h.tmax[tbin]) {

	       if (prop_avail) {
		 sum_levels(sta.observ[prop_indx[i]], de, hdr.nobs, grid[tbin][row][col].prop[i], grid[tbin][row][col].count[i]);
	       } /* End prop_avail check. */
	     } /* End tbin sorting check. */
	   } /* End for loop over all tbins. */
         } /* End for loop over i = all output/requested properties. */
       } /* End if check on get_indices. */
       
       /* Free up space. */
       for (i = 0; i < MAXPROP; ++i) {
	 if (sta.observ[i] != NULL) {
	   free((void *) sta.observ[i]);
	   sta.observ[i] = (double *) NULL;
	 }
       }
       
     }  /* End while looping over all stations in a file. */ 
     if (error > 0) {
       report_status(error, stderr);
       exit(1);
     }
     
   NEXTFILE2:
     close(infile);
     
   } while (curfile++ < nfiles ); /* End looping over input files. */
   
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
   
   fprintf(stderr,"\n  filling pycnostads ...");
   n_filled = 0;
   
   /* Check each isopycnally averaged gridnode.
    * If depth level is flagged as missing, check
    * the isobarically averaged data to
    * differentiate a vertical datagap from a
    * pycnostad. */
   for (tbin = 0; tbin < ntbins; ++tbin) {
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
   } /* end for tbin*/
   
   fprintf(stderr,"\n  %d levels were identified and filled", n_filled);
   fprintf(stderr,"\n End of %s\n", argv[0]);
   
   exit(0);

} /* end main */

/****************************************************************************/

void print_usage(char *program) {
   fprintf(stderr,"\nComputes an average profile for each property specified");
   fprintf(stderr,"\nat each lat/lon gridnode from randomly spaced stations.");
   fprintf(stderr,"\nInput profiles are weighted by distance from gridnode center");
   fprintf(stderr,"\nunless weighting is turned off with -N.");
   fprintf(stderr,"\nOutputs a netCDF file of hydro properties gridded as a ");
   fprintf(stderr,"\nfunction of lat, lon, and depth.  The depth dimension ");
   fprintf(stderr,"\nconsists of a series of standard depths, optionally ");
   fprintf(stderr,"\nspecified with -Z. Averaging is done on isopycnal surfaces");
   fprintf(stderr,"\nand interpolated back onto the depth levels.  The user");
   fprintf(stderr,"\nmust specify these density surfaces with -S<pref><filename>.");
   fprintf(stderr,"\nThe sigma values should be customized to reflect the");
   fprintf(stderr,"\nlocal vertical density stratification.");
   fprintf(stderr,"\n");   
   fprintf(stderr,"\nUsage:  %s filename_root(s)  -B/west/east/south/north -C<cdf_output_file> -I<gridspacing> -P<list_of_properties> -S<ref_id>/<sigma_series_file>  [-D<dirname>] [-E<file_extent>]  [-L<lengthscale>]  [-M<mixed_layer_def>] [-T<tmin>/<tmax>] [-W<window>[/<w_incr>]]  [-Z<std_depth_file>]  \n", program);
   fprintf(stderr,"\n -B   specifies grid bounds");
   fprintf(stderr,"\n -C   name of netCDF output file.");
   fprintf(stderr,"\n -I   specifies grid increment in degrees;  ex: -I0.5");
   fprintf(stderr,"\n          OR specify separate x,y increments with a slash");
   fprintf(stderr,"\n          to delimit xincr/yincr;  ex: -I2.0/0.5\n");
   fprintf(stderr,"\n -P   list of properties to project onto surface;");
   fprintf(stderr,"\n          ex:  -Ppr/th/sa/ox/ht");
   fprintf(stderr,"\n               -P (by itself) produces a list of available properties\n");
   fprintf(stderr,"\n -S   ref_id = 0,1, 2,3, or 4 ");
   fprintf(stderr,"\n          Specify EACH reference pressure with its own -S argument");
   fprintf(stderr,"\n          filename = file containing sigma series definitions.");
   fprintf(stderr,"\n          Each line in file contains sigmin, sigmax, incr");
   fprintf(stderr,"\n          (sigmin and sigmax are INCLUSIVE in generating series).");
   fprintf(stderr,"\n          Values MUST be monotonically INCREASING.");
   fprintf(stderr,"\n          If providing path_and_file referenced to root,");
   fprintf(stderr,"\n          be certain to include an / between the -S<ref_id>");
   fprintf(stderr,"\n          and the root-referenced pathname.");
   fprintf(stderr,"\n\n OPTIONS:");
   fprintf(stderr,"\n -D  specifies directory for input data files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/\n ");
   fprintf(stderr,"\n -E  specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat \n");
   fprintf(stderr,"\n -L lengthscale (L, in km) for distance weighting =  e ^[-(dist/L)^2]");
   fprintf(stderr,"\n     -L0 turns off distance weighting (all points contribute equally to mean)" );
   fprintf(stderr,"\n          default is L =%.1f  km", lengthscale);
   fprintf(stderr,"\n     The weight is computed based on each observation relative to the center");
   fprintf(stderr,"\n     of the horizontal bin that holds the observation.  Observations from");
   fprintf(stderr,"\n     neighboring bins do not contribute - see hb_smooth3d for that.");
   fprintf(stderr,"\n -M  specify the definition of the mixed layer as ");
   fprintf(stderr,"\n          sigma(sea surface) + this value.  default: %f ", MIX_LAYER_DEF);
   fprintf(stderr,"\n          ex:  -M.02");
   fprintf(stderr,"\n -T  optional minyear/maxyear to constrain time interval of observations.");
   fprintf(stderr,"\n -W  Specifies pressure window (db) for computing ");
   fprintf(stderr,"\n          gradient properties (bvfreq, pv...) This ");
   fprintf(stderr,"\n          constitutes the range over which observations");
   fprintf(stderr,"\n          are incorporated into the gradient computation.");
   fprintf(stderr,"\n          The window is centered around each pressure ");
   fprintf(stderr,"\n          level being evaluated.");
   fprintf(stderr,"\n          w_incr specifies how finely to subdivide the");
   fprintf(stderr,"\n          window into increments(db) to create");
   fprintf(stderr,"\n          an evenly spaced pressure series over the window.");
   fprintf(stderr,"\n          defaults: -W100/10");
   fprintf(stderr,"\n -Z  file containing list of standard depths.");
   fprintf(stderr,"\n          Values MUST be monotonically INCREASING.\n");
   fprintf(stderr,"\n          No interpolation over data gaps of more than:");
   fprintf(stderr,"\n          %d m above %d m",GAP_SHALLOW,GAP_CHANGE_DEPTH);
   fprintf(stderr,"\n          %d m below %d m",GAP_DEEP,GAP_CHANGE_DEPTH);
   fprintf(stderr,"\n -h  help ... prints this message.");
   fprintf(stderr,"\n\n");  
   return;
} /* End print_usage() */

/****************************************************************************/
void parse_s_option(char *str, FILE **file_ptr) {

  /* get the sigma values defining the isopycnal surfaces */

  int j;

  /* Check which reference level the file is for. */
  switch (*str) {

    case '0':
      j = 0;
      break;

    case '1':
      j = 1;
      break;

    case '2':
      j = 2;
      break;

    case '3':
      j = 3;
      break;

    case '4':
      j = 4;
      break;

    default:
      fprintf(stderr,"\n Error parsing -S option\n");
      fprintf(stderr,"USAGE: -S<ref_lev_id>/<file_name>");
      fprintf(stderr,"\n ex: -S4/sig4levs.list\n");
      fprintf(stderr,"\n Exiting, end of hb_grid3d.");
      exit(1);
      break;

  }  /* End switch */

  /* Move the pointer to the string one character forward
   * if there is a slash since the syntax for the option
   * is -S<0|1|2|3|4>/<file_name>.  Note that for file
   * names referenced to root, this will remove the root
   * and cause a error in reading if there is no extra
   * "/" between the -S0 and the full root-referenced path.*/
  if (*(++str) == '/') {
    ++str;
  }

  file_ptr[j] = fopen(str,"r");
  if (file_ptr[j] == NULL) {
    fprintf(stderr,"\n Error opening %s \n", str);
    fprintf(stderr,"\n Exiting, end of hb_grid3d.");
    exit(1);
  }

  return;

}   /* End parse_s_option() */

/*****************************************************************************/
int parse_p_option(char *st, int *prop_indx) {

  /* This routine reads the command line list
   * of requested properties and returns the
   * number of properties requested as well
   * as loading the prop IDs into the
   * prop_indx array. */

  char prop[4];
  int n;
  double reflev;

  if (*st == '\0') {
    print_prop_menu();
    exit(0);
  }

  n = 0;
  do {
    if (*st == '/')
      ++st;
    sscanf(st,"%2s", prop);
    prop_indx[n] = get_prop_indx(prop);
    if (prop_indx[n] < 0)  {
      fprintf(stderr,"\n Unknown property '%.2s' specified in -P%s\n", prop, st);
      exit(1);
    }
          
    ++st;
    ++st;

    /* !**!  Special cases for properties */ 
    if ((prop_indx[n] == (int)S_)  || (prop_indx[n] == (int)PE) || (prop_indx[n] == (int)HT) ) {
      if (sscanf(st, "%lf", &reflev) != 1) {
	fprintf(stderr, "\n Specify a ref pressure for  %.2s", prop);
	fprintf(stderr,"\n   ex: -P%.2s1500/th/sa\n", prop);
	exit(1);
      }
      switch ((enum property) prop_indx[n]) {
        case S_:
	  s_pref = reflev;
	  break;

        case PE:
	  pe_pref = reflev;
	  break;

        case HT:
	  ht_pref = reflev;
	  break;

        default:
	  ;
      }  /* end switch */
      while (!(*st=='/' || *st==0 || *st==' '))
	++st;
    }
    
    if ((prop_indx[n] == (int) GE) || (prop_indx[n] == (int) GN)) {
      fprintf(stderr,"\n WARNING: It is not appropriate to average");
      fprintf(stderr,"\n neutral density (gamma-n) isopycnally.");
      fprintf(stderr,"\n Compute this quantity from averaged");
      fprintf(stderr,"\n pr,te,sa output from hb_grid3d.\n");
    }
     
    if ((prop_indx[n] != (int) GE) && (prop_indx[n] != (int) GN)
        &&(prop_indx[n] != (int) DE))  {
      /* depth is automatically done so don't count it here */
      ++n;
    }

  } while (*st == '/'); /* End of looping on command line characters. */
   
  /* Return the number of properties requested. */
  return (n);

}  /* End parse_p_option() */

/*****************************************************************************/
struct gridnode **alloc_grid(int nx, int ny, int nz, int nprops) {

  /* Given the number of columns (nx), rows (ny),
   * depth levels (nz), and requested properties
   * (nprops), return a 2D array of arrays. */

  int i, j, k, n;      /* Local counters. */
  struct gridnode **g; /* Local holder of 2D array of arrays grid,
                        * will be the returned value. */

  /* Meridional allocation. */
  g = (struct gridnode **) malloc(ny * sizeof(struct gridnode *));

  /* Check that allocation was sucessful. */
  if (g == NULL) {
    fprintf(stderr,"\nUnable to allocate memory for grid.\n");
    exit(1);
  }

  /* For each meridional location,
   * (1) do the zonal allocation of data.
   * (2) check for sucessful allocation of each zonal strip.
   * (3) allocate space for each property/depth
   *     (storage of values, weightsum, count,
   *      depths, depth weightsum, nobs). */
  for (i = 0; i < ny; ++i) {
    g[i] = (struct gridnode *) malloc(nx * sizeof(struct gridnode));
    if (g[i] == NULL) {
      fprintf(stderr,"\nUnable to allocate memory for grid[%d]\n", i);
      exit(1);
    }

    for (j = 0; j < nx; ++j) {

      /*=========Property-dimension allocations===========*/
      /* prop, weightsum, and count */

      g[i][j].prop = (double **) malloc(nprops * sizeof(double *));
      if (g[i][j].prop == NULL) {
	fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].prop\n", i, j);
	exit(1);
      }

      g[i][j].weightsum = (double **) malloc(nprops * sizeof(double *));
      if (g[i][j].weightsum == NULL) {
	fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].weightsum\n", i, j);
	exit(1);
      }
      
      g[i][j].count = (UI **) malloc(nprops * sizeof(UI *));
      if (g[i][j].count == NULL) {
	fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].count\n", i, j);
	exit(1);
      }

      /*===========Depth-dimension allocations============*/
      /* d, dweightsum, nobs, mix_layer, deepest */
      g[i][j].d = (double *) malloc(nz * sizeof(double));
      if (g[i][j].d == NULL) {
	fprintf(stderr,"\nUnable to allocate memory for g[%d][%d].d\n", i, j);
	exit(1);
      }

      g[i][j].dweightsum = (double *) malloc(nz * sizeof(double));
      if (g[i][j].dweightsum == NULL) {
	fprintf(stderr,"\nUnable to allocate memory for g[%d][%d].dweightsum\n", i, j);
	exit(1);
      }

      g[i][j].nobs = (UI *) malloc(nz * sizeof(UI));
      if (g[i][j].nobs == NULL) {
	fprintf(stderr,"\nUnable to allocate memory for g[%d][%d].nobs\n", i, j);
	exit(1);
      }
      
      g[i][j].mix_layer = (struct surfrec *) NULL;
      g[i][j].deepest = (struct deepestrec *)NULL;
      
      /* Initialize depth, weightsum, obs count arrays */
      for (n = 0; n < nz; ++n) {
	g[i][j].d[n] = 0.0;
	g[i][j].dweightsum[n] = 0.0;
	g[i][j].nobs[n] = 0;
      }

      /*=========Property storage allocations========*/
      /* prop, weightsum, count 
       * These three variables are 4D arrays (x,y,z,prop). */
      for (k = 0; k < nprops; ++k) {

	g[i][j].prop[k] = (double *) malloc(nz * sizeof(double));
	if (g[i][j].prop[k] == NULL) {
	  fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].prop[%d]\n", i, j, k);
	  exit(1);
	}

	g[i][j].weightsum[k] = (double *) malloc(nz * sizeof(double));
	if (g[i][j].weightsum[k] == NULL) {
	  fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].weightsum[%d]\n", i, j, k);
	  exit(1);
	}

	g[i][j].count[k] = (UI *) malloc(nz * sizeof(UI));
	if (g[i][j].count[k] == NULL) {
	  fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].count[%d]\n", i, j, k);
	  exit(1);
	}
              
	/* zero each prop, count */
	for (n = 0; n < nz; ++n) {     
	  g[i][j].prop[k][n] = 0.0;
	  g[i][j].weightsum[k][n] = 0.0;
	  g[i][j].count[k][n] = 0;
	}
      } /* End for loop over all properties. */
    } /* End for loop over all zonal points. */
  }  /* End for loop over all meridional points. */

  return g;
  
} /* End alloc_grid */

/*****************************************************************************/
void free_grid(int nx, int ny, int nprops, struct gridnode **gridptr) {

  /* Frees all the memory for ny * nx
   * struct gridnodes.  This does the
   * inverse of alloc_grid. */

  int row, col, k;
  struct gridnode *g;
   
  for (row = 0; row < ny; ++row) {
    for (col = 0; col < nx; ++col) {
      g = &gridptr[row][col];
      for (k = 0; k < nprops; ++k) {
	if (g->prop[k] != NULL)
	  free((void *)g->prop[k]);
	if (g->weightsum[k] != NULL)
	  free((void *)g->weightsum[k]);
	if (g->count[k] != NULL)
	  free((void *)g->count[k]);
      } /* End for each property. */
      free((void *)g->prop);
      free((void *)g->weightsum);
      free((void *)g->count);
      free((void *)g->nobs);
      free((void *)g->d);
      free((void *)g->dweightsum);
      delete_surfrec_list(g->mix_layer);
      delete_list(g->deepest);
    } /* End for each col. */
    
    free((void *) gridptr[row]);
  } /* End for each row. */
  
  free((void *) gridptr);
  
  return;

}  /* End free_grid() */

/*****************************************************************************/
int get_sig_series(int ref_id, FILE *fptr, double *siglist) {
  
  /* The file will be read for (min, max, incr)
   * triplets from which a sigma series will be
   * generated and inserted at siglist.
   * The number of sigma values is returned. 
   *  
   * Arguments:
   *
   *  ref_id:   defines the sigma variable: (int) enum property 
   *    fptr:   pointer to file containing list of sigma values OR nil 
   * siglist:   starting addr of array to insert values (the address
   *            we start writing to, not the head of the array).
   */

  double  min, max, incr;
  double  next;
  int i, count;
  int bigcount;

  bigcount = nsiglevs;

  /* Initialize a counter. */
  count = 0;

  /* Scan each line of the file (triplets). */
  while( fscanf(fptr,"%lf %lf %lf", &min, &max, &incr) == 3) {

    /* Check that we have not loaded too many levels.*/
    if (++bigcount >= MAXSIGLEVS) { /* to avoid a SEG FAULT */
      return (bigcount);
    }

    /* Write the min value of the line to the siglist. */
    *siglist = min;

    /* Augment the number of sigma values added to the list. */
    ++count;

    /* Compute all the intermediate steps, separated by incr. */
    while ( ( next = *siglist + incr) <= max) {

      /* Augment the number of sigma values added to the list. */
      ++count;

      /* Check that we have not loaded too many levels.*/
      if (++bigcount >= MAXSIGLEVS) { /* to avoid a SEG FAULT */
	return (bigcount);
      }

      /* Move the address pointer to the list forward
       * to the next element as well as writing the
       * next sigma value to the list. */
      *(++siglist) = next;
    } /* End of while loop over each newly written value. */

    /* Augment the counter to the next empty position. */
    ++siglist;
  } /* End of while loop over all lines in input sigfile. */
  
  /* Close the input sigfile. */
  fclose(fptr);

  /* Return the number of values added to the siglist. */
  return (count);
  
} /* End get_sig_series */

/*****************************************************************************/
int get_time_bins(int minyr, int maxyr, struct CDF_HDR *hptr) {

  /* Sets time bin min/max values in the CDF header */
  hptr->tmin = (int *) malloc( sizeof(int));
  hptr->tmax = (int *) malloc( sizeof(int));

  hptr->tmin[0] = minyr;        
  hptr->tmax[0] = maxyr;

  return (1);
 } /* End get_time_bins */

/*****************************************************************************/
double get_weight(float lat1, float lon1, float lat2, float lon2, float L) {

  /* Returns weight as a function of
   * distance between points:
   * weight = e ^[-(d/L)^2]  */

   double dx, dy, dist;
   
   /* Find the distance between the given grid point
    * and the location of the station. */
   dx = (lon2 - lon1)* RADperDEG * cos(RADperDEG*.5 *(lat1+lat2)) * EarthRadius ;
   dy = (lat2 - lat1) * RADperDEG * EarthRadius;
   dist = sqrt(dx*dx + dy*dy);

   /* Convert distance between points to the
    * ratio of the distance relative to the
    * user-specified length scale. */
   dist = dist /L;   
   dist = dist * dist;
   
   return (exp( -dist));
} /* End get_weight() */	
 
/*****************************************************************************/
void insert_data(double *y, int ny, double *ysig, double *wghtsum, double wght, UI *count, int test_outcrop) {

  /* For the y property at a station with
   * ny observation levels, interpolate
   * to find the yvals at each level in
   * the sigma series,  weight each yval,
   * add it to its appropriate sum, sum the
   * weights and increment the counter.
   * Checks for vertical data gaps using
   * the globally defined de array and does
   * not interpolate over gaps which exceed
   * 200 m in the thermocline (upper 1000 m)
   * or 1000 m elsewhere.  The sigma series
   * is specified in the global array, siglevs.
   * The observed sigma values at this station
   * are subdivided by ref levels, and stored
   * at global addresses pointed to by sigptr. 
   *  
   *            y:    array of y observed values (vector in depth)
   *           ny:    dimension of y (scalar number of depth surfaces)
   *         ysig:    sum of distance-weighted yvalues on each sigma surface
   *      wghtsum:    sum of weights
   *         wght:    weight for this profile (scalar)
   *        count:    counts # of observations on each surface 
   * test_outcrop:    0 or 1: if set, tests for outcropping surfaces
   */

  int  i, j, datagap, n;
  int  refindx;
  double *xtmp[NREFLEVS], *ytmp, *dtmp, z;
  double reflev;
   
  for (i = 0; i < NREFLEVS; ++i ) {
    xtmp[i] = (double *) malloc((UI) (ny * sizeof(double)));
  }
  ytmp = (double *) malloc((UI) (ny * sizeof(double)));
  dtmp = (double *) malloc((UI) (ny * sizeof(double)));

  if (dtmp == NULL) {
    fprintf(stderr,"\nUnable to allocate memory in insert_data()\n");
    exit(1);
  }

  /* Ensure continuous (no missing values) x, y, and depth arrays */
  n = 0; /* Initialize observation count. */

  /* Loop over all station points. */
  for (i = 0; i < ny; ++i) {
    /* Test for empty value. */
    if ( y[i] > -8.9)  {
      /* We have good values, so for each
       * density surface, store the density
       * surface value. */
      for (j = 0; j < NREFLEVS; ++j ) {
	xtmp[j][n] = sigptr[j][i];
      } /* End of loop over all reference levels. */

      /* Store the good observation and its depth. */
      ytmp[n] = y[i];
      dtmp[n++] = de[i];
    } /* End of empty value check. */ 
  } /* End of loop over all station observations. */
   
  /* Check whether we have enough values. */
  if (n <= 1) {
    /* Not enough to average! Free memory and exit. */
    for (i = 0; i < NREFLEVS; ++i) {
      free((void *)xtmp[i]);
    }
    free((void *)ytmp);
    free((void *)dtmp);
    return;
  }

  /* Loop over all density surfaces. */
  for (i = 0; i < nsiglevs; ++i) {

    /* Determine which ref lev corresponds to this sigseries element ... */
    refindx = 0;      
    j = nsig[0];
    while (i > j) {
      if (++refindx == NREFLEVS-1)
	break;
      j += nsig[refindx];
    }

    /* Compute reference pressure. */
    reflev = refindx * 1000.;
    
    /* Test for surface outcropping. */
    if ((test_outcrop) && (siglevs[i] < xtmp[refindx][0])) {   
      /* NOTE:  this assumes x (density) increases with depth!! */           
      ysig[i] += -1 * wght;     /* use depth= -1 for outcropping surfaces */
      wghtsum[i] += wght;
      ++count[i];
    }

    /* See if associated depth exists at this station. */
    else if ((z = interpolate(siglevs[i], xtmp[refindx], dtmp, n)) > -998.) {

      /* Now check for datagaps around this depth...*/
      j = 0;
      while ((++j < n) &&( dtmp[j] < z)  ) 
	;
      if (j == n)
	--j;
      
      if ((dtmp[j-1] == z) || (dtmp[j] == z) )
	datagap = 0;
      else if (z < GAP_CHANGE_DEPTH)
	datagap = (dtmp[j] - dtmp[j-1]) > GAP_SHALLOW;
      else
	datagap = (dtmp[j] - dtmp[j-1]) > GAP_DEEP;
      
      /* Exclude observations which are
       * more than 1500 m above the ref
       * pressure. This prevents alot of
       * bogus values from being added. */
      datagap = datagap || (z < reflev-1500);
      
      /* For values that are not affected by
       * data gaps, compute the running weighted
       * sum of values, the running sum of
       * weights, and count the number of obs
       * on this surface. */
      if (!datagap) {
	if ((z = interpolate(siglevs[i], xtmp[refindx], ytmp, n)) > -998.) {
	  ysig[i] += z * wght;
	  wghtsum[i] += wght;
	  ++count[i];
	}
      } /* End of datagap check. */
    } /* Check if station has depth at this surface - if not, do nothing. */
  } /* End of looping over density surfaces. */

  /* Free tmp memory. */
  for (i = 0; i < NREFLEVS; ++i) {
    free((void *)xtmp[i]);
  }
  free((void *)ytmp);
  free((void *)dtmp);
  return;
} /* End insert_data() */

/*****************************************************************************/
void insert_deepest_vals(struct gridnode *g_ptr, int *prop_indx, int nprops) {

  /* Maintains a linked list of info regarding
   * the deepest observed level at a grid point. 
   * All observations within the deepest 100
   * meters will be retained.  The list will
   * actively grow and contract as each new
   * station is read. The station data are 
   * accessed through the global variable:
   * 
   * struct HYDRO_DATA sta 
   *
   * Inputs:
   *
   *       g_ptr:   ptr to current gridnode  
   *   prop_indx:   contains list of indices to props being gridded 
   *      nprops:   number of properties being gridded 
   */

  int i, j, reflev;
  int density_greater, depth_greater, depth_diff;
  double key;
  struct deepestrec  *dr_ptr, *r1, *r2;
 
  dr_ptr = g_ptr->deepest; 
  j = hdr.nobs - 1;
 
  /*  Empty list case ...*/ 
  if (dr_ptr == NULL) {
    dr_ptr = create_deepestrec(nprops);
    
    dr_ptr->depth = sta.observ[(int)DE][j];
    dr_ptr->pressure = sta.observ[(int)PR][j];
    dr_ptr->sig_0 = sta.observ[(int)S0][j];
    dr_ptr->sig_1 = sta.observ[(int)S1][j];
    dr_ptr->sig_2 = sta.observ[(int)S2][j];
    dr_ptr->sig_3 = sta.observ[(int)S3][j];
    dr_ptr->sig_4 = sta.observ[(int)S4][j];
      
    for (i=0; i< nprops; ++i) {
      dr_ptr->prop[i] = (double) HBEMPTY;
      if (sta.observ[prop_indx[i]] != NULL) {
	if (sta.observ[prop_indx[i]][j] > -8.9) 
	  dr_ptr->prop[i] = sta.observ[prop_indx[i]][j];
      }
    }
    dr_ptr->next = NULL;
    g_ptr->deepest = dr_ptr;
    return;
  }
  
  /* Determine appropriate sigma reference level. */
  reflev = 0;
  for (i = 1; i < NREFLEVS; ++i) {
    if (sta.observ[(int)PR][j] > zmin[i] && sta.observ[(int)PR][j] <= zmax[i])
      reflev = i;
  }
      
  /* Compare to existing deepest record. */   
  dr_ptr = g_ptr->deepest;
  depth_greater = sta.observ[(int)DE][j] > dr_ptr->depth ? 1 : 0;
  depth_diff =  sta.observ[(int)DE][j] - dr_ptr->depth;
  switch (reflev){
      case 0:
	density_greater = sta.observ[(int)S0][j] > dr_ptr->sig_0 ? 1 : 0;
	break;

      case 1:
	density_greater = sta.observ[(int)S1][j] > dr_ptr->sig_1 ? 1 : 0;
	break;

      case 2:
	density_greater = sta.observ[(int)S2][j] > dr_ptr->sig_2 ? 1 : 0;
	break;

      case 3:
	density_greater = sta.observ[(int)S3][j] > dr_ptr->sig_3 ? 1 : 0;
	break;

      case 4:
	density_greater = sta.observ[(int)S4][j] > dr_ptr->sig_4 ? 1 : 0;
	break;

      default:
	fprintf(stderr,"FATAL ERROR in insert_deepest_vals():  reflev is %d", reflev);
	exit(1); 
  }

  /* It's denser, but considerably more
   * shallow (can happen in Greenland Sea),
   * then skip it. */
  if (density_greater && (depth_diff <= -500)) 
    return;
  
  /* If station is denser than previous densest:
   * or much deeper than previous densest
   * insert new rec at beginning of list
   * check remainder of list and delete
   * any recs > 100 m shallower */
  if (density_greater  || (depth_diff > 1000) ) {
    dr_ptr = create_deepestrec(nprops);
    dr_ptr->depth = sta.observ[(int)DE][j];
    dr_ptr->pressure = sta.observ[(int)PR][j];
    
    dr_ptr->sig_0 = sta.observ[(int)S0][j];
    dr_ptr->sig_1 = sta.observ[(int)S1][j];
    dr_ptr->sig_2 = sta.observ[(int)S2][j];
    dr_ptr->sig_3 = sta.observ[(int)S3][j];
    dr_ptr->sig_4 = sta.observ[(int)S4][j];
    for (i=0; i< nprops; ++i) {
      dr_ptr->prop[i] = (double) HBEMPTY;
      if (sta.observ[prop_indx[i]] != NULL) {
	if (sta.observ[prop_indx[i]][j] > -8.9)
	  dr_ptr->prop[i] = sta.observ[prop_indx[i]][j];
      }
        
      if ((prop_indx[i] == (int)PR)  && !depth_greater)
	dr_ptr->prop[i] =g_ptr->deepest->pressure;
    }
    
    if ( !depth_greater) {
      dr_ptr->depth = g_ptr->deepest->depth;
      dr_ptr->pressure = g_ptr->deepest->pressure;
    }
      
      
    dr_ptr->next = g_ptr->deepest;
    
    key = dr_ptr->depth - 100.;
    r1 = dr_ptr;
    r2 = dr_ptr->next;
    while (r2 != NULL) {
      if (r2->depth < key) {
	delete_list(r2->next);
	free((void *) r2->prop);
	free((void *)r2);
	r2 = NULL;
	r1->next = NULL;
      }
      else {
	r1 = r2;
	r2 = r2->next;
      }
    }
    g_ptr->deepest = dr_ptr;
    return;
  }
   
  /* It's not denser, but if it is deeper
   * than the densest observation so far, 
   * replace the deepest depth/pressure and
   * then go on to insert a new record 
   * into the linked list. */
  if (depth_greater) {
    r1 = g_ptr->deepest;
    r1->depth = sta.observ[(int)DE][j];
    r1->pressure = sta.observ[(int)PR][j];
  }
  
  /* It is not the densest, if it is
   * shallower than 100 m above the
   * densest, we're done. */
  key = g_ptr->deepest->depth -100.;

  if (sta.observ[(int)DE][j]  < key)  
    return;
      
  /*  It is within 100 m of densest observation, insert it into linked list*/  
  dr_ptr = create_deepestrec(nprops);
  dr_ptr->depth = sta.observ[(int)DE][j];
  dr_ptr->pressure = sta.observ[(int)PR][j];
  dr_ptr->sig_0 = sta.observ[(int)S0][j];
  dr_ptr->sig_1 = sta.observ[(int)S1][j];
  dr_ptr->sig_2 = sta.observ[(int)S2][j];
  dr_ptr->sig_3 = sta.observ[(int)S3][j];
  dr_ptr->sig_4 = sta.observ[(int)S4][j];

  for (i=0; i< nprops; ++i) {
    dr_ptr->prop[i] = (double) HBEMPTY;
    if (sta.observ[prop_indx[i]] != NULL) {
      if (sta.observ[prop_indx[i]][j] > -8.9) 
	dr_ptr->prop[i] = sta.observ[prop_indx[i]][j];
    }
  }
  
  r1 = g_ptr->deepest;
  r2 = r1->next;
  do {
    if (r2 == NULL) {   /* insert at end of list */
      r1->next = dr_ptr;
      dr_ptr->next = NULL;
      return;
    }
    
    if (dr_ptr->sig_2 > r2->sig_2) { /*insert before r2 */
      r1->next = dr_ptr;
      dr_ptr->next = r2;
      
      /* and check remainder of list against key */
      
      while (r2 != NULL) {
	if (r2->depth < key) {
	  delete_list(r2->next);
	  free((void *) r2->prop);
	  free((void *)r2);
	  r2 = NULL;
	  r1->next = NULL;
	}
	else {
	  r1 = r2;
	  r2 = r2->next;
	}
      }
      return;
    }
    
    r1 = r2;
    r2 = r2->next;
    
  } while (1);    
} /* End insert_deepest_vals () */

/*****************************************************************************/
struct deepestrec * create_deepestrec(int nprops)
{
   struct deepestrec *r1;
   
   r1 = (struct deepestrec *) malloc(sizeof(struct deepestrec));
   if (r1 == NULL) {
       fprintf(stderr,"\nOut of memory in create_deepestrec() \n");
       exit(1);
   }
   r1->prop = (double *) calloc((size_t)nprops, sizeof(double));
   if (r1->prop == NULL) {
       fprintf(stderr,"\nOut of memory in create_deepestrec() \n");
       exit(1);
   }
  
   return(r1);
   
} /* end create_deepestrec() */

/*****************************************************************************/
void delete_list(struct deepestrec *rptr) {

  /* Recursively traverses a linked list
   * of records and frees up all the memory
   * allocated to those records. */

  /* end of list... */
  if (rptr == NULL)
      return;
  
  /* recursive part... */
  delete_list(rptr->next);
  free((void *) rptr->prop);
  free((void *) rptr);
  return;
  
}  /* End delete_list() */

/*****************************************************************************/
void delete_surfrec_list(struct surfrec *rptr) {

  /* Recursively traverses a linked list
   * of records and frees up all the memory
   * allocated to those records. */
  
  /* end of list... */
  if (rptr == NULL)
    return;
      
  /* recursive part... */
  delete_surfrec_list(rptr->next);
  free((void *) rptr->prop);
  free((void *) rptr->wghtsum);
  free((void *) rptr);
  return;
  
}  /* End delete_surfrec_list() */

/*****************************************************************************/
void insert_surf_vals(struct gridnode *g, int *prop_indx, int nprops)
   /* defines the depth of a mixed layer at the sea surface, determines 
      the average values of properties within it, and inserts the information
      into a linked list of records sorted by density. The station data are 
      accessed through the global variable: struct HYDRO_DATA sta 
              g:   ptr to gridnode  
      prop_indx:   contains index to properties being gridded 
         nprops:   number of properties being gridded 
   */
{
   int i, j, n, weight;
   double x, v, vprev, dprev;
   double dens, depth;
   struct surfrec *m;

   if (de[0] > 10.0)
       return;

   /* round off surface density to nearest tenth */

   dens =  (double) (NINT(sta.observ[(int)S0][0] * 10.)) / 10.0;

  /* bottom of mixed layer is defined where 
    density = surface density + mix_layer_def */  
     
   depth = interpolate((sta.observ[(int)S0][0]+ mix_layer_def), sta.observ[(int)S0], de, hdr.nobs);

   if (depth < -999.)  /* all observations are within defined range of surface density */
       depth = de[hdr.nobs-1];

   m = get_surfrec(g->mix_layer, dens, nprops);

   if (g->mix_layer == NULL ) {   /* empty list */           
      g->mix_layer = m;
   }

   /* add depth info to this record */

   m->depth += depth;
   ++ m->n;

  /* compute average properties for mixed layer and add to appropriate fields.
     The average value is computed by summing the observed values weighted by
     the depth between observations.*/

   for (i = 0; i < nprops; ++i) {
      if (sta.observ[prop_indx[i]] != NULL) {
         j = 0;
         x = sta.observ[prop_indx[i]][0];
         n = 1;               /* weight the observation at sea surface */
         if (x < -8.9 ) {     /* unless it is a missing value */
            x = 0.0;
            n = 0;
         }
         dprev = 0.0;
         vprev = sta.observ[prop_indx[i]][0];
         while ( (j < hdr.nobs) && (de[j] <= depth)) {
            if ( (v = sta.observ[prop_indx[i]][j]) > -8.9) {
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
            m->prop[i] += (x / (float) n);
	    m->wghtsum[i] = m->wghtsum[i] + 1.0;   /* this counts number of profiles in this mixed layer */
         }
      }
   }
   return;

}  /* end insert_surf_vals() */

/*****************************************************************************/
struct surfrec *get_surfrec(struct surfrec *rptr, double d, int nprops)

   /* Recursively searches a linked list of surfrecs to:
        a) find an existing record for the specified density;
     or b) create a new record and insert it into the ordered list.

      Returns a pointer to the appropriate record.
      
  arguments: 
       rptr:    pointer to start of list 
          d:    density to key on 
     nprops:    # of properties to allocate space for
   */
{
    struct surfrec *r1ptr;
    double *tempd, *tempd2;

    if (rptr == NULL) {         /* empty list */
       r1ptr = create_surfrec(nprops);
       r1ptr->density = d;
       return(r1ptr);
    }

    if (NINT(d * 10) == NINT(rptr->density * 10)) {  /* current rec */
        return (rptr);
    }

    if (d < (rptr->density - .00001)) {   /* insert before the current rec */

       r1ptr = create_surfrec(nprops);
       tempd = r1ptr->prop;
       tempd2 = r1ptr->wghtsum;

         /* copy all fields from rptr into r1ptr */
       r1ptr->density = rptr->density;
       r1ptr->depth = rptr->depth;
       r1ptr->prop = rptr->prop;
       r1ptr->wghtsum = rptr->wghtsum;
       r1ptr->n = rptr->n;
       r1ptr->next = rptr->next;

        /* initialize the fields of rptr and link it to r1ptr */
       rptr->density = d;
       rptr->depth = 0;
       rptr->prop = tempd;
       rptr->wghtsum = tempd2;
       rptr->n = 0;
       rptr->next = r1ptr;
       
       return(rptr);
    }

    r1ptr = get_surfrec(rptr->next, d, nprops);  /* search rest of list */
    if (rptr->next == NULL)
          rptr->next = r1ptr;
    return (r1ptr);

}   /* end get_surfrec() */

/*****************************************************************************/
struct surfrec *create_surfrec(int nprops)

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
   r->density = 0;
   r->n = 0;
   r->next = NULL;
   r->wghtsum = (double *) malloc(nprops * sizeof(double));
   r->prop = (double *) malloc(nprops * sizeof(double));

   if (r->prop == NULL) {
      fprintf(stderr,"\nUnable to allocate memory in create_surfrec()\n");
      exit(1);
   }

   for (i = 0; i < nprops; ++i) {
      r->prop[i] = 0.0;
      r->wghtsum[i] = 0.0;
   }   

   return (r);
}   /* end create_surfrec() */

/*****************************************************************************/
void compute_avg(double *sum, double *weightsum, int nlevs) {

  /* Arguments:
   *       sum:    array containing summed, weighted values 
   * weightsum:    array containing sum of weights for each level
   *     nlevs:    # of elements in array 
   */

   int j;

   /* Loop over each level. */
   for (j = 0; j < nlevs; ++j) {
     /* Check that we have a non-zero weightsum
      * and if that is the case, divide the
      * weighted sum by the sum of the weights.
      * Otherwise, insert a fill value. */
     sum[j] = ( weightsum[j] > 0) ?  sum[j] / weightsum[j] : (double) HBEMPTY;
   }
   return;

}  /* End compute_avg() */

/*****************************************************************************/
void define_sigma_transitions(struct gridnode *g, int nprops)

 /*  determines where the averaged depth values cross the depth ranges (zmin
    and zmax) for each sigma reference level and zeros the counters for the 
    accumulated sums of any level outside that range. 
          g:    ptr to a grid node 
     nprops:    number of properties stored at each gridnode 
  */
{
   int i, j, refstart, refstop, iref;

   refstop = 0;
   for (iref = 0; iref < NREFLEVS; ++iref) {
      refstart = refstop;
      refstop = refstart + nsig[iref];

      for (i = refstart; i < refstop; ++i) {
         if (g->nobs[i] > 0) {
            if (g->d[i] < zmin[iref] || g->d[i] >= zmax[iref]) {
                g->nobs[i] = 0;
                for (j = 0; j < nprops; ++j) {
                   g->count[j][i] = 0;
                }
            }
         }
      }

   }
   return;

}  /* end define_sigma_transitions() */
/*****************************************************************************/
struct surfrec *define_avg_mixed_layer(struct surfrec *listptr, int nprops)

   /* Traverses a linked list of surfrecs sorted by density and computes 
      the average density at the sea surface, then searches the linked list
      for that density or 2 density records bracketing the average density and
      returns a pointer to a surfrec containing the property values
      associated with the average density. 
      
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
      x += r1->density ; 
      ++n;   
      r1 = r1->next;    
   }
   key = NINT((x / (double)n) * 10.);     /* average density * 10  */
   
   /* Now search list for the density value ... */

   r1 = get_density_rec(key, listptr, nprops);

   return (r1);

}  /* end define_avg_mixed_layer() */

/***********************************************************************/
struct surfrec *get_density_rec(int key, struct surfrec *rptr, int nprops)
    /* Recursively searches a linked list sorted by density for the record
       corresponding to key or for two records bracketing that density.
       Returns a record containing averaged property values for that density. 
       
          key:  density  being searched for *10  
         rptr:  ptr to element of linked list 
       nprops:  # of properties stored in each surfrec 
 
    */
{
   struct surfrec *r1, *r2, *new;
   double x[2], y[2];
   UI n;
   int  key1, key2, i;

   if (rptr == NULL) {           /* YIKES! */
       fprintf(stderr,"\nError in get_density_rec(): ");
       fprintf(stderr,"\nEnd of linked list encountered before finding avg density!");
       fprintf(stderr,"\nError in program logic....exiting.\n");
       exit(1);
   }

   r1 = rptr;
   r2 = rptr->next;

   if (r2 == NULL) {             /* end of list, so r1 must be the density */
      compute_avg(r1->prop, r1->wghtsum, nprops);
      r1->depth /= (double)r1->n;
      return (r1);
   }

   key2 = NINT(r2->density * 10);

   if (key > key2)                /* recursive part */
       return( get_density_rec(key, r2, nprops));
   

   if (key == key2) {             /* exact match! */
       compute_avg(r2->prop, r2->wghtsum, nprops);
       r2->depth /= (double)r2->n;
       return(r2);
   }
     
   key1 = NINT(r1->density * 10);

   if (key < key1) {              /* YIKES! */
       fprintf(stderr,"\nError in get_density_rec(): ");
       fprintf(stderr,"\nAvg density is less than first density in list!");
       fprintf(stderr,"\nThis is a bug in the program logic....exiting.\n");
       exit(1);
   }
   if (key == key1)  {            /* exact match! */
       compute_avg(r1->prop, r1->wghtsum, nprops);
       r1->depth /= (double)r1->n;
       return(r1);
   }

/* if we get this far, then the key must lie between r1 & r2... so 
      choose r1 or r2-- whichever is closest.   */


   if ((key - key1) < (key2 - key) ) {
       compute_avg(r1->prop, r1->wghtsum, nprops);
       r1->depth /= (double)r1->n;
       return(r1);
   }
   else {
       compute_avg(r2->prop, r2->wghtsum, nprops);
       r2->depth /= (double)r2->n;
       return(r2);
   }

/******************************************
*  This code is currently commented out!!! It was used to construct a new   
*  record with values interpolated between r1 and r2.... 
*        compute_avg(r1->prop, r1->wghtsum, nprops);  
*        compute_avg(r2->prop, r2->wghtsum, nprops);
*        x[0] = r1->density;
*        x[1] = r2->density;
 *  
*        new = create_surfrec(nprops);
*        new->density = (double) key / 10.0;
*        new->n = (r1->n + r2->n) >> 1;      
 *       y[0] = r1->depth / (double) r1->n;
*        y[1] = r2->depth / (double) r2->n;
*        new->depth = interpolate(new->density, x, y, 2);
*   
*        for (i = 0; i < nprops; ++i) {
*         if (! (r1->count[i] && r2->wghtsum[i])) {    
*             new->wghtsum[i] = 0;                    
*             new->prop[i] = (double)HBEMPTY;
*             continue;                         
*         }
*
*         new->wghtsum[i] = (r1->wghtsum[i] + r2->wghtsum[i]) >> 1;
*         y[0] = r1->prop[i];
*         y[1] = r2->prop[i];
*         new->prop[i] = interpolate(new->density, x, y, 2);
*        }
*        return(new);
**************/

}  /* end get_density_rec() */

/***********************************************************************/
int do_std_depth(struct gridnode *g, int npts, int *prop_id, int nprops, int ncols, int icol, float **dataptr, short **countptr, float *bottom)

/*  Takes the info for a gridnode and computes the value of each property 
    interpolated onto the std depth levels.  The arrays at dataptr and countptr
    are assumed to have dimensions data[nprops][nstdlevs*ncols].  This
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
      prop_id:   array of property identifiers  
       nprops:   # of properties in prop array  
        ncols:   # of cols in data and count arrays  
         icol:   column # associated with this gridpt  
      dataptr:   ptr to start of data arrays   
     countptr:   ptr to start of count arrays  
       bottom:   returned depth of deepest observation 

*/
{
   int i, ii, j, jj, k, i_te, i_sa, datagap, room, is_denser;
   double *x, *y, *w, *sigval, *sig, z, r;
   double tref, pref, sig_bml, sig_deepest;
   float  *d;
   short  *n;
   int size, start;
   int deepestlev, reflev;
   struct surfrec *m, *m2;
   struct deepestrec *btmptr, *b1;
   extern double std_depth[];

   sigval = (double *) malloc(npts * sizeof(double));
   sig = (double *) malloc(npts * sizeof(double));
   size = monotonic(g->d, g->nobs, sigval, g->prop, g->count, nprops, npts);

   if (size <= 1) {
      for (i = 0; i < nprops; ++i) {
           d = &dataptr[i][icol * NSTDLEVS];
           n = &countptr[i][icol * NSTDLEVS];
           for (j = 0; j < NSTDLEVS; ++j) {
               *(d++) = (float) HBEMPTY;
               *(n++) = 0;
           }
      }
      *bottom = (float) HBEMPTY;
      free((void *)sigval);
      return(NSTDLEVS);
   }

/* construct temporary arrays, x, y & w to store depth, y-property and 
   w-# of obs for that property.  Enlarge the size of the arrays 
   to accommodate the mixed layer and deepest observations.  */

   room = 50;
   x = (double *) malloc((size+room) * sizeof(double));
   y = (double *) malloc((size+room) * sizeof(double));
   w = (double *) malloc((size+room) * sizeof(double));
   if (w == NULL) {
     fprintf(stderr,"\nUnable to allocate memory in do_std_depth()\n");
     exit(1);
   }

   m = define_avg_mixed_layer(g->mix_layer, nprops);
   
/* identify index for temp and salt */

   i = i_te = i_sa = -1;
   while (++i < nprops) {
      if (prop_id[i] == (int)TE)
          i_te = i;
      if (prop_id[i] == (int) SA)
          i_sa = i;
   }

/* find the deepest sigma_level for which there is data */

   deepestlev = 0;
   for (k = 0; k < size; ++k) {
      if (g->nobs[k] > 0)
           deepestlev = k;
   }


/* Define the starting point in the sigma series which is heavier than
    the density at the bottom of the average mixed layer. Ensure that 
    the depth of this sigma level is deeper than the depth of the mixed 
    layer. */
    
   start = 0;
   if (m != NULL) {
   
      if (m->depth > g->d[deepestlev])
          m->depth = g->d[deepestlev];

      pref = 0.0;
      if (m->depth >= zmin[4])
         pref = 4000;
      else if (m->depth >= zmin[3])
         pref = 3000;
      else if (m->depth >= zmin[2])
         pref = 2000;
      else if (m->depth >= zmin[1]) 
         pref = 1000;
	 
      tref = hb_theta(m->prop[i_sa], m->prop[i_te], m->depth, pref);
      hb_svan(m->prop[i_sa], tref, pref, &sig_bml);  /* sigma at bottom of mixed layer */
       
       while ((sigval[start] < sig_bml) && (start < size))
           ++start;

      while ((m->depth > g->d[start]) && (start < size))
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
            if (m->wghtsum[j] > 0) {   
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

           if (z > (HBEMPTY+10.0)) {                  /* check for vertical datagaps */
              jj = 0;
              while (x[++jj] < std_depth[i]) 
                   ;

              if (((int)x[jj-1] == (int)std_depth[i]) 
                  || ((int)x[jj] == (int)std_depth[i])) 
                  datagap = 0;

              else if (std_depth[i] < GAP_CHANGE_DEPTH) 
                  datagap = (x[jj] - x[jj-1]) > GAP_SHALLOW;

              else
                  datagap = (x[jj] - x[jj-1]) > GAP_DEEP;

              if (datagap && (ABS(sig[jj] - sig[jj-1]) >  0.04)) {      /* check for pycnostad */
                  z = (double) HBEMPTY;
                  if (j == 0 && report_gaps)
                    fprintf(stderr," datagap:  %.1lf  %.1lf\n", x[jj-1], x[jj]);
              }
              else {
                   r = interpolate(std_depth[i], x, w, npts);
                   if (r < 0)
                       k = 0;  /* this should never happen */
                   else {
                       k = (r - (int)r) > 0 ? (int) r + 1 : (int) r;
                   }
              }
	     
           }
            if ((z > (HBEMPTY+10.0)) && (k == 0))
	         k = 1;      /* non-empty node gets at least 1 observation */      
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
   delete_surfrec_list(g->mix_layer);
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
         w[i] = m->wghtsum[iprop];
         x[i] = i * 100.0;           
   }         
   w[npts-1] = m->wghtsum[iprop];
   x[npts-1] = m->depth;

/*  !**! Special cases for individual properties... */

   switch ((enum property) prop_id) {
       case PR :                   
                     for (i = 0; i < npts; ++i) {
                        y[i] = x[i];   
                     }         
                     break;
       case HT:
       case PE:
                     if (i_te < 0 || i_sa < 0) {
                        for (i = 0; i < npts; ++i)
                          y[i] = 0.0;
                        break;
                     }
                     t = (double *) malloc(npts * sizeof(double));
                     s = (double *) malloc(npts * sizeof(double));
                     for (i = 0; i < npts; ++i) {
                        t[i] = m->prop[i_te];
                        s[i] = m->prop[i_sa];
                     }
                     if (prop_id == (int)HT) 
                       compute_height(npts, x, t, s, ht_pref, y);
                     else
                       compute_energy(npts, x, t, s, pe_pref, y);  
                     free(t);
                     free(s);
                     break;
       case S0:
                     if (i_te < 0 || i_sa < 0) {
                        for (i = 0; i < npts; ++i)
                          y[i] = m->prop[iprop];
                        break;
                     }

                        t = (double *) malloc(npts * sizeof(double));
                        s = (double *) malloc(npts * sizeof(double));
                        for (i = 0; i < npts; ++i) {
                           t[i] = m->prop[i_te];
                           s[i] = m->prop[i_sa];
                        }
                        pref = 0.0;
                        compute_sigma(pref, npts, y, x, t, s);
                        free(t);
                        free(s);
                        break;
       case S1:
                     if (i_te < 0 || i_sa < 0) {
                        for (i = 0; i < npts; ++i)
                          y[i] = m->prop[iprop];
                        break;
                     }
                        t = (double *) malloc(npts * sizeof(double));
                        s = (double *) malloc(npts * sizeof(double));
                        for (i = 0; i < npts; ++i) {
                           t[i] = m->prop[i_te];
                           s[i] = m->prop[i_sa];
                        }
                        pref = 1000.0;
                        compute_sigma(pref, npts, y, x, t, s);
                        free(t);
                        free(s);
                        break;
       case S2:
                     if (i_te < 0 || i_sa < 0) {
                        for (i = 0; i < npts; ++i)
                          y[i] = m->prop[iprop];
                        break;
                     }
                        t = (double *) malloc(npts * sizeof(double));
                        s = (double *) malloc(npts * sizeof(double));
                        for (i = 0; i < npts; ++i) {
                           t[i] = m->prop[i_te];
                           s[i] = m->prop[i_sa];
                        }
                        pref = 2000.0;
                        compute_sigma(pref, npts, y, x, t, s);
                        free(t);
                        free(s);
                        break;
       case S3:
                     if (i_te < 0 || i_sa < 0) {
                        for (i = 0; i < npts; ++i)
                          y[i] = m->prop[iprop];
                        break;
                     }
                        t = (double *) malloc(npts * sizeof(double));
                        s = (double *) malloc(npts * sizeof(double));
                        for (i = 0; i < npts; ++i) {
                           t[i] = m->prop[i_te];
                           s[i] = m->prop[i_sa];
                        }
                        pref = 3000.0;
                        compute_sigma(pref, npts, y, x, t, s);
                        free(t);
                        free(s);
                        break;
       case S4:
                     if (i_te < 0 || i_sa < 0) {
                        for (i = 0; i < npts; ++i)
                          y[i] = m->prop[iprop];
                        break;
                     }
                        t = (double *) malloc(npts * sizeof(double));
                        s = (double *) malloc(npts * sizeof(double));
                        for (i = 0; i < npts; ++i) {
                           t[i] = m->prop[i_te];
                           s[i] = m->prop[i_sa];
                        }
                        pref = 4000.0;
                        compute_sigma(pref, npts, y, x, t, s);
                        free(t);
                        free(s);
                        break;
       case S_:
                     if (i_te < 0 || i_sa < 0) {
                        for (i = 0; i < npts; ++i)
                          y[i] = m->prop[iprop];
                        break;
                     }
                        t = (double *) malloc(npts * sizeof(double));
                        s = (double *) malloc(npts * sizeof(double));
                        for (i = 0; i < npts; ++i) {
                           t[i] = m->prop[i_te];
                           s[i] = m->prop[i_sa];
                        }
                        compute_sigma(s_pref, npts, y, x, t, s);
                        free(t);
                        free(s);
                        break;
       case VA:
                     if (i_te < 0 || i_sa < 0) {
                        for (i = 0; i < npts; ++i)
                          y[i] = m->prop[iprop];
                        break;
                     }
                        t = (double *) malloc(npts * sizeof(double));
                        s = (double *) malloc(npts * sizeof(double));
                        for (i = 0; i < npts; ++i) {
                           t[i] = m->prop[i_te];
                           s[i] = m->prop[i_sa];
                        }
                        compute_svan(npts, y, x, t, s);
                        free(t);
                        free(s);
                       break;
       case SV:
                     if (i_te < 0 || i_sa < 0) {
                        for (i = 0; i < npts; ++i)
                          y[i] = m->prop[iprop];
                        break;
                     }
                        t = (double *) malloc(npts * sizeof(double));
                        s = (double *) malloc(npts * sizeof(double));
                        for (i = 0; i < npts; ++i) {
                           t[i] = m->prop[i_te];
                           s[i] = m->prop[i_sa];
                        }
                        compute_sp_vol(npts, y, x, t, s);
                        free(t);
                        free(s);
                       break;
       default:
                        for (i = 0; i < npts; ++i)
                          y[i] = m->prop[iprop];
                        break;
   }  /* end switch */

   return (npts);

}  /* end mixed_layer_vals() */
/****************************************************************************/
int monotonic(double *d, UI  *nobs, double *sigmas, double **xx, UI  **count, int nprops, int npts)

/* sorts the arrays into increasing order by depth, removing levels with no obs and 
   ensuring there are no depth inversions.  Returns the 
   # of depth levels containing data. 
   
           d    starting addr of depth  
        nobs    starting addr of nobs array  
      sigmas    array to put sigma values 
          xx    starting addr of other properties  
       count    starting addr of nobs arrays for other properties  
      nprops    # of rows (other properties) in xx  
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

double interpolate(double xval, double *x, double *y, int nypts) {

  /* Performs a linear interpolation to find
   * the position of xval in array, x, and
   * returns the corresponding value in the
   * array, y.  If xval does not appear in
   * array x, the value HBEMPTY is returned.
   * This routine assumes that the x array
   * is monotonic and continuous (no missing
   * values); and it assumes the y array is
   * continuous. */

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
    
    /* Perform the linear interpolation. */
    return ( y[k] + (y[k+1] - y[k]) * v1 / (x[k+1] - x[k]) );

  } /* End of loop over all points in base array. */

  /* If we did not find a suitable place
   * perform the interpolation, then we
   * return a fill value. */
  return ((double)HBEMPTY);

}   /* End interpolate() */

/****************************************************************************/

void sum_levels(double *y, double *d, int ny, double *ysum, UI *count) {

  /* Interpolates the y and d arrays with
   * ny observation levels onto the stddepth
   * levels (globally defined).  Checks for
   * vertical data gaps and does not interpolate
   * over gaps which exceed 200 m in the
   * thermocline (upper 1000 m) or 1000 m 
   * elsewhere. The interpolated values are
   * added to the appropriate levels of the
   * array defined by ysum and the count array
   * is incremented. */ 
  
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
      while (j < n && dtmp[j] < std_depth[i]) 
	++j ;
      
      datagap = 0;
      if ((dtmp[j] == std_depth[i]) )
	datagap = 0;
      else if ( j > 0) {
	
	if  (dtmp[j-1] == std_depth[i]) 
	  datagap = 0;
	else if ( std_depth[i] < GAP_CHANGE_DEPTH)
	  datagap = (dtmp[j] - dtmp[j-1]) > GAP_SHALLOW;
	else
	  datagap = (dtmp[j] - dtmp[j-1]) > GAP_DEEP;
      }              
      
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

/***************************************************************************/
