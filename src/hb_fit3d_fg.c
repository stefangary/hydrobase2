/*  hb_fit3d_fg.c

................................................................................
                          *******  HydroBase 2 *******
                             Ruth Curry
                             Woods Hole Oceanographic Institution
                             Dec 2000
...................................................
*
*  Interpolates for missing values in 3-dimensional grids of 
*  hydrographic properties by propagating observed anomalies relative to an initial background
*  field along isopycnal surfaces.
*  Reads lat/lon/depth gridded fields representing an initial field of properties in HydroBase cdf format.
*  Reads a second cdf file containing the gridded fields that are to be interpolated 
*  and a binary topography file.   Points that are missing (determined from 
*  the topography data) are filled in by first estimating an isopycnal value for the missing
* point from the initial
*  field and surrounding observed points, then using property difference between the
* initial field and the observed field for surrounding points to determine by interpolation of nearest points
*  the anomaly for the gridpoint in question. The anomaly is then applied to initial field value to
*  fill in the missing node with a new value.
*   Masked nodes are not incorporated into the fit.
*  A x-/y- Search radius parameter determines the maximum acceptable distance an
*  interpolated value can be from an observed value. If a masked node is
*  encountered while searching along radius, the search is suspended in 
*  that direction.  If no observation is within that distance, the node is 
*  flagged as "empty". Depths below the seafloor depth are flagged as "masked".
*  
*  Output is a HydroBase cdf file.  

*/

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include "hydrobase.h"
#include "hydro_cdf.h"
#include "hb_fit.h"
#include "hb_memory.h"
#include "hb_paths.h"

/* input_file pathnames */

#define    EXTENT   ""
#define    DIR      ""
#define    FG_EXTENT   ""
#define    FG_DIR      ""

#define   RADperDEG 0.017453292             /* pi/180 */
#define  EarthRadius  6371.0    /* in km */

#define    PREFILL      0        /* turn off prefilling of cdf variables */

/* global variables */

float HB_f_mask;          /* float values are used since the output */
float HB_f_empty;         /* cdf values are of this type */
 

/* needed by print_usage() */
    
double xradius, yradius;
double lengthscale;

/******  define local data structures  **********/

struct MIXLAYER {
    double density;
    double depth;
    int nobs;
};

/* prototypes for locally defined functions */	

void print_usage(char *);
int parse_p_option(char *, int *);
int cdf_construct(struct GRID_INFO *, int, int *, char *, int, char **);
void get_seafloor(FILE *, double, double, short *, struct GRID_INFO *);
void get_x_profile(float, float, int, int, char **, char *, char *, double *, short *, float);
int define_bottom(float *, short);
void define_mixed_layer( double *, double **,  short *, double *, struct GRID_INFO *, int, int, float, struct MIXLAYER *);
void get_sig_profiles(float, float, float *, short, int, char **, char *, char *, double *, double *, double *, short *, short *, short *, double *, double *, double *, double *, double *);
double find_bottom_isopycnal(int, int, int, double **, float *, float, double);
double get_neighbors_anom(double, double,   double *, int , struct GRID_INFO *, double **, double **,  float *, double **, double **, float *);
void get_neighbors( double, double, double *,  int , double **, double **,  float *);
void get_weights(double *, float *, double, struct GRID_INFO *);
double get_xval(double *, double *, double, float, double);
void x_to_stddepth(double *, short *, float *, double *, double *, float, int );


main(int argc, char **argv)
{
  int i, j, n, nprops;
  int nfiles, nfiles_fg, start_list_fg;
  int nz, nsqout, nwsq1, nwsq2;
  int row, col, sq, sq_out, error, tbin;
  int wsq1, wsq2, wsq_mid;
  int wrow, wcol, wcol1, wrow1;
  int do_fit, gotprops, is_pr, is_te, is_sa;
  int aflag, bflag, iflag, oflag, uflag, pflag;
  int infile, outfile, print_msg; 
  int xrad, yrad; 
  int pixel_grid;           /* pixel or gridnode registration */
  int n_empty, n_mask, n_set;
  int xcount;
  int mask_it;
  int *prop_indx, *prop_req, *prop_avail;
  float wlat, wlon, lat, lon, wlat2, wlon2;
  double xoffset, yoffset;
  double topo_xinc, topo_yinc;
  double isopyc_above, curr_depth;
  char modifier;
  char *st;
  char *outfile_name;
  char *dir, *extent;
  char *dir_fg, *extent_fg;
  FILE	*topofile;
  struct CDF_HDR *cdfin;
  struct GRID_INFO hwork1, hout, hwork2;
  short *seafloor;
  double *xsurf, *weights;
  double xbase, xanom, xwghtavg, xstddev;      
  short *new_count;
  short **tcount1, **scount1, **pcount1, **xcount1;
  short **tcount1_fg, **scount1_fg, **pcount1_fg, **xcount1_fg;
  float *bottom_de, *bdwork2, *bdwork1; 
  float *bdwork1_fg, *bdwork2_fg;
  float  *xout, *xtmp_f;
  double *xfit, *xtmp_d;
  double **xwork1, **pwork1, **swork1, **twork1;
  double **xwork1_fg, **pwork1_fg, **swork1_fg, **twork1_fg;
  double **sig0work1, **sig1work1,  **sig2work1, **sig3work1, **sig4work1;
  double  **sig0work1_fg, **sig1work1_fg, **sig2work1_fg, **sig3work1_fg, **sig4work1_fg;
  double **sigptr, **sigptr_fg, **xwork2, **xwork2_fg;
  double **sig0work2, **sig1work2, **sig2work2, **sig3work2, **sig4work2;
  double **sig0work2_fg, **sig1work2_fg, **sig2work2_fg, **sig3work2_fg, **sig4work2_fg;
  double **isopyc, **sfit;
  double **pfit, **tfit, **dfit;
  struct MIXLAYER **mlptr;

/* Set these default values. */

  HB_f_mask = (float) HBMASK;      /* HydroBase cdf file mask flag */
  HB_f_empty = (float) HBEMPTY;    /*  and empty flag */
    
  error = 0;
  aflag = bflag = iflag = oflag = uflag = pflag =  0;
  xradius = 2.0;          /* default search radius in degrees */
  yradius = 2.0;
  pixel_grid = 1;          /* default pixel registration */
  topofile = (FILE *) NULL; 
  topo_xinc = topo_yinc = 0.1;
  dir = DIR;
  extent = EXTENT;
  dir_fg = FG_DIR;
  extent_fg = FG_EXTENT;
  nfiles = nfiles_fg = 0;
  start_list_fg = 0;
  print_msg = 1;
  tbin = 0;
  lengthscale = 100.0;  /* for distance weighting */
  
/*----------------------------------------*/  
  
/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
/*----------------------------------------*/  
/* parse command line arguments */

  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
         case 'A':  /* list of first-guess files follows. */
             aflag = 1;
	     start_list_fg = i; 
	     break;
        case 'B':  /* get grid bounds */
	  bflag = 1;
          st = &argv[i][2];
          if (*st == '/')
             ++st;
          error = (sscanf(st,"%lf", &hout.x_min) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%lf", &hout.x_max) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%lf", &hout.y_min) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%lf", &hout.y_max) != 1);
                        	     
	  if (&hout.x_min > &hout.x_max) {
	    fprintf(stderr,"\nWest bound must be numerically <= east bound");
	    exit(1);
	  }
	  
	  if (&hout.y_min > &hout.y_max) {
	    fprintf(stderr,"\nNorth bound must be numerically <= south bound");
	    exit(1);
	  }
          break;
	  
         case 'C':                   /* input dir for first-guess files*/
          dir_fg = &argv[i][2];
          break;
        case 'D':                   /* get input dir */
          dir = &argv[i][2];
          break;
        case 'E':                    /* get file extent */
          extent = &argv[i][2];
          break;
        case 'F':                    /* extent for first-guess files */
          extent_fg = &argv[i][2];
          break;

  
        case 'G':        /* set to gridnode registration */
	  pixel_grid = 0;
	  break;
	  
        case 'I':
          iflag = 1;
          error = (sscanf(&argv[i][2],"%lf", &hout.x_inc) == 1) ? 0 : 1;
          hout.y_inc = hout.x_inc;
          st = &argv[i][2];
          while (*(st++) != '\0') {
            if (*st == '/') {
               ++st;
               error = (sscanf(st,"%lf", &hout.y_inc) == 1) ? 0 : 1;
               break;
            }
          }
                        
          break;
        case 'L':
          error = (sscanf(&argv[i][2],"%lf", &lengthscale) == 1) ? 0 : 1;
          break;
	  
        case 'O':
	  oflag = 1;
          outfile_name = &argv[i][2];
          break;
	  	      
        case 'P':
          pflag = 1;
	  prop_req = (int *) calloc(MAXPROP, sizeof(int));
          nprops = parse_p_option(&argv[i][2], prop_req);
          break;


        case 'S':
          error = sscanf(&argv[i][2],"%lf", &xradius) != 1;
          yradius = xradius;
          st = &argv[i][2];
          while (*(st++) != '\0') {
            if (*st == '/') {
               ++st;
               error = (sscanf(st,"%lf", &yradius) == 1) ? 0 : 1;
               break;
            }
          }
          break;
	  
	case 'h':  
	   print_usage(argv[0]);
	   exit(0);
	   
        default:
          error = TRUE;
 
      } /* end switch */
      
       
      if (error ) {
         fprintf(stderr,"\nError parsing command line args.\n");
         fprintf(stderr,"     in particular: '%s'\n", argv[i]);
         exit(1);
      }
    }  /* end if */
    
    else {
      if (aflag)
         ++nfiles_fg;
      else
         ++nfiles;
    }  /* end else */
    
  }  /* end for */
  
  
/*--------------------------------------------*/    
/*  Check syntax of options */ 

   
   error = 0;
   
    if (!nfiles ) {
       fprintf(stderr,"\nYou must specify input cdf_files as first argument(s).");
       ++error;    
    }
    if (!bflag ) {
       fprintf(stderr,"\nYou must specify bounds with -B<w/e/s/n> ");
       ++error;    
    }
    if (!iflag ) {
       fprintf(stderr,"\nYou must specify grid spacing with -I<xincr>[/yincr] ");
      ++error;
    }
    if (!oflag ) {
       fprintf(stderr,"\nYou must specify -O<output_cdf_file> ");
      ++error;
    }
                                 
    if (!pflag ) {
       fprintf(stderr,"\nYou must specify a list of properties for output: -Ppr/te/sa/ox ");
      ++error;
    }
    
    if (!aflag) {
       nfiles_fg = nfiles;
       fprintf(stderr,"\nUsing list of input root_file_names for first-guess with appended  extent: %s ", extent_fg);
    }
    
    if (dir_fg == "") {
        dir_fg = dir;
    }
    
    if (topofile == NULL ) {
       topofile = fopen(BATHPATH,"r");
       if (topofile == NULL) {
          fprintf(stderr,"\nUnable to open default topography file\n", BATHPATH);
          exit(1);
      }
          fprintf(stderr,"Opened %s ", BATHPATH); 
    }
    if (!uflag ) {
      fprintf(stderr,"Increment for topofile: [%.2lf/%.2lf] \n", topo_xinc, topo_yinc);
    }
   if (xradius < 0. || yradius < 0.) {
      fprintf (stderr, "SYNTAX ERROR -S option.  Must specify a positive search radius.\n");
      error++;
   }
    
   if (error) {
    fprintf(stderr,"\nUse -h for complete usage info. \n");
    exit(1);
   }
 

/*-------------------------------------------*/   
/* Get array of standard depths from first cdf file in list */
   
   i = 1;
   infile = -1; 
   while (infile < 0 && i <= nfiles) {
     infile = cdf_open(dir, argv[i], extent, FALSE);
     ++i;
   }
   if (infile < 0) {
      fprintf(stderr,"Unable to open any cdf files.");
      fprintf(stderr,"Directory: %s \nRoot: %s \n Extent: %s", dir, argv[i-1], extent);
      exit(1);
    }
       
   cdfin = (struct CDF_HDR *) calloc(1, sizeof(struct CDF_HDR));
   
   if (error = read_cdf_hdr(infile, cdfin)) 
         exit (1);
   
   xtmp_f = (float *) calloc((size_t)cdfin->nz, sizeof(float));
   NSTDLEVS = read_cdf_depths(infile, xtmp_f);
   for (i = 0; i < NSTDLEVS; ++i) 
       std_depth[i] = (double) xtmp_f[i];
       
   std_depth_initialized = 1;
   free(xtmp_f);
   cdf_close(infile);
   nz = NSTDLEVS - 1;   /* index to bottom depth */   
   
/*-------------------------------------------*/   
/* Check for requested properties availability.    After this, prop_indx[] and
    nprops will hold the property info for the output cdf file. */
   
   prop_indx = (int *) calloc((size_t)nprops, sizeof(int)); 
   prop_avail = (int *) calloc((size_t)MAXPROP, sizeof(int)); 
   for (i = 0; i < cdfin->nprops; ++i)    
      prop_avail[get_prop_indx(cdfin->prop_id[i])] = 1;
      
   n = 0;   
   for (i = 0; i < nprops; ++i) {
      if ( prop_avail[prop_req[i]] ) {
        prop_indx[n++] = prop_req[i];
      }
      else {
        if (prop_req[i] == (int) PR || prop_req[i] == (int) TE || prop_req[i] == (int) SA)
          fprintf(stderr,"\n FATAL ERROR!! Property %.2s not available in cdf file (pr te sa are mandatory).\n", get_prop_mne(prop_req[i]));       
	else
          fprintf(stderr,"\n WARNING!! Property %.2s not available in cdf file and will not be output.\n", get_prop_mne(prop_req[i]));       
      }
   }
   nprops = n;
   free(prop_avail);    
   free(prop_req);
  
/*-------------------------------------------*/   
/*   set up CDF_HDR and output file */

   hout.nx = (int) NINT((hout.x_max - hout.x_min)/hout.x_inc);
   hout.ny = (int) NINT((hout.y_max - hout.y_min)/hout.y_inc);
   
   xoffset = 0.5*hout.x_inc;
   yoffset = 0.5*hout.y_inc;

   hout.node_offset = pixel_grid;   
   if (!pixel_grid) {
      xoffset = 0.0;
      yoffset = 0.0;
      ++hout.nx;
      ++hout.ny;
   }
   nsqout = hout.nx * hout.ny;

   outfile = cdf_construct(&hout, nprops, prop_indx, outfile_name, argc, argv);
/*-------------------------------------------*/   
/* Adjust order of properties so that sa is fitted before te -- since 
   salinity is needed to convert potential temperature back to in situ temperature. 
*/
   prop_indx[1] = (int)SA;
   prop_indx[2] = (int)TE;
 
/*-------------------------------------------*/   
/*  Adjust xmin/ymin to be lower/left gridnode
    and xmax/ymax to upper/right gridnode depending upon the value of 
    node_offset.   
*/
    
   hout.x_min += xoffset;                         /* west gridnode */
   hout.x_max = hout.x_min + (hout.nx -1) * hout.x_inc;   /* east  */
   hout.y_min = hout.y_min + yoffset;            /* south gridnode */
   hout.y_max = hout.y_min + (hout.ny -1) * hout.y_inc;  /* north  */
   hout.node_offset = 0; /* now normalized to gridnode registration */   

/* save these values in gridnode units */

   xrad = NINT(xradius / hout.x_inc);  
   yrad = NINT(yradius / hout.y_inc);
           
  if (pixel_grid)
     fprintf (stderr, "\nUsing %s registration", "pixel");
  else
     fprintf (stderr, "\nUsing %s registration", "gridnode");
   fprintf (stderr, "\nGrid dimensions are nx = %d, ny = %d", hout.nx, hout.ny);
   fprintf (stderr, "\nOutput xmin/ymin gridnode: %8.3lf/%8.3lf ", hout.x_min, hout.y_min);
   
      
   fprintf (stderr, "\nNumber of gridnodes in search radii: %d/%d", xrad, yrad);
   xradius = xrad * hout.x_inc;
   yradius = yrad * hout.y_inc;
   fprintf (stderr, "\n Search radius must be a multiple of grid increment -- using: %.3lf/%.3lf ", xradius, yradius);
/*-------------------------------------------*/   
/* Set up big work grid -- output grid + extended bounds  */

   hwork1.x_min = hout.x_min - xradius;
   hwork1.x_max = hout.x_max + xradius;
   hwork1.y_min = hout.y_min - yradius;
   hwork1.y_max = hout.y_max + yradius;
   hwork1.x_inc = hout.x_inc;
   hwork1.y_inc = hout.y_inc;
   hwork1.node_offset = 0;     /* normalized to gridnode registration */
   
   hwork1.nx = hout.nx + 2 * xrad;
   hwork1.ny = hout.ny + 2 * yrad;
   nwsq1 =  hwork1.nx * hwork1.ny;  


   xwork1 = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   pwork1 = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   twork1 = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   swork1 = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   sig0work1 = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   sig1work1 = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   sig2work1 = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   sig3work1 = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   sig4work1 = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   pcount1 = (short **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(short *));
   tcount1 = (short **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(short *));
   scount1 = (short **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(short *));
   xcount1 = (short **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(short *));

   xwork1_fg = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   pwork1_fg = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   twork1_fg = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   swork1_fg = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   sig0work1_fg = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   sig1work1_fg = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   sig2work1_fg = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   sig3work1_fg = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   sig4work1_fg = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   pcount1_fg = (short **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(short *));
   tcount1_fg = (short **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(short *));
   scount1_fg = (short **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(short *));
   xcount1_fg = (short **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(short *));
  

   for (sq = 0; sq < nwsq1; ++sq) {
     pwork1[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     twork1[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     swork1[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     sig0work1[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     sig1work1[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     sig2work1[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     sig3work1[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     sig4work1[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     pcount1[sq] = (short *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(short));
     tcount1[sq] = (short *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(short));
     scount1[sq] = (short *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(short));
 
     pwork1_fg[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     twork1_fg[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     swork1_fg[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     sig0work1_fg[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     sig1work1_fg[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     sig2work1_fg[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     sig3work1_fg[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     sig4work1_fg[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     pcount1_fg[sq] = (short *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(short));
     tcount1_fg[sq] = (short *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(short));
     scount1_fg[sq] = (short *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(short));
   }

   sfit = (double **) get_memory((void *)NULL, (size_t)nsqout, sizeof(double *));
   tfit = (double **) get_memory((void *)NULL, (size_t)nsqout, sizeof(double *));
   pfit = (double **) get_memory((void *)NULL, (size_t)nsqout, sizeof(double *));
   dfit = (double **) get_memory((void *)NULL, (size_t)nsqout, sizeof(double *));
   isopyc = (double **) get_memory((void *)NULL, (size_t)nsqout, sizeof(double *));
   mlptr = (struct MIXLAYER **) get_memory((void *)NULL, (size_t)nsqout, sizeof(struct MIXLAYER *));
   for (sq = 0; sq < nsqout; ++sq) {
      isopyc[sq] = (double *) get_memory((void *)NULL, (size_t)NSTDLEVS, sizeof(double));
      mlptr[sq] = (struct MIXLAYER *) get_memory((void *)NULL, (size_t)1, sizeof(struct MIXLAYER ) );
   }   
   xfit = (double *) get_memory((void *)NULL, (size_t)NSTDLEVS, sizeof(double));


/*-------------------------------------------*/   
/* Set up small work grid to store isopycnal surface sent to interp2d */

   hwork2.nx = xrad * 2 + 1;
   hwork2.ny = yrad * 2 + 1;
   hwork2.x_inc = hout.x_inc;
   hwork2.y_inc = hout.y_inc;
   hwork2.node_offset = 0;  /* normalized to gridnode registration */
   nwsq2 = hwork2.nx * hwork2.ny;
   wsq_mid  = yrad * hwork2.nx + xrad;  /* by definition the middle square of work2 */
   
   xwork2 = (double **) get_memory((void *)NULL,(size_t)nwsq2, sizeof(double *));
   sig0work2 = (double **) get_memory((void *)NULL,(size_t)nwsq2, sizeof(double *));
   sig1work2 = (double **) get_memory((void *)NULL,(size_t)nwsq2, sizeof(double *));
   sig2work2 = (double **) get_memory((void *)NULL,(size_t)nwsq2, sizeof(double *));
   sig3work2 = (double **) get_memory((void *)NULL,(size_t)nwsq2, sizeof(double *));
   sig4work2 = (double **) get_memory((void *)NULL,(size_t)nwsq2, sizeof(double *));

   xwork2_fg = (double **) get_memory((void *)NULL,(size_t)nwsq2, sizeof(double *));
   sig0work2_fg = (double **) get_memory((void *)NULL,(size_t)nwsq2, sizeof(double *));
   sig1work2_fg = (double **) get_memory((void *)NULL,(size_t)nwsq2, sizeof(double *));
   sig2work2_fg = (double **) get_memory((void *)NULL,(size_t)nwsq2, sizeof(double *));
   sig3work2_fg = (double **) get_memory((void *)NULL,(size_t)nwsq2, sizeof(double *));
   sig4work2_fg = (double **) get_memory((void *)NULL,(size_t)nwsq2, sizeof(double *));

   weights = (double *) get_memory((void *) NULL, (size_t)nwsq2, sizeof(double));
 /*-------------------------------------------*/   
/*-------------------------------------------*/   
  /* Read in topography values */
  
   seafloor = (short *) get_memory((void *)NULL, (size_t) nwsq1, sizeof(short));
   get_seafloor(topofile, topo_xinc, topo_yinc, seafloor, &hwork1);
/*-------------------------------------------*/ 

  /* set up array to store bottom depths for output and work */
  
   bottom_de = (float *) get_memory((void *)NULL, (size_t) nsqout, sizeof(float));
   bdwork1 = (float *) get_memory((void *)NULL,(size_t)nwsq1, sizeof(float));
   bdwork2 = (float *) get_memory((void *)NULL,(size_t)nwsq2, sizeof(float));
   bdwork1_fg = (float *) get_memory((void *)NULL,(size_t)nwsq1, sizeof(float));
   bdwork2_fg = (float *) get_memory((void *)NULL,(size_t)nwsq2, sizeof(float));
/*-------------------------------------------*/   
/* Visit each gridpoint in big workspace, read in profiles, first-guess profiles, and compute sigmas. */

   fprintf(stderr,"\nReading in profiles and first-guess fields...");
   
   for (row = 0; row < hwork1.ny; ++row ) {
      lat = hwork1.y_min + row * hwork1.y_inc;
      
      for (col = 0; col < hwork1.nx; ++col) {
         lon = hwork1.x_min + col * hwork1.x_inc;
  
	 sq = row * hwork1.nx + col;
	 
	 get_sig_profiles(lat, lon, &bdwork1_fg[sq], seafloor[sq], nfiles_fg, &argv[start_list_fg], dir_fg, extent_fg, pwork1_fg[sq], twork1_fg[sq], swork1_fg[sq], pcount1_fg[sq], tcount1_fg[sq], scount1_fg[sq], sig0work1_fg[sq], sig1work1_fg[sq], sig2work1_fg[sq], sig3work1_fg[sq], sig4work1_fg[sq]);
	 
	 /* Reset seafloor value if first guess profile is deeper */
	 if (bdwork1_fg[sq] > seafloor[sq])
	    seafloor[sq] = bdwork1_fg[sq];
	    
	 get_sig_profiles(lat, lon, &bdwork1[sq], seafloor[sq], nfiles, argv, dir, extent, pwork1[sq], twork1[sq], swork1[sq], pcount1[sq], tcount1[sq], scount1[sq], sig0work1[sq], sig1work1[sq], sig2work1[sq], sig3work1[sq], sig4work1[sq]);
	   
	  /* Convert in situ temperatures to potential temperatures for fitting */

         for (j = 0; j < NSTDLEVS; ++j) {
           if (twork1[sq][j] > -3.0 && twork1[sq][j] < 100.)
	      twork1[sq][j] = hb_theta(swork1[sq][j], twork1[sq][j], pwork1[sq][j], 0.0);
	      
           if (twork1_fg[sq][j] > -3.0 && twork1_fg[sq][j] < 100.)
	      twork1_fg[sq][j] = hb_theta(swork1_fg[sq][j], twork1_fg[sq][j], pwork1_fg[sq][j], 0.0);
	 }
	 
      } /* end for col */
   } /* end for row */
   
   free((void *) seafloor);

   fprintf(stderr,"\nNow fitting data...  ");


/* Loop for each property including pressure 
   The order of props have been pre-set to 1)pr 2)sa  3)te
   Pressure must be first and the output pressure profiles must be
   converted to depth values and saved in memory
   in order to interpolate all properties back onto std-depths */

   for (i = 0; i < nprops; ++i) {

       is_pr = (prop_indx[i] == (int)PR);
       is_sa = (prop_indx[i] == (int)SA);
       is_te = (prop_indx[i] == (int)TE);
       
       fprintf(stderr,"\n%2s ", get_prop_mne(prop_indx[i]));
      
       switch ((enum property) prop_indx[i]) {
       
 	  case PR :          /* just set the pointers */
	     xwork1 = pwork1;
	     xcount1 = pcount1;
	     xwork1_fg = pwork1_fg;
	     xcount1_fg = pcount1_fg;
	     break;
 	  case TE :          /* just set the pointers */
	     xwork1 = twork1;
	     xcount1 = tcount1;
	     xwork1_fg = twork1_fg;
	     xcount1_fg = tcount1_fg;
	     break;
	  case SA :
	     xwork1 = swork1;
	     xcount1 = scount1;
	     xwork1_fg = swork1_fg;
	     xcount1_fg = scount1_fg;
	     break;
	  
	  default:
	  
	   /* Read property from input files */
	   
	   for ( sq = 0; sq < nwsq1; ++sq) {
	      wrow = sq / hwork1.nx;
	      wcol = sq - wrow * hwork1.nx;
	      wlat = hwork1.y_min + wrow * hwork1.y_inc;
	      wlon = hwork1.x_min + wcol * hwork1.x_inc;
              xwork1[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
              xcount1[sq] = (short *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(short));
              xwork1_fg[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
              xcount1_fg[sq] = (short *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(short));
	      get_x_profile(wlat, wlon, prop_indx[i], nfiles, argv, dir, extent, xwork1[sq], xcount1[sq], bdwork1[sq]);
	      get_x_profile(wlat, wlon, prop_indx[i], nfiles_fg, &argv[start_list_fg], dir_fg, extent_fg, xwork1_fg[sq], xcount1_fg[sq], bdwork1_fg[sq]);
	   }
       }  /* end switch */

       /* free up some memory that never gets used */       
       for ( sq = 0; sq < nwsq1; ++sq) {
            free(xcount1_fg[sq]);
	    xcount1_fg[sq] = NULL;
       }
       	           
        /* Visit each gridpoint in output file and determine whether 
	gaps need to be filled in by fitting.  Produce output profile 
	that includes mask and empty flags where appropriate.  
	Use lat/lon (not row/col) to translate between grids for 
	output/input, seafloor, and work because each has its own 
	dimension and order. */  


        for (row = 0; row < hout.ny; ++row ) {
           lat = hout.y_max - row * hout.y_inc;
           fprintf(stderr,".");
      
           for (col = 0; col < hout.nx; ++col) {
              lon = hout.x_min + col * hout.x_inc;
	 
	      sq_out = row * hout.nx + col;
	 
	  /* reset these at each gridnode*/
	  
	      hwork2.x_min = lon - xrad * hwork2.x_inc;  
	      hwork2.x_max = lon + xrad * hwork2.x_inc; 
	      hwork2.y_min = lat - yrad * hwork2.y_inc;
	      hwork2.y_max = lat + yrad * hwork2.y_inc; 
	 
	 /* set pointers to define small work grids */
	 
	      for (wrow = 0; wrow < hwork2.ny; ++wrow) {
	         for (wcol = 0; wcol < hwork2.nx; ++wcol) {
	    
	           wsq2 = wrow * hwork2.nx + wcol; /* index to small work grids */
                   wlat2 = hwork2.y_min + wrow * hwork2.y_inc;
	           wlon2 = hwork2.x_min + wcol * hwork2.x_inc;
		   
	           wcol1 = NINT((wlon2 - hwork1.x_min + 0.0001) / hwork1.x_inc);
	           wrow1 = NINT((wlat2 - hwork1.y_min + 0.0001) / hwork1.y_inc);
	           wsq1 = wrow1 * hwork1.nx + wcol1;   /* index to big work grids */
	      
		   xwork2[wsq2] = xwork1[wsq1];
	           sig0work2[wsq2] = sig0work1[wsq1];
	           sig1work2[wsq2] = sig1work1[wsq1];
	           sig2work2[wsq2] = sig2work1[wsq1];
	           sig3work2[wsq2] = sig3work1[wsq1];
	           sig4work2[wsq2] = sig4work1[wsq1];
	           bdwork2[wsq2] = bdwork1[wsq1];  /* store value, not ptr */
	    
		   xwork2_fg[wsq2] = xwork1_fg[wsq1];
	           sig0work2_fg[wsq2] = sig0work1_fg[wsq1];
	           sig1work2_fg[wsq2] = sig1work1_fg[wsq1];
	           sig2work2_fg[wsq2] = sig2work1_fg[wsq1];
	           sig3work2_fg[wsq2] = sig3work1_fg[wsq1];
	           sig4work2_fg[wsq2] = sig4work1_fg[wsq1];
	           bdwork2_fg[wsq2] = bdwork1_fg[wsq1];  /* store value, not ptr */
		   
		   if (wsq2 == wsq_mid) {
		     new_count = xcount1[wsq1];
		   }
	      
	         } /* end for wcol */
	      }  /*end for wrow */
	      
	      bottom_de[sq_out] = bdwork2[wsq_mid];
	      
	      /* The first time each output square is visited,
	      construct array of isopycnals corresponding to std depth levels
	      and define a mixed-layer density and depth */
	      
	      if (is_pr) {

                  get_weights(weights, bdwork2_fg, lengthscale, &hwork2);
		  define_mixed_layer(sig0work2_fg[wsq_mid], sig0work2, new_count, weights, &hwork2, xrad, yrad,  bottom_de[sq_out], mlptr[sq_out]);

	          isopyc_above = 0.0;
		  n = 0;
		  if (mlptr[sq_out]->depth > -1.0 && mlptr[sq_out]->depth < 100000) {
		     while (std_depth[n] <= mlptr[sq_out]->depth) {
		        isopyc[sq_out][n] = mlptr[sq_out]->density;
		        ++n;
		     }
		     isopyc_above = mlptr[sq_out]->density + 0.02;  /* by definition, 
		                                      density at bottom of mixed layer */
		  }
	          for (j = n; j < NSTDLEVS; ++j) {
		  
	              curr_depth = std_depth[j];
		      if (j == nz) {
		         curr_depth = bottom_de[sq_out];
		       } 
		       
		      if (curr_depth <= 500 ) {
			      sigptr = sig0work2;
			      sigptr_fg = sig0work2_fg;
		       }
		       else if (curr_depth > 500 && curr_depth <= 1500) {
		              sigptr = sig1work2;
		              sigptr_fg = sig1work2_fg;
		       }
		       else if (curr_depth > 1500 && curr_depth <= 2500) {
		              sigptr = sig2work2;
			      sigptr_fg = sig2work2_fg;
			}
		        else if (curr_depth > 2500 && curr_depth <= 3500) {
		              sigptr = sig3work2;
			      sigptr_fg = sig3work2_fg;
			}
		        else {
		              sigptr = sig4work2;
			      sigptr_fg = sig4work2_fg;
			}
			 
			if (std_depth[j] >= bottom_de[sq_out]) {
			     isopyc[sq_out][j] = (double) HB_f_mask;
			}    
			else if ( sigptr[wsq_mid][j] > 0.0 &&  sigptr[wsq_mid][j] < 100.0 ) {
			     isopyc[sq_out][j] =  sigptr[wsq_mid][j];
			 }			
			else  {  /* estimate an isopycnal for this level*/
			     isopyc[sq_out][j] = sigptr_fg[wsq_mid][j];
		             if (j == nz)  /* special case for bottom */
			         isopyc[sq_out][j] = find_bottom_isopycnal(xrad, yrad, hwork2.nx, sigptr, bdwork2, bottom_de[sq_out], sigptr_fg[wsq_mid][j]);
			  
			}
			
			if ( isopyc[sq_out][j] > 0 && isopyc[sq_out][j] < 100) {
			   /* check for density inversion */
			    if (isopyc[sq_out][j] < isopyc_above - 0.001)    
		                 isopyc[sq_out][j] = isopyc_above;
		            
		            isopyc_above = isopyc[sq_out][j];
			 }
	          } /* end for j */
		  
	     } /*end if is_pr */

             do_fit = 0;
	     mask_it = is_flagged(bottom_de[sq_out], HB_f_mask);     
	     for (j = 0; j < NSTDLEVS; ++j) {
                if (mask_it) { 
		  xfit[j] = (double) HB_f_mask;
		  new_count[j] = 0;
		 }
		else {
	           xfit[j] =  xwork2[wsq_mid][j];
		   if (is_flagged( (float) xfit[j], HB_f_empty) ) {
		      ++do_fit;
		      new_count[j] = 0;
		   }
		}
	     }

	     if (do_fit) {
	     
               if ( ! is_pr ) 
	             get_weights(weights, bdwork2_fg, lengthscale, &hwork2);

            /* Loop for each depth level including bottom ...
	     If level needs to be filled, use isopycnal value estimated for this level
	     and interpolate to find property value on that isopycnal at all the neighboring
	     gridnodes of first-guess and grid being interpolated.  Construct array of
	      anomalies along the isopycnal surface
	     (value - first_guess_value) at each of the neighboring gridnodes and send 
	     array to interp2d().  Add value returned to first-guess value to obtain the 
	     fitted value to be output. */
	     	       
                for (j = 0; j < NSTDLEVS; ++j ) {
	          if (std_depth[j] >= bottom_de[sq_out]) {
		      xfit[j] = (double) HB_f_mask;
		      new_count[j] = 0;
		  }
		  else if ( ! is_flagged((float) xwork2[wsq_mid][j], HB_f_empty) ) {
		      xfit[j] =  xwork2[wsq_mid][j];
		   }
		   else  if ( is_flagged((float) isopyc[sq_out][j], HB_f_empty)) {
		      xfit[j] = (double) HB_f_empty;
		      new_count[j] = 0;
		   }
		   else {    /* isopycnally fit this level */
		       curr_depth = std_depth[j];
	               if (j == nz)
		           curr_depth = bottom_de[sq_out];
			   
		      if (curr_depth <= 500  || curr_depth <= mlptr[sq_out]->depth) {
			      sigptr = sig0work2;
			      sigptr_fg = sig0work2_fg;
		       }
		       else if (curr_depth > 500 && curr_depth <= 1500) {
		              sigptr = sig1work2;
		              sigptr_fg = sig1work2_fg;
		       }
		       else if (curr_depth > 1500 && curr_depth <= 2500) {
		              sigptr = sig2work2;
			      sigptr_fg = sig2work2_fg;
			}
		        else if (curr_depth > 2500 && curr_depth <= 3500) {
		              sigptr = sig3work2;
			      sigptr_fg = sig3work2_fg;
			}
		        else {
		              sigptr = sig4work2;
			      sigptr_fg = sig4work2_fg;
			}
			
                       xsurf = (double *) get_memory((void *)NULL, (size_t) nwsq2, sizeof(double));
		   
                       if ( is_pr && (curr_depth <= mlptr[sq_out]->depth)) {  
			    xfit[j] = hb_p80(curr_depth, (double) lat);
			    if (curr_depth < 5.0)  /* explicitly set surface pressure to 0 */
			       xfit[j] = 0.0;  
			    new_count[j] = -(mlptr[sq_out]->nobs);
		       }
		       else {
 		          xbase = get_neighbors_anom(isopyc[sq_out][j], curr_depth, xsurf, nwsq2, &hwork2, xwork2, sigptr, bdwork2, xwork2_fg, sigptr_fg, bdwork2_fg);
			 
			   if  (is_flagged((float) xbase, HB_f_empty) || is_flagged( (float)xbase, HB_f_mask) ) {
			       /* isopyc not present in first-guess profile, load neighbor values instead of anomalies into xsurf */
			       get_neighbors(isopyc[sq_out][j], curr_depth, xsurf, nwsq2,  xwork2, sigptr, bdwork2);
		               interp2d(xsurf, weights, (double) HBEMPTY, (double) HBMASK, xrad, yrad, &hwork2, &xfit[j], &xstddev, &xwghtavg, &xcount);
			       
		               if (xcount > 0) {
			           if (xcount == 1) ++xcount;  /* don't let this be -1, which is a flag for first-guess value*/
			           new_count[j] =  -xcount;  /* mark interpolated values with negative counts */
		               }
			       else {
			        /* Situation can arise where mixed layer sigma < sig0_fg and the only other 
			neighboring profiles are not along the 8-directions that are searched in interp2d.  
			In this case, revert isopycnal to sig0_fg and estimate value from anomalies */
			         if ((curr_depth <= mlptr[sq_out]->depth)  && (isopyc[sq_out][j] < sigptr_fg[wsq_mid][0]) ) {
				     isopyc[sq_out][j] = sigptr_fg[wsq_mid][0];
 		                     xbase = get_neighbors_anom(isopyc[sq_out][j], curr_depth, xsurf, nwsq2, &hwork2, xwork2, sigptr, bdwork2, xwork2_fg, sigptr_fg, bdwork2_fg);
		                     interp2d(xsurf, weights, (double) HBEMPTY, (double) HBMASK, xrad, yrad, &hwork2, &xanom, &xstddev, &xwghtavg, &xcount);
		                     if (xcount > 0) {
		                        xfit[j] = xbase + xanom;
			                new_count[j] =  - (xcount +1);  /* mark interpolated values with negative counts */
		                     }
		                     else {
		                        xfit[j] = xbase;
			                new_count[j] = -1;   /* mark squares filled only with first guess value as -1 */
		                     }
				  }
			       }
		            }
		            else {
		              interp2d(xsurf, weights, (double) HBEMPTY, (double) HBMASK, xrad, yrad, &hwork2, &xanom, &xstddev, &xwghtavg, &xcount);
			      
		              if (xcount > 0) {
		                xfit[j] = xbase + xanom;
			        new_count[j] =  - (xcount +1);  /* mark interpolated values with negative counts */
		              }
		              else {
		                 xfit[j] = xbase;
			         new_count[j] = -1;   /* mark squares filled only with first guess value as -1 */
		              }
			      
			      if (is_pr)  {  /* check for estimated pr > bottom pr */
			         if ( xfit[j] >  hb_p80((double)bottom_de[sq_out], (double)lat) ) {
			             xfit[j] =  hb_p80((double) bottom_de[sq_out], (double)lat);
			         }
			     }
	 /* insert code to determine stddev of xvalues (not anomaly values) for
         isopycnal and depth surfaces.  Use get_neighbors() to fill xsurf
	 with xvalues on each type of surface.  Then call interp2d() to get stddev */
	 
		         }
		     }	
		     	  
		     if (is_pr &&  (j == nz))  { /* explicitly set bottom pressure */
			    xfit[j] = hb_p80(curr_depth, (double) lat);
		     }    
		     free((void *)xsurf);
	         } /* end else isopycnally fit this level */
	      } /* end for j */
	   } /* end if do_fit */
	   
           if (is_pr) {
	           /* save pressure profile with fitted values.  These will be used 
		to convert potential temperatures back to in situ values. 
		Compute depths associated with fitted profile for interpolating
		all property profiles back onto standard depth levels. */
	   
	         pfit[sq_out] = (double *) get_memory((void *)NULL, (size_t)NSTDLEVS, sizeof(double));
	         dfit[sq_out] = (double *) get_memory((void *)NULL, (size_t)NSTDLEVS, sizeof(double));
	         for (j = 0; j < nz; ++j) {
		 
		     pfit[sq_out][j] =  xfit[j];
		     dfit[sq_out][j] = xfit[j];
		     if ((xfit[j] > -1.0) && (xfit[j] < 20000))
		        dfit[sq_out][j] = hb_depth(pfit[sq_out][j], (double) lat);
		 } /* for j */
		 
		 /* explicitly set bottom pressure and bottom depth */
		 
		 if ((xfit[nz] > -1.0) && xfit[nz] < 20000.) {
                    dfit[sq_out][nz] = (double) bottom_de[sq_out] + 5;  /* add a bit for interpolating*/
	            pfit[sq_out][nz] =  hb_p80(dfit[sq_out][nz] , (double)lat);
	         }
	   }	   
	   if (is_sa) {  
	      
	           /* save salinity profile with fitted values 
		These will be used to convert potential
		   temperatures back to in situ values. */
		      
	         sfit[sq_out] = (double *) get_memory((void *)NULL, (size_t)NSTDLEVS, sizeof(double));
		      
	         for (j = 0; j <= nz; ++j) {
		     sfit[sq_out][j] =  xfit[j];
		 }
	   }
	      
	   if (is_te) {  /* adiabatically adjust potential temperature to pressure */
	         xtmp_d = xfit;
	         xfit = (double *) get_memory((void *)NULL, (size_t)NSTDLEVS, sizeof(double));
		 
	         for (j = 0; j < NSTDLEVS; ++j) {
	    
	             if (xtmp_d[j] > -3.0 && xtmp_d[j] < 100.0) {
	                xfit[j] = hb_theta(sfit[sq_out][j], xtmp_d[j], 0.0, pfit[sq_out][j]);
			if ( ! (xfit[j] > -3.0 && xfit[j] < 100.0) ) {
			   xfit[j] = (double) HB_f_empty;
			}
		    }
		     else if (xtmp_d[j] < -9.0) {
	                xfit[j] = (double) HB_f_empty;
	             }
	  	     else {
	                xfit[j] = (double) HB_f_mask;
                     }
	       
	         }  /* end for j */
	    
	         free((void *) xtmp_d);
		 free((void *) sfit[sq_out]);
		 free((void *) pfit[sq_out]);
	    
	   }  /* end if is_te */


          /* interpolate fitted properties in xfit back onto standard depths and 
	  store in a float array for output to nedcdf files.*/

          xout = (float *) get_memory((void *)NULL, (size_t)NSTDLEVS, sizeof(float));  
	 
	  x_to_stddepth (xfit, new_count, xout, dfit[sq_out], std_depth, bottom_de[sq_out], NSTDLEVS);
	  	   
	   write_prop_cdf(outfile, xout, get_prop_mne(prop_indx[i]), row, col, tbin, 0, 1,1,1,NSTDLEVS);
	   write_prop_count_cdf(outfile, new_count, get_prop_mne(prop_indx[i]), row, col, tbin, 0, 1,1,1,NSTDLEVS);

	   free((void *)xout);
	   
         }  /* end for col */
       } /* end for row */
     
       if (! is_pr ) {  /* pr surf is used for interp2d() of all props in upper ocean  */
           for (sq = 0; sq < nwsq1; ++sq) {
	      free((void *) xwork1[sq]);
	      free((void *) xwork1_fg[sq]);
	      free((void *) xcount1[sq]);
	   }
       }
      
    } /* end for i*/
  
   write_bottom_depth_cdf(outfile, 0, 0, 0, hout.ny, hout.nx, 1, bottom_de);	 
   
   cdf_close(outfile); 
   
   fprintf(stderr,"\nEnd of %s.\n", argv[0]);
   exit(0);
	
}  /* end main */


/************************************************************************/
void print_usage(char *program)
{

  fprintf(stderr,"\nIsopycnally interpolates missing gridpoints in HydroBase cdf files ");
  fprintf(stderr,"\nusing a first-guess grid to define a background density structure ");
  fprintf(stderr,"\nand property difference fields.  Whether or not a gridnode is missing ");
  fprintf(stderr,"\nis determined with real seafloor topography. For each missing value, ");
  fprintf(stderr,"\nits isopycnal is determined from the first-guess field at that gridnode.  ");
  fprintf(stderr,"\nNext, property differences between first-guess field and ");
  fprintf(stderr,"\ngrid-being-interpolated are computed along that isopycnal surface ");
  fprintf(stderr,"\nfor surrounding gridnodes.  A difference value is computed for the missing");
  fprintf(stderr,"\ngridnode using a distance-weighted (gaussian) interpolation. The interpolated ");
  fprintf(stderr,"\ndifference value is added to the first-guess property value to estimate a new value");
  fprintf(stderr,"\nat the missing gridnode.  A search ellipse (xradius, yradius) parameter ");
  fprintf(stderr,"\ndetermines the maximum acceptable distance an interpolated value ");
  fprintf(stderr,"\ncan be from an observed value in each direction searched (8 compass pts). ");
  fprintf(stderr,"\nIf the isopycnal runs into topography, the search is suspended in that direction. ");
  fprintf(stderr,"\nThe input domain is expanded relative to the output domain ");
  fprintf(stderr,"\nso that data around borders of the grid can be incorporated.");
  fprintf(stderr,"\nAll input cdf files must possess the same depth dimensions");
  fprintf(stderr,"\nand grid spacing.  The output cdf file will have the");
  fprintf(stderr,"\nsame depth dimension, but the lat/lon dimensions specified");
  fprintf(stderr,"\nwith -B<w/e/s/n> and -I<x/yincr>");
  
  fprintf(stderr,"\n\nUSAGE:  %s cdf_file_root_name(s) -B<w/e/s/n> -I<x_inc>[/<y_inc>] -O<outfile> -P<properties>  [-A] [-C<dir_for_first_guess_files>] [-D<input_dir>] [-E<input_file_extent>]  [-E<extent_for_first_guess_files>]  [-G] [-L<lengthscale>]  [-S<xradius>[/<yradius>]] [-h]\n\n", program);
  fprintf(stderr," -B  sets the output grid bounds w/e/s/n.\n");
  fprintf(stderr," -I  sets the output grid spacing for x/y dimensions.\n");
  fprintf(stderr," -O  name of output cdf file.\n");
  fprintf(stderr," -P  list of properties to include in output file\n");
  fprintf(stderr,"        ex:  -Ppr/th/sa/ox/ht\n");
  fprintf(stderr,"       -P (by itself) produces a list of available properties\n");
  fprintf(stderr, "\n\tOPTIONS:\n");
  fprintf(stderr,"-A list of first-guess file_root_names follows (leave a blank space before list).\n");
  fprintf(stderr,"     If not specified, cdf_file_roots of previous list will be used with extent specifed by -F \n");
  fprintf(stderr,"-C different directory for first-guess cdf files.\n");
  fprintf(stderr,"     If not specified, directory  specifed by -D will be used \n");
  fprintf(stderr,"-D directory for input cdf files.\n");
  fprintf(stderr,"-E file extent for input cdf files.\n");
  fprintf(stderr,"-F file extent for first-guess cdf files.\n");
  fprintf(stderr,"-G force gridnode registration for ouput (default is pixel registration).\n");
  fprintf(stderr,"-L lengthscale (in km) for distance weighting =  e ^[-(dist/L)^2]\n"); 
   fprintf(stderr,"          default is L =%.1lf  km\n", lengthscale);
  fprintf(stderr,"-S search for a non-empty gridnode within <xradius/yradius> of an empty node\n");
  fprintf(stderr,"      in the initial grid.  Estimate a value only if a non-empty node is found.\n");
  fprintf(stderr,"      radius is specified in xy-units (degrees) .\n");
  fprintf(stderr,"      <radius> = 0 interpolates no gridnodes but puts in masking info.\n"); 
  fprintf(stderr,"      Default is [%.1lf/%.1lf]\n", xradius, yradius);
  fprintf(stderr,"-h help...... prints this message. \n");
  return;
  
} /* end print_usage() */
/*****************************************************************************/
int parse_p_option(char *st, int *prop_indx)
{
  char prop[4];
  int n;

  if (*st == '\0') {
         print_prop_menu();
         exit(0);
  }
  
  
  /* Explicitly set the pr, te, sa because they are mandatory */
  
  prop_indx[0] = (int)PR;
  prop_indx[1] = (int)TE;
  prop_indx[2] = (int)SA;
  
  n = 3;
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

     /* !**!  Special cases for properties ...
     
       No need for pref specs for props like s_ , ht and pe.  No computation of 
       additional props will be done. Any property 
       requested must already exist in the cdf file.  */
     
         /* de, pr, te, sa are automatically done so don't count them here */

     if ( !( (prop_indx[n] == (int)DE) || (prop_indx[n] == (int)PR) || (prop_indx[n] == (int)TE) || (prop_indx[n] == (int)SA)))  
            ++n;
     

   } while (*st == '/');
   return (n);
}  /* end parse_p_option() */

/*****************************************************************************/
int cdf_construct(struct GRID_INFO *hptr, int nprops, int *prop_indx,  char *filename, int nargs, char **arglist)

   /* Opens a cdf output file and writes an appropriate header,
      standard depths, and time bin info.  Returns the id associated
      with the open file */
{
   struct CDF_HDR cdfhdr;
   int i, error, cdfid;
  
   cdfhdr.xmax = (float) hptr->x_max;
   cdfhdr.xmin = (float) hptr->x_min;
   cdfhdr.ymax = (float) hptr->y_max;
   cdfhdr.ymin = (float) hptr->y_min;
   cdfhdr.xincr = (float) hptr->x_inc;
   cdfhdr.yincr = (float) hptr->y_inc;
   
   cdfhdr.nx = (int) NINT((hptr->x_max - hptr->x_min)/hptr->x_inc);
   cdfhdr.ny = (int) NINT((hptr->y_max - hptr->y_min)/hptr->y_inc);
   cdfhdr.node_offset = hptr->node_offset;
   if (! hptr->node_offset) {  
      ++cdfhdr.nx;     
      ++cdfhdr.ny;
   }
   
   cdfhdr.nz = NSTDLEVS;
   cdfhdr.nt = 1;
   cdfhdr.tmin = (int *) malloc(sizeof(int));
   cdfhdr.tmax = (int *) malloc(sizeof(int));
   cdfhdr.tmin[0] = 0;
   cdfhdr.tmax[0] = 9999;
   cdfhdr.counts_included = 1;
   cdfhdr.fill_value =  HB_f_empty;
   cdfhdr.mask_value =  HB_f_mask;
   strncpy(cdfhdr.x_units, "degrees", 8);
   strncpy(cdfhdr.y_units, "degrees", 8);
   strncpy(cdfhdr.z_units, "meters", 7);
   strncpy(cdfhdr.title,"HydroBase", 10);
   strcpy(cdfhdr.command, *arglist);
   i = 1;
   while ( (i < nargs) && (strlen(arglist[i]) + strlen(cdfhdr.command) < 900) ) {
          strncat(cdfhdr.command, " ", 1);
          strcat(cdfhdr.command, arglist[i]);
          ++i;	  
   }
   
   cdfhdr.nprops = nprops;
   cdfhdr.prop_id = (char **) malloc(nprops * sizeof(char *));
   cdfhdr.prop_units = (char **) malloc(nprops * sizeof(char *));
   for (i = 0; i < nprops; ++i) {
      cdfhdr.prop_id[i] = (char *) malloc(3);
      cdfhdr.prop_units[i] = (char *) malloc(50);
      strncpy(cdfhdr.prop_id[i], get_prop_mne(prop_indx[i]),3);
      strcpy(cdfhdr.prop_units[i], get_prop_units(prop_indx[i]));
   }
   
/* Open output file and write out some info so we can free up some memory... */
   
   cdfid = cdf_init(filename);   
   error = cdf_define(cdfid, &cdfhdr, PREFILL, cdfhdr.counts_included);
   if (error)  exit(1);
   
   error = write_std_depths_cdf(cdfid, &cdfhdr);
   error = write_time_bins_cdf(cdfid, &cdfhdr);
   
   for (i = 0; i < nprops; ++i) {
      free((void *)cdfhdr.prop_id[i]);
      free((void *)cdfhdr.prop_units[i]);
   }
   free((void *) cdfhdr.prop_id);
   free((void *) cdfhdr.prop_units);
  
   return(cdfid); 
} /*end cdf_construct() */ 

/************************************************************************/
void get_seafloor(FILE *fptr, double dx, double dy, short *ztopo, struct GRID_INFO *hptr)
   /* reads in topography values from an already open file.  Returns the
      region correponding to the global variables xmin/xmax/ymin/ymax in
      ztopo. 
      An error causes a message to be printed and an exit.   */
{
  int row, col, nrowsin, ncolsin, n, ix, iy;
  int  lon0to360;
  double lat, lon;
  short **z;

  nrowsin = NINT(180/dy) + 1;  /* by definition topography domain is global */
  ncolsin = NINT(360/dx) + 1;
  
   /* Allocate memory */
   
   z = (short **) malloc(nrowsin * sizeof(short *));
   if (z == NULL) {
      fprintf(stderr,"\nError allocating memory in get_seafloor().\n");
      exit(1);
   }
   for (row = 0; row < nrowsin; ++row) {
     z[row] = (short *) malloc(ncolsin * sizeof(short));
     if (z[row] == NULL) {
         fprintf(stderr,"\nError allocating memory in get_seafloor().\n");
         exit(1);
      }
   }
   
   fprintf(stderr,"\nReading in topography values... ");
   for (row = 0; row < nrowsin; ++row) {
     n = fread((short *)z[row], sizeof(short), ncolsin, fptr);
     if (n != ncolsin) {
        fprintf(stderr,"\nError reading the topofile at row %d\n", row);
        fprintf(stderr,"Not enough 2-byte integer values in file.\n");
        exit(1);
     }
   }
   
   fprintf(stderr,"finished.\n");
   fclose(fptr);
   
   /* extract values at gridnodes corresponding to 3D grid area
     and negate to make seafloor values positive. Store in row order to
     match other work grids:  (sq = row * ncols + col) */
   n = 0;
   for (iy = 0; iy < hptr->ny; ++iy) {
      lat = hptr->y_min + iy * hptr->y_inc;
      row = (int) (NINT((lat - (-90.)) / dy ) + 0.0001); 

      if (row < 0)
          row = 0;
      if (row >= nrowsin)
           row = nrowsin-1;   
	     
      for (ix = 0; ix < hptr->nx; ++ix) {
         lon = hptr->x_min + ix * hptr->x_inc;
         lon0to360 = 0;
         if (lon < 0)
           lon0to360 = 360;
         col = (int) (NINT((lon + lon0to360) / dx) + .0001);
         ztopo[n++] =  - z[row][col];
      }
   }
        
   for (row = 0; row < nrowsin; ++row) 
    free((void *)z[row]);
    
   free((void *) z);
   
   return;
  
   
}  /* end get_seafloor() */
/************************************************************************/

void get_x_profile(float lat, float lon, int pindex, int nfiles, char **argv, char *dir, char *ext, double *x, short *count, float bdepth )

 /*  Reads each cdf file in list and extracts depth profile of property
     indexed by pindex.  Returned array contains either a value, HB_f_empty
     of HB_f_mask.  If a count variable is present the info is 
     stored in the count array, otherwise the count array is assigned
     either 0 or 1.    */

{
   struct CDF_HDR cdf;
   int cdfid, curfile, print_msg = 0;
   int lon0to360;
   float *xin;
   short *nobs;
   char *mne;
   int error, i,  row, col;
   int tbin = 0, lev, new_bd, bindex;
   float bd;
   
   

   lon0to360 = 1;
   if (lon < 0)
      lon0to360 = 0;

/* Loop for each input file */

   curfile = 1;
   
   do {
      cdfid = cdf_open(dir, argv[curfile], ext, print_msg);
      if (cdfid < 0)
         goto NEXTFILE;
	 
      if (error = read_cdf_hdr(cdfid, &cdf)) 
         exit (1);
	 
/* compare bounds of file to the lat/lon specified.  Try
   to match the convention for longitude -- usually W is neg, E is pos --
   but we cannot be certain it is always specified like this .*/
	 
      if (lon0to360) {    /* target lon is positive */
         if (cdf.xmax < 0) {
	   cdf.xmax += 360.0;
           if (cdf.xmin < 0)
	       cdf.xmin += 360.0;
         }
      }
      else {              /* target lon is negative */
           if (cdf.xmin > 0) {
	      cdf.xmin -= 360.0;
              if (cdf.xmax > 0)
	         cdf.xmax -= 360.0;
	   }
      } 

      if ((cdf.xmin > lon) || (cdf.xmax < lon) 
      || (cdf.ymin > lat) || (cdf.ymax < lat)) {
         cdf_close(cdfid);
	 goto NEXTFILE;
      }


      error = get_indices(&cdf, lat, lon, &row, &col);
      
      if (row < 0 || row >= cdf.ny || col < 0 || col >= cdf.nx) {
         cdf_close(cdfid);
         goto NEXTFILE;
      }
      
/* check that cdf file has same number of standard depths  */
     
      if (cdf.nz != NSTDLEVS) {
         fprintf(stderr, "\nFATAL ERROR:  Mismatch of standard depths in cdf files.\n");  
         exit(1);
      }


/* determine if property is available in this file */

      i = 0;
      while (i < cdf.nprops && !(pindex == get_prop_indx(cdf.prop_id[i])))
         ++i;
	 
      	 
      if (i >= cdf.nprops) {
      
            fprintf(stderr, "\nWARNING: Property [%2s] not available in this file.", get_prop_mne(pindex));
            cdf_close(cdfid);
	    goto NEXTFILE;
      } 
	 
      
/* allocate space for property arrays; */

      xin = (float *) get_memory((void *)NULL, (size_t)cdf.nz, sizeof(float));
      nobs = (short *) get_memory((void *)NULL, (size_t)cdf.nz, sizeof(nobs));
      

      mne = get_prop_mne(pindex);
      error = read_cdf_prop(cdfid, mne, xin, row, col, tbin, 0, cdf.nz);
      if (error > 0) {
            fprintf(stderr,"\nError attempting to read %.2s at row,col =  %d,%d from cdf file.", mne, row, col);
            exit(1);
      }
	   
     if (cdf.counts_included) {
        error = read_cdf_prop_count(cdfid,mne,nobs,row,col,tbin,0,cdf.nz);
        if (error > 0) {
            fprintf(stderr,"\nError attempting to read %.2s_cnt at row,col =  %d,%d from cdf file.", mne, row, col);
            exit(1);
        }

     }
     
     read_cdf_bottom_depth(cdfid, &bd, row, col, tbin);
     new_bd = 0;
     if (bd < (bdepth-10.0))
          new_bd = 1;    /* this bottom depth is not the real bottom */
	
     /* check for flagged levels.  in case of multiple entries,
        sum up the data at each level.  */ 
	
     bindex = cdf.nz -1;	   
     for (lev = 0; lev < bindex; ++lev) {
       if ( ! is_flagged(xin[lev], cdf.fill_value) && ! is_flagged(xin[lev], HB_f_mask)) {
       
	   if (!cdf.counts_included || nobs[lev] == 0) 
	      nobs[lev] = 1;
	   
	   x[lev] += (double) (ABS(nobs[lev]) * xin[lev]);
	   if (count[lev] * nobs[lev] < 0)  /* test for difference in sign*/
	       count[lev] = ABS(count[lev]) + ABS(nobs[lev]);
	    else
	       count[lev] += nobs[lev];
	}
     } /* end for lev */

/* add bottom value only if it is a true bottom value */     

     if (!new_bd) {
         if ( ! is_flagged(xin[bindex], cdf.fill_value)  && ! is_flagged(xin[bindex], HB_f_mask)) {
	   if (!cdf.counts_included || nobs[bindex] == 0) 
	      nobs[bindex] = 1;
	      
	   x[bindex] += (double) (ABS(nobs[bindex]) * xin[bindex]);
	   if (count[bindex] * nobs[bindex] < 0)  /* test for difference in sign*/
	       count[bindex] = ABS(count[bindex]) + ABS(nobs[bindex]);
	    else
	      count[bindex] += nobs[bindex];
	 }
     }
     
     free((void *) xin);
     free((void *) nobs);
     cdf_close(cdfid);
     
NEXTFILE:
      ;
      	 
   }  while (curfile++ < nfiles); 
   
   
   /* flag levels with zero count */

   for (lev = 0; lev < cdf.nz; ++lev) {
   
      if (count[lev]) {
         x[lev] /= ABS(count[lev]);
      }
      else {
         x[lev] = HB_f_empty;
        if ( !is_flagged(bdepth, HB_f_empty) && !is_flagged(bdepth, HB_f_mask)) {	 
           if (std_depth[lev] >= bdepth)
                 x[lev] = (double) HB_f_mask; 
	} 
      }
   
   } /* end for lev */
   
      
   return;

}  /* end get_x_profile() */

/***************************************************/
int define_bottom(float *bd_addr, short seafloor)
  /* Determines an appropriate bottom depth based on
      observed data  plus seafloor topography and sets *bd_addr.
      Returns 1 if the bottom depth was changed or
      0 if not changed.
  */
{

  int new_bd;
   
  /* Compare seafloor depth and bottom depth.  
     IF bottom_de is defined: 
         Seafloor negative?
	   do nothing;
         Seafloor positive:
	    Compare seafloor/bottom_de:  
	      Seafloor deeper:
	        set bottom_de to seafloor
	      bottom_de deeper:
	        do nothing 
     ELSE bottom_de is undefined:
         Seafloor negative?  set bottom_de to mask value.  
	 Seafloor positive:  set bottom_de to seafloor-10.
  */

     
  new_bd = 0; 
      
  if (!is_flagged(*bd_addr,HB_f_empty) && !is_flagged(*bd_addr, HB_f_mask)) {  
        
	/* bottom_depth is defined */
     
     if (seafloor > 0){           /* negative seafloor means land */
          if ( ((float)seafloor - *bd_addr) > 100) {
             *bd_addr = seafloor - 10;
	     new_bd = 1;
	  }
	  
     }
  }  
  else {      /* bottom_de not defined */
     
    if (seafloor < 1) 
      *bd_addr = HB_f_mask;
    else 
      *bd_addr =  seafloor;

  }
  return(new_bd);       
  
}  /* end define_bottom() */	
/************************************************************************/

void define_mixed_layer( double *sig_profile_fg, double **sigwork,  short *count, double *weights, struct GRID_INFO *hptr, int xrad, int yrad, float bdepth, struct MIXLAYER *mixlayptr)
     /* Determines density of mixed layer from first guess sigma profile and surrounding sigwork
       profiles. Returns density and depth of mixed layer in the structure. 
    */     
{
   double delta, maxsigma, dummy, sigma, *xsurf, *mldepths ;
   int j, nwsq, sq, sqmid, n;
   
   sqmid = yrad * hptr->nx + xrad;
   delta = 0.02;   /* arbitrary definition of mixed layer sigma range */
   mixlayptr->density = (double) HB_f_empty;
   mixlayptr->depth =  (double) HB_f_empty;
   mixlayptr->nobs = 0;
      
 /* If sigma profile in work array already exists, fill in struct with density/depth values 
    and return */   
    j = 0;
   if (sigwork[sqmid][j] > 0.0 && sigwork[sqmid][j] < 100.0) {  
      mixlayptr->density = sigwork[sqmid][j];
      mixlayptr->depth = std_depth[j];
      mixlayptr->nobs = count[j];
      
      maxsigma = mixlayptr->density +delta;
      while ( (++j < NSTDLEVS) && (sigwork[sqmid][j] <= maxsigma) && (sigwork[sqmid][j] > 0.0)) {
         mixlayptr->depth = std_depth[j];
      }
      if (! (is_flagged ( (float) sigwork[sqmid][j], HB_f_empty) || is_flagged((float)sigwork[sqmid][j], HB_f_mask)) )
          return;
   }

   /* Otherwise, estimate mixed layer params from work arrays */  
    
   nwsq = hptr->ny * hptr->nx;
   xsurf = (double *) get_memory((void *) NULL, (size_t)nwsq, sizeof(double));
   mldepths = (double *) get_memory((void *) NULL, (size_t)nwsq, sizeof(double));
   n = 0;
    
   for (sq = 0; sq < nwsq; ++sq) {
 
      j = 0;
      sigma = sigwork[sq][j];    
      xsurf[sq] = sigma;    
      mldepths[sq] = sigma;   /* in case it is flagged */
      
      /* find depth range over which density is within delta sigma units of surface */  
      if (sigma > 0.0 && sigma < 100.0) {
        maxsigma = sigma + delta;
	
	++j;
        while ((j < NSTDLEVS) && (sigwork[sq][j] <= maxsigma) && (sigwork[sq][j] > 0.0))
	    ++j;
	    
	 mldepths[sq] = std_depth[--j]; 
	 ++n;
      }
      
   }  /* end for sq */
   
   /* add first_guess profile */
   
   sigma = sig_profile_fg[0];
   if ( (sigma > 0.0) && (sigma < 100.0)) {
      xsurf[sqmid] = sigma;
      maxsigma = sigma + delta;
      j = 1;
      while ((j < NSTDLEVS) && (sig_profile_fg[j] <= maxsigma) && (sig_profile_fg[j] > 0.0))
	    ++j;
	    
      mldepths[sqmid] = std_depth[--j];
      ++n;
    }
    
       interp2d(mldepths, weights, (double)HB_f_empty, (double)HB_f_mask, xrad, yrad, hptr, &mixlayptr->depth, &dummy, &dummy, &mixlayptr->nobs);
      interp2d(xsurf, weights, (double)HB_f_empty, (double)HB_f_mask, xrad, yrad, hptr, &mixlayptr->density, &dummy, &dummy, &mixlayptr->nobs);
     
     return;
   
}  /* end define_mixed_layer() */
/************************************************************************/
void get_sig_profiles(float lat, float lon, float *bdepth_ptr, short seafloor, int nfiles, char **argv, char *dir, char * ext, double *ppro, double *tpro, double *spro, short *pcount, short *tcount, short *scount, double *s0pro, double *s1pro, double *s2pro, double *s3pro, double *s4pro)

/*  Search cdf files for lon/lat, read in pr, te, sa, and compute all sigmas.  
    Read in bottom depth into the pointer and compare to seafloor to define
    a realistic bottom depth.  Mark masked levels with HB_f_mask.  
    Missing values are flagged with HB_f_empty.  Returns 1 if a new bottom depth 
    was defined or 0 if bottom_depth remains same at this square.

*/
{
   struct CDF_HDR cdf;
   int cdfid, curfile, print_msg = 0;
   int lon0to360, new_bd;
   float deepest, bd;
   float *xin;
   short *nobs;
   int error, i,  row, col;
   int tbin = 0, lev;
   int maxlev, bindx;

  /* initialize these */
  
   lon0to360 = 1;
   if (lon < 0)
      lon0to360 = 0;

   bindx = NSTDLEVS - 1;
   bd = HB_f_empty;
   deepest = HB_f_empty;           
   
   xin = (float *) get_memory((void *)NULL, (size_t)NSTDLEVS, sizeof(float));
   nobs = (short *) get_memory((void *)NULL, (size_t)NSTDLEVS, sizeof(short));

/* Loop for each input file */

   curfile = 1;
   
   do {
      cdfid = cdf_open(dir, argv[curfile], ext, print_msg);
      if (cdfid < 0)
         goto NEXTFILE2;
	 
      if (error = read_cdf_hdr(cdfid, &cdf)) 
         exit (1);
	 
/* compare bounds of file to the lat/lon specified.  Try
   to match the convention for longitude -- usually W is neg, E is pos --
   but we cannot be certain it is always specified like this .*/
	 
      if (lon0to360) {    /* target lon is positive */
        if (cdf.xmax < 0)  {
	   cdf.xmax += 360.0;
           if (cdf.xmin < 0) 
	      cdf.xmin += 360.0;
        }
      }
      else {              /* target lon is negative */
           if (cdf.xmin > 0)  {
	      cdf.xmin -= 360.0;
              if (cdf.xmax > 0)
	         cdf.xmax -= 360.0;
	   }	 
      } 

      if ((cdf.xmin > lon) || (cdf.xmax < lon) 
      || (cdf.ymin > lat) || (cdf.ymax < lat)) {
         cdf_close(cdfid);
	 goto NEXTFILE2;
      }


      error = get_indices(&cdf, lat, lon, &row, &col);
      
      if (row < 0 || row >= cdf.ny || col < 0 || col >= cdf.nx) {
         cdf_close(cdfid);
         goto NEXTFILE2;
      }
      
/* check that cdf file has same number of standard depths  */
     
      if (cdf.nz != NSTDLEVS) {
         fprintf(stderr, "\nFATAL ERROR:  Mismatch of standard depths in cdf files.\n");  
         exit(1);
      }

   read_cdf_bottom_depth(cdfid, &bd, row, col, tbin);
   if (is_flagged(bd, cdf.fill_value) ) 
       bd = HB_f_empty;
       
   new_bd = define_bottom(&bd, seafloor);

   if (!is_flagged(bd, HB_f_mask) && (bd > deepest)) { 
      deepest = bd;        /* this profile is the deepest */   
      ppro[bindx] = 0.0;   /* zero previous bottom obs to handle the case */
      pcount[bindx] = 0;   /* of multiple entries for same square */
      tpro[bindx] = 0.0;
      tcount[bindx] = 0;  
      spro[bindx] = 0.0;
      scount[bindx] = 0;  
   }
   
/******************/      
/* extract pr, te, sa at all depth levels.  For now, set masked or missing values
       to zero  */
    

   
     error = read_cdf_prop(cdfid, "te", xin, row, col, tbin, 0, cdf.nz);
     if (error) {
       fprintf(stderr, "\n No te variable in cdf file:  %s\n", argv[curfile]);
       exit(1);
     }
     
     if (cdf.counts_included) {
        error = read_cdf_prop_count(cdfid,"te",nobs,row,col,tbin,0,cdf.nz);
        if (error > 0) {
            fprintf(stderr,"\nError attempting to read te_cnt from cdf file.");
            exit(1);
        }

     }
     
     if (new_bd) {        /* cancel bottom observation */
       xin[bindx] = cdf.fill_value;
       nobs[bindx] = 0;
     }
     
     for (lev = 0; lev < NSTDLEVS; ++lev) {
     
       if (is_flagged(xin[lev], cdf.fill_value) || is_flagged(xin[lev], HB_f_mask) ) {
          xin[lev] = 0.0;
       }  
       else {
          if (cdf.counts_included){
	     if (tcount[lev] * nobs[lev] < 0 )  /* test for different sign */
	        tcount[lev] = ABS(tcount[lev]) + ABS(nobs[lev]);
	     else
	        tcount[lev] += nobs[lev];
	  }
	  else {
             ++tcount[lev];
	     nobs[lev] = 1;
	  }
          tpro[lev] += (double) (xin[lev] * ABS(nobs[lev]));
       }
     }
  
     error = read_cdf_prop(cdfid,"sa", xin, row, col, tbin, 0, cdf.nz);
     if (error) {
       fprintf(stderr, "\n No sa variable in cdf file:  %s\n", argv[curfile]);
       exit(1);
     }
     if (cdf.counts_included) {
        error = read_cdf_prop_count(cdfid,"sa",nobs,row,col,tbin,0,cdf.nz);
        if (error > 0) {
            fprintf(stderr,"\nError attempting to read sa_cnt from cdf file.");
            exit(1);
        }
     }

     if (new_bd) {        /* cancel bottom observation */
       xin[bindx] = cdf.fill_value;
       nobs[bindx] = 0;
     }
     
     for (lev = 0; lev < NSTDLEVS; ++lev) {
     
       if (is_flagged(xin[lev], cdf.fill_value) || is_flagged(xin[lev], HB_f_mask) ) {
          xin[lev] = 0.0;
       }  
       else {
          if (cdf.counts_included) {
	     if (scount[lev] * nobs[lev] < 0 )  /* test for different sign */
	        scount[lev] = ABS(scount[lev]) + ABS(nobs[lev]);
	     else
	        scount[lev] += nobs[lev];
	  }
	  else {
             ++scount[lev];
	     nobs[lev] = 1;
	  }
          spro[lev] += (double) (xin[lev] * ABS(nobs[lev]));
       }

     }
     error = read_cdf_prop(cdfid,"pr", xin, row, col, tbin, 0, cdf.nz);
     if (error) {
       fprintf(stderr, "\n No pr variable in cdf file:  %s\n", argv[curfile]);
       exit(1);
     }
     if (cdf.counts_included) {
        error = read_cdf_prop_count(cdfid,"pr",nobs,row,col,tbin,0,cdf.nz);
        if (error > 0) {
            fprintf(stderr,"\nError attempting to read pr_cnt from cdf file.");
            exit(1);
        }
     }
     
     if (new_bd) {        /* cancel bottom observation */
       xin[bindx] = cdf.fill_value;
       nobs[bindx] = 0;
     }
     
     for (lev = 0; lev < NSTDLEVS; ++lev) {
     
       if (is_flagged(xin[lev], cdf.fill_value) || is_flagged(xin[lev], HB_f_mask) ) {
          xin[lev] = 0.0;
       }  
       else {
          if (cdf.counts_included)  {
	     if (pcount[lev] * nobs[lev] < 0 )  /* test for different sign */
	        pcount[lev] = ABS(pcount[lev]) + ABS(nobs[lev]);
	     else
	        pcount[lev] += nobs[lev];
	  }
	  else {
             ++pcount[lev];
	     nobs[lev] = 1;
	  }
          ppro[lev] += (double) (xin[lev] * ABS(nobs[lev]));
       }
    }
  
     
     cdf_close(cdfid);

NEXTFILE2:
      ;	 
   }  while (curfile++ < nfiles); 
   
/*****************************************************/ 
/* find average in case there were multiple entries */

   for (lev = 0; lev < NSTDLEVS; ++lev) {
   
     if (pcount[lev])
       ppro[lev] /= ABS(pcount[lev]);
     else 
       ppro[lev] = HB_f_empty;
   
     if (tcount[lev])
       tpro[lev] /=  ABS(tcount[lev]);
     else 
       tpro[lev] = HB_f_empty;
   
     if (scount[lev])
       spro[lev] /= ABS(scount[lev]);
     else 
       spro[lev] = HB_f_empty;
   
   }

   free((void *)xin);
   free((void *)nobs);
     
  /* gridpoint is on land, mask all arrays */
  
   if ( is_flagged(bd, HB_f_mask) && deepest < 0.0) {
   
      *bdepth_ptr = HB_f_mask;
   
       for (lev = 0; lev <= bindx; ++lev) {
          s0pro[lev] = HB_f_mask;
          s1pro[lev] = HB_f_mask;
          s2pro[lev] = HB_f_mask;
          s3pro[lev] = HB_f_mask;
          s4pro[lev] = HB_f_mask;
          ppro[lev] = HB_f_mask;
          tpro[lev] = HB_f_mask;
          spro[lev] = HB_f_mask;
      }
      return;
   }
   
   
  /* Not on land ...compute sigmas */
  
   compute_sigma(0.0, cdf.nz, s0pro, ppro, tpro, spro);
   compute_sigma(1000.0, cdf.nz, s1pro, ppro, tpro, spro);
   compute_sigma(2000.0, cdf.nz, s2pro, ppro, tpro, spro);
   compute_sigma(3000.0, cdf.nz, s3pro, ppro, tpro, spro);
   compute_sigma(4000.0, cdf.nz, s4pro, ppro, tpro, spro);

/* replace missing values with empty flag */

   for (lev = 0; lev <= bindx; ++lev) {
        if (s0pro[lev] < 0) {
             s0pro[lev] = HB_f_empty;
             s1pro[lev] = HB_f_empty;
             s2pro[lev] = HB_f_empty;
             s3pro[lev] = HB_f_empty;
             s4pro[lev] = HB_f_empty;
             ppro[lev] = HB_f_empty;
             tpro[lev] = HB_f_empty;
             spro[lev] = HB_f_empty;
        }
   } /* end for lev */

    if (deepest > 0.0)
         *bdepth_ptr = deepest;   /* this value gets returned:  a true bottom depth  */
    else 
         *bdepth_ptr = bd;
	       
   if (deepest < 0.) {        /* is bottom missing? */
       s0pro[bindx] = HB_f_empty;
       s1pro[bindx] = HB_f_empty;
       s2pro[bindx] = HB_f_empty;
       s3pro[bindx] = HB_f_empty;
       s4pro[bindx] = HB_f_empty;
       ppro[bindx] = HB_f_empty;
       tpro[bindx] = HB_f_empty;
       spro[bindx] = HB_f_empty;
       
       return;       
   }

  /* flag values below the seafloor depth with HB_f_mask */

   for (lev = 0; lev < bindx; ++lev) {
        if (std_depth[lev] >= *bdepth_ptr) {
             s0pro[lev] = HB_f_mask;
             s1pro[lev] = HB_f_mask;
             s2pro[lev] = HB_f_mask;
             s3pro[lev] = HB_f_mask;
             s4pro[lev] = HB_f_mask;
             ppro[lev] = HB_f_mask;
             tpro[lev] = HB_f_mask;
             spro[lev] = HB_f_mask;
        }
   } /* end for lev */
    
 /* All profiles should now contain either a value, HB_f_empty or HB_f_mask */
         
   return;
    
} /* end get_sig_profiles() */
/*******************************************/

/****************************************************/
double find_bottom_isopycnal(int xrad, int yrad, int ncols, double **sigptr, float *bdwork, float bdepth, double sig_fg)
/*
   Search in 8 directions to find profiles with data near this
   depth .  Suspend search in a direction if run into a masked value
   or if xrad/yrad is reached.  Return weighted average of first-guess and 
   up to 8 closest sigma values (weighted by distance from central gridnode).  
   Return HB_f_empty if no value could be computed.
*/

{
   int i, row, col, stop, sq, xradius, yradius, zindex, wflag;
   double sig, ssum, wght, wsum;
   float bd;

   
   
   ssum = 0.0;
   wsum = 0.0;
   wflag = 0;
   zindex = NSTDLEVS-1;

   
   /* search west */
   
   row = yrad;  /* by definition, row/col corresponding to gridnode being fit */
   col = xrad;
   xradius = 0;
   wght = 2.0; 
   stop = 0; 
   while (++xradius <= xrad && !stop)  {
      --col;
      wght = wght * 0.5;
      sq = row * ncols + col;
      bd = bdwork[sq];
      
      sig = sigptr[sq][zindex];  /* first assume neighbor is close to same bottom depth */

      if (ABS(bd - bdepth) > 151.) {
         if (bd < bdepth) 
	    sig = (double) HB_f_mask;  /* neighbor is too shallow */
	    
	 else {
	    /* find 2 stddepths that bracket bdepth and interpolate sig value */
	    i = 0;
	    while (i < zindex && (bdepth >= std_depth[i]))
	      ++i;
	      
	    if (i == 0)
	      	sig = (double) HB_f_empty;  /* bdepth < first std_depth  */
		
            else {   /* i and i-1 bracket the value ... */

                if (sigptr[sq][i-1] < 0.0 && sigptr[sq][i] < 0.0) {
		   sig = (double) HB_f_empty;  /* neighbor has no value  */
		}
		else { 
		  if (sigptr[sq][i-1] < 0.0)  /* upper level empty */ 
		     sig = sigptr[sq][i];
		  else if (sigptr[sq][i] < 0.0 )  /* lower level empty */
		     sig = sigptr[sq][i-1];
		  else 	                          /* interpolate */
	             sig = sigptr[sq][i-1] + (sigptr[sq][i] - sigptr[sq][i-1]) * (bdepth - std_depth[i-1]) / (std_depth[i] - std_depth[i-1]);
		  
		}
	    }
	 }
      }
      
      if (sig > 0 ) {  /* either a sig value or a mask */
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  +=  wght;
	    ++wflag;
	}
        stop = 1;
      }
   }
   
   /* search east */
   row = yrad;
   col = xrad;
   xradius = 0;
   wght = 2.0; 
   stop = 0; 
   while (++xradius <= xrad && !stop)  {
      ++col;
      wght = wght * 0.5;
      sq = row * ncols + col;
      bd = bdwork[sq];
      
      sig = sigptr[sq][zindex];  /* first assume neighbor is close to same bottom depth */
      
      if (ABS(bd - bdepth) > 201.) {

         if (bd < bdepth) 
	    sig = (double) HB_f_mask;  /* neighbor is too shallow */
	    
	 else {
	    /* find 2 stddepths that bracket bdepth and interpolate sig value */
	    i = 0;
	    while (i < zindex && (bdepth >= std_depth[i]))
	      ++i;
	      
	    if (i == 0)
	      	sig = (double) HB_f_empty;  /* bdepth < first std_depth  */
		
            else {   /* i and i-1 bracket the value ... */

                if (sigptr[sq][i-1] < 0.0 && sigptr[sq][i] < 0.0) {
		   sig = (double) HB_f_empty;  /* neighbor has no value  */
		}
		else { 
		  if (sigptr[sq][i-1] < 0.0)  /* upper level empty */ 
		     sig = sigptr[sq][i];
		  else if (sigptr[sq][i] < 0.0 )  /* lower level empty */
		     sig = sigptr[sq][i-1];
		  else 	                          /* interpolate */
	             sig = sigptr[sq][i-1] + (sigptr[sq][i] - sigptr[sq][i-1]) * (bdepth - std_depth[i-1]) / (std_depth[i] - std_depth[i-1]);
		  
		}
	    }
	 }
      }
      
      if (sig > 0 ) {
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  +=  wght;
	    ++wflag;
	}
        stop = 1;
      }
   }
  
   /* search north */
   
   row = yrad;
   col = xrad;
   yradius = 0;
   wght = 2.0; 
   stop = 0; 
   while (++yradius <= yrad && !stop)  {
      ++row;
      wght = wght * 0.5;
      sq = row * ncols + col;
      bd = bdwork[sq];
      sig = sigptr[sq][zindex];  /* first assume neighbor is close to same bottom depth */

      if (ABS(bd - bdepth) > 201.) {
         if (bd < bdepth) 
	    sig = (double) HB_f_mask;  /* neighbor is too shallow */
	    
	 else {
	    /* find 2 stddepths that bracket bdepth and interpolate sig value */
	    i = 0;
	    while (i < zindex && (bdepth >= std_depth[i]))
	      ++i;
	      
	    if (i == 0)
	      	sig = (double) HB_f_empty;  /* bdepth < first std_depth  */
		
            else {   /* i and i-1 bracket the value ... */

                if (sigptr[sq][i-1] < 0.0 && sigptr[sq][i] < 0.0) {
		   sig = (double) HB_f_empty;  /* neighbor has no value  */
		}
		else { 
		  if (sigptr[sq][i-1] < 0.0)  /* upper level empty */ 
		     sig = sigptr[sq][i];
		  else if (sigptr[sq][i] < 0.0 )  /* lower level empty */
		     sig = sigptr[sq][i-1];
		  else 	                          /* interpolate */
	             sig = sigptr[sq][i-1] + (sigptr[sq][i] - sigptr[sq][i-1]) * (bdepth - std_depth[i-1]) / (std_depth[i] - std_depth[i-1]);
		  
		}
	    }
	 }
      }
      
      
      if (sig > 0 ) {
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  += wght;
	    ++wflag;
	}
        stop = 1;
      }
   }
  
  
   /* search south */
   
   row = yrad;
   col = xrad;
   yradius = 0;
   wght = 2.0; 
   stop = 0; 
   while (++yradius <= yrad && !stop)  {
      --row;
      wght = wght * 0.5;
      sq = row * ncols + col;
      bd = bdwork[sq];
      sig = sigptr[sq][zindex];  /* first assume neighbor is close to same bottom depth */

      if (ABS(bd - bdepth) > 201.) {
         if (bd < bdepth) 
	    sig = (double) HB_f_mask;  /* neighbor is too shallow */
	    
	 else {
	    /* find 2 stddepths that bracket bdepth and interpolate sig value */
	    i = 0;
	    while (i < zindex && (bdepth >= std_depth[i]))
	      ++i;
	      
	    if (i == 0)
	      	sig = (double) HB_f_empty;  /* bdepth < first std_depth  */
		
            else {   /* i and i-1 bracket the value ... */

                if (sigptr[sq][i-1] < 0.0 && sigptr[sq][i] < 0.0) {
		   sig = (double) HB_f_empty;  /* neighbor has no value  */
		}
		else { 
		  if (sigptr[sq][i-1] < 0.0)  /* upper level empty */ 
		     sig = sigptr[sq][i];
		  else if (sigptr[sq][i] < 0.0 )  /* lower level empty */
		     sig = sigptr[sq][i-1];
		  else 	                          /* interpolate */
	             sig = sigptr[sq][i-1] + (sigptr[sq][i] - sigptr[sq][i-1]) * (bdepth - std_depth[i-1]) / (std_depth[i] - std_depth[i-1]);
		  
		}
	    }
	 }
      }
      
      
      if (sig > 0 ) {
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  +=  wght;
	    ++wflag;
	}
        stop = 1;
      }
   }

   /* search northeast */
   row = yrad;
   col = xrad;
   xradius = 0;
   yradius = 0;
   wght = 1.5; 
   stop = 0; 
   while ((++xradius < xrad) && (++yradius < yrad) && !stop)  {
      ++col;
      ++row;
      wght = wght * 0.5;
      sq = row * ncols + col;
      bd = bdwork[sq];
      
      sig = sigptr[sq][zindex];  /* first assume neighbor is close to same bottom depth */
      
      if (ABS(bd - bdepth) > 201.) {

         if (bd < bdepth) 
	    sig = (double) HB_f_mask;  /* neighbor is too shallow */
	    
	 else {
	    /* find 2 stddepths that bracket bdepth and interpolate sig value */
	    i = 0;
	    while (i < zindex && (bdepth >= std_depth[i]))
	      ++i;
	      
	    if (i == 0)
	      	sig = (double) HB_f_empty;  /* bdepth < first std_depth  */
		
            else {   /* i and i-1 bracket the value ... */

                if (sigptr[sq][i-1] < 0.0 && sigptr[sq][i] < 0.0) {
		   sig = (double) HB_f_empty;  /* neighbor has no value  */
		}
		else { 
		  if (sigptr[sq][i-1] < 0.0)  /* upper level empty */ 
		     sig = sigptr[sq][i];
		  else if (sigptr[sq][i] < 0.0 )  /* lower level empty */
		     sig = sigptr[sq][i-1];
		  else 	                          /* interpolate */
	             sig = sigptr[sq][i-1] + (sigptr[sq][i] - sigptr[sq][i-1]) * (bdepth - std_depth[i-1]) / (std_depth[i] - std_depth[i-1]);
		  
		}
	    }
	 }
      }
      
      if (sig > 0 ) {
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  +=  wght;
	    ++wflag;
	}
        stop = 1;
      }
   }
  
    /* search southeast */
   row = yrad;
   col = xrad;
   xradius = 0;
   yradius = 0;
   wght = 1.5; 
   stop = 0; 
   while ((++xradius < xrad) && (++yradius < yrad) && !stop)  {
      ++col;
      --row;
      wght = wght * 0.5;
      sq = row * ncols + col;
      bd = bdwork[sq];
      
      sig = sigptr[sq][zindex];  /* first assume neighbor is close to same bottom depth */
      
      if (ABS(bd - bdepth) > 201.) {

         if (bd < bdepth) 
	    sig = (double) HB_f_mask;  /* neighbor is too shallow */
	    
	 else {
	    /* find 2 stddepths that bracket bdepth and interpolate sig value */
	    i = 0;
	    while (i < zindex && (bdepth >= std_depth[i]))
	      ++i;
	      
	    if (i == 0)
	      	sig = (double) HB_f_empty;  /* bdepth < first std_depth  */
		
            else {   /* i and i-1 bracket the value ... */

                if (sigptr[sq][i-1] < 0.0 && sigptr[sq][i] < 0.0) {
		   sig = (double) HB_f_empty;  /* neighbor has no value  */
		}
		else { 
		  if (sigptr[sq][i-1] < 0.0)  /* upper level empty */ 
		     sig = sigptr[sq][i];
		  else if (sigptr[sq][i] < 0.0 )  /* lower level empty */
		     sig = sigptr[sq][i-1];
		  else 	                          /* interpolate */
	             sig = sigptr[sq][i-1] + (sigptr[sq][i] - sigptr[sq][i-1]) * (bdepth - std_depth[i-1]) / (std_depth[i] - std_depth[i-1]);
		  
		}
	    }
	 }
      }
      
      if (sig > 0 ) {
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  +=  wght;
	    ++wflag;
	}
        stop = 1;
      }
   }
  
   /* search southwest */
   row = yrad;
   col = xrad;
   xradius = 0;
   yradius = 0;
   wght = 1.5; 
   stop = 0; 
   while ((++xradius < xrad) && (++yradius < yrad) && !stop)  {
      --col;
      --row;
      wght = wght * 0.5;
      sq = row * ncols + col;
      bd = bdwork[sq];
      
      sig = sigptr[sq][zindex];  /* first assume neighbor is close to same bottom depth */
      
      if (ABS(bd - bdepth) > 201.) {

         if (bd < bdepth) 
	    sig = (double) HB_f_mask;  /* neighbor is too shallow */
	    
	 else {
	    /* find 2 stddepths that bracket bdepth and interpolate sig value */
	    i = 0;
	    while (i < zindex && (bdepth >= std_depth[i]))
	      ++i;
	      
	    if (i == 0)
	      	sig = (double) HB_f_empty;  /* bdepth < first std_depth  */
		
            else {   /* i and i-1 bracket the value ... */

                if (sigptr[sq][i-1] < 0.0 && sigptr[sq][i] < 0.0) {
		   sig = (double) HB_f_empty;  /* neighbor has no value  */
		}
		else { 
		  if (sigptr[sq][i-1] < 0.0)  /* upper level empty */ 
		     sig = sigptr[sq][i];
		  else if (sigptr[sq][i] < 0.0 )  /* lower level empty */
		     sig = sigptr[sq][i-1];
		  else 	                          /* interpolate */
	             sig = sigptr[sq][i-1] + (sigptr[sq][i] - sigptr[sq][i-1]) * (bdepth - std_depth[i-1]) / (std_depth[i] - std_depth[i-1]);
		  
		}
	    }
	 }
      }
      
      if (sig > 0 ) {
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  +=  wght;
	    ++wflag;
	}
        stop = 1;
      }
   }
  
    /* search northwest */
   row = yrad;
   col = xrad;
   xradius = 0;
   yradius = 0;
   wght = 1.5; 
   stop = 0; 
   while ((++xradius < xrad) && (++yradius < yrad) && !stop)  {
      --col;
      ++row;
      wght = wght * 0.5;
      sq = row * ncols + col;
      bd = bdwork[sq];
      
      sig = sigptr[sq][zindex];  /* first assume neighbor is close to same bottom depth */
      
      if (ABS(bd - bdepth) > 201.) {

         if (bd < bdepth) 
	    sig = (double) HB_f_mask;  /* neighbor is too shallow */
	    
	 else {
	    /* find 2 stddepths that bracket bdepth and interpolate sig value */
	    i = 0;
	    while (i < zindex && (bdepth >= std_depth[i]))
	      ++i;
	      
	    if (i == 0)
	      	sig = (double) HB_f_empty;  /* bdepth < first std_depth  */
		
            else {   /* i and i-1 bracket the value ... */

                if (sigptr[sq][i-1] < 0.0 && sigptr[sq][i] < 0.0) {
		   sig = (double) HB_f_empty;  /* neighbor has no value  */
		}
		else { 
		  if (sigptr[sq][i-1] < 0.0)  /* upper level empty */ 
		     sig = sigptr[sq][i];
		  else if (sigptr[sq][i] < 0.0 )  /* lower level empty */
		     sig = sigptr[sq][i-1];
		  else 	                          /* interpolate */
	             sig = sigptr[sq][i-1] + (sigptr[sq][i] - sigptr[sq][i-1]) * (bdepth - std_depth[i-1]) / (std_depth[i] - std_depth[i-1]);
		  
		}
	    }
	 }
      }
      
      if (sig > 0 ) {
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  +=  wght;
	    ++wflag;
	}
        stop = 1;
      }
   }
  
/* Compare to first guess isopycnal if it is an actual value */

   if (sig_fg > 0.0 && sig_fg < 100.0) {
            
      if (wsum > 0) {
            sig = ssum / wsum;
/*            fprintf(stderr, "Difference at bottom: sig - sig_fg = %.3lf \n", sig - sig_fg); */
     }
      
      
      ssum += sig_fg * 1.0;
      wsum += 1.0;
      ++wflag ;
   }
   
 
   if (!wflag) 
      return((double) HB_f_empty);
      
   return (ssum / wsum);
   
}  /* end find_bottom_isopycnal() */

/****************************************************/
double get_neighbors_anom( double isopyc, double depth_being_interpolated, double *xsurf,  int nsq, struct GRID_INFO *hptr, double **xwork2, double **sigwork2,  float *bdpth, double **xwork2_fg, double **sigwork2_fg,  float *bdpth_fg)

/*
   Traverse each gridnode in work arrays. The bdpth array has the same
   dimension as the work grids and x grid -- it contains the seafloor
   depth at each gridnode. Interpolate to find prop value 
   associated with specified isopycnal in first guess and xwork2 arrays.
   Store difference (xwork2 - first_guess) in xsurf
   If isopycnal runs into bottom or outcrops, mark xsurf[sq] with HB_MASK.
   If no anomaly value can be computed, mark with HB_EMPTY.
   Returns xvalue in first-guess field associated with middle square.
*/

{
   double  xval, xbase;
   int  sq, sqmid;

   sqmid = nsq / 2;    /* by definition, the square to be fit -- just label it empty */

   if (depth_being_interpolated < 0) {
       fprintf(stderr,"FATAL ERROR:  depth_being_interpolated is negative in get_neighbors_anom()\n");  
    } 
    
   for (sq = 0; sq < nsq; ++sq) {
      if (sq == sqmid) {
         xsurf[sq] = (double) HB_f_empty;
	 xbase = get_xval(xwork2_fg[sq], sigwork2_fg[sq], isopyc, bdpth_fg[sq], depth_being_interpolated);
      }
      else {
         /* get xval for this square in grid being fitted */
         xsurf[sq] = get_xval(xwork2[sq], sigwork2[sq], isopyc, bdpth[sq], depth_being_interpolated);
	 
         if (!  (is_flagged((float)xsurf[sq], HB_f_empty) || is_flagged((float)xsurf[sq], HB_f_mask ) )) {
	         /* get xval in first-guess field */
           xval =  get_xval(xwork2_fg[sq], sigwork2_fg[sq], isopyc, bdpth_fg[sq], depth_being_interpolated);
	  
	   if ( (is_flagged((float)xval, HB_f_empty) || is_flagged((float)xval, HB_f_mask ) ) )
	     xsurf[sq] = xval;
	   else 
	     xsurf[sq] -= xval;   /* store difference */
	   
         }
      }
   }  /* end for sq */
   return (xbase);
}  /* end get_neighbors_anom() */
/****************************************************/
void get_neighbors( double isopyc, double depth_being_interpolated, double *xsurf,  int nsq, double **xwork2, double **sigwork2,  float *bdpth)
/*
   Traverse each gridnode in work arrays. The bdpth array has the same
   dimension as the work grids and x grid -- it contains the seafloor
   depth at each gridnode. Interpolate to find prop value 
   associated with specified isopycnal in xwork2 arrays -- Store in xsurf.
   If isopycnal runs into bottom or outcrops, mark xsurf[sq] with HB_MASK.
   If no  value can be computed, mark with HB_EMPTY.
*/
{
   double  xval;
   int  sq;

   if (depth_being_interpolated < 0) {
       fprintf(stderr,"FATAL ERROR:  depth_being_interpolated is negative in get_neighbors()\n");  
    } 
    
   for (sq = 0; sq < nsq; ++sq) {
          xsurf[sq] = get_xval(xwork2[sq], sigwork2[sq], isopyc, bdpth[sq], depth_being_interpolated);
   }  /* end for sq */
   return;
}  /* end get_neighbors() */
/***************************************************/
void get_weights(double *weights, float *bdepth, double L, struct GRID_INFO *hptr)
  /* sets weights at each gridpoint in weights array
    based on distance from central gridnode 
	   e^-([dist/L]^2) 
	   
	 ****  plus added weight for bathymetry gradient (if < 100m/deg):
	     1.0 - bathgrad/100 
 */
{
   int row, col, sq, imid, jmid;
  double dx, dy, dist;
  double *bath, *x, xav, yav, bathgrad;
  float clat, clon, lat, lon;
  
 /* determine central node's lat/lon */ 
  imid = hptr->nx / 2;
  jmid = hptr->ny / 2;
  ij2xy(hptr, imid, jmid, &clon, &clat);
 
 /* compute distance and weight at each point in array */ 
  for (row = 0; row < hptr->ny; ++row) {
     for (col = 0; col < hptr->nx; ++col) {
        sq = row * hptr->nx + col;
        ij2xy(hptr, col, row, &lon, &lat);
        dx = (lon - clon)*RADperDEG * cos(RADperDEG*.5*(lat+clat)) * EarthRadius;
        dy = (lat - clat) * RADperDEG * EarthRadius;
        dist = sqrt(dx*dx + dy*dy);
	dist = dist / L;
	weights[sq] = exp(-(dist*dist));
     }
   }
   
   return;

}  /* end get_weights() */
/****************************************************/

double get_xval(double *xptr, double *sigptr, double isopyc, float btmd, double depth_being_interpolated)
/* interpolates to find value of x in array xptr corresponding to isopyc in array sigptr.  
        btmd :  bottom depth, depth values are stored in global variable std_depth. 
	depth_being_interpolated:  std_depth of level being interpolated.
    Returns value of x, HB_f_empty or HB_f_mask.
	 */
{
   double *d, *x, *sig, xval;
   double diff, lastdiff, zval, best_zval;
   double emptyflag;
   int n, nz, lev;
   int no_bottom, j, datagap;

/* Allocate memory to hold profiles */

   d = (double *) get_memory((void *)NULL, (size_t) NSTDLEVS, sizeof(double));
   x = (double *) get_memory((void *)NULL, (size_t) NSTDLEVS, sizeof(double));
   sig = (double *) get_memory((void *)NULL, (size_t) NSTDLEVS, sizeof(double));

   nz = NSTDLEVS - 1;
   emptyflag  = HB_f_empty + 10;
   
      /* load  work profiles into local profile array and weed out
         missing values for interpolation function */
      n = 0;
      for (lev = 0; lev < nz; ++lev) {
        if (sigptr[lev] > emptyflag && sigptr[lev] < 100000.0 && xptr[lev] > emptyflag) {
	  d[n] = std_depth[lev];
	  x[n] = xptr[lev];
	  sig[n] = sigptr[lev];
	  
	  ++n;
	}
      }
      
      /* check seafloor for valid observation */
       
      no_bottom = 1;     
      if (sigptr[nz] > emptyflag && sigptr[nz] < 100000.0 && xptr[nz] > emptyflag) {
          d[n] = btmd;  /* insert seafloor depth for this square */
	  x[n] = xptr[nz];
	  sig[n] = sigptr[nz];
	  ++n;
	  no_bottom = 0;
      }
      
      if (n == 0) {   /* no values at this grid node */
             free(d);
             free(x);
             free(sig);
         if (is_flagged(btmd, HB_f_mask))  {  /* node is masked */
	   return ((double) HB_f_mask);
	 }
	 return ((double) HB_f_empty);
      }
      
      else  {   /* check for all crossings and use closest to this stddepth */
           j = 0;
	   lastdiff = (double)  HB_f_mask;  /* initialize with very large number */ 
	   best_zval =  (double) HB_f_empty;  /* initialize with very negative number */ 
	   do {
               zval = hb_linterp(isopyc, &sig[j], &d[j], n-j);
	       if (zval > -9990.) {
	        
	           if ( (diff = ABS(zval - depth_being_interpolated)) < lastdiff ) {
		            best_zval = zval;
			    lastdiff = diff;
		    }
		    
		    while ((d[j] <= zval) && (j < n))
	               ++j;
	       }
           } while (zval > -99990. && (j < n));
	 
	   if (best_zval >= 0 ) {   /* found an appropriate crossing */
	   
	       /* if this is at the top of a pycnostad, get as close as possible 
	       to depth being interpolated */
	       
	       if (lastdiff > 300) {
	          zval = best_zval;
	          diff =  depth_being_interpolated - best_zval;
		  if (diff > 0)  {   /* zval is located > 300 meters above the level being interpolated for*/
		     j = 0;
		    while ((d[j] <= zval) && (j < n))
	               ++j;
		     /* now j points to level just below best_zval */
		     
		     while ((j < n) && (d[j] < depth_being_interpolated) && ABS((sig[j] - isopyc) < 0.005) )
		         ++j;
			 
			 /*now j points to level closest to depth_being_interpolated that is 
	           within 0.01 sigma units of the isopycnal value being sought */
		       
		       if (j == n)  /* should never happen, but check anyway */
		             --j;
			 
		       best_zval = d[j];
		  }/* end if diff > 0 */
		  
		/*  else  best_zval is below depth being interpolated  and so can't get any closer  
		(assumes that density is monotonically increasing) */
		  
	       } /* end if lastdiff > 300 */
	       
	       zval = best_zval;
	       
                /* get xvalue and check for datagaps around this depth...*/
	    
	           xval = hb_linterp(zval, d, x, n);
	    	    
                   j = 0;
	           while (d[j] <= zval && j < n)
	            ++j ;
		 
	           if  (j == n  && d[j] == zval)
	                datagap = 0;
		    else if (j == 0)
		        datagap = 0;
	           else if ((d[j-1] == zval) || (d[j] == zval))
	               datagap = 0;
	           else if (zval < 1001)
	               datagap = (d[j] - d[j-1]) > GAP_SHALLOW;
                   else
                       datagap = (d[j] - d[j-1]) > GAP_DEEP;
		
                   free(d);
                   free(x);
                   free(sig);
		
	           if (!datagap) {
	              return(xval);
	           }
		   return((double) HB_f_empty);
		
	  }
	 else {  /* no crossing -- does it outcrop?*/
	    
             free(d);
             free(x);
             free(sig);
	     if (!no_bottom && (isopyc > sig[n-1])) {
	         return((double)HB_f_mask);   /* isopycnal runs into sea floor */
	       }
	       
	    return((double) HB_f_empty);
	 }
      } /* end else check all crossings */
      
} /* end get_xval() */


/****************************************************/

void x_to_stddepth(double *xin, short *xcnt, float *xout, double *din, double *stdd, float btmd, int nlevs)
/* interpolates x,depth arrays (xin, din) onto depth levels in stdd.  
        xcnt :  < 0 for values that were fitted, > 0 for levels that did not need fitting,
	         = 0 for levels that are masked or missing. 
        btmd :  bottom depth. 
    Uses xcnt to determine which values need to be interpolated:  if xcnt >= 0 the xin
    value is copied to xout.
    Returns array  xout ready to be written to netcdf file: depths below btmd are assigned HB_f_mask.
    Missing values are assigned HB_f_empty.
*/
{
   double *dtmp, *xtmp,  xval;
   double emptyflag, prevdepth;
   int n, nz, lev, j;
    
/*  case where there is no ocean */
   if ( is_flagged(btmd, HB_f_mask) ) {
      for (lev = 0; lev < nlevs; ++lev)  {
          if (xcnt[lev] != 0) {
	    fprintf(stderr," >> SOFTWARE BUG:  masked square has non-zero count for fitted values.\n");
	  }
          xout[lev] =  (float) xin[lev];
      }	  
      return;
   }
   
/* Allocate memory to hold profiles */

   dtmp = (double *) get_memory((void *)NULL, (size_t) nlevs, sizeof(double));
   xtmp = (double *) get_memory((void *)NULL, (size_t) nlevs, sizeof(double));

   nz = nlevs - 1;    /* index to deepest level */
   emptyflag  = HB_f_empty + 10;
   
      /* load  work profiles into local profile array and weed out
         missing values for interpolation function */
      prevdepth = -10.0; 
      n = 0;
      for (lev = 0; lev < nz; ++lev) {
        if (xcnt[lev] != 0 ) {
	  dtmp[n] = din[lev];
	  if (dtmp[n] < prevdepth) {
/*	      fprintf(stderr," >> MESSAGE from x_to_stddepth(): depth array not monotonically increasing. prevdepth: %.1lf   thisdepth: %.1lf.\n",  prevdepth, dtmp[n]); */
	  }
	  else {
	     prevdepth = dtmp[n];
	     xtmp[n] = xin[lev];
	     ++n;
	  }
	}
      }
      if (xcnt[nz] != 0) {     /* add bottom level */
         dtmp[n] =  din[nz];
	  if (dtmp[n] < prevdepth) {
/*	      fprintf(stderr," >> MESSAGE from x_to_stddepth(): bottom depth not monotonically increasing.  prevdepth: %.1lf   thisdepth: %.1lf\n", prevdepth, dtmp[n]); */
	  }
	 xtmp[n] = xin[nz];
	 ++n;
      }
           
      if (n < 2) {  /* not enough xvalues to interpolate, just copy input array to output array */
          for (lev = 0; lev <= nz; ++lev)  
	      xout[lev] =  (float) xin[lev];
	      
	  free(dtmp);
	  free(xtmp);
	  return;
      } 
      
      /* interpolate fitted values (where xcnt < 0) onto stddepths. */
      
      for (lev = 0; lev < nz; ++lev) {
      
         if (xcnt[lev] >= 0)  /* level was not interpolated */
	      xout[lev] = (float) xin[lev];
	 else {
	      xval = hb_linterp(stdd[lev], dtmp, xtmp, n);
	   
	      xout[lev] = (float) xval;
	      if (xval < -9999.0){
/*	          fprintf(stderr," >> SOFTWARE BUG in x_to_stddepth():  level was fitted but cannot be interpolated.   sdepth= %.1lf lev=%2d\n", stdd[lev], lev); */
	      
	          xout[lev] = HB_f_empty;
	      }
	 }
       }
       
       xout[nz] = (float) xin[nz];   /* explicity set, don't interpolate, deepest value */
      free(dtmp);
      free(xtmp);
      return;

} /* end x_to_stddepth */

