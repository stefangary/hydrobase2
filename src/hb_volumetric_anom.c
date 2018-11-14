/*  hb_volumetric_anom.c

................................................................................
                          *******  HydroBase2 *******
               ...................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
                             2004
................................................................................
____________________________________________________________________________
  USAGE:  

 hb_volumetric_anom filename(_roots) -1mne/value[t] -2mne/value[b]  -B<w/e/s/n> -Ixinc[/yinc] -Mjan<filename> -Mfeb<filename>   -Mmar<filename> -Mapr<filename> -Mmay<filename> -Mjun<filename>   -Mjul<filename> -Maug<filename>  -Msep<filename>   -Moct<filename> -Mnov<filename> -Mdec<filename> -S<std_depth_filename> [-D<dirname>] [-E<file_extent>] [-Z<zmin/zmax>]

  list of filenames MUST be first argument

 -1 : property and value at upper level;
 -2 : property and value at lower level;
 -B : grid bounds (west/east/south/north)
 -I : grid increments (xgrid, ygrid)
 -M : profiles for each month to 

 -D : specifies directory for input files (default is current dir)  
      ex: -D/home/pete/data

 -E : specifies input_file extent (default is no extent) ex: -E.sum
 -Z : specifies depth limits
____________________________________________________________________________
Computes volumetric heat content anomaly and fresh water anomaly for a vertical layer defined by two surfaces over a given region.  Deviations of Theta and Salinity for measured profiles are computed relative to a monthly average profile (i.e. the seasonal cycle is removed).  These monthly averages are provided in filenames specified by -M arguments -- one for each month.

Foreach gridpoint, deltaT and deltaS are accumulated over the specified layer --weighting in the vertical is by dz. Rectangular area surrounding each gridpoint is dx*dy (xinc/yinc) at each profile's latitude.

Outputs  lon, lat, delta_theta, delta_salt, avg_theta(per unit volume), avg_salt, avg_rho, avg_Cp, dz, dxdy, heat content anomaly (delta-Q J/m^2:  delta_theta*avg_rho*avg_Cp*dz) and fresh water anomaly (delta-F in meters/unit area: dz*delta_salt/avg_salt).
The sums of property content at all gridpoints are accumulated (weighted by area of each gridpoint) and on the last lines are listed the following:
   Total area (km^2)
   Avg thickness of layer (area weighted) (km)
   Total heat content anomaly (J/m^2) 
   Total fresh water anomaly (m)
   Avg Theta/unit volume
   Avg delta_theta / unit volume
   Avg Salinity/unit volume 
   Avg delta_salt /unit volume
   Avg Density /unit volume
   Avg Heat Content/unit vol
   
____________________________________________________________________________
*/


#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include "hydrobase.h"
#include "hb_gamma.h"
#include "hb_paths.h"

/* input file pathnames */

#define    EXTENT   ""
#define    DIR      ""
#define    MAXOBS  6000
#define    NMONTHS  13

#define  PI     3.141592654
#define  RADperDEG 0.017453292    /* pi/180 */
#define  EarthRadius  3437.747    /* in nm */
#define  KMperNM  1.853248652

      
				    
  /*  define structures to store profile information at each gridnode */

struct PROFILE {
        int nobs;
	double *d, *t, *s;         
};				    

struct GRIDNODE {
        int    count;
	double dxdy;  /* area (km) (function of lon,lat) at gridnode */        
	double dTdz;          
	double dSdz;          
	double rho_dz;          
	double cp_dz;          
	double dz;   
	double dz_diff; 
	double dz_perfect;      
	double Tdz;          
	double Sdz;          
};

struct PROFILE **stdprofile;				    
struct GRIDNODE *grid;

enum MONTH { CLIM, JAN, FEB, MAR, APR, MAY, JUN, JUL, AUG, SEP, OCT, NOV, DEC};
int  mfile[NMONTHS];   /* files containing monthly and climatological average profiles */

   /* global variables to store station ... */

struct HYDRO_DATA station;
struct HYDRO_HDR hdr;
double x1, x2;
double pref1, pref2;    /* ref pressures for top/bottom of layer */
int prop1, prop2;
int nrows, ncols, nsq;  /* for grids */
int xdateline;

struct GAMMA_NC ginfo;   /* used for neutral density */

int warnflag;
int iflag, zflag;
int check_bottom_outcrop;
int check_top_outcrop;
double zmin, zmax;
double xmin, xmax, ymin, ymax;
double xinc, yinc;

   /* prototypes for locally defined functions */

void print_usage(char *);
int open_monthly_file(char *);
int get_monthly_data(int, struct PROFILE *);
double get_top_of_layer(int);
double get_bottom_of_layer(double, int, double * );
void get_hydro_data(int);
void get_weighted_ts_anom(double, double, double *, double *, double *, struct PROFILE *);
double get_weighted_avg(double, double, double *,double *, double *, double *, double *, double *);
void compute_prop(int, double **, double *, double *, double *, int, double, double);
int vector(double, double, double, double,  double *);
void output_data(FILE *);

int main (int argc, char **argv)
{
   short   layer1, layer2;
   int     nprops, i, n, index;
   int     curfile = 1, nfiles = 0;
   int     bflag, sflag, mflag; 
   int     print_msg = 1;       /* print message in file open routines*/
   int     error, prompt;
   char    *extent, *dir, *st;
   char    id[3];
   int     infile;
   

/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
 
/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    layer1 = layer2 = 0;
    error = 0;
    nprops = 2;
    zmin = 0.0;
    zmax = 10000.0;
    infile = STDIN;
    iflag = zflag = 0;
    bflag = sflag = mflag = 0;
    warnflag = 0;           /* set to 1 after warning is printed */
    check_bottom_outcrop = 0;
    check_top_outcrop = 0;
     xdateline = 0;
   
    
    
/* initialize these ... */

   for (i = 0; i < MAXPROP; ++i) {
      station.observ[i] = (double *) NULL;
   }
   hdr.prop_id = (int *) NULL;

   for (i = 0; i < NMONTHS; ++i)    /* file descriptors */
      mfile[i] = -1;

/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case '1':
                        layer1 = 1;
                        error = (sscanf(&argv[i][2],"%2s", id) == 1) ? 0 : 1;
                        if ((prop1 = get_prop_indx(id)) < 0) {
                           fprintf(stderr,"\n%.2s is not an appropriate property.", id);
                           ++error;
                        }
                        
                        if (prop1 != (int)PR  && prop1 != (int)DE  ) {
                           fprintf(stderr,"\nZlayer must be pr or de\n%.2s is not an appropriate property.", id);
                           ++error;
                        }
                        else {
			 st = &argv[i][4];
			 if (*st == '/')
			   ++st;
                         if (sscanf(st,"%lf", &x1) != 1)
                            ++error;
                        }
			n = strlen(argv[i]);
			if (argv[i][n-1] == 't')   /* check for optional t */
			    check_top_outcrop = 1;
                        break;

               case '2':
                        layer2 = 1;
                        error = (sscanf(&argv[i][2],"%2s", id) == 1) ? 0 : 1;
                        if ((prop2 = get_prop_indx(id)) < 0) {
                           fprintf(stderr,"\n%.2s is not an appropriate property.", id);
                           ++error;
                        }
                        
                        if (prop2 != (int)PR  && prop2 != (int)DE  ) {
                           fprintf(stderr,"\nZlayer must be pr or de\n%.2s is not an appropriate property.", id);
                           ++error;
                        }
                        else {
			 st = &argv[i][4];
			 if (*st == '/')
			   ++st;
                         if (sscanf(st,"%lf", &x2) != 1)
                            ++error;
                        }
			n = strlen(argv[i]);
			if (argv[i][n-1] == 'b')   /* check for optional b */
			    check_bottom_outcrop = 1;
                        break;
               case 'D':                   /* get input dir */
                        dir = &argv[i][2];
                        break;
               case 'E':                    /* get file extent */
                        extent = &argv[i][2];
                        break;
               case 'B':
                        bflag = 1;
                        error = (sscanf(&argv[i][2],"%lf/%lf/%lf/%lf", &xmin, &xmax, &ymin, &ymax) == 4) ? 0 : 1;
                        
                        if (xmin > 0 && xmax < 0)
                           xmax += 360.;
                        if (xmax > 180.)
                           xdateline = 1;
			if (ymin > ymax) {
			   fprintf(stderr,"Southern grid boundary exceeds northern boundary.\n");
			   error = 1;
			}
			if (xmin > xmax) {
			   fprintf(stderr,"Western grid boundary exceeds eastern boundary.\n");
			   error = 1;
			}
                        break;

               case 'I':
                        iflag = 1;
                        error = (sscanf(&argv[i][2],"%lf", &xinc) == 1) ? 0 : 1;
                        yinc = xinc;
			
                        st = &argv[i][2];
                        while (*(st++) != '\0') {
                            if (*st == '/') {
                              ++st;
                              error = (sscanf(st,"%lf", &yinc) == 1) ? 0 : 1;
                              break;
                            }
                        }
                        break;
			
	       case 'M':
	               error = open_monthly_file(&argv[i][2]);
	               if (!error) 
		          ++mflag;
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

       else  {
           ++nfiles;
       }
   }  /* end for */

   if (!layer1 || !layer2 ) {
       fprintf(stderr,"\nYou must define the layer with -1 and -2.\n");
       exit(1);
   }
   if (  !(iflag && bflag)) {
       fprintf(stderr,"\nYou must specify grid bounds and spacing with -B and -I \n");
        exit(1);
  }
   
   if (  !mflag) {
       fprintf(stderr,"\nYou must specify at least one file containing average profiles to compute delta-T and delta-S. Use a 3-char string to specify month (jan..dec) or (cli) for overall climatology to be used if a month is not  available.\n");
        exit(1);
  }
   if (zflag) {
       fprintf(stderr,"\nUsing depth limits %.2lf - %.2lf db", zmin, zmax);
   }
   
     
   
/* compute dimensions of grid and allocate space ...*/

   nrows = NINT( ceil(((ymax - ymin) / yinc)));
   ncols = NINT(ceil(((xmax - xmin) / xinc)));
   nsq = nrows * ncols;
   grid = (struct GRIDNODE *) calloc((size_t)nsq, sizeof(struct GRIDNODE));
   stdprofile = (struct PROFILE **) calloc((size_t)NMONTHS, sizeof(struct PROFILE *));
   
/* read in monthly average profiles at each gridpoint */

   for (i = 0; i < NMONTHS; ++i) {
      if (mfile[i] >= 0) {
          stdprofile[i] = (struct PROFILE *)calloc((size_t)nsq, sizeof(struct PROFILE));
          if (stdprofile[i] == NULL) {
            fprintf(stderr,"\n Unable to allocate memory for monthly profiles\n");
            exit(1);
          }
          get_monthly_data(mfile[i], stdprofile[i]);
	  close(mfile[i]);
      }
   }

/* loop for each input file */

   do {

    if (!nfiles) {
       fprintf(stderr,"\nExpecting input from stdin ... \n");
    }
    else {
       infile = open_hydro_file(dir, argv[curfile], extent, print_msg);
       if (infile < 0)
          goto NEXTFILE;
    }

            /* read each file completely */

    get_hydro_data(infile);

NEXTFILE:

         close(infile);

   } while (curfile++ < nfiles );
   

   output_data(stdout);
   

   fprintf(stderr,"\n\nEnd of %s.\n", argv[0]);
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\n%s computes anomalies of heat content and fresh water anomaly from depth-averaged temperature and salinity over a layer specified by two surfaces", program);
   fprintf(stderr,"\nInput can be either observed profiles or gridded profiles");
   fprintf(stderr,"\nOutput values are gridded.");
   
   fprintf(stderr,"\nDelta-Theta and Delta-Salt are determined relative to monthly average profiles that are read in from files (HydroBase station format -- can be derived from CDF gridded files using hb_cdf2asc).  These files are specified with -Mmonth$filename");     
   fprintf(stderr,"\ndTdz and dSdz are estimated for a series of vertical zbins and their weighted sum provides an estimate of the anomaly at each gridnode.");
    fprintf(stderr,"\nThe zlayer can be strictly adhered to (zlayer must exist in entirety) or relaxed to utilize whatever portion of the layer exists.");
   fprintf(stderr,"\nFor each gridpoint, outputs a line containing:");
   fprintf(stderr,"\nlon, lat, dtheta, dsalt, avg_theta, avg_salt, avg_rho, avg_Cp, dz(layer_thickness), dx*dy, delta-Q, delta-F ");
   fprintf(stderr,"\nThe sum of properties for all grid squares is accumulated and on the last lines are listed the following:");
   fprintf(stderr,"\nTotal area (km^2)");
   fprintf(stderr,"\nAvg thickness of layer (area weighted) (km)");
   fprintf(stderr,"\nAvg theta/unit volume (deg C)");
   fprintf(stderr,"\nAvg salinity/unit volume");
   fprintf(stderr,"\nAvg density (rho)/unit volume (kg/m^3");
   fprintf(stderr,"\nHeat content change (deltaQ) (J/m^2)");
   fprintf(stderr,"\nFreshwater anomaly per unit area (deltaF) (meters)");
   
   
   fprintf(stderr,"\n\nUsage:  %s filename_root(s)  -1<mne[pref]/value> -2<mne[pref]/value[b]> -Ixinc[/yinc] [-D<dirname>] [-E<file_extent>][-Z<zmin/zmax>]", program);

   fprintf(stderr,"\n\n  List of filenames MUST be first argument");
   fprintf(stderr,"\n  or input is expected to come from stdin...");
   fprintf(stderr,"\n    -1  : surface 1 z-property/value (property = pr or de");
   fprintf(stderr,"\n          Append 't' to start with the first observation level available");
   fprintf(stderr,"\n          in the event  it is deeper than value at the top of z-layer.   ");
   fprintf(stderr,"\n    -2  : surface 2 z-property/value (property = de or pr)");
   fprintf(stderr,"\n          Append b to integrate to the bottom if this zvalue");
   fprintf(stderr,"\n          is deeper than the bottom.  ");
   fprintf(stderr,"\n    -B  : specify grid bounds (w/e/n/s) in degrees");
   fprintf(stderr,"\n    -I  : specify grid spacing (xinc/yinc) in degrees to compute dx and dy");
   fprintf(stderr,"\n    -M  : filename for monthly average profiles: -M<month><filename>");
   fprintf(stderr,"\n           month is 3-char identifier (jan, feb, mar....dec,and  cli for climatology)");
   fprintf(stderr,"\n           Specify one or more month files.  If a needed month file is not supplied, the cli file will be used to calculate delta-T and delta-S");
   fprintf(stderr,"\n           At least one -M file must be specified.");
   fprintf(stderr,"\n           If neither the correct month nor a cli file is provided, no values for that month will be computed in the volumetric summation");
   
   fprintf(stderr,"\n   [-D] : specifies dirname for input files (default is current");
   fprintf(stderr,"\n          directory)  ex: -D../data/ ");
   fprintf(stderr,"\n   [-E] : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n   [-h] : help... prints this message.");
   fprintf(stderr,"\n\n");  
       
   return;
}


/*****************************************************************************/
int open_monthly_file(char *st)
  /* Opens the specified file for reading and stores file descriptor in the global array mfile. 
     EXITS with an error message if file open fails*/
{
   int index, error;
   
   error = 0;

   switch (*st) {
	case 'a':
	    switch (*(st+1)) {
		case 'p':
			 index = (int) APR;
			 break;
		case  'u':
			index = (int) AUG;
			 break;
	        default:
		        error = 1;
	    }
	    break;
	case 'c':
	    index = (int) CLIM;
	    break;
	case 'd':
	    index = (int) DEC;
	    break;
	case 'f':
	    index = (int) FEB;
	    break;
	case 'j':
	    if (*(st+1) == 'a') {
	    	index = (int) JAN;
	    }
	    else {
	       switch (*(st+2)) {
		   case 'n':
			index = (int) JUN;
			break;
		   case 'l':
			index = (int) JUL;
			break;
	           default:
		        error = 1;
	       }
	    }
	    break;
	case 'm':
	    switch (*(st+2)) {
		case 'r':
			index = (int) MAR;
			break;
		case 'y':
			index = (int) MAY;
			break;
	        default:
		        error = 1;
	    }
	    break;
	case 'n':
	    index = (int) NOV;
	    break;
	case 'o':
	    index = (int) OCT;
	    break;
	case 's':
	    index = (int) SEP;
	    break;
	default:
	   fprintf(stderr,"\n Unable to parse month in -M option.");
	   error = 1;  
			
   }  /* end switch */
   if (error) {
	fprintf(stderr,"\n -M%s ", st);
	fprintf(stderr,"\n Correct usage : -M(3_letter_month)FILEPATH \nex: -Mjan/Data/ASC/NordicSeas.jan.asc ");
        exit(1);
   }
   
   if ((mfile[index] = open_hydro_file("",st+3, "", 1)) < 0) {
	   fprintf(stderr,"\n Unable to open file for reading [%s]", st+3);
	   exit(1);
   
   }
   return(0);
	       

} /* end open_monthly_file() */
/*****************************************************************************/
double get_top_of_layer(int main_props_avail)
/* get depth associated with surface at top of layer. Return depth or -9999 if no
   depth can be determined.  */  
{
int j;
double z1;
  
      if (station.observ[prop1] == NULL)
         return(-9999.);

      z1 = hb_linterp(x1, station.observ[prop1], station.observ[(int)DE], hdr.nobs);
           
      if ((z1 < -8.) && (check_top_outcrop))  {  /* didn't find surf 1 */
           
        switch ((enum property) prop1) {
            case S0:
            case S1:   
            case S2:  
            case S3:
            case S4:
            case S_:
            case GN:
	    case PR:
	    case DE:
             if (station.observ[prop1][0] > x1) {
                z1 = station.observ[(int)DE][0];
             }
             break;
          default:
             ;
        }  /* end switch */
      
      }
            
      return(z1);
        
} /* end get_top_of_layer() */

           
/*****************************************************************************/
double get_bottom_of_layer(double z1, int main_props_avail, double *zperfect_addr)
/* get depth associated with surface at bottom of layer. Return depth or -9999 if no
   depth can be determined. zperfect_addr is the depth of the top of the layer upon entry 
   into this subroutine and the thickness of the layer for a complete water column upon exit 
 */  
{
int j;
double z2, zbotm;  
      
      
      if (station.observ[prop2] == NULL)
         return(-9999.);
           
      zbotm = x2;
      z2 = hb_linterp(x2, station.observ[prop2], station.observ[(int)DE], hdr.nobs);
      if ((z2 < 0) && check_bottom_outcrop) {
         zbotm =  hdr.pdr;
         if ((z1 >= 0) && (z2 < 0)) {  
           
           switch ((enum property) prop2) {
            case S0:
            case S1:   
            case S2:  
            case S3:
            case S4:
            case S_:
            case GN:
	    case PR:
	    case DE:
             if ( x2 > station.observ[prop2][hdr.nobs -1] ) {
                z2 = station.observ[(int)DE][hdr.nobs -1];
             }
             break;
             default:
                ;
           }  
         }
      }
      *zperfect_addr = zbotm - *zperfect_addr;
      return(z2);     
} /* end get_bottom_of_layer() */       
/********************************************/

int get_monthly_data(int file, struct PROFILE *data) 
/*  Reads each profile 
   in a HydroBase file and stores pr, th. sa. 
*/
{
   int error, i, j, row, col, sq;
   int main_props_avail;
   double  p1, p2, sbar, tbar, zperfect;
   
   
/* read each station in file ... */

    while ((error = get_station(file, &hdr, &station)) == 0) {

       /* is station within bounds ?   */
       
       if (xdateline && hdr.lon < 0)
          hdr.lon += 360.;
	         
       row = NINT((hdr.lat + 0.0001 - ymin) / yinc - 0.5);
       col = NINT((hdr.lon + 0.0001 - xmin) / xinc - 0.5);
       
       if ((row >= 0) && (row < nrows) && (col >= 0) && (col < ncols)) {

	  sq = row * ncols + col;
       /* ensure that pr, de, te, and sa are available to compute theta ... */
       
	  data[sq].nobs = 0;
          main_props_avail = 1;
          if (!(available(PR, &hdr) && available(DE, &hdr) && available(TE, &hdr) 
              && available(SA, &hdr))) {
            main_props_avail = 0;
          }
           zperfect = x1;                                       /* depth at top of layer requested  */
           p1 = get_top_of_layer(main_props_avail);
           p2 = get_bottom_of_layer(p1, main_props_avail, &zperfect);
	 
	   if (p1 >= 0. && p2 >= 0. && p2 > p1) {
	         if (!available(TH, &hdr) && main_props_avail) {
	     
                    free_and_alloc(&station.observ[(int)TH], hdr.nobs);
	            compute_theta(hdr.nobs, station.observ[(int)TH], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA]);
			     }  
		   data[sq].nobs =  hdr.nobs;
		   data[sq].d = (double *) calloc(hdr.nobs, sizeof(double));
		   data[sq].t = (double *) calloc(hdr.nobs, sizeof(double));
		   data[sq].s = (double *) calloc(hdr.nobs, sizeof(double));
		   
		   for ( i = 0; i < hdr.nobs; ++i) {
		      data[sq].d[i] = station.observ[(int)DE][i];
		      data[sq].t[i] = station.observ[(int)TH][i];
		      data[sq].s[i] = station.observ[(int)SA][i];
		   }
           }
      } /* end if in bounds */
      
      /* clean up...*/  
          
      for (i = 0; i < MAXPROP; ++i) {
         if (station.observ[i] != NULL) {
          free(station.observ[i]);
          station.observ[i] = NULL;
         }
      }    
   }  /* end while */

   if (error > 0)
      report_status(error, stderr);
      
   return;
} /* end get_monthly_data() */
/*****************************************************************************/
void get_hydro_data(int file)
   /*  Reads each station in a HydroBase file and computes anomaly of theta and salt  
   over specified layer.  All other necessary output information is computed and stored 
   in the global GRIDNODE array. */
{
   int error, i, j, row, col, sq;
   int main_props_avail;
   double  p1, p2, dir, pref=0.0;
   double rhobar, tbar, sbar, cpbar;
   double dx1, dx2, dy, dz, dTdz, dSdz, dzdiff, dz_std, zperfect;
   double *sp_heat = NULL;
   struct PROFILE *pptr;

/* read each station in file ... */

    while ((error = get_station(file, &hdr, &station)) == 0) {	  

       /* is station within bounds ?   */
       
       if (xdateline && hdr.lon < 0)
          hdr.lon += 360.;
	         
       row = NINT((hdr.lat + 0.0001 - ymin) / yinc - 0.5);
       col = NINT((hdr.lon + 0.0001 - xmin) / xinc - 0.5);
       sq = row * ncols + col;
       
       if ((row < 0) || (row >= nrows) || (col < 0) || (col >= ncols)) 
          continue;

       /* ensure that pr, de, te, and sa are available to compute rho ... */
       
       main_props_avail = 1;
       if (!(available(PR, &hdr) && available(DE, &hdr) && available(TE, &hdr) 
           && available(SA, &hdr))) {
         main_props_avail = 0;
       }
        zperfect = x1;   /* depth at top of layer requested */
        p1 = get_top_of_layer(main_props_avail);
        p2 = get_bottom_of_layer(p1, main_props_avail, &zperfect);
 	
	
	/* compute rho from specific volume */
	
        free_and_alloc(&station.observ[(int)SV], hdr.nobs);
	compute_sp_vol(hdr.nobs, station.observ[(int)SV], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA]);
	
	for (i=0; i < hdr.nobs; ++i) {
	   station.observ[(int)SV][i] = 1.0e8 / station.observ[(int)SV][i];
	}
	
	/* compute specific heat of seawater  */
	
        for (i = 0; i < hdr.nobs; ++i) {
             free_and_alloc(&sp_heat, hdr.nobs);
	     sp_heat[i] = hb_cpsw(station.observ[(int)SA][i], station.observ[(int)TE][i], station.observ[(int)PR][i]);
	}
	
	/* compute potential temperature */
       if (!available(TH, &hdr)) {
                free_and_alloc(&station.observ[(int)TH], hdr.nobs);
                
                compute_prop((int)TH, &station.observ[(int)TH], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA], hdr.nobs, pref, (double) hdr.lat);
       }
       
       tbar = get_weighted_avg(p1, p2, station.observ[(int)TH], station.observ[(int)SV], sp_heat, &rhobar, &cpbar, &sbar);
	  
       if (tbar > -3.0 ) {
       
            /* get climatological values */
	    
	 pptr = NULL;
	    
	 if (hdr.month < 1 || hdr.month > 12) {
	    fprintf(stderr, "\nWARNING: Unknown month encountered ...");
	    write_hydro_hdr(STDERR, &hdr);
	 }
	 else {
	    if (stdprofile[hdr.month] != NULL )
	       pptr = &stdprofile[hdr.month][sq];   /* use monthly profile if available */
	    else if (stdprofile[(int)CLIM] != NULL)
	       pptr = &stdprofile[(int)CLIM][sq];   /* use general climatological profiles */
	 }   
	    
	 if (pptr != NULL) {
	 
	 
            get_weighted_ts_anom(p1, p2, &dTdz, &dSdz, &dz_std, pptr);
	    if (dz_std > 0) {	 
	        dz = (p2 - p1)/1000.0;  /* change from m to km */
	        dz_std /= 1000.0;
	        zperfect /= 1000.0;
	    
	        ++grid[sq].count;
	    
	        if (grid[sq].dxdy <= 0.0) { 
	        /* compute area of grid square (trapezoid) in km^2 at this lat/lon */
	           hdr.lon = xmin + (col + 0.5) * xinc;
	           hdr.lat = ymin + (row + 0.5) * yinc;
	           dir = vector((double)hdr.lat-(0.5*yinc), (double)hdr.lon, (double)hdr.lat+(0.5*yinc), (double)hdr.lon, &dy);
	           dir = vector((double)hdr.lat-(0.5*yinc), (double)hdr.lon-(0.5*xinc), (double)hdr.lat-(0.5*yinc), (double)hdr.lon+(0.5*xinc), &dx1);
	           dir = vector((double)hdr.lat+(0.5*yinc), (double)hdr.lon-(0.5*xinc), (double)hdr.lat+(0.5*yinc), (double)hdr.lon+(0.5*xinc), &dx2);

                   grid[sq].dxdy = (dx1 + dx2) * 0.5 * dy;
	        }
                grid[sq].Tdz += tbar * dz;
                grid[sq].Sdz += sbar * dz;
                grid[sq].rho_dz += rhobar * dz;
                grid[sq].cp_dz += cpbar * dz;
	        grid[sq].dTdz += dTdz * dz;
	        grid[sq].dSdz += dSdz * dz;
	        grid[sq].dz_diff += dz - dz_std;
                grid[sq].dz += dz;
	        grid[sq].dz_perfect = zperfect;
	    }
	 }/* end if pptr */
	  	   
      } /* end if tbar */
         

  /* clean up...*/  
          
      for (i = 0; i < MAXPROP; ++i) {
         if (station.observ[i] != NULL) {
          free(station.observ[i]);
          station.observ[i] = NULL;
         }
      }    
   }  /* end while */

   if (error > 0)
      report_status(error, stderr);

   return;

} /* end get_hydro_data() */
/****************************************************************************/
void get_weighted_ts_anom(double p1, double p2, double *dTdz_addr, double *dSdz_addr, double *dz_addr, struct PROFILE *stdprofile)
 /*  computes weighted property anomalies for the depth interval
    specified by p1 and p2. */
{
   int i, n, start, end;
   double t1, t2, s1, s2;
   double tavg, savg, t1std, t2std, s1std, s2std;
   double dz, dz_sum, dTdz_sum, dSdz_sum;
   double *ttmp, *dtmp, *stmp;
 
   ttmp = (double *) calloc ((size_t)hdr.nobs, sizeof(double));
   stmp = (double *) calloc ((size_t)hdr.nobs, sizeof(double));
   dtmp = (double *) calloc ((size_t)hdr.nobs, sizeof(double));

   n = 0;
   for (i = 0; i < hdr.nobs; ++i) {
      if ((station.observ[(int)TH][i] > -8.0) && (station.observ[(int)SA][i] > 0.0) ) {
         ttmp[n] = station.observ[(int)TH][i];
         stmp[n] = station.observ[(int)SA][i];
         dtmp[n] = station.observ[(int)DE][i];
         ++n;
      }
   }
   if (n == 0) {
      free(dtmp);
      free(ttmp);
      free(stmp);
      *dTdz_addr = -999.;
      *dSdz_addr = -999.;
      *dz_addr = 0;
      return;
   }   
      
   t1 = hb_linterp(p1, dtmp, ttmp, n);
   s1 = hb_linterp(p1, dtmp, stmp, n);
   t2 = hb_linterp(p2, dtmp, ttmp, n);
   s2 = hb_linterp(p2, dtmp, stmp, n);
  
   if ( (t1 < -8.0) || (t2 < -8.0) || (t2 < -8.0) || (t2 < -8.0) ) {
      /* this shouldn't happen, but if it does, write a message */
      fprintf(stderr,"\n WARNING:  Program bug occurred in get_weighted_ts_anom()");
      free(dtmp);
      free(ttmp);
      free(stmp);
      return;
   }

   /* find index of first level greater than depth at top of layer */ 
     
   start = 0;
   while (p1 > dtmp[start])
     ++start;
   /*  find largest index of depth array less than depth at bottom of layer */
   end = n-1;
   while ((p2 <= dtmp[end]) && (end >= start))
     --end;

      
   /* now find average t and s for consecutive layers between top and bottom of big layer -- and associated t and s for same depth layers in standard profile.  Subtract to compute the anomaly and sum the anomalies weighted by their respective layer thickness... */
   
   dz_sum = 0.0;
   dTdz_sum = 0.0;
   dSdz_sum = 0.0;
    t1std = hb_linterp(p1, stdprofile->d, stdprofile->t, stdprofile->nobs);
    s1std = hb_linterp(p1, stdprofile->d, stdprofile->s, stdprofile->nobs); 
   
   for (i = start; i <= end; ++i) {

        t2std  = hb_linterp(dtmp[i], stdprofile->d, stdprofile->t, stdprofile->nobs);
        s2std = hb_linterp(dtmp[i], stdprofile->d, stdprofile->s, stdprofile->nobs); 
	
	/* if for some reason this layer doesn't exist in the std profile, skip the layer */

        if (t1std > -9.0 && t2std > -9.0   && s1std > -9.0 && s2std > -9.0)  {
	        
            dz = dtmp[i] - p1;
            dz_sum += dz;
	    
            tavg = (ttmp[i] + t1) * 0.5 ;
            savg = (stmp[i] + s1) * 0.5;
	    dTdz_sum +=  ( tavg - (t1std + t2std)  * 0.5 ) * dz;
	    dSdz_sum +=  ( savg - (s1std + s2std)  * 0.5 ) * dz;
	  }
	  
	 p1 = dtmp[i];
         t1 = ttmp[i];
         s1 = stmp[i];
	 t1std = t2std;
	 s1std = s2std;
	
  }
   
   /* add last depth interval ... */
   
   t2std  = hb_linterp(p2, stdprofile->d, stdprofile->t, stdprofile->nobs);
   s2std = hb_linterp(p2, stdprofile->d, stdprofile->s, stdprofile->nobs); 
   if (t1std > -9.0 && t2std > -9.0   && s1std > -9.0 && s2std > -9.0)  {
            dz = p2 - p1;
            dz_sum += dz;
            tavg = (t2 + t1) * 0.5 ;
            savg = (s2 + s1) * 0.5;
            dTdz_sum += ( tavg - (t1std + t2std)  * 0.5 ) * dz;
	    dSdz_sum +=  ( savg - (s1std + s2std)  * 0.5 ) * dz;
    }
   free(dtmp);
   free(ttmp);
   free(stmp);
   
   if (dz_sum <= 0) {
      *dTdz_addr = -999.;
      *dSdz_addr = -999.;
      *dz_addr = dz_sum;
      return;
   }
        
   *dTdz_addr = dTdz_sum / dz_sum;
   *dSdz_addr = dSdz_sum / dz_sum;
   *dz_addr = dz_sum;
   return;
} /* end get_weighted_ts_anom() */

/****************************************************************************/
double get_weighted_avg(double p1, double p2, double *xptr, double *rho_ptr,double *Cp_ptr, double *rho_avg_ptr, double *Cp_avg_ptr, double *sa_avg_ptr)

/*  computes weighted property averages over the depth interval
    specified by p1 and p2.  Returns the average property, salt, density, and spec.heat or -9999
    if no value can be computed. */
{
   int i, n, start, end;
   double xp1, xp2, rp1,rp2, Cp1, Cp2, sp1, sp2;
   double Cp_sum, rho_sum;
   double dz, dz_sum, theta_sum, salt_sum;
   double *xtmp, *dtmp, *rtmp, *Cptmp, *stmp;
  
   if (xptr == NULL) 
      return (-999.0);
   
   xtmp = (double *) calloc ((size_t)hdr.nobs, sizeof(double));
   stmp = (double *) calloc ((size_t)hdr.nobs, sizeof(double));
   dtmp = (double *) calloc ((size_t)hdr.nobs, sizeof(double));
   rtmp = (double *) calloc ((size_t)hdr.nobs, sizeof(double));
   Cptmp = (double *) calloc ((size_t)hdr.nobs, sizeof(double));
   n = 0;
   for (i = 0; i < hdr.nobs; ++i) {
      if ((xptr[i] > -8.0) && (station.observ[(int)SA][i] > 0.0) ) {
         xtmp[n] = xptr[i];
         stmp[n] = station.observ[(int)SA][i];
         dtmp[n] = station.observ[(int)DE][i];
	 rtmp[n] = rho_ptr[i];
	 Cptmp[n] = Cp_ptr[i];
         ++n;
      }
   }
   if (n == 0){
      free(dtmp);
      free(xtmp);
      free(stmp);
      free(rtmp);
      free(Cptmp);
      return (-999.);
   }   
      
   xp1 = hb_linterp(p1, dtmp, xtmp, n);
   xp2 = hb_linterp(p2, dtmp, xtmp, n);
   sp1 = hb_linterp(p1, dtmp, stmp, n);
   sp2 = hb_linterp(p2, dtmp, stmp, n);
   rp1 = hb_linterp(p1, dtmp, rtmp, n);
   rp2 = hb_linterp(p2, dtmp, rtmp, n);
   Cp1 = hb_linterp(p1, dtmp, Cptmp, n);
   Cp2 = hb_linterp(p2, dtmp, Cptmp, n);
  
   if ((xp1 < -8.0) || (xp2 < -8.0)) {
      free(dtmp);
      free(xtmp);
      free(stmp);
      free(rtmp);
      free(Cptmp);
      return (-999.);
   }

   /* find index of first level greater than depth at top of layer */   
   start = 0;
   while (p1 > dtmp[start])
     ++start;
   /*  find largest index of depth array less than depth at bottom of layer */
   end = n-1;
   while ((p2 <= dtmp[end]))
     --end;
   
   /* now do weighted average... */
   
   dz_sum = 0.0;
   rho_sum = 0.0;
   Cp_sum = 0.0;
   theta_sum = 0.0;
   salt_sum = 0.0;
   
   for (i = start; i <= end; ++i) {
         dz = dtmp[i] - p1;
         dz_sum += dz;
	 rho_sum += (rtmp[i] + rp1)  * 0.5 * dz;
	 Cp_sum += (Cptmp[i] + Cp1)  * 0.5 * dz;
         theta_sum += (xtmp[i] + xp1) * 0.5 * dz;
         salt_sum += (stmp[i] + sp1) * 0.5 * dz;
         p1 = dtmp[i];
         xp1 = xtmp[i];
         sp1 = stmp[i];
         rp1 = rtmp[i];
         Cp1 = Cptmp[i];
  }
   
   /* add last depth interval ... */
   
   dz = p2 - p1;
   dz_sum += dz;
   rho_sum += (rp2 + rp1)  * 0.5 * dz;
   Cp_sum += (Cp2 + Cp1)  * 0.5 * dz;
   theta_sum += (xp1 + xp2) *0.5 * dz;
   salt_sum += (sp1 +sp2) *0.5 * dz;
   
   free(dtmp);
   free(xtmp);
   free(stmp);
   free(rtmp);
   free(Cptmp);
   
   if (dz_sum <= 0)
     return(-9999.);
        
   *rho_avg_ptr = rho_sum / dz_sum;
   *Cp_avg_ptr = Cp_sum / dz_sum;
   *sa_avg_ptr = salt_sum / dz_sum;
   return (theta_sum / dz_sum);
   
}  /* end get_weighted_avg() */
/****************************************************************************/
void compute_prop(int i, double **xaddr, double *pptr, double *tptr, double *sptr, int n, double ref_val, double dlat)

/*
     int i;         index to enum property 
     double **xaddr address of array already allocated 
     double *pptr   pressure array 
     double *tptr   temperature array 
     double *sptr   salinity array 
     int n          number of obs 
     double ref_val ref pressure for s_, ht, or pe 
     double dlat    latitude 
 */
{

int j;
double *xptr, *xxptr;

   xptr = *xaddr;

/* !**! Special cases for individual properties... */
   switch ((enum property) i) {
       case TH:
               compute_theta(n, xptr, pptr, tptr, sptr);
               break;
       case S0:
               compute_sigma(0., n, xptr, pptr, tptr, sptr);
               break;
       case S1:
               compute_sigma(1000., n, xptr, pptr, tptr, sptr );
               break;
       case S2:
               compute_sigma(2000., n, xptr, pptr, tptr, sptr);
               break;
       case S3:
               compute_sigma(3000., n, xptr, pptr, tptr, sptr);
               break;
       case S4:
               compute_sigma(4000., n, xptr, pptr, tptr, sptr);
               break;
       case S_:
               compute_sigma(ref_val, n, xptr, pptr, tptr, sptr);
               break;
       case HT:
               compute_height(n, pptr, tptr, sptr, ref_val, xptr);
               break;
       case PE:
               compute_energy(n, pptr, tptr, sptr, ref_val, xptr);
               break;
       case SV:
               compute_sp_vol( n, xptr, pptr, tptr, sptr);
               break;

       case VA:
               compute_svan( n, xptr, pptr, tptr, sptr);
               break;

       case GN:
               compute_gamma_n(&ginfo, n, xptr, pptr, tptr, sptr, (double) hdr.lon, dlat);
	       break;
       case GE:
               xxptr = (double *) calloc((size_t)n, sizeof(double));
               compute_gamma_nerr(&ginfo, n, xxptr, xptr, pptr, tptr, sptr, (double) hdr.lon, dlat);
	       free((void *)xxptr);
	       break;
       default: 
               free(xptr);
               xptr = NULL;
               break;
  } /* end switch */
          
  return;

}   
/****************************************************************************/
int vector(double lat1, double lon1, double lat2, double lon2, double *dist_addr)
  /* Returns the direction from point1->point2 in degrees from north.
     The distance in km is returned at dist_addr */
{
   double dx, dy, phi;
   
   dx = (lon2 - lon1)* RADperDEG * cos(RADperDEG*.5 *(lat1+lat2)) * EarthRadius ;
   dy = (lat2 - lat1) * RADperDEG * EarthRadius;
   
   if (dy == 0) {  /* course is zonal */
     *dist_addr = dx * KMperNM;
     if (dx < 0) {
       *dist_addr = -dx * KMperNM;
       return(270);
     }
     return (90);
   }
   
   if (dx == 0.) {   /* meridional */
     *dist_addr = dy * KMperNM;
     if (dy < 0) {
       *dist_addr = -dy * KMperNM;
       return(180);
     }
     return(0);
   }

   phi = atan(dx/dy);
   *dist_addr = sqrt(dx*dx + dy*dy) * KMperNM;
   
   if (phi > 0) {   
      if (dx > 0)
         return ((int)((rint) (phi / RADperDEG)));   /* 0 -> 90 */
      
      return (180 + (int)((rint) (phi / RADperDEG))); /* 180 -> 270 */
   }
   
   if (dx > 0) 
     return (180 + (int)((rint) (phi / RADperDEG)));  /* 90 -> 180 */
     
   return (360 + (int)((rint) (phi / RADperDEG)));  /* 270 -> 360 */

}  /* end vector() */
      
/****************************************************************************/
void output_data(FILE *fptr)
 /*  Traverse the grid, weight and sum the values and output appropriate information */
{
  int i, row, col, sq;   
  double q_sum, f_sum, s_sum, t_sum, rho_sum, vol_sum, area_sum;
  double  dz_avg, deltaQ, deltaF;
  double lat, lon;
  struct GRIDNODE *curptr;

  q_sum = 0.0;
  f_sum = 0.0;
  s_sum = 0.0;
  t_sum = 0.0;
  rho_sum = 0.0;
  vol_sum = 0.0;
  area_sum = 0.0;

  sq = -1;   
  for (row = 0; row < nrows; ++row) {
     for (col = 0; col < ncols; ++col) {
        ++sq;  
        curptr = &grid[sq];
        if (curptr->count > 0 && curptr->dSdz > -9.) {
	
	   /* first divide each field by sum of weights */
	   
	   curptr->dTdz /= curptr->dz;
	   curptr->dSdz /= curptr->dz;
	   curptr->Tdz /= curptr->dz;
	   curptr->Sdz /= curptr->dz;
	   curptr->rho_dz /= curptr->dz;
	   curptr->cp_dz /= curptr->dz;
	   /* determine average thickness and compute heat and fresh water content changes */
           curptr->dz = curptr->dz / (double)curptr->count;  /* average thickness of layer in km*/
	   curptr->dz_diff = curptr->dz_diff /  (double)curptr->count;  /* difference of layer thickness for which anomaly was computed */
           deltaQ = curptr->dTdz * curptr->dz * 1000.0 * curptr->rho_dz * curptr->cp_dz;
	   deltaF = -(curptr->dSdz) * curptr->dz * 1000.0 / curptr->Sdz;
	   
           area_sum += curptr->dxdy; 
           vol_sum += curptr->dz * curptr->dxdy;
           q_sum += deltaQ * curptr->dxdy;
	   f_sum += deltaF * curptr->dxdy;
           t_sum += curptr->Tdz * curptr->dz * curptr->dxdy;
           s_sum += curptr->Sdz * curptr->dz * curptr->dxdy;
           rho_sum += curptr->rho_dz * curptr->dz * curptr->dxdy;
	   
	   lon = xmin + (col + 0.5) * xinc;
	   lat = ymin + (row + 0.5) * yinc;
     
           fprintf(fptr, "%8.3lf %8.3lf (%d) %.4e %.4e %8.4lf %8.4lf %9.4lf %6.4lf %9.3lf  [%.1lf] %.0lf%%  %.4e %.4e %.4lf \n",  lon, lat, curptr->count, curptr->dTdz, curptr->dSdz, curptr->Tdz, curptr->Sdz, curptr->rho_dz, curptr->cp_dz, curptr->dz, curptr->dz_diff,  (100. * curptr->dz/curptr->dz_perfect), curptr->dxdy, deltaQ, deltaF);
      
        } 	/* end if */
	
     }  /* end for col */
   } /* end for row*/
   
   dz_avg = vol_sum / area_sum;

   fprintf(fptr, "\nTotal area: %.4e km^2 \nAverage thickness of layer %9.3lf km \nAvg Theta per unit volume: %10.4lf (degC) \nAvg Salinity per unit volume: %8.4lf  \nAverage Density (rho): %10.4lf (kg/m^3) \nTemperature content  (deltaQ): %.8e (J/m^2) \nFresh water anomaly per unit area %.4lf (meters)\n", area_sum, dz_avg, t_sum/vol_sum, s_sum/vol_sum, rho_sum/vol_sum, q_sum/area_sum, f_sum/area_sum); 
  
  return;
   
} /* end output_data() */
