/*  hb_layerinteg.c

................................................................................
                          *******  HydroBase2 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             1996
			     updated to ANSI standard Dec 2001
................................................................................
____________________________________________________________________________
  USAGE:  

 hb_layerinteg filename(_roots) -1mne/value -2mne/value -P<prop> [-Wwindow/incr]  [-B<endlat/lon>] [-A<startlat/lon>] [-M] [-L] [-Y] [-I] [-D<dirname>] [-E<file_extent>]

  list of filenames MUST be first argument

 -1 : property and value at upper level;
 -2 : property and value at lower level;
 -P : property to integrate
 -A : starting lat/lon
 -B : bounds -- lonmin/lonmax/latmin/latmax
 -I : output id field

 -D : specifies directory for input files (default is current dir)  
      ex: -D/home/pete/data

 -E : specifies input_file extent (default is no extent) ex: -E.sum

____________________________________________________________________________
Computes an area-weighted integral of some property for a vertical interval 
defined by two surfaces. Weighting in the vertical is by rho*dz and in the 
horizontal by dist*dz.  
Outputs  a specified id (like year or latitude) avg_prop, avg_rho, thickness of layer .
If a start and end position are supplied then the stations are sorted by distance
from the start position.    Only stations lying between the two endpoints are 
included in the output.                             ____________________________________________________________________________
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

#define  PI     3.141592654
#define  RADperDEG 0.017453292    /* pi/180 */
#define  EarthRadius  3437.747    /* in nm */
#define  KMperNM  1.853248652


    /* flags to indicate which properties are to be output and which are
       needed for computation */
      
int prop_req;                  /* index to prop requested for output */
int prop_needed[MAXPROP];      /* set of props requested || defining a
                                    surface || needed for computation */
				    
  /*  define a structure to store layer information */

struct LAYER_INFO {
        int year;
	int month;
	int id;
	double lat;
	double lon;
	double dist;          /* distance (km) from starting point to sort list */
	double avg_prop;      /* weighted average over layer */    
	double avg_rho;       /* average density over layer */
	double dz;            /* thickness of layer */
        struct LAYER_INFO *next;
};				    
struct LAYER_INFO *list_ptr;

   /* global variables to store station ... */

struct HYDRO_DATA station;
struct HYDRO_HDR hdr;
double x1, x2;
double pref1, pref2;    /* ref pressures for top/bottom of layer */
int prop1, prop2;
int window, w_incr;     /* used to specify pr window for gradient properties*/

struct GAMMA_NC ginfo;   /* used for neutral density */

int warnflag, firsta;
int iflag, aflag, zflag;
int check_bottom_outcrop;
int xdateline;
double zmin, zmax;
double start_lat, start_lon;
double latmin, latmax, lonmin, lonmax;
char *identifier;

   /* prototypes for locally defined functions */

void    print_usage(char *);
void    get_hydro_data(int);
int     parse_prop_list(char *);
double get_weighted_avg(double, double, double *, double *, double *);
void compute_prop(int, double **, double *, double *, double *, int, double, double);
int 	vector(double, double, double, double,  double *);
struct LAYER_INFO *insert_rec(struct LAYER_INFO *, struct LAYER_INFO *);
void output_data(FILE *);

int main (int argc, char **argv)
{
   short   popt, opt1, opt2;
   int     nprops, i, n;
   int     curfile = 1, nfiles = 0; 
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
    popt = opt1 = opt2 = 0;
    error = 0;
    nprops = 0;
    window = 50;
    w_incr = 10;
    zmin = 0.0;
    zmax = 10000.0;
    infile = STDIN;
    iflag = zflag = 0;
    aflag = 0;
    warnflag = 0;           /* set to 1 after warning is printed */
    check_bottom_outcrop = 0;
    firsta = 1;
    list_ptr = (struct LAYER_INFO *) NULL;
    latmin = -90.0;
    latmax = 90.0;
    lonmin = -360.;
    lonmax = 360.;
    xdateline = 0;

/* initialize these ... */

   for (i = 0; i < MAXPROP; ++i) {
      station.observ[i] = (double *) NULL;
   }
   hdr.prop_id = (int *) NULL;



/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case '1':
                        opt1 = 1;
                        error = (sscanf(&argv[i][2],"%2s", id) == 1) ? 0 : 1;
                        if ((prop1 = get_prop_indx(id)) < 0) {
                           fprintf(stderr,"\n%.2s is not an appropriate property.", id);
                           ++error;
                        }
                        
                        if (prop1 == (int)S_ || prop1 == (int)HT || prop1 == (int)PE ) {
                         if (sscanf(&argv[i][4],"%lf/%lf", &pref1, &x1) != 2)
                            ++error;
                        }
                        else {
                         if (sscanf(&argv[i][5],"%lf", &x1) != 1)
                            ++error;
                        }
                        break;

               case '2':
                        opt2 = 1;
                        error = (sscanf(&argv[i][2],"%2s", id) == 1) ? 0 : 1;
                        if ((prop2 = get_prop_indx(id)) < 0) {
                           fprintf(stderr,"\n%.2s is not an appropriate property.", id);
                           ++error;
                        }
                        
                        if (prop2 == (int)S_ || prop2 == (int)HT || prop2 == (int)PE ) {
                         if (sscanf(&argv[i][4],"%lf/%lf", &pref2, &x2) != 2)
                            ++error;
                        }
                        else {
                         if (sscanf(&argv[i][5],"%lf", &x2) != 1)
                            ++error;
                        }
			n = strlen(argv[i]);
			if (argv[i][n-1] == 'b')   /* check for optional b */
			    check_bottom_outcrop = 1;
                        break;
               case 'A':
                        aflag = 1;
                        error = (sscanf(&argv[i][2],"%lf/%lf", &start_lat, &start_lon) == 2) ? 0 : 1;
                        break;
               case 'B':
                        error = (sscanf(&argv[i][2],"%lf/%lf/%lf/%lf", &lonmin, &lonmax, &latmin, &latmax) == 4) ? 0 : 1;

                        if (lonmin > lonmax) {
                           fprintf(stderr,"\nERROR specifying bounds:  lonmin must be <= lonmax");
			   exit(1);
			}			
                        if (latmin > latmax) {
                           fprintf(stderr,"\nERROR specifying bounds:  latmin must be <= latmax");
			   exit(1);
			}			
                        if (lonmax > 180)
                           xdateline = 1;
                        break;
               case 'D':                   /* get input dir */
                        dir = &argv[i][2];
                        break;
               case 'E':                    /* get file extent */
                        extent = &argv[i][2];
                        break;

               case 'I':
                        iflag = 1;
			identifier = &argv[i][2];
                        break;
	  
               case 'P':
                        popt = 1;
                        if ( argv[i][2] == '\0') {
                               print_prop_menu();
                               exit(0);
                        }
                        nprops = parse_prop_list(&argv[i][2]);
                        if (nprops <= 0)
                             error = 1;
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
                        
               case 'Z':
	                zflag = 1;
                        error = (sscanf(&argv[i][2],"%lf", &zmin) == 1) ? 0 : 1;
                        st = &argv[i][2];
                        while (*(st++) != '\0') {
                            if (*st == '/') {
                              ++st;
                              error = (sscanf(st,"%lf", &zmax) == 1) ? 0 : 1;
                              break;
                            }
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
             exit(1);
          }

       }  /* end if */

       else  {
           ++nfiles;
       }
   }  /* end for */

   if (!popt ) {
       fprintf(stderr,"\nYou must specify an output property with -P.\n");
       exit(1);
   }
   if (!opt1 || !opt2 ) {
       fprintf(stderr,"\nYou must define the layer with -1 and -2.\n");
       exit(1);
   }
   if (  !iflag) {
       fprintf(stderr,"\nWARNING: You have not specified any identifier}.  Only properties will be output.\n");
   }
   
   if (zflag) {
       fprintf(stderr,"\nUsing depth limits %.2lf - %.2lf db", zmin, zmax);
   }
   
  
   if ((prop1 == (int)GN) || (prop2 == (int)GN)) {
     gamma_nc_init((char *) GAMMA_NC_PATH, &ginfo);
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
   fprintf(stderr,"\n%s integrates a property over a layer specified by two surfaces", program);
   fprintf(stderr,"\nThe property is weighted in the vertical by rho and dz.");
   fprintf(stderr,"\nThe layer can be further limited to a pressure range using -Z<zmin/zmax>");     
   fprintf(stderr,"\nFor each station, outputs a line containing some combination of year, month, lat, lon, station_id  ");
   fprintf(stderr,"plus avg_prop, avg_rho, layer_thickness, distance_between_stations, avg_prop*avg_rho*layer_thickness*distance, cumulative_sum\n ");
   
    fprintf(stderr,"\n\nUsage:  %s filename_root(s)  -1<mne[pref]/value> -2<mne[pref]/value[b]> -P<list_of_properties> [-A<startlat/lon>] [-B<minlon/maxlon/minlat/maxlat>] [-D<dirname>] [-E<file_extent>] [-I] [-W<window>[/<w_incr>] [-Z<zmin/zmax>]", program);

   fprintf(stderr,"\n\n  List of filenames MUST be first argument");
   fprintf(stderr,"\n  or input is expected to come from stdin...");
   fprintf(stderr,"\n    -1  : surface 1 property/value (if property = s_");
   fprintf(stderr,"\n          you need to specify pref after s_");
   fprintf(stderr,"\n    -2  : surface 2 property/value (if property = s_");
   fprintf(stderr,"\n          you need to specify pref after s_");
   fprintf(stderr,"\n          Append b to to integrate to the bottom if this surface");
   fprintf(stderr,"\n          is deeper than the bottom.  This will only occur if the");
   fprintf(stderr,"\n          deepest observation is within 100 meters of the seafloor depth.");
   fprintf(stderr,"\n    -P  : property to integrate (only one allowed): ex: -Psa");
   fprintf(stderr,"\n          by itself, -P will list available properties");


   fprintf(stderr,"\n   [-A] : specify starting position (lat/lon) and output distance between stations");
   fprintf(stderr,"\n          The stations will be sorted by increasing distance from this point.");
   fprintf(stderr,"\n          If no starting point is specified, the stations will not be sorted");
   fprintf(stderr,"\n         but will be taken in the order they appear in the file.");
   fprintf(stderr,"\n   [-B] : specify geographic bounds.  Stations that do not ");
   fprintf(stderr,"\n          fall within bounds will not be included in the output. ");
   fprintf(stderr,"\n            ex: -B-75/-10/20/30 ");
 
   fprintf(stderr,"\n   [-D] : specifies dirname for input files (default is current");
   fprintf(stderr,"\n          directory)  ex: -D../data/ ");
   fprintf(stderr,"\n   [-E] : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n   [-I] : specify an identifier (char string) to be output (e.g. year, lat, lon");
   fprintf(stderr,"\n            Be sure to enclose in double quotes if string contains white space");
   fprintf(stderr,"\n            ex: -I1999  or -I24   or -IKN164  ");
   fprintf(stderr,"\n   [-W] : pressure window in db for computing gradient");
   fprintf(stderr,"\n          properties (bvfreq, pv...) The window is applied");
   fprintf(stderr,"\n          both above and below each observation level and");
   fprintf(stderr,"\n          constitutes the range over which observations");
   fprintf(stderr,"\n          are incorporated into the gradient computation.");
   fprintf(stderr,"\n          The optional w_incr value specifies an interval");
   fprintf(stderr,"\n          (db) into which the pressure range is subdivided");
   fprintf(stderr,"\n          to create an evenly spaced pr series over the");
   fprintf(stderr,"\n          window.  Default values: [50/10]");
   fprintf(stderr,"\n   [-Z] : specify pressure limits for layer.");
   fprintf(stderr,"\n          Exclude pressures that fall outside of these limits from weighted average.");
   fprintf(stderr,"\n   [-h] : help... prints this message.");
   fprintf(stderr,"\n\n");  
       
   return;
}


/*****************************************************************************/
/*****************************************************************************/
int parse_prop_list(char *st)
    /*  Parses a list of property mnemonics and sets global flags accordingly.
        Returns the number of properties.  An error will cause an exit. */
{
   int index;
   char prop[4];
   
      if (*st == '/')
         ++st;
      sscanf(st,"%2s", prop);
      index = get_prop_indx(prop);
      if (index < 0) {
         fprintf(stderr,"\n Unknown property requested: %s\n", prop);
         exit(1);
      }
      prop_req = index;
      prop_needed[index] = 1;
      ++st;
      ++st;

      if (*st != '\0') {
        fprintf(stderr,"\nOnly one property can be integrated in a run of this program.");
      }
      
      /* !**!  Special cases for properties */

      switch ((enum property)index) {
          case OX:
	  case O2:
	  case N2:
	  case N3:
	  case P4:
 	  case SI:
 	  case HE:
	  case SA: 
 	  case F1:
 	  case F2:
 	  case TU:
	  case TH:
	  case TE: 
	     /* these properties are sensible to integrate */
	             break;  
	  default:
	     fprintf(stderr,"\n Cannot compute volumetric quantity for %s ", get_prop_mne(index));
	     exit;
	  
      } /* end switch */

     /* end of special cases */

  return (1);
}  /* end parse_prop_list() */
/*****************************************************************************/
void get_hydro_data(int file)
   /*  Reads each station in a HydroBase file and integrates property 
   over specified layer.  All potentially necessary output information are stored 
   in a linked list of records. */
{
   int error, i, j, nreq;
   int main_props_avail;
   double  p1, p2, dir, pref=0.0;
   struct LAYER_INFO *rec_ptr, *r1ptr;


/* read each station in file ... */

    while ((error = get_station(file, &hdr, &station)) == 0) {

       if (xdateline && hdr.lon < 0)
          hdr.lon += 360;


     /* check for and skip out of bounds stations */  
          
    if ((hdr.lon <= lonmax) && (hdr.lon >= lonmin) && (hdr.lat <= latmax) && (hdr.lat >= lonmin))  {
	  
	  

       /* ensure that pr, de, te, and sa are available to compute rho ... */
       main_props_avail = 1;
       if (!(available(PR, &hdr) && available(DE, &hdr) && available(TE, &hdr) 
           && available(SA, &hdr))) {
         main_props_avail = 0;
       }

/* get depth associated with 1st surface...*/  
  
      if (!(available((enum property) prop1, &hdr) ) && main_props_avail ) {
         free_and_alloc(&station.observ[prop1], hdr.nobs);
	 if (prop1 == (int)OX) {
            if (available(O2, &hdr)) {
               for (j=0; j < hdr.nobs; ++j) {
                  station.observ[prop1][j] = ox_kg2l(station.observ[(int)O2][j], station.observ[(int)PR][j], station.observ[(int)TE][j],station.observ[(int)SA][j]);
               }
            }
	 }
	 else if (prop1 == (int)O2) {
            if (available(OX, &hdr)) {
              for (j=0; j < hdr.nobs; ++j) {
                 station.observ[prop1][j] = ox_l2kg(station.observ[(int)OX][j], station.observ[(int)PR][j], station.observ[(int)TE][j],station.observ[(int)SA][j]);
              }
            }
	 
	 }
         else 
	    compute_prop(prop1, &station.observ[prop1], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA], hdr.nobs, pref1, (double) hdr.lat); 
      }
       
      p1 = hb_linterp(x1, station.observ[prop1], station.observ[(int)DE],
           hdr.nobs);
           
      if (p1 < -8.) {  /* didn't find surf 1 */
           
        switch ((enum property) prop1) {
         case S0:
         case S1:   /* if surface density obs > density we are */
         case S2:   /* seeking, set p1 to first depth in array */
         case S3:
         case S4:
         case S_:
         case GN:
             if ((station.observ[prop1][0] > x1) 
             && (station.observ[(int)DE][0] < 150.)) {
                p1 = station.observ[(int)DE][0];
             }
             break;
          default:
             ;
        }  /* end switch */
      
      }
      
      if (zflag && (p1 >= 0.)) {   /* check depth limits */
          if (p1 < zmin)
	     p1 = zmin;
	  if (p1 > zmax)
	     p1 = -9999.;
      } 
           
/* get depth associated with 2nd surface...*/    
 
      if ((station.observ[prop2] == NULL) && (main_props_avail)) {
            free_and_alloc(&station.observ[prop2], hdr.nobs);
	    if (prop2 == (int)OX) {
               if (available(O2, &hdr)) {
                  for (j=0; j < hdr.nobs; ++j) {
                     station.observ[prop2][j] = ox_kg2l(station.observ[(int)O2][j], station.observ[(int)PR][j], station.observ[(int)TE][j],station.observ[(int)SA][j]);
                  }
               }
	    }
	    else if (prop2 == (int)O2) {
               if (available(OX, &hdr)) {
                 for (j=0; j < hdr.nobs; ++j) {
                    station.observ[prop2][j] = ox_l2kg(station.observ[(int)OX][j], station.observ[(int)PR][j], station.observ[(int)TE][j],station.observ[(int)SA][j]);
                 }
               }
	    
	    }
            else 
               compute_prop(prop2, &station.observ[prop2], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA], hdr.nobs, pref2, (double) hdr.lat);
      }
           
      p2 = hb_linterp(x2, station.observ[prop2], station.observ[(int)DE], hdr.nobs);

      if (zflag && (p2 >= 0.)) {   /* check depth limits */
          if (p2 < zmin)
	     p2 = -9999.;
	  if (p2 > zmax)
	     p2 = zmax;
      } 

      if (check_bottom_outcrop) {
         if ((p1 > -8) && (p2 < -8.)) {  
           
           switch ((enum property) prop2) {
            case S0:
            case S1:   
            case S2:  
            case S3:
            case S4:
            case S_:
	    case PR:
	    case DE:
             if ((x2 > station.observ[prop2][hdr.nobs -1] )
             && (hdr.pdr != 0) 
             && (station.observ[(int)DE][hdr.nobs -1] > (hdr.pdr - 100))) {
                p2 = station.observ[(int)DE][hdr.nobs -1];
             }
             break;
             default:
                ;
           }  
         }
      }

      if ((p1 > -1) && (p2 > -1) && (p2 >= p1)) {
      
        if (firsta) {
	  firsta = 0;
	  if (!aflag) {
	     start_lat = (double) hdr.lat;
	     start_lon = (double) hdr.lon;
	  }
	}
	
        rec_ptr = (struct LAYER_INFO *) calloc(1, sizeof(struct LAYER_INFO));
	if (rec_ptr == NULL) {
	   fprintf(stderr,"\nUnable to allocate memory in get_hydro_data()\n");
	   exit(1);
	}
	rec_ptr->year = hdr.year;
	rec_ptr->month = hdr.month;
	rec_ptr->lat = (double) hdr.lat;
	rec_ptr->lon = (double) hdr.lon;
	rec_ptr->id = hdr.station;
	
	/* compute rho from specific volume */
	
        free_and_alloc(&station.observ[(int)SV], hdr.nobs);
	compute_sp_vol(hdr.nobs, station.observ[(int)SV], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA]);
	for (i=0; i < hdr.nobs; ++i) {
	   station.observ[(int)SV][i] = 1.0e8 / station.observ[(int)SV][i];
	}
	
        if (!available((enum property) prop_req, &hdr)) {
                free_and_alloc(&station.observ[prop_req], hdr.nobs);
                
                compute_prop(prop_req, &station.observ[prop_req], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA], hdr.nobs, pref, (double) hdr.lat);
          }
          rec_ptr->avg_prop = get_weighted_avg(p1, p2, station.observ[prop_req], station.observ[(int)SV], &rec_ptr->avg_rho);
	  rec_ptr->dz = p2 - p1;
	  dir = vector(start_lat, start_lon, rec_ptr->lat, rec_ptr->lon, &rec_ptr->dist);
	   
	  if (aflag)                        /* sort the list by distance */
             list_ptr = insert_rec(list_ptr, rec_ptr);
	     
	  else {
	     if (list_ptr == NULL) {        /* empty list */
	       rec_ptr->next = NULL;
	       list_ptr = rec_ptr;
	     }
	     else {  /* traverse the linked list to the end and append record */
	        r1ptr = list_ptr;
		while (r1ptr->next != NULL)
		    r1ptr = r1ptr->next;
		    
	        rec_ptr->next = NULL;
		r1ptr->next = rec_ptr;
	     }
	  }
      }
         
    } /* end if !skipsta */

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
struct LAYER_INFO *insert_rec(struct LAYER_INFO *startptr, struct LAYER_INFO *rptr)
   /* inserts a record into a linked list sorted by increasing distance and returns 
      a pointer to the start of the list */
{

   struct LAYER_INFO *restoflist;

   /* case 1: empty list */
      
   if (startptr == NULL)  {       /* insert at end */
      rptr->next = NULL;
      return(rptr);
   }
   
   /* not an empty list */
   
   
   if (rptr->dist < startptr->dist) {  /* insert before start of list */
       rptr->next = startptr;
       return(rptr);
   }
   
   /* recursively search remainder of list */
   
   restoflist = insert_rec(startptr->next, rptr);
   startptr->next = restoflist;
   return (startptr);
   
} /* end insert_rec() */
/****************************************************************************/
double get_weighted_avg(double p1, double p2, double *xptr, double *rho_ptr, double *rho_avg_ptr)

/*  computes a weighted property average over the depth interval
    specified by p1 and p2.  Returns the average  prop and average density or -9999
    if no value can be computed. */
{
   int i, n, start, end;
   double xp1, xp2, rp1,rp2;
   double rho_dz, rho_dz_sum;
   double dz, dz_sum, layer_sum;
   double *xtmp, *dtmp, *rtmp;
  
   if (xptr == NULL) 
      return (-999.0);
   
   xtmp = (double *) calloc ((size_t)hdr.nobs, sizeof(double));
   dtmp = (double *) calloc ((size_t)hdr.nobs, sizeof(double));
   rtmp = (double *) calloc ((size_t)hdr.nobs, sizeof(double));
   n = 0;
   for (i = 0; i < hdr.nobs; ++i) {
      if (xptr[i] >= 0.0) {
         xtmp[n] = xptr[i];
         dtmp[n] = station.observ[(int)DE][i];
	 rtmp[n] = rho_ptr[i];
         ++n;
      }
   }
   if (n == 0){
      free(dtmp);
      free(xtmp);
      free(rtmp);
      return (-999.);
   }   
      
   xp1 = hb_linterp(p1, dtmp, xtmp, n);
   xp2 = hb_linterp(p2, dtmp, xtmp, n);
   rp1 = hb_linterp(p1, dtmp, rtmp, n);
   rp2 = hb_linterp(p2, dtmp, rtmp, n);
   
   if ((xp1 < 0.0) || (xp2 < 0.0)) {
      free(dtmp);
      free(xtmp);
      free(rtmp);
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
   rho_dz_sum = 0.0;
   layer_sum = 0.0;
   
   for (i = start; i <= end; ++i) {
         dz = dtmp[i] - p1;
	 rho_dz = (rtmp[i] + rp1)  * 0.5 * dz;
         dz_sum += dz;
	 rho_dz_sum += rho_dz;
         layer_sum += (xtmp[i] + xp1) * 0.5 * rho_dz;
         p1 = dtmp[i];
         xp1 = xtmp[i];
         rp1 = rtmp[i];
  }
   
   /* add last depth interval ... */
   
   dz = p2 - p1;
   dz_sum += dz;
   rho_dz = (rp2 + rp1)  * 0.5 * dz;
   rho_dz_sum += rho_dz;
   layer_sum += (xp1 + xp2) *0.5 * rho_dz;
   
   free(dtmp);
   free(xtmp);
   free(rtmp);
   
   *rho_avg_ptr = rho_dz_sum / dz_sum;
   return (layer_sum / rho_dz_sum);
   
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
       case TP:
               compute_dtdp(xptr,pptr,tptr,sptr,n,window,w_incr);
               break;
       case TZ:
	       compute_dtdz(xptr,pptr,tptr,sptr,n,window,w_incr,dlat);
               break;
       case HC:
               compute_heat_capacity(n,xptr,pptr,tptr,sptr);
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

       case PV: 
               buoy_freq(xptr, pptr, tptr, sptr, n, window, w_incr);
               for (j = 0; j < n; ++j) {
                 if (xptr[j] > -99.) 
                   xptr[j] *= xptr[j];
               }
               po_vort(xptr, xptr, n, dlat);

               break;

       case RR:
	 fprintf(stderr,"\nhb_layerinteg: ERROR: Approx. Rossby radius not relevant here!");
               free(xptr);
               xptr = NULL;
               break;
		 
       case BF:
               buoy_freq(xptr, pptr, tptr, sptr, n, window, w_incr);
               for (j = 0; j < n; ++j) {
                    if (xptr[j] > -9998.)   /* convert the units */
                        xptr[j] *= 1.0e5;
               }
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
 /*  Traverse the linked list, compute distance between stations 
      convert everything to meters, and output appropriate information */
{
      
  struct LAYER_INFO *curptr, *nextptr;
  double dist_prev, dist_next, ddist;
  double area, x_sum, rho_sum, area_sum, ddist_sum;
  int azimuth;
     
  curptr = list_ptr;
  nextptr = list_ptr->next;
  dist_next = 0.0;
  x_sum = 0.0;
  rho_sum = 0.0;
  area_sum = 0.0;
  ddist_sum = 0.0;
   
   while (curptr != NULL) {
   
     dist_prev = dist_next;
     if (nextptr == NULL)
         dist_next = 0.0;
     else
         azimuth = vector(curptr->lat, curptr->lon, nextptr->lat, nextptr->lon, &dist_next);
	 
     
     
     ddist = (dist_prev + dist_next) * 0.5;    

     if (ddist == 0.0)  /*there is a single station to evaluate; set distance to unit */
        ddist = 1.0;
	
	
     /* No need to correct units that are in km to m:
         rho = kg/m^3   
	  dz = m        
       ddist = km  since this is just a weight  */
     	
     area = curptr->dz * ddist;
     area_sum += area; 
     ddist_sum += ddist;
     x_sum += curptr->avg_prop * area;
     rho_sum += curptr->avg_rho * area;
     
     
     curptr = nextptr;
     nextptr = nextptr->next;
	
   } /* end while */
   
	
   if (iflag)
        fprintf(fptr, "%s", identifier);

   fprintf(fptr, "  %10.4lf %10.4lf %10.4lf %10.4lf\n",  x_sum/area_sum, rho_sum/area_sum, area_sum/ddist_sum,  x_sum/area_sum * rho_sum/area_sum);
   
   
} /* end output_data() */
