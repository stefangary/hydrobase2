/*   hb_slice.c
................................................................................
                              *  HydroBase2 *
................................................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             May 2001 
................................................................................
/*  For a list of positions (lon/lat pairs) and set of HydroBase gridded 
    cdf files, produces a HydroBase station file containing "stations" 
    that lie along the line connecting the positions.  In essence the
    positions define a pathway along which a vertical slice is cut
    through the 3D gridded property fields.

*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "hydrobase.h"
#include "hydro_cdf.h"
#include "hb_memory.h"
#include "hb_gamma.h"
#include "hb_paths.h"


/* input_file pathnames */

#define    EXTENT   ""
#define    DIR      ""

#define PI     3.141592654
#define  RADperDEG 0.017453292   /* pi/180 */
#define  EarthRadius  3437.747    /* in nm */

struct POS_REC {
   double lat, lon, dist;
   int type;
};


/* globally defined variables */

int prop_req[MAXPROP];         /* set of props requested for output */
int prop_needed[MAXPROP];      /* set of props requested || defining a
                                    surface || needed for computation */
double s_pref;           /* ref lev for computing a non-standard sigma level */
double pe_pref;           /* ref lev for computing potential energy */
double ht_pref;           /* ref lev for computing dynamic height */
int window, w_incr;     /* used to specify pr window for gradient properties*/

int xdateline;   

struct GAMMA_NC ginfo;   /* used for neutral density */
  
/* prototypes for locally defined functions */	

void print_usage(char *);
int parse_p_option(char *);
int vector(double, double, double, double, double *, double *);
int compare_points(const void *, const void *);
double getlat(double, double, double, double, double *);
double getlon(double, double, double, double, double *);
int check_bounds(struct CDF_HDR *, float, float);
int load_properties(int, struct CDF_HDR *, struct HYDRO_DATA *, float, float);
int get_profile(struct POS_REC *, int, char **, char *, char *, struct HYDRO_HDR *, struct HYDRO_DATA *);
void compute_props(struct HYDRO_HDR *, struct HYDRO_DATA *);


main(int argc, char **argv) {
  int i, j,  n, npoints;
  int nfiles, outfile, nalloc;
  int error, nz;
  int nprops;
  int pixel_grid, quadrant;
  int iflag, pflag, lflag;
  int hdg, zonal, meridional;
  int startpoint;
  char *dir, *extent;
  char *st;
  float xinc, yinc;
  float xoffset, yoffset;
  double phi, dist, startdist;
  double lat0, lon0, lonmin, lonmax;
  double startlat, endlat, startlon, endlon;
  FILE *posfile;
  struct POS_REC *point;
  struct HYDRO_DATA data;
  struct HYDRO_HDR hdr;
   
  /* Set default values. */
  error = 0;
  dir = DIR;
  extent = EXTENT;
  nfiles = 0;
  outfile = STDOUT;
  pixel_grid = 1;
  iflag = pflag = lflag = 0;
  posfile = stdin;
  window = 100;
  w_incr = 10;

  /* Initialize variables ... */
  hdr.prop_id = NULL;
  for (i = 0; i < MAXPROP; ++i) 
    data.observ[i] = (double *) NULL;
      
  /* Set up header for creating HydroBase stations... */
  strncpy(hdr.country, "XX", 3);
  strncpy(hdr.ship, "XX", 3);
  hdr.cruise = 999;
  hdr.station = 999;
  hdr.year = 9999;
  hdr.month = 99;
  hdr.day = 99;
  hdr.qual[0] = '0';
  hdr.qual[1] = '0';
  hdr.qual[2] = '0';
  hdr.qual[3] = '0';
  hdr.instrument = 'u';
  hdr.origin = '0';
   
  /*----------------------------------------*/  
  
  /* Are there command line arguments? */
  if (argc < 1) {
      print_usage(argv[0]);
      exit(1);
   }
   
  /*----------------------------------------*/ 
 
  /* Parse command line arguments */

  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
        case 'D':                   /* get input dir */
          dir = &argv[i][2];
          break;
	  
        case 'E':                    /* get file extent */
          extent = &argv[i][2];
          break;
	  
        case 'G':        /* set to gridnode registration */
	  pixel_grid = 0;
	  break;
	  
        case 'I':
          iflag = 1;
          error = (sscanf(&argv[i][2],"%f", &xinc) == 1) ? 0 : 1;
          yinc = xinc;
          st = &argv[i][2];
          while (*(st++) != '\0') {
            if (*st == '/') {
               ++st;
               error = (sscanf(st,"%f", &yinc) == 1) ? 0 : 1;
               break;
            }
          }
                        
          break;
	  
         case 'L':
	  lflag = 1;
	  posfile = fopen(&argv[i][2],"r");
	  if (posfile == NULL) {
	    fprintf(stderr,"\nUnable to open %s for input.\n", &argv[i][2]);
	    exit(1);
	  }
	  break;
	  
       case 'O':
          outfile = create_hydro_file(&argv[i][2], OVERWRITE);
	  if (outfile < 0) {
	    fprintf(stderr,"\nUnable to open %s for output.", &argv[i][2]);
	    exit(1);
	  }
          break;
	  
        case 'P':
          pflag = 1;
          nprops = parse_p_option(&argv[i][2]);
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
      ++nfiles;
    }  /* end else */
    
  }  /* end for */

  /*--------------------------------------------*/    
  /*  Check syntax of arguments */ 
  error = 0;
   
  if (!nfiles ) {
    fprintf(stderr,"\nYou must specify input cdf_files as first argument(s).");
    ++error;    
  }

  if (!iflag ) {
    fprintf(stderr,"\nYou must specify grid spacing with -I<xincr>[/yincr] ");
    ++error;
  }
    
  if (!lflag ) {
    fprintf(stderr,"\nExpecting list of points from stdin: \n");
  }
  
  if (!pflag ) {
    fprintf(stderr,"\nYou must specify properties for output: -Ppr/te/sa/ox ");
    ++error;
  }
    
  if (error) {
    fprintf(stderr,"\nUse -h for complete usage info. \n");
    exit(1);
  }

  xoffset = yoffset = 0.0;  /* if gridlines on whole degrees */
   
  if (pixel_grid) {         /* if not...*/
    xoffset = 0.5 * xinc;
    yoffset = 0.5 * yinc;
  }

  /*--------------------------------------------*/
    
  /* Read in first  lat/lon point */

  nalloc = 10000;

  point = (struct POS_REC *) get_memory((void *)NULL, (size_t)nalloc, sizeof(struct POS_REC));
   
  n = fscanf(posfile, "%lf%lf",  &startlon, &startlat);
  if (n != 2) {
    fprintf(stderr,"\nError reading startlon, startlat\n");
    exit(1);
  }
   
  point[0].lat = startlat;
  point[0].lon = startlon;
  point[0].dist = 0.0;
  point[0].type = 2;
  
  lonmin = 999.0;
  lonmax = -999.0;
  startdist = 0.0;
  
  if (startlon > lonmax) lonmax = startlon;
  if (startlon < lonmin) lonmin = startlon;

  fprintf(stderr,"hb_slice got first point in slice path.\n");

  /*--------------------------------------------*/    
  /* Get consecutive points, construct an array
   * of gridline crossings between each pair of
   * points. */
          
  npoints = 1;
  startpoint = 0;  
  while  ((n = fscanf(posfile, "%lf%lf",  &endlon, &endlat)) == 2) { 
    if (endlon > lonmax) lonmax = endlon;
    if (endlon < lonmin) lonmin = endlon;
    hdg = vector(startlat, startlon, endlat, endlon, &phi, &dist);
    dist += startdist;
    
    if (hdg >= 360) 
      hdg = 360 - hdg;
    
    meridional = zonal = 0;
    
    if (hdg == 180 || hdg == 0 )
      meridional = 1;
    if (hdg == 90 || hdg == 270)
        zonal = 1;
    
    if (hdg < 90)
      quadrant = 0;
    else if (hdg < 180)
      quadrant = 1;
    else if (hdg < 270)
      quadrant = 2;
    else 
      quadrant = 3;
    
    /*--------------------------------------------*/    
    if (!zonal) {
      
      /* Get latitudinal gridline xings */
      
      lat0 = (double) NINT(startlat);
      if (quadrant == 0 || quadrant == 3) {
	
	/* Heading north */
	
	lat0 -=  1;
	lat0 += yoffset;      /* find nearest y-gridline */
	
	/* find first gridline past startlat */  
	while ((lat0 += yinc) <= startlat)  
	  ;  
	
	/* find all lat gridline xings until endlat */ 
	while (lat0 < endlat) {
	  point[npoints].lat = lat0;
	  point[npoints].lon = getlon(lat0, startlon, startlat, phi, &point[npoints].dist);
	  point[npoints].dist += startdist;
	  point[npoints].type =  1;
	  if (++npoints > nalloc) {
	    nalloc += 100;
	    point = (struct POS_REC *) get_memory((void *)point, (size_t) nalloc, sizeof(struct POS_REC));
	  }
	  
	  lat0 += yinc;
	}
	
      }  /* end if */
      
      
      else {                 /* heading south */
	
	/* find nearest whole degree */
	lat0 +=  1;
	lat0 -= yoffset;      /* find nearest y-gridline */
	
	/* find first gridline past startlat */  
	while ((lat0 -= yinc) >= startlat)  
	  ;  
	
	/* find all lat gridline xings until endlat */ 
	    while (lat0 > endlat) {
	      point[npoints].lat = lat0;
	      point[npoints].lon = getlon(lat0, startlon, startlat, phi, &point[npoints].dist);
	      point[npoints].dist += startdist;
	      point[npoints].type =  1;
	      if (++npoints > nalloc) {
		nalloc += 100;
		point = (struct POS_REC *) get_memory((void *)point, (size_t) nalloc, sizeof(struct POS_REC));
	      }
	      
	      lat0 -= yinc;
	    }
	    
      } /* end else */ 
      
    } /* end if !zonal */
    
    
    /*--------------------------------------------*/    
    if (!meridional) { 
      
      /* get longitudinal gridline xings */
      
      lon0 = (double) NINT(startlon);
      if (quadrant == 0 || quadrant == 1) {
	
	/* heading east */
	
	lon0 -=  1;
	lon0 += xoffset;      /* find nearest x-gridline */
	
	/* find first gridline past startlon */  
	while ((lon0 += xinc) <= startlon)  
	  ;  
	
	/* find all lon gridline xings until endlon */ 
	while (lon0 < endlon) {
	  point[npoints].lon = lon0;
	  point[npoints].lat = getlat(lon0, startlon, startlat, phi, &point[npoints].dist);
	  point[npoints].dist += startdist;
	  point[npoints].type =  0;
	  if (++npoints > nalloc) {
	    nalloc += 100;
	    point = (struct POS_REC *) get_memory((void *)point, (size_t) nalloc, sizeof(struct POS_REC));
	  }
	  
	  lon0 += xinc;
	}
	
      }  /* end if */
      
      
      else {                 /* heading west */
	
	lon0 +=  1;
	lon0 -= xoffset;      /* find nearest x-gridline */
	
	/* find first gridline past startlon */  
	while ((lon0 -= xinc) >= startlon)  
	  ;  
	
	/* find all lon gridline xings until endlon */ 
	while (lon0 > endlon) {
	  point[npoints].lon = lon0;
	  point[npoints].lat = getlat(lon0, startlon, startlat, phi, &point[npoints].dist);
	  point[npoints].dist += startdist;
	  point[npoints].type =  0;
	  if (++npoints > nalloc) {
	    nalloc += 100;
	    point = (struct POS_REC *) get_memory((void *)point, (size_t) nalloc, sizeof(struct POS_REC));
	  }
	  
	  lon0 -= xinc;
	}
	
      } /* end else */ 
      
      
    } /* end if !meridional */
    
    /*--------------------------------------------*/    
    /* Add the endpoint to the array */
    
    point[npoints].lat = endlat;
    point[npoints].lon = endlon;
    point[npoints].dist = dist;
    point[npoints].type = 2;
    ++npoints;
    
    fprintf(stderr,"hb_slice got next point in slice path.\n");

    /* sort the array of points in order by increasing dist from startpoint */
    
    qsort((void *) &point[startpoint], npoints - startpoint, sizeof(struct POS_REC), compare_points);
    
    
    startlat = endlat;
    startlon = endlon;
    startpoint = npoints;
    startdist = dist;
    
  } /* end while sscanf(posfile) */

  if (lflag)
    fclose(posfile);
  
  if (n != EOF) {
    fprintf(stderr,"\n ERROR reading lon/lat pair.... continuing.\n");
  }
  
  if (npoints < 2) {
    fprintf(stderr,"\n You must provide at least 2 lon/lat end points.\n");
    fprintf(stderr,"\nUse -h for complete usage info. \n");
    exit(1);
  }
  
  point = get_memory((void *)point, (size_t)npoints, sizeof(struct POS_REC));
  
  fprintf(stderr,"hb_slice done creating pathway.\n");

  /*--------------------------------------------*/    
  /*        end of defining pathway */
  /*--------------------------------------------*/    
  
  xdateline = 0;
  xdateline = (lonmin < 180. && lonmax > 180.);
  if (xdateline) {
    fprintf(stderr,"Dateline is crossed:  converting longitudes 0 -> 360");
  }
  
  /*  For each point along pathway, obtain a profile from cdf files */
  fprintf(stderr,"hb_slice starting to get data from cdf file...\n");
  
  for (i = 0; i < npoints; ++i) {

    iflag = get_profile(&point[i], nfiles, argv, dir, extent, &hdr, &data);

    if ( iflag ) {

      /* Got valid data so write it to output file. */
      if (hdr.nobs > 0) {
	compute_props(&hdr, &data);
	write_hydro_station(outfile, &hdr, &data);
      }

      /* Deallocate memory */
      free((void *)hdr.prop_id);
      for (n = 0; n < MAXPROP; ++n) {
	if (data.observ[n] != NULL) {
	  free((void *) data.observ[n]);
	  data.observ[n] = NULL;
	}
      }
    }
  } /* end for i */

  fprintf(stderr,"\nhb_slice done looping over each profile point.");
  
  fprintf(stderr,"\n\nEnd of %s\n", argv[0]);
  exit(0);
}  /* end main */

/************************************************************************/
void print_usage(char *program)
{
   fprintf(stderr,"\nUsage:  %s cdf_file(s) -I<delta_x[/delta_y]>  -P<list_of_properties> [-O<outfile>] [-G] [-L<position_file>] [-D<dirname>] [-E<file_extent>] -W<window>[/<w_incr>]] ", program);

   fprintf(stderr,"\n\n  List of filenames must be first argument ");
   fprintf(stderr,"\n -I  : xgrid[/ygrid] increment in degrees.  ");
   fprintf(stderr,"\n -P  : list of properties for output;");
   fprintf(stderr,"\n       ex:  -Ppr/th/sa/ox/ht");
   fprintf(stderr,"\n            -P (by itself) produces a list of available properties");
   fprintf(stderr,"\n\n   OPTIONS:  ");
   fprintf(stderr,"\n[-D] : directory for input cdf files  ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n[-E] : file extent for input cdf files");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n[-G] : Gridnodes fall on whole degrees in cdf files");
   fprintf(stderr,"\n       default is gridlines offset by 0.5 * gridinc (pixel registration)");
   fprintf(stderr,"\n[-L] : file containing list of lon/lat positions.");
   fprintf(stderr,"\n       If not specified, this list will be read from stdin ");
   fprintf(stderr,"\n[-O] : name of output file");
   fprintf(stderr,"\n[-W] : Specifies pressure window (db) for computing ");
   fprintf(stderr,"\n       gradient properties (bvfreq, pv...) This ");
   fprintf(stderr,"\n       constitutes the range over which observations");
   fprintf(stderr,"\n       are incorporated into the gradient computation.");
   fprintf(stderr,"\n       The window is centered around each pressure ");
   fprintf(stderr,"\n       level being evaluated.");
   fprintf(stderr,"\n       w_incr specifies how finely to subdivide the");
   fprintf(stderr,"\n       window into increments(db) to create");
   fprintf(stderr,"\n       an evenly spaced pressure series over the window.");
   fprintf(stderr,"\n       defaults: -W%d/%d", window, w_incr);
   fprintf(stderr,"\n[-h] : help -- prints this message.");
   fprintf(stderr,"\n\n");  
   return;
} /* end print_usage() */

/*****************************************************************************/
int parse_p_option(char *st)
    /*  Parses a list of property mnemonics and sets global flags accordingly.
        Returns the number of properties.  An error will cause an exit. */
{
   int index, nprops;
   char prop[4];
   double ref_val;

   nprops = 0;

   do {
      if (*st == '/')
         ++st;
      sscanf(st,"%2s", prop);
      index = get_prop_indx(prop);
      if (index < 0) {
         fprintf(stderr,"\n Unknown property requested: %s\n", prop);
         exit(1);
      }
      prop_req[index] = 1;
      prop_needed[index] = 1;
      ++nprops;
      ++st;
      ++st;

     /* !**!  Special cases for properties */

      if (((enum property)index == S_ ) || ((enum property) index == HT) || ((enum property) index == PE)) {
         if (sscanf(st, "%lf", &ref_val) != 1) {
             fprintf(stderr, "\n Specify a ref pressure for  %.2s", prop);
             fprintf(stderr,"\n   ex: -P%.2s1500/th/sa\n", prop);
             exit(1);
         }
         while (!(*st=='/' || *st=='\0' || *st==' '))
            ++st;

        switch ((enum property) index) {
           case S_:
              s_pref = ref_val;
              break;
           case PE:
              pe_pref = ref_val;
              break;
           case HT:
              ht_pref = ref_val;
              break;
           default:
              ;
        }  /* end switch */
      }

   /* end of special cases */

  } while (*st == '/');
  
  if (prop_req[(int)GE] || prop_req[(int)GN]) 
     gamma_nc_init((char *) GAMMA_NC_PATH, &ginfo);
  

  return (nprops);
}  /* end parse_p_option() */
/***************************************************************************/
int compare_points(const void *pt1, const void *pt2)
   /* Routine for qsort to sort structure by its .dist field.
      Returns -1, 0, or 1 for pt1 < pt2, pt1 == pt2, or pt1 > pt2 */
{
   struct POS_REC *p1, *p2;
   double dist1, dist2;
   
   p1 = (struct POS_REC *)pt1;
   p2 = (struct POS_REC *)pt2;
   dist1 = p1->dist;
   dist2 = p2->dist;
   
   if  (ABS(dist1 - dist2) <= 1) 
      return(0);
   
   if (dist1 > dist2)
      return(1);
    
    return(-1);
   
}
/**************************************************************************/
int vector(double lat1, double lon1, double lat2, double lon2, double *phi_addr, double *dist_addr)
  /* Returns the direction from point1->point2 in degrees from north.
     The arctan of this direction in radians  (domain is -pi/2 to +pi/2)
     is returned at phi_addr. The distance in nautical miles is returned
     at dist_addr */
{
   int hdg;
   double dx, dy;
   
   dx = (lon2 - lon1)* RADperDEG * cos(RADperDEG*.5 *(lat1+lat2)) * EarthRadius ;
   dy = (lat2 - lat1) * RADperDEG * EarthRadius;
   
   if (dy == 0) {  /* course is zonal */
     *dist_addr = dx;
     if (dx < 0) {
       *phi_addr = -PI/2;
       *dist_addr = -dx;
       return(270);
     }
     *phi_addr = PI /2;
     return (90);
   }
   
   if (dx == 0.) {   /* meridional */
     *dist_addr = dy;
     *phi_addr = 0.0;
     if (dy < 0) {
       *dist_addr = -dy;
       return(180);
     }
     return(0);
   }

   *phi_addr = atan(dx/dy);
   *dist_addr = sqrt(dx*dx + dy*dy);
   
   if (*phi_addr > 0) {   
      if (dx > 0){
         hdg = (int) NINT(*phi_addr / RADperDEG);   /* 0 -> 90 */
	 return(hdg);
      }
      
      return ((180 + NINT(*phi_addr / RADperDEG))); /* 180 -> 270 */
   }
   
   if (dx > 0) 
     return (180 + NINT (*phi_addr / RADperDEG));  /* 90 -> 180 */
     
   return (360 + NINT(*phi_addr / RADperDEG));  /* 270 -> 360 */

}  /* end vector() */
      
/****************************************************************************/
double getlat(double lon2, double lon1, double lat1, double phi, double *dist_addr)
  /* Returns the latitude that intersects at lon2 the ray whose vertex
     is lon1,lat1 with angle phi in radians and domain -pi/2 to pi/2.
     The distance in naut miles between those points is also returned.  */
{
   double dx, dy, lat2;
   int i;
   
   lat2 = lat1;
   
   /* iterate a few times to converge on lat2 */
   
   for (i = 0; i < 10; ++i) {
     dx = (lon2 - lon1)* RADperDEG * cos(RADperDEG * 0.5*(lat2+lat1)) * EarthRadius;
     dy = dx / tan(phi);
     lat2 = dy /(RADperDEG * EarthRadius) + lat1;
   }
   
   *dist_addr = sqrt(dx*dx + dy*dy);
   return( dy /(RADperDEG * EarthRadius) + lat1);
   
}

/****************************************************************************/
double getlon(double lat2, double lon1, double lat1, double phi, double *dist_addr)
  /* Returns the longitude that intersects at lat2 the ray whose vertex
     is lon1,lat1 with angle phi in radians and domain -pi/2 to pi/2.
      The distance in naut miles between those points is also returned.    */
{
   double dx, dy;
   
   dy = (lat2 - lat1) * RADperDEG * EarthRadius;
   dx =  dy * tan(phi);
   
   *dist_addr = sqrt(dx*dx + dy*dy);
   
   return(dx / (EarthRadius * cos(RADperDEG*.5 *(lat1+lat2)) * RADperDEG) + lon1);
  
}
/****************************************************************************/
int get_profile(struct POS_REC *pos, int nfiles, char **arglist, char *dir, char *extent, struct HYDRO_HDR *hptr, struct HYDRO_DATA *dptr)

   /* Returns 1 if successful, 0 if not */
{
   int i, j, n, cdfid, error, row, col;
   int hdg, ip, nz, np2;
   int found, dcl;
   struct CDF_HDR cdf;
   float lat, lon;
   double phi, dist, dist2;
   double flag, weight1, weight2;
   struct HYDRO_DATA *dptr2;
     
   int print_msg = 0;
   int curfile = 1;
   
   
   flag = (double) -8.9;

   if (xdateline) {      /* ensure longitude range 0 -> 360 */
      if (pos->lon < 0)
         pos->lon += 360.0;
   }
   else {                  /* ensure longitude range -180 -> 180 */
      if (pos->lon > 180.)
         pos->lon -= 360.0;
   }

   hptr->lat = (float) pos->lat;
   hptr->lon = (float) pos->lon;
   hptr->ms10 = ms10(hptr->lat, hptr->lon, &hptr->ms1);

   do {
   
      cdfid = cdf_open(dir, arglist[curfile], extent, print_msg);
      if (error = read_cdf_hdr(cdfid, &cdf)) 
         exit (1);
	 
      found = check_bounds(&cdf, hptr->lat, hptr->lon);   
	 
      if (!found) 
         cdf_close(cdfid);
      
   }  while (!found && curfile++ < nfiles);

   if (!found) {
     /* SFG added brackets and fprintf. */
     fprintf(stderr,"Profile not found, skipping!");
     return (0);
   }
 /*  utilize this point ...	  */
 
   hptr->nprops = load_properties(cdfid, &cdf, dptr, hptr->lat, hptr->lon);
   
   hptr->prop_id = (int *) calloc((size_t)hptr->nprops, sizeof(int));
   
   /* Now that the prop_id vector is allocated, it must
    * be populated.  First initialize a counter for the
    * number of active properties, n.*/
   n = 0;
   for (i = 0; i < MAXPROP; ++i) {
     if (dptr->observ[i] != NULL){
       hptr->prop_id[n++] = i;
	/* For the case where there is space allocated for this property,
	 * simultaneously augment the number of active properties
	 * and set the property id in the header. For locations where
	 * there is no data then only depth and pressure can be set.
	 * SFG added brackets after the if() above and below the fprints
	 * below and this changes the execution behavior! Otherwise, this
	 * loop counted n up to MAXPROP, but it should not have.*/

	 /* Commented out for less verbosity.*/
	 fprintf(stderr,"\nDebug: setting prop_id = %d",i);
	 fprintf(stderr,"\nDebug: setting n count = %d",n);
     }
   } 

   /* Spit out the hptr->prop_id to verify correct ordering. */
   for ( i = 0; i < hptr->nprops; ++i) {
     fprintf(stderr,"\nIndex %d prop_id %d",i,hptr->prop_id[i]);
   }

   
   nz = cdf.nz; 
    
   /* Is this gridnode at the same position as lat/lon requested? */
   
   get_indices(&cdf, hptr->lat, hptr->lon, &row, &col);
   get_lat_lon(&cdf, row, col, &lat, &lon);
   hdg = vector((double)lat, (double)lon, pos->lat, pos->lon, &phi, &dist);

   if (dist > 2 ) {   /* pt does not coincide with gridnode  */
   
	/* find another gridpt and interpolate the two */
	
	switch (pos->type) {
	   case 0:         /* longitude stays the same */
	      if (hdg > 90 && hdg < 180 ) {
		   lat = lat - cdf.yincr;
		   break;
              }
	      lat = lat + cdf.yincr;
	      break;
	   
	   case 1:        /* latitude stays the same */
	      if (hdg > 180) {
		   lon = lon - cdf.xincr;
		   break;
              }
	      lon = lon + cdf.xincr;
	      break;
	   
	   default:   /* unlikely case where endpt doesn't fall 
	                 on a gridnode.  */
	      
	      if ((hdg > 45 && hdg < 135) || (hdg > 225 && hdg < 315)) {
	         if (hdg > 180) {
		   lon = lon - cdf.xincr;
		   break;
                 }
	         lon = lon + cdf.xincr;
	         break;
              }
	      
	      if (hdg > 90 && hdg < 180 ) {
		   lat = lat - cdf.yincr;
		   break;
              }
	      lat = lat + cdf.yincr;
	   
	} /* end switch */
	
	/* first check currently open cdf file for match */
	
        found = check_bounds(&cdf, lat, lon);
	
	curfile = 0;
	while (!found && ++curfile <= nfiles) {
           cdf_close(cdfid);
           cdfid = cdf_open(dir, arglist[curfile], extent, print_msg);
           if (error = read_cdf_hdr(cdfid, &cdf)) 
               exit (1);
           found = check_bounds(&cdf, lat, lon);
	}   
	
	if (found) {

          if (cdf.nz != nz) {
	     fprintf(stderr,"\nFATAL ERROR:  mismatch of standard levels in cdf files.\n");
	     exit(1);
	  }	
	  dptr2 = (struct HYDRO_DATA *) get_memory((void *)NULL, 1, sizeof(struct HYDRO_DATA));
	
          np2 = load_properties(cdfid, &cdf, dptr2, lat, lon);

          /* assign relative weights to each profile inversely proportional
	     to distance */
	  	   
          get_indices(&cdf, lat, lon, &row, &col);
          get_lat_lon(&cdf, row, col, &lat, &lon);
          hdg = vector((double)lat, (double)lon, pos->lat, pos->lon, &phi, &dist2);

          weight1 = dist2 / (dist + dist2);
	  weight2 = dist / (dist + dist2);
	  
	  /* find deepest common depth and stdlevel */

          	  
	  dcl = (int) dptr->observ[(int)DE][nz-1];
	  np2 = 0;   /* keep track of which profile */
	  if (dcl < 0) {
	     dcl = dptr2->observ[(int)DE][nz-1];
	     np2 = 1;
	  }
          else {
	     if (dptr2->observ[(int)DE][nz-1] > 0)
	        if (dcl > dptr2->observ[(int)DE][nz-1]) {
	           dcl = (int) dptr2->observ[(int)DE][nz-1];
		   np2 = 1;
		}
	  }

          if (dcl < 0)
	    np2 = 2;	  
	  i = 0;
	  while ((i < nz-1)  && (std_depth[i] < dcl)) {
	     ++i;
	  }
	  
	  dcl = i;
	   
	   /* interpolate each level except bottom to produce a single profile */
	   
          for (i = 0; i < hptr->nprops; ++i) {
	    ip = hptr->prop_id[i];
	    if (dptr2->observ[ip] != NULL ) {
	    
               for (j = 0; j < nz-1; ++j) {
	         if (dptr->observ[ip][j] <= flag) 
		     dptr->observ[ip][j] = dptr2->observ[ip][j];
		 else if (dptr2->observ[ip][j] <= flag)
		     ;
		 else 
		    dptr->observ[ip][j] = weight1 * dptr->observ[ip][j] + weight2 * dptr2->observ[ip][j];
	    
	       }
	       
	       /* set levels below deepest common level to empty */
	       for (j = dcl; j < nz-1; ++j) {
	             dptr->observ[ip][j] = (double) HBEMPTY;
	       }
	       
	       if (np2 == 1)
		  dptr->observ[ip][nz-1] = dptr2->observ[ip][nz-1];
	       
	       if (np2 == 2)
		  dptr->observ[ip][nz-1] = HBEMPTY;
	    }
	  } 

          for (i = 0; i < MAXPROP; ++i) {
	      if (dptr2->observ[i] != NULL)
	         free((void *) dptr2->observ[i]);
	  }	   
	  free((void *) dptr2);
	}
	 
   }  /* end if dist > 2 */

   cdf_close(cdfid);

   /* Prepare station for output.
    * Define an index property to check for
    * missing data.  We don't want it to be
    * pressure or depth. */
   i = 0;
   ip = hptr->prop_id[i];

   /*
   fprintf(stderr,"Debug: ip = %d\n",ip);
   fprintf(stderr,"Debug: hptr->nprops = %d\n",hptr->nprops);
   */
   
   /* Depth and pressure are at index values 0 and 1,
    * respectively, so this while loop will skip over
    * either of those two with the continue so that
    * a non-pressure or non-depth variable can be found.
    * The problem with this loop is that if only depth
    * and pressure are available in this profile
    * (empty profile) then the loop will overrun and
    * eventually pick up a spurious value for ip in
    * whatever happens to be in memory at
    * hptr->prop_id[i]. Instead,
    * we need to stop the while loop explicitly.*/
   while ( (ip == (int)DE) || (ip == (int)PR) ) {
     fprintf(stderr,"Debug: i = %d\n",i);
     
     /* Increment a counter and test that it is equal
      * to the number of properties in this station.
      * If this happens, we should exit the loop, but
      * here it just continues the loop to the next
      * iteration.  So the fix is here.  I comment
      * out the continue and replace it with a check
      * for whether anything other than DE and PR
      * has been found.
      * If it is not equal to the number of properties,
      * skip to top of the loop.*/
     if (++i == hptr->nprops) {
       /* SFG
       continue;
       */

       /* We should end this loop as we're
	* about to overrun the data otherwise.
	* Get whatever index we have. */
       ip = hptr->prop_id[i];

       /* Check if that index is still DE or PR
	* and if it is, then we return(0) to
	* mark an empty profile.*/
       if ( (ip == (int)DE) || (ip == (int)PR) ) {
	 fprintf(stderr,"Cannot find any variables other than DE and PR. Skipping!");
	 return(0);
       } else {
	 break;
       }  
     }
     /* SFG added the brackets around the if() {} above.
      * Counter should not overrun nprops. The goal is
      * to just skip over DE and PR and get to the first
      * active property.*/
     
      
     /* We still have some properties to scan through,
      * set ip to the next property id in the list.
      * We are left with the ID of the last property
      * in the station.  The current problem is that
      * this value is being set to 49, which is beyond
      * the maximum allowed number of properties, so
      * below we throw a seg. fault because there's no
      * allocated space.*/
     ip = hptr->prop_id[i];
     
     fprintf(stderr,"Debug: ip = %d\n",ip);
   }

   /* Remove flagged levels from data arrays */    
   n = 0;

   fprintf(stderr,"\n Debug: n = %d\n",n);
   fprintf(stderr,"Debug: nz = %d\n",nz);
   fprintf(stderr,"Debug: ip = %d\n",ip);
   
    for (j = 0; j < nz; ++j) {

      /* Check over each depth level for a bad value. */
      if (dptr->observ[ip][j] > flag) {

	/* For good values, copy over the good data and
	 * augment a the good data counter (here, n). */
	for (i = 0; i < hptr->nprops; ++i) {
	  dptr->observ[hptr->prop_id[i]][n] = dptr->observ[hptr->prop_id[i]][j];
	}
	++n;
      }

      /* Bad values are simply skipped. */
    }

    hptr->nobs = n;
    hptr->pdr = dptr->observ[(int)DE][n-1] + 10;  
    dptr->nobs = n;
    dptr->nprops = hptr->nprops;

    return (1);
} /* end get_profile() */

/****************************************************************************/
int check_bounds(struct CDF_HDR *cdfptr, float lat, float lon)

      /* Returns 1 if this point falls into the cdf bounds
         or 0 if not */
      
{
   float xoff, yoff;
      

      xoff = yoff = 0.0;      
      if (!cdfptr->node_offset) {
         xoff = 0.5 * cdfptr->xincr;
         yoff = 0.5 * cdfptr->yincr;
      }

     /* within bounds? */


      if ((lat < cdfptr->ymin -yoff) || (lat >= cdfptr->ymax) )   /* lat NOT in bounds */
          return (0);    
	            
      if ((lon >= (cdfptr->xmin - xoff))
        && (lon < (cdfptr->xmax + xoff)) )
        return (1);                                          /* lat,lon both in bounds */
	
	/* make sure longitude ranges are comparable */
	
       if (lon >= 0.) {          /* positive longitude */
      
          if (cdfptr->xmax < 0) {  /*lon range in cdf file is negative, make it pos */
             cdfptr->xmax += 360;
             cdfptr->xmin += 360;
	  }
	  else {   /* cdf file crosses greenwich, make lon neg and compare */
	      lon -= 360;
	  }
	  
         if ((lon >= (cdfptr->xmin - xoff)) && (lon < (cdfptr->xmax + xoff)) )
               return (1);                                          /* lat,lon both in bounds */
	  
	  return (0);  
      }
      
      /* negative longitude */
         if (cdfptr->xmin >= 0) {  /* lon range in cdf file is positive, make it neg */
             cdfptr->xmin -= 360;
             cdfptr->xmax -= 360;
	 }    
          else {   /* make lon positive and compare */
	     lon += 360;
	  }
	  
         if ((lon >= (cdfptr->xmin - xoff)) && (lon < (cdfptr->xmax + xoff)) )
               return (1);                                          /* lat,lon both in bounds */
	  
	  return (0);  
        

}  /* end check_bounds() */
/****************************************************************************/

int load_properties(int cdfid, struct CDF_HDR *cdfptr, struct HYDRO_DATA *dptr, float lat, float lon)

/*  Loads properties from cdf_file into dptr for gridpt nearest lat/lon.
   Returns number of properties, including depth. */

{
   int *prop_avail, error;
   int i, j, n, nprops;
   float *x;
   int row, col;
   
   
  prop_avail = (int *) calloc((size_t)MAXPROP, sizeof(int));	

  /*
  fprintf(stderr,"\nDebug: cdfptr->nprops = %d",cdfptr->nprops);
  */
  
  for (i = 0; i < cdfptr->nprops; ++i) {
    prop_avail[get_prop_indx(cdfptr->prop_id[i])] = 1;

    /*
    fprintf(stderr,"\nDebug: cdfptr->prop_id[i] = %s",cdfptr->prop_id[i]);
    */
  }
  
  get_indices(cdfptr, lat, lon, &row, &col);
  
  n = cdfptr->nz;
  x = (float *) get_memory((void *)NULL, n, sizeof(float));
  free_and_alloc(&dptr->observ[(int)DE], n);
  
  error = read_cdf_depths(cdfid, x);
  error = read_cdf_bottom_depth(cdfid, &x[n-1], row, col, 0);       
  for (j = 0; j < n; ++j) {
     dptr->observ[(int)DE][j] = (double) x[j];
     std_depth[j] = (double) x[j];
  }
  
  nprops = 1;
  for (i=0; i < MAXPROP; ++i) {
    if (prop_avail[i] && (prop_req[i] || prop_needed[i])) {
       free_and_alloc(&dptr->observ[i], n);
       
       error = read_cdf_prop(cdfid, get_prop_mne(i), x, row, col, 0, 0, n);
       for (j = 0; j < n; ++j) {
           if (is_flagged(x[j], cdfptr->fill_value) || is_flagged(x[j], (float) HBMASK))
	        x[j] = (float) HBEMPTY;
           dptr->observ[i][j] = (double) x[j];

	   fprintf(stderr,"\n%s data value = %f",cdfptr->prop_id[nprops],dptr->observ[i][j]);
       } 
       ++nprops;      
    }
  }
  
 free(prop_avail);
 free(x);
 return(nprops);
 
} /* end load_properties()*/

/**********************************************************************************/
void compute_props(struct HYDRO_HDR *hptr, struct HYDRO_DATA *dptr)

   /* For all properties in the global array prop_req, computes it 
     if necessary, then frees up any unnecessary property. */
{
   int i, j, n, ratio_done;
   int main_props_avail;
   double dlat;
   double *e;
   
/* determine if pr, de, te, and sa are available ... */

    main_props_avail = 1;
    if (!(available(PR, hptr) && available(DE, hptr) && available(TE, hptr) 
              && available(SA, hptr))) {
         fprintf(stderr,"\n>>>>> WARNING!!!  ");
         fprintf(stderr,"Station does not include pr, de, te, and sa.");
         main_props_avail = 0;
    }
   
   ratio_done = 0;
   
/* !**! Special cases for individual properties... */
   
   for (i = 0; i < MAXPROP; ++i) {
      if (prop_req[i] && !available((enum property)i, hptr)  && main_props_avail) {
      
           switch ((enum property) i) {
             case OX:
                 if (available(O2, hptr)) {
                    free_and_alloc(&dptr->observ[i], hptr->nobs);
                    for (j=0; j < hptr->nobs; ++j) {
                       dptr->observ[i][j] = ox_kg2l(dptr->observ[(int)O2][j], dptr->observ[(int)PR][j], dptr->observ[(int)TE][j],dptr->observ[(int)SA][j]);
                    }
                 }
               break;
             case O2:
                 if (available(OX, hptr)) {
                    free_and_alloc(&dptr->observ[i], hptr->nobs);
                    for (j=0; j < hptr->nobs; ++j) {
                       dptr->observ[i][j] = ox_l2kg(dptr->observ[(int)OX][j], dptr->observ[(int)PR][j], dptr->observ[(int)TE][j],dptr->observ[(int)SA][j]);
                    }
                 }
               break;
             case TH:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_theta(hptr->nobs, dptr->observ[i], dptr->observ[(int)PR], dptr->observ[(int)TE],dptr->observ[(int)SA] );
               break;
             case TP:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_dtdp(dptr->observ[i], dptr->observ[(int)PR], dptr->observ[(int)TE],dptr->observ[(int)SA], hptr->nobs, window, w_incr);
               break;
             case TZ:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
		 dlat = (double) hptr->lat;
                 compute_dtdz(dptr->observ[i], dptr->observ[(int)PR], dptr->observ[(int)TE], dptr->observ[(int)SA], hptr->nobs,window,w_incr,dlat);
               break;
             case HC:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_heat_capacity(hptr->nobs,dptr->observ[i],dptr->observ[(int)PR],dptr->observ[(int)TE],dptr->observ[(int)SA] );
               break;
             case S0:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_sigma(0., hptr->nobs, dptr->observ[i], dptr->observ[(int)PR], dptr->observ[(int)TE],dptr->observ[(int)SA]);
               break;
             case S1:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_sigma(1000., hptr->nobs, dptr->observ[i], dptr->observ[(int)PR], dptr->observ[(int)TE],dptr->observ[(int)SA]);
               break;
             case S2:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_sigma(2000., hptr->nobs, dptr->observ[i], dptr->observ[(int)PR], dptr->observ[(int)TE],dptr->observ[(int)SA]);
               break;
             case S3:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_sigma(3000., hptr->nobs, dptr->observ[i], dptr->observ[(int)PR], dptr->observ[(int)TE],dptr->observ[(int)SA]);
               break;
             case S4:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_sigma(4000., hptr->nobs, dptr->observ[i], dptr->observ[(int)PR], dptr->observ[(int)TE],dptr->observ[(int)SA]);
               break;
             case S_:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_sigma(s_pref, hptr->nobs, dptr->observ[i], dptr->observ[(int)PR], dptr->observ[(int)TE],dptr->observ[(int)SA]);
               break;
             case HT:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_height(hptr->nobs, dptr->observ[(int)PR], dptr->observ[(int)TE],dptr->observ[(int)SA], ht_pref, dptr->observ[i]);
               break;
             case PE:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_energy(hptr->nobs, dptr->observ[(int)PR], dptr->observ[(int)TE],dptr->observ[(int)SA], pe_pref, dptr->observ[i]);
               break;
             case SV:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_sp_vol( hptr->nobs, dptr->observ[i], dptr->observ[(int)PR], dptr->observ[(int)TE],dptr->observ[(int)SA]);
               break;

             case VA:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_svan( hptr->nobs, dptr->observ[i], dptr->observ[(int)PR], dptr->observ[(int)TE],dptr->observ[(int)SA]);
               break;
             case DR:
	     case AL:
	     case BE:
	         if (! ratio_done) {
                    free_and_alloc(&dptr->observ[(int)DR], hptr->nobs);
                    free_and_alloc(&dptr->observ[(int)AL], hptr->nobs);
                    free_and_alloc(&dptr->observ[(int)BE], hptr->nobs);
                    compute_ratio( hptr->nobs, dptr->observ[(int)DR], dptr->observ[(int)PR], dptr->observ[(int)TE],dptr->observ[(int)SA], dptr->observ[(int)AL], dptr->observ[(int)BE]);
		    ratio_done = 1;
		  }
               break;

             case VS:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_sound_vel( dptr->observ[i], dptr->observ[(int)PR], dptr->observ[(int)TE],dptr->observ[(int)SA], hptr->nobs);
               break;
             case PV:
               free_and_alloc(&dptr->observ[i], hptr->nobs);
               buoy_freq(dptr->observ[i], dptr->observ[(int)PR], dptr->observ[(int)TE],dptr->observ[(int)SA], hptr->nobs, window, w_incr);
	       
	       dlat = (double) hptr->lat;
	       
	       e = (double *) calloc(hptr->nobs, sizeof(double));
	       
	       for (j = 0; j < hptr->nobs; ++j) {
	         e[j] = dptr->observ[i][j];
	         if (e[j] > -999.0)
		     e[j] *= e[j];
	       }
	       
               po_vort(dptr->observ[i], e, hptr->nobs, dlat);
	       free((void *)e);
               break;

             case RR:
               free_and_alloc(&dptr->observ[i], hptr->nobs);
	       dlat = (double) hptr->lat;
	       compute_approx_rossby_radius(dptr->observ[i], hptr->nobs, hptr->pdr, dptr->observ[(int)DE], dptr->observ[(int)PR], dptr->observ[(int)TE], dptr->observ[(int)SA], dlat, window, w_incr);
               break;

	       
             case BF:
               free_and_alloc(&dptr->observ[i], hptr->nobs);
               buoy_freq(dptr->observ[i],  dptr->observ[(int)PR], dptr->observ[(int)TE],dptr->observ[(int)SA], hptr->nobs, window, w_incr);
               break;

             case GN:
	         if (!prop_req[(int)GE]) {
                   free_and_alloc(&dptr->observ[GN], hptr->nobs);
	           compute_gamma_n(&ginfo, hptr->nobs, dptr->observ[GN], 
		       dptr->observ[(int)PR], dptr->observ[(int)TE],                                   dptr->observ[(int)SA], (double) hptr->lon, 
		       (double) hptr->lat);
		 }
	         break;
	       
             case GE:
                 free_and_alloc(&dptr->observ[(int)GE], hptr->nobs);
                 free_and_alloc(&dptr->observ[(int)GN], hptr->nobs);
	         compute_gamma_nerr(&ginfo, hptr->nobs, dptr->observ[(int)GN],   		     dptr->observ[(int)GE],  dptr->observ[(int)PR],
 		     dptr->observ[(int)TE], dptr->observ[(int)SA], 
		     (double) hptr->lon, (double) hptr->lat);
	         break;
             default:
               break;
          } /* end switch */
      }  /* end if */
   } /* end for i */
   
   /* Now free up any properties that were not requested. 
      Adjust the header appropriately. */
   
   for (i = 0; i < MAXPROP; ++i) {
       if (dptr->observ[i] != NULL) {
         if (!prop_req[i]) {
	    free((void *)dptr->observ[i]);
	    dptr->observ[i] = NULL;
	 }
       }
   } 

   n = 0;   
   for (i = 0; i < MAXPROP; ++i) {
       if (dptr->observ[i] != NULL) {
          ++n;
       }
   } 

   hptr->nprops = dptr->nprops = n;
   free((void *) hptr->prop_id);
   hptr->prop_id = (int *) calloc(n, sizeof(int));
   n = 0;
   for (i = 0; i < MAXPROP; ++i) {
       if (dptr->observ[i] != NULL) {
            hptr->prop_id[n++] = i;
       }
   }
   
   return;
}  /* end compute_props() */




