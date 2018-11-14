/* hb_toposlice.c
................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                    author:  Ruth Curry
                             Woods Hole Oceanographic Institution
                             May 2001
................................................................................
..................................................................  
.  Reads binary gridded bathymetry values consisting of 2-byte integers
.  and writes out lon/lat/z triplets  corresponding to points along
.  pathway. Distance is an optional 4th column. For a pathway that is not 
.  zonal or meridional, the output will not be equally spaced, but instead 
.  correspond to points at the intersection of the pathway and the gridlines.
..................................................................  
*/ 
#include <stdio.h> 
#include <stdlib.h>   /*sfg added to deal with exit and malloc warnings*/
#include <math.h>
#include <string.h>
#include "hb_memory.h"
#include "hb_paths.h"
#include "hydrobase.h"

#define   NINT(x)       (((x) < 0) ? ((int)((x) - 0.5)) : ((int)((x) + 0.5)))
#define    ABS(x)       (((x) < 0) ? -(x) : (x))

#define PI     3.141592654
#define  RADperDEG 0.017453292   /* pi/180 */
#define  EarthRadius  3437.747    /* in nm */
#define  KMperNM  1.853248652

struct POS_REC {
   double lat, lon, dist;
   int type;
};


/* globally defined variables */


double xmin = 0.0;        /* input grid bounds assumed */
double xmax = 360.0;
double ymin = -90.0;
double ymax = 90.0;
double xincr, yincr;
int nrowsin, ncolsin;  
int lon0to360, xgreenwich;
int  depth_is_neg; 
   
/*  prototypes for locally defined functions */

void print_usage(char *);
short **get_topo(FILE *);
int vector(double, double, double, double, double *, double *);
int compare_points(const void *, const void *);
double getlat(double, double, double, double, double *);
double getlon(double, double, double, double, double *);
double get_dist(double, double, double, double);
double get_depth(struct POS_REC *, short **, double *, double *);
           
int main (int argc, char **argv)
{ 
   int i, j, n, npoints;
   int i_flag, t_flag, dflag, kflag, lflag, pflag;
   int mflag;
   int error, nalloc;
   int depth_is_neg;   
   int hdg, zonal, meridional;
   int startpoint, quadrant;
   int xparm;
   short int **zin;
   double zout, lat, lon, x;
   double phi, dist;
   double lat0, lon0, lonmin, lonmax;
   double startlat, endlat, startlon, endlon;
   double prevlat, prevlon, cumdist;
   char *s;
   FILE *outfile;
   FILE *topofile;
   FILE *posfile;
   struct POS_REC *point;

   if (argc < 1 ){
      print_usage(argv[0]);
   }
   
/* initialize these ... */

   t_flag = 0;
   i_flag = 0;
   pflag = 0;
   dflag = kflag = lflag = mflag = 0;
   depth_is_neg = 1; 
   xincr = yincr = 0.1;
   error = 0;  
   outfile = stdout;
   posfile = stdin;

   
/* parse the command line arguments */

   for (i = 1; i < argc; i++) { 
      if (argv[i][0] == '-') {
         s = &argv[i][1]; 
         switch (*s) { 

            case 'D':
	       dflag = 1;
	       ++s;
	       if (*s == 'k' || *s == 'K')
	          kflag = 1;
	       break;
	       
            case 'I':
               i_flag = 1;
               ++s;
               if (*s == '/')
                  ++s; 
               error += (sscanf(s,"%lf", &xincr) != 1); 
               s = strchr(s,'/'); /* check for another delimiter*/ 
               yincr = xincr; 
               if (s != NULL) { 
                  ++s; /* move past delimiter */ 
                  error += (sscanf(s,"%lf", &yincr) != 1); 
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
	  
           case 'M':
	       mflag = 1;
               ++s;
               error += (sscanf(s,"%d", &xparm) != 1);
	       error += (xparm < 1) || (xparm > 3); 
	       break;
              
           case 'N' :
               depth_is_neg = 0;
               break;
	       
           case 'O':
               outfile = fopen(&argv[i][2],"w");
               if (outfile == NULL) {
                  fprintf(stderr,"\nError opening %s for output\n", &argv[i][2]);
                  exit(1);
               }
               break;
           case 'P':
	       pflag = 1;
	       break;
               
           case 'T':
               topofile = fopen(&argv[i][2],"r");
               if (topofile == NULL) {
                  fprintf(stderr,"\nError opening %s for input\n", &argv[i][2]);
                  exit(1);
               }
               t_flag = 1;
               fprintf(stderr,"\nOpened %s \n", &argv[i][2]);
               break;
               
            case 'h': 
               print_usage(argv[0]);
               exit(0);
               
            default:  
               error = 1;
               
         } /* end switch */
             
      }
      else  {
        error = 1;
      }
      
      if (error ) { 
           fprintf(stderr,"\nError parsing command line args.\n");
           fprintf(stderr," in particular:  '%s'\n", argv[i]); 
           fprintf(stderr,"\nType %s -h for help\n", argv[0]); 
           exit(1); 
      }
   } /* end for */

   if (! t_flag ) {
   
      topofile = fopen(BATHPATH,"r");
      if (topofile == NULL) {
          fprintf(stderr,"\nUnable to open default topography file\n", BATHPATH);
          fprintf(stderr,"You can specify the name of a binary topography file for input with [-T]\n"); 
          exit(1);
      }
    
   } 
   if (! i_flag) {
      fprintf(stderr,"Increment for topofile: [%.2lf/%.2lf] \n", xincr, yincr);
   } 
   
    if (! lflag ) {
      fprintf(stderr,"Enter longitude latitude pairs ....\n  "); 
   } 

   nrowsin = (int) (NINT((ymax - ymin) / yincr)) + 1;
   ncolsin = (int) (NINT((xmax - xmin)/ xincr)) + 1; 
   
  

/*--------------------------------------------*/    
/* Read in first  lat/lon point */

   nalloc = 10000;

   point = (struct POS_REC *) get_memory((void *)NULL, (size_t)nalloc, sizeof(struct POS_REC));
   
   if (lflag)
       n = fscanf(posfile, "%lf%lf",  &startlon, &startlat);
   else 
       n = fscanf(stdin, "%lf%lf\n",  &startlon, &startlat);
       
   
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
   
   if (startlon > lonmax) lonmax = startlon;
   if (startlon < lonmin) lonmin = startlon;

/*--------------------------------------------*/    
/* Get consecutive points, construct an array of gridline crossings between
   each pair of points */
          
   npoints = 1;
   startpoint = 0; 
   if (lflag)
      n = fscanf(posfile, "%lf%lf",  &endlon, &endlat);
   else
      n = fscanf(stdin, "%lf%lf\n",  &endlon, &endlat);
       
   while  (n == 2) {
      
      if (endlon > lonmax) lonmax = endlon;
      if (endlon < lonmin) lonmin = endlon;
     hdg = vector(startlat, startlon, endlat, endlon, &phi, &dist);
     
     
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
     
       /* get latitudinal gridline xings */
       
       lat0 = (double) NINT(startlat);
       if (quadrant == 0 || quadrant == 3) {
       
                /*heading north */
	 
	    if (lat0 > startlat)  /* find nearest whole degree */
	      lat0 -=  1;
	      
	      
	       /* find first gridline past startlat */  
	    while ((lat0 += yincr) <= startlat)  
	        ;  
		
		/* find all lat gridline xings until endlat */ 
	    while (lat0 < endlat) {
	       point[npoints].lat = lat0;
	       point[npoints].lon = getlon(lat0, startlon, startlat,  phi, &point[npoints].dist);
	       point[npoints].type =  1;
	       if (++npoints > nalloc) {
	          nalloc += 100;
		  point = (struct POS_REC *) get_memory((void *)point, (size_t) nalloc, sizeof(struct POS_REC));
	       }
	       
	       lat0 += yincr;
	    }
	    
	}  /* end if */
	
	
	else {                 /* heading south */
	 
	    if (lat0 < startlat)  /* find nearest whole degree */
	      lat0 +=  1;
	      
	      
	       /* find first gridline past startlat */  
	    while ((lat0 -= yincr) >= startlat)  
	        ;  
		
		/* find all lat gridline xings until endlat */ 
	    while (lat0 > endlat) {
	       point[npoints].lat = lat0;
	       point[npoints].lon = getlon(lat0, startlon, startlat, phi, &point[npoints].dist);
	       point[npoints].type =  1;
	       if (++npoints > nalloc) {
	          nalloc += 100;
		  point = (struct POS_REC *) get_memory((void *)point, (size_t) nalloc, sizeof(struct POS_REC));
	       }
	       
	       lat0 -= yincr;
	    }
	   
       } /* end else */ 
     
     } /* end if !zonal */
     
     
/*--------------------------------------------*/    
     if (!meridional) { 
     
       /* get longitudinal gridline xings */
     
       lon0 = (double) NINT(startlon);
       if (quadrant == 0 || quadrant == 1) {
       
                /* heading east */
	 
	    if (lon0 > startlon)  /* find nearest whole degree */
	      lon0 -=  1;
	      
	      
	       /* find first gridline past startlon */  
	    while ((lon0 += xincr) <= startlon)  
	        ;  
		
		/* find all lon gridline xings until endlon */ 
	    while (lon0 < endlon) {
	       point[npoints].lon = lon0;
	       point[npoints].lat = getlat(lon0, startlon, startlat, phi, &point[npoints].dist);
	       point[npoints].type =  0;
	       if (++npoints > nalloc) {
	          nalloc += 100;
		  point = (struct POS_REC *) get_memory((void *)point, (size_t) nalloc, sizeof(struct POS_REC));
	       }
	       
	       lon0 += xincr;
	    }
	    
	}  /* end if */
	
	
	else {                 /* heading west */
	 
	    if (lon0 < startlon)  /* find nearest whole degree */
	      lon0 +=  1;
	      	      
	       /* find first gridline past startlon */  
	    while ((lon0 -= xincr) >= startlon)  
	        ;  
		
		/* find all lon gridline xings until endlon */ 
	    while (lon0 > endlon) {
	       point[npoints].lon = lon0;
	       point[npoints].lat = getlat(lon0, startlon, startlat, phi, &point[npoints].dist);
	       point[npoints].type =  0;
	       if (++npoints > nalloc) {
	          nalloc += 100;
		  point = (struct POS_REC *) get_memory((void *)point, (size_t) nalloc, sizeof(struct POS_REC));
	       }
	       
	       lon0 -= xincr;
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
      
      /* sort the array of points in order by increasing dist from startpoint */
            
      qsort((void *) &point[startpoint], npoints - startpoint, sizeof(struct POS_REC), compare_points);
      
      
      startlat = endlat;
      startlon = endlon;
      startpoint = npoints;
      
      if (lflag)
         n = fscanf(posfile, "%lf%lf",  &endlon, &endlat);
      else
         n = fscanf(stdin, "%lf%lf\n",  &endlon, &endlat);
	 
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
   
   point = (struct POS_REC *) get_memory((void *)point, (size_t)npoints, sizeof(struct POS_REC));

/*--------------------------------------------*/    
/*        end of defining pathway */
/*--------------------------------------------*/    

   lon0to360 = 1;
   if (lonmin < 0)
      lon0to360 = 0;
   xgreenwich = 0;
   xgreenwich = (!lon0to360 && lonmax >= 0.);
   
 /* read in topography file ... */

    zin = get_topo(topofile);
    fclose(topofile);
  
      
    fprintf(stderr,"\nWriting out values.....\n");

    prevlat = point[0].lat;
    prevlon = point[0].lon;
    cumdist = 0.0;
    if (xparm == 3) 
       dflag = 1;
       
    for (i = 0; i < npoints; ++i) {
      zout = get_depth(&point[i], zin, &lat, &lon);
      if (depth_is_neg)
	zout = -zout;
	   
      if (pflag) 
	zout = hb_p80( zout, lat);
	     
      if (dflag)  {
	dist = get_dist(prevlat, prevlon, lat, lon);
	if (kflag)
	  dist = dist * KMperNM;
	cumdist += dist;
	prevlat = lat;
	prevlon = lon;
      } 
      if (mflag) {
	
	x = lon;
	if (xparm == 2)
	  x = lat;
	if (xparm == 3)
	  x = cumdist;
	
	if (i == 0) {
	  if ( depth_is_neg == 0 ) {
	    fprintf(outfile,">I\n%8.3lf -10000\n", x);
	  }
	  else {
	    fprintf(outfile,">I\n%8.3lf  10000\n", x);
	  }
	}
	fprintf(outfile,"%8.3lf %8.0lf\n", x, zout);
	
      }
      else {
	fprintf(outfile,"%8.3lf %8.3lf %8.0lf", lon, lat, zout);
	if (dflag)  
	  fprintf(outfile," %10.2lf", cumdist);
	fprintf(outfile,"\n");
      }
    }
    
    if (mflag) {
      x = lon;
      if (xparm == 2)
	x = lat;
      if (xparm == 3)
	x = cumdist;
      
      if ( depth_is_neg == 0 ) {
	fprintf(outfile,"%8.3lf -10000\n", x);
      }
      else {
	fprintf(outfile,"%8.3lf  10000\n", x);
      }
    }
     
    fprintf(stderr,"\n End of %s.\n", argv[0]);
    exit(0);
    
} /* end main */


/****************************************************************************/

void print_usage(char *program) 
{ 
   fprintf(stderr,"\n%s outputs bathymetry along a pathway specified", program);
   fprintf(stderr,"\nas a list of lon/lat positions. The default bathymetry database");
   fprintf(stderr,"\nis %s ", BATHPATH);
   fprintf(stderr,"\nAn alternative file may be specified with the -T option. ");
   fprintf(stderr,"\nIf the input grid spacing is other than %.2lf, specify grid increment",xincr);
   fprintf(stderr,"\nwith -I.  Pathway must be a minimum of 2 points, but can be");  
   fprintf(stderr,"\nmultiple points. -L specifies filename containing these points.");  
   fprintf(stderr,"\nIf no posfile is specified, the points are expected from stdin.");
   fprintf(stderr,"\nLon/lat/z values are written to outfile or stdout.");
   fprintf(stderr,"\n -D causes the distance along the pathway to be written"); 
   fprintf(stderr,"\nout as a fourth column.");
   fprintf(stderr,"\nTo generate an xy mask file for use with hb_fitsection,");
   fprintf(stderr,"\nspecify -M<xval> where xval is 1, 2, or 3 (lon, lat or distance)");
   fprintf(stderr,"\n\nUsage:  %s  [-T<global_topo_file>] [-I<xincr[/yincr]>][-D[k]] [-L<posfile>] [-M<xval>] [-N] [-O<outfile>] [-P]", program);
   fprintf(stderr,"\n\n   OPTIONS:"); 
   fprintf(stderr,"\n[-D] : include distance (nautical miles) in output. "); 
   fprintf(stderr,"\n       Append k for distance in km [-Dk]. "); 
   fprintf(stderr,"\n[-I] : specify grid increments of input topo file. "); 
   fprintf(stderr,"\n        ex: -I.1/.5   Default increment is  [%.2lf/%.2lf]", xincr, yincr); 
   fprintf(stderr,"\n[-L] : position file contains lon/lat pairs (one per line) to define pathway"); 
   fprintf(stderr,"\n[-M] : outputs xy values for use as a mask file with hb_fitsection."); 
   fprintf(stderr,"\n       specify xval = 1 for longitude"); 
   fprintf(stderr,"\n                    = 2 for latitude"); 
   fprintf(stderr,"\n                    = 3 for distance"); 
   fprintf(stderr,"\n       Info is added at beginning and end for polygon masking mode."); 
   
   fprintf(stderr,"\n[-N] : make seafloor values negative. ");
   fprintf(stderr,"\n          [default is seafloor positive] ");
   fprintf(stderr,"\n[-O] : output file. If not specified output stdout"); 
   fprintf(stderr,"\n[-P] : output pressure instead of depth"); 
   fprintf(stderr,"\n[-T] : specify name of binary topography file. Default[%s]", BATHPATH);
   fprintf(stderr,"\n[-h] : help -- prints this message.");
   fprintf(stderr,"\n\n");
   return;
   
} /*end print_usage() */
   
/****************************************************************************/
/****************************************************************************/

short **get_topo(FILE *fptr)
   /* allocates memory and reads in topography values.  Returns a pointer to
      the start of the memory block. An error causes an error message to be
      printed and an exit.   */
{

   int row, n;
   short **z;
   
   /* Allocate memory */
   
   z = (short **) malloc(nrowsin * sizeof(short *));
   if (z == NULL) {
      fprintf(stderr,"\nError allocating memory.\n");
      exit(1);
   }
   
   for (row = 0; row < nrowsin; ++row) {
     z[row] = (short *) malloc(ncolsin * sizeof(short));
     if (z[row] == NULL) {
         fprintf(stderr,"\nError allocating memory.\n");
         exit(1);
      }
   }
   
   fprintf(stderr,"\nReading in topography values ");
   
   for (row = 0; row < nrowsin; ++row) {
     n = fread((short *)z[row], sizeof(short), ncolsin, fptr);
     if (n != ncolsin) {
         fprintf(stderr,"\nError reading the topofile at row %d\n", row);
         exit(1);
     }
     if ((row % 10) == 0)   /* signal progress */
        fprintf(stderr,".");
   }
   
   fprintf(stderr,"\nFinished reading topofile.\n");
   
   return(z);

} /* end get_topo() */
/****************************************************************************/
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
     dx = (lon2 - lon1)* RADperDEG * cos(RADperDEG * 0.5*(lat2+lat1) ) * EarthRadius;
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
   
   return(dx / (EarthRadius * cos(RADperDEG*0.5*(lat1+lat2)) * RADperDEG) + lon1);
  
}
/****************************************************************************/
double get_depth(struct POS_REC *pos, short **zin, double *latptr, double *lonptr)	 {
   int i, n, row, col;
   int hdg;
   double lon0, lat, lon;
   double phi, dist, dist2, zout, z2;
   double weight1, weight2;

   lon0 = pos->lon;
   if (!lon0to360) {          /* make lon0 be in range 0 -> 360 */
      if (pos->lon < 0)
         lon0 += 360;
   }
   
   /* find nearest point in zin grid */
   
   row = NINT((pos->lat + .0001 - ymin ) / yincr);
   col = NINT((lon0 + .0001 - xmin) / xincr);
   zout = (double) zin[row][col];
   lat = ymin + row * yincr;
   lon = xmin + col * xincr;
   
   if (xgreenwich) {    /* make lon range -180 -> 180 */
       if (lon > 180)
           lon -= 360.;
       if (lon0 > 180)
           lon0 -= 360.;
   }
   
   /* Does the gridnode coincide with the lat/lon point?*/   
   
   hdg = vector(lat, lon, pos->lat, lon0, &phi, &dist);
   
   if (dist > 1 ) {   /* pt does not coincide with gridnode  */
   
	/* find neighboring gridpt and interpolate the two */
	
	switch (pos->type) {
	   case 0:         /* longitude stays the same */
	      if (hdg > 90 && hdg < 180 ) {
	           --row;
		   lat -=  yincr;
		   break;
              }
	      ++row;
	      lat += yincr;
	      break;
	   
	   case 1:        /* latitude stays the same */
	      if (hdg > 180) {
	           --col;
		   lon -= xincr;
		   break;
              }
	      ++col;
	      lon += xincr;
	      break;
	   
 	   default:   /*  endpt doesn't fall  on a gridnode.  */
	      if ((hdg > 45 && hdg < 135) || (hdg > 225 && hdg < 315)) {
	         if (hdg > 180) {
		   lon -= xincr;
		   --col;
		   break;
                 }
	         lon += xincr;
		 ++col;
	         break;
              }
	      
	      if (hdg > 90 && hdg < 180 ) {
		   --row;
		   lat -= yincr;
		   break;
              }
	      lat += yincr;
	      ++row;
	   
	} /* end switch */
	      
        if ((row >= 0) && (row < nrowsin) && (col >= 0) && (col < ncolsin)) {
	
          hdg = vector((double)lat, (double)lon, pos->lat, lon0, &phi, &dist2);

          /* weight inversely proportional to distance from point */
          weight1 = dist2 / (dist + dist2);
	  weight2 =  dist / (dist + dist2);
          z2 = (double) zin[row][col];
	  zout = weight1 *zout + weight2 *z2;
	}
   }
   
   *latptr = pos->lat;
   *lonptr = pos->lon;
   return(zout);
}  
/****************************************************************************/
double get_dist(double prev_lat, double prev_lon, double lat, double lon)
  /* Returns distance between points in nautical miles */

{
   double dx, dy;
   double dist;
   
   
   dy = (double) ABS(lat - prev_lat);
   dx = cos((lat + prev_lat) *.5 * RADperDEG) * ABS(lon - prev_lon);
   dist = RADperDEG * EarthRadius * sqrt(dx * dx + dy * dy);  
   return(dist);

} /* end get_dist() */
