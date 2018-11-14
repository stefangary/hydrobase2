/*  hb_volumetric.c

................................................................................
                          *******  HydroBase2 *******
               ...................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
                             2004
................................................................................
____________________________________________________________________________
  USAGE:  

 hb_volumetric filename(_roots) -1mne/value[t] -2mne/value[b]  -Ixinc[/yinc] [-D<dirname>] [-E<file_extent>] [-Z<zmin/zmax>]

  list of filenames MUST be first argument

 -1 : property and value at upper level;
 -2 : property and value at lower level;
 -I : grid increments (xgrid, ygrid)

 -D : specifies directory for input files (default is current dir)  
      ex: -D/home/pete/data

 -E : specifies input_file extent (default is no extent) ex: -E.sum
 -Z : specifies depth limits
____________________________________________________________________________
Computes volumetric heat and salt content for a vertical layer defined by two surfaces. 
Weighting in the vertical is by dz. Rectangular area surrounding each profile 
is dx*dy for the specified gridspacing (xinc/yinc) at each profile's latitude.
Outputs  lon, lat, avg_theta, avg_salt, avg_rho, avg_Cp, dz, dxdy and 
heat content (avg_theta*avg_rho*avg_Cp*dx*dy*dz.
The sum of property content at all gridpoints is accumulated and on the last lines are
listed the following:
   Total area (km^2)
   Avg thickness of layer (area weighted) (km)
   Total heat content
   Avg Theta/unit volume
   Avg Salinity/unit volume 
   Avg Density 
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

#define  PI     3.141592654
#define  RADperDEG 0.017453292    /* pi/180 */
#define  EarthRadius  3437.747    /* in nm */
#define  KMperNM  1.853248652

      
				    
  /*  define a structure to store layer information */

struct LAYER_INFO {
	double lat;
	double lon;
	double avg_theta;     /* depth-averaged potential temperature for layer */    
	double avg_salt;      /* depth-averaged salinity for layer */    
	double avg_rho;       /* depth-averaged in situ density for layer */
	double avg_Cp;        /* depth-averaged specific heat for  layer */
	double dz;            /* thickness of layer */
	double dxdy;          /* area (km) surrounding each grid point (profile) */
	double Qdzdxdy;       /* volumetric heat content (theta*rho*Cp*km^3) */
        struct LAYER_INFO *next;
};				    
struct LAYER_INFO *list_ptr;

   /* global variables to store station ... */

struct HYDRO_DATA station;
struct HYDRO_HDR hdr;
double x1, x2;
double pref1, pref2;    /* ref pressures for top/bottom of layer */
int prop1, prop2;
double sp_heat[MAXOBS];   /* store specific heat for computing heat content */

struct GAMMA_NC ginfo;   /* used for neutral density */

int warnflag;
int iflag, zflag;
int check_bottom_outcrop;
int check_top_outcrop;
double zmin, zmax;
double xinc, yinc;

   /* prototypes for locally defined functions */

void    print_usage(char *);
void    get_hydro_data(int);
double get_weighted_avg(double, double, double *,double *, double *, double *, double *, double *);
void compute_prop(int, double **, double *, double *, double *, int, double, double);
int 	vector(double, double, double, double,  double *);
void output_data(FILE *);

int main (int argc, char **argv)
{
   short   opt1, opt2;
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
    opt1 = opt2 = 0;
    error = 0;
    nprops = 2;
    zmin = 0.0;
    zmax = 10000.0;
    infile = STDIN;
    iflag = zflag = 0;
    warnflag = 0;           /* set to 1 after warning is printed */
    check_bottom_outcrop = 0;
    check_top_outcrop = 0;
    list_ptr = (struct LAYER_INFO *) NULL;

    
    
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
			 st = &argv[i][4];
			 if (*st == '/')
			   ++st;
			 
                         if (sscanf(st,"%lf/%lf", &pref1, &x1) != 2)
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
                        opt2 = 1;
                        error = (sscanf(&argv[i][2],"%2s", id) == 1) ? 0 : 1;
                        if ((prop2 = get_prop_indx(id)) < 0) {
                           fprintf(stderr,"\n%.2s is not an appropriate property.", id);
                           ++error;
                        }
                        
                        if (prop2 == (int)S_ || prop2 == (int)HT || prop2 == (int)PE ) {
			 st = &argv[i][4];
			 if (*st == '/')
			   ++st;
                         if (sscanf(st,"%lf/%lf", &pref2, &x2) != 2)
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

   if (!opt1 || !opt2 ) {
       fprintf(stderr,"\nYou must define the layer with -1 and -2.\n");
       exit(1);
   }
   if (  !iflag) {
       fprintf(stderr,"\nYou must specify grid spacing (in degrees: xinc/yinc}. \n");
        exit(1);
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
   fprintf(stderr,"\n%s computes volumetric heat and salinity content over a layer specified by two surfaces", program);
   fprintf(stderr,"\nThe layer can be strictly adhered to or relaxed to incorporate values from sea surface (if surface1 outcrops above depth=0) down to the seafloor (if surface2 > value at seafloor)");
   fprintf(stderr,"\nThe layer can be further limited to a depth range using -Z<zmin/zmax>");     
   fprintf(stderr,"\nFor each gridpoint, outputs a line containing:");
   fprintf(stderr,"\nlon, lat, avg_theta, avg_salt, avg_rho, avg_Cp, dz(layer_thickness), dx*dy, heat_content ");
   fprintf(stderr,"\nThe sum of properties is accumulated and on the last lines are listed the following:");
   fprintf(stderr,"\nTotal area (km^2)");
   fprintf(stderr,"\nAvg thickness of layer (area weighted) (km)");
   fprintf(stderr,"\nTotal heat content (J)");
   fprintf(stderr,"\nAvg theta/unit volume (deg C)");
   fprintf(stderr,"\nAvg salinity/unit volume");
   fprintf(stderr,"\nAvg density/unit volume (kg/m^3");
   
   fprintf(stderr,"\n\nUsage:  %s filename_root(s)  -1<mne[pref]/value> -2<mne[pref]/value[b]> -Ixinc[/yinc] [-D<dirname>] [-E<file_extent>][-Z<zmin/zmax>]", program);

   fprintf(stderr,"\n\n  List of filenames MUST be first argument");
   fprintf(stderr,"\n  or input is expected to come from stdin...");
   fprintf(stderr,"\n    -1  : surface 1 property/value (if property = s_");
   fprintf(stderr,"\n          you need to specify pref after s_");
   fprintf(stderr,"\n          Append 't' to start with the first observation level available");
   fprintf(stderr,"\n          in the event that it is greater than the property value at the top of the layer specified.   ");
   fprintf(stderr,"\n    -2  : surface 2 property/value (if property = s_");
   fprintf(stderr,"\n          you need to specify pref after s_");
   fprintf(stderr,"\n          Append b to integrate to the bottom if this surface");
   fprintf(stderr,"\n          is deeper than the bottom.  ");
   fprintf(stderr,"\n    -I  : specify grid spacing (xinc/yinc) in degrees to compute dx and dy");
   fprintf(stderr,"\n   [-D] : specifies dirname for input files (default is current");
   fprintf(stderr,"\n          directory)  ex: -D../data/ ");
   fprintf(stderr,"\n   [-E] : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n   [-Z] : specify depth limits for layer.");
   fprintf(stderr,"\n          Exclude depths that fall outside of these limits from weighted average.");
   fprintf(stderr,"\n   [-h] : help... prints this message.");
   fprintf(stderr,"\n\n");  
       
   return;
}


/*****************************************************************************/
/*****************************************************************************/
void get_hydro_data(int file)
   /*  Reads each station in a HydroBase file and integrates property 
   over specified layer.  All potentially necessary output information are stored 
   in a linked list of records. */
{
   int error, i, j, nreq;
   int main_props_avail;
   double  p1, p2, dir, pref=0.0;
   double dx1, dx2, dy;
   struct LAYER_INFO *rec_ptr, *r1ptr;


/* read each station in file ... */

    while ((error = get_station(file, &hdr, &station)) == 0) {	  

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
           
      if ((p1 < -8.) && (check_top_outcrop))  {  /* didn't find surf 1 */
           
        switch ((enum property) prop1) {
         case S0:
         case S1:   /* if shallowest observation > value we are */
         case S2:   /* seeking, set p1 to first depth in array */
         case S3:
         case S4:
         case S_:
         case GN:
	 case PR:
	 case DE:
             if (station.observ[prop1][0] > x1) {
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
      
      if ((station.observ[prop2] == NULL) && (main_props_avail) ) {
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

      
      if ((p2 < 0) && check_bottom_outcrop) {
         if ((p1 >= 0) && (p2 < 0)) {  
           
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
             if ((x2 > station.observ[prop2][hdr.nobs -1] )
             && (station.observ[(int)DE][hdr.nobs -1] > (hdr.pdr - 100))) {
                p2 = station.observ[(int)DE][hdr.nobs -1];
             }
             break;
             default:
                ;
           }  
         }
      }
      
      if (zflag && (p2 >= 0.)) {   /* check depth limits */
          if (p2 < zmin)
	     p2 = -9999.;
	  if (p2 > zmax)
	     p2 = zmax;
      } 

      if ((p1 > -1) && (p2 > -1) && (p2 > p1)) {
      
 	
        rec_ptr = (struct LAYER_INFO *) calloc(1, sizeof(struct LAYER_INFO));
	if (rec_ptr == NULL) {
	   fprintf(stderr,"\nUnable to allocate memory in get_hydro_data()\n");
	   exit(1);
	}
	rec_ptr->lat = (double) hdr.lat;
	rec_ptr->lon = (double) hdr.lon;
	
	/* compute rho from specific volume */
	
        free_and_alloc(&station.observ[(int)SV], hdr.nobs);
	compute_sp_vol(hdr.nobs, station.observ[(int)SV], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA]);
	
	for (i=0; i < hdr.nobs; ++i) {
	   station.observ[(int)SV][i] = 1.0e8 / station.observ[(int)SV][i];
	}
	
	/* compute specific heat of seawater  */
	
        for (i = 0; i < hdr.nobs; ++i) {
	     sp_heat[i] = hb_cpsw(station.observ[(int)SA][i], station.observ[(int)TE][i], station.observ[(int)PR][i]);
	}
	
	/* compute potential temperature */
       if (!available(TH, &hdr)) {
                free_and_alloc(&station.observ[(int)TH], hdr.nobs);
                
                compute_prop((int)TH, &station.observ[(int)TH], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA], hdr.nobs, pref, (double) hdr.lat);
       }
       
          rec_ptr->avg_theta = get_weighted_avg(p1, p2, station.observ[(int)TH],  station.observ[(int)SV], sp_heat, &rec_ptr->avg_rho, &rec_ptr->avg_Cp, &rec_ptr->avg_salt);
	  
	  if (rec_ptr->avg_theta < -9.0) {
	    free(rec_ptr);
	  }
	  else {
	  
	    rec_ptr->dz = (p2 - p1)/1000.0;  /* in km */
	  
	  /* compute area of grid square (trapezoid) in km^2 at this lat/lon */
	    dir = vector(rec_ptr->lat-(0.5*yinc), rec_ptr->lon, rec_ptr->lat+(0.5*yinc), rec_ptr->lon, &dy);
	    dir = vector(rec_ptr->lat-(0.5*yinc), rec_ptr->lon-(0.5*xinc), rec_ptr->lat-(0.5*yinc), rec_ptr->lon+(0.5*xinc), &dx1);
	    dir = vector(rec_ptr->lat+(0.5*yinc), rec_ptr->lon-(0.5*xinc), rec_ptr->lat+(0.5*yinc), rec_ptr->lon+(0.5*xinc), &dx2);

            rec_ptr->dxdy = (dx1 + dx2) * 0.5 * dy;
            rec_ptr->Qdzdxdy =  rec_ptr->avg_theta * rec_ptr->avg_rho * 1.0e9 * rec_ptr->avg_Cp * rec_ptr->dz *  rec_ptr->dxdy ;
          
	  
	  	   
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
 /*  Traverse the linked list, sum the property and area weights,  
       and output appropriate information */
{
      
  struct LAYER_INFO *curptr, *nextptr;
  double s_sum, t_sum, q_sum, rho_sum, vol_sum, area_sum, cp_sum, dz_avg;
  int n;
     
  curptr = list_ptr;
  if (curptr == NULL) {
     fprintf(stderr,"\nNo values to print out\n");
     return;
  }
  nextptr = list_ptr->next;
  q_sum = 0.0;
  s_sum = 0.0;
  t_sum = 0.0;
  rho_sum = 0.0;
  vol_sum = 0.0;
  area_sum = 0.0;
   
   while (curptr != NULL) {
   
     area_sum += curptr->dxdy; 
     vol_sum += curptr->dz * curptr->dxdy;
     q_sum += curptr->Qdzdxdy;
     t_sum += curptr->avg_theta * curptr->dz * curptr->dxdy;
     s_sum += curptr->avg_salt * curptr->dz * curptr->dxdy;
     rho_sum += curptr->avg_rho * curptr->dz * curptr->dxdy;
     cp_sum += curptr->avg_Cp *  curptr->dz * curptr->dxdy;
     
      fprintf(fptr, "%8.3lf %8.3lf  %10.4lf %9.4lf %9.4lf %6.4lf %.4e %.4e %.4e \n",  curptr->lon, curptr->lat,  curptr->avg_theta, curptr->avg_salt, curptr->avg_rho, curptr->avg_Cp, curptr->dz, curptr->dxdy, curptr->Qdzdxdy);
      
      curptr = nextptr;
      if (nextptr != NULL)  nextptr = nextptr->next;
		
   } /* end while */
   
   dz_avg = vol_sum / area_sum;

   fprintf(fptr, "\nTotal area: %.4e km^2 \nAverage thickness of layer %.4e km \nHeat content per unit area: %.8e (J/m^2) \nTheta per unit volume: %10.4lf (degC/km^3) \nSalinity per unit volume: %8.4lf  \nAverage Density (rho): %10.4lf (kg/m^3)\nAverage Cp (sp heat): %10.4lf (J/kg*degC)\n", area_sum, dz_avg, q_sum/area_sum * 1.0e-06, t_sum/vol_sum, s_sum/vol_sum, rho_sum/vol_sum, cp_sum/vol_sum); 
  
 
   
} /* end output_data() */
