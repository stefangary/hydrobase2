/*  hb_layerav.c

................................................................................
                          *******  HydroBase2 *******
               ...................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
                             1996
			     updated to ANSI standard Dec 2001
................................................................................
____________________________________________________________________________
  USAGE:  

 hb_layerav filename(_roots) -1mne/value -2mne/value -P<proplist> [-Wwindow/incr]  [-D<dirname>] [-E<file_extent>]

  list of filenames MUST be first argument

 -1 : property and value at upper level;
 -2 : property and value at lower level;

 -D : specifies directory for input files (default is current dir)  
      ex: -D/home/pete/data

 -E : specifies input_file extent (default is no extent) ex: -E.sum

____________________________________________________________________________
hb_layerav computes average value of properties over a pressure interval
defined by two surfaces. Each observation contributing to the average
is weighted by the pressure over which it was observed.  Outputs 
year, month, lat, lon, prop1 .. propN .                                                  ____________________________________________________________________________
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

    /* flags to indicate which properties are to be output and which are
       needed for computation */
      
short prop_req[MAXPROP];         /* set of props requested for output */
short prop_needed[MAXPROP];      /* set of props requested || defining a
                                    surface || needed for computation */

   /* global variables to store station ... */

struct HYDRO_DATA station;
struct HYDRO_HDR hdr;
double s_pref;           /* ref lev for computing a non-standard sigma level */
double pe_pref;           /* ref lev for computing potential energy */
double ht_pref;           /* ref lev for computing dynamic height */
double x1, x2;
double pref1, pref2;
int prop1, prop2;
int window, w_incr;     /* used to specify pr window for gradient properties*/

struct GAMMA_NC ginfo;   /* used for neutral density */
int warnflag;
int yflag, mflag, lflag, iflag, tflag, zflag;
double zmin, zmax;
int check_bottom_outcrop;

   /* prototypes for locally defined functions */

void    print_usage(char *);
void    get_hydro_data(int);
int     parse_prop_list(char *);
double get_weighted_av(double, double, double *);
void compute_prop(int, double **, double *, double *, double *, int, double, double);

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
    yflag = mflag = lflag = iflag = tflag = zflag = 0;
    warnflag = 0;           /* set to 1 after warning is printed */
    check_bottom_outcrop = 0;

/* initialize these ... */

   for (i = 0; i < MAXPROP; ++i) {
      station.observ[i] = (double *) NULL;
   }
   hdr.prop_id = (int *) NULL;



/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'D':                   /* get input dir */
                        dir = &argv[i][2];
                        break;
               case 'E':                    /* get file extent */
                        extent = &argv[i][2];
                        break;

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
               case 'M':
                        mflag = 1;
                        break;
               case 'Y':
                        yflag = 1;
                        break;
               case 'L':
                        lflag = 1;
                        break;
               case 'I':
                        iflag = 1;
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
                        
               case 'T':
                        tflag = 1;
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
       fprintf(stderr,"\nYou must specify output properties with -P.\n");
       exit(1);
   }
   if (!opt1 || !opt2 ) {
       fprintf(stderr,"\nYou must define the layer with -1 and -2.\n");
       exit(1);
   }
   if (  !(mflag || yflag || lflag || iflag)) {
       fprintf(stderr,"\nWARNING: You have not specified any output parameter from {month|year|position|station_id}.  Only properties will be output.\n");
   }
   
   if (zflag) {
       fprintf(stderr,"\nUsing pressure limits %.2lf - %.2lf db", zmin, zmax);
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

   fprintf(stderr,"\n\nEnd of layerav.\n");
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nhb_layerav computes average value of properties over a pressure interval");
   fprintf(stderr,"\ndefined by two surfaces. Each observation contributing to the average");
   fprintf(stderr,"\nis weighted by the pressure over which it was observed.");      fprintf(stderr,"\nThe layer can be further limited to a pressure range using -Z<zmin/zmax>");     fprintf(stderr,"\nOutputs year, month, lat, lon, prop1 .. propN. ");
   
    fprintf(stderr,"\n\nUsage:  %s filename_root(s)  -1<mne[pref]/value> -2<mne[pref]/value[b]> -P<list_of_properties> [-D<dirname>] [-E<file_extent>] [-I] [-L] [-M] [-T] [-Y] [-W<window>[/<w_incr>] [-Z<zmin/zmax>]", program);

   fprintf(stderr,"\n\n  List of filenames MUST be first argument");
   fprintf(stderr,"\n  or input is expected to come from stdin...");
   fprintf(stderr,"\n    -1  : surface 1 property/value (if property = s_");
   fprintf(stderr,"\n          you need to specify pref after s_");
   fprintf(stderr,"\n    -2  : surface 2 property/value (if property = s_");
   fprintf(stderr,"\n          you need to specify pref after s_");
   fprintf(stderr,"\n          Append b to average down to the bottom if this");
   fprintf(stderr,"\n          surface is deeper than the bottom.  This will");
   fprintf(stderr,"\n          only occur if deepest observation is within 100m of the seafloor.");
   fprintf(stderr,"\n    -P  : list of properties: ex: -Pth/sa");
   fprintf(stderr,"\n          by itself, -P will list available properties");
   fprintf(stderr,"\n   [-D] : specifies dirname for input files (default is current");
   fprintf(stderr,"\n          directory)  ex: -D../data/ ");
   fprintf(stderr,"\n   [-E] : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n   [-I] : include station ID in output listing");
   fprintf(stderr,"\n   [-L] : include lat/lon in output listing");
   fprintf(stderr,"\n   [-M] : include month in output listing");
   fprintf(stderr,"\n   [-Y] : include year in output listing");
   fprintf(stderr,"\n   [-T] : include thickness of layer in output listing");
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

      if (((enum property)index == S_ ) || ((enum property) index == PE) || ((enum property) index == HT) ) {
         if (sscanf(st, "%lf", &ref_val) != 1) {
             fprintf(stderr, "\n Specify a ref pressure for  %.2s", prop);
             fprintf(stderr,"\n   ex: -P%.2s1500/th/sa\n", prop);
             exit(1);
         }
         while (!(*st=='/' || *st=='\0' || *st==' '))
            ++st;

         if (s_pref >= 0) {   /* has this already been requested? */
                     if ( NINT( s_pref) != NINT(ref_val)) {
                         fprintf(stderr,"Sorry. You have requested the property s_ with multiple reference levels.\n");
                         fprintf(stderr,"You can only use one of those prefs");
                         fprintf(stderr,"  You will probably have to do multiple runs for each different pref you want to associate with s_\n\n");
                        exit(1);
                     }
          }
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
  
 /* !**!  Special cases for properties */
  
  if (prop_req[(int)GE] && !prop_req[(int)GN]) {
     prop_req[(int)GN] = 1;
     prop_needed[(int)GN] = 1;
     ++nprops;
  }
  
  if (prop_req[(int)GE] || prop_req[(int)GN]) 
     gamma_nc_init((char *) GAMMA_NC_PATH, &ginfo);


 /* end of special cases */

  return (nprops);
}  /* end parse_prop_list() */
/*****************************************************************************/
void get_hydro_data(int file)
   /*  Reads each station in a HydroBase file and computes property values
       to find value at each surface.  This module requires that the HydroBase
       file contains a minimum of pr, de, te, sa observations.  */
{
   int error, i, j, nreq;
   int main_props_avail;
   double  p1, p2, avg, pref;


/* read each station in file ... */

    while ((error = get_station(file, &hdr, &station)) == 0) {

  /* ensure that pr, de, te, and sa are available ... */
       main_props_avail = 1;
       if (!(available(PR, &hdr) && available(DE, &hdr) && available(TE, &hdr) 
           && available(SA, &hdr))) {
         main_props_avail = 0;
       }

/* get pressure associated with 1st surface...*/  
  
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
       
      p1 = hb_linterp(x1, station.observ[prop1], station.observ[(int)PR],
           hdr.nobs);
           
      if (p1 < -8.) {  /* didn't find surf 1 */
           
        switch ((enum property) prop1) {
         case S0:
         case S1:   /* if surface density obs > density we are */
         case S2:   /* seeking, set p1 to first pressure in array */
         case S3:
         case S4:
         case S_:
         case GN:
             if ((station.observ[prop1][0] > x1) 
             && (station.observ[(int)PR][0] < 150.)) {
                p1 = station.observ[(int)PR][0];
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
           
/* get pressure associated with 2nd surface...*/    
 
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
           
      p2 = hb_linterp(x2, station.observ[prop2], station.observ[(int)PR], hdr.nobs);

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
                   p2 = station.observ[(int)PR][hdr.nobs -1];
                }
                break;
             default:
                ;
           }  
         }
      }

      if ((p1 > -1) && (p2 > -1) && (p2 >= p1)) {
        if (yflag)
           fprintf(stdout,"%5d ", hdr.year);
        if (mflag)
           fprintf(stdout,"%2d ", hdr.month);
        if (lflag)
           fprintf(stdout,"%8.3f %8.3f ", hdr.lon, hdr.lat);
        if (iflag)
          fprintf(stdout,"%5d ", hdr.station);
     
        if (tflag)
	  fprintf(stdout, "%10.3lf ",  p2-p1);
       
        for (i = 0; i < MAXPROP; ++i) {
          if (prop_req[i]) {
             if (!available((enum property) i, &hdr)) {
                free_and_alloc(&station.observ[i], hdr.nobs);
                switch ((enum property) i) {
                  case S_:
                      pref = s_pref;
                      break;
                  case PE:
                      pref = pe_pref;
                      break;
                  case HT:
                      pref = ht_pref;
                      break;
                  default:
                      pref = 0;
                }  /* end switch */
                
                compute_prop(i, &station.observ[i], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA], hdr.nobs, pref, (double) hdr.lat);
             }
             avg = get_weighted_av(p1, p2, station.observ[i]);
             fprintf(stdout, " %10.4lf", avg);
          }
        }
        fprintf(stdout, "\n");
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
double get_weighted_av(double p1, double p2, double *xptr)

/*  computes a weighted property average over the pressure interval
    specified by p1 and p2.  Returns the average or -999. for no value.  */
{
   int i, n, start, end;
   double x1, x2;
   double weight, w, sum;
   double *xtmp, *ptmp;
  
   if (xptr == NULL) 
      return (-999.0);
   
   xtmp = (double *) calloc ((size_t)hdr.nobs, sizeof(double));
   ptmp = (double *) calloc ((size_t)hdr.nobs, sizeof(double));
   n = 0;
   for (i = 0; i < hdr.nobs; ++i) {
      if (xptr[i] >= 0.0) {
         xtmp[n] = xptr[i];
         ptmp[n] = station.observ[(int)PR][i];
         ++n;
      }
   }
   if (n == 0){
      free(ptmp);
      free(xtmp);
      return (-999.);
   }   
      
   x1 = hb_linterp(p1, ptmp, xtmp, n);
   x2 = hb_linterp(p2, ptmp, xtmp, n);
   
   if ((x1 < 0.0) || (x2 < 0.0)) {
      free(ptmp);
      free(xtmp);
      return (-999.);
   }

   /* find index of first pressure greater than pressure at top of layer */   
   start = 0;
   while (p1 > ptmp[start])
     ++start;
   /*  find largest index of pressure array less than pressure at bottom of layer */
   end = n-1;
   while ((p2 <= ptmp[end]))
     --end;
   
   /* now do weighted average... */
   
   w = 0.0;
   sum = 0.0;
   for (i = start; i <= end; ++i) {
         weight = ptmp[i] - p1;
         w += weight;
         sum += (xtmp[i] + x1) * 0.5 * weight;
         p1 = ptmp[i];
         x1 = xtmp[i];
   }
   
   /* add last pressure interval ... */
   
   weight = p2 - p1;
   w += weight;
   sum += (x1 + x2) *.5 * weight;
   
   free(ptmp);
   free(xtmp);
   
   return (sum / w);
   
}  /* end get_weighted_av() */
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
	       
        case DR:
               for ( j= 0; j< n; ++j) {
	          xptr[j] = hb_ratio(sptr[j], tptr[j], pptr[j]);
	       }
               break;
        case AL:
               for (j= 0; j< n; ++j) {
	          xptr[j] = hb_alpha(sptr[j], tptr[j], pptr[j]);
	       }
               break;
        case BE:
               for (j= 0; j< n; ++j) {
	          xptr[j] = hb_beta(sptr[j], tptr[j], pptr[j]);
	       }
               break;
       case VS:
               compute_sound_vel( xptr, pptr, tptr, sptr, n);
               break;

       case PV: 
               buoy_freq(xptr, pptr, tptr, sptr, n, window, w_incr);
               for (j = 0; j < n; ++j) {
                 if (xptr[j] > -99.) 
                   xptr[j] *= xptr[j];
               }
               po_vort(xptr, xptr, n, dlat);

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
