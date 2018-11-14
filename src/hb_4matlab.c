/*  hb_4matlab.c

................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Institution
                             1995
                             updated to ANSI-C Feb 2000
................................................................................
____________________________________________________________________________
hb_4matlab computes properties at each pressure level in a station and 
outputs an ascii listing of each observed level containing all the properties
specified with the -P option plus user-specified info: { lat/lon year month station/cruise}
____________________________________________________________________________

  USAGE:  

 hb_4matlab filename(_roots) -P<property_list> [-L] [-Y] [-M] [-S] [-W<window>[/<w_incr>] [-D<dirname>] [-E<file_extent>]

  list of filenames must be first argument or input is expected from stdin
  
 -P : list of properties to evaluate; ex: -Ppr/te/th/sa/ox
      -P (by itself) will print a list of the available properties;
      
 -L : output lat/lon of observed level

 -Y : output year of observed level

 -M : output month of observed level

 -S : output station#/cruise# of observed level

 -D : specifies directory for input files (default is current dir)  
      ex: -D/home/pete/data

 -E : specifies input_file extent (default is no extent) ex: -E.sum

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

    /* flags to indicate which properties are to be output and which are
       needed for computation */
      
short prop_req[MAXPROP];         /* set of props requested for output */
short prop_needed[MAXPROP];      /* set of props requested
                                     or needed for computation */

   /* global variables to store station ... */

struct HYDRO_DATA station;
struct HYDRO_HDR hdr;
double s_pref;           /* ref lev for computing a non-standard sigma level */
double pe_pref;           /* ref lev for computing potential energy */
double ht_pref;           /* ref lev for computing dynamic height */
int window, w_incr;     /* used to specify pr window for gradient properties*/

int yopt, lopt, sopt, mopt, dopt;   /* option flags */

struct GAMMA_NC ginfo;


  /* prototypes for locally defined functions */
  
void print_usage(char *);
int parse_prop_list(char *);
void get_hydro_data(int);

main (int argc, char **argv)
{
   short   popt;
   int     index, nprops = 0;
   int     i;
   int     curfile = 1, nfiles = 0; 
   int     print_msg = 1;       /* print message in file open routines*/
   int     error, prompt;
   char    *extent, *dir, *st;
   int     infile;


/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
 
/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    popt = yopt = lopt = sopt = mopt = dopt = 0;
    error = 0;
    s_pref = -1;
    window = 100;
    w_incr = 10;

/* initialize these ... */

   for (i = 0; i < MAXPROP; ++i) {
      prop_req[i] = 0;
      prop_needed[i] = 0;
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

               case 'L':
                        lopt = 1;
                        break;
                        
               case 'Y':
                        yopt = 1;
                        break;
                        
               case 'M':
                        mopt = 1;
			if (argv[i][2] == 'd')
			    dopt = 1;
                        break;
                        
               case 'S':
                        sopt = 1;
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

   if ( !popt) {
       fprintf(stderr,"\nYou must specify a list of properties to output.\n");
       exit(1);
   }
   if (!nfiles) {
       fprintf(stderr,"\nExpecting input from stdin ... ");
       infile = STDIN;
   }

   if (! ( lopt || mopt || sopt || yopt)) {
       fprintf(stderr,"\nWARNING! no header info (year, month, position, cruise or station id) \nhas been specified for output.\n");
   }

/* loop for each input file */

   do {

     if (nfiles > 0) {
       infile = open_hydro_file(dir, argv[curfile], extent, print_msg);
       if (infile < 0)
          goto NEXTFILE;

     }
            /* read each file completely */

 
     get_hydro_data(infile);

NEXTFILE:

     if (nfiles > 0) 
         close(infile);

   } while (curfile++ < nfiles );

   fprintf(stderr,"\n\nEnd of hb_4matlab.\n");
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{  
 fprintf(stderr,"\n%s computes properties specified with -P option at each depth level in a station \n and outputs an ascii listing of each level plus the lat/lon/year/month/station/cruise of the observation.\n", program);

   fprintf(stderr,"\nUsage:  %s filename_root(s)  -P<list_of_properties>  [-D<dirname>] [-E<file_extent>] [-W<window>/<w_incr>]", program);

   fprintf(stderr,"\n\n  List of filenames MUST be first argument!");
   fprintf(stderr,"\n    -P  : list of properties to project onto surface;");
   fprintf(stderr,"\n          ex:  -Ppr/th/sa/ox/ht");
   fprintf(stderr,"\n               -P (by itself) produces a list of available properties");
   fprintf(stderr,"\n   [-L] : option to output lat/lon ");
   fprintf(stderr,"\n   [-Y] : option to output year ");
   fprintf(stderr,"\n   [-M] : option to output month.  Append d to optionally output day: -Md ");
   fprintf(stderr,"\n   [-S] : option to output station#/cruise# ");
   fprintf(stderr,"\n   [-D] : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n   [-E] : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n   [-W] : Specifies pressure window (db) for computing ");
   fprintf(stderr,"\n          gradient properties (bvfreq, pv...) This ");
   fprintf(stderr,"\n          constitutes the range over which observations");
   fprintf(stderr,"\n          are incorporated into the gradient computation.");
   fprintf(stderr,"\n          The window is centered around each pressure ");
   fprintf(stderr,"\n          level being evaluated.");
   fprintf(stderr,"\n          w_incr specifies how finely to subdivide the");
   fprintf(stderr,"\n          window into increments(db) to create");
   fprintf(stderr,"\n          an evenly spaced pressure series over the window.");
   fprintf(stderr,"\n          defaults: -W100/10");
   fprintf(stderr,"\n    [-h] help...... prints this message. \n");

   fprintf(stderr,"\n\n");  
   return;
}

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

  /* some special cases ... */

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
  
 
 /* !**!  Special cases for properties */
  
  if (prop_req[(int)GE] && !prop_req[(int)GN]) {
     prop_req[(int)GN] = 1;
     prop_needed[(int)GN] = 1;
  }
  
  if (prop_req[(int)GE] || prop_req[(int)GN]) 
     gamma_nc_init((char *) GAMMA_NC_PATH, &ginfo);
     
 /* end of special cases */
  
  return (nprops);
}  /* end parse_prop_list() */
 

/*****************************************************************************/
void get_hydro_data(int file)

   /*  Reads each station in a HydroBase file and computes property values
       at each standard level.    */
{
   int error, i, j, nreq, nbytes, pts_avail, ratio_done;
   double dlon, dlat, deltap;
   double *scan;
   char str[40], *st;


/* read each station in file ... */

    while ((error = get_station(file, &hdr, &station)) == 0) {

  /* check if basic properties pr, te, and sa are available ... */

    pts_avail = (available(PR, &hdr) && available(TE, &hdr) 
                 && available(SA, &hdr));

    ratio_done = 0;
    
 /* compute appropriate properties at each level in station ... */

/* !**! Special cases for individual properties... */

    for (i = 0; i < MAXPROP; ++i) {
       if (prop_needed[i] && !available((enum property)i, &hdr)) {
          switch ((enum property) i) {
             case OX:
                 if (available(O2, &hdr) && pts_avail) {
                    free_and_alloc(&station.observ[i], hdr.nobs);
                    for (j=0; j < hdr.nobs; ++j) {
                       station.observ[i][j] = ox_kg2l(station.observ[(int)O2][j], station.observ[(int)PR][j], station.observ[(int)TE][j],station.observ[(int)SA][j]);
                    }
                 }
               break;
             case O2:
                 if (available(OX, &hdr) && pts_avail) {
                    free_and_alloc(&station.observ[i], hdr.nobs);
                    for (j=0; j < hdr.nobs; ++j) {
                       station.observ[i][j] = ox_l2kg(station.observ[(int)OX][j], station.observ[(int)PR][j], station.observ[(int)TE][j],station.observ[(int)SA][j]);
                    }
                 }
               break;
             case TH:
	       if (pts_avail) {
                  free_and_alloc(&station.observ[i], hdr.nobs);
                  compute_theta(hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA] );
	       }
               break;
             case TP:
	       if (pts_avail) {
		 free_and_alloc(&station.observ[i], hdr.nobs);
		 compute_dtdp(station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], hdr.nobs, window, w_incr);
	       }
	       break;
             case TZ:
	       if (pts_avail) {
		 free_and_alloc(&station.observ[i], hdr.nobs);
		 dlat = (double) hdr.lat;
		 compute_dtdz(station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA],hdr.nobs,window,w_incr,dlat);
	       }
               break;
             case HC:
	       if (pts_avail) {
		 free_and_alloc(&station.observ[i], hdr.nobs);
		 compute_heat_capacity(hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA] );
	       }
               break;
             case S0:
	       if (pts_avail) {
                  free_and_alloc(&station.observ[i], hdr.nobs);
                  compute_sigma(0., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               }
	       break;
             case S1:
	       if (pts_avail) {
                  free_and_alloc(&station.observ[i], hdr.nobs);
                  compute_sigma(1000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               }
	       break;
             case S2:
	       if (pts_avail) {
                  free_and_alloc(&station.observ[i], hdr.nobs);
                  compute_sigma(2000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               }
	       break;
             case S3:
	       if (pts_avail) {
                  free_and_alloc(&station.observ[i], hdr.nobs);
                  compute_sigma(3000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               }
	       break;
             case S4:
	       if (pts_avail) {
                  free_and_alloc(&station.observ[i], hdr.nobs);
                  compute_sigma(4000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               }
	       break;
             case S_:
	       if (pts_avail) {
                  free_and_alloc(&station.observ[i], hdr.nobs);
                  compute_sigma(s_pref, hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               }
	       break;
             case HT:
	       if (pts_avail) {
                  free_and_alloc(&station.observ[i], hdr.nobs);
                  compute_height(hdr.nobs, station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], ht_pref, station.observ[i]);
               }
	       break;
             case PE:
	       if (pts_avail) {
                  free_and_alloc(&station.observ[i], hdr.nobs);
                  compute_energy(hdr.nobs, station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], pe_pref, station.observ[i]);
               }
	       break;
             case SV:
	       if (pts_avail) {
                  free_and_alloc(&station.observ[i], hdr.nobs);
                  compute_sp_vol( hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               }
	       break;
             case DR:
	     case AL:
	     case BE:
	       if (pts_avail && !ratio_done) {
                  free_and_alloc(&station.observ[DR], hdr.nobs);
                  free_and_alloc(&station.observ[AL], hdr.nobs);
                  free_and_alloc(&station.observ[BE], hdr.nobs);
                  compute_ratio( hdr.nobs, station.observ[DR], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], station.observ[(int)AL], station.observ[(int)BE]);
		  ratio_done = 1;
               }
	       break;

             case VA:
	       if (pts_avail) {
                  free_and_alloc(&station.observ[i], hdr.nobs);
                  compute_svan( hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               }
	       break;
	       
             case VS:
	       if (pts_avail) {
                  free_and_alloc(&station.observ[i], hdr.nobs);
                  compute_sound_vel( station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], hdr.nobs);
               }
               break;

             case PV: 
	       if (pts_avail) {
                  free_and_alloc(&station.observ[i], hdr.nobs);
                  dlat = (double) hdr.lat;
                  buoy_freq(station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], hdr.nobs, window, w_incr);
                  for (j = 0; j < hdr.nobs; ++j) {
                    if (station.observ[i][j] > -99.)  
                       station.observ[i][j] *= station.observ[i][j];
                  }
                  po_vort(station.observ[i],station.observ[i], hdr.nobs, dlat);
               }
               break;
             case BF:
	       if (pts_avail) {
                  free_and_alloc(&station.observ[i], hdr.nobs);
                  buoy_freq(station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], hdr.nobs, window, w_incr);
                  for (j = 0; j < hdr.nobs; ++j) {
                    if (station.observ[i][j] > -99.)   /* convert the units */
                        station.observ[i][j] *= 1.0e5;
                  }
               }
	       break;

             case GN:
	       if (pts_avail && !prop_req[(int)GE]) {
                  free_and_alloc(&station.observ[i], hdr.nobs);
                  dlat = (double) hdr.lat;
                  dlon = (double) hdr.lon;
		  compute_gamma_n(&ginfo, hdr.nobs, station.observ[i], station.observ[(int)PR],station.observ[(int)TE],station.observ[(int)SA], dlon, dlat);
               }
	       break;
	       
             case GE:
	       if (pts_avail) {
                  free_and_alloc(&station.observ[(int)GE], hdr.nobs);
                  free_and_alloc(&station.observ[(int)GN], hdr.nobs);
                  dlat = (double) hdr.lat;
                  dlon = (double) hdr.lon;
		  compute_gamma_nerr(&ginfo, hdr.nobs, station.observ[(int)GN],   
		     station.observ[(int)GE], station.observ[(int)PR],
		     station.observ[(int)TE], station.observ[(int)SA], dlon,
		     dlat);
               }
	       break;
	       
             default:
               break;
          } /* end switch */
       } /* end if */
     } /* end for */
     
 /* determine availability of requested properties... */
    nreq = 0;
    for (i = 0; i < MAXPROP; ++i) {
       if (prop_req[i] && station.observ[i] != NULL )
          ++nreq;
    }
    
 /* output the station */
    hdr.nprops = station.nprops = nreq;
    free((void *)hdr.prop_id);
    hdr.prop_id = (int *) malloc(nreq * sizeof(int));
    j = 0;
    for (i = 0; i < MAXPROP; ++i) {
       if (! prop_req[i] && station.observ[i] != NULL) {
           free((void *)station.observ[i]);
           station.observ[i] = NULL;
       }
       else if (station.observ[i] != NULL)  {
         hdr.prop_id[j++] = i;
       }
    }
    
    nbytes = 0;
    st = &str[0];
    if (lopt) { 
       sprintf(st,"%9.3f %9.3f ", hdr.lon, hdr.lat);
       nbytes += 20;
    }
    
    st = &str[nbytes];
    if (yopt) {
       sprintf(st,"%4d ", hdr.year);
       nbytes += 5;
    }
   
    st = &str[nbytes];
    if (mopt) { 
       sprintf(st,"%2d ", hdr.month);
       nbytes += 3;
    }
    st = &str[nbytes];
    if (dopt) { 
       sprintf(st,"%2d ", hdr.day);
       nbytes += 3;
    }
   
    st = &str[nbytes];
    if (sopt) {
       sprintf(st,"%5d %4d ", hdr.cruise, hdr.station);
       nbytes += 11;
    }
   
    scan = (double *) malloc(hdr.nprops * sizeof(double));
    for (i = 0; i < hdr.nobs; ++i) {
       for (j = 0; j < hdr.nprops; ++j) {
          scan[j] = station.observ[hdr.prop_id[j]][i];
       }
       write(STDOUT, str, nbytes);
       write_hydro_scan(STDOUT, scan, hdr.nprops, hdr.prop_id);
    }
    free(scan);

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

   
