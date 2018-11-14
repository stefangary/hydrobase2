/*  hb_strfun.c

_................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Institution
                             1993
                             Updated to ANSI standards Feb 2000
................................................................................
____________________________________________________________________________
Computes Montgomery Stream Function for a specified density surface.  
A reference surface is optional (default is sea surface).
     msf =  ht1 - (p1 - pbar1) * sva1 - [ht2 - (p2 - pbar2) *sva2]
___________________________________________________________________________

 */


#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include "hydrobase.h"
#include "hydro_cdf.h"
#include "hb_gamma.h"
#include "hb_paths.h"


/* input file pathnames */

#define    EXTENT   ""
#define    DIR      ""


/* set up data structures to store  data */

struct point {
          double  pr1, sv1, ht1;    
          double  pr2, sv2, ht2;    
	  double  lat, lon;
    struct point *next;    
};

double pbar1, pbar2;
double np1, np2;
int reflev_is_pr;
double val1, val2;
double pref1, pref2;
int index1, index2; 

   /* global variables to store station ... */

struct HYDRO_DATA station;
struct HYDRO_HDR hdr;

float HB_f_mask;        /* flag for masked nodes in .cdf file */

struct GAMMA_NC ginfo;   /* used for neutral density */
int warnflag;            

     /* box boundaries  */

float   xmin, xmax, ymin, ymax;  
int 	xdateline;   

   /* prototypes for functions defined locally */
   
void print_usage(char *);
int parse_prop_list(char *, int *, double *, double *, int *);
void get_hydro_data(int, struct point *);
void get_cdf_data(int, struct point *);
void insertdata(struct point *);

main (int argc, char **argv)
{
   short   cdf_flag;
   int     i, blah;
   int     curfile = 1, nfiles = 0; 
   int     print_msg = 1;       /* print message in file open routines*/
   int     error;
   char    *extent, *dir, *st;
   int     infile;
   struct point *listptr;


/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
 
/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    cdf_flag = 0;
    warnflag = 0;          /* set to 1 after warning is printed */
    error = 0;
    xdateline = 0;
    HB_f_mask = HBMASK;
    reflev_is_pr = 1;
    val2 = 0.0;
    pref2 = 0.0;
    index2 = (int)PR;
    sopt = 0;
    ymin = -90;
    ymax = 90;
    xmin = -360;
    xmax = 360;
    
    

/* initialize these ... */

   for (i = 0; i < MAXPROP; ++i) {
      station.observ[i] = (double *) NULL;
   }
   hdr.prop_id = (int *) NULL;



/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'B':                    /* optional box bounds */
                        bopt = 1;
                        st = &argv[i][2];
                           if (*st == '/')
                               ++st;
                        error = (sscanf(st,"%f", &xmin) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%f", &xmax) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%f", &ymin) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%f", &ymax) != 1);
                        
                        if (xmin > 0 && xmax < 0)
                            xmax += 360;
                        if (xmax > 180)
                           xdateline = 1;
                        break;

               case 'C':
                        cdf_flag = 1;
                        break;
               case 'D':                   /* get input dir */
                        dir = &argv[i][2];
                        break;
               case 'E':                    /* get file extent */
                        extent = &argv[i][2];
                        break;
               case 'R':                    /* define reference surface */
                        error = parse_prop_list(&argv[i][2], &index2, &pref2, &val2, &reflev_is_pr);
                        break;

               case 'S':                    /* define surface */
                        sopt = 1;
                        error = parse_prop_list(&argv[i][2], &index1, &pref1, &val1, &blah);	
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

   if (!sopt) {
       fprintf(stderr,"\nYou must specify a surface with -S<prop_mne>[ref_lev]/<value> \n");
       exit(1);
   }

   if (!nfiles) {
       if (cdf_flag) {
          fprintf(stderr,"\nCDF files cannot be piped in via stdin");
	  fprintf(stderr,"\nYou must specify a list. ");
           exit(1);
       }
       fprintf(stderr,"\nExpecting input from stdin ... ");
       infile = STDIN;
   }
   fprintf(stderr,"\n\n  bounds: %.3f %.3f %.3f %.3f\n", ymin,ymax,xmin,xmax );


 

/* loop for each input file */

   do {

               /* netCDF file format */

     if (cdf_flag) {                         

       infile = cdf_open(dir, argv[curfile], extent, print_msg);
       if (infile < 0 )
          goto NEXTFILE;
     }
              /* HydroBase file format */
     else {                    
       if (nfiles > 0) {
         infile = open_hydro_file(dir, argv[curfile], extent, print_msg);
         if (infile < 0)
            goto NEXTFILE;
       }
     }

     
            /* read each file completely */

     if (cdf_flag) 
         get_cdf_data(infile, list) ;
     else 
         get_hydro_data(infile, list);

NEXTFILE:

     if (cdf_flag) 
         cdf_close(infile);
     else {
       if (nfiles > 0)
         close(infile);
     }

   } while (curfile++ < nfiles );


/*  compute stream function and write out results ... */


   fprintf(stderr,"\n    writing to outfile ...\n");

   do_stats(list);

   fprintf(stderr,"\n\nEnd of %s.\n", argv[0]);
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nUsage:  %s filename_root(s)  [-C] -B/west/east/south/north -I<delta_x[/delta_y]> -P<list_of_properties> [-S<surface_def_file>] [-W<window>[/<w_incr>]] [-D<dirname>] [-E<file_extent>] [-Z]", program);

   fprintf(stderr,"\n\n  List of filenames MUST be first argument or ");
   fprintf(stderr,"\nstation data is expected from stdin.");
   fprintf(stderr,"\n -S  : defines surface for stream function:  prop/value.");
   fprintf(stderr,"\n       <prop> can be one of the following density types:");
   fprintf(stderr,"\n       s0  s1  s2  s3  s4  s_<pref>   gn ");
   fprintf(stderr,"\n\n    OPTIONS:");
   fprintf(stderr,"\n[-B]  : specifies grid bounds");
   fprintf(stderr,"\n[-C] : input files are cdf format files.");
   fprintf(stderr,"\n       (default is HydroBase format)");
   fprintf(stderr,"\n[-D] : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n       ex: -D../data/ ");
   fprintf(stderr,"\n[-E] : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n[-R]  : defines reference level:  prop/value.");
   fprintf(stderr,"\n       <prop> can be a density or pressure surface:");
   fprintf(stderr,"\n       (default is pr/0.0)");
   fprintf(stderr,"\n[-h]  :  help -- prints this message.");
   fprintf(stderr,"\n\n");  
   return;
}

/*****************************************************************************/
int parse_prop_list(char *st, int *index_ptr, double *pref_ptr, double *value_ptr, int *is_pressure_ptr)
    /*  Parses a string for property mnemonic and value
        and returns appropriate information.
          Returns 0 if successful, or 1 for an error. */
{
   int index;
   char prop[4];


      if (*st == '/')
         ++st;
      sscanf(st,"%2s", prop);
      *index_ptr = get_prop_indx(prop);
      if (*index_ptr < 0) {
         fprintf(stderr,"\n Unknown property requested: %s\n", prop);
         return(1);
      }
      if (*index_ptr == (int)PR)
        *is_pressure_ptr = 1;
	
      ++st;
      ++st;
      if (*st == '/')
         ++st;


      if ((enum property)index == GN )  {
         gamma_nc_init((char *) GAMMA_NC_PATH, &ginfo);
      }
      
      if ((enum property)index == S_ )  {
         if (sscanf(st, "%lf", pref_ptr) != 1) {
             fprintf(stderr, "\n Specify a ref pressure for  %.2s", prop);
             fprintf(stderr,"\n   ex: -S%.2s1500/34.65\n", prop);
             exit(1);
         }
         while (!(*st=='/' || *st=='\0' || *st==' '))
            ++st;

      }

      if (*st == '/')
         ++st;
	 
      if (sscanf(st, "%lf", value_ptr) != 1) {
             fprintf(stderr,"\nError parsing property value." );
             return(1);
      }

      return(0);
}  /* end parse_prop_list() */
/****************************************************************************/
void get_hydro_data(int file, struct surface *listptr)

   /*  Reads each profile in a HydroBase file and adds the property information
       to a linked list.   This module requires that the HydroBase
       file contains a minimum of pr,  te, sa observations.  */
{
   int error, i;
   int main_props_avail;


/* read each station in file ... */

    while ((error = get_station(file, &hdr, &station)) == 0) {

  /* check for and skip over out of bounds stations   */
  
       if (xdateline && hdr.lon < 0)
          hdr.lon += 360;
          

       if ((hdr.lon > xmax) || (hdr.lon < xmin) || (hdr.lat > ymax) || (hdr.lat < ymin))  
       continue;
         
         
/* ensure that pr, te, and sa are available ... */


       if (!(available(PR, &hdr) &&  available(TE, &hdr) && available(SA, &hdr))) 
         continue;


 /* compute appropriate properties at each level in station */

/* !**! Special cases for individual properties... */

       
        free_and_alloc(&station.observ[(int)HT], hdr.nobs);
        compute_height(hdr.nobs, station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], 0.0, station.observ[HT]);
	   
        free_and_alloc(&station.observ[(int)VA], hdr.nobs);
        compute_svan( hdr.nobs, station.observ[VA], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);

        i = index1;
	if (!available((enum property)i, &hdr)) {
           switch (i) {
             case TH:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_theta(hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA] );
	       ++hdr.nprops;
	       hdr.prop_id = (int *) realloc (hdr.prop_id, hdr.nprops);
	       hdr.prop_id[hdr.nprops-1] = i;
               break;
             case S0:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(0., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
	       hdr.prop_id = (int *) realloc (hdr.prop_id, hdr.nprops);
	       hdr.prop_id[hdr.nprops-1] = i;
               break;
             case S1:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(1000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
	       hdr.prop_id = (int *) realloc (hdr.prop_id, hdr.nprops);
	       hdr.prop_id[hdr.nprops-1] = i;
               break;
             case S2:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(2000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
	       hdr.prop_id = (int *) realloc (hdr.prop_id, hdr.nprops);
	       hdr.prop_id[hdr.nprops-1] = i;
               break;
             case S3:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(3000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
	       hdr.prop_id = (int *) realloc (hdr.prop_id, hdr.nprops);
	       hdr.prop_id[hdr.nprops-1] = i;
               break;
             case S4:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(4000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
	       hdr.prop_id = (int *) realloc (hdr.prop_id, hdr.nprops);
	       hdr.prop_id[hdr.nprops-1] = i;
               break;
             case S_:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(s_pref, hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
	       hdr.prop_id = (int *) realloc (hdr.prop_id, hdr.nprops);
	       hdr.prop_id[hdr.nprops-1] = i;
               break;
             case GN:
               free_and_alloc(&station.observ[i], hdr.nobs);
	       compute_gamma_n(&ginfo, hdr.nobs, station.observ[i],
	           station.observ[(int)PR],station.observ[(int)TE],
		   station.observ[(int)SA], (double)hdr.lon, (double)hdr.lat);
	       hdr.prop_id = (int *) realloc (hdr.prop_id, hdr.nprops);
	       hdr.prop_id[hdr.nprops-1] = i;
	       break;
	       
               default:
	         fprintf(stderr,"\nProperty [%s] is not an appropriate surface for stream function.", get_prop_mne(i));
                 exit(1);
            } /* end switch */
	  } /* end if !available*/
	  
	  i =  index2;
	  if (!available((enum property)i, &hdr) && (index2 != index1)) {
            switch (i) {
              case TH:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_theta(hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA] );
               break;
             case S0:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(0., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               break;
             case S1:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(1000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               break;
             case S2:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(2000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               break;
             case S3:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(3000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               break;
             case S4:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(4000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               break;
             case S_:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(s_pref, hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               break;
             case GN:
               free_and_alloc(&station.observ[i], hdr.nobs);
	       compute_gamma_n(&ginfo, hdr.nobs, station.observ[i],
	          station.observ[(int)PR],station.observ[(int)TE],
		  station.observ[(int)SA], (double)hdr.lon, (double)hdr.lat);
	       break;
	       
             default:
	         fprintf(stderr,"\nProperty [%s] is not an appropriate surface for stream function.", get_prop_mne(i));
                 exit(1);
             } /* end switch */
	  } /* end if !available*/
	  
      /* add a point to the linked list & interpolate data onto each 
         surface */

	    
         listptr = insertdata(listptr);


      /* free up station storage */
       
      for (i = 0; i < MAXPROP; ++i) {
          if (station.observ[i] != NULL) {
             free((void *) station.observ[i]);
             station.observ[i] = NULL; 
          }    
      }

   }  /* end while */

   if (error > 0)
      report_status(error, stderr);

   return;

} /* end get_hydro_data() */

/****************************************************************************/

struct point *insertdata(struct point *start_ptr)

   /*   Projects property values onto surface, creates a new record
        and adds values to the appropriate structure field.   
	Inserts record at start of list and 
	returns a  pointer to new start of list. 
        
        arguments:
            start_ptr:   pointer to start of linked list
   */
{
   int i, datagap;
   
   
   

}  /* end insertdata() */

/****************************************************************************/
double project_prop(double *obs, int nobs, double p0)

   /*   Projects property values onto a specified surface...Returns the value
        of the property at the pressure level corresponding to this surface,
        or -9999 if the property cannot be projected.  The global pressure arrays 
        are used.  
        
            obs:   array of observations for property
           nobs:   number of observations
             p0:   interpolated pressure of surface
    */
{
   double val;
   double x[2], y[2];
   int    j;
   int   prevdepth, curdepth;


   prevdepth = curdepth = 0;
   j = 0;

   /* find first valid observation level. */

  while (obs[j] < -8.9) {
     if (++j == nobs)
          return (-9999.);
  }
  x[0] = station.observ[(int)PR][j];
  y[0] = obs[j];

 /* Now, check successive pairs of datapoints until the 
    surface is found or the last observation is encountered.   */

  while (++j < nobs) {
    if (obs[j]  > -8.9) {
        x[1] = station.observ[(int)PR][j]; 
        y[1] = obs[j];
        if ((val = hb_linterp(p0, x, y, 2)) > -9998.) {
              return (val);
        }
        x[0] = x[1];
        y[0] = y[1];
     }
  }  /* end while */

  return (-9999.);
} /* end project_prop() */

/****************************************************************************/
               
void get_cdf_data(int cdfid, struct surface *listptr)
{
   struct CDF_HDR cdf;
   struct surface *surf;
   int startrow, startcol;
   int endrow, endcol;
   float minlat, minlon;
   float maxlat, maxlon;
   float xoffset, yoffset;
   float lat, lon;
   double dlat, dlon;
   float *z, *x;
   short **count;
   char *mne;
   int error, i, j, k, sq, row, col, npts;
   int index_prop;
   short prop_avail[MAXPROP];     /* set of props in .cdf file */
   short compute_flag[MAXPROP];   /* props NOT in .cdf which can be computed*/


   if (error = read_cdf_hdr(cdfid, &cdf)) 
         exit (1);

/* compare bounds of file to these grid bounds to ensure overlap */

   if (xdateline && cdf.xmin < 0)
      cdf.xmin += 360;
   if (xdateline && cdf.xmax < 0)
      cdf.xmax += 360;
      
   if ((cdf.xmin > xmax) || (cdf.xmax < xmin) 
      || (cdf.ymin > ymax) || (cdf.ymax < ymin)) {
       fprintf(stderr, "\nNo data in cdf file is within grid bounds.\n");
       return;
   }

   xoffset = yoffset = 0.0;
   
   if (cdf.node_offset) {   /* for pixel gridded cdf file */
      xoffset = 0.5 * cdf.xincr;
      yoffset = 0.5 * cdf.yincr;
   }
    
   minlon = ((cdf.xmin + xoffset) < xmin) ? xmin : cdf.xmin + xoffset;
   maxlon = ((cdf.xmax - xoffset) > xmax) ? xmax : cdf.xmax - xoffset;
   minlat = ((cdf.ymin + yoffset) < ymin) ? ymin : cdf.ymin + yoffset;
   maxlat = ((cdf.ymax - yoffset) > ymax) ? ymax : cdf.ymax - yoffset;
   
   
   error = get_indices(&cdf, maxlat, minlon, &startrow, &startcol);
   error = get_indices(&cdf, minlat, maxlon, &endrow, &endcol);
   
   if (startrow < 0 ) startrow = 0;
   if (endrow >= cdf.ny) endrow = cdf.ny-1;
   if (startcol < 0 ) startcol = 0;
   if (endcol >= cdf.nx) endcol = cdf.nx-1;
   
/* determine which properties are available in the cdf file */

   for (i = 0; i < MAXPROP; ++i) {  /* initialize these flags */
      prop_avail[i] = 0;
      compute_flag[i] = 0;
   }

   prop_avail[(int)DE] = 1;        /* depth is an index variable */
   prop_needed[(int)DE] = 1;
   prop_needed[(int)PR] = 1;      /* must always be in a cdf file */

   for (i = 0; i < cdf.nprops; ++i)    
      prop_avail[get_prop_indx(cdf.prop_id[i])] = 1;
   
/* determine which properties can and should be computed ...   */

/*  !**! Special cases for individual properties */

   for (i = 0; i < MAXPROP; ++i) {
      if (prop_needed[i] && !prop_avail[i]) {
         switch ((enum property) i) {
            case OX: 
                compute_flag[i] = prop_avail[(int)O2] && prop_avail[(int)PR] &&
                                  prop_avail[(int)TE] && prop_avail[(int)SA];
                if (compute_flag[i] == 0) {
                   fprintf(stderr,"No oxygens available in this file. \n"); 
                   prop_needed[i] = 0;
                   prop_req[i] = 0;
                }
                else {
                   prop_needed[(int)PR] = 1;
                   prop_needed[(int)TE] = 1;
                   prop_needed[(int)SA] = 1;
                   prop_needed[(int)O2] = 1;
                   fprintf(stderr,"\n OX not available in cdf file, but will be computed from o2, pr, te, sa values");
                }
                break; 
            case O2:
                compute_flag[i] = prop_avail[(int)OX] && prop_avail[(int)PR] &&
                                  prop_avail[(int)TE] && prop_avail[(int)SA];
                if (compute_flag[i] == 0) {
                   fprintf(stderr,"No oxygens available in this file. \n"); 
                   prop_needed[i] = 0;
                   prop_req[i] = 0;
               }
                else {
                   prop_needed[(int)PR] = 1;
                   prop_needed[(int)TE] = 1;
                   prop_needed[(int)SA] = 1;
                   prop_needed[(int)OX] = 1;
                   fprintf(stderr,"\n O2 not available in cdf file, but will be computed from ox, pr, te, sa values");
               }
               break; 
	       
            case S0: 
            case S1: 
            case S2:    /* fall through */
            case S3: 
            case S4: 
            case S_:
            case TH:
            case GN:
            case GE:
            case HT:
            case PE:
            case BF:
            case PV:
            case SV:
            case VA:
                compute_flag[i] = prop_avail[(int)PR] && prop_avail[(int) TE] 
                         && prop_avail[(int) SA];

                if (compute_flag[i] == 0) {
                  fprintf(stderr,"\n\n FATAL ERROR!");
                  fprintf(stderr,"\n Unable to compute %.2s for this file.\n", 
                          get_prop_mne(i));
                  fprintf(stderr,"The .cdf file must include pr, te, sa \n"); 
                  exit (0);
                }

                prop_needed[(int)PR] = 1;
                prop_needed[(int)TE] = 1;
                prop_needed[(int)SA] = 1;

                fprintf(stderr,"\n %.2s not available in cdf file, but will be computed from averaged p,t,s values.", get_prop_mne(i));
                break;

            default:
                fprintf(stderr,"\n\n FATAL ERROR!");
                fprintf(stderr,"\n  Property <%.2s> not available in this file.\n", get_prop_mne(i));
                exit(0);

         }  /* end switch */
      }
   }   /* end for */

/* get depth values */

   z = (float *) malloc ((size_t) cdf.nz * sizeof(float));

   if (read_cdf_depths(cdfid, z) < 0) {
      exit(1);
   }

/* allocate space for property and count arrays; initialize count arrays to 0 */

   for (i = 0; i < MAXPROP; ++i) {
      if (prop_needed[i]) {
         free_and_alloc(&station.observ[i], cdf.nz);
      }
   }

   x = (float *) malloc ((size_t) cdf.nz * sizeof(float));
   
   if (cdf.counts_included) {
      count = (short **) malloc ((size_t) MAXPROP * sizeof(short *));
      if (count == NULL) {
          fprintf(stderr, "\nInsufficient memory for call to malloc()");
          exit(1);
      }
      for (i = 0; i < MAXPROP; ++i) {
         count[i] = NULL;
         if (prop_needed[i]) {
            count[i] = (short *) malloc((size_t) cdf.nz * sizeof(short));
            if (count[i] == NULL) {
              fprintf(stderr, "\nInsufficient memory for call to malloc()");
              exit(1);
            }
            for (j = 0; j < cdf.nz; ++j)
               count[i][j] = 0;
         }
      }
   }
   
/* visit each appropriate gridpt in the cdf file  */

   for (row = startrow; row <= endrow; ++row) {
       for (col = startcol; col <= endcol; ++col) {

          
          /* convert lat/lon of this data from cdf file to the grid position
             for the surfaces. */

           error = get_lat_lon(&cdf, row, col, &lat, &lon);
           if ((lat > ymax) || (lat < ymin) || (lon > xmax) || (lon < xmin))
               continue;

           dlat = (double) lat;
	   dlon = (double) lon; 
	   
	   /* this assumes we are using pixel grid registration */
	     
           sq =  NINT((lat + .0001 - ymin) / delta_y - 0.5) * ncols + 
                 NINT((lon + .0001 - xmin) / delta_x - 0.5);

           /* get bottom depth for this square ...*/

           read_cdf_bottom_depth(cdfid, &z[cdf.nz-1], row, col, tbin);

           for (i = 0; i < MAXPROP; ++i) {

              /* extract available properties ... */

              if (prop_needed[i] && prop_avail[i]) {

                 mne = get_prop_mne(i);
                 error = read_cdf_prop(cdfid, mne, x, row, col, tbin, 0, cdf.nz);
                 if (error > 0) {
                    fprintf(stderr,"\nError attempting to read %.2s at row,col =  %d,%d from cdf file.", mne, row, col);
                    exit(1);
                 }

                 if (cdf.counts_included) {
                    error = read_cdf_prop_count(cdfid, mne, count[i],
                             row, col, tbin, 0, cdf.nz);
                    if (error > 0) {
                       fprintf(stderr,"\nError attempting to read %.2s_cnt at row,col =  %d,%d from cdf file.", mne, row, col);
                       exit(1);
                    }

                 }

            /* place extracted data into appropriate column of station.observ */

                 for (j = 0; j < cdf.nz; ++j)
                    station.observ[i][j] = (double) x[j];
              }

           } /* end for  */

           /* choose the best property to determine whether any data is
              available at a given depth in .cdf file ... */

           if (station.observ[(int)TE] != NULL)
               index_prop = (int)TE;
           else if (station.observ[(int)SA] != NULL)
               index_prop = (int)SA;
           else if (station.observ[(int)PR] != NULL)
               index_prop = (int)PR;
           else 
               index_prop = get_prop_indx(cdf.prop_id[0]);

           /* eliminate levels with no data */

           k = 0;
           for (i = 0; i < cdf.nz; ++i) {
              if (! (is_flagged((float)station.observ[index_prop][i], cdf.fill_value) 
	       || is_flagged((float)station.observ[index_prop][i], HB_f_mask)) ) {
                 for (j = 0; j < MAXPROP; ++j) {
                   if (prop_needed[j] && prop_avail[j] && (j != (int)DE)) {
                      station.observ[j][k] = station.observ[j][i];
		      if (is_flagged((float)station.observ[j][k], cdf.fill_value) || is_flagged((float)station.observ[j][k], HB_f_mask)) {
		        station.observ[j][k] = HB_MISSING;
		      }
                      if (cdf.counts_included )
                         count[j][k] = count[j][i];
                   }
                 }
                 station.observ[(int)DE][k] = (double) z[i];      
                 if (cdf.counts_included)
                     count[(int)DE][k] = count[index_prop][i];
                 ++k;
              }
           }
           if ((npts = k) <= 0)
                continue;

           /* next, compute properties that were flagged. Set count for
              each computed property   */
      
/*  !**! Special cases for individual properties */
           for (i = 0; i < MAXPROP; ++i) {
              if (compute_flag[i]) {
                  switch ((enum property) i) {
                     case OX:
                        for (j=0; j < npts; ++j) {
                           station.observ[i][j] = ox_kg2l(station.observ[(int)O2][j], station.observ[(int)PR][j], station.observ[(int)TE][j],station.observ[(int)SA][j]);
                        }
                        if (cdf.counts_included) {
                            for (j = 0; j < npts; ++j) {
                               count[i][j] = count[(int)O2][j];
                            }
                        }
                        break;
                      case O2:
                        for (j=0; j < npts; ++j) {
                           station.observ[i][j] = ox_l2kg(station.observ[(int)OX][j], station.observ[(int)PR][j], station.observ[(int)TE][j],station.observ[(int)SA][j]);
                        }
                        if (cdf.counts_included) {
                            for (j = 0; j < npts; ++j) {
                               count[i][j] = count[(int)OX][j];
                            }
                        }
                        break;
                   case S0: 
                         compute_sigma(0., npts, station.observ[(int)S0], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA]);
                        if (cdf.counts_included) {
                            for (j = 0; j < npts; ++j) {
                               count[i][j] = count[(int)TE][j];
                            }
                        }
                         break;
                     case S1: 
                         compute_sigma(1000., npts, station.observ[(int)S1], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA]);
                        if (cdf.counts_included) {
                            for (j = 0; j < npts; ++j) {
                               count[i][j] = count[(int)TE][j];
                            }
                        }
                         break;
                     case S2:   
                         compute_sigma(2000., npts, station.observ[(int)S2], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA]);
                        if (cdf.counts_included) {
                            for (j = 0; j < npts; ++j) {
                               count[i][j] = count[(int)TE][j];
                            }
                        }
                         break;
                     case S3: 
                         compute_sigma(3000., npts, station.observ[(int)S3], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA]);
                        if (cdf.counts_included) {
                            for (j = 0; j < npts; ++j) {
                               count[i][j] = count[(int)TE][j];
                            }
                        }
                         break;
                     case S4: 
                         compute_sigma(4000., npts, station.observ[(int)S4], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA]);
                        if (cdf.counts_included) {
                            for (j = 0; j < npts; ++j) {
                               count[i][j] = count[(int)TE][j];
                            }
                        }
                         break;
                     case S_: 
                         compute_sigma(s_pref, npts, station.observ[(int)S_], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA]);
                        if (cdf.counts_included) {
                            for (j = 0; j < npts; ++j) {
                               count[i][j] = count[(int)TE][j];
                            }
                        }
                         break;
                     case TH:
                         compute_theta(npts, station.observ[(int)TH], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA]);
                        if (cdf.counts_included) {
                            for (j = 0; j < npts; ++j) {
                               count[i][j] = count[(int)TE][j];
                            }
                        }
                         break;
                     case HT:
                         compute_height(npts, station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA], ht_pref, station.observ[(int)HT]);
                        if (cdf.counts_included) {
                            for (j = 0; j < npts; ++j) {
                               count[i][j] = count[(int)TE][j];
                            }
                        }
                         break;
                     case PE:
                         compute_energy(npts, station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA], pe_pref,  station.observ[(int)PE]);
                        if (cdf.counts_included) {
                            for (j = 0; j < npts; ++j) {
                               count[i][j] = count[(int)TE][j];
                            }
                        }
                         break;
                     case SV: 
                         compute_sp_vol(npts, station.observ[(int)SV], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA]);
                        if (cdf.counts_included) {
                            for (j = 0; j < npts; ++j) {
                               count[i][j] = count[(int)TE][j];
                            }
                        }
                         break;
                     case VA: 
                         compute_svan( npts, station.observ[(int)VA], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA]);
                        if (cdf.counts_included) {
                            for (j = 0; j < npts; ++j) {
                               count[i][j] = count[(int)TE][j];
                            }
                        }
                         break;

                     case GN: 
		         if (!compute_flag[(int)GE]) { 
                            compute_gamma_n( &ginfo, npts, station.observ[(int)GN], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA], dlon, dlat);
                           if (cdf.counts_included) {
                             for (j = 0; j < npts; ++j) {
                               count[i][j] = count[(int)TE][j];
                             }
                           }
			 }
                         break;

                     case GE:
	                compute_gamma_nerr(&ginfo, npts, station.observ[(int)GN],   
		          station.observ[(int)GE], station.observ[(int)PR],
		          station.observ[(int)TE], station.observ[(int)SA], dlon, dlat);
                        if (cdf.counts_included) {
                            for (j = 0; j < npts; ++j) {
                               count[(int)GE][j] = count[(int)TE][j];
                               count[(int)GN][j] = count[(int)TE][j];
                            }
                        }
	                break;

                     default:
                         break;
              
                  }  /* end switch */

              }
           }

         /* traverse the linked list of surfaces & interpolate data onto each 
            surface */

           surf = listptr;
           while (surf != NULL) {

                 insertdata(surf, sq, npts, dlat, prop_avail);

              surf = surf->next;
           }  /* end while */

       }  /* end for col */
   }  /* end for row */

   free((void *) z);
   free((void *) x);
   for (i = 0; i < MAXPROP; ++i) {
       if (station.observ[i] != NULL) {
          free((void *)station.observ[i]);
          station.observ[i] = NULL;
       }
       if (cdf.counts_included && (count[i] != NULL)) {
          free((void *)count[i]);
          count[i] = NULL;
       }
   }
   if (cdf.counts_included)
      free((void *)count);

   return;
}  /* end get_cdf_data() */
/****************************************************************************/
/***********************************************************************/
double msf(p, ht, sv, pbar)
double p;    /* pressure in dbars */
double ht;   /* dynamic height in dynamic meters =  m**2/s**2 * 10^-1 */
double sv;   /* specific volume anomaly:  m**3/kg * 10^8 */
double pbar; /* average pressure for entire surface in dbars  */
{
  /* convert dbars to Pa; sv to m**3/kg; and dynamic meters to m**2/s**s */

   return( -(1.0e-4 * (p - pbar) * sv - 10.0 * ht));
}

