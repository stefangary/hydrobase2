/*  hb_linew_avg_section.c
................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             2006
               ...................................................
..........................................................................
. Isopycnally averages properties onto a regularly spaced grid
. as a function of distance (x) and pressure or depth (y), aligning the Gulf Stream
. onto a specified grid coordinate in all sections.  The y
. dimension consists of a series of standard depths or pressures determined
.  by the y-increment specified (with -I).  The averaging is done 
. on isopycnal (sigma) surfaces (definitions specified with -S) and interpolated back onto the standard levels. 
. The section is divided into 2 segments:  DWBC and Gulf Stream . 
.  Observed profiles are manually separated and input as 2 separate files:
   -1:  DWBC ( usually standard profiles 9001 - 9014). Profiles between Gulf Stream north wall and 300 km have 
      upper water column stripped off (usually at LSW/LNADW interface gamma-n=27.983)
   -2:  profiles located at the Gulf Stream's north wall and southward.
   A third file gives the position of the Gulf Stream north wall for each cruise.  This must be created
   separately and contains:
   -3:   ship cruise lat lon (separated by spaces, one line per synoptic crossing of Gulf Stream)
.  DWBC profiles are gridded as a function of distance from  40.3N, 70.2W.
.   For the Gulf Stream portion all profiles are gridded as a function
.  of distance from the  GS north wall in each synoptic section. For the default averaged section, 
.  the DWBC spans the range 0-300 km, and the Gulf Stream starts at 240 km
.  Note that these segments overlap, because the deep Gulf Stream is offset relative to the
   surface flow.  In manually creating the DWBC and Gulf Stream files (1 & 2 above) the deep layers
   under the Gulf Stream should be separated and added to the DWBC file, (vice versa for the portion of
   the DWBC underlying the Gulf Stream. These segment lengths can be alternately specified with -B<xmin1//xmax1/xmin2/xmax2>
..........................................................................
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hydrobase.h"
#include "hb_gamma.h"
#include "hb_paths.h"

#define   RADperDEG 0.017453292             /* pi/180 */
#define  EarthRadius  6371.0    /* in km */
#define PI  3.14159

#define    PRINT_MSG  1          /* turn on message printing to stderr */
#define    MIX_LAYER_DEF 0.02    /* default sigma range for mixed layer */
#define    HBEMPTY -0.9e035   /* define here as large negative value */


struct GAMMA_NC ginfo;   /* used for neutral density */

/***************************************************************/
/*  define the standard sigma levels on which data will be averaged;
     sigma series is divided by reference levels ...*/

#define  MAXSIGLEVS 800    /* max size of sigma series -- arbitrary number */
#define  NREFLEVS  5       /* # of ref levs defining sigma values */

double siglevs[MAXSIGLEVS];   /* stores the sigma series  */
int  nsiglevs;                /* number of elements in full sigma series */

int ref_prop_id[NREFLEVS] = {(int)S0,(int)S1,(int)S2,(int)S3,(int)S4};
int nsig[NREFLEVS];       /* number of elements of sigseries for each ref lev */
int zmin[NREFLEVS] = {-1,500,1500,2500,3500}; /* depth ranges for ref levs */
int zmax[NREFLEVS] = {500,1500,2500,3500,99999};
double *sigptr[NREFLEVS];   /* points to station arrays of referenced sigma values */
FILE *sigfile[NREFLEVS];

/***************************************************************/
/* define standard depth levels on which data will be averaged for determining variance. */

#define MAX_DEPTH_LEVS  600
double std_depth[MAX_DEPTH_LEVS];
int nstddepths;          /* depends upon yinc specified by user */

int isobaric_avg;
int dflag;    /* depth, not pressure, is y-dimension */
FILE *topofile;


/***************************************************************/

/* globally referenced variables for input station data */

struct HYDRO_DATA sta;
struct HYDRO_HDR hdr;
double *pr, *de, *te, *sa;
double s_pref;
double pe_pref;
double ht_pref;

double mix_layer_def;      /* option to specify mixed layer definition */

/* globally referenced grid units */

#define  NSEGS  2      /* # of segments comprising entire section */
#define NFILES 3     /* # of input files containing profiles */
double xinc, yinc;
double xlength;          /*  length scale for weight function = e^[-(d/L)^2] */
double xmin[NSEGS] = {0, 240};
double xmax[NSEGS] = {300, 460};
int nx[NSEGS];
double *xdist;    /* vector of distance values at gridpts along x-axis */
double *xtopo;   /* vector of seafloor depth or pressure along x-axis */
double lat0 = 40.3;
double lon0 = -70.2;  /* lat/lon of start of section */



/************ structures to represent grid node ****************/

struct surfrec {
          double  density, depth;
          double *prop;
	  double *prop_sqsum;
	  double *wghtsum;
              int  *count;
	      int   n;
  struct surfrec *next;
};

struct deepestrec {
          double  depth, pressure;
	  double sig_0, sig_1, sig_2, sig_3, sig_4;
	  int count;
          struct deepestrec *next;
};

struct gridnode {
         double **prop;
	 double **prop_sqsum;
	 double **weightsum;
         double  *d;
	 double  *dweightsum;
	 double **d_prop;
	 double **d_prop_sqsum;
	 double **d_weightsum;
             int **count;
	     int **d_count;
             int  *nobs;
 struct deepestrec *deepest;
 struct surfrec  *mix_layer;
};

struct profile_loc {
        char ship[3];
	int cruise;
	float lat, lon;
	struct profile_loc *next;
};
/***************************************************************/
/*  prototypes for functions defined within this module... */

struct gridnode *alloc_grid(int, int, int);

int get_sig_series(int, FILE *,  double *);
int parse_p_option(char *, int *);
int do_std_depth(struct gridnode *, int, int *, int, int);
int mixed_layer_vals(int, int, struct surfrec *, double *, double *, double *, double *,  double *,  int, int, int);
int monotonic(double *, int  *, double *, double **, double **, double **, int  **, int, int);
int interpolate_count(double, double *, int *, int);

void print_usage(char *);
void free_grid(int, int, struct gridnode *);
void parse_s_option(char *);
void get_topo(FILE *, int);
void insert_data(double *,int, double *, double *, double *, double, int *, int);
void insert_data_on_depthlevs(double *, int, double *, double *, double *, double, int *);
void insert_surf_vals(struct gridnode *,int *,int, double);
void insert_deepest_vals(struct gridnode *);
void compute_avg(double *,double *, double *, int *, int);
void define_bottom(struct gridnode *, int);
void define_avg_mixed_layer(struct gridnode *, int, int);
void do_middle_segment(struct gridnode *, int, int);
void define_sigma_transitions(struct gridnode *, int);
void delete_list(struct deepestrec *);
void delete_surfrec_list(struct surfrec *);

struct surfrec *get_surfrec(struct surfrec *, double, int);
struct surfrec *create_surfrec(int);
struct surfrec *get_density_rec(int, struct surfrec *, int);
struct deepestrec *create_deepestrec();
struct profile_loc *insert_GS_rec(struct profile_loc *, FILE *);
struct profile_loc *find_GS_rec(struct profile_loc *, struct HYDRO_HDR *);

double interpolate(double, double *, double *, int);
double get_weight(double, double); 


/***************************************************************/

int main (int argc, char **argv)
{
   short   bopt, iopt, popt;
   int     curfile = 1, nfiles; 
   int     print_mess; 
   int     nlevs, nprops, *prop_indx;
   int     i, j, error, test_outcrop, prop_avail;
   int     seg_index, nx_section, iseg, iend, istart, in_bounds;
   FILE   *outfile, *GSpos_file;
   char   *st;
   int     infile[NFILES];
   double dx, dy, dist;
   double weight;
   struct gridnode *grid[NSEGS], *section;
   struct profile_loc *GS_listptr, *GS_north;


 
/*  set these default values */

    nfiles = 0;
    nprops = 0;
    popt =  0;
    error = 0;
    dflag = 0;
    print_mess = PRINT_MSG;
    mix_layer_def = MIX_LAYER_DEF;
    xlength = 30.0;   /* L in km for distance weighting = e^- PI *(d/L)^2 */
    xinc = 20.0;  /* km */
    yinc = 10.0;  /* meters or decibars*/
    isobaric_avg = 1;
    topofile = NULL;
    xtopo = NULL;

    GS_listptr = NULL;
    for (i = 0; i < NFILES; ++i)
       infile[i] = 0;
    
    for (j = 0; j < NREFLEVS; ++j)
       sigfile[j] = NULL;
       
     prop_indx = (int *) calloc(MAXPROP, sizeof(int));    
    
/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
   
/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
	          
               case '1':
	           if ( (infile[0] = open_hydro_file("", &argv[i][2], "", print_mess)) < 0) {
		      fprintf(stderr,"Unable to open %s for input\n", &argv[i][2]);
		      error = 1;
		   }
		   ++nfiles;
	           break;
	       case '2':
	           if ( (infile[1] = open_hydro_file("", &argv[i][2], "", print_mess)) < 0) {
		      fprintf(stderr,"Unable to open %s for input\n", &argv[i][2]);
		      error = 1;
		   }
		   ++nfiles;
	           break;
	       case '3':
	           if ( (GSpos_file = fopen(&argv[i][2], "r") ) == NULL) {
		      fprintf(stderr,"Unable to open %s for input\n", &argv[i][2]);
		      exit(1);
		   }
		   GS_listptr = insert_GS_rec(GS_listptr, GSpos_file);

	           break;
               case 'B':                    /* get bounds for section segments */
	       
	                error = (sscanf(&argv[i][2], "%lf/%lf/%lf/%lf", &xmin[0], &xmax[0], &xmin[1], &xmax[1]) == 4) ? 0 : 1;
                        break;
                case 'D':    
		        dflag = 1;              
                        break;

               case 'I':
                        error = (sscanf(&argv[i][2],"%lf", &xinc) == 1) ? 0 : 1;
                        st = &argv[i][2];
                        while (*(st++) != '\0') {
                            if (*st == '/') {
                              ++st;
                              error = (sscanf(st,"%lf", &yinc) == 1) ? 0 : 1;
                              break;
                            }
                        }
                        
                        break;
               case 'L':
                        error = (sscanf(&argv[i][2],"%lf", &xlength) == 1) ? 0 : 1;
                        break;
                        
               case 'M':
                        error = (sscanf(&argv[i][2],"%lf", &mix_layer_def) == 1) ? 0 : 1;
                        break;

                case 'N':                   /* turn off isobaric averaging output */
		        isobaric_avg = 0 ;              
                        break;
                        

               case 'P':
                        popt = 1;
                        nprops = parse_p_option(&argv[i][2], prop_indx);
                        break;

               case 'S':
                        parse_s_option(&argv[i][2]);
                        break;
                 
               case 'T':
                        topofile = fopen(&argv[i][2], "r");
			if (topofile == NULL) {
			   fprintf(stderr,"Unable to open %s for reading.\n", &argv[i][2]);
			   exit(1);
			}
			fprintf(stderr,"Opened %s \n", &argv[i][2]);
                        break;
               case 'X':                    /* set lat/lon for start of section*/
	                error = (sscanf(&argv[i][2], "%lf/%lf", &lat0, &lon0) == 2) ? 0 : 1;
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
             print_usage(argv[0]);
             exit(1);
          }

       }  /* end if */
       else {
             fprintf(stderr,"\nError parsing command line args.\n");
             fprintf(stderr,"     in particular: '%s'\n", argv[i]);
             print_usage(argv[0]);
             exit(1);

       }
   }  /* end for */

   if ( !nfiles || !popt || sigfile[0] == NULL) {
       fprintf(stderr,"\nYou must specify input file(s), properties and siglist filename\n");
       print_usage(argv[0]);
       exit(1);
   }
   if (nfiles > 1 && GS_listptr == NULL) {
       fprintf(stderr,"\n A file containing the position of the Gulf Stream north wall is necessary \n");
       print_usage(argv[0]);
       exit(1);
   }
     
/* define values for sigma series ... */

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
   fprintf(stderr,"\n total number of sigma levels: %d\n", nsiglevs);

  /*  define std depth levels from yinc */
  
     nstddepths = NINT(6000.00001 / yinc);
     for (i = 0; i < nstddepths; ++i) {
        std_depth[i] = i * yinc;
     }
 
/* compute dimensions of grid and allocate space for computation ...*/
    
   for (i = 0; i < NSEGS; ++i) {
      nx[i] = NINT((xmax[i] + 0.0001 - xmin[i]) / xinc);
   }
   nx_section =  NINT( (xmax[1] + 0.0001 - xmin[0]) / xinc);
   section = alloc_grid(nx_section, nsiglevs, nprops);
   grid[0] = &section[0];                 /* set pointers to start of each section division, these are overlapping */
   grid[1] = &section[nx_section - nx[1]];

   xdist = (double *) calloc(nx_section, sizeof(double));   
   for (i = 0; i < nx_section; ++i)     /* store vector of x (distance) values */
      xdist[i] = xmin[0] + i * xinc;
      
   get_topo(topofile, nx_section);  
   
   fprintf(stderr,"xlength for distance weighting [= e^-[ pi *(d/L)^2]  is %.1lf km\n", xlength);
   
/* initialize these... */

   for (i = 0; i < MAXPROP; ++i)
      sta.observ[i] = (double *) NULL;
   hdr.prop_id = (int *) NULL;

/* loop for each input_file */

   for (curfile = 0; curfile < NFILES; ++curfile) {
       if (infile[curfile] == 0) 
          continue;
	  
	  
     /* loop for each station */
    while ((error = get_station(infile[curfile], &hdr, &sta)) == 0) {
    
       /* determine which gridnodes are within xlength of this profile.
       Compute offset distance from start of each segment */
       in_bounds = 1;
       switch (curfile) {
          case 0:               /* Profiles in DWBC */
	      iseg = 0;
	      dx = (hdr.lon - lon0)*RADperDEG * cos(RADperDEG*0.5 * (hdr.lat + lat0)) * EarthRadius;
	      dy = (hdr.lat - lat0) * RADperDEG * EarthRadius;
	      dist = sqrt(dx*dx + dy*dy);
	      seg_index = NINT((dist + 0.0001) / xinc);         /* offset from start of segment */
	      istart = NINT((dist - xlength + .0001) / xinc);
	      if (istart < 0)
	           istart = 0;
	      if (istart >= nx[0])
		      in_bounds = 0;    
	      iend = NINT((dist + xlength + .0001) / xinc);
	      if (iend >= nx[0] )
	           iend = nx[0] - 1;
	       if (iend < 0)
	            in_bounds = 0;
	     break;
          case 1:              /* Profiles in GULF STREAM and southward */
	      iseg = 1;
	      GS_north = find_GS_rec(GS_listptr, &hdr);
	      if (GS_north == NULL) {
	         in_bounds = 0;
	      }
	      else {
	         dx = (hdr.lon - GS_north->lon)*RADperDEG * cos(RADperDEG*0.5 * (hdr.lat + GS_north->lat)) * EarthRadius;
	         dy = (hdr.lat - GS_north->lat) * RADperDEG * EarthRadius;
	         dist = sqrt(dx*dx + dy*dy);
		 seg_index = NINT((dist + .0001) / xinc);        /* offset from start of segment */
	         istart = NINT((dist - xlength + .0001) / xinc);
	         if (istart < 0)
	              istart = 0;
		  if (istart >= nx[iseg])
		      in_bounds = 0;    
	        iend = NINT((dist + xlength + .0001) / xinc);
	         if (iend >= nx[iseg] )
	             iend = nx[iseg] - 1;
	         if (iend < 0)
	            in_bounds = 0;
	      }
	     break;
      
       } /* end switch */
        
       if (in_bounds) {
         pr = sta.observ[(int)PR];   /* set frequently used pointers */
         de = sta.observ[(int)DE];
         te = sta.observ[(int)TE];
         sa = sta.observ[(int)SA];
         if (pr == NULL || te == NULL || sa == NULL || de == NULL) {
           fprintf(stderr, "\nThe requisite parameters: pr, te, sa, de are not available at this station.\n");
           write_hydro_hdr(STDERR, &hdr);
           continue;
         }

        /* explicitly set the first observation to pressure = depth = 0 if it is within 15 meters
	 of the sea surface */
         if (pr[0] > 0.0 && pr[0] < 15.0)  {
	       pr[0] = 0.0; 
	       de[0] = 0.0;  
	  } 
         /* compute referenced sigmas for this station... */

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

         /* compute other appropriate properties for this station ... */

         for (i = 0; i < nprops; ++i) {
            test_outcrop = 0;
            prop_avail = 1;

         /* !**! Special cases for individual properties... */
            switch ((enum property) prop_indx[i]) {
               case PR:
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
                        sta.observ[O2][j] = ox_l2kg(sta.observ[(int)OX][j], pr[j], te[j], sa[j]);
                      }
                      prop_avail = 1;
                    }
                  }
                  break;
               case TH:
                  free_and_alloc(&sta.observ[(int)TH], hdr.nobs);
                  compute_theta(hdr.nobs, sta.observ[(int)TH], pr, te, sa);
                  break;

               case PV:
                  free_and_alloc(&sta.observ[(int)PV], hdr.nobs);
                  buoy_freq(sta.observ[(int)PV],pr,te,sa,hdr.nobs,100,10);
                  for (j = 0; j < hdr.nobs; ++j) 
                    sta.observ[(int)PV][j] *= sta.observ[(int)PV][j];
                  po_vort(sta.observ[(int)PV],sta.observ[(int)PV], hdr.nobs, (double)hdr.lat);
                  break;
               default:
                   prop_avail = available((enum property) prop_indx[i], &hdr);
                   break;
            } /* end switch */
	    
            for ( j = istart; j <= iend; ++j) {  
	       /* insert data into each gridnode that falls within xlength distance */
	       weight = get_weight(ABS(j * xinc  - dist), xlength);

               if (prop_avail) {
                  insert_data(sta.observ[prop_indx[i]], hdr.nobs, grid[iseg][j].prop[i],  grid[iseg][j].prop_sqsum[i], grid[iseg][j].weightsum[i], weight, grid[iseg][j].count[i], test_outcrop);
                  insert_data_on_depthlevs(sta.observ[prop_indx[i]], hdr.nobs, grid[iseg][j].d_prop[i],  grid[iseg][j].d_prop_sqsum[i], grid[iseg][j].d_weightsum[i], weight, grid[iseg][j].d_count[i] ) ;
	       }
            } /* end for j */
         } /* end for i < nprops*/

      /* insert depth values at each gridnode within xlength... */

          for ( j = istart; j <= iend; ++j) {   /* insert data into each gridnode that falls within xlength scale */
	       weight = get_weight(ABS( j * xinc - dist), xlength);
                test_outcrop = (de[0] <= 51.) ? 1 : 0;
                insert_data( de, hdr.nobs, grid[iseg][j].d, (double *)NULL,  grid[iseg][j].dweightsum, weight,
                            grid[iseg][j].nobs, test_outcrop);
	        insert_surf_vals(&grid[iseg][j], prop_indx, nprops, weight);
          } /* end for j */

     /* add deepest vals only for the gridnode in which this profile is located  */	  
          if (seg_index >= 0 && seg_index < nx[iseg])
                 insert_deepest_vals(&grid[iseg][seg_index]);
	  
       } /* end if in_bounds */
       
       /* free up space... */

       for (i = 0; i < MAXPROP; ++i) {
          if (sta.observ[i] != NULL) {
             free((void *) sta.observ[i]);
             sta.observ[i] = (double *) NULL;
          }
       }

     }  /*end while !eof */ 
     if (error > 0) {
         report_status(error, stderr);
         exit(1);
   }

     close(infile);

   } /* end for curfile */             

   fprintf(stderr,"  computing averages ...\n");
   
   define_bottom(section, nx_section);
   define_avg_mixed_layer(section, nx_section, nprops);

/* for each gridnode, compute means at all sigma levels and std depth levels 
    label each level with neutral density   */
    
    
     gamma_nc_init((char *) GAMMA_NC_PATH, &ginfo);

   for (j = 0; j < nx_section; ++j) {
            for (i = 0; i < nprops; ++i) {
               compute_avg(section[j].prop[i], section[j].prop_sqsum[i], section[j].weightsum[i], section[j].count[i], nsiglevs);
               compute_avg(section[j].d_prop[i], section[j].d_prop_sqsum[i], section[j].d_weightsum[i], section[j].d_count[i], nstddepths);
            }
            compute_avg(section[j].d, (double *)NULL, section[j].dweightsum, section[j].nobs,  nsiglevs);
            define_sigma_transitions(&section[j], nprops);
            nlevs = do_std_depth(&section[j], nsiglevs, prop_indx,  nprops, j);
    }
   exit(0);
} /* end main */


/****************************************************************************/

void print_usage(char *program)
{
   
   fprintf(stderr,"\nUsage:  %s filename_root(s)  -1<dwbc_file>  -2<gulfstream_south_file>  -3<gulfstream_northwall_positions> -P<list_of_properties> -S<sigseries_root_name>  [-B<xmin0/xmax0/xmin1/xmax1> ] [-I<xinc>/<yinc>]   [-D]   [-L<xlength>]  [-M<mixed_layer_def>] [-N]  [-T<topofile>]  [-X<lat0/lon0>]  \n", program);
    fprintf(stderr,"\n -1   File containing DWBC stations;");
    fprintf(stderr,"\n -2   File containing stations south of Gulf Stream north wall;");
    fprintf(stderr,"\n -3   File identifying position of Gulf Stream north wall for each synoptic section.");
    fprintf(stderr,"\n       it contains one line per section  listing ship cruise lat lon ");
  
   fprintf(stderr,"\n -P   list of properties to project onto surface;");
   fprintf(stderr,"\n          ex:  -Ppr/th/sa/ox");
   fprintf(stderr,"\n               -P (by itself) produces a list of available properties\n");
   fprintf(stderr,"\n -S   rootname (relative to current directory) of 5 sigma series files ");
   fprintf(stderr,"\n          filenames are expected to be <root>.sig0list, <root>.sig1list, etc");
   fprintf(stderr,"\n          for sigma 0 , 1, 2, 3 and 4.");
   fprintf(stderr,"\n          Each line in file contains sigmin, sigmax, incr");
   fprintf(stderr,"\n          (sigmin and sigmax are INCLUSIVE in generating series).");
   fprintf(stderr,"\n          Values MUST be monotonically INCREASING.");
   fprintf(stderr,"\n\n OPTIONS:");
   fprintf(stderr,"\n -B   specifies x bounds (in km) for DWBC and Gulf Stream segments.");
   fprintf(stderr,"\n        Defaults:  segment 1 = %.1lf/%.1lf", xmin[0], xmax[0]);
   fprintf(stderr,"\n                      segment 2 = %.1lf/%.1lf", xmin[1], xmax[1]);
   fprintf(stderr,"\n -D  y-dimension for output is standard depth levels. (default is pressure levels.");
   fprintf(stderr,"\n -I   specifies x/y grid increments.  Defaults: %.1lf/%.1lf", xinc, yinc);
   fprintf(stderr,"\n -L specify decorrelation scale (L, in km) for distance weighting =  e ^[-pi * (dist/L)^2]");
   fprintf(stderr,"\n    weight is 1 at distance = 0,  ~0.1 at this lengthscale and drops off rapidly beyond this point" );
   fprintf(stderr,"\n     -L0 turns off distance weighting (all points contribute equally to mean)" );
   fprintf(stderr,"\n          default is L =%.1lf  km", xlength);
   fprintf(stderr,"\n -M  specify the definition of the mixed layer as ");
   fprintf(stderr,"\n          sigma(sea surface) + this value.  default: %f ", MIX_LAYER_DEF);
   fprintf(stderr,"\n          ex:  -M.02");
   fprintf(stderr,"\n -N  do NOT output isobarically averaged values");
   fprintf(stderr,"\n -T specifiy file with distance, topography (either depth or pressure) values to define seafloor along section");
   fprintf(stderr,"\n -X  set lat/lon for start of section. Defaults: %.3lf/%.3lf", lat0, lon0);
   fprintf(stderr,"\n -h  help ... prints this message.");
   fprintf(stderr,"\n\n");  
   return;
} /* end print_usage() */
/****************************************************************************/

void parse_s_option(char *root)
   /* get the sigma list files defining the isopycnal surfaces */
{
  int j;
  char fname[1000];

   for (j = 0; j <= 4; ++j) {
   strcpy(fname, root);
    switch (j) {
       case 0: 
          strcat(fname, ".sig0list");
          break;
       case 1: 
          strcat(fname, ".sig1list");
          break;
       case 2: 
          strcat(fname, ".sig2list");
          break;
       case 3: 
          strcat(fname, ".sig3list");
          break;
       case 4: 
         strcat(fname, ".sig4list");
          break;
    }
    sigfile[j] = fopen(fname,"r");
    if (sigfile[j] == NULL) {
      fprintf(stderr,"\nError opening %s \n", fname);
      exit(1);
    }
  } /* end for */
   return;
}   /* end parse_s_option() */

/*****************************************************************************/
int parse_p_option(char *st, int *prop_indx)
{
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
        fprintf(stderr,"\nWARNING neutral density (gamma-n) is not an appropriate property to average isopycnally.");
        fprintf(stderr,"\nCompute from averaged pr,te,sa output by hb_grid3d.\n");
     }
     
     if ((prop_indx[n] != (int) GE) && (prop_indx[n] != (int) GN))  {
            ++n;
     }

   } while (*st == '/');
   
   
   return (n);
}  /* end parse_p_option() */

/*****************************************************************************/
struct gridnode *alloc_grid(int nx, int nz, int nprops)
{  /* allocates memory for averaging properties on sigma-levels and standard depth levels  */
   int i, j, k, n;
   struct gridnode *g;

         g = (struct gridnode *) malloc(nx * sizeof(struct gridnode));
         if (g == NULL) {
           fprintf(stderr,"\nUnable to allocate memory for grid[%d]\n", i);
           exit(1);
         }
         for (j = 0; j < nx; ++j) {

            g[j].prop = (double **) malloc(nprops * sizeof(double *));
            if (g[j].prop == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].prop\n", i, j);
                exit(1);
            }
            g[j].prop_sqsum = (double **) malloc(nprops * sizeof(double *));
            if (g[j].prop_sqsum == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].prop\n", i, j);
                exit(1);
            }
            g[j].weightsum = (double **) malloc(nprops * sizeof(double *));
            if (g[j].weightsum == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].weightsum\n", i, j);
                exit(1);
            }

            g[j].count = (int **) malloc(nprops * sizeof(int *));
            if (g[j].count == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].count\n", i, j);
                exit(1);
            }
            g[j].d_prop = (double **) malloc(nprops * sizeof(double *));
            if (g[j].prop == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].prop\n", i, j);
                exit(1);
            }
            g[j].d_prop_sqsum = (double **) malloc(nprops * sizeof(double *));
            if (g[j].prop_sqsum == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].prop\n", i, j);
                exit(1);
            }
            g[j].d_weightsum = (double **) malloc(nprops * sizeof(double *));
            if (g[j].weightsum == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].weightsum\n", i, j);
                exit(1);
            }

            g[j].d_count = (int **) malloc(nprops * sizeof(int *));
            if (g[j].count == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].count\n", i, j);
                exit(1);
            }

            g[j].d = (double *) malloc(nz * sizeof(double));
            if (g[j].d == NULL) {
               fprintf(stderr,"\nUnable to allocate memory for g[%d][%d].d\n", i, j);
               exit(1);
            }
            g[j].dweightsum = (double *) malloc(nz * sizeof(double));
            if (g[j].dweightsum == NULL) {
               fprintf(stderr,"\nUnable to allocate memory for g[%d][%d].dweightsum\n", i, j);
               exit(1);
            }
            g[j].nobs = (int *) malloc(nz * sizeof(int));
            if (g[j].nobs == NULL) {
               fprintf(stderr,"\nUnable to allocate memory for g[%d][%d].nobs\n", i, j);
               exit(1);
            }

            g[j].mix_layer = (struct surfrec *) NULL;
            g[j].deepest = (struct deepestrec *)NULL;

            for (n = 0; n < nz; ++n) {       /* initialize depth,count array */
                 g[j].d[n] = 0.0;
                 g[j].dweightsum[n] = 0.0;
                 g[j].nobs[n] = 0;
            }

            for (k = 0; k < nprops; ++k) {

              g[j].prop[k] = (double *) malloc(nz * sizeof(double));
              if (g[j].prop[k] == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].prop[%d]\n", i, j, k);
                exit(1);
              }
              g[j].prop_sqsum[k] = (double *) malloc(nz * sizeof(double));
              if (g[j].prop_sqsum[k] == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].prop[%d]\n", i, j, k);
                exit(1);
              }
              g[j].weightsum[k] = (double *) malloc(nz * sizeof(double));
              if (g[j].weightsum[k] == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].weightsum[%d]\n", i, j, k);
                exit(1);
              }

              g[j].count[k] = (int *) malloc(nz * sizeof(int));
              if (g[j].count[k] == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].count[%d]\n", i, j, k);
                exit(1);
              }
              
               g[j].d_prop[k] = (double *) malloc(nstddepths * sizeof(double));
              if (g[j].prop[k] == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].prop[%d]\n", i, j, k);
                exit(1);
              }
              g[j].d_prop_sqsum[k] = (double *) malloc(nstddepths * sizeof(double));
              if (g[j].prop_sqsum[k] == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].prop[%d]\n", i, j, k);
                exit(1);
              }
              g[j].d_weightsum[k] = (double *) malloc(nstddepths * sizeof(double));
              if (g[j].weightsum[k] == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].weightsum[%d]\n", i, j, k);
                exit(1);
              }

              g[j].d_count[k] = (int *) malloc(nstddepths * sizeof(int));
              if (g[j].count[k] == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].count[%d]\n", i, j, k);
                exit(1);
              }
                /* zero each prop, count */
              for (n = 0; n < nz; ++n) {        /* for averaging along isopycnals */
                 g[j].prop[k][n] = 0.0;
                 g[j].prop_sqsum[k][n] = 0.0;
                 g[j].weightsum[k][n] = 0.0;
                 g[j].count[k][n] = 0;
              }
	      
	      for (n = 0; n < nstddepths; ++n) {  /* for averaging on depth surfaces */
                 g[j].d_prop[k][n] = 0.0;
                 g[j].d_prop_sqsum[k][n] = 0.0;
                 g[j].d_weightsum[k][n] = 0.0;
                 g[j].d_count[k][n] = 0;
	      }
            } /* end for k*/
         } /* end for j*/

   return g;
} /* end alloc_grid */
/*****************************************************************************/
void free_grid(int nx, int nprops, struct gridnode *gridptr)
  /* frees all the memory for ny * nx struct gridnodes */
{
   int col, k;
   struct gridnode *g;
   
      for (col = 0; col < nx; ++col) {
         g = &gridptr[col];
         for (k = 0; k < nprops; ++k) {
           if (g->prop[k] != NULL)
	     free((void *)g->prop[k]);
           if (g->prop_sqsum[k] != NULL)
	     free((void *)g->prop_sqsum[k]);
           if (g->weightsum[k] != NULL)
	     free((void *)g->weightsum[k]);
           if (g->count[k] != NULL)
	     free((void *)g->count[k]);

           if (g->d_prop[k] != NULL)
	     free((void *)g->d_prop[k]);
           if (g->d_prop_sqsum[k] != NULL)
	     free((void *)g->d_prop_sqsum[k]);
           if (g->d_weightsum[k] != NULL)
	     free((void *)g->d_weightsum[k]);
           if (g->d_count[k] != NULL)
	     free((void *)g->d_count[k]);
         }
	 free((void *)g->prop);
	 free((void *)g->prop_sqsum);
	 free((void *)g->weightsum);
	 free((void *)g->count);
	 free((void *)g->d_prop);
	 free((void *)g->d_prop_sqsum);
	 free((void *)g->d_weightsum);
	 free((void *)g->d_count);
	 free((void *)g->nobs);
	 free((void *)g->d);
	 free((void *)g->dweightsum);
	 delete_surfrec_list(g->mix_layer);
	 delete_list(g->deepest);
      } /* end for col */
         
   free((void *) gridptr);
   
   return;
}  /* end free_grid() */
/*****************************************************************************/
int get_sig_series(int ref_id, FILE *fptr, double *siglist)
  
    /*  the file will be read for (min, max, incr) triplets
       from which a sigma series will be generated and inserted at siglist.
       The number of sigma values is returned. 
       
       arguments:
       ref_id:   defines the sigma variable: (int) enum property 
         fptr:   pointer to file containing list of sigma values OR nil 
      siglist:   starting addr of array to insert values 
   */
{
   double  min, max, incr;
   double  next;
   int i, count;
   int bigcount;


   bigcount = nsiglevs;

     count = 0;
     while( fscanf(fptr,"%lf %lf %lf", &min, &max, &incr) == 3) {
         if (++bigcount >= MAXSIGLEVS) {
           return (bigcount);
         }
         *siglist = min;
         ++count;
         while ( ( next = *siglist + incr) <= max) {
           ++count;
           if (++bigcount >= MAXSIGLEVS) { /* avoid a SEG FAULT */
              return (bigcount);
           }
           *(++siglist) = next; 
         }
         ++siglist;
     }

     fclose(fptr);
     return (count);

} /* end get_sig_series */
/*****************************************************************************/
double get_weight( double dist, double L)
   /* Returns weight as a function of distance between points:
         weight = e ^[- PI *(d/L)^2]
  */
{  
   dist = dist /L;
   dist = dist * dist;
   
   return (exp( -PI * dist));
}	
/*****************************************************************************/
void get_topo(FILE *topofile, int nx)
/* Reads an already opened file containing distance, depth values and interpolates
   the bottom depth for each value in the global xdist vector */
{
   char buffer[1000];
   double dist[1000], dep[1000];
   int i, npts, nempty;
   
   if (topofile == NULL)
      return;
      
   if (xdist == NULL) {
      fprintf(stderr, "FATAL ERROR in get_topo():  No xdist vector defined\n");
       return;
   }
   
   if (nx <= 0) {
       fprintf(stderr, "FATAL ERROR in get_topo():  number of nx values is <= 0\n");
       return;
   }
   
   xtopo = (double *) calloc(nx, sizeof(double));
   npts = 0;
   while ( fscanf(topofile, "%[^\n]", buffer) == 1) {
     if (strchr(buffer, '>') == NULL) {
        if ( sscanf(buffer, "%lf %lf", &dist[npts], &dep[npts]) == 2) 
	   ++npts;
     }
     getc(topofile);  /* move past LF */
   }
   fclose(topofile);
   
   nempty = 0;
   for (i = 0; i < nx; ++i) {
      xtopo[i] = hb_linterp(xdist[i], dist, dep, npts);
      if (xtopo[i] < -9998.0) {
        ++nempty;
      }
   }
   
   if (nempty > 0) {
      fprintf(stderr,"WARNING in get_topo(): Unable to interpolate %d topo value(s) along the section.", nempty);
   }
   
   return;

} /* end get_topo() */
/*****************************************************************************/
struct profile_loc *insert_GS_rec(struct profile_loc *GS_listptr, FILE *gsfile)
{
   struct profile_loc *recptr, *r1ptr;
   char  buffer[1000];
   int i, error;

   i = 0;
   while (fscanf(gsfile,"%[^\n]", buffer) != EOF) {
       getc(gsfile);
       ++i;
       recptr = (struct profile_loc *) calloc( 1, sizeof(struct profile_loc));
       error = sscanf(buffer, "%2s %d %f %f", &recptr->ship, &recptr->cruise, &recptr->lat, &recptr->lon) != 4;
       if (error) {
          fprintf(stderr,"Error reading Gulf Stream north wall position at line %d\n", i);
	  exit(1);
       }
      recptr->next = (struct profile_loc *) NULL;

      if (GS_listptr == NULL)
         GS_listptr = recptr;
      else
          r1ptr->next = recptr;
	  
      r1ptr = recptr;
       
   } /* end while */
   
   fclose(gsfile);
   return(GS_listptr);
      
}/* end insert_GS_rec() */
/*****************************************************************************/
struct profile_loc *find_GS_rec(struct profile_loc *GS_listptr, struct HYDRO_HDR *hptr)
/* recursively searches a linked list for a match of ship and cruise and returns
    a pointer to the record.  Returns NULL if no match can be found */
{
   int found;
   
   if (GS_listptr == NULL)
      return ((struct profile_loc *)NULL);
      
   found = (GS_listptr->cruise == hptr->cruise) ? 1 : 0;
   
   if (found) {
      found =  (strncmp(GS_listptr->ship, hptr->ship, 2) == 0);
   }
   
   if (! found)  {   /* recursive part */
 
      GS_listptr = find_GS_rec(GS_listptr->next, hptr);
   }
          
   return (GS_listptr);

}  /* end find_GS_rec() */

/*****************************************************************************/
void insert_data(double * y, int ny, double *ysig, double *ysig_sqsum, double *wghtsum, double wght, int *count, int test_outcrop)

/* For the y property at a station with ny observation levels, interpolate
   to find the yvals at each level in the sigma series,  weight each
   yval, add it to its appropriate sum, and squared sum, sum the weights and increment the counter.  Checks for vertical
   data gaps using the globally defined de array and does not interpolate over
   gaps which exceed 200 m in the thermocline (upper 1000 m) or 1000 m 
   elsewhere.  The sigma series is specified in the global array, siglevs.
   The observed sigma values at this station are subdivided by ref levels, 
   and stored at global addresses pointed to by sigptr. 
   
           y:    array of y observed values
          ny:    dimension of y 
        ysig:    sum of distance-weighted yvalues on each sigma surface
  ysig_sqsum:    sum of distance-weighted y_squared values on each sigma surface
     wghtsum:    sum of weights
       wght:     weight for this profile
       count:    counts # of observations on each surface 
test_outcrop:    0 or 1: if set, tests for outcropping surfaces
*/
{
   int  i, j, datagap, n;
   int  refindx;
   double *xtmp[NREFLEVS], *ytmp, *dtmp, z;
   double reflev;
   
   for (i = 0; i < NREFLEVS; ++i )
        xtmp[i] = (double *) malloc((int) (ny * sizeof(double)));
   ytmp = (double *) malloc((int) (ny * sizeof(double)));
   dtmp = (double *) malloc((int) (ny * sizeof(double)));
   if (dtmp == NULL) {
       fprintf(stderr,"\nUnable to allocate memory in insert_data()\n");
       exit(1);
   }

/* Ensure continuous (no missing values) x,  y, and depth arrays */
   n = 0;
   for (i = 0; i < ny; ++i) {
         if ( y[i] > -8.9)  {
            for (j = 0; j < NREFLEVS; ++j )
               xtmp[j][n] = sigptr[j][i];
            ytmp[n] = y[i];
            dtmp[n++] = de[i];
         }      
   }
   
   if (n <= 1) {               /* not enough values */
     for (i = 0; i < NREFLEVS; ++i)
        free((void *)xtmp[i]);
     free((void *)ytmp);
     free((void *)dtmp);
     return;
   }

   for (i = 0; i < nsiglevs; ++i) {

    /* determine which ref lev corresponds to this sigseries element ... */

      refindx = 0;      
      j = nsig[0];
      while (i > j) {
         if (++refindx == NREFLEVS-1)
             break;
         j += nsig[refindx];
      }

      reflev = refindx * 1000.;
      
      if ((test_outcrop) && (siglevs[i] < xtmp[refindx][0])) {   
        /* note:  this assumes x (density) increases with depth!! */           
             ysig[i] += -1 * wght;     /* use depth= -1 for outcropping surfaces, no need to add to y-squared sum */
	     wghtsum[i] += wght;
             ++count[i];
      }

        /* see if associated depth exists at this station ... */

      else if ((z = interpolate(siglevs[i], xtmp[refindx], dtmp, n)) > -998.) {

            /* now check for datagaps around this depth...*/
            j = 0;
            while ((++j < n) &&( dtmp[j] < z)  ) 
                ;
	    if (j == n)
	        --j;
		
            if ((dtmp[j-1] == z) || (dtmp[j] == z) )
                datagap = 0;
            else if (z < 1001)
                datagap = (dtmp[j] - dtmp[j-1]) > GAP_SHALLOW;
            else
                datagap = (dtmp[j] - dtmp[j-1]) > GAP_DEEP;
                
        /* exclude observations which are more than 1500 m above the reference pressure.
           This prevents alot of bogus values from being added. */
           
           datagap = datagap || (z < reflev-1500);
  
            if (!datagap) {
               if ((z = interpolate(siglevs[i], xtmp[refindx], ytmp, n)) > -998.) {
                  ysig[i] += z * wght;
		  if (ysig_sqsum != NULL)
		            ysig_sqsum[i] += z * z * wght;
		  wghtsum[i] += wght;
                  ++count[i];
               }
            }
      }
   }

   for (i = 0; i < NREFLEVS; ++i)
      free((void *)xtmp[i]);
   free((void *)ytmp);
   free((void *)dtmp);
   return;
} /* end insert_data() */

/*****************************************************************************/
void insert_data_on_depthlevs(double *y, int ny, double *ysum, double *ysqsum, double *wghtsum, double wght, int *count)

/* For the y property at a station with ny observation levels, interpolate
   to find the yvals at each level in the stddepth series,  weight each
   y-value, add it to its appropriate sum, and squared sum, sum the weights and increment the counter.  
   Checks for vertical data gaps using the globally defined de array and does not interpolate over
   gaps which exceed 200 m in the thermocline (upper 1000 m) or 1000 m 
   elsewhere.  The depth series is specified in the global array, std_depth.
   
           y:    array of y observed values
          ny:    dimension of y 
        ysum:    sum of distance-weighted yvalues on each sigma surface
      ysqsum:    sum of distance-weighted y_squared values on each sigma surface
     wghtsum:    sum of weights
       wght:     weight for this profile
       count:    counts # of observations on each surface 
*/
{
   int  i, j, datagap, n;
   double  *ytmp, *dtmp, z;
   
   ytmp = (double *) calloc(ny, sizeof(double));
   dtmp = (double *) calloc(ny, sizeof(double));
   if (dtmp == NULL) {
       fprintf(stderr,"\nUnable to allocate memory in insert_data_on_depthlevs()\n");
       exit(1);
   }

/* Ensure continuous (no missing values) y, and depth arrays */
   n = 0;
   for (i = 0; i < ny; ++i) {
         if ( y[i] > -8.9)  {
            ytmp[n] = y[i];
            dtmp[n++] = de[i];
	 }
   }
   
   if (n <= 1) {               /* not enough values */
     free((void *)ytmp);
     free((void *)dtmp);
     return;
   }

   for (i = 0; i < nstddepths; ++i) {
        z = interpolate(std_depth[i], dtmp, ytmp, n) ; 
        if (z > -998.)  {

            /* now check for datagaps around this depth...*/
            j = 0;
            while ((++j < n) &&( dtmp[j] < std_depth[i])  ) 
                ;
	    if (j == n)
	        --j;
		
            if ((dtmp[j-1] == std_depth[i]) || (dtmp[j] == std_depth[i]) )
                datagap = 0;
            else if (std_depth[i] < 1001)
                datagap = (dtmp[j] - dtmp[j-1]) > GAP_SHALLOW;
            else
                datagap = (dtmp[j] - dtmp[j-1]) > GAP_DEEP;
                  
            if (!datagap) {
                  ysum[i] += z * wght;
		  if (ysqsum != NULL)
		        ysqsum[i] += z * z * wght;
		  wghtsum[i] += wght;
                  ++count[i];
            }
      } /* end if z */
      else if (std_depth[i] < 10.0)  {  /*  get a surface value from 10 m instead */
      
          z = interpolate(10.0, dtmp, ytmp, n) ; 
          if (z > -998.)  {
                  ysum[i] += z * wght;
		  if (ysqsum != NULL) {
		        ysqsum[i] += z * z * wght;
		  }
		  wghtsum[i] += wght;
                  ++count[i];
            }
      }
   } /* end for */

   free((void *)ytmp);
   free((void *)dtmp);
   return;
} /* end insert_data_on_depthlevs() */

/*****************************************************************************/
void insert_deepest_vals(struct gridnode *g_ptr)
   /* Maintains a linked list of info regarding the deepest observed 
      level. 
      All observations within the deepest 100 meters will be retained.
      The list will actively grow and contract as each new station is
      read. The station data are 
      accessed through the global variable: struct HYDRO_DATA sta 
       
         g_ptr:   ptr to current gridnode  
      prop_indx:   contains index to properties being gridded 
         nprops:   number of properties being gridded 
   */
{
   int i, j, reflev, index;
   int density_greater, density_diff, depth_greater, depth_diff;
   double key;
   double *sptr;
   struct deepestrec  *dr_ptr, *r1, *r2;
 
   dr_ptr = g_ptr->deepest; 
   j = hdr.nobs - 1;
 
 /*  empty list case ...*/ 
     
   if (dr_ptr == NULL) {
      dr_ptr = create_deepestrec();
      
      dr_ptr->depth = sta.observ[(int)DE][j];
      dr_ptr->pressure = sta.observ[(int)PR][j];

     /* determine appropriate sigma reference level */

     reflev = 0;
      for (i = 1; i < NREFLEVS; ++i) {
           if (sta.observ[(int)PR][j] > zmin[i] && sta.observ[(int)PR][j] <= zmax[i])
               reflev = i;
      }
      switch (reflev) {
          case 0:
	       sptr = sta.observ[(int)S0];
	       break;
          case 1:
	       sptr = sta.observ[(int)S1];
	       break;
          case 2:
	       sptr = sta.observ[(int)S2];
	       break;
          case 3:
	       sptr = sta.observ[(int)S3];
	       break;
          case 4:
	       sptr = sta.observ[(int)S4];
	       break;
          default:
             fprintf(stderr,"FATAL ERROR in insert_deepest_vals():  reflev is %d", reflev);
	     exit(1); 
	  
      } /* end switch */

    /* find densest level in profile (which is not always the bottom level!) */      
      key = sptr[0];
      index = 0;
      for ( i = 1; i < hdr.nobs; ++i) {
         if (sptr[i] > key) {
            index = i;
	    key = sptr[i];
	 }
      }

      dr_ptr->sig_0 = sta.observ[(int)S0][index];
      dr_ptr->sig_1 = sta.observ[(int)S1][index];
      dr_ptr->sig_2 = sta.observ[(int)S2][index];
      dr_ptr->sig_3 = sta.observ[(int)S3][index];
      dr_ptr->sig_4 = sta.observ[(int)S4][index];
      
      dr_ptr->next = NULL;
      g_ptr->deepest = dr_ptr;
      return;
   }
 
 
   /* list is not empty.... */  
   /* determine appropriate sigma reference level */

   reflev = 0;
   for (i = 1; i < NREFLEVS; ++i) {
      if (sta.observ[(int)PR][j] > zmin[i] && sta.observ[(int)PR][j] <= zmax[i])
          reflev = i;
   }
      
       switch (reflev) {
          case 0:
	       sptr = sta.observ[(int)S0];
	       break;
          case 1:
	       sptr = sta.observ[(int)S1];
	       break;
          case 2:
	       sptr = sta.observ[(int)S2];
	       break;
          case 3:
	       sptr = sta.observ[(int)S3];
	       break;
          case 4:
	       sptr = sta.observ[(int)S4];
	       break;
          default:
             fprintf(stderr,"FATAL ERROR in insert_deepest_vals():  reflev is %d", reflev);
	     exit(1); 
	  
      } /* end switch */

    /* find densest level in profile (which is not always the bottom level!) */      
      key = sptr[0];
      index = 0;
      for ( i = 1; i < hdr.nobs; ++i) {
         if (sptr[i] > key) {
            index = i;
	    key = sptr[i];
	 }
      }
  
/* Compare to existing deepest record.... */   
   
   dr_ptr = g_ptr->deepest;
   depth_greater = sta.observ[(int)DE][j] > dr_ptr->depth ? 1 : 0;
   depth_diff =  sta.observ[(int)DE][j] - dr_ptr->depth;
   switch (reflev){
      case 0:
         density_greater = key > dr_ptr->sig_0 ? 1 : 0;
	 density_diff = ABS(key - dr_ptr->sig_0);
         break;
      case 1:
         density_greater = key > dr_ptr->sig_1 ? 1 : 0;
	 density_diff = ABS(key - dr_ptr->sig_1);
	 break;
      case 2:
         density_greater = key > dr_ptr->sig_2 ? 1 : 0;
	 density_diff = ABS(key - dr_ptr->sig_2);
         break;
      case 3:
         density_greater = key > dr_ptr->sig_3 ? 1 : 0;
	 density_diff = ABS(key - dr_ptr->sig_3);
         break;
      case 4:
         density_greater = key > dr_ptr->sig_4 ? 1 : 0;
	 density_diff = ABS(key - dr_ptr->sig_4);
         break;
      default:
         fprintf(stderr,"FATAL ERROR in insert_deepest_vals():  reflev is %d", reflev);
	 exit(1); 
   }

   
   if (density_greater && (depth_diff <= -500)) 
      return;
     
/* If station is denser than previous densest:
    or much deeper than previous densest
    insert new rec at beginning of list
    check remainder of list and delete any recs > 100 m shallower */
     
   if (density_greater  || (depth_diff > 1000) ) {
      dr_ptr = create_deepestrec();
      dr_ptr->depth = sta.observ[(int)DE][j];
      dr_ptr->pressure = sta.observ[(int)PR][j];
      
      dr_ptr->sig_0 = sta.observ[(int)S0][index];
      dr_ptr->sig_1 = sta.observ[(int)S1][index];
      dr_ptr->sig_2 = sta.observ[(int)S2][index];
      dr_ptr->sig_3 = sta.observ[(int)S3][index];
      dr_ptr->sig_4 = sta.observ[(int)S4][index];
      
/* COMMENTED OUT... for now retain depth of densest profile

     if ( !depth_greater) {
         dr_ptr->depth = g_ptr->deepest->depth;
         dr_ptr->pressure = g_ptr->deepest->pressure;
      }
*/
      
      
      dr_ptr->next = g_ptr->deepest;
      
      key = dr_ptr->depth - 100.;
      r1 = dr_ptr;
      r2 = dr_ptr->next;
      while (r2 != NULL) {
         if (r2->depth < key) {
            delete_list(r2->next);
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
   } /* end if density-greater */
   
   
 /*  it's not the densest, if it is shallower than 100 m above the
     densest, we're done.....*/

   key = g_ptr->deepest->depth -100.;

   if (sta.observ[(int)DE][j]  < key)  
      return;
      
   /* It's not denser, but if it is deeper than the densest observation so far, 
        OR It is within 100 m of densest observation, insert it into linked list*/
        
   dr_ptr = create_deepestrec();
   dr_ptr->depth = sta.observ[(int)DE][j];
   dr_ptr->pressure = sta.observ[(int)PR][j];
   
      dr_ptr->sig_0 = sta.observ[(int)S0][index];
      dr_ptr->sig_1 = sta.observ[(int)S1][index];
      dr_ptr->sig_2 = sta.observ[(int)S2][index];
      dr_ptr->sig_3 = sta.observ[(int)S3][index];
      dr_ptr->sig_4 = sta.observ[(int)S4][index];

   
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
} /* end insert_deepest_vals () */
/*****************************************************************************/
struct deepestrec * create_deepestrec()
{
   struct deepestrec *r1;
   
   r1 = (struct deepestrec *) malloc(sizeof(struct deepestrec));
   if (r1 == NULL) {
       fprintf(stderr,"\nOut of memory in create_deepestrec() \n");
       exit(1);
   }
   return(r1);
   
} /* end create_deepestrec() */
/*****************************************************************************/
void delete_list(struct deepestrec *rptr)

  /* Recursively traverses a linked list of records and frees
     up all the memory allocated to those records */
{
  /* end of list... */
   if (rptr == NULL)
      return;
      
  /* recursive part... */
   delete_list(rptr->next);
   free((void *) rptr);
   return;

}  /* end delete_list() */
/*****************************************************************************/
/*****************************************************************************/
void delete_surfrec_list(struct surfrec *rptr)

  /* Recursively traverses a linked list of records and frees
     up all the memory allocated to those records */
{
  /* end of list... */
   if (rptr == NULL)
      return;
      
  /* recursive part... */
   delete_surfrec_list(rptr->next);
   free((void *) rptr->prop);
   free((void *) rptr->prop_sqsum);
   free((void *) rptr->wghtsum);
   free((void *) rptr->count);
   free((void *) rptr);
   return;

}  /* end delete_surfrec_list() */
/*****************************************************************************/
/*****************************************************************************/
void insert_surf_vals(struct gridnode *g, int *prop_indx, int nprops, double wght)
   /* defines the depth of a mixed layer at the sea surface, determines 
      the average values of properties within it, and inserts the information
      into a linked list of records sorted by density. The station data are 
      accessed through the global variable: struct HYDRO_DATA sta 
              g:   ptr to gridnode  
      prop_indx:   contains index to properties being gridded 
         nprops:   number of properties being gridded 
	 weight:   distance weighting for this observation
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
            m->prop[i] +=  wght * x / (double) n ;
	    m->prop_sqsum[i] += x/n * x/n * wght;
	    m->wghtsum[i] +=  wght;   
	    ++(m->count[i]);        /* counts number of profiles for each property */
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
    int * tempi;
    struct surfrec *r1ptr;
    double *tempd, *tempd2, *tempd3;

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
       tempd3 = r1ptr->prop_sqsum;
       tempi = r1ptr->count;

         /* copy all fields from rptr into r1ptr */
       r1ptr->density = rptr->density;
       r1ptr->depth = rptr->depth;
       r1ptr->n = rptr->n;
       r1ptr->prop = rptr->prop;
       r1ptr->prop_sqsum = rptr->prop_sqsum;
       r1ptr->wghtsum = rptr->wghtsum;
       r1ptr->count = rptr->count;
       r1ptr->next = rptr->next;

        /* initialize the fields of rptr and link it to r1ptr */
       rptr->density = d;
       rptr->depth = 0;
       rptr->n = 0;
       rptr->prop = tempd;
       rptr->prop_sqsum = tempd3;
       rptr->wghtsum = tempd2;
       rptr->count = tempi;
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
   r->count = (int *) malloc(nprops * sizeof(int));
   r->next = NULL;
   r->wghtsum = (double *) malloc(nprops * sizeof(double));
   r->prop = (double *) malloc(nprops * sizeof(double));
   r->prop_sqsum = (double *) malloc(nprops * sizeof(double));

   if (r->prop == NULL) {
      fprintf(stderr,"\nUnable to allocate memory in create_surfrec()\n");
      exit(1);
   }

   for (i = 0; i < nprops; ++i) {
      r->prop[i] = 0.0;
      r->prop_sqsum[i] = 0.0;
      r->wghtsum[i] = 0.0;
      r->count[i] = 0;
   }   

   return (r);
}   /* end create_surfrec() */
/*****************************************************************************/
void compute_avg(double *sum, double *sqsum, double *weightsum, int *count, int nlevs)

   /*   sum:    array containing summed, weighted values 
  weightsum:    array containing sum of weights for each level
      nlevs:    # of elements in array 
   */
{
   int j;

   for (j = 0; j < nlevs; ++j) {
      sum[j] = ( weightsum[j] > 0) ?  sum[j] / weightsum[j] : (double) HBEMPTY;
      if (sqsum != NULL)
         sqsum[j] = ( weightsum[j] > 0) ?  sqsum[j] / weightsum[j] : (double) HBEMPTY;
	 
      if (count[j] > 0)	 
         weightsum[j] /= count[j];
   }
   return;

}  /* end compute_avg() */
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
void define_avg_mixed_layer(struct gridnode *section, int nx, int nprops)

   /* Visits each gridnode along the section, determines the average density at the sea surface, depth of mixed
     layer, and average properties for that mixed layer. If no surface info is available,
     appropriate values are determined from neighboring gridnodes. Replaces the linked
      list of surfrecs with a record containing averaged mix layer values. 
      
   */ 
{
   struct surfrec *r1;
   double x, *ml_depth, *ml_density;
   double   **ml_prop, **ml_prop_sqsum, **ml_wghtsum;
   int **ml_count;
   int i, index1, index2, count, n;
   int key;
 
   
   ml_depth = (double *) calloc(nx, sizeof(double));
   ml_density = (double *) calloc(nx, sizeof(double));
   ml_prop = (double **) calloc(nx, sizeof(double *));
   ml_prop_sqsum = (double **) calloc(nx, sizeof(double *));
   ml_wghtsum = (double **) calloc(nx, sizeof(double *));
   ml_count = (int **) calloc(nx, sizeof(int *));
   for (i = 0; i < nx; ++i) {
        ml_prop[i] = (double *) calloc(nprops,  sizeof(double));
        ml_prop_sqsum[i] = (double *) calloc(nprops,  sizeof(double));
        ml_wghtsum[i] = (double *) calloc(nprops,  sizeof(double));
        ml_count[i] = (int *) calloc(nprops,  sizeof(int));
   }	

   for (i = 0; i < nx; ++i) {
 
    /* find average mixed layer density */  
       r1 = section[i].mix_layer;
       x = 0.0;
       count = 0;
       while (r1 != NULL) {
           x += r1->density * r1->n ; 
          count += r1->n;   
          r1 = r1->next;    
        }
	
	if (count > 0) {
            key = NINT((x / (double)count) * 10.);     /* average density * 10  */
   
        /*  search list for the density value and replace list with a single record ... */
   
             r1 = get_density_rec(key, section[i].mix_layer, nprops);
	     
	     ml_depth[i] = r1->depth;
	     ml_density[i] = r1->density;
	     for (n = 0; n < nprops; ++n) {
	         ml_prop[i][n] = r1->prop[n];
	         ml_prop_sqsum[i][n] = r1->prop_sqsum[n];
	         ml_wghtsum[i][n] = r1->wghtsum[n];
	         ml_count[i][n] = r1->count[n];
	     }
   
             delete_surfrec_list(section[i].mix_layer );
	     r1 = create_surfrec(nprops);
	     r1->depth = ml_depth[i];
	     r1->density = ml_density[i];
	     for (n = 0; n < nprops; ++n) {
	          r1->prop[n] = ml_prop[i][n];
	         r1->prop_sqsum[n] = ml_prop_sqsum[i][n] ;
	         r1->wghtsum[n] = ml_wghtsum[i][n];
	          r1->count[n] = ml_count[i][n];
	     }
	     r1->n = r1->count[0];
	     r1->next = NULL;
	     section[i].mix_layer = r1;
	     
	 } /* end if count */
   } /* end for i */

   /* go along array and interpolate where no surface values exist*/
   
   for (i= 0; i < nx; ++i) {
      if (ml_density[i] <= 0) {
         index1 = i;
	 while ((--index1 >= 0) && (ml_density[index1] <= 0))
	      ;
 	   if (index1 < 0 )  
	      index1 = i;
	      
	    index2 = i;
	    while ((++index2 < nx) && (ml_density[index2] <= 0)) /* find neighbor */
	       ;
	    if (index2 == nx)
	       index2 = i;   
	       
	    if (index1 == i)  {              /* only one endpoint for interpolating */
                   ml_depth[i] = ml_depth[index2];
		   ml_density[i] = ml_density[index2];
		   for (n= 0; n < nprops; ++n) {
	   	      ml_prop[i][n] = ml_prop[index2][n];
	    	      ml_prop_sqsum[i][n] = ml_prop_sqsum[index2][n];
	     	      ml_wghtsum[i][n] = ml_wghtsum[index2][n];
	     	      ml_count[i][n] = ml_count[index2][n] ;
		   }
            }
	    else if (index2 == i) {
                   ml_depth[i] = ml_depth[index1];
		   ml_density[i] = ml_density[index1];
		   for (n= 0; n < nprops; ++n) {
	   	      ml_prop[i][n] = ml_prop[index1][n];
	    	      ml_prop_sqsum[i][n] = ml_prop_sqsum[index1][n];
	     	      ml_wghtsum[i][n] = ml_wghtsum[index1][n];
	     	      ml_count[i][n] = ml_count[index1][n] ;
		   }
	    }
	    else {
	        ml_depth[i] = ml_depth[index1] + (i - index1) * (ml_depth[index2] - ml_depth[index1]) / (index2 - index1);
	        ml_density[i] = ml_density[index1] + (i - index1) * (ml_density[index2] - ml_density[index1]) / (index2 - index1);
		   for (n= 0; n < nprops; ++n) {
	   	      ml_prop[i][n] = ml_prop[i][index1] + (i - index1) * (ml_prop[i][index2] - ml_prop[i][index1]) / (index2 - index1);
	    	      ml_prop_sqsum[i][n] =ml_prop_sqsum[i][index1] + (i - index1) * (ml_prop_sqsum[i][index2] - ml_prop_sqsum[i][index1]) / (index2 - index1);
	     	      ml_wghtsum[i][n] =ml_wghtsum[i][index1] + (i - index1) * (ml_wghtsum[i][index2] - ml_wghtsum[i][index1]) / (index2 - index1);
	     	      ml_count[i][n] = ml_count[index1][n];
	           }
	    }
	    
	    r1 = create_surfrec(nprops);
	    r1->density = ml_density[i];
	    r1->depth = ml_depth[i];
	    r1->next = NULL;
	    for (n = 0; n < nprops; ++n) {
	         r1->prop[n] = ml_prop[i][n];
	         r1->prop_sqsum[n] = ml_prop_sqsum[i][n];
	         r1->wghtsum[n] = ml_wghtsum[i][n];
	         r1->count[n] = ml_count[i][n];
	     }
            if (section[i].mix_layer != NULL)
	           fprintf(stderr,"in define_avg_mixed_layer()  section[i].mix_layer should be NULL but isn't\n");	     
	     section[i].mix_layer = r1;
	    
      } /* end if ml_density < 0*/
   
   }  /* end for i */
    
   for (i = 0; i < nx; ++i)  {
      free((void *)ml_prop[i]);
      free((void *)ml_prop_sqsum[i]);
      free((void *)ml_wghtsum[i]);
      free((void *)ml_count[i]);
   }
   
   free((void *)ml_depth);
   free((void *)ml_density);
   free((void *)ml_prop );
   free((void *)ml_prop_sqsum );
   free((void *)ml_wghtsum);
   free((void *)ml_count);

   return;
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
   int n;
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
      compute_avg(r1->prop, r1->prop_sqsum, r1->wghtsum, r1->count, nprops);
      r1->depth /= (double)r1->n;
      return (r1);
   }

   key2 = NINT(r2->density * 10);

   if (key > key2)                /* recursive part */
       return( get_density_rec(key, r2, nprops));
   

   if (key == key2) {             /* exact match! */
       compute_avg(r2->prop, r2->prop_sqsum, r2->wghtsum, r2->count, nprops);
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
       compute_avg(r1->prop, r1->prop_sqsum, r1->wghtsum,  r1->count, nprops);
       r1->depth /= (double)r1->n;
       return(r1);
   }

/* if we get this far, then the key must lie between r1 & r2... so 
      choose r1 or r2-- whichever is closest.   */


   if ((key - key1) < (key2 - key) ) {
       compute_avg(r1->prop, r1->prop_sqsum, r1->wghtsum,  r1->count, nprops);
       r1->depth /= (double)r1->n;
       return(r1);
   }
   else {
       compute_avg(r2->prop, r2->prop_sqsum, r2->wghtsum, r2->count, nprops);
       r2->depth /= (double)r2->n;
       return(r2);
   }


}  /* end get_density_rec() */
/***********************************************************************/
void define_bottom(struct gridnode *section, int nx)
   /* Visits each gridnode along section and determines average bottom depth and density 
     or interpolates from surrounding gridnodes, if no bottom info available */
{
   int i, index1, index2, count;
   double *bd, *bpr, **sigma;
   struct deepestrec *recptr;

/* allocate space to collect bottom info */

   bd = (double *) calloc(nx, sizeof (double));
   bpr = (double *) calloc(nx, sizeof (double));
   sigma = (double **) calloc(NREFLEVS, sizeof (double *));
   for (i = 0; i < NREFLEVS; ++i)
        sigma[i]  = (double *) calloc(nx, sizeof (double));
  
    /* find average depth and density at each gridnode and replace linked list with a single record */
   for ( i = 0; i < nx; ++i) {
   
      count = 0;
      recptr = section[i].deepest;
      while (recptr != NULL) {
          bd[i] += recptr->depth;
	  bpr[i] += recptr->pressure;
	  sigma[0][i] += recptr->sig_0;
	  sigma[1][i] +=  recptr->sig_1;
	  sigma[2][i] +=  recptr->sig_2;
	  sigma[3][i] +=  recptr->sig_3;
	  sigma[4][i] +=  recptr->sig_4;
	  recptr = recptr->next;
	  ++count;
      }
      
      if (count > 0) {
         bd[i] /= count;
	 bpr[i] /= count;
	  sigma[0][i] /= count;
	  sigma[1][i] /= count;
	  sigma[2][i] /= count;
	  sigma[3][i] /= count;
	  sigma[4][i] /= count;
          delete_list(section[i].deepest);
          recptr = create_deepestrec();
	     recptr->depth = bd[i];
	     recptr->pressure = bpr[i];
	     recptr->sig_0 = sigma[0][i];
	     recptr->sig_1 = sigma[1][i];
	     recptr->sig_2 = sigma[2][i];
	     recptr->sig_3 = sigma[3][i];
	     recptr->sig_4 = sigma[4][i];
	     recptr->count = count;
	     recptr->next = NULL;
          section[i].deepest = recptr;
      } /* end if */
   }  /* end for i  */
   
   
   /* go along array and interpolate where no bottom info exists */
   
   for (i= 0; i < nx; ++i) {
         
         if (bd[i] <= 0) {
	   index1 = i;
	   while ((--index1 >= 0) && (bd[index1] <= 0))  /* find neighboring bottom info */
	      ;
	   if (index1 < 0 )  
	      index1 = i;
	      
	    index2 = i;
	    while ((++index2 < nx) && (bd[index2] <= 0)) /* find neighboring bottom info */
	       ;
	    if (index2 == nx)
	       index2 = i;   
	       
	    if (index1 == i)  {              /* only one endpoint for interpolating */
	       bd[i] = bd[index2];
	       bpr[i] = bpr[index2];
	       sigma[0][i] = sigma[0][index2];
	       sigma[1][i] = sigma[1][index2];
	       sigma[2][i] = sigma[2][index2]; 
	       sigma[3][i] = sigma[3][index2];
	       sigma[4][i] = sigma[4][index2]; 
	     }
	     else if (index2 == i) {
	       bd[i] = bd[index1];
	       bpr[i] = bpr[index1];
	       sigma[0][i] = sigma[0][index1];
	       sigma[1][i] = sigma[1][index1];
	       sigma[2][i] = sigma[2][index1]; 
	       sigma[3][i] = sigma[3][index1];
	       sigma[4][i] = sigma[4][index1]; 
	     
	     }
	     else {
	       bd[i] = bd[index1] + (i - index1) * (bd[index2] - bd[index1]) / (index2 - index1);
	       bpr[i] = bpr[index1] + (i - index1) * (bpr[index2] - bpr[index1]) / (index2 - index1);
	       sigma[0][i] = sigma[0][index1] + (i - index1) * (sigma[0][index2] - sigma[0][index1]) / (index2 - index1);
	       sigma[1][i] = sigma[1][index1] + (i - index1) * (sigma[1][index2] - sigma[1][index1]) / (index2 - index1);
	       sigma[2][i] = sigma[2][index1] + (i - index1) * (sigma[2][index2] - sigma[2][index1]) / (index2 - index1);
	       sigma[3][i] = sigma[3][index1] + (i - index1) * (sigma[3][index2] - sigma[3][index1]) / (index2 - index1);
	       sigma[4][i] = sigma[4][index1] + (i - index1) * (sigma[4][index2] - sigma[4][index1]) / (index2 - index1);
	     }
	     
	     recptr = create_deepestrec();
	     recptr->depth = bd[i];
	     recptr->pressure = bpr[i];
	     recptr->sig_0 = sigma[0][i];
	     recptr->sig_1 = sigma[1][i];
	     recptr->sig_2 = sigma[2][i];
	     recptr->sig_3 = sigma[3][i];
	     recptr->sig_4 = sigma[4][i];
	     recptr->count = 1;
	     recptr->next = NULL;
	     section[i].deepest = recptr;
	      
	 } /* end if bd[i] */
   }  /* end for i */
   
   for  (i = 0; i < NREFLEVS; ++i)
        free((void *)sigma[i]);
   free((void *)bd);
   free((void *)bpr);
   free((void *)sigma);
    return;
} /* end define_bottom() */
/***********************************************************************/
void do_middle_segment(struct gridnode *section, int nx, int nprops)
   /* copies the isopycnally summed values from the central node in this segment 
      to the rest of the gridnodes. Mixed layer and deepest vals handled 
       in define_avg_mixed_layer() and define_bottom() */
     
{
   int imid, i, iprop, n;
   struct gridnode *node1, *node2;
     
   if (nx <= 1)  
       return;
       
    imid = nx / 2;
    node1 = &section[imid];
    
    for (i = 0; i < nx; ++i) {
       if (i != imid) {
          node2 = &section[i];
	  for (n = 0; n < nsiglevs; ++ n) {
               for (iprop = 0; iprop < nprops; ++iprop) {
	          node2->prop[iprop][n] = node1->prop[iprop][n];
	          node2->prop_sqsum[iprop][n] = node1->prop_sqsum[iprop][n];
	          node2->weightsum[iprop][n] = node1->weightsum[iprop][n];
	          node2->count[iprop][n] = node1->count[iprop][n];
	      }
	      node2->d[n] = node1->d[n];
	      node2->dweightsum[n] = node1->dweightsum[n];
	      node2->nobs[n] = node1->nobs[n];
	  }
	  
	  for (n = 0; n < nstddepths; ++n) {
               for (iprop = 0; iprop < nprops; ++iprop) {
	             node2->d_prop[iprop][n] = node1->d_prop[iprop][n];
	             node2->d_prop_sqsum[iprop][n] = node1->d_prop_sqsum[iprop][n];
	             node2->d_weightsum[iprop][n] = node1->d_weightsum[iprop][n];
	             node2->d_count[iprop][n] = node1->d_count[iprop][n];
               }
	   }
       }
    }
   return;
}  /* end do_middle_segment() */
/***********************************************************************/
int do_std_depth(struct gridnode *g, int npts, int *prop_id, int nprops, int ix)

/*  Takes the info for a gridnode and computes the mean, stddev, and avg weight of each property 
    interpolated onto standard pressure levels.  std_depth and nstddepths are defined globally.
   
    The function returns the # of std levels output  . 
    
            g:   address of this gridnode  
         npts:   size of arrays for this gridnode  
      prop_id:   array of property identifiers  
       nprops:   # of properties in prop array  
           ix:   index of this point in xdist and xtopo arrays

*/

{
   int i, ii, j, jj, k, i_te, i_sa, i_pr, datagap, room, is_denser;
   double z, r, stddev, bottom_val;
   double tref, pref, sig_bml, sig_deepest;
   double  *d, *var, *wght;
   double *x, *y, *v, *w, *c, *sigval, *sig; 
   int  *n;
   int size, start;
   int deepestlev, reflev;
   int **cptr, **d_cptr;
   double **dptr, **vptr, **wptr;
   double **d_dptr,  **d_vptr,  **d_wptr;
   struct surfrec *m, *m2;
   struct deepestrec *btmptr;


   sigval = (double *) calloc(npts,  sizeof(double));
      
   size = monotonic(g->d, g->nobs, sigval, g->prop, g->prop_sqsum, g->weightsum, g->count, nprops, npts);
   if (size <= 1) {
      free((void *)sigval);
      return(0);
   }

/* construct temporary arrays, x, y, v,  w and c to store depth, y-property, variance, weight and 
   # of obs for that property.  Enlarge the size of the arrays 
   to accommodate the mixed layer and deepest observations.  */

   room = 50;
   x = (double *) malloc((size+room) * sizeof(double));
   y = (double *) malloc((size+room) * sizeof(double));
   v = (double *) malloc((size+room) * sizeof(double));
   w = (double *) malloc((size+room) * sizeof(double));
   c = (double *) malloc((size+room) * sizeof(double));
   if (c == NULL) {
     fprintf(stderr,"\nUnable to allocate memory in do_std_depth()\n");
     exit(1);
   }

   sig = (double *) malloc(npts * sizeof(double));
   cptr = (int **) malloc(nprops *  sizeof(int *) );
   vptr = (double **) malloc(nprops * sizeof(double *) );
   wptr = (double **) malloc(nprops * sizeof(double *) );
   dptr = (double **) malloc(nprops * sizeof(double *) );

   for (i = 0; i < nprops; ++i) {
      dptr[i] = (double *) calloc(nstddepths, sizeof(double));
      vptr[i] = (double *) calloc(nstddepths, sizeof(double));
      wptr[i] = (double *) calloc(nstddepths, sizeof(double));
      cptr[i] = (int *) calloc(nstddepths, sizeof(int));
   }
	   
/* identify index for pr, temp and salt */

   i = i_te = i_sa = i_pr = -1;
   while (++i < nprops) {
      if (prop_id[i] == (int)PR)
          i_pr = i;
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
   
/* define y-value corresponding to seafloor */   
   btmptr = g->deepest;
   if (xtopo != NULL)
       bottom_val = xtopo[ix];
    else {
       bottom_val = btmptr->pressure;
       if (dflag)
           bottom_val = btmptr->depth;
    }
      
/* Define the starting point in the sigma series which is heavier than
    the density at the bottom of the average mixed layer. Ensure that 
    the depth of this sigma level is deeper than the depth of the mixed 
    layer. */
  
   m = g->mix_layer;  
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
	   
      if (dflag) {
          while ((m->depth > g->d[start]) && (start < size))
           ++start;
      }
      else {
         while ((m->depth > g->prop[i_pr][start]) && (start < size))
           ++start;
      }
   }

/* Eliminate missing values from each y-property. Interpolate the property
   onto std y-levels .  Check for vertical datagaps and flag a standard
   level as missing if the gap is too large.  Also interpolate to approximate 
   the # of observations and variance at each std level.  The interpolation routine returns 
   HBEMPTY when appropriate, so that each level is assigned a value
   whether or not any data is present.  The value of count at a level with 
   no data is assigned zero. */


   for (j = 0; j < nprops; ++j) {
         npts = 0;
         if (m != NULL) {
            if (m->wghtsum[j] > 0) {   
              npts = mixed_layer_vals(j, prop_id[j], m, x, y, v, w, c, i_pr, i_te, i_sa);
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
	       if (dflag)
                  x[npts] = g->d[k];
	       else
                  x[npts] = g->prop[i_pr][k];
		  
               y[npts] = g->prop[j][k];
	       v[npts] = g->prop_sqsum[j][k];
	       w[npts] = g->weightsum[j][k];
               c[npts] = (double) g->count[j][k];
               sig[npts] = sigval[k];
               ++npts;
            }  
         } 
	 
         d = dptr[j];
	 var= vptr[j];
	 wght = wptr[j];
         n = cptr[j];
         
        for (i = 0; i < nstddepths; ++i) {   /* interpolate onto std levels */
	
           if (npts <= 1) {
               z = (double) HBEMPTY;
               k = 0;
           }
	   else if (std_depth[i] > bottom_val) { /* skip surfaces deeper than seafloor */
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

              else if (std_depth[i] < 1001) 
                  datagap = (x[jj] - x[jj-1]) > GAP_SHALLOW;

              else
                  datagap = (x[jj] - x[jj-1]) > GAP_DEEP;

              if (datagap && (ABS(sig[jj] - sig[jj-1]) >  0.04)) {      /* check for pycnostad */
                  z = (double) HB_MISSING;
		  k = 0;
              }
              else {
                   r = interpolate(std_depth[i], x, c, npts);
                   if (r < 0)
                       k = 0;  /* this should never happen */
                   else {
                       k = (r - (int)r) > 0 ? (int) r + 1 : (int) r;
		       if (k == 0)
		          ++k;
                   }
              }
	     
              *(d++) =  z;
              *(n++) =  k;
	      if (k > 0) {
	          *(var++) = interpolate(std_depth[i], x, v, npts);
	          *(wght++) = interpolate(std_depth[i], x, w, npts);
	       }
	       else {
	         *(var++) =  HB_MISSING;
	         *(wght++) =  HB_MISSING;
	       }
	   } 
	   else {
	         *(d++) =  HB_MISSING;
	         *(n++) =  0;
	         *(var++) =  HB_MISSING;
	         *(wght++) =  HB_MISSING;
	   }  /* end if z */
        } /* end for i */
    } /* end for j */ 
    
   free((void *)sigval);
   free((void *)sig);
   free ((void *)x);
   free ((void *)y);
   free ((void *)v);
   free ((void *)w);
   free ((void *)c);

    
    /* now interpolate depth-averaged data onto standard  levels */
 
   size = nstddepths + 50;
   x = (double *) calloc((size), sizeof(double));
   y = (double *) calloc((size), sizeof(double));
   v = (double *) calloc((size), sizeof(double));
   w = (double *) calloc((size), sizeof(double));
   c = (double *) calloc((size), sizeof(double));
   if (c == NULL) {
     fprintf(stderr,"\nUnable to allocate memory in do_std_depth()\n");
     exit(1);
   }
    
         d_dptr = (double **) calloc(nprops, sizeof(double *));
         d_vptr = (double **) calloc(nprops, sizeof(double *));
         d_wptr = (double **) calloc(nprops, sizeof(double *));
         d_cptr = (int **) calloc(nprops, sizeof(int *));
   
         for (i = 0; i < nprops; ++i) {
            d_dptr[i] = (double *) calloc(nstddepths, sizeof(double));
            d_vptr[i] = (double *) calloc(nstddepths, sizeof(double));
            d_wptr[i] = (double *) calloc(nstddepths, sizeof(double));
            d_cptr[i] = (int *) calloc(nstddepths, sizeof(int));
         }
 
    for (j = 0; j < nprops; ++j) {
         npts = 0;
         for (k = 0; k < nstddepths; ++k) {
            if (g->d_count[j][k] > 0) {
	       if (dflag) 
                   x[npts] = std_depth[k];
	       else
                    x[npts] = g->d_prop[i_pr][k];
               y[npts] = g->d_prop[j][k];
	       v[npts] = g->d_prop_sqsum[j][k];
	       w[npts] = g->d_weightsum[j][k];
               c[npts] = (double) g->d_count[j][k];
               ++npts;
            }  
         } 
	 
         d = d_dptr[j];
	 var= d_vptr[j];
	 wght = d_wptr[j];
         n = d_cptr[j];
         
        for (i = 0; i < nstddepths; ++i) {   /* interpolate onto std levels */
           if (npts <= 1) {
               z = (double) HBEMPTY;
               k = 0;
           }
	   else if (std_depth[i] > bottom_val) { /* skip surfaces deeper than seafloor */
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

              else if (std_depth[i] < 1001) 
                  datagap = (x[jj] - x[jj-1]) > GAP_SHALLOW;

              else
                  datagap = (x[jj] - x[jj-1]) > GAP_DEEP;

              if (datagap) {     
                  z = (double) HB_MISSING;
		  k = 0;
              }
              else {
                   r = interpolate(std_depth[i], x, c, npts);
                   if (r < 0)
                       k = 0;  /* this should never happen */
                   else {
                       k = (r - (int)r) > 0 ? (int) r + 1 : (int) r;
                   }
              }
	      
             *(d++) =  z;
             *(n++) =  k;
	      if (k > 0) {
	         *(var++) = interpolate(std_depth[i], x, v, npts);
	         *(wght++) = interpolate(std_depth[i], x, w, npts);
	       }
	       else {
	       	   *(var++) = HB_MISSING;
	           *(wght++) = HB_MISSING;

	       }
	    } /* end if z */
	    else {
	       	   *(d++) = HB_MISSING;
	           *(n++) = 0;
	       	   *(var++) = HB_MISSING;
	           *(wght++) = HB_MISSING;
	    }  /* end else */
         }  /* end for i */
	 	 
    } /* end for j */ 

 
   /* write averages to stdout */ 

    k = 0;    
    for  (i = 0; i < nstddepths; ++i) {
    
	 /* check for pycnostads */
        if (cptr[i_pr][i] == 0 && d_cptr[i_pr][i] > 0 )  {
	    
	    for (j = 0; j < nprops; ++j) {
	         dptr[j][i] = d_dptr[j][i];
	         vptr[j][i] = d_vptr[j][i];
	         wptr[j][i] = d_wptr[j][i];
	         cptr[j][i] = d_cptr[j][i];
	    }
	}

        if (cptr[i_pr][i] > 0) {
	
	      compute_gamma_n(&ginfo, 1, x, &dptr[i_pr][i], &dptr[i_te][i], &dptr[i_sa][i], (double) -69.0, (double) 38.0);
		    
              fprintf(stdout, "%.1lf  %.1lf  %.3lf  ", xdist[ix], std_depth[i], x[0]);
	
	   for (j = 0; j < nprops; ++j) {
	      if (cptr[j][i] <= 0) {
	         dptr[j][i] = HB_MISSING;
	      }
	      stddev = HB_MISSING;
	      if (cptr[j][i] > 1)
	            stddev = sqrt (ABS(vptr[j][i] - dptr[j][i] * dptr[j][i]) * cptr[j][i] / (cptr[j][i] - 1));
		    
	      fprintf(stdout, " %.3lf %.4lf %.4lf %6d  ", dptr[j][i], stddev, wptr[j][i], cptr[j][i]);
	   }
	   
           if (isobaric_avg) {	   
	      if (d_cptr[i_pr][i] <= 0) {
                    fprintf(stdout, "   %.3lf  ", HB_MISSING);
	      }
	      else {
	           compute_gamma_n(&ginfo, 1, x, &d_dptr[i_pr][i], &d_dptr[i_te][i], &d_dptr[i_sa][i], (double) -69.0, (double) 38.0);
                    fprintf(stdout, "   %.3lf  ", x[0]);
 	      }
	      
	      for (j = 0; j < nprops; ++j) {

	         if (d_cptr[j][i] <= 0) {
	            d_dptr[j][i] = HB_MISSING;
	         }
	         stddev = HB_MISSING;
	         if (d_cptr[j][i] > 1)
	            stddev = sqrt ( ABS(d_vptr[j][i] - d_dptr[j][i] * d_dptr[j][i]) * d_cptr[j][i] / (d_cptr[j][i] - 1));
		    
	         fprintf(stdout, " %.3lf %.4lf %.4lf %3d  ", d_dptr[j][i], stddev, d_wptr[j][i], d_cptr[j][i]);
	      } /* end for j */
	   } /* end if isobaric_avg */
           fprintf(stdout, "\n");
           ++k;
       }  /* end if cptr[ipr][i] > 0 */
    } /* end for i */
    
/* clean up space no longer needed ... */

   free ((void *)x);
   free ((void *)y);
   free ((void *)v);
   free ((void *)w);
   free ((void *)c);

    for (i = 0; i < nprops; ++i) {
            free(d_dptr[i] );
            free(d_vptr[i] );
            free(d_wptr[i]);
            free(d_cptr[i]);
            free(dptr[i] );
            free(vptr[i] );
            free(wptr[i]);
            free(cptr[i]);
     }
     free(d_dptr ) ;
     free(d_vptr );
     free(d_wptr );
     free(d_cptr );
     free(dptr ) ;
     free(vptr );
     free(wptr );
     free(cptr );
    return(k);

} /* end do_std_depth() */
/****************************************************************************/
int mixed_layer_vals(int iprop, int prop_id, struct surfrec *m, double *x, double *y, double *v, double *w, double *c, int i_pr, int i_te, int i_sa)

   /* Creates arrays of depth (= pressure), property, variance, avg weights and counts at
      100 m intervals for the requested property.  Returns
      the number of points in each array. 
      
         iprop:          indicates the ith property of nprops 
       prop_id:       indicates the ith property of MAXPROPS 
             m:      info defining the mixed layer 
             x:       depth array 
             y:       property array 
             v:       variance array 
             w:       avg weight array 
             c:       count array 
    i_pr, i_te, i_sa:       index to pressure, temperature, salt (or -1) 
 
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
	 c[i] = (double) m->count[iprop];
	 v[i] = m->prop_sqsum[iprop];          
   }         
   w[npts-1] = m->wghtsum[iprop];
   x[npts-1] = m->depth;
   c[npts-1] = (double) m->n;
   v[npts-1] = m->prop_sqsum[iprop];          

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
int monotonic(double *d, int  *nobs, double *sigmas, double **xx, double **xxvar, double **xxwgt, int  **count, int nprops, int npts)

/* sorts the arrays into increasing order by depth, removing levels with no obs and 
   ensuring there are no depth inversions.  Returns the 
   # of depth levels containing data. 
   
           d    starting addr of depth  
        nobs    starting addr of nobs array  
      sigmas    array to put sigma values 
          xx    starting addr of other properties  
       xxvar    starting addr of variance arrays
       xxwgt    starting addr of average weight arrays
       count    starting addr of nobs arrays for other properties  
      nprops    # of rows (other properties) in xx  
        npts    # of cols (points) in each property array 

   */
     
{
   int i, j, k, mono;
   int size;
   double  *diff;


   diff = (double *) calloc(npts, sizeof(double));
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
	   xxvar[j][k] = xxvar[j][i];
	   xxwgt[j][k] = xxwgt[j][i];
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
	     xxvar[j][k] = xxvar[j][i];
	     xxwgt[j][k] = xxwgt[j][i];
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
                   xxvar[j][k] = xxvar[j][i];
                   xxwgt[j][k] = xxwgt[j][i];
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
} /* end monotonic() */ 

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/* Performs a linear interpolation to find the position of xval in array, x,
   and returns the corresponding value in the array, y.  If xval does not
   appear in array x, the value HBEMPTY is returned.  This routine assumes 
   that the x array is monotonic and continuous (no missing values); and it
   assumes the y array is continuous.    */

double interpolate(double xval, double *x, double *y, int nypts)
{
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

      return ( y[k] + (y[k+1] - y[k]) * v1 / (x[k+1] - x[k]) );
   }

   return ((double)HBEMPTY);

}   /* end interpolate() */

/****************************************************************************/
/* Performs a linear interpolation to find the position of xval in array, x,
   and returns the corresponding value in the array, y.  If xval does not
   appear in array x, the value HBEMPTY is returned.  This routine assumes 
   that the x array is monotonic and continuous (no missing values); and it
   assumes the y array is continuous.    */

int interpolate_count(double xval, double *x, int *y, int nypts)
{
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

      return ( NINT (y[k] + (y[k+1] - y[k]) * v1 / (x[k+1] - x[k]) ) );
   }

   return (0);

}   /* end interpolate() */
