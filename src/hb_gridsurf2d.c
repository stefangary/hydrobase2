/*  hb_gridsurf2d.c

_................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Institution
                             1993
                             Updated to ANSI standards Feb 2000
................................................................................
____________________________________________________________________________
hb_gridsurf2d projects hydrographic properties onto one or more surfaces
and computes mean and standard deviation values of those properties
for each node in a grid. The grid spacing and bounds are specified by 
the user. All points within a square are used to compute the mean and 
are weighted equally. Squares do not overlap.  The centers of each 
square form the x-y loci of elements of the matrix representing each
surface.

The surfaces may be defined by some value of any property supported
by HydroBase (depth, pressure, temperature, density, etc.)
The properties at each station are linearly interpolated onto the
surface; and the averaging, which produces the grid, is done on that
surface.  That is, if the surface is a 100 m depth surface, the
interpolated properties will be averaged on that depth surface.  
IT IS STRONGLY RECOMMENDED that properties be averaged on density surfaces
to avoid problems of mixing water masses, especially in frontal zones.
If the surface being gridded with gridsurf is of a type other than
density, we recommend that the data be gridded (averaged) first with
the program hb_grid3d, and use that product as input to hb_gridsurf.

For each surface, an output file is generated containing :
     
     nrows, ncols, grid spacing, grid bounds         (1st line)
     list_of_properties                              (2nd line)
     lat lon   n1 p1 p1dev  n2 p2 p2dev ...  ( 1 line for each gridpt)
      where for each property:      n = # of obs
                                    p = mean 
                                 pdev = std dev    
                                                      
____________________________________________________________________________
___________________________________________________________________________
  USAGE:  

 hb_gridsurf filename(_roots) [-C] -B/west/east/south/north -I<gridspacing>  -S<surface_def_file> -P<property_list> [-D<dirname>] [-E<file_extent>] [-Z]

  list of filenames MUST be first argument or input is expected to come
  from stdin.

 -C : input files are cdf format (default is HydroBase format);

 -B : specifies grid bounds

 -I : specifies grid increment in degrees;  ex: -I0.5

 -S : file containing surface definitions. If this is not specified
      these values are expected to come from stdin .

 -P : list of properties to evaluate; ex: -Ppr/te/th/sa/ox
      -P (by itself) will print a list of the available properties;

 -D : specifies directory for input files (default is current dir)  
      ex: -D/home/pete/data

 -E : specifies input_file extent (default is no extent) ex: -E.sum

 -Z : use deepest occurance of surface (for cases where property defining
      surface is not monotonic).

 -K : if present, a blank file will be output.  This is a list of
      all lat/lon where there are no observations in the cells.

 -T : if present and using -C option (netcdf input), then a list of
      cells that are masked out by topography are included.  Same
      format as the blank file.  However, this allows one to
      isolate the nodes that are blank b/c of topography rather than
      missing values.  NOT ADDED!

***************************************************************/

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

/* Data structure primitive to
 * store gridded data */
struct surface {
  double  value;    /* value of property on this surface */
  double  pref;
  int   data_ind;   /* index ID of property */
  char  density_flag ;  /* set if this is a density surface */
  double *prop[MAXPROP]; /* ptrs to grid of prop values */
  double *propsq[MAXPROP];
  unsigned *count[MAXPROP];
  FILE *fptr;    /* output file */
  FILE *bfptr;   /* output blank file */
  FILE *tfptr;   /* output topo blank file */
  struct surface *next;    
};

/* Flags to indicate which properties
 * are to be output and which are
 * needed for computation */
short prop_req[MAXPROP];         /* set of props requested for output */
short prop_needed[MAXPROP];      /* set of props requested || defining a
                                    surface || needed for computation */
/* Global variables to store stations */
struct HYDRO_DATA station;
struct HYDRO_HDR hdr;


double s_pref;          /* ref lev for computing a non-standard sigma level */
double pe_pref;         /* ref lev for computing potential energy */
double ht_pref;         /* ref lev for computing dynamic height */
double ih_pref;         /* ref lev for de-integrated height */
int window, w_incr;     /* used to specify pr window for gradient properties*/
int use_deepest;        /* flag for non-monotonic properties */
float HB_f_mask;        /* flag for masked nodes in .cdf file */

struct GAMMA_NC ginfo;   /* used for neutral density */
int warnflag;            

/* Boundaries for grid */
float   xmin, xmax, ymin, ymax, delta_x, delta_y;  
int 	lon0to360;    /* =  1 for all positive longitudes,
                       * = -1 for all negative,
                       * =  0 for crossing Greenwich */
int *merid_range;     /* required by is_in_range() */
int  ncols, nrows;
int  tbin;

/* Prototypes for locally defined functions */
void print_usage(char *);
int parse_prop_list(char *);
struct surface *add_surf(FILE *, int, int, int);
void alloc_grids(struct surface *, unsigned int);
void get_hydro_data(int, struct surface *);
void insertdata(struct surface *, int, int, double, short *);
double project_prop(double *, int, double );
double project_deriv_prop(int, double, int, double);
void get_cdf_data(int, struct surface *);
void insert_counted_data(struct surface *, int, int, short **, double, short *);
void do_stats(struct surface *, int, int);

main (int argc, char **argv)
{
  short   bopt, iopt, popt, zopt, kopt, topt, copt;  /* Presence of option flags. */
  short   cdf_flag;                /* Flag for hb ascii or cdf input. */
  int     blank_flag;              /* Flag for creating blank file. */
  int     topob_flag;              /* Flag for creating a topo blank file.*/
  int     index, nprops = 0;
  int     i;
  int     curfile = 1, nfiles = 0; /* File counters. */
  int     print_msg = 1;           /* Print message in file opening. */
  int     error, prompt;
  unsigned gridsize;
  FILE    *def_file;               /* Link to surface definition file. */
  char    *extent, *dir, *st;      /* Characters for reading command line. */
  int     infile;                  /* Input file ID. */
  struct surface *list = NULL, *surf;


  /* Are there command line arguments? */
  if (argc < 2) {
    print_usage(argv[0]);
    exit(1);
  }
  
  /*  Set default values */
  dir = DIR;
  extent = EXTENT;
  cdf_flag = 0;
  blank_flag = 0;
  topob_flag = 0;
  def_file = (FILE *) stdin;
  warnflag = 0;          /* set to 1 after warning is printed */
  tbin = 0;
  prompt = 1;            /* for surface definitions */
  bopt = iopt = popt  = kopt = copt = topt = 0;
  use_deepest = 0;
  error = 0;
  s_pref = -1;
  lon0to360 = 1;    /* all positive longitudes */
  window = 100;
  w_incr = 10;
  HB_f_mask = HBMASK;
  merid_range = (int *) NULL;
  
  /* Initialize data needed flags and
   * data storage structure. */
  for (i = 0; i < MAXPROP; ++i) {
    prop_req[i] = 0;
    prop_needed[i] = 0;
    station.observ[i] = (double *) NULL;
  }
  hdr.prop_id = (int *) NULL;
   
  /* Parse command line arguments */
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
               case 'C':                   /* Check if cdf input. */
		        copt = 1;
		        cdf_flag = 1;
                        break;

               case 'D':                   /* Get input dir */
                        dir = &argv[i][2];
                        break;

               case 'E':                   /* Get file extent */
                        extent = &argv[i][2];
                        break;

               case 'B':                   /* Get grid bounds */
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
                        
                        if (xmin >= 0 && xmax < 0)
                            xmax += 360;
                        
			if (xmin < 0 ) {
			   lon0to360 = 0;
			   if (xmax <= 0)
			        lon0to360 = -1;
			}
			
                        break;

               case 'K':
		        kopt = 1;
		        blank_flag = 1;
			break;

               case 'T':
		        topt = 1;
		        topob_flag = 1;
		        break;
               case 'I':                   /* Get grid increments. */
                        iopt = 1;
                        error = (sscanf(&argv[i][2],"%f", &delta_x) == 1) ? 0 : 1;
                        delta_y = delta_x;
                        st = strchr(&argv[i][2],'/');
                        if (st != NULL) {
                          sscanf(++st,"%f", &delta_y);
                        }
                        break;

               case 'S':                   /* Get surface definition file. */
                        def_file = fopen(&argv[i][2],"r");
                        prompt = 0;              /* turn off prompt flag */
                        if (def_file == NULL) {
                           fprintf(stderr,"\nError opening %s.\n",&argv[i][2]);
                           exit(1);
                        }
                        break;

               case 'P':                   /* Get output properties. */
                        popt = 1;
                        if ( argv[i][2] == '\0') {
                               print_prop_menu();
                               exit(0);
                        }
                        nprops = parse_prop_list(&argv[i][2]);
                        if (nprops <= 0)
                             error = 1;
                        break;

               case 'W':                   /* Get window bounds. */
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

               case 'Z':                   /* Get deepest preference. */
                        use_deepest = 1;
                        break;
			
	       case 'h':
	                print_usage(argv[0]);
			exit(0);
               default:
                        error = 1;

          }    /* End switch */

          if (error ) {
             fprintf(stderr,"\nError parsing command line args.\n");
             fprintf(stderr,"     in particular: '%s'\n", argv[i]);
             exit(1);
          }

       }  /* End if */

       else  {
           ++nfiles;
       }
   }  /* End for loop parsing command line args. */


   /* Check for any mistakes on the command line. */
   if (!bopt || !iopt || !popt) {
       fprintf(stderr,"\nYou must specify bounds, properties and gridspacing! \n");
       exit(1);
   }

   if (topt && !copt) {
       fprintf(stderr,"\nYou can only use -T when using -C (reading from netCDF)! \n");
       exit(1);     
   }

   if (!nfiles) {
       if (prompt) {
          fprintf(stderr,"\nYou must specify a surface definition file with -S");
	  fprintf(stderr,"\nwhen station files are input via stdin. ");
	  exit(1);
       }
       
       if (cdf_flag) {
	 fprintf(stderr,"\nCDF files cannot be piped in via stdin");
	 fprintf(stderr,"\nYou must specify a list. ");
	 exit(1);
       }
       fprintf(stderr,"\nExpecting input from stdin ... ");
       infile = STDIN;
   }
   fprintf(stderr,"\n\n  bounds: %.3f %.3f %.3f %.3f\n", ymin,ymax,xmin,xmax );
   /* End of checking for command line errors. */


   /* Compute dimensions of matrix formed by grid.
    * Note that grids are actually stored as a
    * single array rather than a "2D" array of
    * arrays. */
   nrows = (int) (ceil((double)((ymax - ymin) / delta_y)) + .0001);
   ncols = (int) (ceil((double)((xmax - xmin) / delta_x)) + .0001);
   gridsize = nrows * ncols;

   fprintf(stderr,"\n Grid will be:  %d rows by %d cols\n",  nrows, ncols);

   /* Get info on each surface, 
    * add surface to the linked list, 
    * allocate space for computation,
    * write heading to each output file ... */
   while ( (surf = add_surf(def_file,prompt,blank_flag,topob_flag))  != NULL) {
     surf->next = list;
     list = surf;
     alloc_grids(list, gridsize);
     fprintf(list->fptr,"%d %d %.2f %.2f %.3f %.3f %.3f %.3f %2d\n", 
             nrows, ncols, delta_x, delta_y, xmin, ymin, xmax, ymax, nprops);
     for (i = 0; i < MAXPROP; ++i ) {
       if (prop_req[i]) {
	 fprintf(list->fptr, "%2s ", get_prop_mne(i));
       }
     }
     fprintf(list->fptr,"\n");
   } /* End while loop over all surfaces. Done adding surfaces. */

   /* Check for neutral density. */
   if (prop_needed[(int)GN] || prop_needed[(int)GE]) 
     gamma_nc_init((char *) GAMMA_NC_PATH, &ginfo);

   /* Loop for each input file to load input data. */
   do {
     if (cdf_flag) {      /* For netCDF file format */

       infile = cdf_open(dir, argv[curfile], extent, print_msg);
       /* Check for bad file/error and if
	* present, skip to the next file.*/
       if (infile < 0 )
	 goto NEXTFILE;
     }
     else {               /* For hydrobase file format */
       if (nfiles > 0) {
         infile = open_hydro_file(dir, argv[curfile], extent, print_msg);
	 /* Check for bad file/error and if
	  * present, skip to the next file.*/
         if (infile < 0)
	   goto NEXTFILE;
       }
     }
     
     /* Read each file completely and
      * close the file once read. */
     if (cdf_flag) {
       get_cdf_data(infile, list) ;
     }
     else {
       get_hydro_data(infile, list);
     }

   NEXTFILE:

     if (cdf_flag) 
         cdf_close(infile);
     else {
       if (nfiles > 0)
         close(infile);
     }

   } while (curfile++ < nfiles );


   /* Compute means and std deviations and write statistics to outfile ... */
   fprintf(stderr,"\n    writing stats to outfile ...\n");

   do_stats(list,blank_flag,topob_flag);

   fprintf(stderr,"\n\nEnd of %s.\n", argv[0]);
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nUsage:  %s filename_root(s)  [-C] -B/west/east/south/north -I<delta_x[/delta_y]> [-K] -P<list_of_properties> [-S<surface_def_file>] [-W<window>[/<w_incr>]] [-D<dirname>] [-E<file_extent>] [-Z]", program);

   fprintf(stderr,"\n\n  List of filenames MUST be first argument or ");
   fprintf(stderr,"\nstation data is expected from stdin.");
   fprintf(stderr,"\n -B  : specifies grid bounds");
   fprintf(stderr,"\n -I  : specifies grid increment in degrees.  ex: -I0.5");
   fprintf(stderr,"\n       If x and y increments differ, separate with");
   fprintf(stderr,"\n       a /       ex: -I0.5/1.0");
   fprintf(stderr,"\n -P  : list of properties to project onto surface;");
   fprintf(stderr,"\n       ex:  -Ppr/th/sa/ox/ht1000");
   fprintf(stderr,"\n            -P (by itself) produces a list of available properties");
   fprintf(stderr,"\n -S  : file containing surface definitions.");
   fprintf(stderr,"\n       If this is not specified these values are expected");
   fprintf(stderr,"\n       to come from the standard input device.");
   fprintf(stderr,"\n\n    OPTIONS:");
   fprintf(stderr,"\n[-C] : input files are cdf format files.");
   fprintf(stderr,"\n       (default is HydroBase format)");
   fprintf(stderr,"\n[-D] : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n       ex: -D../data/ ");
   fprintf(stderr,"\n[-E] : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n[-K] : specifies whether a blankfile should be output.");
   fprintf(stderr,"\n       The blankfile is a list of lat/lon of zero value");
   fprintf(stderr,"\n       cells, named <surf_root>.blk, for use in masking.");
   fprintf(stderr,"\n[-T] : specifies whether a topo blankfile should be output.");
   fprintf(stderr,"\n       This is a list of nodes that are blank due to");
   fprintf(stderr,"\n       topography and not just missing values.  Only works");
   fprintf(stderr,"\n       when reading from HB2 netCDF files (-C option).");
   fprintf(stderr,"\n       THIS OPTION IS NOT YET FUNCTIONAL!");
   fprintf(stderr,"\n[-W] : Specifies pressure window (db) for computing ");
   fprintf(stderr,"\n       gradient properties {bvfreq, pv}. This ");
   fprintf(stderr,"\n       constitutes the range over which observations");
   fprintf(stderr,"\n       are incorporated into the gradient computation.");
   fprintf(stderr,"\n       The window is centered around each pressure ");
   fprintf(stderr,"\n       level being evaluated.");
   fprintf(stderr,"\n       w_incr specifies how finely to subdivide the");
   fprintf(stderr,"\n       window into increments(db) to create");
   fprintf(stderr,"\n       an evenly spaced pressure series over the window.");
   fprintf(stderr,"\n          defaults: -W100/10");
   fprintf(stderr,"\n[-Z] : use deepest occurrence of surface in each profile.");
   fprintf(stderr,"\n       The default is to use the first occurrence of surface.");
   fprintf(stderr,"\n[-h]  :  help -- prints this message.");
   fprintf(stderr,"\n\n");  
   return;
}

/*****************************************************************************/
int parse_prop_list(char *st) {
  /* Parses a list of property mnemonics
   * and sets global flags accordingly.
   * Returns the number of properties.
   * An error will cause an exit. */

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
     
     /* Special cases for properties that
      * require user-specified ref values. */
     if (((enum property)index == S_ ) || ((enum property) index == HT) || ((enum property) index == IH) || ((enum property) index == PE)) {
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
	      fprintf(stderr,"\n HT reference is: %f",ht_pref);
              break;

           case IH:
	      ih_pref = ref_val;
	      fprintf(stderr,"\n IH reference is: %f",ih_pref);
	      break;

           default:
              ;
       }  /* end switch */
     }
     
     /* End of special cases */

   } while (*st == '/');

   return (nprops);
}  /* end parse_prop_list() */

/*****************************************************************************/

struct surface *add_surf(FILE *infile,int  prompt, int bflag, int tflag) {
  /*  -  Allocates storage space for a record of type <struct surface>;
   *  -  queries the user for information about the surface (via stdin);
   *  -  sets the fields of the struct accordingly;
   *  -  opens the output file;
   *  -  returns a ptr to the struct or NULL if no struct was added.
   */
  
  static n = 0;
  char id[4];
  int  index;
  struct surface *surf;
  char    fname[80];   

  /* Commented out prompt because onerous.
  if (prompt) {
    fprintf(stderr,"\nDefine each projection surface ...\n");
    fprintf(stderr,"\n    -----  Surface property options  -------\n");
    print_prop_menu();
    fprintf(stderr,"\nbot:  (bottom observation)");
    fprintf(stderr,"\n\nend:  (end list)");
    fprintf(stderr,"\n    ----------------------------------------"); 
    fprintf(stderr,"\nChoose a property type for surface #%d: ", ++n);
  }
  */
  if ( fscanf(infile, "%s", id) != 1) {
    fprintf(stderr,"\nError reading from surface definition file\n");
    exit(1);
  }
 
  if (strncmp(id, "end", 3) == 0) 
    return(NULL);          /* exit function */
  
  /* Special cases for bottom values request. */
  if (strncmp(id, "bot", 3) == 0) {
    surf = (struct surface *) calloc(1, sizeof(struct surface));
    surf->data_ind = -1;
    if (prompt) 
      fprintf(stderr,"Enter path/name of outfile for this surface: ");
    
    if (fscanf(infile, "%s", fname) != 1) {
      fprintf(stderr,"\nError reading surf.filename from surface definition file\n");
      exit(1);
    }else{
      fprintf(stderr,"\nGot path of: %s",fname);
    }
    
    if ((surf->fptr = fopen(fname, "w")) == NULL) {
      fprintf(stderr,"\nError opening %s for output.",fname);
      exit(1);
    }
    
    if ( bflag ) {
      /* Create a blankfile name and open that file too. */
      strncat(fname,".blk",4);
	if((surf->bfptr = fopen(fname, "w")) == NULL) {
	  fprintf(stderr,"\nError opening %s for output.",fname);
	  exit(1);
	}
    }

    if ( tflag ) {
      /* Create a topo blankfile name and open that file too. */
      strncat(fname,".t",2);
	if((surf->tfptr = fopen(fname, "w")) == NULL) {
	  fprintf(stderr,"\nError opening %s for output.",fname);
	  exit(1);
	}
    }

    return (surf);

  } /* End of special setup for bottom values request. */


  /* Regular request - if there is a bottom values request,
   * this function would already have exited. */
  if ((index = get_prop_indx(id)) < 0) {       /* is choice appropriate? */
    fprintf(stderr,"\n%s is not an option!\n\n", id);
    exit(1);    
  }  /* end if */ 

  surf = (struct surface *) calloc(1, sizeof(struct surface));
  
  surf->data_ind = index;

  /*  and reference pressure, if appropriate... */
  if ( index == (int)S_  ){
    if (prompt) {
      fprintf(stderr,"enter reference Pr:  ");
    }
    if (fscanf(infile, "%lf", &surf->pref) != 1) {
      fprintf(stderr,"\nError reading pref for property %s from surface definition file\n", id);
      exit(1);
    }
  }

  /* get value of surface ... */
  if (prompt) {
    fprintf(stderr,"enter %s value", get_prop_descrip(index));
    fprintf(stderr,": ");
  }
  if (fscanf(infile, "%lf", &surf->value) != 1) {
    fprintf(stderr,"\nError reading surf.value from surface definition file\n");
    exit(1);
  }else{
    fprintf(stderr,"\nGot value of %f\n",surf->value);
  }
  
  /* !**!  Special cases for individual properties... */
  switch ((enum property) index)  {
    case GN: 
      surf->pref = 0.;
      surf->density_flag = 1;
      prop_needed[index] = 1;     
      break;
		 
    case S_:
      surf->density_flag = 1;
      if (s_pref >= 0) {   /* has this already been requested? */
	if ( NINT( s_pref) != NINT(surf->pref)) {
	  fprintf(stderr,"\nSorry. You have requested the property\n");
	  fprintf(stderr,"s_ with multiple reference levels.\n");
	  fprintf(stderr,"You can only use one of those prefs.\n");
	  fprintf(stderr,"You will probably have to do multiple\n");
	  fprintf(stderr,"runs for each different pref you want\n");
	  fprintf(stderr,"to associate with s_\n\n");
	  exit(1);
	}
      }
      s_pref = surf->pref;
      prop_needed[index] = 1;     
      break;

    case S0: 
      surf->pref = 0.;
      surf->density_flag = 1;
      prop_needed[index] = 1;     
      break;

    case S1: 
      surf->pref = 1000.;
      surf->density_flag = 1;
      prop_needed[index] = 1;     
      break;

    case S2: 
      surf->pref = 2000.;
      surf->density_flag = 1;
      prop_needed[index] = 1;     
      break;
      
    case S3: 
      surf->pref = 3000;      
      surf->density_flag = 1;
      prop_needed[index] = 1;     
      break;
      
    case S4: 
      surf->pref = 4000.;
      surf->density_flag = 1;
      prop_needed[index] = 1;     
      break;

    default:
      surf->pref = 0.;       
      surf->density_flag = 0;
      prop_needed[index] = 1; 
      break;   
      
  }  /* End switch */

  /* Open output file */

  if (prompt) 
    fprintf(stderr,"Enter path/name of outfile for this surface: ");
  
  if (fscanf(infile, "%s", fname) != 1) {
    fprintf(stderr,"\nError reading surf.filename from surface definition file\n");
    exit(1);
  }else{
    fprintf(stderr,"\nGot path of: %s",fname);
  }
   
  if ((surf->fptr = fopen(fname, "w")) == NULL) {
    fprintf(stderr,"\nError opening %s for output.",fname);
    exit(1);
  }

  if ( bflag ) {
    /* Create a blankfile name and open that file too. */
    strncat(fname,".blk",4);
      if((surf->bfptr = fopen(fname, "w")) == NULL) {
	fprintf(stderr,"\nError opening %s for output.",fname);
	exit(1);
      }
  }

  if ( tflag ) {
    /* Create a topo blankfile name and open that file too. */
    strncat(fname,".t",2);
    if((surf->tfptr = fopen(fname, "w")) == NULL) {
      fprintf(stderr,"\nError opening %s for output.",fname);
      exit(1);
    }
  }

  return (surf);
  
} /* End add_surf() */

/*****************************************************************************/

void alloc_grids(struct surface *s,unsigned int n) {

  /* Allocates memory for each appropriate array in the struct * passed as 
   * an argument and initializes all values to 0.
   *
   * The second argument contains the amount of memory to allocate.
   *
   * The count array for PR is always initialized since it is the benchmark
   * which determines whether there are any stations at a gridpt.   */

  int i, j;
  
  /* Allocate memory of a much as all
   * possible properties. */
  for (i = 0; i < MAXPROP; ++i ) {
   
    s->count[i] = NULL;   /* Counter grid */
    s->prop[i] = NULL;    /* Property grid */
    s->propsq[i] = NULL;  /* Squared value grid */

    if (prop_req[i]) {
      if ((s->count[i] = (unsigned *) malloc(n * sizeof(unsigned))) == NULL) {
	fprintf(stderr,"\n\nUnable to allocate memory for s->count");
	exit(1);
      }
      if ((s->prop[i] = (double *) malloc(n * sizeof(double))) == NULL) {
	fprintf(stderr,"\n\nUnable to allocate memory for s->prop");
	exit(1);
      }
      if ((s->propsq[i] = (double *) malloc(n * sizeof(double))) == NULL) {
	fprintf(stderr,"\n\nUnable to allocate memory for s->propsq");
	exit(1);
      }

      /* Initialize array-grids to zero. */
      for (j = 0; j < n; ++j) {
	s->count[i][j] = 0;
	s->prop[i][j] = 0.;
	s->propsq[i][j] = 0.;
      }
    }
    else if ((enum property) i == PR) {
      s->count[i] = (unsigned *) malloc(n * sizeof(unsigned));
      for (j = 0; j < n; ++j) {
	s->count[i][j] = 0;
      }
    }
    
  } /* End for loop over all properties. */

  return;

}  /* End alloc_grids() */

/****************************************************************************/

void get_hydro_data(int file, struct surface *listptr) {

  /* Reads each station in a HydroBase file and adds property values
   * to the appropriate surfaces.   This module requires that the HydroBase
   * file contains a minimum of pr, de, te, sa observations.  */

  int error, i, j;
  int row, col, sq;
  int main_props_avail, ratio_done;
  double dlat;
  short derived[MAXPROP];
  struct surface  *surf;
  
  /* Initialize this so all derived
   * properties get computed in insertdata() */
  
  for (i = 0; i < MAXPROP; ++i) {
    derived[i] = 0;
    if (station.observ[i] != NULL) {
      free((void *) station.observ[i]);
      station.observ[i] = NULL; 
    }    
  }

  /* Read each station in file ... */
  while ((error = get_station(file, &hdr, &station)) == 0) {

    /* Is station within gridbounds? */
    if (hdr.lat <= ymax && hdr.lat >= ymin && is_in_range(hdr.lon, xmin, xmax, merid_range, lon0to360) ) {

      if (lon0to360 == 1) {
	if (hdr.lon < 0)
	  hdr.lon += 360.0;
      }  
      else if (lon0to360 < 0) {
	if (hdr.lon >= 0 )
	  hdr.lon -= 360.0;
      }
      else {   /* range crosses Greenwich */
	if (hdr.lon < 0)   /* first make it positive */
	  hdr.lon += 360.;
	
	if (hdr.lon > xmax ) /* make it neg if not in interval 0->xmax */
	  hdr.lon -= 360.;
      }
      
      /* The extra smidgeon puts stations
       * that fall on the border into
       * the proper grid box. */
       row = NINT( (hdr.lat + .0001 - ymin)/ delta_y - 0.5);
       col = NINT( (hdr.lon + .0001 - xmin)/ delta_x - 0.5);

       /* Check again for  out of bounds stations   */
       if ((row < 0) || (row >= nrows) || (col < 0) || (col >= ncols))  {
	 fprintf(stderr,"\n Software bug in get_hydro_data():  is_in_range() returned TRUE, but row,col out_of_bounds");
	 continue;
       }
       sq = row * ncols + col;
           
       /* Ensure that pr, de, te, and sa are available ... */
       main_props_avail = 1;
       ratio_done = 0;
       
       if (!(available(PR, &hdr) && available(DE, &hdr) && available(TE, &hdr) 
	     && available(SA, &hdr))) {
         main_props_avail = 0;
       }

       /* Compute appropriate properties at each
        * level in station.  Derivative properties
	* (pot.vorticity, buoyancy...) get
        * computed later in insertdata(). */

       /* !**! Special cases for individual properties... */
       for (i = 0; i < MAXPROP; ++i) {
	 if (prop_needed[i] && main_props_avail && !available((enum property)i, &hdr)) {
	   switch ((enum property) i) {
	     case OX:
	       if (available(O2, &hdr)) {
		 free_and_alloc(&station.observ[i], hdr.nobs);
		 for (j=0; j < hdr.nobs; ++j) {
		   station.observ[i][j] = ox_kg2l(station.observ[(int)O2][j], station.observ[(int)PR][j], station.observ[(int)TE][j],station.observ[(int)SA][j]);
		 }
               }
               break;

             case O2:
               if (available(OX, &hdr)) {
		 free_and_alloc(&station.observ[i], hdr.nobs);
		 for (j=0; j < hdr.nobs; ++j) {
		   station.observ[i][j] = ox_l2kg(station.observ[(int)OX][j], station.observ[(int)PR][j], station.observ[(int)TE][j],station.observ[(int)SA][j]);
		 }
               }
               break;

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

             case HT:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_height(hdr.nobs, station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], ht_pref, station.observ[i]);
               break;

	     case IH:
	       free_and_alloc(&station.observ[i], hdr.nobs);
               compute_htdz_over_f(hdr.nobs, station.observ[(int)DE], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], ih_pref, hdr.lat, station.observ[i]);
	       break;

             case PE:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_energy(hdr.nobs, station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], pe_pref, station.observ[i]);
               break;

             case SV:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sp_vol( hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               break;

	     case DR: /* Fall through */

             case AL:

             case BE:
	       if ( ! ratio_done) {
                  free_and_alloc(&station.observ[(int)DR], hdr.nobs);
                  free_and_alloc(&station.observ[(int)AL], hdr.nobs);
                  free_and_alloc(&station.observ[(int)BE], hdr.nobs);
                  compute_ratio( hdr.nobs, station.observ[DR], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], station.observ[(int)AL], station.observ[(int)BE]);
		  ratio_done = 1;
	       }
               break;

             case VA:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_svan( hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               break;
	       
             case VS:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sound_vel( station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], hdr.nobs);
               break;

             case PV:   /* Fall through */
             case BF:
                  prop_needed[(int)PR] = 1;
                  break;

             case GN:
	       if (!prop_needed[(int)GE]) {
		 free_and_alloc(&station.observ[i], hdr.nobs);
		 compute_gamma_n(&ginfo, hdr.nobs, station.observ[i],
				 station.observ[(int)PR],station.observ[(int)TE],
				 station.observ[(int)SA], (double)hdr.lon, (double)hdr.lat);
	       }
	       break;
	       
             case GE:
               free_and_alloc(&station.observ[(int)GE], hdr.nobs);
               free_and_alloc(&station.observ[(int)GN], hdr.nobs);
	       compute_gamma_nerr(&ginfo, hdr.nobs, station.observ[(int)GN],   
				  station.observ[(int)GE], station.observ[(int)PR],
				  station.observ[(int)TE], station.observ[(int)SA], 
				  (double)hdr.lon, (double)hdr.lat);
	       break;

             default:
               break;

	   } /* end switch */
	 } /* end if */
       }  /* end for */

       /* Traverse the linked list of surfaces & interpolate data onto each 
	  surface */

       dlat = (double) hdr.lat;
       surf = listptr;
       while (surf != NULL) {
	 insertdata(surf, sq, hdr.nobs, dlat, derived);
	 surf = surf->next;
       }  /* end while */

    } /* end if is_in_range*/
       
    /* Clean up storage */
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

void insertdata(struct surface *sptr, int sq, int nobs, double lat, short *already_deriv)

   /*   Projects property values onto a surface and adds each value to the 
        appropriate grid element.   The station data arrays are global. 
        
        arguments:
            sptr:   defines projection surface
              sq:   grid offset for this station 
            nobs:   number of observation levels in this station 
             lat:   latitude for this station 
   already_deriv:   flags derivative props which are already computed 
   */
{
   int i, datagap;
   double p0, val, last_val, *d;
   double sns[1], tns[1], pns[1], dummy;
   double dsns[1], dtns[1], dpns[1], glevels[1];

   /*****************************************************************/
   /* return BOTTOM OBSERVATION ... */
   /*****************************************************************/
   if (sptr->data_ind == -1) {     
   
        for (i = 0; i < MAXPROP; ++i ) {

          if (prop_req[i]) {   /* was this property requested? */

/* !**! Special cases for individual properties... */

             switch ((enum property) i) {

               case BF:    /* derivative properties */
               case PV:
                    if (already_deriv[i]) 
                       val = station.observ[i][nobs-1];
                    else
                       val = project_deriv_prop(i, station.observ[(int)PR][nobs-1], nobs, lat);
                    break;

               case O2:    /* optional observed properties */
               case OX:    
               case N2:    
               case N3:
               case P4:
               case SI:
               case F1:
               case F2:
               case F3:
               case HE:
               case TU:
               case VE:
               case VN:
                   if (station.observ[i] == NULL) {
                       val = -9999.;
                       break;
                    }
                    /* if available fall through to default*/

               default:    /* all other properties */ 
                    val = station.observ[i][nobs-1];
                    break;
             } /* end switch */

             if (val > -9998.) {
               sptr->prop[i][sq] += val;
               sptr->propsq[i][sq] += (val * val);
               ++(sptr->count[i][sq]);
             }
          } 
        }  /* end for */
        return;
   }
   
   /*****************************************************************/
   /* SURFACE is a neutral surface (gamma-n) */
   /*****************************************************************/
   if (sptr->data_ind == (int)GN) {
   
      if (station.observ[(int)GN] == NULL) {
         fprintf(stderr,"\nWARNING: No pr,te,sa info to compute neutral surface.");
         fprintf(stderr,"   skipping...");
         return;
      }
   
      glevels[0] = sptr->value;
      
      neutral_surfaces(station.observ[(int)SA], station.observ[(int)TE], station.observ[(int)PR], station.observ[(int)GN], nobs, glevels, 1, sns, tns, pns, dsns, dtns, dpns);
      
      if (pns[0] < 0.0) {   /* neutral surface not found */
      
         /* If it outcrops at sea surface, increment pressure and depth counters
	    and return */
	 
         if (sptr->value < station.observ[(int)GN][0] && station.observ[(int)PR][0] < 21) {
	    ++(sptr->count[(int)PR][sq]);
	    if (prop_req[(int)DE])
	       ++(sptr->count[(int)DE][sq]);
	 }
	 
	 return;
      }  /* end neutral surface not found */
      
      
      /* check for multiple xings of the neutral surface */
      
      if (dpns[0] != 0.0 && !warnflag) {
          fprintf(stderr,"\nWARNING: multiple crossings of gamma-n. Median of the surfaces will be used. Check <ns-multiples.dat> for additonal info.\n");
	  warnflag = 1;
      }
      
      p0 = pns[0];
      
 /*   check for vertical datagaps.  An unacceptable datagap is defined
       as :       > GAP_SHALLOW for surfaces in the upper 1000 m
               > GAP_DEEP for all other surfaces   */
	       
      d = station.observ[(int)PR];  /* set this ptr for easier referencing */
      
        /* after this loop, d[i-1] and d[i] will bracket the pr value
          associated with the projection surface ... */
      i = 0;
      while (d[++i] < p0 ) 
            ;
	    
	    
      if ((p0 - d[i-1] < 10) || (d[i] - p0 < 10) )
          datagap = 0;
      else if (p0 < GAP_CHANGE_DEPTH)
          datagap = (d[i] - d[i-1]) > GAP_SHALLOW;
      else
          datagap = (d[i] - d[i-1]) > GAP_DEEP;

      if (datagap)  return;
     
     
    /* !**! Special cases for individual properties... */
      
      for (i = 0; i < MAXPROP; ++i) {
      
         if (prop_req[i]) {
	 
	    switch ((enum property) i)  {
	    
	       case DE: 
	          val = hb_depth(p0,lat);
		  break;
	    
	       case TE: 
	          val = tns[0];
		  break;
	       case SA: 
	          val = sns[0];
		  break;
	       case PR: 
	          val = pns[0];
		  break;
	       case GN: 
	          val = sptr->value;
		  break;
	       
	       case TH: 
	          val = hb_theta(sns[0], tns[0], pns[0], 0.0);
		  break;
	      case DR:
	          val = hb_ratio(sns[0], tns[0], pns[0]);
		  break;
	      case AL:
	          val = hb_alpha(sns[0], tns[0], pns[0]);
		  break;
	      case BE:
	          val = hb_beta(sns[0], tns[0], pns[0]);
		  break;
	       case S0: 
	          dummy = hb_svan(sns[0], hb_theta(sns[0],tns[0],pns[0],0.0), 0.0, &val);
		  break;
		  
	       case S1: 
	          dummy = hb_svan(sns[0], hb_theta(sns[0],tns[0],pns[0],1000.0), 1000.0, &val);
		  break;
		  
	       case S2: 
	          dummy = hb_svan(sns[0], hb_theta(sns[0],tns[0],pns[0],2000.0), 2000.0, &val);
		  break;
		  
	       case S3: 
	          dummy = hb_svan(sns[0], hb_theta(sns[0],tns[0],pns[0],3000.0), 3000.0, &val);
		  break;
		  
	       case S4: 
	          dummy = hb_svan(sns[0], hb_theta(sns[0],tns[0],pns[0],4000.0), 4000.0, &val);
		  break;
		  
	       case S_: 
	          dummy = hb_svan(sns[0], hb_theta(sns[0],tns[0],pns[0],s_pref), s_pref, &val);
		  break;
		  
	       case SV: 
	          dummy = hb_svan(sns[0], tns[0], pns[0], &val);
		  val = 1.0e8 / (val + 1000.);
		  break;
		  
	       case VA: 
	          val = hb_svan(sns[0], tns[0], pns[0], &dummy);
		  break;
	       
		  
               case BF:    /* both derivative properties */
               case PV:
                   val = project_deriv_prop(i, pns[0], nobs, lat);
                   break;
		   
		   
               case OX:     
               case O2: 
               case N2:    /* optional observed properties */
               case N3:
               case P4:
               case SI:
               case F1:
               case F2:
               case F3:
               case HE:
               case TU:
               case VE:
               case VN:
               case VS:
                    if (station.observ[i] == NULL ) {
                       val = -9999.;
                       break;
                    }
                    /* if available fall through to default*/

               default:    /* all other properties interpolate onto level of neutral surface */ 
                    val = project_prop(station.observ[i], nobs, p0);
                    break;
  	    
	    }  /* end switch */
	 
            if (val > -9998.) {
               sptr->prop[i][sq] += val;
               sptr->propsq[i][sq] += (val * val);
               ++(sptr->count[i][sq]);
            }
	 }
      }
   
      return;
   }   

   /*****************************************************************/
   /* SURFACE IS NOT THE BOTTOM OR A NEUTRAL SURFACE... */
   /*****************************************************************/
   
   if (station.observ[sptr->data_ind] == NULL) {
         fprintf(stderr,"\nWARNING: No pr,te,sa info to compute %2s surface", get_prop_mne(sptr->data_ind));
         return;
   }
   
   d = station.observ[(int)PR];  /* set this ptr for easier referencing */

 /* determine if the surface exists at this station
    and check for vertical datagaps.  An unacceptable datagap is defined
    as :       > GAP_SHALLOW for surfaces in the upper 1000 m
               > GAP_DEEP for all other surfaces   */

   if ((p0 = hb_linterp(sptr->value, station.observ[sptr->data_ind], d, nobs)) > -9998.) {
   
       /* after this loop, d[i-1] and d[i] will bracket the de value
          associated with the projection surface ... */
     i = 0;
     while (d[++i] < p0 ) 
            ;

     if (use_deepest && (i < (nobs-1))) {   /* check for deeper surface */
        last_val = p0;
        while ( (p0 = hb_linterp(sptr->value, &station.observ[sptr->data_ind][i], &d[i], nobs-i)) > -9998.) {
              last_val = p0;
              while ( d[++i] < p0 ) 
                  ;
        }
        p0 = last_val;
        i = 0;
        while (d[++i] < p0 )   /* find  d[i] associated with this deepest p0 */
            ;
     }
    

     if ((p0 - d[i-1] < 10) || (d[i] - p0 < 10) )
          datagap = 0;
     else if (p0 < GAP_CHANGE_DEPTH)
          datagap = (d[i] - d[i-1]) > GAP_SHALLOW;
     else
          datagap = (d[i] - d[i-1]) > GAP_DEEP;


      /* Exclude observations which are more than 1500 m above the 
         reference pressure, if one is specified.
         This prevents alot of bogus values from being added.
         For the case when we have a pr surface (pref = 0), then
         this will always result in a negative number = 0 - 1500,
         and p0 < neg will always result in false.*/

      datagap = datagap || (p0 < sptr->pref-1500);

     if (!datagap) {
     
       --i;     /* i still points to index of surface in pressure array */
       
       /* note that use of i changes here! */

        for (i = 0; i < MAXPROP; ++i ) {  

          if (prop_req[i]) {   /* was this property requested? */

/* !**! Special cases for individual properties... */

             switch ((enum property) i) {

               case BF:    /* derivative properties */
               case PV:
                    if (already_deriv[i]) 
                       val = project_prop(station.observ[i], nobs, p0);
                    else
                       val = project_deriv_prop(i, p0, nobs, lat);
                    break;

               case OX:     
               case O2: 
               case N2:    /* optional observed properties */
               case N3:
               case P4:
               case SI:
               case F1:
               case F2:
               case F3:
               case VE:
               case VN:
               case HE:
               case TU:
                    if (station.observ[i] == NULL ) {
                       val = -9999.;
                       break;
                    }
                    /* if available fall through to default*/

               default:    /* all other properties */
                    val = project_prop(station.observ[i], nobs, p0);
                    break;
             } /* end switch */

             if (val > -9998.) {
               sptr->prop[i][sq] += val;
               sptr->propsq[i][sq] += (val * val);
               ++(sptr->count[i][sq]);
             }
          } 
          else {                          /* prop not requested */
           if ((enum property) i == PR)   /* even if pr is not being stored */
             ++(sptr->count[i][sq]);      /* the counter is always used */
          }
        }  /* end for */
      }   /* end if  !datagap */
       
   }
   else {    /* surface not at this station */

     /* if this is a density surface, check if it outcrops at the sea surface. 
        Outcropping surfaces get a zero added to the pressure and depth
        arrays (operationally their counters just get incremented). 
        First be sure that the station isn't missing too many surface
        observations ... */

     if ((sptr->density_flag) && (d[0] < 21.)) {
       if (sptr->value < station.observ[sptr->data_ind][0]) {
           ++(sptr->count[(int) PR][sq]);    
           if (prop_req[(int) DE])
             ++(sptr->count[(int)DE][sq]); 
       }   
     }

     /* For any surface, check if it is below the bottom.  If it is,
      * put a negative counter for this square. */

   }

   return;

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
double project_deriv_prop(int index, double p0, int nobs, double lat)

   /* Computes the value of a derivative property at the specified pressure 
      level from the global station data arrays. 
      
       index:   identifies property to be projected 
          p0:   pressure at surface on which property will be projected 
        nobs:   number of p,t,s observation levels at this station 
         lat:   latitude of this station 

       */
{
   double y;

/*  !**!  Special cases for individual properties... */

   switch ((enum property)index) {
      case BF:
          y = buoyancy(p0, station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], nobs, window, w_incr);
          if (y > -9998.)
             y = y * 1.0e5;
          break;
      case PV:
          y = buoyancy(p0, station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], nobs, window, w_incr);
          if (y < -9998.0)
              break;
          y = potvort((y*y), lat);
          break;
      default:
          y = -9999.0;
          break;

   }  /* end switch */

   return (y);

}  /* end project_deriv_prop() */
/****************************************************************************/
               
void get_cdf_data(int cdfid, struct surface *listptr)
{
   struct CDF_HDR cdf;
   struct surface *surf;
   int startrow, startcol;
   int endrow, endcol;
   int cdf0to360, split, ratio_done, jj;
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

/* compare bounds of file to range of surface grid to ensure overlap */


   cdf0to360 = 1;            /* cdf longitudes are all positive */
   if (cdf.xmin < 0) {
       cdf0to360 = 0;       /* cdf longitudes are mixed sign */
       if (cdf.xmax <= 0)
          cdf0to360 = -1;  /* cdf longitudes are all negative */
   }

   lon = cdf.xmin;
   while ((lon < cdf.xmax) && ! (i = is_in_range(lon, xmin, xmax, merid_range, lon0to360))) {
       lon += cdf.xincr;
   }

   if ((i == 0) || (cdf.ymin > ymax) || (cdf.ymax < ymin)) {
       fprintf(stderr, "\nNo data in cdf file is within grid bounds.\n");
       return;
   }

   /* adjust longitudes if necessary to match sign conventions */
   
   split = 0;
   if (cdf0to360 != lon0to360) {

      if (lon0to360 > 0) {   /* for surface grid lons all positive */
     
          if (cdf0to360 < 0) {  /*cdf lons all neg */
	       cdf.xmin += 360.0;
	       cdf.xmax += 360.0;
	       cdf0to360 = 1;
	  }
	  if (cdf0to360 == 0) 	{ /*cdf lons are mixed */
	       if (cdf.xmax < xmin) {
		   cdf.xmin += 360.0;
	           cdf.xmax += 360.0;
		   cdf0to360 = 1;
	       }
	       else {
		   if ( (xmax > (cdf.xmin+360)) && (xmax < (cdf.xmax+360)))
	              split = 1;
	       }
	  }       
      }
      else if (lon0to360 < 0) { /* for surface grid lons all negative */
          if (cdf0to360 > 0) {  /*cdf lons all pos */
	       cdf.xmin -= 360.0;
	       cdf.xmax -= 360.0;
	       cdf0to360 = -1;
	  }	       
	  if (cdf0to360 == 0) 	{ /*cdf lons are mixed */
	       if (cdf.xmin > xmax) {
		   cdf.xmin -= 360.0;
	           cdf.xmax -= 360.0;
	           cdf0to360 = -1;
	       }
	      else {
	       /* test for 2 split pieces of cdf file overlapping with grid */
	   	 if ((xmin > (cdf.xmin-360)) && (xmin < (cdf.xmax-360))) 
	           split = 1;
	     }
	  }       
      }
      else {  /* lon0to360 == 0*/
          if (cdf0to360 > 0) {  /* cdf lons are are pos */
	     if (cdf.xmin > xmax) {
	        cdf.xmin -= 360;
		cdf.xmax -= 360;
		cdf0to360 = -1;
	     }
	     else {
	       /* test for 2 split pieces of cdf file overlapping with grid */
	   	if ((xmin > (cdf.xmin-360)) && (xmin < (cdf.xmax-360))) 
	           split = 1;
	     }
	  }
	  else {  /* cdf lons are are neg */
	     if (cdf.xmax < xmin) {
	        cdf.xmin += 360;
		cdf.xmax += 360;
		cdf0to360 = 1;
	     }
	     else {
	       /* test for 2 split pieces of cdf file overlapping with grid */
		  if ( (xmax > (cdf.xmin+360)) && (xmax < (cdf.xmax+360)))
	              split = 1;
	     }
	  }
      }
   }
   
   
   xoffset = yoffset = 0.0;
   
   if (cdf.node_offset) {   /* for pixel gridded cdf file */
      xoffset = 0.5 * cdf.xincr;
      yoffset = 0.5 * cdf.yincr;
   }
   
for (jj = 0; jj <= split; jj++) {

   if (jj == 1) {
      if (cdf0to360 > 0) {
         cdf.xmin -= 360;
	 cdf.xmax -= 360;
	 cdf0to360 = -1;
      }
      else {
         cdf.xmin += 360;
	 cdf.xmax += 360;
	 cdf0to360 =  1;
      
      }
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
	    case IH:
            case PE:
            case BF:
            case PV:
            case SV:
            case VS:
            case VA:
            case DR:
	    case AL:
	    case BE:
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
 
           ratio_done = 0;
      
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
		     case IH:
		       compute_htdz_over_f(npts, station.observ[(int)DE], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA], ih_pref, dlat, station.observ[(int)IH]);
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
                     case DR: 
		     case AL:
		     case BE:
		         if ( ! ratio_done ) {
                            compute_ratio( npts, station.observ[(int)DR], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA], station.observ[(int)AL], station.observ[(int)BE]);
			 
                           if (cdf.counts_included) {
                               for (j = 0; j < npts; ++j) {
			          if (prop_needed[(int)DR] )
                                        count[(int) DR][j] = count[(int)TE][j];
			          if (prop_needed[(int)AL] )
                                      count[(int) AL][j] = count[(int)TE][j];
			          if (prop_needed[(int)BE] )
                                      count[(int) BE][j] = count[(int)TE][j];
                               }
                           }
			   ratio_done = 1;
			 }
                         break;

                     case VS: 
                         compute_sound_vel( station.observ[(int)VS], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA], npts);
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

              if (cdf.counts_included) 
                 insert_counted_data(surf, sq, npts, count, dlat, prop_avail);
              else 
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
      
} /* end for jj */
   return;
}  /* end get_cdf_data() */
/****************************************************************************/

void insert_counted_data(struct surface *sptr, int sq, int nlevs, short **count, double lat, short *already_deriv)

   /*   Projects property values onto each surface and adds  values to the 
        appropriate grid element. The property arrays are globally defined,
        the address of the count arrays is passed as an argument. 
    */
{
   double  p0, val, last_val, r;
   double *x, *d;
   double sns[1], tns[1], pns[1], dummy;
   double dsns[1], dtns[1], dpns[1], glevels[1];
   int    nobs, datagap;                   
   int    i, j;

/*  bottom observations were requested... */

   if (sptr->data_ind == -1) {     
        for (i = 0; i < MAXPROP; ++i ) {

          if (prop_req[i]) {   /* was this property requested? */

/* !**! Special cases for individual properties... */

             switch ((enum property) i) {

               case BF:    /* derivative properties */
               case PV:
                    if (already_deriv[i]) 
                       val = station.observ[i][nlevs-1];
                    else
                       val = project_deriv_prop(i, station.observ[(int)PR][nlevs-1], nlevs, lat);
                    break;

               case OX:    /* optional observed properties */
               case O2:    
               case N2:
               case N3:
               case P4:
               case SI:
               case F1:
               case F2:
               case F3:
               case VE:
               case VN:
               case HE:
               case TU:
                    if (station.observ[i] == NULL) {
                       val = -9999.;
                       break;
                    }
                    /* if available fall through to default*/

               default:    /* all other properties */ 
                    val = station.observ[i][nlevs-1];
                    break;
             } /* end switch */

             if (val > -9998.) {
 	       nobs = ABS(count[i][nlevs-1]);
	       if (nobs == 0)     /* in case bottom has a zero count*/
	          nobs = 1;
               sptr->prop[i][sq] += nobs * val;
               sptr->propsq[i][sq] += (nobs * val * val);
               sptr->count[i][sq] += nobs;
             }
          } 
        }  /* end for */
        return;
   }
   
/* SURFACE is a neutral surface (gamma-n) */


   if (sptr->data_ind == (int)GN) {
   
      glevels[0] = sptr->value;
      
      neutral_surfaces(station.observ[(int)SA], station.observ[(int)TE], station.observ[(int)PR], station.observ[(int)GN], nlevs, glevels, 1, sns, tns, pns, dsns, dtns, dpns);
      
      if (pns[0] < 0.0) {   /* neutral surface not found */
      
         /* If it outcrops at sea surface, increment pressure and depth counters (effectively
	    add zero to the sums) and return */
	 
         if (sptr->value < station.observ[(int)GN][0] && station.observ[(int)PR][0] < 21) {
	    ++(sptr->count[(int)PR][sq]);
	    if (prop_req[(int)DE])
	       ++(sptr->count[(int)DE][sq]);
	 }
	 
	 return;
      }  /* end neutral surface not found */
      
      
      /* check for multiple xings of the neutral surface */
      
      if (dpns[0] != 0.0 && !warnflag) {
          fprintf(stderr,"\nWARNING: multiple crossings of gamma-n. The median will be used. Check <ns-multiples.dat> for additonal info.\n");
	  warnflag = 1;
      }
      
 /* find depth associated with pressure of neutral surface */
      
      p0 = pns[0];
      x = (double *) malloc(nlevs * sizeof(double));
      
    /* !**! Special cases for individual properties... */
      
      for (i = 0; i < MAXPROP; ++i) {
      
         if (prop_req[i]) {
	 
	    switch ((enum property) i)  {
	    
	       case DE: 
	          val = hb_depth(p0, lat);
		  break;
	    
	       case TE: 
	          val = tns[0];
		  break;
	       case SA: 
	          val = sns[0];
		  break;
	       case PR: 
	          val = pns[0];
		  break;
	       case GN: 
	          val = sptr->value;
		  break;
	       
	       case TH: 
	          val = hb_theta(sns[0], tns[0], pns[0], 0.0);
		  break;
	    
	       case S0: 
	          dummy = hb_svan(sns[0], hb_theta(sns[0],tns[0],pns[0],0.0), 0.0, &val);
		  break;
		  
	       case S1: 
	          dummy = hb_svan(sns[0], hb_theta(sns[0],tns[0],pns[0],1000.0), 1000.0, &val);
		  break;
		  
	       case S2: 
	          dummy = hb_svan(sns[0], hb_theta(sns[0],tns[0],pns[0],2000.0), 2000.0, &val);
		  break;
		  
	       case S3: 
	          dummy = hb_svan(sns[0], hb_theta(sns[0],tns[0],pns[0],3000.0), 3000.0, &val);
		  break;
		  
	       case S4: 
	          dummy = hb_svan(sns[0], hb_theta(sns[0],tns[0],pns[0],4000.0), 4000.0, &val);
		  break;
		  
	       case S_: 
	          dummy = hb_svan(sns[0], hb_theta(sns[0],tns[0],pns[0],s_pref), s_pref, &val);
		  break;
		  
	       case SV: 
	          dummy = hb_svan(sns[0], tns[0], pns[0], &val);
		  val = 1.0e8 / (val + 1000.);
		  break;
		  
	       case VA: 
	          val = hb_svan(sns[0], tns[0], pns[0], &dummy);
		  break;
	       
	       case DR: 
	          val = hb_ratio(sns[0], tns[0], pns[0]);
		  break;
	       case AL: 
	          val = hb_alpha(sns[0], tns[0], pns[0]);
		  break;
	       case BE: 
	          val = hb_beta(sns[0], tns[0], pns[0]);
		  break;
		  
               case BF:    /* both derivative properties */
               case PV:
                   val = project_deriv_prop(i, pns[0], nlevs, lat);
		   for (j = 0; j < nlevs; ++j) {
		      count[i][j] = 1;
		   }
                   break;
		   
		   
               case OX:     
               case O2: 
               case N2:    /* optional observed properties */
               case N3:
               case P4:
               case SI:
               case F1:
               case F2:
               case F3:
               case VE:
               case VN:
               case HE:
               case TU:
               case VS:
                    if (station.observ[i] == NULL ) {
                       val = -9999.;
                       break;
                    }
                    /* if available fall through to default*/

               default:   /* all other properties interpolate onto level of neutral surface */ 
                    val = project_prop(station.observ[i], nlevs, p0);
                    break;
  	    
	    }  /* end switch */
	 
            if (val > -9998.) {
	       /* determine number of obs which contributed to mean value */
	       
	       for  (j = 0; j < nlevs; ++j) {
	          x[j] = (double) ABS(count[i][j]);
	       }
	       d = station.observ[(int)PR];
	       r = hb_linterp(p0, d, x, nlevs);
	       nobs = (r - (int)r) > 0 ? (int)r + 1 : (int)r;
	       
               sptr->prop[i][sq] += (nobs * val);
               sptr->propsq[i][sq] += (val * val * nobs);
               sptr->count[i][sq] += nobs;
            }
	 } /* end if prop_req[i] */
	 
	 else {
	    if ((enum property) i == PR)
	      ++(sptr->count[i][sq]);
	 }
      }
      
      free((void *)x);
   
      return;
   }   


/* surface is not the bottom or gamma_n  ... */

   d = station.observ[(int)PR];     /* set this pointer for easier references */

 /* determine if the surface exists at this station
    and check for vertical datagaps.  An unacceptable datagap is defined
    in hydrobase.h:  GAP_SHALLOW
                     GAP_DEEP   */


   if ((p0 = hb_linterp(sptr->value, station.observ[sptr->data_ind], d, nlevs)) > -9998.) {


       /* after this loop, d[i-1] and d[i] will bracket the pressure value
          associated with the projection surface ... */

     i = 0;
     while (d[++i] < p0 ) 
            ;

     if (use_deepest && (i < (nlevs-1))) {   /* check for deeper surface */
        last_val = p0;
        while ( (p0 = hb_linterp(sptr->value, &station.observ[sptr->data_ind][i], &d[i], nlevs-i)) > -9998.) {
              last_val = p0;
              while ( d[++i] < p0 ) 
                  ;
        }
        p0 = last_val;
        i = 0;
        while (d[++i] < p0 )   /* find  d[i] associated with this deepest p0 */
            ;
     }
    
     if ((p0 - d[i-1] < 10) || (d[i] - p0 < 10) )
          datagap = 0;
     else if (p0 < GAP_CHANGE_DEPTH)
          datagap = (d[i] - d[i-1]) > GAP_SHALLOW;
     else
          datagap = (d[i] - d[i-1]) > GAP_DEEP;


      /* Exclude observations which are more than 1500 m above the 
         reference pressure, if one is specified.
         This prevents alot of bogus values from being added. */

      datagap = datagap || (p0 < sptr->pref-1500);

     if (!datagap) {

       --i;     /* i still points to index of surface in depth array */
       
        x = (double *) malloc(nlevs * sizeof(double));

       /* note that use of i changes here! */
        for (i = 0; i < MAXPROP; ++i ) {

          if (prop_req[i]) {   /* was this property requested? */
             switch ((enum property) i) {
                case BF:
                case PV:
                    if (already_deriv[i]) 
                       val = project_prop(station.observ[i], nlevs, p0);
                    else {
                       val = project_deriv_prop(i, p0, nlevs, lat);
		       for (j = 0; j < nlevs; ++j)
		          count[i][j] = 1;
		    }
                    break;
                case DE:
                    val = hb_depth(p0, lat);
                    break;
                case PR:
                    val = p0;
                    break;
                default:
                    val = project_prop(station.observ[i], nlevs, p0);
                    break;
             } /* end switch */

             if (val > -9998.) {
                 /* determine number of obs which contributed to mean value */
                  for (j = 0; j < nlevs; ++j) {
                      x[j] = (double) ABS(count[i][j]);
                  }
                  r = hb_linterp(p0, d, x, nlevs);
                  nobs = (r - (int) r) > 0 ? (int) r + 1: (int) r;

                  sptr->prop[i][sq] += val * nobs;
                  sptr->propsq[i][sq] += val * val * nobs;
                  sptr->count[i][sq]+= nobs;
             }

          } 
          else {                          /* prop not requested */
           if ((enum property) i == PR)   /* even if pr is not being stored */
             ++(sptr->count[i][sq]);      /* the counter is always used */
          }
        }  /* end for */

        free((void *)x);
     } /* end if !datagap */
       
   }
   else {    /* surface not at this station */

     /* if this is a density surface, check if it outcrops at the sea surface.  
        Outcropping surfaces get a zero added to the pressure and depth arrays.
        First be sure that it isn't missing too many surface observations ... */

     if ((sptr->density_flag) && (d[0] < 21.)) {
       if (sptr->value < station.observ[sptr->data_ind][0]) {
           ++(sptr->count[(int) PR][sq]);    
           if (prop_req[(int) DE])
             ++(sptr->count[(int)DE][sq]);   
       }   
     }

   }

   return;

}  /* end insert_counted_data() */

/****************************************************************************/

void do_stats(struct surface *listptr, int bflag, int tflag) {

  struct surface *surf;
  int    sq, row, col, i, n;
  int    print_it;
  int    print_blank;
  int    print_topob;
  float  lat, lon;
  float  mean, stddev;
  double var;
  
  surf = listptr;
  while (surf != NULL) {     /* traverse linked list of surfaces */
    sq = -1;      /* initialize, then simply visit each square*/
    
    for (row = 0; row < nrows; ++row) {
      for (col = 0; col < ncols; ++col) {
	lat = (row  + .5) * delta_y + ymin ;
	lon = (col + .5) * delta_x + xmin;
	++sq;
	
	/* Reset valid location, blank, and topo
	 * blank flags for each grid node. */
	print_it = 0;
	print_blank = 0;
	print_topob = 0;
	for (i = 0; i < MAXPROP; ++i ) {
	  if ( prop_req[i] ) {
	    /* Changed != to > so that negative surface
	     * counts can mean topography interactions,
	     * zero surface counts are blank nodes, and
	     * positive surface counts are valid nodes. */
	    if ( surf->count[i][sq] > 0 ) {
	      ++print_it;
	    }
	    /* Added extra clause to check for negative
	     * counts to identify when we have hit the
	     * topography/land. Note that this clause is
	     * only activated when tflag is engaged (-T)
	     * so that if -T is not specified, then the
	     * blank file behaves as it did without -T,
	     * i.e. logging all unoccupied nodes, even
	     * the nodes affected by topography. */
	    else if ( (surf->count[i][sq] < 0) && tflag ) {
	      ++print_topob;
	    }
	    else {
	      ++print_blank;
	    }
	  } 
	}
	     
	if (print_it) {	  
	  fprintf(surf->fptr,"%7.3f %8.3f", lat, lon);
	  for (i = 0; i < MAXPROP; ++i ) {
	    if (prop_req[i]) {
	      mean = surf->prop[i][sq];
	      stddev = 0.;
	      if ( (n = ABS(surf->count[i][sq]) ) > 1) {    /*count can be negative */
		mean = (float) surf->prop[i][sq] / (float) n;
		/* Compute variance by subtracting square of mean
		 * from mean(squared value). */
		var = (surf->propsq[i][sq] - (double) (mean * mean * n)) / (double) (n-1);
		stddev = (float) sqrt (ABS(var));
	      }
	      fprintf(surf->fptr,"  %8.3f %8.3f %4d", mean, stddev, n);
	    } /* End of prop required check. */
	  } /* End for loop over all possible props on surface. */
	  fprintf(surf->fptr,"\n");
	} /* End print_it check. */
	
	if (print_blank && bflag) {
	  fprintf(surf->bfptr,"%8.3f %7.3f\n", lon, lat);
	} /* End of blank output check. */

	if (print_topob && tflag) {
	  fprintf(surf->tfptr,"%8.3f %7.3f\n", lon, lat);
	} /* End of topo blank output check. */

      } /* End for loop over columns. */
    } /* End for loop over rows. */

    /* Close main output file */
    fclose(surf->fptr);

    /* Close blank file */
    if(bflag) {
      fclose(surf->bfptr);
    }

    /* Close topo blank file*/
    if(tflag) {
      fclose(surf->tfptr);
    }

    /* Move onto the next surface (if any) */
    surf = surf->next;

  }  /* End while loop over all surfaces. */

  return;

}  /* End do_stats(). */

/**************************************************************************/
