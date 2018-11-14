/*  hb_surf.c

................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth  Curry
                             Woods Hole Oceanographic Inst
                             1993
................................................................................
 
surf reads a HydroBase station file, projects hydrographic properties 
onto one or more surfaces and outputs the values of those properties together 
with optional year, month, lat, and lon for each station in file.

The surfaces may be defined by some value of any property supported
by HydroBase (depth, pressure, temperature, density, etc.)
The properties at each station are linearly interpolated onto the
surface;

For each surface, an output file is generated containing :
     
        one or more
     {year month lat lon station_id}   p1 ... pn  ( 1 line for each station)
                                                      
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

/* Input file pathnames */
#define    EXTENT   ""
#define    DIR      ""

/* Set up data structure to store info about each surface */
struct surface {
  double  value;        /* Value of property defining this surface */
  double  pref;
  int     data_ind;     /* Index ID of property defining the surf. */
  char    density_flag; /* Set if a density surface */
  FILE   *fptr;         /* Output file */
  struct surface *next; /* Surfaces are a linked list. */   
};

/* Flags to indicate which properties
 * are to be output and which are
 * needed for computation */   
short prop_req[MAXPROP];         /* set of props requested for output */
short prop_needed[MAXPROP];      /* set of props requested || defining a
                                    surface || needed for computation */

/* Global variables to store station */
struct HYDRO_DATA station;
struct HYDRO_HDR hdr;
double s_pref;           /* ref lev for computing a non-standard sigma level */
double ht_pref;          /* ref lev for computing dynamic height */
double pe_pref;          /* ref lev for computing potential energy */
int window, w_incr;      /* used to specify pr window for gradient properties*/
struct GAMMA_NC ginfo;   /* used for neutral density */
int warnflag;            

/* Boundaries  */
float   xmin, xmax, ymin, ymax;     
int  xdateline;
int yflag, mflag, lflag, iflag, dopt;

/* Prototypes for locally defined functions */
struct surface *add_surf(FILE *, int);
void    print_usage(char *);
void    get_hydro_data(int, struct surface *);
int     parse_prop_list(char *);

main (int argc, char **argv) {

  short   bopt, popt;              /* Flag presence of command line opts. */
  int     index, nprops = 0;
  int     i;
  int     curfile = 1, nfiles = 0; /* File counters */
  int     print_msg = 1;           /* Print message while opening files */
  int     error, prompt;
  FILE    *def_file;
  char    *extent, *dir, *st;      /* Strings for reading command line */
  int     infile;                  /* File ID for hydrobase input file */
  struct surface *list = NULL, *surf; /* The surface list */

  /* Are there command line arguments? */
  if (argc < 2) {
    print_usage(argv[0]);
    exit(1);
  }
 
  /*  Set default values */
  dir = DIR;
  extent = EXTENT;
  def_file = (FILE *) stdin;
  prompt = 1;            /* for surface definitions */
  bopt  = popt = dopt = 0;
  error = 0;
  xmin = -360;
  xmax = 360;
  ymin = -90;
  ymax = 90;
  s_pref = -1;
  xdateline = 0;
  window = 100;
  w_incr = 10;
  yflag = mflag = lflag = iflag = 0;
  warnflag = 0;           /* set to 1 after warning is printed */

  /* Initialize these ... */
  for (i = 0; i < MAXPROP; ++i) {
    prop_req[i] = 0;
    prop_needed[i] = 0;
    station.observ[i] = (double *)NULL; 
  }
  hdr.prop_id = (int *) NULL;
  
  /* Are there command line arguments? */
  if (argc <= 1) {
    print_usage(argv[0]);
    exit(1);
  }
 
  /*  Parse the command line arguments */
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
               case 'D':                   /* Get input dir. */
                        dir = &argv[i][2];
                        break;

               case 'E':                    /* Get file extent. */
                        extent = &argv[i][2];
                        break;

               case 'B':                    /* Get grid bounds. */
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

               case 'P':                    /* Get output properties. */
                        popt = 1;
                        if ( argv[i][2] == '\0') {
                               print_prop_menu();
                               exit(0);
                        }
                        nprops = parse_prop_list(&argv[i][2]);
                        if (nprops <= 0)
                             error = 1;
                        break;

               case 'S':                    /* Get surface file. */
                        def_file = fopen(&argv[i][2],"r");
                        prompt = 0;              /* turn off prompt flag */
                        if (def_file == NULL) {
                           fprintf(stderr,"\nError opening %s.\n",&argv[i][2]);
                           exit(1);
                        }
                        break;

               case 'W':                    /* Get pressure window. */
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

               case 'M':                    /* Include month in listing. */
                        mflag = 1;
			if (argv[i][2] == 'd')
			    dopt = 1;
                        break;

               case 'Y':                    /* Include year in listing. */
                        yflag = 1;
                        break;

               case 'L':                    /* Include lon/lat in listing. */
                        lflag = 1;
                        break;

               case 'I':                    /* Include station number in listing. */
                        iflag = 1;
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

   if ( !popt ) {
     fprintf(stderr,"\nYou must specify properties");
     fprintf(stderr,"\nwith the -P flag.\n");
     exit(1);
   }

   if ( !(mflag || yflag || lflag || iflag) ) {
     fprintf(stderr,"\n WARNING: You have not specified");
     fprintf(stderr,"\n any output parameter from");
     fprintf(stderr,"\n {month|year|position|station_id}.");
     fprintf(stderr,"\n Only properties will be output.\n");
   }
   
   if ( !nfiles ) {
     if ( prompt ) {
       fprintf(stderr,"\nYou must specify a surface definition file with -S");
       fprintf(stderr,"\nwhen station files are input via stdin. ");
       exit(1);
     } 
     fprintf(stderr,"\nExpecting input from stdin ... ");
   }

   fprintf(stderr,"\n\n  bounds: %.3f %.3f %.3f %.3f\n", ymin, ymax, xmin, xmax );

   /* Get info on each surface,
    * add surface to the linked list, 
    * allocate space for computation,
    * write heading to each output file. */
   while ( (surf = add_surf(def_file, prompt))  != NULL) {
     surf->next = list;
     list = surf;
   }

   fclose(def_file);
   
   if (prop_needed[(int)GN] || prop_needed[(int)GE]) 
     gamma_nc_init((char *) GAMMA_NC_PATH, &ginfo);

   infile = STDIN;

   /* Loop for each input file */
   do {

     if (nfiles > 0) 
       infile = open_hydro_file(dir, argv[curfile], extent, print_msg);

     if (infile >= 0)
       get_hydro_data(infile, list);

     if (nfiles > 0) 
       close(infile);

   } while (curfile++ < nfiles );

   fprintf(stderr,"\n\nEnd of %s.\n", argv[0]);
   exit(0);

} /* End main */

/****************************************************************************/

void print_usage(char *program) {
  fprintf(stderr,"\n%s projects hydrographic properties onto one or more surfaces", program);
  fprintf(stderr,"\nand outputs the values of those properties together with");
  fprintf(stderr,"\nsome combination of year, month, lat, and lon for each");
  fprintf(stderr,"\nprofile in a HydroBase station file.");
  fprintf(stderr,"\nSurfaces may be defined by some value of any property supported");
  fprintf(stderr,"\n(e.g. depth, pressure, temperature, density, etc.)");
 
  fprintf(stderr,"\n\nUsage:  %s filename_root(s) -P<list_of_properties>",program);
  fprintf(stderr,"\n       [-B/west/east/south/north] [-S<surface_def_file>]");
  fprintf(stderr,"\n       [-W<window>[/<w_incr>]] [-D<dirname>]");
  fprintf(stderr,"\n       [-E<file_extent>] [-Y] [-M[d]] [-L] [-I]");

  fprintf(stderr,"\n\nList of filenames MUST be first argument.");
  fprintf(stderr,"\n  If no files are named, station data expected from stdin.");
  fprintf(stderr,"\n   -P  list of properties to project onto surface;");
  fprintf(stderr,"\n          ex:  -Ppr/th/sa/ox/ht");
  fprintf(stderr,"\n               -P (by itself) produces a list of available properties");
  fprintf(stderr,"\n   OPTIONS:");
  fprintf(stderr,"\n   -B  specifies grid bounds: w/e/s/n");
  fprintf(stderr,"\n   -D  specifies dirname for input files (default is current directory) ");
  fprintf(stderr,"\n            ex: -D../data/ ");
  fprintf(stderr,"\n   -E  specifies input_file extent (default is no extent)");  
  fprintf(stderr,"\n            ex: -E.dat ");
  fprintf(stderr,"\n   -I  include station ID in output listing");
  fprintf(stderr,"\n   -L  include lat/lon in output listing");
  fprintf(stderr,"\n   -M  include month  in output listing.  Append 'd' to also output day. ex: -Md");
  fprintf(stderr,"\n   -Y  include year in output listing");
  fprintf(stderr,"\n   -S  file containing surface definitions.");
  fprintf(stderr,"\n       If this is not specified these values are expected");
  fprintf(stderr,"\n       to come from the standard input device.  This file");
  fprintf(stderr,"\n       must be specified if station data are input via stdin.");
  fprintf(stderr,"\n   -W  Specifies pressure window (db) for computing ");
  fprintf(stderr,"\n       gradient properties (bvfreq, pv...) This ");
  fprintf(stderr,"\n       constitutes the range over which observations");
  fprintf(stderr,"\n       are incorporated into the gradient computation.");
  fprintf(stderr,"\n       The window is centered around each pressure ");
  fprintf(stderr,"\n       level being evaluated.");
  fprintf(stderr,"\n       w_incr specifies how finely to subdivide the");
  fprintf(stderr,"\n       window into increments(db) to create");
  fprintf(stderr,"\n       an evenly spaced pressure series over the window.");
  fprintf(stderr,"\n          defaults: -W100/10 ");
  fprintf(stderr,"\n   -h  help (prints this message)");  
  
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

      /* !**!  Special cases for properties */

      if (((enum property)index == S_ )|| ((enum property) index == HT) || ((enum property) index == PE)) {
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


   /* End of special cases */

  } while (*st == '/');

  return (nprops);

} /* End parse_prop_list() */

/*****************************************************************************/
struct surface *add_surf(FILE *infile, int prompt)
/*  -  Allocates storage space for a record of type <struct surface>;
    -  queries the user for information about the surface (via stdin);
    -  sets the fields of the struct accordingly;
    -  opens the output file;
    -  returns a ptr to the struct  or NULL if no struct was added.
*/

{
   static n = 0;
   char id[4];
   int  index;
   struct surface *surf;
   char    fname[80];   


   if (!n && prompt) {
   fprintf(stderr,"\nDefine each projection surface ...\n");

   fprintf(stderr,"\n    -----  Surface property options  -------\n");
   print_prop_menu();
   fprintf(stderr,"\nbot:  (bottom observation)");
   fprintf(stderr,"\n\nend:  (end list)");
   fprintf(stderr,"\n    ----------------------------------------"); 

   }

   if (prompt)
         fprintf(stderr,"\nChoose a property type for surface #%d: ", ++n);
   if ( fscanf(infile, "%s", id) != 1) {
         fprintf(stderr,"\nError reading from surface definition file\n");
         exit(1);
   }
 
   if (strncmp(id, "end", 3) == 0) 
               return(NULL);          /* exit function */

   if (strncmp(id, "bot", 3) == 0) {
         surf = (struct surface *) malloc(sizeof(struct surface));
         surf->data_ind = -1;
         if (prompt) 
            fprintf(stderr,"Enter path/name of outfile for this surface: ");

         if (fscanf(infile, "%s", fname) != 1) {
            fprintf(stderr,"\nError reading surf.filename from surface definition file\n");
            exit(1);
         }
   
         if ((surf->fptr = fopen(fname, "w")) == NULL) {
             fprintf(stderr,"\nError opening %s for output.");
             exit(1);
         }

         return (surf);
   }

      if ((index = get_prop_indx(id)) < 0) {       /* is choice appropriate? */
          fprintf(stderr,"\n%s is not an option!\n\n", id);
          exit(1);    
      }  /* end if */ 


       surf = (struct surface *) malloc(sizeof(struct surface));

       surf->data_ind = index;


  /*  and reference pressure, if appropriate... */

      if ( (index== (int)S_ ) || (index == (int)HT) || (index == (int)PE)){
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
         fprintf(stderr,"enter %s value", get_prop_mne(index));
         fprintf(stderr,": ");
      }
      if (fscanf(infile, "%lf", &surf->value) != 1) {
        fprintf(stderr,"\nError reading surf.value from surface definition file\n");
        exit(1);
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
                         fprintf(stderr,"Sorry. You have requested the property s_ with multiple reference levels.\n");
                         fprintf(stderr,"You can only use one of those prefs");
                         fprintf(stderr,"  You will probably have to do multiple runs for each different pref you want to associate with s_\n\n");
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

          case HT: 
                 ht_pref = surf->pref;       
                 surf->density_flag = 0;
                 prop_needed[index] = 1; 
                 break;   
	         
          case PE: 
                 pe_pref = surf->pref;       
                 surf->density_flag = 0;
                 prop_needed[index] = 1; 
                 break;   
	         
          default:
                 surf->pref = 0.;       
                 surf->density_flag = 0;
                 prop_needed[index] = 1; 
                 break;   

       }  /* end switch */



/* open output file */

   if (prompt) 
      fprintf(stderr,"Enter path/name of outfile for this surface: ");

   if (fscanf(infile, "%s", fname) != 1) {
         fprintf(stderr,"\nError reading surf.filename from surface definition file\n");
         exit(1);
   }
   
   if ((surf->fptr = fopen(fname, "w")) == NULL) {
      fprintf(stderr,"\nError opening %s for output.");
      exit(1);
   }

   return (surf);

}  /* end add_surf() */

/*****************************************************************************/

void get_hydro_data(int file, struct surface *listptr) {

  /* Reads each station in a HydroBase file
   * and adds property values to the appropriate
   *surfaces.  This module requires that the
   * HydroBase file contains a minimum of pr,
   * de, te, sa observations.  */

  int error, i, main_props_avail, ratio_done, j;
  double dlat;
  short derived[MAXPROP];
  struct surface  *surf;
  void outputdata(struct surface *, short *);

  /* Initialize this to flag all derived properties. */
  for (i = 0; i < MAXPROP; ++i) 
    derived[i] = 0;

  /* Read each station in file. */

  while ((error = get_station(file, &hdr, &station)) == 0) {

    /* Check for and skip over out of bounds stations. */
    if (xdateline && hdr.lon < 0)
      hdr.lon += 360;
    
    if ((hdr.lat > ymax) || (hdr.lat < ymin) 
        || (hdr.lon < xmin) || (hdr.lon > xmax))  
      continue;
   
    /* Ensure that pr, de, te, and sa are available. */
    main_props_avail = 1;
    if (!(available(PR, &hdr) && available(DE, &hdr) && available(TE, &hdr) 
	  && available(SA, &hdr))) {
      main_props_avail = 0;
    }
    ratio_done = 0;

    /* Compute appropriate properties at each level
     * in station.  Derivative properties (pot.vorticity,
     * buoyancy...) get computed later . */

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

  	  case PE:
	    free_and_alloc(&station.observ[i], hdr.nobs);
	    compute_energy(hdr.nobs, station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], pe_pref, station.observ[i]);
	    break;

	  case SV:
	    free_and_alloc(&station.observ[i], hdr.nobs);
	    compute_sp_vol( hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
	    break;

	  case DR: /* Fall though */
	  case AL:
	  case BE:
	    if (! ratio_done) {
	      free_and_alloc(&station.observ[(int)DR], hdr.nobs);
	      free_and_alloc(&station.observ[(int)AL], hdr.nobs);
	      free_and_alloc(&station.observ[(int)BE], hdr.nobs);
	      compute_ratio( hdr.nobs, station.observ[(int) DR], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], station.observ[(int)AL], station.observ[(int)BE]);
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
	    
	  case PV: /* Fall through */
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
	} /* End switch. */
      } /* End of prop needed check. */
    } /* End of loop over all properties. */

    /* Traverse the linked list of surfaces and
     * interpolate data onto each surface */
    surf = listptr;
    while (surf != NULL) {
      outputdata(surf, derived);
      surf = surf->next;
    }  /* end while */

    /* Once all stations have been incorporated
     * into the surfaces, clear input data. */
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

void outputdata(struct surface *sptr, short *already_deriv) {

  /* Projects property values onto a surface and
   * writes values to output file. Note that the
   * station data arrays are global. 
   *
   * Arguments:
   *             sptr: defines projection surface 
   *    already_deriv: flags derivative props which
   *                   are already computed 
   */

  int i, datagap;
  char field_string[20];
  double p0, d0, val, dlat, *d;
  double sns[1], tns[1], pns[1], dummy;
  double dsns[1], dtns[1], dpns[1], glevels[1];
  double project_deriv_prop(int, double, int, double);
  double project_prop(double *, int, struct surface *);
  double project_on_depth(double *, int, double);

  /* Return bottom observations for
   * surface property index value = -1. */
  if (sptr->data_ind == -1) {

    /* Print time info. */
    if (yflag)
      fprintf(sptr->fptr,"%5d ", hdr.year);
    if (mflag)
      fprintf(sptr->fptr,"%2d ", hdr.month);
    if (dopt)
      fprintf(sptr->fptr,"%2d ", hdr.day);
    
    /* Print position info. */
    if (lflag)
      fprintf(sptr->fptr,"%8.3f %8.3f ", hdr.lon, hdr.lat);

    /* Print station ID info. */
    if (iflag)
      fprintf(sptr->fptr,"%5d ", hdr.station);
    
    /* Loop over each possible property. */
    for (i = 0; i < MAXPROP; ++i ) {
      
      /* Check that this property was requested. */
      if (prop_req[i]) {
	
	/* !**! Special cases for individual properties... */
	switch ((enum property) i) {
	  
	  case BF:    /* derivative properties */
	  case PV:
	    if (already_deriv[i]) {
	      val = station.observ[i][hdr.nobs-1];
	    }
	    else{
	      dlat = (double) hdr.lat;
	      val = project_deriv_prop(i, station.observ[(int)PR][hdr.nobs-1], hdr.nobs, dlat);
	    }
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
	  case VE:
	  case VN:
	  case VS:
	  case HE:
	  case TU:
	    if (station.observ[i] == NULL) {
	      val = -9999.;
	      break;
	    }
	    /* if available fall through to default*/
	  default:    /* all other properties */ 
	    val = station.observ[i][hdr.nobs-1];
	    break;

	} /* End switch */

	/* Decide on output field parameters. */
	sprintf(field_string," %c%d.%dlf", '%', get_field_width(i), get_field_precis(i));

	/* Print the value to output. */
	fprintf(sptr->fptr, field_string, val);
      } /* End of property requested check. */
    }  /* End for loop over all possible properties. */

    /* Done printing any possible values from
     * this station.  Add a newline to the
     * output before proceding to the next station.*/
    fprintf(sptr->fptr, "\n");
    return;
  } /* End of bottom values check. */
   
  /* SURFACE is a neutral surface (gamma-n) */
  if (sptr->data_ind == (int)GN) {
   
    /* Check that there is sufficient info to
     * compute the local neutral surface for
     * this station. */
    if (station.observ[(int)GN] == NULL) {
      fprintf(stderr,"\nWARNING: No pr,te,sa info to compute neutral surface.");
      fprintf(stderr,"   skipping...");
      return;
    }
   
    glevels[0] = sptr->value;
      
    neutral_surfaces(station.observ[(int)SA], station.observ[(int)TE], station.observ[(int)PR], station.observ[(int)GN], hdr.nobs, glevels, 1, sns, tns, pns, dsns, dtns, dpns);
      
    if (pns[0] < 0.0) {   /* neutral surface not found */
      return;
    }  
      
    /* Check for multiple xings of the neutral surface */
    if (dpns[0] != 0.0 && !warnflag) {
      fprintf(stderr,"\nWARNING: multiple crossings of gamma-n.");
      fprintf(stderr,"\nWill use median of surfaces.");
      fprintf(stderr,"\nCheck <ns-multiples.dat> for additonal info.\n");
      warnflag = 1;
    }
      
    /* Compute the depth of neutral surface given
     * its pressure location. */
    d0 = hb_linterp(pns[0], station.observ[(int)PR], station.observ[(int)DE], hdr.nobs);
      
    /* Check for vertical datagaps. An unacceptable
     * datagap is defined as :
     *   > GAP_SHALLOW for surfaces in the upper 1000 m
     *   > GAP_DEEP for all other surfaces
     */
	       
    d = station.observ[(int)DE];  /* set this ptr for easier referencing */
      
    /* After this loop, d[i-1] and d[i] will
     * bracket the de value associated with
     * the projection surface ... */
    i = 0;
    while (d[++i] < d0 ) 
      ;
    
    if ((d0 - d[i-1] < 10) || (d[i] - d0 < 10) )
      datagap = 0;
    else if (d0 < GAP_CHANGE_DEPTH)
      datagap = (d[i] - d[i-1]) > GAP_SHALLOW;
    else
      datagap = (d[i] - d[i-1]) > GAP_DEEP;
    
    if (datagap)  return;
    
    /* If we don't have a data gap, proceed
     * to print out time, location, cruise ID
     * to output file. */
    if (yflag)
      fprintf(sptr->fptr,"%5d ", hdr.year);
    if (mflag)
      fprintf(sptr->fptr,"%2d ", hdr.month);
    if (dopt)
      fprintf(sptr->fptr,"%2d ", hdr.day);
    if (lflag)
      fprintf(sptr->fptr,"%8.3f %8.3f ", hdr.lon, hdr.lat);
    if (iflag)
      fprintf(sptr->fptr,"%5d ", hdr.station);
    
    /* Special cases for values of individual
     * properties. Loop over every possible prop. */
    for (i = 0; i < MAXPROP; ++i) {
      
      /* Check that we want this property. */
      if (prop_req[i]) {
	
	/* Assign the value of the current
	 * property to val. */
	switch ((enum property) i)  {
	  
	  case DE: 
	    val = d0;
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
	    
	  case BF:    /* Fall through: both derivative properties */
	  case PV:    /* Note that lat calc is kept with val calc.*/
	    dlat = (double) hdr.lat;
	    val = project_deriv_prop(i, pns[0], hdr.nobs, dlat);
	    break;
		   	   
	  case OX:     
	  case O2: 
	  case N2:    /* Optional observed properties */
	  case N3:    /* all fall through. */
	  case P4:
	  case SI:
	  case F1:
	  case F2:
	  case F3:
	  case VE:
	  case VN:
	  case VS:
	  case HE:
	  case TU:
	    if (station.observ[i] == NULL ) {
	      val = -9999.;
	      break;
	    }
	    /* if available fall through to default*/
	    
 	  default:    /* all other properties interpolate onto depth of neutral surface */ 
	    val = project_on_depth(station.observ[i], hdr.nobs, d0);
	    break;
  	    
	}  /* end switch */
	    
	/* Choose data output field size. */    
	sprintf(field_string," %c%d.%dlf", '%', get_field_width(i), get_field_precis(i));

	/* Write data to output file. */
	fprintf(sptr->fptr, field_string, val);
      }
    } /* End for loop over all possible properties. */
      
    /* Done printing any possible values from
     * this station.  Add a newline to the
     * output before proceding to the next station.*/
    fprintf(sptr->fptr, "\n");
    return;
  } /* End of neutral surface check. */ 

  /* Surface is not the bottom or a neutral
   * so work with a GENERALIZED SURFACE (pr,
   * density, salinity, etc.). */

  /* Set ptr to list of station obs depths
   * for easier referencing */
  d = station.observ[(int)DE];

  /* Determine if the surface exists at
   * this station and check for vertical
   * datagaps.  An unacceptable datagap
   * is defined as :
   *   > GAP_SHALLOW for surfaces in the upper 1000 m
   *   > GAP_DEEP  for all other surfaces
   *
   * The way to check if the surface exists
   * at this station is to check if there is
   * a valid depth interpolation of the surface
   * value (sptr->value) within the values of
   * surface coordinates (station.observ[sptr->data_ind]).
   * For example, if we seek a density surface,
   * then sptr->value is the density value of
   * the surface and station.observ[sptr->data_ind]
   * is the density profile.  The depth of the
   * surface is d0.
   */
  if ((d0 = hb_linterp(sptr->value, station.observ[sptr->data_ind], d, hdr.nobs)) > -9998.) {

    /* Check that pr is required for
     * derivative props and if it is
     * required, interpolate the station
     * pressure to the surface. */
    if (prop_needed[(int)PR])
      p0 = hb_linterp(sptr->value, station.observ[sptr->data_ind], station.observ[(int)PR], hdr.nobs);

    /* Find the depth indices that lie
     * above and below the surface in
     * the station.  After this loop,
     * d[i-1] and d[i] will bracket the
     * de value associated with the
     * projection surface. */
    i = 0;
    while (d[++i] < d0 ) 
      ;

    if ((d[i-1] == d0) || (d[i] == d0) )
      datagap = 0;
    else if (d0 < GAP_CHANGE_DEPTH)
      datagap = (d[i] - d[i-1]) > GAP_SHALLOW;
    else
      datagap = (d[i] - d[i-1]) > GAP_DEEP;
    
    /* Exclude observations which are
     * more than 1500 m above the
     * reference pressure, if one is
     * specified. This prevents alot
     * of bogus values from being added,
     * but may cause results to be
     * different than expected if user
     * is unaware of this. */
    datagap = datagap || (d0 < sptr->pref-1500);

    if (!datagap) {
      /* For a good value with a small gap,
       * print the coordinates to output file. */
      if (yflag)
	fprintf(sptr->fptr,"%5d ", hdr.year);
      if (mflag)
	fprintf(sptr->fptr,"%2d ", hdr.month);
      if (dopt)
	fprintf(sptr->fptr,"%2d ", hdr.day);
      if (lflag)
	fprintf(sptr->fptr,"%8.3f %8.3f ", hdr.lon, hdr.lat);
      if (iflag)
	fprintf(sptr->fptr,"%5d ", hdr.station);

      /* Print any requested values to output file. */
      for (i = 0; i < MAXPROP; ++i ) {

	/* Was this property requested? */
	if (prop_req[i]) {

	  /* Special cases for individual properties. */
	  switch ((enum property) i) {

	    case BF:    /* Both derivative properties */
	    case PV:    /* fall through. */
	      if (already_deriv[i]) {
		val = project_prop(station.observ[i], hdr.nobs, sptr);
	      }
	      else {
		dlat = (double) hdr.lat;
		val = project_deriv_prop(i, p0, hdr.nobs, dlat);
	      }
	      break;

	    case OX:    /* Optional observed properties */
	    case O2:    /* all fall through. */
	    case N2:
	    case N3:
	    case P4:
	    case SI:
	    case F1:
	    case F2:
	    case F3:
	    case VE:
	    case VN:
	    case VS:
	    case HE:
	    case TU:
	      /* If available fall through to default*/
	      if (station.observ[i] == NULL ) {
		val = -9999.;
		break;
	      }

	    default:    /* All other properties */ 
	      val = project_prop(station.observ[i], hdr.nobs, sptr);
	      break;
	  } /* End switch selecting property. */

	  /* Get info about output formatting. */
	  sprintf(field_string," %c%d.%dlf", '%', get_field_width(i), get_field_precis(i));

	  /* Print the actual data value to output file. */
	  fprintf(sptr->fptr, field_string, val);

	} /* End property requested check. */
      }  /* End for loop over all possible properties. */

      /* Done with all possible values on the surface
       * at this station.  Terminate line with newline. */
      fprintf(sptr->fptr, "\n");

    } /* End datagap check. */
  } /* Presence of surface check. */
  
  return;

}  /* End outputdata() */

/****************************************************************************/

double project_prop(double *obs, int nobs, struct surface *sptr) {

  /* Projects property values onto a
   * specified surface.  Returns the
   * value of the property at the
   * pressure level corresponding to
   * this surface, or -9999 if the
   * property cannot be projected.
   * Uses global arrays pr & de.
   *
   * arguments:
   *    obs:    array of observations for property 
   *    nobs:   number of observations 
   *    sptr:   specifies the surface on which
   *            to project the prop. Note that
   *            thoughout this routine, sptr->data_ind
   *            specifies what type of surface we
   *            are looking for.
   */

  double val;
  double x[2], y[2];
  int    j;
  int   prevdepth, curdepth;

  /* Initialize depth counters. */
  prevdepth = curdepth = 0;
  j = 0;

  /* Find first valid observation level.
   * If we overrun the profile without
   * finding the first level, quit with
   * error flag. */
  while ((obs[j] < -8.9)||(station.observ[sptr->data_ind][j] < -8.9)) {
    if (++j == nobs)
      return (-9999.);
  }

  /* Initialize the "previous" values
   * or "initial" values. */
  x[0] = station.observ[sptr->data_ind][j];
  y[0] = obs[j];

  /* Check successive pairs of datapoints
   * until the  surface is found or the
   * last observation is encountered. */
  while (++j < nobs) {
    if ((obs[j] > -8.9) && (station.observ[sptr->data_ind][j] > -8.9)) {

      /* Update the "current" or "next" values. */
      x[1] = station.observ[sptr->data_ind][j];
      y[1] = obs[j];

      /* If we perform a valid interpolation
       * between the initial and next values,
       * then return the result. Several non-
       * valid interpolations will be performed
       * in the process of searching for the
       * final value. */
      if ((val = hb_linterp(sptr->value, x, y, 2)) > -9998.) {
	return (val);
      }

      /* Update the "intial" values for next loop. */
      x[0] = x[1];
      y[0] = y[1];
    }
  }  /* End while loop over all observations. */
  
  /* If we did not exit earlier with a valid
   * interpolation, then we must quit with
   * and error flag value. */
  return (-9999.);

} /* End project_prop() */

/****************************************************************************/

double project_deriv_prop(int index, double p0, int nobs, double lat) {

  /* Computes the value of a derivative
   * property at the specified pressure 
   * level from the global station data. 
   *
   * arguments:
   *  index:   identifies property to be projected 
   *  p0:      pressure at surface on which property will be projected 
   *  nobs:    number of p,t,s observation levels at this station 
   *  lat:     latitude of this station 
   */

   double y;

   /* Special cases for individual properties. */

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
double project_on_depth(double *obs, int nobs, double d0) {

  /* Projects property values onto the depth
   * surface (d0)...Returns the value of the
   * property at the level corresponding to
   * this surface, or -9999 if the property
   * cannot be projected.  The global arrays 
   * station.observ are used.  
   *     
   *         obs:   array of observations for property
   *        nobs:   number of observations
   *          d0:   interpolated depth of surface
   */

  double val;
  double x[2], y[2];
  int    j;

  j = 0;

  /* Find first valid observation level.
   * If we overrun the profile without
   * finding the first level, quit with
   * error flag. */
  while (obs[j] < -8.9) {
    if (++j == nobs)
      return (-9999.);
  }

  /* Initialize the "previous" values
   * or "initial" values. */
  x[0] = station.observ[(int)DE][j];
  y[0] = obs[j];

  /* Check successive pairs of datapoints
   * until the  surface is found or the
   * last observation is encountered. */
  while (++j < nobs) {
    if (obs[j]  > -8.9) {

      /* Update the "current" or "next" values. */
      x[1] = station.observ[(int)DE][j]; 
      y[1] = obs[j];

      /* If we perform a valid interpolation
       * between the initial and next values,
       * then return the result. Several non-
       * valid interpolations will be performed
       * in the process of searching for the
       * final value. */
      if ((val = hb_linterp(d0, x, y, 2)) > -9998.) {
	return (val);
      }

      /* Update the "initial" values for
       * interpolation happening one depth
       * level below. */
      x[0] = x[1];
      y[0] = y[1];
    }
  }  /* End while loop over all observations. */

  /* If we did not exit early with a valid
   * interpolation, then quit w/ a error flag. */
  return (-9999.);

} /* End project_on_depth() */

/****************************************************************************/
