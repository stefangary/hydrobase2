/*  hb_cdfmerge.c

................................................................................
                          *******  HydroBase 2 *******
                             Ruth Curry
                             Woods Hole Oceanographic Institution
                             Dec 2000
...................................................
*
*  Merges any number of HydroBase2 cdf files into one file with appropriate bounds.
*  Obviously grid increments must be compatible amongst the files.
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "hydrobase.h"
#include "hydro_cdf.h"
#include "hb_memory.h"
#include "hb_paths.h"

/* input_file pathnames */

#define    EXTENT   ""
#define    DIR      ""

#define    PREFILL      1       /* turn on prefilling of cdf variables in cdf_construct() */

/* global variables for topography grid*/	

int lon0to360, xgreenwich;
float  *topo_lat, *topo_lon;
short **ztopo;
int 	 depth_is_neg;   
int     ncols, nrows;
double   xmin, xmax, ymin, ymax, delta_x, delta_y;

float HB_f_mask;          /* float values are used since the output */
float HB_f_empty;         /* cdf values are of this type */

/* prototypes for locally defined functions */	

void print_usage(char *);
int parse_p_option(char *, int *);
int cdf_construct(struct CDF_HDR *, int, int *, char *, int, char **);
void get_bounds(struct CDF_HDR *, int, char **, char *, char * );
void get_profiles( int, int, struct CDF_HDR *,  int ,  int *);
short **get_topo(FILE *);
float find_nearest_topo_val(float, float, short **);
void add_bottom_mask(int, struct CDF_HDR *);

main(int argc, char **argv)
{
  int i, j, n, nprops;
  int nfiles, curfile;
  int  error;
  int bflag, oflag;
  int infile, outfile, print_msg; 
  int *prop_req;
  char *st, *dir, *extent;
  char *outfile_name;
  float *xtmp;
  struct CDF_HDR *cdfin, hout;
   FILE   *topofile;
 
/* Set these default values. */
    
  error = 0;
  dir = DIR;
  extent = EXTENT;
  nfiles = 0;
  print_msg = 1;
  bflag = oflag = nprops = 0;
  lon0to360 = 0;
  xgreenwich = 0;
  depth_is_neg = 1;
  xmin = 0;
  xmax = 360;
  ymin = -90.0;
  ymax = 90.0;
  delta_x = delta_y = 0.1;
  HB_f_mask = (float)HBMASK;
  HB_f_empty = (float)HBEMPTY;
 
/*----------------------------------------*/  
  
/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
/*----------------------------------------*/  
/* parse command line arguments */

  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      
        case 'B':  /* get grid bounds */
	  bflag = 1;
          st = &argv[i][2];
          if (*st == '/')
             ++st;
          error = (sscanf(st,"%f", &hout.xmin) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%f", &hout.xmax) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%f", &hout.ymin) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%f", &hout.ymax) != 1);
                        	     
	  if (hout.xmin > hout.xmax) {
	    fprintf(stderr,"\nWest bound must be numerically <= east bound");
	    exit(1);
	  }
	  
	  if (hout.ymin > hout.ymax) {
	    fprintf(stderr,"\nNorth bound must be numerically >= south bound");
	    exit(1);
	  }
          break;
 
        case 'D':                   /* get input dir */
          dir = &argv[i][2];
          break;
	  
        case 'E':                    /* get file extent */
          extent = &argv[i][2];
          break;
  	  
        case 'L':        /* force sign of longitude  */
	  switch (argv[i][2]) {
	      case '+':
	          lon0to360 = 1;
		  break;
	      case '-':
	          lon0to360 = -1;
	      case '0':
	          lon0to360 = 0;
		  xgreenwich = TRUE;
		  break;
	      default:
	          fprintf(stderr,"\nUse -L+  or -L- to force sign of  longitude");
	          fprintf(stderr,"\nUse -L0 (zero) to use mixed sign longitudes (implies crossing Greenwich ");
	          error = 1;
	  }
	  break;
	  
        case 'O':
	  oflag = 1;
          outfile_name = &argv[i][2];
          break;
	  	      
        case 'P':
	  prop_req = (int *) calloc(MAXPROP, sizeof(int));
          nprops = parse_p_option(&argv[i][2], prop_req);
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
/*  Check syntax of options */ 

   
   error = 0;
   
    if (!nfiles ) {
       fprintf(stderr,"\nYou must specify input cdf_files as first argument(s).");
       ++error;    
    }
    if (!oflag ) {
       fprintf(stderr,"\nYou must specify -O<output_cdf_file> ");
      ++error;
    }
    if (!nprops ) {
       fprintf(stderr,"\nYou must specify output properties -P<property_list> ");
      ++error;
    }
                                     
   if (error) {
    fprintf(stderr,"\nUse -h for complete usage info. \n");
    exit(1);
   }
 
    if (bflag ) {
       fprintf(stderr,"\nMerged file bounds will be  [west:%8.3f ]  [east:%8.3f ]  [south:%8.3f]  [north:%8.3f] ", hout.xmin, hout.xmax, hout.ymin, hout.ymax);
       
       if (hout.xmin >= 0)
           lon0to360 = 1;
       else  {
               lon0to360 = -1;
	       if (hout.xmax > 0) {
	          xgreenwich = TRUE;
		  lon0to360 = 0;
		}  
	}
    }

/* Get standard depths from first cdf file in list */
   
   i = 1;
   infile = -1; 
   while (infile < 0 && i <= nfiles) {
     infile = cdf_open(dir, argv[i], extent, FALSE);
     ++i;
   }
   if (infile < 0) {
      fprintf(stderr,"Unable to open any cdf files.");
      exit(1);
    }
       
   cdfin = (struct CDF_HDR *) calloc(1, sizeof(struct CDF_HDR));
   
   if (error = read_cdf_hdr(infile, cdfin)) 
         exit (1);
   
   xtmp = (float *) calloc((size_t)cdfin->nz, sizeof(float));
   NSTDLEVS = read_cdf_depths(infile, xtmp);
   for (i = 0; i < NSTDLEVS; ++i) 
       std_depth[i] = (double) xtmp[i];
       
   std_depth_initialized = 1;
   free(xtmp);
   
   /* set these other grid specs */
   
   hout.yincr = cdfin->yincr;
   hout.xincr = cdfin->xincr;
   hout.node_offset = cdfin->node_offset;
   hout.counts_included = cdfin->counts_included;
   
   cdf_close(infile);
   free((void *) cdfin);
      
/*-------------------------------------------*/   
/*   set up CDF_HDR and output file */
/* Get min/max bounds of all files */

   if (!bflag) {
        get_bounds(&hout, nfiles, argv, dir, extent );
	fprintf(stderr,"\nOutput file bounds will be  [west:%8.3f ]  [east:%8.3f ]  [south:%8.3f]  [north:%8.3f] ", hout.xmin, hout.xmax, hout.ymin, hout.ymax);
   }	
   hout.nz = NSTDLEVS;
   hout.nx = (int) NINT((hout.xmax - hout.xmin) / hout.xincr);
   hout.ny = (int) NINT((hout.ymax - hout.ymin) / hout.yincr);
   if (hout.node_offset == 0) {  /* not pixel grid registration */
     ++hout.nx;
     ++hout.ny;
   }
   outfile = cdf_construct(&hout, nprops, prop_req, outfile_name, argc, argv);

  if (hout.node_offset)
     fprintf (stderr, "\nUsing %s registration", "pixel");
  else
     fprintf (stderr, "\nUsing %s registration", "gridnode");
     
     
   fprintf (stderr, "\nGrid dimensions are nx = %d, ny = %d", hout.nx, hout.ny);
      
/* read in topography file ... */

   topofile = fopen(BATHPATH,"r");
   if (topofile == NULL) {
          fprintf(stderr,"\nUnable to open default topography file\n", BATHPATH);
          exit(1);
   }
    ztopo = get_topo(topofile);
    fclose(topofile);
    
/* add seafloor/land masking info to output file ... */

   fprintf(stderr,"\nAdding topography masking info...");
    add_bottom_mask(outfile, &hout);
    
/* Traverse list of input files and merge info into new output file */

   fprintf(stderr,"\nReading in profiles...");
   
   curfile = 1;
   do {
       infile = cdf_open(dir, argv[curfile], extent, print_msg);
      if (infile < 0)
         goto NEXTFILE;
	 
      get_profiles(infile, outfile, &hout, nprops, prop_req);
      cdf_close(infile); 
      
NEXTFILE:
      ;   
   } while (curfile++ < nfiles);
    
   cdf_close(outfile); 
   
   fprintf(stderr,"\nEnd of %s.\n", argv[0]);
   exit(0);
	
}  /* end main */


/************************************************************************/
void print_usage(char *program)
{

  fprintf(stderr,"\nMerges multiple HydroBase gridded cdf files into one file with appropriate bounds");
  fprintf(stderr,"\nGrid increments and standard depths must be compatible among all files. ");
  fprintf(stderr,"\nAny gridnodes not covered in the input files are flagged as missing in output file. ");
  
  fprintf(stderr,"\n\nUSAGE:  %s cdf_file(s) -O<outfile> -P<properties> [-B<w/e/s/n>]  [-D<input_dir>] [-E<input_file_extent>]  [-L<+|-|0>]  [-h] [-:]\n\n", program);
  fprintf(stderr," -O  name of output cdf file.\n");
  fprintf(stderr," -P  list of properties to include in output file\n");
  fprintf(stderr,"        ex:  -Ppr/th/sa/ox/ht\n");
  fprintf(stderr,"       -P (by itself) produces a list of available properties\n");
  fprintf(stderr, "\n\tOPTIONS:\n");
  fprintf(stderr,"-B  sets the output grid bounds w/e/s/n.\n");
  fprintf(stderr,"-D directory for input cdf files.\n");
  fprintf(stderr,"-E file extent for input cdf files.\n");
  fprintf(stderr,"-L specify sign of longtitudes to be positive (0 to 360), negative (-360 to 0), \n     or mixed (western bound negative) \n");
  fprintf(stderr,"        ex:  -L+  [for 0 to 360]\n");
  fprintf(stderr,"        OR -L- [ lon all negative]  \n");
  fprintf(stderr,"        OR  -L0 [cross Greenwich Meridian] \n");
  fprintf(stderr,"-h help...... prints this message. \n");
  return;
  
} /* end print_usage() */
/*****************************************************************************/
int parse_p_option(char *st, int *prop_indx)
{
  char prop[4];
  int n;

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

     /* !**!  Special cases for properties ...
     
       No need for pref specs for props like s_ , ht and pe.  No computation of 
       additional props will be done. Any property 
       requested must already exist in the cdf file.  */
     
         /* de is automatically done so don't count it here */

     if ( !(prop_indx[n] == (int)DE ) )  
            ++n;

   } while (*st == '/');
   return (n);
}  /* end parse_p_option() */

/*****************************************************************************/
int cdf_construct(struct CDF_HDR *hptr, int nprops, int *prop_indx,  char *filename, int nargs, char **arglist)

   /* Opens a cdf output file and writes an appropriate header,
      standard depths, and time bin info.  Returns the id associated
      with the open file */
{
   int i, error, cdfid;
  
   hptr->nt = 1;
   hptr->tmin = (int *) malloc(sizeof(int));
   hptr->tmax = (int *) malloc(sizeof(int));
   hptr->tmin[0] = 0;
   hptr->tmax[0] = 9999;
   hptr->fill_value =  (float) HBEMPTY;
   hptr->mask_value = (float)  HBMASK;
   strncpy(hptr->x_units, "degrees", 8);
   strncpy(hptr->y_units, "degrees", 8);
   strncpy(hptr->z_units, "meters", 7);
   strncpy(hptr->title,"HydroBase", 10);
   strcpy(hptr->command, *arglist);
   for (i = 1; i < nargs; ++i) {
      strncat(hptr->command, " ", 1);
      strcat(hptr->command, arglist[i]);
   }
   
   
   hptr->nprops = nprops;
   hptr->prop_id = (char **) malloc(nprops * sizeof(char *));
   hptr->prop_units = (char **) malloc(nprops * sizeof(char *));
   for (i = 0; i < nprops; ++i) {
      hptr->prop_id[i] = (char *) malloc(3);
      hptr->prop_units[i] = (char *) malloc(50);
      strncpy(hptr->prop_id[i], get_prop_mne(prop_indx[i]),3);
      strcpy(hptr->prop_units[i], get_prop_units(prop_indx[i]));
   }
   
/* Open output file and write out some info so we can free up some memory... */
   
   cdfid = cdf_init(filename);   
   error = cdf_define(cdfid, hptr, PREFILL, hptr->counts_included);
   if (error)  exit(1);
   
   error = write_std_depths_cdf(cdfid, hptr);
   error = write_time_bins_cdf(cdfid, hptr);
     
   return(cdfid); 
} /*end cdf_construct() */ 

/************************************************************************/
void get_bounds(struct CDF_HDR *hptr, int nfiles, char **argv, char *dir, char *ext )

 /*  Reads each cdf file in list and determines boundaries of output grid and gridnode spacing */

{
   struct CDF_HDR cdf;
   int i, index1, index2, error;
   int merid[360];
   int cdfid, curfile, print_msg = 0;

/*initialize these */
  for (i = 0; i < 360; ++i)
     merid[i] = 0;
     
  hptr->ymax = -90;
  hptr->ymin = 90;   
  hptr->xmax = 0;
  hptr->xmin = 360;   

  
/* Loop for each input file */

   curfile = 1;
   
   do {
      cdfid = cdf_open(dir, argv[curfile], ext, print_msg);
      if (cdfid < 0)
         goto NEXTFILE1;
	 
      if (error = read_cdf_hdr(cdfid, &cdf)) 
         exit (1);
	
/* find min/max boundaries. */
	 
        if (cdf.xmin < 0)   {     /* initially set all longitudes positive */
	   cdf.xmin += 360.0;
	   cdf.xmax += 360.0;
       }
       index1 = NINT (cdf.xmin);
       index2 =  NINT(cdf.xmax - 1);
       if (index2 < 360) {
           for (i = index1; i <= index2; ++i)
              merid[i] = 1;
       }
       else {      /* crossed greenwich */
           for (i = index1; i < 360; ++i)
              merid[i] = 1;
	   cdf.xmax -= 360.0;
	   index2 -= 360;
           for (i =0; i <= index2; ++i)
              merid[i] = 1;
       }

      if ( cdf.ymin < hptr->ymin)
          hptr->ymin = cdf.ymin;
      if ( cdf.ymax > hptr->ymax)
          hptr->ymax = cdf.ymax;


       cdf_close(cdfid);
       
NEXTFILE1:
      ;
      	 
   }  while (curfile++ < nfiles); 

    if (merid[0] && merid[359]) {   /* either crossed Greenwich, or bounds are entire 360 of longitude */
      /*search backward from Greenwich to find western bound and make it negative*/
       i = 359;  
       while ( (i >= 0) && (merid[i] > 0))
           --i;
	   
	if (i < 0 ) {   /* entire 360 of lon */
	    hptr->xmin = 0.0;   
	    hptr->xmax = 360.0;
	    lon0to360 = 1;
	    xgreenwich = FALSE;
	    return;
	 }  
	  
         hptr->xmin = (float) (i + 1) - 360.;
	    
	      /* continue searching backward for east boundary and make it positive */
	  while (i > 0 && merid[i] == 0)
	      --i;
	  hptr->xmax = (float) (i+1);
	  fprintf(stderr, "\nCrossed Greenwich ");
	  xgreenwich = TRUE;
	  lon0to360 = 0;
	  return; 
     }
     
     /* search forward to find western bound */
     i = 0;
     while ( (i < 360) && (merid[i] == 0))
        ++i;

    if (i == 360) {
       fprintf(stderr,"\nNo longitude bounds could be determined:  SOFTWARE BUG 1\n");
       exit(1);
    }
     
     hptr->xmin = (float) i;
     while ((i < 360) && (merid[i] == 1) )  
        ++i;
	
     hptr->xmax = (float) i;
     
     if (hptr->xmax <= hptr->xmin) {
        fprintf(stderr,"\nError determining longitude bounds:  SOFTWARE BUG 2\n");
       exit(1);
     }
     
     if (lon0to360 < 0) {  /* force negative longitudes */
        hptr->xmax -= 360;
        hptr->xmin -= 360;
        return;
     }	  	         
     if (lon0to360 == 0) {  /* if lons are all positive, set lon0to360 to 1) */
        if ((hptr->xmin >= 0) && (hptr->xmax > 0) && (hptr->xmax <= 360) )
	   lon0to360 = 1;
     }
     return;
      
} /* end get_bounds() */
/************************************************************************/

void get_profiles( int cdfid, int outfile, struct CDF_HDR *hptr,  int nprop_out,  int *propindx_out)

 /*  Reads  cdf file, reads each property profile and counts and writes to output file   */

{
   struct CDF_HDR cdf;
   float *xin_ptr;
   short *nobs_ptr;
   char *mne;
   int error, i,  j, row, col;
   int rowout, colout, index;
   int tbin = 0;
   float bdepth;
   float lat, lon;
   
   if (error = read_cdf_hdr(cdfid, &cdf)) 
         exit (1);
	 
	 
      if (lon0to360 > 0) {    /* longitude range is positive */
        if (cdf.xmin < 0)
	   cdf.xmin += 360.0;
        if (cdf.xmax < 0)
	   cdf.xmax += 360.0;
     
      }
      else if (lon0to360 < 0) {     /* longitude range is negative */
           if (cdf.xmin > 0)
	      cdf.xmin -= 360.0;
           if (cdf.xmax > 0)
	      cdf.xmax -= 360.0;
      }
      else {
          if (xgreenwich) {    /* mixed range with a negative western bound */
	     if (cdf.xmin > hptr->xmax) {
	          cdf.xmin -= 360.0;
	          cdf.xmax -= 360.0;
	     }
	  }
      }

   /* check that cdf file has same number of standard depths  */
     
      if (cdf.nz != NSTDLEVS) {
         fprintf(stderr, "\nFATAL ERROR:  Mismatch of standard depths in cdf files.\n");  
         exit(1);
      }

/* allocate space for property arrays; */
      xin_ptr = (float *) get_memory((void *)NULL, (size_t)cdf.nz, sizeof(float));
      nobs_ptr = (short *) get_memory((void *)NULL, (size_t)cdf.nz, sizeof(short *));
     
/* for each property, visit every row,col in input file.  Read profile in, 
       find  output row,col, 
      determine row/col in output file and write it */

    for (i = 0; i < nprop_out; ++i ) {
	  index = propindx_out[i];
	    /* is requested property available? */
	  j = 0;
	  while (j < cdf.nprops && ! (index == get_prop_indx(cdf.prop_id[j])))
	          ++j;
	  if (j >= cdf.nprops) 
	          fprintf(stderr,"\nRequested property [%2s] not available ", get_prop_mne(index) );
	  
	  if ( j < cdf.nprops) {
	         mne = cdf.prop_id[j];
	  
                 for (row = 0; row < cdf.ny; ++row ) {
                     for (col = 0; col < cdf.nx; ++col) {
	                if (get_lat_lon(&cdf, row, col, &lat, &lon) < 0 ) {
	                    fprintf(stderr,"\nFATAL ERROR translating row, col [%d,%d] into lat,lon\n", row, col);
		            exit(1);
	                }
			if (lon > 359.9999)
			    lon -= 360.0;
			    
	                if (get_indices(hptr, lat, lon, &rowout, &colout) < 0) {
	                    fprintf(stderr,"\nFATAL ERROR translating lat,lon [%8.3f,%8.3f] into row, col\n", lat, lon);
		            exit(1);
	                }
	                if (rowout >= hptr->ny || rowout < 0 || colout >= hptr->nx || colout < 0) {
	                    fprintf(stderr,"\n WARNING: lat,lon  [%8.3f,%8.3f] not within output bounds w/e/s/n [%8.3f/%8.3f/%8.3f/%8.3f] \n", lat, lon, hptr->xmin, hptr->xmax, hptr->ymin, hptr->ymax);
			    continue;
	                }
                        error = read_cdf_prop(cdfid, mne, xin_ptr, row, col, tbin, 0, cdf.nz);
                        if (error > 0) {
                              fprintf(stderr,"\nError attempting to read %.2s at row,col =  %d,%d from cdf file.", mne, row, col);
                              exit(1);
                        }
	       
	                write_prop_cdf(outfile, xin_ptr, mne, rowout, colout, tbin, 0, 1, 1, 1, NSTDLEVS);
	   
                       if (cdf.counts_included) {
                          error = read_cdf_prop_count(cdfid,mne,nobs_ptr,row,col,tbin,0,cdf.nz);
                          if (error > 0) {
                              fprintf(stderr,"\nError attempting to read %.2s_cnt at row,col =  %d,%d from cdf file.", mne, row, col);
                              exit(1);
                          }
		          write_prop_count_cdf(outfile, nobs_ptr, mne, rowout, colout, tbin, 0, 1, 1, 1, cdf.nz);
                       }
	      
	           }  /*end for col */
               } /* end for row */   
	 }  /* end if */     
      } /* for i */

     
     free((void *) xin_ptr);
     free((void *) nobs_ptr);
     
     /* Read in bottom depths and write to output file */
   
     for (row = 0; row < cdf.ny; ++row ) {
         for (col = 0; col < cdf.nx; ++col) {
	      if (get_lat_lon(&cdf, row, col, &lat, &lon) < 0 ) {
	          fprintf(stderr,"\nFATAL ERROR translating row, col [%d,%d] into lat,lon\n", row, col);
	          exit(1);
	      }
	      if (lon > 359.9999)
	         lon -= 360.0;
		 
	      if (get_indices(hptr, lat, lon, &rowout, &colout) < 0) {
	          fprintf(stderr,"\nFATAL ERROR translating lat,lon [%8.3f,%8.3f] into row, col\n", lat, lon);
		  exit(1);
	      }
	      if (rowout >= hptr->ny || rowout < 0 || colout >= hptr->nx || colout < 0) {
	          fprintf(stderr,"\n WARNING: lat,lon  [%8.3f,%8.3f] not within output bounds w/e/s/n [%8.3f/%8.3f/%8.3f/%8.3f] \n", lat, lon, hptr->xmin, hptr->xmax, hptr->ymin, hptr->ymax);
	          continue;
	       }
    
 	       error = read_cdf_bottom_depth(cdfid, &bdepth, row, col, tbin);
               if (error < 0) {
                     fprintf(stderr,"\nError attempting to read bottom depth at row,col =  %d,%d from cdf file.",  row, col);
                      exit(1);
               }
	       write_bottom_depth_cdf(outfile, rowout, colout, tbin, 1, 1, 1, &bdepth);
	 }  /*end for col */
     } /* end for row */   
     
   return;

}  /* end get_profiles() */

/****************************************************************************/
short **get_topo(FILE *fptr)
   /* allocates memory and reads in topography values.  Returns a pointer to
      the start of the memory block. An error causes an error message to be
      printed and an exit.   */
{

   int row, n;
   short **z;
   
   /* Set some globally defined variables */
   
   nrows = (int)NINT((ymax - ymin) / delta_y) + 1;
   ncols = (int)NINT((xmax - xmin) / delta_x) + 1;
   
   topo_lat = (float *) malloc(nrows * sizeof(float));
   topo_lon = (float *) malloc(ncols * sizeof(float));
   
   for (n = 0; n < nrows; ++n) {
     topo_lat[n] = (float) ((double) n * delta_y + ymin);
   }
   
   for (n = 0; n < ncols; ++n) {
     topo_lon[n] = (float) ((double) n * delta_x + xmin);
   }
   
   z = (short **) malloc(nrows * sizeof(short *));
   if (z == NULL) {
      fprintf(stderr,"\nError allocating memory.\n");
      exit(1);
   }
   
   for (row = 0; row < nrows; ++row) {
     z[row] = (short *) malloc(ncols * sizeof(short));
     if (z[row] == NULL) {
         fprintf(stderr,"\nError allocating memory.\n");
         exit(1);
      }
   }
   
   fprintf(stderr,"\nReading in topography values ");
   
   for (row = 0; row < nrows; ++row) {
     n = fread((short *)z[row], sizeof(short), ncols, fptr);
     if (n != ncols) {
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

float find_nearest_topo_val(float lat, float lon, short **topo)
   /* interpolates nearest values in topo matrix and returns bottom
      depth for lat/lon specified.  Depth returned is positive below 
      sea level/ negative above. */
      
{
  int row, col, row2, col2;
  short z;
  double zn, zs, dlat, dlon;
  
  dlat = (double) lat;

  if (lon0to360 <= 0) {  
     if (lon < 0) lon += 360.;
  }
  dlon = (double) lon;
  
  row = (int) ((dlat - ymin) / delta_y);
  col = (int) ((dlon - xmin) / delta_x);
  
  /* row and row+1 should bracket the lat; col and col+1 the lon */
  
  if ((row < 0) || (row >= nrows)) {
    fprintf(stderr,"\nRow (%d) computed for latitude (%7.3f) does not exist in topo file.\n", row, lat);
    exit(1);
  }
  
   if ((col < 0) || (col >= ncols)) {
    fprintf(stderr,"\nColumn (%d) computed for longitude (%8.3f) does not exist in topo file.\n", col, lon);
    exit(1);
  }
  
  col2 = col + 1;
  row2 = row + 1;
  
  if (col2 >= ncols) {
     fprintf(stderr,"\nDid not expect to have longitude right on xmax -- report bug 1 in find_nearest_topo_val()\n");
     exit(1);
  }
  
  zs = (double)topo[row][col] + (double)(topo[row][col2] - topo[row][col]) * (dlon - topo_lon[col]) / (double)(topo_lon[col2] - topo_lon[col]);
  
  if (row2 >= nrows) {  /* Must be at the north pole */
      z = (short) zs;
      if (depth_is_neg)
         z = -z;
      return(z);
   }
      
  zn = topo[row2][col] + (short) ((double)(topo[row2][col2] - topo[row2][col]) * (dlon - (double)topo_lon[col]) / (double)(topo_lon[col2] - topo_lon[col]) );
  
  z = (short) (zs + (zn - zs) * (dlat - (double)topo_lat[row]) / (double)(topo_lat[row2] - topo_lat[row]));
  
  if (depth_is_neg)
    z = -z;
      
  return( (float) z);
 

}  /* end find_nearest_topo_val() */
/**************************************************/
void add_bottom_mask(int outfile, struct CDF_HDR *hptr)

/* traverse output file, adding bottom mask values where appropriate */

{
    float *x_ptr;
   char *mne;
   int error, i,  np,  row, col;
   int tbin = 0;
   float bdepth;
   float lat, lon;
  
/* allocate space for property arrays; */
      x_ptr = (float *) get_memory((void *)NULL, (size_t)hptr->nz, sizeof(float));
   
/* visit every row,col in output file, find bottom depth for gridnode, and fill in empty value 
   -- or mask value everywhere beneath bottom depth */

     for (row = 0; row < hptr->ny; ++row ) {
         for (col = 0; col < hptr->nx; ++col) {
	       if (get_lat_lon(hptr, row, col, &lat, &lon) < 0 ) {
	           fprintf(stderr,"\nFATAL ERROR translating row, col [%d,%d] into lat,lon\n", row, col);
		   exit(1);
	       }
   	       
	       bdepth =  find_nearest_topo_val(lat, lon, ztopo);
	       
	       for (i = 0; i < hptr->nz; ++i ) {
	           x_ptr[i] = HB_f_empty;
	           if (std_depth[i] > bdepth )
		        x_ptr[i] = HB_f_mask;
	       }
	       /* explicitly mask deepest level only if seafloor is land */
	        x_ptr[hptr->nz - 1] = HB_f_empty;
	       if (bdepth <= 0)
	           x_ptr[hptr->nz - 1] = HB_f_mask;
		   
	       for (np = 0; np < hptr->nprops; ++np)  {
	          mne = hptr->prop_id[np];
	          write_prop_cdf(outfile, x_ptr, mne, row, col, tbin, 0, 1, 1, 1, NSTDLEVS);
	       }
	       write_bottom_depth_cdf(outfile, row, col, tbin, 1, 1, 1, &bdepth);
	      
	  }  /*end for col */
      } /* end for row */   
     return;
} /* end add_bottom_mask () */
