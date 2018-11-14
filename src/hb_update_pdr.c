/* hb_update_pdr.c
................................................................................
                         *******  HydroBase2 *******
................................................................................
                          
                    author:  Ruth  Curry
                             Woods Hole Oceanographic Inst
                             July 2000 
................................................................................
................................................................................
.  Reads in a file of global topography. Then for each HydroBase station file
.  specified, the pdr field of the station headers are checked.  If the pdr
.  value is <= 0 OR it is > 0 but < deepest observed depth, that pdr value is
.  assigned EITHER the depth from the topography database corresponding to 
.  the station lat/lon OR the deepest observed depth -- whichever is deeper.  
.
.  The updated station file is named by adding the extent *.up to the original
.  filename.
.    
.
................................................................................
................................................................................

*/
#include <stdio.h>
#include <stdlib.h>   /*sfg added to deal with exit and malloc warnings*/
#include <ctype.h>
#include <string.h>
#include "hydrobase.h"
#include "hb_paths.h"

#define    DIR     ""
#define    EXTENT   ""
#define    PRINT_MSG 1


     /* boundaries for grid */

double   xmin, xmax, ymin, ymax, delta_x, delta_y;
float  *topo_lat, *topo_lon;
int 	lon0to360, depth_is_neg;   
int     ncols, nrows;


  /* prototypes for functions declared locally */
  
void print_usage(char *);
short find_nearest_topo_val(float, float, short **, int);
short **get_topo(FILE *, int);

main (int argc, char **argv)
{
   int     curfile = 1, nfiles = 0; 
   int     i, infile, outfile;
   int     i_flag, b_flag, t_flag, v_flag;
   int     status, force_it;
   int     error = 0;
   short **ztopo, z;
   struct HYDRO_HDR h;
   struct HYDRO_DATA data;
   char   *dir, *extent, *st;
   char    fname[200];
   FILE   *topofile;
   int     nstation;
   int     nedited;


/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    i_flag = b_flag = t_flag = v_flag = 0;
    lon0to360 = 1;
    depth_is_neg = 1;
    xmin = 0;
    xmax = 360;
    ymin = -90.0;
    ymax = 90.0;
    delta_x = delta_y = 0.1;
    force_it = 0;
    nstation = 0;
    nedited = 0;

/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'B':                    /* get grid bounds */
                        b_flag = 1;
                        st = &argv[i][2];
                           if (*st == '/')
                               ++st;
                        error = (sscanf(st,"%lf", &xmin) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%lf", &xmax) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%lf", &ymin) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%lf", &ymax) != 1);

                        if (xmin > xmax) {
                          fprintf(stderr,"\nW bound cannot exceed E bound.\n");
			  exit(1);
			}     
                        if (ymin > ymax) {
                          fprintf(stderr,"\nS bound cannot exceed N bound.\n");
			  exit(1);
			}
                        if (xmin < 0)
                           lon0to360 = 0; 
                        break;

               case 'D':
                        dir = &argv[i][2];
                        break;

               case 'E':
                        extent = &argv[i][2];
                        break;

               case 'F' :
                        force_it = 1;
                        break;

               case 'I' :
                        i_flag = 1;
                        error = (sscanf(&argv[i][2],"%lf", &delta_x) == 1) ? 0 : 1;
                        delta_y = delta_x;
                        st = strchr(&argv[i][2],'/');
                        if (st != NULL) {
                          sscanf(++st,"%lf", &delta_y);
                        }
                        break;

               case 'N' :
                        depth_is_neg = 0;
                        break;

               case 'T':
                        topofile = fopen(&argv[i][2],"r");
                        if (topofile == NULL) {
                          fprintf(stderr,"\nError opening %s for input\n", &argv[i][2]);
                          exit(1);
                        }
                        t_flag = 1;
                        break;

	       case 'V':
		        v_flag = 1;
			break;

               case 'h':
                        print_usage(argv[0]);
                        exit(0);
               default  :
                        fprintf(stderr,"\nError parsing command line");
                        fprintf(stderr,"\n in particular: %s\n", argv[i]);
                        exit(1);
            }  /* end switch */
       }  /* end if */
       else  {
           ++nfiles;
       }
   }  /* end for */

   if ( !t_flag ) {   
      topofile = fopen(BATHPATH,"r");
      if (topofile == NULL) {
          fprintf(stderr,"\nUnable to open default topography file\n", BATHPATH);
          fprintf(stderr,"You can specify the name of a binary topography file for input with [-T]\n"); 
          exit(1);
      }
    
   } 
   if ( (!i_flag) && v_flag ) {
      fprintf(stderr,"\nIncrement for topofile: [%.2lf/%.2lf]", delta_x, delta_y);
   } 
   
   if ( (!nfiles) && v_flag ) {
      fprintf(stderr,"\n Expecting station files from stdin.... ");
      infile = STDIN;
      outfile = STDOUT;
   }

/* initialize these... */
   for (i = 0; i < MAXPROP; ++i)
      data.observ[i] = (double *) NULL;
   h.prop_id = (int *) NULL;
   
/* read in topography file ... */
   ztopo = get_topo(topofile, v_flag);
    fclose(topofile);

 /* loop for each input file */
   do {

     if (nfiles) {
       if ( v_flag ) {
	 infile = open_hydro_file(dir, argv[curfile], extent, PRINT_MSG);
       }
       else {
	 infile = open_hydro_file(dir, argv[curfile], extent, 0);
       }
       if (infile  < 0) 
	 goto NEXTFILE;

       strcpy(fname, dir);
       if ((i = strlen(dir)) != 0) {
	 if (dir[i-1] != '/')
	   strncat(fname,"/",1);
       }
       strcat(fname, argv[curfile]);
       strcat(fname,extent);
       strncat(fname,".up",3);
       outfile = create_hydro_file(fname, OVERWRITE);
       if (outfile < 0)
	 exit(1);

     }
     
     /* loop for each station */

     while ((status = get_station(infile, &h, &data)) == 0) { 

       nstation++;

        if ((force_it) || (h.pdr <= 0)|| (h.pdr < data.observ[(int)DE][h.nobs-1])) {
	  nedited++;
	  z = find_nearest_topo_val(h.lat, h.lon, ztopo, v_flag);
           h.pdr = (int) z;
           if ((h.pdr == 0) || (h.pdr < data.observ[(int)DE][h.nobs-1])) 
              h.pdr = (int) data.observ[(int)DE][h.nobs-1] + 10;
        } 

        error = write_hydro_station(outfile, &h, &data);
        if (error) {
          fprintf(stderr,"\nError in write_hydro_station(): error code = %d.\n", error);
          exit(1);
        }

     }  /* end while */

     report_status(status, stderr);
     close(infile);
     close(outfile);

NEXTFILE:
     ;
   } while (curfile++ < nfiles);

   if ( v_flag ) {
     fprintf(stderr,"\n\n %d stations read in and %d stations edited.\n",nstation,nedited);

     fprintf(stderr,"\n\nEnd of %s.\n\n", argv[0]);
   }
   exit(0);

} /* end main() */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\n%s fills in missing pdr fields for HydroBase station files.", program);
   fprintf(stderr,"\nIf no station file is specified, stdin and stdout ");
   fprintf(stderr,"\nare used for input/output.  If a filename is");
   fprintf(stderr,"\nspecified, an output file is created with extent .up");
   fprintf(stderr,"\nappended to the original filename.");
   fprintf(stderr,"\nThe default topography file is %s", BATHPATH);
   fprintf(stderr,"\nFor an alternate topography file, use -T to specify its name,");
   fprintf(stderr,"\n-B to specify its bounds, and -I for its grid increment");
   fprintf(stderr,"\nUse -N if the alternate topography has positive seafloor depths.");
   
   fprintf(stderr,"\n\nUsage:  %s [filename_root(s)]", program);
   fprintf(stderr,"\n   [-Ddirname] [-Eextent] [-F] [-T<topofile>]");
   fprintf(stderr,"\n   [-B<west/east/south/north>] [-I<xincr/yincr>]");
   fprintf(stderr,"\n   [-N] [-V] [-h]");
   fprintf(stderr,"\n");
   
   fprintf(stderr,"\tOPTIONS:");
   fprintf(stderr,"\n[-B] : specifies bounds of topography grid. Default [0/360/-90/90]");
   fprintf(stderr,"\n[-D] : directory of input station files ");
   fprintf(stderr,"\n        ex: -D../data/ ");
   fprintf(stderr,"\n[-E] : specifies file extent (default is no extent)");  
   fprintf(stderr,"\n        ex: -E.dat ");
   fprintf(stderr,"\n[-F] : force pdr field to be updated according to topography.");
   fprintf(stderr,"\n         Default only updates if pdr field <= 0 or < bottom observed depth");
   fprintf(stderr,"\n[-I] : specify x-increment [/y-increment] for topography grid");
   fprintf(stderr,"\n         Default [%.2lf/%.2lf]", delta_x, delta_y);
   fprintf(stderr,"\n[-T] : name of topography file. Default is [%s]", BATHPATH);
   fprintf(stderr,"\n[-N] : specify this flag if topography file has positive seafloor depths.");
   fprintf(stderr,"\n[-V] : will run in verbose mode (update messages to stderr).");
   fprintf(stderr,"\n[-h] : help -- prints this message");
   fprintf(stderr,"\n\n");  
   return;
}

/****************************************************************************/

short **get_topo(FILE *fptr, int v_flag)
   /* allocates memory and reads in topography values.  Returns a pointer to
      the start of the memory block. An error causes an error message to be
      printed and an exit.  If the v_flag is set to 1, then this routine
      will print out update messages.  Otherwise, no messages. */
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
   
   if ( v_flag ) {
     fprintf(stderr,"\nReading in topography values ");
   }

   for (row = 0; row < nrows; ++row) {
     n = fread((short *)z[row], sizeof(short), ncols, fptr);
     if (n != ncols) {
       fprintf(stderr,"\nError reading the topofile at row %d\n", row);
       exit(1);
     } /* End of good read check. */
     if ( v_flag ) {
       if ((row % 10) == 0) {  /* signal progress */
	 fprintf(stderr,".");
       } /* End of progress check. */
     } /* End of verbose check. */
   } /* End of loop over all rows. */

   if ( v_flag ) {
     fprintf(stderr,"\nFinished reading topofile.\n");
   }

   return(z);

} /* End get_topo() */
/****************************************************************************/

short find_nearest_topo_val(float lat, float lon, short **topo, int v_flag)
   /* interpolates nearest values in topo matrix and returns bottom
      depth for lat/lon specified.  Depth returned is positive below 
      sea level/ negative above. */
      
{
  int row, col, row2, col2;
  short z;
  double zn, zs, dlat, dlon;
  
  /* Compute double preciesion lat (deg N) and lon (deg E). */
  dlat = (double) lat;
  
  if (lon0to360) {  
    if (lon < 0) lon += 360.;
  }
  dlon = (double) lon;
  
  /* Find the index values of the domain that
   * are closest to the current location. Note:
   *  row and row+1 should bracket the lat
   *  col and col+1 should bracket the lon
   * since the int typecast will truncate any
   * figures to the right of the decimal pt. */
  row = (int) ((dlat - ymin) / delta_y);
  col = (int) ((dlon - xmin) / delta_x);
  
  if ((row < 0) || (row >= nrows)) {
    fprintf(stderr,"\nRow (%d) computed for latitude (%8.3f) does not exist in topo file.\n", row, lat);
    exit(1);
  }
  
   if ((col < 0) || (col >= ncols)) {
    fprintf(stderr,"\nColumn (%d) computed for longitude (%8.3f) does not exist in topo file.\n", col, lon);
    exit(1);
  }
  
  col2 = col + 1;
  row2 = row + 1;
  
  if (col2 >= ncols) {
    fprintf(stderr,"\nWARNING: Station location is at or beyond edge of topo map.");
    fprintf(stderr,"Limitation of find_nearest_topo_val()\n");
    exit(1);
  }
  
  /* Zonal linear interpolation along the southern edge. */
  zs = (double)topo[row][col] + ((double)(topo[row][col2] - topo[row][col])) * (dlon - (double)topo_lon[col]) / ((double)(topo_lon[col2] - topo_lon[col]));
  
  /* Check for north pole along northern edge.
   * If we encouter the pole, return southern
   * edge value and exit early. */
  if (row2 >= nrows) {
    z = (short) zs;
    if (depth_is_neg) {
      z = -z;
    }

    if (z <= 0) {
      if ( v_flag ) {
	fprintf(stderr,"\nWARNING: Interpolated value above sea level from topo file at: ");
      }
      fprintf(stderr,"%8.3fN %8.3fE\n", lat, lon);
    } /* End of above SL check. */
    
    return(z);
  }
     
  /* Zonal linear interpolation along the northern edge. */
  zn = (double)topo[row2][col] + ((double)(topo[row2][col2] - topo[row2][col])) * (dlon - (double)topo_lon[col]) / ((double)(topo_lon[col2] - topo_lon[col]));

  /* Original line - inconsistent typecasting? */
  /*zn = topo[row2][col] + (short) ((double)(topo[row2][col2] - topo[row2][col]) * (dlon - (double)topo_lon[col]) / (double)(topo_lon[col2] - topo_lon[col]) );*/
  
  /* Meridional linear interpolation between
   * northern and southern edges. */
  z = (short) (zs + (zn - zs) * (dlat - (double)topo_lat[row])/((double)(topo_lat[row2] - topo_lat[row])));
  
  if (depth_is_neg) {
    z = -z;
  } /* End of conversion of negative depths to positive depths. */

  if (z <= 0) {
    if ( v_flag ) {
      fprintf(stderr,"\nWARNING: Interpolated value above sea level from topo file at: ");
    }
    fprintf(stderr,"%8.3fN %8.3fE\n", lat, lon);
  } /* End of above SL check. */
  
  return(z);
 
}  /* end find_nearest_topo_val() */

/*************************************************/

