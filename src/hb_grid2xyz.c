/*  hb_grid2xyz.c
................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             1993
			     updated to ANSI Dec 1999
................................................................................
_______________________________________________________________________________

   USAGE:  hb_grid2xyz <input_filename> -O<output_file_root> -P<property_list>
             [-R<smoothing_radius>] [-N] [-S] [-B<blanking_file>]
_______________________________________________________________________________
.  Separates a HydroBase *.grid file (output by hb_gridsurf) into individual
.  *.xyz files containing lon, lat, property triplets.  Optionally smooths
.  each point with information from surrounding gridpoints (radius specified
.  with -R) weighted by a distance factor = 1/(d*d), where d = # of gridpts
.  away from the central node.  The number of observations
.  is included in the output if -N option is specified.  -S creates files
.  with std deviation information.
.
.  The program outputs up to 4 files for each property specified: 
.     1) property info:   lon, lat, av_value [, n]
.     2) blank gridpts:   lon, lat
.     3) stddev info  :   lon, lat, stddev [, n]
.     4) stddev blanks:   lon, lat
.  The stddev information is only output if the -S option is specified.
.
.   Output files are named:     <root><property_id>.xyz
.                                  "       "       .blk
.                                  "       "       _var.xyz
.                                  "       "       _var.blk
.   where root is supplied by the -O option and property_id is the same as
.   the property_list_ids in the -P option.
.
.   The -B option permits you to supply a file containing the lon,lat of
.   gridpoints to be blanked.  Its purpose is to mask out bathymetry,
.   coastlines, etc...  This list is checked before computing each new
.   gridpoint.
.
.   Stefan Gary added the -M option which allows for smoothing based
.   on the number of contributing observations at each grid point.
.   This was done to be in line with the methodology of Lozier, Owens,
.   and Curry (1995).
.
*/


#include <stdio.h>
#include <stdlib.h>   /*sfg added to deal with exit and malloc warnings*/
#include <math.h>
#include <string.h>
#include "hydrobase.h"


/* Data structures to represent gridpoints within
 * each grid matrix.
 * (1) In the old_grid_pt, mean, var, and n will
 *     be converted to sumpts, sumofsquares, and
 *     total number of points in grid box.
 * (2) The new_grid_pt holds running sums for
 *     computing the mean (sum), and stddev (sumsq).
 *     These are weighted averages of the form:
 *     weighted_avg = sum_over_i(weight_i * value_i)/sum_over_i( weight_i)
 *     where if the weight_i were all equal, we'd
 *     recover the pure mean.  The value wghts thus
 *     is the running sum of the weights and n is
 *     the total count of observations being included
 *     into this point. */
struct old_grid_pt {
  float mean;
  double var;
  int n;
};

struct new_grid_pt {
  float sum;
  double sumsq;
  int n;      
  float wghts;
};

/* Function primitives. */
void  add_wghted_layer(int, double, struct new_grid_pt *, int, int, struct old_grid_pt *, int, int);
void  print_usage(char *);

int main (int argc, char **argv) {
  /* Counters:
   * i = counts all possible properties.
   * j = property index for memonic
   * row = current row in matrix
   * col = current column in matrix
   * sq = square number in vectorized matrix
   *      (the matrices are actually stored
   *       in a long, "unwrapped" vector.) */
  int    row, col, n, i, j, sq;

  /* Current (search) smoothing
   *    radius and upper limit. */
  int    rad,  maxrad;

  double weight_factor; /* Smoothing weighting factor. */
  int    maxpts2smo;    /* Max number of points to smooth. */
  int    ncols, nrows;  /* Upper limits to col, row. */
  int    gridsize;      /* Upper limit to sq. */

  int    error, index, nprops;
  int    prop_index[MAXPROP];     
  char   prop_request[MAXPROP], prop_avail[MAXPROP];                      

  int    bflag, nflag, sflag; /* Flags for presence of */
  int    pflag, oflag, mflag; /* options on command line.*/

  /* id = holds the property mnemonic. */
  char   *id, field_descrip[12];
  double var, mean, stddev;
  double propsum, propsqsum;
  float  lat, lon;

  /* Grid information. */
  float  xmin, ymin, xmax, ymax, xspacing, yspacing;

  /* File links. */
  FILE  *infile, *outfile[4], *blankfile;

  char  *fname, *root, *s;
  char  *blank;

  /* _pt variables are single points.
   * _matrix variables are whole grids,
   *  one level for each possible variable. */
  struct old_grid_pt *gridpt, *old_matrix[MAXPROP];
  struct new_grid_pt *newpt, *new_matrix[MAXPROP];
  
  /* Check for command line arguments. */
  if (argc < 2) {
    print_usage(argv[0]);
    exit(1);
  }

  /* Initialize flags and defaults. */
  for (i = 0; i < MAXPROP; ++i) {
    prop_request[i] = prop_avail[i] = 0;
  }
  id = (char *) malloc(3);
  maxrad = 0;
  weight_factor = 0.5;
  maxpts2smo = 9999;
  oflag = pflag = mflag = 0;
  bflag = nflag = sflag = 0;
  error = 0;
  
  /* Parse the command line arguments */
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
        case 'O':                    /* get output file root_name */
	  oflag = 1;
	  root = &argv[i][2];
	  break;

        case 'P':                    /* get list of properties */
	  pflag = 1;
	  s = &argv[i][2];
	  if (*s == '\0') {
	    print_prop_menu();
	    exit(0);
	  }
	  do {
	    if (*s == '/')
	      ++s;
	    sscanf(s,"%2s", id);
	    index = get_prop_indx(id);
	    if (error = (index < 0) ? 1 : 0)
	      break;
	    prop_request[index] = '1';
	    ++s;
	    ++s;
	  } while (*s == '/');
	  break;

        case 'N':
	  nflag = 1;
	  break;

        case 'R':
	  error = (sscanf(&argv[i][2],"%d", &maxrad) == 1) ? 0 : 1;
	  break;
	  
        case 'W':
	  error = (sscanf(&argv[i][2],"%f", &weight_factor) == 1) ? 0 : 1;
	  if( weight_factor < 0 || weight_factor > 1 ) {
	    fprintf(stderr,"\nWeight factor must be in range [0,1].");
	    exit(1);
	  }
	  break;

        case 'M':
	  mflag = 1;
	  error = (sscanf(&argv[i][2],"%d", &maxpts2smo) == 1) ? 0 : 1;
	  break;

        case 'S':
	  sflag = 1;
	  break;

        case 'B':
	  bflag = 1;
	  if ((blankfile = fopen(&argv[i][2], "r")) == NULL) {
	    error = 1;
	    fprintf(stderr,"\nError opening %s\n", &argv[i][2]);
	  }
	  break;

        case 'h':
	  print_usage(argv[0]);
	  exit(0);

        default:
	  error = 1;
	  
      }    /* End switch */
      if (error ) {
	print_usage(argv[0]);
	fprintf(stderr,"\nError parsing command line arg:\n");
	fprintf(stderr,"      '%s'\n", argv[i]);
	exit(1);
      }
    }  /* end if */
    else  {
      infile = fopen(argv[i],"r");
      if (infile == NULL) {
	fprintf(stderr,"Unable to open %s for reading.\n", argv[i]);
	exit(1);
      }
      fprintf(stderr, "\nOpened %s ... \n", argv[i]);
    }  /* end else */
  }  /* end for */

  if (!(oflag &&  pflag)) {
       fprintf(stderr,"\nYou must specify -O and -P arguments.\n");
       exit(1);
  }

  /* Get row/column dimensions from infile */
  if (fscanf(infile,"%d%d%f%f%f%f%f%f%d", &nrows, &ncols, &xspacing, &yspacing, &xmin,
	     &ymin, &xmax, &ymax, &nprops )  != 9) {
    fprintf(stderr,"\nError reading heading from infile\n\n");
    exit(1);
  }

  fprintf(stderr, "\nInput grid is %d rows by %d cols", nrows, ncols);

  gridsize = nrows * ncols;

  fprintf(stderr,"\n Properties available: ");
  for (i = 0; i < nprops; ++i) {
    fscanf(infile,"%2s", id);
    prop_index[i] = get_prop_indx(id);
    prop_avail[prop_index[i]] = '1';
    fprintf(stderr," %2s", id);
  }
  fprintf(stderr,"\n");

  /* Check that requested properties are
   * available, allocate space for matrices,
   * and initialize values. */
  for (i = 0; i < MAXPROP; ++i ) {
    if (prop_request[i] && !prop_avail[i]) {
      fprintf(stderr,"\nRequested property: %2s not available.", get_prop_mne(i));
      prop_request[i] = 0;
    }

    /* Memory allocation for new grid. */
    new_matrix[i] = NULL;
    if (prop_request[i]) {
      new_matrix[i] = (struct new_grid_pt *) malloc(sizeof(struct new_grid_pt) * gridsize);
      if (new_matrix[i] == NULL) {
	fprintf(stderr,"\nUnable to allocate memory for new matrix.");
	exit(1);
      }

      /* For each requested property, zero out new grid. */
      for (sq= 0; sq < gridsize; ++sq) {
	newpt = &new_matrix[i][sq];
	newpt->sum = 0.;
	newpt->sumsq = 0.;
	newpt->n = 0;
	newpt->wghts = 0.;
      } /* End of loop over all grid points. */
    } /* End prop available check. */

    /* Memory allocation for old grid. */
    old_matrix[i] = NULL;
    if (prop_avail[i]) {
      old_matrix[i] = (struct old_grid_pt *) malloc(sizeof(struct old_grid_pt) * gridsize);
      if (old_matrix[i] == NULL) {
	fprintf(stderr,"\nUnable to allocate memory for old matrix.");
	exit(1);
      }

      /* For each available property, zero out new grid. */
      for (sq= 0; sq < gridsize; ++sq) {
	gridpt = &old_matrix[i][sq];
	gridpt->mean = 0.;
	gridpt->var = 0.;
	gridpt->n = 0;
      } /* End of loop over all grid points. */
    } /* End prop available check. */
  } /* End for loop over all possible properties. */

  /* Allocate blanking matrix. */
  blank = (char *) malloc(gridsize);
  if (blank == NULL) {
    fprintf(stderr,"\nUnable to allocate memory for blanking matrix.");
    exit(1);
  }

  /* Initialize blanking matrix. */
  for (sq = 0; sq < gridsize; ++sq) {
    blank[sq] = 0;
  }

  /* Read in blanking matrix. */ 
  if (bflag) {

    /* For each line in the file, */
    while (fscanf(blankfile,"%f %f", &lon, &lat) != EOF) {

      /* Compute the row and column indices of the point. */
      row = (int) (.0001 + (lat - ymin) / yspacing);
      col = (int) (.0001 + (lon - xmin) / xspacing);

      /* Check that the supplied points are within the bounds. */
      if ((row >= 0) && (row < nrows) && (col >= 0) && (col < ncols)) {

	/* Convert row, col indices to grid-vector index
	 * and assign all nodes to be blank with a 1. */
	blank[row * ncols + col] = '1';
      }
    } /* End of while loop over blankfile lines. */
  } /* End of blankfile present check. */

  /* Read data from infile and store in matrices
   * by row, col.  Multiply mean by nobs to produce
   * raw sums; square the std dev before multiplying
   * by (nobs - 1) and adding (sumx * sumx / n) to
   * produce the raw sumofsquares. */
  fprintf(stderr,"\nReading input file...\n");

  while ( fscanf(infile, "%f%f", &lat, &lon) != EOF ) {

    /* Convert lat and lon to matrix indices. */
    row = (int) (.0001 + (lat - ymin) / yspacing);
    col = (int) (.0001 + (lon - xmin) / xspacing);
    
    /* For each property available, */
    for (i = 0; i < nprops; ++i) {
      if (fscanf( infile,"%lf%lf%d", &mean, &stddev, &n) != 3) {
	fprintf(stderr,"\nError reading input file at ");
	fprintf(stderr,"values of (row, col) = ( %d , %d)\n", row, col);
      }
      
      /* Compute raw sums and raw sum of squares since:
       * (1) mean = propsum/n.
       * (2) mean = the only value for n = 1 case
       * (3) stddev = sqrt( (mean - prop)/(n-1) )*/
      if (n > 0) {                   
	propsum = mean * n;          /* (1) */
	propsqsum = mean * mean;     /* (2) */
	if (n > 1) {                 /* raw sum of squares */
	  propsqsum = stddev * stddev * (n-1) + mean * mean * n ;
	}
      }
      
      /* Find the menomic index. */
      j = prop_index[i];

      /* Find index to old_matrix. */
      sq = row * ncols + col;

      /* Get addr of gridpoint. */
      gridpt = &old_matrix[j][sq];

      /* Store info at that addr. */
      gridpt->mean = propsum;
      gridpt->var = propsqsum;
      gridpt->n = n;
    }

  }  /* End while looping over all lines in input file. */

  /*  For each pt in the new grid, sum the values
   * from surrounding squares in old grid, weighted
   * by distance. Move successively outward until 
   * search radius is reached.  Compute means and
   * standard deviations from weighted sums. */
  fprintf(stderr,"\nComputing new matrix ...\n");
  fname = (char *)malloc(50);
  free(id);

  /* Loop over each possible property. */
  for (i = 0; i < MAXPROP; ++i) {

    /* Check for requested (and also available) properties. */
    if (prop_request[i]) {

      /* Get prop info. */
      id = get_prop_mne(i);    
      sprintf(field_descrip, " %c%d.%dlf", '%',get_field_width(i), get_field_precis(i)); 
      fprintf(stderr,"\n%s ",id);

      /* Open output file for prop.
       * <root><property_id>.xyz */
      fname = strcpy(fname,root);
      fname = strncat(fname,id,2);
      fname = strncat(fname,".xyz",4);
      outfile[0] = fopen(fname,"w");

      /* Open blanking file for prop.
       * <root><property_id>.blk */
      fname = strcpy(fname,root);
      fname = strncat(fname,id,2);
      fname = strncat(fname,".blk",4);
      outfile[1] = fopen(fname,"w");       

      if (sflag) {
	/* Open stddev file for prop.
	 * <root><property_id>_var.xyz */
	fname = strcpy(fname,root);
	fname = strncat(fname,id,2);
	fname = strncat(fname,"_var.xyz",8);
	outfile[2] = fopen(fname,"w");

	/* Open stddev blanking file.
	 * <root><property_id>_var.blk */	
	fname = strcpy(fname,root);
	fname = strncat(fname,id,2);
	fname = strncat(fname,"_var.blk",8);
	outfile[3] = fopen(fname,"w");
      } /* End of sflag check. */

      /* Initialize the gridpoint counter. */
      sq = -1;

      /* Cycle through all grid points. */
      for (row=0; row < nrows; ++row) {
	fprintf(stderr,"*");
	for (col=0; col < ncols; ++col) {
	  ++sq;

	  /* Lats and lons are automatically set
	   * to the centers between row and column
	   * lines (Pixel Registration). */
	  lat = ymin + (row + .5) * yspacing;
	  lon = xmin + (col + .5) * xspacing;
	  
	  /* Test whether to blank this gridpt.
	   * and if so, write blank coordinates. */
	  if (blank[sq]) {
	    fprintf(outfile[1],"%.3f %.3f\n", lon, lat);
	    if (sflag) {
	      fprintf(outfile[3],"%.3f %.3f\n", lon, lat);
	    }                  
	  }
	  else { /* Do not blank this gridpt.  Apply the smoother! */
	    /* Initialise the search radius. */
	    rad = 0;

	    /* Initialize new point with old matrix values.
	     * New, init, #points = #points in core/center point.
	     * New, init, sum of values = "mean" = propsum
             * New, init, sum sq values = "var" = propsqsum
	     * New, init, weights = #points in core/center point. */
	    newpt = &new_matrix[i][sq];
	    gridpt = &old_matrix[i][sq];
	    newpt->n = gridpt->n;
	    newpt->sum = gridpt->mean;
	    newpt->sumsq = gridpt->var;
	    newpt->wghts = gridpt->n;
	    
	    /* Slowly increase the size of the search
	     * radius and add information to the running
	     * sums in this grid point from surrounding
	     * points. This will *not* engage if maxrad
	     * stays at its default value of 0 or if we
	     * have more points than maxpts2smo. */
	    while  ( (rad < maxrad) && (newpt->n <= maxpts2smo) ) {
	      ++rad;
	      add_wghted_layer(rad, weight_factor, newpt, row, col,
			       old_matrix[i], nrows, ncols);
	    }
	    
	    /* Write to blanking files if we still have no
	     * contributing values to this node even after
	     * going up to the full search radius/smoothing. */
	    if ((n = newpt->n) == 0) {
	      fprintf(outfile[1],"%.3f %.3f\n", lon, lat);
	      if (sflag) {
	        fprintf(outfile[3],"%.3f %.3f\n", lon, lat);
	      }
	    }
	    else {
	      /* Write valid data to property file. */
	      mean = newpt->sum / newpt->wghts;
	      fprintf(outfile[0],"%.3f %.3f ", lon, lat); 
	      fprintf(outfile[0], field_descrip, mean);
	      if (nflag) {
		fprintf(outfile[0]," %4d",  n);
	      }
	      fprintf(outfile[0],"\n");
	      
	      if (sflag) {  
		/* Write to stddev file */       
		stddev = 0.0;
		if (n > 1) {
		  var = (newpt->sumsq - newpt->sum * newpt->sum /
			 (double) newpt->wghts ) / (double) (newpt->wghts - 1);
		  stddev = sqrt(ABS(var));
		  fprintf(outfile[2],"%.3f %.3f %.7lf", lon, lat, stddev);
		  if (nflag) {
		    fprintf(outfile[2]," %4d",  n);
		  }
		  fprintf(outfile[2],"\n");
		}
		else {         /* Write to stddev blanking file */
		  fprintf(outfile[3],"%.3f %.3f\n", lon, lat);
		}
	      } /* End standard deviation output check. */
	    } /* End else (for grid points with non-zero nobs). */
	  } /* End else (non-blanked grid points). */
	} /* End for loop over all columns. */
      } /* End for loop over all rows. */

      /* Close all the outfiles associated
       * with this property. */
      fclose(outfile[0]);
      fclose(outfile[1]);
      if (sflag) {         
	fclose(outfile[2]);
	fclose(outfile[3]);
      }

    } /* End check if property available & requested. */
  }  /* End for loop over all possible properties. */
  
  fprintf(stderr,"\nEnd of %s.\n", argv[0]);
  exit(0);
} /* End main */

/****************************************************************************/
void print_usage(char *program) {
  fprintf(stderr,"\n Separates a HydroBase *.grid file (output by");
  fprintf(stderr,"\n hb_gridsurf) into individual *.xyz files");
  fprintf(stderr,"\n containing lon, lat, property triplets.");
  fprintf(stderr,"\n Optionally smooths each point with information");
  fprintf(stderr,"\n from surrounding gridpoints (radius specified");
  fprintf(stderr,"\n with -R) weighted by a distance factor = W/d");
  fprintf(stderr,"\n (W specifed by -W) and where d = # of gridpts");
  fprintf(stderr,"\n away from the central node.  Also, smoothing" );
  fprintf(stderr,"\n can be selectively applied to gridpoints with");
  fprintf(stderr,"\n a small number of points (max number of points");
  fprintf(stderr,"\n to be smoothed is specified by -M).");
  fprintf(stderr,"\n");
  fprintf(stderr,"\n The program outputs up to 4 files for each");
  fprintf(stderr,"\n property specified:"); 
  fprintf(stderr,"\n   1) property info:   lon, lat, av_value [, n]");
  fprintf(stderr,"\n   2) blank gridpts:   lon, lat");
  fprintf(stderr,"\n   3) stddev info  :   lon, lat, stddev [, n]");
  fprintf(stderr,"\n   4) stddev blanks:   lon, lat");
  fprintf(stderr,"\n Note that stddev and number of points info is");
  fprintf(stderr,"\n output only if the -S and -N options, respectively,");
  fprintf(stderr,"\n are specified.");
  fprintf(stderr,"\n");
  fprintf(stderr,"\n Output files are named: <root><property_id>.xyz");
  fprintf(stderr,"\n                         <root><property_id>.blk");
  fprintf(stderr,"\n                         <root><property_id>_var.xyz");
  fprintf(stderr,"\n                         <root><property_id>_var.blk");
  fprintf(stderr,"\n where root is supplied by the -O option and");
  fprintf(stderr,"\n property_id is the same as the property_list_ids");
  fprintf(stderr,"\n in the -P option.");
  fprintf(stderr,"\n");
  fprintf(stderr,"\n The -B option permits you to supply a file");
  fprintf(stderr,"\n containing the lon,lat of gridpoints to be blanked.");
  fprintf(stderr,"\n Its purpose is to mask out bathymetry, coastlines,");
  fprintf(stderr,"\n etc.  This list is checked before computing each new");
  fprintf(stderr,"\n gridpoint.");
  fprintf(stderr,"\n");
  fprintf(stderr,"\nUsage:  %s input_filename -O<output_file_root>",program);
  fprintf(stderr,"\n                          -P<property_list>");
  fprintf(stderr,"\n                         [-R<smoothing_radius>]");
  fprintf(stderr,"\n                         [-W<weighting_factor>]");
  fprintf(stderr,"\n                         [-M<max_pts_2_smooth>]");
  fprintf(stderr,"\n                         [-N] [-S] [-B<blanking_info_file> [-h]");
  fprintf(stderr,"\n\n   specify input file_name (HydroBase *.grid file ");
  fprintf(stderr,"\n   output by hb_gridsurf) as first argument");
  fprintf(stderr,"\n    -O  : specifies root name for output files");
  fprintf(stderr,"\n          ex: -Opr200.1970_75.");
  fprintf(stderr,"\n    -P  : specifies properties to project onto surface");
  fprintf(stderr,"\n          Enter a list of your choices (separated by slashes)...");
  fprintf(stderr,"\n     ex: -Ppr/th/sa/ox/de/ht");
  fprintf(stderr,"\n         -P (by itself) will print a menu of available properties");
  fprintf(stderr,"\n   [-R] : smoothing radius (# of gridpts away from central node)");
  fprintf(stderr,"\n            ex: -R0 provides no smoothing [default]"); 
  fprintf(stderr,"\n                -R1 gives 3*gridincr smoothing");
  fprintf(stderr,"\n                -R2 gives 5*gridincr smoothing");
  fprintf(stderr,"\n   [-W] : weight factor [0,1] for surrounding points when smoothing.");
  fprintf(stderr,"\n            ex: -W0.5 gives (default) weighting of 0.5/d.");
  fprintf(stderr,"\n   [-M] : max number of points to smooth selectively limits smoothing");
  fprintf(stderr,"\n          to regions of low data density.  If not specified, all grid");
  fprintf(stderr,"\n          nodes are smoothed equally [default].");
  fprintf(stderr,"\n            ex: -M9 applies smoothing to regions with 9 and fewer");
  fprintf(stderr,"\n                contributing measurements, as determined by hb_gridsurf2d.");
  fprintf(stderr,"\n   [-N] : output number of obs incorporated into each gridpt.");
  fprintf(stderr,"\n   [-S] : output std dev file for each property.");
  fprintf(stderr,"\n   [-B] : specifies name of file with lon/lat of pts");
  fprintf(stderr,"\n          to be blanked.  This list is checked BEFORE");
  fprintf(stderr,"\n          computing each new gridpoint.");
  fprintf(stderr,"\n   [-h] : help -- prints this message");
  fprintf(stderr,"\n\n");  
  return;
}

/*******************************************************************************/
void add_wghted_layer(int radius, double wfactor, struct new_grid_pt *gridpt, int i, int j, struct old_grid_pt *oldprop_ptr, int nrows, int ncols) {

  /* add_wghted_layer() identifies pts in
   * an existing matrix (size: nrows * ncols)
   * that are radius # of points away from a
   * point centered at (i,j)  The values of
   * each property at each point thus
   * identified, are weighted by (1/radius)*
   * wfactor, and agglomerated into the
   * appropriate gridpt in the new matrix.
   * The same is done for the squared sums
   * and the number of pts.
   *
   * NOTE: This routine just touches on the
   * "rim" of grid points defined by the radius
   * and does not include points on the interior
   * (if we have a large radius).  Thus, if we
   * have a radius larger than 1, this routine
   * needs be run for all interior radii as well
   * as the outer radius.
   *
   * radius    smoothing radius (in # of gridpts )
   * wfactor   base value for weighting factor
   * gridpt    address of current gridpt -> target for calculations. 
   * i, j      row, col of currrent gridpt in matrix 
   * oldprop_ptr   starting addr of prop in old matrix 
   * nrows, ncols  dimensions of  matrix (same for both old and new).
   */

  int top, bot, left, right;
  int row, col, layer, side;
  float  weight;
  struct old_grid_pt *oldgridpt;
  
  /* Return if we have a zero search radius. */
  if (radius <= 0) {
    return ;
  }
  
  /* Define weight for values in this layer */
  weight = wfactor/(2.0 * radius);

  /* Define row,col of lower left corner of first layer */
  row = i;
  col = j;
  layer = 1;
  
  /* Move to lower left corner of layer corresponding to radius */
  while (layer++ < radius) {   
    --row;
    --col;
  }  /* End while */
  
  /* Length of a side of the square defined
   * by locus of pts in layer within the radius. */
  side = 2 * radius + 1; 
  
  /* Define corners of square
   *         (right,top)
   *     +-------+
   *     |       |
   *     |   +   |
   *     |       |
   *     +-------+
   * (left,bot)
   *
   * Need to subtract 1 when adding by
   * side number of indices, ex. (1,1)
   * and (9,5) for a box 9x5.
   */
  bot = row;
  left = col;
  top = row + side - 1;
  right = col + side - 1;

  /* Move up left side of square by first
   * testing are we within column bounds
   * of the matrix? */
  if ((col >= 0) && (col < ncols)) {
    while ( row < top) {

      /* Test are we within row bounds? */
      if ((row >= 0) && (row < nrows)) {

	/* Retrieve the old grid information. */
	oldgridpt = &oldprop_ptr[row * ncols + col];

	/* Test if there are any obs at this pt. */
	if (oldgridpt->n > 0) {

	  /* If there are observations, add values
	   * to the running sums. */
	  gridpt->sum += weight * oldgridpt->mean;
	  gridpt->sumsq += weight * oldgridpt->var;
	  gridpt->n += oldgridpt->n;
	  gridpt->wghts += weight * oldgridpt->n;
	}
      }
      ++row;
    } /* end while */
  }  /* end if */
  
  /* Set row, col to point to top, left corner */
  row = top;
  col = left;
  
  /* Move across top of square */
  if ((row >= 0) && (row < nrows)) {     /* within row bounds? */
    while ( col < right) {
      if ((col >= 0) && (col < ncols)) {      /* within column bounds? */
	oldgridpt = &oldprop_ptr[row * ncols + col];
	if (oldgridpt->n > 0) {             /* any obs at this pt? */
	  gridpt->sum += weight * oldgridpt->mean;
	  gridpt->sumsq += weight * oldgridpt->var;
	  gridpt->n += oldgridpt->n;
	  gridpt->wghts += weight * oldgridpt->n;
	}
      }
      ++col;
    } /* end while */
  }  /* end if */
  
  /* set row, col to point to top, right corner */  
  row = top;
  col = right;
  
  /*    move down right side of square */
  
  if ((col >= 0) && (col < ncols)) {     /* within col bounds? */
    while ( row > bot) {
      if ((row >= 0) && (row < nrows)) {      /* within row bounds? */
	oldgridpt = &oldprop_ptr[row * ncols + col];
	if (oldgridpt->n > 0) {             /* any obs at this pt? */
	  gridpt->sum += weight * oldgridpt->mean;
	  gridpt->sumsq += weight * oldgridpt->var;
	  gridpt->n += oldgridpt->n;
	  gridpt->wghts += weight * oldgridpt->n;
	}
      }
      --row;
    } /* end while */
  }  /* end if */
  
  /* set row, col to point to bottom, right corner */
  
  row = bot;
  col = right;
  
  /*    move across bottom of square */
  
  if ((row >= 0) && (row < nrows)) {     /* within row bounds? */
    while ( col > left) {
      if ((col >= 0) && (col < ncols)) {      /* within col bounds? */
	oldgridpt = &oldprop_ptr[row * ncols + col];
	if (oldgridpt->n > 0) {             /* any obs at this pt? */
	  gridpt->sum += weight * oldgridpt->mean;
	  gridpt->sumsq += weight * oldgridpt->var;
	  gridpt->n += oldgridpt->n;
	  gridpt->wghts += weight * oldgridpt->n;
	}
      }
      --col;
    } /* end while */
  }  /* end if */
  
  return;
} /* End add_wghted_layer() */
/*****************************************************************/

