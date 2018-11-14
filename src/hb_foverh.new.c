/* hb_foverh.c
................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                    author:  Ruth Curry
                             Woods Hole Oceanographic Institution
                             May 2001
................................................................................
..................................................................  
.  Reads binary gridded bathymetry values consisting of 2-byte integers
.  and writes out lon lat f/h triplets

..................................................................  
*/ 
#include <stdio.h> 
#include <math.h>
#include <string.h>
#include "hydrobase.h"
#include "hb_memory.h"
#include "hb_paths.h"

#define   NINT(x)       (((x) < 0) ? ((int)((x) - 0.5)) : ((int)((x) + 0.5)))
#define    ABS(x)       (((x) < 0) ? -(x) : (x))

struct GRID_INFO {
	int nx;			/* Number of columns */
	int ny;			/* Number of rows */
	int node_offset;	/* 0 for node grids, 1 for pixel grids */
	double x_min;		/* Minimum x coordinate */
	double x_max;		/* Maximum x coordinate */
	double y_min;		/* Minimum y coordinate */
	double y_max;		/* Maximum y coordinate */
	double x_inc;		/* x increment */
	double y_inc;		/* y increment */

};

/* globally defined variables */


double xmin = 0.0;        /* input grid bounds assumed */
double xmax = 360.0;
double ymin = -90.0;
double ymax = 90.0;
double xincr, yincr;
int  depth_is_neg; 
int nrowsin, ncolsin;  

/*  prototypes for locally defined functions */

void print_usage(char *);
short **get_topo(FILE *);

           
int main (int argc, char **argv)
{ 
   int i, j, k, row, col;
   int i_flag, b_flag, t_flag;
   int nrowsout, ncolsout;
   int startcol, startcol2;
   int endcol, endcol2;
   int startrow, endrow;
   int error, latfirst, neglon;
   int lon0to360, depth_is_neg, xgreenwich;   
   int dummy[4] = {0,0,0,0};
   short int **zin;
   float zout, lat, lon, tmp;
   double fcor;
   char *s;
   FILE *outfile;
   FILE *topofile;
   FILE *ufile, *lfile;
   struct GRID_INFO h;

   if (argc < 1 ){
      print_usage(argv[0]);
   }
   
/* initialize these ... */

   h.x_min = xmin;
   h.x_max = xmax;
   h.y_min = ymin;
   h.y_max = ymax;
   t_flag = 0;
   i_flag = 0;
   b_flag = 0;
   latfirst = 0;
   startcol = endcol = startcol2 = endcol2 = 0;
   startrow = endrow = 0;
   xincr = yincr = 0.1;
   depth_is_neg = 1; 
   error = 0;  
   outfile = stdout;
   ufile = lfile = topofile = (FILE *)NULL;

   
/* parse the command line arguments */

   for (i = 1; i < argc; i++) { 
      if (argv[i][0] == '-') {
         s = &argv[i][1]; 
         switch (*s) { 
            case 'B':                    /* get output grid bounds */
                s = &argv[i][2];
                if (*s == '/')
                      ++s;
                error = (sscanf(s,"%lf", &h.x_min) != 1);
                while (*(s++) != '/')
                    ;  
                error += (sscanf(s,"%lf", &h.x_max) != 1);
                while (*(s++) != '/')
                   ;  
                error += (sscanf(s,"%lf", &h.y_min) != 1);
                while (*(s++) != '/')
                   ;  
                error += (sscanf(s,"%lf", &h.y_max) != 1);

                if (h.x_min > h.x_max)  {
                  fprintf(stderr,"\nW bound cannot exceed E bound.\n");
                  error = 1;
                }
                                               
                if (h.y_min > h.y_max) { 
                  fprintf(stderr,"\nS bound cannot exceed N bound.\n");
                  error = 1;
                }
                b_flag = 1;          
                break;

            case 'I':
               i_flag = 1;
               ++s;
               if (*s == '/')
                  ++s; 
               error += (sscanf(s,"%lf", &xincr) != 1); 
               s = strchr(s,'/'); /* check for another delimiter*/ 
               yincr = xincr; 
               if (s != NULL) { 
                  ++s; /* move past delimiter */ 
                  error += (sscanf(s,"%lf", &yincr) != 1); 
               } 
               break;
               
           case 'L':
               lfile = fopen(&argv[i][2],"r");
               if (lfile == NULL) {
                  fprintf(stderr,"\nError opening %s for input\n", &argv[i][2]);
                  exit(1);
               }
               break;
           case 'O':
               outfile = fopen(&argv[i][2],"w");
               if (outfile == NULL) {
                  fprintf(stderr,"\nError opening %s for output\n", &argv[i][2]);
                  exit(1);
               }
               break;
               
           case 'T':
               topofile = fopen(&argv[i][2],"r");
               if (topofile == NULL) {
                  fprintf(stderr,"\nError opening %s for input\n", &argv[i][2]);
                  exit(1);
               }
               t_flag = 1;
               break;
               
           case 'U':
               ufile = fopen(&argv[i][2],"r");
               if (ufile == NULL) {
                  fprintf(stderr,"\nError opening %s for input\n", &argv[i][2]);
                  exit(1);
               }
               break;
	       
             case 'h': 
               print_usage(argv[0]);
               exit(0);
               
            case ':':
               latfirst = 1;
               break;
	       
            default:  
               error = 1;
               
         } /* end switch */
             
      }
      else  {
        error = 1;
      }
      
      if (error ) { 
           fprintf(stderr,"\nError parsing command line args.\n");
           fprintf(stderr," in particular:  '%s'\n", argv[i]); 
           fprintf(stderr,"\nType %s -h for help\n", argv[0]); 
           exit(1); 
      }
   } /* end for */
   
    if (!ufile && !b_flag) {
          fprintf(stderr,"\nMust specify either -B OR -U arguments!\n");
          fprintf(stderr,"\nType %s -h for help\n", argv[0]); 
          exit(1);
    }

    if (ufile && lfile && !b_flag) {
          fprintf(stderr,"\nMust specify -B if using both -U and -L arguments!\n");
          fprintf(stderr,"\nType %s -h for help\n", argv[0]); 
          exit(1);
    }
    
    if (lfile && !ufile) {
        fprintf(stderr,"\nMust specify -U if using -L argument!\n");
        fprintf(stderr,"\nType %s -h for help\n", argv[0]); 
        exit(1);
    }
    
    if ( !lfile && !t_flag) {
        topofile = fopen(BATHPATH,"r");
        if (topofile == NULL) {
          fprintf(stderr,"\nUnable to open default topography file\n", BATHPATH);
          fprintf(stderr,"You can specify the name of a binary topography file for input with [-T]\n"); 
          exit(1);
        }
   } 
   if (topofile && !i_flag) {
      fprintf(stderr,"Increment for topofile: [%.2lf/%.2lf] \n", xincr, yincr);
   } 
   

   if (topofile) {
      nrowsin = (int) (NINT((ymax - ymin) / yincr)) + 1;
      ncolsin = (int) (NINT((xmax - xmin)/ xincr)) + 1; 
   
      /* read in topography file ... */

       zin = get_topo(topofile);
       fclose(topofile);
    }
    else {
      nrowsin = (int) (NINT((h.y_max - h.y_min) / h.y_inc)) + 1;
      ncolsin = (int) (NINT((h.x_max - h.x_min)/ h.x_inc)) + 1; 
    
    }

    
    
/*  BRANCH to output values by gridbounds */

  if (!ufile) {
      startrow = (int) (NINT(( h.y_min - ymin) / yincr) + 0.0001);
      endrow = (int) (NINT(( h.y_max - ymin) / yincr) + 0.0001);
   

     /* find starting/ending columns.  Allow for crossing Greenwich meridian */
   
      lon0to360 = 1;   
      if (h.x_min < 0)
         lon0to360 = 0;
                             
      xgreenwich = 0;
      if (!lon0to360)
         if (h.x_max >= 0) 
             xgreenwich = 1;
      
   
      if (lon0to360) {   /* the simple case where lon bounds are all positive */
         startcol = (int)( NINT(h.x_min / xincr) + 0.0001);
         endcol = (int)( NINT(h.x_max  / xincr) + 0.0001);
      }
   
      if (!lon0to360) {  /* requested at least one negative longitude */
   
        startcol = (int) (NINT((h.x_min + 360.) / xincr) + 0.0001);
        endcol = (int) (NINT((h.x_max + 360.) / xincr) + 0.0001);
     
        if (xgreenwich) {
           endcol = ncolsin - 1;  /* include greenwich */
           startcol2 = 1;     /* don't repeat greenwich */
           endcol2 = (int) (NINT(h.x_max / xincr) + .0001);
        }
      }
   
    /* adjust output grid bounds to mirror input gridnodes */

      if (h.x_min >= 0) 
          h.x_min = startcol * xincr + xmin;
      else
          h.x_min = startcol * xincr + xmin - 360.0;
      
      if (xgreenwich)
           h.x_max = endcol2 * xincr + xmin;
      else if (h.x_max >= 0)
           h.x_max = endcol * xincr + xmin;
      else     
           h.x_max = endcol * xincr + xmin - 360.0;
       
      h.y_min = startrow * yincr + ymin;
      h.y_max = endrow * yincr + ymin;

      fprintf(stderr,"\nOutput grid bounds: %8.1lf/%8.1lf/%8.1lf/%8.1lf", h.x_min, h.x_max, h.y_min, h.y_max); 
      fprintf(stderr,"\nGrid increment: %6.3lf/%6.3lf\n", xincr, yincr); 
  
  
       h.x_inc = xincr;
       h.y_inc = yincr;
       h.node_offset = 0;
      
      /* allocate memory and initialize arrays...  */
  
   
       nrowsout = endrow - startrow + 1;  
       ncolsout = (int) (NINT((h.x_max - h.x_min) / h.x_inc) + 1); 
   
      fprintf(stderr,"\n nrows: %4d   ncols: %4d\n", nrowsout, ncolsout); 
   
    /* subsample the input array and output values ... */
      
    fprintf(stderr,"\nWriting xyz file with seafloor as positive values.\n");
    
    for (row = startrow; row <= endrow; ++row) {
    
      for (col = startcol; col <= endcol; ++col) {
      
        zout = (float) zin[row][col];
        
        if (depth_is_neg)
           zout = -zout;
	   
	lat = ymin + row * yincr; 
	lon = xmin + col * xincr;
	if (!lon0to360 && lon > 180)
	    lon -=  360.;
	    
	fcor = hb_coriol((double)lat);
        if (lat < 0) fcor = -fcor;
	zout = (float) (fcor /(double) zout);
	   
	fprintf(outfile,"%8.3f %8.3f %.4e\n", lon, lat, zout);
      }
      
      if (startcol2 > 0) {
        for (col = startcol2; col <= endcol2; ++col) {
           zout = (float) zin[row][col];
           if (depth_is_neg)
              zout = -zout;
	   lat = ymin + row * yincr; 
	   lon = xmin + col * xincr;
	   if (!lon0to360 && lon > 180)
	       lon -= - 360.;
	   fcor = hb_coriol((double)lat);
	   zout = (float) (fcor /(double) zout);
	   fprintf(outfile,"%8.3f %8.3f %.4e\n", lon, lat, zout);
        }
      }
    } 
         
    fprintf(stderr,"\n End of %s.\n", argv[0]);
    exit(0);
 
 }  /* end if !ufile */
 
/*  BRANCH to use input file and subtract depths */
 
 while ((error = fscanf(ufile, "%f %f %f", &lon, &lat, &zout)) == 3) {
    if (latfirst) {
       tmp = lat;
       lat = lon;
       lon = tmp;
    }

    /* find index to topo grid from lat/lon */
    
    row = NINT(((double)lat - ymin ) / yincr);
    if (row >= nrowsin) 
          fprintf(stderr,"\nError computing gridnode for latitude %.3f\n", lat);
    
    neglon = 0;
    if (lon < 0) {
        lon += 360.0;
	neglon = 1; 
    }
    
    col = NINT(((double)lon - xmin) / xincr);
    
    if (col >= ncolsin) 
          fprintf(stderr,"\nError computing gridnode for longitude %.3f\n", lon);


    fcor = hb_coriol((double)lat);
    if (lat < 0) fcor = -fcor;

    tmp = (float) zin[row][col];
    if (depth_is_neg)
       tmp = -tmp;
     
    if (neglon) lon = lon - 360. ; 
    fprintf(outfile,"%8.3f %8.3f %.4e\n", lon, lat, (fcor / (double)(tmp - zout)));
        

 }  /* end while */
 
 if (error != EOF){
   fprintf(stderr,"\nError reading H file \n");
   exit(1);
 }
 
 fprintf(stderr,"\n End of %s.\n", argv[0]);
 exit(0);

} /* end main */


/****************************************************************************/

void print_usage(char *program) 
{ 
   fprintf(stderr,"\n%s reads gridded bathymetry values from a binary", program);
   fprintf(stderr,"\ntopography file and computes f/H values.");
   fprintf(stderr,"\nThe user must specify the grid increment of");
   fprintf(stderr,"\nthe binary file with -I.  Use -B to specify the");
   fprintf(stderr,"\nbounds of the area to extract bathymetry values."); 
   fprintf(stderr,"\nLon/lat/z values are written to outfile or stdout.");
   fprintf(stderr,"\n\nUsage:  %s  -T<topo.onetenthdeg.dat> [-B<w/e/s/n>] [-I<xincr[/yincr]>] [-N] [-O<outfile>]", program);
   fprintf(stderr,"\n\n   OPTIONS:"); 
   fprintf(stderr,"\n[-B] : specifies output bounds, if not 0/360/-90/90");
   fprintf(stderr,"\n[-H] : file defining lon/lat/depth of upper surface");
   fprintf(stderr,"\n        to be used.  This depth will be subtracted");     fprintf(stderr,"\n        from seafloor depth to give H."); 

   fprintf(stderr,"\n[-I] : specify grid increments of input topo file. "); 
   fprintf(stderr,"\n        ex: -I.1/.5    default increment is .1 "); 
   fprintf(stderr,"\n[-O] : output file. If not specified output stdout"); 
   fprintf(stderr,"\n[-T]: specify full pathname of binary topography file.  ");
   fprintf(stderr,"\n       Default is [%s]", BATHPATH);
   fprintf(stderr,"\n[-h] : help -- prints this message.");
   fprintf(stderr,"\n\n");
   return;
   
} /*end print_usage() */
   
/****************************************************************************/
/****************************************************************************/

short **get_topo(FILE *fptr)
   /* allocates memory and reads in topography values.  Returns a pointer to
      the start of the memory block. An error causes an error message to be
      printed and an exit.   */
{

   int row, n;
   short **z;
   
   /* Allocate memory */
   
   z = (short **) malloc(nrowsin * sizeof(short *));
   if (z == NULL) {
      fprintf(stderr,"\nError allocating memory.\n");
      exit(1);
   }
   
   for (row = 0; row < nrowsin; ++row) {
     z[row] = (short *) malloc(ncolsin * sizeof(short));
     if (z[row] == NULL) {
         fprintf(stderr,"\nError allocating memory.\n");
         exit(1);
      }
   }
   
   fprintf(stderr,"\nReading in topography values ");
   
   for (row = 0; row < nrowsin; ++row) {
     n = fread((short *)z[row], sizeof(short), ncolsin, fptr);
     if (n != ncolsin) {
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
