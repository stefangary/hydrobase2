/* hb_fit.h
................................................................................
                          *******  HydroBase2  *******
               ...................................................
                          
                    author:  Ruth  Curry
                             Woods Hole Oceanographic Institution
                             2000
................................................................................
*/


#ifndef HB_FIT_PARMS
#define HB_FIT_PARMS 1


#ifndef TRUE
#define TRUE  1
#endif

#ifndef  FALSE
#define  FALSE 0
#endif

#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))     /* min and max value macros */
#endif


#ifndef TOOBIG
#define TOOBIG 0.9e+35     /* value larger than any oceanographic property */
#endif                     /* but small enough to fit in a float variable */

#define ZGRID_EMPTY 1.0e+35
#define MEM_CHUNK  2000
#define SSIZE       500    /* arbitrary work array size */


 
#ifndef NINT
#define   NINT(x)       (((x) < 0) ? ((int)((x) - 0.5)) : ((int)((x) + 0.5)))
#endif

#ifndef ABS
#define    ABS(x)       (((x) < 0) ? -(x) : (x))
#endif

typedef int BOOLEAN;                    /* used for logical variables */

struct GRID_INFO {
	int nx;				/* Number of columns */
	int ny;				/* Number of rows */
	int node_offset;		/* 0 for node grids, 1 for pixel grids */
	double x_min;			/* Minimum x coordinate */
	double x_max;			/* Maximum x coordinate */
	double y_min;			/* Minimum y coordinate */
	double y_max;			/* Maximum y coordinate */
	double x_inc;			/* x increment */
	double y_inc;			/* y increment */

};

struct POINT {               /* Stores input data */
  double XP, YP, ZP;
} ;


extern int ij2xy(struct GRID_INFO *, int, int, float *, float *);
extern int xy2ij(struct GRID_INFO *, float, float, int *, int *);
extern void interp2d(double *, double *,  double, double, int , int , struct GRID_INFO *, double *, double *, double *, int *);
extern double * zgrid (double, int, int, int, int, double,  struct GRID_INFO *, int *, int *, double *, double *, struct POINT *, int *, int *, int *);
extern int inside (double, double, double *, double *, int);
extern int get_mask(FILE *, struct GRID_INFO *, char *, BOOLEAN);

#endif  /* ifndef HB_FIT_PARMS */


/*
-----------------------------------------------------------------------------------------
 	Notes on node_offset:

	Assume x_min = y_min = 0 and x_max = y_max = 10 and x_inc = y_inc = 1.
	For a normal node grid we have:
		(1) nx = (x_max - x_min) / x_inc + 1 = 11
		    ny = (y_max - y_min) / y_inc + 1 = 1
		(2) node # 0 is at (x,y) = (x_min, y_max) = (0,10) and represents the surface
		    value in a box with dimensions (1,1) centered on the node.
	For a pixel grid we have:
		(1) nx = (x_max - x_min) / x_inc = 10
		    ny = (y_max - y_min) / y_inc = 10
		(2) node # 0 is at (x,y) = (x_min + 0.5*x_inc, y_max - 0.5*y_inc) = (0.5, 9.5)
		    and represents the surface value in a box with dimensions (1,1)
		    centered on the node.
-------------------------------------------------------------------------------------------*/
