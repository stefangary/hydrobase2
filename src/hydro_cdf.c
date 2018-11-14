/* hydro_cdf.c
................................................................................
                          *******  HydroBase2  *******
               ...................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
                             2000
................................................................................
.
.  A set of routines to facilitate the creation and reading of netCDF files 
.  in conjunction with HydroBase2 software.  These functions utilize the UCAR
.  netCDF C-interface routines version 3+.
.

.  int cdf_init(char *filename)    
	Opens a new CDF file for writing.   

.  int cdf_open(char *dir, char *file_root, char *extent, int print_mess)
           Open a cdf file for reading.

.  int cdf_update(char *dir, char *file_root, char *extent, int print_mess)
           Open a cdf file for reading and writing.
	   
.  void cdf_close(int cdfid)
         Close a cdf file.

.  int cdf_define(int cdfid, struct CDF_HDR *hptr, int prefill, int add_count_vars)
        Defines all variables and attributes for cdf file.

.  int write_prop_cdf(int cdfid, float *data_ptr, char *prop_mne, int row, int col, int tbin, int zlev, int nrows, int ncols, int ntbins, int nlevs)
     Writes values of a property at a lat/lon node.

.  int write_prop_count_cdf(int cdfid, short *cnt_ptr, char *prop_mne, int row, int col,int tbin, int zlev, int nrows, int ncols, int ntbins, int nlevs)
	Writes # of obs for a property at a lat/lon node.

.  int write_std_depths_cdf(int cdfid, struct CDF_HDR *hptr)
	 Writes standard depths to de variable.

.  int write_time_bins_cdf(cdfid, hptr)
	Writes min/max values defining time bins.

.  int write_bottom_depth_cdf(int cdfid, int row, int col, int tbin, int nrows, int ncols, int ntbins, float *data_ptr)
	Writes the bottom depth values to bottom variable.

.  int read_cdf_hdr(int cdfid, struct CDF_HDR *haddr)
      Reads header info from a cdf file.

.  int read_cdf_prop(int cdfid, char *varname, float *dataptr, int row, int col, int tbin, int zlev, int npts)
      Reads all values of a property at a lat/lon node.

.  int read_cdf_prop_count(int cdfid, char *prop_mne, short *cnt_ptr, int row, int col, int tbin, int zlev, int npts)
	 Reads # of obs for a property at a lat/lon node.

.  int read_cdf_depths(int cdfid, float *d_ptr)
	Reads values of depth variable from cdf file.

.  int read_cdf_bottom_depth(int cdfid, float *d_ptr, int row, int col, int tbin)
	 Reads values of bottom depth from cdf file.

.  int get_indices(struct CDF_HDR *hptr, float lat, float lon, int *row_ptr, int *col_ptr)
       Converts lat/lon to row/col for a cdf file.

.  int get_lat_lon(struct CDF_HDR *hptr, int row, int col, float *lat_ptr, float *lon_ptr)
       Converts row/col to lat/lon.

.  int std_depth_init(FILE *depth_file)
     Initializes globally defined std_depth array.

.  int d2stdlev(double depth)
  	returns stddepth index associated with depth

.  double stdlev2depth(int index)
	 returns depth associated with stddepth index.
  
.  int is_flagged(float x, float flag)
        Determines whether x is a missing value based on the value of flag.

 */
  
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "netcdf.h"
#include "hydrobase.h"
#include "hydro_cdf.h"

#define PROP_DIM 4        /* # of dimensions describing each property */
#define NDEF_DEPTHS  76  /* dimension of def_depth array */

double def_depth[NDEF_DEPTHS] = 
        { 0,     10,   20,   30,   50,   75,  100,  125,  150, 200,
         250,   300,  350,  400,  450,  500,  550,  600,  650, 700,
	 750,   800,  850,  900,  950, 
	 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900,
	 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900,
	 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900,
	 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5400, 5600, 5800,
	 6000, 6200, 6400, 6600, 6800, 7000, 7500, 8000, 8500, 9000,
          -99
 };

/*********************************************************************/ 

int cdf_init(char *filename)
/* Creates a new netCDF file and returns the netCDF id for that file */
{
   int cdfid, error;

   error = nc_create(filename, NC_CLOBBER, &cdfid);
   if (error != NC_NOERR) {
      fprintf(stderr, "\nUnable to create cdf_file: %s", filename);
      fprintf(stderr, "\n Error: \n%s", nc_strerror(error));
      exit(1);
   }
   return (cdfid);
   
}  /* end cdf_init() */

/***************************************************************************/

int cdf_open(char *dir, char *file_root, char *extent, int print_mess)

   /* opens an existing cdf file for reading and returns the file descriptor 
      to be used in future read operations. If file does not exist, a message
      is printed on stderr and -1 is returned.
      arguments:
      
      dir;           directory of file  OR "" 
      file_root;     root of file name  OR entire file_path 
      extent;        file extent        OR "" 
      print_mess;    1/0 to print/suppress message on stderr 
      
   */
{
   int cdfid;
   int i, error;
   char fname[500];

/* construct filename from arguments */
        strcpy(fname, dir);
        if ((i = strlen(dir)) != 0) {
           if (dir[i-1] != '/')
              strncat(fname,"/",1);
        }

        strcat(fname, file_root);

        if ((strlen(extent) > 0) && (extent[0] != '.'))
           strncat(fname,".",1);
        strcat(fname, extent);

   error = nc_open(fname, NC_NOWRITE, &cdfid); 
   
   if (error == NC_NOERR) { 
      if (print_mess) 
          fprintf(stderr,"\nOpened %s ...", fname);
      return(cdfid);
   }
   else {
      if (print_mess) {
         fprintf(stderr,"\n Unable to open %s.\n", fname);
         fprintf(stderr, "\n Error: \n%s", nc_strerror(error));
      }
      return(-1);
   }


} /* end cdf_open() */

/***************************************************************************/

int cdf_update(char *dir, char *file_root, char *extent, int print_mess)

   /* opens an existing cdf file for reading/writing and returns the 
      file descriptor to be used in future read/write operations. 
      If file does not exist, a message is printed on stderr and -1 
      is returned.
      arguments:
      
      dir;           directory of file  OR "" 
      file_root;     root of file name  OR entire file_path 
      extent;        file extent        OR "" 
      print_mess;    1/0 to print/suppress message on stderr 
      
   */
{
   int cdfid;
   int i, error;
   char fname[500];

/* construct filename from arguments */
        strcpy(fname, dir);
        if ((i = strlen(dir)) != 0) {
           if (dir[i-1] != '/')
              strncat(fname,"/",1);
        }

        strcat(fname, file_root);

        if ((strlen(extent) > 0) && (extent[0] != '.'))
           strncat(fname,".",1);
        strcat(fname, extent);

   if (print_mess) {
      if ((error = nc_open(fname, NC_WRITE, &cdfid)) != NC_NOERR) {
         fprintf(stderr,"\n Unable to open %s.\n", fname);
         fprintf(stderr, "\n Error: \n%s", nc_strerror(error));
     }
      else {
          fprintf(stderr,"\nOpened %s ...", fname);
      }
   }

   return (cdfid);

} /* end cdf_update() */

/***************************************************************************/

void cdf_close(int cdfid)
{
   nc_close(cdfid);
   fprintf(stderr,"\nNetcdf file closed.");
   return;

} /* end cdf_close() */

/*********************************************************************/ 

int cdf_define(int cdfid, struct CDF_HDR *hptr, int prefill, int add_count_vars)

/* Defines all variables and attributes for the CDF file based on the info
   passed in struct CDF_HDR. If the prefill flag is set, each variable will
   be prefilled with the fill-value assigned in the header info.  By not
   prefilling, the write operations will be more efficent, however each point
   must be filled by the program calling the writing routines.  If 
   add_count_vars is set, variables will be defined to store the number of
   observations which determined the value of each property at each point.
   By convention, these will be named by appending "_cnt" to the property
   mnemonic.  
  
   After writing this info to the cdf_file, it takes the open file out of
   define mode.  Returns -1 if an error occurs.
   
   arguments:
   int cdfid:               id of already open cdf file 
   struct CDF_HDR *hptr:    header defined in hydro_cdf.h 
   int prefill:             0 = NOT prefill; or 1 = prefill 
   int add_count_vars:      1/0 = include/not_include variables counting # obs for each prop
*/ 
{
   int lat_dim, lon_dim, depth_dim, time_dim;
   int depth_id, bottom_id;
   int minyear_id, maxyear_id;
   int prop_id[MAXPROP * 2];             /* accommodate all props + counts */
   int dims[PROP_DIM], nchars, error;
   int oldfillmode;
   char *descrip;
   char varname[10];
   int i, nvars;
   short e_code;
   short short_val = 0;                   /* fill value for count variables */

/* define dimensions */

   error = nc_def_dim(cdfid, "lat", (size_t) hptr->ny, &lat_dim);
   if (error != NC_NOERR) {
      fprintf(stderr, "\n Error defining lat_dim: \n%s", nc_strerror(error));
      return (-1);
   }
   error = nc_def_dim(cdfid, "lon", (size_t) hptr->nx, &lon_dim);
   error = nc_def_dim(cdfid, "de", (size_t) hptr->nz, &depth_dim);
   error = nc_def_dim(cdfid, "time", (size_t) hptr->nt, &time_dim);

/* define variables */

   dims[0] = depth_dim;
   error = nc_def_var(cdfid, "de", NC_FLOAT, 1, dims, &depth_id);
   if (error != NC_NOERR) {
      fprintf(stderr, "\nCan't define variable:\n%s", nc_strerror(error));
      return (-1);
   }
   error = nc_put_att_text(cdfid, depth_id, "units", 7, "meters");
   error = nc_put_att_text(cdfid, depth_id, "long_name", 11, "DEPTH  (M)");
   error = nc_put_att_text(cdfid, depth_id, "generic_name", 6, "depth");
   e_code = 3;
   error = nc_put_att_short(cdfid, depth_id, "epic_code", NC_SHORT, 1, &e_code);

   dims[0] = time_dim;
   error = nc_def_var(cdfid, "minyear", NC_INT, 1, dims, &minyear_id);
   error = nc_put_att_text(cdfid, minyear_id, "units", 6, "years");

   error = nc_def_var(cdfid, "maxyear", NC_INT, 1, dims, &maxyear_id);
   error = nc_put_att_text(cdfid, maxyear_id, "units", 6, "years");

   dims[0] = time_dim;
   dims[1] = lat_dim;
   dims[2] = lon_dim;

   e_code = 0;
   error = nc_def_var(cdfid,"bottom", NC_FLOAT, 3, dims, &bottom_id);
   error = nc_put_att_text(cdfid, bottom_id, "units", 7, "meters");
   error = nc_put_att_text(cdfid, bottom_id, "long_name", 13, "Bottom Depth");
   error = nc_put_att_short(cdfid, bottom_id, "epic_code", NC_SHORT, 1, &e_code);

   dims[0] = time_dim;
   dims[1] = lat_dim;
   dims[2] = lon_dim;
   dims[3] = depth_dim;
   for (i = 0; i < hptr->nprops; ++i) {
      error = nc_def_var(cdfid, hptr->prop_id[i], NC_FLOAT, 4, dims, &prop_id[i]);
      nchars = strlen(hptr->prop_units[i]) + 1;
      error = nc_put_att_text(cdfid, prop_id[i], "units", nchars, hptr->prop_units[i]);
      error = nc_put_att_float(cdfid, prop_id[i], "_FillValue", NC_FLOAT, 1, &hptr->fill_value);
      error = nc_put_att_float(cdfid, prop_id[i], "_MaskValue", NC_FLOAT, 1, &hptr->mask_value);
      descrip = get_prop_descrip(get_prop_indx(hptr->prop_id[i]));
      nchars = strlen(descrip) + 1;
      error = nc_put_att_text(cdfid, prop_id[i], "long_name", nchars, descrip);
   }
   nvars = i;
   if (add_count_vars) {   /* add a count var for each prop + one for depth */
      for (i= 0; i < hptr->nprops; ++i) {
         strncpy(varname, hptr->prop_id[i], 3);
         strcat (varname, COUNT_VAR_SUFFIX);
         error = nc_def_var(cdfid, varname, NC_SHORT, 4, dims, &prop_id[nvars]);
         error = nc_put_att_short(cdfid, prop_id[nvars], "_FillValue", NC_SHORT, 1, &short_val);
         ++nvars;
      }
         error = nc_def_var(cdfid, "de_cnt", NC_SHORT, 4, dims, &prop_id[nvars]);
         error = nc_put_att_short(cdfid, prop_id[nvars], "_FillValue", NC_SHORT, 1, &short_val);
   }

/* assign global attributes */ 

   error = nc_put_att_float(cdfid, NC_GLOBAL, "latmin", NC_FLOAT, 1, &hptr->ymin);
   error = nc_put_att_float(cdfid, NC_GLOBAL, "latmax", NC_FLOAT, 1, &hptr->ymax);
   error = nc_put_att_float(cdfid, NC_GLOBAL, "latincr", NC_FLOAT, 1, &hptr->yincr);
   error = nc_put_att_float(cdfid, NC_GLOBAL, "lonmin", NC_FLOAT, 1, &hptr->xmin);
   error = nc_put_att_float(cdfid, NC_GLOBAL, "lonmax", NC_FLOAT, 1, &hptr->xmax);
   error = nc_put_att_float(cdfid, NC_GLOBAL, "lonincr", NC_FLOAT, 1, &hptr->xincr);
   error = nc_put_att_int(cdfid, NC_GLOBAL, "node_offset", NC_INT, 1, &hptr->node_offset);
   error = nc_put_att_int(cdfid, NC_GLOBAL, "nprops", NC_INT, 1, &hptr->nprops);
   error = nc_put_att_int(cdfid, NC_GLOBAL, "counts_included", NC_INT, 1, &add_count_vars);
   nchars = strlen(hptr->title) + 1;
   error = nc_put_att_text(cdfid, NC_GLOBAL, "title", nchars, hptr->title);
   nchars = strlen(hptr->command) + 1;
   error = nc_put_att_text(cdfid, NC_GLOBAL, "command", nchars, hptr->command);
   error = nc_put_att_text(cdfid, NC_GLOBAL, "compliance",32,"HydroBase PMEL/EPIC WHOI/OARS");

   if (!prefill) {
      error = nc_set_fill(cdfid, NC_NOFILL, &oldfillmode);
   }
   else {
       error = nc_set_fill(cdfid, NC_FILL, &oldfillmode);
   }
   error = nc_enddef(cdfid);       /* leave define mode */

   return 0;

}  /* end cdf_define() */

/*********************************************************************/ 


int write_prop_cdf(int cdfid, float *data_ptr, char *prop_mne, int row, int col, int tbin, int zlev, int nrows, int ncols, int ntbins, int nlevs)

/* Writes the data representing a single property at a single gridnode for 
   n time bins and all depths to the already opened netcdf file.  Returns 0  
   for a successful write or prints an error message and exits in case of an 
   error.
   
   arguments:
   int cdfid;             cdf file identifier 
   float *data_ptr;       start of data array
   char *prop_mne;        2-char property mnemonic
   int row, col, tbin, zlev;    grid position to start writing data 
   int nrows, ncols;      # of contiguous rows and cols stored in data array 
   int ntbins, nlevs;     # of contiguous time bins and zlevs to output


*/
{
   size_t start[PROP_DIM], count[PROP_DIM];
   int error, varid;
   
   start[0] = tbin;      
   start[1] = row;
   start[2] = col;
   start[3] = zlev;

   count[0] = ntbins;
   count[1] = nrows;
   count[2] = ncols;
   count[3] = nlevs;

   error = nc_inq_varid(cdfid, prop_mne, &varid);
   error = nc_put_vara_float(cdfid, varid, start, count, data_ptr);
   if (error != NC_NOERR) {
      fprintf(stderr,"\nError writing property %.2s to cdf_file.\n", prop_mne);
      fprintf(stderr, "\n %s", nc_strerror(error));
      fprintf(stderr,"\n\n       write_prop_cdf() exiting.\n\n");
   }
   return 0;

}  /* end write_prop_cdf() */

/*********************************************************************/ 


int write_prop_count_cdf(int cdfid, short *cnt_ptr, char *prop_mne, int row, int col,int tbin, int zlev, int nrows, int ncols, int ntbins, int nlevs)

/* Writes the data representing the number of observations at each depth for a  
   property at a lat/lon position to an already opened netcdf file.
   Returns 0 for a successful write or prints an error message
   and exits in case of an error.

   arguments:
   cdfid;           cdf file identifier  
   cnt_ptr;         start of count array  
   prop_mne;        2-char property mnemonic  
   row, col, zlev;  grid position to start writing data  
   tbin;             "                              "    
   nrows, ncols;    # of contiguous rows and cols stored in count array  
   ntbins;          # of contiguous time bins in count array  
   nlevs;           # of zlevs to output  


*/
{
   size_t start[PROP_DIM], count[PROP_DIM];
   char varname[20];
   int error, varid;
   
   start[0] = tbin;      
   start[1] = row;
   start[2] = col;
   start[3] = zlev;

   count[0] = ntbins;
   count[1] = nrows;
   count[2] = ncols;
   count[3] = nlevs;

   strncpy(varname, prop_mne, 3);
   strcat(varname, COUNT_VAR_SUFFIX);
   error = nc_inq_varid(cdfid, varname, &varid);
   error = nc_put_vara_short(cdfid, varid, start, count, cnt_ptr);
   if (error != NC_NOERR) {
      fprintf(stderr,"\nError writing property %.2s to cdf_file.\n", prop_mne);
      fprintf(stderr, "\n %s", nc_strerror(error));
      fprintf(stderr,"\n\n       write_prop_count_cdf() exiting.\n\n");
   }
   return 0;

}  /* end write_prop_count_cdf() */

/***************************************************************************/
int write_std_depths_cdf(int cdfid, struct CDF_HDR *hptr)

/* Write the standard depth values contained in global variable std_depth[]
    to an already opened cdf file. */
{
   int error, i;
   size_t start[1], count[1];
   int depth_id;
   float *depth_ptr;


   start[0] = 0;
   count[0] = hptr->nz;

   if (!std_depth_initialized) {
       fprintf(stderr, "\n standard depths not initialized! : write_std_depths_cdf()\n");
       exit(1);
   }

 /*  They must be of type float; so first create an array of appropriate type*/

   depth_ptr = (float *) malloc (hptr->nz * sizeof(float));
   for (i = 0; i < hptr->nz; ++i)
      depth_ptr[i] = (float) std_depth[i];

   error = nc_inq_varid(cdfid, "de", &depth_id);
   if (error != NC_NOERR) {
      fprintf(stderr,"\nDepth undefined in the cdf file");
      fprintf(stderr, "\n %s", nc_strerror(error));
       exit(1);
   }

   error = nc_put_vara_float(cdfid, depth_id, start, count, depth_ptr);
   free((void *)depth_ptr);
   if (error != NC_NOERR) {
      fprintf(stderr, "\nError writing standard depths to cdf file.\n");
      fprintf(stderr, "\n %s", nc_strerror(error));
      exit(1);
   }
   return(0);
} /* end write_std_depths_cdf() */
/***************************************************************************/
int write_bottom_depth_cdf(int cdfid, int row, int col, int tbin, int nrows, int ncols, int ntbins, float *data_ptr)

   /* arguments:
       cdfid;                   cdf file identifier 
       row, col, tbin;          grid position to start writing data 
       nrows, ncols, ntbins;   # of datapts in each dimension to output 
       data_ptr;              start of data array 
   */
{
   int  bottom_id, error;
   size_t start[3], count[3];

   start[0] = tbin;
   start[1] = row;
   start[2] = col;

   count[0] = ntbins;
   count[1] = nrows;
   count[2] = ncols;

   error = nc_inq_varid(cdfid, "bottom", &bottom_id);
   error = nc_put_vara_float(cdfid, bottom_id, start, count, data_ptr);
   if (error != NC_NOERR) {
      fprintf(stderr, "\nError writing bottom depths to cdf file.\n");
      fprintf(stderr, "\n %s", nc_strerror(error));
      exit(1);
   }
   return (error);

} /* end write_bottom_depth_cdf() */
/***************************************************************************/
int write_time_bins_cdf(int cdfid, struct CDF_HDR *hptr)
/* Write the min/max values defining timebins to an already opened cdf file. */
{
   int error;
   size_t start[1], count[1];
   int minyear_id, maxyear_id;


   start[0] = 0;
   count[0] = hptr->nt;


   error = nc_inq_varid(cdfid, "minyear", &minyear_id);
   error = nc_inq_varid(cdfid, "maxyear", &maxyear_id);

   error = nc_put_vara_int(cdfid, minyear_id, start, count, hptr->tmin);
   error = nc_put_vara_int(cdfid, maxyear_id, start, count, hptr->tmax);
   if (error != NC_NOERR) {
       fprintf(stderr, "\nError writing time bins to cdf file.\n");
      fprintf(stderr, "\n %s", nc_strerror(error));
      exit(1);
   }
   return(0);
} /* end write_std_depths_cdf() */
/***************************************************************************/


int read_cdf_hdr(int cdfid, struct CDF_HDR *haddr)

   /* Reads an already opened cdf file written with the HydroBase utilities.
      Returns 0 for a successful read.  In case of an error, returns a 
      number > 0.   
               error codes:   1 :  file does not exist;
                              2 :  not a HydroBase cdf file;
   */
{
   int ndims, nvars, ngatts, recdim, natts;
   int lat_dim, lon_dim, depth_dim, time_dim;
   int  dim[MAX_VAR_DIMS];
   char varname[MAX_NC_NAME];
   int error, i, varid, np, n, found;
   char *mne;
   float fval;
   int lval;
   size_t start[1], count[1], len;


   /* get # of properties stored in file */ 

   error = nc_get_att_int(cdfid, NC_GLOBAL, "nprops", &lval);
   np = lval;

   if ((error = nc_inq(cdfid, &ndims, &nvars, &ngatts, &recdim)) != NC_NOERR){
      fprintf(stderr,"\nError reading cdf file.\n");
      return (1);
   }
       

   if ((error = nc_inq_dimid(cdfid, "lat", &lat_dim)) != NC_NOERR ) {
      fprintf(stderr,"\nNot a HydroBase cdf file.\n");
       return (2);
   }
   error = nc_inq_dimlen(cdfid, lat_dim, &len );
   haddr->ny = (int) len;

   if (( error = nc_inq_dimid(cdfid, "lon", &lon_dim)) != NC_NOERR ) {
      fprintf(stderr,"\nNot a HydroBase cdf file.\n");
       return (2);
   }
   error = nc_inq_dimlen(cdfid, lon_dim, &len );
   haddr->nx = (int) len;
   
   if ((error = nc_inq_dimid(cdfid, "de", &depth_dim)) != NC_NOERR ) {
      fprintf(stderr,"\nNot a HydroBase cdf file.\n");
       return (2);
   }
   error = nc_inq_dimlen(cdfid, depth_dim, &len );
   haddr->nz = (int) len;

   if ((error = nc_inq_dimid(cdfid, "time", &time_dim)) != NC_NOERR ) {
      fprintf(stderr,"\nNot a HydroBase cdf file.\n");
       return (2);
   }
   error = nc_inq_dimlen(cdfid, time_dim, &len);
   haddr->nt = (int) len;

   /* Check each cdf variable .
      Get mnemonic and units for each property, and ranges for time bins.  */
   
   haddr->prop_id = (char **) malloc(np * sizeof(char *));
   haddr->prop_units = (char **) malloc(np * sizeof(char *));
   haddr->tmin = (int *) malloc(haddr->nt * sizeof(int));
   haddr->tmax = (int *) malloc(haddr->nt * sizeof(int));

   n = 0;
   for (varid = 0; varid < nvars; ++varid) {

      error = nc_inq_varname(cdfid, varid, varname);

      found = 0;
      i = 0;
      while ((!found) && (i < MAXPROP)) {   /* Is it a property? */
          mne = get_prop_mne(i++);
          found = ! (strncmp(varname, mne, 3));  
      }

      if (found) {           

         if (strncmp(mne, "de", 2) == 0) {       /* for depth, just get units */
            error = nc_get_att_text(cdfid, varid, "units", haddr->z_units);
         }
         else {
            haddr->prop_id[n] = (char *) malloc(3);
            strncpy(haddr->prop_id[n], mne, 3);
            error = nc_inq_attlen(cdfid, varid, "units", &len);
            haddr->prop_units[n] = (char *) malloc(len);
            error = nc_get_att_text(cdfid, varid, "units", haddr->prop_units[n]);
            ++n;
         }
      }
      else {
        start[0] = 0;
        count[0] = haddr->nt;
        if (strncmp(varname, "minyear", 7) == 0) {  /* get the time bin mins */
           error = nc_get_vara_int(cdfid, varid, start, count, haddr->tmin);
        }
        if (strncmp(varname, "maxyear", 7) == 0) {  /* get the time bin max */
           error = nc_get_vara_int(cdfid, varid, start, count, haddr->tmax);
        }
      }

   }

   if (n != np) {
     fprintf(stderr,"\nWARNING!");
     fprintf(stderr,"\n # of properties found in file [%d] does not match", n);
     fprintf(stderr,"\n  the # of props in file attribute information [%d]\n", np);
   }
   haddr->nprops = n; 



   strncpy(haddr->x_units, "degrees", 8);    /* longitude */
   strncpy(haddr->y_units, "degrees", 8);    /* latitude */
   strncpy(haddr->t_units, "years", 6);      /* time */


   /* get the fill value by checking the first variable's attribute */
       
   error = nc_inq_varid(cdfid, haddr->prop_id[0], &varid);
   error = nc_get_att_float(cdfid, varid, "_MissingValue", &fval );
   
   if (error != NC_NOERR) {  /* for backward compatibility: it used to be called _FillValue */
       error = nc_get_att_float(cdfid, varid, "_FillValue", &fval );
   }
   haddr->fill_value = fval;
   
   /* get the mask value by checking the first variable's attribute */
       
   error = nc_inq_varid(cdfid, haddr->prop_id[0], &varid);
   error = nc_get_att_float(cdfid, varid, "_MaskValue", &fval );
   
   if (error != NC_NOERR) {  /* for backward compatibility: this attribute wasn't always defined*/
       fval = (float) HBMASK;
   }
   haddr->mask_value = fval;

   /* get global attributes */

   error = nc_get_att_float(cdfid, NC_GLOBAL, "latmin", &fval);
   haddr->ymin = fval;
   error = nc_get_att_float(cdfid, NC_GLOBAL, "latmax", &fval);
   haddr->ymax = fval;
   error = nc_get_att_float(cdfid, NC_GLOBAL, "latincr", &fval);
   haddr->yincr = fval;
   error = nc_get_att_float(cdfid, NC_GLOBAL, "lonmin", &fval);
   haddr->xmin = fval;
   error = nc_get_att_float(cdfid, NC_GLOBAL, "lonmax", &fval);
   haddr->xmax = fval;
   error = nc_get_att_float(cdfid, NC_GLOBAL, "lonincr", &fval);
   haddr->xincr = fval;
   error = nc_get_att_text(cdfid, NC_GLOBAL, "title", haddr->title);
   error = nc_get_att_int(cdfid, NC_GLOBAL, "node_offset", &lval);
   haddr->node_offset = lval;
   error = nc_get_att_int(cdfid, NC_GLOBAL, "counts_included", &lval);
   haddr->counts_included = lval;
   error = nc_get_att_text(cdfid, NC_GLOBAL, "command", &haddr->command[0]);

   return(0);
}  /* end read_cdf_hdr() */

/***************************************************************************/


int read_cdf_depths(int cdfid, float *d_ptr)

   /* Retrieves the standard depth values used as a coordinate variable. 
      Returns the number of depth values or -1 in case of an error.
      
      arguments:
      int cdfid:      id of cdf file already opened for reading
      float *d_ptr:   address to return depth values 
   */
{
   int depth_dim, varid, error;
   size_t  start[1], count[1];

   if ((error = nc_inq_dimid(cdfid, "de", &depth_dim)) != NC_NOERR) {
       fprintf(stderr,"\n Unable to retrieve depth dimension from cdf file\n");
       return (-1);
   }

   if (nc_inq_dimlen(cdfid, depth_dim, &count[0]) != NC_NOERR)
       return (-1);

   if ( nc_inq_varid(cdfid, "de", &varid) != NC_NOERR)
       return (-1);

   start[0] = 0;
   if ((error = nc_get_vara_float(cdfid, varid, start, count, d_ptr)) != NC_NOERR){
       fprintf(stderr,"\n Error retrieving depths from cdf file:");
       fprintf(stderr,"\n%s\n", nc_strerror(error));
       return (-1);
   }
   return ((int) count[0]);

} /* end read_cdf_depths() */
      
/***************************************************************************/
int read_cdf_bottom_depth(int cdfid, float *d_ptr, int row, int col, int tbin)

   /* Retrieves a single bottom depth value at row, col from a cdf file. 
      Returns 1 for a successful read, or -1 in case of an error. 
      
      arguments: 
       cdfid;             id of cdf file already opened for reading 
       d_ptr;             address to return depth value 
       row, col, tbin;    offset into grid 
 */
{
   int varid, error;
   size_t start[3], count[3];

   start[0] = tbin;
   start[1] = row;
   start[2] = col;

   count[0] = 1;
   count[1] = 1;
   count[2] = 1;

   if ((error = nc_inq_varid(cdfid, "bottom", &varid)) != NC_NOERR) {
       fprintf(stderr,"\n Error retrieving bottom depth id from cdf file:");
       fprintf(stderr,"\n%s\n", nc_strerror(error));
       return (-1);
   }

   if ((error = nc_get_vara_float(cdfid, varid, start, count, d_ptr)) != NC_NOERR) {
       fprintf(stderr,"\n Error retrieving bottom depth value from cdf file:");
       fprintf(stderr,"\n%s\n", nc_strerror(error));
       return (-1);
   }
   
   return(1);

}  /* read_cdf_bottom_depth() */
/***************************************************************************/

int read_cdf_prop(int cdfid, char *varname, float *dataptr, int row, int col, int tbin, int zlev, int npts)

   /* Reads an already opened cdf file written with the HydroBase utilities.
      Dataptr is the address of an array with enough space to hold npts
      amount of data. 
      Returns 0 for a successful read.  In case of an error, returns a 
      number > 0.
         
        error codes:   1 :  file not open;
	       
	arguments:
	       
	cdfid:     id of cdf file already opened for reading 
	varname:   property mnemonic 
	dataptr:   address of array into which values will be read 
	row:       row corresponding to lat of gridnode to read 
	col:       column corresponding to lon of gridnode to read 
	tbin:      index to time bin 
	zlev:      depth level at which to begin 
	npts:      # of values to read 
   */
{
   int varid, error;
   size_t start[PROP_DIM], count[PROP_DIM];

   error = nc_inq_varid(cdfid, varname, &varid);
   if (error != NC_NOERR ) {
     fprintf(stderr,"\nError finding property (%s)\n", varname);
     fprintf(stderr,"\n%s\n", nc_strerror(error));
     return(1);
   }

   start[0] = tbin;     
   start[1] = row;
   start[2] = col;
   start[3] = zlev;

   count[0] =  1;
   count[1] =  1;
   count[2] =  1;
   count[3] =  npts;

   error = nc_get_vara_float(cdfid, varid, start, count, dataptr);
   if (error != NC_NOERR) {
     fprintf(stderr,"\nError reading values for property (%s)\n", varname);
     fprintf(stderr,"\n%s\n", nc_strerror(error));
      return (1);
   }

   return 0;

} /* read_cdf_prop() */

/***************************************************************************/


int read_cdf_prop_count(int cdfid, char *prop_mne, short *cnt_ptr, int row, int col, int tbin, int zlev, int npts)

   /* Reads an already opened cdf file written with the HydroBase utilities.
      cnt_ptr is the address of an array with enough space to hold npts
      amount of data. 
      Returns 0 for a successful read.  In case of an error, returns a 
      number > 0.   
               error codes:   1 :  file not open;


	arguments:
	cdfid     id of cdf file already opened for reading 
	prop_mne  2-character property mnemonic 
	cnt_ptr   address of array into which values will be read
	row       lat of gridnode to read
	col       lon of gridnode to read 
	tbin      index of time bin to read 
	zlev      depth level at which to begin 
	npts      # of values to read 

   */
{
   int varid, error;
   char varname[10];
   size_t start[PROP_DIM], count[PROP_DIM];

   strncpy(varname, prop_mne, 3);
   strcat(varname, COUNT_VAR_SUFFIX);
   error = nc_inq_varid(cdfid, varname, &varid);
   if ( error != NC_NOERR) {
     fprintf(stderr,"\nError finding property (%s)\n", varname);
     fprintf(stderr,"\n%s\n", nc_strerror(error));
     return(1);
   }

   start[0] = tbin;    
   start[1] = row;
   start[2] = col;
   start[3] = zlev;

   count[0] =  1;
   count[1] =  1;
   count[2] =  1;
   count[3] =  npts;

   error = nc_get_vara_short(cdfid, varid, start, count, cnt_ptr);
   if (error != NC_NOERR) {
     fprintf(stderr,"\nError reading property values(%s)\n", varname);
     fprintf(stderr,"\n%s\n", nc_strerror(error));
     return (1);
   }

   return 0;

} /* read_cdf_prop_count() */

/***************************************************************************/

int get_indices(struct CDF_HDR *hptr, float lat, float lon, int *row_ptr, int *col_ptr)

   /* Translates a latitude, longitude into a corresponding (row,col)
      from the info stored in the CDF_HDR. The lat/lon does
      not necessarily have to fall exactly on a grid node: i.e. each
      grid node has an associated area = yincr * xincr, whose definition
      depends on whether it is a node grid or a pixel grid. 
      Returns 0 for successful conversion within bounds; or a -1 if the
      row, col does not fall within the grid bounds.
      
   */
{  
   if (hptr->node_offset == 0)  {                         /* for node grids */ 

   /* the extra smidgeon puts stations that fall right on a border into the 
      proper box */


     *row_ptr = NINT((hptr->ymax - (lat + .0001)) / hptr->yincr);
     *col_ptr = NINT((lon + .0001 - hptr->xmin) / hptr->xincr);
   }

   else {                                             /* for pixel grids */

   /* the extra smidgeon puts stations that fall right on a border into the 
      proper box */

     *row_ptr = NINT((hptr->ymax - (lat + .0001)) / hptr->yincr - 0.5);
     *col_ptr = NINT((lon + .0001 - hptr->xmin ) / hptr->xincr - 0.5);
   }

   if ((*row_ptr < 0) || (*row_ptr >= hptr->ny))
      return(-1);
   if ((*col_ptr < 0) || (*col_ptr >= hptr->nx))
      return(-1);
      
   return 0;
 
} /* end get_indices() */

/***************************************************************************/

int get_lat_lon(struct CDF_HDR *hptr, int row, int col, float *lat_ptr, float *lon_ptr)

    /* Translates the row,col into the lat,lon of its associated gridnode.
       Returns 0 for a successful conversion or -1 if the point is out
       of the gridbounds.
    */ 
{
   if ((row >= hptr->ny) || (row < 0) || (col >= hptr->nx) || (col < 0))
      return (-1);


   if (hptr->node_offset == 0) {                          /* for node grids */
      *lat_ptr = hptr->ymax - row * hptr->yincr;
      *lon_ptr = hptr->xmin + col * hptr->xincr;
   }
   else {                                             /* for pixel grids */
      *lat_ptr =  hptr->ymax - (row + .5) * hptr->yincr;
      *lon_ptr =  hptr->xmin + (col + .5) * hptr->xincr;
   }

  return (0);

}  /* end get_lat_lon() */
/**********************************************************************/
int std_depth_init(FILE *depth_file)

   /* Initializes the global variable std_depth and NSTDLEVS using values in 
      depth_file or the default values if depth_file points to NULL. Returns 
      the number of std depths */
{
   int i;

/*  use default standard depths ... */
   if (depth_file == NULL) {
      for (i = 0; i < NDEF_DEPTHS; ++i)
         std_depth[i] = def_depth[i];

      std_depth_initialized = 1;
      NSTDLEVS = NDEF_DEPTHS;

      return(NDEF_DEPTHS);
   }

/*  get standard depths from file ... */
   i = 0;
   while (fscanf(depth_file,"%lf", &std_depth[i] ) == 1) {
     if (++i >= MAXSTDLEVS) {
        fprintf(stderr,"\nOnly able to use %d standard depths from file.\n", i);
        break;
     }
  
   }
   fclose (depth_file);
   NSTDLEVS = i;
   std_depth_initialized = 1;
   return (i);

}  /* end std_depth_init() */

/**********************************************************************/
int d2stdlev(double depth)

   /* Returns the standard level index corresponding to depth.  Or a
      negative number if some error occurs:
                 -1 :  not a standard depth 
                 -2 :  standard depths have not been defined yet.
                 -99:  depth exceeds greatest standard depth.
      The global variables std_depth[] and NSTDLEVS are necessary and 
      should have been initialized prior to calling this function.
   */
{

   int i, d, sd;
   
/*  make sure standard depths have been defined ... */

   if (!std_depth_initialized)
       return (-2);

/*  check each standard depth and compare integer values to avoid problems 
    of floating point storage ... */

   d = (int) (depth +.00001);
   i = -1;
   while (++i < NSTDLEVS) {

      if (d < (sd = (int) (std_depth[i] +.00001)))
          return (-1);

      if (d == sd)
          return(i);
   }

/* if this statement is reached, the depth exceeds greatest standard depth 
     defined */
   return (-99);

}
/**********************************************************************/
double stdlev2depth(int index)

   /* Returns the depth corresponding to stddepth index .  Or a
      -99.0 if the index is not a stddepth index.
   */
{
   if ((index < 0) || (index >= NSTDLEVS))
       return (-99.0);

   return (std_depth[index]);

} /* end of stdlev2depth */
/***************************************************/

int is_flagged(float x, float flag)

  /* Determines whether x is a missing value based on
  the value of flag.  If flag is negative, values less
  than or equal to flag are considered missing.  For a 
  positive flag, values greater than or equal to flag 
  are interpreted as missing. 
  Returns 1 if missing, 0 if not.  */
{
  if (flag < 0) {
     if (x <= (flag + 1.0))
        return(1);
	
     return(0);
  }
  
  if (x >= (flag - 1.0))
    return (1);
    
  return (0);
}  /* end is_flagged() */

