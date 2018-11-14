/*  hydro_utils.c
................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
                             1993
                             updated to ANSI C Nov 1999
................................................................................
.
.  A set of routines to facilitate reading, writing files in HydroBase
.  format.
................................................................................
................................................................................
 int  open_hydro_file(char *dir, char *root, char *extent, int print_msg):
       Opens an existing hydro file for reading.
 int  get_station(int file, struct HYDRO_HDR *h_addr, struct HYDRO_DATA *d_ptr):
       Reads a complete station.
 int  read_hydro_hdr(int file, struct HYDRO_HDR *haddr):    
       Gets station header info.
 void report_status(int status, file *file)
       Reports any errors reported by read_hydro_hdr().
 int  get_data_scan(int file, double *scan, int n, int *prop_order) :  
       Reads the next scan in an open hydro file
 int  get_separator(int file) :  
       Advances file beyond the station separator.
 int  create_hydro_file(char *fname, int mode):  
       Opens a HydroBase file for writing with specified mode.
 int  write_hydro_station(int file, struct HYDRO_HDR *hptr, struct HYDRO_DATA *dptr):
       Writes an entire station to an open hydro file .  
 int  write_hydro_hdr(int file, struct HYDRO_HDR *hptr) : 
       Writes a header to an open hydro file 
 int  write_hydro_scan(int file, double *scan, int n, int *prop_order) : 
       Writes an individual data scan 
 int  write_separator(int file) : 
       Writes a separator line.
 int  ms10(float lat, float lon, int *ms1_ptr) : 
       Computes the 10-degree & 1-degree Marsden Square # for a lat,lon pair.
 int  available(enum property p, struct HYDRO_HDR *hptr) :
       Returns a 1 if the specified property is available for a station. 
 int  is_in_range(float lon, float xmin, float xmax, int *merid, int lon0to360)
       Returns a 1 if specified lon is in the range xmin -> xmax 
 void free_and_alloc(double **dptr, int n)
       Frees up space pointed to by dptr and allocates space for n elements.
 void list_origin(FILE *fptr);
       Lists origination code table
 void list_instrument(FILE *fptr);
       Lists instrument code table
................................................................................
................................................................................
*/

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include "hydrobase.h"
/*sfg added include*/
#include <string.h>

#define NORIG 11
char orig_code[NORIG] = {   '0',  /* miscellaneous */
                            '1',  /* NODC / WOD01 */
                            '2',  /* ICES */
                            '3',  /* NSIDC */
                            '4',  /* WHPO */
                            '5',  /* BBSR */
                            '6',  /* WHOI */
                            '7',  /* SIO */
                            '8',  /* WGHC */
                            '9',   /* individual PI */
                            ' '   /* unknown */
 };
                         
char *orig_name[NORIG] = {  "miscellaneous",
                            "NODC/WOD2001" ,
                            "ICES",
                            "NSIDC",
                            "WHPO",
                            "BBSR",
                            "WHOI",
                            "SIO",
                            "WOCE Global Hydrographic Climatology (Gouretski & Koltermann, 2004)",
                            "individual PI",
                            "unknown"
}; 

#define NINSTR 7
char instr_code[NINSTR] = {'b',  /* bottle */
                           'c',  /* ctd */
                           'f',  /* profiling float */
                           's',  /* seasoar */
                           'm',  /* moored profiler */
                           'u',  /* unknown */
                           ' '   /* unspecified */
};
                   
char *instr_name[NINSTR] = {"bottle",
                            "ctd",
                            "float",
                            "seasoar",
                            "moored_profiler",
                            "unknown",
                            "unspecified"
};
      
/**********************************************************************/

int open_hydro_file(char *dir, char *root, char *extent, int print_mess)

  /* opens an existing hydro file for READING only. Returns a file
   descriptor if successful or -1 if an error occurs. 
    
     char *dir       name of directory   or "" 
     char *root      root of filename    or full path to file 
     char *extent    extent for file     or "" 
     int  print_mess 0 to suppress printing messages on stderr device 
  */
{ 

   char fname[800];
   int  file, i;

    strcpy(fname, dir);
     if ((i = strlen(dir)) != 0) {
        if (dir[i-1] != '/')
           strncat(fname,"/",1);
     }
     strcat(fname, root);
     strcat(fname, extent);

     file = open(fname, O_RDONLY);

     if (print_mess) {
        if (file < 0) {
          fprintf (stderr,"\nunable to open %s for input\n\n", fname);
        }
        else {
          fprintf(stderr,"\nOpened %s ...\n", fname);
        }
     }

     return (file);
   
}  /* end open_hydro_file() */

/**********************************************************************/

int get_station(int file, struct HYDRO_HDR *h_addr, struct HYDRO_DATA *d_ptr)

  /* Retrieves an entire station and places the data in struct HYDRO_DATA. 
     It advances beyond the station separator at the end of the current station. 
     Returns 0 for a successful read, -1 for end_of_file,  
     and an error code if an error occurs.

      error code =  1 :  error parsing header.
                    2 :  error parsing datascan.
                    3 :  error reading station separator.
                    4 :  unexpected end of file occurred.
      arguments:              
      file:    already opened file handle 
      h_addr:  station header info (already allocated) 
      d_ptr:   pointer to station data (already allocated)
  */
{
   int i, j, iprop;
   char line[NBSEP];
   double *scan;
 
   if ((i = read_hydro_hdr(file, h_addr)) < 0) 
      return (-1);
   else if (i > 0)
      return (1);
   else
      ;

/* assign values at d_ptr and allocate memory for data scans */

   d_ptr->nobs = h_addr->nobs;
   d_ptr->nprops = h_addr->nprops;
   for (i = 0; i < d_ptr->nprops; ++i) {
      iprop = h_addr->prop_id[i]; 
      if (d_ptr->observ[iprop] != NULL) {
           free((void *)d_ptr->observ[iprop]);
           d_ptr->observ[iprop] = NULL;
      }
      d_ptr->observ[iprop] = (double *) calloc((size_t)h_addr->nobs, sizeof(double));
      if (d_ptr->observ[iprop] == NULL) {
         fprintf(stderr,"\nUnable to allocate memory in get_station()\n");
         exit(1);
      }
   }

   scan = (double *) calloc((size_t)h_addr->nprops, sizeof(double));

/* read each scan from file */
   
   for (i = 0; i < h_addr->nobs; ++i) {
      if (get_data_scan(file, scan, h_addr->nprops, h_addr->prop_id) > 0) {
           free((void *)scan);
           return (2);
      }
      for (j = 0; j < d_ptr->nprops; ++j) {
         iprop = h_addr->prop_id[j]; 
         d_ptr->observ[iprop][i] = scan[j];
      }
   }  /* end for */

   free((void *)scan);

   if (read(file, line, NBSEP) != NBSEP) {   /* move past station separator */
       return (3);
   }

   return (0);
}  /* end get_station() */

/**********************************************************************/

int read_hydro_hdr(int file, struct HYDRO_HDR *haddr)

  /* reads a header from an already opened hydro file, parses the
    information and stores it at haddr.  Returns 0 for a successful
    read, -1 for end-of-file, or a code > 0 to indicate an error.

          error codes returned :
                1 : error parsing header.
                2 : error parsing property codes
  */
{
    char *s, line[NBHEAD];
    int i, nbytes;

    if ( read(file, line, NBHEAD) != NBHEAD)
         return (-1);
    if (sscanf(&line[0],"%s", &haddr->country[0]) != 1){
         fprintf(stderr,"%s",line);
         return (1);
     }
    if (sscanf(&line[3],"%s", &haddr->ship[0]) != 1){
         fprintf(stderr,"%s",line);
          return (1);
     }
    if (sscanf(&line[6],"%d", &haddr->cruise) != 1){
         fprintf(stderr,"%s",line);
         return (1);
     }
    if (sscanf(&line[12],"%d", &haddr->station) != 1){
         fprintf(stderr,"%s",line);
         return (1);
     }
    if (sscanf(&line[17],"%d %d %d",&haddr->year,&haddr->month,&haddr->day) != 3){
         fprintf(stderr,"%s",line);
         return (1);
     }
    if (haddr->year < 1000)
       haddr->year += 1900;
       
    sscanf(&line[68],"%d", &haddr->ms10);  /* offset of ms10 is diagnostic of hdr */
    
    if (haddr->ms10 >= 1000) {               /*  hydrobase2 hdr */
    
      if (sscanf(&line[28],"%c%c", &haddr->origin, &haddr->instrument) != 2)
         return (1);
      if (sscanf(&line[31],"%f %f %d", &haddr->lat, &haddr->lon, &haddr->pdr) != 3)
         return (1);
      if (sscanf(&line[59],"%d %d %d %d", &haddr->nobs, &haddr->nprops, &haddr->ms10, &haddr->ms1) != 4)
         return (1);
         
      for (i = 0; i < 4; ++i)   {    /* load quality codes this way in case of blanks */
         haddr->qual[i] = line[54+i];  
         if (haddr->qual[i] == ' ')
           haddr->qual[i] = '0';
      } 
     
    }
    else {                                  /* old hydrobase hdr */
    
       haddr->origin = '0';
       haddr->instrument = 'u';
       if (sscanf(&line[30],"%f %f %d %d %d %d %d", &haddr->lat, &haddr->lon, &haddr->pdr, &haddr->nobs, &haddr->nprops, &haddr->ms10, &haddr->ms1) != 7)
         return (1);
     }     
   
    
    nbytes = 3 * haddr->nprops;  /* 3 bytes per prop incl trailing whitespace */
    if (read(file, line, nbytes) != nbytes)
        return (2);

    if ( line[nbytes-1] != LF )  /* check that last char read was a LF */
        return (2);

/* allocate space to cross-reference properties for this station */

    if (haddr->prop_id != NULL)
         free((void *)haddr->prop_id);

    haddr->prop_id = (int *) calloc((size_t)haddr->nprops, sizeof(int));

/* now get the properties */
    s = line;
    for (i = 0; i < haddr->nprops; ++i) {
        if ((haddr->prop_id[i] = get_prop_indx(s)) < 0) {
           fprintf(stderr,"\n Unknown property listed in file: %.2s\n", s);
           return (2);
        }
        s += 3;   /* advance to next property */
    }
    return(0);

} /* end read_hydro_hdr() */

/**********************************************************************/

void report_status(int status, FILE *file)
  /* Translates and reports errors during read operations in get_station() 
  */
{
  switch (status) {
       case -1:
            break;
       case 0:
            break;
       case 1:
            fprintf(file,"\nError parsing header in get_station(). \n");
            break;
       case 2:
            fprintf(file,"\nError parsing datascan in get_station(). \n");
            break;
       case 3:
            fprintf(file,"\nError reading station separator in get_station(). \n");
            break;
       case 4:
            fprintf(file,"\nUnexpected eof encountered in get_station(). \n");
            break;
       default:
            fprintf(file,"\n Unknown error code returned by get_station() : %d\n", status);

  } /* end switch */
  return;
}  /* end report_status() */

/**********************************************************************/

int get_data_scan(int file, double *scan,int n, int *prop_order)

  /* Reads and parses one data scan.  Returns 0 for a successful
     read  or an error code if an error occurs.

      error code =  2 :  error parsing datascan.
                    4 :  unexpected end of file occurred.
                    
     arguments:
       file:         already opened file 
       scan:         array to store data values 
       n:            number of properties per scan 
       prop_order:   index to enumerated property list 
  */
{
   char *s, *line;
   int nbytes, i;

/* determine number of bytes in a scan allowing for the LF at the end... */

   nbytes = 0;
   for (i = 0; i < n; ++i) {             
     nbytes += get_field_width(prop_order[i]) + 1;  /* add a space to field */
   }
   ++nbytes;                                        /* allow for LF */
   line = (char *) calloc((size_t)nbytes, sizeof(char));
 
   if (read(file, line, nbytes) != nbytes ) {
          free((void *)line);
          return(4);  
   }

   if (line[nbytes-1] != LF) {
     free((void *)line);
     return(2);
   }

/* parse the line into property fields... */
   s = line;
   for (i = 0; i < n; ++i) {             
      if (sscanf(s,"%lf", &scan[i] ) != 1) {
          free((void *)line);
          return(2);
      }
      s += get_field_width(prop_order[i]) + 1;
   }

   free((void *)line);
   return (0);

} /* get_data_scan() */
/**********************************************************************/

int get_separator(int file)

   /* Advances a file beyond the station separator.  Returns 0 for
      a successful read or 1 if an error occurs.  */
{
   char line[NBSEP];

   if (read(file, line, NBSEP) != NBSEP) {   /* move past station separator */
       return (1);
   }

   return (0);

}  /* end get_separator() */

/**********************************************************************/
int create_hydro_file(char *fname, int mode)

  /* opens a file for writing in specified mode. An error will occur
    if noclobber is specified and the file already exists. Returns a file  
    descriptor for a successful open,  or -1 for an error.
    
    arguments:
       fname: string containing name of file to be opened 
       mode:  0=overwrite if file exists, 1=noclobber, 2=append 
   */
{
   int file;

   switch (mode) {
       case APPEND:
               if ((file = open(fname, 1)) >= 0) {
                  lseek(file,0L,2);  /* position pointer at end of file */
                  return(file);
               }
               /* else fall through */
       case OVERWRITE:
          if ((file = creat(fname, 0666)) < 0) { 
              return(-1);
          }
          return(file);

       case NOCLOBBER:  /* fall through */
       default:
          if ((file = open(fname, 1)) >= 0) {  /* check if file exists already*/
              return (-1);                    /* return an error if it does */
          }

          if ((file = creat(fname, 0666)) < 0) { /* create if it doesn't exist*/
              return(-1);
          }
          return(file);

   } /* end switch */
   
}  /* end create_hydro_file() */
/**********************************************************************/
int write_hydro_station(int file, struct HYDRO_HDR *hptr, struct HYDRO_DATA *dptr)

   /* writes a complete station to an open file.  Returns 0 for
      a successful write or a number > 0 if an error occurs.
   */
{
   int i, j;
   double *scan;

   /* Write the station header and check for good write. */
   if (write_hydro_hdr(file, hptr) > 0)
      return (1);

   /* Allocate space for a single scan based on the number
    * of properties in this station. */
   scan = (double *) calloc((size_t) hptr->nprops, sizeof(double));

   /* For each depth level (each observation)... */
   for (i = 0; i < hptr->nobs; ++i) {

     /* For each property, load a value into the scan line. */
     for (j = 0; j < hptr->nprops; ++j) {
       scan[j] = dptr->observ[hptr->prop_id[j]][i];
     }

     /* Write the scan line and check for good write. */
     if (write_hydro_scan(file, scan, hptr->nprops, hptr->prop_id) > 0) {
       free((void *)scan);
       return(2);
     }
   }
   free((void *)scan);

   write_separator(file);

   return(0);

} /* write_hydro_station() */

/**********************************************************************/
int write_hydro_station_no_missing(int file, struct HYDRO_HDR *hptr, struct HYDRO_DATA *dptr)

   /* writes a station to an open file, skipping any scans
    * that have HB_MISSING in them.  Returns 0 for
    * a successful write or a number > 0 if an error occurs.
    * Kept separate from write_hydro_station because the
    * error checking slows things down.
    *
    * CAREFUL: DO NOT USE YET -> WILL NOT CHANGE hdr.nobs
    * so HB files will become corrupt.
   */
{
  int i, j, scan_ok;
   double *scan;

   /* Write the station header and check for good write. */
   if (write_hydro_hdr(file, hptr) > 0)
      return (1);

   /* Allocate space for a single scan based on the number
    * of properties in this station. */
   scan = (double *) calloc((size_t) hptr->nprops, sizeof(double));

   /* For each depth level (each observation)... */
   for (i = 0; i < hptr->nobs; ++i) {

     /* Set scan checker flag */
     scan_ok = 1;

     /* For each property, load a value into the scan line. */
     for (j = 0; j < hptr->nprops; ++j) {
       scan[j] = dptr->observ[hptr->prop_id[j]][i];

       /* Check for HB_MISSING */
       if ( scan[j] == HB_MISSING ) {
	 scan_ok = 0;
       }
     }

     if ( scan_ok ) {
       /* Write the scan line and check for good write. */
       if (write_hydro_scan(file, scan, hptr->nprops, hptr->prop_id) > 0) {
	 free((void *)scan);
	 return(2);
       }
     }
   }
   free((void *)scan);

   write_separator(file);

   return(0);

} /* write_hydro_station_no_missing() */
 
/**********************************************************************/

int write_hydro_hdr(int file, struct HYDRO_HDR *hptr)

   /* writes a station header. Returns 0 for a successful operation or
      number > 0 if an error occurs. 
               error code : 2 :  error writing line to file. */
{

   char line[NBHEAD], *s;
   int nbytes, i;

/* check that year is upgraded to 4 digits, not 2 */

   if (hptr->year < 1000)
     hptr->year += 1900;
           
      
/* make sure quality bytes are not blank */

   for (i = 0; i < NQUAL; ++i) {
      if (hptr->qual[i] < '0' || hptr->qual[i] > '9')
         hptr->qual[i] = '0';
   }
      
/* make sure origin and instrument are not blank */

   if (hptr->origin == '\0') 
      hptr->origin = ' ';
      

   if (hptr->instrument == '\0') 
      hptr->instrument = ' ';


/* compose first header line and write to file ...*/

    sprintf(line,"%.2s %.2s %5d %4d %4d %2d %2d %c%c %7.3f %8.3f %5d %c%c%c%c %4d %3d %4d %2d", hptr->country, hptr->ship, hptr->cruise, hptr->station, hptr->year, hptr->month, hptr->day, hptr->origin, hptr->instrument, hptr->lat, hptr->lon, hptr->pdr, hptr->qual[0], hptr->qual[1], hptr->qual[2], hptr->qual[3], hptr->nobs, hptr->nprops, hptr->ms10, hptr->ms1); 

   line[NBHEAD-1] = LF;
   if (write(file, line, NBHEAD) != NBHEAD)
     return (2);

/*  compose line of property descriptors and write to file ... */

   nbytes = hptr->nprops * 3;     /* # of bytes needed to write prop descriptors */
   s = line;
   for (i = 0; i < hptr->nprops; ++i) {
      sprintf(s,"%.2s ", get_prop_mne(hptr->prop_id[i]));
      s += 3;
   }
   line[nbytes-1] = LF;

   if (write(file, line, nbytes) != nbytes)
      return (2);

   return(0);
} /* end write_hydro_hdr() */

/**********************************************************************/

int write_hydro_scan(int file, double *scan, int n, int *prop_id)

   /* writes a single data scan to an open hydro file.  Returns 0 for
      a successful write or a number > 0 if an error occurs.
           error code =   2 :  error writing to output file.
    arguments:
         file:   file descriptor of file opened for writing 
         scan:   array of property values in appropriate order 
         n:      number of properties 
        prop_id: index to enum property for each value in scan 

   */
{

   char *line, *s, str[10];
   int nbytes, i, jj, kk, k;

/* determine number of bytes needed for scan */

   nbytes = 0;
   for (i = 0; i < n; ++i) {
      nbytes += get_field_width(prop_id[i]) + 1;
   }
   ++nbytes;         /* allow for LF */
   line = (char *) calloc((size_t)nbytes, sizeof(char));
   s = line;

/* now write value in format appropriate to each property */

   for (i = 0; i < n; ++i ) {
      if (prop_id[i] < 0 || prop_id[i] >= MAXPROP) {
          fprintf(stderr,"FATAL ERROR: prop_id index passed to write_hydro_scan() is out of range!!");
          exit(1);
      }
      sprintf(str," %c%d.%dlf", '%', get_field_width(prop_id[i]), 
              get_field_precis(prop_id[i]));
              
      /* some  values are too large to fit in the field width: */
      
      kk = get_field_width(prop_id[i])- get_field_precis(prop_id[i]) - 1;
           
      jj = 1;
      for (k = 2; k < kk; ++k) {
           jj *= 10;
      }
      if (scan[i] > (99.999 * jj) || scan[i] < (-9.999 * jj)) {
        scan[i] = HB_MISSING;
      }
        
      sprintf(s, str, scan[i]);
      s += get_field_width(prop_id[i]) + 1;
   } /* end for */

   line[nbytes-1] = LF;
   if (write(file, line, nbytes) != nbytes) {
        free((void *)line);
        return (2);
   }

   free((void *)line);
   return (0);

}  /* end write_hydro_scan() */

/**********************************************************************/
int write_separator(int file)

   /* writes a separator line to a hydro file.  Returns 0 for a successful
      write or a code > 0 if an error occurs.
   */
{
   char line[NBSEP];

   line[0] = '*';
   line[1] = '*';
   line[2] = LF;
   if (write(file, line, NBSEP) != NBSEP)
       return (2);

   return (0);     

}  /* end write_separator() */

/**********************************************************************/
int ms10(float lat, float lon, int *ms1_ptr)

   /*  Computes the 10-degree and 1-degree Marsden Square designations for a 
       particular lat/lon combination. */
{
   int  quadrant[2][2] = { 7000, 1000,
                           5000, 3000 };
   int  k, kk;
   int  ilat, ilon;


   k  = (lat < 0.) ? 1 : 0;   /* determine earth quadrant */
   kk = (lon < 0.) ? 0 : 1;
   if (lon >= 180.) {
      lon = lon - 360;
      kk = 0;
   }
   
   /* This smidgeon handles borderline
    * cases correctly in all hemispheres.
    * Integer conversion truncates. */
   ilat = (int) (lat + .00001);  
   ilon = (int) (lon + .00001);  

   if (ilat < 0) ilat = -ilat;
   if (ilon < 0) ilon = -ilon;

   if (ilat >= 90) ilat = 89;    /* special case at the poles */

   *ms1_ptr = (ilat % 10) * 10 + (ilon % 10);

   return (quadrant[k][kk] + (ilat / 10) * 100 + (ilon / 10));

}  /* end ms10() */

/**********************************************************************/
int available(enum property prop, struct HYDRO_HDR *hptr)

   /*  Determines if the specified property is available at a particular
       station.  Returns 0 if not, 1 if available.  */
{
   int i, found;

   found = 0;
   i = -1;
   while (++i < hptr->nprops) {
         if (hptr->prop_id[i] == (int) prop)
             return (1);
   }
   return (0);

}  /* end available() */
/****************************************************************************/
int  is_in_range(float lon, float xmin, float xmax, int *merid, int lon0to360)

    /*  Returns 1 if lon falls within the longitide range xmin -> xmax
       or 0 if it does not
      for any combination of negative/positive longitude values.
      lon0to360 must be set upon entry into this function:
             = 1  if longitude range is all positive values
	     = 0   if mixed and crosses Greenwich meridian
	     = -1  if all negative values
       merid is an array used to quickly check whether a lon is in the interval
          by assigning a 0 or 1 values to each of its 360 elements. If merid is NULL
	  upon entry into the function, the array is created and its values are set
	  based upon the values of xmin,xmax and lon0to360.  Subsequent entries use
	  this array
	  
   */	   
{
   int i, index1, index2;
   
   if (merid == NULL) {
     merid = (int *) calloc(360, sizeof(int));
     
     if (lon0to360 > 0) {              /*lon range is all positive */
        index1 = NINT(xmin);
	index2 = NINT(xmax -1);
	for (i = index1; i <= index2; ++i) 
	    merid[i] = 1;
	
     }
     else if (lon0to360 < 0) {  /*lon range is all negative */
        index1 = NINT(xmin + 360.0);
	index2 = NINT(xmax + 360.0 - 1);
	for (i = index1; i <= index2; ++i) 
	    merid[i] = 1;
     
     }
     else {   /*  range crosses Greenwich */
          index1 = NINT(xmin + 360.0);
	  for (i = index1; i < 360; ++i) 
	     merid[i] = 1;
	  index2 = NINT(xmax - 1);
	  for (i = 0; i <= index2; ++i)
	     merid[i] = 1;
     }
   }  /* end if merid == NULL */
   
   i = NINT(lon);
   if (lon < 0)
       i += 360;
       
   if (i >= 360) 
       i = 359; 
         
   if (merid[i] ) 
        return (1);
   return (0);

}  /* end is_in_range() */
/****************************************************************************/

void free_and_alloc(double **dptr, int n)

   /* Frees up the space pointed to by dptr.  It MUST have been allocated
      using malloc()!!  Then mallocs space for n double values and sets
      dptr pointing to it.  If insufficient memory, a message is written to
      stderr and an exit() is made. */
{

   if (*dptr != NULL) {
      free((char *)*dptr);
      *dptr = NULL;
   }

   *dptr = (double *) malloc(n * sizeof(double));

   if (*dptr == NULL) {
      fprintf(stderr, "\nInsufficient memory for call to malloc()\n");
      exit(1);
   }

   return;
}  /* end free_and_alloc() */
/****************************************************************************/
void list_origin(FILE *fptr)
  /* List out the origination code to file */
{
 int i;
 
 for (i = 0; i < NORIG; ++i ){
   fprintf(fptr,"  %c  %s\n", orig_code[i], orig_name[i]);
 }
 return;
}
/****************************************************************************/
void list_instrument(FILE *fptr)
  /* List out the instrument code to file */
{
 int i;
 
 for (i = 0; i < NINSTR; ++i ){
   fprintf(fptr,"  %c  %s\n", instr_code[i], instr_name[i]);
 }
 return;
}
/****************************************************************************/
