/* hydrobase.h 

................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth G Curry
                             Woods Hole Oceanographic Instition
                             1993  
                             updated 1999 to ANSI standards
................................................................................
*/


#ifndef    HYDROBASE2
#define    HYDROBASE2  1

/** Definitions for HydroBase structures and supported properties **/
/** Additional property definitions (such as names, units, and field
 ** writing information are in prop_subs.c **/
#define  MAXPROP  42
enum property { PR, DE, TE, TH, TP, TZ, HC, SA, OX, O2, N2, N3, P4, SI, HT, IH, PE, S0, S1, S2, S3, S4, S_, SD, BF, PV, SV, VA, F1, F2, F3, HE, TU, GN, GE, VN, VE, VS, DR, AL, BE, RR };

#define    NBHEAD    76      /* # of bytes in header including LF */
#define    NBSEP      3      /* # of bytes separating stations incl LF */
#define    NQUAL      4      /* # of bytes for quality control id */

#define   OVERWRITE  0       /* define modes for output files */
#define   NOCLOBBER  1
#define   APPEND     2 

/* No interpolation for vertical data gaps exceeding GAP_SHALLOW
 * and GAP_DEEP, with the boundary at GAP_CHANGE_DEPTH.
 * Added . after each number to make them floating point since
 * compared to pressure. */
#define   GAP_SHALLOW  210.     /* datagap in thermocline > 210 m , actual value was 260 for 250m, but will reset to 210, sfg*/
#define   GAP_DEEP   610.       /* datagap below thermocline > 610 m */
#define   GAP_CHANGE_DEPTH 1001. /* transition depth for shallow and deep gap*/

#define   HB_MISSING  -9.0     /* denotes missing observation */

struct HYDRO_HDR {
       char country[3], ship[3];      /* 2 char code + end-of-string char */
       char origin, instrument;            /* 1 char codes */
       int  cruise, station, year, month, day;
       float lat, lon;
       int  pdr, nobs, nprops, ms10, ms1;
       int *prop_id;               /* index to observ_prop identifiers */
       char qual[NQUAL];           /* quality control bytes */ 
};

struct HYDRO_DATA {
       int nobs, nprops;
       double *observ[MAXPROP];
};

/* frequently used macros */
/* Makes use of the C ternary operator. First one is the
 * same as if { x < 0 } then {(int) x-0.5 } else {(int) x+0.5 } */
#define   NINT(x)       (((x) < 0) ? ((int)((x) - 0.5)) : ((int)((x) + 0.5)))
#define    ABS(x)       (((x) < 0) ? -(x) : (x))


/* frequently used definitions */

#define   STDIN    0
#define   STDOUT   1
#define   STDERR   2
#define   LF  0x0a       /* ASCII code for linefeed */

#ifndef TRUE
#define TRUE  1
#endif

#ifndef  FALSE
#define  FALSE 0
#endif

#endif /*ifndef HYDROBASE2 */

/* prototypes for functions in hydro_utils.c  */

#ifndef HYDRO_UTILS
#define HYDRO_UTILS 1
extern  int  open_hydro_file(char *, char *, char *, int );
extern  int  get_station(int, struct HYDRO_HDR *, struct HYDRO_DATA *);
extern  int  read_hydro_hdr(int, struct HYDRO_HDR *);    
extern  void report_status(int, FILE *);
extern  int  get_data_scan(int, double *, int, int *); 
extern  int  get_separator(int );  
extern  int  create_hydro_file(char *, int );  
extern  int  write_hydro_station(int, struct HYDRO_HDR *, struct HYDRO_DATA *);
extern  int  write_hydro_hdr(int, struct HYDRO_HDR *); 
extern  int  write_hydro_scan(int, double *, int, int *); 
extern  int  write_separator(int); 
extern  int  ms10(float, float, int *); 
extern  int  available(enum property, struct HYDRO_HDR *);
extern  int  is_in_range(float, float, float, int *, int );
extern  void free_and_alloc(double **, int);
extern  void list_origin(FILE *);
extern  void list_instrument(FILE *);
extern  void *get_memory(void *, size_t, size_t);
#endif  /*ifndef HYDRO_UTILS*/

/* prototypes for functions in phyprops.c  */

#ifndef PHYPROPS
#define PHYPROPS 1
  extern double hb_sal78(double,double,double,int); 
  extern double hb_svan(double, double, double, double *); 
  extern double hb_depth(double, double); 
  extern double hb_tf(double, double);
  extern double hb_cpsw(double, double, double);
  extern double hb_atg(double, double, double);
  extern double hb_theta(double, double, double, double);
  extern double hb_svel(double, double, double);
  extern double hb_bvfrq(double *, double *, double *, int, double *, double *);
  extern double hb_bvgrdts(double *, double *, double *, double *,int, double *);
  extern double hb_grady(double *, double *,int, double *, double *);
  extern double hb_p80(double, double);
  extern double hb_gravity(double );
  extern double hb_coriol(double );
  extern double hb_linterp(double, double *, double *,int);
#endif  /*ifndef PHYPROPS*/

/* prototypes for functions in eos80.c  */
#ifndef EOS80
#define EOS80 1
  extern double hb_eos80d(double, double, double, double **);
  extern double hb_alpha(double, double, double);
  extern double hb_beta(double, double, double);
  extern double hb_gamma(double, double, double);
  extern double hb_ratio(double, double, double);
#endif  /*ifndef EOS80*/

/* prototypes for functions in prop_subs.c  */
#ifndef PROP_SUBS
#define PROP_SUBS 1
 extern int get_prop_indx(char *);
 extern void print_prop_menu(); 
 extern char *get_prop_mne(int);
 extern char *get_prop_descrip(int);
 extern char *get_prop_units(int);
 extern int get_field_width(int);
 extern int get_field_precis(int);
 extern void compute_sigma(double, int, double *, double *, double *, double *);
 extern void compute_ratio( int , double *, double *, double *, double *, double *, double *);
 extern void compute_svan(int, double *, double *, double *, double *);
 extern void compute_sp_vol(int, double *, double *,  double *, double *);
 extern void compute_height(int, double *, double *, double *, double, double *);
 extern void compute_energy(int, double *, double *, double *, double, double *);
 extern void compute_theta(int, double *, double *, double *, double *);
 extern void buoy_freq(double *, double *, double *, double *, int, int, int);
 extern void po_vort(double *, double *, int, double);
 extern void compute_sound_vel(double *, double *,  double *, double *, int);
 extern double buoyancy(double , double *, double *, double *, int,int, int);
 extern double potvort(double, double);
 extern double ox_kg2l(double, double, double, double);
 extern double ox_l2kg(double, double, double, double);
#endif  /*ifndef PROP_SUBS*/
