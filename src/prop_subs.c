/* prop_subs.c
................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             1993
                             Updated to ANSI standards:  Nov 1999
................................................................................
.
.  Functions which facilitate use of properties listed in "hydrobase.h"
.
.  int get_prop_indx(char *s) :
.           translates character string into an integer
.           representing an index to the list of properties.
.  void print_prop_menu() : 
.           prints a list of available properties on the stderr device.
.
.  char *get_prop_mne(int i)  :
.           returns the mnemonic associated with property index i.
.
.  char *get_prop_descrip(int i) :
.           Returns a char description of property[i].
.
.  char *get_prop_units(int i) :
.           Returns the units associated with property[i].
.
.  int get_field_width(int i) :
.           Returns the field width associated with property[i].
.
.  int get_field_precis(int i) :
.           Returns the decimal precision associated with property[i].
.
.  void compute_sigma(double pref, int nobs, double *sigma, double *p, 
.                     double *t, double *s):
.           Computes an array of sigma values for each level in p,t,s arrays.
.
.  void compute_dref_sigma(double *d, int nobs, double *sigma, double *p,
.                     double *t, double *s):
.           Sigma at s0, s2, and s4 depending on depth for each level
.
.  void compute_heat_capacity(int nobs, double *heat_cap, double *p,
.                     double *t, double *s):
.            Heat capacity, ref. to surface for each level
.
.  void compute_ratio( int nobs, double *ratio, double *p, double *t, double *s)
.  	 Computes density ratio = -ds/dt from  arrays of p, t, s at each observation
.    	level.     
.	
.  void compute_svan(int nobs, double *sva, double *p, 
.                      double *t, double *s).
.           Computes an array of specific volume anomaly values for each level 
.           in p,t,s arrays.
.
.  void compute_sp_vol( int nobs, double *sv, double *p,  double *t, double *s)
.           Computes an array of spec. volume values (in situ) for each level  
.            in p,t,s arrays.
.
.  void compute_height(int nobs, double *p, double *t, double *s, double pref, double *ht) :
.           Computes dynamic height relative to pref for each level in the
.           p, t, s arrays.
.
.  void compute_htdz_over_f(int nobs, double *p, double *t, double *s, double pref, double *htdz) :
.           Computes depth-integrated dynamic height over f relative to pref for each level in the
.           p, t, s arrays.
.
.  void compute_energy(int nobs, double *p, double *t, double *s, double pref, double *chi)
.           Computes potential energy relative to reference level pref for each
.           level in the p,t,s arrays.
.
.  void compute_theta(int nobs, double *th, double *p, double *t, double *s);
.           Computes potential temp rel to the sea surface for each level 
.           in the p, t, s arrays.
.
. void compute_dtdp(double *tp, double *p, double *t, double *s, 
.                 int  nobs,  int window, int w_incr);
.
.           Computes dT/dp, the vertical gradient of temperature in pressure
.           coordinates.
.
. void compute_dtdz(double *tz, double *p, double *t, double *s, 
.                 int  nobs,  int window, int w_incr, double lat);
.
.           Computes dT/dz, the vertical gradient of temperature in z-level
.           (depth) coordinates.
.
. void buoy_freq(double *bf, double *p, double *t, double *s, 
.                 int nobs,  int window, int w_incr);
.           Computes brunt-vaisala frequency in radians/sec at each pressure 
.           level, for a set of p,t,s observations. Window and w_incr together 
.           specify the pressure interval and how that interval is divided up
.           for computing the specific vol anomaly gradient near each observed
.           pressure. 
.
. void po_vort(double *pv, double *e, int nobs, double lat)
.           Computes potential vorticity for each value of the stability 
.           parameter (n^2 in (radians/sec)^2) and latitude.
.
.  void compute_sound_vel( double *sv, double *p,  double *t, double *s, int nobs)
.           Computes an array of sound velocity values (in m/sec) for each level  
.            in p,t,s arrays.
.
. double buoyancy(double p0, double *p, double *t, double *s, 
.                 int  nobs,  int window, int w_incr);
.           Computes buoyancy one value at a time.
.           Computes brunt-vaisala frequency in radians/sec at pressure, 
.           p0 given a set of p,t,s observations. The window specifies 
.           the pressure interval surrounding p0 in which observations 
.           will be considered; w_incr is the increment by which the window
.           is subdivided into smaller pieces to approximate the gradient.
.
. double potvort(double e, double lat)
.           Computes PV one value at a time.
.           Returns potential vorticity for the specified value of stability 
.           parameter,e, (= N^2 in (radians/sec)^2) and latitude.
.
. double ox_kg2l(double oxkg, double pr, double te, double sa)
.           Returns oxygen in ml/l from ox in umole/kg.
.
. double ox_l2kg(double oxl, double pr, double te, double sa)
     Returns oxygen in umole/kg from ox in ml/l.
.
*/
/*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PROP_SUBS 1
#include "hydrobase.h"
#define PI 3.141592654

/*************************************************************************/
/*************************************************************************/

/* character string mnemonics which correspond to properties */

char *prop_mne[MAXPROP] =  { "pr",   
                             "de",  
                             "te",   
                             "th",
			     "tp",
			     "tz",
                             "hc",   
                             "sa",   
                             "ox",   
                             "o2",   
                             "n2",
                             "n3",
                             "p4",
                             "si",  
                             "ht",
			     "ih",
                             "pe",  
                             "s0",  
                             "s1",  
                             "s2", 
                             "s3", 
                             "s4",
                             "s_",
			     "sd",
                             "bf",  
                             "pv",
                             "sv",
                             "va",   
                             "f1",
                             "f2",
                             "f3",
                             "he",
                             "tu",
                             "gn",
                             "ge",
                             "vn",
                             "ve",
			     "vs",
			     "dr",
			     "al",
			     "be",
			     "rr"
};

char *prop_descrip[MAXPROP] = { "pressure",
                                "depth",
                                "in situ temperature",
                                "potential temperature: pref=0.",
				"vert. grad. pot. temp. pref=0.",
				"vert. grad. pot. temp. pref=0.",
                                "heat capacity: pref = 0.",
                                "salinity",
                                "oxygen",
                                "oxygen",
                                "nitrite",
                                "nitrate",
                                "phosphate",
                                "silicate",
                                "dynamic height",
				"de-integrated dynamic height/f",
                                "potential energy anomaly",
                                "potential density: pref = 0.",
                                "potential density: pref = 1000.",
                                "potential density: pref = 2000.",
                                "potential density: pref = 3000.",
                                "potential density: pref = 4000.",
                                "potential density: pref = ?.",
				"potential density: s0, s2, s4.",
                                "buoyancy frequency",
                                "potential vorticity",
                                "specific volume",
                                "specific volume anomaly",
                                "cfc-11",
                                "cfc-12",
                                "cfc-113",
                                "helium",
                                "tritium",
                                "neutral density (gamma-n)",
                                "gamma-n errorbar",
				"velocity north",
				"velocity east",
				"sound velocity",
				"density ratio",
				"thermal expansion",
				"haline contraction",
				"APPROX. Rossby radius of def."
} ;

char *prop_units[MAXPROP] = { "dbars",
                              "meters",
                              "degrees C",
                              "degrees C",
			      "degrees C/dbars",
			      "degrees C/meters",
                              "J/(kg degC)",
                              "psu",
                              "ml/liter",
                              "micromole/kg",
                              "micromole/kg",
                              "micromole/kg",
                              "micromole/kg",
                              "micromole/kg",
                              "dyn. meters (= *10 m**2/s**2)",
			      "Sverdrups (= *10**6 m**3/s)",
                              "10**6 g/sec**2 (ergs/cm**2)",
                              "kg/m**3",
                              "kg/m**3",
                              "kg/m**3",
                              "kg/m**3",
                              "kg/m**3",
                              "kg/m**3",
                              "kg/m**3",
                              "* 1.e-5 radians/sec",
                              "* 1.e-12 m^-1 sec^-1",
                              "* 1.e-8 m**3/kg",
                              "* 1.e-8 m**3/kg",
                              "picomole/kg",
                              "picomole/kg",
                              "picomole/kg",
                              "micromole/kg",
                              "tritium units",
                              "kg/m**3",
                              "kg/m**3",
                              "m/sec",
                              "m/sec",
                              "m/sec",
			      "",
			      "10**7 alpha",
			      "10**7 beta",
			      "meters"
};

int field_width[MAXPROP] = {     
                                 8, /* pr */
                                 8, /* de */
                                 8, /* te */
                                 8, /* th */
				10, /* tp */
				10, /* tz */
                                10, /* hc */
                                 8, /* sa */
                                 8, /* ox */
                                 8, /* o2 */
                                 8, /* n2 */
                                 8, /* n3 */
                                 8, /* p4 */
                                 8, /* si */
                                 8, /* ht */
				10, /* ih */
                                10, /* pe */
                                 8, /* s0 */
                                 8, /* s1 */
                                 8, /* s2 */
                                 8, /* s3 */
                                 8, /* s4 */
                                 8, /* s_ */
				 8, /* sd */
                                10, /* bf */
                                10, /* pv */
                                10, /* sv */
                                 8, /* va */
                                 8, /* f1 */
                                 8, /* f2 */
                                 8, /* f3 */
                                 8, /* he */
                                 8, /* tu */
                                 8, /* gn */
                                 8, /* ge */
                                 9, /* vn */
                                 9, /* ve */
                                 9,  /* vs */
				 8,  /* dr */
				 9, /* al */
				 9, /* be */
				 9  /* rr */
};

int field_precis[MAXPROP] =  {   
                                 1, /* pr */
                                 1, /* de */
                                 4, /* te */
                                 4, /* th */
                                 4, /* tp */
                                 4, /* tz */
				 3, /* hc */
                                 4, /* sa */
                                 3, /* ox */
                                 1, /* o2 */
                                 3, /* n2 */
                                 3, /* n3 */
                                 3, /* p4 */
                                 3, /* si */
                                 3, /* ht */
				 4, /* ih */
                                 1, /* pe */
                                 4, /* s0 */
                                 4, /* s1 */
                                 4, /* s2 */
                                 4, /* s3 */
                                 4, /* s4 */
                                 4, /* s_ */
				 4, /* sd */
                                 4, /* bf */
                                 3, /* pv */
                                 3, /* sv */
                                 4, /* va */
                                 3, /* f1 */
                                 3, /* f2 */
                                 3, /* f3 */
                                 3, /* he */
                                 3, /* tu */
                                 4, /* gn */
                                 5, /* ge */
				 4, /* vn */
				 4, /* ve */
				 3, /* vs */
				 4, /* dr */
				 1, /* al */
				 1, /* be */
				 1  /* rr */
};

/******************************************************/

void print_prop_menu()
{
   int i;

   for (i=0; i < MAXPROP; ++i) {
     fprintf(stderr,"\n%s : %s [%s]", prop_mne[i], prop_descrip[i], prop_units[i]);
   }
   fprintf(stderr,"\n");
   return;
} /* end print_prop_menu() */

/******************************************************/

int get_prop_indx(char *str)
{
   int error = -1;
   char *s;

   s = str;
   switch (*s) {
      case 'a':
      case 'A':
            switch (*(++s)) {
               case 'L':
               case 'l':
                   return ((int) AL);  
               default:      
                   return (error);
             } 
            break;
      
      case 'b':
      case 'B':
            switch (*(++s)) {
               case 'E':
               case 'e':
                   return ((int) BE);  
               case 'F':
               case 'f':
                   return ((int) BF);  
               default:      
                   return (error);
             } 
            break;
      case 'd':
      case 'D':
            switch (*(++s)) {
               case 'E':
               case 'e':
                   return ((int) DE);  
               case 'R':
               case 'r':
                   return ((int) DR);  
               default:      
                   return (error);
             } 
            break;
      case 'F':
      case 'f':
            switch (*(++s)) {
               case '1':
                   return ((int) F1);
               case '2' :
                   return ((int) F2);  
               case '3' :
                   return ((int) F3);  
               default:      
                   return (error);
             } 
            break;
      case 'G':
      case 'g':
            switch (*(++s)) {
               case 'E':
               case 'e' :
                   return ((int) GE);  
               case 'N':
               case 'n' :
                   return ((int) GN);  
               default:      
                   return (error);
             } 
            break;
      case 'H':
      case 'h':
            switch (*(++s)) {
               case 'C':
               case 'c':
                   return ((int) HC);  
               case 'E':
               case 'e':
                   return ((int) HE);  
               case 'T':
               case 't':
                   return ((int) HT);  
               default:      
                   return (error);
             } 
            break;
      case 'I':
      case 'i':
	    switch (*(++s)) {
	       case 'H':
	       case 'h':
	           return ((int) IH);
	       default:
	           return (error);
	    }
            break;
      case 'N':
      case 'n':
            switch (*(++s)) {
              case '2':
                 return ((int) N2);
              case '3':
                 return ((int) N3);
              default:
                 return error;
            }
            break;
      case 'O':
      case 'o':
            switch (*(++s)) {
               case 'X':
               case 'x':
                   return ((int) OX);  
               case '2':
                   return ((int) O2);  
              default:      
                   return (error);
             } 
            break;
      case 'P':
      case 'p':
            switch (*(++s)) {
               case 'E':
               case 'e':
                   return ((int) PE);
               case 'R':
               case 'r':
                   return ((int) PR);
               case 'V':
               case 'v':
                   return ((int) PV);  
               case '4':
                   return ((int) P4);  
               default:      
                   return (error);
             } 
            break;
      case 'R':
      case 'r':
            switch (*(++s)) {
               case 'R':
               case 'r':
                   return ((int) RR);  
              default:      
                   return (error);
             } 
            break;
      case 'S':
      case 's':
            switch (*(++s)) {
              case 'A':
              case 'a':
                 return ((int) SA);
	      case 'D':
	      case 'd':
		 return ((int) SD);
              case 'I':
              case 'i':
                 return ((int) SI);
              case 'V':
              case 'v':
                 return ((int) SV);
              case '_':
                 return ((int) S_);
              case '0':
                 return ((int) S0);
              case '1':
                 return ((int) S1);
              case '2':
                 return ((int) S2);
              case '3':
                 return ((int) S3);
              case '4':
                 return ((int) S4);
              default:
		return (error);
            } /* end switch */
            break;
      case 'T':
      case 't':
            switch (*(++s)) {
              case 'E' :
              case 'e' :
                 return ((int) TE);
              case 'H' :
              case 'h' :
                 return ((int) TH);
	      case 'P' :
	      case 'p' :
		 return ((int) TP);
              case 'U' :
              case 'u' :
                 return ((int) TU);
	      case 'Z' :
	      case 'z' :
		 return ((int) TZ);
              default:
		return (error);
            } /* end switch */
            break;
      case 'V':
      case 'v':
            switch (*(++s)) {
              case 'A':
              case 'a':
                 return ((int) VA);
               case 'E':
               case 'e' :
                   return ((int) VE);  
               case 'N':
               case 'n' :
                   return ((int) VN);  
               case 'S':
               case 's' :
                   return ((int) VS);  
              default:
		return(error);
            } /* end switch */
            break;

      default:
	return(error);
   } /* end switch */
} /* end get_prop_indx() */

/*****************************************************************/
char *get_prop_mne(int i)
/*  returns the ith property mnemonic */
{
   return prop_mne[i];
}
/*****************************************************************/
char *get_prop_descrip(int i)
   /* returns the ith property description */
{
   return prop_descrip[i];
}
/*****************************************************************/
int get_field_width(int i)
   /* returns the ith property field width */
{
   return field_width[i];
}
/*****************************************************************/
int get_field_precis(int i)
   /* returns the decimal precision used to describe property i */
{
   return field_precis[i];
}
/*****************************************************************/
char *get_prop_units(int i)
    /* returns the units for property[i] */
{
   return prop_units[i];
}
/*****************************************************************/
void compute_sigma(double pref, int nobs, double *sigma, double *p, double *t, double *s)
  /* Computes sigma values from  arrays of p, t, s at each observation
    level.     */
{
   int    j;
   double tref, sv;

      for (j = 0; j < nobs; ++j) {
         if (s[j] < 0 || p[j] < 0)
	    sigma[j] = HB_MISSING;
	 else {
            tref = hb_theta(s[j], t[j], p[j], pref);
            sv = hb_svan(s[j], tref, pref, &sigma[j]);
	 }
      }

      return;
}  /* end compute_sigma() */
/*****************************************************************/

void compute_dref_sigma(double *d, int nobs, double *sigma, double *p, double *t, double *s) {
  /* Computes sigma values from  arrays of d, p, t, s
   * at each observation level.  The depth referenced
   * sigma is a sigma profile computed using sigma_0,
   * sigma_2, and sigma_4 above 1000m, between 1000m
   * and 3000m and below 3000m, respectively.  This
   * is the same sigma calculation as is done for
   * sorting scans into sigma bins during the quality
   * control routines (hb_statfit_ts and hb_statchk_ts). */
  int    j;
  double tref, pref, sv;

  for (j = 0; j < nobs; ++j) {
    if (s[j] < 0.0 || t[j] < -8.0 || p[j] < 0.0)
      sigma[j] = HB_MISSING;
    else {

      /* Determine ref pressure:
       * Depths less than 1km = 0 dbar
       * Depths 1km to 3km = 2000 dbar
       * Depths greater 3km = 4000 dbar */
      pref = 2000.0;
      if (d[j] <= 1000.0)
	pref = 0.0;
      if (d[j] > 3000.0)
	pref = 4000.0;
      
      tref = hb_theta(s[j], t[j], p[j], pref);
      sv = hb_svan(s[j], tref, pref, &sigma[j]);
    }
  }
  return;
}  /* end compute_dref_sigma() */

/*****************************************************************/

void compute_heat_capacity( int nobs, double *heat_cap, double *p, double *t, double *s) {
  /* Compute specific heat capacity referenced at
   * the surface given pressure, insitu temp, and salinity
   */

   int    j;
   double tref;
   double pref = 0.0;

   for (j = 0; j < nobs; ++j) {
     if (s[j] < 0.0 || t[j] < -8.0 || p[j] < 0.0)
       heat_cap[j] = HB_MISSING;
     else {
       tref = hb_theta(s[j], t[j], p[j], pref);
       heat_cap[j] = hb_cpsw(s[j], tref, pref);
     }
   }
   return;
} /* End of compute_heat_capacity */

/*****************************************************************/

void compute_ratio( int nobs, double *ratio, double *p, double *t, double *s, double *alpha, double *beta)
  /* Computes density ratio = -ds/dt, thermal expansion (alpha) and haline contraction (beta) coefficients  from  arrays of p, t, s at each observation level. Set ratio, alpha, or beta to NULL if either of these quantities are not returned   */
{
   int    i, j;
   double **DRV;
   double SV;
   
   DRV = (double **) calloc(3, sizeof(double *));
   for (i= 0; i < 3; ++i)
       DRV[i] = (double *) calloc(8, sizeof(double));
       
      for (j = 0; j < nobs; ++j) {
         if (s[j] < 0 || p[j] < 0) {
	    if (alpha != NULL)
	         alpha[j] = HB_MISSING;
	    if (beta != NULL)	 
	         beta[j] = HB_MISSING;
	    if (ratio != NULL)	 
	        ratio[j] = HB_MISSING;
	 }
	 else {
	 
	    SV = hb_eos80d(s[j], t[j], p[j], DRV);
	    
	    if (alpha != NULL)
	         alpha[j] = -1.0e7 * DRV[1][2] / (DRV[0][2] + 1000.0);
	    if (beta != NULL)	 
	         beta[j] = 1.0e7 * DRV[0][7] / (DRV[0][2] + 1000.0);
	    if (ratio != NULL)	 
	        ratio[j] = -DRV[1][2] / DRV[0][7]  ;
	 }
      }

   for (i = 0; i < 3; ++i)       
       free( DRV[i]); 
   free(DRV);
   return;
}  /* end compute_ratio() */

/****************************************************************************/
void compute_sp_vol(int nobs, double *spvol, double *p, double *t, double *s)
  /*  Computes in situ Specific Volume from p, t, s at each observation
    level.     */
{
   int    j;
   double  sigma, sv;

   /* Units are  * 1e-8 m**3/kg 
      Convert the sigma value (kg/m**3) to a rho value by adding 1000.
      Specific volume is 1/rho, but is multiplied here by 10^8 to 
      convert to same units as sp.vol. anomaly. */
 
      for (j = 0; j < nobs; ++j) {
         if (s[j] < 0 || p[j] < 0)
	    spvol[j] = HB_MISSING;
	 else {
           sv = hb_svan(s[j], t[j], p[j], &sigma);
           spvol[j] = 1.0e8 / (sigma + 1000.);
	 }
      }

      return;

}  /* end compute_sp_vol() */
/****************************************************************************/
void compute_svan(int nobs, double *sva, double *p, double *t, double *s)
  /* Computes specific volume anomaly ( = spvol(p,t,s) - spvol(p,0,35)
     from  arrays of p, t, s at each observation level.  
     Units are in 1.0e-8 m**3/kg 
               or 1.0e-5 cm**3/g  */
{
   int    j;
   double sigma;

      for (j = 0; j < nobs; ++j) {
         if (s[j] < 0 || p[j] < 0)
	    sva[j] = HB_MISSING;
	 else
            sva[j] = hb_svan(s[j], t[j], p[j], &sigma);
      }

      return;
}  /* end compute_svan() */
/****************************************************************************/
void compute_height(int nobs, double *p, double *t, double *s, double pref, double *h)

/* computes dynamic height relative to pref at each level in pressure
   array.  The units are :  
             dynamic height  : dyn meters = 1/10 * m**2/s**2 .  
           spec vol anomaly  :  1e-8 * m**3/kg
                   pressure  :  dbars = 1e4 N/m**2 = 1e4 kg/m s**2

   If a vertical datagap is encountered, 
    no height is computed beneath that level

SFG changed -99999.9 to HB_MISSING and generalized the gap detection
to work with the GAP_SHALLOW and GAP_DEEP global variables.
 */
{
  int j, start, i, datagap;
  double sref, tref;
  double sva1, sva0, sig, last_h, last_p;
  

/* first find the reference pressure level ... */

  /* ...by finding the temperature, tref at pref...*/
   if ((tref = hb_linterp(pref, p, t, nobs)) < -8.) {
      for (i = 0; i < nobs; ++i) {
         h[i] = HB_MISSING;
      }
      return;
   }

   /* ...and by finding the salinity, sref, at pref...*/
   sref = hb_linterp(pref, p, s, nobs);

   start = 0;
   while (p[start] < pref)
      ++start;
      
   if (start > 0)
      --start;  /* start now points to the first pr level above the ref pr */

   /* check for datagap around the ref level */
      
   datagap = (p[start+1] - p[start]) > GAP_DEEP;
   if (datagap) {
     for (i = 0; i < nobs; ++i) {
        h[i] = HB_MISSING;   /* can't compute height for this station */
     }
     return;
   }
   

/* determine height between the ref lev and the first pr level above...*/   
   
   sva0 = hb_svan(sref, tref, pref, &sig);
   h[start] = 0.0;
   if (start >= 0) {
      sva1 = hb_svan(s[start], t[start], p[start], &sig);
      h[start] = (sva0 + sva1) *.5e-5 * (pref - p[start]);
      sva0 = sva1;
   }
   last_h = h[start];
   last_p = p[start];

/* now integrate upward through the station, check for missing values in p,t,s
   and for vertical datagaps ... */

   for (j = start-1; j >= 0; --j) {
      if (s[j] < -8. || t[j] < -8. || p[j] < -8) {
         h[j] = HB_MISSING;
         continue;
      }
      if (p[j] < GAP_CHANGE_DEPTH)
        datagap = ( last_p - p[j]) > GAP_SHALLOW;
      else
        datagap = ( last_p - p[j]) > GAP_DEEP;

      if (datagap) {
        for (i= j; i >= 0; --i) {
          h[i] = HB_MISSING;
        }
        j = 0;  /* don't bother integrating upward any farther */
      }
      else {
        sva1 = hb_svan(s[j], t[j], p[j], &sig);

/* the 1e-5 term corrects the units : 10^-8  * 10^4 * 10^-1 
                                      (svan)   Pa/db   dyn meters   */
        h[j] = last_h + ((sva0 + sva1) * 0.5e-5 * (last_p - p[j]));
        last_h = h[j];
        last_p = p[j];
        sva0 = sva1;
      }
   } /* end for */

   
 /* find ht between ref level and the first observation beneath that level... */
  
   ++start;   /* start now points to first level at or below ref pr */
   sva0 = hb_svan(sref, tref, pref, &sig);
   h[start] = 0.0;
   if (start < nobs) {
      sva1 = hb_svan(s[start], t[start], p[start], &sig);
      h[start] = (sva0 + sva1) * 0.5e-5 * (pref - p[start]) ;
      sva0 = sva1;
   }
   last_h = h[start];
   last_p = p[start];
  
/* now integrate downward through the station, check for missing values 
    in p,t,s and for vertical datagaps ... */
    
   for (j = start+1; j < nobs; ++j) {
      if (s[j] < -8. || t[j] < -8. || p[j] < -8) {
         h[j] = HB_MISSING;
         continue;
      }

      if (p[j] < GAP_CHANGE_DEPTH )
        datagap = (  p[j] - last_p) > GAP_SHALLOW;
      else
        datagap = (  p[j] - last_p) > GAP_DEEP;

      if (datagap) {
        for (i= j; i < nobs; ++i) {
          h[i] = HB_MISSING;
        }
        return;  /* don't bother integrating downward any farther */
      }
      
      sva1 = hb_svan(s[j], t[j], p[j], &sig);

      h[j] = last_h + ((sva0 + sva1) *.5e-5 * (last_p - p[j]) );
      last_h = h[j];
      last_p = p[j];
      sva0 = sva1;
   }

   return;
      
} /* end compute_height() */
/****************************************************************************/
void compute_htdz_over_f(int nobs, double *d, double *p, double *t, double *s, double pref, double lat, double *htdz)

/* computes depth integrated dynamic height over f relative to pref at each level in pressure
   array.  The units are :  
             dynamic height  : dyn meters = 1/10 * m**2/s**2 .  
           spec vol anomaly  :  1e-8 * m**3/kg
                   pressure  :  dbars = 1e4 N/m**2 = 1e4 kg/m s**2

   Copied from compute_height since dynamic height is needed, and then we need to integrate that.  Difference is that we integrate by depth, not pressure.

   If a vertical datagap is encountered, 
    no height is computed beneath that level
 */
{
  int j, start, i, datagap;
  double sref, tref, href, dref, f;
  double sva1, sva0, sig, last_h, last_htdz, last_d, last_p;
  double *h;

  /*===============================================
   * First find the reference pressure level.
   * If the pressure cannot be found in the profile,
   * flag all output values to -99999.9 and quit.*/
   if ((tref = hb_linterp(pref, p, t, nobs)) < -8.) {
      for (i = 0; i < nobs; ++i) {
         htdz[i] = HB_MISSING;
      }
      return;
   }

   /* Find salinity at the reference pressure.*/
   sref = hb_linterp(pref, p, s, nobs);

  /*================================================
   * Use compute_height to find the dynamic height first time.*/
   h = (double *)malloc(sizeof(double)*nobs); 
   compute_height(nobs,p,t,s,pref,h);

  /* Find the dynamic height at the reference pressure.
   * Check value, should be zero. Find the depth at the
   * reference pressure too since we have to integrate
   * by dz, not dp. */
   href = hb_linterp(pref, p, h, nobs);
   dref = hb_linterp(pref, p, d, nobs);

  /*=================================================
   * Set up the integration */
   start = 0;
   while (p[start] < pref)
     ++start;
      
   if (start > 0)
     --start;  /* start now points to the first pr level above the ref pr */

   /* Check for datagap around the ref level
    * If data gap too big, flag all output to
    * -99999.9 and quit. This is usually due to
    * the fact that the reference level can be
    * much deeper than stations over shallow
    * topography.*/   
   datagap = (p[start+1] - p[start]) > GAP_DEEP ;
   if (datagap) {
     for (i = 0; i < nobs; ++i) {
       htdz[i] = HB_MISSING;   /* can't compute height for this station */
       /*fprintf(stderr,"\ncompute_htdz_over_f WARNING: gap at ref. level!");*/
     }
     return;
   }
   
/* determine height between the ref lev and the first pr level above.
 * The 0.5 factor is a 1/2 to average href and h[start]
 * (necessary for the trapezoidal integration).*/   
   htdz[start] = 0.0;
   if (start >= 0) {
      htdz[start] = (href + h[start])*0.5*(dref - d[start]);
      last_h = h[start];
   }
   last_htdz = htdz[start];
   last_d = d[start];
   last_p = p[start];

/* now integrate upward through the station, check for missing values
   and for vertical datagaps ... */

   for (j = start-1; j >= 0; --j) {
      if (s[j] < -8. || t[j] < -8. || p[j] < -8. || h[j] < -8. || d[j] < -8. ) {
         htdz[j] = HB_MISSING;
	 /*fprintf(stderr,"\ncompute_htdz_over_f WARNING: non-sane values!");*/
         continue;
      }

      if (p[j] < GAP_CHANGE_DEPTH)
        datagap = ( last_p - p[j]) > GAP_SHALLOW;
      else
        datagap = ( last_p - p[j]) > GAP_DEEP;

      if (datagap) {
        for (i= j; i >= 0; --i) {
          htdz[i] = HB_MISSING;
	  fprintf(stderr,"\ncompute_htdz_over_f WARNING: gap detected while integrating up!");
        }
        j = 0;  /* don't bother integrating upward any farther */
      }
      else {
        htdz[j] = last_htdz + ((last_h + h[j])*0.5*(last_d - d[j]));
        last_h = h[j];
	last_htdz = htdz[j];
        last_d = d[j];
	last_p = p[j];
      }
   } /* end for */

   /*==========================================================
    * Done integrating upward, now need to integrate downward.*/
   
   /* find ht between ref level and the first observation beneath that level... */
  
   ++start;   /* start now points to first level at or below ref pr */
   htdz[start] = 0.0;
   if (start < nobs) {
      htdz[start] = (href + h[start])*0.5*(dref - d[start]) ;
      last_h = h[start];
   }
   last_htdz = htdz[start];
   last_p = p[start];
   last_d = d[start];

/* now integrate downward through the station, check for missing values 
    in p,t,s and for vertical datagaps ... */
    
   for (j = start+1; j < nobs; ++j) {
      if (s[j] < -8. || t[j] < -8. || p[j] < -8. || h[j] < -8. || d[j] < -8. ) {
         htdz[j] = HB_MISSING;
	 /*fprintf(stderr,"\ncompute_htdz_over_f WARNING: non-sane values!");*/
         continue;
      }

      if (p[j] < GAP_CHANGE_DEPTH)
        datagap = (  p[j] - last_p) > GAP_SHALLOW;
      else
        datagap = (  p[j] - last_p) > GAP_DEEP;

      if (datagap) {
        for (i= j; i < nobs; ++i) {
          htdz[i] = HB_MISSING;
	  fprintf(stderr,"\ncompute_htdz_over_f WARNING: gap detected while integrating down!");
        }
        return;  /* don't bother integrating downward any farther */
      }
      htdz[j] = last_htdz + ((last_h + h[j])*0.5*(last_d - d[j]) );
      last_htdz = htdz[j];
      last_p = p[j];
      last_h = h[j];
      last_d = d[j];
   }

   /* Final step: divide all values by Coriolis parameter.
    * Can change f to 1 or other orders of magnitude for
    * checking the accuracy of the integration.  Divide by
    * 1e6 to get Sv as well as multiply out the extra 1/10 in
    * the dynamic meters by default.*/
   f = hb_coriol(lat);
   /*fprintf(stderr,"\n Coriolis parameter is %20.10e\n",f);*/
   for (i = 0; i < nobs; ++i) {
     if ( htdz[i] != HB_MISSING ) {
       htdz[i] = htdz[i]/(1e5*f);
       /*htdz[i] = htdz[i];*/
       /*fprintf(stderr,"\nData dump: %d %f",i,htdz[i]);*/
     }
   }

   return;
      
} /* end compute_htdz_over_f() */

/****************************************************************************/
void compute_approx_rossby_radius(double *ro, int nobs, int pdr, double *d, double *p, double *t, double *s, double lat, int window, int w_incr)

/* computes the Rossby radius of deformation using equations 2.2 and 2.3ab
 * in Chelton et al. (1998), which is an APPROXIMATE method and, on average,
 * biased high by about 6.5% according to their results.  This is
 * copied from compute_htdz_over_f since
 * the core of this calculation is depth integrated buoyancy frequency.
 * 
 * Note that a key difference in the integration here
 * compared to dynamic height integrations is that we
 * are only interested in integrating from the bottom
 * to the surface, not from a specified reference depth.
 * This is why there is no pref in the function.
 *
 * The units are :  
 * Rossby radius, ro: meters
 * spec vol anomaly  :  1e-8 * m**3/kg
 * pressure  :  dbars = 1e4 N/m**2 = 1e4 kg/m s**2
 *
 * If a vertical datagap is encountered, 
 * no height is computed beneath that level
 */
{
  int j, start, i, datagap;
  double f, beta;
  double sva1, sva0, sig, last_n, last_ro, last_d, last_p;
  double *n;

  /* Check that pdr is deeper than the deepest measurement. */
  if (d[nobs-1] > ((double) pdr)) {
    fprintf(stderr,"\ncompute_approx_rossby_radius ERROR: Deepest obs deeper than seafloor.");
    fprintf(stderr,"\ncompute_approx_rossby_radius ERROR: hb_update_pdr needs to be run first.");
    /* Return missing values. */    
    for (i = 0; i < nobs; ++i) {
      ro[i] = HB_MISSING;
    }
    return;
  }
  
  /*===============================================
   * First find the depth of the water column and
   * the depth of the deepest measurement.  If the
   * discrepancy is bigger than a data gap, we return
   * a bunch of missing values.  An alternative could
   * be to project uniform properties to the bottom,
   * but we hold off on that extrapolation for now.
   * Or, we could just ignore any data gaps and allow
   * for the integration to continue, regardless.
   * However, ignoring data gaps can result in
   * a few really bad values (i.e. a shallow profile
   * in deep water that happens to have very high
   * stratification.  The buoyacy frequency that will be
   * used in this case is based only on the upper
   * waters, not the waters at depth, which will be
   * completely unrepresentative of the integrated
   * value along the whole profile.
   */

  if (d[nobs-1] < GAP_CHANGE_DEPTH){
    datagap = ( pdr - d[nobs-1]) > GAP_SHALLOW;
    fprintf(stderr,"\ncompute_approx_rossby_radius WARNING: GAP_SHALLOW detected at bottom!");
  }
  else{
    datagap = ( pdr - d[nobs-1]) > GAP_DEEP;
    fprintf(stderr,"\ncompute_approx_rossby_radius WARNING: GAP_DEEP detected at bottom!");
  }

  if (datagap) {
    for (i = 0; i < nobs; ++i) {
      ro[i] = HB_MISSING;
    }
    return;
  }   
  
  /*================================================
   * Compute the buoyancy frequency
   */
   n = (double *)malloc(sizeof(double)*nobs); 
   buoy_freq(n,p,t,s,nobs,window,w_incr);

   /* Force the removal of any density inversions
    * by changing any negative values to zero.
    * Note, if we just check for < 0, it will
    * lose memory of any values set to HB_MISSING,
    * so keep those. */
   for (j = start-1; j >= 0; --j) {
     if ( n[j] < 0.0 ) {
       if ( n[j] == HB_MISSING ) {
       /* Leave the HB_MISSING value */
       }
       else{
	 n[j] = 0.0;
       }
     }
   }

   /* If the bottom value of the profile is missing,
    * we cannot start the integration, so exit. */
   if ( n[nobs-1] == HB_MISSING ) {
     fprintf(stderr,"\ncompute_approx_rossby_radius WARNING: HB_MISSING detected at deepest observation!");
      /* Return missing values. */    
      for (i = 0; i < nobs; ++i) {
        ro[i] = HB_MISSING;
      }
      return;
   }
   
  /*=================================================
   * Set up the integration
   */

   /* We always start from the bottom.*/
   start = nobs-1;
   
   /* Determine height between the bottom and the first pr level above.
    * The 0.5 factor is a 1/2 to average the value at the first
    * observation level and the bottom.  We assume N = 0 at the
    * bottom.  This is trapezoidal integration.*/   
   ro[start] = 0.0;
   if (start >= 0) {
      ro[start] = (0.0 + n[start])*0.5*(pdr - d[start]);
   }
   last_ro = ro[start];
   last_n = n[start];
   last_d = d[start];
   last_p = p[start];

   /* now integrate upward through the station, check for missing values
      and for vertical datagaps ... */

   for (j = start-1; j >= 0; --j) {
     /* Any negative values of n[j] here are HB_MISSING,
      * not density inversions as this was explicitly checked above.
      * Note that if this HB_MISSING value is inserted, the last_
      * values are NOT updated meaning that the integration will
      * just jump over this missing value, provided that the data
      * gap is not too big.*/
      if (s[j] < -8. || t[j] < -8. || p[j] < -8. || n[j] < -8. || d[j] < -8. ) {
         ro[j] = HB_MISSING;
	 /*fprintf(stderr,"\ncompute_approx_rossby_radius WARNING: non-sane value!");*/
         continue;
      }

      if (p[j] < GAP_CHANGE_DEPTH)
        datagap = ( last_p - p[j]) > GAP_SHALLOW;
      else
        datagap = ( last_p - p[j]) > GAP_DEEP;

      if (datagap) {
        /* Allow INTERNAL data gaps for now. */
        ro[j] = HB_MISSING;

	/* REMOVE THE ABOVE LINE AND ACTIVATE THE
	 * lines below for data gap filtering.
        for (i= j; i >= 0; --i) {
          ro[i] = HB_MISSING;
	  fprintf(stderr,"\ncompute_approx_rossby_radius WARNING: gap detected while integrating up!");
        }
        j = 0;  /* don't bother integrating upward any farther */
      }
      else {
        ro[j] = last_ro + ((last_n + n[j])*0.5*(last_d - d[j]));
        last_n = n[j];
	last_ro = ro[j];
        last_d = d[j];
	last_p = p[j];
      }
   } /* end for */

   /*==========================================================
    * Done integrating upward, no need to integrate downward.
    */
   
   /*==============================================================
    * Sanity check - if surface value is less than zero, there was
    * a massive density inversion.  Report missing values.
    */
   if ( ro[0] < 0 ) {
     fprintf(stderr,"\ncompute_approx_rossby_radius WARNING: Massive density inversion detected!");
     for (i = 0; i < nobs; ++i) {
       ro[i] = HB_MISSING;
     }
     return;
   }
   
   /*==========================================================
    * Final step: divide all values by Coriolis parameter
    * or beta plane, depending on the latitude.
    * Can manually change f to 1 or other orders of magnitude for
    * checking the accuracy of the integration.
    */
   if ( lat > 5. ) {
     /* Use Eq. 2.3a */
     f = hb_coriol(lat);
     /*fprintf(stderr,"\n Coriolis parameter is %20.10e\n",f);*/

     for (i=0; i < nobs; ++i) {
       if ( ro[i] != HB_MISSING ) {
	 ro[i] = ro[i]/(f*PI);
       }
     }
     
   }
   else {
     /* Use Eq. 2.3b */
     beta = hb_beta_plane(lat);
     /*fprintf(stderr,"\n Beta plane parameter is %20.10e\n",beta);*/

     for (i=0; i < nobs; ++i) {
       if ( ro[i] != HB_MISSING ) {
	 ro[i] = sqrt(ro[i]/(2*beta));
       }
     }
   }

   return;
      
} /* end compute_approx_rossby_radius() */

/****************************************************************************/
void compute_energy(int nobs, double *p, double *t, double *s, double pref, double *chi)

/* computes potential energy anomaly relative to the specified ref level at 
   each level in pressure array.  
        chi = 1/g * Integral of (Pr * sp_vol_anom * delta-P)
   The units are :  
           potential energy  : * 1e6 ergs/cm**2 (or g/sec**2) = 1e3 J/m**2
           spec vol anomaly  : * 1e-5 cm**3/g
                   pressure  :  dbars = 1e5 dyn/cm**2 = 1e5 g/cm/s**2
                          g  : 980 cm/s**2
                          
          
      value of constant comes from :
              correcting spvolanom units: 1e-5
              getting mean spvolanom:    .5
              correcting delta-pr units: 1e5
              getting mean pr:           .5
              correcting pr units:       1e5
              1/g:                       .00102040816
              adjust decimal place:      1e-6

    If a vertical datagap is encountered, the integration is stopped and
    no potential energy is computed above (beneath) that level.
 */
{
   int j, start, i, datagap;
   double sref, tref;
   double sva1, sva0, sig, last_chi, last_p;
   double C = .255102e-4;  /*constant defined above */

/* first find the reference pressure level ... */

   if ((tref = hb_linterp(pref, p, t, nobs)) < -8.) {
      for (i = 0; i < nobs; ++i) {
         chi[i] = -99999.9;
      }
      return;
   }
   
   sref = hb_linterp(pref, p, s, nobs);
   
   start = 0;
   while (p[start] < pref)
      ++start;
      
   --start;  /* start now points to the first pr level above the ref pr */

   /* check for datagap around the ref level */
      
   datagap = (p[start+1] - p[start]) > 600 ;
   if (datagap) {
     for (i = 0; i < nobs; ++i) {
        chi[i] = -99999.9;   /* can't compute energy for this station */
     }
     return;
   }
   
/* determine chi between the ref lev and the first pr level above...*/   
   
   sva0 = hb_svan(sref, tref, pref, &sig);
   chi[start] = 0.0;
   if (start >= 0) {
      sva1 = hb_svan(s[start], t[start], p[start], &sig);
      chi[start] = (sva0 + sva1) * (pref - p[start]) * (pref + p[start]) * C;
      sva0 = sva1;
   }
   last_chi = chi[start];
   last_p = p[start];

/* now integrate upward through the station, check for missing values in p,t,s
   and for vertical datagaps ... */

   for (j = start-1; j >= 0; --j) {
      if (s[j] < -8. || t[j] < -8. || p[j] < -8) {
         chi[j] = -99999.9;
         continue;
      }

      if (p[j] < 1001.)
        datagap = ( last_p - p[j]) > 250;
      else
        datagap = ( last_p - p[j]) > 600;

      if (datagap) {
        for (i= j; i >= 0; --i) {
          chi[i] = -99999.9;
        }
        j = 0;  /* don't bother integrating upward any farther */
      }
      else {
        sva1 = hb_svan(s[j], t[j], p[j], &sig);

            /* C is a constant defined above */
        chi[j] = last_chi + ((sva0 + sva1) * (last_p - p[j])*(p[j]+ last_p) * C);
        last_chi = chi[j];
        last_p = p[j];
        sva0 = sva1;
      }
   }
   
 /* find chi between ref level and the first observation beneath that level... */
  
   ++start;   /* start now points to first level at or below ref pr */
   sva0 = hb_svan(sref, tref, pref, &sig);
   chi[start] = 0.0;
   if (start < nobs) {
      sva1 = hb_svan(s[start], t[start], p[start], &sig);
      chi[start] = (sva0 + sva1) * (pref - p[start]) * (pref + p[start]) * C;
      sva0 = sva1;
   }
   last_chi = chi[start];
   last_p = p[start];
  
/* now integrate downward through the station, check for missing values 
    in p,t,s and for vertical datagaps ... */
    
   for (j = start+1; j < nobs; ++j) {
      if (s[j] < -8. || t[j] < -8. || p[j] < -8) {
         chi[j] = -99999.9;
         continue;
      }

      if (p[j] < 1001.)
        datagap = (  p[j] - last_p) > 250;
      else
        datagap = (  p[j] - last_p) > 600;

      if (datagap) {
        for (i= j; i < nobs; ++i) {
          chi[i] = -99999.9;
        }
        return;  /* don't bother integrating downward any farther */
      }
      
      sva1 = hb_svan(s[j], t[j], p[j], &sig);

            /* C is a constant defined above */
      chi[j] = last_chi + ((sva0 + sva1) * (last_p - p[j])*(p[j]+ last_p) * C);
      last_chi = chi[j];
      last_p = p[j];
      sva0 = sva1;
   }
   return;
      
} /* end compute_energy() */

/***************************************************************************/
void compute_theta(int n, double *th, double *p, double *t, double *s)
/* int n;       # of levels in arrays 
   double *th;  pointer to array of theta values 
   double *p;   pointer to array of pressure 
   double *t;   pointer to array of in situ temperature 
   double *s;   pointer to array of salinity 
*/
{
   int i;
   double pref = 0.0;

   for (i = 0; i < n; ++i) {
       if (s[i] < 0.0 || p[i] < 0.0)
         th[i] = HB_MISSING;
       else
         th[i] = hb_theta(s[i], t[i], p[i], pref);
   }
   return;
}
/***************************************************************************/
void compute_dtdp(double *tp, double *p, double *t, double *s, int nobs, int window, int w_incr)
/*  double *tp;         array to hold returned vertical temperature
                        gradient values.
    double *p, *t, *s;  pressure and (insitu) temperature values 
    int nobs;           number of observations in p, t, s, tp arrays
    int window;         size of  pressure window (db) into which
                        observations will be incorporated for estimating t
                        gradients
    int w_incr;         subdivide the window into increments (db) for
                          approximating the gradients 
*/

/* Computes the vertical temperature gradient in deg. C/decibar using gradients 
 * near each observed pressure. A series of p,t points is generated at w_incr
 * increments from the shallowest observed pressure to the deepest. A temperature
 * gradient is then estimated over the specified pressure window at each
 * increment of the pressure series.  These values of temperature gradient are 
 * interpolated to obtain a value at each pressure level of the original
 * observations.
 * NOTE: Input temperature is assumed to be insitu temperature, not potential
 * temperature, which is consistent with all of the other hydrobase conventions.
 */
{
  double *pr_win, *te_win, *sa_win, *tp_win, *th_win, p0;
   double  e, bflast;
   int mid, npts_win, npts;
   int i, j, k, n;

/* Determine maximum size of arrays to hold pressure
 * series for a window.  The first step is to test
 * if the given w_incr is bigger than window.  Then,
 * we find the total pressure range for this profile,
 * (p[nobs-1] - p[0]), and divide that by the pressure
 * increment to get the number of points in the 
 * profile (e). Finally, allocate space.*/

   window = window / 2;
   if (w_incr > window)
       w_incr = window;
   
   e = (p[nobs-1] - p[0]) / (double) w_incr;
   npts = (int) (e + .00001);
   ++npts; /* add one for top point */
   ++npts; /* add one for bottom point */
   
   pr_win = (double *) malloc(npts * sizeof(double));
   te_win = (double *) malloc(npts * sizeof(double));
   sa_win = (double *) malloc(npts * sizeof(double));
   th_win = (double *) malloc(npts * sizeof(double));
   tp_win = (double *) malloc(npts * sizeof(double));
   if (pr_win == NULL) {
     fprintf(stderr, "\nUnable to allocate memory in compute_dtdp()\n");
     exit(2);
   }
   if (te_win == NULL) {
     fprintf(stderr, "\nUnable to allocate memory in compute_dtdp()\n");
     exit(2);
   }
   if (sa_win == NULL) {
     fprintf(stderr, "\nUnable to allocate memory in compute_dtdp()\n");
     exit(2);
   }   
   if (th_win == NULL) {
     fprintf(stderr, "\nUnable to allocate memory in compute_dtdp()\n");
     exit(2);
   }
   if (tp_win == NULL) {
     fprintf(stderr, "\nUnable to allocate memory in compute_dtdp()\n");
     exit(2);
   }

/* Construct a pressure series of t at the
 * specified increment.  Note that these profiles
 * are linearly interpolated from the observations.
 * The whole profile is dividied up into pressure
 * increments.  A pressue and temperature is 
 * associated with each increment.  If the
 * increments are consistently smaller than the
 * data spacing, then we are not really adding any
 * extra information because of the simple linear
 * interpolation between data points.
 */

   pr_win[0] = p[0];
   te_win[0] = t[0];
   sa_win[0] = s[0];
   th_win[0] = hb_theta(s[0],t[0],p[0],0.0);

   --npts;  /* to avoid subtracting each time test is performed in the for loop */
   for (i = 1; i < npts; ++i) {
      /* Incrementally step down in pressure. */
      pr_win[i] = pr_win[i-1] + w_incr;

      /* Linearly interpolate t to this pressure level. */
      te_win[i] = hb_linterp(pr_win[i], p, t, nobs);

      /* Linearly interpolate s to this pressure level. */
      sa_win[i] = hb_linterp(pr_win[i], p, s, nobs);
 
      /* Get the potential temperature at this level. */
      th_win[i] = hb_theta(sa_win[i],te_win[i],pr_win[i],0.0);
   }
   
   /* Set the final points at the bottom. */
   pr_win[npts] = p[nobs-1];
   te_win[npts] = t[nobs-1];
   sa_win[npts] = s[nobs-1];
   th_win[npts] = hb_theta(sa_win[npts],te_win[npts],pr_win[npts],0.0);

   ++npts;  /* Back to correct npts */

/* Compute index of midpoint of window. For the default
 * -W100/10, then mid = 50/10 = 5.  Recall that window
 * has been halved, above. Window and w_incr are integers,
 * so this is integer division! The biggest w_incr can be
 * is window, so mid is at least 1 point large. */
   mid = window / w_incr;
   npts_win = mid * 2 + 1; /* define a window of 5 points width minimum */
   
/* Compute temperature gradient at each point in the
 * pressure series sliding the pressure window 
 * incrementally through.  Adjust the window at the top 
 * of the profile to use all the points that are available.
 * This will decrease to a minimum of 2 points at the top
 * (where i = 0). Do not do the 
 * same at the bottom because it leads to unusually high 
 * gradients.
 */

/* Unclear why n needs to be initialized here. 
 * Perhaps we set the minimum value just in
 * case it is not set below. */ 
   n = 2;

/* Loop 1: Go from the surface to the middle of the
 * uppermost possible window.  Specify the index of the
 * middle of the uppermost possible window as j, the
 * loop limit.  Also verify that if the middle of this
 * window is more than the number of available points
 * in the profile, set to the number of available points
 * in the profile.  So this first loop is for just the
 * surface points of a tall profile or it takes care of
 * all the points in a very short profile.  In this
 * case, tall and short are defined as relative to the
 * size of the window over which the calculations take
 * place.  Therefore, this loop is a special case loop
 * taking care of the surface boundary condition on the
 * gradient calculation.
 */
   j = mid;
   if (mid > npts)
     j = npts;
   for (i = 0; i < j; ++i) {
     /* Set number of points to use. */
     n = i * 2 + 1;

     /* Top of profile, so set window size to minimum of 2. */
     if (i == 0)
        n = 2;

     /* Check that number of points to use does
      * not exceed the number of (interpolated)
      * positions in the profile.*/
     if (n >= npts)
        n = npts;

     /* Compute the vertical gradient starting from 0
      * (the first element in the p,t arrays) to only
      * the first n points of these arrays. The gradient
      * is stored in o_win and the avg pressure over the
      * interval is in p0 and average temp is in e. */
     tp_win[i] = hb_grady(&th_win[0],&pr_win[0],n,&p0,&e);
   }

/* Loop 2: For all the windows centered below the upper
 * most possible window, simply slide the window down and
 * compute the vertical temperature gradient over each
 * window.  Note that in this case, npts_win does not 
 * change.  Also, the loop goes from the midpoint of the
 * uppermost window (which was the lower limit of Loop 1,
 * above) to the midpoint of the lowermost window 
 * (j = npts-mid).  Since we are computing the gradient
 * within all the FULL windows, the data that are included
 * in the gradient calculation start at index k=0 in 
 * _win[k] and as the window slides deeper, k is augmented
 * so the shallower data are not included.  Therefore,
 * this loop is the main working loop in that most of the
 * data should be processed here.
 */
   n = npts_win;
   j = npts - mid;
   k = 0;
   for (i = mid; i < j; ++i) {
        tp_win[i] = hb_grady(&th_win[k],&pr_win[k],n,&p0,&e);
        ++k;
   }
   
   /* Set bottom of window to last computed o_win value */
   e = tp_win[i-1];   
   if (j > mid) {
     for (i = j; i < npts; ++i) 
        tp_win[i] = e;
   }
   
   /* Interpolate from the window uniform spacing profile
    * to the original observations depths. */
   for (i = 0; i < nobs; ++i) {
      tp[i] = hb_linterp(p[i],pr_win,tp_win,npts);
      /* No need for units changes */
   }
 
   free((void *)th_win);
   free((void *)te_win);
   free((void *)pr_win);
   free((void *)sa_win);
   free((void *)tp_win);
   return;

} /* end compute_dtdp() */

/***************************************************************************/

void compute_dtdz(double *tz, double *p, double *t, double *s, int nobs, int window, int w_incr, double lat)
/*  double *tz;         array to hold returned vertical temperature
                        gradient values.
    double *p, *t, *s;  pressure and (insitu) temperature values 
    int nobs;           number of observations in p, t, s, tz arrays
    int window;         size of  pressure window (db) into which
                        observations will be incorporated for estimating t
                        gradients
    int w_incr;         subdivide the window into increments (db) for
                          approximating the gradients 
    double lat;         latitude, degrees, necessary for depth calculations
*/

/* Computes the vertical temperature gradient in deg. C/meter using gradients 
 * near each observed pressure.
 * NOTE: Input temperature is assumed to be insitu temperature, not potential
 * temperature, which is consistent with all of the other hydrobase conventions.
 * Abbreviated commenting here - check out compute_dtdp for more details.
 */
{
  double *pr_win, *te_win, *sa_win, *tz_win, *th_win, *de_win, p0;
   double  e, bflast;
   int mid, npts_win, npts;
   int i, j, k, n;

/* Determine maximum size of arrays to hold pressure
 * series for a window.
 */
   window = window / 2;
   if (w_incr > window)
       w_incr = window;
   
   e = (p[nobs-1] - p[0]) / (double) w_incr;
   npts = (int) (e + .00001);
   ++npts; /* add one for top point */
   ++npts; /* add one for bottom point */
   
   pr_win = (double *) malloc(npts * sizeof(double));
   de_win = (double *) malloc(npts * sizeof(double));
   te_win = (double *) malloc(npts * sizeof(double));
   sa_win = (double *) malloc(npts * sizeof(double));
   th_win = (double *) malloc(npts * sizeof(double));
   tz_win = (double *) malloc(npts * sizeof(double));
   if (pr_win == NULL) {
     fprintf(stderr, "\nUnable to allocate memory in compute_dtdz()\n");
     exit(2);
   }
   if (de_win == NULL) {
     fprintf(stderr, "\nUnable to allocate memory in compute_dtdz()\n");
     exit(2);
   }
   if (te_win == NULL) {
     fprintf(stderr, "\nUnable to allocate memory in compute_dtdz()\n");
     exit(2);
   }
   if (sa_win == NULL) {
     fprintf(stderr, "\nUnable to allocate memory in compute_dtdz()\n");
     exit(2);
   }   
   if (th_win == NULL) {
     fprintf(stderr, "\nUnable to allocate memory in compute_dtdz()\n");
     exit(2);
   }
   if (tz_win == NULL) {
     fprintf(stderr, "\nUnable to allocate memory in compute_dtdz()\n");
     exit(2);
   }

/* Construct a pressure series at the
 * specified increment.
 */
   pr_win[0] = p[0];
   te_win[0] = t[0];
   sa_win[0] = s[0];
   th_win[0] = hb_theta(s[0],t[0],p[0],0.0);
   de_win[0] = hb_depth(p[0],lat);

   --npts;  /* to avoid subtracting each time test is performed in the for loop */
   for (i = 1; i < npts; ++i) {
      /* Incrementally step down in pressure. */
      pr_win[i] = pr_win[i-1] + w_incr;

      /* Linearly interpolate t to this pressure level. */
      te_win[i] = hb_linterp(pr_win[i], p, t, nobs);

      /* Linearly interpolate s to this pressure level. */
      sa_win[i] = hb_linterp(pr_win[i], p, s, nobs);
 
      /* Get the potential temperature at this level. */
      th_win[i] = hb_theta(sa_win[i],te_win[i],pr_win[i],0.0);

      /* Get the depth at this pressure level. */
      de_win[i] = hb_depth(pr_win[i],lat);
   }
   
   /* Set the final points at the bottom. */
   pr_win[npts] = p[nobs-1];
   te_win[npts] = t[nobs-1];
   sa_win[npts] = s[nobs-1];
   th_win[npts] = hb_theta(sa_win[npts],te_win[npts],pr_win[npts],0.0);
   de_win[npts] = hb_depth(pr_win[npts],lat);

   ++npts;  /* Back to correct npts */

/* Compute index of midpoint of window.*/
   mid = window / w_incr;
   npts_win = mid * 2 + 1; /* define a window of 5 points width minimum */
   
/* Compute temperature gradient at each point in the
 * pressure series sliding the pressure window 
 * incrementally through.
 */

/* Loop 1: Go from the surface to the middle of the
 * uppermost possible window.
 */
   n = 2;
   j = mid;
   if (mid > npts)
     j = npts;
   for (i = 0; i < j; ++i) {
     /* Set number of points to use. */
     n = i * 2 + 1;

     /* Top of profile, so set window size to minimum of 2. */
     if (i == 0)
        n = 2;

     /* Check that number of points to use does
      * not exceed the number of (interpolated)
      * positions in the profile.*/
     if (n >= npts)
        n = npts;

     /* Compute the vertical gradient. */
     tz_win[i] = hb_grady(&th_win[0],&de_win[0],n,&p0,&e);
   }

/* Loop 2: For all the windows centered below the upper
 * most possible window, simply slide the window down and
 * compute the vertical temperature gradient over each
 * window.
 */
   n = npts_win;
   j = npts - mid;
   k = 0;
   for (i = mid; i < j; ++i) {
        tz_win[i] = hb_grady(&th_win[k],&de_win[k],n,&p0,&e);
        ++k;
   }
   
   /* Set bottom of window to last computed o_win value */
   e = tz_win[i-1];   
   if (j > mid) {
     for (i = j; i < npts; ++i) 
        tz_win[i] = e;
   }
   
   /* Interpolate from the window uniform spacing profile
    * to the original observations depths. */
   for (i = 0; i < nobs; ++i) {
      tz[i] = hb_linterp(p[i],pr_win,tz_win,npts);
      /* No need for units changes */
   }
 
   free((void *)th_win);
   free((void *)te_win);
   free((void *)pr_win);
   free((void *)sa_win);
   free((void *)tz_win);
   free((void *)de_win);
   return;

} /* end compute_dtdz() */

/***************************************************************************/

void buoy_freq(double *bf, double *p, double *t, double *s,int nobs,int window, int w_incr)
/*  double *bf;         array to hold returned buoyancy values
    double *p, *t, *s;  observed pressure, temperature, salinity values 
    int nobs;           number of observations in p,t,s arrays
    int window;         size of  pressure window (db) into which
                        observations will be incorporated for estimating t,s
                        gradients
    int w_incr;         subdivide the window into increments (db) for
                          approximating the gradients 
*/

/* Computes buoyancy frequency in radians/sec using gradients of 
   specific vol anomaly near each observed pressure. A pressure series of 
   p,t,s points is generated at w_incr increments from the shallowest observed
   pressure to the deepest. A sp.vol.anom. gradient is then estimated over 
   the specified pressure window and buoyancy is computed at each
   increment of the pressure series.  These values of buoyancy are interpolated
   to obtain a value at each pressure level of the original observations.
*/
{
   double *p_win, *t_win, *s_win, *b_win, p0;
   double  e, bflast;
   double cph2rps = 0.001745329;   /* (2*pi)/3600 */
   int mid, npts_win, npts;
   int i, j, k, n;

/* Determine maximum size of arrays to hold pressure
 * series for a window.  The first step is to test
 * if the given w_incr is bigger than window.  Then,
 * we find the total pressure range for this profile,
 * (p[nobs-1] - p[0]), and divide that by the pressure
 * increment to get the number of points in the profile.
 * Finally, allocate space.*/

   window = window / 2;
   if (w_incr > window)
       w_incr = window;
   
   e = (p[nobs-1] - p[0]) / (double) w_incr;
   npts = (int) (e + .00001);
   ++npts;
   ++npts; /* add one more for bottom point */
   
   p_win = (double *) malloc(npts * sizeof(double));
   t_win = (double *) malloc(npts * sizeof(double));
   s_win = (double *) malloc(npts * sizeof(double));
   b_win = (double *) malloc(npts * sizeof(double));
   if (b_win == NULL) {
     fprintf(stderr, "\nUnable to allocate memory in buoy_freq()\n");
     exit(2);
   }

/* Construct a pressure series of t and s at the
 * specified increment.  Note that these profiles
 * are linearly interpolated from the observations. */

   p_win[0] = p[0];
   t_win[0] = t[0];
   s_win[0] = s[0];
   --npts;  /* to avoid subtracting everytime test is performed, in the for loop () */
   for (i = 1; i < npts; ++i) {
      /* Incrementally step down in pressure. */
      p_win[i] = p_win[i-1] + w_incr;

      /* Linearly interpolate t and s to these pressure levels. */
      t_win[i] = hb_linterp(p_win[i], p, t, nobs);
      s_win[i] = hb_linterp(p_win[i], p, s, nobs);
   }
   
   /* Set the final points at the bottom. */
   p_win[npts] = p[nobs-1];
   t_win[npts] = t[nobs-1];
   s_win[npts] = s[nobs-1];
   
   ++npts;  /* back to correct npts */

/* Compute bvfreq at each point in the pressure series sliding the pressure
   window incrementally through.  Adjust the window at the top 
   of the profile to use all the points that are available.  This will decrease
   to 2 points at the top. Do not do the same at the bottom because it
   leads to unusually high pv values  */

   mid = window / w_incr;  /* index of midpoint of window */
   npts_win = mid * 2 + 1; /* define a window of 5 points width minimum */
   
   n = 2;
   j = mid;
   if (mid > npts)
     j = npts;
   for (i = 0; i < j; ++i) {
     n = i * 2 + 1;
     if (i == 0)
        n = 2;
     if (n >= npts)
        n = npts;
     b_win[i] = hb_bvfrq(&s_win[0], &t_win[0], &p_win[0], n, &p0, &e);
   }
   n = npts_win;
   j = npts - mid;
   k = 0;
   for (i = mid; i < j; ++i) {
        b_win[i] = hb_bvfrq(&s_win[k], &t_win[k], &p_win[k], n, &p0, &e);
        ++k;
   }
   
   /* Set bottom of window to last computed b_win value */
   e = b_win[i-1];   
   if (j > mid) {
     for (i = j; i < npts; ++i) 
        b_win[i] = e;
   }
   
   /* Interpolate from the window uniform spacing profile
    * to the observations depths. */
   for (i = 0; i < nobs; ++i) {
      bf[i] = hb_linterp(p[i], p_win, b_win, npts);
      if (bf[i] > -999.)
         bf[i] *= cph2rps;   /* correct the units */
   }
 
   free((void *)b_win);
   free((void *)p_win);
   free((void *)t_win);
   free((void *)s_win);
   return;

} /* end buoy_freq() */

/****************************************************************************/
void po_vort( double *pv, double *e, int nobs, double lat)
/* Returns potential vorticity for the specified value of stability 
   parameter (n-squared) and latitude. 
    pv;    pointer to array of returned povo values  
    e;     pointer to array of n-squared (buoyancy freq)values in (rad/sec)^2  
    nobs;  dimension of pv and e arrays 
    lat;   latitude in degrees 
*/
{
   int i;
   for (i = 0; i < nobs; ++i) {
      
      if (e[i] < -9990.)       /* check for missing value flag ... */
        pv[i] =  -999.;
      else {
	/* Multiply by 1e14 = 1e12 * 100
	 * for output units of 1e-12 (ms)^-1
	 * for output units of
	 * AND divide by extra gravity term because
	 * e = g/rho * drho/dz but pv = f/rho * drho/dz */
        pv[i] = e[i] * 1.0e14  * hb_coriol(lat) / hb_gravity(lat);
     }
   }
   return ;

} /* end potvort() */


/***************************************************************************/
void compute_sound_vel(double *svel, double *p, double *t, double *s, int nobs)
  /*  Computes sound velocity from p, t, s at each observation
    level.     */
{
   int    j;
 
      for (j = 0; j < nobs; ++j) {
         if (s[j] < 0 || p[j] < 0)
	    svel[j] = HB_MISSING;
	 else 
           svel[j] = hb_svel(s[j], t[j], p[j]);
      }

      return;

}  /* end compute_sp_vol() */
/****************************************************************************/
double buoyancy(double p0, double *p, double *t, double *s, int nobs, int window, int w_incr)

/* Computes a single value of buoyancy frequency (radians/sec) using gradients of 
   t and s near the pressure designated by p0.  The t,s gradients are 
   estimated by linearly interpolating t and s as a function of p from 
   the observations which which fall within a pressure range of window
   on either side of p0.  A pressure series of p,t,s points is generated in
   this manner at w_incr intervals from p0.  If the p0 specified does 
   not exist in the pressure series represented by p, the value -9999 is
   returned.
   
   double p0;          pressure level at which buoyancy is computed 
   double *p, *t, *s;  observed pressure, temperature, salinity values 
   int nobs;           number of observations in p,t,s arrays 
   int window;         size of pressure interval (db) on either side of p0 in which
                      observations will be incorporated into estimating t,s
                      gradients
   int w_incr;         subdivide the window into increments (db) for
                      approximating the gradients 
*/
{
   double *p_win, *t_win, *s_win;
   double bv, e;
   double cph2rps = 0.001745329;   /* (2*pi)/3600 */
   int mid, npts_win;
   int i, n, end, start;

/* determine maximum size of arrays to hold pressure series for a window */

   window = window / 2;
   if (w_incr > window)
       w_incr = window;
   mid = window / w_incr;  /* index of midpoint of array */
   npts_win = mid * 2 + 1;
   
   p_win = (double *) malloc(npts_win * sizeof(double));
   t_win = (double *) malloc(npts_win * sizeof(double));
   s_win = (double *) malloc(npts_win * sizeof(double));

/* determine whether specified pressure level exists at this station */

   if (hb_linterp(p0, p, t, nobs) < -9990.) {
           return (-9999.);
   }

/* set up the window of pressure around the specifed pressure level...
   Account for cases where window juts above top observed pressure or 
   below bottom observed pressure ...  
   start represents the index at the top of window, end is
   the index at the bottom.  Keep the top and bottom as equidistant as
   possible from the midpoint because bvf will be computed at the mid-pressure
   of the window.*/

   p_win[mid] = p0;
   
   start = 0;
   i = mid;
   while (--i >= start) {
     if ((p_win[i] = p_win[i+1] - w_incr) < p[0] )
       start = i+1;
   }
   
   end = mid + mid - start;
   if (end == mid)
      ++end;
      
   i = mid;
   n = nobs-1;
   while (++i <= end) {
      if ((p_win[i] = p_win[i-1] + w_incr) > p[n])
        end = i-1;
   }
   
   if (start < mid)   
      start = mid - (end - mid);
   
   if (end == mid)
      start = mid - 1;
      
   npts_win = (end - start) + 1;

   for (i = start; i <= end; ++i) {
      t_win[i] = hb_linterp(p_win[i], p, t, nobs);
      s_win[i] = hb_linterp(p_win[i], p, s, nobs);
   }
   
   bv = hb_bvfrq(&s_win[start], &t_win[start], &p_win[start], npts_win, &p0, &e);

   free((void *)p_win);
   free((void *)t_win);
   free((void *)s_win);
   return (bv * cph2rps);

} /* end buoyancy() */
/****************************************************************************/
double potvort(double e, double lat)

/* Returns a single value of potential vorticity for the specified value of 
stability  parameter, e, (= buoyancy frequency squared) and latitude.
   
    e :    n-squared (buoyancy freq) in (rad/sec)^2 
   lat:   latitude in degrees
*/
{

   if (e < -9990.)       /* check for missing value flag ... */
       return (-99.);

    return (e * 1.0e14  * hb_coriol(lat) / hb_gravity(lat) );

} /* end potvort() */

/****************************************************************************/
double ox_kg2l(double ox, double pr, double te, double sa)
/*  Returns oxygen in units of ml/l 
   double ox:   oxygen value in micromoles/kg 
   double pr, te, sa:    in situ properties 
*/
{
    double x, pt;
    double sig, zero = 0.0;
    
    if (ox <= 0) 
       return ox;
    
    pt = hb_theta(sa, te, pr, zero);
    x = hb_svan(sa, pt, zero, &sig);
    return (ox * 0.022392 * (1 + sig /1000.));

}
/****************************************************************************/
double ox_l2kg(double ox, double pr, double te, double sa)
/*  Returns oxygen in units of umole/kg 
   double ox:   oxygen value in ml/l 
   double pr, te, sa:    in situ properties 
*/
{
    double x, pt;
    double sig, zero = 0.0;
    
    if (ox <= 0) 
       return ox;
    pt = hb_theta(sa, te, pr, zero);
    x = hb_svan(sa, pt, zero, &sig);
    return (ox / (0.022392 * (1 + sig /1000.)));

}
/****************************************************************************/
/****************************************************************************/
