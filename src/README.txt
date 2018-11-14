List of major changes compared to original HB2:

- fixed bug in hb_interpolate.  HB3 has replaced this
  with hb_prseries whose interpolation core has been
  rewritten.

- Increased the size of MAXSTDLEVS and also required
  that std_depths be initialized to MAXSTDLEVS.  Before,
  the two variables were not linked.  All changes in
  hydro_cdf.h

- Added comments, hello debugging to hb_fit3d.
  Only the depth levels reading in was added.

- Added IH variable (17/12/2015)
  > hydrobase.h  MAX_PROP = 41 (from 40)
                 IH added in list after HT.
  > prop_subs.c  added function compute_htdz_over_f
    		 changed -99999.9 to HB_MISSING in function compute_height
		 added GAP_CHANGE_DEPTH, GAP_SHALLOW, GAP_DEEP to
		 compute_height and compute_htdz_over_f.
  > hb_propcalc.c added ih_pref variable and case for IH using the
    		  function added in prop_subs.c.
  > hb_gridsurf2d.c added ih_pref variable and cases for IH.

