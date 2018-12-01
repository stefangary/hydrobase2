The scripts in this folder do work together as
a cohesive QC pipeline, however they are primarily
included here as a template as each particular QC
application may require fine tuning.

=================================================
DEPENDENCIES
=================================================

1) In order to use state binning in a time efficient
way, you will want to install GNU parallel, i.e.

sudo apt-get install parallel

...if you have access to more than 2 cores.

2) In order to automatically generate T-S plots
of the qc range and raw data points during the
QC process, you will need the Generic Mapping Tools v5+.

=================================================

The first two passes are during the
QC stage:
1) raw data are range checked
2) then state-binned
3) then ts-checked

The data in this main directory and not
in the archive_pass directories are the
state-binned QC'ed data that were then
used to construct the final climatologies.

The current script is parallelized, there
are two levels - the fe_ front end for setting
up the environment and the wk_ worker script
for the individual parallelized steps.  GNU
parallel seems to run fine on small runs,
but will crash if not enough processesors are
reserved for system use.

Performance increases are critical
now as the error bars estimate are based on leaving
out random profiles which will have an impact on the
state-binning process.  State binning is probably
the single most expensive part the climatology
construction.

=================================================

Notes on the resolution of fill values in hb_bin3d:

Use of hb_bin3d will result in lots of fill values,
especially temperatures.  I'm not clear as to why
the temperatures are filled but the salinities
are not.  However, looking at the locations of these
fills, they are of two broad categories:
1) along the bottom - due to the fact that hb_bin3d does
   not fill in the bottom observation in the HB-netcdf
   file but hb_nc2asc does attempt to pull that value
   out, so many nodes at the very bottom of the profile
   have temperature with HB_FILL.
2) HB_FILL in each profile - again, only te, might be due to
   vertical data gaps that don't allow interpolation.
The worrying observation is that many times, when there
is a temperature HB_FILL in the middle of the column,
all salinity and temperature values below that fill
value are nearly constant.  This means that there
are potentially EXTRAPOLATED data below the fill
value that are entirely bogus.

Might be an issue with using TE instead of T90.
Try it.

Even with standard approach - input TE and output
TE, data look OK.  The uniform salinities at the
bottom with slightly changing temperatures look
like they match the data - could be a precision
limit of HB for deep salinities.

Fill values in the raw data, usually in salinity
due to bubbles near the surface or glider apogee,
seem to be gracefully ignored by the binning.

The mystery is why is there a stretch of flagged
temperatures along a patch of good data from T~9.4
to T~8?  The salinities are OK (uniformly spaced),

YES! This was an issue with TE instead of T90.
When specifying output as T90, the spurious fill
values go away and a continuous profile is generated.
Looks like it is OK for input files to be TE, just
not the output.

Also, note that about 1/15 ~ 7% of scans differ
by exactly =/- 0.0001 oC (the limit of HB precision
and instruments as well) if input is TE or T90.
This must be due to a storage precision issue ->
if data are precomputed and stored from TE to T90
before running hb_bin3d, then the T90 values are
only to 0.0001.  When TE is converted to T90
inside hb_bin3d during the automatic conversion
routine, the T90 values are stored directly to
memory without being saved in a file and they have
full double precision, not just 0.0001.  Hence
the difference, exactly at the limit of precision.

Finally, noted that hb_bin3d doesn't deal gracefully
with input fill values (HB_MISSING).  Make certain
that HB_MISSING are filtered out with remove_fill_values.sh.

