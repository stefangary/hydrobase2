#!/bin/csh -f
#-----------------------------------------
# This script uses hydrobase and gmt
# routines to make a ts plot for the
# file named in $1 and writes the pdf
# output file to $2.
#
# Stefan Gary, 2018
#-----------------------------------------
# This software is distributed under the
# terms of the GNU LGPL v3 or later.
#-----------------------------------------


# Set a temporary filename for the
# postscript output.
set outps = out.ps

# Set the temp and salt ranges on the plot
set trange = -2/30
set tinc = 0.5
set srange = 30.0/40.0
set sinc = 0.1

# Make sigma basemap(s).
hb_sigbasemap -Osig0basemap.xyz -T/${trange}/${tinc} -S/${srange}/${sinc} -P0.0
hb_sigbasemap -Osig2basemap.xyz -T/${trange}/${tinc} -S/${srange}/${sinc} -P2000.0
hb_sigbasemap -Osig4basemap.xyz -T/${trange}/${tinc} -S/${srange}/${sinc} -P4000.0

# Define contour intervals.
makecpt -Crainbow -V -T0/50/0.5 > sig.cpt

# Other preferences
gmtset LABEL_FONT_SIZE 12
gmtset ANNOT_FONT_SIZE_PRIMARY 12

# Contour the basemap(s).
pscontour sig0basemap.xyz -Csig.cpt -R${srange}/${trange} -A+ap -JX7i/7i -P -Ba0.5g0.1f0.1:'Salinity [PSU]':/a1g1f0.1:'Potential Temperature [deg C]':WeSn -V -Wthin,green -X1i -Y1i -K -Gn1/2i > $outps

pscontour sig2basemap.xyz -Csig.cpt -R -A+ap -J -P -B -V -Wthin,red -O -K -Gn1/2i >> $outps

pscontour sig4basemap.xyz -Csig.cpt -R -A+ap -J -P -B -V -Wthin,cyan -O -K -Gn1/2i >> $outps

# Compute potential temperature
# from the input data and then
# write theta-S data to output
# file (xy) for plotting.  Use
# the -M option in hb_xyprop to
# place a ">" between profiles
# so GMT makes separate lines.
hb_propcalc $1 -Pth/sa | hb_xyprop -M -Xsa -Yth > ts.xy

# Overlay plot the T-S data to the
# basemap to make a T-S plot.
#psxy ts.xy -A -J -R -B -Wthick,gray -O -K -P -V -m >> $outps
psxy ts.xy -A -J -R -B -Wblack -Sp0.05 -O -P -V -m >> $outps

# Convert and clean up.
set fname = (`basename $outps .ps`)
ps2pdf ${fname}.ps
mv ${fname}.pdf $2
rm ${fname}.ps
rm sig*basemap.xyz
rm ts.xy
rm *.cpt
