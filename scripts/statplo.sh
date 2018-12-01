#!/bin/csh -f
#-----------------------------------------
# This script uses hydrobase and gmt
# routines to make a ts plot for the
# file named in $1, plots its statistical
# envelope defined by the root $2 and
# writes the pdf output file to $3.
#
# Stefan Gary, 2018.
#
# Mods for GMT5 - removed -m option
# in psxy.
# Autorange T-S plot domain
#-----------------------------------------
# This software is distributed under the terms
# of the GNU LGPL v3 or later.
#-----------------------------------------

# Set a temporary filename for the
# postscript output.
set outps = out.ps

# Set the temp and salt ranges on the plot

# New dynamic ranges depending on input data
set smin = `awk '$4 == "c" {print $1}' $1_ts.dat | gmtmath STDIN LOWER -Sl =`
set smax = `awk '$4 == "c" {print $1}' $1_ts.dat | gmtmath STDIN UPPER -Sl =`
set srange = ${smin}/${smax}
set sinc = `gmtmath -Q $smax $smin SUB 50 DIV =`
set tmin = `awk '$4 == "c" {print $2}' $1_ts.dat | gmtmath STDIN LOWER -Sl =`
set tmax = `awk '$4 == "c" {print $2}' $1_ts.dat | gmtmath STDIN UPPER -Sl =`
set trange = ${tmin}/${tmax}
set tinc = `gmtmath -Q $tmax $tmin SUB 50 DIV =`

# Old static ranges
#set trange = -2/30
#set tinc = 0.5
#set srange = 30.0/40.0
#set sinc = 0.1

# Make sigma basemap(s).
hb_sigbasemap -Osig0basemap.xyz -T/${trange}/${tinc} -S/${srange}/${sinc} -P0.0
hb_sigbasemap -Osig2basemap.xyz -T/${trange}/${tinc} -S/${srange}/${sinc} -P2000.0
hb_sigbasemap -Osig4basemap.xyz -T/${trange}/${tinc} -S/${srange}/${sinc} -P4000.0

# Define contour intervals.
makecpt -Crainbow -T18.0/26.5/0.5 > sig0.cpt
makecpt -Crainbow -T35.5/37.0/0.5 > sig2.cpt
makecpt -Crainbow -T46.0/50.0/0.5 > sig4.cpt
makecpt -Crainbow -T0/6000/100 > depth.cpt

# Other preferences
#gmtset LABEL_FONT_SIZE 12
#gmtset ANNOT_FONT_SIZE_PRIMARY 12

# Contour the basemap(s).
pscontour sig0basemap.xyz -Csig0.cpt -R${srange}/${trange} -A+ap -JX7i/7i -P -Ba0.2g0.1f0.05:'Salinity [PSU]':/a1g1f0.1:'Potential Temperature [deg C]':WeSn -Wthin,black -X1i -Y3i -K -Gn1/2i > $outps

pscontour sig2basemap.xyz -Csig2.cpt -R -A+ap -J -P -B -Wthin,black -O -K -Gn1/2i >> $outps

pscontour sig4basemap.xyz -Csig4.cpt -R -A+ap -J -P -B -Wthin,black -O -K -Gn1/2i >> $outps

# Overlay plot the T-S data to the
# basemap to make a T-S plot (includes
# both rejected and accepted scans).
# The first two columns are (x,y) = (s,t)
# while the next column is pressure
# which is mapped to the color axis.  The
# final column is a symbol designation for
# either accepted (c = circle) or rejected
# (+ = plus sign) scans.  The -S flag's value
# specifies the symbol size.
psxy $1_ts.dat -A -J -R -B -Cdepth.cpt -S0.08 -O -K -P >> $outps

# Overlay plot the mean T-S curve (black circles)
psxy ${2}.mean.ts -A -J -R -B -Sc0.1i -O -K -P >> $outps

# Overlay plot the lines tangent T-S curve (red)
psxy ${2}.slope.ts -A -J -R -B -Wthick,red -O -K -P >> $outps

# Overlay plot the envelope about the T-S curve (green)
psxy ${2}.stddev.ts -A -J -R -B -Wthick,green -O -K -P >> $outps

# Overlay the depth colorbar
psscale -D3.5i/-1i/6i/0.5ih -Cdepth.cpt -Ba1000f500g500 -O >> $outps

# Convert and clean up.
set fname = (`basename $outps .ps`)
ps2pdf ${fname}.ps
mv ${fname}.pdf $3
rm ${fname}.ps
rm sig*basemap.xyz
rm *.cpt
