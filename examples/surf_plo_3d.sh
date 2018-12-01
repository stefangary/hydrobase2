#!/bin/csh -f
#-----------------------------------------
# This script uses hydrobase and gmt
# routines to read a hydrobase data set
# and plot one of:
#  (1) pressure
#  (2) potential vorticity
# on a density surface along with
# distributions of the number of points
# and standard deviations.
#
# The name of the hydrobase dataset to
# use for the plot is $1 and the name
# of the property to plot (in hydrobase
# syntax) is $2.  The output plots will
# be named with $2 as a root name.
#
# This routine is specifically tailored
# to run with hydrobase 3D CDF files
# rather than station files.  Since the
# 3D files are already fully interpolated,
# the blockmean and surface fitting
# routines are skipped.
#
# Stefan Gary, 2018
#-----------------------------------------
# This software is distributed under the
# terms of the GNU LGPL v3 or later.
#-----------------------------------------

set infile = $1
set pname = $2

set topofile = topo.grd.keep
set topolevs = topo.cpt.keep

# Define clipping contour line
echo 2000 c > contint.txt

# Define colorscales for the plots
set pvmax = 30
set prmax = 300
if( $pname == 'pr' ) then
    set vmax = ${prmax}
endif

if( $pname == 'pv' ) then
    set vmax = ${pvmax}
endif
echo $vmax

#------------sigma_2 = 36.95-------------
#makecpt -Crainbow -T1000/2200/50 > pr1.cpt
#makecpt -Crainbow -T2225/2500/25 > pr2.cpt
#------------sigma_2 = 36.98-------------
makecpt -Crainbow -T1000/2400/50 > pr1.cpt
makecpt -Crainbow -T2425/2700/25 > pr2.cpt
makecpt -Crainbow -T0/200/10 > pr_std.cpt
makecpt -Crainbow -T0/20/1 > pv1.cpt
makecpt -Crainbow -T30/100/10 > pv2.cpt
makecpt -Crainbow -T0/20/1 > pv_std.cpt

# Set a temporary filename for the
# postscript output.
set out_mean = mean.ps
set out_std = std.ps
set out_num = num.ps

# Define grid of the output plot
@ xmin = -85
@ xmax = 0
@ xinc = 1

@ ymin = 0
@ ymax = 65
@ yinc = 1

# Other preferences
gmtset LABEL_FONT_SIZE 12
gmtset ANNOT_FONT_SIZE_PRIMARY 12

# Interpolate properties vertically to surface.
# Colums in output file are:
# $1   $2   $3
# lon  lat  prop
hb_gridsurf2d $infile -C -P${pname} -B${xmin}/${xmax}/${ymin}/${ymax} -I1.0 -W100/10 <<ENDLIST
s2 36.98 surf.grid
end
ENDLIST

# Convert surface file .grid to .xyz and
# do some smoothing in the process.
hb_grid2xyz surf.grid -Osurf -P${pname} -BNA_s2_3695.blk.keep -S -R1 -W0.25 -M9 

# Remove points that have less than ten measurments.
#awk '{if ($4 > 9) {print $1, $2, $3}}' surf${pname}.xyz > tmp.xyz
#awk '{if ($4 > 9) {print $1, $2, $3}}' surf${pname}_var.xyz > tmp_var.xyz
#mv tmp.xyz surf${pname}.xyz
#mv tmp_var.xyz surf${pname}_var.xyz

# Truncate lines that have HUGE standard deviations.
awk -v vmax="$vmax" '{if($3 <= vmax) {print $1, $2, $3} else {print $1, $2, vmax}}' surf${pname}_var.xyz > tmp.xyz
mv tmp.xyz surf${pname}_var.xyz

# One further averaging:
blockmean surf${pname}.xyz -R${xmin}/${xmax}/${ymin}/${ymax} -I2.0 -F -V > surf${pname}.bkm
blockmean surf${pname}_var.xyz -R${xmin}/${xmax}/${ymin}/${ymax} -I2.0 -F -V > surf${pname}_var.bkm
awk '{print $1, $2}' surf${pname}.bkm > pos.xy

# Interpolate missing features and convert to grd.
#---------------USING SURFACE------------------
surface surf${pname}.bkm -Gmean.grd -I2.0 -R${xmin}/${xmax}/${ymin}/${ymax} -T0.5i0b -Vl

surface surf${pname}_var.bkm -Gstd.grd -I1.0 -R${xmin}/${xmax}/${ymin}/${ymax} -T0.5i0b -Vl

#=============================================
# Contour plot mean field
#=============================================
# Lay down deep bathymetry shading
grdimage $topofile -C${topolevs} -JM6i -Ba10/a10WeSn -R${xmin}/${xmax}/${ymin}/${ymax} -X1i -Y2i -V -K -P > $out_mean

# Define clipping paths
grdcontour $topofile -Ccontint.txt -J -Dclip -m > cliplog.ps

# Convert clipping contours from .xyz to .xy
awk '{print $1, $2}' clip > clip.xy
rm clip

# Apply the clipping paths
psclip clip.xy -J -B -R -X0i -Y0i -B -m -P -O -K -V >> $out_mean

# Contour plot data level 1.
grdcontour mean.grd -C${pname}1.cpt -J -B -R -Wthin,black -P -Q -V -O -K -Gn1/2i -A+ap >> $out_mean

# Contour plot data level 2.
grdcontour mean.grd -C${pname}2.cpt -J -B -R -Wthin,black,-- -P -Q -V -O -K -Gn1/2i -A+ap >> $out_mean

# Plot data point locations
psxy pos.xy -J -R -B -Gred -Wred -O -K -P -Sp0.05 -X0i -Y0i -V >> $out_mean

# Turn off plot clipping
psclip -C -O -K -X0i -Y0i >> $out_mean

# Lay down intermediate, shallow bathymetry shading.
#grdimage $topofile -C${topoul} -J -B -R -X0i -Y0i -V -O -K -P -G:B- >> $out_mean

# Plot coastlines
pscoast -V -R -J -B -P -Gblack -O >> $out_mean

# Plot colorbar
#psscale -D3.5i/-1i/6i/0.5ih -C${pname}.cpt -O >> $out_mean

# Convert to pdf.
set name_mean = (`basename $out_mean .ps`)
ps2pdf ${name_mean}.ps
mv ${name_mean}.pdf ${pname}_mean.pdf

#=============================================
# Shade plot std field with mean countours
#=============================================
# Lay down bathymetry shading
grdimage $topofile -C${topolevs} -JM6i -Ba10g1/a10g1WeSn -R${xmin}/${xmax}/${ymin}/${ymax} -X1i -Y2i -V -K -P > $out_std

# Apply the clipping paths
psclip clip.xy -J -B -R -X0i -Y0i -B -m -P -O -K -V >> $out_std

# Shade plot std
grdimage std.grd -C${pname}_std.cpt -J -B -R -P -Q -X0i -Y0i -V -O -K >> $out_std

# Contour plot data level 1.
grdcontour mean.grd -C${pname}1.cpt -J -B -R -Wthin,black -P -Q -V -O -K -S5 -Gn1/2i -A+ap >> $out_std

# Contour plot data level 2.
grdcontour mean.grd -C${pname}2.cpt -J -B -R -Wthin,black,-- -P -Q -V -O -K -S5 -Gn1/2i -A+ap >> $out_std

# Shade plot data point locations
psxy pos.xy -J -R -B -Gblack -O -K -P -Sp -X0i -Y0i -V >> $out_std

# Lay down intermediate, shallow bathymetry shading.
#grdimage $topofile -C${topoul} -J -B -R -X0i -Y0i -V -O -K -P >> $out_std

# Turn off plot clipping
psclip -C -O -K -X0i -Y0i >> $out_std

# Plot coastlines
pscoast -V -R -J -B -P -Gblack -O -K >> $out_std

# Plot colorbar
psscale -D3.5i/-1i/6i/0.5ih -C${pname}_std.cpt -O >> $out_std

# Convert to pdf.
set name_std = (`basename $out_std .ps`)
ps2pdf ${name_std}.ps
mv ${name_std}.pdf ${pname}_std.pdf

# Clean up
rm *.ps
rm *.cpt
rm *.xyz
rm *.bkm
rm *.blk
rm *.grd
rm *.xy
rm contint.txt
rm *.grid
