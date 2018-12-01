#!/bin/csh -f
#---------------------------------------------
# This script will read a 3D HB CDF format
# climatology file ($1) and produce section
# plots of the the requested property ($2).
#
# Stefan Gary, 2018
#---------------------------------------------
# This software is distributed under the terms
# of the GNU LGPL v3 or later.
#---------------------------------------------

# Define the hydrobase cdf file to input.
set infile = $1

# Define the property we wish to plot
set pname = $2

# Define the topography file for plotting.
set topomask = ../lib/topo.onetenthdeg.swap.dat
set blkfile  = NA_s2_3695.blk.keep

# Define colorscales for the plots
makecpt -Crainbow -T1000/2200/50 > pr1.cpt
makecpt -Crainbow -T2225/2500/25 > pr2.cpt
makecpt -Crainbow -T0/20/1 > pv1.cpt
makecpt -Crainbow -T30/100/10 > pv2.cpt
makecpt -Crainbow -T36.90/37.00/0.01 > s2.cpt

# Set a temporary filename for the
# postscript output.
set out_zonal = zonal.ps
set out_merid = merid.ps

# Define grid of the output plot
@ xmin = -85
@ xmax = 0
@ xinc = 1

@ ymin = 0
@ ymax = 65
@ yinc = 1

# Define the MERIDIONAL section path
@ mlon = -60
@ mlati = 20
@ mlatf = 45
echo $mlon $mlati >> merid.sct
echo $mlon $mlatf >> merid.sct

# Define the ZONAL section path
@ zlat = 38
@ zloni = -75
@ zlonf = -5
echo $zloni $zlat >> zonal.sct
echo $zlonf $zlat >> zonal.sct

# Get section masking information
hb_toposlice -T${topomask} -I0.1 -M1 -Olonmask.xy -Lzonal.sct
hb_toposlice -T${topomask} -I0.1 -M2 -Olatmask.xy -Lmerid.sct

# Other preferences
gmtset LABEL_FONT_SIZE 12
gmtset ANNOT_FONT_SIZE_PRIMARY 12

#=============================================
# Contour plot MERIDIONAL section
#=============================================

# Interpolate properties vertically and horizontally
# to create a section.

if ( $pname == pr ) then
hb_slice $infile -I1 -Pde/${pname}/s2 -O${pname}.hb.tmp -Lmerid.sct
endif

if ( $pname == pv ) then
hb_slice $infile -I1 -Pde/pr/${pname}/s2 -O${pname}.hb.tmp -Lmerid.sct
endif

# Convert slice data to section data.
hb_section ${pname}.hb.tmp -Xla -Yde -Z${pname} > sect${pname}.xyz
hb_section ${pname}.hb.tmp -Xla -Yde -Zs2 > sects2.xyz

# Blockmean the data to 2 degree columns, 200 dbar intervals
blockmean sect${pname}.xyz -I2/200 -R${mlati}/${mlatf}/0/6000 -V > sect${pname}.bkm
blockmean sects2.xyz -I2/200 -R${mlati}/${mlatf}/0/6000 -V > sects2.bkm

# Get data point locations
awk '{print $1, $2}' sect${pname}.bkm > pos.xy

# Interpolate missing values 
surface sect${pname}.bkm -Gsect${pname}.grd -I1/100 -R${mlati}/${mlatf}/0/6000 -T0.5i0b -Vl
surface sects2.bkm -Gsects2.grd -I1/100 -R${mlati}/${mlatf}/0/6000 -T0.5i0b -Vl


# Contour the meridional section
grdcontour sect${pname}.grd -C${pname}1.cpt -JX6i/-6i -P -R${mlati}/${mlatf}/0/6000 -Ba5f1:'Latitude':/a500f100:'Depth [m]':WeSn -V -Wthin,black -X1i -Y1i -K -S5 -Gn1/2i > $out_merid

grdcontour sects2.grd -Cs2.cpt -J -P -R -B -V -Wthin,red -X0i -Y0i -O -K -S5 -Gn1/2i >> $out_merid

# Add data point locations
psxy pos.xy -J -R -B -Sp -O -K -P -V >> $out_merid

# Adding bathymetry polygon
psxy latmask.xy -A -J -R -B -Gblack -Wthick -L -O -P -V -m >> $out_merid

# Convert to pdf.
set name_merid = (`basename $out_merid .ps`)
ps2pdf ${name_merid}.ps
mv ${name_merid}.pdf ${pname}_merid.pdf

#=============================================
# Contour plot ZONAL section
#=============================================

# Interpolate properties vertically and horizontally
# to create a section.

if ( $pname == pr ) then
hb_slice $infile -I1 -Pde/${pname}/s2 -O${pname}.hb.tmp -Lzonal.sct
endif

if ( $pname == pv ) then
hb_slice $infile -I1 -Pde/pr/${pname}/s2 -O${pname}.hb.tmp -Lzonal.sct
endif

# Convert slice data to section data.
hb_section ${pname}.hb.tmp -Xlo -Yde -Z${pname} > sect${pname}.xyz
hb_section ${pname}.hb.tmp -Xlo -Yde -Zs2 > sects2.xyz

# Blockmean the data to 2 degree columns, 200 dbar intervals
blockmean sect${pname}.xyz -I2/200 -R${zloni}/${zlonf}/0/6000 -V > sect${pname}.bkm
blockmean sects2.xyz -I2/200 -R${zloni}/${zlonf}/0/6000 -V > sects2.bkm

# Get data point locations
awk '{print $1, $2}' sect${pname}.bkm > pos.xy

# Interpolate missing values 
surface sect${pname}.bkm -Gsect${pname}.grd -I1/100 -R${zloni}/${zlonf}/0/6000 -T0.5i0b -Vl
surface sects2.bkm -Gsects2.grd -I1/100 -R${zloni}/${zlonf}/0/6000 -T0.5i0b -Vl

# Contour the meridional section
grdcontour sect${pname}.grd -C${pname}1.cpt -JX6i/-6i -P -R${zloni}/${zlonf}/0/6000 -Ba5f1:'Longitude':/a500f100:'Depth [m]':WeSn -V -Wthin,black -X1i -Y1i -K -Gn1/2i > $out_zonal

grdcontour sects2.grd -Cs2.cpt -J -P -R -B -V -Wthin,red -X0i -Y0i -O -K -Gn1/2i >> $out_zonal

# Add data point locations
psxy pos.xy -J -R -B -Sp -O -K -P -V >> $out_zonal

# Adding bathymetry polygon
psxy lonmask.xy -A -J -R -B -Gblack -Wthick -L -O -P -V -m >> $out_zonal

# Convert to pdf.
set name_zonal = (`basename $out_zonal .ps`)
ps2pdf ${name_zonal}.ps
mv ${name_zonal}.pdf ${pname}_zonal.pdf

# Clean up
rm sect${pname}.xyz
rm *.cpt
rm latmask.xy
rm lonmask.xy
rm *.ps
rm *.bkm
rm *.sct
rm *.hb.tmp
rm *.blk
rm *.grd
