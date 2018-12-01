#!/bin/tcsh -f
#---------------------------------------------
# This script will scan through a series of
# hydrobase stations, infiles, and use the
# suite of hb_statfit_ts and hb_statchk_ts
# to remove data points that do not fit
# within the statistical variability of the
# data at each density level.
#
# Stefan Gary, 2018
#
# This rountine is different from the first
# version of stat_check_basin (e.g. GLBB2012)
# because the previous version just used all
# ship and argo float data and QC'ed all at
# once in a data-density adaptive manner.
# This routine first uses state-binned
# isopycnal-averaged data to compute the
# statistical envelope of acceptable data
# so that gliders, ships, or argo floats density
# in a particular year/season don't overly
# bias the mean and acceptable ranges.
#
# Then, all data is QC'ed on a 1 degree by
# 1 degree grid (not data-density adaptive)
# because there's more or less enough data
# to do this.
#
# Data are not divided by seasons - all seasons
# and months QC'd together.
#---------------------------------------------
# This software is distributed under the terms
# of the GNU LGPL v3 or later.
#---------------------------------------------

# Define the hydrobase files to input

##############################
########## PASS 1 ############
##############################
# Raw data to be filtered
#set rawfile = ../step2a_range_check/archive_pass1/all.espna.nf.1950to2015.ge200.hb.rchk1_schk0

# Data to set the envelope of variability
#set smofile = ../step2b_state_bin/archive_pass1/all.espna.nf.1950to2015.ge200.hb.rchk1_schk0.stb

# Define the output filename
#set outfile = all.espna.nf.1950to2015.ge200.hb.rchk1_schk1
#set badfile = all.espna.nf.1950to2015.ge200.hb.rchk1_schk1.bad

##############################
########## PASS 2 ############
##############################
# Raw data to be filtered
set rawfile = ../step2a_range_check/archive_pass2/all.espna.nf.1950to2015.ge200.hb.rchk2_schk1

# Data to set the envelope of variability
set smofile = ../step2b_state_bin/archive_pass2/all.espna.nf.1950to2015.ge200.hb.rchk2_schk1.stb

# Define the output filename
set outfile = all.espna.nf.1950to2015.ge200.hb.rchk2_schk2
set badfile = all.espna.nf.1950to2015.ge200.hb.rchk2_schk2.bad

##############################
##############################
##############################

if ( -e $outfile ) then
    echo ERROR: $outfile already exists.  Exiting.
    exit
endif

if ( -e $badfile ) then
    echo ERROR: $badfile already exists.  Exiting.
    exit
endif

# Define the minimum number of stations
# required to compute mean ts-curves.
@ minsta = 10

# Define the sigma bin input file w/ format:
# sig_bin1_min sig_bin1_max 
#cat << ENDLIST > sigbins.txt
#14.00 14.50 14.50 15.00 15.00 15.50 15.50 16.00 16.00 16.50 16.50 17.00 17.00 17.50 17.50 18.00 18.00 18.50 18.50 19.00 19.00 19.50 19.50 20.00 20.00 20.50 20.50 21.00 21.00 21.50 21.50 22.00 22.00 22.50 22.50 23.00 23.00 23.50 23.50 24.00 24.00 24.50 24.50 25.00 25.00 25.50 25.50 25.60 25.60 25.70 25.70 25.80 25.80 25.90 25.90 26.00 26.00 26.10 26.10 26.20 26.20 26.30 26.30 26.40 26.40 26.50 26.50 26.60 26.60 26.70 26.70 26.80 26.80 26.90 26.90 27.00 27.00 27.05 27.05 27.10 27.10 27.15 27.15 27.20 27.20 27.25 27.25 27.30 27.30 27.35 27.35 27.40 27.40 27.45 27.45 27.50 27.50 27.55 27.55 27.60 27.60 27.65 27.65 27.70 27.70 27.75 27.75 27.80 27.80 27.85 27.85 27.90 27.90 27.95 27.95 28.00 28.00 28.05 28.05 28.10 28.10 28.15 28.15 28.20 28.20 28.25 28.25 28.30 28.30 28.35 28.35 28.40 28.40 28.45 28.45 28.50 28.50 28.55 28.55 28.60 28.60 28.65 28.65 28.70 28.70 28.75 28.75 28.80 28.80 28.85 28.85 28.90 28.90 28.95 28.95 29.00 29.00 29.05 29.05 29.10 29.10 29.15 29.15 29.20 29.20 29.25 29.25 29.30 36.00 36.05 36.05 36.10 36.10 36.15 36.15 36.20 36.20 36.25 36.25 36.30 36.30 36.35 36.35 36.40 36.40 36.45 36.45 36.50 36.50 36.55 36.55 36.60 36.60 36.65 36.65 36.70 36.70 36.75 36.75 36.80 36.80 36.85 36.85 36.90 36.90 36.95 36.95 37.00 37.00 37.05 37.05 37.10 37.10 37.15 37.15 37.20 37.20 37.25 37.25 37.30 37.30 37.35 37.35 37.40 37.40 37.45 37.45 37.50 37.50 37.55 37.55 37.60 37.60 37.65 37.65 37.70 37.70 37.75 37.75 37.80 37.80 37.85 37.85 37.90 37.90 37.95 45.70 45.75 45.75 45.80 45.80 45.85 45.85 45.90 45.90 45.95 45.95 46.00 46.00 46.05 46.05 46.10 46.10 46.15 46.15 46.20 46.20 46.25 46.25 46.30
#ENDLIST

# Create directories to hold ms1 files.
set msdir = ms01_tmp
if ( -e $msdir ) then
    echo ERROR: $msdir already exists.  Exiting.
    exit
endif
mkdir ${msdir}

# Split the input files into 1 degree squares
hb_mssort $smofile -S1 -O${msdir} -N.smo.hb > log_ms5.smo
rm -f ./${msdir}/msextra.dat
hb_mssort $rawfile -S1 -O${msdir} -N.raw.hb > log_ms5.raw
rm -f ./${msdir}/msextra.dat

# Check that the same number of squares are represented.
@ num_ms_smo = `ls -1 ${msdir}/*.smo.hb | wc -l`
@ num_ms_raw = `ls -1 ${msdir}/*.raw.hb | wc -l`
if ( $num_ms_smo != $num_ms_raw ) then
    echo '-----------------------------------------------'
    echo 'WARNING WARNING WARNING WARNING WARNING WARNING'
    echo '-----------------------------------------------'
    echo Number of MS for envelope is different from raw
    echo data.  Specifically:
    echo Env: $num_ms_smo
    echo Raw: $num_ms_raw
    echo '-----------------------------------------------'
    # There is almost always a mismatch and this is due to
    # the fact that the state-binning geographic grids can
    # pull in data from different Marsden Squares and so the
    # data on the edges of the domain do not necessarily have
    # a representation in the other parts.  VERIFY THIS.
endif

# Count the number of stations in each square.
foreach msfile ( ./${msdir}/*.smo.hb )

    # Get the basename of this ms file.
    set msnum = `basename $msfile .smo.hb`

    @ nstations_smo = `hb_countsta $msfile`
    @ nstations_raw = `hb_countsta ./${msdir}/${msnum}.raw.hb`

    echo '-----------------------------------------------'
    echo Square ${msnum} has ${nstations_smo} smo, ${nstations_raw} raw sta.
    echo '-----------------------------------------------'

    # Check that there are enough stations
    # to use this square in the climatology.
    if ( $nstations_smo < $minsta ) then
	echo '-----------------------------------------------'
	echo 'WARNING WARNING WARNING WARNING WARNING WARNING'
	echo '-----------------------------------------------'
    	echo Skipping ${msnum} b/c only ${nstations_smo} stations
	echo '-----------------------------------------------'
    else
	# Compute the statistical envelope around the TS curve.
	hb_statfit_ts $msfile -Ostats_tmp.txt -Ssigbins.txt -Pfit_tmp

	# Filter scans based on statistical envelope.
	hb_statchk_ts -I./${msdir}/${msnum}.raw.hb -Oout_tmp.hb -Sstats_tmp.txt -N2.3 -Q20.0 -Pchk_tmp -Bbad_tmp.hb

	# Use GMT to make diagnostic/documentation plots.
	statplo.sh chk_tmp fit_tmp ts_${msnum}.pdf

	# Add the filtered output to main output file.
	cat out_tmp.hb >> $outfile
	cat bad_tmp.hb >> $badfile

	# Clean up
	rm -f fit_tmp*
	rm -f chk_tmp*
	rm -f out_tmp.hb
	rm -f bad_tmp.hb
	rm -f stats_tmp.txt
	rm -f chk_tmp_ts.dat

    endif  # end of check for sufficient number of stations
end  # of loop over each msfile

# Clean up
rm -rf ${msdir}
