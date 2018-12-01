#!/bin/tcsh -f
#---------------------------------------------
# This script will scan through a series of
# hydrobase stations, all in one input file
# ($1) and creates a single output file, in
# hydrobase format (${1}.stb) that is state-
# binned.
#
# This fe_ version is the front end for a
# parallelized implementation.  The front
# end will set up the environment for the
# work and will then call the parallel work,
# each parallel work is in the wk_ version. 
#
# State-binning is used to compensate for the
# fact that glider data is much more dense
# than ship based data so we don't want the
# final averages to be biased by multiple
# glider profiles in on area with fewer
# ship profiles.
#
# The general premise is that the ocean
# state is on a time scale of one month
# and a length scale of 0.25 degrees, so
# all profiles that fit in each bin in
# time and space will be averaged together
# and then the average values will be
# used as a merged dataset.
#
# The _monthly variant explicitly has
# monthly subdivions built in to allow
# for more streamlined division and then
# recombination of data.
#
# NOTE: hb_bin3d, which is at the core of this
# routine, does not deal gracefully with
# 1) outputting te - only output t90 in wk_ and
#    all t90 are converted back to te at the
#    end of the fe_ script which recombines all
#    subfiles into a single output file.
# 2) inputting HB_MISSING - check for fill values
#    on data before running.  I could include a
#    fill value filter on here before processing,
#    but I want to keep a close eye on the QC
#    process and for now manually want to check
#    for fill values.
#
# Stefan Gary, 2018
#---------------------------------------------
# This software is distributed under the terms
# of the GNU LGPL v3 or later.
#---------------------------------------------

echo start: `date` > time.log

#=============================================
# Some basic HB setup params:
#=============================================

# Define location of the HB installations
# Using HB2 only for hb_getpos because sfg
# added the -T option to HB2 hb_getpos and
# not present in HB3.
set hb3dir = '/home/sa03sg/src/HB3/bin'
set hb2dir = '/home/sa03sg/src/HB2/bin'

#=============================================
# Define grid parameters for output.
#=============================================

# List the month names in a look up table
set month_nums = ( 01 02 03 04 05 06 07 08 09 10 11 12 )

#=============================================
# Inputs and outputs
#=============================================
# Define the hydrobase station file to input
set infile = $1
set in_bn = `basename $infile`

# Define the output filename
set outfile = ${in_bn}.stb

# Check if this output file exists
if ( -e $outfile ) then
    # Does exist, tell user and quit.
    echo ERROR: $outfile exists, exiting...
    exit
endif
# Initialize a blank output file.
touch $outfile

#=============================================
# Subdomains
#=============================================
# All the spatial binning is split into 10 deg
# Marsden Squares because the density bins for
# each area of the ocean vary slightly from
# location to location and HB provides a very
# good series of suggested density bins for
# each 10 degree square.

# Create a directory to hold subsquares.
set ms10dir = ms10_tmp
if ( -e ${ms10dir} ) then
    echo ERROR: $ms10dir already exists.  Exiting...
    exit
endif
mkdir ${ms10dir} 

# Split the input files into 10 degree squares
${hb3dir}/hb_mssort ${infile} -S10 -O./${ms10dir} -N.hb > log_ms10.txt.tmp

# Split the input 10 degree squares further into months and years
foreach ms10file ( ./${ms10dir}/*.hb )

    # Get the basename of this ms10 file.
    set ms10rt = `basename $ms10file .hb`

    ${hb3dir}/hb_mysort ${ms10rt} -A -D./${ms10dir} -E.hb -N.hb -O./${ms10dir}

    # Clean up - remove this MS10 file so we just have monthly
    # split files in the ms10dir and none of the original spatial
    # splits.
    rm -f $ms10file

end   # End of loop over MS10 files

echo sorted: `date` >> time.log

#=============================================================
# Operate
#=============================================================

# All files in this directory should be sorted by
# MS10, month, year with filenames MS10.yyyy.mm.hb.
#foreach file ( ./${ms10dir}/*.hb )
    # Test that the parallel subset works as expected
    # before parallelizing.
#    wk_state_bin_monthly.sh $file

#end   # End of loop over all files.

# First attempt at parallelization based on
# examples in man parallel.  --jobs option
# limits the number of jobs at any given
# time.  14 on this machine will leave one
# core available for parallelizing and other
# system functions and one core for the user.
set curdir = `pwd`
ls ./${ms10dir}/*.hb | parallel --jobs 14 ${curdir}/wk_state_bin_monthly.sh {} ${curdir}

# Then need to recombine the ${1}.stb files into a single .stb
foreach file ( ./${ms10dir}/*.stb )
    cat $file >> $outfile
end

# Final clean up
rm -rf ${ms10dir}
rm -f log_ms10.txt.tmp

echo done `date` >> time.log

