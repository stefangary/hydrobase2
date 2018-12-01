#!/bin/tcsh -f
#---------------------------------------------
# This script will scan through a series of
# hydrobase stations, all in one input file
# ($1) and creates a single output file, in
# hydrobase format (${1}.stb) that is state-
# binned.
#
# The fe_ version is the front end for a
# parallelized implementation.  The front
# end will set up the environment for the
# work and will then call the parallel work,
# each parallel work is in the wk_ version. 
#
# Call this wk_ version with >& dump.log
# to redirect both stdout and stderr to
# the dump.log.
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
# and a length scale of ~0.2 degrees, so
# all profiles that fit in each bin in
# time and space will be averaged together
# and then the average values will be
# used as a merged dataset.
#
# NOTE: OUTPUT IS ALWAYS IN T90 FROM THIS
# ROUTINE DUE TO BUGS IN hb_bin3d.
#
# NOTE: MAKE CERTAIN YOU CHECK FOR FILL
# VALUES IN DATA BEFORE RUNNING THIS ROUTINE.
#
# Stefan Gary, August 2018
#---------------------------------------------
# This software is distributed under the terms
# of the GNU LGPL v3 or later.
#---------------------------------------------

# The input file name is $1.
# Useage of parallel will start a new bash shell
# in my home directory, not in the current working
# directory.  So, force it to go to the calling
# directory with $2.

# Set HB3 installation
set hb3dir = /home/sa03sg/src/HB3/bin
set sigdir = /home/sa03sg/src/HB3/lists

# Set binning increment [deg]
set xinc_bin = 0.20
set yinc_bin = 0.20

# Get MS
set ms10rt = `basename $1 .hb | awk -F. '{print $1}'`

# Get month
set month = `basename $1 .hb | awk -F. '{print $3}'`

# Get year
set year = `basename $1 .hb | awk -F. '{print $2}'`

echo Working on MS ${ms10rt} month ${month} year ${year}...

echo Moving to $2
cd $2

#-----------------------------------------------
# STATE BIN- AND ISOPYCNAL- AVERAGE PROFILES
#-----------------------------------------------
# When using HB3, monthly time binning is built
# into the hb_bin3d.  However, monthly binning
# will always attempt to reconstruct the deep
# waters (so anything below the mixed layer and
# the seasonally heated layer) assuming no
# seasonality at depth.  This means that
# the data AT DEPTH are not binned by month,
# rather they are binned by year.  With -U0 and
# -X0, this means that all the monthly files
# are actually identical!  So, do not use
# the built-in monthly binning.
#
# Turn off mixed layer and seasonally heated
# layer special cases with -U0 and -X0
# Specify a 10 m depth spacing for reprojecting
# data onto z-levels with -Zzlev_10m.txt
# Do not specify working on one year's data only with -Y
#   because of year sorting above.
# No distance weighting, -L0
# Create an output file for all possible months, -M1/2/.../11/12
# Specify the root name for the output files, -N
# Force t90 output since te output has lots of spurious
# masking in it.  Convert back to te separately if needed.
# It's OK to input te as it gets automatically converted
# to t90, just don't output it.
${hb3dir}/hb_bin3d ${1} -B`${hb3dir}/hb_msq10bounds ${ms10rt}` -I${xinc_bin}/${yinc_bin} -Pde/pr/t90/sa -Sd${sigdir}/ -Sr${ms10rt} -L0 -Zzlev_10m.txt -M1/2/3/4/5/6/7/8/9/10/11/12 -N${ms10rt}.${year}.${month}.nc.tmp -U0 -X0
#-----------------------------------------------

#-----------------------------------------------
# Dump state-binned, isopycnal-averaged profiles
# back into hb format for use later.  Also,
# update the header on each state-binned file since
# we know the year and month.  Not updating the
# header will cause hb_bin3d to crash if it is
# ever used to operate on state-binned profiles.
${hb3dir}/hb_nc2asc ${ms10rt}.${year}.${month}.nc.tmp | ${hb3dir}/hb_updatehdr -M${month} -Y${year} -F | ${hb3dir}/hb_propcalc -Pde/pr/te/sa > ${1}.stb

# Clean up
rm -f ${ms10rt}.${year}.${month}.nc.tmp
