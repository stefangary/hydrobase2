#!/bin/csh -f
#---------------------------------------------
# This script will scan through a series of
# hydrobase stations, infiles, and use the
# hb_rangechk_ts option to remove the fill
# values before other filters are applied.
#
# Operate on $1.
#
# Stefan Gary, 2018
#---------------------------------------------
# This software is distributed under the terms
# of the GNU LGPL v3 or later.
#---------------------------------------------

# Define the hydrobase files to input
set infile = $1
set bn = `basename $infile`

# Define the ending for all the files stored in this run.
set runid = nofill

# Define the output file
set outfile = ${bn}.${runid}
set badfile = ${bn}.${runid}.bad
set logfile = ${bn}.${runid}.log

# No need to run this based on geography
# since the fill values should be uniformly
# -9 (large neg numbers outside normal obs)

# Define the acceptable T-S ranges with depth
# dep_min dep_max t_min t_max s_min s_max
set rangefile = ts_range.txt
cat << ENDLIST > ${rangefile}
   0.0  500.0 -2.0 40.0  0.1 39.0
 500.0 1000.0 -2.0 19.0 34.0 39.0
1000.0 2000.0 -2.0 13.0 34.0 39.0
2000.0 3000.0 -2.0 13.0 34.0 39.0
3000.0 9000.0 -2.0  5.0 34.0 36.0
ENDLIST

# Run the rangecheck
hb_rangechk_ts $infile -B${badfile} -R${rangefile} -O${outfile} -L${logfile} -Q100.0

# Final clean up
rm ${rangefile}
