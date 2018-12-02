#!/bin/csh -f
#---------------------------------------------
# This script will scan through a series of
# hydrobase stations, infiles, and use the
# hb_rangechk_ts option to remove the bad
# scans of T and S based on the range file
# (created herein) as well as an automated
# calculation.
#
# Stefan Gary, 2018
#---------------------------------------------
# This software is distributed under the terms
# of the GNU LGPL v3 or later version.
#---------------------------------------------

# Define location of hydrobase install
set hb2_dir = /usr/local/hb2/bin

# Define the hydrobase files to input
########################
# Pass 1:
#set infile = ../step1_import/all.espna.nf.1950to2015.ge200.hb

########################
# Pass 2:
#ln -sv ../step2c_ts_check/archive_pass1/all.espna.nf.1950to2015.ge200.hb.rchk1_schk1 ./all.espna.nf.1950to2015.ge200.hb
#set infile = ./all.espna.nf.1950to2015.ge200.hb

#######################
# IN GENERAL:
set infile = $1
set bn = `basename $infile`

# Define the ending for all the files stored in this run.
#set runid = rchk2_schk1
set runid = $2

# Define the output file
set outfile = ${bn}.${runid}
set badfile = ${bn}.${runid}.bad
set logfile = ${bn}.${runid}.log

# Create and define a directory to hold intermediate subdivisions.
set ms01dir = ms01_${runid}
if ( -e $ms01dir ) then
    echo ERROR: $ms01dir already exists.  Exiting.
    exit
endif
mkdir $ms01dir

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

# Break up the input files into ms01.
${hb2_dir}/hb_mssort $infile -S1 -O${ms01dir}/ -N.hb > ms01.log.${runid}

# For each ms01 file,
foreach msfile ( ${ms01dir}/*.hb )

    # Get the basename
    set msbn = `basename $msfile .hb`

    # Run the rangecheck
    # Use the 20% criteria as in Lozier et al., 1995
    ${hb2_dir}/hb_rangechk_ts $msfile -B${msbn}.bad.tmp -R${rangefile} -O${msbn}.chk.tmp -L${msbn}.log.tmp -A3.0 -Q20.0

    # Concatenate information
    cat ${msbn}.chk.tmp >> $outfile
    cat ${msbn}.bad.tmp >> $badfile
    cat ${msbn}.log.tmp >> $logfile

    # Incremental clean up
    rm -f ${msbn}.chk.tmp
    rm -f ${msbn}.bad.tmp
    rm -f ${msbn}.log.tmp

end

# Final clean up
rm -rf ${ms01dir}
rm ${rangefile}
