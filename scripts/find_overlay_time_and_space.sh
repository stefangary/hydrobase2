#!/bin/tcsh -f
#--------------------------------------
# This script is designed to count the
# number of stations in each month in
# each grid node in a climatology.  All
# stations already have the same coordinates, so
# this is a matter of finding the number
# of stations with the same coordinates.
#
# This _time_and_space variant is designed
# to compute the number of samples in each
# space bin as well as, separately, the
# number of years represented by the samples
# in each bin.
# 
# Stefan Gary, 2018
#
# Modified to just get all stations
# for each grid node (i.e. for the time mean
# maps).
#---------------------------------------
# This software is distributed under the
# terms of the GNU LGPL v3 or later version.
#---------------------------------------

# PART 1 - MONTHLY - disable with if statement.
if ( 1 == 0 ) then

    # Get all station positions and times.
    hb_getpos $1 -T > out.txy.tmp

    @ ii = `wc -l out.txy.tmp | awk '{print $1}'`

    echo Found $ii total stations.

    # Loop over each month
    @ mm = 1
    @ jj = 0
    @ kk = 0
    while ( $mm <= 12 )

	# Extract data for this month only, positions only.
	# Duplicates in space are stored together and unique
	# values are counted.
	awk -v mon=$mm '$2 == mon {print $4,$5}' out.txy.tmp | sort -nr | uniq -c > ${1}.${mm}.cxy

	# Extract data for this month only, positions
	# AND years are stored. Sort stations by the positions
	# and years, duplicates stored together.
	# Count up the number of duplicates, so each line is
	# the number of stations per year per position, x, y, year.
	# Then, strip away the number of stations per year and year
	# and count again to the get the number of unique years
	# sampled at each x,y point.
	#awk -v mon=$mm '$2 == mon {print $4,$5,$1}' out.txy.tmp | sort -nr | uniq | sort -nr | awk '{print $1,$2}' | uniq -c > ${1}.${mm}.cxy.time
	# Duh - since this loop is for each month, and the state
	# binning has already restricted one station per month per year per
	# position, the above line is already computing the number of
	# unique years that are sampling this month at this position.

	@ jj = `gmtmath -C0 ${1}.${mm}.cxy SUM -Sl = | awk '{print $1}'`

	echo Found $jj stations for month $mm

	@ kk = $kk + $jj

	# Move to the next month
	@ mm = $mm + 1

    end

    echo In total, found $kk stations.

    # Clean up
    rm -f out.xyt.tmp
endif

# PART 2 - SEASONAL
# This part is not necessarily the number of unique years since the
# state binning doesn't bin by season, only by month.
if ( 1 == 0 ) then

    # Define seasons and season numbers (subsets of seasons)
    #set season_names = ( DJF MAM JJA SON )
    #set season_nums = ( 1/2/12 3/4/5 6/7/8 9/10/11 )
    #set season_names = ( NDJ FMA MJJ ASO )
    #set season_nums = ( 1/11/12 2/3/4 5/6/7 8/9/10 )

    # All seasons in one go
    #set season_names = ( DJF MAM JJA SON JFM AMJ JAS OND NDJ FMA MJJ ASO )
    #set season_nums = ( 1/2/12 3/4/5 6/7/8 9/10/11 1/2/3 4/5/6 7/8/9 10/11/12 1/11/12 2/3/4 5/6/7 8/9/10 )

    # Two "seasons" warm and cold
    set season_names = ( SONDJF MAMJJA )
    set season_nums = ( 1/2/9/10/11/12 3/4/5/6/7/8 )

    # Loop over each season
    @ ss = 1
    @ jj = 0
    @ kk = 0
    foreach season ( $season_names )

	hb_extract $1 -Tm${season_nums[$ss]} | hb_getpos -T > out.txy.tmp

	# Sort all the positions so that duplicates are
	# together.  Use uniq to count up the duplicates.
	awk '{print $4,$5}' out.txy.tmp | sort -nr | uniq -c > ${1}.${season}.cxy.space

	# Count up all the uniq years sampling this season.
	awk '{print $4,$5,$1}' out.txy.tmp | sort -nr | uniq | sort -nr | awk '{print $1,$2}' | uniq -c > ${1}.${season}.cxy.time


	@ jj = `gmtmath -C0 ${1}.${season}.cxy.space SUM -Sl = | awk '{print $1}'`

	echo Found $jj stations for season $season

	@ kk = $kk + $jj

	# Move to the next season
	rm -f out.txy.tmp
	@ ss = $ss + 1

    end 

endif

# PART 3 - ALL
# Just get all the profiles for each grid node, regardless
# of the time step.
if ( 1 == 1 ) then

    # Get positions
    hb_getpos $1 -T > out.txy.tmp

    # Sort all the positions so that duplicates are
    # together.  Use uniq to count up the duplicates.
    awk '{print $4,$5}' out.txy.tmp | sort -nr | uniq -c > ${1}.ALL.cxy.space

    # Double check that the result is the same as the number of input stations.
    set input = `hb_countsta $1`
    set output = `gmtmath -C0 ${1}.ALL.cxy.space SUM -Sl = | awk '{print $1}'`
    echo Started with $input and ended with $output

endif


