#!/bin/tcsh -f
# Quick script to count stations in Rockall Trough.
# Stefan Gary 2018
# This software is distributed under the terms of
# the GNU LGPL version 3 or later.

# Set the seasons:
set season_names = ( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec)
set season_nums = ( 1 2 3 4 5 6 7 8 9 10 11 12 )

#set season_names = ( DJF MAM JJA SON )
#set season_nums = ( 1/2/12 3/4/5 6/7/8 9/10/11 )

# Set domain to count stations
set xmin = -15
set xmax = 8
set ymin = 56.5
set ymax = 58.5
set domain = ${xmin}/${xmax}/${ymin}/${ymax}

# HB3
set hb3dir = ~/src/HB3/bin/

#=============================================================
# Operate
#============================================================= 

# For each season,
@ ss = 1
foreach season ( $season_names )

    # Extract the data in time and space
    ${hb3dir}/hb_extract $1 -Tm${season_nums[$ss]} -Tg/${domain} > tmp.hb

    # Count the stations
    set numsta = `hb_countsta tmp.hb`

    # Clean up
    rm -f tmp.hb

    # Print results to stdout
    echo $season $season_nums[$ss] $numsta

    # Move season counter to next season
    @ ss = $ss + 1

end
