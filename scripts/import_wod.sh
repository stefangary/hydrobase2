#!/bin/tcsh -f
#----------------------------------------
# Script to import the WOD13 download.
# Stefan Gary, 2018
#----------------------------------------
# This software is distributed under the
# terms of the GNU LGPL v3 and later versions.
#----------------------------------------

# Specify HB2 installation
set hb2 = /usr/local/hb2/bin

# Specify the WOD data directory
set dat = ./

# Specify the NODC ship codes file
set ship_codes = ship_codes.txt

# Comment out or use which ever
# relevant instrument types apply.
#${hb2}/wod13_convert ${dat}/*CTD* -Tc -Bbad.ctd.hb -Ogood.ctd.hb -S$ship_codes
${hb2}/wod13_convert ${dat}/*OSD* -Tc -Bbad.osd.hb -Ogood.osd.hb -S$ship_codes
${hb2}/wod13_convert ${dat}/*PFL* -Tc -Bbad.pfl.hb -Ogood.pfl.hb -S$ship_codes

