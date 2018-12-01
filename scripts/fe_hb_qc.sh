#!/bin/tcsh -f
#===========================
# Front end script to controll
# all the HB QC pipeline as
# contained in the scripts here.
#
# Stefan Gary, 2018
#
#===========================
# This software is distributed
# under the terms of the GNU
# LGPL v3 or later.
#===========================

#===========================
# Step 1: Import data to HB
#===========================

# This step may change if
# the input data format
# is different.  In this case,
# we will use World Ocean Database
# (WOD) data as that is
# very common.

echo Importing data...
cd step1_import
ls
cd ..

#===========================
# Step 2: Range check data
#===========================

echo Range checking data...
cd step2_range_check
ls
cd ..

#===========================
# Step 3: State bin data
#===========================

echo State binning data...
cd step3_state_bin
ls
cd ..

#===========================
# Step 4: TS check data
#===========================

echo TS checking data...
cd step4_ts_check
ls
cd ..

#===========================
# Step 5: loop steps 2-4
#===========================

echo Looping over range check, state binning, and ts check...

cd step2_range_check
ls
cd ..

cd step3_state_bin
ls
cd ..

cd step4_ts_check
ls
cd ..

#===========================
# Done!
#===========================
