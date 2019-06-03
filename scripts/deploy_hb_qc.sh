#!/bin/tcsh -f
#==========================
# Script to create and
# populate the HB data QC
# directories.  Scripts
# are *copied* into those
# directories as a way to
# keep any local changes
# and modifications with
# the data that are changed
# as well.
#
# The user must specify where
# the put the HB QC directories.
#
# Stefan Gary, 2018
#===========================
# This software is distributed
# under the terms of the GNU
# LGPL v3 or later.
#===========================

# Specify the location of
# where to deploy the HB QC
set base_hb_qc_dir = $1

echo Setting up HB QC directories in ${base_hb_qc_dir}...

# Get directory of hydrobase installation
# assuming we are launching this script
# from the scripts directory of the
# hydrobase2 folder cloned from github.
cd ..
set hb_dir = `pwd`
cd scripts

# Make directories.
# Note that state binning depends on HB3 (not just HB2)
# but, if skipped, the rest of the pipeline resides in
# just HB2.
mkdir -p ${base_hb_qc_dir}/step1_import
mkdir -p ${base_hb_qc_dir}/step2_range_check
mkdir -p ${base_hb_qc_dir}/step3_state_bin
mkdir -p ${base_hb_qc_dir}/step4_ts_check

# Copy scripts into the directories.
# Any modifications to the local scripts
# will be stored locally in the directories
# with the data they impact.

# Step 1 converts World Ocean Database data
# to hydrobase format.
cp import_wod.sh ${base_hb_qc_dir}/step1_import
ln -s ${hb_dir}/lists/platformlist_wod13.txt ${base_hb_qc_dir}/step1_import/ship_codes.txt

# Step 2 has two scripts because the removal of
# fill values uses hb_rangechk_ts in a very simple
# way and then the automated range check uses
# hb_rangechk_ts in a more complicated way.
cp remove_fill_values.sh ${base_hb_qc_dir}/step2_range_check
cp range_check_basin.sh ${base_hb_qc_dir}/step2_range_check
ln -s ${hb_dir}/lists/zlev_10m.txt ${base_hb_qc_dir}/step2_range_check/zlev_10m.txt

# Step 3 state bins the data using Gary et al., (2018)
# approach to ensuring that data distribution does not
# overly bias the quality control ranges.  There are two
# scripts, one a front end (fe) controlling multiple
# parallel instances of the worker script (wk).
cp fe_state_bin_monthly.sh ${base_hb_qc_dir}/step3_state_bin
cp wk_state_bin_monthly.sh ${base_hb_qc_dir}/step3_state_bin
ln -s ${hb_dir}/lists/zlev_10m.txt ${base_hb_qc_dir}/step3_state_bin/zlev_10m.txt
ln -s ${hb_dir}/lists ${base_hb_qc_dir}/step3_state_bin/lists

# Step 4 checks the data in temperature-salinity
# space.
cp stat_check_basin.sh ${base_hb_qc_dir}/step4_ts_check
ln -s ${hb_dir}/lists/sigbins.txt ${base_hb_qc_dir}/step4_ts_check/sigbins.txt
ln -s ${hb_dir}/scripts/statplo.sh ${base_hb_qc_dir}/step4_ts_check/statplo.sh

# Copy the main QC front end script
cp fe_hb_qc.sh ${base_hb_qc_dir}/

# Done!
echo Hydrobase QC directories and scripts deployed.
echo You will want to copy or link the input data
echo into the directory ${base_hb_qc_dir}/step1_import
echo and then you can run fe_hb_qc.sh.

# And, include these lines if you want to use
# the sample data included with the HB distribution.
gunzip -c ${hb_dir}/examples/ocldb1543698377.11445.OSD.gz > ${base_hb_qc_dir}/step1_import/OSD.wod
gunzip -c ${hb_dir}/examples/ocldb1543698377.11445.PFL.gz > ${base_hb_qc_dir}/step1_import/PFL.wod

