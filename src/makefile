#                 ****************************************
#                               HydroBase 2 
#                 ****************************************
#                    Make file for HydroBase 2 Software
#
#                               Ruth Curry
#                      Woods Hole Oceanographic Institution
#                           Woods Hole, MA  02543
#
#
#-----------------------------------------------------------------------------
#  define these paths before making anything ...
#-----------------------------------------------------------------------------
INCLUDE= .
NETCDFLIB = /usr/local/netcdf3/lib
NETCDFINC = /usr/local/netcdf3/include
#BINDIR = ../bin
BINDIR = /usr/local/hb2/bin
LIBDIR = /usr/local/hb2/lib
#
#-----------------------------------------------------------------------------
DEBUG=  # set to -g to turn on debugging and set OPTIMIZE to a blank
#
OPTIMIZE = +O3  #  optimization flags
#
# these define the compiler commands, flags and system libraries needed ...
#   appropriate options for some systems have been included ...
#   UNcomment the set of options appropriate to your machine
#
#
CC= gcc                       # invoke your C compiler
CFLAGS=-c -I$(INCLUDE)        # compile, include  
LINK= gcc                     # linker 
LDFLAGS =                     # any linker flags
LIBS= -lc -lm                 # math, c libs
#
#-----------------------------------------------------------------------------
#       !! STOP EDITING HERE, THE REST IS FIXED !!
#-----------------------------------------------------------------------------
#
CDF_H = $(NETCDFINC)/netcdf.h
#CDF_L = -L$(NETCDFLIB) -lnetcdf 
CDF_L = $(NETCDFLIB)/libnetcdf.a
FILTERS_H = $(INCLUDE)/hb_filters.h
FIT_H = $(INCLUDE)/hb_fit.h
GAMMA_H = $(INCLUDE)/hb_gamma.h
HYDRO_H = $(INCLUDE)/hydrobase.h
HYDRO_CDF_H = $(INCLUDE)/hydro_cdf.h
MEMORY_H = $(INCLUDE)/hb_memory.h
PATHS_H = $(INCLUDE)/hb_paths.h
SIGBINS_H = $(INCLUDE)/hb_sigbins.h
#
PROG_SRC =  hb_cdf2asc.c hb_cdfinfo.c  hb_cdfmerge.c hb_coriolis.c hb_decimate.c hb_distdepth.c hb_extract.c hb_thicksta.c hb_dupcheck.c hb_findblanks.c hb_fit3d.c hb_fit3d_fg.c hb_fitsurf2d.c hb_foverh.c  hb_getbottom.c hb_getdist.c hb_getpos.c hb_gettrack.c hb_gridsection.c hb_grid2xyz.c  hb_grid3d.c hb_gridsurf2d.c hb_interpolate.c hb_layerav.c  hb_layerinteg.c hb_4matlab.c hb_msq10bounds.c hb_mssort.c hb_propcalc.c hb_section.c hb_siftlevs.c hb_siftsta.c hb_sigma.c hb_sigma_eval.c hb_sigma024_eval.c hb_slice.c hb_smooth2d.c hb_smooth3d.c hb_smoothsection.c  hb_stationsort.c hb_stats1d.c hb_surf.c hb_surfdiff.c hb_theta.c hb_topo2xyz.c hb_toposlice.c hb_tssig.c hb_update_pdr.c hb_read_topo.c hb_updatehdr.c hb_volumetric.c  hb_volumetric_anom.c hb_xyprop.c ices_convert.c sac_convert.c wghc_convert.c  wghc_grid_convert.c woce_btl2hb.c woce_btl2hb_qc.c   woce_btl_qualchk.c  woce_ctd2hb.c sgnc_convert.c wod13_convert.c wod09_convert.c wod05_convert.c wod01_convert.c  wod98_convert.c woa01_convert.c woa98_convert.c woa05_convert.c hb_msnum
PROG_OBJ = ${PROG_SRC:.c=.o}
PROGRAMS = ${PROG_SRC:.c=}
UTILS_SRC = phyprops.c eos80.c gamma_subs.c hydro_utils.c prop_subs.c hydro_cdf.c filters.c fit_subs.c mem_subs.c julconv.c
UTILS_OBJ = ${UTILS_SRC:.c=.o}
OTHER_SRC =  apnc_convert.c awi_convert.c a20_ctdconvert.c bbsr_ctdconvert.c coastal_sta.c ctda_convert.c gage2_btl_convert.c gage2_dcc_convert.c hb_add_ladcp.c hb_saltshift.c hb_linew_avg_section.c hb_monochk.c hb_vrotate.c ifremer_convert.c igor_ctdconvert.c gamma_ex.c kiel_ctdconvert.c mmp_convert.c  npi_ctdconvert.c odf_ctdconvert.c seabird_btl_convert.c seabird_dcc_convert.c primer_ctdconvert.c woce_btl2hb_a20.c yearsmooth.c
OTHER_OBJ = ${OTHER_SRC:.c=.o}
OTHER = ${OTHER_SRC:.c=}
QC_SRC = hb_statchk_ts.c hb_statfit_ts.c hb_rangechk_ts.c hb_statfit_to.c hb_statchk_to.c hb_histsta.c hb_countsta.c hb_sigbasemap.c
QC_OBJ = ${QC_SRC:.c=.o}
QC = ${QC_SRC:.c=}
#-----------------------------------------------------------------------------
#
hydrobase:	$(PROGRAMS)
#
qc: 		$(QC)
#
bathymetry: 	$(BATH)
#
other:          $(OTHER)
#
install:	hydrobase 
	chmod 755 $(PROGRAMS)
	mv $(PROGRAMS) $(BINDIR)
	mkdir -p $(LIBDIR)
	cp ../lib/topo.onetenthdeg.* $(LIBDIR)
	cp ../lib/gamma.nc $(LIBDIR)
	cp ../lib/shipcodes*.asc $(LIBDIR)
#
install_qc:	qc 
	chmod 755 $(QC)
	mv $(QC) $(BINDIR)
#
install_bath:	bathymetry 
	chmod 755 $(BATH)
	mv $(BATH) $(BINDIR)
#
install_other:	other 
	chmod 755 $(OTHER)
	mv $(OTHER) $(BINDIR)
#
clean:
	rm  -f $(PROG_OBJ) $(UTILS_OBJ) $(PROGRAMS) $(QC_OBJ) $(QC) $(OTHER_OBJ) $(OTHER)
clean_all:
	rm  -f $(PROG_OBJ) $(UTILS_OBJ) $(PROGRAMS) $(BATH_OBJ) $(BATH) $(QC_OBJ) $(QC) $(OTHER_OBJ) $(OTHER)
#-----------------------------------------------------------------------------
#
#
#   PROGRAM DEPENDENCIES
#
hb_cdf2asc:	hb_cdf2asc.o hydro_utils.o hydro_cdf.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_cdf2asc.o hydro_utils.o hydro_cdf.o prop_subs.o phyprops.o  eos80.o $(CDF_L) -o hb_cdf2asc $(LIBS)
#
hb_cdfinfo:	hb_cdfinfo.o hydro_cdf.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_cdfinfo.o hydro_cdf.o prop_subs.o phyprops.o  eos80.o $(CDF_L) -o hb_cdfinfo $(LIBS)
#
hb_cdfmerge:	hb_cdfmerge.o hydro_utils.o hydro_cdf.o prop_subs.o phyprops.o eos80.o  mem_subs.o
	$(LINK) $(LDFLAGS) hb_cdfmerge.o hydro_utils.o hydro_cdf.o prop_subs.o phyprops.o  eos80.o  mem_subs.o $(CDF_L) -o hb_cdfmerge $(LIBS)
#
hb_coriolis: 	hb_coriolis.o phyprops.o  eos80.o
	$(LINK) $(LDFLAGS)  hb_coriolis.o phyprops.o eos80.o $(CDF_L) -o hb_coriolis $(LIBS)
#
hb_decimate: 	hb_decimate.o phyprops.o prop_subs.o hydro_utils.o eos80.o
	$(LINK) $(LDFLAGS)  hb_decimate.o phyprops.o prop_subs.o hydro_utils.o eos80.o   $(CDF_L)  -o hb_decimate $(LIBS)
#
hb_distdepth:  	hb_distdepth.o 
	$(LINK) $(LDFLAGS) hb_distdepth.o  -o hb_distdepth $(LIBS)
#
hb_extract:	hb_extract.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_extract.o  hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_extract $(LIBS)
#
hb_thicksta:	hb_thicksta.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_thicksta.o  hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_thicksta $(LIBS)
#
hb_dupcheck:	hb_dupcheck.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_dupcheck.o  hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_dupcheck $(LIBS)
#
hb_findblanks:		hb_findblanks.o 
	$(LINK) $(LDFLAGS) hb_findblanks.o $(LIBS) -o hb_findblanks 
#
hb_fit3d:	hb_fit3d.o fit_subs.o hydro_cdf.o prop_subs.o hydro_utils.o phyprops.o eos80.o mem_subs.o
	$(LINK) $(LDFLAGS) hb_fit3d.o hydro_cdf.o prop_subs.o phyprops.o eos80.o hydro_utils.o fit_subs.o  mem_subs.o $(CDF_L) $(LIBS) -o hb_fit3d 
#
hb_fit3d_fg:	hb_fit3d_fg.o fit_subs.o hydro_cdf.o prop_subs.o hydro_utils.o phyprops.o eos80.o mem_subs.o
	$(LINK) $(LDFLAGS) hb_fit3d_fg.o hydro_cdf.o prop_subs.o phyprops.o eos80.o hydro_utils.o fit_subs.o  mem_subs.o $(CDF_L) $(LIBS) -o hb_fit3d_fg
#
hb_fitsurf2d:	hb_fitsurf2d.o fit_subs.o mem_subs.o
	$(LINK) $(LDFLAGS) hb_fitsurf2d.o fit_subs.o mem_subs.o $(LIBS) -o hb_fitsurf2d 
#
hb_foverh:	hb_foverh.o mem_subs.o phyprops.o eos80.o 
	$(LINK) $(LDFLAGS) hb_foverh.o  phyprops.o eos80.o  mem_subs.o $(LIBS) -o hb_foverh 
#

hb_getbottom: 	hb_getbottom.o phyprops.o prop_subs.o hydro_utils.o gamma_subs.o  eos80.o
	$(LINK) $(LDFLAGS)  hb_getbottom.o phyprops.o prop_subs.o hydro_utils.o eos80.o  gamma_subs.o  $(CDF_L)  -o hb_getbottom $(LIBS)
#
hb_getdist:	hb_getdist.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_getdist.o  hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_getdist $(LIBS)
#
hb_getpos:	hb_getpos.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_getpos.o  hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_getpos $(LIBS)
#
hb_gettrack:  	hb_gettrack.o 
	$(LINK) $(LDFLAGS) hb_gettrack.o  -o hb_gettrack $(LIBS)
#
hb_grid2xyz:	hb_grid2xyz.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_grid2xyz.o  hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_grid2xyz $(LIBS)
#
hb_grid3d:	hb_grid3d.o hydro_cdf.o prop_subs.o hydro_utils.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_grid3d.o hydro_cdf.o prop_subs.o phyprops.o eos80.o hydro_utils.o   $(CDF_L) -o hb_grid3d $(LIBS) 
#
hb_gridsection:	hb_gridsection.o fit_subs.o mem_subs.o
	$(LINK) $(LDFLAGS) hb_gridsection.o fit_subs.o  mem_subs.o $(LIBS) -o hb_gridsection 
#
hb_gridsurf2d: 	hb_gridsurf2d.o phyprops.o prop_subs.o hydro_cdf.o hydro_utils.o eos80.o gamma_subs.o
	$(LINK) $(LDFLAGS)  hb_gridsurf2d.o phyprops.o prop_subs.o hydro_cdf.o hydro_utils.o eos80.o  gamma_subs.o $(CDF_L)  -o hb_gridsurf2d $(LIBS)
#
hb_interpolate: 	hb_interpolate.o phyprops.o prop_subs.o hydro_utils.o eos80.o
	$(LINK) $(LDFLAGS)  hb_interpolate.o phyprops.o prop_subs.o hydro_utils.o eos80.o   $(CDF_L)  -o hb_interpolate $(LIBS)
#
hb_layerav: 	hb_layerav.o phyprops.o prop_subs.o hydro_utils.o eos80.o gamma_subs.o
	$(LINK) $(LDFLAGS)  hb_layerav.o phyprops.o prop_subs.o  hydro_utils.o eos80.o gamma_subs.o $(CDF_L) $(LIBS)  -o hb_layerav
#
hb_layerinteg: 	hb_layerinteg.o phyprops.o prop_subs.o hydro_utils.o eos80.o gamma_subs.o
	$(LINK) $(LDFLAGS)  hb_layerinteg.o phyprops.o prop_subs.o  hydro_utils.o eos80.o gamma_subs.o $(CDF_L) $(LIBS)  -o hb_layerinteg
#
hb_4matlab: 	hb_4matlab.o phyprops.o prop_subs.o hydro_utils.o gamma_subs.o eos80.o
	$(LINK) $(LDFLAGS)  hb_4matlab.o phyprops.o prop_subs.o hydro_utils.o gamma_subs.o eos80.o $(CDF_L) -o hb_4matlab $(LIBS)
#
hb_msq10bounds:	hb_msq10bounds.o
	$(LINK) $(LDFLAGS) hb_msq10bounds.o -o hb_msq10bounds $(LIBS)
#
hb_mssort:	hb_mssort.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_mssort.o  hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_mssort $(LIBS)
#
hb_propcalc: 	hb_propcalc.o phyprops.o prop_subs.o hydro_utils.o eos80.o gamma_subs.o 
	$(LINK) $(LDFLAGS)  hb_propcalc.o phyprops.o prop_subs.o hydro_utils.o eos80.o  gamma_subs.o  $(CDF_L) -o hb_propcalc $(LIBS)
#
hb_section: 	hb_section.o phyprops.o prop_subs.o hydro_utils.o eos80.o gamma_subs.o 
	$(LINK) $(LDFLAGS)  hb_section.o phyprops.o prop_subs.o hydro_utils.o eos80.o  gamma_subs.o $(CDF_L)  -o hb_section $(LIBS)
#
hb_siftlevs:	hb_siftlevs.o hydro_utils.o prop_subs.o phyprops.o eos80.o gamma_subs.o
	$(LINK) $(LDFLAGS) hb_siftlevs.o hydro_utils.o prop_subs.o phyprops.o eos80.o  gamma_subs.o  $(CDF_L) -o hb_siftlevs $(LIBS)
#
hb_siftsta:	hb_siftsta.o hydro_utils.o prop_subs.o phyprops.o eos80.o gamma_subs.o
	$(LINK) $(LDFLAGS) hb_siftsta.o hydro_utils.o prop_subs.o phyprops.o eos80.o  gamma_subs.o  $(CDF_L) -o hb_siftsta $(LIBS)
#
hb_sigma:	hb_sigma.o  phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_sigma.o phyprops.o eos80.o -o hb_sigma $(LIBS)
#
hb_sigma_eval:	hb_sigma_eval.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_sigma_eval.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_sigma_eval $(LIBS)
#
hb_sigma024_eval:	hb_sigma024_eval.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_sigma024_eval.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_sigma024_eval $(LIBS)
#
hb_slice: 	hb_slice.o mem_subs.o phyprops.o prop_subs.o hydro_cdf.o hydro_utils.o gamma_subs.o  eos80.o
	$(LINK) $(LDFLAGS)  hb_slice.o mem_subs.o phyprops.o prop_subs.o hydro_cdf.o hydro_utils.o eos80.o  gamma_subs.o   $(CDF_L)  -o hb_slice $(LIBS)
#
hb_smooth2d:	hb_smooth2d.o filters.o  fit_subs.o mem_subs.o
	$(LINK) $(LDFLAGS) hb_smooth2d.o filters.o fit_subs.o mem_subs.o -o hb_smooth2d $(LIBS)
#
hb_smooth3d:	hb_smooth3d.o filters.o hydro_utils.o hydro_cdf.o prop_subs.o phyprops.o eos80.o mem_subs.o
	$(LINK) $(LDFLAGS) hb_smooth3d.o filters.o  hydro_utils.o hydro_cdf.o  prop_subs.o phyprops.o eos80.o mem_subs.o $(CDF_L)  -o hb_smooth3d $(LIBS)
#
hb_smoothsection:	hb_smoothsection.o filters.o fit_subs.o mem_subs.o
	$(LINK) $(LDFLAGS) hb_smoothsection.o filters.o  fit_subs.o mem_subs.o $(CDF_L)  -o hb_smoothsection $(LIBS)
#
hb_stationsort:	hb_stationsort.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_stationsort.o  hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_stationsort $(LIBS)
#
hb_stats1d:	hb_stats1d.o 
	$(LINK) $(LDFLAGS)  hb_stats1d.o  -o hb_stats1d $(LIBS)
#
hb_surf: 	hb_surf.o phyprops.o prop_subs.o hydro_utils.o eos80.o gamma_subs.o
	$(LINK) $(LDFLAGS)  hb_surf.o phyprops.o prop_subs.o  hydro_utils.o eos80.o gamma_subs.o $(CDF_L) $(LIBS)  -o hb_surf
#
hb_surfdiff:	hb_surfdiff.o 
	$(LINK) $(LDFLAGS) hb_surfdiff.o $(LIBS) -o hb_surfdiff 
#
hb_theta:	hb_theta.o  phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_theta.o phyprops.o eos80.o -o hb_theta $(LIBS)
#
hb_topo2xyz:	hb_topo2xyz.o mem_subs.o
	$(LINK) $(LDFLAGS) hb_topo2xyz.o  mem_subs.o $(LIBS) -o hb_topo2xyz 
#
hb_toposlice:	hb_toposlice.o mem_subs.o phyprops.o prop_subs.o hydro_utils.o eos80.o
	$(LINK) $(LDFLAGS) hb_toposlice.o  mem_subs.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(LIBS) -o hb_toposlice 
#
hb_tssig: 	hb_tssig.o phyprops.o prop_subs.o hydro_utils.o eos80.o
	$(LINK) $(LDFLAGS)  hb_tssig.o phyprops.o prop_subs.o  hydro_utils.o eos80.o  $(CDF_L) $(LIBS)  -o hb_tssig
#
hb_update_pdr:	hb_update_pdr.o phyprops.o prop_subs.o hydro_utils.o eos80.o
	$(LINK) $(LDFLAGS)  hb_update_pdr.o phyprops.o prop_subs.o  hydro_utils.o eos80.o  -o hb_update_pdr $(LIBS)
#
hb_read_topo:	hb_read_topo.o phyprops.o prop_subs.o hydro_utils.o eos80.o
	$(LINK) $(LDFLAGS)  hb_read_topo.o phyprops.o prop_subs.o  hydro_utils.o eos80.o  -o hb_read_topo $(LIBS)
#
hb_updatehdr:	hb_updatehdr.o phyprops.o prop_subs.o hydro_utils.o eos80.o
	$(LINK) $(LDFLAGS)  hb_updatehdr.o phyprops.o prop_subs.o  hydro_utils.o eos80.o  -o hb_updatehdr $(LIBS)
#
hb_volumetric: 	hb_volumetric.o phyprops.o prop_subs.o hydro_utils.o eos80.o gamma_subs.o
	$(LINK) $(LDFLAGS)  hb_volumetric.o phyprops.o prop_subs.o  hydro_utils.o eos80.o gamma_subs.o $(CDF_L) $(LIBS)  -o hb_volumetric
#
hb_volumetric_anom: 	hb_volumetric_anom.o phyprops.o prop_subs.o hydro_utils.o eos80.o gamma_subs.o
	$(LINK) $(LDFLAGS)  hb_volumetric_anom.o phyprops.o prop_subs.o  hydro_utils.o eos80.o gamma_subs.o $(CDF_L) $(LIBS)  -o hb_volumetric_anom
#
hb_xyprop:	hb_xyprop.o phyprops.o prop_subs.o hydro_utils.o eos80.o
	$(LINK) $(LDFLAGS)  hb_xyprop.o phyprops.o prop_subs.o  hydro_utils.o eos80.o -o hb_xyprop $(LIBS)
#
ices_convert:	ices_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o filters.o mem_subs.o
	$(LINK) $(LDFLAGS)  ices_convert.o  phyprops.o prop_subs.o hydro_utils.o eos80.o filters.o mem_subs.o $(CDF_L) -o ices_convert $(LIBS)
#
sac_convert:	sac_convert.o hydro_cdf.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  sac_convert.o  hydro_cdf.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o sac_convert $(LIBS)
#
seabird_btl_convert:	seabird_btl_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  seabird_btl_convert.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o seabird_btl_convert $(LIBS)
#
seabird_dcc_convert:	seabird_dcc_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  seabird_dcc_convert.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o seabird_dcc_convert $(LIBS)
#
wghc_convert:	wghc_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  wghc_convert.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o wghc_convert $(LIBS)
#
wghc_grid_convert:	wghc_grid_convert.o hydro_cdf.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  wghc_grid_convert.o  hydro_cdf.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o wghc_grid_convert $(LIBS)
#
woce_btl2hb:	woce_btl2hb.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  woce_btl2hb.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o woce_btl2hb $(LIBS)
#
woce_btl2hb_qc:	woce_btl2hb_qc.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  woce_btl2hb_qc.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o woce_btl2hb_qc $(LIBS)
#
woce_btl_qualchk:	woce_btl_qualchk.o 
	$(LINK) $(LDFLAGS)  woce_btl_qualchk.o $(CDF_L) -o woce_btl_qualchk $(LIBS)
#
woce_ctd2hb:	woce_ctd2hb.o hydro_utils.o prop_subs.o phyprops.o eos80.o filters.o mem_subs.o
	$(LINK) $(LDFLAGS)  woce_ctd2hb.o phyprops.o prop_subs.o hydro_utils.o eos80.o filters.o mem_subs.o  $(CDF_L) -o woce_ctd2hb $(LIBS)
#
woa01_convert:	woa01_convert.o hydro_cdf.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  woa01_convert.o  hydro_cdf.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o woa01_convert $(LIBS)
#
woa05_convert:	woa05_convert.o hydro_cdf.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  woa05_convert.o  hydro_cdf.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o woa05_convert $(LIBS)
#
woa98_convert:	woa98_convert.o hydro_cdf.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  woa98_convert.o  hydro_cdf.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o woa98_convert $(LIBS)
#
sgnc_convert:	sgnc_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o hydro_cdf.o
	$(LINK) $(LDFLAGS)  sgnc_convert.o phyprops.o prop_subs.o hydro_utils.o eos80.o hydro_cdf.o $(CDF_L) -o sgnc_convert $(LIBS)
#
wod13_convert:	wod13_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  wod13_convert.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o wod13_convert $(LIBS)
#
wod09_convert:	wod09_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  wod09_convert.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o wod09_convert $(LIBS)
#
wod05_convert:	wod05_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  wod05_convert.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o wod05_convert $(LIBS)
#
wod01_convert:	wod01_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  wod01_convert.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o wod01_convert $(LIBS)
#
wod98_convert:	wod98_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  wod98_convert.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o wod98_convert $(LIBS)
#
hb_msnum:	hb_msnum.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  hb_msnum.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o hb_msnum $(LIBS)
#
#
#
#   OBJECT DEPENDENCIES
#
eos80.o:	eos80.c $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) eos80.c
#
filters.o:	filters.c $(FILTERS_H) $(MEMORY_H)
	$(CC)  $(CFLAGS) $(DEBUG) filters.c
#
fit_subs.o:	fit_subs.c $(FIT_H) $(MEMORY_H)
	$(CC)  $(CFLAGS) $(DEBUG) fit_subs.c
#
gamma_subs.o:	gamma_subs.c $(GAMMA_H) $(HYDRO_H) $(CDF_H)
	$(CC) $(CFLAGS)  -I$(NETCDFINC) $(DEBUG) gamma_subs.c
#
hydro_cdf.o:	hydro_cdf.c  $(HYDRO_H) $(HYDRO_CDF_H) $(CDF_H)
	$(CC) $(CFLAGS)  -I$(NETCDFINC) $(DEBUG) hydro_cdf.c
#
hydro_utils.o:	hydro_utils.c $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) hydro_utils.c
#
mem_subs.o:	mem_subs.c $(MEMORY_H)
	$(CC) $(CFLAGS) $(DEBUG) mem_subs.c
#
phyprops.o:	phyprops.c $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) phyprops.c
#
prop_subs.o:	prop_subs.c $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) prop_subs.c
#
# 
#
hb_cdf2asc.o:       hb_cdf2asc.c $(HYDRO_CDF_H) $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_cdf2asc.c
#
hb_cdfinfo.o:       hb_cdfinfo.c $(HYDRO_CDF_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_cdfinfo.c
#
hb_cdfmerge.o:       hb_cdfmerge.c $(HYDRO_CDF_H) $(HYDRO_H) $(MEMORY_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_cdfmerge.c
#
hb_coriolis.o:       hb_coriolis.c $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_coriolis.c
#
hb_decimate.o:	hb_decimate.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_decimate.c
#
hb_distdepth.o:	hb_distdepth.c $(MEMORY_H) $(PATHS_H) 
	$(CC)  $(CFLAGS) $(DEBUG) hb_distdepth.c
#
hb_extract.o:	hb_extract.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_extract.c
#
hb_thicksta.o:	hb_thicksta.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_thicksta.c
#
hb_dupcheck.o:	hb_dupcheck.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_dupcheck.c
#
hb_findblanks.o:	hb_findblanks.c 
	$(CC)  $(CFLAGS) $(DEBUG) hb_findblanks.c
#
hb_fit3d.o:	hb_fit3d.c $(FIT_H) $(HYDRO_H) $(HYDRO_CDF_H) $(MEMORY_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_fit3d.c
#
hb_fit3d_fg.o:	hb_fit3d_fg.c $(FIT_H) $(HYDRO_H) $(HYDRO_CDF_H) $(MEMORY_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_fit3d_fg.c
#
hb_fitsurf2d.o:	hb_fitsurf2d.c $(FIT_H) $(MEMORY_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_fitsurf2d.c
#
hb_foverh.o:	hb_foverh.c  $(MEMORY_H) $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_foverh.c
#
hb_getbottom.o:	hb_getbottom.c $(HYDRO_H) $(PATHS_H) $(GAMMA_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_getbottom.c
#
hb_getdist.o:	hb_getdist.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_getdist.c
#
hb_getpos.o:	hb_getpos.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_getpos.c
#
hb_gettrack.o:	hb_gettrack.c $(FIT_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_gettrack.c
#
hb_grid2xyz.o:	hb_grid2xyz.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_grid2xyz.c
#
hb_grid3d.o:	hb_grid3d.c  $(HYDRO_H) $(HYDRO_CDF_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_grid3d.c
#
hb_gridsection.o:	hb_gridsection.c $(FIT_H) $(MEMORY_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_gridsection.c
#
hb_gridsurf2d.o:	hb_gridsurf2d.c   $(HYDRO_H) $(HYDRO_CDF_H) $(PATHS_H) $(GAMMA_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_gridsurf2d.c
#
hb_interpolate.o:	hb_interpolate.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_interpolate.c
#
hb_layerav.o:	hb_layerav.c   $(HYDRO_H) $(PATHS_H) $(GAMMA_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_layerav.c
#
hb_layerinteg.o:	hb_layerinteg.c   $(HYDRO_H) $(PATHS_H) $(GAMMA_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_layerinteg.c
#
hb_4matlab.o:	hb_4matlab.c $(HYDRO_H) $(PATHS_H) $(GAMMA_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_4matlab.c
#
hb_msq10bounds.o:	hb_msq10bounds.c
	$(CC)  $(CFLAGS) $(DEBUG) hb_msq10bounds.c
#
hb_mssort.o:	hb_mssort.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_mssort.c
#
hb_propcalc.o:	hb_propcalc.c $(HYDRO_H) $(PATHS_H) $(GAMMA_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_propcalc.c
#
hb_section.o:	hb_section.c $(HYDRO_H) $(PATHS_H) $(GAMMA_H)  
	$(CC)  $(CFLAGS) $(DEBUG) hb_section.c
#
hb_siftlevs.o:	hb_siftlevs.c  $(HYDRO_H) $(PATHS_H) $(GAMMA_H) 
	$(CC) $(CFLAGS) $(DEBUG) hb_siftlevs.c
#
hb_siftsta.o:	hb_siftsta.c  $(HYDRO_H) $(PATHS_H) $(GAMMA_H) 
	$(CC) $(CFLAGS) $(DEBUG) hb_siftsta.c
#
hb_sigma.o:	hb_sigma.c  $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_sigma.c
#
hb_sigma_eval.o:	hb_sigma_eval.c  $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_sigma_eval.c
#
hb_sigma024_eval.o:	hb_sigma024_eval.c  $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_sigma024_eval.c
#
hb_slice.o:	hb_slice.c   $(HYDRO_H) $(HYDRO_CDF_H) $(MEMORY_H) $(PATHS_H) $(GAMMA_H)  
	$(CC) $(CFLAGS) $(DEBUG) hb_slice.c
#
hb_smooth2d.o:	hb_smooth2d.c $(MEMORY_H) $(FILTERS_H) $(FIT_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_smooth2d.c
#
hb_smooth3d.o:	hb_smooth3d.c $(HYDRO_H) $(HYDRO_CDF_H) $(FILTERS_H) $(MEMORY_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_smooth3d.c
#
hb_smoothsection.o:	hb_smoothsection.c  $(FILTERS_H) $(MEMORY_H) $(FIT_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_smoothsection.c
#
hb_stationsort.o:	hb_stationsort.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_stationsort.c
#
hb_stats1d.o:	hb_stats1d.c   
	$(CC) $(CFLAGS) $(DEBUG) hb_stats1d.c
#
hb_surf.o:	hb_surf.c   $(HYDRO_H) $(PATHS_H) $(GAMMA_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_surf.c
#
hb_surfdiff.o:	hb_surfdiff.c 
	$(CC)  $(CFLAGS) $(DEBUG) hb_surfdiff.c
#
hb_theta.o:	hb_theta.c  $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_theta.c
#
hb_topo2xyz.o:	hb_topo2xyz.c  $(MEMORY_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_topo2xyz.c
#
hb_toposlice.o:	hb_toposlice.c  $(MEMORY_H) $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_toposlice.c
#
hb_tssig.o:	hb_tssig.c   
	$(CC) $(CFLAGS) $(DEBUG) hb_tssig.c
#
hb_update_pdr.o:	hb_update_pdr.c $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_update_pdr.c
#
hb_read_topo.o:	hb_read_topo.c $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_read_topo.c
#
hb_updatehdr.o:	hb_updatehdr.c $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_updatehdr.c
#
hb_volumetric.o:	hb_volumetric.c   $(HYDRO_H) $(PATHS_H) $(GAMMA_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_volumetric.c
#
hb_volumetric_anom.o:	hb_volumetric_anom.c   $(HYDRO_H) $(PATHS_H) $(GAMMA_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_volumetric_anom.c
#
hb_xyprop.o:	hb_xyprop.c   $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_xyprop.c
#
ices_convert.o:	ices_convert.c $(HYDRO_H) $(FILTERS_H) $(MEMORY_H)
	$(CC)  $(CFLAGS) $(DEBUG) ices_convert.c
#
sac_convert.o:	sac_convert.c $(HYDRO_H) $(HYDRO_CDF_H)
	$(CC)  $(CFLAGS) $(DEBUG) sac_convert.c
#
seabird_btl_convert.o:	seabird_btl_convert.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) seabird_btl_convert.c
#
seabird_dcc_convert.o:	seabird_dcc_convert.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) seabird_dcc_convert.c
#
wghc_convert.o:	wghc_convert.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) wghc_convert.c
#
wghc_grid_convert.o:	wghc_grid_convert.c $(HYDRO_H) $(HYDRO_CDF_H)
	$(CC)  $(CFLAGS) $(DEBUG) wghc_grid_convert.c
#
woce_btl2hb.o:	woce_btl2hb.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) woce_btl2hb.c
#
woce_btl2hb_qc.o:	woce_btl2hb_qc.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) woce_btl2hb_qc.c
#
woce_btl_qualchk.o:	woce_btl_qualchk.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) woce_btl_qualchk.c
#
woce_ctd2hb.o:	woce_ctd2hb.c $(HYDRO_H) $(FILTERS_H)
	$(CC)  $(CFLAGS) $(DEBUG) woce_ctd2hb.c
#
woa01_convert.o:	woa01_convert.c $(HYDRO_H) $(HYDRO_CDF_H)
	$(CC)  $(CFLAGS) $(DEBUG) woa01_convert.c
#
woa05_convert.o:	woa05_convert.c $(HYDRO_H) $(HYDRO_CDF_H)
	$(CC)  $(CFLAGS) $(DEBUG) woa05_convert.c
#
woa98_convert.o:	woa98_convert.c $(HYDRO_H) $(HYDRO_CDF_H)
	$(CC)  $(CFLAGS) $(DEBUG) woa98_convert.c
#
sgnc_convert.o:		sgnc_convert.c $(HYDRO_H) $(CDF_H)
	$(CC)  $(CFLAGS) $(DEBUG) sgnc_convert.c
#
wod13_convert.o:	wod13_convert.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) wod13_convert.c
#
wod09_convert.o:	wod09_convert.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) wod09_convert.c
#
wod05_convert.o:	wod05_convert.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) wod05_convert.c
#
wod01_convert.o:	wod01_convert.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) wod01_convert.c
#
wod98_convert.o:	wod98_convert.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) wod98_convert.c
#
hb_msnum.o:	hb_msnum.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_msnum.c
#
#
#
#    OTHER PROGRAMS
#
apnc_convert:	apnc_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o  julconv.o
	$(LINK) $(LDFLAGS)  apnc_convert.o phyprops.o prop_subs.o hydro_utils.o eos80.o  julconv.o  -o apnc_convert $(CDF_L) $(LIBS)
#
awi_convert:	awi_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o  mem_subs.o
	$(LINK) $(LDFLAGS)  awi_convert.o phyprops.o prop_subs.o hydro_utils.o eos80.o mem_subs.o  -o awi_convert $(LIBS)
#
a20_ctdconvert:     a20_ctdconvert.o hydro_utils.o phyprops.o eos80.o prop_subs.o mem_subs.o
	$(LINK) $(LDFLAGS) a20_ctdconvert.o hydro_utils.o prop_subs.o phyprops.o eos80.o mem_subs.o -o a20_ctdconvert  $(LIBS)
#
bbsr_ctdconvert:     bbsr_ctdconvert.o hydro_utils.o phyprops.o eos80.o filters.o prop_subs.o mem_subs.o
	$(LINK) $(LDFLAGS) bbsr_ctdconvert.o hydro_utils.o prop_subs.o phyprops.o eos80.o filters.o mem_subs.o -o bbsr_ctdconvert  $(LIBS)
#
coastal_sta:	coastal_sta.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) coastal_sta.o  hydro_utils.o prop_subs.o phyprops.o eos80.o -o coastal_sta $(LIBS)
#
ctda_convert:	ctda_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o filters.o mem_subs.o
	$(LINK) $(LDFLAGS)  ctda_convert.o phyprops.o prop_subs.o hydro_utils.o eos80.o filters.o mem_subs.o  -o ctda_convert $(LIBS)
#
gage2_btl_convert:	gage2_btl_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  gage2_btl_convert.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o gage2_btl_convert $(LIBS)
#
gage2_dcc_convert:	gage2_dcc_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  gage2_dcc_convert.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o gage2_dcc_convert $(LIBS)
#
gamma_ex:	gamma_ex.o gamma_subs.o phyprops.o eos80.o 
	$(LINK) $(LDFLAGS) gamma_ex.o gamma_subs.o phyprops.o eos80.o  $(CDF_L) -o gamma_ex $(LIBS)
#
hb_4mike:	hb_4mike.o  hydro_utils.o prop_subs.o phyprops.o eos80.o mem_subs.o
	$(LINK) $(LDFLAGS)  hb_4mike.o phyprops.o prop_subs.o hydro_utils.o eos80.o mem_subs.o  -o hb_4mike $(LIBS)
#
hb_add_ladcp:   hb_add_ladcp.o  hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_add_ladcp.o  hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_add_ladcp  $(LIBS)
#
hb_linew_avg_section:	hb_linew_avg_section.o hydro_utils.o prop_subs.o phyprops.o eos80.o gamma_subs.o
	$(LINK) $(LDFLAGS) hb_linew_avg_section.o  hydro_utils.o prop_subs.o phyprops.o eos80.o gamma_subs.o $(CDF_L)  -o hb_linew_avg_section $(LIBS)
#
hb_monochk: 	hb_monochk.o phyprops.o prop_subs.o hydro_utils.o eos80.o
	$(LINK) $(LDFLAGS)  hb_monochk.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o hb_monochk $(LIBS)
#
hb_saltshift: 	hb_saltshift.o phyprops.o prop_subs.o hydro_utils.o eos80.o
	$(LINK) $(LDFLAGS)  hb_saltshift.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o hb_saltshift $(LIBS)
#
hb_vrotate:	hb_vrotate.o  hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_vrotate.o  hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_vrotate  $(LIBS)
#
ifremer_convert:     ifremer_convert.o hydro_utils.o phyprops.o eos80.o prop_subs.o julconv.o
	$(LINK) $(LDFLAGS) ifremer_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o  julconv.o $(CDF_L) -o ifremer_convert  $(LIBS)
#
igor_ctdconvert:     igor_ctdconvert.o hydro_utils.o phyprops.o eos80.o prop_subs.o mem_subs.o
	$(LINK) $(LDFLAGS) igor_ctdconvert.o hydro_utils.o prop_subs.o phyprops.o eos80.o mem_subs.o -o igor_ctdconvert  $(LIBS)
#
kiel_ctdconvert:	kiel_ctdconvert.o  hydro_utils.o prop_subs.o phyprops.o eos80.o 
	$(LINK) $(LDFLAGS)  kiel_ctdconvert.o phyprops.o prop_subs.o hydro_utils.o eos80.o -o kiel_ctdconvert $(LIBS)
#
mmp_convert:	mmp_convert.o  hydro_utils.o prop_subs.o phyprops.o eos80.o mem_subs.o
	$(LINK) $(LDFLAGS)  mmp_convert.o phyprops.o prop_subs.o hydro_utils.o eos80.o mem_subs.o  -o mmp_convert $(LIBS)
#
npi_ctdconvert:	npi_ctdconvert.o hydro_utils.o prop_subs.o phyprops.o eos80.o filters.o mem_subs.o
	$(LINK) $(LDFLAGS)  npi_ctdconvert.o phyprops.o prop_subs.o hydro_utils.o eos80.o filters.o mem_subs.o  -o npi_ctdconvert $(LIBS)
#
odf_ctdconvert:	odf_ctdconvert.o hydro_utils.o prop_subs.o phyprops.o eos80.o filters.o mem_subs.o
	$(LINK) $(LDFLAGS)  odf_ctdconvert.o phyprops.o prop_subs.o hydro_utils.o eos80.o filters.o mem_subs.o  -o odf_ctdconvert $(LIBS)
#
primer_ctdconvert:	primer_ctdconvert.o hydro_utils.o prop_subs.o phyprops.o eos80.o filters.o mem_subs.o
	$(LINK) $(LDFLAGS)  primer_ctdconvert.o phyprops.o prop_subs.o hydro_utils.o eos80.o filters.o mem_subs.o  -o primer_ctdconvert $(LIBS)
#
woce_btl2hb_a20:	woce_btl2hb_a20.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  woce_btl2hb_a20.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o woce_btl2hb_a20 $(LIBS)
#
yearsmooth:	yearsmooth.o
	$(CC) $(DEBUG) yearsmooth.o -o yearsmooth -lm
#
#
#    OTHER OBJECT DEPENDENCIES
#
apnc_convert.o:	apnc_convert.c $(HYDRO_H) 
	$(CC) $(CFLAGS)  -I$(NETCDFINC)  $(DEBUG) apnc_convert.c
#
julconv.o:      julconv.c
	$(CC) $(CFLAGS) $(DEBUG) julconv.c
#
awi_convert.o:	awi_convert.c $(HYDRO_H) $(MEMORY_H)
	$(CC) $(CFLAGS) $(DEBUG) awi_convert.c
#
a20_ctdconvert.o:	a20_ctdconvert.c $(HYDRO_H) 
	$(CC) $(CFLAGS) $(DEBUG) a20_ctdconvert.c
#
bbsr_ctdconvert.o:	bbsr_ctdconvert.c $(HYDRO_H) $(FILTERS_H)
	$(CC) $(CFLAGS) $(DEBUG) bbsr_ctdconvert.c
#
coastal_sta.o:	coastal_sta.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) coastal_sta.c
#
ctda_convert.o:	ctda_convert.c $(HYDRO_H) $(FILTERS_H) $(MEMORY_H)
	$(CC) $(CFLAGS) $(DEBUG) ctda_convert.c
#
gage2_btl_convert.o:	gage2_btl_convert.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) gage2_btl_convert.c
#
gage2_dcc_convert.o:	gage2_dcc_convert.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) gage2_dcc_convert.c
#
gamma_ex.o:	gamma_ex.c $(GAMMA_H) $(PATHS_H)
	$(CC)  $(CFLAGS) $(DEBUG) gamma_ex.c
#
hb_4mike.o:	hb_4mike.c $(HYDRO_H) $(MEMORY_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_4mike.c
#
hb_add_ladcp.o:  hb_add_ladcp.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_add_ladcp.c
#
hb_linew_avg_section.o:	hb_linew_avg_section.c $(HYDRO_H) $(GAMMA_H)  $(PATHS_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_linew_avg_section.c
#
hb_monochk.o:	hb_monochk.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_monochk.c
#
hb_saltshift.o:	hb_saltshift.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_saltshift.c
#
hb_vrotate.o:	 hb_vrotate.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_vrotate.c
#
ifremer_convert.o:	ifremer_convert.c $(HYDRO_H) $(CDF_H) 
	$(CC) $(CFLAGS) $(DEBUG) -I$(NETCDFINC) ifremer_convert.c
#
igor_ctdconvert.o:	igor_ctdconvert.c $(HYDRO_H) $(FILTERS_H)
	$(CC) $(CFLAGS) $(DEBUG) igor_ctdconvert.c
#
kiel_ctdconvert.o:	kiel_ctdconvert.c $(HYDRO_H) 
	$(CC) $(CFLAGS) $(DEBUG) kiel_ctdconvert.c
#
mmp_convert.o:  mmp_convert.c  $(HYDRO_H)  $(MEMORY_H)
	$(CC) $(CFLAGS) $(DEBUG) mmp_convert.c
#
npi_ctdconvert.o:	npi_ctdconvert.c $(HYDRO_H) $(FILTERS_H) $(MEMORY_H)
	$(CC) $(CFLAGS) $(DEBUG) npi_ctdconvert.c
#
odf_ctdconvert.o:	odf_ctdconvert.c $(HYDRO_H) $(FILTERS_H) $(MEMORY_H)
	$(CC) $(CFLAGS) $(DEBUG) odf_ctdconvert.c
#
primer_ctdconvert.o:	primer_ctdconvert.c $(HYDRO_H) $(FILTERS_H) $(MEMORY_H)
	$(CC) $(CFLAGS) $(DEBUG) primer_ctdconvert.c
#
woce_btl2hb_a20.o:	woce_btl2hb_a20.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) woce_btl2hb_a20.c
#
yearsmooth.o:	yearsmooth.c
	$(CC) $(CFLAGS) $(DEBUG) yearsmooth.c
#
#
#  QC PROGRAMS
#
hb_rangechk_ts:	hb_rangechk_ts.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_rangechk_ts.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_rangechk_ts $(LIBS)
#
hb_histsta:	hb_histsta.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_histsta.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_histsta $(LIBS)
#
hb_countsta:	hb_countsta.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_countsta.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_countsta $(LIBS)
#
hb_sigbasemap:	hb_sigbasemap.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_sigbasemap.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_sigbasemap $(LIBS)
#
hb_statchk_ts:	hb_statchk_ts.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_statchk_ts.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_statchk_ts $(LIBS)
#
hb_statfit_ts:	hb_statfit_ts.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_statfit_ts.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_statfit_ts $(LIBS)
#
hb_rangechk_to:	hb_rangechk_to.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_rangechk_to.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_rangechk_to $(LIBS)
#
hb_statchk_to:	hb_statchk_to.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_statchk_to.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_statchk_to $(LIBS)
#
hb_statfit_to:	hb_statfit_to.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_statfit_to.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_statfit_to $(LIBS)
#
#
#  QC OBJECT DEPENDENCIES
#
hb_rangechk_ts.o:	hb_rangechk_ts.c $(HYDRO_H) 
	$(CC)  $(CFLAGS) $(DEBUG) hb_rangechk_ts.c
#
hb_histsta.o:	hb_histsta.c $(HYDRO_H) 
	$(CC)  $(CFLAGS) $(DEBUG) hb_histsta.c
#
hb_countsta.o:	hb_countsta.c $(HYDRO_H) 
	$(CC)  $(CFLAGS) $(DEBUG) hb_countsta.c
#
hb_sigbasemap.o:	hb_sigbasemap.c $(HYDRO_H) 
	$(CC)  $(CFLAGS) $(DEBUG) hb_sigbasemap.c
#
hb_statchk_ts.o:	hb_statchk_ts.c $(HYDRO_H) $(SIGBINS_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_statchk_ts.c
#
hb_statfit_ts.o:	hb_statfit_ts.c $(HYDRO_H) $(SIGBINS_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_statfit_ts.c
#
hb_rangechk_to.o:	hb_rangechk_to.c $(HYDRO_H) 
	$(CC)  $(CFLAGS) $(DEBUG) hb_rangechk_to.c
#
hb_statchk_to.o:	hb_statchk_to.c $(HYDRO_H) $(SIGBINS_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_statchk_to.c
#
hb_statfit_to.o:	hb_statfit_to.c $(HYDRO_H) $(SIGBINS_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_statfit_to.c
