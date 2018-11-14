/* hb_paths.h
................................................................................
                          *******  HydroBase 2 *******
                             Ruth Curry
                             Woods Hole Oceanographic Institution
                             Mar 2001
...................................................
*/

/*  Replace /usr/local/hydrobase  with the directory path to your local HydroBase2 installation  in each define statement below.  

For BATHPATH, use the file topo.onetenthdeg.swap.dat if you are installing on a PC (Linux) or SGI platform.  For other platforms (HP UNIX, SUN ...) use the file  topo.onetenthdeg.dat    The difference is in the way bytes are stored on different platforms.

*/

#define GAMMA_NC_PATH "/home/sa03sg/src/HB2/lib/gamma.nc"
#define BATHPATH "/home/sa03sg/src/HB2/lib/topo.onetenthdeg.swap.dat"
