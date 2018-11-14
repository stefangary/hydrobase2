/*  hb_rangechk_ts.c
************************************************************************
*  Usage:  hb_rangechk_ts infile_list  -R<rangefile>
*                                      -O<outfile>
*                                      -B<badfile>
*                                     [-L<logfile>]
*                                     [-A<number_std_deviations>]
*                                     [-D<directory>]
*                                     [-E<extension>]
*                                     [-Q<percent_bad_limit>]
*                                     [-h]
*
*  Reads HydroBase data files and checks that each 
*  temperature and salt value falls within an
*  appropriate range.
*
*  Theta and S  are specified for depth ranges in 
*  the rangefile, specified by the -R option.  The
*  format of all rows in the range file is:
*
*  depth_min depth_max temp_min temp_max salt_min salt_max
*
*  For each bad point, a counter is incremented for
*  that depth range and the scan is eliminated from
*  that profile.
* 
*  The updated station is written to the output file,
*  -O, the bad scans are written to the badfile, -B.
*
*  After the entire file has been searched, the
*  numbers of out-of-range values are written to the
*  logfile, -L.
*
*  An additional option is to run in automatic mode
*  (-A) where all the input stations are binned by
*  depth (bin edges are defined by standard levels)
*  and any scans (data points) in depth bins are
*  eliminated if they lie a certain number of standard
*  deviations from the mean (specified by the number
*  after the -A flag).  This range check is done
*  AFTER the user-specified range check so that any
*  really bad values do not throw off the depth bin
*  means.
*
*  In order to add more control to the quality
*  check process, we also define the quota option
*  (-Q) where we specify a percentage of scans.  If
*  the number of out of range scans in a station
*  exceeds this quota, then we reject the whole
*  station.
*
**************************************************************************
*/ 
#include <stdio.h>
#include <stdlib.h>   /*sfg added to deal with exit and malloc warnings*/
#include <math.h>
#include "hydrobase.h"

#define  PRINT_MSG   1 
#define  MAXRANGES  200 
#define    DIR     ""
#define    EXTENT   ""

/* Define arrays to store min/max depth, temperature,
 * and salinity at each depth bin - NORMAL MODE. */
double dmin[MAXRANGES], dmax[MAXRANGES];
double tmin[MAXRANGES], tmax[MAXRANGES];
double smin[MAXRANGES], smax[MAXRANGES];

/* Define arrays to store min/max depth, mean,
 * and standard deviation at each depth bin.
 * zbin = depth bins
 * _avg = bin averages
 * _std = bin standard deviations
 * _npt = bin number of points */
double zbin[85] = {0,10,20,30,50,75,100,125,150,200, 
		 250,300,350,400,450,500,550,600,
		 650,700,750,800,850,900,950,
		 1000,1100,1200,1300,1400,1500,
		 1600,1700,1800,1900, 
		 2000,2100,2200,2300,2400,2500,
		 2600,2700,2800,2900,
		 3000,3100,3200,3300,3400,3500,
		 3600,3700,3800,3900,
		 4000,4100,4200,4300,4400,4500,
		 4600,4700,4800,4900,
		 5000,5100,5200,5300,5400,5500,
		 5600,5700,5800,5900,
		 6000,6200,6400,6600,6800,
		 7000,7500,8000,9000,99999};
double tavg[84], tstd[84];
double savg[84], sstd[84];
int tnpt[84], snpt[84];
int nbin = 84;

/* Prototypes for locally defined functions */
int  load_ranges(FILE *);      /* Load rangefile */
int  find_index(double, int);  /* Find index in depth arrays dmin/dmax */
void print_usage(char *, double);      /* Print usage if user error */

main(int argc, char **argv) {

  /* Local counter variables:
   * ngood    = number of good scans per station
   * nbad     = number of bad scans per station
   * nranges  = number of depth bins in rangefile
   * stacount = number of stations read (over all files)
   * outcount = number of good stations
   * badcount = number of bad stations
   * indx     = temp holder of index in user-spec'd dep bins
   * binx     = temp holder of index in auto-mode dep bins
   */
   int    ngood, nbad, nranges, stacount, outcount, badcount;
   int    status, indx, binx, badflag;
   int    n, i, j;

   /* nfiles = number of files on command line
    * curfile = current file number */
   int    nfiles, curfile;

   /* File IDs of hydrobase format files */
   int    infile, outfile, badfile;

   /* Arrays to store the number of bad
    * scans in both normal and automatic modes. */
   int    *tbad, *sbad;
   int    *tbada, *sbada;

   double s, t, *theta;
   FILE   *logfile, *rangefile;
   char   *dir, *extent;
   
   /* Input, good, and bad data storage */
   struct HYDRO_HDR hdr;   
   struct HYDRO_DATA data, newdata, baddata;
   
   /* Command line flag for automatic mode option. */
   int a_flag;
   double nstd, delta;
   char *tmpstr;

   /* Command line flag for quota check option. */
   int q_flag, qtrump;
   double quota;
   
   int v_flag;

   fprintf(stderr,"\n Starting hb_rangechk_ts.");

   /* Are there command line arguments? */
   if (argc < 2) {
     print_usage(argv[0],quota);
      exit(1);
   }

   /* Set these default values... */
   dir = DIR;
   extent = EXTENT;
   logfile = (FILE *)NULL;
   rangefile = (FILE *)NULL;
   infile = outfile = badfile = 0;
   stacount = outcount = badcount = 0; /* Initialize station counts */
   nfiles = 0;
   curfile = 1;

   a_flag = 0;
   nstd = 4.0;

   q_flag = 0;
   qtrump = 0;
   quota = 100.0;
   
   v_flag = 0;

   /* Initialize input, good, and bad data storage */
   for (i = 0; i < MAXPROP; ++i) {
     data.observ[i] = (double *) NULL;
     newdata.observ[i] = (double *) NULL;
     baddata.observ[i] = (double *) NULL;
   }
   hdr.prop_id = (int *)NULL;
   theta = (double *)NULL;
 

   /* Parse command line arguments... */
   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {

               case 'D':
		 /* Set the input file directory. */
		 dir = &argv[i][2];
		 break;

               case 'E':
		 /* Set the input file extension. */
		 extent = &argv[i][2];
		 break;

               case 'B':
		 /* Create the bad data output file = badfile */
		 badfile = create_hydro_file(&argv[i][2],OVERWRITE);
		 if (badfile < 0) {
		   fprintf(stderr,"\nUnable to open %s for output. \n", &argv[i][2]);
		   exit(1);
		 }
		  
		 fprintf(stderr,"\nOpened %s for output.\n",&argv[i][2]); 		  
		 break;

               case 'O':
		 /* Create the good data output file = outfile */
		 outfile = create_hydro_file(&argv[i][2],OVERWRITE);
		 if (outfile < 0) {
		   fprintf(stderr,"\nUnable to open %s for output. \n", &argv[i][2]);
		   exit(1);
		 }

		 fprintf(stderr,"\nOpened %s for output.\n",&argv[i][2]);
		 break;

               case 'L':
		 /* Create the logfile. */
		 logfile = fopen(&argv[i][2],"a");
		 if (logfile == NULL) {
		   fprintf(stderr,"\nUnable to open %s for writing. \n", &argv[i][2]);
		   exit(1);
		 }
		 
		 fprintf(stderr,"\nOpened %s in append mode.\n",&argv[i][2]);
		 break;

               case 'R':
		 /* Open the rangefile. */
		 rangefile = fopen(&argv[i][2],"r");
		 if (rangefile == NULL) {
		   fprintf(stderr,"\nUnable to open %s for reading. \n", &argv[i][2]);
		   exit(1);
		 }
		 fprintf(stderr,"\nOpened %s \n",&argv[i][2]); 		         
		 break;

	       case 'A':
		 /* Automatic mode specified.
		  * Read in the number of std dev. */
		 a_flag = 1;
		 tmpstr = &argv[i][2];
		 if( sscanf(tmpstr,"%lf", &nstd) != 1 ) {
		   fprintf(stderr,"\nError reading -A option %lf\n",nstd);
		   exit(1);
		 }
		 break;

	       case 'Q':
		 /* Quota check is specified. */
		 q_flag = 1;
		 tmpstr = &argv[i][2];
		 if ( sscanf(tmpstr,"%lf", &quota) != 1 ) {
		   fprintf(stderr,"\nError reading -A option %lf\n",nstd);
		   exit(1);
		 }
		 if ( quota < 0.0 || quota > 100.0 ) {
		   fprintf(stderr,"\nQuota -Q%lf is not within 0 to 100!",quota);
		 }
		 break;

	       case 'h':
		 print_usage(argv[0],quota);
		  exit(0);

	       case 'V':
		 v_flag = 1;
		 break;

               default :
                   fprintf(stderr,"\nError parsing command line");
                   fprintf(stderr,"\n in particular: %s\n", argv[i]);
                   fprintf(stderr,"Use -h for help\n");
                   exit(1);
            }  /* end switch */
       }  /* end if */
       else {
	 /* If none of the switches apply, then the
	  * command line argument must be one of the
	  * list of infiles.  Augment the input file
	  * counter. */
          ++nfiles;
       }
   }  /* end for */

   /* Check that the user has given enough information on command line. */
   if ( rangefile == NULL  ) {
     print_usage(argv[0],quota);
       fprintf(stderr,"\n\nYou must specify a range file with -R.\n");
       exit(1);
   }

   if ( !outfile  ) {
     print_usage(argv[0],quota);
     fprintf(stderr,"\n You must specify an output file for good scans with -O.\n");
     exit(1);
   }

   if ( !badfile  ) {
     print_usage(argv[0],quota);
     fprintf(stderr,"\n You must specify an output file for bad scans with -B.\n");
     exit(1);
   }
   
   if ( v_flag == 1 ) {
     fprintf(stderr,"\nReading rangefile...");
   }
   /* Read rangefile and allocate memory for keeping
    * track of numbers of bad and good values in
    * each depth range.*/
   if ((nranges = load_ranges(rangefile)) < 0) {
     /* Returned value, nranges, is the number of
      * depth bins/lines in the range file. */
     fprintf(stderr,"\n No range data loaded from rangefile, -R.\n");
     exit(1);
   }
   tbad = (int *) calloc(nranges,  sizeof(int));
   sbad = (int *) calloc(nranges,  sizeof(int));
   tbada = (int *) calloc(nbin,  sizeof(int));
   sbada = (int *) calloc(nbin,  sizeof(int));
   if ( v_flag == 1 ) {
     fprintf(stderr,"\nDone reading rangefile.");
   }

   /* Initialize the profile statistics. */
   for (i=0; i < nbin; i++) {
     tavg[i] = 0.0;
     savg[i] = 0.0;
     tstd[i] = 0.0;
     sstd[i] = 0.0;
     tnpt[i] = 0;
     snpt[i] = 0;
     tbada[i] = 0;
     sbada[i] = 0;
   }

   /*==================FIRST PASS ON DATA==================
    * If the user requests an automatic mode check,
    * then we need to first open all the files/stations
    * and compute the profile mean and standard deviation.
    * The next pass is the one where we actually edit the
    * data and write to file. Note that we must
    * do two passes through the data because we only
    * have the final value for the mean/std once we've
    * cycled through all the data at least one complete time. */
   if( a_flag ) {

     /* Loop for each input file. */
     do {
       if ( !nfiles) {
	 /* If there are no infiles, look to stdin. */
	 infile = STDIN;
	 fprintf(stderr,"\n Expecting data from stdin....  ");
       }
       else {
	 /* Open the current infile. */
	 infile = open_hydro_file(dir, argv[curfile], extent, PRINT_MSG);
	 if (infile  < 0)
	   /* We do not have a valid infile.x */
	   goto NEXTFILE;
       }
       /* Loop for each station */
       while ((status = get_station(infile, &hdr, &data)) == 0) {
	 /* get_station returns 0 for sucessful read. */
	 /* Input data are in structures hdr and data. */

	 /* Augment the number of stations opened. */
	 ++stacount;

	 /* Initialize good and bad scan counters. */
	 ngood = nbad = 0;

	 /* Compute potential temperature from given T,S,P.
	  * theta is a local scope variable only. */
	 free_and_alloc(&theta, hdr.nobs);
	 compute_theta(hdr.nobs, theta, data.observ[(int)PR], data.observ[(int)TE], data.observ[(int)SA]);

	 /* Loop for each data scan */
	 for (i = 0; i < hdr.nobs; ++i) {
	   /* For each depth in the profile, */

	   /* Find the depth bin (from the rangefile) that bounds this depth */
	   indx = find_index(data.observ[(int)DE][i], nranges);
	   badflag = 0; /* Initialize the badflag = "so far, so good" */

	   /* Extract the pot. temp. value at this level */
	   t = theta[i];

	   /* Check that data are within pot. temp. range at this depth. */
	   if ((t < tmin[indx]) || (t > tmax[indx])) {
	     badflag = 2;  /* Flag the scan. */
	     tbad[indx]++;
	   }

	   /* Extract the salinity value at this level */
	   s = data.observ[(int)SA][i];

	   /* Check that the salinity at this level is in good range. */
	   if ((s < smin[indx]) || (s > smax[indx])) {
	     ++badflag;  /* Flag the scan */
	     sbad[indx]++;
	   }

	   /* Buffer this scan if a flag has been thrown. */
	   if (badflag) {
	     /* Do nothing here since this pass is
	      * for computing profile statistics. */
	   }
	   else {
	     /* The data are good. */
	     
	     /* Include these good data into the
	      * running means and std dev. */

	     /* (1)---Find the bin to which the scan belongs.--- */
	     binx = find_bin(data.observ[(int)DE][i],nbin);

	     /* (2)---Compute statistics using online algorithm.--- */
	     
	     /* Running mean and (unnormalized) variance. */
	     delta = t - tavg[binx];
	     tavg[binx] = (tnpt[binx]*tavg[binx] + t)/(tnpt[binx] + 1);
	     tnpt[binx]++;
	     tstd[binx] = tstd[binx] + delta*(t - tavg[binx]);
	     
	     delta = s - savg[binx];
	     savg[binx] = (snpt[binx]*savg[binx] + s)/(snpt[binx] + 1);
	     snpt[binx]++;
	     sstd[binx] = sstd[binx] + delta*(s - savg[binx]);

	   } /* End of data flag check */
	 } /* End for loop over each scan */
       }  /* End of while loop over all stations */
       
       report_status(status, stderr);
       if (nfiles) close(infile);
       
     NEXTFILE:
       ;
     } while (curfile++ < nfiles); /* End of loop over all input files */

     /* Finalize computation of std dev.
      * from un-normalized variances. Also,
      * Multiply in the number of standard
      * deviations at this point to reduce
      * the need for computing this many,
      * repeated time later. */
     for (i=0; i<nbin; i++) {
       if( tnpt[i] > 1 ) {
	 tstd[i] = nstd * sqrt(fabs(tstd[i]/(tnpt[i]-1)));
       }
       else {
	 tstd[i] = 0.0;
       }

       if( snpt[i] > 1 ) {
	 sstd[i] = nstd * sqrt(fabs(sstd[i]/(snpt[i]-1)));
       }
       else {
	 sstd[i] = 0.0;
       }
     } /* End of for loop over all auto-mode bins. */
   } /* End automatic mode check. */

   /*===================SECOND PASS================
    * Once we know the mean and standard deviation
    * of the profile, we can proceed with sorting
    * the good and bad data points. */
   
   /* Loop for each input file. */
   curfile = 1;
   stacount = 0;
   do {
      if ( !nfiles) {
	/* If there are no infiles, look to stdin. */
	infile = STDIN;
	fprintf(stderr,"\n Expecting data from stdin....  ");
      }
      else {
	/* Open the current infile. */
	infile = open_hydro_file(dir, argv[curfile], extent, PRINT_MSG);
	if (infile  < 0)
	  /* We do not have a valid infile.x */
          goto NEXTFILE2;
      }
      /* Loop for each station */
      while ((status = get_station(infile, &hdr, &data)) == 0) {
	/* get_station returns 0 for sucessful read. */
	/* Input data are in structures hdr and data. */

	/* Augment the number of stations opened. */
	++stacount;

	/* Initialize good and bad scan counters. */
	ngood = nbad = 0;

	/* Create space to store each variable in both good and bad data struct */
	for (i = 0; i < hdr.nprops; ++i) {
	  newdata.observ[hdr.prop_id[i]] = (double *) calloc((size_t)hdr.nobs, sizeof(double));
	  baddata.observ[hdr.prop_id[i]] = (double *) calloc((size_t)hdr.nobs, sizeof(double));
	}

	/* Compute potential temperature from given T,S,P.
	 * theta is a local scope variable only. */
	free_and_alloc(&theta, hdr.nobs);
	/* WARNING! NO CHECK THAT PR IS AVAILABLE OR NOT!!!  Will simply segfault with no
	 * additional information if PR is not available. */
	compute_theta(hdr.nobs, theta, data.observ[(int)PR], data.observ[(int)TE], data.observ[(int)SA]);
	/* Loop for each data scan */
	for (i = 0; i < hdr.nobs; ++i) {
	  /* For each depth in the profile, */

	  /* Find the depth bins (from rangefile and auto-mode
	   * bins) that bound this depth */
	  indx = find_index(data.observ[(int)DE][i], nranges);
	  binx = find_bin(data.observ[(int)DE][i],nbin);
	  badflag = 0; /* Initialize the badflag = "so far, so good" */

	  /* Extract the pot. temp. value at this level */
	  t = theta[i];

	  /* Check that data are within usr-specified
	   * pot. temp. range at this depth. */
	  if ( (t < tmin[indx]) || (t > tmax[indx]) ) {
	    badflag = 2;  /* Flag the scan. */
	  }

	  /* Check that data are within auto-mode
	   * pot. temp. range at this depth. */
	  if ( a_flag ) {
	    if (( t<(tavg[binx]-tstd[binx]) )||( t>(tavg[binx]+tstd[binx]) )) {
	      badflag = 2;  /* Flag the scan. */
	    }
	  }

	  /* Extract the salinity value at this level */
	  s = data.observ[(int)SA][i];
	  
	  /* Check that the salinity at this level is in good range. */
	  if ( ( s < smin[indx] ) || ( s > smax[indx] ) ) {
	    ++badflag;  /* Flag the scan */
	  }

	  if ( a_flag ) {
	    if (( s<(savg[binx]-sstd[binx]) )||( s>(savg[binx]+sstd[binx]) )) {
	      if ( badflag == 1 || badflag == 3 ) {
		/* No need to flag the scan b/c it
		* has been salt-only flagged in previous check.*/
	      }
	      else {
		++badflag;  /* Flag the scan */
	      }
	    }
	  }

	  /* Buffer this scan if a flag has been thrown. */
	  if (badflag) {

	    /* Copy all properties at this level to the bad data. */
	    for (j = 0; j < hdr.nprops; ++j) {
	      baddata.observ[hdr.prop_id[j]][nbad] = data.observ[hdr.prop_id[j]][i];
	    }

	    /* Increment the number of bad scans */
	    ++nbad;

	    if (badflag == 1) {
	      /* Salinity only was bad, augment bad salt counter at this level. */
	      /*++sbad[indx];*/
	      ++sbada[binx];
	    }
	    if (badflag >= 2) {
	      /* Temperature was bad, augment bad temp counter at this level. */
	      /*++tbad[indx];*/
	      ++tbada[binx];
	    }
	    if (badflag > 2) {
	      /* Saltinity was also bad (along with temp), augment counter. */
	      /*++sbad[indx];*/
	      ++sbada[binx];
	    }
	  }
	  else {
	    /* The data are good. */

	    /* Copy scan to good data struct. */
	    for (j = 0; j < hdr.nprops; ++j) {
	      newdata.observ[hdr.prop_id[j]][ngood] = data.observ[hdr.prop_id[j]][i];
	    }

	    /* Augment good scans counter */
	    ++ngood;

	  } /* End of data flag check */
	} /* End for loop over each scan */

	/* Set quality control bytes. */
	hdr.qual[3] = '1';

	/* Perform the quota check */
	qtrump = 0;
	if ( q_flag ) {
	  if ( 100.0*((double)nbad)/((double)(hdr.nobs)) > quota ) {

	    /* We have exceeded the quota for this station.
	     * Setting ngood to zero will prevent this
	     * station from being written to the outfile. */
	    ngood = 0;

	    /* Then, we raise a flag so that we know to write
	     * the whole station to the bad file (not just its
	     * bad scans. */
	    qtrump = 1;
	  }
	} /* End quota check. */

	/* If there are any good scans, */
	if (ngood) {

	  /* Set the number of good scans in good data.
           * Note the similar header as the original
           * input data. */
	  hdr.nobs = newdata.nobs = ngood;

	  /* Keep the number of props the same. */
	  newdata.nprops = data.nprops;

	  /* Write only the good scans to the outfile. */
	  write_hydro_station(outfile, &hdr, &newdata);

	  /* Increment the number of profiles written to outfile. */
	  ++outcount;
	}
	
	/* If there are any bad scans, */
	if (nbad) {

	  if (qtrump == 0) {
	    /* We write only the bad scans to the badfile. */

	    /* Set the number of bad scans in bad data.
	     * Note the similar header as the original
	     * input data. */
	    hdr.nobs = baddata.nobs = nbad;

	    /* Keep the number of props the same. */
	    baddata.nprops = data.nprops;
	    
	    /* Write only the bad scans to the badfile. */
	    write_hydro_station(badfile, &hdr, &baddata);
	    
	  }
	  else {
	    /* We write the whole station to the badfile. */
	    write_hydro_station(badfile, &hdr, &data);
	  }

	  /* Increment the number of profiles written to the badfile. */
	  ++badcount;
	  
	  /* Reset the quota trumping flag. */
	  qtrump = 0;
	}
	
	/* Free memory for data storage. */
	for (i = 0; i < hdr.nprops; ++i) {
          free(newdata.observ[hdr.prop_id[i]]);
          free(baddata.observ[hdr.prop_id[i]]);
	}
      }  /* End of while loop over all stations. */

      report_status(status, stderr);
      if (nfiles) close(infile);
      
   NEXTFILE2:
      ;
   } while (curfile++ < nfiles); /* End of loop over all input files. */

   /* Write to logfile if user desires log info. */   
   if (logfile) {
     fprintf(logfile,"\n\nstations in: %7d  retained stations : %7d  edited stations : %7d\n", stacount, outcount, badcount);
     fprintf(logfile,"\nDepth range    #Tbad    #Sbad \n");
     for (i = 0; i < nranges; ++i) {
       fprintf(logfile,"%7.1lf %7.1lf %7d %7d\n", dmin[i], dmax[i], tbad[i], sbad[i]);
     }
     if ( a_flag ) {
       fprintf(logfile,"\n\nDepth (auto)   #Tbad    #Sbad    tavg     tstd     savg     sstd      tnpt     snpt     \n");
       for (i = 0; i < nbin; i++) {
	 fprintf(logfile,"%7.1lf %7.1lf %7d %7d %7.3f %7.3f %7.3f %7.3f %7d %7d \n",
		 zbin[i], zbin[i+1], tbada[i], sbada[i], 
		 tavg[i], tstd[i], savg[i], sstd[i],
		 tnpt[i],snpt[i]);
       }
     }
   } /* End of logfile check. */
   
   fprintf(stderr,"\nEnd of %s\n", argv[0]);
   exit(0);

}  /* End main. */

/****************************************************************************/

void print_usage(char *program, double q) {
  fprintf(stderr,"\nUsage:  %s list_of_infile(s) -B<badfile>",program);
  fprintf(stderr,"\n                             -R<rangefile>");
  fprintf(stderr,"\n                             -O<outfile>");
  fprintf(stderr,"\n                            [-L<logfile>]");
  fprintf(stderr,"\n                            [-A<num_std_dev>]");
  fprintf(stderr,"\n                            [-D<directory>]");
  fprintf(stderr,"\n                            [-E<extension>]");
  fprintf(stderr,"\n                            [-Q<percent>]");
  fprintf(stderr,"\n                            [-h]");
  fprintf(stderr," ");
  fprintf(stderr,"\n List of infiles should be first argument(s) ");
  fprintf(stderr,"\n   -R  : specifies file containing min/max values for");
  fprintf(stderr,"\n         acceptable depth, theta, salt");
  fprintf(stderr,"\n        ex: -R7102.ranges ");
  fprintf(stderr,"\n   -B  : specifies output file for bad scans.");  
  fprintf(stderr,"\n        ex: -I7102.rbad ");
  fprintf(stderr,"\n   -O  : specifies output file for good scans.  ");
  fprintf(stderr,"\n        ex: -O7102.rchk ");
  fprintf(stderr,"\n   -D  : specifies dir for input files (default is ./) ");
  fprintf(stderr,"\n        ex: -D../data/ ");
  fprintf(stderr,"\n   -E  : specifies file extent for input files (default is no extent)");  
  fprintf(stderr,"\n        ex: -E.dat ");
  fprintf(stderr,"\n  [-L] : specifies file to which log info is written ");
  fprintf(stderr,"\n        ex: -L7102.rchklog ");
  fprintf(stderr,"\n  [-A] : specifies automatic mode where the mean and");
  fprintf(stderr,"\n         standard deviation of the salt and temp");
  fprintf(stderr,"\n         profiles (with depth) are computed and all");
  fprintf(stderr,"\n         values outside the specified number of ");
  fprintf(stderr,"\n         standard deviations are rejected as bad scans.");
  fprintf(stderr,"\n  [-Q] : specifies the max percentage of bad scans");
  fprintf(stderr,"\n         for a station to be accepted as good.  If");
  fprintf(stderr,"\n         the percentage of bad scans exceeds this value");
  fprintf(stderr,"\n         the whole station is rejected. [%lf]",q);
  fprintf(stderr,"\n        ex: -Q100.0 will allow any number of scans to be dropped.");
  fprintf(stderr,"\n            -Q0.0   will reject the station even if one scan is bad.");
  fprintf(stderr,"\n  [-h] : help -- prints this message. ");
  fprintf(stderr,"\n  [-V] : Verbose mode.");
  fprintf(stderr,"\n\n WARNING! NO ERROR CHECKING IF PR, TE, OR SA IS MISSING");  
  fprintf(stderr,"\n THIS PROGRAM WILL SIMPLY SEGFAULT.");  
  fprintf(stderr,"\n\n");  
  return;
} /* end print_usage() */
/****************************************************************************/


int load_ranges(FILE *rangefile)

/* Reads the rangefile and loads the range values into the global 
 * variables.  Returns the number of range values read or -1 if an
 * error occurs.  In case of an error, an appropriate message is 
 * written to the stderr device.  The format of the range file is:
 *
 * dens_min dens_max temp_min temp_max salt_min salt_max
 *
 **/
{
   int n, i;   /* n = file column counter, i = file line counter */
   double r;
   char buffer[200];

   i = 0;
   
   /* Read in each line, up to the newline termination */
   while (fscanf(rangefile,"%[^\n]", buffer) != EOF) {
   
     getc(rangefile);  /* move past LF */
   
     /* Scan each line. sscanf returns n, the number of columns. */
     if ((n = sscanf(buffer,"%lf %lf %lf %lf %lf %lf", &dmin[i], &dmax[i], &tmin[i], &tmax[i], &smin[i], &smax[i])) != 6) 
       fprintf(stderr, "\nWARNING: Error reading rangefile : only %1d columns were read at line %3d", n, i);
     else
       /* Check that we have not exceeded acceptable number of
        * density bins.  Also, note that we augment the line
        * counter, i. */
       if (++i >= MAXRANGES) {
	 fprintf(stderr,"\nExceeded array size for number of ranges.");
	 fprintf(stderr,"\nChange #define MAXRANGES and recompile.\n");
	 exit(1);
       }
   }
   
   close(rangefile);
   
   if (i == 0) {
     /* No ranges read from rangefile, return -1. */
      fprintf(stderr,"\nNo ranges were read in.");
      i = -1;
   }

   /* Return the number of lines read to invoking program. */
   return (i);

}  /* end load_ranges() */

/****************************************************************************/

int find_index(double depth, int n) {

/* Determines the index corresponding to the dmin/dmax pair in
 * which the specified depth falls. Returns the index number, or
 * -1 if the depth is not in any range. */

   int i;

   i = -1;
   while (++i < n) {
     if (depth >= dmin[i] && depth < dmax[i])
        return (i);
   }
   return (-1);

}  /* End find_index() */

/****************************************************************************/

int find_bin(double depth, int n) {

/* Determines the index corresponding to the
 * zbin edges in which the specified depth
 * falls. Returns the index number, or
 * -1 if the depth is not in any range. */
   int i;

   i = -1;
   while (++i < n) {
     if (depth >= zbin[i] && depth < zbin[i+1])
        return (i);
   }
   return (-1);

}  /* End find_bin() */

/*********************************************************************/
