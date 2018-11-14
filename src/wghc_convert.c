/*  wghc_convert.c
................................................................................
.   Reads WOCE global hydrographic climatology (Gouretski) observed level hydrographic data files
.   extracts :    header info
.                 p,d,t,s,ox,n2,n3,si,p4
.   
.   USAGE: wghc_convert infile_list -B<badfile> -O<outfile> -T<instr_type> -Sship_lookup_file -Ccountry_code_file  [-D<dir>] [-E<extent>]
...............................................................................
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "hydrobase.h"

struct S_REC {
    char shipname[30];
    char shipcode[3];
    char ccode[3];
    struct S_REC *next;
};

struct C_REC {
    char name[30];
    char code[3];
    struct S_REC *shipptr;
    struct C_REC *next;
};
  /*  HydroBase globally defined variables */
  
struct HYDRO_HDR  hdr, bad_hdr;
struct HYDRO_DATA data, bad_data;

char *buffer;
int n_sa, n_te, n_ox, n_si, n_p4, n_n3;
int bopt;
struct C_REC *table;

  /*  prototypes for locally defined functions */

void print_usage(char *);
FILE *openfile(char *, char *, char *);
int compare(char *, char *) ;
struct C_REC *build_cc_lookup(char *);
void build_shipcode_lookup(char *, struct C_REC *);
struct C_REC * find_country(struct C_REC *);
struct S_REC * find_ship(struct S_REC *);
int wghc_read(FILE *);
char *get_code(char *, char *);
void check_sta(int, int);

main (int  argc, char **argv)
{

   int error, nfiles, curfile, status;
   int i, j,  staread, staout, stabad;
   int outfile, badfile;
   char *shipfilename, *country_filename;
   char  *dir, *extent, *st;
   FILE *fp;
   struct C_REC  *country;

/* check for command line arguments */

   if (argc < 2) {
     print_usage(argv[0]);
     exit(1);
   }

/*  set these default values */

    dir = "";
    extent = "";
    bopt = 0;
    error = 0;
    nfiles = 0;
    curfile = 1;
    staout = staread = stabad = 0;
    outfile = STDOUT;
    shipfilename = country_filename = NULL;
    buffer = NULL;
    for (i = 0; i < MAXPROP; ++i) {
       data.observ[i] = NULL;
       bad_data.observ[i] = NULL;
   }
    hdr.prop_id = NULL;
    bad_hdr.prop_id = NULL;
 
/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'B':                    /* open  file for output */
                        badfile = create_hydro_file(&argv[i][2], OVERWRITE);
                        if (badfile < 1) {
                           fprintf(stderr,"\nError opening output file: %s\n", 
                                  &argv[i][2]);
                          exit(1);
                        }
                        fprintf(stderr,"\nBad data will be written to: %s", &argv[i][2]);

			bopt = 1;
                        break;
               case 'C':                    /* get country code file  */
                        country_filename = &argv[i][2];
                         break;
               case 'D':                   /* get input dir */
                        dir = &argv[i][2];
                        break;
               case 'E':                    /* get file extent */
                        extent = &argv[i][2];
                        break;
               case 'O':                    /* get output file  */
                        outfile = create_hydro_file(&argv[i][2], OVERWRITE);
                        if (outfile < 1) {
                           fprintf(stderr,"\nError opening output file: %s\n", 
                                  &argv[i][2]);
                           exit(1);
                        }
                        fprintf(stderr,"\nOpening or overwriting: %s", &argv[i][2]);
                        break;
               case 'S':                    /* get shipcode file  */
                        shipfilename = &argv[i][2];
                         break;
      	       case 'h':  
	           print_usage(argv[0]);
	           exit(0);

             default:
                        error = 1;

          }    /* end switch */

          if (error ) {
             fprintf(stderr,"\nError parsing command line args.\n");
             fprintf(stderr,"     in particular: '%s'\n", argv[i]);
             exit(1);
          }

       }  /* end if */

       else  {
           ++nfiles;
       }
   }  /* end for */
   
   if (shipfilename == NULL) {
       fprintf(stderr,"\nYou must specify the path to OCL/NODC shipcode lookup file [ship_lookup]. \n");
       exit(1);
   }
   if (country_filename == NULL) {
       fprintf(stderr,"\nYou must specify the path to OCL/NODC country code file [country_codes]. \n");
       exit(1);
   }
   if (! bopt) {
       fprintf(stderr,"\nBad stations will not be saved in a separate file\n");
   }
   
   if (! nfiles) {
       fprintf(stderr,"\nExpecting data to come from STDIN...");
       fp = stdin;
   }
   
/*  Read in files to create a lookup table for country/ship codes */

   table = build_cc_lookup(country_filename);
   build_shipcode_lookup(shipfilename, table );
   
    
 /*  loop for each input file */
 
 
  do {
  
     if (nfiles) {
  	fp = openfile(dir, argv[curfile], extent);
  	if (fp == NULL) 
        	goto NEXTFILE;
     }
     /* loop for each station */
     
      do  {
     
            status = wghc_read(fp);
            ++staread;
	    
	   if (hdr.nobs > 0) {        
              write_hydro_station(outfile, &hdr, &data);
              ++staout;
	   }
	   
	   if (bopt && (bad_hdr.nobs > 0)) {
	      write_hydro_station(badfile, &bad_hdr, &bad_data);
              ++stabad;
	   }
        
           if (hdr.prop_id != NULL) {
	       free(hdr.prop_id);
	       hdr.prop_id = NULL;
	    }
	          
     } while (status != EOF );
       
NEXTFILE:
     if (nfiles) 
       fclose(fp);
  } while   (curfile++ < nfiles);       

   fprintf(stderr,"\nEnd of conversion.");
   fprintf(stderr,"\n  %d stations read in", staread);
   fprintf(stderr,"\n  %d stations accepted", staout);
   fprintf(stderr,"\n  %d stations rejected ", staread-staout);
   fprintf(stderr,"\n  %d stations contained bad scans \n\n", stabad);
   exit(0);

} /* end main() */ 
  
/*****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"%s converts QC'd observed profile data from the WOCE Global Hydrographic Climatology ", program);
   fprintf(stderr,"((Gouretski & Koltermann, 2004) into HydroBase2 station format. \n");
   fprintf(stderr,"\n\nUSAGE:  %s infile_list -C<country_code_file> -S<shipcode_file> [-O<outfile>]  [-B<bad_file>]  [-D<dir>] [-E<extent>] ");
   fprintf(stderr,"\n List of filenames must be first arguments");
   fprintf(stderr,"\n     -C:  supply country code filename ");
   fprintf(stderr,"\n     -S:  specify ship code filename ");
   fprintf(stderr,"\n OPTIONS:     ");
   fprintf(stderr,"\n    [-B]:  output file for bad observations ");
   fprintf(stderr,"\n    [-D]:  directory for input files ");
   fprintf(stderr,"\n    [-E]:  extent for input files ");
   fprintf(stderr,"\n    [-O]:  specify output file for HydroBase station files ");
   fprintf(stderr,"\n               --or stations are output to STDOUT");
   fprintf(stderr," \n");
   return;
}
   
/*****************************************************************************/
/*****************************************************************************/
FILE *openfile(char *dir, char *root, char *extent)
{
   char st[80];
   int i;
   FILE *infile;
   
   strcpy(st, dir);
   if ((i=strlen(dir)) != 0) {
      if (dir[i-1] != '/')
         strncat(st,"/",1);
   }
   strcat(st, root);
   strcat(st, extent);
   infile = fopen(st,"r");
   if (infile == NULL)
         fprintf(stderr,"\n Unable to open %s \n", st);
   else
         fprintf(stderr,"\n   opened %s ... ", st);
   
   return(infile);
   
}  /* end openfile() */
/*****************************************************************************/
struct C_REC * build_cc_lookup(char *filename)

   /*  Uses information in file to build an alphetical 
      linked list of countries, forming the
      basis for a lookup table of country and ship codes.
      Returns a pointer to the start of list or NULL */
{
   int n, count;
   FILE *fptr;
   char  *line, *st;
   struct C_REC *listptr, *cptr, *rptr, *prevptr;
   
   
   fptr = openfile("", filename, "");
   if (fptr == NULL) {
      exit (1);
   }
   line = (char *) calloc(100, sizeof(char));
   listptr = (struct C_REC *) NULL;
   count = 0;
   while (( n = fscanf(fptr, "%[^\n]", line ) ) == 1) {
       fgetc(fptr);   /* move past LF character */
       ++count;
        cptr = (struct C_REC *)calloc(1, sizeof(struct C_REC));

        st = line;
        strncpy(cptr->code, st, 2);
        ++st;
        ++st;
     
        while (*st == ' ')
             ++st;
	  
        strcpy(cptr->name, st );
        cptr->shipptr = NULL;
     
     /*insert record into linked list */
     
        if (listptr == NULL) {    /* empty list */
            cptr->next = NULL;
	    listptr = cptr;
        }
        else {
	 
	    rptr = listptr;
	    prevptr = NULL;
	 
	 /* traverse linked list to find place to insert this record */
	 
	    while  (  (rptr != NULL) && (compare(cptr->name, rptr->name) > 0)  ) {
	        prevptr = rptr;
	         rptr = rptr->next;
	     }
	  
	     if (rptr == NULL) {       /*insert at end of list */
	           prevptr->next = cptr;
		   cptr->next = NULL;
	     }
	     else if (prevptr == NULL ) {  /* insert record at beginning of list */
	            cptr->next = listptr;
		    listptr = cptr;
	     }    
	     else {                          /*insert before rptr */
	          cptr->next = rptr;
		  prevptr->next = cptr;
	     }
        } /*end else */
   
   } /* end while */
   
   if ( n == EOF) {
          fclose(fptr);
	  free(line);
          return (listptr);
   }
	  
   fprintf(stderr,"\nError reading %s  ", filename);
   fprintf(stderr,"\ncount = %d   ", count);
   if (count)
      fprintf(stderr," last country read is [%s] \n ", prevptr->name);
   exit(1);

} /* end build_cc_lookup() */
/*****************************************************************************/

int compare(char *key, char *str) 
      /*  Compares 2 strings recursively.  
         Returns:
           1 if key comes alphabetically after str
	  -1 if key comes before
	   0  if they match
    Individual characters are compared until they differ or end-of-string is encountered 
    in either.  When one is at the end-of-string, if all previous chars matched, then
    0 (match) is returned.
	  */
{
   char *k, *s;

   k = key;
   s = str;
   
   if (*k == '\0' || *s == '\0')    /* if we're at the end of string, they match */
        return 0;
	
   if ( (int) *k > (int) *s )
       return 1;
   if ( (int) *k < (int) *s )
       return -1;
       
    return (compare(++k, ++s));   	
   
}   /* end compare () */

/*****************************************************************************/
void build_shipcode_lookup(char *fname, struct C_REC *start_of_table)
   /*  Opens ship_lookup file and adds info to linked list at appropriate country.
   */
{
   FILE * fptr;
   int n, match;
   char  *line, *st, *xx;
   struct C_REC *cptr;
   struct S_REC  *sptr, *rptr, *prevptr;
 
  fptr = openfile("", fname, "");
   if (fptr == NULL) {
      exit (1);
   }
   line = (char *) calloc(200, sizeof(char));
   
   while (( n = fscanf(fptr, "%[^\n]", line )) == 1) {
       fgetc(fptr);   /* move past LF character */
   
      if (line[0] == '*' && line[1] == '*') {    /* country */
      
           st = &line[2];   /* move past asterisks  and white space*/
	   while (*st == ' ')
	        ++st;  
		
	   cptr = start_of_table;
	   while (!(match = (compare(st, cptr->name) == 0))  && (cptr->next != NULL))
	        cptr = cptr->next;
	    
	    if (!match) {
	       fprintf(stderr, " \n Unrecognized country [%s]\n", st );
	       exit(1);
	    }
	   	      
      }
      else  {      /*ship */
           sptr = (struct S_REC *) calloc(1, sizeof(struct S_REC));
	   strncpy(sptr->ccode, line, 2);
	   strncpy(sptr->shipcode, &line[2], 2);
	   st = &line[5];
	   while (*st == ' ')
	        ++st;
		
	   xx = st;	  /* remove any parentheses, commas, periods */
	   while (*xx != '\0') {
	      if (*xx == '(' || *xx == ',' || *xx == '.')
	         *xx = '\0';
	      else
	         ++xx;
	   }
	   
	   /* xx should now be a null (endOfString) value. 
	  Remove any trailing white space */
	   --xx;    
	   while  ( *xx == ' ') {   
	      *xx = '\0';
	      --xx;
	   }
	   
	   strcpy(sptr->shipname, st);
	    
	    /*insert this record alphabetically into the linked list */
	    
	    if (cptr->shipptr == NULL) {  /* empty list */
	       sptr->next = NULL;
               cptr->shipptr = sptr;
	    }
	    else {
	        
	         rptr = cptr->shipptr;
	         prevptr = NULL;
	 /* traverse linked list to find place to insert this record */
	         while ( (rptr != NULL) && (compare(sptr->shipname, rptr->shipname) > 0 )) {
	            prevptr = rptr;
	            rptr = rptr->next;
	         }
		 
		 if (rptr == NULL) {  /* insert at end of list */
		     prevptr->next = sptr;
		     sptr->next = NULL;
		 }
		 else if (prevptr == NULL) {   /* insert at beginning of list */
		     sptr->next = cptr->shipptr;
                     cptr->shipptr = sptr;
		 }
	         else {                          /*insert before rptr */
	              sptr->next = rptr;
		      prevptr->next = sptr;
	         }
	    
	    } /* end else */
      } /* else ship */
   
   } /* end while */
   
   if ( n == EOF) {
          fclose(fptr);
	  free(line);
          return;
   }
   
   /* if we get here, an error occurred... */
    fprintf(stderr,"\nError reading %s  ", fname);
    fprintf(stderr,"\nlast line read is [%s] \n ", line);
    exit(1);

}  /*end build_shipcode_lookup() */
/*****************************************************************************/
char *get_code(char *ship, char *nation)
   /* Looks up nation and ship in global table, and returns a 4-char string
      containing the OCL codes.  Returns NULL if country cannot be found.
      Returns country code and shipcode = XX if only the ship cannot be found.
    */
{
  struct C_REC *cptr;
  struct S_REC *sptr;
  int found, status;
  char *st;
  
  if ((cptr = table) == NULL) {
      fprintf(stderr,"\nUnable to access lookup table\n");
      exit(1);
  }
  
  /* search for matching nation */
  found = 0;
  while (!found) {
     status = compare(nation, cptr->name);
     if (status == 0)  {
	  sptr = cptr->shipptr;
          found = 1;
     }
     else if (status <  0 || cptr->next == NULL) {
          fprintf(stderr,"Unknown country: [%s]\n", nation);
	  return NULL;
     }
     else
        cptr = cptr->next;
  }  /* end while */
  
  /* store country code  */
  
      st = (char *) calloc(5, sizeof(char));
      *st =  cptr->code[0];
      *(st+1) = cptr->code[1];

  if (sptr == NULL) {
      fprintf(stderr,"\nFound country(%2s), but unable to access ship list.\n", cptr->code);
       return(st);
  }
   /* search for matching shipname */
 
  while (sptr != NULL) {
     status = compare(ship, sptr->shipname);
     if (status == 0) {
          *(st+2) = sptr->shipcode[0];
          *(st+3) = sptr->shipcode[1];
           return(st);
     }
     if (status < 0 || sptr->next == NULL) {
          fprintf(stderr,"Unknown ship: [%s]\n", ship);
         *(st+2) = 'X';
         *(st+3) = 'X';
           return(st);
     }
     sptr = sptr->next;
  } /* end while */
  
  fprintf(stderr,"\n ERROR in get_code():  should never get this far.\n");
  return NULL;

}  /* end get_code() */

/*****************************************************************************/
int wghc_read(FILE *fptr)
   /*  Upon entry, buffer is NULL if first attempt to read this file.
      Otherwise, buffer holds last line read.
      Reads one entire station profile  (until new station or end_of_file encountered)
       with all its properties and flags into global Hydrobase structures 
      Returns status = 1 for successful read, 
                       EOF for end of file 
	An error reading or parsing causes an error message and exit from program.
 */
{
   char *ship, *country, *st;
   char err_str1[10], err_str2[10], err_str3[10];
   int end_of_sta, status, wghc_id, wghc_id2;
   int   i, count, badcount, propindex;
   float blah, bias_sa, bias_ox, bias_si, bias_n3, bias_p4;

   n_sa = n_te = n_ox = n_si = n_p4 = n_n3 = 0;
   
/* If buffer is empty, this is first line of file */   

   if (buffer == NULL) {
      buffer = (char *) calloc( 200, sizeof(char));
      status =  fscanf(fptr, "%[^\n]", buffer );
      if (status == EOF)
          return(EOF);
      if (status != 1) {
         fprintf(stderr,"Error reading first line of  wghc file.");
         exit (1);
      } 
      fgetc(fptr);   /* move past LF character */
    }  


   /* alloc space in global structures */
    
    hdr.prop_id = (int *) calloc(8, sizeof(int));
    hdr.nprops = 8;
    hdr.prop_id[0] = (int) PR; 
     hdr.prop_id[1] = (int)DE;
     hdr.prop_id[2] = (int)TE;
     hdr.prop_id[3] = (int)SA;
     hdr.prop_id[4] = (int)OX;
     hdr.prop_id[5] = (int)N3;
     hdr.prop_id[6] = (int)P4;
     hdr.prop_id[7] = (int)SI;
     
     for (i = 0; i < hdr.nprops; ++i) {
        free_and_alloc(&data.observ[hdr.prop_id[i]], 3000);
	if (bopt)
            free_and_alloc(&bad_data.observ[hdr.prop_id[i]], 3000);
     }
    count = 0;  
    badcount = 0;
   end_of_sta = 0;
   
   while (!end_of_sta) {
      
      switch (buffer[0] )  {
   
          case 'A':
	           ship  = &buffer[13];
		       /* remove trailing blanks */
		   st = &ship[30];
		   while (*st == ' ' && st >= ship) {
		      *st = '\0';
		      --st;
		   }
		   
		   country = &buffer[44];
		          /* remove trailing blanks */
		   st = &country[15];
		   while (*st == ' ' && st >= country) {
		      *st = '\0';
		      --st;
		   }
		   hdr.country[0] = 'X';
		   hdr.country[1] = 'X';
		   hdr.country[2] = '\0';
		   hdr.ship[0] = 'X';
		   hdr.ship[1] = 'X';
		   hdr.ship[2] = '\0';
		   
		   st = get_code(ship, country);
		   if (st != NULL) {
		        strncpy(hdr.country, st, 2);
		        strncpy(hdr.ship, &st[2], 2);
			free((char *)st);
		   }
	          status = sscanf(&buffer[2], "%d", &wghc_id);
		   
	      break;
	      
          case 'B': 
	          status = sscanf(&buffer[2], "%d", &wghc_id2);
		  if (status != 1) {
                       fprintf(stderr,"Error parsing wghc cruise id.");
		       fprintf(stderr, "\n%s\n", buffer);
                       exit (1);
		  }
		  
	          status = sscanf(&buffer[13], "%d", &hdr.station);
		  if (status != 1) {
                       fprintf(stderr,"Error parsing station id.");
		       fprintf(stderr, "\n%s\n", buffer);
                       exit (1);
		  }
		  if (hdr.station > 9999) {
		     hdr.station = hdr.station % 1000;
		  }
		  if (hdr.station <= 0) {
		     hdr.station = wghc_id % 1000;
		  }
		  switch (buffer[27]) {
		      case 'N':
		      case 'R':
		      case 'X':
		             hdr.instrument = 'b';
			     break;
	              case 'C':
		             hdr.instrument = 'c';
			     break;
		      default:
		             hdr.instrument = ' ';    
			     
		  } /* end switch */
		  
		  hdr.origin = '8';
		  
	          status = sscanf(&buffer[35], "%d", &hdr.cruise);
		  if (status != 1) {
                       fprintf(stderr,"Error parsing cruise id.");
		       fprintf(stderr, "\n%s\n", buffer);
                       exit (1);
		  }
		  if (hdr.cruise > 99999) {
		     hdr.cruise = hdr.cruise % 10000;
		  }
		  
		  if (hdr.cruise <= 0) {
		       hdr.cruise = wghc_id2 % 10000;
		  }
		  
	      break;
	      
          case 'C':
	          status = sscanf(&buffer[2], "%d %d %d", &hdr.year, &hdr.month, &hdr.day);
		  if (status != 3) {
                       fprintf(stderr,"Error parsing year, month, day");
		       fprintf(stderr, "\n%s\n", buffer);
                       exit (1);
		  }
		  
	          status = sscanf(&buffer[22], "%f %f", &hdr.lat, &hdr.lon);
		  if (status != 2) {
                       fprintf(stderr,"Error parsing lat, lon");
		       fprintf(stderr, "\n%s\n", buffer);
                       exit (1);
		  }
		   if (hdr.cruise <= 0) {
		     hdr.cruise =  (hdr.year - (hdr.year % 100)) + hdr.month;
		   }
	      break;
	      
          case 'D': 
		  
	          status = sscanf(&buffer[36], "%d", &hdr.pdr);
		  if (status != 1) {
                       fprintf(stderr,"Error parsing bottom depth");
		       fprintf(stderr, "\n%s\n", buffer);
                       exit (1);
		  }
		  
	      break;
	      
          case 'E':
		  status = sscanf(&buffer[2], "%f %f %f %f %f %f %f %f %f %f %f %f %f", &blah, &blah, &blah, &blah, &bias_sa, &blah, &blah, &blah, &blah, &bias_ox, &bias_si, &bias_n3, &bias_p4);
		  if (status != 13) {
                       fprintf(stderr,"Error parsing biases");
		       fprintf(stderr, "\n%s\n", buffer);
                       exit (1);
		  }
	      break;
	      
          case '/':     /* skip over comment */
	      break;
          default:
	        /* data line*/
		
		  status = sscanf(buffer, "%lf %lf %lf %f %lf %f %f %f %f %lf %lf %lf %lf %s %s %s", &data.observ[(int)DE][count], &data.observ[(int)PR][count], &data.observ[(int)TE][count], &blah, &data.observ[(int)SA][count], &blah, &blah, &blah, &blah, &data.observ[(int)OX][count], &data.observ[(int)SI][count], &data.observ[(int)N3][count], &data.observ[(int)P4][count], err_str1, err_str2, err_str3);
		  if (status != 16) {
                       fprintf(stderr,"Error parsing properties");
		       fprintf(stderr, "\n%s\n", buffer);
                       exit (1);
		  }
		  
		  st = err_str1;
		  if (st[0] == '1' ) data.observ[(int)DE ][count] = HB_MISSING;
		  if (st[1] == '1' ) data.observ[(int)TE ][count] = HB_MISSING;
		  if (st[2] == '1' ) data.observ[(int)SA ][count] = HB_MISSING;
		  if (st[3] == '1' ) data.observ[(int)OX ][count] = HB_MISSING;
		  if (st[4] == '1' ) data.observ[(int)SI ][count] = HB_MISSING;
		  if (st[5] == '1' ) data.observ[(int)N3 ][count] = HB_MISSING;
 		  if (st[6] == '1' ) data.observ[(int)P4 ][count] = HB_MISSING;
		  st = err_str2;
		  if (st[0] == '1' ) data.observ[(int)DE ][count] = HB_MISSING;
		  if (st[1] == '1' ) data.observ[(int)TE ][count] = HB_MISSING;
		  if (st[2] == '1' ) data.observ[(int)SA ][count] = HB_MISSING;
		  if (st[3] == '1' ) data.observ[(int)OX ][count] = HB_MISSING;
		  if (st[4] == '1' ) data.observ[(int)SI ][count] = HB_MISSING;
		  if (st[5] == '1' ) data.observ[(int)N3 ][count] = HB_MISSING;
 		  if (st[6] == '1' ) data.observ[(int)P4 ][count] = HB_MISSING;

   
		  st = err_str3;
		  if (st[0] == '1' ) data.observ[(int)DE ][count] = HB_MISSING;
		  if (st[1] == '1' ) data.observ[(int)TE ][count] = HB_MISSING;
		  if (st[2] == '1' ) data.observ[(int)SA ][count] = HB_MISSING;
		  if (st[3] == '1' ) data.observ[(int)OX ][count] = HB_MISSING;
		  if (st[4] == '1' ) data.observ[(int)SI ][count] = HB_MISSING;
		  if (st[5] == '1' ) data.observ[(int)N3 ][count] = HB_MISSING;
 		  if (st[6] == '1' ) data.observ[(int)P4 ][count] = HB_MISSING;
		  
                  if (bias_sa > -8.0 ) 
		      data.observ[(int)SA ][count] -= bias_sa;
                  if (bias_si > 0. && data.observ[(int)SI ][count] > 0)
		      data.observ[(int)SI ][count] *= bias_si;
                  if (bias_ox > 0. && data.observ[(int)OX ][count] > 0)
		      data.observ[(int)OX ][count] *= bias_ox;
                  if (bias_p4 > 0. && data.observ[(int)P4 ][count] > 0)
		      data.observ[(int)P4 ][count] *= bias_p4;
                  if (bias_n3 > 0. && data.observ[(int)N3 ][count] > 0)
		      data.observ[(int)N3 ][count] *= bias_n3;
			   
		   if ((data.observ[(int)DE ][count] >= 0.0) &&  (data.observ[(int)SA ][count] > 0.0) && (data.observ[(int)TE ][count] > -3.0) ) {
		         ++n_sa;   
		         ++n_te;   
		      if (data.observ[(int)OX ][count] > 0.0)
		         ++n_ox;   
		      if (data.observ[(int)P4 ][count] >= 0.0)
		         ++n_p4;   
		      if (data.observ[(int)SI ][count] >= 0.0)
		         ++n_si;   
		      if (data.observ[(int)N3 ][count] >= 0.0)
		         ++n_n3; 
		      ++count;  
		   }	 
		    else {
		       
		       if (bopt) {
                             status = sscanf(buffer, "%lf %lf %lf %f %lf %f %f %f %f %lf %lf %lf %lf", &bad_data.observ[(int)DE][badcount], &bad_data.observ[(int)PR][badcount], &bad_data.observ[(int)TE][badcount], &blah, &bad_data.observ[(int)SA][badcount], &blah, &blah, &blah, &blah, &bad_data.observ[(int)OX][badcount], &bad_data.observ[(int)SI][badcount], &bad_data.observ[(int)N3][badcount], &bad_data.observ[(int)P4][badcount]);
			    ++badcount;
		       }
		    }
      } /* end switch */
     
      status =  fscanf(fptr, "%[^\n]", buffer );
     
      if (status == EOF) {
          free(buffer);
	  buffer = (char *) NULL;
	  check_sta(count, badcount);
	  return(EOF);
      }
      if (status != 1) {    /* error */
         fprintf(stderr,"Error reading wghc file: \n%s\n", buffer);
         exit (1);
      }
            
      fgetc(fptr);  /* buffer contains a new line, move past LF */
      end_of_sta = buffer[0] == 'A' ;
   
   }  /* end while !end_of_sta */
   
   check_sta(count, badcount);
    return(1);
}  /* end wghc_read() */
/*****************************************************************************/
void check_sta(int count, int badcount)
   /* Gets station ready for output.  Determines what properties have valid observations
      and adjusts hdr (and bad_hdr) as appropriate.  */
{
   int i, npropsout;
   
   hdr.nobs = count;
   
   if (count) {    /* pr, de, te, sa -- so far */
       npropsout = 4;
       if ( n_ox) {
          hdr.prop_id[npropsout] = (int)OX;
           ++npropsout;
       }
       if ( n_n3) {
          hdr.prop_id[npropsout] = (int)N3;
           ++npropsout;
       }
        if ( n_p4) {
          hdr.prop_id[npropsout] = (int)P4;
           ++npropsout;
       }
       if ( n_si) {
          hdr.prop_id[npropsout] = (int)SI;
           ++npropsout;
       }
       
       hdr.nprops = npropsout;
       data.nprops = npropsout;
       data.nobs = hdr.nobs;
 
       hdr.qual[0] = '0';
       hdr.qual[1] = '0';
       hdr.qual[2] = '0';
       hdr.qual[3] = '1';
       if (n_n3 || n_si || n_p4)
          hdr.qual[1] = '1';
       if (n_ox)
          hdr.qual[2] = '1';
	  
       hdr.ms10 = ms10(hdr.lat, hdr.lon, &hdr.ms1);
   }  

   if (!bopt)  return;
   
   /* check and adjust info for writing bad observations */
      
   bad_hdr.nobs = badcount;
   if (badcount) {
       bad_hdr.prop_id = hdr.prop_id;
       bad_hdr.nprops = hdr.nprops;
       bad_hdr.ms10 = ms10(hdr.lat, hdr.lon, &bad_hdr.ms1);
       bad_hdr.pdr = hdr.pdr;
       strncpy(bad_hdr.country, hdr.country, 2);
       strncpy(bad_hdr.ship, hdr.ship, 2);
       bad_hdr.origin = hdr.origin;
       bad_hdr.instrument = hdr.instrument;
       bad_hdr.cruise = hdr.cruise;
       bad_hdr.station = hdr.station;
       bad_hdr.year = hdr.year;
       bad_hdr.month = hdr.month;
       bad_hdr.day = hdr.day;
       bad_hdr.lat = hdr.lat;
       bad_hdr.lon = hdr.lon;
       bad_data.nobs = bad_hdr.nobs;
       bad_data.nprops = bad_hdr.nprops;
   }
   return;

} /* end check_sta */

