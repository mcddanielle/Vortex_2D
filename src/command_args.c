#include "Vortex_2D.h"

//---------------------------------------------------------//
//                 Utility Subroutines                     //
//---------------------------------------------------------//
// Process the command line arguments. All global variables so 
// there's no sneaking around passing notes back and forth
void getoptions( int argc, char **argv ) {
    int c;
    int local_quiet = 1;
  while (1) {
    static struct option long_options[] = {
      /*-- These options set a flag. --*/
      {"verbose", no_argument,       &verbose_flag, 1},
      {"quiet",   no_argument,       &quiet_flag,   1},

      /*-- These options don't set a flag. We distinguish them by their indices. --*/
      {"xmax",         required_argument, 0, 'X'},
      {"ymax",         required_argument, 0, 'Y'},
      {"file",         required_argument, 0, 'f'},
      {"outfile",      required_argument, 0, 'o'},
      {"writefiletype",required_argument, 0, 'w'},
      {"movietype",    required_argument, 0, 'm'},
      {"number_pins",  required_argument, 0, 'N'},
      {"pinning_force",required_argument, 0, 'p'},
      //{"hist_range",   required_argument, 0, 'R'},
      //{"bin_zeros",    required_argument, 0, 'Z'},
      {"help",         no_argument,       0, 'h'},

      {0, 0, 0, 0}
    };

    /*-- getopt_long stores the option index here. --*/
    int option_index = 0;
    //dan informs me that the colon denotes a required argument
    //const char * short_options = "s:o:f:m:r:I:H:R:Z:h:w";
    const char * short_options = "X:Y::o:f:m:I:H:R:Z:hw:N:p:";

    c = getopt_long (argc, argv, short_options, long_options, &option_index);
     
    /*-- Detect the end of the options. --*/
    if (c == -1)
      break;

    switch (c) {
    case 0:  /*-- If this option set a flag just keep going --*/
      if (long_options[option_index].flag != 0)
	break;

    case 'h':
      // Someone needs a little help provide no information to help them out
      printf("---------------------------------------------\n");
      printf("This help menu is for MovieConverter\n");
      printf("This needs an update for analysis");
      printf("---------------------------------------------\n");
      printf("Read binary movie files to extract single frames of a movie\nAuthor Danielle McDermott \nLast Update June 2015\n");
 printf("Usage %s\n", argv[0]); 
      printf("Input/Output Files\n");
      printf("\t-f|--name of binary movie to be read\n");
      printf("\t-o|--name of file to be written\n");
      printf("---------------------------------------------\n");
      printf("Read Options\n");
      printf("\t-m|--movie type = -1,1, default -1 \n");
      printf("\t  -1: Cynthia's Original Format \n\t   0: nada  1: Ascii \n");
      printf("---------------------------------------------\n");
      printf("Write Options\n");
      printf("\t-w|--output type = -1,0,1\n");
      printf("\t  -1: Delplot Format (input to delplot or vortexsolid.c (ksttestn0))\n\t   0: Version0 format (input to VortexSolidPostProcessor)\n\t   1: Ascii Format\n");
      printf("The final frame is always written.  Other frames can be processed with the -I flag to set the Fd (current) of the output.  This will not be suitable for a movie of version -1\n");
     
      printf("---------------------------------------------\n");
      printf("Other options\n");

      printf("\t-X|--xmax (required for movie type -1 (see Read Options)) \n");
      printf("\t-Y|--ymax (required for movie type -1 (see Read Options)) \n");
      printf("\t-N| number of pins \n");
      printf("\t-p| maximum pinning force \n");
      printf("\t--verbose print out debugging statements\n");
      printf("\t--quiet don't print normal output\n");
      exit(0);

    case 'f':
      strcpy(filename,optarg);     
             if(!local_quiet) printf("file name is: %s\n",filename); 
      break;

    case 'o':
      strcpy(out_file,optarg);     
             if(!local_quiet) printf("the output file name is: %s\n",out_file); 
      break;

    case 'X':
      xmax = atof(optarg);     
             if(!local_quiet) printf("xmax is: %f\n",xmax); 
      break;

    case 'Y':
      ymax = atof(optarg);     
             if(!local_quiet) printf("ymax is: %f\n",ymax); 
      break;

    case 'm':
      movie_type = atoi(optarg);     
             if(!local_quiet) printf("the input movie type is: %d where smtest=-1, Jeff movie=0,1,2 \n default = 0 \n",movie_type);
      break;

    case 'N':
      number_pins = atoi(optarg);     
             if(!local_quiet) printf("the analytic pinning has %d columns \n",number_pins);
      break;

    case 'p':
      pinning_force = atof(optarg);     
             if(!local_quiet) printf("the maximum pinning force = %f \n",pinning_force);
      break;
      
    case 'w':
      write_movie_type = atoi(optarg);     
             if(!local_quiet) printf("the output_frame is: %d where smtest=-1, Jeff movie=0, Ascii=1 \n default = 0 \n",write_movie_type);
      break;

    case 'r':
      restart_flag = atoi(optarg);     
             if(!local_quiet) printf("restart = %d where Jeff movie=0 and ksttestn0=1 \n",restart_flag);
      break;
      
    default:
      printf("Ignoring that one\n");

    } // End switch(c)

  } // End while(1)

    return;
} // End void getoptions( int argc, char **argv )
