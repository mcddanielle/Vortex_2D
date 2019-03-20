# Vortex_2D

### Author Danielle McDermott, November 1, 2013
 
compile with Makefile

2019: it appears I had some trouble linking certain libraries, so I cut them and created basic_src to avoid those calls.

 
use --help flag to view following docs:
-------------------------------------------------

---------------------------------------------
This help menu is for MovieConverter
This needs an update for analysis---------------------------------------------
Read binary movie files to extract single frames of a movie
Author Danielle McDermott 
Last Update June 2015
Usage basic_vortex_geometry
Input/Output Files
	-f|--name of binary movie to be read
	-o|--name of file to be written
---------------------------------------------
Read Options
	-m|--movie type = -1,1, default -1 
	  -1: Cynthia's Original Format 
	   0: nada  1: Ascii 
---------------------------------------------
Write Options
	-w|--output type = -1,0,1
	  -1: Delplot Format (input to delplot or vortexsolid.c (ksttestn0))
	   0: Version0 format (input to VortexSolidPostProcessor)
	   1: Ascii Format
The final frame is always written.  Other frames can be processed with the -I flag to set the Fd (current) of the output.  This will not be suitable for a movie of version -1
---------------------------------------------
Other options
	-X|--xmax (required for movie type -1 (see Read Options)) 
	-Y|--ymax (required for movie type -1 (see Read Options)) 
	-N| number of pins 
	-p| maximum pinning force 
	--verbose print out debugging statements
	--quiet don't print normal output
-------------------------------------------------

 
//-----------------------------------------------------------------------
 *Software requires an input movie name to run and will assume a 
 *system size of 36X36 unless told otherwise using the command line
 *argument -s
 *
 *given the argument -m 1, rather that convert between smtest and a Jeff-movie, 

*the code will simply read in a Jeff-movie and write its final frame to a file
 *
 *TODO: Include a "write ascii" option
 *TODO: Use strtok() to parse the filenames, right now it will only
 *      accept filenames in the given directory
 *TODO: use the new subroutine "print_frame" to do all frame printing
 */

//TO RUN
// movie_convert -f N4.3_drive-down -I 0.15 -H 1 -s 72.0 -m 1
//I = current
//H =1 turns on histograms
//m =1 implies a Jeff movie, which goes with the histograms

//-----------------------------------------------------------------
LENGTH SCALES OF VORTEX GEOMETRY
//-----------------------------------------------------------------
1) Pinning period - command line read
2) xi = 1.0  attractive length scale
3) not long-range cutoff, but roughly r=8.0