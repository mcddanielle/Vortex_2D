# Vortex_2D
ABOUT
/*Author Danielle McDermott, November 1, 2013
 */

 *
 *compile with Makefile
 *

 *BIG GOAL: to combine movie reading of MovieConverter
 *      with the geometry processing of Nanoparticle
 *

 *from Jeff's source code:*-----------------------------
   SMTEST movie FRAME formats:
  -1: has ids, positions
 *----------------------------------

 /*----------------------November 1st, 2013 GOAL-----------------------*/
 1) Working Code!  defined by: 
 e) then the ensuing sorting, nn, g(r), q6, etc may begin
 f) linked list + malloc the bt[C] array in vortex_list
 
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