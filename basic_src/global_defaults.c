#include "Vortex_2D.h"

void assign_global_defaults(){

// Global variables -- Default values given here
// getoptions() assigns runtime values

  //flag -X and -Y//
  xmin = 0.0;   
  ymin = 0.0;
  xmax = 36.5;
  ymax = 36.5;
  xsi  = 1.0;       // hardcoded as length scale, 
                    //but xsi = 1.0, kappa = \xsi/\lambda
  bondLim = 1.75;   //the zero pinning case roughly 1.2*r where F(r) = 0.0

  number_pins = 0;
  pinning_force = 0.0;
  //flag -I //
  frame_print_current = 0.0;  //frame printing + 2D histogramming

  verbose_flag = 0;   // Print out lots of debugging statements if set
  quiet_flag = 0;     // DO NOT! print normal program output if set
                      //DM Note: I never use this flag, so it isn't programmed
  //flag -r//
  restart_flag = 0;
  //flag -m//
  movie_type = -1;  //default is a jeff v0
  write_movie_type = -1;  //default is a jeff v0

  out_file[0] = 0; //if you don't give a file name, code will know from 
                   //if(!out_file[0])
//---------------------------------end globals

  return;
}


void check_globals_for_errors(){

  //error checking// -- this is a subroutine waiting to happen
  if(!xmax && movie_type == -1){
     printf("Cannot process this movie type without size!, xmax=%lf",xmax);
     exit (-1);
    }


  if(write_movie_type == -1) restart_flag = 1;//make a delplot file

   if(movie_type == 1){
     printf("new ascii converter!\n");
     //     exit(1);
   }

  /*------------------------------------------------------------*/
   if( ! *out_file ) strcpy(out_file,"frame_"); //seems silly now

  if(verbose_flag){
    printf("Processing movie:\n\t %s \n",filename);
    printf("Writing to:\n\t %s \n",out_file);
    printf("System size: \n\t %f X %f \n",xmax,ymax);
  }

  //subroutine adds pinning force to particle-particle force
  //and finds the first minima... NOT LINEAR!!!
  
  float temp = estimate_neighbor_distance(pinning_force);
  /*
  if (temp < 1.0){
      bondLim = pinning_force / 2.0 * temp;
      }
  */

  if(!quiet_flag) printf("VERY rough estimate hard wall repulsion limit = %f\n",temp);

  return;
}
