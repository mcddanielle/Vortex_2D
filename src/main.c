#include "Vortex_2D.h"

int main(int argc, char **argv){

  /*-------------------configure globals ---------------------- */
  assign_global_defaults();
  getoptions( argc, argv );   //command line arg processor
  check_globals_for_errors(); 

  /*------------------set pin lattice ------------------------- */
  float pin_location[number_pins];  //declare here for variable scope
  if(number_pins){
    find_pinning_minima(&pin_location);  
    //printf("pin loc: %f\n",pin_location[0]);
  }
  else{
    //number_pins = 0 and it is a NULL array
    //set the null pin in the center of the system
    pin_location[0] = xmax/2;
    //look on either side to the edges
    pinning_period = xmax;  //globals
    }//else may prove completely unnecessary

  /*-----------------------Print information--------------------*/
  //won't want this in a true production run
  if (quiet_flag) printf("#1 F_p    2 N_p    3 r0     4 gr_min   5 gr_norm_factor   (6 <d_avg> +/- 7 sigma_{d_avg}) \t 8 max_neigh_dist \t (9 <d_pin> +/- 10 sigma_{d_pin})  (11 rho_digital  12 percent_edge_bin 13 percent_center 14 percent_edge_particle 15 first_max  16 last_max ) 17 N_v  18 bondLim \n ");//quiet flag
  if (quiet_flag) printf("%f\t%d\t",pinning_force,number_pins);


  /*------------------- read movie SUBROUTINE-------------------*/

  struct vor_bond *vhead = NULL;  //pass ptr to this from load_coordinates
  int num_vor, i;
  //struct vortex *vortex_list;

  vhead = load_movie_ptr(&num_vor, pin_location);//, &vortex_list);
  //vortex_list comes back blank... you don't know what you're doing
  /*
  printf("\n vortex list: size = %d\n", sizeof(vortex_list) );
  printf("\n vortex size = %d\n", sizeof(struct vortex) );
  printf("vor list x[0]=%f\n",vortex_list[0].x[1]);
  */

  if (quiet_flag) printf("%d\t%f\t",num_vor,bondLim);

  //ideally pass back the atom pointer here and then run angular_correlation
  //the problem with vor_list[N] is the malloc... 

  //save head
  struct vor_bond *v1 = vhead;

  //calculate the g(r) data
  int bondcount = num_vor * (num_vor-1 ) / 2;  //unnecesary?
  if(!quiet_flag) printf("num vortices = %d, total neighbors = %d", num_vor, bondcount);
  //pair_correlation( num_vor, vhead );

  //reset head for q6 calc
  vhead = v1;


   //----------if you haven't printed yet, do the final frame---------//
  /*
   if(frame_printed_flag == 2){
          printf("printing final frame at current=%f\n",current);

          if(write_movie_type == 1){
              print_ascii_frame(num_vor,temp,current,x,y);
	  }
	  else print_frame(num_vor,temp,current,x,y);
   }//--------------------end print if you haven't yet----------------//

  
  */
if (quiet_flag) printf("\n");//quiet flag
  return 0;

}//main function ends here!

//-------------------------------------------------------------------------//
void find_pinning_minima(float *pin_location){

  pinning_period = xmax / (float) number_pins;  //globals!
  if(!quiet_flag) printf("pinning period T = %f\n",pinning_period);
  int i;
  //  float pin_location[number_pins];

  //minimum attractive well located at \theta = 3/2*PI
  //since F_p \sim cos() and U_p \sim sin()
  //i.e. sin(2 PI x / T) is minimized at x = 3/4
  //each period has one minima 3/4T from it's left edge
  for(i=0;i<number_pins;i++){
    pin_location[i] = (i + 0.75) * pinning_period;
    //printf("location of attractive wells: %f\n",pin_location[i]);
  }
  return;
}

//-----------------------------------------------------------------// 
//subroutine to find the pin location simply based on x location
//does take periodicity into account
//if x<T/4, it should wrap back to the final pin
//which is located at xmax - T/4
//TODO: should pass back pin_dist by reference while I'm at it
int find_pin(float x, float *pin_location){

  int pin_id = 0;
  //printf("x=%f\n",x);
  //case where x is nearly zero and pbc wraps to max pin number
  if (x < (float)pinning_period/4.0){
    pin_id = number_pins - 1;
  }
  else{
    //integer division: floor: 0 to number_pins-1
    pin_id = (int) x/pinning_period;  
    if(DEBUG) printf("\n pin_id = %d\n",pin_id);

    //make sure that is the closet pin via period/2
    if ( x-pin_location[pin_id] > pinning_period/2 ){
      pin_id++;
    }
    if ( x-pin_location[pin_id] < -pinning_period/2 ){
      pin_id--;
    }       

    //check the periodic condition and wrap back to zero if necessary
    if (!abs(pin_id - number_pins)) pin_id = 0;
    if (pin_id == -1) pin_id = number_pins - 1;

    if(DEBUG) printf("update to pin_id = %d\n",pin_id);
  }

  
  return pin_id;
}
