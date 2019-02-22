#include "Vortex_2D.h"

/*TO DO
 *add PBC
 */

//if T/2 ~ \xsi, just bin by T/2... 
//or
//if particle arrangement is 1D, just bin by T/2

int density2D(int N, struct vortex *vortex_list, float *pin_location){

  int i,j;
  float system_density = ( (float) N ) / (xmax * ymax) ;
  int pin_bin[number_pins];
  int pin_density[number_pins];
  int number_NOT_pinned = 0;
  int number_pinned = 0;

  for(i=0;i<number_pins;i++) pin_bin[i] = 0;

  if(number_pins){
    /* We know the pin locations from main.c*/                          

    //populate the bins, where each pin gets a bin... 
    //this is sufficient for 1D and quasi 1D systems
    for(i=0;i<N;i++){
      //printf("here? %d %d", i, vortex_list[i].pin_number);
        pin_bin[vortex_list[i].pin_number]++;
    }
  }//end if(number_pins)
  else{
    if (!quiet_flag) printf("system density = %f\n",system_density);
  }

  //fine-grained density where \kappa = 1.0
  int bins_per_period = (int) ( (pinning_period) / xsi );

  //want an odd number of periods/pin s.t.
  //a pin may center on x_Np
  if (bins_per_period < 4)  //if bins = 0,1,2,3 too small to split up
    bins_per_period = 1;
  else if (bins_per_period % 2)  //if even, make odd
    bins_per_period += 1;
  
  int middle_bin = (int) bins_per_period/2;  
  //bins = 0,1,...bins_per_period-1
  //therefore 0th bin accounted for 
  //i.e. if n_bins==1, everything goes in bin[0]
  //-----if n_bins==5, middle bin is int(5/2) = 2, of 0,1,2,3,4, middle!!!

  if(!quiet_flag) printf("middle bin = %d\n", middle_bin);
  if(!quiet_flag) printf("distance between pins = %f",pinning_period);

  //open file to print to... either fine grained, or just by pin
  if((local_density_file = fopen("LocalDensity.txt","w"))==NULL) {
      printf("Cannot open file.\n");
      exit(1);
    }
  //print top line of file, either way
  fprintf(local_density_file,"#x_pos \t y_pos \t num_vor \t n_v \n");

  //if there are no pins, just print the homogeneous density
  if(!number_pins){
    fprintf(local_density_file,"%f %d %f\n",0.0, N, system_density );
  }
  else if(bins_per_period > 1){

    //warning, scoped within if
    //make a local density histogram based on xsi lengthscale
    int histogram_fine_grain[number_pins][bins_per_period];

    //zero the elements
    for(i=0;i<number_pins;i++){ 
      for(j=0;j<bins_per_period;j++){
	histogram_fine_grain[i][j] = 0;   
      }
    }  
  
    for(i=0;i<N;i++){
      //assign to bin
      float xloc = vortex_list[i].pin_dist;

      //fg_bin will be positive/negative about zero
      //xsi doesn't work very well for clusters
      int fg_bin = middle_bin + (int) xloc/xsi;  

      if(abs(fg_bin) < (bins_per_period) ){
	histogram_fine_grain[vortex_list[i].pin_number][fg_bin]++;
	number_pinned++;
	//printf("I am vortex: %d, my location is (%f,%f)",i,vortex_list[i].x[0],vortex_list[i].x[1]);
	//printf("x_loc = %f",xloc);
	//printf(" #%d, %d\n",vortex_list[i].pin_number,fg_bin);
      }
      //edge effect i think... probably should work on this
      else{
  	    number_NOT_pinned++;
	    //printf("\nvortex %d is not pinned, total not pinned = %d \n"
	    //,i,number_NOT_pinned);
	    //printf("\n x_loc = %f",xloc);
	    //exit(-1);
      }
      //printf("\nMax Pins = %d, bins/period = %d, 
      //I belong in bin %d \n", number_pins, bins_per_period, fg_bin);

	  
    }
    //print the fine grained histogram to LocalDensity.txt
    for(i=0;i<number_pins;i++){       
      for(j=0;j<bins_per_period;j++){
	float fg_bin_loc = pin_location[i] + xsi*(j-middle_bin);

	if (fg_bin_loc > xmax) fg_bin_loc -= xmax;
	if (fg_bin_loc < 0) fg_bin_loc += xmax;

	//add x-y values to Local Density
	fprintf(local_density_file,"%f \t 0.0000 %d \t %f\n",fg_bin_loc,
	       histogram_fine_grain[i][j],
	       (float) histogram_fine_grain[i][j]/(ymax*xsi));
	//fprintf(local_density_file,"%d %d \n",i,j);
      }
    }
    fprintf(local_density_file,"#number not pinned = %d \n", number_NOT_pinned);
    fprintf(local_density_file,"#density not pinned = %f \n", (float) number_NOT_pinned/(ymax * xsi) );
    fprintf(local_density_file,"#number pinned = %d \n", number_pinned);
    fprintf(local_density_file,"#density pinned = %f \n", (float) number_pinned/(ymax * pinning_period/2) );

  }//end if(bins_per_period > 1)
  else{
    //make a local density file?
      for(i=0;i<number_pins;i++){ 

	fprintf(local_density_file,"%f \t 0.0000 \t %d \t %f\n",\
               pin_location[i], \
	       pin_bin[i],pin_bin[i]/(ymax * pinning_period/2) );

      }//end for(i=0... loop
  }//end if (there is only one bin) else

  fclose(local_density_file);

  //--------------------------------------------------------------------//
  //--------------------------- Y histogram ----------------------------//
  //--------------------------------------------------------------------//

  //TODO pin the y bins by xsi-ish
  //open file to print to... either fine grained, or just by pin
  

  //hardcoding y-binning length
  int number_y_bins = number_pins;                            
  if(number_y_bins < 5) number_y_bins = 20;
  float y_bin_width = xsi*ymax/(float)number_y_bins;   
  float x_bin_width = xsi*ymax/(float)number_y_bins;   

  //printf("y_width=%f, xsi=%f, ymax=%f", y_bin_width, xsi, ymax);

  int y_histogram[number_y_bins]; //[number_x_bins]--- 2D
  int x_histogram[number_y_bins]; //[number_x_bins]--- 2D
  
  //try a stamp
  float xy_hist[number_y_bins][number_y_bins];  //start coarse

    for(i=0;i<number_y_bins;i++){ 
      y_histogram[i] = 0;
      x_histogram[i] = 0;
      for(j=0;j<number_y_bins;j++){
	xy_hist[i][j] = 0.0;
      }
    }  

    for(i=0;i<N;i++){
      int ybin = (int) vortex_list[i].x[1] / y_bin_width;
      int xbin = (int) vortex_list[i].x[0] / x_bin_width;
      y_histogram[ybin]++;
      x_histogram[xbin]++;

      //update the xy grid
      //stamp(xy_hist,i,j,number_y_bins);
      xy_hist[xbin][ybin]++;

      // printf("i=%d, ybin=%d, \n",i,ybin);
    }

  //Reuse pointer... why not?
  if((local_density_file = fopen("XY_LocalDensity.txt","w"))==NULL) {
      printf("Cannot open file.\n");
      exit(1);
    }
      fprintf(local_density_file,"\n");
      fprintf(local_density_file,"\n#x\t    y \t    N \t    N/N_v\n");

    for(i=0;i<number_y_bins;i++){ 
      float y_pos = ((float) i + 1.0)/number_y_bins * ymax;
      //add both x- and y- to LocalDensity.txt
      fprintf(local_density_file,"%10.3f %10.3f %6d %6.3f\n",xmax, y_pos, y_histogram[i], y_histogram[i] / (float) N);
    }

      fprintf(local_density_file,"\n");
      fprintf(local_density_file,"\n#x\t    y \t    N \t    N/N_v\n");

    for(i=0;i<number_y_bins;i++){ 
      float x_pos = ((float) i+ 1.0)/number_y_bins * xmax;
      //add both x- and y- to LocalDensity.txt
      fprintf(local_density_file,"%10.3f %10.3f %6d %6.3f\n",x_pos, ymax, x_histogram[i], x_histogram[i] / (float) N);
    }

      fprintf(local_density_file,"\n");
      fprintf(local_density_file,"\n#x\t    y \t    N \t    N/N_v\n");

    for(i=0;i<number_y_bins;i++){ 
      float x_pos = ((float) i+ 1.0)/number_y_bins * xmax;

      for(j=0;j<number_y_bins;j++){
	float y_pos = ((float) j+ 1.0)/number_y_bins * ymax;

	fprintf(local_density_file,"%10.3f %10.3f %6d %6.3f\n", x_pos, y_pos, xy_hist[i][j], xy_hist[i][j] / (float) N);
      }
    }

  
  fclose(local_density_file);

  return 0;
}

float avg_dist_neighbor(int N, struct vortex *vortex_list, float *std_dev, float *d_neighbor_max){
  int i;
  float avg_distance = 0.0;
  float max_d_local = 0.0;
  float sum_deviation = 0.0;
  float avg_x_dist = 0.0;
  float avg_y_dist = 0.0; 
  float x_dev = 0.0;
  float y_dev = 0.0;
  float std_dev_x;
  float std_dev_y;

  for(i=0;i<N;i++){
    avg_distance += vortex_list[i].d_neighbor;  //distance, NOT vector
    avg_x_dist += vortex_list[i].x_min;  //x-component of dmin
    avg_y_dist += vortex_list[i].y_min;  //y-component

    if( vortex_list[i].d_neighbor > max_d_local){
      max_d_local = vortex_list[i].d_neighbor;
    }
    //printf("length=%f,v_l.d=%f \n ", avg_distance, vortex_list[i].d_neighbor);
    //printf("v_l.d=%f \n ", vortex_list[i].d_neighbor);
  }

  //after summing all distances, divide by total particle number
  avg_distance /= ((float) N);
  avg_x_dist /= ((float) N);
  avg_y_dist /= ((float) N);

  *d_neighbor_max = max_d_local;

  //calculate the standard deviation for each quantity
  //(1) sum the squares
    for(i=0; i<N;i++)
      sum_deviation += (vortex_list[i].d_neighbor - avg_distance)*(vortex_list[i].d_neighbor - avg_distance);
      x_dev += (vortex_list[i].x_min - avg_x_dist)*(vortex_list[i].x_min - avg_x_dist);
      y_dev += (vortex_list[i].y_min - avg_y_dist)*(vortex_list[i].y_min - avg_y_dist);
    //printf("\nsun_deviation=%f\n",sum_deviation);

    *std_dev = sqrt(sum_deviation/(float)N);
    std_dev_x = sqrt(x_dev/(float)N);
    std_dev_y = sqrt(y_dev/(float)N);

    if(quiet_flag) printf("%f\t%f\t%f\t",avg_distance, *std_dev, max_d_local); 
   //these would go well in a global info file
   if(!quiet_flag) printf("\n average distance to neighbor = %f +/- %f\n",avg_distance,*std_dev);
   if(!quiet_flag) printf("\n MAX distance to neighbor = %f\n",max_d_local);

  return avg_distance;
}
//--------------------------------------------------------------------------//
float estimate_neighbor_distance(float fp){

  float r, force_r, min_r, max_r;
  int i;
  min_r = 0.05;
  max_r = 3.0;
  float dr = 0.01;
  float r0;
  force_r = 1000000;  //big number

  while (r < max_r){
    force_r = -fp + 1/(r*r) - 2*exp(-r);
    if( force_r < 0.01 ){  //just grab it as soon as it is close to zero
      r0 = r;
      //if(quiet_flag) printf("%10.6f",r0);  //NOT quiet_flagged!!!
      break;
    }
    r += dr;
  }
  if(!r0) printf("problem!\n");

  return r0;  //arbitrary... double it and add a fudge factor
}
//----------------------------------------------------------------------------//
//take the average ABSOLUTE distance from a pinning center
//and the standard deviation
//examine as a function of N_pins, N_vortices, and Force_pins
int check_1D(int N, struct vortex *vor){

  float avg_x_distance = 0.0;
  float sum_deviation = 0.0;
  float std_dev = 0.0;
  int system_dimension = 2;
  int i;
  float x;
  for(i=0;i<N;i++){
    //x = fabs(vor[i].pin_dist);
    //printf("x=%f\n",x);
    avg_x_distance += fabs(vor[i].pin_dist);
  }  

  avg_x_distance /= ((float) N);
    
  for(i=0; i<N;i++)
    sum_deviation += (fabs(vor[i].pin_dist) - avg_x_distance)*(fabs(vor[i].pin_dist) - avg_x_distance);

  std_dev = sqrt(sum_deviation/(float)N);
  
  if ( fabs(avg_x_distance) < 2*std_dev && fabs(avg_x_distance)<0.1 ){
    system_dimension = 1;
    }
  if(!quiet_flag) printf("average distance to pinning center= %f +/- %f\n",avg_x_distance,std_dev);
  if(quiet_flag) printf("%f\t%f\t", avg_x_distance,std_dev);

  return system_dimension;
}

//----------------------------------------------------------------------
//omg, omg, omg
// i dynamically allocated an array of floats within this
// subroutine!!! this is powerful and can potentially
// help me recode the "vortex.bonds_to" static declaration crap

void find_nearest_y_neighbors(int N, struct vortex *vor, float *pin_location){
  
  int i,j;

  struct pin_site pin[number_pins];  //malloc dynamically!
  int size = 2*N/number_pins;       //roughly twice number of vortices/pin

  for(i=0;i<number_pins;i++){
    pin[i].pin_id = i;       //see, unnecessary
    pin[i].x_pin = pin_location[i]; //1-to-1 mapping
    pin[i].y_pos = (float*) malloc(size * sizeof(float) );
    pin[i].num_vor = 0;
  }

  for(i=0;i<N;i++){
    int array_index = pin[vor[i].pin_number].num_vor;
    if(array_index < size){
      pin[vor[i].pin_number].y_pos[array_index] = vor[i].x[1];
    }
    else{
      float *tmp = realloc(pin[vor[i].pin_number].y_pos, 2*size*sizeof(float));
      if (tmp == NULL){
        //Error
      }
      else{
	pin[vor[i].pin_number].y_pos = tmp;
      }
    }//end realloc-ing else

    pin[vor[i].pin_number].num_vor++;
  }  

  //okay, now gotta sort the damn thing, then print
  //qsort!
  float dist1, dist2;

  for(i=0;i<number_pins;i++){

    qsort(pin[i].y_pos, pin[i].num_vor, sizeof (float), compare_floats);

    printf("#Pin %d ", pin[i].pin_id);
    printf("has %d vortices, neighbor distances:\n", pin[i].num_vor);

    //start one up from zero
    //end one down from max
    //pbc would account for the 0th and Nth vortices
    for(j=1;j<pin[i].num_vor-1;j++){

      dist1 = pin[i].y_pos[j+1]-pin[i].y_pos[j]; //just lengths
      dist2 = pin[i].y_pos[j]-pin[i].y_pos[j-1];
      //okay the difference between distance 1 and 2
      //is actually interesting, 
      //second order differences will reveal grains...
      printf("%f %f %f \n",pin[i].x_pin,pin[i].y_pos[j],fabs(dist1-dist2) );
    }

    printf("\n");
  }

  //if you extend this past the 1D systems, you can use this to detect
  //bunches in the y-direction in a single pin

  return;
}


int compare_floats (const void *a, const void *b){
       const float *da = (const float *) a;
       const float *db = (const float *) b;
     
       return (*da > *db) - (*da < *db);
     }

//-------------------------------------------------------------------
void pin_centered_y_neighbors(int N, struct vortex *vor, float *pin_location){
  
  int i,j;

  struct pin_site pin[number_pins];  //malloc dynamically!
  int size = 2*N/number_pins;       //roughly twice number of vortices/pin

  for(i=0;i<number_pins;i++){
    pin[i].pin_id = i;       //see, unnecessary
    pin[i].x_pin = pin_location[i]; //1-to-1 mapping
    pin[i].y_pos = (float*) malloc(size * sizeof(float) );
    pin[i].num_vor = 0;
  }

  for(i=0;i<N;i++){
    if( fabs(vor[i].pin_dist) < 0.5 ){
      int array_index = pin[vor[i].pin_number].num_vor;
      if(array_index < size){
	pin[vor[i].pin_number].y_pos[array_index] = vor[i].x[1];
      }
      else{
	float *tmp = realloc(pin[vor[i].pin_number].y_pos, 2*size*sizeof(float));
	if (tmp == NULL){
	  //Error
	}
	else{
	  pin[vor[i].pin_number].y_pos = tmp;
	}
      }//end realloc-ing else

      pin[vor[i].pin_number].num_vor++;
    }//end if close enough
  }  

  //okay, now gotta sort the damn thing, then print
  //qsort!
  float dist1, dist2;

  for(i=0;i<number_pins;i++){

    qsort(pin[i].y_pos, pin[i].num_vor, sizeof (float), compare_floats);

    printf("#Pin %d ", pin[i].pin_id);
    printf("has %d central vortices, neighbor distances:\n", pin[i].num_vor);

    //start one up from zero
    //end one down from max
    //pbc would account for the 0th and Nth vortices
    for(j=1;j<pin[i].num_vor-1;j++){

      dist1 = pin[i].y_pos[j+1]-pin[i].y_pos[j]; //just lengths
      dist2 = pin[i].y_pos[j]-pin[i].y_pos[j-1];
      //okay the difference between distance 1 and 2
      //is actually interesting, 
      //second order differences will reveal grains...
      printf("%f %f %f \n",pin[i].x_pin,pin[i].y_pos[j],fabs(dist1-dist2) );
    }

    printf("\n");
  }

  //if you extend this past the 1D systems, you can use this to detect
  //bunches in the y-direction in a single pin

  return;
}

