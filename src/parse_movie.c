#include "Vortex_2D.h"

//--------------------------parse movie-----------------------------------//
/*Each frame has a header containing the number of vortices, the time step,
 *int i will correspond to the num_vor OR integers 0 to (num_vor-1)
 */

  /* We know the pin locations from main.c                            
   *for(i=0;i<number_pins;i++) 
   *  printf("location of attractive wells: %f\n",pin_location[i]);
   */

struct vor_bond *load_movie_ptr(int *num_vor, float *pin_location){//, struct vortex *vor_list){
//, float *xp){  //works
//, struct vortex **vortex_list){ //doesn't work

  /*-----------------flags to control printing------------------*/

  int frame_printed_flag = 0;  //refers to single frame, not used 
  
 //---------------------ascii... not coded-----------------------//
  if(movie_type == 1){ 
    //read_single_ascii();  
    printf("Not yet coded, use a binary smtest.\n");
    return 0;
  }

  /*-----------------------declare movie variables-------------------*/
  int N, time, i, j;      
  //use a while loop to read movie, but increment i 
  i = 0;           

  /*--------------------Open the movie for Reading----------------------*/

  if ((open_movie = fopen(filename, "rb")) == NULL){
         printf("Can't open %s.\n",filename);
         exit(1);
  }

  //read header to get number of vortices
  fread(&N,sizeof(int),1,open_movie);
  fread(&time,sizeof(int),1,open_movie);

  //return num_vor to main.c
  *num_vor = N;
  
  float x[*num_vor], y[*num_vor];

  //-----------------------------------------------------------------//
  //---------------------the main loop starts here!------------------//
  //-----------------------------------------------------------------//
   

  do{//--------------------------------DO{ --- } while(not end of movie)   
    if(ferror(open_movie)){
      printf("Error reading movie\n");
      exit(-1);  
    }

    /*----Read first entry as i, check whether it matches num_vor.---*/
    fread(&i,sizeof(int),1,open_movie);  //on first loop i = id = 0
    /*---------  it will either be vortex id or num_vor  ------------*/

    //---------------------------------------------------------------//
    //------------------------- NEW FRAME!!! ------------------------//
    //---------------------------------------------------------------//

     if(i==*num_vor){                      
       /*If this condition is met, read time*/                             
       fread(&time,sizeof(int),1,open_movie);
      
      /*THIS WORKS BECAUSE each line starts with an int (id)
       *that can always be checked against num_vor,
       *since id runs from 0 to num_vor-1, this condition is special
       *
       *DM: 11/1/2013 can you tell I went through
       *several existential crisis related to the comment above?
       */ 
     }
     //-----------------------  NOT NEW ---------------------------//
     //read in individual vortex data
     else{               
       fread(&x[i],sizeof(float),1,open_movie);
       fread(&y[i],sizeof(float),1,open_movie);
       if(DEBUG) printf("%10d %10f %10f\n", i, x[i], y[i] );
     }//end else----- where if = if(i==num_vor)

     //--------------------- End of Frame!!! -------------------------//

  }while (!feof(open_movie));  // end do until EOF

  fclose (open_movie);

   /*keep assignments separate from the load,
    *potentially handle an entire movie rather than a single frame
    *make a linked list of num_vor length
     struct vortex tempv = malloc(2*sizeof(struct vortex));  //1-element array
    */

   //-----------------------------------------------------------------//
   //-----------------------Assign Structs----------------------------//
   //-----------------------------------------------------------------//
   struct vortex vortex_list[N];  //malloc dynamically!

   //new 11/22/2013
   //struct vortex *vortex_list, *vl_head;

   //temporary distances between two particles... written over throughout
   float length, vector_IJ[3];   
   int pin_id;

   float avg_d, std_dev, max_neighbor_distance;    //<d_avg> measure

   //get every neighbor interaction.  
   //i.e. N(N-1)/2
   struct vor_bond *bond, *vb_head;

   //vb_head = head of linked list of bonds
   vb_head = NULL; 

   //new 11/22/2013
   //vl_head = NULL;

   for (i=0;i<(*num_vor);i++){

     //zero these things
     vortex_list[i].bondcount = 0;         //number of neighbors
     vortex_list[i].pin_number = 0;        //closest pin number
     vortex_list[i].pin_dist = 0.0;        //closest pin distance in x

     //arbitrarily big, xmax = 36.5
     vortex_list[i].d_neighbor = xmax;     //closest neighbor distance
     vortex_list[i].x_min = xmax;          //closest x distance (1D check)
     vortex_list[i].y_min = xmax;          //closest x distance (1D check)
     vortex_list[i].y_min2 = xmax;          //closest x distance (1D check)

     //bond_counter = 0;
     vortex_list[i].id = i;                //from movie 0 to num_vor
     vortex_list[i].x[0] =  x[i];        
     vortex_list[i].x[1] =  y[i];
     vortex_list[i].x[2] =  0.0;           //z to work with 3D routines

     //an experiment that i'm struggling with
     //1-element array
     vortex_list[i].bt_new = (struct vortex*) malloc(2*sizeof(struct vortex));  

     /*-------------------  Find Closet Pin ----------------------*/
     //TODO include &pin_dist in find_pin() subroutine
     if (number_pins > 1){
       pin_id = find_pin(vortex_list[i].x[0], pin_location);
     }
     else{
       pin_id = 0;
     }
     //printf("pin_id=%d",pin_id);

     vortex_list[i].pin_number = pin_id;
     vortex_list[i].pin_dist = vortex_list[i].x[0] - pin_location[pin_id];
     //periodic wrap here
     if (x[i] < (float)pinning_period/4.0){
       vortex_list[i].pin_dist += xmax;
     }

     if(DEBUG){//------------------------------------------------------//
       if(fabs(vortex_list[i].pin_dist) > pinning_period){
	 printf("pin=%d ",vortex_list[i].pin_number);
	 printf("pin distance =%f,pin location= %f \n vortex location = %f,%f\n", \
		vortex_list[i].pin_dist, pin_location[pin_id],		\
		vortex_list[i].x[0],vortex_list[i].x[1]);
       }
     }//end if(debug)//-----------------------------------------------//

     /*-----------------------------------------------------------*/

   for(j=0;j<i;j++){
     length = 0.0; //

     /*now here we must be very careful about PBC using CIJOR methods */
     length = periodic_bondlength( vortex_list[j].x , vortex_list[i].x );

     /*calculate the vector quantity of length*/
     //should write the periodic part into vecDiff, but this is
     //a macro, NOT a subroutine.
     vecDiff(vector_IJ, vortex_list[i].x, vortex_list[j].x);

     if(vector_IJ[0] > xmax/2) vector_IJ[0] -= xmax;
     if(vector_IJ[1] > ymax/2) vector_IJ[1] -= ymax;
     if(vector_IJ[0] < -xmax/2) vector_IJ[0] += xmax;
     if(vector_IJ[1] < -ymax/2) vector_IJ[1] += ymax;
	    

     /*if the neighbor distance is shorter than previously encountered, reset*/
     /*simultaneously set xmin and ymin as vector components*/
     if(length < vortex_list[i].d_neighbor){
	 vortex_list[i].d_neighbor = length;
	 vortex_list[i].x_min = fabs(vector_IJ[0]);
	 vortex_list[i].y_min = fabs(vector_IJ[1]);
     }
     if(length < vortex_list[j].d_neighbor){
	 vortex_list[j].d_neighbor = length;
	 vortex_list[j].x_min = fabs(vector_IJ[0]);
	 vortex_list[j].y_min = fabs(vector_IJ[1]);
     }       

     //allocate mem for next item, part of the 'bond' LINKED list
     bond = (struct vor_bond *)malloc(sizeof(struct vor_bond)); 
     //ptrs to bonding atoms... 
     bond->vor1 = &vortex_list[i];       
     bond->vor2 = &vortex_list[j];

     //store periodic length in struct
     bond->length = length;        

     //store periodic distance vector in struct
     bond->vec_x[0] = vector_IJ[0];
     bond->vec_x[1] = vector_IJ[1];
     bond->vec_x[2] = vector_IJ[2];

     //0 = same pin, 1 nn pins, 2 nnn pins, etc
     bond->bond_type = abs(vortex_list[i].pin_number-vortex_list[j].pin_number);

     /*******************************************************************/
     /*Commenting out 1D pin calculation... it is ugly*/
     //find the closest x distance as a test for 1D systems
     //i.e. everybody is in a line in a pin
     //okay, OF COURSE this value is arbitrarily close
     //to zero as written... how can you apply a smart cut to this???
     /*
     if(  bond->bond_type == 1 ){
	 if( fabs(vector_IJ[0]) < fabs(vortex_list[j].x_min) )
	   vortex_list[j].x_min = vector_IJ[0];
	 if( fabs(vector_IJ[0]) < fabs(vortex_list[i].x_min) )
	   vortex_list[i].x_min = vector_IJ[0];
     }
     */
     //--------------------------------------------------------------
     //NOT GREAT BECAUSE you want to detect the +/- for
     //each vortex, even if there is a large void to the next one

     //these aren't consistently getting set... is it else if?
     //even if they are... not good enough... need to sort
     //probably better off with a histogram than an average
     //with this measure since there are plenty of y_min~=y_min2 cases
     /*
     if( bond->bond_type == 0){
       if( fabs(vector_IJ[1]) < fabs(vortex_list[j].y_min) )
	   vortex_list[j].y_min = vector_IJ[1];
       else if( fabs(vector_IJ[1]) < fabs(vortex_list[j].y_min2))
		vortex_list[j].y_min2 = vector_IJ[1];
       if( fabs(vector_IJ[1]) < fabs(vortex_list[i].y_min))
	   vortex_list[i].y_min = vector_IJ[1];
       else if( fabs(vector_IJ[1]) < fabs(vortex_list[i].y_min2))
		vortex_list[j].y_min2 = vector_IJ[1];
     }
     */
     /*******************************************************************/

     /*Lower cutoff, 0.1, error prevention, corresponds with wstripe code
      *Upper cutoff, bondLim, is arbitrary... working on that
      */
     if(length < bondLim && length > 0.1){ //periodic issues
         /*---------------------"bonds to"----------------------
	  *internal counter atm[i].bondcount 
	  *runs from 0 to Catom-1
	  *upper bound is C=20... could malloc it instead
	  *atom_i "bt" atom_j
	  */

         vortex_list[i].bt[vortex_list[i].bondcount] = &vortex_list[j]; 
         vortex_list[j].bt[vortex_list[j].bondcount] = &vortex_list[i]; 

	 /*SEE ATTEMPT TO MALLOC BELOW*/
	 
	 //increment bondcount AFTER update "bonds to" array
	 vortex_list[i].bondcount += 1;
	 vortex_list[j].bondcount += 1;	
	 // could always realloc to shorten array
	 // AFTER all ij combinations are checked

       }/*-------------------end "bonds to"------------------*/
     
     bond->next = vb_head;          //set pointer to previous item
     vb_head = bond;                //increment list
     }//close -------------------------j loop
   }//close --------------------------i loop 
   //all ij combinations have been checked... 

   //------------------------------------------------------------------//
   //---------------------measure some stuff---------------------------//
   //------------------------------------------------------------------//

   int layer = layer2D(N,vortex_list);

   //run bond pair correlation, return first shell closure
   //returns nothing, but resets global variable bondLim
   pair_correlation( *num_vor, vb_head );

   //What is the average nearest-neighbor distance?
   avg_d = avg_dist_neighbor(N,vortex_list,&std_dev,&max_neighbor_distance);

   //Are all of the vortices within some limit of pinning centers?
   //if AVG distance to pin is SMALL... keep an eye on this 
   int dimension = check_1D(N,vortex_list);

   if(DEBUG){//--------------------------------------------------------
     if(dimension == 1 && number_pins){
       //make it a void temporarily
       find_nearest_y_neighbors(N,vortex_list,pin_location);     
     }
     else if(number_pins){
       pin_centered_y_neighbors(N,vortex_list,pin_location); 
     }
   }//end if debug-----------------------------------------------------

   //so far dim==1 says all vortices are tightly pinned
   //what about pin interactions?

   if(!quiet_flag) printf("This is effectively a %d-D system\n",dimension);
 
   /*NOT sure what this was meant to do...
   if(dimension == 1 && pinning_period < avg_d){
     //do nothing
   }
   else{
     dimension = 2;
   }
   */

   //-----------------------------------------------------------------//


   //Can you use the above (and fp) to set length considerations?

   //Local Density
   //   if(number_pins == 1){
   //DM seg faulting, not priority right now
   //commmeting out
   //density2D(N,vortex_list,pin_location);
     //}
   stamp_vortices(N, vortex_list);

  //BOO parameters... cut off neighbors at next-nearest for hex lattice
   angle2D_correlation(N, vortex_list, bondLim);

  /*to use the above subroutine, define nn by strength of attractive potential*/


  /*
  for(i=0;i<N;i++){
    int N_nb = vortex_list[i].bondcount;    //make a shorthand for loops
    for(j=0;j<N_nb;j++){
      printf("old: %d ", vortex_list[i].bt[j]->id);
      //printf("new: %d ", vortex_list[i].bt_new[j].id);

    }
  }
  */

   //printf("\n before memset vortex list = %d\n",sizeof(*vortex_list));
   //vor_list = (struct vortex*) malloc( N * sizeof(struct vortex) );

   //   printf("after memset size of vor_list = %d\n",sizeof(vor_list));
   //*vor_list = *vortex_list;
   //printf("size of vor_list = %d\n",sizeof(vor_list));

   return vb_head;  
}


//--------------------------------------------------------------------------
void print_frame( int num_vor, float temp, float current, float *x, float *y){

   float vel1=0.0;
   float vel2=0.0;
   int time=0;
   int i;

   //convert current to a string and make a file name-------------------
   char str_current[20];
   sprintf(str_current,"%f",current); 

   if(! out_file[0] ){
     strcpy(out_file,"frame_current_");
   }

   strcat(out_file,str_current);     //cat into "frame_current_0.1"
   //printf("out_file=%s\n",out_file);
   //------------------------------------------------file NAMED!

   if ((kmovie = fopen(out_file, "wb")) == NULL){
         printf("Can't open %s.\n",kmovie);
         exit(1);
   }

   fwrite(&num_vor,sizeof (int),1,kmovie);
   fwrite(&time,sizeof (int),1,kmovie); 
  
   //If the restart flag is zero, write a jeff movie
   //this is deprecated... should be if(movie_type != -1)

   if(!restart_flag){
     fwrite(&temp,sizeof (float),1,kmovie);
     fwrite(&current,sizeof (float),1,kmovie); 

     fwrite(&xmin,sizeof (float),1,kmovie);
     fwrite(&xmax,sizeof (float),1,kmovie); 
     fwrite(&xmin,sizeof (float),1,kmovie);  //NOT redundant xmin=ymin
     fwrite(&xmax,sizeof (float),1,kmovie);  //xmax=ymax
   }

   for(i=0;i<num_vor;i++){

    fwrite(&i,sizeof (int),1,kmovie);
    fwrite(&x[i],sizeof (float),1,kmovie);
    fwrite(&y[i],sizeof (float),1,kmovie);

   //If the restart flag is zero, write a jeff movie
    if(!restart_flag){
      fwrite(&vel1,sizeof (float),1,kmovie);
      fwrite(&vel2,sizeof (float),1,kmovie); 
    }

   }

    fclose (kmovie);

  return;
}

//--------------------------------------------------------------------------
void print_ascii_frame( int num_vor, float temp, float current, float *x, float *y){

   int i;

   char str_current[20];
   sprintf(str_current,"%f",current); //convert current to a string
   char out_file[100] = "ascii_";    //beginning of out_file name
   strcat(out_file,str_current);     //cat into "frame_current_0.1"

   printf("out_file=%s\n",out_file);

   if ((kmovie = fopen(out_file, "w")) == NULL){
         printf("Can't open %s.\n",kmovie);
         exit(1);
   }

   fprintf(kmovie,"#Number of Vortices: %d\n",num_vor);
   fprintf(kmovie,"#X Domain: %f,%f\n",xmin,xmax); 
   fprintf(kmovie,"#Y Domain: %f,%f\n",ymin,ymax); 

   for(i=0;i<num_vor;i++){
     fprintf(kmovie,"%d %f %f\n", i, x[i],y[i] ); 
   }

    fclose (kmovie);

  return;
}

/*---------------  malloc a linked list within the vortex_list--------------*/
/*ATTEMPT TO MALLOC*/
	 /*
	 vortex_list[i]->bt_new[vortex_list[i].bondcount] = vortex_list[j];
	 vortex_list[j]->bt_new[vortex_list[j].bondcount] = vortex_list[i];
	 
	 void *temp_i = realloc(vortex_list[i]->bt_new, (2*vortex_list[i].bondcount)*sizeof(struct vortex));
	 void *temp_j = realloc(vortex_list[j]->bt_new, (2*vortex_list[j].bondcount)*sizeof(struct vortex));
	 if ( temp_i != NULL ){ //realloc was successful
	   vortex_list[i]->bt_new = temp_i;
	 }
	 else{ //there was an error
	   free(vortex_list[i]->bt_new);
	   printf("Error allocating memory!\n");
	   return 1;
	 }

	 if ( temp_j != NULL ){ //realloc was successful
	   vortex_list[j]->bt_new = temp_j;
	 }
	 else{ //there was an error
	   free(vortex_list[j]->bt_new);
	   printf("Error allocating memory!\n");
	   return 1;
	 }
	 */
	 //vortex_list->next->next->next
	 //where the linked list 
	 //is a more efficient version of the array
	 //pointer that will require some
	 //recoding in angular_parameters()
