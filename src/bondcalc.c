#include "Vortex_2D.h"

/*
 *5/16/2014 
 *Implementing the Gaussian blip histogram requires a careful tuning of 
 *sigma and dr for a particular density system.  
 *TODO: re-implement the config/parameter file for this
 */

float bondlength(float *x1, float *x2){

  float l12 = sqrt( pow( (x1[0] - x2[0]) , 2 ) 
                   + pow( (x1[1] - x2[1]) , 2 ) 
                   + pow( (x1[2] - x2[2]) , 2 ) );
  return (l12);
} 

//END: double bondlength(float * x1,float * x2)

float periodic_bondlength(float *x1, float *x2){

  float distx = fabs(x1[0] - x2[0]);
  float disty = fabs(x1[1] - x2[1]);
  //printf("before distx: %f, fabs disty: %f", fabs(distx), fabs(disty));
  if (distx > xmax/2) distx -= xmax;
  if (disty > ymax/2) disty -= ymax;
  //ignore z coordinate here
  //printf("after distx: %f, fabs disty: %f ", fabs(distx), fabs(disty));

  float l12 = sqrt( distx*distx + disty*disty ); 
  //inefficient to use pow() AND to keep z = xn[2]

  return (l12);
} 

//VECTOR simple pythagoran theorem subroutine
float vec_length(float *x1){

  float L = sqrt( pow(x1[0] , 2) +
                   pow(x1[1] , 2) +
                   pow(x1[2] , 2));

  //printf("L = %f x1 = %f %f %f\n",L,x1[0],x1[1],x1[2]);
  return (L);
}
//-------------------------------------------end vec_length

/*
 *DM 5/16/2014 FINALLY FINALLY FINALLY UNDERSTAND NORMALIZATION,
 DON'T CHANGE THE PARAMETERS OR THE PREFACTOR OF G(R)
 */

//------------------------------------------------------
int pair_correlation( int N, struct vor_bond *b1 ){

  struct bond_data *bhead = b1;
  FILE *file_pair;
  FILE *get_params;

  /*and here lies the new coding task... figure out
   *what the attractive length scale should be
   *and a good guestimate of the parameters
   *to put into gr_paramter read
   ****** super bonus for using dan's library
   *also gotta get the period in there
   */

  //define and get g(r) parameters here:
  double sigma, RMAX, RMIN, dr;

      //DM 5/1/2014
      RMIN = 0.5;    //Hardwire!!!
      RMAX = 10.0;
      sigma = 0.001;  //DM 5/16/2014 KEEP ratio sigma/dr = 1/4-1/5 or smaller
      dr = 0.04;

      //Nv=400,Np=0 sigma = 0.0005, dr = 0.01

      //want the height of the Gaussian to equal 1
      //want Gaussian sufficiently narrow compared with 

      /*TODO: look at how other people handle this... including Dan!
  if((get_params = fopen("gr_parameters.txt","r"))==NULL){
      printf("Cannot open file, using default parameters\n");
      return 1;
      RMIN = 0.2;
      RMAX = 36.5/2;
      sigma = 0.001;
      //dr = 0.005;
      //seems to be an error when no gr_params is provided
  }
  else{
    //fscanf parameters here
    if( fscanf(get_params,"%lf%lf%lf",&sigma,&RMAX,&RMIN) == EOF){
      printf("improper parameter read!\n");
      exit(-1);
      }
          fclose(get_params);
    }//end else
      */

      //      dr = sigma*10.0;

    if (!quiet_flag) printf("\nparameters are hardcoded: %lf %lf %lf %lf\n", \
                                      sigma,RMAX,RMIN,dr);

    if(RMAX < RMIN){
      printf("(RMAX < RMIN) those parameters are wrong!\n");
      exit(-1);
    }

    //name the file... kinda ugly!
  char *sigma_str[6];
  sprintf(sigma_str,"%0.3f",sigma);
  
  char output_str[80];
  strcpy (output_str,"gr_sigma");
  strcat (output_str,sigma_str);
  strcat (output_str,".dat");

  if (!quiet_flag) printf("write g(r) to: %s\n",output_str);
  //end name the file...

  if((file_pair = fopen(output_str, "w"))==NULL){
      printf("Cannot open file.\n");
      exit(1);
  }

  //write a simple header to g(r) data file
  fprintf(file_pair, "#r      total \n");

  /*for each r, calculate a Gaussian blip that describes the correlation
   *between that r value and all the lengths to describe the cluster
   *technically, should start at r = 0.0, but this contribution is nil
   */

  int i, j;
  double r = RMIN;        //initialize r
  double g_r_temp = 0.0;  //g_r for each step
  double g_r_new;
  double norm_factor = 0.0;
  //prefactor

  //5/1/2014 removed 1/(sqrt(2pi)*sigma)
  double Cg_r = 36.5*36.5/(N*N*2*PI*dr) ; ///(2*PI*dr);

  if(0){
    Cg_r *= 36.5*36.5/N*N;
  }

  double new_prefactor = 36.5/(2*PI*dr * N * N );// missing r and \rho

  double g_r, maxg_r, rmax, norm_g_r, gr_min1, r_min1;

  int flag_r_max = 0;
  int flag_r_min = 0;

  double grX = 0.0;
  double grY = 0.0;
  double grX_temp = 0.0;
  double grY_temp = 0.0;
  
  maxg_r   = 0.0; //first peak 
  rmax     = 0.0; //bondlength of first peak g(r_max)
  norm_g_r = 0.0; //divide by this at the end
  gr_min1  = 0.0;
  r_min1   = 0.0;

  //each r value has a unique g(r) to calculate
  //these are printed, not stored
  while(r<RMAX){     //increment: r+=dr;

    //zero g_r for every new r value
    g_r = 0.0;  
    grX = 0.0;
    grY = 0.0;
    g_r_new = 0.0;

    //reset to head of bond list for each r 
    b1 = bhead;        

    //loop through every bond and count up the contributions
    
    while(b1){
      //printf("%d %d %f\n",b1->vor1->id,b1->vor2->id,b1->length);
      //no r^{-2} value,to speed up calculations
      g_r_temp = exp(-pow((r-(b1->length) ),2) / (2*sigma*sigma) );
      g_r += g_r_temp;                  

      //where dr for continuous is too small for the discrete measure
      if( fabs( r - b1->length ) <= dr ) g_r_new++;

      if(r>0.6){
	grX_temp = exp(-pow((r-fabs(b1->vec_x[0]) ),2) / (2*sigma*sigma) );
	grY_temp = exp(-pow((r-fabs(b1->vec_x[1]) ),2) / (2*sigma*sigma) );
      //add up all length contributions

	grX += grX_temp;
	grY += grY_temp;
      }
      //cute, the bond->type allows us to effectively
      //use a hash to assign the value
      //rather than if/case statements
      g_r += g_r_temp;
      
      b1 = b1->next;
     }  //end while(b1) 
     

    /*calculate the r^{-2} value before determining maxima
     *otherwise the algorithm for finding maxima tends
     *to miss the first peak (if peaks are close)
     */

    g_r *= Cg_r * pow(r,-1); /// r;           //2D
    grX *= Cg_r * pow(r,-2); /// r;
    grY *= Cg_r * pow(r,-2); /// r;
    g_r_new *= new_prefactor / r;

    norm_factor += g_r*dr;    //NOT zeroed for new r-values
    //check if g_r is a maxima at this r value

    //    
    fprintf(file_pair, "%2.4f  %2.4f  %2.4f \n", r, g_r, g_r_new);

    if( (g_r >= maxg_r) && !flag_r_max){
	 maxg_r = g_r;
	 rmax = r;
    }
    else{
            //only turn on the flag for nonzero g(r)
	    if(g_r)  flag_r_max = 1;
    }
      
    if(flag_r_max && !flag_r_min && !gr_min1 && g_r < 0.1){
	 gr_min1 = g_r;
	 r_min1  = r;
	 //printf("\nfirst shell:%f,%f\n",r,r_min1);
    }
    else{
      if(gr_min1){
	flag_r_min = 1;
      }
    }

    //end maxima check
     
    //-------------print the value to the file-------------//
       //       fprintf(file_pair, "%2.4f  %2.4f %2.4f %2.4f\n", r,	\
      //         g_r, grX, grY);

    //----------------end printing section-----------------//

    r += dr;
  }  //------------------------------------end r<RMAX while loop

  //printf("\nfirst shell closes at %f,%f???\n",r_min1,bondLim);	
  //redeclare global
  if(flag_r_min && r_min1 < bondLim) bondLim = r_min1;  //only reset if r_min1 != 0.00

      if(!quiet_flag){
	printf("\nfirst shell closes at %f\n",bondLim);
        printf("\n normalization factor %f\n", norm_factor);
      }
      if(quiet_flag) printf("%f\t%f\t%f\t ",rmax,maxg_r,norm_factor);  


      fprintf(file_pair,"\n# set arrow from %f,0 to %f,%f nohead\n",
                                                   rmax,
                                                   rmax,
                                                   maxg_r);

  fclose(file_pair);
        return 0;
}

//-------------------------------------------------------------
//-------------------------------------------------------------
//-------------------------------------------------------------

