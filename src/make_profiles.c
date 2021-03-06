/*  DM 2015 working on id-ing edge and center particles*/
/*  DM 2013 this code has been made into a 2D contour map 
    that ids edge and center of the contour*/
/*  FROM Cynthia Reichhardt */
/*  This program takes a file containing the x and y positions of vortices
    and creates a grayscale bitmap whose intensity is proportional to the
    local flux density.  */ 
/*  The form of the vortex field used is 
    Bz = (1+z)/[2*Pi*(R^2 + (1+z)^2)^(3/2)], 
    which is the z-component of the field, a distance z above the surface,
    and at a radius r from the origin (r as in cylindrical coordinates),
    due to a monopole a distance 1 (i.e. lambda) below the surface.
    To get real units, multiply by Phi_0/lambda^2. */
/*  The number P (below) determines how far out the vortex field is calculated. 
    If the cutoff radius is rc, we have Bz(rc) = P*Bz(0).  Thus P is the 
    fraction by which Bz has dropped from its r = 0 value, when you are at 
    r = rc.
    In terms of P, rc is given as rc = (1+z)*sqrt( P^(-2/3) - 1). */
/*  A final problem is this.  If one integrates Bz over the entire plane at 
    some fixed z, you get the total flux to be 1.  However, if we integrate 
    only out to rc, the	flux is then F = 1 - P^(1/3).  To make the total flux 
    come out correctly, if we approximate Bz with a cutoff of rc, we must 
    rescale it by 1/(1-P^(1/3)) to make the total flux = 1, as it must be.  
    Otherwise, the average field would be too low. */
/* Revisions log:
 * 6.23.08 Digging this code, which originates from Stuart and perhaps also
           Jared, out of mothballs and trying to use it on the new vortex 
	   avalanche project.  Cleaning the code a bit.
*/


#include "Vortex_2D.h"

#define TWOPI	6.283185308

/* X-size of image array in pixels (Y-size is computed) */
#define XSIZE 128  /* 128 */ 

/* X-size of vortex system in lambda */
#define SX  36.5
/* Y-size of vortex system in lambda */
#define SY  36.5
/* Height above surface plus one */
#define Z1  1  		
/* Vortex field cutoff when Bz(r) = P*Bz(0). */
#define P	0.01	/* .005 */
/* Take every SKIP'th frame (SKIP = 1 is every frame) */
#define SKIP	3       
#ifndef SEEK_END
#define SEEK_END 2
#endif

	
void stamp_vortices(int N, struct vortex *vor){

    float **matrix();  

    int i,j,k; 
    int YSIZE,pixelradius,skip;
    int p,q,ix,iy;
    long ii,jj,n,length,star,nextStar,total;
    int nVortices,time;
    float	**image,delta,**stamp,max;
    float	r2,Bz;
    /* Array which contains final Bean profile */
    /*DM including a y profile as well*/
    float	profile_y[XSIZE], profile_x[XSIZE];			

    /*DM Binary binning for clustering*/
    int binary_measure[XSIZE][XSIZE];
    int binary_edge_array[XSIZE][XSIZE];
    int binary_center_array[XSIZE][XSIZE];
    int binary_edge = 0;  //number of '1' grid points without 4nn
    int binary_center = 0; //number of '0' grid points with 4nn
    int total_on = 0;       
    int total_grid_points = XSIZE*XSIZE;

    int total_edge_particles = 0; //count of actual number of particles, not grid

    //binary_density = total_on/total_grid_points

    /* DM: Positions of vortices in vor.x[3] array*/

    /* DM Make a histogram of the image[i][j] values
     * bin 0 to 100 in steps of 2
     */
    int Nh = 50;
    int histogram_vortex_amplitude[Nh];
    float first_nonzero;
    int first_max = 0;        //not doing much with first/last max 
    int last_max[2] = {0,0};  //i,hist
    for(i=0;i<Nh;i++){
      histogram_vortex_amplitude[i]=0;
    }

    /* Cutoff radius for field from a vortex */
    float	radiuscutoff;						
    char file_name[80];
    /* in = input, op = output (all frames), fskip = output */
    FILE 	*in, *out, *fskip, *fedge, *fcenter;
    //there is a MUCH more efficient way to do this!
  
    YSIZE = (SY * (long) XSIZE)/ SX;  //should be square if SY==SX
	
    /* DM Array containing minima (-1), maxima(1) and nothing*/
    int min_max[XSIZE][YSIZE];

    /* Make array to hold integer bitmap */
    image = matrix(0,XSIZE-1,0,YSIZE-1);
	
    /*  size of one pixel, in lambda */
    delta = (float) SX/XSIZE;		
  
    /* radiuscutoff = cutoff (in lambda) */
    radiuscutoff = Z1*sqrt(1/pow(P,2.0/3.0) - 1);		
    if(!quiet_flag) printf("\nVortex field radius cuttoff = %g\n",radiuscutoff);
  
    /* pixelradius = radius in pixels of field stamp */
    pixelradius = radiuscutoff/delta;							
    /* Make array which stamps vortex field */
    stamp = matrix(0,2*pixelradius,0,2*pixelradius);			
  
    /* Fill in vortex field stamp with correct field */
    for (i = 0; i <= 2*pixelradius; i++) {
      for (j = 0; j <= 2*pixelradius; j++) {
	ii = i - pixelradius; 
	jj = j - pixelradius;
	r2 = (ii*ii + jj*jj)*delta*delta;
	Bz = Z1/pow(r2 + Z1*Z1,1.5);
	Bz /= TWOPI;
	/* Rescale to get total flux right */
	Bz /= 1 - pow(P,1.0/3.0);			
	stamp[i][j] = Bz;
      }
    }
  
    /*  Open output files  */
    //out = fopen("output","w");
    if( (out = fopen("output","w")) == NULL) {
      printf("Can't open output file!\n");
      exit(-1);
    }

    //    fskip = fopen("xy_output","w");
    if( (fskip  = fopen("xy_output","w")) == NULL) {
      printf("Can't open output file!\n");
      exit(-1);
    }
  
    skip = 1;
  
  /* Vortex positions stored in vor*/

  /* Zero out image array */  
    for (i = 0; i < XSIZE; i++)
      for (j = 0; j < YSIZE; j++)
	image[i][j] = 0;
    
    /* Now stamp each vortex position. */
    for (k = 0; k < N; k++){	
      p = vor[k].x[0]/delta;
      q = vor[k].x[1]/delta;
      //set all edge_particle variables in vortex struct equal to zero here
      vor[k].edge_particle = 0;
      
      for (i = 0; i <= 2*pixelradius; i++){
	for (j = 0; j <= 2*pixelradius; j++){
	  if (p < 0 || p > XSIZE || q < 0 || q > YSIZE)	
	    /* Data off border. */
	    continue;
	  ix = i - pixelradius + p;
	  iy = j - pixelradius + q;
	  
	  /* Periodic BC's */
	  if (ix < 0)			
	    ix += XSIZE;
	  if (ix >= XSIZE)
	    ix -= XSIZE;
	  
	  if (iy < 0)
	    iy += YSIZE;	
	  if (iy >= YSIZE)
	    iy -= YSIZE;

	  // CIJOL No antivortices. 6.23.08
	  //if (v[k].sign == 1)
	  //image[ix][iy] += stamp[i][j];
	  //else
	  //image[ix][iy] -= stamp[i][j];
	  image[ix][iy] += stamp[i][j];
	  
	}
      }//end for i loop

    }//end k loop
    //-----------------------------------------------------------------------
    /* Array is complete; Look for maxima*/

    /*now average over y */
    
    /* Zero out x and y profile arrays. */
    for (i = 0; i < XSIZE; i++)		
      profile_y[i] = 0;

    for (i = 0; i < YSIZE; i++)		
      profile_x[i] = 0;
    
    //----------------------------------------------------------//
    /* Average image over y, put in profile_y. */
    for (i = 0; i < XSIZE; i++){		
      for (j = 0; j < YSIZE; j++){
	profile_y[i] += image[i][j];

	//create amplitude histogram
	histogram_vortex_amplitude[(int) (Nh*image[i][j])]++;
      }
    }
    //----------------------------------------------------------//

    /*histogram populated, find first nonzero and maxima*/
    for(i=0;i<Nh;i++){
      if(histogram_vortex_amplitude[i]){
	first_nonzero = (float)i/((float)Nh);
      }
      //will only work well for a sharp first minima
      if(!first_max && i &&					\
          histogram_vortex_amplitude[i] > histogram_vortex_amplitude[i-1] &&  \
          histogram_vortex_amplitude[i] > histogram_vortex_amplitude[i+1]){
	first_max = i;
	break;
      }
    }//end for(i histogram loop with first nonzero and first max

    //will only work well for a dominate last maxima
    //NOT a good identifier
    //UGLY and not terribly accurate
    for(i=Nh-2;i>first_max;i--){
      if( !last_max[1] && histogram_vortex_amplitude[i] > last_max[1] && histogram_vortex_amplitude[i] > histogram_vortex_amplitude[i-1] &&  histogram_vortex_amplitude[i] > histogram_vortex_amplitude[i+1]){
	last_max[0] = i;  //will be written over
	last_max[1] = histogram_vortex_amplitude[i];  //will be written over
      }//end if statement
    }//end for i loop
        

    if(!quiet_flag) printf("first nonzero amp at %f\n",first_nonzero);
    if(!quiet_flag) printf("first max at %f\n",(float)first_max/((float)Nh));
    if(!quiet_flag) printf("last max at %f,%d\n", (float)last_max[0]/((float)Nh),last_max[1]);
    


    /* Average image over x, put in profile_x. */
    for (j = 0; j < YSIZE; j++){		
      for (i = 0; i < XSIZE; i++){
	profile_x[j] += image[i][j];	
      }
    }

    /* Print out profile array. */
    /* DM cheat, XSIZE=YSIZE */
    for (i = 0; i < XSIZE; i++)	{	
      fprintf(out, "%g\t%g\t%g\n",(i+0.5)*delta, \
                                  profile_y[i]/YSIZE, \
                                  profile_x[i]/XSIZE);

      //print xy_data to "xy_output" yay!
      for(j=0; j< YSIZE; j++){

	binary_measure[i][j]      = 0;  //assume zero
	binary_edge_array[i][j]   = 0;  //assume zero
	binary_center_array[i][j] = 0;  //assume zero

	/*arbitrary cutoff to make a 1/0 grid of the vortex density
	 *populated_density = total_on/total_grid_points
	 *should be compared with 2d density \rho = N / SX*SY
	 
	 *a particle with cutoff = 0.13 is 7x7, 
	 *and should be counted at part of binary measure,
	 *however, edge vs. center needs some work!!!
	 
	 */

	if(image[i][j] > 0.13){  //if above cutoff
	  binary_measure[i][j] = 1;    //turn to 'on'
	  total_on++;            //track total ones in system
	}
	fprintf(fskip, "%g\t%g\t%g\t%d\n",(i+0.5)*delta, (j+0.5)*delta, \
                                           image[i][j], binary_measure[i][j]);
      }//end for(j=0 loop

      fprintf(out,"\n");    //these blanks lines make contour plotting possible
    fprintf(fskip,"\n");
    }
    fprintf(out,"\n");

    fclose(out);
    fclose(fskip);

    //------------------------------------------------------------------
    //find minima/maxima...
    //------------------------------------------------------------------
    //ungodly, but it works!!!
    //returns many values for messy systems

    for (i=0; i<XSIZE; i++)	
      for(j=0; j<YSIZE; j++)
	min_max[i][j]=0;
    //-----------------------//end zeroing for loops

    fskip = fopen("maxima_xy_output","w");

    if (fskip == NULL) {
      printf("Can't open output file!\n");
      exit(-1);
    }

    //open the edge files, but consolidate these soon!!!
    if( (fedge  = fopen("xy_edge","w")) == NULL) {
      printf("Can't open output file!\n");
      exit(-1);
    }

    if( (fcenter  = fopen("xy_center","w")) == NULL) {
      printf("Can't open output file!\n");
      exit(-1);
    }
    //----------------------------------------------all files open

    for (i=1; i<XSIZE-1; i++)	{	
      for(j=1; j<YSIZE-1; j++){

	/*Identify min/max of density contour*/
	if( (image[i][j] > image[i+1][j]) &&		\
	     (image[i][j] > image[i-1][j]) &&		\
	      (image[i][j] > image[i+1][j+1]) &&	\
	       (image[i][j] > image[i-1][j-1]) &&	\
		(image[i][j] > image[i+1][j-1]) &&	\
		 (image[i][j] > image[i-1][j+1]) &&	\
		  (image[i][j] > image[i][j+1]) &&	\
		   (image[i][j] > image[i][j-1]) ){
	  min_max[i][j] = 1;
	  fprintf(fskip,"%g\t%g\t%g\n",(i+0.5)*delta, (j+0.5)*delta, image[i][j]);
	}

	//now determine edge/center, even print it for trouble shooting
	//this makes beautiful plots that clearly define the actual edges
	//however you made need to put this back into a particle space 
	//to get more accurate results
	//since you want to eventually combine with a p6 measure 
	//that will ultimately be a good thing

	//if you simply print this data, how will you compare to 
	//individual particle positions?
	//unfortunately you are going to have to save it
	if(binary_measure[i][j]){
	  int binary_check = binary_measure[i+1][j]	\
	    + binary_measure[i-1][j]			\
	    + binary_measure[i][j+1]			\
	    + binary_measure[i][j-1];	  
	  if (binary_check < 4 && binary_check){  //don't take the zero particles
	    //edge particle
	    binary_edge++;
	    binary_edge_array[i][j]   = 1;  
	    fprintf(fedge, "%g\t%g\n",(i+0.5)*delta, (j+0.5)*delta );
	  }
	  else if(binary_check == 4) {
	    //inner particle
	    binary_center++;
	    binary_center_array[i][j] = 1;  
	    fprintf(fcenter,  "%g\t%g\n",(i+0.5)*delta, (j+0.5)*delta );
	  }
	}//end binary - edge/center determination

      }
    }
    fclose(fskip);
    fclose(fedge);
    fclose(fcenter);
    //--------------------------------------------------------------------

    /*DM 7.18.2015 
     *Here is where we connect particle positions 
     *with edge/center determination
     *so far too much hard coding 
     *works well of isolated particles and low density system
     *fails for high density

     *ADD PERIODIC BOUNDARY CONDITIONS
     *MEASURE LOCAL DENSITY TO0 - DMIN?
     *detect edge dislocations?
     *could you do a relative measure of edge?  
     *problem is to determine what grid to surround the particles with
     *with is local density dependent.
     
     *if you know xmin, you could search in a regime just greater than that

     *is there a way to normalize from binary_edge to particle edge?
     */

    //Loop through the particles:
    for (k = 0; k < N; k++){	

      //the gridded particle positions
      //you should be able to reduce these to integers 
      //and use the look up in the binary_edge_array
      //to see how close they are.
      //hard code, but if they are <4 away from an edge, they're an edge

      if(0){
	printf("%f,%f,%f\n", vor[k].x[0], vor[k].x[1],delta);
      }

      p = vor[k].x[0]/delta;
      q = vor[k].x[1]/delta;

      int num_edge=0;

      //check this more fully, but we can't access array elements
      //that we don't have

      int gb=6;
      
      for ( i=(p-gb); i<(p+gb+1); i++ ){ 
	for ( j=(q-gb); j<(q+gb+1); j++ ){ 

	  if (i < 0 || i > XSIZE || j < 0 || j > YSIZE){
	    continue;
	  }
	  else{

	    if(binary_edge_array[i][j]){
	      num_edge++;
	    }
			       
	  }//end else we have the case we want, 

	}//end i loop
      }//end j loop

      if(0){
	printf("%d,%d,%d\n", p, q, num_edge);
	fflush(stdout);
      }

      if (num_edge>5){
	vor[k].edge_particle = 1;
	total_edge_particles++;
      }

    }//end k loop

    /*------------------------------------------------------------------*/


    //--------------------------------------------------------------------
    fskip = fopen("minima_xy_output","w");

    if (fskip == NULL) {
      printf("Can't open output file!\n");
      exit(-1);
    }

    for (i=1; i<XSIZE-1; i++)	{	
      for(j=1; j<YSIZE-1; j++){

	if( (image[i][j] < image[i+1][j]) &&		\
	     (image[i][j] < image[i-1][j]) &&		\
	      (image[i][j] < image[i+1][j+1]) &&	\
	       (image[i][j] < image[i-1][j-1]) &&	\
		(image[i][j] < image[i+1][j-1]) &&	\
		 (image[i][j] < image[i-1][j+1]) &&	\
		  (image[i][j] < image[i][j+1]) &&	\
		   (image[i][j] < image[i][j-1]) ){

	  min_max[i][j] = -1;	  
	  fprintf(fskip,"%g\t%g\t%g\n",(i+0.5)*delta, (j+0.5)*delta, image[i][j]);
	}
      }
    }

    fclose(fskip);
    //--------------------------------------------------------------------
    //Debug Max/Min Array
    //--------------------------------------------------------------------
    /*
    printf("\n max array printing...\n");
    fflush(NULL);
    for (i=0; i<XSIZE; i++){	
      for(j=0; j<YSIZE; j++){
	printf("%d  ",min_max[i][j]);
      }
      printf("\n");
    }
    */
      //-----------------------------------------------------------

    /*Write the histogram*/
    out = fopen("histogram_vortex_amplitude","w");

    if (out == NULL) {
      printf("Can't open output file!\n");
      exit(-1);
    }

    for(i=0;i<Nh;i++){
      //if(histogram_vortex_amplitude[i]) first_nonzero = i;
      //print pinning force for contour plots of histograms
      fprintf(out,"%f\t%f\t%d\n",pinning_force, ((float)i)/((float)Nh),histogram_vortex_amplitude[i]);
    }

    fprintf(out,"#first nonzero amp at %f\n",first_nonzero);
    fprintf(out,"#first max at %f\n",(float)first_max/((float)Nh));
    fprintf(out,"#last max at %f,%d\n", (float)last_max[0]/((float)Nh),last_max[1]);
    fprintf(out,"\n");

    fclose(out);

  if(!quiet_flag) printf("\nAll Done with Vortex Contour!\n");

  //populated_density = total_on/total_grid_points
  if(!quiet_flag) printf("local density roughly %f\n",(float)total_on/total_grid_points);
  if(quiet_flag) printf("%f\t%f\t%f\t%f\t%f\t",(float)total_on/total_grid_points,\
			(float)binary_edge/total_grid_points,\
			(float)binary_center/total_grid_points,	\
			(float)total_edge_particles/(float)N,	\
			(float)first_max/((float)Nh),\
			(float)last_max[0]/((float)Nh) );
  //
  //should be compared with 2d density \rho = N / SX*SY
  //


}

/*-------------------- matrix() --------------------------------------*/


float **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
  void nrerror();
  int i;
  float **m;
    
  m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m -= nrl;
    
  for(i=nrl;i<=nrh;i++) {
    m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
    if (!m[i]) nrerror("allocation failure 2 in matrix()");
    m[i] -= ncl;
  }
  return m;
}

/*-------------------- nrerror() --------------------------------------*/

void nrerror(error_text)
char error_text[];
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}
    
