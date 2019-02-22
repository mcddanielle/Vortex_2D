#include "Vortex_2D.h"

/*Layer Order Parameter:
 *\Psi_{layer,n_l} = 
 *|\frac{1}{n_{bin} \sum_{j=1}^{n_{bin}} exp{i \frac{2\pi(n_l-1)}{Ly} y_j}|
 *
 *re-write to be pin-specific
 */

int layer2D(int N, struct vortex *vor){

  //set an x bin / y bin size, loop through testing for x/y in range, 
  //add to particular sum

  //start with xbin = pinning_period; ybin = xbin
  int num_bin = 10;

  float xbin = xmax / (float) num_bin; 
  float ybin = xbin;
 
  int hist_xbin[num_bin];  //keep track of n_bin
  int hist_ybin[num_bin];

  gsl_complex psi_x[num_bin];
  gsl_complex psi_y[num_bin];

  float Ly = ymax; //simply to match notation

  int i,nl;

  float avg_psi_x, avg_psi_y;

  gsl_complex prefactor_nl, psi_jx, psi_jy, old_psi_x, old_psi_y; 

  if((layer_file = fopen("LayerOP_2D.txt","w"))==NULL) {
    printf("Cannot open file.\n");
    exit(1);
  }

  //start at the nl-1=1 case
  for(nl=3;nl<7;nl++){ 
    prefactor_nl = gsl_complex_polar(1.0, 2*PI*(nl-1)/Ly );
    
    //zero for every value of nl
    for(i=0;i<num_bin;i++){
      hist_xbin[i]=0;
      hist_ybin[i]=0;    
      psi_x[i] = gsl_complex_polar(0.0, 0.0 );
      psi_y[i] = gsl_complex_polar(0.0, 0.0 );
    }  

    //printf("prefactor=%f\n",gsl_complex_abs(prefactor_nl));

    for(i=0;i<N;i++){
      //printf("i,nl=%d,%d\t",i,nl);
      //fflush(NULL);
      int x_bin_loc = vor[i].x[0] / xbin;
      int y_bin_loc = vor[i].x[1] / ybin;

      //printf("x,y=%f,%d,%f,%d\t",vor[i].x[0],x_bin_loc,vor[i].x[1],y_bin_loc);
      //fflush(NULL);

      if (x_bin_loc < num_bin && x_bin_loc > -1 )
	hist_xbin[ x_bin_loc ]++;
      else 
	return -1;
      if (y_bin_loc < num_bin && y_bin_loc > -1)
	hist_ybin[ y_bin_loc ]++;
      else 
	return -1;

      double x = vor[i].x[0];
      double y = vor[i].x[1];

      psi_jx = gsl_complex_pow_real(prefactor_nl,x );
      psi_jy = gsl_complex_pow_real(prefactor_nl,y );

      //printf("psi_jx,psi_jy=%f,%f\t",gsl_complex_abs(psi_jx),gsl_complex_abs(psi_jy));
      //fflush(NULL);

      psi_x[x_bin_loc] = gsl_complex_add(psi_x[x_bin_loc],psi_jx);
      psi_y[y_bin_loc] = gsl_complex_add(psi_y[y_bin_loc],psi_jy);

      //printf("psix,psiy=%f,%f\n",gsl_complex_abs(psi_x[x_bin_loc]),gsl_complex_abs(psi_y[y_bin_loc]));
      //fflush(NULL);

    }//all vortices accounted for, print!

    fprintf(layer_file,"#n nl = %d \n",nl);   
    for(i=0;i<num_bin;i++){
      avg_psi_x = gsl_complex_abs( gsl_complex_div_real(psi_x[i],hist_xbin[i]));
      avg_psi_y = gsl_complex_abs( gsl_complex_div_real(psi_y[i],hist_ybin[i]));
      fprintf(layer_file,"%f\t%f\t%f\t%f\n",i*xbin,avg_psi_x,i*xbin,avg_psi_y);
    }
    fprintf(layer_file,"\n",nl);   

  }

     fclose(layer_file);


  return 0;
}
