#include "Vortex_2D.h"

/*
 * For uncorrelated systems, the order parameter Q6 = 1/√12 = 0.289.
 *
 *include Tetragonality Metric   
 *q = 1 - 3/8 \sigma_{j=1}^3 \sigma_{k=j+1}^4 (Cos(theta_jik)+1/3)^2

"This definition is a rescaled version of the parameter introduced in ref. 4, in such a way that the average value of q varies between 0 (in an ideal gas) and 1 (in a perfect tetrahedral network)."

AU  - Errington, Jeffrey R., Debenedetti, Pablo G.
TI  - Relationship between structural order and the anomalies of liquid water

JA  - Nature
PY  - 2001/01/18/print
VL  - 409
IS  - 6818
SP  - 318
EP  - 321
SN  - 0028-0836
UR  - http://dx.doi.org/10.1038/35053024
 */

/*
 * struct atom_data contains an array of ptrs: atm.bt[C]
 *
 * creates an atom centered bond network, each atom has an array
 * of pointers to each bonded atom 
 *
 * First calculate the local order parameter ql[i] for each atom
 *
 * Next calculate the global order parameter Ql for N atoms 
 */

FILE *ang_file;

int angle_coor( int N, struct vortex *vor )
{
  int i,ip,j;

  int l,m;        //lp counts 0,1,2 while l counts 4,6,8
  int lp = 0;     //lp is the array variable, l is the Y_lm index

  double ql[4][N];   //used in the counting loops for summing
  double q2, q4, q6, q8; //local order parameter for each atom  l = 4,6,8

  double Ql[4] = {0.0, 0.0, 0.0, 0.0}; //global order parameter 
  double Q2, Q4, Q6, Q8;

  int N_nb;   //bonds for each atoms
  double Nt;  //total bonds to normalize by

    if((ang_file = fopen("BondOrder.txt","w"))==NULL) 
    {printf("Cannot open file.\n");
      exit(1);
    }

    fprintf(ang_file,"#i     q2[i]    q4[i]    q6[i]    q8[i]    N_nb   \n");

  //----------------------histogram piece----------------//
    int hist_max = 1000;
  int hist_value_q4, hist_value_q6, hist_value_q8;
  int hist_q4[hist_max];
  int hist_q6[hist_max];
  int hist_q8[hist_max];

  for(i=0;i<hist_max;i++){
    hist_q4[i]=0;
    hist_q6[i]=0;
    hist_q8[i]=0;
  }
  //----------------------end histogram piece------------//

  /* The gsl_sf library contains needed error handling.
   * The defined struct is used in the calculation of 
   * the value (result1.val) and error (result1.err)
   * in the calculation of the legendre polynomial P_lm
   */

  gsl_sf_result result1, result2; 

  /*----------gsl_complex library to handle P_lm properly-----------*/

  //struct gsl_complex z = gsl_complex_polar (double r, double theta)
  gsl_complex ql_m, qlm, temp1, temp2, Q_c[4][17];

  /*Q_c[lp][17] is a way to keep track of the global order
   * have to map for(m=0;m<(l+1)... (Plm calculations)
   * to m = -l...0...l = -8...0...8
   *   m' =  0....2l+1 =  0...8...16  =17 elements
   *  -m' = 8 - m
   *   m' = 8 + m
   *notice that the Q_c[0] and Q_c[1] will have null terms
   *since Q4 and Q6 only have 9 and 13 elements respectively
   */
  //zero all global elements (just in case)

  for(j=0;j<4;j++){
  for(i=0;i<17;i++){

    Q_c[j][i] = gsl_complex_polar(0.0,0.0);

     }//end j
   }//end i

  //--------------------------------------------------//

  /*The following variables are used to re-center the origin at 
   *atom i, and define a vector to atom j w/out changing the 
   *overall orientation of the cartesian coordinate system.
   */

    double xj[3];                 //local vector magnitude r_ij
    double r_j, cos_th_j, phi_j;  //local length, theta, phi values

  for(i=0;i<N;i++){
    
    N_nb = vor[i].bondcount;    //make a shorthand for loops
    Nt += N_nb;                 //running total for normalization of Ql

    ql[0][i] = 0.0;    //individual sum for l=2
    ql[1][i] = 0.0;    //individual sum for l=4
    ql[2][i] = 0.0;    //individual sum for l=6
    ql[3][i] = 0.0;    //8

    /*the j loop calculates the local parameters 
     *for each vor[i].bt[j < N_nb] 
     *where *bt ="bonds to" points to a different atom struct
     */

    /*-----------------ang_p = angular parameters-----------------
     *declare angular parameters for each bond_ij in the loop
     *i.e. atom_i bonds to N_nb atom_j's (nb = number of bonds)
    */

    double ang_p[N_nb][2]; 

    //j loop to calculate local coordinates stored 
    //ang_p[N_nb][2] = cos(theta_ij),phi_ij

 for(j=0;j<N_nb;j++){

           vecDiff(xj, vor[i].bt[j]->x, vor[i].x); //Macro in pytha.h
           r_j = vec_length(xj);                   //subroutine in bondcalc.c
           cos_th_j = xj[2] / r_j;        //cos(th)=z/r
           phi_j = atan2(xj[1] , xj[0]);  //phi=arctan(y/x)

//if(!i) printf("%d,cos(th)=%lf,phi=%lf\n",j+1,cos_th_j,phi_j);

	     ang_p[j][1] = cos_th_j;  //assign to an array so 
	     ang_p[j][2] = phi_j;     //values not recalculated
                                      //for each j loop (inner loop)
 }//end j loop

      //---------------end coordinate calculations-----------------//

     /*
      *angles have been calculated, do the spherical harmonics 
      */

      lp = 0;      //could do a for loop for lp also, but overkill

      for(l=2;l<9;l+=2){    

         for(m=0;m<(l+1);m++){

           //zero the local values in z = re^(i*theta) form
	   qlm = gsl_complex_polar(0.0,0.0);
	   ql_m = gsl_complex_polar(0.0,0.0);

	   for(j=0; j<N_nb; j++){    
       
	     /*-----------------P_lm_error handling:----------------
	      *http://www.gnu.org/software/gsl/manual/html_node/
               Special-Function-Usage.html
             
               *the gsl_function_e takes l,m,cos(theta) and stores
	       *in result1.value and result1.error
 	      -------------------------------------------------------*/

             //Plm given cos(theta)
	     gsl_sf_legendre_sphPlm_e(l, m, ang_p[j][1],&result1);
             
             //temp2 = Plm * e^(i*m*phi)
	     temp2 = gsl_complex_polar( result1.val, m*ang_p[j][2]);
	   
	     //---------temp1 is the negative m values--------------//
	     //---------------see notes above on Q_c----------------//

             if(m!=0){

             //temp1 = (-1)^m * Plm * e^(-i*m*phi)
	     temp1 = gsl_complex_polar(pow(-1.0,m)*result1.val,-m*ang_p[j][2]);

             //ql_m += temp1... 
	     //difference with Q_c is not zeroed like qlm
	      ql_m = gsl_complex_add(ql_m,temp1);
 	      Q_c[lp][8-m] = gsl_complex_add(Q_c[lp][8-m],temp1);

	     }//end if(m!=0) for negative m values

	     //--------positive m values----------------------------//

             //qlm += temp2
	     qlm = gsl_complex_add(qlm,temp2);
             Q_c[lp][8+m] = gsl_complex_add(Q_c[lp][8+m],temp2);

	     }//end (j loop)



           /*
	    *update qlm when completing j loop (all ij angles)
	    *Q_c[l][m'] will be averaged when i loop terminates

	    *normalize qlm and ql_m 
	    *by by dividing by the number of bonding pairs N_nb
	    */
	  if(N_nb){
	   ql_m = gsl_complex_div_real(ql_m, N_nb);
	    qlm = gsl_complex_div_real( qlm, N_nb);
	     }

	/*Careful:  
	 *ql[lp][i] will become q2[atom i] when all bonds processed
	 *now it is only sums 
	 *without multiplication by prefactors
	 *or the final squareroot
	 */
	    //complex_abs2 gives (psi*psi)
	   ql[lp][i] += gsl_complex_abs2(ql_m) + gsl_complex_abs2(qlm);

      } //end---------------------for(m=-l...)

      lp++;

    }//end---------------------------------for(l=0...l+=2)

      //report the local parameter for each i value
     q2 = sqrt(2*PI/ 5 * ql[0][i]);	      
     q4 = sqrt(4*PI/ 9 * ql[1][i]); 
     q6 = sqrt(4*PI/13 * ql[2][i]);
     q8 = sqrt(4*PI/17 * ql[3][i]);

     //--------------histogram increment -------------------------//
     /*
     hist_value_q4 = hist_max*q4;
     hist_value_q6 = hist_max*q6;
     hist_value_q8 = hist_max*q8;

     hist_q4[hist_value_q4] += 1;
     hist_q6[hist_value_q6] += 1;
     hist_q8[hist_value_q8] += 1;
     */
     //-------------end histogram increment-----------------------//

     fprintf(ang_file,"%2d %9.6lf %9.6lf %9.6lf %9.6lf %2d\n",
       i+1, q2, q4, q6, q8, N_nb);

  }//end-------------------------for(i=0;i<N;i++)

  //------------------ Global Q_c array ----------------------//

  for(lp=0;lp<4;lp++){
    for(m=0;m<17;m++){
  
     //sum over m: Ql+=(psi*psi) of local Qlm
     Ql[lp] += gsl_complex_abs2(Q_c[lp][m]);
 
    }//m loop
  }//lp loop

  //-------------------normalize and squareroot--------------//

     Q2 = sqrt(4*PI/ 9 * Ql[0]) /Nt; 
     Q4 = sqrt(4*PI/ 9 * Ql[1]) /Nt; 
     Q6 = sqrt(4*PI/13 * Ql[2]) /Nt;
     Q8 = sqrt(4*PI/17 * Ql[3]) /Nt;	  

     printf("         Q2 = %lf Q4 = %lf  Q6 = %lf  Q8 = %lf Nt = %f\n",Q2,Q4,Q6,Q8,Nt);
     fprintf(ang_file,"\n #Q2 = %lf Q4 = %lf  Q6 = %lf  Q8 = %lf Nt = %f\n",Q2,Q4,Q6,Q8,Nt);

     fclose(ang_file);

  return 0;

     FILE *file_ql;

  if((file_ql = fopen("ql_data.txt","w"))==NULL){
      printf("Cannot open file.\n");
      exit(1);
    }

     fprintf(file_ql, "   i     q4    q6   q8\n");
     for(i=0;i<hist_max;i++){
       fprintf(file_ql,"%6d %6d %6d %6d\n",i,hist_q4[i],hist_q6[i],hist_q8[i]);
     }

     fclose(file_ql);

  return 0;
}


