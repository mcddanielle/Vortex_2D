#include "Vortex_2D.h"

void print_movie( int num_vor, int time, float temp, float current, float *x, float *y, float *vx, float *vy){

   int i;

   fwrite(&num_vor,sizeof (int),1,trans_movie);
   fwrite(&time,sizeof (int),1,trans_movie); 
  
   //If the restart flag is zero, write a jeff movie
   if(!restart_flag){
     fwrite(&temp,sizeof (float),1,trans_movie);
     fwrite(&current,sizeof (float),1,trans_movie); 

     fwrite(&xmin,sizeof (float),1,trans_movie);
     fwrite(&xmax,sizeof (float),1,trans_movie); 
     fwrite(&xmin,sizeof (float),1,trans_movie);  //NOT redundant xmin=ymin
     fwrite(&xmax,sizeof (float),1,trans_movie);  //xmax=ymax
   }

   for(i=0;i<num_vor;i++){

    fwrite(&i,sizeof (int),1,trans_movie);
    fwrite(&x[i],sizeof (float),1,trans_movie);
    fwrite(&y[i],sizeof (float),1,trans_movie);

   //If the restart flag is zero, write a jeff movie
    if(!restart_flag){
      fwrite(&vx[i],sizeof (float),1,trans_movie);
      fwrite(&vy[i],sizeof (float),1,trans_movie); 
    }

   }

  return;
}
