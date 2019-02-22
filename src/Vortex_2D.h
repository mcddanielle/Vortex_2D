/*Author Danielle McDermott, November 1, 2013
 *stee ABOUT file
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <unistd.h>
//command line arguments
#include <getopt.h>
#include <errno.h>
#include <fcntl.h>
//math libraries
#include <gsl/gsl_sf_legendre.h>  
#include <gsl/gsl_sf_result.h>
#include <complex.h>
#include <gsl_complex.h>
#include <gsl_complex_math.h>
//Dan's library
//#include <parameters.h>

#ifndef PYTH_H_GUARD   
#define PYTH_H_GUARD

#define DEBUG 0
#define PI 3.14159265359
#define C 30

//file pointers
FILE *open_movie;  //open the binary movie
FILE *kmovie;      //write the binary final frame (for voronoi)
FILE *trans_movie; //translate entire binary movie
FILE *ang_file;
FILE *ang_hist_file;
FILE *local_density_file;
FILE *layer_file;
//----------------------struct definitions----------------------//

//linked vortex list
struct vortex {
  int id;
  float x[3];   //x,y,(z?)->x[3]
  int pin_number;
  float pin_dist;
  int bondcount;
  float d_neighbor;
  float x_min;
  float y_min;
  float y_min2;
  int edge_particle;      //basically a bool that is set to 1 if edge particle
  struct vortex *bt_new;
  struct vortex *bt[C];   //C is chosen artifically-INEFFICIENT!
                          //the array works better for one
                          //reason, you don't want a second linked
                          //list, BUT you can dynamically malloc
                          //this array
  //struct vortex *next;    //this allows for linked list, not too scary
};


//linked neighbor list for g(r) calculations
struct vor_bond {
  struct vortex *vor1, *vor2;
  float length;
  float vec_x[3];      //the vector! 
  int bond_type;       //where 0=in trap
                       //1 nn trap
                       //2 nnn trap  
  struct vor_bond * next;
};

struct pin_site {
  float x_pin;
  int pin_id;   //unnecessary

  //pointer i.e. an array containing y points
  float *y_pos;  
  //limits this to use in the 1D case
  //would be better with a linked list id,x,y
  int num_vor;
};


//------------------- Function prototypes----------------------------//
void getoptions( int, char ** );
void print_frame( int num_vor, float temp, float current, float *x, float *y);
void print_movie( int num_vor, int time, float temp, float current, float *x, float *y, float *vx, float *vy);
void print_ascii_frame( int num_vor, float temp, float current, float *x, float *y);
void read_single_ascii();

void get_global_defaults();
void check_globals_for_errors();

void find_pinning_minima(float *pin_minima);

void stamp_vortices(int N, struct vortex *vor);

//Geometry calculations
float bondlength(float * x1, float * x2);
float periodic_bondlength(float * x1, float * x2);
float vec_length(float *x1);
float estimate_neighbor_distance(float fp);

//NOT explicitly listing args
//admittedly because ptrs are hard for me... should be fine now
int pair_correlation();         //(int N, struct vor_bond *bond);
int angle3D_correlation();      //(int N, struct vortex[N])
int density2D();                //(int N, struct vortex[N])
int layer2D(int N, struct vortex[N]);
float avg_dist_neighbor(int N, struct vortex *vor, float *std_dev, float *max_dist);
void find_nearest_y_neighbors(int N, struct vortex *vor, float *pin_location);
void pin_centered_y_neighbors(int N, struct vortex *vor, float *pin_location);

int check_1D(int N, struct vortex *vor);
int find_pin(float x, float *pin_location);
//function definition!!!
struct vor_bond *load_movie_ptr(int *N, float *pin_location);//, struct vortex *vortex_list); 
//, float *x);//, struct vortex **vortex_list);
void angle2D_correlation(int N, struct vortex *vor, float bondlength);

//comparison for qsort in density.c
int compare_floats (const void *a, const void *b);

/*------------------ Global variables defined here--------------------//
 * getoptions() gets runtime values
 * assign_global_defaults() and check_globals_for_errors() do the rest
 */ 

  float xmin, xmax, ymin, ymax;
  float frame_print_current;
  float pinning_period;
  float pinning_force;
  float xsi;          //interaction length scale
  float bondLim;
  int number_pins;

  int verbose_flag;   // Print out lots of debugging statements if set
  int quiet_flag;     // DO NOT! print normal program output if set

  int restart_flag;     //repurposed.  defined in top Pof main.c
  int movie_type;
  int write_movie_type;

  char filename[100];
  char out_file[120];


//----------------------------------------------------//
//------------- Vector Macros-------------------------//
//----------------------------------------------------//
#define sign2(x) (( x > 0 ) - ( x < 0 ));

#define dotProduct(a,b,c) \
  a = (b)[0] * (c)[0] + (b)[1] * (c)[1] + (b)[2] * (c)[2];

#define crossProduct(a,b,c) \
    (a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2]; \
    (a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0]; \
    (a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1];
#define vecSum(a,b,c) \
    (a)[0] = (b)[0] + (c)[0]; \
    (a)[1] = (b)[1] + (c)[1]; \
    (a)[2] = (b)[2] + (c)[2];
#define vecDiff(a,b,c) \
    (a)[0] = (b)[0] - (c)[0]; \
    (a)[1] = (b)[1] - (c)[1]; \
    (a)[2] = (b)[2] - (c)[2];
#define vecNeg(a,b) \
    (a)[0] = -(b)[0]; \
    (a)[1] = -(b)[1]; \
    (a)[2] = -(b)[2];

#endif           
//--------------------------------end the PYTH_H_GUARD def


/*  HISTOGRAM STUFF
//void print_histogram_data(int num_vor, int *vx, int *vy, int *speed, float current, int binNum, int high, int low, int normalize);
//int bin_data(float binMe, float low, float high, int binNum);
//int int_cmp(const void *a, const void *b);  //for sorting

  float hist_range;  //max speed of particles
  int hist_flag;
  int hist_dimension;
  int bin_zeros ;

  char histogram_file_name[50];
  char histogram_peak_file[50];
  char histogram_peak_vy_file[50];

*/
