#ifndef SS_LIB_H
#define SS_LIB_H

#include "block_dims.h"

#define TRUE 1
#define FALSE 0

#define CUCALL( call )               \
{                                       \
cudaError_t result = call;              \
if ( cudaSuccess != result )            \
    fprintf(stderr, "CUDA error  %s \n", cudaGetErrorString( result ) ); \
}

//Sphere radius
#define R 64/2


typedef real bl_array[Bz][By][Bx];
struct block {
    bl_array u;
    bl_array u1;
    char k[Bz][By][Bx];
    int up; //z direction
    int down;
    int left; //x direction
    int right;
    int fore; //y direction.
    int aft;
};

typedef struct {
    //contains all data for a aos based representation of a room.
    int X;
    int Y;
    int Z;
    char* k; //array of k values (can be 0-6)
    int blocks_in; //total number of internal blocks
    struct block* aos; //the main array
    int* index_of_struct; //Array containing index of struct based on blocks coords.
} ss_t ;

ss_t create_aos(int X, int Y, int Z, char* is_in_sphere, int* blocks_in, struct block** aos_pp, int** index_of_struct);
int get_coords(int coords[3], int X, int Y, int Z, int* index_of_struct, int* io_block_ind, int arrind[3]);
int free_ss(ss_t data);
real *hanning_window(int big_n);

#endif

