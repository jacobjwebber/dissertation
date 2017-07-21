
#define TRUE 1
#define FALSE 0
#define real double

#define CUCALL( call )               \
{                                       \
cudaError_t result = call;              \
if ( cudaSuccess != result )            \
    fprintf(stderr, "CUDA error  %s \n", cudaGetErrorString( result ) ); \
}

//Block size
#define Bz 2
#define By 2
#define Bx 32

//Sphere radius
#define R 64/2


typedef real bl_array[Bz][By][Bx];
struct block {
    int up; //z direction
    int down;
    int left; //x direction
    int right;
    int fore; //y direction.
    int aft;
    bl_array u;
    bl_array u1;
    char k[Bz][By][Bx];
};

int create_aos(int X, int Y, int Z, char* is_in_sphere, int* blocks_in, struct block** aos_pp, int** index_of_struct);
int get_coords(int coords[3], int X, int Y, int Z, int* index_of_struct, int* io_block_ind, int arrind[3]);
