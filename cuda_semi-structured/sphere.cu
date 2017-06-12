#include <stdio.h>
#include <math.h>

#define TRUE 1
#define FALSE 0

#define real double

//Block size
#define Bl 8
#define Bm 8
#define Bp 4

//Bounding box dimensions
#define L 8 //64
#define M 8 //64
#define P 8 //64

//Sphere radius
#define R 4 //64/2

struct block {
    int up; //z direction
    int down;
    int left; //x direction
    int right;
    int fore; //y direction.
    int aft;
    float u[Bl][Bm][Bp];
    float u1[Bl][Bm][Bp];
    char k[Bl][Bm][Bp];
    int coords[3];
};

__global__ void perform_stencil_internal(struct block *u, struct block *u1, real l2, real l, real g);

//DUPLICATE
char is_block_in(int x, int y, int z, char is_in_sphere[L][M][P]) {

    int i,j,k;
    for (i = 0; i < Bl; i++) {
        for (j = 0; j < Bm; j++) {
            for (k = 0; k < Bp; k++) {
                if (is_in_sphere[i][j][k]) {return TRUE;}
            }
        }
    }
    return FALSE;
}

int coords_to_index(int x, int y, int z) {
    int area = M*P;
    return x*area + y*P + x;
}

char is_block_internal(int x, int y, int z, char *array, int xmax, int ymax, int zmax) {
    //This function takes coordinates and returns true if any points within a 
    // blocksize starting on that point are inside the room.
    int i, j, k;
    for (i = x; i < x+Bl; i++) {
        for (j = y; j < y+Bm; j++) {
            for (k = z; k < z+Bp; k++) {
                if ( array[i*ymax*zmax + j*zmax + k]) {
                    //printf("boop , %i %i %i \n", x, y, z);
                    return TRUE;
                 }
            }
        }
    }

    return FALSE;
}

   


int main() {

    //Coefficients
    //courant number and courant number squared.
    real l2 = 1.0/3.0;
    real l = sqrt(l2);
    //Wall reflection coefficient.
    real r = 0.9;
    real g = (1-r)/(1+r);
    //grid spacing (m)
    real h = 0.1;
    //Speed of sound
    real c = 343;

    //Error checking - Bounding grid must be divisible by block.
    if(L%Bl != 0 || M%Bm != 0 || P%Bp != 0) {
        printf("Not divisible by block\n");
        return 0;
    }

    //Define array of grid points within sphere.
    char is_in_sphere[L][M][P];
    int i,j,k;
    for (i = 0; i < L; i++) {
        for (j = 0; j < M; j++) {
            for (k = 0; k < P; k++) {
                is_in_sphere[i][j][k] = sqrt( (i-R)*(i-R) + (j-R)*(j-R) + (k-R)*(k-R)) < (float) R;
            }
        }
    }

    /*
    //Print slice in middle
    for (i = 0; i < L; i++) {
        for (j = 0; j < M; j++) {
            printf("%i ", is_in_sphere[i][j][R]);
        }
        printf("\n");
    }
    // */

    //Assign a block for all volumes containing points
    //by first assigning enough blocks for the whole lot.

    int num_blocks_l = L/Bl;
    int num_blocks_m = M/Bm;
    int num_blocks_p = P/Bp;
    int total_blocks = num_blocks_l*num_blocks_m*num_blocks_p;

    //Create an array storing the location of each block.
    int *index_of_struct = (int*) calloc(total_blocks, sizeof(int));

    struct block *array = (struct block *) calloc(total_blocks, sizeof(struct block));

    //total number of internal blocks.
    int blocks_in = 0;

    for (i = 0; i < num_blocks_l; i++) {
        for (j = 0; j < num_blocks_m; j++) {
            for (k = 0; k < num_blocks_p; k++) {
                if ( is_block_internal(i * Bl, j * Bm, k * Bp, &(is_in_sphere[0][0][0]), L,M,P) ) {
                    index_of_struct[i*num_blocks_m*num_blocks_p + j*num_blocks_p + k] = blocks_in;
                    blocks_in++;
                }
            }
        }
    }

    printf("%i blocks are internal out of %i blocks in total\n", blocks_in, total_blocks);

}

__global__ void perform_stencil_internal(struct block *u, real l2, real l, real g) {
    //Do stencil.
    //internal elements
    //Block size should be dimensions of struct blocks -2 in each direction.

    //The following statement effectively caches the block.
    //FIX THIS using SHARED.???
    struct block bl = u[blockIdx.x];
    
    int x = threadIdx.x + 1;
    int y = threadIdx.y + 1;
    int z = threadIdx.z + 1;
    int k = bl.k[x][y][z];

    bl.u[x][y][z] = ((2 - l2 * k) * bl.u1[x][y][z] +
                    l2 * ( bl.u1[x][y][z+1] + 
                           bl.u1[x][y][z-1] + 
                           bl.u1[x][y+1][z] + 
                           bl.u1[x][y-1][z] + 
                           bl.u1[x+1][y][z] +
                           bl.u1[x-1][y][z] ) +
                           (0.5 * l * g * (6 - k) - 1) * bl.u[x][y][z])/(1 + 0.5 * l * g * (6 - k));

    if (k == 0) {
        bl.u[x][y][z] = 0;
    }
    
}

__global__ void perform_stencil_surfaces_LR(struct block *u, real l2, real l, real g) {
    //Do stencil.
    //Edge elements left and right.
    struct block bl = u[blockIdx.x];
    struct block bl_r = u[bl.right];
    struct block bl_l = u[bl.left];
    
    real left, right;

    int y = threadIdx.y + 1;
    int z = threadIdx.z + 1;
    
    if ( bl.left != -1 ) {
        left = u[bl.left].u1[Bl-1][y][z];
    } else {
        left = 0;
    }
    
    int k = bl.k[0][y][z];
    
    
    bl.u[0][y][z] = ((2 - l2 * k) * bl.u1[0][y][z] +
                    l2 * ( bl.u1[0][y][z+1] + 
                           bl.u1[0][y][z-1] + 
                           bl.u1[0][y+1][z] + 
                           bl.u1[0][y-1][z] + 
                           bl.u1[0+1][y][z] +
                           left ) +
                           (0.5 * l * g * (6 - k) - 1) * bl.u[0][y][z])/(1 + 0.5 * l * g * (6 - k));
    
   

    if (k == 0) {
        bl.u[0][y][z] = 0;
    }
    
}
