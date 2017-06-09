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
#define L 64
#define M 64
#define P 64

//Sphere radius
#define R 64/2

struct block {
    int up;
    int down;
    int left;
    int right;
    int fore;
    int aft;
    float u[Bl][Bm][Bp];
    float u1[Bl][Bm][Bp];
    char inside[Bl][Bm][Bp];
    int coords[3];
};

__global__ void perform_stencil_internal(struct block *u, struct block *u1, real scale);


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

int main() {

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
    */

    //Assign a block for all volumes containing points
    //by first assigning enough blocks for the whole lot.

    int num_blocks_l = L/Bl;
    int num_blocks_m = M/Bm;
    int num_blocks_p = P/Bp;

    struct block *array = (struct block *) calloc(num_blocks_l*num_blocks_m*num_blocks_p, sizeof(struct block));
    int blocks_in = 0;

    for (i = 0; i < num_blocks_l; i++) {
        for (j = 0; j < num_blocks_m; j++) {
            for (k = 0; k < num_blocks_p; k++) {
            }
        }
    }

}


__global__ void perform_stencil_internal(struct block *u, struct block *u1, real scale) {
    //Do stencil.
    //internal elements
    struct block bl = u[blockIdx.x];
    
    int x = threadIdx.x + 1;
    int y = threadIdx.y + 1;
    int z = threadIdx.z + 1;
    
    bl.u[x][y][z] = scale*( bl.u1[x][y][z+1] + 
                            bl.u1[x][y][z-1] + 
                            bl.u1[x][y+1][z] + 
                            bl.u1[x][y-1][z] + 
                            bl.u1[x+1][y][z] +
                            bl.u1[x-1][y][z]) - bl.u[x][y][z];
    
}
